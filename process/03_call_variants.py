#!/usr/bin/env python3
"""
Call variants across the CSP callable window.

freebayes runs in pooled-continuous mode, which reports allele frequency without
assuming a ploidy. A malaria infection can carry several clones, so a fixed ploidy
model would collapse the frequency information the clone analysis depends on.

freebayes output is treated as a candidate list only. Alleles are then recounted
from the BAM files at every candidate site in every sample, which removes the
difference between "no variant" and "not enough depth", and a plate wide filter
drops the recurrent low level calls typical of ONT systematic error.

Outputs under <out_dir>/05_variants:
  per_sample/<sample>.freebayes.vcf.gz
  candidate_sites_raw.tsv.gz
  sites_filtered.tsv        one row per site, PASS or ARTEFACT, with annotation
  per_sample_calls.tsv      one row per sample and site above threshold
  allele_depth_matrix.tsv   A, C, G, T counts per sample and site

Author: Khadim Gueye
"""

import argparse
import gzip
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

CODON_TABLE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def call_one(job):
    bam, cfg = job
    s = bam.name.replace(".ontarget.bam", "")
    out_vcf = Path(cfg["out_dir"]) / "05_variants" / "per_sample" / f"{s}.freebayes.vcf.gz"
    log = Path(cfg["out_dir"]) / "logs" / f"{s}.freebayes.log"
    cmd = (
        f"freebayes -f {cfg['reference']} -r {cfg['region']} --pooled-continuous "
        f"--haplotype-length 0 --min-alternate-fraction {cfg['min_af']} "
        f"--min-alternate-count {cfg['min_ac']} --min-base-quality {cfg['min_bq']} "
        f"--min-mapping-quality {cfg['min_mapq']} --limit-coverage 5000 {bam} "
        f"| bcftools norm -f {cfg['reference']} -m -any "
        f"| bcftools view -v snps -Oz -o {out_vcf} && bcftools index -t {out_vcf}"
    )
    with open(log, "wb") as fh:
        p = subprocess.run(["bash", "-c", f"set -o pipefail; {cmd}"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        fh.write(p.stdout or b"")
    return s, p.returncode


def load_candidates(out_dir, samples):
    rows = []
    for s in samples:
        vcf = out_dir / "05_variants" / "per_sample" / f"{s}.freebayes.vcf.gz"
        if not vcf.exists():
            continue
        q = subprocess.run(
            ["bcftools", "query", "-f",
             "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AO]\t[%DP]\t%SAF\t%SAR\n", str(vcf)],
            capture_output=True, text=True).stdout
        for line in q.splitlines():
            c, pos, ref, alt, qual, ao, dp, saf, sar = line.split("\t")
            rows.append((s, c, int(pos), ref, alt, float(qual),
                         int(ao), int(dp), int(saf), int(sar)))
    return pd.DataFrame(rows, columns=["sample", "chrom", "pos", "ref", "alt",
                                       "qual", "ao", "dp", "saf", "sar"])


def count_alleles(job):
    bam_path, sites, min_bq = job
    s = Path(bam_path).name.replace(".ontarget.bam", "")
    bam = pysam.AlignmentFile(bam_path)
    out = []
    for chrom, pos in sites:
        cov = bam.count_coverage(chrom, pos - 1, pos, quality_threshold=min_bq,
                                 read_callback="all")
        out.append((s, chrom, pos, *[int(x[0]) for x in cov]))
    bam.close()
    return out


def annotate(codon_map, fa, chrom, pos, alt):
    """Return gene level consequence of one base substitution using the codon map."""
    hit = codon_map.get(pos)
    if hit is None:
        return "", "", "non-coding"
    n, aa_ref, positions, strand = hit
    codon = [fa.fetch(chrom, p - 1, p).upper() for p in positions]
    idx = positions.index(pos)
    mutated = codon.copy()
    mutated[idx] = alt
    ref_codon, alt_codon = "".join(codon), "".join(mutated)
    if strand == "-":
        ref_codon = ref_codon.translate(COMPLEMENT)
        alt_codon = alt_codon.translate(COMPLEMENT)
    aa_alt = CODON_TABLE.get(alt_codon, "X")
    effect = "synonymous" if aa_ref == aa_alt else "missense"
    if aa_alt == "*":
        effect = "stop_gained"
    return n, f"{aa_ref}{n}{aa_alt}", effect


def main():
    ap = argparse.ArgumentParser(description="Call CSP variants in pooled-continuous mode.")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--codon-map", required=True, type=Path)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--callable-start", required=True, type=int)
    ap.add_argument("--callable-end", required=True, type=int)
    ap.add_argument("--min-alt-fraction", type=float, default=0.03)
    ap.add_argument("--min-alt-count", type=int, default=5)
    ap.add_argument("--min-base-qual", type=int, default=10)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--min-allele-freq", type=float, default=0.10)
    ap.add_argument("--min-allele-reads", type=int, default=10)
    ap.add_argument("--min-variant-qual", type=float, default=0.0)
    ap.add_argument("--strand-filter", choices=["auto", "on", "off"], default="auto")
    ap.add_argument("--ntc-pattern", default="NTC")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    var_dir = args.out_dir / "05_variants"
    (var_dir / "per_sample").mkdir(parents=True, exist_ok=True)
    (args.out_dir / "logs").mkdir(parents=True, exist_ok=True)

    bams = sorted((args.out_dir / "03_alignment").glob("*.ontarget.bam"))
    if not bams:
        sys.exit(f"[ERROR] no ontarget bam in {args.out_dir / '03_alignment'}")
    samples = [b.name.replace(".ontarget.bam", "") for b in bams]
    region = f"{args.chrom}:{args.callable_start}-{args.callable_end}"

    cfg = {"out_dir": str(args.out_dir), "reference": str(args.reference), "region": region,
           "min_af": args.min_alt_fraction, "min_ac": args.min_alt_count,
           "min_bq": args.min_base_qual, "min_mapq": args.min_mapq}

    print(f"[INFO] freebayes pooled-continuous on {len(bams)} samples, region {region}")
    with ProcessPoolExecutor(args.jobs) as ex:
        results = list(ex.map(call_one, [(b, cfg) for b in bams]))
    bad = [s for s, rc in results if rc != 0]
    if bad:
        print(f"[WARN] freebayes failed on {len(bad)} samples: {', '.join(bad[:5])}")

    cand = load_candidates(args.out_dir, samples)
    cand.to_csv(var_dir / "candidate_sites_raw.tsv.gz", sep="\t", index=False)
    if cand.empty:
        sys.exit("[ERROR] freebayes produced no candidate variant")

    cand["af"] = cand.ao / cand.dp
    keep = ((cand.qual >= args.min_variant_qual) & (cand.af >= args.min_allele_freq)
            & (cand.ao >= args.min_allele_reads))

    # A strand balance filter only makes sense when both orientations are sequenced.
    # This library is single orientation, so the filter would reject every site.
    total_saf, total_sar = int(cand.saf.sum()), int(cand.sar.sum())
    minor_strand = min(total_saf, total_sar) / max(1, total_saf + total_sar)
    if args.strand_filter == "on" or (args.strand_filter == "auto" and minor_strand >= 0.05):
        keep &= (np.minimum(cand.saf, cand.sar) / cand.ao) >= 0.1
        print(f"[INFO] strand balance filter applied, minor strand {minor_strand:.1%}")
    else:
        print(f"[INFO] strand balance filter skipped, minor strand {minor_strand:.1%}")
    kept = cand[keep]
    sites = kept[["chrom", "pos"]].drop_duplicates().sort_values("pos")
    site_list = list(sites.itertuples(index=False, name=None))
    if not site_list:
        sys.exit("[ERROR] no candidate site after filtering. Check MIN_ALLELE_FREQ and "
                 "STRAND_FILTER in config.sh")
    print(f"[INFO] candidate sites: {len(site_list)}")

    with ProcessPoolExecutor(args.jobs) as ex:
        counted = [r for batch in ex.map(
            count_alleles, [(str(b), site_list, args.min_base_qual) for b in bams]) for r in batch]
    A = pd.DataFrame(counted, columns=["sample", "chrom", "pos", "A", "C", "G", "T"])
    A["depth"] = A[list("ACGT")].sum(axis=1)
    A.to_csv(var_dir / "allele_depth_matrix.tsv", sep="\t", index=False)

    is_ntc = A["sample"].str.contains(args.ntc_pattern, case=False, regex=True)
    field = A[~is_ntc]

    fa = pysam.FastaFile(str(args.reference))
    cm = pd.read_csv(args.codon_map, sep="\t")
    codon_map = {}
    for r in cm.itertuples(index=False):
        for p in (r.pos1, r.pos2, r.pos3):
            codon_map[p] = (r.aa_position, r.aa_3d7, [r.pos1, r.pos2, r.pos3], r.strand)

    site_rows = []
    for (chrom, pos), grp in field.groupby(["chrom", "pos"]):
        ref_base = fa.fetch(chrom, pos - 1, pos).upper()
        grp = grp[grp.depth >= args.min_depth]
        if grp.empty:
            continue
        alts = sorted(set(kept[(kept.chrom == chrom) & (kept.pos == pos)].alt))
        for alt in alts:
            if alt not in "ACGT" or alt == ref_base:
                continue
            f = grp[alt] / grp.depth
            minor = np.minimum(f, 1 - f)
            background = float(np.percentile(minor, 25))
            n_dominant = int((f >= 0.5).sum())
            floor = max(args.min_allele_freq, 3 * background)
            n_minor = int(((f >= floor) & (f < 0.5)).sum())
            artefact = n_dominant == 0 and (background >= 0.03 or n_minor < 2)
            aa_pos, aa_change, effect = annotate(codon_map, fa, chrom, pos, alt)
            site_rows.append({
                "chrom": chrom, "pos": pos, "ref": ref_base, "alt": alt,
                "aa_position": aa_pos, "aa_change": aa_change, "effect": effect,
                "n_samples_evaluated": len(grp), "background_af_q25": round(background, 4),
                "n_samples_alt_dominant": n_dominant, "n_samples_alt_minor": n_minor,
                "max_af": round(float(f.max()), 4), "site_floor": round(floor, 4),
                "status": "ARTEFACT_recurrent_low_level" if artefact else "PASS",
            })
    if not site_rows:
        sys.exit("[ERROR] no site passed the plate level filter. Check MIN_ALLELE_FREQ, "
                 "MIN_SAMPLE_DEPTH and STRAND_FILTER in config.sh")
    S = pd.DataFrame(site_rows).sort_values("pos")
    S.to_csv(var_dir / "sites_filtered.tsv", sep="\t", index=False)

    ntc_depth = A[is_ntc].groupby(["chrom", "pos"]).depth.max() if is_ntc.any() else None
    passing = S[S.status == "PASS"]
    idx = A.set_index(["sample", "chrom", "pos"])
    calls = []
    for s in samples:
        for r in passing.itertuples(index=False):
            key = (s, r.chrom, r.pos)
            if key not in idx.index:
                continue
            x = idx.loc[key]
            depth = int(x.depth)
            if depth < args.min_depth:
                continue
            floor = r.site_floor
            if ntc_depth is not None and (r.chrom, r.pos) in ntc_depth.index:
                floor = max(floor, 2 * float(ntc_depth[(r.chrom, r.pos)]) / depth)
            alt_n, ref_n = int(x[r.alt]), int(x[r.ref])
            alt_f, ref_f = alt_n / depth, ref_n / depth
            if alt_n < args.min_allele_reads or alt_f < floor:
                continue
            mixed = ref_n >= args.min_allele_reads and ref_f >= floor
            calls.append({
                "sample": s, "chrom": r.chrom, "pos": r.pos, "ref": r.ref, "alt": r.alt,
                "aa_position": r.aa_position, "aa_change": r.aa_change, "effect": r.effect,
                "depth": depth, "ref_reads": ref_n, "alt_reads": alt_n,
                "alt_freq": round(alt_f, 4), "ref_freq": round(ref_f, 4),
                "genotype": "MIXED" if mixed else "ALT", "floor_used": round(floor, 4),
            })
    C = pd.DataFrame(calls)
    C.to_csv(var_dir / "per_sample_calls.tsv", sep="\t", index=False)

    print(f"[INFO] sites: {S.status.value_counts().to_dict()}")
    print(f"[INFO] sample level calls: {len(C)}, of which mixed: "
          f"{int((C.genotype == 'MIXED').sum()) if len(C) else 0}")
    print(f"[INFO] output: {var_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
