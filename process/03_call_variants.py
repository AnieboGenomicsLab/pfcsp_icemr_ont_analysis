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

Indels are read straight from the CIGAR strings and reported in indel_qc.tsv before
any filtering, so excluding them stays an auditable decision rather than a silent one.
Only substitutions are carried into the consensus and protein steps.

Clair3 can be run as an orthogonal check. It is a deep learning caller trained on ONT
R10.4.1 data, so agreement between the two is evidence that the calls do not depend on
the freebayes error model. Its genotype field assumes a ploidy and is not used.

Outputs under <out_dir>/05_variants:
  per_sample/<sample>.freebayes.vcf.gz
  candidate_sites_raw.tsv.gz
  sites_filtered.tsv        one row per site allele, with its status:
                              PASS                          dominant somewhere, or a
                                                            credible minor in 2 samples
                              SINGLETON_unconfirmed         seen in one sample only, kept
                                                            but resting on that sample
                              ARTEFACT_recurrent_low_level  never dominant and present as
                                                            background across the plate
  per_sample_calls.tsv      one row per sample and site above threshold
  allele_depth_matrix.tsv   A, C, G, T counts per sample and site
  indel_qc.tsv              indels seen in the BAM, with homopolymer context
  caller_concordance.tsv    freebayes against Clair3, when enabled

Author: Khadim Gueye
"""

import argparse
import gzip
import os
import subprocess
import sys
from collections import Counter
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
        f"| bcftools view {'' if cfg['keep_indels'] else '-v snps '}-Oz -o {out_vcf} "
        f"&& bcftools index -t {out_vcf}"
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


def scan_indels(job):
    """Count indels per reference position from the CIGAR, independently of any caller.

    An insertion is anchored on the preceding reference base, as in VCF. The frequency
    of an indel is taken over the reads that cover its position, which for a deletion
    is the aligned base count plus the deleted reads.
    """
    bam_path, chrom, lo, hi = job
    sample = Path(bam_path).name.replace(".ontarget.bam", "")
    bam = pysam.AlignmentFile(bam_path)
    counts, n_read = Counter(), 0
    for read in bam.fetch(chrom, lo - 1, hi):
        if read.is_secondary or read.is_supplementary:
            continue
        n_read += 1
        ref_pos = read.reference_start + 1
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                ref_pos += length
            elif op == 2:
                if lo <= ref_pos <= hi:
                    counts[(ref_pos, "DEL", length)] += 1
                ref_pos += length
            elif op == 1:
                if lo <= ref_pos - 1 <= hi:
                    counts[(ref_pos - 1, "INS", length)] += 1
            elif op == 3:
                ref_pos += length
    if n_read < 100:
        bam.close()
        return []

    cov = bam.count_coverage(chrom, lo - 1, hi, quality_threshold=0,
                             read_callback="all")
    bam.close()
    aligned = {lo + i: sum(int(c[i]) for c in cov) for i in range(hi - lo + 1)}
    deleted = Counter()
    for (pos, kind, _), n in counts.items():
        if kind == "DEL":
            deleted[pos] += n

    out = []
    for (pos, kind, ln), n in counts.items():
        depth = aligned.get(pos, 0) + deleted.get(pos, 0)
        if depth == 0:
            continue
        out.append((sample, pos, kind, ln, n, depth, n / depth))
    return out


def homopolymer_run(seq, lo, pos):
    """Length of the homopolymer run the reference base at pos belongs to."""
    i = pos - lo
    if not 0 <= i < len(seq):
        return 0
    base, n = seq[i], 1
    j = i - 1
    while j >= 0 and seq[j] == base:
        n += 1
        j -= 1
    j = i + 1
    while j < len(seq) and seq[j] == base:
        n += 1
        j += 1
    return n


def indel_report(var_dir, bams, reference, chrom, lo, hi, jobs, min_freq, hp_min):
    """Classify every indel seen in the BAM files.

    CANDIDATE_REAL        in frame, or reaching 30 percent in at least one sample.
                          Either could be a real allele, so it is reported as a warning.
    HOMOPOLYMER_ARTEFACT  inside a homopolymer run, out of frame, never clone level.
                          This is the signature of ONT basecalling, not of biology.
    LOW_LEVEL_NOISE       everything else, which is residual alignment error.
    """
    with ProcessPoolExecutor(jobs) as ex:
        rows = [r for batch in ex.map(
            scan_indels, [(str(b), chrom, lo, hi) for b in bams]) for r in batch]
    D = pd.DataFrame(rows, columns=["sample", "pos", "type", "len",
                                    "reads_indel", "depth_at_pos", "af"])
    if D.empty:
        print("[INFO] indel QC: no indel observed in the callable window")
        return
    seq = pysam.FastaFile(str(reference)).fetch(chrom, lo - 1, hi).upper()
    g = (D[D.af >= min_freq].groupby(["pos", "type", "len"])
         .agg(n_samples=("sample", "size"), max_af=("af", "max"),
              median_af_detected=("af", "median")).reset_index())
    if g.empty:
        print(f"[INFO] indel QC: no indel reaches {min_freq:.0%} in any sample")
        return
    g["n_samples_total"] = D["sample"].nunique()
    g["homopolymer"] = [homopolymer_run(seq, lo, p) for p in g.pos]
    g["in_frame"] = g.len % 3 == 0
    g["status"] = np.where(
        g.in_frame | (g.max_af >= 0.30), "CANDIDATE_REAL",
        np.where(g.homopolymer >= hp_min, "HOMOPOLYMER_ARTEFACT", "LOW_LEVEL_NOISE"))
    g = g.sort_values(["status", "n_samples"], ascending=[True, False])
    g.to_csv(var_dir / "indel_qc.tsv", sep="\t", index=False)
    print(f"[INFO] indel QC: {len(g)} indels above {min_freq:.0%}, "
          f"{g.status.value_counts().to_dict()}")
    real = g[g.status == "CANDIDATE_REAL"]
    if len(real):
        print("[WARN] indels that are not explained by homopolymer error:")
        print(real.to_string(index=False))
        print("[WARN] the consensus and protein steps model substitutions only. "
              "Review indel_qc.tsv before using the protein output.")


def clair3_one(job):
    """Run Clair3 on one sample. Used as an orthogonal check, never as the caller."""
    bam, cfg = job
    sample = Path(bam).name.replace(".ontarget.bam", "")
    wd = Path(cfg["out_dir"]) / "05_variants" / "clair3" / sample
    log = Path(cfg["out_dir"]) / "logs" / f"{sample}.clair3.log"
    cmd = (f"run_clair3.sh --bam_fn={bam} --ref_fn={cfg['reference']} "
           f"--model_path={cfg['model']} --platform=ont --threads={cfg['c3_threads']} "
           f"--bed_fn={cfg['bed']} --ctg_name={cfg['chrom']} --include_all_ctgs "
           f"--no_phasing_for_fa --var_pct_full=1.0 --ref_pct_full=1.0 --output={wd}")
    with open(log, "wb") as fh:
        p = subprocess.run(["bash", "-c", cmd],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        fh.write(p.stdout or b"")
    return sample, p.returncode


def clair3_concordance(var_dir, samples, calls):
    """Join the Clair3 allele frequencies onto the freebayes calls."""
    rows = []
    for sample in samples:
        vcf = var_dir / "clair3" / sample / "merge_output.vcf.gz"
        if not vcf.exists():
            continue
        q = subprocess.run(
            ["bcftools", "query", "-i", 'FILTER="PASS"', "-f",
             "%POS\t%REF\t%ALT\t[%GT]\t[%AF]\n", str(vcf)],
            capture_output=True, text=True).stdout
        for line in q.splitlines():
            pos, ref, alt, gt, af = line.split("\t")
            rows.append((sample, int(pos), ref, alt, gt, float(af)))
    C3 = pd.DataFrame(rows, columns=["sample", "pos", "ref", "alt",
                                     "clair3_gt", "clair3_af"])
    if C3.empty or calls.empty:
        print("[WARN] Clair3 produced no comparable call")
        return
    merged = calls[["sample", "pos", "ref", "alt", "alt_freq", "genotype"]].merge(
        C3, on=["sample", "pos", "ref", "alt"], how="outer", indicator=True)
    merged["support"] = merged._merge.map({"left_only": "freebayes_only",
                                           "right_only": "clair3_only", "both": "both"})
    merged = merged.drop(columns="_merge")
    both = merged[merged.support == "both"]
    merged["af_difference"] = (merged.alt_freq - merged.clair3_af).abs().round(4)
    merged.to_csv(var_dir / "caller_concordance.tsv", sep="\t", index=False)
    print(f"[INFO] caller concordance: {merged.support.value_counts().to_dict()}")
    if len(both):
        d = (both.alt_freq - both.clair3_af).abs()
        print(f"[INFO] allele frequency agreement on shared calls: "
              f"median {d.median():.4f}, max {d.max():.4f}")


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
    ap.add_argument("--artefact-background", type=float, default=0.03)
    ap.add_argument("--keep-indels", action="store_true",
                    help="keep indels in the per sample VCF, reported but not consensed")
    ap.add_argument("--indel-min-freq", type=float, default=0.05)
    ap.add_argument("--indel-homopolymer-min", type=int, default=5)
    ap.add_argument("--use-clair3", action="store_true")
    ap.add_argument("--clair3-model", default="")
    ap.add_argument("--clair3-threads", type=int, default=2)
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
           "min_bq": args.min_base_qual, "min_mapq": args.min_mapq,
           "keep_indels": args.keep_indels}

    print(f"[INFO] freebayes pooled-continuous on {len(bams)} samples, region {region}")
    with ProcessPoolExecutor(args.jobs) as ex:
        results = list(ex.map(call_one, [(b, cfg) for b in bams]))
    bad = [s for s, rc in results if rc != 0]
    if bad:
        print(f"[WARN] freebayes failed on {len(bad)} samples: {', '.join(bad[:5])}")

    indel_report(var_dir, bams, args.reference, args.chrom, args.callable_start,
                 args.callable_end, args.jobs, args.indel_min_freq,
                 args.indel_homopolymer_min)

    if args.use_clair3:
        model = args.clair3_model or os.path.join(
            os.environ.get("CONDA_PREFIX", ""), "bin", "models",
            "r1041_e82_400bps_sup_v500")
        bed = var_dir / "callable.bed"
        bed.write_text(f"{args.chrom}\t{args.callable_start - 1}\t{args.callable_end}\n")
        c3cfg = {"out_dir": str(args.out_dir), "reference": str(args.reference),
                 "model": model, "bed": str(bed), "chrom": args.chrom,
                 "c3_threads": args.clair3_threads}
        print(f"[INFO] Clair3 concordance check, model {Path(model).name}")
        with ProcessPoolExecutor(max(1, args.jobs // 2)) as ex:
            c3 = list(ex.map(clair3_one, [(str(b), c3cfg) for b in bams]))
        bad3 = [x for x, rc in c3 if rc != 0]
        if bad3:
            print(f"[WARN] Clair3 failed on {len(bad3)} samples: {', '.join(bad3[:5])}")

    cand = load_candidates(args.out_dir, samples)
    cand.to_csv(var_dir / "candidate_sites_raw.tsv.gz", sep="\t", index=False)
    if cand.empty:
        sys.exit("[ERROR] freebayes produced no candidate variant")

    # Only substitutions reach the consensus and protein steps. Indels are reported by
    # indel_report, which reads the CIGAR strings directly.
    n_raw = len(cand)
    cand = cand[(cand.ref.str.len() == 1) & (cand.alt.str.len() == 1)
                & cand.alt.isin(list("ACGT"))]
    if n_raw != len(cand):
        print(f"[INFO] {n_raw - len(cand)} indel records set aside, see indel_qc.tsv")
    if cand.empty:
        sys.exit("[ERROR] no substitution candidate after removing indels")

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
    print(f"[INFO] candidate positions: {len(site_list)}")

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
            if len(alt) != 1 or alt not in "ACGT" or alt == ref_base:
                continue
            f = grp[alt] / grp.depth
            minor = np.minimum(f, 1 - f)
            background = float(np.percentile(minor, 25))
            n_dominant = int((f >= 0.5).sum())
            floor = max(args.min_allele_freq, 3 * background)
            n_minor = int(((f >= floor) & (f < 0.5)).sum())
            if n_dominant == 0 and background >= args.artefact_background:
                status = "ARTEFACT_recurrent_low_level"
            elif n_dominant == 0 and n_minor < 2:
                status = "SINGLETON_unconfirmed"
            else:
                status = "PASS"
            aa_pos, aa_change, effect = annotate(codon_map, fa, chrom, pos, alt)
            site_rows.append({
                "chrom": chrom, "pos": pos, "ref": ref_base, "alt": alt,
                "aa_position": aa_pos, "aa_change": aa_change, "effect": effect,
                "n_samples_evaluated": len(grp), "background_af_q25": round(background, 4),
                "n_samples_alt_dominant": n_dominant, "n_samples_alt_minor": n_minor,
                "max_af": round(float(f.max()), 4), "site_floor": round(floor, 4),
                "status": status,
            })
    if not site_rows:
        sys.exit("[ERROR] no site passed the plate level filter. Check MIN_ALLELE_FREQ, "
                 "MIN_SAMPLE_DEPTH and STRAND_FILTER in config.sh")
    S = pd.DataFrame(site_rows).sort_values("pos")
    S.to_csv(var_dir / "sites_filtered.tsv", sep="\t", index=False)

    # Contamination floor. What matters is how many reads of the alt allele the worst
    # NTC well carries, not how deep that well is.
    ntc_alt = (A[is_ntc].groupby(["chrom", "pos"])[["A", "C", "G", "T"]].max()
               if is_ntc.any() else None)
    passing = S[S.status != "ARTEFACT_recurrent_low_level"]
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
            if ntc_alt is not None and (r.chrom, r.pos) in ntc_alt.index:
                floor = max(floor, 2 * float(ntc_alt.at[(r.chrom, r.pos), r.alt]) / depth)
            alt_n, ref_n = int(x[r.alt]), int(x[r.ref])
            alt_f, ref_f = alt_n / depth, ref_n / depth
            if alt_n < args.min_allele_reads or alt_f < floor:
                continue
            mixed = ref_n >= args.min_allele_reads and ref_f >= floor
            calls.append({
                "sample": s, "chrom": r.chrom, "pos": r.pos, "ref": r.ref, "alt": r.alt,
                "aa_position": r.aa_position, "aa_change": r.aa_change, "effect": r.effect,
                "site_status": r.status,
                "depth": depth, "ref_reads": ref_n, "alt_reads": alt_n,
                "alt_freq": round(alt_f, 4), "ref_freq": round(ref_f, 4),
                "genotype": "MIXED" if mixed else "ALT", "floor_used": round(floor, 4),
            })
    C = pd.DataFrame(calls)
    C.to_csv(var_dir / "per_sample_calls.tsv", sep="\t", index=False)

    if args.use_clair3:
        clair3_concordance(var_dir, samples, C)

    print(f"[INFO] site alleles: {S.status.value_counts().to_dict()}")
    kept_sites = S[S.status != "ARTEFACT_recurrent_low_level"]
    print(f"[INFO] effect of the kept alleles: "
          f"{kept_sites.effect.value_counts().to_dict()}")
    print(f"[INFO] sample level calls: {len(C)}, of which mixed: "
          f"{int((C.genotype == 'MIXED').sum()) if len(C) else 0}")
    print(f"[INFO] output: {var_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
