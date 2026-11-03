#!/usr/bin/env python3
"""
Build the per sample consensus and resolve the co-circulating clones.

Two views of the same data:

  consensus    one sequence per sample over the callable window, with an IUPAC
               ambiguity code wherever two alleles pass the frequency floor. A
               mixed infection therefore produces ambiguous positions that cannot
               be scored against 3D7 either way.

  haplotypes   reads carrying the mixed sites are grouped into haplotypes, then
               denoised with the UNOISE3 rule and screened for PCR chimeras. Each
               surviving haplotype is one resolved clone and can be scored on its
               own. MOI is the number of haplotypes that pass.

Alleles the variant step flagged as recurrent ONT error are masked before the mixed
sites are chosen, so a systematic artefact cannot be counted as a clone. An allele seen
in a single sample is not masked: a sample specific allele is the opposite of a
recurrent error, and removing it would discard rare clones. The dominant base of a
position is never masked.

Outputs under <out_dir>:
  06_consensus/consensus_nt.fasta             IUPAC consensus per sample
  06_consensus/consensus_positions.tsv        per position base, depth, ambiguity
  07_haplotypes/haplotypes.tsv                one row per sample and haplotype
  07_haplotypes/haplotypes_nt.fasta           full length haplotype sequences
  07_haplotypes/moi_per_sample.tsv            MOI and clone proportions

Author: Khadim Gueye
"""

import argparse
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import pysam

IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M", frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V", frozenset("ACGT"): "N",
}


def iupac(bases):
    bases = frozenset(bases)
    if len(bases) == 1:
        return next(iter(bases))
    return IUPAC.get(bases, "N")


def denoise(counts, alpha):
    """UNOISE3: fold a haplotype into a more abundant neighbour when it looks like error."""
    centroids = []
    for seq, n in sorted(counts.items(), key=lambda x: -x[1]):
        best = None
        for c in centroids:
            d = sum(1 for a, b in zip(seq, c[0]) if a != b)
            if d > 0 and n / c[1] <= 1 / 2 ** (alpha * d + 1) and (best is None or d < best[0]):
                best = (d, c)
        if best:
            best[1][1] += n
        else:
            centroids.append([seq, n])
    return centroids


def is_bimera(seq, parents):
    """DADA2 style: left part of one parent joined to the right part of another."""
    for a in parents:
        for b in parents:
            if a != b and any(seq[:k] == a[:k] and seq[k:] == b[k:] for k in range(1, len(seq))):
                return True
    return False


def process_one(job):
    sample, cfg = job
    bam_path = Path(cfg["out_dir"]) / "03_alignment" / f"{sample}.ontarget.bam"
    if not bam_path.exists():
        return sample, None, [], []

    chrom, lo, hi = cfg["chrom"], cfg["start"], cfg["end"]
    bam = pysam.AlignmentFile(str(bam_path))

    cons_rows = []
    cons_seq = []
    mixed_sites = []
    for pos in range(lo, hi + 1):
        cov = bam.count_coverage(chrom, pos - 1, pos, quality_threshold=cfg["min_bq"],
                                 read_callback="all")
        counts = {b: int(c[0]) for b, c in zip("ACGT", cov)}
        depth = sum(counts.values())
        if depth < cfg["min_depth"]:
            cons_seq.append("N")
            cons_rows.append((sample, pos, "N", depth, "", 0, "low_depth", ""))
            continue
        present = {b: n for b, n in counts.items()
                   if n >= cfg["min_allele_reads"] and n / depth >= cfg["min_freq"]}
        if not present:
            top = max(counts, key=counts.get)
            present = {top: counts[top]}
        dropped = ""
        if cfg["masked"]:
            major = max(present, key=present.get)
            out = [b for b in present
                   if b != major and (pos, b) in cfg["masked"]]
            if out:
                dropped = ";".join(f"{b}:{present[b]}" for b in out)
                present = {b: n for b, n in present.items() if b not in out}
        base = iupac(present)
        cons_seq.append(base)
        cons_rows.append((sample, pos, base, depth,
                          ";".join(f"{b}:{n}" for b, n in sorted(present.items(),
                                                                 key=lambda x: -x[1])),
                          len(present), "mixed" if len(present) > 1 else "single",
                          dropped))
        if len(present) == 2:
            ordered = sorted(present, key=lambda b: -present[b])
            mixed_sites.append((pos, ordered[0], ordered[1],
                                present[ordered[1]] / depth))

    consensus = "".join(cons_seq)

    hap_rows = []
    n_mixed = len(mixed_sites)
    if mixed_sites:
        # Every read spans the whole amplicon, so any mixed site can be phased. Above the
        # cap, keep the most balanced sites: they separate clones best. Taking them in
        # genomic order would drop one epitope in favour of the other.
        sites = sorted(sorted(mixed_sites, key=lambda x: -x[3])[: cfg["max_sites"]],
                       key=lambda x: x[0])
        strings = []
        for read in bam.fetch(chrom, lo - 1, hi):
            if read.is_secondary or read.is_supplementary:
                continue
            r2q = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
            seq, qual = read.query_sequence, read.query_qualities
            s = ""
            for pos, major, minor, _ in sites:
                q = r2q.get(pos - 1)
                b = seq[q] if (q is not None and qual[q] >= cfg["min_bq_hap"]) else "N"
                s += "0" if b == major else ("1" if b == minor else "N")
            strings.append(s)
        complete = Counter(s for s in strings if "N" not in s)
        n_complete = sum(complete.values())
        if n_complete >= cfg["min_phased_reads"]:
            centroids = sorted(denoise(complete, cfg["alpha"]), key=lambda x: -x[1])
            chimeras = set()
            for i, (s, n) in enumerate(centroids):
                parents = [p for p, pn in centroids[:i] if pn >= n]
                if len(parents) >= 2 and is_bimera(s, parents):
                    chimeras.add(s)
            n_chimeric = sum(n for s, n in centroids if s in chimeras)
            centroids = [c for c in centroids if c[0] not in chimeras]
            for s in strings:
                if "N" in s and s.count("N") <= 0.3 * len(s):
                    match = [c for c in centroids
                             if all(x == "N" or x == y for x, y in zip(s, c[0]))]
                    if len(match) == 1:
                        match[0][1] += 1
            total = sum(n for _, n in centroids)
            keep = [(s, n) for s, n in centroids
                    if n >= cfg["min_allele_reads"] and n / total >= cfg["min_freq"]]
            kept_total = sum(n for _, n in keep)
            site_label = ";".join(f"{p}{a}>{b}" for p, a, b, _ in sites)
            for i, (s, n) in enumerate(sorted(keep, key=lambda x: -x[1]), 1):
                full = list(consensus)
                for (pos, major, minor, _), flag in zip(sites, s):
                    full[pos - lo] = major if flag == "0" else minor
                hap_rows.append({
                    "sample": sample, "haplotype_id": f"{sample}_H{i}",
                    "pattern": s, "sites_used": site_label,
                    "reads": n, "freq": round(n / kept_total, 4),
                    "reads_phased": n_complete, "chimeric_reads_removed": n_chimeric,
                    "mixed_sites_found": n_mixed, "mixed_sites_used": len(sites),
                    "sequence": "".join(full),
                })
    if not hap_rows:
        depth = max((r[3] for r in cons_rows), default=0)
        if depth >= cfg["min_depth"]:
            hap_rows.append({
                "sample": sample, "haplotype_id": f"{sample}_H1",
                "pattern": "", "sites_used": "", "reads": depth, "freq": 1.0,
                "reads_phased": depth, "chimeric_reads_removed": 0,
                "mixed_sites_found": n_mixed, "mixed_sites_used": 0,
                "sequence": consensus,
            })
    bam.close()
    return sample, consensus, cons_rows, hap_rows


def main():
    ap = argparse.ArgumentParser(description="Consensus, haplotypes and MOI for CSP.")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--callable-start", required=True, type=int)
    ap.add_argument("--callable-end", required=True, type=int)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--min-allele-freq", type=float, default=0.10)
    ap.add_argument("--min-allele-reads", type=int, default=10)
    ap.add_argument("--min-phased-reads", type=int, default=50)
    ap.add_argument("--min-base-qual", type=int, default=10)
    ap.add_argument("--denoise-alpha", type=float, default=2.0)
    ap.add_argument("--max-phased-sites", type=int, default=20)
    ap.add_argument("--artefact-sites", type=Path,
                    help="sites_filtered.tsv, alleles not marked PASS are masked")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    bams = sorted((args.out_dir / "03_alignment").glob("*.ontarget.bam"))
    if not bams:
        sys.exit(f"[ERROR] no ontarget bam in {args.out_dir / '03_alignment'}")
    samples = [b.name.replace(".ontarget.bam", "") for b in bams]

    cons_dir = args.out_dir / "06_consensus"
    hap_dir = args.out_dir / "07_haplotypes"
    cons_dir.mkdir(parents=True, exist_ok=True)
    hap_dir.mkdir(parents=True, exist_ok=True)

    masked = set()
    if args.artefact_sites and args.artefact_sites.exists():
        S = pd.read_csv(args.artefact_sites, sep="\t")
        masked = {(int(r.pos), str(r.alt)) for r in S.itertuples(index=False)
                  if str(r.status).startswith("ARTEFACT")}
        print(f"[INFO] masking {len(masked)} artefact alleles from the variant step")

    cfg = {
        "masked": masked,
        "out_dir": str(args.out_dir), "chrom": args.chrom,
        "start": args.callable_start, "end": args.callable_end,
        "min_depth": args.min_depth, "min_freq": args.min_allele_freq,
        "min_allele_reads": args.min_allele_reads, "min_phased_reads": args.min_phased_reads,
        "min_bq": args.min_base_qual, "min_bq_hap": max(7, args.min_base_qual - 3),
        "alpha": args.denoise_alpha, "max_sites": args.max_phased_sites,
    }

    print(f"[INFO] consensus and haplotypes for {len(samples)} samples, "
          f"window {args.chrom}:{args.callable_start}-{args.callable_end}")
    with ProcessPoolExecutor(args.jobs) as ex:
        results = list(ex.map(process_one, [(s, cfg) for s in samples]))

    pos_rows, hap_all = [], []
    with open(cons_dir / "consensus_nt.fasta", "w") as fh:
        for sample, consensus, rows, haps in results:
            pos_rows.extend(rows)
            hap_all.extend(haps)
            if consensus is None:
                continue
            fh.write(f">{sample}\n")
            for i in range(0, len(consensus), 60):
                fh.write(consensus[i:i + 60] + "\n")

    P = pd.DataFrame(pos_rows, columns=["sample", "pos", "base", "depth",
                                        "alleles", "n_alleles", "call", "masked_allele"])
    P.to_csv(cons_dir / "consensus_positions.tsv", sep="\t", index=False)
    hit = P[P.masked_allele != ""]
    if len(hit):
        print(f"[INFO] artefact alleles masked in {hit['sample'].nunique()} samples "
              f"at {hit.pos.nunique()} positions")

    H = pd.DataFrame(hap_all)
    if H.empty:
        sys.exit("[ERROR] no haplotype produced")
    H.to_csv(hap_dir / "haplotypes.tsv", sep="\t", index=False)
    with open(hap_dir / "haplotypes_nt.fasta", "w") as fh:
        for r in H.itertuples(index=False):
            fh.write(f">{r.haplotype_id} reads={r.reads} freq={r.freq}\n")
            for i in range(0, len(r.sequence), 60):
                fh.write(r.sequence[i:i + 60] + "\n")

    amb = P[P.call == "mixed"].groupby("sample").size().rename("ambiguous_positions")
    cov = P[P.call != "low_depth"].groupby("sample").size().rename("positions_covered")
    depth = P.groupby("sample").depth.median().rename("median_depth")
    moi = H.groupby("sample").agg(
        MOI=("haplotype_id", "size"),
        mixed_sites_found=("mixed_sites_found", "max"),
        mixed_sites_used=("mixed_sites_used", "max"),
        reads_phased=("reads_phased", "max"),
        chimeric_reads_removed=("chimeric_reads_removed", "max"),
    )
    moi["clone_proportions"] = H.groupby("sample").freq.apply(
        lambda x: "/".join(f"{v:.2f}" for v in sorted(x, reverse=True)))
    M = moi.join([amb, cov, depth]).fillna({"ambiguous_positions": 0, "positions_covered": 0})
    M["clonality"] = ["polyclonal" if m > 1 else "monoclonal" for m in M.MOI]
    M.to_csv(hap_dir / "moi_per_sample.tsv", sep="\t")
    capped = M[M.mixed_sites_found > M.mixed_sites_used.where(M.mixed_sites_used > 0,
                                                              M.mixed_sites_found)]
    if len(capped):
        print(f"[WARN] {len(capped)} samples carry more than "
              f"{args.max_phased_sites} mixed sites, the most balanced were phased")

    print(f"[INFO] MOI: {M.MOI.value_counts().sort_index().to_dict()}, "
          f"median {M.MOI.median()}")
    print(f"[INFO] clonality: {M.clonality.value_counts().to_dict()}")
    print(f"[INFO] output: {cons_dir} and {hap_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
