#!/usr/bin/env python3
"""
Filter and map the CSP amplicon reads, one sample per input fastq.

Steps per sample: chopper length and quality filter, minimap2 map-ont, sort, index,
then an on target pass that keeps primary reads with MAPQ above the cut off which
overlap the amplicon by at least the required fraction.

With MAP_TARGET=genome the reads are mapped against the whole Pf3D7 genome and
restricted to the amplicon afterwards, so host and off target reads land elsewhere
instead of being forced onto a short CSP reference.

Outputs under <out_dir>:
  02_filtered_reads/<sample>.filtered.fastq.gz
  03_alignment/<sample>.bam, .bam.bai
  03_alignment/<sample>.ontarget.bam, .bam.bai
  01_read_qc/mapping_summary.tsv

Author: Khadim Gueye
"""

import argparse
import json
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pysam


def run(cmd, log_path):
    with open(log_path, "ab") as log:
        p = subprocess.run(cmd, shell=isinstance(cmd, str), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
        log.write(p.stdout or b"")
    return p.returncode


def sample_name(path):
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def process_one(job):
    fq, cfg = job
    s = sample_name(fq)
    out = Path(cfg["out_dir"])
    log = out / "logs" / f"{s}.map.log"
    filt = out / "02_filtered_reads" / f"{s}.filtered.fastq.gz"
    bam = out / "03_alignment" / f"{s}.bam"
    ont = out / "03_alignment" / f"{s}.ontarget.bam"

    chopper = (
        f"zcat -f {fq} | chopper -q {cfg['min_qual']} --minlength {cfg['min_len']} "
        f"--maxlength {cfg['max_len']} --threads 2 | gzip -1 > {filt}"
    )
    if run(["bash", "-c", f"set -o pipefail; {chopper}"], log) != 0:
        return {"sample": s, "status": "chopper_failed"}

    mm2 = (
        f"minimap2 -t {cfg['threads']} -ax map-ont --secondary=no "
        f"-R '@RG\\tID:{s}\\tSM:{s}' {cfg['reference']} {filt} "
        f"| samtools sort -o {bam} - && samtools index {bam}"
    )
    if run(["bash", "-c", f"set -o pipefail; {mm2}"], log) != 0:
        return {"sample": s, "status": "minimap2_failed"}

    chrom, start, end = cfg["chrom"], cfg["start"], cfg["end"]
    span = end - start + 1
    counts = {"reads_mapped": 0, "reads_unmapped": 0, "reads_low_mapq": 0,
              "reads_partial": 0, "reads_on_target": 0}
    lengths = []

    src = pysam.AlignmentFile(str(bam))
    dst = pysam.AlignmentFile(str(ont), "wb", template=src)
    for read in src.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_unmapped:
            counts["reads_unmapped"] += 1
            continue
        counts["reads_mapped"] += 1
        if read.mapping_quality < cfg["min_mapq"]:
            counts["reads_low_mapq"] += 1
            continue
        if read.reference_name != chrom:
            counts["reads_partial"] += 1
            continue
        overlap = min(end, read.reference_end) - max(start - 1, read.reference_start)
        if overlap <= 0 or overlap / span < cfg["min_overlap"]:
            counts["reads_partial"] += 1
            continue
        counts["reads_on_target"] += 1
        lengths.append(read.query_alignment_length)
        dst.write(read)
    dst.close()
    src.close()
    pysam.index(str(ont))

    total = counts["reads_mapped"] + counts["reads_unmapped"]
    rec = {
        "sample": s,
        "reads_after_filter": total,
        **counts,
        "pct_on_target": round(100 * counts["reads_on_target"] / total, 2) if total else 0.0,
        "mean_aligned_length": round(sum(lengths) / len(lengths), 1) if lengths else 0.0,
        "status": "ok",
    }
    (out / "03_alignment" / f"{s}.mapping.json").write_text(json.dumps(rec, indent=2))
    return rec


def main():
    ap = argparse.ArgumentParser(description="Filter and map CSP amplicon reads.")
    ap.add_argument("--input-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", required=True, type=int)
    ap.add_argument("--end", required=True, type=int)
    ap.add_argument("--min-qual", type=int, default=10)
    ap.add_argument("--min-len", type=int, default=250)
    ap.add_argument("--max-len", type=int, default=600)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-overlap", type=float, default=0.8)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    fastqs = sorted(
        p for p in args.input_dir.iterdir()
        if re.search(r"\.(fastq|fq)(\.gz)?$", p.name) and p.is_file()
    )
    if not fastqs:
        sys.exit(f"[ERROR] no fastq found in {args.input_dir}")

    for sub in ("02_filtered_reads", "03_alignment", "01_read_qc", "logs"):
        (args.out_dir / sub).mkdir(parents=True, exist_ok=True)

    cfg = {
        "out_dir": str(args.out_dir), "reference": str(args.reference),
        "chrom": args.chrom, "start": args.start, "end": args.end,
        "min_qual": args.min_qual, "min_len": args.min_len, "max_len": args.max_len,
        "min_mapq": args.min_mapq, "min_overlap": args.min_overlap, "threads": args.threads,
    }

    print(f"[INFO] samples: {len(fastqs)}  reference: {args.reference.name}  jobs: {args.jobs}")
    with ProcessPoolExecutor(args.jobs) as ex:
        records = list(ex.map(process_one, [(f, cfg) for f in fastqs]))

    fields = ["sample", "reads_after_filter", "reads_mapped", "reads_unmapped",
              "reads_low_mapq", "reads_partial", "reads_on_target", "pct_on_target",
              "mean_aligned_length", "status"]
    summary = args.out_dir / "01_read_qc" / "mapping_summary.tsv"
    with open(summary, "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in sorted(records, key=lambda x: x["sample"]):
            fh.write("\t".join(str(r.get(k, "")) for k in fields) + "\n")

    failed = [r["sample"] for r in records if r["status"] != "ok"]
    ok = [r for r in records if r["status"] == "ok"]
    depths = sorted(r["reads_on_target"] for r in ok)
    median = depths[len(depths) // 2] if depths else 0
    print(f"[INFO] mapped {len(ok)}/{len(records)} samples, median on target depth {median}")
    if failed:
        print(f"[WARN] failed: {', '.join(failed)}")
    print(f"[INFO] summary: {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
