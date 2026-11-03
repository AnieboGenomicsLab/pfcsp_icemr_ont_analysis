#!/usr/bin/env python3
"""
Extract the PfCSP reference set from the Pf3D7 genome and GFF.

Writes into the reference directory:
  pfcsp_gene.fasta       gene region, genome strand
  pfcsp_cds.fasta        spliced CDS, gene strand, start to stop
  pfcsp_protein.faa      translated CDS
  pfcsp_amplicon.fasta   amplicon slice, genome strand, used when MAP_TARGET=amplicon
  pfcsp_regions.bed      amplicon, callable window, CDS, TH2R, TH3R in genome coordinates
  pfcsp_codon_map.tsv    one row per codon: aa position, 3D7 amino acid, three genome positions
  pfcsp_reference_summary.txt

The codon map is what links read level base calls back to amino acid positions.
CSP is on the minus strand, so codon position 1 of the protein sits at the highest
genome coordinate of the CDS.

Author: Khadim Gueye
"""

import argparse
import sys
from pathlib import Path

import pysam

CODON_TABLE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3))


def read_cds(gff_path, transcript):
    """Return (chrom, strand, [(start, end), ...]) for the transcript CDS, ascending."""
    chrom = strand = None
    exons = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            if f"Parent={transcript}" not in f[8]:
                continue
            chrom, strand = f[0], f[6]
            exons.append((int(f[3]), int(f[4])))
    if not exons:
        sys.exit(f"[ERROR] no CDS found for {transcript} in {gff_path}")
    return chrom, strand, sorted(exons)


def write_fasta(path, name, seq, width=60):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description="Extract the PfCSP reference set.")
    ap.add_argument("--genome", required=True, type=Path)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--transcript", default="PF3D7_0304600.1")
    ap.add_argument("--gene-id", default="PF3D7_0304600")
    ap.add_argument("--amplicon-chrom", required=True)
    ap.add_argument("--amplicon-start", required=True, type=int)
    ap.add_argument("--amplicon-end", required=True, type=int)
    ap.add_argument("--primer-trim", type=int, default=25)
    ap.add_argument("--th2r", nargs=2, type=int, default=[311, 327], metavar=("START", "END"))
    ap.add_argument("--th3r", nargs=2, type=int, default=[352, 367], metavar=("START", "END"))
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fa = pysam.FastaFile(str(args.genome))

    chrom, strand, exons = read_cds(args.gff, args.transcript)
    if chrom != args.amplicon_chrom:
        sys.exit(f"[ERROR] CDS is on {chrom} but amplicon is on {args.amplicon_chrom}")

    gene_start, gene_end = exons[0][0], exons[-1][1]
    gene_seq = fa.fetch(chrom, gene_start - 1, gene_end).upper()
    write_fasta(args.out_dir / "pfcsp_gene.fasta",
                f"{args.gene_id}_{chrom}:{gene_start}-{gene_end}_strand{strand}", gene_seq)

    cds_positions = [p for s, e in exons for p in range(s, e + 1)]
    if strand == "-":
        cds_positions = cds_positions[::-1]
    cds_seq = "".join(fa.fetch(chrom, p - 1, p) for p in cds_positions).upper()
    if strand == "-":
        cds_seq = revcomp("".join(fa.fetch(chrom, p - 1, p) for p in cds_positions[::-1]).upper())
    write_fasta(args.out_dir / "pfcsp_cds.fasta", args.transcript, cds_seq)

    protein = translate(cds_seq)
    aa_seq = protein.rstrip("*")
    write_fasta(args.out_dir / "pfcsp_protein.faa", args.transcript, aa_seq)

    amp_seq = fa.fetch(chrom, args.amplicon_start - 1, args.amplicon_end).upper()
    amp_name = f"csp_{chrom}_{args.amplicon_start}_{args.amplicon_end}"
    write_fasta(args.out_dir / "pfcsp_amplicon.fasta", amp_name, amp_seq)
    pysam.faidx(str(args.out_dir / "pfcsp_amplicon.fasta"))

    codon_rows = []
    for n in range(1, len(aa_seq) + 1):
        ps = cds_positions[3 * (n - 1):3 * n]
        if len(ps) < 3:
            break
        codon = "".join(fa.fetch(chrom, p - 1, p) for p in ps).upper()
        if strand == "-":
            codon = codon.translate(COMPLEMENT)
        codon_rows.append((n, aa_seq[n - 1], codon, chrom, strand, ps[0], ps[1], ps[2]))

    with open(args.out_dir / "pfcsp_codon_map.tsv", "w") as fh:
        fh.write("aa_position\taa_3d7\tcodon_3d7\tchrom\tstrand\tpos1\tpos2\tpos3\n")
        for r in codon_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    def aa_window(start, end):
        ps = [p for n, _, _, _, _, a, b, c in codon_rows if start <= n <= end for p in (a, b, c)]
        return (min(ps), max(ps)) if ps else (None, None)

    th2r_lo, th2r_hi = aa_window(*args.th2r)
    th3r_lo, th3r_hi = aa_window(*args.th3r)
    callable_start = args.amplicon_start + args.primer_trim
    callable_end = args.amplicon_end - args.primer_trim

    cds_in_amp = [n for n, _, _, _, _, a, b, c in codon_rows
                  if callable_start <= min(a, b, c) and max(a, b, c) <= callable_end]
    cterm_lo, cterm_hi = (min(cds_in_amp), max(cds_in_amp)) if cds_in_amp else (None, None)

    bed = [
        (chrom, args.amplicon_start - 1, args.amplicon_end, "amplicon"),
        (chrom, callable_start - 1, callable_end, "callable"),
        (chrom, gene_start - 1, gene_end, f"CDS_{args.gene_id}"),
        (chrom, th2r_lo - 1, th2r_hi, f"TH2R_aa{args.th2r[0]}-{args.th2r[1]}"),
        (chrom, th3r_lo - 1, th3r_hi, f"TH3R_aa{args.th3r[0]}-{args.th3r[1]}"),
    ]
    with open(args.out_dir / "pfcsp_regions.bed", "w") as fh:
        for c, s, e, name in sorted(bed, key=lambda x: x[1]):
            fh.write(f"{c}\t{s}\t{e}\t{name}\n")

    lines = [
        f"transcript           {args.transcript}",
        f"gene                 {args.gene_id}",
        f"location             {chrom}:{gene_start}-{gene_end} strand {strand}",
        f"CDS length           {len(cds_seq)} bp",
        f"protein length       {len(aa_seq)} aa",
        f"amplicon             {chrom}:{args.amplicon_start}-{args.amplicon_end} ({len(amp_seq)} bp)",
        f"callable window      {chrom}:{callable_start}-{callable_end} "
        f"(primer trim {args.primer_trim} bp per edge)",
        f"codons fully callable aa {cterm_lo} to {cterm_hi}",
        f"TH2R                 aa {args.th2r[0]}-{args.th2r[1]} = "
        f"{aa_seq[args.th2r[0] - 1:args.th2r[1]]}  ({chrom}:{th2r_lo}-{th2r_hi})",
        f"TH3R                 aa {args.th3r[0]}-{args.th3r[1]} = "
        f"{aa_seq[args.th3r[0] - 1:args.th3r[1]]}  ({chrom}:{th3r_lo}-{th3r_hi})",
    ]
    (args.out_dir / "pfcsp_reference_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    if protein.count("*") != 1 or not protein.endswith("*"):
        print(f"[WARN] CDS does not translate to a single terminal stop: {protein.count('*')} stops")
    if cterm_lo is None:
        sys.exit("[ERROR] no codon falls entirely inside the callable window")
    return 0


if __name__ == "__main__":
    sys.exit(main())
