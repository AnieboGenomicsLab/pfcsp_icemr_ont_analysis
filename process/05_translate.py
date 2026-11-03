#!/usr/bin/env python3
"""
Translate the consensus and the haplotypes, then score them against 3D7.

Codons are read through the codon map, so each amino acid keeps its position on the
397 aa PF3D7_0304600.1 protein. CSP sits on the minus strand and the codon map already
holds the three genome positions per amino acid in transcription order.

An ambiguous nucleotide is expanded to the bases it stands for. When every
combination gives the same amino acid the call is kept, otherwise the amino acid is
recorded as X and the position is not scored against 3D7. That is why a mixed
infection often scores as incomplete at consensus level while its individual
haplotypes score cleanly.

Outputs under <out_dir>/08_protein:
  consensus_protein.faa, haplotypes_protein.faa
  substitutions_consensus.tsv, substitutions_haplotypes.tsv
  substitution_prevalence.tsv
  epitope_status_consensus.tsv, epitope_status_haplotypes.tsv
  epitope_alignment_th2r.tsv, epitope_alignment_th3r.tsv

Author: Khadim Gueye
"""

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

CODON_TABLE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
AMBIGUITY = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def translate_codon(bases, strand):
    """Resolve an ambiguous codon. Returns one amino acid letter, X when degenerate."""
    options = [AMBIGUITY.get(b.upper(), "ACGT") for b in bases]
    if any(len(o) > 1 for o in options) and len(set(map(len, options))) and \
            len(options[0]) * len(options[1]) * len(options[2]) > 64:
        return "X"
    amino = set()
    for combo in itertools.product(*options):
        codon = "".join(combo)
        if strand == "-":
            codon = codon.translate(COMPLEMENT)
        amino.add(CODON_TABLE.get(codon, "X"))
        if len(amino) > 1:
            return "X"
    return amino.pop() if amino else "X"


def translate_window(seq, window_start, codons):
    """Return {aa_position: amino acid} for every codon fully inside the window."""
    out = {}
    for n, aa_ref, positions, strand in codons:
        idx = [p - window_start for p in positions]
        if any(i < 0 or i >= len(seq) for i in idx):
            continue
        bases = [seq[i] for i in idx]
        if any(b == "-" for b in bases):
            out[n] = "X"
            continue
        out[n] = translate_codon(bases, strand)
    return out


def protein_string(calls, lo, hi):
    return "".join(calls.get(n, "X") for n in range(lo, hi + 1))


def score_epitope(calls, ref, start, end):
    """Return match, mismatch or incomplete over one epitope window."""
    subs, unknown = [], 0
    for n in range(start, end + 1):
        aa = calls.get(n, "X")
        r = ref.get(n)
        if aa == "X" or r is None:
            unknown += 1
        elif aa != r:
            subs.append(f"{r}{n}{aa}")
    if subs:
        return "mismatch", subs, unknown
    if unknown:
        return "incomplete", subs, unknown
    return "match", subs, unknown


def main():
    ap = argparse.ArgumentParser(description="Translate CSP consensus and haplotypes.")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--codon-map", required=True, type=Path)
    ap.add_argument("--protein", required=True, type=Path)
    ap.add_argument("--callable-start", required=True, type=int)
    ap.add_argument("--callable-end", required=True, type=int)
    ap.add_argument("--th2r", nargs=2, type=int, default=[311, 327], metavar=("START", "END"))
    ap.add_argument("--th3r", nargs=2, type=int, default=[352, 367], metavar=("START", "END"))
    ap.add_argument("--metadata", type=Path, default=None)
    args = ap.parse_args()

    prot_dir = args.out_dir / "08_protein"
    prot_dir.mkdir(parents=True, exist_ok=True)

    cm = pd.read_csv(args.codon_map, sep="\t")
    lo, hi = args.callable_start, args.callable_end
    codons = [(r.aa_position, r.aa_3d7, [r.pos1, r.pos2, r.pos3], r.strand)
              for r in cm.itertuples(index=False)
              if lo <= min(r.pos1, r.pos2, r.pos3) and max(r.pos1, r.pos2, r.pos3) <= hi]
    if not codons:
        sys.exit("[ERROR] no codon falls inside the callable window")
    aa_lo, aa_hi = min(c[0] for c in codons), max(c[0] for c in codons)
    ref_aa = {c[0]: c[1] for c in codons}
    print(f"[INFO] translating amino acids {aa_lo} to {aa_hi}")

    consensus = read_fasta(args.out_dir / "06_consensus" / "consensus_nt.fasta")
    hap_table = pd.read_csv(args.out_dir / "07_haplotypes" / "haplotypes.tsv", sep="\t")

    meta = None
    if args.metadata and args.metadata.exists():
        sep = "\t" if args.metadata.suffix in (".tsv", ".txt") else ","
        meta = pd.read_csv(args.metadata, sep=sep, dtype=str)
        meta.columns = [c.strip().lower() for c in meta.columns]
        key = next((c for c in meta.columns if c in ("sample_id", "sample", "id")), None)
        meta = meta.set_index(key) if key else None

    def meta_field(sample, field):
        if meta is None or sample not in meta.index or field not in meta.columns:
            return ""
        return meta.loc[sample, field]

    sub_rows, status_rows, prot_records = [], [], []
    for sample, seq in sorted(consensus.items()):
        calls = translate_window(seq, lo, codons)
        prot_records.append((sample, protein_string(calls, aa_lo, aa_hi)))
        for n in range(aa_lo, aa_hi + 1):
            aa = calls.get(n, "X")
            if aa != "X" and aa != ref_aa[n]:
                region = ("TH2R" if args.th2r[0] <= n <= args.th2r[1]
                          else "TH3R" if args.th3r[0] <= n <= args.th3r[1] else "outside")
                sub_rows.append({"sample": sample, "aa_position": n, "aa_3d7": ref_aa[n],
                                 "aa_sample": aa, "substitution": f"{ref_aa[n]}{n}{aa}",
                                 "region": region})
        row = {"sample": sample,
               "country": meta_field(sample, "country"), "state": meta_field(sample, "state"),
               "n_ambiguous_aa": sum(1 for n in range(aa_lo, aa_hi + 1)
                                     if calls.get(n, "X") == "X")}
        for label, (start, end) in (("th2r", tuple(args.th2r)), ("th3r", tuple(args.th3r))):
            st, subs, unknown = score_epitope(calls, ref_aa, start, end)
            row[f"{label}_status"] = st
            row[f"{label}_substitutions"] = ";".join(subs)
            row[f"{label}_ambiguous"] = unknown
        both = {row["th2r_status"], row["th3r_status"]}
        row["vaccine_status"] = ("mismatch" if "mismatch" in both
                                 else "incomplete" if "incomplete" in both else "match")
        row["mismatch_localisation"] = (
            "TH2R+TH3R" if row["th2r_status"] == row["th3r_status"] == "mismatch"
            else "TH2R only" if row["th2r_status"] == "mismatch"
            else "TH3R only" if row["th3r_status"] == "mismatch" else "")
        status_rows.append(row)

    with open(prot_dir / "consensus_protein.faa", "w") as fh:
        for sample, prot in prot_records:
            fh.write(f">{sample} aa{aa_lo}-{aa_hi}\n{prot}\n")

    SC = pd.DataFrame(sub_rows)
    SC.to_csv(prot_dir / "substitutions_consensus.tsv", sep="\t", index=False)
    TC = pd.DataFrame(status_rows)
    TC.to_csv(prot_dir / "epitope_status_consensus.tsv", sep="\t", index=False)

    hap_sub, hap_status, hap_prot = [], [], []
    for r in hap_table.itertuples(index=False):
        calls = translate_window(r.sequence, lo, codons)
        prot = protein_string(calls, aa_lo, aa_hi)
        hap_prot.append((r.haplotype_id, prot, r.reads, r.freq))
        for n in range(aa_lo, aa_hi + 1):
            aa = calls.get(n, "X")
            if aa != "X" and aa != ref_aa[n]:
                region = ("TH2R" if args.th2r[0] <= n <= args.th2r[1]
                          else "TH3R" if args.th3r[0] <= n <= args.th3r[1] else "outside")
                hap_sub.append({"sample": r.sample, "haplotype_id": r.haplotype_id,
                                "aa_position": n, "aa_3d7": ref_aa[n], "aa_sample": aa,
                                "substitution": f"{ref_aa[n]}{n}{aa}", "region": region})
        row = {"sample": r.sample, "haplotype_id": r.haplotype_id,
               "reads": r.reads, "freq": r.freq,
               "country": meta_field(r.sample, "country"),
               "state": meta_field(r.sample, "state")}
        for label, (start, end) in (("th2r", tuple(args.th2r)), ("th3r", tuple(args.th3r))):
            st, subs, unknown = score_epitope(calls, ref_aa, start, end)
            row[f"{label}_status"] = st
            row[f"{label}_substitutions"] = ";".join(subs)
        both = {row["th2r_status"], row["th3r_status"]}
        row["vaccine_status"] = ("mismatch" if "mismatch" in both
                                 else "incomplete" if "incomplete" in both else "match")
        hap_status.append(row)

    with open(prot_dir / "haplotypes_protein.faa", "w") as fh:
        for hid, prot, reads, freq in hap_prot:
            fh.write(f">{hid} aa{aa_lo}-{aa_hi} reads={reads} freq={freq}\n{prot}\n")

    SH = pd.DataFrame(hap_sub)
    SH.to_csv(prot_dir / "substitutions_haplotypes.tsv", sep="\t", index=False)
    TH = pd.DataFrame(hap_status)
    TH.to_csv(prot_dir / "epitope_status_haplotypes.tsv", sep="\t", index=False)

    n_samples = len(consensus)
    if not SC.empty:
        prev = (SC.groupby(["aa_position", "aa_3d7", "aa_sample", "substitution", "region"])
                .agg(n_samples=("sample", "nunique")).reset_index())
        prev["pct_samples"] = (100 * prev.n_samples / n_samples).round(1)
        prev = prev.sort_values("n_samples", ascending=False)
        prev.to_csv(prot_dir / "substitution_prevalence.tsv", sep="\t", index=False)

    for label, (start, end) in (("th2r", tuple(args.th2r)), ("th3r", tuple(args.th3r))):
        rows = [{"sequence_id": "3D7_reference",
                 **{str(n): ref_aa[n] for n in range(start, end + 1)}, "status": "reference"}]
        for sample, prot in prot_records:
            rows.append({
                "sequence_id": sample,
                **{str(n): prot[n - aa_lo] if aa_lo <= n <= aa_hi else "X"
                   for n in range(start, end + 1)},
                "status": TC.set_index("sample").loc[sample, f"{label}_status"],
            })
        pd.DataFrame(rows).to_csv(prot_dir / f"epitope_alignment_{label}.tsv",
                                  sep="\t", index=False)

    print(f"[INFO] consensus level vaccine status: {TC.vaccine_status.value_counts().to_dict()}")
    print(f"[INFO] haplotype level vaccine status: {TH.vaccine_status.value_counts().to_dict()}")
    if not SC.empty:
        top = prev.head(5)[["substitution", "region", "pct_samples"]].to_string(index=False)
        print(f"[INFO] most prevalent substitutions:\n{top}")
    print(f"[INFO] output: {prot_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
