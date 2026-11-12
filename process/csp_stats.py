#!/usr/bin/env python3
"""
Reading the pipeline output and the population genetics computed on top of it.

No plotting here and no matplotlib import, so the pipeline can build the plot tables on
a machine with no display. The functions were previously inside the notebook helper,
which meant a figure could only be redrawn by recomputing everything.

Author: Khadim Gueye
"""

import csv
import itertools
import re
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

# A residue or base the pipeline could not resolve.
MISSING_AA = set("X*-?")
VALID_NT = set("ACGT")

# IUPAC codes that mean two bases at once. N means no call, not heterozygous.
HET_CODES = set("RYSWKMBDHV")


# --- reading the pipeline output -------------------------------------------------

def read_fasta(path):
    out, name = OrderedDict(), None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name is not None:
                out[name].append(line.upper())
    return OrderedDict((k, "".join(v)) for k, v in out.items())


def load_metadata(path):
    """sample_id to country and state."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df.columns = [c.strip().lower() for c in df.columns]
    key = next(c for c in ("sample_id", "sample", "id") if c in df.columns)
    return df.rename(columns={key: "sample_id"}).set_index("sample_id")


def sample_of(clone_id):
    return clone_id.rsplit("_H", 1)[0]


def is_control(name, ntc="NTC", positive="3D7|3d7"):
    return bool(re.search(ntc, name) or re.search(positive, name))


def load_clones(out_dir, meta, ntc="NTC", positive="3D7|3d7"):
    """One row per phased clone, with its country and state. Controls are dropped."""
    nt = read_fasta(Path(out_dir) / "07_haplotypes" / "haplotypes_nt.fasta")
    aa = read_fasta(Path(out_dir) / "08_protein" / "haplotypes_protein.faa")
    rows = []
    for clone, seq in nt.items():
        if is_control(clone, ntc, positive):
            continue
        s = sample_of(clone)
        if s not in meta.index:
            continue
        rows.append({
            "clone": clone, "sample": s,
            "country": meta.at[s, "country"], "state": meta.at[s, "state"],
            "nt": seq, "aa": aa.get(clone, ""), "usable": usable(seq),
        })
    return pd.DataFrame(rows)


def count_controls(out_dir, ntc="NTC", positive="3D7|3d7"):
    """Control samples, counted once each. 03_alignment also holds .ontarget copies."""
    stems = {p.name.split(".")[0] for p in Path(out_dir).glob("03_alignment/*.bam")}
    return sum(1 for s in stems if is_control(s, ntc, positive))


def codon_index(codon_map, start, end):
    """Genomic position to amino acid position, inside the window."""
    out = {}
    with open(codon_map) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                trio = [int(row[c]) for c in ("pos1", "pos2", "pos3")]
            except (KeyError, ValueError):
                continue
            for g in trio:
                if start <= g <= end:
                    out[g] = int(row["aa_position"])
    return out


def heterozygous_codons(out_dir, meta, codon_map, start, end):
    """One row per sample and codon carrying an ambiguity code in the consensus.

    Read from the consensus rather than from the per sample calls, because the call
    table also carries the contamination floor and that floor suppresses true calls
    wherever the no template controls are themselves contaminated.
    """
    gpos = codon_index(codon_map, start, end)
    cons = read_fasta(Path(out_dir) / "06_consensus" / "consensus_nt.fasta")
    rows = []
    for name, seq in cons.items():
        if name not in meta.index:
            continue
        hits = {gpos[start + i] for i, ch in enumerate(seq)
                if ch in HET_CODES and (start + i) in gpos}
        rows += [{"sample": name, "country": meta.at[name, "country"],
                  "aa_position": aa} for aa in sorted(hits)]
    return pd.DataFrame(rows, columns=["sample", "country", "aa_position"])


def samples_with_consensus(out_dir, meta):
    """Samples that actually have a consensus, the honest denominator for a fraction."""
    cons = read_fasta(Path(out_dir) / "06_consensus" / "consensus_nt.fasta")
    return [n for n in cons if n in meta.index]


def load_moi(out_dir, meta):
    df = pd.read_csv(Path(out_dir) / "07_haplotypes" / "moi_per_sample.tsv", sep="\t")
    df["MOI"] = df["MOI"].astype(float).round().astype(int)
    df = df[df["sample"].isin(meta.index)].copy()
    df["country"] = df["sample"].map(meta["country"])
    df["state"] = df["sample"].map(meta["state"])
    return df


# --- population genetics ---------------------------------------------------------

def usable(seq):
    """True when every base is resolved. One unresolved clone would otherwise blank the
    columns it is ambiguous at for the whole pool, and those are the informative ones."""
    return bool(seq) and set(seq) <= VALID_NT


def _matrix(seqs):
    """Sequences as a character array, dropping any that is not fully resolved."""
    kept = [s for s in seqs if usable(s)]
    if not kept:
        raise ValueError("no fully resolved sequence in this pool")
    return np.array([list(s) for s in kept])


def nucleotide_diversity(seqs):
    """Mean pairwise difference per site, pi."""
    arr = _matrix(seqs)
    n, L = arr.shape
    if n < 2 or L == 0:
        return float("nan")
    total = 0
    for i, j in itertools.combinations(range(n), 2):
        total += np.count_nonzero(arr[i] != arr[j])
    return total / (n * (n - 1) / 2) / L


def haplotype_diversity(seqs):
    """1 minus the summed squared frequencies, corrected for sample size."""
    seqs = [s for s in seqs if usable(s)]
    n = len(seqs)
    if n < 2:
        return float("nan")
    freqs = np.array(list(Counter(seqs).values())) / n
    return (1 - np.sum(freqs ** 2)) * n / (n - 1)


def _mean_diff(a, b):
    return float(np.mean([np.count_nonzero(x != y) for x in a for y in b]))


def hudson_fst(seqs_a, seqs_b):
    """Hudson FST from mean pairwise differences, 0 when the two pools are one."""
    both = _matrix(list(seqs_a) + list(seqs_b))
    a, b = both[:len(seqs_a)], both[len(seqs_a):]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    within_a = sum(np.count_nonzero(a[i] != a[j])
                   for i, j in itertools.combinations(range(len(a)), 2))
    within_a /= len(a) * (len(a) - 1) / 2
    within_b = sum(np.count_nonzero(b[i] != b[j])
                   for i, j in itertools.combinations(range(len(b)), 2))
    within_b /= len(b) * (len(b) - 1) / 2
    between = _mean_diff(a, b)
    if between == 0:
        return 0.0
    return max(0.0, (between - (within_a + within_b) / 2) / between)


def pcoa(seqs, labels, n_axes=2):
    """Classical multidimensional scaling on the pairwise Hamming distance."""
    pairs = [(s, l) for s, l in zip(seqs, labels) if usable(s)]
    seqs, labels = [p[0] for p in pairs], [p[1] for p in pairs]
    arr = _matrix(seqs)
    n = len(arr)
    d = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d[i, j] = d[j, i] = np.count_nonzero(arr[i] != arr[j])
    d2 = d ** 2
    j_mat = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j_mat @ d2 @ j_mat
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals > 0
    coords = vecs[:, pos][:, :n_axes] * np.sqrt(vals[pos][:n_axes])
    explained = vals[pos][:n_axes] / vals[pos].sum() * 100
    out = pd.DataFrame(coords, columns=[f"PCoA{i + 1}" for i in range(coords.shape[1])])
    out["label"] = list(labels)
    return out, explained


def distance_to_reference(clones, ref_seq):
    """Substitutions from 3D7 per clone, the number behind the PCoA."""
    keep = clones[clones["nt"].map(usable)].copy()
    arr = _matrix(keep["nt"].tolist() + [ref_seq])
    ref = arr[-1]
    keep["nt_differences_vs_3d7"] = [int(np.count_nonzero(r != ref)) for r in arr[:-1]]
    return keep[["clone", "sample", "country", "state", "nt_differences_vs_3d7"]]
