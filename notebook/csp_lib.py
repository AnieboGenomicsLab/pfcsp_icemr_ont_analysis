#!/usr/bin/env python3
"""
Analysis and plotting helpers for the PfCSP notebook.

Everything here reads the pipeline output and adds nothing to it. The pipeline stops at
the protein level, so the population genetics computed below, nucleotide diversity,
haplotype diversity, FST and the PCoA, live here rather than in a pipeline step.

The 3D7 reference is carried through the PCoA as an ordinary point, because the question
the figure answers is how far the field sequences sit from the vaccine sequence.

Author: Khadim Gueye
"""

import csv
import itertools
import re
from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

COUNTRY_COLOR = {"Nigeria": "#2E6DB4", "Cameroon": "#D95F43", "Control": "#9E9E9E"}
REF_COLOR = "#111111"
EPITOPE_SHADE = {"TH2R": "#FBE3DA", "TH3R": "#DCE9F5"}
ACCENT = "#1F4E79"
BAR_INK = "#3F3F3F"
GOLD = "#F2C230"

# A residue or base the pipeline could not resolve.
MISSING_AA = set("X*-?")
VALID_NT = set("ACGT")


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
    """Control samples, counted once each. The alignment folder also holds .ontarget copies."""
    stems = {p.name.split(".")[0] for p in Path(out_dir).glob("03_alignment/*.bam")}
    return sum(1 for s in stems if is_control(s, ntc, positive))


# IUPAC codes that mean two bases at once. N means no call, not heterozygous.
HET_CODES = set("RYSWKMBDHV")


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


# --- figures ---------------------------------------------------------------------

def _stems(ax, values, top, min_show, dx, fontsize=6):
    """One stem per substitution. Labels sharing a position sit side by side.

    The labels are rotated, so stacking them vertically would still overprint. They are
    spread along x instead, from a common baseline above the tallest stem of the group.
    """
    by_pos = {}
    for (pos, sub), v in values.items():
        by_pos.setdefault(pos, []).append((v, sub))
    for pos, items in by_pos.items():
        items.sort(reverse=True)
        base = max(v for v, _ in items) + top * 0.04
        shown = [(v, sub) for v, sub in items if v >= min_show]
        for v, sub in items:
            ax.vlines(pos, 0, v, color="#555555", lw=0.9, zorder=2)
            ax.plot([pos], [v], "o", ms=4, color="#B02418", zorder=3)
        for i, (v, sub) in enumerate(shown):
            x = pos + (i - (len(shown) - 1) / 2) * dx
            ax.text(x, base, sub, rotation=90, ha="center", va="bottom",
                    fontsize=fontsize, fontweight="bold", zorder=4)


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontweight("bold")


def fig_sample_counts(meta, n_controls, title=None):
    counts = meta["country"].value_counts()
    labels = [c for c in ("Nigeria", "Cameroon") if c in counts.index]
    values = [int(counts[c]) for c in labels]
    if n_controls:
        labels, values = labels + ["Control"], values + [n_controls]
    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    bars = ax.bar(labels, values, color=[COUNTRY_COLOR.get(l, BAR_INK) for l in labels],
                  edgecolor="black", linewidth=1.1, width=0.62)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v + max(values) * 0.015, str(v),
                ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of samples", fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(values) * 1.12)
    _style(ax)
    ax.tick_params(labelsize=11)
    for lab in ax.get_xticklabels():
        lab.set_fontweight("bold")
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig


def fig_moi(moi, country_a="Nigeria", country_b="Cameroon"):
    """One violin per state, the second country kept as a single group on the right."""
    groups, colours = [], []
    for state in sorted(moi.loc[moi.country == country_a, "state"].unique()):
        groups.append((state, moi.loc[moi.state == state, "MOI"].values))
        colours.append("#A8C4E0")
    if (moi.country == country_b).any():
        groups.append((country_b, moi.loc[moi.country == country_b, "MOI"].values))
        colours.append("#F2B7AC")

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    parts = ax.violinplot([g[1] for g in groups], showextrema=False, widths=0.85)
    for body, c in zip(parts["bodies"], colours):
        body.set_facecolor(c)
        body.set_edgecolor("black")
        body.set_linewidth(0.9)
        body.set_alpha(1.0)
    for i, (_, vals) in enumerate(groups, start=1):
        ax.scatter([i], [np.mean(vals)], marker="*", s=230, color=GOLD,
                   edgecolor="black", linewidth=0.8, zorder=3)
        ax.scatter([i], [np.median(vals)], marker="o", s=95, color=GOLD,
                   edgecolor="black", linewidth=0.8, zorder=3)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([g[0] for g in groups], rotation=45, ha="right")
    ax.set_xlabel("State", fontsize=11, fontweight="bold")
    ax.set_ylabel("Multiplicity of infection (MOI)", fontsize=11, fontweight="bold")
    ax.legend(handles=[plt.Line2D([], [], marker="*", color="none", markersize=15,
                                  markerfacecolor=GOLD, markeredgecolor="black",
                                  label="Mean"),
                       plt.Line2D([], [], marker="o", color="none", markersize=10,
                                  markerfacecolor=GOLD, markeredgecolor="black",
                                  label="Median")],
              frameon=False, loc="upper left", prop={"weight": "bold", "size": 10})
    _style(ax)
    fig.tight_layout()
    return fig


def fig_diversity(clones, country_a="Nigeria", country_b="Cameroon"):
    """Nucleotide diversity and haplotype diversity, one bar per country."""
    countries = [c for c in (country_a, country_b)
                 if (clones.country == c).any()]
    pis, hds = [], []
    for c in countries:
        seqs = clones.loc[clones.country == c, "nt"].tolist()
        pis.append(nucleotide_diversity(seqs))
        hds.append(haplotype_diversity(seqs))

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6))
    for ax, vals, label, fmt in ((axes[0], pis, "Nucleotide diversity (pi)", "{:.4f}"),
                                 (axes[1], hds, "Haplotype diversity (Hd)", "{:.4f}")):
        bars = ax.bar(countries, vals, color=[COUNTRY_COLOR[c] for c in countries],
                      edgecolor="black", linewidth=1.1, width=0.6)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + max(vals) * 0.02, fmt.format(v),
                    ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.set_ylabel(label, fontsize=11, fontweight="bold")
        ax.set_ylim(0, max(vals) * 1.18)
        _style(ax)
    fig.tight_layout()
    return fig


def fig_substitutions(het, subs, th2r, th3r, n_by_country, window, top_n=15):
    """Three panels. Heterozygosity, the commonest changes per sample, and prevalence.

    het comes from heterozygous_codons, subs from 08_protein/substitutions_consensus.tsv.
    n_by_country gives the denominator per country, window the codon range to span.
    """
    fig = plt.figure(figsize=(13.0, 8.6))
    ax_a = fig.add_axes([0.06, 0.60, 0.38, 0.34])
    ax_b = fig.add_axes([0.54, 0.60, 0.42, 0.34])
    ax_c = fig.add_axes([0.06, 0.08, 0.90, 0.38])

    # Panel A, how often each position is heterozygous within a host.
    span = range(window[0], window[1] + 1)
    for country, colour in (("Nigeria", COUNTRY_COLOR["Nigeria"]),
                            ("Cameroon", COUNTRY_COLOR["Cameroon"])):
        total = n_by_country.get(country, 0)
        if not total:
            continue
        pool = het[het.country == country]
        frac = (pool.groupby("aa_position")["sample"].nunique() / total
                ).reindex(span, fill_value=0.0)
        ax_a.plot(frac.index, frac.values, lw=1.1, color=colour, label=country)
    for name, (lo, hi) in (("TH2R", th2r), ("TH3R", th3r)):
        ax_a.axvspan(lo, hi, color=EPITOPE_SHADE[name], zorder=0)
        ax_a.text((lo + hi) / 2, ax_a.get_ylim()[1], name, ha="center", va="bottom",
                  fontsize=9, fontweight="bold",
                  color="#C0504D" if name == "TH2R" else ACCENT)
    ax_a.set_xlabel("CSP amino acid position", fontsize=10, fontweight="bold")
    ax_a.set_ylabel("Fraction of samples heterozygous", fontsize=10, fontweight="bold")
    ax_a.legend(frameon=False, prop={"weight": "bold", "size": 9})
    _style(ax_a)

    # Panel B, presence or absence of the commonest epitope changes, one column per sample.
    ep = subs[subs.region.isin(("TH2R", "TH3R"))]
    order = [s for s, _ in Counter(ep.substitution).most_common(top_n)]
    samples = sorted(subs["sample"].unique())
    grid = np.zeros((len(order), len(samples)))
    index = {s: i for i, s in enumerate(samples)}
    for sub, grp in ep[ep.substitution.isin(order)].groupby("substitution"):
        for s in grp["sample"].unique():
            grid[order.index(sub), index[s]] = 1
    ax_b.imshow(grid, aspect="auto", cmap="Blues", vmin=0, vmax=1.4,
                interpolation="nearest")
    ax_b.set_yticks(range(len(order)))
    ax_b.set_yticklabels(order, fontsize=8, fontweight="bold")
    ax_b.set_xticks([])
    ax_b.set_xlabel(f"Samples (n = {len(samples)})", fontsize=10, fontweight="bold")
    ax_b.legend(handles=[Patch(facecolor="#08306B", label="Present"),
                         Patch(facecolor="#F7FBFF", edgecolor="black", label="Absent")],
                loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False,
                prop={"weight": "bold", "size": 9})
    for s in ax_b.spines.values():
        s.set_visible(False)

    # Panel C, how many samples carry each change, by position.
    per_pos = subs.groupby(["aa_position", "substitution"])["sample"].nunique()
    for name, (lo, hi) in (("TH2R", th2r), ("TH3R", th3r)):
        ax_c.axvspan(lo, hi, color=EPITOPE_SHADE[name], zorder=0)
    top = per_pos.max()
    _stems(ax_c, per_pos, top, top * 0.05, dx=0.9)
    for name, (lo, hi) in (("TH2R", th2r), ("TH3R", th3r)):
        ax_c.text((lo + hi) / 2, top * 1.06, name, ha="center", va="bottom",
                  fontsize=10, fontweight="bold",
                  color="#C0504D" if name == "TH2R" else ACCENT)
    ax_c.set_xlabel("CSP amino acid position", fontsize=10, fontweight="bold")
    ax_c.set_ylabel("Samples with substitution", fontsize=10, fontweight="bold")
    ax_c.set_ylim(0, top * 1.30)
    _style(ax_c)

    for ax, letter in ((ax_a, "A"), (ax_b, "B"), (ax_c, "C")):
        ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom")
    return fig


def fig_vaccine_match(consensus, clone_status, th2r, th3r, subs_clone):
    """Five panels on how far the cohort sits from the 3D7 antigen."""
    fig = plt.figure(figsize=(12.0, 9.4))
    ax_a = fig.add_axes([0.07, 0.70, 0.36, 0.25])
    ax_b = fig.add_axes([0.58, 0.70, 0.36, 0.25])
    ax_c1 = fig.add_axes([0.07, 0.38, 0.40, 0.22])
    ax_c2 = fig.add_axes([0.55, 0.38, 0.40, 0.22])
    ax_d = fig.add_axes([0.07, 0.06, 0.34, 0.22])
    ax_e = fig.add_axes([0.56, 0.06, 0.32, 0.22])

    ORDER = ["match", "mismatch", "incomplete"]
    FILL = {"match": "#4FA64F", "mismatch": "#B02418", "incomplete": "#9E9E9E"}

    for ax, series, total, label in (
            (ax_a, consensus["vaccine_status"], len(consensus),
             f"Samples (of {len(consensus)})"),
            (ax_b, clone_status["vaccine_status"], len(clone_status),
             f"Clones (of {len(clone_status)})")):
        counts = [int((series == k).sum()) for k in ORDER]
        bars = ax.bar([k.capitalize() for k in ORDER], counts,
                      color=[FILL[k] for k in ORDER], edgecolor="black", linewidth=1.0,
                      width=0.6)
        for b, v in zip(bars, counts):
            ax.text(b.get_x() + b.get_width() / 2, v + max(counts) * 0.02, str(v),
                    ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.set_ylabel(label, fontsize=10, fontweight="bold")
        ax.set_ylim(0, max(counts) * 1.16)
        _style(ax)
    ax_a.set_title("Vaccine match status\n(sample consensus level)", fontsize=10,
                   fontweight="bold")
    ax_b.set_title("Vaccine match status\n(individual clone level)", fontsize=10,
                   fontweight="bold")

    # Panels C, the substitutions responsible, one axis per epitope.
    n_clones = len(clone_status)
    for ax, name, (lo, hi) in ((ax_c1, "TH2R", th2r), (ax_c2, "TH3R", th3r)):
        ax.axvspan(lo, hi, color=EPITOPE_SHADE[name], zorder=0)
        inside = subs_clone[subs_clone.region == name]
        pct = inside.groupby(["aa_position", "substitution"]).size() / n_clones * 100
        _stems(ax, pct, 100.0, 1.0, dx=0.55)
        ax.set_xlim(lo - 2, hi + 3)
        ax.set_ylim(0, 150)
        ax.set_xlabel("CSP amino acid position", fontsize=9, fontweight="bold")
        ax.set_ylabel("Clones (%)", fontsize=9, fontweight="bold")
        ax.set_title(name, fontsize=10, fontweight="bold",
                     color="#C0504D" if name == "TH2R" else ACCENT)
        _style(ax)

    # Panel D, where a mismatching sample mismatches.
    loc = consensus.loc[consensus.vaccine_status == "mismatch", "mismatch_localisation"]
    keys = ["TH2R only", "TH3R only", "TH2R+TH3R"]
    vals = [int((loc == k).sum()) for k in keys]
    bars = ax_d.bar(keys, vals, color=["#B02418", "#7FB3D9", "#B02418"],
                    edgecolor="black", linewidth=1.0, width=0.6)
    for b, v in zip(bars, vals):
        ax_d.text(b.get_x() + b.get_width() / 2, v + max(max(vals), 1) * 0.03, str(v),
                  ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_d.set_ylabel(f"Mismatched samples (of {int(vals and sum(vals))})", fontsize=9,
                    fontweight="bold")
    ax_d.set_ylim(0, max(max(vals), 1) * 1.2)
    ax_d.set_title("Mismatch localisation", fontsize=10, fontweight="bold")
    _style(ax_d)

    # Panel E, the two epitopes scored against each other.
    tab = pd.crosstab(consensus["th2r_status"], consensus["th3r_status"])
    tab = tab.reindex(index=ORDER, columns=ORDER, fill_value=0)
    im = ax_e.imshow(tab.values, cmap="Blues")
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            ax_e.text(j, i, int(tab.values[i, j]), ha="center", va="center", fontsize=9,
                      fontweight="bold",
                      color="white" if tab.values[i, j] > tab.values.max() / 2 else "black")
    ax_e.set_xticks(range(3), [k.capitalize() for k in ORDER], fontweight="bold",
                    fontsize=8, rotation=30, ha="right")
    ax_e.set_yticks(range(3), [k.capitalize() for k in ORDER], fontweight="bold",
                    fontsize=8)
    ax_e.set_xlabel("TH3R status", fontsize=9, fontweight="bold")
    ax_e.set_ylabel("TH2R status", fontsize=9, fontweight="bold")
    ax_e.set_title("TH2R x TH3R match status", fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=ax_e, fraction=0.046, label="Samples")

    for ax, letter in ((ax_a, "A"), (ax_b, "B"), (ax_c1, "C"), (ax_d, "D"), (ax_e, "E")):
        ax.text(-0.14, 1.08, letter, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom")
    return fig


def fig_structure(clones, ref_seq, country_a="Nigeria", country_b="Cameroon",
                  min_per_group=2):
    """FST by state and by country, then a PCoA carrying the 3D7 reference."""
    fig = plt.figure(figsize=(15.0, 4.8))
    ax_a = fig.add_axes([0.05, 0.18, 0.22, 0.70])
    ax_b = fig.add_axes([0.37, 0.18, 0.16, 0.70])
    ax_c = fig.add_axes([0.63, 0.18, 0.24, 0.70])
    for ax in (ax_a, ax_b):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    # Groups are the states of country A plus country B as one group.
    pools = {}
    for state, grp in clones[clones.country == country_a].groupby("state"):
        if len(grp) >= min_per_group:
            pools[state] = grp["nt"].tolist()
    if (clones.country == country_b).any():
        pools[country_b] = clones.loc[clones.country == country_b, "nt"].tolist()
    names = sorted(pools)

    mat = np.full((len(names), len(names)), np.nan)
    for i, j in itertools.combinations(range(len(names)), 2):
        mat[i, j] = mat[j, i] = hudson_fst(pools[names[i]], pools[names[j]])
    masked = np.ma.masked_invalid(np.tril(mat, -1) + np.triu(np.full_like(mat, np.nan)))
    im = ax_a.imshow(masked, cmap="Blues", vmin=0, aspect="auto")
    ax_a.set_xticks(range(len(names)), names, rotation=45, ha="right", fontsize=8,
                    fontweight="bold")
    ax_a.set_yticks(range(len(names)), names, fontsize=8, fontweight="bold")
    ax_a.set_xlabel("State", fontsize=9, fontweight="bold")
    ax_a.set_ylabel("State", fontsize=9, fontweight="bold")
    ax_a.set_title("FST, state level", fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=ax_a, fraction=0.046, label="FST")

    # Country level, a single pair.
    if (clones.country == country_b).any():
        value = hudson_fst(clones.loc[clones.country == country_a, "nt"].tolist(),
                           clones.loc[clones.country == country_b, "nt"].tolist())
        im2 = ax_b.imshow([[np.nan, np.nan], [value, np.nan]], cmap="Blues", vmin=0,
                          vmax=max(value * 10, 0.05), aspect="auto")
        ax_b.text(0, 1, f"{value:.4f}", ha="center", va="center", fontsize=11,
                  fontweight="bold")
        ax_b.set_xticks([0, 1], [country_a, country_b], fontweight="bold", fontsize=9)
        ax_b.set_yticks([0, 1], [country_a, country_b], fontweight="bold", fontsize=9)
        ax_b.set_title("FST, country level", fontsize=10, fontweight="bold")
        fig.colorbar(im2, ax=ax_b, fraction=0.046, label="FST")
    else:
        ax_b.axis("off")

    # PCoA. The reference is one more point, so the distance to the vaccine is readable.
    # Identical clones land on the same spot, so each sequence is plotted once and the
    # marker is sized by how many clones carry it.
    grouped = (clones.groupby(["nt", "country"]).size()
               .reset_index(name="n_clones"))
    seqs = grouped["nt"].tolist() + [ref_seq]
    labels = grouped["country"].tolist() + ["3D7 reference"]
    coords, explained = pcoa(seqs, labels)
    coords["n_clones"] = list(grouped["n_clones"]) + [1]
    # The smaller cohort is drawn last so it is not buried under the larger one.
    present = [c for c in (country_a, country_b) if c in set(labels)]
    present.sort(key=lambda c: -int((clones.country == c).sum()))
    for name in present:
        sub = coords[coords.label == name]
        ax_c.scatter(sub.PCoA1, sub.PCoA2, s=18 + 9 * sub.n_clones,
                     color=COUNTRY_COLOR[name], edgecolor="white", linewidth=0.5,
                     alpha=0.85, label=name, zorder=2)
    ref = coords[coords.label == "3D7 reference"]
    ax_c.scatter(ref.PCoA1, ref.PCoA2, s=260, marker="*", color=REF_COLOR,
                 edgecolor="white", linewidth=1.2, zorder=4, label="3D7 reference")
    ax_c.annotate("3D7", (float(ref.PCoA1.iloc[0]), float(ref.PCoA2.iloc[0])),
                  textcoords="offset points", xytext=(10, 6), fontsize=10,
                  fontweight="bold", color=REF_COLOR, zorder=5)
    ax_c.set_xlabel(f"PCoA1 ({explained[0]:.1f} %)", fontsize=9, fontweight="bold")
    ax_c.set_ylabel(f"PCoA2 ({explained[1]:.1f} %)", fontsize=9, fontweight="bold")
    ax_c.set_title("PCoA with the 3D7 reference", fontsize=10, fontweight="bold")
    ax_c.legend(frameon=False, prop={"weight": "bold", "size": 9},
                loc="center left", bbox_to_anchor=(1.02, 0.5))
    _style(ax_c)

    for ax, letter in ((ax_a, "A"), (ax_b, "B"), (ax_c, "C")):
        ax.text(-0.18, 1.06, letter, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom")
    return fig, coords, explained


def distance_to_reference(clones, ref_seq):
    """Substitutions from 3D7 per clone, the number behind the PCoA."""
    keep = clones[clones["nt"].map(usable)].copy()
    arr = _matrix(keep["nt"].tolist() + [ref_seq])
    ref = arr[-1]
    keep["nt_differences_vs_3d7"] = [int(np.count_nonzero(r != ref)) for r in arr[:-1]]
    return keep[["clone", "sample", "country", "state", "nt_differences_vs_3d7"]]
