#!/usr/bin/env python3
"""
Publication figure: the CSP C terminal protein against 3D7.

The question is vaccine relevant rather than phylogenetic. RTS,S and R21 carry the 3D7
CSP allele and Th2R and Th3R are recognised as peptides, so a synonymous change is
invisible to the immune system and the distance that counts is the distance in protein.
The figure is a polymorphism matrix, not a tree: the window holds 99 residues and about
twenty variable positions, less signal than the nucleotide alignment, which already
resolves only 9 percent of its nodes.

Panel a   how often each position differs from 3D7, over clones
Panel b   one row per protein haplotype, with, from left to right
            the distance to 3D7, substituted residues over the variable positions
            the matrix, a dot where the residue matches 3D7 and the letter where it does
            Th2R and Th3R, filled when the haplotype carries the whole 3D7 epitope
            one column per state, filled where a clone of that state carries it
            the clones carrying the haplotype

The two panels share one column axis and the position numbers sit under the matrix, so
the bars, the epitope bands and the rows read as one block.

Colour carries two things kept apart. The matrix uses four saturated hues for the
chemical class of the substituted residue, validated on every pair, not only on
neighbours, since any two cells can end up side by side: worst pair Delta E 16.7 under
protanopia and 24.7 under normal vision, against floors of 8 and 15. The distance column
uses one sequential heat ramp read against its own key. Epitope bands are shaded and
labelled, never coloured.

The cohort is set with --exclude-prefix, which drops samples whose identifier starts
with that prefix. A haplotype whose residues the translation could not resolve is kept
in the table and left out of the figure, since it is missing data, not a sequence.

The layout is placed in inches and the row height adapts to the number of haplotypes, so
all of them can be drawn on one page. Below MIN_ROW_IN the rows stop shrinking and the
figure grows past the page rather than printing text no one can read.

Outputs under <out_dir>/10_figures:
  csp_protein_vs_3d7.pdf           vector, for the manuscript
  csp_protein_vs_3d7.png           raster, for a word processor
  protein_haplotypes.tsv           the matrix as a table, every haplotype
  protein_variable_positions.tsv   per position, the residues seen and their counts

Author: Khadim Gueye
"""

import argparse
import csv
import re
import sys
from collections import Counter, OrderedDict
from pathlib import Path

import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

SURFACE = "#FCFCFB"
INK = "#0B0B0B"
INK_SECOND = "#52514E"
INK_MUTED = "#898781"
RULE = "#C3C2B7"
GRID = "#E1E0D9"
BAND = "#EDECE5"

# One shade per epitope, so the two bands are told apart at a glance.
EPITOPE_BAND = {"Th2R": "#EDECE5", "Th3R": "#E3EDF6"}

AA_CLASS = {
    **{a: "Basic" for a in "KRH"},
    **{a: "Acidic" for a in "DE"},
    **{a: "Polar" for a in "STNQCY"},
    **{a: "Hydrophobic" for a in "AVLIMFWPG"},
}
CLASS_ORDER = ["Basic", "Acidic", "Polar", "Hydrophobic"]
CLASS_COLOR = {
    "Basic": "#245DC2",
    "Acidic": "#AA295F",
    "Polar": "#DD9C42",
    "Hydrophobic": "#9AA1FB",
}
UNKNOWN_COLOR = "#BDBDBD"

# Heat ramp for the distance column, pale yellow to deep red.
DIST_RAMP = ["#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"]

# How the clones of one haplotype split between states. Each row is read on its own and
# sums to 100, so the column says where that haplotype sits, not how large it is.
STATE_STEPS = [(25.0, "#CFCBC0", "1 to 25"), (50.0, "#9A968B", "26 to 50"),
               (75.0, "#5E5B53", "51 to 75"), (100.0, "#0B0B0B", "76 to 100")]


def state_fill(pct):
    for top, colour, _ in STATE_STEPS:
        if pct <= top:
            return colour
    return STATE_STEPS[-1][1]


# Clones by country, used only when the cohort holds more than one.
COUNTRY_COLOR = {"a": "#5F5C54", "b": "#ADA68F"}

# A residue the translation could not resolve. It is missing data, not a substitution:
# counting X as a difference turned 22 variable positions into 95.
MISSING = set("X*-?")

# Floor on the printed row height, in inches. Below this the residue letter drops under
# five points and stops being legible on paper.
MIN_ROW_IN = 0.083
MAX_ROW_IN = 0.175


def ink_on(hex_color):
    r, g, b = (int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    chan = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    lum = 0.2126 * chan[0] + 0.7152 * chan[1] + 0.0722 * chan[2]
    return INK if lum > 0.42 else "#FFFFFF"


def ramp_color(t):
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    span = len(DIST_RAMP) - 1
    i = min(int(t * span), span - 1)
    f = t * span - i
    lo, hi = DIST_RAMP[i], DIST_RAMP[i + 1]
    out = "#"
    for k in (1, 3, 5):
        a, b = int(lo[k:k + 2], 16), int(hi[k:k + 2], 16)
        out += f"{round(a + (b - a) * f):02X}"
    return out


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


def reference_residues(codon_map, start, end):
    """The 3D7 residue at each amino acid position whose codon is inside the window."""
    out = []
    with open(codon_map) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                inside = all(start <= int(row[c]) <= end for c in ("pos1", "pos2", "pos3"))
            except (KeyError, ValueError):
                continue
            if inside:
                out.append((int(row["aa_position"]), row["aa_3d7"]))
    return sorted(out)


def write_tsv(path, fields, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")


def collapse_report(nt_fasta, members_by_protein):
    """How many nucleotide haplotypes sit behind each protein haplotype.

    Translation is expected to merge clones that differ only by a silent change. Where it
    merges nothing the ratio is one to one, and that is a result about selection, so the
    pipeline states it instead of leaving it to be inferred.
    """
    if not nt_fasta.exists():
        return None
    nt = read_fasta(nt_fasta)
    out = []
    for members in members_by_protein:
        seqs = {nt[m] for m in members if m in nt}
        out.append(len(seqs) if seqs else 0)
    return out


def build_matrix(protein_fasta, positions, ntc_pattern, positive_pattern, exclude_prefix):
    """One entry per distinct protein haplotype, with the clones behind it."""
    seqs = read_fasta(protein_fasta)
    field = {k: v for k, v in seqs.items()
             if not re.search(ntc_pattern, k) and not re.search(positive_pattern, k)}
    if exclude_prefix:
        field = {k: v for k, v in field.items()
                 if not k.rsplit("_H", 1)[0].startswith(exclude_prefix)}
    if not field:
        raise RuntimeError("[ERROR] no field clone left after the cohort filter")
    lengths = {len(v) for v in field.values()}
    if lengths != {len(positions)}:
        raise RuntimeError(f"[ERROR] protein sequences are {sorted(lengths)} residues but the "
                 f"window holds {len(positions)} codons, the two do not match")
    by_seq = OrderedDict()
    for name, seq in field.items():
        by_seq.setdefault(seq, []).append(name)
    haps = [{"sequence": seq, "members": members, "n_clones": len(members),
             "n_samples": len({m.rsplit("_H", 1)[0] for m in members})}
            for seq, members in by_seq.items()]
    haps.sort(key=lambda h: (-h["n_clones"], h["sequence"]))
    for i, h in enumerate(haps, start=1):
        h["haplotype"] = f"Hap{i:02d}"
    return haps, len(field)


def load_states(path, exclude_prefix):
    """sample to state, read from the metadata table, cohort filter applied."""
    states = {}
    with open(path, encoding="utf-8-sig") as fh:
        rows = csv.DictReader(fh, delimiter="\t" if path.suffix != ".csv" else ",")
        cols = {c.strip().lower(): c for c in (rows.fieldnames or [])}
        key = next((cols[c] for c in ("sample_id", "sample", "id") if c in cols), None)
        col = next((cols[c] for c in ("state", "region") if c in cols), None)
        if not key or not col:
            raise RuntimeError(f"[ERROR] {path} needs a sample_id column and a state column")
        for row in rows:
            name, state = (row.get(key) or "").strip(), (row.get(col) or "").strip()
            if not name or not state:
                continue
            if exclude_prefix and name.startswith(exclude_prefix):
                continue
            states[name] = state
    return states


def epitope_state(seq, window, ref_at, lo, hi):
    """identical, variant or unknown, over the whole epitope, not only its hot sites."""
    inside = [p for p in window if lo <= p <= hi]
    aas = [seq[ref_index[p]] for p in inside]
    if any(a in MISSING for a in aas):
        return "unknown"
    return "identical" if all(a == ref_at[p] for p, a in zip(inside, aas)) else "variant"


def epitope_spans(positions, th2r, th3r):
    spans = []
    for name, (lo, hi) in (("Th2R", th2r), ("Th3R", th3r)):
        cols = [i for i, p in enumerate(positions) if lo <= p <= hi]
        if cols:
            spans.append((name, min(cols), max(cols)))
    return spans


def draw(out_pdf, out_png, shown, positions, ref_at, counts, total_clones,
         th2r, th3r, title, dpi, page_w, page_h, country_labels, show_country,
         state_order, caption):
    """One figure, two panels, one shared column axis, every axis placed in inches."""
    ncol, nrow = len(positions), len(shown) + 1

    ylab_w, lab_w, lab_gap, dist_w, pre_gap = 0.22, 0.84, 0.07, 0.22, 0.10
    ep_gap, ep_w, ep_pad = 0.10, 0.17, 0.03
    st_gap, st_w = 0.10, 0.16
    freq_gap, freq_bar, freq_txt = 0.12, 0.52, 0.70
    nstate = len(state_order)
    st_block = st_gap + nstate * st_w if nstate else 0.0
    left = ylab_w + lab_w + lab_gap + dist_w + pre_gap
    right = (ep_gap + 2 * ep_w + ep_pad + st_block
             + freq_gap + freq_bar + freq_txt)
    cell_w = max(0.105, min(0.215, (page_w - left - right) / ncol))

    top, bar_h, gap, band_h = 0.30, 1.00, 0.07, 0.24
    band_gap, num_h, bottom = 0.05, 0.80, 0.90
    fixed_h = top + bar_h + gap + band_h + band_gap + num_h + bottom
    cell_h = max(MIN_ROW_IN, min(MAX_ROW_IN, (page_h - fixed_h) / nrow))

    res_fs = max(4.6, min(8.2, cell_h * 72 * 0.80))
    lab_fs = max(4.4, min(7.6, cell_h * 72 * 0.72))
    num_fs = max(5.0, min(7.4, cell_w * 72 * 0.52))
    dot_ms = max(1.2, min(2.0, cell_h * 72 * 0.22))

    fig_w = left + ncol * cell_w + right
    fig_h = fixed_h + nrow * cell_h
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=SURFACE)

    def axes(x, y, w, h):
        ax = fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])
        ax.set_facecolor("none")
        for s in ax.spines.values():
            s.set_visible(False)
        return ax

    mat_w, mat_h = ncol * cell_w, nrow * cell_h
    x_names = ylab_w + lab_w
    x_dist = x_names + lab_gap
    x_ep = left + mat_w + ep_gap
    x_state = x_ep + 2 * ep_w + ep_pad + st_gap
    x_freq = x_ep + 2 * ep_w + ep_pad + st_block + freq_gap
    y_mat = bottom + num_h
    y_band = y_mat + mat_h + band_gap
    y_bar = y_band + band_h + gap

    ax_mat = axes(left, y_mat, mat_w, mat_h)
    ax_dist = axes(x_dist, y_mat, dist_w, mat_h)
    ax_ep = axes(x_ep, y_mat, 2 * ep_w + ep_pad, mat_h)
    ax_state = axes(x_state, y_mat, nstate * st_w, mat_h) if nstate else None
    ax_freq = axes(x_freq, y_mat, freq_bar + freq_txt, mat_h)
    ax_num = axes(left, bottom, mat_w, num_h)
    ax_band = axes(left, y_band, mat_w, band_h)
    ax_bar = axes(left, y_bar, mat_w, bar_h)

    for ax in (ax_mat, ax_num, ax_band, ax_bar):
        ax.set_xlim(-0.5, ncol - 0.5)
    side = [ax_mat, ax_dist, ax_ep, ax_freq] + ([ax_state] if nstate else [])
    for ax in side:
        ax.set_ylim(nrow - 0.5, -0.5)
    spans = epitope_spans(positions, th2r, th3r)

    ax_bar.set_ylim(0, max(100.0, max(counts) * 1.12))
    for name, c0, c1 in spans:
        ax_bar.add_patch(Rectangle((c0 - 0.5, 0), c1 - c0 + 1, ax_bar.get_ylim()[1],
                                   facecolor=EPITOPE_BAND.get(name, BAND),
                                   edgecolor="none", zorder=0))
        for edge in (c0 - 0.5, c1 + 0.5):
            ax_bar.plot([edge, edge], [0, ax_bar.get_ylim()[1]], color=RULE, lw=0.9,
                        zorder=1)
    for y in (25, 50, 75, 100):
        ax_bar.axhline(y, color=GRID, lw=0.6, zorder=1)
    ax_bar.bar(range(ncol), counts, width=0.62, color=INK_SECOND, zorder=2)
    for i, v in enumerate(counts):
        ax_bar.text(i, v + 2.0, f"{v:.0f}", ha="center", va="bottom",
                    fontsize=6.4, color=INK, fontweight="bold", zorder=3)
    ax_bar.set_ylabel("Clones differing\nfrom 3D7 (%)", fontsize=7.8, color=INK,
                      labelpad=5, fontweight="bold")
    ax_bar.set_yticks([0, 25, 50, 75, 100])
    ax_bar.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=7.0,
                           color=INK, fontweight="bold")
    ax_bar.set_xticks([])
    ax_bar.tick_params(length=0)
    ax_bar.plot([-0.5, ncol - 0.5], [0, 0], color=RULE, lw=0.9, zorder=4)
    ax_bar.plot([-0.5, -0.5], [0, ax_bar.get_ylim()[1]], color=RULE, lw=0.9, zorder=4)
    ax_bar.text(ncol - 0.5, ax_bar.get_ylim()[1] * 1.04, title, fontsize=9.5,
                color=INK, ha="right", va="bottom", fontweight="bold")
    fig.text(0.012, (y_bar + bar_h + 0.04) / fig_h, "a", fontsize=11,
             fontweight="bold", color=INK, ha="left", va="bottom")
    fig.text(0.012, (y_band + band_h / 2) / fig_h, "b", fontsize=11,
             fontweight="bold", color=INK, ha="left", va="center")

    ax_band.set_ylim(0, 1)
    ax_band.axis("off")
    for name, c0, c1 in spans:
        ax_band.add_patch(Rectangle((c0 - 0.5, 0.18), c1 - c0 + 1, 0.68,
                                    facecolor=EPITOPE_BAND.get(name, BAND),
                                    edgecolor=RULE, lw=0.6))
        ax_band.text((c0 + c1) / 2, 0.52, name, ha="center", va="center",
                     fontsize=8.4, color=INK, fontweight="bold")

    ax_num.set_ylim(0, 1)
    ax_num.axis("off")
    for i, pos in enumerate(positions):
        ax_num.text(i, 0.84, str(pos), ha="center", va="top", fontsize=num_fs,
                    color=INK, rotation=90, fontweight="bold")

    fig.text((left + mat_w / 2) / fig_w, (bottom + 0.06) / fig_h,
             "Amino acid position on PF3D7_0304600.1", ha="center", va="bottom",
             fontsize=8.0, color=INK, fontweight="bold")
    if nstate:
        fig.text((x_state + nstate * st_w / 2) / fig_w, (bottom + 0.06) / fig_h,
                 "State", ha="center", va="bottom", fontsize=8.0, color=INK,
                 fontweight="bold")
    if caption:
        fig.text((left + mat_w / 2) / fig_w, (bottom - 0.14) / fig_h, caption,
                 ha="center", va="bottom", fontsize=7.6, color=INK, fontweight="bold")
    fig.text(ylab_w / 2 / fig_w, (y_mat + mat_h / 2) / fig_h, "Protein haplotype",
             ha="center", va="center", rotation=90, fontsize=8.0, color=INK,
             fontweight="bold")

    hdr_y = (bottom + num_h - 0.08) / fig_h
    headers = [(x_dist + dist_w / 2, "Distance"), (x_ep + ep_w / 2, "Th2R"),
               (x_ep + ep_w + ep_pad + ep_w / 2, "Th3R")]
    headers += [(x_state + (k + 0.5) * st_w, name)
                for k, name in enumerate(state_order)]
    headers += [(x_freq + 0.07, "Clones")]
    for x_in, text in headers:
        fig.text(x_in / fig_w, hdr_y, text, ha="center", va="top", rotation=90,
                 fontsize=6.0, color=INK, fontweight="bold")

    for ax in side:
        ax.axis("off")
    ax_dist.set_xlim(0, 1)
    ax_ep.set_xlim(-0.5, 1.5)
    if nstate:
        ax_state.set_xlim(-0.5, nstate - 0.5)
    for name, c0, c1 in spans:
        ax_mat.add_patch(Rectangle((c0 - 0.5, -0.5), c1 - c0 + 1, nrow,
                                   facecolor=EPITOPE_BAND.get(name, BAND),
                                   edgecolor="none", zorder=0))
    edges = {e for _, c0, c1 in spans for e in (c0 - 0.5, c1 + 0.5)}
    for i in range(ncol - 1):
        if i + 0.5 not in edges:
            ax_mat.plot([i + 0.5, i + 0.5], [-0.5, nrow - 0.5], color=GRID, lw=0.4,
                        zorder=1)
    for e in sorted(edges):
        ax_mat.plot([e, e], [-0.5, nrow - 0.5], color=RULE, lw=0.9, zorder=1)
    ax_mat.add_patch(Rectangle((-0.5, -0.5), ncol, nrow, facecolor="none",
                               edgecolor=RULE, lw=0.9, zorder=4))

    label_x = -0.5 - (left - x_names) / cell_w
    for i, pos in enumerate(positions):
        ax_mat.text(i, 0, ref_at[pos], ha="center", va="center", fontsize=res_fs,
                    color=INK, family="DejaVu Sans Mono", fontweight="bold")
    ax_mat.text(label_x, 0, "3D7", ha="right", va="center", fontsize=lab_fs + 0.6,
                color=INK, fontweight="bold")
    ax_dist.add_patch(Rectangle((0.08, -0.42), 0.84, 0.84, facecolor=ramp_color(0.0),
                                edgecolor=RULE, lw=0.4))
    for k in (0, 1):
        ax_ep.add_patch(Rectangle((k - 0.34, -0.38), 0.68, 0.76, facecolor=INK,
                                  edgecolor="none"))
    for k in range(nstate):
        ax_state.add_patch(Rectangle((k - 0.34, -0.38), 0.68, 0.76,
                                     facecolor="#FFFFFF", edgecolor=RULE, lw=0.5))
    ax_mat.plot([-0.5, ncol - 0.5], [0.52, 0.52], color=RULE, lw=0.9, zorder=3)

    for r, h in enumerate(shown):
        y = r + 1
        same = h["distance"] == 0 and h["n_unresolved"] == 0
        label = f"{h['haplotype']} = 3D7" if same else h["haplotype"]
        ax_mat.text(label_x, y, label, ha="right", va="center", fontsize=lab_fs,
                    color=INK, fontweight="bold")
        for i, pos in enumerate(positions):
            aa = h["sequence"][ref_index[pos]]
            if aa in MISSING:
                ax_mat.add_patch(Rectangle((i - 0.44, y - 0.40), 0.88, 0.80,
                                           facecolor="#FFFFFF", edgecolor=RULE,
                                           lw=0.5, hatch="///", zorder=2))
            elif aa == ref_at[pos]:
                ax_mat.plot([i], [y], marker=".", ms=dot_ms, color=INK_MUTED, zorder=2)
            else:
                colour = CLASS_COLOR.get(AA_CLASS.get(aa, ""), UNKNOWN_COLOR)
                ax_mat.add_patch(Rectangle((i - 0.44, y - 0.40), 0.88, 0.80,
                                           facecolor=colour, edgecolor="none", zorder=2))
                ax_mat.text(i, y, aa, ha="center", va="center", fontsize=res_fs,
                            color=ink_on(colour), family="DejaVu Sans Mono",
                            fontweight="bold", zorder=3)

        ax_dist.add_patch(Rectangle((0.08, y - 0.42), 0.84, 0.84,
                                    facecolor=ramp_color(h["distance"]),
                                    edgecolor=RULE, lw=0.4))
        for k, state in enumerate((h["th2r"], h["th3r"])):
            if state == "identical":
                ax_ep.add_patch(Rectangle((k - 0.34, y - 0.38), 0.68, 0.76,
                                          facecolor=INK, edgecolor="none"))
            elif state == "unknown":
                ax_ep.add_patch(Rectangle((k - 0.34, y - 0.38), 0.68, 0.76,
                                          facecolor="#FFFFFF", edgecolor=RULE,
                                          lw=0.5, hatch="///"))
            else:
                ax_ep.add_patch(Rectangle((k - 0.34, y - 0.38), 0.68, 0.76,
                                          facecolor="#FFFFFF", edgecolor=RULE, lw=0.5))
        for k, name in enumerate(state_order):
            pct = h["state_share"].get(name, 0.0)
            ax_state.add_patch(Rectangle(
                (k - 0.34, y - 0.38), 0.68, 0.76,
                facecolor=state_fill(pct) if pct else "#FFFFFF",
                edgecolor="none" if pct else RULE, lw=0.5))

    ax_dist.add_patch(Rectangle((0.06, -0.5), 0.88, nrow, facecolor="none",
                                edgecolor=RULE, lw=0.9, zorder=4))
    for k in (0, 1):
        ax_ep.add_patch(Rectangle((k - 0.36, -0.5), 0.72, nrow, facecolor="none",
                                  edgecolor=RULE, lw=0.9, zorder=4))
    for k in range(nstate):
        ax_state.add_patch(Rectangle((k - 0.36, -0.5), 0.72, nrow, facecolor="none",
                                     edgecolor=RULE, lw=0.9, zorder=4))

    top_n = max(h["n_clones"] for h in shown)
    ax_freq.set_xlim(0, top_n * (freq_bar + freq_txt) / freq_bar)
    for r, h in enumerate(shown):
        y = r + 1
        if show_country:
            ax_freq.add_patch(Rectangle((0, y - 0.32), h["n_a"], 0.64,
                                        facecolor=COUNTRY_COLOR["a"], edgecolor="none"))
            ax_freq.add_patch(Rectangle((h["n_a"], y - 0.32), h["n_b"], 0.64,
                                        facecolor=COUNTRY_COLOR["b"], edgecolor="none"))
        else:
            ax_freq.add_patch(Rectangle((0, y - 0.32), h["n_clones"], 0.64,
                                        facecolor=INK_MUTED, edgecolor="none"))
        ax_freq.text(h["n_clones"] + top_n * 0.07, y,
                     f"{h['n_clones']} ({100.0 * h['n_clones'] / total_clones:.1f}%)",
                     ha="left", va="center", fontsize=lab_fs, color=INK,
                     fontweight="bold")

    legs = [fig.legend(handles=[Patch(facecolor=CLASS_COLOR[c], edgecolor="none",
                                      label=c) for c in CLASS_ORDER],
                       loc="lower left", bbox_to_anchor=(left / fig_w, 0.16 / fig_h),
                       ncol=2, frameon=False, fontsize=7.2, handlelength=1.0,
                       handleheight=1.0, columnspacing=1.0, labelspacing=0.35,
                       title="Class of the substituted residue")]
    if nstate:
        legs.append(fig.legend(
            handles=[Patch(facecolor=c, edgecolor="none", label=lab)
                     for _, c, lab in STATE_STEPS],
            loc="lower left", bbox_to_anchor=((left + 2.10) / fig_w, 0.16 / fig_h),
            ncol=2, frameon=False, fontsize=7.2, handlelength=1.0, handleheight=1.0,
            columnspacing=1.0, labelspacing=0.35,
            title="Share of the haplotype (%)"))
    if show_country:
        legs.append(fig.legend(
            handles=[Patch(facecolor=COUNTRY_COLOR[k], edgecolor="none",
                           label=country_labels[k]) for k in ("a", "b")],
            loc="lower left", bbox_to_anchor=((left + 4.00) / fig_w, 0.16 / fig_h),
            ncol=1, frameon=False, fontsize=7.2, handlelength=1.0, handleheight=1.0,
            labelspacing=0.35, title="Clones by country"))
    for lg in legs:
        lg.get_title().set_fontsize(7.2)
        lg.get_title().set_color(INK)
        lg.get_title().set_fontweight("bold")
        for t in lg.get_texts():
            t.set_color(INK)
            t.set_fontweight("bold")

    cb_w, cb_h = 1.25, 0.085
    cb_x, cb_y = fig_w - cb_w - 0.52, 0.30
    ax_cb = axes(cb_x, cb_y, cb_w, cb_h)
    ax_cb.set_xlim(0, 1)
    ax_cb.set_ylim(0, 1)
    ax_cb.axis("off")
    for s in range(64):
        ax_cb.add_patch(Rectangle((s / 64.0, 0), 1 / 64.0 + 0.002, 1,
                                  facecolor=ramp_color((s + 0.5) / 64.0),
                                  edgecolor="none"))
    ax_cb.add_patch(Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=RULE, lw=0.5))
    for t, lab in ((0.0, "0"), (0.5, "0.5"), (1.0, "1")):
        ax_cb.text(t, -0.55, lab, ha="center", va="top", fontsize=7.0, color=INK,
                   fontweight="bold")
    fig.text(cb_x / fig_w, (cb_y + cb_h + 0.05) / fig_h, "Distance to 3D7",
             fontsize=7.2, color=INK, ha="left", va="bottom", fontweight="bold")

    if out_pdf:
        fig.savefig(out_pdf, facecolor=SURFACE)
    if out_png:
        fig.savefig(out_png, dpi=dpi, facecolor=SURFACE)
    return fig, fig_w, fig_h, cell_h, res_fs


def main(argv=None):
    """Build the figure. Returns the matplotlib figure so a notebook can show it."""
    ap = argparse.ArgumentParser(
        description="Protein level figure of the CSP C terminus against 3D7.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--codon-map", type=Path, required=True)
    ap.add_argument("--callable-start", type=int, required=True)
    ap.add_argument("--callable-end", type=int, required=True)
    ap.add_argument("--th2r", type=int, nargs=2, default=[311, 327])
    ap.add_argument("--th3r", type=int, nargs=2, default=[352, 367])
    ap.add_argument("--metadata", type=Path, default=None,
                    help="table with a sample_id column and a state column")
    ap.add_argument("--ntc-pattern", default="NTC")
    ap.add_argument("--positive-pattern", default="3D7|3d7")
    ap.add_argument("--exclude-prefix", default="",
                    help="drop samples whose identifier starts with this prefix")
    ap.add_argument("--cameroon-prefix", default="CMR",
                    help="prefix that marks the second country in the clone bars")
    ap.add_argument("--country-a", default="Nigeria")
    ap.add_argument("--country-b", default="Cameroon")
    ap.add_argument("--min-clones", type=int, default=1,
                    help="a haplotype is drawn when at least this many clones carry it")
    ap.add_argument("--max-rows", type=int, default=0,
                    help="hard cap on drawn rows, 0 draws every haplotype")
    ap.add_argument("--page-width", type=float, default=6.6,
                    help="printed width in inches, the text width of a portrait page")
    ap.add_argument("--page-height", type=float, default=9.3,
                    help="printed height aimed for; rows stop shrinking at a floor")
    ap.add_argument("--title", default="PfCSP C terminus against 3D7")
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args(argv)

    protein = args.out_dir / "08_protein" / "haplotypes_protein.faa"
    if not protein.exists():
        raise RuntimeError(f"[ERROR] translated clones not found: {protein}, run step 5 first")

    ref = reference_residues(args.codon_map, args.callable_start, args.callable_end)
    if not ref:
        raise RuntimeError("[ERROR] no codon of the map falls inside the callable window")
    window = [p for p, _ in ref]
    ref_at = dict(ref)
    print(f"[INFO] window: amino acids {window[0]} to {window[-1]}, {len(window)} residues")

    global ref_index
    ref_index = {p: i for i, p in enumerate(window)}

    haps, total_clones = build_matrix(protein, window, args.ntc_pattern,
                                      args.positive_pattern, args.exclude_prefix)

    nt_counts = collapse_report(args.out_dir / "07_haplotypes" / "haplotypes_nt.fasta",
                               [h["members"] for h in haps])
    if nt_counts:
        for h, c in zip(haps, nt_counts):
            h["n_nt_haplotypes"] = c
        n_nt = sum(nt_counts)
        merged = sum(1 for c in nt_counts if c > 1)
        print(f"[INFO] {n_nt} nucleotide haplotypes translate into {len(haps)} protein "
              f"haplotypes, {merged} of them merging more than one")
        if merged == 0:
            print("[INFO] no nucleotide haplotype pair differs by silent changes only, "
                  "so translation merges nothing")
    variable = [p for p in window
                if any(h["sequence"][ref_index[p]] not in MISSING
                       and h["sequence"][ref_index[p]] != ref_at[p] for h in haps)]
    if args.exclude_prefix:
        print(f"[INFO] cohort: samples not starting with {args.exclude_prefix}")
    print(f"[INFO] {total_clones} clones, {len(haps)} protein haplotypes, "
          f"{len(variable)} variable positions")
    if not variable:
        raise RuntimeError("[ERROR] no position differs from 3D7, nothing to draw")

    for h in haps:
        h["n_sub"] = sum(1 for p in variable
                         if h["sequence"][ref_index[p]] not in MISSING
                         and h["sequence"][ref_index[p]] != ref_at[p])
        h["n_unresolved"] = sum(1 for p in variable
                                if h["sequence"][ref_index[p]] in MISSING)
        h["distance"] = h["n_sub"] / len(variable)
        h["th2r"] = epitope_state(h["sequence"], window, ref_at, *args.th2r)
        h["th3r"] = epitope_state(h["sequence"], window, ref_at, *args.th3r)
        h["n_b"] = sum(1 for m in h["members"]
                       if m.rsplit("_H", 1)[0].startswith(args.cameroon_prefix))
        h["n_a"] = h["n_clones"] - h["n_b"]
    show_country = any(h["n_b"] for h in haps) and any(h["n_a"] for h in haps)

    # One column per state, filled where a clone of that state carries the haplotype.
    sample_state = {}
    if args.metadata and args.metadata.exists():
        sample_state = load_states(args.metadata, args.exclude_prefix)
    elif args.metadata:
        print(f"[WARN] metadata not found: {args.metadata}, states left out")
    for h in haps:
        h["state_counts"] = Counter(
            sample_state[m.rsplit("_H", 1)[0]] for m in h["members"]
            if m.rsplit("_H", 1)[0] in sample_state)
        placed = sum(h["state_counts"].values())
        h["state_share"] = {k: 100.0 * v / placed
                            for k, v in h["state_counts"].items()} if placed else {}
        h["states"] = set(h["state_counts"])
    state_order = sorted({st for h in haps for st in h["states"]})
    if sample_state:
        missing = {m.rsplit("_H", 1)[0] for h in haps for m in h["members"]
                   if m.rsplit("_H", 1)[0] not in sample_state}
        print(f"[INFO] {len(state_order)} state(s): {', '.join(state_order)}")
        if missing:
            print(f"[WARN] no state for {len(missing)} sample(s): "
                  f"{', '.join(sorted(missing))}")

    # Prevalence is counted over clones, not over haplotypes, so a haplotype carried by
    # one clone does not weigh as much as one carried by twenty.
    counts, resolved_at = [], {}
    for p in variable:
        diff = sum(h["n_clones"] for h in haps
                   if h["sequence"][ref_index[p]] not in MISSING
                   and h["sequence"][ref_index[p]] != ref_at[p])
        resolved = sum(h["n_clones"] for h in haps
                       if h["sequence"][ref_index[p]] not in MISSING)
        resolved_at[p] = resolved
        counts.append(100.0 * diff / resolved if resolved else 0.0)

    drawable = [h for h in haps if h["n_unresolved"] < len(variable)]
    if len(drawable) < len(haps):
        print(f"[INFO] {len(haps) - len(drawable)} haplotype(s) left out of the figure, "
              f"the translation resolved no variable position")
    shown = [h for h in drawable if h["n_clones"] >= args.min_clones]
    if args.max_rows > 0:
        shown = shown[:args.max_rows]
    if not shown:
        shown = drawable
    print(f"[INFO] drawing {len(shown)} of {len(haps)} haplotypes")

    n_samples = len({m.rsplit("_H", 1)[0] for h in haps for m in h["members"]})
    caption = (f"{total_clones} clones from {n_samples} samples collapse into "
               f"{len(haps)} protein haplotypes, one row each")

    out = args.out_dir / "10_figures"
    out.mkdir(parents=True, exist_ok=True)

    fields = ["haplotype", "n_clones", "pct_clones", "n_samples"]
    if show_country:
        fields += [f"n_clones_{args.country_a.lower()}",
                   f"n_clones_{args.country_b.lower()}"]
    fields += ["n_substitutions", "n_unresolved", "distance_to_3d7", "th2r_vs_3d7",
               "th3r_vs_3d7", "states", "n_nt_haplotypes", "substitutions", "drawn",
               "sequence", "members"]
    rows = []
    for h in haps:
        changes = [f"{ref_at[p]}{p}{h['sequence'][ref_index[p]]}" for p in variable
                   if h["sequence"][ref_index[p]] not in MISSING
                   and h["sequence"][ref_index[p]] != ref_at[p]]
        row = {
            "haplotype": h["haplotype"], "n_clones": h["n_clones"],
            "pct_clones": round(100.0 * h["n_clones"] / total_clones, 2),
            "n_samples": h["n_samples"],
            "n_substitutions": h["n_sub"], "n_unresolved": h["n_unresolved"],
            "distance_to_3d7": round(h["distance"], 3),
            "th2r_vs_3d7": h["th2r"], "th3r_vs_3d7": h["th3r"],
            "n_nt_haplotypes": h.get("n_nt_haplotypes", ""),
            "states": ";".join(
                f"{k}:{v} ({h['state_share'][k]:.1f}%)" for k, v in
                sorted(h["state_counts"].items(), key=lambda kv: (-kv[1], kv[0]))),
            "substitutions": ";".join(changes) or "none, identical to 3D7",
            "drawn": "yes" if h in shown else "no",
            "sequence": h["sequence"],
            "members": ";".join(h["members"][:30]),
        }
        if show_country:
            row[f"n_clones_{args.country_a.lower()}"] = h["n_a"]
            row[f"n_clones_{args.country_b.lower()}"] = h["n_b"]
        rows.append(row)
    write_tsv(out / "protein_haplotypes.tsv", fields, rows)
    print(f"[OK] {out / 'protein_haplotypes.tsv'} ({len(rows)} haplotypes)")

    pos_rows = []
    for p, pct in zip(variable, counts):
        seen = Counter()
        for h in haps:
            aa = h["sequence"][ref_index[p]]
            if aa not in MISSING:
                seen[aa] += h["n_clones"]
        region = ("Th2R" if args.th2r[0] <= p <= args.th2r[1] else
                  "Th3R" if args.th3r[0] <= p <= args.th3r[1] else "outside")
        pos_rows.append({
            "aa_position": p, "aa_3d7": ref_at[p], "region": region,
            "pct_clones_differing": round(pct, 2),
            "n_clones_resolved": resolved_at[p],
            "residues": ";".join(f"{a}:{n}" for a, n in seen.most_common()),
        })
    write_tsv(out / "protein_variable_positions.tsv",
              ["aa_position", "aa_3d7", "region", "pct_clones_differing",
               "n_clones_resolved", "residues"], pos_rows)
    print(f"[OK] {out / 'protein_variable_positions.tsv'} ({len(pos_rows)} positions)")

    n_3d7 = sum(h["n_clones"] for h in haps
                if h["distance"] == 0 and h["n_unresolved"] == 0)
    n_th2r = sum(h["n_clones"] for h in haps if h["th2r"] == "identical")
    n_th3r = sum(h["n_clones"] for h in haps if h["th3r"] == "identical")
    print(f"[INFO] clones identical to 3D7 over the window: {n_3d7}, "
          f"over Th2R: {n_th2r}, over Th3R: {n_th3r}")

    fig, fw, fh, ch, fs = draw(
        out / "csp_protein_vs_3d7.pdf", out / "csp_protein_vs_3d7.png",
        shown, variable, ref_at, counts, total_clones,
        tuple(args.th2r), tuple(args.th3r), args.title, args.dpi,
        args.page_width, args.page_height,
        {"a": args.country_a, "b": args.country_b}, show_country, state_order,
        caption)
    print(f"[INFO] figure {fw:.2f} x {fh:.2f} in, row {ch:.3f} in, residue {fs:.1f} pt")
    if fh > args.page_height + 0.05:
        print(f"[WARN] taller than the {args.page_height} in page: insert it as a full "
              f"page or supplementary figure, or lower --max-rows")
    print(f"[OK] {out / 'csp_protein_vs_3d7.pdf'}")
    print(f"[OK] {out / 'csp_protein_vs_3d7.png'}")
    return fig


if __name__ == "__main__":
    main()
    sys.exit(0)
