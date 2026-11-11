#!/usr/bin/env python3
"""
Publication figure: the CSP C terminal protein against 3D7.

The question this answers is vaccine relevant rather than phylogenetic. RTS,S and R21
carry the 3D7 CSP allele and Th2R and Th3R are T cell epitopes recognised as peptides,
so a synonymous change is invisible to the immune system and the distance that counts
is the distance in protein. This step works on the translated clones for that reason.

It is deliberately not a tree. The window holds 99 residues and 23 variable positions,
which is less signal than the nucleotide alignment already had, and that one resolves
9 percent of its nodes. A polymorphism matrix states the same data without implying an
order of descent nothing supports.

The figure has two panels over one shared set of columns, the variable positions only:

  a   how often each position differs from 3D7 across the clones
  b   one row per protein haplotype, a dot where it matches 3D7 and the residue where
      it does not, with the count behind each row on the right

Colour carries one thing, the chemical class of the substituted residue, because a
charge change inside an epitope is the part a reader is looking for. The four class
colours were validated for colour blind separation on every pair, not only on
neighbours, since any two cells of a matrix can end up side by side: worst pair
Delta E 16.7 under protanopia and 24.7 under normal vision, against floors of 8 and 15.
The epitope bands are shaded and labelled rather than coloured, so the colour channel
stays on the data.

Outputs under <out_dir>/10_figures:
  csp_protein_vs_3d7.pdf           vector, for the manuscript
  csp_protein_vs_3d7.png           raster at 600 dpi, for a word processor
  protein_haplotypes.tsv           the matrix as a table, every haplotype
  protein_variable_positions.tsv   per position, the residues seen and their counts

Author: Khadim Gueye
"""

import argparse
import csv
import re
import sys
import textwrap
from collections import Counter, OrderedDict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

# Ink and surface, so every panel sits on one neutral ground.
SURFACE = "#FCFCFB"
INK = "#0B0B0B"
INK_SECOND = "#52514E"
INK_MUTED = "#898781"
RULE = "#C3C2B7"
GRID = "#E1E0D9"
BAND = "#EDECE5"

# Chemical class of a residue. A substitution that changes the class is the one worth
# seeing first, which is what the colour is spent on.
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

# A residue the translation could not resolve, because the codon carried an ambiguity
# code or too little depth. It is missing data, not a substitution: counting X as a
# difference turned 23 variable positions into 95 and would have put a figure in the
# manuscript claiming four times the diversity the data holds.
MISSING = set("X*-?")


def ink_on(hex_color):
    """Black or white for a label sitting inside a fill, by relative luminance."""
    r, g, b = (int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    chan = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    lum = 0.2126 * chan[0] + 0.7152 * chan[1] + 0.0722 * chan[2]
    return INK if lum > 0.42 else "#FFFFFF"


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


def build_matrix(protein_fasta, positions, ntc_pattern, positive_pattern):
    """One entry per distinct protein haplotype, with the clones behind it."""
    seqs = read_fasta(protein_fasta)
    field = {k: v for k, v in seqs.items()
             if not re.search(ntc_pattern, k) and not re.search(positive_pattern, k)}
    if not field:
        sys.exit("[ERROR] no field clone found in the protein fasta")
    lengths = {len(v) for v in field.values()}
    if lengths != {len(positions)}:
        sys.exit(f"[ERROR] protein sequences are {sorted(lengths)} residues but the "
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


def epitope_spans(positions, th2r, th3r):
    """Column index ranges the two epitopes cover, as drawn."""
    spans = []
    for name, (lo, hi) in (("Th2R", th2r), ("Th3R", th3r)):
        cols = [i for i, p in enumerate(positions) if lo <= p <= hi]
        if cols:
            spans.append((name, min(cols), max(cols)))
    return spans


def draw(out_pdf, out_png, haps, shown, positions, ref_at, counts, total_clones,
         th2r, th3r, title, dpi):
    """One figure, two panels, one shared column axis.

    Every axis is placed in inches rather than left to a layout engine. The matrix has
    a fixed cell size, so the figure has to grow with the number of rows, and an
    automatic layout would rescale the cells instead and let the labels collide.
    """
    ncol, nrow = len(positions), len(shown) + 1          # the reference takes a row
    # Drawn at the size it is printed at.
    #
    # A word processor scales a figure down to the text width of the page. At 11.7
    # inches wide that is a factor of 0.55, so an 8.8 point residue reaches the reader
    # at 4.9 points and nothing is legible. Raising the point sizes does not help,
    # because they are scaled by the same factor. The figure is therefore built 6.6
    # inches wide, which is the text width of a portrait page, and inserted at full
    # size: every point size below is the size the reader sees.
    cell_w, cell_h = 0.198, 0.158
    left, right_gap = 0.88, 1.36
    bottom, top = 1.38, 0.40
    bar_h, band_h, gap = 1.08, 0.24, 0.12

    fig_w = left + ncol * cell_w + right_gap
    fig_h = bottom + nrow * cell_h + gap + band_h + gap + bar_h + top
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=SURFACE)

    def axes(x, y, w, h):
        ax = fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])
        ax.set_facecolor(SURFACE)
        for s in ax.spines.values():
            s.set_visible(False)
        return ax

    mat_h = nrow * cell_h
    ax_mat = axes(left, bottom, ncol * cell_w, mat_h)
    ax_band = axes(left, bottom + mat_h + gap, ncol * cell_w, band_h)
    ax_bar = axes(left, bottom + mat_h + gap + band_h + gap, ncol * cell_w, bar_h)
    ax_freq = axes(left + ncol * cell_w + 0.12, bottom, right_gap - 0.9, mat_h)

    for ax in (ax_mat, ax_band, ax_bar):
        ax.set_xlim(-0.5, ncol - 0.5)
    spans = epitope_spans(positions, th2r, th3r)

    # Panel a, how often each position differs from 3D7.
    ax_bar.set_ylim(0, max(100.0, max(counts) * 1.12))
    for name, c0, c1 in spans:
        ax_bar.add_patch(Rectangle((c0 - 0.5, 0), c1 - c0 + 1, ax_bar.get_ylim()[1],
                                   facecolor=BAND, edgecolor="none", zorder=0))
    for y in (25, 50, 75, 100):
        ax_bar.axhline(y, color=GRID, lw=0.6, zorder=1)
    ax_bar.bar(range(ncol), counts, width=0.62, color=INK_SECOND, zorder=2)
    for i, v in enumerate(counts):
        if v >= 10:
            ax_bar.text(i, v + 2.0, f"{v:.0f}", ha="center", va="bottom",
                        fontsize=6.4, color=INK, fontweight="bold", zorder=3)
    ax_bar.set_ylabel("Clones differing\nfrom 3D7 (%)", fontsize=7.6, color=INK,
                      labelpad=5)
    ax_bar.set_yticks([0, 25, 50, 75, 100])
    ax_bar.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=6.8,
                           color=INK_SECOND)
    ax_bar.set_xticks([])
    ax_bar.tick_params(length=0)
    ax_bar.text(ncol - 0.5, ax_bar.get_ylim()[1] * 1.04, title, fontsize=9.5,
                color=INK, ha="right", va="bottom", fontweight="bold")
    letter_x = 0.012
    fig.text(letter_x, (bottom + mat_h + gap + band_h + gap + bar_h + 0.06) / fig_h,
             "a", fontsize=11, fontweight="bold", color=INK, ha="left", va="bottom")
    fig.text(letter_x, (bottom + mat_h + 0.06) / fig_h, "b", fontsize=11,
             fontweight="bold", color=INK, ha="left", va="bottom")

    # The epitope strip between the two panels.
    ax_band.set_ylim(0, 1)
    ax_band.axis("off")
    for name, c0, c1 in spans:
        ax_band.add_patch(Rectangle((c0 - 0.5, 0.30), c1 - c0 + 1, 0.42,
                                    facecolor=BAND, edgecolor=RULE, lw=0.6))
        ax_band.text((c0 + c1) / 2, 0.51, name, ha="center", va="center",
                     fontsize=8.4, color=INK, fontweight="bold")
    outside = [i for i in range(ncol)
               if not any(c0 <= i <= c1 for _, c0, c1 in spans)]
    for i in outside:
        ax_band.plot([i], [0.51], marker="_", ms=5, color=INK_MUTED, lw=1)

    # Panel b, the matrix. Row 0 is 3D7, then the haplotypes.
    ax_mat.set_ylim(nrow - 0.5, -0.5)
    ax_mat.axis("off")
    for name, c0, c1 in spans:
        ax_mat.add_patch(Rectangle((c0 - 0.5, -0.5), c1 - c0 + 1, nrow,
                                   facecolor=BAND, edgecolor="none", zorder=0))
    for i, pos in enumerate(positions):
        ax_mat.text(i, nrow - 0.18, str(pos), ha="center", va="top", fontsize=6.4,
                    color=INK_SECOND, rotation=90)
        ax_mat.text(i, -0.5 - 0.30, ref_at[pos], ha="center", va="center",
                    fontsize=7.4, color=INK, family="DejaVu Sans Mono",
                    fontweight="bold")
    ax_mat.text(-0.9, -0.5 - 0.30, "3D7", ha="right", va="center", fontsize=7.6,
                color=INK, fontweight="bold")
    ax_mat.plot([-0.5, ncol - 0.5], [0.05, 0.05], color=RULE, lw=0.9, zorder=3)

    for r, h in enumerate(shown):
        y = r + 0.6
        identical = all(h["sequence"][ref_index[p]] in MISSING
                        or h["sequence"][ref_index[p]] == ref_at[p] for p in positions)
        label = f"{h['haplotype']} = 3D7" if identical else h["haplotype"]
        ax_mat.text(-0.9, y, label, ha="right", va="center", fontsize=6.8,
                    color=INK, fontweight="bold" if identical else "normal")
        for i, pos in enumerate(positions):
            aa = h["sequence"][ref_index[pos]]
            if aa in MISSING:
                ax_mat.add_patch(Rectangle((i - 0.44, y - 0.40), 0.88, 0.80,
                                           facecolor="#FFFFFF", edgecolor=RULE,
                                           lw=0.5, hatch="///", zorder=2))
                continue
            if aa == ref_at[pos]:
                ax_mat.plot([i], [y], marker=".", ms=2.0, color=INK_MUTED, zorder=2)
                continue
            colour = CLASS_COLOR.get(AA_CLASS.get(aa, ""), UNKNOWN_COLOR)
            ax_mat.add_patch(Rectangle((i - 0.44, y - 0.40), 0.88, 0.80,
                                       facecolor=colour, edgecolor="none", zorder=2))
            ax_mat.text(i, y, aa, ha="center", va="center", fontsize=7.0,
                        color=ink_on(colour), family="DejaVu Sans Mono",
                        fontweight="bold", zorder=3)
    # Counts beside each row.
    ax_freq.set_ylim(nrow - 0.5, -0.5)
    top_n = max(h["n_clones"] for h in shown)
    ax_freq.set_xlim(0, top_n * 1.52)
    ax_freq.axis("off")
    for r, h in enumerate(shown):
        y = r + 0.6
        ax_freq.add_patch(Rectangle((0, y - 0.30), h["n_clones"], 0.60,
                                    facecolor=RULE, edgecolor="none"))
        ax_freq.text(h["n_clones"] + top_n * 0.05, y,
                     f"{h['n_clones']}  ({100.0 * h['n_clones'] / total_clones:.1f}%)",
                     ha="left", va="center", fontsize=6.4, color=INK)
    ax_freq.text(0, -0.5 - 0.30, "Clones", ha="left", va="center", fontsize=7.6,
                 color=INK, fontweight="bold")

    # Legend and the sentence that scopes the figure.
    handles = [Patch(facecolor=CLASS_COLOR[c], edgecolor="none", label=c)
               for c in CLASS_ORDER]
    leg = fig.legend(handles=handles, loc="lower left",
                     bbox_to_anchor=(left / fig_w, 0.62 / fig_h), ncol=4, frameon=False,
                     fontsize=7.4, handlelength=1.0, handleheight=1.0,
                     columnspacing=1.2, title="Class of the substituted residue")
    leg.get_title().set_fontsize(7.4)
    leg.get_title().set_color(INK)
    leg.get_title().set_fontweight("bold")
    for t in leg.get_texts():
        t.set_color(INK)

    hidden = len(haps) - len(shown)
    hidden_clones = sum(h["n_clones"] for h in haps[len(shown):])
    note = (f"A dot marks a residue identical to 3D7, a hatched cell a residue the "
            f"translation could not resolve. {len(haps)} protein haplotypes "
            f"were found among {total_clones} clones; the {len(shown)} shown carry "
            f"{total_clones - hidden_clones} of them "
            f"({100.0 * (total_clones - hidden_clones) / total_clones:.0f}%).")
    if hidden:
        note += (f" The remaining {hidden} are rarer and are listed in "
                 f"protein_haplotypes.tsv.")
    # Wrapped to the drawn width, so the sentence never runs off the page.
    chars = int((fig_w - left - 0.12) / 0.0455)
    fig.text(left / fig_w, 0.14 / fig_h, "\n".join(textwrap.wrap(note, chars)),
             fontsize=6.6, color=INK_SECOND, ha="left", va="bottom", linespacing=1.5)

    fig.savefig(out_pdf, facecolor=SURFACE)
    fig.savefig(out_png, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Protein level figure of the CSP C terminus against 3D7.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--codon-map", type=Path, required=True)
    ap.add_argument("--callable-start", type=int, required=True)
    ap.add_argument("--callable-end", type=int, required=True)
    ap.add_argument("--th2r", type=int, nargs=2, default=[311, 327])
    ap.add_argument("--th3r", type=int, nargs=2, default=[352, 367])
    ap.add_argument("--ntc-pattern", default="NTC")
    ap.add_argument("--positive-pattern", default="3D7|3d7")
    ap.add_argument("--min-clones", type=int, default=3,
                    help="a haplotype is drawn when at least this many clones carry it")
    ap.add_argument("--max-rows", type=int, default=45,
                    help="hard cap on drawn rows, so the figure stays a page")
    ap.add_argument("--title", default="PfCSP C terminus against 3D7")
    ap.add_argument("--dpi", type=int, default=600,
                    help="raster resolution; 600 is what a journal asks for")
    args = ap.parse_args()

    protein = args.out_dir / "08_protein" / "haplotypes_protein.faa"
    if not protein.exists():
        sys.exit(f"[ERROR] translated clones not found: {protein}, run step 5 first")

    ref = reference_residues(args.codon_map, args.callable_start, args.callable_end)
    if not ref:
        sys.exit("[ERROR] no codon of the map falls inside the callable window")
    window = [p for p, _ in ref]
    ref_at = dict(ref)
    print(f"[INFO] window: amino acids {window[0]} to {window[-1]}, {len(window)} residues")

    global ref_index
    ref_index = {p: i for i, p in enumerate(window)}

    haps, total_clones = build_matrix(protein, window, args.ntc_pattern,
                                      args.positive_pattern)
    variable = [p for p in window
                if any(h["sequence"][ref_index[p]] not in MISSING
                       and h["sequence"][ref_index[p]] != ref_at[p] for h in haps)]
    print(f"[INFO] {total_clones} clones, {len(haps)} protein haplotypes, "
          f"{len(variable)} variable positions")
    if not variable:
        sys.exit("[ERROR] no position differs from 3D7, nothing to draw")

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

    shown = [h for h in haps if h["n_clones"] >= args.min_clones][:args.max_rows]
    if not shown:
        shown = haps[:args.max_rows]
    print(f"[INFO] drawing {len(shown)} haplotypes at {args.min_clones} clones or more")

    out = args.out_dir / "10_figures"
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for h in haps:
        changes = [f"{ref_at[p]}{p}{h['sequence'][ref_index[p]]}" for p in variable
                   if h["sequence"][ref_index[p]] not in MISSING
                   and h["sequence"][ref_index[p]] != ref_at[p]]
        unresolved = sum(1 for p in variable
                         if h["sequence"][ref_index[p]] in MISSING)
        rows.append({
            "haplotype": h["haplotype"], "n_clones": h["n_clones"],
            "pct_clones": round(100.0 * h["n_clones"] / total_clones, 2),
            "n_samples": h["n_samples"], "n_substitutions": len(changes),
            "n_unresolved": unresolved,
            "substitutions": ";".join(changes) or "none, identical to 3D7",
            "drawn": "yes" if h in shown else "no",
            "sequence": h["sequence"],
            "members": ";".join(h["members"][:30]),
        })
    write_tsv(out / "protein_haplotypes.tsv",
              ["haplotype", "n_clones", "pct_clones", "n_samples", "n_substitutions",
               "n_unresolved", "substitutions", "drawn", "sequence", "members"], rows)
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

    draw(out / "csp_protein_vs_3d7.pdf", out / "csp_protein_vs_3d7.png",
         haps, shown, variable, ref_at, counts, total_clones,
         tuple(args.th2r), tuple(args.th3r), args.title, args.dpi)
    print(f"[OK] {out / 'csp_protein_vs_3d7.pdf'}")
    print(f"[OK] {out / 'csp_protein_vs_3d7.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
