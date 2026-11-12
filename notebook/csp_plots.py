#!/usr/bin/env python3
"""
Drawing only. Every function takes a table written by process/07_plot_tables.py.

Nothing here reads the pipeline output, computes a statistic or holds a colour. A figure
is therefore reproducible from 10_table_for_plot alone, and a palette is changed by
editing plot_style.tsv and rerunning the table step, not by editing this file.

Author: Khadim Gueye
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def load_tables(table_dir):
    """Every plot table in the folder, keyed by file stem."""
    d = Path(table_dir)
    out = {p.stem: pd.read_csv(p, sep="\t") for p in sorted(d.glob("*.tsv"))}
    missing = [n for n in ("plot_style", "index") if n not in out]
    if missing:
        raise FileNotFoundError(f"{d} is missing {', '.join(missing)}")
    return out


def palette(style, role):
    """The entries of one role, in their declared order."""
    sub = style[style.role == role].sort_values("order")
    return {r.key: r for _, r in sub.iterrows()}


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontweight("bold")


def _bands(ax, style, span_text=False, y=None):
    """The two epitope bands, from the style table."""
    for _, b in style[style.role == "epitope_band"].sort_values("order").iterrows():
        ax.axvspan(float(b.start), float(b.end), color=b.color, zorder=0)
    if not span_text:
        return
    txt = palette(style, "epitope_text")
    for key, b in txt.items():
        ax.text((float(b.start) + float(b.end)) / 2, y, key, ha="center", va="bottom",
                fontsize=9 if span_text == "small" else 10, fontweight="bold",
                color=b.color)


def _stems(ax, df, top, dx, style, fontsize=6):
    """One stem per substitution. Labels sharing a position sit side by side.

    The labels are rotated, so stacking them vertically would still overprint. They are
    spread along x instead, from a common baseline above the tallest stem of the group.
    """
    line = style[(style.role == "stem") & (style.key == "line")].color.iloc[0]
    head = style[(style.role == "stem") & (style.key == "point")].color.iloc[0]
    # A table holding several panels reads x as text when another panel has text there.
    df = df.assign(x=pd.to_numeric(df.x), value=pd.to_numeric(df.value))
    for pos, grp in df.groupby("x"):
        grp = grp.sort_values("value", ascending=False)
        base = float(grp.value.max()) + top * 0.04
        shown = grp[grp.show_label == 1]
        for _, r in grp.iterrows():
            ax.vlines(pos, 0, r.value, color=line, lw=0.9, zorder=2)
            ax.plot([pos], [r.value], "o", ms=4, color=head, zorder=3)
        for i, (_, r) in enumerate(shown.iterrows()):
            x = pos + (i - (len(shown) - 1) / 2) * dx
            ax.text(x, base, r.label, rotation=90, ha="center", va="bottom",
                    fontsize=fontsize, fontweight="bold", zorder=4)


def _letters(pairs, dx=-0.08, dy=1.06, size=13):
    for ax, letter in pairs:
        ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=size,
                fontweight="bold", va="bottom")


# --- the figures -----------------------------------------------------------------

def fig_sample_counts(t, style, title=None):
    d = t[t.panel == "A"].sort_values("order")
    edge = style[(style.role == "ink") & (style.key == "bar_edge")].color.iloc[0]
    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    bars = ax.bar(d.x, d.value, color=list(d.color), edgecolor=edge, linewidth=1.1,
                  width=0.62)
    for b, v, lab in zip(bars, d.value, d.label):
        ax.text(b.get_x() + b.get_width() / 2, v + d.value.max() * 0.015, str(lab),
                ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_ylabel(d.y_axis_label.iloc[0], fontsize=12, fontweight="bold")
    ax.set_ylim(0, d.value.max() * 1.12)
    _style(ax)
    ax.tick_params(labelsize=11)
    for lab in ax.get_xticklabels():
        lab.set_fontweight("bold")
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig


def fig_moi(t, style):
    pts = t[t.panel == "points"]
    summ = t[t.panel == "summary"]
    groups = pts.drop_duplicates("x").sort_values("order")
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    parts = ax.violinplot([pts.loc[pts.x == g, "value"].values for g in groups.x],
                          showextrema=False, widths=0.85)
    for body, c in zip(parts["bodies"], groups.color):
        body.set_facecolor(c)
        body.set_edgecolor("black")
        body.set_linewidth(0.9)
        body.set_alpha(1.0)
    at = {g: i + 1 for i, g in enumerate(groups.x)}
    for _, r in summ.iterrows():
        ax.scatter([at[r.x]], [r.value], marker=r.marker, s=r["size"], color=r.color,
                   edgecolor="black", linewidth=0.8, zorder=3)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels(groups.x, rotation=45, ha="right")
    ax.set_xlabel(pts.x_axis_label.iloc[0], fontsize=11, fontweight="bold")
    ax.set_ylabel(pts.y_axis_label.iloc[0], fontsize=11, fontweight="bold")
    marks = summ.drop_duplicates("series").sort_values("order")
    ax.legend(handles=[plt.Line2D([], [], marker=r.marker, color="none",
                                  markersize=15 if r.marker == "*" else 10,
                                  markerfacecolor=r.color, markeredgecolor="black",
                                  label=r.label) for _, r in marks.iterrows()],
              frameon=False, loc="upper left", prop={"weight": "bold", "size": 10})
    _style(ax)
    fig.tight_layout()
    return fig


def fig_diversity(t, style):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6))
    for ax, panel in zip(axes, ("A", "B")):
        d = t[t.panel == panel].sort_values("order")
        bars = ax.bar(d.x, d.value, color=list(d.color), edgecolor="black",
                      linewidth=1.1, width=0.6)
        for b, v, lab in zip(bars, d.value, d.label):
            ax.text(b.get_x() + b.get_width() / 2, v + d.value.max() * 0.02, str(lab),
                    ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.set_ylabel(d.y_axis_label.iloc[0], fontsize=11, fontweight="bold")
        ax.set_ylim(0, d.value.max() * 1.18)
        _style(ax)
    fig.tight_layout()
    return fig


def fig_substitutions(t, style):
    fig = plt.figure(figsize=(13.0, 8.6))
    ax_a = fig.add_axes([0.06, 0.60, 0.38, 0.34])
    ax_b = fig.add_axes([0.54, 0.60, 0.42, 0.34])
    ax_c = fig.add_axes([0.06, 0.08, 0.90, 0.38])

    a = t[t.panel == "A"]
    for country, grp in a.groupby("series", sort=False):
        grp = grp.sort_values("x")
        ax_a.plot(grp.x, grp.value, lw=1.1, color=grp.color.iloc[0], label=country)
    _bands(ax_a, style)
    _bands(ax_a, style, span_text="small", y=ax_a.get_ylim()[1])
    ax_a.set_xlabel(a.x_axis_label.iloc[0], fontsize=10, fontweight="bold")
    ax_a.set_ylabel("Fraction of samples heterozygous", fontsize=10, fontweight="bold")
    ax_a.legend(frameon=False, prop={"weight": "bold", "size": 9})
    _style(ax_a)

    b = t[t.panel == "B"]
    order = b.drop_duplicates("series").sort_values("order")
    grid = b.pivot_table(index="y", columns="x", values="value").values
    ax_b.imshow(grid, aspect="auto", cmap="Blues", vmin=0, vmax=1.4,
                interpolation="nearest")
    ax_b.set_yticks(range(len(order)))
    ax_b.set_yticklabels(order.series, fontsize=8, fontweight="bold")
    ax_b.set_xticks([])
    ax_b.set_xlabel(f"Samples (n = {int(b.n_samples_in_panel.iloc[0])})", fontsize=10,
                    fontweight="bold")
    g = palette(style, "grid")
    ax_b.legend(handles=[Patch(facecolor=g["present"].color, label=g["present"].label),
                         Patch(facecolor=g["absent"].color, edgecolor="black",
                               label=g["absent"].label)],
                loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2, frameon=False,
                prop={"weight": "bold", "size": 9})
    for s in ax_b.spines.values():
        s.set_visible(False)

    c = t[t.panel == "C"]
    top = float(c.panel_top.iloc[0])
    _bands(ax_c, style)
    _stems(ax_c, c, top, dx=0.9, style=style)
    _bands(ax_c, style, span_text=True, y=top * 1.06)
    ax_c.set_xlabel(c.x_axis_label.iloc[0], fontsize=10, fontweight="bold")
    ax_c.set_ylabel("Samples with substitution", fontsize=10, fontweight="bold")
    ax_c.set_ylim(0, top * 1.30)
    _style(ax_c)

    _letters(((ax_a, "A"), (ax_b, "B"), (ax_c, "C")))
    return fig


def fig_vaccine_match(t, style):
    fig = plt.figure(figsize=(12.0, 9.4))
    ax_a = fig.add_axes([0.07, 0.70, 0.36, 0.25])
    ax_b = fig.add_axes([0.58, 0.70, 0.36, 0.25])
    ax_c1 = fig.add_axes([0.07, 0.38, 0.40, 0.22])
    ax_c2 = fig.add_axes([0.55, 0.38, 0.40, 0.22])
    ax_d = fig.add_axes([0.07, 0.06, 0.34, 0.22])
    ax_e = fig.add_axes([0.56, 0.06, 0.32, 0.22])

    for ax, panel, title in (
            (ax_a, "A", "Vaccine match status\n(sample consensus level)"),
            (ax_b, "B", "Vaccine match status\n(individual clone level)")):
        d = t[t.panel == panel].sort_values("order")
        bars = ax.bar(d.x, d.value, color=list(d.color), edgecolor="black",
                      linewidth=1.0, width=0.6)
        for bar, v, lab in zip(bars, d.value, d.label):
            ax.text(bar.get_x() + bar.get_width() / 2, v + d.value.max() * 0.02,
                    str(lab), ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.set_ylabel(d.y_axis_label.iloc[0], fontsize=10, fontweight="bold")
        ax.set_ylim(0, d.value.max() * 1.16)
        ax.set_title(title, fontsize=10, fontweight="bold")
        _style(ax)

    txt = palette(style, "epitope_text")
    band = palette(style, "epitope_band")
    for ax, name in ((ax_c1, "TH2R"), (ax_c2, "TH3R")):
        d = t[t.panel == f"C_{name}"]
        ax.axvspan(float(band[name].start), float(band[name].end),
                   color=band[name].color, zorder=0)
        _stems(ax, d, float(d.panel_top.iloc[0]), dx=0.55, style=style)
        ax.set_xlim(float(d.x_min.iloc[0]), float(d.x_max.iloc[0]))
        ax.set_ylim(0, 150)
        ax.set_xlabel("CSP amino acid position", fontsize=9, fontweight="bold")
        ax.set_ylabel(d.y_axis_label.iloc[0], fontsize=9, fontweight="bold")
        ax.set_title(name, fontsize=10, fontweight="bold", color=txt[name].color)
        _style(ax)

    d = t[t.panel == "D"].sort_values("order")
    bars = ax_d.bar(d.x, d.value, color=list(d.color), edgecolor="black", linewidth=1.0,
                    width=0.6)
    for bar, v, lab in zip(bars, d.value, d.label):
        ax_d.text(bar.get_x() + bar.get_width() / 2, v + max(d.value.max(), 1) * 0.03,
                  str(lab), ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_d.set_ylabel(d.y_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_d.set_ylim(0, max(d.value.max(), 1) * 1.2)
    ax_d.set_title("Mismatch localisation", fontsize=10, fontweight="bold")
    _style(ax_d)

    e = t[t.panel == "E"].copy()
    for c in ("x", "y", "value"):
        e[c] = pd.to_numeric(e[c])
    n = int(e.x.max()) + 1
    mat = e.pivot_table(index="y", columns="x", values="value").values
    im = ax_e.imshow(mat, cmap="Blues")
    for _, r in e.iterrows():
        ax_e.text(r.x, r.y, str(r.label), ha="center", va="center", fontsize=9,
                  fontweight="bold",
                  color="white" if r.value > mat.max() / 2 else "black")
    cols = e.drop_duplicates("x").sort_values("x").x_label
    rows = e.drop_duplicates("y").sort_values("y").y_label
    ax_e.set_xticks(range(n), cols, fontweight="bold", fontsize=8, rotation=30,
                    ha="right")
    ax_e.set_yticks(range(n), rows, fontweight="bold", fontsize=8)
    ax_e.set_xlabel(e.x_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_e.set_ylabel(e.y_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_e.set_title("TH2R x TH3R match status", fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=ax_e, fraction=0.046, label="Samples")

    _letters(((ax_a, "A"), (ax_b, "B"), (ax_c1, "C"), (ax_d, "D"), (ax_e, "E")),
             dx=-0.14, dy=1.08)
    return fig


def fig_structure(t, style):
    cmap = style[(style.role == "colormap") & (style.key == "fst")].color.iloc[0]
    fig = plt.figure(figsize=(15.0, 4.8))
    ax_a = fig.add_axes([0.05, 0.18, 0.22, 0.70])
    ax_b = fig.add_axes([0.37, 0.18, 0.16, 0.70])
    ax_c = fig.add_axes([0.63, 0.18, 0.24, 0.70])
    for ax in (ax_a, ax_b):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    a = t[t.panel == "A"].copy()
    for c in ("x", "y", "value", "lower_triangle"):
        a[c] = pd.to_numeric(a[c])
    names = a.drop_duplicates("x").sort_values("x").x_label.tolist()
    mat = np.full((len(names), len(names)), np.nan)
    for _, r in a[a.lower_triangle == 1].iterrows():
        mat[int(r.y), int(r.x)] = r.value
    im = ax_a.imshow(np.ma.masked_invalid(mat), cmap=cmap, vmin=0, aspect="auto")
    ax_a.set_xticks(range(len(names)), names, rotation=45, ha="right", fontsize=8,
                    fontweight="bold")
    ax_a.set_yticks(range(len(names)), names, fontsize=8, fontweight="bold")
    ax_a.set_xlabel(a.x_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_a.set_ylabel(a.y_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_a.set_title("FST, state level", fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=ax_a, fraction=0.046, label="FST")

    b = t[t.panel == "B"].copy()
    if len(b):
        b["value"] = pd.to_numeric(b.value)
        r = b.iloc[0]
        im2 = ax_b.imshow([[np.nan, np.nan], [r.value, np.nan]], cmap=cmap, vmin=0,
                          vmax=max(r.value * 10, 0.05), aspect="auto")
        ax_b.text(0, 1, str(r.label), ha="center", va="center", fontsize=11,
                  fontweight="bold")
        ax_b.set_xticks([0, 1], [r.x_label, r.y_label], fontweight="bold", fontsize=9)
        ax_b.set_yticks([0, 1], [r.x_label, r.y_label], fontweight="bold", fontsize=9)
        ax_b.set_title("FST, country level", fontsize=10, fontweight="bold")
        fig.colorbar(im2, ax=ax_b, fraction=0.046, label="FST")
    else:
        ax_b.axis("off")

    c = t[t.panel == "C"].copy()
    for col in ("x", "y", "size", "order"):
        c[col] = pd.to_numeric(c[col])
    for series, grp in sorted(c.groupby("series"), key=lambda kv: kv[1].order.iloc[0]):
        ax_c.scatter(grp.x, grp.y, s=grp["size"], color=grp.color.iloc[0],
                     marker=grp.marker.iloc[0], edgecolor="white",
                     linewidth=1.2 if grp.marker.iloc[0] == "*" else 0.5,
                     alpha=1.0 if grp.marker.iloc[0] == "*" else 0.85,
                     zorder=4 if grp.marker.iloc[0] == "*" else 2, label=series)
    ref = c[c.marker == "*"]
    if len(ref):
        r = ref.iloc[0]
        ax_c.annotate(str(r.label), (r.x, r.y), textcoords="offset points",
                      xytext=(10, 6), fontsize=10, fontweight="bold", color=r.color,
                      zorder=5)
    ax_c.set_xlabel(c.x_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_c.set_ylabel(c.y_axis_label.iloc[0], fontsize=9, fontweight="bold")
    ax_c.set_title("PCoA with the 3D7 reference", fontsize=10, fontweight="bold")
    ax_c.legend(frameon=False, prop={"weight": "bold", "size": 9},
                loc="center left", bbox_to_anchor=(1.02, 0.5))
    _style(ax_c)

    _letters(((ax_a, "A"), (ax_b, "B"), (ax_c, "C")), dx=-0.18)
    return fig
