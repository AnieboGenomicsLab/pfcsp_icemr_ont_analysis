#!/usr/bin/env python3
"""
Write one table per figure, holding everything that figure draws.

The rule this step exists for: a figure is never redrawn by rerunning the pipeline. The
pipeline writes the numbers and the visual encoding into 10_table_for_plot, and the
notebook reads only that folder. A year from now the figures still rebuild, and the same
numbers can be drawn another way without touching the analysis.

Every table shares the leading columns panel, series, x, y, value, n, pct, color, size,
marker, order, label, then adds its own. plot_style.tsv declares each palette once, so a
colour is changed in one place and this step rerun. index.tsv lists what is here.

Tables for a manuscript belong in 12_tables, not here.

Author: Khadim Gueye
"""

import argparse
import itertools
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import csp_stats as S

HEAD = ["panel", "series", "x", "y", "value", "n", "pct",
        "color", "size", "marker", "order", "label"]

COUNTRY_COLOR = {"Nigeria": "#2E6DB4", "Cameroon": "#D95F43", "Control": "#9E9E9E"}
REF_COLOR = "#111111"
GOLD = "#F2C230"
ACCENT = "#1F4E79"
STATUS_COLOR = {"match": "#4FA64F", "mismatch": "#B02418", "incomplete": "#9E9E9E"}
LOCALISATION_COLOR = {"TH2R only": "#B02418", "TH3R only": "#7FB3D9",
                      "TH2R+TH3R": "#B02418"}


def write(path, rows, note=""):
    """One table, the shared columns first and the figure's own columns after."""
    df = pd.DataFrame(rows)
    extra = [c for c in df.columns if c not in HEAD]
    df = df.reindex(columns=[c for c in HEAD if c in df.columns] + extra)
    df.to_csv(path, sep="\t", index=False)
    print(f"[OK] {path.name} ({len(df)} rows){note}")
    return len(df)


def style_table(th2r, th3r):
    """Every palette the figures use, declared once."""
    rows = []

    def add(role, key, color, order, label="", start="", end=""):
        rows.append({"role": role, "key": key, "color": color, "order": order,
                     "label": label or key, "start": start, "end": end})

    for i, (k, c) in enumerate(COUNTRY_COLOR.items(), start=1):
        add("country", k, c, i)
    add("reference", "3D7", REF_COLOR, 1, "3D7 reference")
    add("violin", "state", "#A8C4E0", 1, "Nigerian state")
    add("violin", "country_b", "#F2B7AC", 2, "Second country")
    add("summary_marker", "Mean", GOLD, 1)
    add("summary_marker", "Median", GOLD, 2)
    for i, (k, c) in enumerate(STATUS_COLOR.items(), start=1):
        add("status", k, c, i, k.capitalize())
    for i, (k, c) in enumerate(LOCALISATION_COLOR.items(), start=1):
        add("localisation", k, c, i)
    add("epitope_band", "TH2R", "#FBE3DA", 1, "TH2R", th2r[0], th2r[1])
    add("epitope_band", "TH3R", "#DCE9F5", 2, "TH3R", th3r[0], th3r[1])
    add("epitope_text", "TH2R", "#C0504D", 1, "TH2R", th2r[0], th2r[1])
    add("epitope_text", "TH3R", ACCENT, 2, "TH3R", th3r[0], th3r[1])
    add("stem", "line", "#555555", 1, "stem line")
    add("stem", "point", "#B02418", 2, "stem head")
    add("grid", "present", "#08306B", 1, "Present")
    add("grid", "absent", "#F7FBFF", 2, "Absent")
    add("ink", "bar_edge", "#000000", 1, "bar outline")
    add("ink", "bar_default", "#3F3F3F", 2, "bar with no palette entry")
    add("colormap", "fst", "Blues", 1, "FST heatmaps")
    add("colormap", "counts", "Blues", 2, "count heatmaps")
    add("colormap", "distance_to_3d7", "YlOrRd", 3, "distance to 3D7")
    return rows


# --- one function per figure -----------------------------------------------------

def table_fig1(meta, n_controls):
    """Sample counts, one bar per country plus the controls."""
    counts = meta["country"].value_counts()
    rows = []
    order = 0
    for c in ("Nigeria", "Cameroon"):
        if c in counts.index:
            order += 1
            rows.append({"panel": "A", "series": c, "x": c, "value": int(counts[c]),
                         "n": int(counts[c]), "color": COUNTRY_COLOR[c],
                         "order": order, "label": str(int(counts[c]))})
    if n_controls:
        order += 1
        rows.append({"panel": "A", "series": "Control", "x": "Control",
                     "value": int(n_controls), "n": int(n_controls),
                     "color": COUNTRY_COLOR["Control"], "order": order,
                     "label": str(int(n_controls))})
    for r in rows:
        r["y_axis_label"] = "Number of samples"
    return rows


def table_fig2(moi, country_a, country_b):
    """MOI. One row per sample for the violins, then the mean and median per group."""
    groups = sorted(moi.loc[moi.country == country_a, "state"].unique())
    has_b = bool((moi.country == country_b).any())
    rows = []
    for i, state in enumerate(groups, start=1):
        sub = moi[moi.state == state]
        for _, r in sub.iterrows():
            rows.append({"panel": "points", "series": state, "x": state,
                         "value": int(r.MOI), "color": "#A8C4E0", "order": i,
                         "sample": r["sample"], "country": r.country})
    if has_b:
        i = len(groups) + 1
        for _, r in moi[moi.country == country_b].iterrows():
            rows.append({"panel": "points", "series": country_b, "x": country_b,
                         "value": int(r.MOI), "color": "#F2B7AC", "order": i,
                         "sample": r["sample"], "country": r.country})
    pts = pd.DataFrame(rows)
    for (grp, order), sub in pts.groupby(["x", "order"]):
        for stat, val, marker, size in (("Mean", sub.value.mean(), "*", 230),
                                        ("Median", sub.value.median(), "o", 95)):
            rows.append({"panel": "summary", "series": stat, "x": grp,
                         "value": round(float(val), 4), "n": len(sub), "color": GOLD,
                         "size": size, "marker": marker, "order": order,
                         "label": stat})
    for r in rows:
        r["y_axis_label"] = "Multiplicity of infection (MOI)"
        r["x_axis_label"] = "State"
    return rows


def table_fig3(clones, country_a, country_b):
    """Nucleotide diversity and haplotype diversity, one bar per country."""
    rows = []
    for panel, fn, label in (("A", S.nucleotide_diversity, "Nucleotide diversity (pi)"),
                             ("B", S.haplotype_diversity, "Haplotype diversity (Hd)")):
        order = 0
        for c in (country_a, country_b):
            if not (clones.country == c).any():
                continue
            order += 1
            seqs = clones.loc[clones.country == c, "nt"].tolist()
            v = float(fn(seqs))
            rows.append({"panel": panel, "series": c, "x": c, "value": round(v, 6),
                         "n": len(seqs), "color": COUNTRY_COLOR[c], "order": order,
                         "label": f"{v:.4f}", "y_axis_label": label,
                         "n_haplotypes": int(pd.Series(seqs).nunique())})
    return rows


def table_fig4(het, subs, n_by_country, n_with_consensus, window, top_n=15):
    """Heterozygosity by position, the presence grid, and prevalence by position."""
    rows = []
    span = list(range(window[0], window[1] + 1))

    # Panel A. value keeps the denominator the figure has always used, the metadata
    # count. n_with_consensus is the stricter denominator, carried so the fraction can
    # be recomputed without rerunning anything.
    for country in ("Nigeria", "Cameroon"):
        total = n_by_country.get(country, 0)
        if not total:
            continue
        pool = het[het.country == country]
        frac = (pool.groupby("aa_position")["sample"].nunique()
                .reindex(span, fill_value=0))
        for pos in span:
            k = int(frac.loc[pos])
            rows.append({"panel": "A", "series": country, "x": int(pos),
                         "value": round(k / total, 6), "n": k,
                         "color": COUNTRY_COLOR[country], "order": 1,
                         "n_denominator": int(total),
                         "n_with_consensus": int(n_with_consensus.get(country, 0))})

    # Panel B, the commonest epitope changes, present or absent per sample.
    ep = subs[subs.region.isin(("TH2R", "TH3R"))]
    order = [s for s, _ in Counter(ep.substitution).most_common(top_n)]
    samples = sorted(subs["sample"].unique())
    carriers = {sub: set(g["sample"]) for sub, g in ep.groupby("substitution")}
    for i, sub in enumerate(order):
        for j, sample in enumerate(samples):
            on = sample in carriers.get(sub, ())
            rows.append({"panel": "B", "series": sub, "x": j, "y": i,
                         "value": int(on), "color": "#08306B" if on else "#F7FBFF",
                         "order": i, "label": sub, "sample": sample,
                         "n_samples_in_panel": len(samples)})

    # Panel C, how many samples carry each change, by position.
    per = subs.groupby(["aa_position", "substitution"])["sample"].nunique()
    top = float(per.max())
    n_samples = len(samples)
    for (pos, sub), v in per.items():
        region = subs.loc[subs.substitution == sub, "region"].iloc[0]
        rows.append({"panel": "C", "series": sub, "x": int(pos), "value": int(v),
                     "n": int(v), "pct": round(100 * v / n_samples, 2),
                     "color": "#B02418", "order": int(pos), "label": sub,
                     "region": region, "show_label": int(v >= top * 0.05),
                     "panel_top": top})
    for r in rows:
        r.setdefault("x_axis_label", "CSP amino acid position")
    return rows


def table_fig7(cons, clone_status, subs_clone, th2r, th3r):
    """The five panels of the vaccine comparison."""
    rows = []
    ORDER = ["match", "mismatch", "incomplete"]

    for panel, series, total, axis in (
            ("A", cons["vaccine_status"], len(cons), f"Samples (of {len(cons)})"),
            ("B", clone_status["vaccine_status"], len(clone_status),
             f"Clones (of {len(clone_status)})")):
        for i, k in enumerate(ORDER, start=1):
            v = int((series == k).sum())
            rows.append({"panel": panel, "series": k, "x": k.capitalize(), "value": v,
                         "n": v, "pct": round(100 * v / total, 2) if total else "",
                         "color": STATUS_COLOR[k], "order": i, "label": str(v),
                         "y_axis_label": axis, "n_total": total})

    # Panels C, the substitutions responsible, as a percentage of clones.
    n_clones = len(clone_status)
    for name, (lo, hi) in (("TH2R", th2r), ("TH3R", th3r)):
        inside = subs_clone[subs_clone.region == name]
        per = inside.groupby(["aa_position", "substitution"]).size()
        for (pos, sub), v in per.items():
            pct = 100 * v / n_clones
            rows.append({"panel": f"C_{name}", "series": sub, "x": int(pos),
                         "value": round(pct, 3), "n": int(v), "pct": round(pct, 2),
                         "color": "#B02418", "order": int(pos), "label": sub,
                         "region": name, "show_label": int(pct >= 1.0),
                         "panel_top": 100.0, "x_min": lo - 2, "x_max": hi + 3,
                         "n_total": n_clones, "y_axis_label": "Clones (%)"})

    # Panel D, where a mismatching sample mismatches.
    loc = cons.loc[cons.vaccine_status == "mismatch", "mismatch_localisation"]
    keys = ["TH2R only", "TH3R only", "TH2R+TH3R"]
    vals = [int((loc == k).sum()) for k in keys]
    for i, (k, v) in enumerate(zip(keys, vals), start=1):
        rows.append({"panel": "D", "series": k, "x": k, "value": v, "n": v,
                     "color": LOCALISATION_COLOR[k], "order": i, "label": str(v),
                     "n_total": sum(vals),
                     "y_axis_label": f"Mismatched samples (of {sum(vals)})"})

    # Panel E, the two epitopes scored against each other.
    tab = pd.crosstab(cons["th2r_status"], cons["th3r_status"])
    tab = tab.reindex(index=ORDER, columns=ORDER, fill_value=0)
    for i, r in enumerate(ORDER):
        for j, c in enumerate(ORDER):
            rows.append({"panel": "E", "series": f"{r} x {c}", "x": j, "y": i,
                         "value": int(tab.values[i, j]), "n": int(tab.values[i, j]),
                         "order": i * 3 + j, "label": str(int(tab.values[i, j])),
                         "x_label": c.capitalize(), "y_label": r.capitalize(),
                         "x_axis_label": "TH3R status", "y_axis_label": "TH2R status"})
    return rows


def table_fig8(clones, ref_seq, country_a, country_b, min_per_group=2):
    """FST by state, FST by country, and the PCoA carrying 3D7."""
    rows = []

    pools = {}
    for state, grp in clones[clones.country == country_a].groupby("state"):
        if len(grp) >= min_per_group:
            pools[state] = grp["nt"].tolist()
    if (clones.country == country_b).any():
        pools[country_b] = clones.loc[clones.country == country_b, "nt"].tolist()
    names = sorted(pools)

    for i, j in itertools.combinations(range(len(names)), 2):
        v = S.hudson_fst(pools[names[i]], pools[names[j]])
        for a, b in ((i, j), (j, i)):
            rows.append({"panel": "A", "series": f"{names[a]} vs {names[b]}",
                         "x": b, "y": a, "value": round(float(v), 6),
                         "order": a * len(names) + b, "x_label": names[b],
                         "y_label": names[a], "n_x": len(pools[names[b]]),
                         "n_y": len(pools[names[a]]), "lower_triangle": int(a > b),
                         "x_axis_label": "State", "y_axis_label": "State"})

    if (clones.country == country_b).any():
        v = S.hudson_fst(clones.loc[clones.country == country_a, "nt"].tolist(),
                         clones.loc[clones.country == country_b, "nt"].tolist())
        rows.append({"panel": "B", "series": f"{country_a} vs {country_b}",
                     "x": 0, "y": 1, "value": round(float(v), 6), "order": 1,
                     "label": f"{v:.4f}", "x_label": country_a, "y_label": country_b,
                     "n_x": int((clones.country == country_a).sum()),
                     "n_y": int((clones.country == country_b).sum())})

    # Identical clones land on one spot, so each sequence is a single point sized by
    # how many clones carry it. The reference is one more point.
    grouped = clones.groupby(["nt", "country"]).size().reset_index(name="n_clones")
    seqs = grouped["nt"].tolist() + [ref_seq]
    labels = grouped["country"].tolist() + ["3D7 reference"]
    coords, explained = S.pcoa(seqs, labels)
    coords["n_clones"] = list(grouped["n_clones"]) + [1]

    # The smaller cohort is drawn last so it is not buried under the larger one.
    present = [c for c in (country_a, country_b) if c in set(labels)]
    present.sort(key=lambda c: -int((clones.country == c).sum()))
    draw_order = {c: i + 1 for i, c in enumerate(present)}
    for _, r in coords.iterrows():
        is_ref = r.label == "3D7 reference"
        rows.append({
            "panel": "C", "series": r.label,
            "x": round(float(r.PCoA1), 6), "y": round(float(r.PCoA2), 6),
            "n": int(r.n_clones),
            "color": REF_COLOR if is_ref else COUNTRY_COLOR.get(r.label, "#3F3F3F"),
            "size": 260 if is_ref else 18 + 9 * int(r.n_clones),
            "marker": "*" if is_ref else "o",
            "order": 99 if is_ref else draw_order.get(r.label, 50),
            "label": "3D7" if is_ref else "",
            "x_axis_label": f"PCoA1 ({explained[0]:.1f} %)",
            "y_axis_label": f"PCoA2 ({explained[1]:.1f} %)"})
    for i, v in enumerate(explained, start=1):
        rows.append({"panel": "C_meta", "series": f"PCoA{i}", "x": i,
                     "value": round(float(v), 4), "order": i,
                     "label": "percent of variation explained"})
    return rows


def table_fig10(phylo_dir):
    """The iTOL annotation files, in the order they are loaded onto the tree."""
    order = [
        ("tree", "nigeria.rooted.nwk", "the tree itself"),
        ("labels", "itol_labelcolors_origin.txt", "label colour, no ring"),
        ("ring 1", "itol_ringlabel_origin.txt", "ring title"),
        ("ring 1", "itol_colorstrip_origin.txt", "origin of the haplotype"),
        ("ring 2", "itol_ringlabel_our_clones.txt", "ring title"),
        ("ring 2", "itol_symbol_our_clones.txt", "clones in our samples, circles"),
        ("ring 3", "itol_ringlabel_region.txt", "ring title"),
        ("ring 3", "itol_piechart_by_region.txt", "region of the Pf8 samples"),
        ("ring 4", "itol_ringlabel_country.txt", "ring title"),
        ("ring 4", "itol_piechart_by_country.txt", "country of the Pf8 samples"),
        ("alternative", "itol_bar_our_clones.txt", "the same counts as bars"),
    ]
    rows = []
    for i, (ring, name, what) in enumerate(order, start=1):
        p = Path(phylo_dir) / name
        rows.append({"panel": ring, "series": name, "order": i, "label": what,
                     "file": str(p), "present": int(p.exists())})
    return rows


def main():
    ap = argparse.ArgumentParser(description="Write the plot tables for every figure.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--reference-dir", type=Path, required=True)
    ap.add_argument("--codon-map", type=Path, required=True)
    ap.add_argument("--amplicon-start", type=int, required=True)
    ap.add_argument("--callable-start", type=int, required=True)
    ap.add_argument("--callable-end", type=int, required=True)
    ap.add_argument("--th2r", type=int, nargs=2, required=True)
    ap.add_argument("--th3r", type=int, nargs=2, required=True)
    ap.add_argument("--window", type=int, nargs=2, default=None,
                    help="amino acid range panel A of figure 4 spans, "
                         "by default the range the substitutions cover")
    ap.add_argument("--country-a", default="Nigeria")
    ap.add_argument("--country-b", default="Cameroon")
    ap.add_argument("--tables-dir-name", default="10_table_for_plot")
    args = ap.parse_args()

    out = args.out_dir / args.tables_dir_name
    out.mkdir(parents=True, exist_ok=True)
    th2r, th3r = tuple(args.th2r), tuple(args.th3r)

    meta = S.load_metadata(args.metadata)
    clones_all = S.load_clones(args.out_dir, meta)
    dropped = clones_all[~clones_all.usable]
    clones = clones_all[clones_all.usable].copy()
    n_controls = S.count_controls(args.out_dir)
    print(f"[INFO] {len(meta)} samples, {len(clones_all)} clones, "
          f"{n_controls} controls")
    if len(dropped):
        print(f"[INFO] dropped for unresolved bases: {len(dropped)} "
              f"({', '.join(dropped.clone)})")

    amplicon = list(S.read_fasta(args.reference_dir / "pfcsp_amplicon.fasta").values())[0]
    off = args.callable_start - args.amplicon_start
    ref_nt = amplicon[off:off + (args.callable_end - args.callable_start + 1)]
    if len(ref_nt) != len(clones.nt.iloc[0]):
        raise RuntimeError(f"reference window is {len(ref_nt)} bp but the clones are "
                           f"{len(clones.nt.iloc[0])} bp")

    moi = S.load_moi(args.out_dir, meta)
    het = S.heterozygous_codons(args.out_dir, meta, args.codon_map,
                                args.callable_start, args.callable_end)
    with_cons = Counter(meta.at[s, "country"]
                        for s in S.samples_with_consensus(args.out_dir, meta))

    subs = pd.read_csv(args.out_dir / "08_protein" / "substitutions_consensus.tsv",
                       sep="\t")
    subs = subs[subs["sample"].isin(meta.index)].copy()
    subs["country"] = subs["sample"].map(meta["country"])

    window = (tuple(args.window) if args.window else
              (int(subs.aa_position.min()), int(subs.aa_position.max())))
    print(f"[INFO] figure 4 panel A spans amino acids {window[0]} to {window[1]}")

    cons = pd.read_csv(args.out_dir / "08_protein" / "epitope_status_consensus.tsv",
                       sep="\t")
    cons = cons[cons["sample"].isin(meta.index)].copy()
    cons["country"] = cons["sample"].map(meta["country"])
    cons = cons[cons.country == args.country_a].fillna("")

    eps = pd.read_csv(args.out_dir / "08_protein" / "epitope_status_haplotypes.tsv",
                      sep="\t")
    eps["country"] = eps["sample"].map(meta["country"])
    eps = eps[eps.country == args.country_a]

    sub_h = pd.read_csv(args.out_dir / "08_protein" / "substitutions_haplotypes.tsv",
                        sep="\t")
    sub_h["country"] = sub_h["sample"].map(meta["country"])
    sub_h = sub_h[sub_h.country == args.country_a]

    made = []
    made.append(("fig1_sample_counts.tsv", "Samples per country and controls",
                 write(out / "fig1_sample_counts.tsv",
                       table_fig1(meta, n_controls))))
    made.append(("fig2_moi.tsv", "MOI per sample by state, with mean and median",
                 write(out / "fig2_moi.tsv",
                       table_fig2(moi, args.country_a, args.country_b))))
    made.append(("fig3_diversity.tsv", "Nucleotide and haplotype diversity",
                 write(out / "fig3_diversity.tsv",
                       table_fig3(clones, args.country_a, args.country_b))))
    made.append(("fig4_substitutions.tsv",
                 "Heterozygosity, presence grid and prevalence by position",
                 write(out / "fig4_substitutions.tsv",
                       table_fig4(het, subs, meta["country"].value_counts().to_dict(),
                                  with_cons, window))))
    made.append(("fig7_vaccine_match.tsv", "Vaccine match, five panels",
                 write(out / "fig7_vaccine_match.tsv",
                       table_fig7(cons, eps, sub_h, th2r, th3r))))
    made.append(("fig8_structure.tsv", "FST by state, FST by country, PCoA with 3D7",
                 write(out / "fig8_structure.tsv",
                       table_fig8(clones, ref_nt, args.country_a, args.country_b))))
    made.append(("fig10_tree_annotations.tsv", "iTOL files and their load order",
                 write(out / "fig10_tree_annotations.tsv",
                       table_fig10(args.out_dir / "09_phylogeny" / "nigeria"))))

    dist = S.distance_to_reference(clones, ref_nt)
    dist.to_csv(out / "clone_distance_to_3d7.tsv", sep="\t", index=False)
    print(f"[OK] clone_distance_to_3d7.tsv ({len(dist)} rows)")
    made.append(("clone_distance_to_3d7.tsv",
                 "Nucleotide differences from 3D7 per clone", len(dist)))

    st = style_table(th2r, th3r)
    pd.DataFrame(st).to_csv(out / "plot_style.tsv", sep="\t", index=False)
    print(f"[OK] plot_style.tsv ({len(st)} rows)")

    # Figure 9 writes its own two tables, it owns the protein matrix.
    for name, what in (("fig9_protein_haplotypes.tsv", "one row per protein haplotype"),
                       ("fig9_protein_positions.tsv", "one row per variable position"),
                       ("fig9_protein_cells.tsv", "one row per drawn residue cell")):
        p = out / name
        made.append((name, what, len(pd.read_csv(p, sep="\t")) if p.exists() else 0))

    index = [{"figure": n.split("_")[0], "file": n, "rows": k, "describes": w,
              "written_by": "process/07_plot_tables.py" if not n.startswith("fig9")
                            else "process/csp_protein_matrix.py"}
             for n, w, k in made]
    index.append({"figure": "all", "file": "plot_style.tsv", "rows": len(st),
                  "describes": "every palette, declared once",
                  "written_by": "process/07_plot_tables.py"})
    pd.DataFrame(index).to_csv(out / "index.tsv", sep="\t", index=False)
    print(f"[OK] index.tsv ({len(index)} rows)")
    print(f"[OK] plot tables in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
