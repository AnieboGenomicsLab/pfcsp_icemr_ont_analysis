#!/usr/bin/env python3
"""
Phylogeny of the CSP C terminal window, our samples against MalariaGEN Pf8.

No aligner is run. Every upstream step writes the window in 3D7 coordinates and the
Pf8 genotypes are called against the same reference, so position i of one sequence is
already homologous to position i of every other. Realigning 299 bp of a hypervariable
gene would invent gaps the data does not support.

Two datasets per run, sharing the same Pf8 context and the same reference so that the
trees differ only by the sequences under study:

  nigeria   our Nigerian samples, Cameroon removed
  all       every field sample we sequenced

Sequences come from the resolved clones, not the sample consensus. A polyclonal
infection has an IUPAC consensus that no parasite carries, and most of our samples are
polyclonal. Pf8 cannot be phased from genotypes alone, so its side is restricted to
monoclonal samples instead, which is the matching choice rather than the same one.

Identical sequences are collapsed to a single tip. 299 bp cannot place a thousand near
identical tips, duplicated sequences carry zero length branches, and the count behind
each tip is kept in the metadata where a figure can read it.

A position carrying a Pf8 SNP record that did not pass Pf8 filters is blanked in every
sequence. Left alone it would be written as reference on the Pf8 side only, and our
samples would separate from Pf8 at a site where Pf8 was never callable.

Outputs under <out_dir>/09_phylogeny, with one directory per dataset:

  alignment_sequences.fasta        one record per clone and per Pf8 sample
  alignment_haplotypes.fasta       unique sequences, the tree input
  metadata_sequences.tsv           one row per record
  metadata_haplotypes.tsv          one row per tip, with its countries
  haplotype_country_matrix.tsv     tip by country counts
  colour_key.tsv                   every colour this step assigns
  <dataset>.treefile               IQ-TREE maximum likelihood tree
  <dataset>.rooted.nwk             the same tree rooted on 3D7
  <dataset>.support_summary.tsv    how much of the tree the alignment resolves
  <dataset>.popart_by_group.nex    PopART network input, sliced by country
  <dataset>.popart_by_region.nex   the same, sliced by region
  itol_colorstrip_origin.txt       ring 1, ours, shared with Pf8, or Pf8 only
  itol_piechart_by_region.txt      ring 2
  itol_piechart_by_country.txt     ring 3
  itol_ringlabel_origin.txt        names ring 1 on the tree
  itol_ringlabel_region.txt        names ring 2 on the tree
  itol_ringlabel_country.txt       names ring 3 on the tree
  itol_labelcolors_origin.txt      leaf label font colour, adds no ring

and, beside those directories:

  pf8_selected_samples.tsv         which Pf8 samples entered, and why
  control_check.tsv                what the NTC and 3D7 wells carry
  masked_positions.tsv             positions blanked in every sequence

Author: Khadim Gueye
"""

import argparse
import re
import subprocess
import sys
import unicodedata
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

# MalariaGEN population codes to a name a figure legend can carry.
REGION_OF_POPULATION = {
    "AF-W": "West Africa",
    "AF-C": "Central Africa",
    "AF-E": "East Africa",
    "AF-NE": "Northeast Africa",
    "AS-S-E": "South Asia (east)",
    "AS-S-FE": "South Asia (far east)",
    "AS-SE-W": "Southeast Asia (west)",
    "AS-SE-E": "Southeast Asia (east)",
    "OC-NG": "Oceania",
    "SA": "South America",
}

REFERENCE_NAME = "Pf3D7_ref"
UNAMBIGUOUS = set("ACGT")

OUR_COLOR = "#E7298A"
REFERENCE_COLOR = "#000000"
UNKNOWN_COLOR = "#BDBDBD"

# Country and region palettes.
#
# Shading each country inside a region hue put the twelve West African countries within
# a Delta E of 2.6 of one another in OKLab, where 15 is the floor at which a reader with
# full colour vision still separates two fills, so Ghana, Guinea and Benin arrived as
# one brown. These two lists were searched instead: candidates on several lightness
# bands at high chroma, then farthest point sampling for the set, then a 2-opt pass
# maximising the weakest neighbouring pair.
#
# Only neighbours touch, because a pie writes its slices in one fixed order, so the
# chain is what the search maximises and the order is part of the result. Our own
# samples and the reference are fixed colours sitting at their own alphabetical place,
# and the search held them in position: dropping the pink in afterwards had put it
# beside Nigeria at a Delta E of 3.4, which is the one pair the figure exists to
# separate.
#
#   country   worst adjacent pair Delta E 15.3 over all 33 slices, 17.3 normal vision
#   region    worst adjacent pair Delta E 19.0
#
# Colours are assigned to labels alphabetically, which is the order the pies and the
# legends use, so a change to the set of countries invalidates the ordering and the
# search has to be run again.
COUNTRY_PALETTE = [
    "#CE5A95", "#96B95A", "#C85CA4", "#CA6E00", "#903798", "#D95B52", "#9F2F7E",
    "#2EBBE8", "#007F67", "#00C1D4", "#A52B6F", "#ED86A7", "#AA295F", "#009EC6",
    "#75C075", "#009AD1", "#0BA458", "#0072AB", "#00C4BB", "#5B9D2C", "#C05FB2",
    "#AF283D", "#767BE7", "#DD9C42", "#9AA1FB", "#D65875", "#245DC2", "#AD274F",
    "#7C41AE", "#A34300", "#D35985",
]
REGION_PALETTE = [
    "#D65875", "#9AA1FB", "#A34300", "#767BE7", "#DD9C42",
    "#AA295F", "#ED86A7", "#245DC2", "#CA6E00", "#C05FB2",
]

# A coarser grouping, written to the metadata and to colour_key.tsv rather than drawn
# as a ring. Thirty one countries are too many for a legend that has to be read at a
# glance, and at this level the panel still answers the question the figure asks: where
# our sequences sit relative to Nigeria, to Cameroon and to the wider population. Its
# eight colours were validated on the same adjacency rule, worst adjacent pair Delta E
# 9.2 under deuteranopia and 22.9 under normal vision.
GROUP_SLOTS = [
    ("Our data", "#E7298A"),
    ("Nigeria", "#008300"),
    ("Cameroon", "#2A78D6"),
    ("Other West Africa", "#EB6834"),
    ("Central Africa", "#1BAF7A"),
    ("East Africa", "#6250D6"),
    ("Northeast Africa", "#EDA100"),
    ("Outside Africa", "#E34948"),
]
GROUP_ORDER = [name for name, _ in GROUP_SLOTS] + ["Reference"]
GROUP_COLORS = dict(GROUP_SLOTS)
GROUP_COLORS["Reference"] = REFERENCE_COLOR

OUTSIDE_AFRICA = {"South Asia (east)", "South Asia (far east)",
                  "Southeast Asia (west)", "Southeast Asia (east)",
                  "Oceania", "South America"}


def country_group(record):
    """Which of the eight groups a record falls in."""
    if record["source"] == "reference":
        return "Reference"
    if record["source"] == "our_data":
        return "Our data"
    country, region = record["country"], record["region"]
    if country in ("Nigeria", "Cameroon"):
        return country
    if region in OUTSIDE_AFRICA:
        return "Outside Africa"
    return "Other West Africa" if region == "West Africa" else region


def read_fasta(path):
    """Name to sequence, upper case, first word of the header as the name."""
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


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def write_tsv(path, fields, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")


def ambiguous_fraction(seq):
    if not seq:
        return 1.0
    return sum(1 for c in seq if c not in UNAMBIGUOUS) / len(seq)


def reference_window(genome, chrom, start, end):
    """The 3D7 bases of the window, as samtools reports them, 1 based inclusive."""
    res = subprocess.run(["samtools", "faidx", str(genome), f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True, check=True)
    return "".join(res.stdout.splitlines()[1:]).upper()


# --- our sequences -------------------------------------------------------------

def classify_well(sample_id, ntc_pattern, positive_pattern):
    if re.search(ntc_pattern, sample_id):
        return "negative_control"
    if re.search(positive_pattern, sample_id):
        return "positive_control"
    return "field"


def country_of(sample_id, cameroon_prefix):
    """Our plates carry two countries and the identifier says which."""
    return "Cameroon" if sample_id.startswith(cameroon_prefix) else "Nigeria"


def load_our_sequences(fasta, ntc_pattern, positive_pattern, cameroon_prefix):
    """One record per sequence, with the sample it came from.

    A haplotype is named <sample>_H<n>, so the clone suffix is stripped only when it
    looks like one. Splitting on the last underscore unconditionally renamed the
    control well c1_3d7 to the sample c1, which matches neither control pattern, so
    every control entered the consensus datasets as a field sample.
    """
    rows = []
    for name, seq in read_fasta(fasta).items():
        m = re.match(r"^(.+)_(H\d+)$", name)
        sample, clone = (m.group(1), m.group(2)) if m else (name, "")
        rows.append({
            "sequence_id": name,
            "sample_id": sample,
            "clone": clone,
            "well_type": classify_well(sample, ntc_pattern, positive_pattern),
            "country": country_of(sample, cameroon_prefix),
            "sequence": seq,
        })
    return rows


# --- Pf8 -----------------------------------------------------------------------

def load_pf8_metadata(path):
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        need = ["Sample", "Country", "Year", "Population", "QC pass", "Study"]
        missing = [n for n in need if n not in idx]
        if missing:
            sys.exit(f"[ERROR] Pf8 metadata is missing columns: {', '.join(missing)}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            rows[f[idx["Sample"]]] = {
                "country": f[idx["Country"]],
                "year": f[idx["Year"]],
                "population": f[idx["Population"]],
                "qc_pass": f[idx["QC pass"]],
                "study": f[idx["Study"]],
            }
    return rows


def load_pf8_fws(path):
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            try:
                out[f[0]] = float(f[1])
            except ValueError:
                continue
    return out


# How many candidates are genotyped per country before the cap is applied. Picking by
# Fws alone brought in Nigerian samples that were 23 percent missing across the window,
# which is a row of N in the alignment rather than a sequence. Candidates are built
# first, then ranked on how complete they came out, so the cap is filled with sequences
# that exist rather than with the first names off the Fws list.
PF8_OVERSAMPLE = 5


def pf8_candidates(meta, fws, vcf_samples, focus, context, max_focus, max_context,
                   min_fws):
    """Pf8 samples worth genotyping, before completeness is known.

    QC pass and a Fws at or above the floor, because an unphased mixed infection
    becomes a consensus of two clones and would sit on a branch that belongs to no
    parasite.
    """
    in_vcf = set(vcf_samples)
    wanted = {c: max_focus for c in focus}
    wanted.update({c: max_context for c in context if c not in wanted})
    by_country = defaultdict(list)
    for sample, m in meta.items():
        if sample not in in_vcf:
            continue
        if m["qc_pass"] != "True" or m["country"] not in wanted:
            continue
        f = fws.get(sample)
        if f is None or f < min_fws:
            continue
        by_country[m["country"]].append((f, sample))
    out = []
    for country, cap in wanted.items():
        ranked = sorted(by_country.get(country, []), reverse=True)
        for f, sample in ranked[:cap * PF8_OVERSAMPLE]:
            m = meta[sample]
            out.append({
                "sample_id": sample, "country": country, "year": m["year"],
                "population": m["population"], "study": m["study"],
                "fws": round(f, 4), "n_available": len(ranked), "cap": cap,
            })
    return out


def pf8_keep_most_complete(candidates, max_ambiguous):
    """Fill each country cap with its most complete sequences."""
    by_country = defaultdict(list)
    for c in candidates:
        by_country[c["country"]].append(c)
    for country, rows in by_country.items():
        rows.sort(key=lambda r: (r["ambiguous_fraction"], -r["fws"]))
        cap = rows[0]["cap"]
        taken = 0
        for r in rows:
            if taken < cap and r["ambiguous_fraction"] <= max_ambiguous:
                r["used"] = "yes"
                taken += 1
            else:
                r["used"] = ("no, above the ambiguity floor"
                             if r["ambiguous_fraction"] > max_ambiguous
                             else "no, country cap already filled")
    return candidates


def pf8_window_sequences(vcf, chrom, start, end, samples, ref_seq, mask_non_pass):
    """Rebuild each Pf8 sample over the window from its genotypes.

    Only sites Pf8 itself passed are applied. A heterozygous call, a missing call and
    a spanning deletion all become N: the first cannot be phased from a genotype, and
    the other two are an absence of evidence, not a reference base.
    """
    sample_list = [s["sample_id"] for s in samples]
    region = f"{chrom}:{start}-{end}"
    proc = subprocess.run(
        ["bcftools", "query", "-r", region, "-s", ",".join(sample_list),
         "-f", "%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n", str(vcf)],
        capture_output=True, text=True, check=True)

    seqs = {s: list(ref_seq) for s in sample_list}
    masked, applied = set(), 0
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        pos, ref, alt, filt = int(f[0]), f[1], f[2], f[3]
        if not (start <= pos <= end):
            continue
        alts = alt.split(",")
        real = [a for a in alts if a != "*"]
        # SNP record: one reference base and every true alternative one base too.
        if len(ref) != 1 or not real or any(len(a) != 1 for a in real):
            continue
        col = pos - start
        if filt != "PASS":
            if mask_non_pass:
                masked.add(pos)
            continue
        applied += 1
        alleles = [ref] + alts
        for sample, gt in zip(sample_list, f[4:]):
            call = gt.replace("|", "/").split("/")
            if len(call) < 1 or call[0] == "." or len(set(call)) != 1:
                seqs[sample][col] = "N"
                continue
            try:
                a = alleles[int(call[0])]
            except (ValueError, IndexError):
                seqs[sample][col] = "N"
                continue
            seqs[sample][col] = "N" if a == "*" else a
    return {s: "".join(v) for s, v in seqs.items()}, sorted(masked), applied


# --- assembly ------------------------------------------------------------------

def apply_mask(seq, columns):
    if not columns:
        return seq
    out = list(seq)
    for c in columns:
        out[c] = "N"
    return "".join(out)


def collapse(records):
    """Identical sequences become one tip. The reference names the tip it falls in."""
    groups = OrderedDict()
    for r in records:
        groups.setdefault(r["sequence"], []).append(r)
    tips, n = [], 0
    for seq, members in groups.items():
        has_ref = any(m["source"] == "reference" for m in members)
        if has_ref:
            name = REFERENCE_NAME
        else:
            n += 1
            name = f"H{n:03d}"
        tips.append({"tip": name, "sequence": seq, "members": members})
    return tips


def build_palette(records):
    """A colour per country and a colour per region, keyed by the label a file carries.

    Our own samples and the reference are fixed, the rest take the searched palettes in
    alphabetical order, which is the order the slices are written in.
    """
    fixed = {"Our data": OUR_COLOR, REFERENCE_NAME: REFERENCE_COLOR,
             "Reference": REFERENCE_COLOR}
    group, region = dict(fixed), dict(fixed)
    for key, field, palette in (("group", "group", COUNTRY_PALETTE),
                                ("region", "region", REGION_PALETTE)):
        labels = sorted({r[field] for r in records if r[field] and r[field] not in fixed})
        target = group if key == "group" else region
        for i, label in enumerate(labels):
            if i >= len(palette):
                print(f"[WARN] {len(labels)} {key} labels but {len(palette)} colours, "
                      f"{label} and the rest fall back to grey")
                target[label] = UNKNOWN_COLOR
            else:
                target[label] = palette[i]
    return group, region


def write_itol_piechart(path, tips, trait_key, colours, title):
    """One pie per tip, sliced by country or by region, sized by how many sequences
    the tip holds. This is the shape a collapsed haplotype tree needs: a tip is not one
    country, it is every country that carries that sequence."""
    present = {(m.get(trait_key) or "unknown") for t in tips for m in t["members"]}
    values = sorted(present)
    cols = [colours.get(v, UNKNOWN_COLOR) for v in values]
    with open(path, "w") as out:
        out.write("DATASET_PIECHART\nSEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{title}\n")
        out.write(f"COLOR\t{OUR_COLOR}\n")
        out.write("FIELD_COLORS\t" + "\t".join(cols) + "\n")
        out.write("FIELD_LABELS\t" + "\t".join(v.replace(" ", "_") for v in values) + "\n")
        out.write(f"LEGEND_TITLE\t{title}\n")
        out.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in values) + "\n")
        out.write("LEGEND_COLORS\t" + "\t".join(cols) + "\n")
        out.write("LEGEND_LABELS\t" + "\t".join(values) + "\n")
        out.write("DATA\n")
        for t in tips:
            counts = Counter((m.get(trait_key) or "unknown") for m in t["members"])
            total = sum(counts.values())
            # Area, not radius, carries the count, so a tip of 36 does not read as six
            # times a tip of 1.
            radius = max(4, min(30, round(4 + 5 * total ** 0.5)))
            out.write(f"{t['tip']}\t-1\t{radius}\t"
                      + "\t".join(str(counts.get(v, 0)) for v in values) + "\n")


# Tip label text, coloured by the same three categories as the origin strip. The ring
# says where a haplotype was seen, but the tip names are what a reader follows when
# tracing one, and a name in plain black carries none of that.
ORIGIN_COLORS = {
    "Our data only": OUR_COLOR,
    "Shared with Pf8": "#7570B3",
    "Pf8 only": "#9E9E9E",
    "Reference": REFERENCE_COLOR,
}


def origin_of(tip):
    ours = any(m["source"] == "our_data" for m in tip["members"])
    pf8 = any(m["source"] == "Pf8" for m in tip["members"])
    if tip["tip"] == REFERENCE_NAME:
        return "Reference"
    if ours and pf8:
        return "Shared with Pf8"
    return "Our data only" if ours else "Pf8 only"


def write_itol_label_colors(path, tips):
    """Leaf label font colour and weight, by origin.

    Our own haplotypes are set bold as well as coloured. Colour alone separates them
    on screen and stops separating them in a printed or a colour blind reading, and
    these are the tips the figure exists to point at.
    """
    with open(path, "w") as out:
        out.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for t in tips:
            kind = origin_of(t)
            style = "bold" if kind in ("Our data only", "Shared with Pf8") else "normal"
            out.write(f"{t['tip']}\tlabel\t{ORIGIN_COLORS[kind]}\t{style}\t1\n")


def write_itol_ring_label(path, tips, label):
    """A label printed beside a ring on the tree itself.

    DATASET_PIECHART has no SHOW_LABELS, so a pie ring cannot name itself and three
    concentric rings reach the reader unlabelled. DATASET_COLORSTRIP does have it, so
    this writes a one pixel white strip whose only job is to carry the label. Load it
    immediately before the ring it names.
    """
    with open(path, "w") as out:
        out.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{label}\n")
        out.write("COLOR\t#FFFFFF\n")
        out.write("COLOR_BRANCHES\t0\n")
        out.write("STRIP_WIDTH\t1\nMARGIN\t40\nBORDER_WIDTH\t0\n")
        out.write("SHOW_LABELS\t1\nLABEL_ROTATION\t0\nLABEL_SHIFT\t0\n")
        out.write("LABEL_ALIGN_TO_TREE\t0\n")
        out.write("DATA\n")
        for t in tips:
            out.write(f"{t['tip']}\t#FFFFFF\n")


def write_itol_colorstrip(path, tips):
    """Where our sequences are: a tip we carry, a tip we share with Pf8, or neither."""
    kinds = ORIGIN_COLORS
    with open(path, "w") as out:
        out.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        out.write("DATASET_LABEL\tOrigin\n")
        out.write(f"COLOR\t{OUR_COLOR}\n")
        out.write("COLOR_BRANCHES\t0\nSTRIP_WIDTH\t28\n")
        out.write("SHOW_LABELS\t0\n")
        out.write("LEGEND_TITLE\tOrigin of the haplotype\n")
        out.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in kinds) + "\n")
        out.write("LEGEND_COLORS\t" + "\t".join(kinds.values()) + "\n")
        out.write("LEGEND_LABELS\t" + "\t".join(kinds) + "\n")
        out.write("DATA\n")
        for t in tips:
            kind = origin_of(t)
            out.write(f"{t['tip']}\t{kinds[kind]}\t{kind}\n")


def write_colour_key(path, group_colours, region_colours):
    """The same colours as a table, for a figure drawn anywhere other than iTOL."""
    rows = [{"level": "grouped", "label": k, "colour": GROUP_COLORS[k]}
            for k in GROUP_ORDER]
    rows += [{"level": "country", "label": k, "colour": v}
             for k, v in sorted(group_colours.items())]
    rows += [{"level": "region", "label": k, "colour": v}
             for k, v in sorted(region_colours.items())]
    write_tsv(path, ["level", "label", "colour"], rows)


def nexus_label(name):
    """PopART reads a trait label up to the first space, so spaces become underscores
    and punctuation is dropped rather than left to split a country in two.

    Accents are folded to their base letter first. Dropping them outright turned
    Cote d'Ivoire into Cte_dIvoire on the legend of the figure.
    """
    folded = unicodedata.normalize("NFKD", name)
    folded = "".join(c for c in folded if not unicodedata.combining(c))
    return re.sub(r"[^A-Za-z0-9_]", "", folded.replace(" ", "_").replace("-", "_"))


def write_popart_nexus(path, tips, trait_key, missing="?"):
    """The haplotype alignment with a TRAITS block, which PopART turns into a network
    whose node areas are the counts and whose slices are the groups.

    A network rather than only a tree because this window carries few informative
    sites. A bifurcating tree has to commit to an order of descent that 31 informative
    sites cannot support, while a network shows the same sequences and leaves the
    unresolved part visibly unresolved.
    """
    names = [t["tip"] for t in tips]
    seqs = {t["tip"]: t["sequence"].replace("N", missing) for t in tips}
    traits = []
    for t in tips:
        for m in t["members"]:
            v = m.get(trait_key) or "unknown"
            if v not in traits:
                traits.append(v)
    traits.sort()
    labels = [nexus_label(v) for v in traits]
    nchar = len(next(iter(seqs.values()))) if seqs else 0
    with open(path, "w") as out:
        out.write("#NEXUS\n\n")
        out.write("BEGIN TAXA;\n")
        out.write(f"DIMENSIONS NTAX={len(names)};\n")
        out.write("TAXLABELS\n")
        for name in names:
            out.write(f"{name}\n")
        out.write(";\nEND;\n\n")
        out.write("BEGIN CHARACTERS;\n")
        out.write(f"DIMENSIONS NCHAR={nchar};\n")
        out.write("FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. ;\n")
        out.write("MATRIX\n\n")
        for name in names:
            out.write(f"{name} {seqs[name]}\n")
        out.write(";\nEND;\n\n")
        out.write("BEGIN TRAITS;\n")
        out.write(f"  Dimensions NTRAITS={len(traits)};\n")
        out.write("  Format labels=yes missing=? separator=Comma;\n")
        out.write("  TraitLabels " + " ".join(labels) + ";\n")
        out.write("  Matrix\n")
        for t in tips:
            counts = Counter((m.get(trait_key) or "unknown") for m in t["members"])
            out.write(f"{t['tip']} " + ",".join(str(counts.get(v, 0)) for v in traits) + "\n")
        out.write(";\nEND;\n")


def support_summary(treefile):
    """How much of the tree the alignment actually resolves.

    IQ-TREE writes SH-aLRT/UFBoot on each internal branch. A node is reported as
    resolved only when it carries both, which is the pair the IQ-TREE manual asks for:
    either one alone accepts nodes the other rejects.
    """
    text = Path(treefile).read_text()
    pairs = re.findall(r"\)([0-9.]+)/([0-9.]+):", text)
    if not pairs:
        return {}
    alrt = [float(a) for a, _ in pairs]
    boot = [float(b) for _, b in pairs]
    both = sum(1 for a, b in zip(alrt, boot) if a >= 80 and b >= 95)
    return {
        "internal_nodes": len(pairs),
        "ufboot_ge_95": sum(1 for b in boot if b >= 95),
        "shalrt_ge_80": sum(1 for a in alrt if a >= 80),
        "resolved_both": both,
        "resolved_pct": round(100.0 * both / len(pairs), 1),
    }


def run_iqtree(alignment, prefix, outgroup, bootstrap, alrt, threads, model, log):
    cmd = ["iqtree", "-s", str(alignment), "-m", model, "--prefix", str(prefix),
           "-T", str(threads), "-o", outgroup, "--redo", "--quiet"]
    if bootstrap:
        cmd += ["-B", str(bootstrap)]
    if alrt:
        cmd += ["--alrt", str(alrt)]
    print("[INFO] " + " ".join(cmd))
    with open(log, "w") as fh:
        res = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return res.returncode


def root_tree(treefile, out_path, outgroup):
    """Root on the reference and write the result beside the IQ-TREE output.

    IQ-TREE places the outgroup first but writes the tree unrooted. Rerooting here
    means the file a plotting tool opens already hangs from 3D7 and no one has to
    remember to do it in the viewer.
    """
    try:
        from Bio import Phylo
    except ImportError:
        print("[WARN] biopython not available, the rooted copy was not written")
        return False
    tree = Phylo.read(str(treefile), "newick")
    names = [t.name for t in tree.get_terminals()]
    if outgroup not in names:
        print(f"[WARN] {outgroup} is not a tip of {treefile}, left unrooted")
        return False
    tree.root_with_outgroup(outgroup)
    Phylo.write(tree, str(out_path), "newick")
    return True


def main():
    ap = argparse.ArgumentParser(description="CSP phylogeny against MalariaGEN Pf8.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--callable-start", type=int, required=True)
    ap.add_argument("--callable-end", type=int, required=True)
    ap.add_argument("--sequence-source", choices=["haplotypes", "consensus"],
                    default="haplotypes")
    ap.add_argument("--max-ambiguous", type=float, default=0.05,
                    help="drop a sequence with more than this fraction of non ACGT")
    ap.add_argument("--ntc-pattern", default="NTC")
    ap.add_argument("--positive-pattern", default="3D7|3d7")
    ap.add_argument("--cameroon-prefix", default="CMR")
    ap.add_argument("--our-year", default="2025")
    ap.add_argument("--pf8-vcf", type=Path)
    ap.add_argument("--pf8-metadata", type=Path)
    ap.add_argument("--pf8-fws", type=Path)
    ap.add_argument("--pf8-focus-countries", default="Nigeria,Cameroon")
    ap.add_argument("--pf8-context-countries", default="")
    ap.add_argument("--pf8-max-focus", type=int, default=80)
    ap.add_argument("--pf8-max-context", type=int, default=20)
    ap.add_argument("--pf8-min-fws", type=float, default=0.95)
    ap.add_argument("--mask-pf8-non-pass", type=int, default=1)
    ap.add_argument("--model", default="MFP")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--alrt", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--datasets", default="nigeria,all")
    ap.add_argument("--skip-tree", action="store_true",
                    help="rebuild the alignments, metadata and Nexus without IQ-TREE")
    args = ap.parse_args()

    phylo = args.out_dir / "09_phylogeny"
    phylo.mkdir(parents=True, exist_ok=True)
    start, end = args.callable_start, args.callable_end

    ref_seq = reference_window(args.genome, args.chrom, start, end)
    print(f"[INFO] reference window {args.chrom}:{start}-{end}, {len(ref_seq)} bp")

    source = ("07_haplotypes/haplotypes_nt.fasta" if args.sequence_source == "haplotypes"
              else "06_consensus/consensus_nt.fasta")
    src_path = args.out_dir / source
    if not src_path.exists():
        sys.exit(f"[ERROR] sequences not found: {src_path}")
    ours = load_our_sequences(src_path, args.ntc_pattern, args.positive_pattern,
                              args.cameroon_prefix)
    bad_len = {r["sequence_id"] for r in ours if len(r["sequence"]) != len(ref_seq)}
    if bad_len:
        sys.exit(f"[ERROR] {len(bad_len)} sequences are not {len(ref_seq)} bp, "
                 f"the window does not match: {sorted(bad_len)[:5]}")

    # What the control wells carry, written whether or not they enter a tree.
    controls = [r for r in ours if r["well_type"] != "field"]
    rows = []
    for r in controls:
        diff = [i + start for i, (a, b) in enumerate(zip(ref_seq, r["sequence"]))
                if b in UNAMBIGUOUS and a != b]
        rows.append({
            "sequence_id": r["sequence_id"], "well_type": r["well_type"],
            "n_differences_vs_3D7": len(diff),
            "ambiguous_fraction": round(ambiguous_fraction(r["sequence"]), 4),
            "positions": ";".join(str(d) for d in diff),
        })
    write_tsv(phylo / "control_check.tsv",
              ["sequence_id", "well_type", "n_differences_vs_3D7",
               "ambiguous_fraction", "positions"], rows)
    print(f"[OK] control_check.tsv ({len(rows)} control sequences)")

    field = [r for r in ours if r["well_type"] == "field"]
    kept, dropped = [], []
    for r in field:
        (dropped if ambiguous_fraction(r["sequence"]) > args.max_ambiguous
         else kept).append(r)
    print(f"[INFO] our sequences: {len(field)} clones, {len(kept)} kept, "
          f"{len(dropped)} above {args.max_ambiguous:.0%} ambiguous")
    if dropped:
        print("[INFO] dropped: " + ", ".join(r["sequence_id"] for r in dropped[:10]))
    for r in kept:
        r.update({"source": "our_data", "year": args.our_year, "population": "",
                  "fws": "", "study": ""})

    # Pf8 context.
    pf8_records, pf8_chosen, masked_pos = [], [], []
    if args.pf8_vcf and args.pf8_metadata and args.pf8_fws:
        meta = load_pf8_metadata(args.pf8_metadata)
        fws = load_pf8_fws(args.pf8_fws)
        listed = subprocess.run(["bcftools", "query", "-l", str(args.pf8_vcf)],
                                capture_output=True, text=True, check=True)
        vcf_samples = listed.stdout.split()
        focus = [c.strip() for c in args.pf8_focus_countries.split(",") if c.strip()]
        context = [c.strip() for c in args.pf8_context_countries.split(",") if c.strip()]
        pf8_chosen = pf8_candidates(meta, fws, vcf_samples, focus, context,
                                    args.pf8_max_focus, args.pf8_max_context,
                                    args.pf8_min_fws)
        print(f"[INFO] Pf8: {len(pf8_chosen)} candidates over "
              f"{len({c['country'] for c in pf8_chosen})} countries")
        if pf8_chosen:
            seqs, masked_pos, applied = pf8_window_sequences(
                args.pf8_vcf, args.chrom, start, end, pf8_chosen, ref_seq,
                bool(args.mask_pf8_non_pass))
            print(f"[INFO] Pf8: {applied} passing SNP sites applied, "
                  f"{len(masked_pos)} sites blanked in every sequence")
            for c in pf8_chosen:
                c["ambiguous_fraction"] = round(
                    ambiguous_fraction(seqs[c["sample_id"]]), 4)
            pf8_keep_most_complete(pf8_chosen, args.max_ambiguous)
            for c in pf8_chosen:
                if c["used"] != "yes":
                    continue
                pf8_records.append({
                    "sequence_id": c["sample_id"], "sample_id": c["sample_id"],
                    "clone": "", "well_type": "field", "country": c["country"],
                    "sequence": seqs[c["sample_id"]], "source": "Pf8",
                    "year": c["year"], "population": c["population"],
                    "fws": c["fws"], "study": c["study"],
                })
            write_tsv(phylo / "pf8_selected_samples.tsv",
                      ["sample_id", "country", "population", "year", "study", "fws",
                       "n_available", "cap", "ambiguous_fraction", "used"], pf8_chosen)
            per_country = Counter(c["country"] for c in pf8_chosen if c["used"] == "yes")
            print(f"[OK] pf8_selected_samples.tsv ({len(pf8_records)} usable: "
                  + ", ".join(f"{k} {v}" for k, v in sorted(per_country.items())) + ")")
    else:
        print("[WARN] no Pf8 input given, the trees carry our samples and 3D7 only")

    write_tsv(phylo / "masked_positions.tsv", ["position", "reason"],
              [{"position": p, "reason": "Pf8 SNP record did not pass Pf8 filters"}
               for p in masked_pos])
    mask_cols = [p - start for p in masked_pos]

    reference = {
        "sequence_id": REFERENCE_NAME, "sample_id": REFERENCE_NAME, "clone": "",
        "well_type": "reference", "country": "Reference", "sequence": ref_seq,
        "source": "reference", "year": "", "population": "", "fws": "", "study": "",
    }

    for ds in [d.strip() for d in args.datasets.split(",") if d.strip()]:
        if ds == "nigeria":
            mine = [r for r in kept if r["country"] == "Nigeria"]
        elif ds == "all":
            mine = list(kept)
        else:
            print(f"[WARN] unknown dataset {ds}, skipped")
            continue
        records = [dict(r) for r in mine + pf8_records + [reference]]
        for r in records:
            r["sequence"] = apply_mask(r["sequence"], mask_cols)
            r["region"] = REGION_OF_POPULATION.get(
                r["population"], "West Africa" if r["source"] == "our_data"
                else ("Reference" if r["source"] == "reference" else ""))
            r["group"] = ("Our data" if r["source"] == "our_data"
                          else (REFERENCE_NAME if r["source"] == "reference"
                                else r["country"]))
            r["group_detail"] = (f"Our data - {r['country']}"
                                 if r["source"] == "our_data" else r["group"])
        for r in records:
            r["group8"] = country_group(r)

        d = phylo / ds
        d.mkdir(parents=True, exist_ok=True)
        write_fasta(d / "alignment_sequences.fasta",
                    [(r["sequence_id"], r["sequence"]) for r in records])
        write_tsv(d / "metadata_sequences.tsv",
                  ["sequence_id", "sample_id", "clone", "source", "group8", "group",
                   "group_detail", "country", "region", "population", "year", "fws",
                   "study"], records)

        tips = collapse(records)
        write_fasta(d / "alignment_haplotypes.fasta",
                    [(t["tip"], t["sequence"]) for t in tips])

        countries = sorted({r["group"] for r in records})
        tip_rows, matrix_rows = [], []
        for t in tips:
            m = t["members"]
            counts = Counter(x["group"] for x in m)
            tip_rows.append({
                "tip": t["tip"], "n_sequences": len(m),
                "n_our_clones": sum(1 for x in m if x["source"] == "our_data"),
                "n_our_samples": len({x["sample_id"] for x in m
                                      if x["source"] == "our_data"}),
                "n_pf8_samples": sum(1 for x in m if x["source"] == "Pf8"),
                "contains_reference": "yes" if t["tip"] == REFERENCE_NAME else "no",
                "groups": ";".join(f"{k}:{v}" for k, v in counts.most_common()),
                "regions": ";".join(sorted({x["region"] for x in m if x["region"]})),
                "members": ";".join(x["sequence_id"] for x in m[:20]),
            })
            row = {"tip": t["tip"]}
            row.update({c: counts.get(c, 0) for c in countries})
            matrix_rows.append(row)
        write_tsv(d / "metadata_haplotypes.tsv",
                  ["tip", "n_sequences", "n_our_clones", "n_our_samples",
                   "n_pf8_samples", "contains_reference", "groups", "regions",
                   "members"], tip_rows)
        write_tsv(d / "haplotype_country_matrix.tsv", ["tip"] + countries, matrix_rows)
        write_popart_nexus(d / f"{ds}.popart_by_group.nex", tips, "group")
        write_popart_nexus(d / f"{ds}.popart_by_region.nex", tips, "region")
        group_colours, region_colours = build_palette(records)
        write_itol_piechart(d / "itol_piechart_by_country.txt", tips, "group",
                            group_colours, "Country")
        write_itol_piechart(d / "itol_piechart_by_region.txt", tips, "region",
                            region_colours, "Region")
        write_itol_colorstrip(d / "itol_colorstrip_origin.txt", tips)
        write_itol_ring_label(d / "itol_ringlabel_origin.txt", tips, "Origin")
        write_itol_ring_label(d / "itol_ringlabel_region.txt", tips, "Region")
        write_itol_ring_label(d / "itol_ringlabel_country.txt", tips, "Country")
        write_itol_label_colors(d / "itol_labelcolors_origin.txt", tips)
        write_colour_key(d / "colour_key.tsv", group_colours, region_colours)

        print(f"[OK] {ds}: {len(records)} sequences, {len(tips)} unique haplotypes "
              f"({len(mine)} ours, {len(pf8_records)} Pf8, 1 reference)")

        # --skip-tree rebuilds the alignments, the metadata and the Nexus files
        # without spending another run of IQ-TREE. An existing tree is still
        # summarised, so the flag adds files rather than dropping them.
        if args.skip_tree:
            if not (d / f"{ds}.treefile").exists():
                print(f"[INFO] {ds}: no tree and IQ-TREE skipped on request")
                continue
            print(f"[INFO] {ds}: IQ-TREE skipped, the existing tree is summarised")
        else:
            if len(tips) < 4:
                print(f"[WARN] {ds} has {len(tips)} tips, too few for a tree")
                continue
            code = run_iqtree(d / "alignment_haplotypes.fasta", d / ds, REFERENCE_NAME,
                              args.bootstrap, args.alrt, args.threads, args.model,
                              d / f"{ds}.iqtree.log")
            if code != 0 or not (d / f"{ds}.treefile").exists():
                print(f"[ERROR] IQ-TREE failed for {ds}, see {d / ds}.log")
                continue
        root_tree(d / f"{ds}.treefile", d / f"{ds}.rooted.nwk", REFERENCE_NAME)
        model_line, informative = "", ""
        for line in (d / f"{ds}.iqtree").read_text().splitlines():
            if line.startswith("Model of substitution:"):
                model_line = line.split(":", 1)[1].strip()
            elif line.startswith("Number of parsimony informative sites"):
                informative = line.split(":", 1)[1].strip()
        sup = support_summary(d / f"{ds}.treefile")
        write_tsv(d / f"{ds}.support_summary.tsv",
                  ["dataset", "model", "parsimony_informative_sites", "internal_nodes",
                   "ufboot_ge_95", "shalrt_ge_80", "resolved_both", "resolved_pct"],
                  [{"dataset": ds, "model": model_line,
                    "parsimony_informative_sites": informative, **sup}])
        print(f"[OK] {ds}: tree written, model {model_line or 'unknown'}, "
              f"{informative} informative sites, "
              f"{sup.get('resolved_both', 0)}/{sup.get('internal_nodes', 0)} nodes "
              f"resolved ({sup.get('resolved_pct', 0)}%)")

    print(f"[OK] phylogeny in {phylo}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
