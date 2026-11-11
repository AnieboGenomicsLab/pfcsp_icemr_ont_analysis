#!/usr/bin/env python3
"""
Build the cohort metadata table the pipeline reads through METADATA.

One row per sample, with the Nigerian state the sample was collected in. The state is
not in the run folder, it lives in the collection tables of the other studies, so it is
looked up by sample identifier in the sources given with --source, first match wins.

Cameroonian samples carry the country in the state column and nothing finer, because no
state was recorded for them. Controls are left out, they are not field samples.

Six ICEMR identifiers are absent from every source. They are filled with the state that
every other identifier of the same series carries, and the state_source column says so
on each of those rows, so an inferred value is never mistaken for a read one.

Author: Khadim Gueye
"""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

ID_KEYS = ("sample_id", "sample id", "sample", "id")
STATE_KEYS = ("state", "region")

# Spelling kept when the same state is written differently by two sources.
CANONICAL = {"kogi": "Kogi", "federal capital territory": "FCT", "fct": "FCT"}


def read_source(path):
    """sample to state, from any table holding an identifier column and a state column."""
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = csv.DictReader(fh, delimiter=sep)
        cols = {(c or "").strip().lower(): c for c in (rows.fieldnames or [])}
        key = next((cols[c] for c in ID_KEYS if c in cols), None)
        col = next((cols[c] for c in STATE_KEYS if c in cols), None)
        if not key or not col:
            sys.exit(f"[ERROR] {path} has no identifier column or no state column")
        out = {}
        for row in rows:
            name = (row.get(key) or "").strip()
            state = (row.get(col) or "").strip()
            if name and state:
                out.setdefault(name, CANONICAL.get(state.lower(), state))
    return out


def sample_ids(input_dir, pattern):
    return sorted({re.sub(pattern, "", p.name) for p in input_dir.glob("*.fastq.gz")})


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", type=Path, required=True,
                    help="folder of <sample>.fastq.gz, the cohort to describe")
    ap.add_argument("--source", type=Path, action="append", required=True,
                    help="table with an identifier column and a state column, repeatable, "
                         "first match wins")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cameroon-prefix", default="CMR")
    ap.add_argument("--cameroon-label", default="Cameroon")
    ap.add_argument("--country-a", default="Nigeria")
    ap.add_argument("--control-pattern", default=r"^c\d+_(3[dD]7|NTC)$",
                    help="identifiers matching this are controls and are left out")
    ap.add_argument("--fallback-series", default=r"^NGA-\d+-",
                    help="identifiers matching this take the state of their own series "
                         "when no source holds them")
    args = ap.parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"[ERROR] input folder not found: {args.input_dir}")

    lookup = {}
    for src in args.source:
        if not src.exists():
            print(f"[WARN] source not found, skipped: {src}")
            continue
        found = read_source(src)
        new = {k: (v, src.name) for k, v in found.items() if k not in lookup}
        lookup.update(new)
        print(f"[OK] {src.name}: {len(found)} entries, {len(new)} new")
    if not lookup:
        sys.exit("[ERROR] no usable source, nothing to look up")

    ids = sample_ids(args.input_dir, r"\.fastq\.gz$")
    if not ids:
        sys.exit(f"[ERROR] no fastq found in {args.input_dir}")
    control = re.compile(args.control_pattern)
    series = re.compile(args.fallback_series)

    # The state every member of a series carries, used only when a member is unlisted.
    by_series = {}
    for name, (state, _) in lookup.items():
        m = series.match(name)
        if m:
            by_series.setdefault(m.group(0), Counter())[state] += 1

    # Cameroon is held back and written after Nigeria, so the cohort of the figure
    # reads first and the country without a state sits at the end.
    rows, tail, skipped, inferred, unknown = [], [], [], [], []
    for name in ids:
        if control.match(name):
            skipped.append(name)
            continue
        if name.startswith(args.cameroon_prefix):
            tail.append((name, args.cameroon_label, args.cameroon_label,
                         "country only, no state recorded"))
            continue
        if name in lookup:
            state, origin = lookup[name]
            rows.append((name, args.country_a, state, origin))
            continue
        base = re.sub(r"[a-z]$", "", name)
        if base != name and base in lookup:
            state, origin = lookup[base]
            rows.append((name, args.country_a, state,
                         f"{origin}, replicate of {base}"))
            continue
        m = series.match(name)
        counts = by_series.get(m.group(0)) if m else None
        if counts and len(counts) == 1:
            state = next(iter(counts))
            rows.append((name, args.country_a, state,
                         f"inferred, every {m.group(0)} sample is {state}"))
            inferred.append(name)
            continue
        rows.append((name, args.country_a, "", "not found"))
        unknown.append(name)

    rows += tail
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        fh.write("sample_id\tcountry\tstate\tstate_source\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    states = Counter(s for _, c, s, _ in rows if c == args.country_a and s)
    print(f"[OK] {args.out} ({len(rows)} samples, {len(skipped)} controls left out)")
    print(f"[INFO] {args.country_a}: {sum(states.values())} samples, "
          f"{len(states)} states")
    for state, n in sorted(states.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"       {state:<12}{n}")
    if inferred:
        print(f"[WARN] state inferred from the series for {len(inferred)} sample(s): "
              f"{', '.join(inferred)}")
    if unknown:
        print(f"[WARN] no state for {len(unknown)} sample(s): {', '.join(unknown)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
