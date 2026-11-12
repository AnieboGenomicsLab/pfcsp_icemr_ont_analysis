# Intermediate processes

One-off work that is not part of the pipeline and is not run by `main.sh`.

A script belongs here when it prepares an input by hand rather than analysing data. The
usual case is pulling a few columns out of a spreadsheet, a REDCap export or a
collaborator's table and writing them into the metadata the pipeline reads. That is a
decision about one dataset, not a general step, so putting it in `main.sh` would make
every future run depend on a file that only ever mattered once.

Run them by hand, from the repository root, and record what was run in the commit
message. Each script writes its output somewhere the pipeline reads, never in place.

| Script | What it does |
|---|---|
| | |

See `process/00_metadata.py` for the shape these take.
