# ampprot — AMP proteomics harmonization

Builds a single harmonized HDF5 artifact from heterogeneous AMP proteomics deliveries — Olink
proximity-extension panels, DIA and TMT mass spectrometry, SomaScan aptamer arrays, SRM — each
arriving with its own sample sheet, its own feature identifiers, its own transform history and its
own idea of what a missing value means.

**This repository holds the code and the documentation. It holds no data.** Every source delivery
and every built artifact is controlled-access and stays outside it.

## Layout

```
src/ampprot/
    config.py       every stratum declared: paths, transforms, platform, evidence grade
    readers.py      one reader per delivery format; the format-specific knowledge lives here
    harmonize.py    transform reconciliation, identifier resolution, feature alignment
    writer.py       the HDF5 writer and its schema
    verify.py       the assertion audit - every claim the build makes about itself
    manifest.py     provenance: what was read, its checksum, what came out
    build.py        the entry point
    readme.py       generates the artifact's own README from the build

HARMONIZATION_PLAN.md    the design, stratum by stratum
REQUIREMENTS.md          what the artifact must satisfy

docs/ARTIFACT_README.md  generated documentation of a built artifact
```

## Running it

```
pip install -r requirements.txt
python -m ampprot.build
```

Source paths are declared per stratum in `config.py` and resolve against a local data tree that is
not part of this repository.

## Two design commitments worth knowing before reading the code

**Evidence is graded, not asserted.** Every platform and pipeline field carries an evidence grade —
`A` read from the delivery itself, `B` transcribed from a deposit README, `C` inferred. A field that
cannot be evidenced says so rather than being filled with something plausible.

**A distinction the delivery does not carry is not invented.** Where a product merges two
fractions, or where one sample token covers what were two runs, the build records the merge and
declines to split it. `config.py` carries the reasoning inline at each such stratum.

## What must never enter this repository

Source deliveries, built artifacts, sample sheets, clinical tables, donor identifiers — and
internal project correspondence: review threads, audit replies, issue discussion. What belongs here
is what describes the software. The
`.gitignore` is deny-by-default: everything at the root is ignored and files are re-included by
name. Add to that list deliberately, never by relaxing the rule.
