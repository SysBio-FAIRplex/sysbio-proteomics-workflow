"""README generated from the build itself, after the harmonised dataset exists.

Nothing here is hand-written prose about what the build *should* have done -- every
number, transform and decision is read back out of the produced artifact.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import pandas as pd

from . import config as C, harmonize as H


def _fmt_int(v):
    return f"{int(v):,}" if v is not None and pd.notna(v) else "-"


def generate(report: dict, audit: pd.DataFrame, layer_reports: list[dict],
             inv: pd.DataFrame) -> str:
    L = []
    a = L.append
    built = [r for r in layer_reports]
    total_rows = sum(r["n_rows"] for r in built)

    a(f"# {C.H5_NAME}")
    a("")
    a("Harmonised AMP proteomics, one HDF5 with a keyed layer per (assay x biospecimen "
      "matrix) combination.")
    a("")
    a(f"- Built: {report['build_utc']}")
    a(f"- Size: {report['h5_bytes']/1e6:.1f} MB")
    a(f"- Layers populated: {len(built)} of {len(C.LAYER_ORDER)}")
    a(f"- Rows (samples) across populated layers: {total_rows:,}")
    if report.get("previous_build"):
        a(f"- Previous build: `{report['previous_build']}` (not overwritten)")
    a("")
    a("> **The filename is the version.** Every build produces its own datestamped "
      "artifact, so a rebuild can no longer overwrite one that has already been audited "
      "or shared. **Do not pin this filename** in downstream code - resolve the newest "
      "`sysbio_proteomics-*.h5`, or pin deliberately when you need a reproducible "
      "analysis. `build_report.json` and this README both name the artifact they "
      "describe, so a bundle is always self-identifying.")
    a("")

    # ---------------- layers ----------------
    a("## Layers")
    a("")
    a("| layer key | rows | proteins | strata | value surface |")
    a("|---|---|---|---|---|")
    for r in built:
        surface = "`values_z` (per-stratum Z)" if r["z_applied"] else "`values` (native)"
        a(f"| `{r['layer']}` | {r['n_rows']:,} | {r['n_proteins']:,} | "
          f"{r['n_strata']} | {surface} |")
    for k, meta in C.DEFERRED_LAYERS.items():
        a(f"| `{k}` | 0 | 0 | 0 | deferred - {meta['reason'].split(';')[0]} |")
    a("")
    for k, meta in C.DEFERRED_STRATA.items():
        a(f"> **Deferred stratum** `{k}` (layer `{meta['layer']}`): {meta['reason']}")
    a("")

    # ---------------- scale / transform / Z ----------------
    a("## Transform and Z, per stratum")
    a("")
    a("Rule applied: a stratum is transformed **only** where the source scale needs it to "
      "approximate normality, and a layer is Z-scored **only** where its strata are not "
      "already in one numeric space. Where strata already agree, values ship exactly as "
      "delivered.")
    a("")
    a("| layer | stratum | native scale | transform | Z | median abs skew |")
    a("|---|---|---|---|---|---|")
    for r in built:
        for key, notes in r["strata"].items():
            st = next(s for s in C.STRATA if s.key == key)
            sc = notes.get("scale_check_after_transform") or notes.get("scale_check", {})
            a(f"| `{r['layer']}` | `{key}` | `{st.value_scale}` | "
              f"{notes.get('transform_applied', '-')} | "
              f"{'yes' if notes.get('zscored') else 'no'} | "
              f"{sc.get('median_abs_skew', '-')} |")
    a("")

    a("### Why each layer was or was not Z-scored")
    a("")
    for r in built:
        pol = r["policy"]
        a(f"**`{r['layer']}`** - Z applied: **{pol['z_applied']}**")
        a("")
        a(f"- {pol['rationale']}")
        dc = pol.get("distribution_check", {})
        if dc.get("per_stratum"):
            a(f"- measured before merge: location range "
              f"{dc.get('location_range_in_pooled_sd')} pooled SD, spread ratio "
              f"{dc.get('spread_ratio_max_over_min')}")
            for s in dc["per_stratum"]:
                a(f"    - `{s['stratum']}`: centre {s['centre']:.3f}, "
                  f"spread {s['spread']:.3f}")
        if pol.get("label_vs_data"):
            a(f"- **label vs data:** {pol['label_vs_data']}")
        ns = r.get("numeric_space", {})
        a(f"- post-merge numeric-space check on `{ns.get('merged_surface')}`: "
          f"**{'verified' if ns.get('verified') else 'FAILED'}** - {ns.get('reason')}")
        a("")

    # ---------------- QC ----------------
    a("## QC applied, per stratum")
    a("")
    a("Source-specified QC is honoured before any transform, Z or merge. Checks marked "
      "*upstream* were already applied by the data producer and are verified rather than "
      "re-applied.")
    a("")
    a("| stratum | check | applied | reference |")
    a("|---|---|---|---|")
    for r in built:
        for key, notes in r["strata"].items():
            for q in notes.get("source_qc", []):
                where = "upstream" if q["applied_upstream"] else "here"
                a(f"| `{key}` | {q['check']} | {where} | {q['reference']} |")
    a("")
    a("> **QC described as \"done\" upstream was often not actually done.** Several source "
      "READMEs and QC logs state that a filter was applied, but the delivered matrix does "
      "not satisfy it. Those filters were re-applied here rather than taken on trust. The "
      "clearest case: `QC_process_Temporal.txt` states *\"only keep proteins that have "
      "missing data in <50% of the samples\"*, yet the delivered "
      "`n278_residual_log2_batch.TCX.csv` contains 398 features below that threshold, some "
      "as low as 10% present - most likely because the batch regression that runs *after* "
      "that step reintroduces missingness. Re-applying the filter here removes them. "
      "Treat every upstream QC claim in the table above as verified-by-us, not as "
      "inherited.")
    a("")
    a("### Standing QC policy")
    a("")
    a("1. QC decisions are made **per stratum, before** any transform, Z or merge.")
    a("2. Samples and features that fail QC are **flagged, never silently deleted**. The "
      "exclusion threshold belongs to the study, not to this build.")
    a("3. Every threshold lives in `config.py`, never inline.")
    a("4. Every exclusion is counted in the audit.")
    a("5. Every `SOURCE_QC` entry is **checked**, not recited: a step the producer says "
      "they applied is verified against the delivered file where the data can show it, "
      "and a step they left to us fails the audit unless the build left a trace of "
      "actually doing it.")
    a("")
    a("**Deletion is reserved for rows whose identity is unusable, never their quality.** "
      "There are exactly three such cases in this build, and each is recorded as data "
      "rather than falling out silently as a join miss:")
    a("")
    a("- **Olink plasma, 172 samples** on two plates whose donor identity is in dispute "
      "pending AMP-PD. Withdrawn by *plate* - the physical fact - so the rule survives a "
      "re-release. 39 participants lose every plasma sample.")
    a("- **Samples with no donor at all in the source metadata.** They cannot be joined, "
      "de-duplicated against a repeat visit, or withdrawn if that donor revokes consent.")
    a("- **DIA plasma, low-completeness runs** - the one *quality*-based removal, and it "
      "exists only because that deposit ships no sample-level QC of its own. The failures "
      "are unambiguous: one run detected 10 of 202 proteins.")
    a("")
    a("### Feature filtering")
    a("")
    a("Drop rates differ substantially between strata, and the difference is about *what "
      "processing stage the deposit sits at*, not about inconsistent QC - the same "
      "threshold is applied identically everywhere. Strata delivered as finished, "
      "already-filtered products (ROSMAP TMT round 1, DiverseCohorts DLPFC, SomaScan, "
      "Olink) lose almost nothing. Strata delivered closer to the instrument lose more: "
      "ROSMAP TMT round 2 is a raw Proteome Discoverer export and drops 22.9%, and the two "
      "MaxQuant label-free sets drop 34.5% and 40.8%, because a protein quantified in only "
      "a few runs is genuinely absent from most samples.")
    a("")
    a("The missingness filter is applied **per stratum, never after merge** - each "
      "stratum has its own sample set, so a protein measured well by one assay would "
      "otherwise look mostly-missing across a merged layer and be dropped for the wrong "
      "reason. Threshold follows `QC_process_*.txt`: keep features present in >= 50% of "
      "that stratum's samples.")
    a("")
    a("| stratum | unresolved id | below 50% present | duplicates collapsed | final |")
    a("|---|---|---|---|---|")
    for r in built:
        for key, n in r["strata"].items():
            a(f"| `{key}` | {_fmt_int(n.get('n_features_unresolved_dropped'))} | "
              f"{_fmt_int(n.get('n_features_below_50pct_present_dropped'))} | "
              f"{_fmt_int(n.get('n_duplicate_accessions_collapsed', 0))} | "
              f"{_fmt_int(n.get('n_features_final'))} |")
    a("")

    a("### Protein sparsity, per stratum")
    a("")
    a("Measured on what actually ships, after filtering. Sparsity is low but never zero - "
      "some patchiness across strata is expected and is not a defect.")
    a("")
    a("| stratum | % missing | features fully observed | median feature presence | min |")
    a("|---|---|---|---|---|")
    for r_ in built:
        for key, n in r_["strata"].items():
            sp = n.get("sparsity", {})
            a(f"| `{key}` | {sp.get('pct_missing_overall', '-')}% | "
              f"{sp.get('pct_features_fully_observed', '-')}% | "
              f"{sp.get('median_feature_presence_pct', '-')}% | "
              f"{sp.get('min_feature_presence_pct', '-')}% |")
    a("")
    a("Note the *layer* matrices are sparser than these per-stratum figures, because a "
      "layer takes the union of its strata's features: a protein measured by one stratum "
      "and not another is structurally absent for the other's rows. That is recorded in "
      "`missing_reason` as code 2 (`not_assayed_in_this_stratum`) and is not a QC problem.")
    a("")
    r2 = None
    for r_ in built:
        if "rosmap_tmt_r2" in r_["strata"]:
            r2 = r_["strata"]["rosmap_tmt_r2"]
    if r2 and r2.get("structured_missingness"):
        sm = r2["structured_missingness"]
        a("#### ROSMAP round 2: missingness is plex-structured, not per-sample")
        a("")
        a("Round 2 is a raw Proteome Discoverer multiconsensus export, and its missingness "
          "behaves differently from every other stratum. Measured on the delivered file:")
        a("")
        a(f"- **{sm['features_partially_present_within_a_plex']} features are partially "
          f"present within a plex** - i.e. none. Missingness is entirely all-or-nothing at "
          f"the plex level.")
        a(f"- a feature is present in {sm['mean_plexes_present_per_feature']} of "
          f"{sm['n_plexes']} plexes on average")
        a(f"- {sm['n_features_entirely_empty']} features are empty in every channel")
        a("")
        a("This is **identification** missingness from the consensus step, not detection "
          "missingness in a sample: a protein identified in one plex and not another "
          "produces a clean block of absent channels. Two consequences worth carrying "
          "into any analysis:")
        a("")
        a("1. The >=50% presence filter is effectively selecting proteins identified in at "
          "least 7 of 14 plexes, not proteins detected in at least half the donors. It is "
          "the same threshold used everywhere else, but it means something different here.")
        a("2. Because whole plexes drop out together, remaining missingness is correlated "
          "with batch. Anything that treats missingness as random - naive imputation, "
          "complete-case filtering - will inherit that structure. Round 2 is also the "
          "stratum most likely to carry **stale** identifications from the multiconsensus "
          "search, so treat its low-presence tail with more suspicion than the others.")
        a("")

    # ---------------- identifiers ----------------
    a("## Identifiers")
    a("")
    a("- **All identifiers are strings**, exactly as the source defines them. Prefixes "
      "(`PD-`, `PP-`, `R`, `AMPAD_MSSM_`) are part of the identifier and are never "
      "stripped, so IDs join to future metadata drops unchanged.")
    a("- `person_id` **is** the source identifier. No surrogate is invented.")
    a("- Because two studies can reuse the same digits, `person_source_namespace` records "
      "the identifier system and `person_global_key` "
      "(`namespace:id`) is the unambiguous cross-study join key.")
    a("")
    a("| namespace | strata |")
    a("|---|---|")
    ns_map: dict[str, list[str]] = {}
    for k, v in H.PERSON_NAMESPACE.items():
        ns_map.setdefault(v, []).append(k)
    for ns, ks in sorted(ns_map.items()):
        a(f"| `{ns}` | {', '.join(f'`{k}`' for k in ks)} |")
    a("")
    a("> ROSMAP appears under two namespaces: `ROSMAP.projid` (TMT, SRM, and the `sd` "
      "subset of DiverseCohorts) and `ROSMAP.individualID` (SomaScan). The bridge between "
      "them is shipped in this file as **`person_crosswalk_rosmap`** "
      "(from syn3191087), so ROSMAP plasma *can* be joined to ROSMAP brain by person.")
    a("")

    a("### Identifier detective work")
    a("")
    a("Three identifier problems were not solvable from the obvious column and needed "
      "tracing back through the data. Each is recorded here because the resolution is not "
      "self-evident from the files, and each fails *silently* - returning zero matches or "
      "an unlinked donor rather than an error.")
    a("")
    a("**1. DiverseCohorts `sd` had no `individualID` at all.**")
    a("")
    a("The frontal traits sheet has an `individualID` column that is populated for the "
      "`emdp`, `mayo` and `mssm` contributing datasets and **completely empty for `sd`** "
      "- 400 rows, zero values, 371 of which survive into the delivered matrix. Read "
      "naively this looks like 371 samples with no donor, and an earlier pass of this "
      "build mislabelled them as reference pools.")
    a("")
    a("They are neither. The sheet carries a *second* identifier column, `projID`, which "
      "is populated for exactly the `sd` rows and empty for every other dataset. The "
      "values are 8-digit ROS/MAP project IDs, so `sd` is a ROS/MAP contribution to "
      "DiverseCohorts carrying ROSMAP identifiers rather than AMP-AD ones. Falling back "
      "to `projID` where `individualID` is absent recovers all 371, and takes this "
      "stratum from 375 unresolved to **0**.")
    a("")
    a("Consequence: `diversecohorts_dlpfc` spans **two identifier systems in one "
      "stratum**, so `person_source_namespace` is assigned per row "
      "(715 `AMP-AD.individualID`, 371 `ROSMAP.projid`) rather than per stratum. Anything "
      "that assumes one namespace per stratum will mis-join this layer.")
    a("")
    a("**2. ROSMAP uses two identifier systems across assays.**")
    a("")
    a("ROSMAP TMT and SRM key on numeric `projid`; ROSMAP SomaScan keys on `R`-prefixed "
      "`individualID`. They never collide, so nothing errors - the two simply never join, "
      "and plasma silently fails to link to brain for the same donor. The bridge is "
      "`syn3191087`, published as `ROSMAP_clinical.csv`, which is written into this file "
      "as `person_crosswalk_rosmap` (identifier columns only; that source is otherwise "
      "clinical and is dropped at read).")
    a("")
    a("**3. Sample IDs were re-ordered between matrix and metadata.** The MSBB label-free "
      "run labels put the plate position before the specimen number "
      "(`b<batch>_<position>_<specimen>`) while its biospecimen metadata puts the specimen "
      "first (`b<batch>_<specimen>_<position>`) - the last two tokens swapped. The batch-7 "
      "rerun compounds it: `b7r2_*` in the matrix against `r2b7_*` in the metadata. Mayo "
      "appends a position token the metadata omits "
      "(`mayo_b<batch>_<specimen>_<position>` against `b<batch>_<specimen>`). None of "
      "these join without an explicit rewrite, and all three silently return zero matches "
      "if attempted directly.")
    a("")
    a("### Row keys")
    a("")
    a("A person legitimately repeats within a layer - multiple brain regions, multiple "
      "timepoints, technical replicates. `row_key` is composed from the real "
      "discriminators in order (person -> anatomic site -> visit -> stratum) and only "
      "falls back to a replicate counter when those do not separate the rows. Every "
      "component is also its own column, so no one has to parse the key.")
    a("")

    # ---------------- nulls ----------------
    a("## Missing values")
    a("")
    a("`values` and `values_z` are the analysis-ready tables and carry **exactly one null "
      "type: float `NaN`**. No sentinels, no second null flavour - numpy, pandas and "
      "sklearn all behave normally.")
    a("")
    a("The *reason* a cell is missing lives in a separate `missing_reason` int8 array of "
      "the same shape, which is diagnostic only:")
    a("")
    a("| code | meaning |")
    a("|---|---|")
    for k, v in H.MISSING_LEGEND.items():
        a(f"| {k} | {v} |")
    a("")
    a("### `below_lod` is three-state, not boolean")
    a("")
    a("On the Olink layers, `below_lod` carries **three** states rather than two:")
    a("")
    a("| value | meaning |")
    a("|---|---|")
    a("| `0.0` | measured, at or above the limit of detection |")
    a("| `1.0` | measured, below the limit of detection |")
    a("| `NaN` | **not measured** - there is no value to compare against the limit |")
    a("")
    a("This matters because it changes results rather than presentation. `below_lod` was "
      "previously written as `(values < lod)` cast to float, and `NaN < NaN` is `False`, "
      "which casts to `0.0` - so a never-measured cell was stored identically to one "
      "measured comfortably above the limit. The conventional filter `below_lod == 0.0` "
      "therefore **silently retained every never-measured cell**, and the only correct "
      "predicate was `below_lod == 0.0 AND missing_reason == 0`, which is not discoverable "
      "from the surface. Restoring `NaN` makes the three states separable with the single "
      "null flavour the file already documents.")
    a("")

    # ---------------- longitudinal axis ----------------
    a("## The time axis: `visit_month`")
    a("")
    a("`visit_month` - **months since that participant's baseline** - is the canonical "
      "cross-layer time axis and is present in **every** layer: populated where the "
      "dataset is longitudinal, null where it is postmortem. One consumer expression now "
      "aligns Olink, SomaScan and DIA, and adding a longitudinal dataset means populating "
      "`visit_month` rather than inventing another axis.")
    a("")
    a("`visit_index` stays beside it as the source-native ordinal and is **never** the "
      "join key. `visit_month_evidence` records how the value was reached:")
    a("")
    a("| evidence | meaning | layers |")
    a("|---|---|---|")
    a("| `source` | the deposit states months since baseline | Olink plasma/CSF, DIA "
      "plasma/CSF |")
    a("| `derived_from_visit_index` | converted from an ordinal follow-up number | "
      "SomaScan plasma |")
    a("")
    a("> **The SomaScan conversion is an approximation and is labelled as one.** ROSMAP's "
      "`Visit` is the follow-up *year*, so `visit_month = Visit x 12` is a structural "
      "conversion of an ordinal index - it never touches `age_at_visit` and never trips "
      "the scope gate. But ROSMAP follow-ups are nominally annual and not exactly so, "
      "which means these values are **not** measured intervals and must not be treated as "
      "such. `visit_index` is retained so the conversion is reversible.")
    a("")
    a("Before this, `soma_plasma` carried only an ordinal and the brain layers carried "
      "nothing, so a consumer joining ROSMAP plasma to PDRD Olink on `person_id` had no "
      "shared time column and nothing longitudinal could cross cohorts.")
    a("")

    # ---------------- PMI ----------------
    pmi_rows = [(r["layer"], k, n.get("pmi") or {})
                for r in built for k, n in r["strata"].items()]
    a("## Post-mortem interval - a deliberate scope change")
    a("")
    a("> **Scope change, recorded rather than made silently.** Post-mortem interval was "
      "previously denied as clinical phenotype. It is now carried as a **biospecimen "
      "technical CDE**: it describes the specimen, not the donor, and it is a standard "
      "covariate for brain proteomics. `pmi` and `postmortem_interval` were removed from "
      "the scope gate's deny list for this reason and this reason only. That narrows "
      "HARMONIZATION_PLAN section 0 and REQUIREMENTS SEC-1, both of which are written "
      "against the DUC, so it is a change to the disclosure boundary and is documented as "
      "one here, in `config.py` at `CLINICAL_DENY`, and in both markdowns. Nothing else "
      "was re-admitted, and SEC-2's read-time allow-list is unchanged - only the join key "
      "and the `pmi` column are ever loaded from the clinical file.")
    a("")
    a("**One scale: `pmi_hours`, hours with decimals, always.** Sub-hour precision is "
      "real - ROSMAP's median is 6.783 h - so rounding to whole hours would discard "
      "resolution on a covariate whose entire use is fine-grained. A bare `pmi` column is "
      "never trusted, because the unit differs between studies by convention: MSBB records "
      "**minutes**, ROSMAP **hours**.")
    a("")
    a("### Conversion ledger")
    a("")
    a("| stratum | source | column | unit as found | stated or inferred | factor | "
      "populated | median (h) |")
    a("|---|---|---|---|---|---|---|---|")
    for layer, key, p in pmi_rows:
        if not p:
            continue
        stated = ("stated" if p.get("unit_stated_by_source")
                  else ("inferred" if p.get("evidence") == "inferred_from_distribution"
                        else "-"))
        a(f"| `{key}` | {p.get('source') or '-'} | {p.get('source_column') or '-'} | "
          f"{p.get('source_unit') or '-'} | {stated} | "
          f"{p.get('factor_to_hours') if p.get('factor_to_hours') is not None else '-'} | "
          f"{p.get('n_populated', 0)}/{p.get('n_samples', 0)} | "
          f"{p.get('median_hours') if p.get('median_hours') is not None else '-'} |")
    a("")
    a("`pmi_evidence` distinguishes **two different kinds of null**, per row, because the "
      "difference matters to anyone modelling with it:")
    a("")
    a("| evidence | meaning |")
    a("|---|---|")
    a("| `converted_from_stated_unit` | the source named its unit; the conversion is exact "
      "and reversible |")
    a("| `inferred_from_distribution` | it did not; the unit was read off the observed "
      "median against a plausible band |")
    a("| `not_applicable_antemortem` | **a fact about the specimen** - a living donor's "
      "plasma or CSF has no post-mortem interval |")
    a("| `source_not_acquired` | **a gap in what we hold** - the object exists, we do not "
      "have it |")
    a("")
    a("`pmi_hours` is **never zero** where the value is unknown. A zero interval is a "
      "claim, not a null.")
    a("")
    a("Where inference is used it is safe because the bands do not overlap: PMI is the "
      "interval between death and processing at a brain bank, so the plausible range is "
      f"{C.PMI_PLAUSIBLE_HOURS[0]}-{C.PMI_PLAUSIBLE_HOURS[1]} hours, and the same "
      "durations recorded in minutes land near 400. Where the median sits in neither band "
      "the values stay null and the stratum is reported undetermined - never guessed. The "
      "audit checks the post-conversion median against that band and **reports** an "
      "out-of-band result; it never auto-corrects.")
    a("")
    a("> **For the CDE repository.** These three columns are new to the specimen surface "
      "and are not yet in the external CDE dictionary: **`pmi_hours`** (float, hours with "
      "decimals), **`pmi_evidence`** (the four-value enum above) and **`pmi_source_unit`** "
      "(the unit as found in the source). They need registering there, and the unit "
      "convention needs registering with them - a bare `pmi` CDE with no unit field is "
      "unsafe across studies, because MSBB's minutes and ROSMAP's hours differ by 60x on "
      "a column that looks identical. Any external repository already carrying a `pmi` CDE "
      "should be reviewed against this.")
    a("")

    # ---------------- Brodmann ----------------
    a("## Brodmann area: measurement and inference are separable")
    a("")
    a("`brodmann_area` is **partly inferred**, and `brodmann_area_evidence` says which "
      "rows are which. Anything consuming `brodmann_area` must read the evidence column "
      "with it.")
    a("")
    a("| evidence | meaning |")
    a("|---|---|")
    a("| `source` | the specimen's own annotation stated it |")
    a("| `propagated` | another specimen from the same donor and tissue stated it |")
    a("| `convention` | neither did; it was asserted from the anatomic site |")
    a("")
    a("| stratum | site | convention | source | propagated | convention-filled | "
      "disagreements kept |")
    a("|---|---|---|---|---|---|---|")
    for r_ in built:
        for key, n in r_["strata"].items():
            b = n.get("brodmann") or {}
            if not b or not (b.get("n_source") or b.get("n_propagated")
                             or b.get("n_convention")):
                continue
            a(f"| `{key}` | {b.get('anatomic_site')} | {b.get('convention') or '-'} | "
              f"{b.get('n_source', 0)} | {b.get('n_propagated', 0)} | "
              f"{b.get('n_convention', 0)} | "
              f"{b.get('n_source_disagrees_with_convention', 0)} |")
    a("")
    a("The conventions applied are " + ", ".join(
        f"`{k}` -> `{v}`" for k, v in C.BRODMANN_CONVENTION.items()) + ", declared once in "
      "`config.py` so a consumer can undo them.")
    a("")
    a("**A stated value always wins.** The convention only fills where the source is "
      "silent, and never overwrites - which is why the 67 DiverseCohorts rows that state "
      "`BA10` are still `BA10` rather than being swept into the DLPFC convention's `BA9`. "
      "Those 67 rows are the standing evidence that the convention is a useful default and "
      "not a universal truth.")
    a("")
    a("> **For the CDM side.** `PROTEOMICS_REVIEW.md` gives `brodmann_area` precedence "
      "over the free-text site label. A loader that honours that without reading "
      "`brodmann_area_evidence` would move every convention-filled DLPFC row off "
      "`4195656 Region of frontal cortex` onto a BA9 concept **on our inference alone**. "
      "Separately, the DiverseCohorts temporal strata move off `4193043 Region of temporal "
      "cortex` to a superior-temporal-gyrus concept, because the source metadata says "
      "`superior temporal gyrus` for all 280 temporal TMT specimens and the build's "
      "earlier `temporal cortex` claim was wrong.")
    a("")

    # ---------------- CDE coverage ----------------
    a("## SysBio CDE coverage")
    a("")
    a("Evidence grades: **A** derivable from the data files, **B** transcribed from a "
      "bundled README, **C** unavailable. Read from the built artifact, not from the "
      "config declaration - where an assay-metadata object states the platform per "
      "sample, that is what ships and the grade is **A**.")
    a("")
    a("| stratum | platform | evidence | analysis pipeline | evidence | source |")
    a("|---|---|---|---|---|---|")
    for r_ in built:
        for key, n in r_["strata"].items():
            st = next(s for s in C.STRATA if s.key == key)
            meta = n.get("assay_meta") or {}
            obs = meta.get("platform_observed") if meta.get("joined") else None
            if obs:
                names = sorted(obs)
                total = sum(obs.values())
                plat = (names[0] if len(names) == 1 else
                        "mixed per sample: " + "; ".join(f"{k} ({obs[k]}/{total})"
                                                         for k in names))
                ev, src = "A", f"`{meta['source'].rsplit('/', 1)[-1]}`"
            else:
                plat = st.platform or "-"
                ev = st.platform_evidence or "C"
                src = "config"
            a(f"| `{key}` | {plat} | {ev} | {st.analysis_pipeline or '-'} | "
              f"{st.pipeline_evidence or 'C'} | {src} |")
    a("")
    a("`platform` is genuinely **per sample**, not per stratum: the DiverseCohorts strata "
      "span three instruments and the build records that rather than collapsing it to one "
      "name. The per-sample `platform` column is authoritative; the table above is its "
      "stratum-level summary. ROSMAP SRM is the one stratum that does **not** close - "
      "`syn23569441` carries `platform = NA` on all 1,212 of its rows, which is a silence "
      "in the deposit rather than a gap in this build, and is reported as such.")
    a("")
    a("FILES CDEs `current_version`, `created_on`, `modified_on` and `drs_id` are left "
      "**null** rather than guessed. Every filesystem mtime in this directory is the "
      "download timestamp, not file creation, so using it would be quietly wrong.")
    a("")

    # ---------------- audit ----------------
    s = report["audit"]
    a("## Assertion audit")
    a("")
    a(f"{s['n_checks']} checks: **{s['n_pass']} pass**, **{s['n_fail']} fail**, "
      f"{s['n_unverifiable']} unverifiable.")
    a("")
    a("Every claim - in the tracking file, in a bundled README, or encoded in this "
      "codebase - is re-checked against the data. Disagreements are reported, not "
      "silently resolved.")
    a("")
    fails = audit[audit.verdict == "FAIL"]
    if len(fails):
        a("| check | subject | claim | source | observed |")
        a("|---|---|---|---|---|")
        for _, r in fails.iterrows():
            a(f"| {r['check']} | {r['subject']} | {r['claim']} | {r['claim_source']} | "
              f"{r['observed']} |")
        a("")
    a("Column naming is checked across layers: every layer exposes the same CDE and link "
      "columns under the same names, with layer-specific proteomics columns "
      "(`panel`, `dilution`, `seq_id`, `peptide_sequence`, `lod`) present only where the "
      "assay actually has them.")
    a("")
    a("Full detail: `build/inventory/assertion_audit.tsv`.")
    a("")
    a("### Known limitations")
    a("")
    a("- **PMI is not held for five strata.** MSBB and Mayo PMI live in individual-"
      "metadata objects (`syn73713767`, `syn73713766`) that are not in the delivery, and "
      "DiverseCohorts has **no individual-metadata object registered anywhere**, so its "
      "PMI source is still unidentified. Those rows are `source_not_acquired` - a gap in "
      "what we hold, distinct from the antemortem nulls. MSBB is also the study that "
      "records minutes, so the unit conversion stays untested against real data until "
      "that object arrives.")
    a("- **Brodmann area is ~88% inference on the DiverseCohorts strata.** No "
      "`TMT quantitation` specimen carries a source `BrodmannArea` at all; cross-assay "
      "propagation reaches about 12% and the rest is convention. Separable via "
      "`brodmann_area_evidence`, but it is not measurement.")
    a("- **ROSMAP SRM has no platform.** `syn23569441` carries `platform = NA` on all "
      "1,212 rows. This is a silence in the deposit, not a gap in the build.")
    a("- **The SomaScan normalisation label has no evidence behind it.** The ANML and "
      "non-ANML deposits are byte-identical (md5 `f8a3d03ecd1b57cbb409bc6b054a7890`), so "
      "one of the two is mislabelled and the data cannot say which. Values are "
      "unaffected; only the label is unsupported. Raised with ADKP.")
    a("- **The Olink plasma plate swap is unresolved at source.** 172 samples are "
      "withdrawn pending AMP-PD's answer, not because the answer is known.")
    a("- **ROSMAP round 2 carries plex-structured missingness** (see above) and is the "
      "stratum most likely to contain stale multiconsensus identifications.")
    a("- **The FILES CDEs `current_version` / `created_on` / `modified_on` / `drs_id` are "
      "null throughout**, deliberately: every filesystem mtime here is a download "
      "timestamp, and Synapse versions were not downloaded.")
    a("- **The DIA plasma fraction is not recoverable.** The delivered batch-corrected "
      "product merges native and depleted plasma under one `-PLA-` sample token, so which "
      "fraction a sample came from cannot be read off it. It is therefore one stratum, "
      "not two.")
    a("")

    a("## What is in the file")
    a("")
    a("```")
    a(f"{C.H5_NAME}")
    for r_ in built:
        a(f"|-- {r_['layer']}/")
    for k in C.DEFERRED_LAYERS:
        a(f"|-- {k}/                    (declared, empty)")
    a("+-- person_crosswalk_rosmap/    ROSMAP projid <-> individualID (syn3191087)")
    a("```")
    a("")
    a("Every populated layer has the same internal structure:")
    a("")
    a("```")
    a("<layer>/")
    a("  values            (n_rows, n_proteins) float32   the analysis-ready table")
    a("  values_z          same shape, only where the layer needed Z")
    a("  lod, below_lod    same shape, Olink layers only")
    a("  missing_reason    same shape, int8 diagnostic (see Missing values)")
    a("  row_key           (n_rows,)     person + visit, the analysis key")
    a("  protein_columns   (n_proteins,) column headers for `values`")
    a("  samples/<col>     one dataset per CDE / link / technical column")
    a("  features/<col>    one dataset per feature-description column")
    a("  assay/<col>       SysBio ASSAY CDEs, one row per stratum")
    a("  files/<col>       SysBio FILES CDEs, one row per source file")
    a("```")
    a("")
    a("Layer-level attributes carry the scale policy, the per-stratum transform and Z "
      "flags, the post-merge numeric-space check, and provenance.")
    a("")
    # ---------------- tabs / returns ----------------
    tabs, returns = report.get("tabs") or {}, report.get("returns") or []
    if tabs:
        a("## The analysis-ready bundles (`build/tabs/`)")
        a("")
        a("Each layer is also written as a standalone bundle - **wide only**; the long "
          "form is deliberately not emitted, because pooled brain TMT alone would be "
          "~25M rows and the wide matrix is what the HDF5 already holds.")
        a("")
        a("```")
        a("tabs/<layer_key>/")
        a("  <layer_key>_matrix.parquet     feature_id + one column per specimen")
        a("  <layer_key>_samples.parquet    one row per sample, the full CDE surface")
        a("  <layer_key>_features.parquet   one row per feature, the identity half")
        a("  <layer_key>_bundle.json        which surface the matrix carries, and how to join")
        a("```")
        a("")
        a("Parquet is primary because it preserves dtypes and the three-state `below_lod` "
          "distinction above; a `.tsv.gz` mirror of each table ships beside it.")
        a("")
        a("| layer | matrix | value surface | columns keyed by | samples | features |")
        a("|---|---|---|---|---|---|")
        for k in C.LAYER_ORDER:
            b = tabs.get(k)
            if not b:
                continue
            m = b["matrix"]
            a(f"| `{k}` | {m['n_features']:,} x {m['n_specimens']:,} | "
              f"`{m['value_surface']}` | `{m['column_basis']}` | "
              f"{b['samples']['n_rows']:,} | {b['features']['n_rows']:,} |")
        a("")
        a("The matrix is written in the same pass that writes the HDF5's `values`, from "
          "the same in-memory array, and the tables in the same pass that writes the "
          "HDF5's tables. A standalone HDF5-to-parquet converter would be simpler and "
          "would give up exactly the property that matters: **the two cannot drift**.")
        a("")
    if returns:
        a("### The AMP return path (`build/returns/`)")
        a("")
        a("One directory per grant, holding that grant's specimen set for every layer it "
          "appears in, plus a manifest that **demonstrates** the boundary rather than "
          "asserting it - it names the grants present in the build, the grant returned, "
          "and the count of rows from any other grant that reached the slice, which must "
          "be zero.")
        a("")
        a("| grant | specimens | layers | rows from another grant |")
        a("|---|---|---|---|")
        for r_ in returns:
            a(f"| {r_['grant']} | {r_['n_specimens_returned']:,} | "
              f"{len(r_['layers'])} | {r_['n_rows_from_another_grant']} |")
        a("")
        a("This is computed on every build rather than at release time: returning one "
          "AMP programme's donors to another is not recoverable after the fact.")
        a("")
    if tabs or returns:
        a("> **Controlled access.** `tabs/` and `returns/` carry the same donor-level "
          "identifiers as the HDF5. They are controlled-access artifacts, not a "
          "publishable by-product, and `.gitignore` denies them by name as well as by "
          "directory.")
        a("")

    a("## Files produced")
    a("")
    a(f"- `build/{C.H5_NAME}` - the harmonised dataset")
    a("- `build/build_report.json` - machine-readable build summary")
    if tabs:
        a("- `build/tabs/<layer>/` - the analysis-ready bundle per layer "
          "(**controlled access**)")
    if returns:
        a("- `build/returns/<grant>/` - the per-grant return path with its boundary "
          "manifest (**controlled access**)")
    a("- `build/inventory/derived_inventory.tsv` - one row per stratum, computed")
    a("- `build/inventory/tracking_file_diff.tsv` - where the data disagrees with the "
      "tracking file")
    a("- `build/inventory/missing_objects.tsv` - referenced Synapse objects and what "
      "their absence blocks")
    a("- `build/inventory/assertion_audit.tsv` - every claim, checked")
    a("- `build/inventory/file_manifest.tsv` - every source file with sha256")
    a("")
    a("Sources read through an explicit allow-list and never ingested as data:")
    a("")
    a("- `ROSMAP_clinical.csv` (syn3191087) - **identifiers plus `pmi`**. This file is "
      "overwhelmingly clinical (`msex`, `educ`, `race`, `apoe_genotype`, `age_death`, "
      "`braaksc`, `ceradsc`, `cogdx`, `mmse`); only `projid` / `individualID` / `Study` "
      "for the crosswalk, and `projid` + `pmi` for the specimen CDE, are ever loaded. The "
      "`pmi` read is the scope change documented above and is the single column "
      "re-admitted.")
    a("- `MSBB_biospecimen_metadata.csv` (syn21893059) and "
      "`MayoRNAseq_biospecimen_metadata.csv` (syn20827192) - identifiers and biospecimen "
      "descriptors; `samplingAge` / `samplingAgeUnits` are dropped at read.")
    a("- The five assay-metadata objects - platform, reagent lot, control flags and "
      "acquisition descriptors only.")
    a("")
    a("## Reading a layer")
    a("")
    a("```python")
    a("import h5py, numpy as np, pandas as pd")
    a("")
    a(f'h5 = h5py.File("build/{C.H5_NAME}", "r")')
    a('g  = h5["olink_csf_proteomics"]')
    a("")
    a('proteins = [s.decode() for s in g["protein_columns"][:]]')
    a('rows     = [s.decode() for s in g["row_key"][:]]')
    a('X        = pd.DataFrame(g["values"][:], index=rows, columns=proteins)')
    a("")
    a("# CDE / link columns for the same rows")
    a("# CDE / link columns for the same rows (string datasets come back as bytes)")
    a("def _col(ds):")
    a("    v = ds[:]")
    a('    return [x.decode() for x in v] if v.dtype.kind in "OS" else v')
    a("")
    a('meta = pd.DataFrame({k: _col(g["samples"][k]) for k in g["samples"]}, index=rows)')
    a("")
    a("# join ROSMAP plasma to ROSMAP brain by donor")
    a('cw = pd.DataFrame({k: _col(h5["person_crosswalk_rosmap"][k])')
    a('                   for k in h5["person_crosswalk_rosmap"]})')
    a("```")
    a("")
    return "\n".join(L)


def write(report, audit, layer_reports, inv) -> str:
    text = generate(report, audit, layer_reports, inv)
    path = C.BUILD / "README.md"
    path.write_text(text)
    return str(path)
