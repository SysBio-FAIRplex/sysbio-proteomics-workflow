# sysbio_proteomics-09042026.h5

Harmonised AMP proteomics, one HDF5 with a keyed layer per (assay x biospecimen matrix) combination.

- Built: 2026-09-04T16:35:33.428381+00:00
- Size: 255.2 MB
- Layers populated: 8 of 8
- Rows (samples) across populated layers: 13,336
- Previous build: `sysbio_proteomics-08222026.h5` (not overwritten)

> **The filename is the version.** Every build produces its own datestamped artifact, so a rebuild can no longer overwrite one that has already been audited or shared. **Do not pin this filename** in downstream code - resolve the newest `sysbio_proteomics-*.h5`, or pin deliberately when you need a reproducible analysis. `build_report.json` and this README both name the artifact they describe, so a bundle is always self-identifying.

## Layers

| layer key | rows | proteins | strata | value surface |
|---|---|---|---|---|
| `tmt_brain_proteomics` | 2,160 | 14,916 | 5 | `values_z` (per-stratum Z) |
| `lfq_brain_proteomics` | 536 | 4,437 | 2 | `values` (native) |
| `srm_brain_proteomics` | 2,704 | 185 | 2 | `values` (native) |
| `soma_plasma_proteomics` | 970 | 7,289 | 1 | `values` (native) |
| `olink_plasma_proteomics` | 1,380 | 1,463 | 1 | `values` (native) |
| `olink_csf_proteomics` | 1,320 | 1,463 | 1 | `values` (native) |
| `dia_plasma_proteomics` | 1,459 | 202 | 1 | `values` (native) |
| `dia_csf_proteomics` | 2,807 | 249 | 1 | `values` (native) |


## Transform and Z, per stratum

Rule applied: a stratum is transformed **only** where the source scale needs it to approximate normality, and a layer is Z-scored **only** where its strata are not already in one numeric space. Where strata already agree, values ship exactly as delivered.

| layer | stratum | native scale | transform | Z | median abs skew |
|---|---|---|---|---|---|
| `tmt_brain_proteomics` | `diversecohorts_dlpfc` | `log2_ratio_total_batch_residual` | none | yes | 0.529 |
| `tmt_brain_proteomics` | `diversecohorts_temporal` | `log2_ratio_total_batch_residual` | none | yes | 0.47 |
| `tmt_brain_proteomics` | `rosmap_tmt_r1` | `log2_rel_batch_median` | none | yes | 0.637 |
| `tmt_brain_proteomics` | `msbb_tmt_phg` | `log2_reporter_intensity` | log2 | yes | 0.337 |
| `tmt_brain_proteomics` | `rosmap_tmt_r2` | `log2_rel_batch_median` | log2 | yes | 0.628 |
| `lfq_brain_proteomics` | `mayo_lfq_tcx` | `log2_lfq_intensity` | log2 | no | 0.472 |
| `lfq_brain_proteomics` | `msbb_lfq_pfc` | `log2_lfq_intensity` | log2 | no | 0.365 |
| `srm_brain_proteomics` | `rosmap_srm_panel1` | `log2_ratio` | none | no | 0.691 |
| `srm_brain_proteomics` | `rosmap_srm_panel2` | `log2_ratio` | none | no | 0.692 |
| `soma_plasma_proteomics` | `rosmap_soma_plasma` | `log10_rfu_anml` | none | no | 2.536 |
| `olink_plasma_proteomics` | `pdrd_olink_plasma` | `npx_log2` | none | no | 0.448 |
| `olink_csf_proteomics` | `pdrd_olink_csf` | `npx_log2` | none | no | 0.511 |
| `dia_plasma_proteomics` | `pdrd_dia_plasma` | `log2_intensity_batch_corrected` | none | no | 1.157 |
| `dia_csf_proteomics` | `pdrd_dia_csf` | `log2_intensity_batch_corrected` | none | no | 0.578 |

### Why each layer was or was not Z-scored

**`tmt_brain_proteomics`** - Z applied: **True**

- measured distributions differ across strata (strata differ materially (location 113.62 pooled SD, spread ratio 1.61)), so each stratum is Z-scored on its own samples before merging
- measured before merge: location range 113.623 pooled SD, spread ratio 1.606
    - `diversecohorts_dlpfc`: centre 0.002, spread 0.150
    - `diversecohorts_temporal`: centre 0.002, spread 0.156
    - `rosmap_tmt_r1`: centre 0.010, spread 0.125
    - `msbb_tmt_phg`: centre 17.104, spread 0.202
    - `rosmap_tmt_r2`: centre -0.001, spread 0.151
- post-merge numeric-space check on `values_z`: **verified** - all strata occupy one numeric space at merge

**`lfq_brain_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- measured before merge: location range 0.298 pooled SD, spread ratio 1.725
    - `mayo_lfq_tcx`: centre 27.887, spread 0.423
    - `msbb_lfq_pfc`: centre 27.716, spread 0.730
- post-merge numeric-space check on `values`: **verified** - all strata occupy one numeric space at merge

**`srm_brain_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- measured before merge: location range 0.0 pooled SD, spread ratio 1.25
    - `rosmap_srm_panel1`: centre 0.000, spread 0.266
    - `rosmap_srm_panel2`: centre 0.000, spread 0.333
- post-merge numeric-space check on `values`: **verified** - all strata occupy one numeric space at merge

**`soma_plasma_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- post-merge numeric-space check on `values`: **verified** - single stratum, trivially one space

**`olink_plasma_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- post-merge numeric-space check on `values`: **verified** - single stratum, trivially one space

**`olink_csf_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- post-merge numeric-space check on `values`: **verified** - single stratum, trivially one space

**`dia_plasma_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- post-merge numeric-space check on `values`: **verified** - single stratum, trivially one space

**`dia_csf_proteomics`** - Z applied: **False**

- measured centres and spreads are comparable across strata, so the native values are kept exactly as delivered
- post-merge numeric-space check on `values`: **verified** - single stratum, trivially one space

## QC applied, per stratum

Source-specified QC is honoured before any transform, Z or merge. Checks marked *upstream* were already applied by the data producer and are verified rather than re-applied.

| stratum | check | applied | reference |
|---|---|---|---|
| `diversecohorts_dlpfc` | GIS/non-GIS separation | upstream | QC_process_Frontal.txt step 1 |
| `diversecohorts_dlpfc` | features <50% missing | upstream | QC_process_Frontal.txt step 2 |
| `diversecohorts_dlpfc` | ratio to total abundance, log2 | upstream | QC_process_Frontal.txt step 2 |
| `diversecohorts_dlpfc` | iterative PCA outlier removal (19 removed) | upstream | QC_process_Frontal.txt step 3 |
| `diversecohorts_dlpfc` | batch effect regressed out | upstream | QC_process_Frontal.txt step 4 |
| `diversecohorts_temporal` | GIS/non-GIS separation | upstream | QC_process_Temporal.txt step 1 |
| `diversecohorts_temporal` | features <50% missing | upstream | QC_process_Temporal.txt step 2 |
| `diversecohorts_temporal` | ratio to total abundance, log2 | upstream | QC_process_Temporal.txt step 2 |
| `diversecohorts_temporal` | iterative PCA outlier removal (2 removed) | upstream | QC_process_Temporal.txt step 3 |
| `diversecohorts_temporal` | batch effect regressed out | upstream | QC_process_Temporal.txt step 4 |
| `rosmap_tmt_r1` | GIS channels excluded | upstream | delivered C2 matrix is 8 channels/batch |
| `rosmap_tmt_r1` | TAMPOR median-polish batch correction | upstream | TAMPOR.R + batch correction PDF |
| `msbb_tmt_phg` | 8 low-quality samples removed (198 -> 190) | upstream | xlsx header note |
| `rosmap_tmt_r2` | 126 GIS reference channel excluded | upstream | PD export reports ratios vs 126, so the reference channel is not a sample |
| `rosmap_tmt_r2` | ratio to internal reference standard | upstream | PD Abundance Ratio (ch / 126) |
| `rosmap_tmt_r2` | log2 of the reference ratio | here | standard TMT practice |
| `rosmap_tmt_r2` | per-protein centring on median of per-plex medians | here | matches the round-1 product construction, so both ROSMAP rounds are comparable |
| `rosmap_tmt_r2` | features <50% missing | here | same threshold as every other stratum |
| `mayo_lfq_tcx` | reverse / contaminant / site-only rows removed | here | MaxQuant convention |
| `msbb_lfq_pfc` | reverse / contaminant / site-only rows removed | here | MaxQuant convention |
| `rosmap_srm_panel1` | control wells flagged | here | isControl / sample.type |
| `rosmap_srm_panel2` | control wells flagged | here | isControl / sample.type |
| `rosmap_soma_plasma` | ANML normalisation | upstream | deposit README (SomaLogic recommendation) |
| `rosmap_soma_plasma` | ColCheck per-analyte QC flagged | here | SomaScan protein metadata |
| `pdrd_olink_plasma` | D01<->D02 bridging | upstream | AMP-PD README |
| `pdrd_olink_plasma` | QC_Warning / Distribution_QC / Outliers_QC flagged | here | AMP-PD README |
| `pdrd_olink_plasma` | below-LOD fraction computed per assay | here | Olink recommendation: consider excluding assays with 25-50%+ below LOD |
| `pdrd_olink_csf` | D01<->D02 bridging | upstream | AMP-PD README |
| `pdrd_olink_csf` | QC_Warning / Distribution_QC / Outliers_QC flagged | here | AMP-PD README |
| `pdrd_olink_csf` | below-LOD fraction computed per assay | here | Olink recommendation: consider excluding assays with 25-50%+ below LOD |

> **QC described as "done" upstream was often not actually done.** Several source READMEs and QC logs state that a filter was applied, but the delivered matrix does not satisfy it. Those filters were re-applied here rather than taken on trust. The clearest case: `QC_process_Temporal.txt` states *"only keep proteins that have missing data in <50% of the samples"*, yet the delivered `n278_residual_log2_batch.TCX.csv` contains 398 features below that threshold, some as low as 10% present - most likely because the batch regression that runs *after* that step reintroduces missingness. Re-applying the filter here removes them. Treat every upstream QC claim in the table above as verified-by-us, not as inherited.

### Standing QC policy

1. QC decisions are made **per stratum, before** any transform, Z or merge.
2. Samples and features that fail QC are **flagged, never silently deleted**. The exclusion threshold belongs to the study, not to this build.
3. Every threshold lives in `config.py`, never inline.
4. Every exclusion is counted in the audit.
5. Every `SOURCE_QC` entry is **checked**, not recited: a step the producer says they applied is verified against the delivered file where the data can show it, and a step they left to us fails the audit unless the build left a trace of actually doing it.

**Deletion is reserved for rows whose identity is unusable, never their quality.** There are exactly three such cases in this build, and each is recorded as data rather than falling out silently as a join miss:

- **Olink plasma, 172 samples** on two plates whose donor identity is in dispute pending AMP-PD. Withdrawn by *plate* - the physical fact - so the rule survives a re-release. 39 participants lose every plasma sample.
- **Samples with no donor at all in the source metadata.** They cannot be joined, de-duplicated against a repeat visit, or withdrawn if that donor revokes consent.
- **DIA plasma, low-completeness runs** - the one *quality*-based removal, and it exists only because that deposit ships no sample-level QC of its own. The failures are unambiguous: one run detected 10 of 202 proteins.

### Feature filtering

Drop rates differ substantially between strata, and the difference is about *what processing stage the deposit sits at*, not about inconsistent QC - the same threshold is applied identically everywhere. Strata delivered as finished, already-filtered products (ROSMAP TMT round 1, DiverseCohorts DLPFC, SomaScan, Olink) lose almost nothing. Strata delivered closer to the instrument lose more: ROSMAP TMT round 2 is a raw Proteome Discoverer export and drops 22.9%, and the two MaxQuant label-free sets drop 34.5% and 40.8%, because a protein quantified in only a few runs is genuinely absent from most samples.

The missingness filter is applied **per stratum, never after merge** - each stratum has its own sample set, so a protein measured well by one assay would otherwise look mostly-missing across a merged layer and be dropped for the wrong reason. Threshold follows `QC_process_*.txt`: keep features present in >= 50% of that stratum's samples.

| stratum | unresolved id | below 50% present | duplicates collapsed | final |
|---|---|---|---|---|
| `diversecohorts_dlpfc` | 0 | 2 | 0 | 9,178 |
| `diversecohorts_temporal` | 0 | 398 | 0 | 9,336 |
| `rosmap_tmt_r1` | 0 | 1 | 0 | 8,816 |
| `msbb_tmt_phg` | 0 | 0 | 1,205 | 10,942 |
| `rosmap_tmt_r2` | 0 | 2,452 | 0 | 8,237 |
| `mayo_lfq_tcx` | 2 | 2,111 | 71 | 3,927 |
| `msbb_lfq_pfc` | 1 | 2,351 | 50 | 3,360 |
| `rosmap_srm_panel1` | 0 | 0 | 0 | 119 |
| `rosmap_srm_panel2` | 0 | 0 | 0 | 67 |
| `rosmap_soma_plasma` | 0 | 0 | 0 | 7,289 |
| `pdrd_olink_plasma` | 0 | 0 | 0 | 1,463 |
| `pdrd_olink_csf` | 0 | 0 | 0 | 1,463 |
| `pdrd_dia_plasma` | 0 | 0 | 0 | 202 |
| `pdrd_dia_csf` | 0 | 0 | 0 | 249 |

### Protein sparsity, per stratum

Measured on what actually ships, after filtering. Sparsity is low but never zero - some patchiness across strata is expected and is not a defect.

| stratum | % missing | features fully observed | median feature presence | min |
|---|---|---|---|---|
| `diversecohorts_dlpfc` | 5.95% | 63.2% | 100.0% | 50.1% |
| `diversecohorts_temporal` | 4.92% | 76.5% | 100.0% | 51.4% |
| `rosmap_tmt_r1` | 7.29% | 59.1% | 100.0% | 51.7% |
| `msbb_tmt_phg` | 0.0% | 100.0% | 100.0% | 100.0% |
| `rosmap_tmt_r2` | 6.8% | 71.3% | 100.0% | 50.0% |
| `mayo_lfq_tcx` | 7.93% | 48.2% | 99.6% | 50.0% |
| `msbb_lfq_pfc` | 8.12% | 45.8% | 99.7% | 50.0% |
| `rosmap_srm_panel1` | 0.1% | 38.7% | 99.9% | 99.0% |
| `rosmap_srm_panel2` | 0.3% | 79.1% | 100.0% | 84.3% |
| `rosmap_soma_plasma` | 0.0% | 100.0% | 100.0% | 100.0% |
| `pdrd_olink_plasma` | 0.07% | 97.1% | 100.0% | 53.8% |
| `pdrd_olink_csf` | 0.0% | 98.0% | 100.0% | 99.8% |
| `pdrd_dia_plasma` | 11.68% | 21.8% | 97.9% | 52.2% |
| `pdrd_dia_csf` | 0.03% | 77.1% | 100.0% | 99.1% |

Note the *layer* matrices are sparser than these per-stratum figures, because a layer takes the union of its strata's features: a protein measured by one stratum and not another is structurally absent for the other's rows. That is recorded in `missing_reason` as code 2 (`not_assayed_in_this_stratum`) and is not a QC problem.

#### ROSMAP round 2: missingness is plex-structured, not per-sample

Round 2 is a raw Proteome Discoverer multiconsensus export, and its missingness behaves differently from every other stratum. Measured on the delivered file:

- **0 features are partially present within a plex** - i.e. none. Missingness is entirely all-or-nothing at the plex level.
- a feature is present in 10.67 of 14 plexes on average
- 166 features are empty in every channel

This is **identification** missingness from the consensus step, not detection missingness in a sample: a protein identified in one plex and not another produces a clean block of absent channels. Two consequences worth carrying into any analysis:

1. The >=50% presence filter is effectively selecting proteins identified in at least 7 of 14 plexes, not proteins detected in at least half the donors. It is the same threshold used everywhere else, but it means something different here.
2. Because whole plexes drop out together, remaining missingness is correlated with batch. Anything that treats missingness as random - naive imputation, complete-case filtering - will inherit that structure. Round 2 is also the stratum most likely to carry **stale** identifications from the multiconsensus search, so treat its low-presence tail with more suspicion than the others.

## Identifiers

- **All identifiers are strings**, exactly as the source defines them. Prefixes (`PD-`, `PP-`, `R`, `AMPAD_MSSM_`) are part of the identifier and are never stripped, so IDs join to future metadata drops unchanged.
- `person_id` **is** the source identifier. No surrogate is invented.
- Because two studies can reuse the same digits, `person_source_namespace` records the identifier system and `person_global_key` (`namespace:id`) is the unambiguous cross-study join key.

| namespace | strata |
|---|---|
| `AMP-AD.individualID` | `diversecohorts_dlpfc`, `diversecohorts_temporal`, `msbb_tmt_phg`, `mayo_lfq_tcx`, `msbb_lfq_pfc` |
| `AMP-PD.participant_id` | `pdrd_olink_plasma`, `pdrd_olink_csf` |
| `ROSMAP.individualID` | `rosmap_soma_plasma` |
| `ROSMAP.projid` | `rosmap_tmt_r1`, `rosmap_tmt_r2`, `rosmap_srm_panel1`, `rosmap_srm_panel2` |

> ROSMAP appears under two namespaces: `ROSMAP.projid` (TMT, SRM, and the `sd` subset of DiverseCohorts) and `ROSMAP.individualID` (SomaScan). The bridge between them is shipped in this file as **`person_crosswalk_rosmap`** (from syn3191087), so ROSMAP plasma *can* be joined to ROSMAP brain by person.

### Identifier detective work

Three identifier problems were not solvable from the obvious column and needed tracing back through the data. Each is recorded here because the resolution is not self-evident from the files, and each fails *silently* - returning zero matches or an unlinked donor rather than an error.

**1. DiverseCohorts `sd` had no `individualID` at all.**

The frontal traits sheet has an `individualID` column that is populated for the `emdp`, `mayo` and `mssm` contributing datasets and **completely empty for `sd`** - 400 rows, zero values, 371 of which survive into the delivered matrix. Read naively this looks like 371 samples with no donor, and an earlier pass of this build mislabelled them as reference pools.

They are neither. The sheet carries a *second* identifier column, `projID`, which is populated for exactly the `sd` rows and empty for every other dataset. The values are 8-digit ROS/MAP project IDs, so `sd` is a ROS/MAP contribution to DiverseCohorts carrying ROSMAP identifiers rather than AMP-AD ones. Falling back to `projID` where `individualID` is absent recovers all 371, and takes this stratum from 375 unresolved to **0**.

Consequence: `diversecohorts_dlpfc` spans **two identifier systems in one stratum**, so `person_source_namespace` is assigned per row (715 `AMP-AD.individualID`, 371 `ROSMAP.projid`) rather than per stratum. Anything that assumes one namespace per stratum will mis-join this layer.

**2. ROSMAP uses two identifier systems across assays.**

ROSMAP TMT and SRM key on numeric `projid`; ROSMAP SomaScan keys on `R`-prefixed `individualID`. They never collide, so nothing errors - the two simply never join, and plasma silently fails to link to brain for the same donor. The bridge is `syn3191087`, published as `ROSMAP_clinical.csv`, which is written into this file as `person_crosswalk_rosmap` (identifier columns only; that source is otherwise clinical and is dropped at read).

**3. Sample IDs were re-ordered between matrix and metadata.** The MSBB label-free run labels put the plate position before the specimen number (`b<batch>_<position>_<specimen>`) while its biospecimen metadata puts the specimen first (`b<batch>_<specimen>_<position>`) - the last two tokens swapped. The batch-7 rerun compounds it: `b7r2_*` in the matrix against `r2b7_*` in the metadata. Mayo appends a position token the metadata omits (`mayo_b<batch>_<specimen>_<position>` against `b<batch>_<specimen>`). None of these join without an explicit rewrite, and all three silently return zero matches if attempted directly.

### Row keys

A person legitimately repeats within a layer - multiple brain regions, multiple timepoints, technical replicates. `row_key` is composed from the real discriminators in order (person -> anatomic site -> visit -> stratum) and only falls back to a replicate counter when those do not separate the rows. Every component is also its own column, so no one has to parse the key.

## Missing values

`values` and `values_z` are the analysis-ready tables and carry **exactly one null type: float `NaN`**. No sentinels, no second null flavour - numpy, pandas and sklearn all behave normally.

The *reason* a cell is missing lives in a separate `missing_reason` int8 array of the same shape, which is diagnostic only:

| code | meaning |
|---|---|
| 0 | observed |
| 1 | not_detected_in_sample (assayed in this stratum) |
| 2 | not_assayed_in_this_stratum (column comes from another stratum) |
| 3 | failed_qc_in_this_stratum (present in source, dropped by QC) |

### `below_lod` is three-state, not boolean

On the Olink layers, `below_lod` carries **three** states rather than two:

| value | meaning |
|---|---|
| `0.0` | measured, at or above the limit of detection |
| `1.0` | measured, below the limit of detection |
| `NaN` | **not measured** - there is no value to compare against the limit |

This matters because it changes results rather than presentation. `below_lod` was previously written as `(values < lod)` cast to float, and `NaN < NaN` is `False`, which casts to `0.0` - so a never-measured cell was stored identically to one measured comfortably above the limit. The conventional filter `below_lod == 0.0` therefore **silently retained every never-measured cell**, and the only correct predicate was `below_lod == 0.0 AND missing_reason == 0`, which is not discoverable from the surface. Restoring `NaN` makes the three states separable with the single null flavour the file already documents.

## The time axis: `visit_month`

`visit_month` - **months since that participant's baseline** - is the canonical cross-layer time axis and is present in **every** layer: populated where the dataset is longitudinal, null where it is postmortem. One consumer expression now aligns Olink, SomaScan and DIA, and adding a longitudinal dataset means populating `visit_month` rather than inventing another axis.

`visit_index` stays beside it as the source-native ordinal and is **never** the join key. `visit_month_evidence` records how the value was reached:

| evidence | meaning | layers |
|---|---|---|
| `source` | the deposit states months since baseline | Olink plasma/CSF, DIA plasma/CSF |
| `derived_from_visit_index` | converted from an ordinal follow-up number | SomaScan plasma |

> **The SomaScan conversion is an approximation and is labelled as one.** ROSMAP's `Visit` is the follow-up *year*, so `visit_month = Visit x 12` is a structural conversion of an ordinal index - it never touches `age_at_visit` and never trips the scope gate. But ROSMAP follow-ups are nominally annual and not exactly so, which means these values are **not** measured intervals and must not be treated as such. `visit_index` is retained so the conversion is reversible.

Before this, `soma_plasma` carried only an ordinal and the brain layers carried nothing, so a consumer joining ROSMAP plasma to PDRD Olink on `person_id` had no shared time column and nothing longitudinal could cross cohorts.

## Post-mortem interval - a deliberate scope change

> **Scope change, recorded rather than made silently.** Post-mortem interval was previously denied as clinical phenotype. It is now carried as a **biospecimen technical CDE**: it describes the specimen, not the donor, and it is a standard covariate for brain proteomics. `pmi` and `postmortem_interval` were removed from the scope gate's deny list for this reason and this reason only. That narrows HARMONIZATION_PLAN section 0 and REQUIREMENTS SEC-1, both of which are written against the DUC, so it is a change to the disclosure boundary and is documented as one here, in `config.py` at `CLINICAL_DENY`, and in both markdowns. Nothing else was re-admitted, and SEC-2's read-time allow-list is unchanged - only the join key and the `pmi` column are ever loaded from the clinical file.

**One scale: `pmi_hours`, hours with decimals, always.** Sub-hour precision is real - ROSMAP's median is 6.783 h - so rounding to whole hours would discard resolution on a covariate whose entire use is fine-grained. A bare `pmi` column is never trusted, because the unit differs between studies by convention: MSBB records **minutes**, ROSMAP **hours**.

### Conversion ledger

| stratum | source | column | unit as found | stated or inferred | factor | populated | median (h) |
|---|---|---|---|---|---|---|---|
| `diversecohorts_dlpfc` | - | - | - | - | - | 0/1082 | - |
| `diversecohorts_temporal` | - | - | - | - | - | 0/278 | - |
| `rosmap_tmt_r1` | ROSMAP_clinical.csv | pmi | hours | stated | 1.0 | 378/400 | 6.5 |
| `msbb_tmt_phg` | - | - | - | - | - | 0/190 | - |
| `rosmap_tmt_r2` | ROSMAP_clinical.csv | pmi | hours | stated | 1.0 | 197/210 | 6.55 |
| `mayo_lfq_tcx` | - | - | - | - | - | 0/230 | - |
| `msbb_lfq_pfc` | - | - | - | - | - | 0/306 | - |
| `rosmap_srm_panel1` | ROSMAP_clinical.csv | pmi | hours | stated | 1.0 | 1186/1364 | 6.583 |
| `rosmap_srm_panel2` | ROSMAP_clinical.csv | pmi | hours | stated | 1.0 | 1186/1340 | 6.583 |
| `rosmap_soma_plasma` | - | - | - | - | - | 0/970 | - |
| `pdrd_olink_plasma` | - | - | - | - | - | 0/1380 | - |
| `pdrd_olink_csf` | - | - | - | - | - | 0/1320 | - |
| `pdrd_dia_plasma` | - | - | - | - | - | 0/1459 | - |
| `pdrd_dia_csf` | - | - | - | - | - | 0/2807 | - |

`pmi_evidence` distinguishes **two different kinds of null**, per row, because the difference matters to anyone modelling with it:

| evidence | meaning |
|---|---|
| `converted_from_stated_unit` | the source named its unit; the conversion is exact and reversible |
| `inferred_from_distribution` | it did not; the unit was read off the observed median against a plausible band |
| `not_applicable_antemortem` | **a fact about the specimen** - a living donor's plasma or CSF has no post-mortem interval |
| `source_not_acquired` | **a gap in what we hold** - the object exists, we do not have it |

`pmi_hours` is **never zero** where the value is unknown. A zero interval is a claim, not a null.

Where inference is used it is safe because the bands do not overlap: PMI is the interval between death and processing at a brain bank, so the plausible range is 0.5-120.0 hours, and the same durations recorded in minutes land near 400. Where the median sits in neither band the values stay null and the stratum is reported undetermined - never guessed. The audit checks the post-conversion median against that band and **reports** an out-of-band result; it never auto-corrects.

> **For the CDE repository.** These three columns are new to the specimen surface and are not yet in the external CDE dictionary: **`pmi_hours`** (float, hours with decimals), **`pmi_evidence`** (the four-value enum above) and **`pmi_source_unit`** (the unit as found in the source). They need registering there, and the unit convention needs registering with them - a bare `pmi` CDE with no unit field is unsafe across studies, because MSBB's minutes and ROSMAP's hours differ by 60x on a column that looks identical. Any external repository already carrying a `pmi` CDE should be reviewed against this.

## Brodmann area: measurement and inference are separable

`brodmann_area` is **partly inferred**, and `brodmann_area_evidence` says which rows are which. Anything consuming `brodmann_area` must read the evidence column with it.

| evidence | meaning |
|---|---|
| `source` | the specimen's own annotation stated it |
| `propagated` | another specimen from the same donor and tissue stated it |
| `convention` | neither did; it was asserted from the anatomic site |

| stratum | site | convention | source | propagated | convention-filled | disagreements kept |
|---|---|---|---|---|---|---|
| `diversecohorts_dlpfc` | dorsolateral prefrontal cortex | BA9 | 0 | 190 | 892 | 67 |
| `diversecohorts_temporal` | superior temporal gyrus | BA22 | 0 | 0 | 278 | 0 |
| `rosmap_tmt_r1` | dorsolateral prefrontal cortex | BA9 | 0 | 0 | 400 | 0 |
| `msbb_tmt_phg` | parahippocampal gyrus | - | 190 | 0 | 0 | 0 |
| `rosmap_tmt_r2` | dorsolateral prefrontal cortex | BA9 | 0 | 0 | 210 | 0 |
| `msbb_lfq_pfc` | prefrontal cortex | - | 284 | 0 | 0 | 0 |
| `rosmap_srm_panel1` | dorsolateral prefrontal cortex | BA9 | 0 | 0 | 1364 | 0 |
| `rosmap_srm_panel2` | dorsolateral prefrontal cortex | BA9 | 0 | 0 | 1340 | 0 |

The conventions applied are `dorsolateral prefrontal cortex` -> `BA9`, `superior temporal gyrus` -> `BA22`, declared once in `config.py` so a consumer can undo them.

**A stated value always wins.** The convention only fills where the source is silent, and never overwrites - which is why the 67 DiverseCohorts rows that state `BA10` are still `BA10` rather than being swept into the DLPFC convention's `BA9`. Those 67 rows are the standing evidence that the convention is a useful default and not a universal truth.

> **For the CDM side.** `PROTEOMICS_REVIEW.md` gives `brodmann_area` precedence over the free-text site label. A loader that honours that without reading `brodmann_area_evidence` would move every convention-filled DLPFC row off `4195656 Region of frontal cortex` onto a BA9 concept **on our inference alone**. Separately, the DiverseCohorts temporal strata move off `4193043 Region of temporal cortex` to a superior-temporal-gyrus concept, because the source metadata says `superior temporal gyrus` for all 280 temporal TMT specimens and the build's earlier `temporal cortex` claim was wrong.

## SysBio CDE coverage

Evidence grades: **A** derivable from the data files, **B** transcribed from a bundled README, **C** unavailable. Read from the built artifact, not from the config declaration - where an assay-metadata object states the platform per sample, that is what ships and the grade is **A**.

| stratum | platform | evidence | analysis pipeline | evidence | source |
|---|---|---|---|---|---|
| `diversecohorts_dlpfc` | mixed per sample: Eclipse (205/1086); Exploris 240 (382/1086); Lumos (499/1086) | A | FragPipe (TMT 18-plex) | A | `AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata_260622.csv.gz` |
| `diversecohorts_temporal` | Lumos | A | FragPipe (TMT 18-plex) | A | `AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata_260622.csv.gz` |
| `rosmap_tmt_r1` | OrbiTrap Fusion | A | Proteome Discoverer 2.3.0.522 + TAMPOR | A | `ROSMAP_assay_proteomics_TMTquantitation_metadata.csv.gz` |
| `msbb_tmt_phg` | Q Exactive HF | A | - | C | `MSBB_assay_TMT_metadata.csv.gz` |
| `rosmap_tmt_r2` | Q Exactive HF-X Orbitrap | A | Proteome Discoverer | A | config |
| `mayo_lfq_tcx` | Q Exactive Plus | A | MaxQuant (version not recorded) | A | `MayoRNAseq_assay_proteomics_metadata.csv.gz` |
| `msbb_lfq_pfc` | Q Extrative Plus | A | MaxQuant (version not recorded) | A | `MSBB_assay_proteomics_metadata.csv.gz` |
| `rosmap_srm_panel1` | triple quadrupole (instrument alias 'Smeagol') | partial | - | C | config |
| `rosmap_srm_panel2` | triple quadrupole (instrument alias 'Smeagol') | partial | - | C | config |
| `rosmap_soma_plasma` | SOMAmer aptamer array | B | SomaLogic ANML normalisation | B | config |
| `pdrd_olink_plasma` | PEA with NGS readout | B | Olink NPX + Olink Analyze bridging | B | config |
| `pdrd_olink_csf` | PEA with NGS readout | B | Olink NPX + Olink Analyze bridging | B | config |
| `pdrd_dia_plasma` | Orbitrap Exploris 480 and TripleTOF 6600 (fractions merged upstream) | C | fragment -> peptide -> protein roll-up; batch correction applied upstream by AMP-PD | B | config |
| `pdrd_dia_csf` | Orbitrap Exploris 480 | B | fragment -> peptide -> protein roll-up; batch correction applied upstream by AMP-PD | B | config |

`platform` is genuinely **per sample**, not per stratum: the DiverseCohorts strata span three instruments and the build records that rather than collapsing it to one name. The per-sample `platform` column is authoritative; the table above is its stratum-level summary. ROSMAP SRM is the one stratum that does **not** close - `syn23569441` carries `platform = NA` on all 1,212 of its rows, which is a silence in the deposit rather than a gap in this build, and is reported as such.

FILES CDEs `current_version`, `created_on`, `modified_on` and `drs_id` are left **null** rather than guessed. Every filesystem mtime in this directory is the download timestamp, not file creation, so using it would be quietly wrong.

## Assertion audit

313 checks: **242 pass**, **18 fail**, 53 unverifiable.

Every claim - in the tracking file, in a bundled README, or encoded in this codebase - is re-checked against the data. Disagreements are reported, not silently resolved.

| check | subject | claim | source | observed |
|---|---|---|---|---|
| unresolved_donor_dropped | diversecohorts_dlpfc | a sample the source cannot resolve to a donor is removed | computed | 4 sample(s) dropped. no donor identifier in the source metadata; identity unusable. Rule: non-pool row whose person_source_value is null or an absent-value sentinel, after the source annotation has been joined. Specimen ids are withheld from this artifact (SEC-8) - they embed a participant identifier, and the rule reproduces the set exactly |
| source_qc | diversecohorts_dlpfc | upstream: features <50% missing | QC_process_Frontal.txt step 2 | 2 features in the delivered file are below 50% present, so the claimed upstream filter does not hold on what was shipped; re-applied here |
| source_qc | diversecohorts_temporal | upstream: features <50% missing | QC_process_Temporal.txt step 2 | 398 features in the delivered file are below 50% present, so the claimed upstream filter does not hold on what was shipped; re-applied here |
| cde_conflict | rosmap_tmt_r1 | config declares platform = 'Orbitrap FTMS; HCD MS2 (model not recorded)' | config | assay metadata says {'OrbiTrap Fusion': 400} - reported, not resolved |
| unresolved_donor_dropped | rosmap_soma_plasma | a sample the source cannot resolve to a donor is removed | computed | 3 sample(s) dropped. no donor identifier in the source metadata; identity unusable. Rule: non-pool row whose person_source_value is null or an absent-value sentinel, after the source annotation has been joined. Specimen ids are withheld from this artifact (SEC-8) - they embed a participant identifier, and the rule reproduces the set exactly |
| samples_withdrawn | pdrd_olink_plasma | samples withdrawn for disputed identity | computed | 172 of 1552 samples (11.1%) withdrawn; 39 participants lose every sample in this stratum. Rule: PlateID in BIOREP_pl1_Samplesheet, BIOREP_pl2_Samplesheet. donor identity in dispute pending AMP-PD (plate swap, D-D) |
| qc_sample_completeness | pdrd_dia_plasma | samples detecting < 50% of the layer's features are dropped | config | 12 sample(s) dropped; they detected [10, 70, 72, 72, 72] of 202 features (lowest first) - sample ids withheld, they embed a participant identifier (SEC-8) |
| tracking | AMP AD / TMT LC/MS / AMP-AD_DiverseCohorts | Number Individuals = 850? | tracking file | derived 1218 (diversecohorts_dlpfc=976; diversecohorts_temporal=242) |
| tracking | AMP AD / SomaScan / ROSMAP | Number Individuals = 610? | tracking file | derived 887 (rosmap_soma_plasma=887) |
| tracking | AMP AD / SomaScan / ROSMAP | Biospecimen source = Brain - dorsolateral prefrontal cortex | tracking file | derived blood (derived from the biospecimen metadata, not from the matrix filename) |
| tracking | AMP AD / LC-MS / MayoRNAseq | Number Individuals = 230 | tracking file | derived 199 (mayo_lfq_tcx=199) |
| tracking | AMP AD / LC-MS / MSBB | Number Individuals = 308 | tracking file | derived 266 (msbb_lfq_pfc=266) |
| tracking | AMP AD / TMT LC/MS / MSBB | Number Individuals = 800 | tracking file | derived 185 (msbb_tmt_phg=185) |
| tracking | AMP PDRD / Olink / Unified Cohorts | Number Individuals = 375? | tracking file | derived 772 (pdrd_olink_plasma=374; pdrd_olink_csf=398) |
| tracking | AMP PDRD / Olink / Unified Cohorts | Biospecimen source = Plasma, CSF | tracking file | derived blood, cerebrospinal fluid (derived from the biospecimen metadata, not from the matrix filename) |
| tracking | AMP PDRD / LC-MS / Unified Cohorts | Number Individuals = 800? | tracking file | derived 937 (pdrd_dia_csf=621; pdrd_dia_plasma=316) |
| tracking | AMP PDRD / LC-MS / Unified Cohorts | Biospecimen source = Plasma, CSF | tracking file | derived blood, cerebrospinal fluid (derived from the biospecimen metadata, not from the matrix filename) |
| tracking | all rows | metadata syn IDs = 15 objects referenced | tracking file | derived 13 present (presence detected by content hash, not by path) |

Column naming is checked across layers: every layer exposes the same CDE and link columns under the same names, with layer-specific proteomics columns (`panel`, `dilution`, `seq_id`, `peptide_sequence`, `lod`) present only where the assay actually has them.

Full detail: `build/inventory/assertion_audit.tsv`.

### Known limitations

- **PMI is not held for five strata.** MSBB and Mayo PMI live in individual-metadata objects (`syn73713767`, `syn73713766`) that are not in the delivery, and DiverseCohorts has **no individual-metadata object registered anywhere**, so its PMI source is still unidentified. Those rows are `source_not_acquired` - a gap in what we hold, distinct from the antemortem nulls. MSBB is also the study that records minutes, so the unit conversion stays untested against real data until that object arrives.
- **Brodmann area is ~88% inference on the DiverseCohorts strata.** No `TMT quantitation` specimen carries a source `BrodmannArea` at all; cross-assay propagation reaches about 12% and the rest is convention. Separable via `brodmann_area_evidence`, but it is not measurement.
- **ROSMAP SRM has no platform.** `syn23569441` carries `platform = NA` on all 1,212 rows. This is a silence in the deposit, not a gap in the build.
- **The SomaScan normalisation label has no evidence behind it.** The ANML and non-ANML deposits are byte-identical (md5 `f8a3d03ecd1b57cbb409bc6b054a7890`), so one of the two is mislabelled and the data cannot say which. Values are unaffected; only the label is unsupported. Raised with ADKP.
- **The Olink plasma plate swap is unresolved at source.** 172 samples are withdrawn pending AMP-PD's answer, not because the answer is known.
- **ROSMAP round 2 carries plex-structured missingness** (see above) and is the stratum most likely to contain stale multiconsensus identifications.
- **The FILES CDEs `current_version` / `created_on` / `modified_on` / `drs_id` are null throughout**, deliberately: every filesystem mtime here is a download timestamp, and Synapse versions were not downloaded.
- **The DIA plasma fraction is not recoverable.** The delivered batch-corrected product merges native and depleted plasma under one `-PLA-` sample token, so which fraction a sample came from cannot be read off it. It is therefore one stratum, not two.

## What is in the file

```
sysbio_proteomics-09042026.h5
|-- tmt_brain_proteomics/
|-- lfq_brain_proteomics/
|-- srm_brain_proteomics/
|-- soma_plasma_proteomics/
|-- olink_plasma_proteomics/
|-- olink_csf_proteomics/
|-- dia_plasma_proteomics/
|-- dia_csf_proteomics/
+-- person_crosswalk_rosmap/    ROSMAP projid <-> individualID (syn3191087)
```

Every populated layer has the same internal structure:

```
<layer>/
  values            (n_rows, n_proteins) float32   the analysis-ready table
  values_z          same shape, only where the layer needed Z
  lod, below_lod    same shape, Olink layers only
  missing_reason    same shape, int8 diagnostic (see Missing values)
  row_key           (n_rows,)     person + visit, the analysis key
  protein_columns   (n_proteins,) column headers for `values`
  samples/<col>     one dataset per CDE / link / technical column
  features/<col>    one dataset per feature-description column
  assay/<col>       SysBio ASSAY CDEs, one row per stratum
  files/<col>       SysBio FILES CDEs, one row per source file
```

Layer-level attributes carry the scale policy, the per-stratum transform and Z flags, the post-merge numeric-space check, and provenance.

## The analysis-ready bundles (`build/tabs/`)

Each layer is also written as a standalone bundle - **wide only**; the long form is deliberately not emitted, because pooled brain TMT alone would be ~25M rows and the wide matrix is what the HDF5 already holds.

```
tabs/<layer_key>/
  <layer_key>_matrix.parquet     feature_id + one column per specimen
  <layer_key>_samples.parquet    one row per sample, the full CDE surface
  <layer_key>_features.parquet   one row per feature, the identity half
  <layer_key>_bundle.json        which surface the matrix carries, and how to join
```

Parquet is primary because it preserves dtypes and the three-state `below_lod` distinction above; a `.tsv.gz` mirror of each table ships beside it.

| layer | matrix | value surface | columns keyed by | samples | features |
|---|---|---|---|---|---|
| `tmt_brain_proteomics` | 14,916 x 2,160 | `values_z` | `specimen_id` | 2,160 | 14,916 |
| `lfq_brain_proteomics` | 4,437 x 536 | `values` | `specimen_id` | 536 | 4,437 |
| `srm_brain_proteomics` | 185 x 2,704 | `values` | `specimen_id` | 2,704 | 185 |
| `soma_plasma_proteomics` | 7,289 x 970 | `values` | `specimen_id` | 970 | 7,289 |
| `olink_plasma_proteomics` | 1,463 x 1,380 | `values` | `specimen_id` | 1,380 | 1,463 |
| `olink_csf_proteomics` | 1,463 x 1,320 | `values` | `specimen_id` | 1,320 | 1,463 |
| `dia_plasma_proteomics` | 202 x 1,459 | `values` | `specimen_id` | 1,459 | 202 |
| `dia_csf_proteomics` | 249 x 2,807 | `values` | `specimen_id` | 2,807 | 249 |

The matrix is written in the same pass that writes the HDF5's `values`, from the same in-memory array, and the tables in the same pass that writes the HDF5's tables. A standalone HDF5-to-parquet converter would be simpler and would give up exactly the property that matters: **the two cannot drift**.

### The AMP return path (`build/returns/`)

One directory per grant, holding that grant's specimen set for every layer it appears in, plus a manifest that **demonstrates** the boundary rather than asserting it - it names the grants present in the build, the grant returned, and the count of rows from any other grant that reached the slice, which must be zero.

| grant | specimens | layers | rows from another grant |
|---|---|---|---|
| AMP AD | 6,370 | 4 | 0 |
| AMP PDRD | 6,966 | 4 | 0 |

This is computed on every build rather than at release time: returning one AMP programme's donors to another is not recoverable after the fact.

> **Controlled access.** `tabs/` and `returns/` carry the same donor-level identifiers as the HDF5. They are controlled-access artifacts, not a publishable by-product, and `.gitignore` denies them by name as well as by directory.

## Files produced

- `build/sysbio_proteomics-09042026.h5` - the harmonised dataset
- `build/build_report.json` - machine-readable build summary
- `build/tabs/<layer>/` - the analysis-ready bundle per layer (**controlled access**)
- `build/returns/<grant>/` - the per-grant return path with its boundary manifest (**controlled access**)
- `build/inventory/derived_inventory.tsv` - one row per stratum, computed
- `build/inventory/tracking_file_diff.tsv` - where the data disagrees with the tracking file
- `build/inventory/missing_objects.tsv` - referenced Synapse objects and what their absence blocks
- `build/inventory/assertion_audit.tsv` - every claim, checked
- `build/inventory/file_manifest.tsv` - every source file with sha256

Sources read through an explicit allow-list and never ingested as data:

- `ROSMAP_clinical.csv` (syn3191087) - **identifiers plus `pmi`**. This file is overwhelmingly clinical (`msex`, `educ`, `race`, `apoe_genotype`, `age_death`, `braaksc`, `ceradsc`, `cogdx`, `mmse`); only `projid` / `individualID` / `Study` for the crosswalk, and `projid` + `pmi` for the specimen CDE, are ever loaded. The `pmi` read is the scope change documented above and is the single column re-admitted.
- `MSBB_biospecimen_metadata.csv` (syn21893059) and `MayoRNAseq_biospecimen_metadata.csv` (syn20827192) - identifiers and biospecimen descriptors; `samplingAge` / `samplingAgeUnits` are dropped at read.
- The five assay-metadata objects - platform, reagent lot, control flags and acquisition descriptors only.

## Reading a layer

```python
import h5py, numpy as np, pandas as pd

h5 = h5py.File("build/sysbio_proteomics-09042026.h5", "r")
g  = h5["olink_csf_proteomics"]

proteins = [s.decode() for s in g["protein_columns"][:]]
rows     = [s.decode() for s in g["row_key"][:]]
X        = pd.DataFrame(g["values"][:], index=rows, columns=proteins)

# CDE / link columns for the same rows
# CDE / link columns for the same rows (string datasets come back as bytes)
def _col(ds):
    v = ds[:]
    return [x.decode() for x in v] if v.dtype.kind in "OS" else v

meta = pd.DataFrame({k: _col(g["samples"][k]) for k in g["samples"]}, index=rows)

# join ROSMAP plasma to ROSMAP brain by donor
cw = pd.DataFrame({k: _col(h5["person_crosswalk_rosmap"][k])
                   for k in h5["person_crosswalk_rosmap"]})
```
