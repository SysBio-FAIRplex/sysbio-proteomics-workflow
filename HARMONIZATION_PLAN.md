# AMP Proteomics Harmonization Plan (Priority-2 datasets)

Three layers — **brain**, **plasma**, **CSF** — each holding the proteomics data files for
every assay run on that matrix, described by the **SysBio extension CDEs** (ASSAY, FILES,
ASSAY_TO_SPECIMEN, ASSAY_INPUT_FILE — the green tables of Figure 2), linked through
**`person_id`** as the central key, with **visit** carrying the longitudinal axis.

Everything else links from elsewhere.

Status: **built.** This document is the original plan, written pre-implementation; the build has
since run (2026-08-22) and been externally audited (2026-09-02). Where this document and the build
disagree, [AUDIT_RESPONSE_PLAN.md](AUDIT_RESPONSE_PLAN.md) is current and
[REQUIREMENTS.md](REQUIREMENTS.md) states what must hold. Sections corrected against the built
artifact are marked **[updated]** inline; §9's blockers and decisions are all resolved and their
outcomes are recorded there.

---

## 0. Scope boundary

These are proteomics artifacts, not phenotype artifacts.

**In:** the measurements; the SysBio extension CDEs; `person_id`; visit; biospecimen and
assay descriptors; technical covariates; provenance.

**Out:** clinical and demographic attributes — diagnosis, cognition, Braak/CERAD, APOE,
education, age, sex, race, ethnicity. These join from the OMOP side under whatever access
controls apply.

**[updated]** **Post-mortem interval is IN, as a CDE** (decision D-J). PMI describes the specimen,
not the donor, and is a standard technical covariate for brain proteomics, so it is carried as a
biospecimen attribute and removed from the clinical deny-list. This narrows the boundary above
deliberately; see REQUIREMENTS SEC-1 and FR-15.

It is reported on **one scale — `pmi_hours`, hours with decimals** — because units differ by study
(MSBB minutes, ROSMAP hours). Values reach that scale by conversion from a stated unit or by
inference from the distribution, and which one applied is recorded per row. The fluid layers
(SomaScan plasma, PDRD Olink CSF/plasma, PDRD DIA) are antemortem samples from living donors, so
PMI is *not applicable* there rather than missing — the two nulls are encoded differently. One
consequence: `syn73713766` and `syn73713767`, listed in §1.2 as off the critical path because they
are clinical, are now on it **for `pmi` / `pmiUnits` only**.

Several annotation files we depend on for `person_id` resolution also carry clinical payload.
We take the identifier and technical columns and discard the rest **at read time**, so clinical
values never enter an in-memory harmonized object:

| File | Keep | Drop |
|---|---|---|
| ROSMAP R2 `..._Traits_FINAL.xlsx` | `ProjID`, `Channel`, `Plex`, `Lot #` | `cogdx_3lev`, `msex`, `braaksc`, `ceradsc`, `pmi` |
| ROSMAP R1 `rosmap_50batch_specimen_metadata...csv` | `SampleID`, `batch.channel`, `projid`, `IndividualID`, `SpecimenID`, `Batch` | `Emory.Coded.Dx.2017`, `DxNumeric.Emory2017`, `EmoryStrictDx.2019` |
| SomaScan `..._sample_metadata.csv` | `projid_visit`, `projid`, `Visit`, `study` | `msex`, `age_at_visit`, `Diagnosis`, `cogn_global`, `apoe_genotype`, `educ`, `age_death`, `Storage_days` |
| DiverseCohorts `final_traits_*.xlsx` | `ID`, `Dataset`, `Batch`, `Channel`, `individualID`, reagent lot | — none present |
| PDRD Olink `proteomics_*_samples.csv` | `participant_id`, `sample_id`, `visit_month` | — none present |

`Dataset` / `study` are retained — which brain bank contributed a specimen is provenance,
not phenotype.

**Enforcement.** A build gate, not a convention: every emitted column name is checked against
an allow-list and a clinical deny-list (`dx`, `diagnosis`, `cogdx`, `cogn`, `braak`, `cerad`,
`apoe`, `educ`, `age`, `sex`, `msex`, `race`, `ethnicity`, `mmse`, `updrs`, `moca`). A hit
**fails the build**. The result is recorded in the release report.

---

## 1. Inventory: what is on disk

Working dir: repository root — 78 GB, 703 files.

| # | Layer | Assay | Study | Site | Primary file | Shape (feat × samp) | Sample key | `person_id` from |
|---|---|---|---|---|---|---|---|---|
| 1 | brain | TMT LC-MS/MS | DiverseCohorts frontal | DLPFC | `syn53185479/.../n1086_residual_log2_batch.csv.gz` | 9,180 × 1,086 | `<dataset>_b<batch>.<channel>` | ✅ traits xlsx |
| 2 | brain | TMT LC-MS/MS | DiverseCohorts temporal | TCX/STG | `syn53185479/.../n278_residual_log2_batch.TCX.csv.gz` | 9,734 × 278 | `b<batch>.<channel>` | ✅ traits xlsx |
| 3 | brain | TMT LC-MS/MS | ROSMAP round 1 | DLPFC | `syn17015098/.../C2.median_polish_corrected_log2(...).csv.gz` | 8,817 × 400 | `b<batch>.<channel>` | ✅ specimen metadata |
| 4 | brain | TMT LC-MS/MS | ROSMAP round 2 | DLPFC | `syn17015098/syn30390636/.../rosmaptmt_r2_14batch_..._Proteins.txt.gz` | PD wide, 224 ch | `F<plex>.<channel>` | ✅ traits xlsx |
| 5 | brain | TMT LC-MS/MS | MSBB 19-batch | PHG | `syn21347564/.../msbb_19batch_normalized_tmt_matrix.xlsx.gz` | 12,147 × 190 | `TMT<batch>_<group>_<specimen>` | ❌ needs `syn21893059` |
| 6 | brain | LFQ LC-MS/MS | MayoRNAseq | TCX | `syn7431760/.../Mayo_Proteomics_TC_proteinoutput.txt.gz` | 6,585 × 230 (30 pool) | `mayo_b<batch>_<specimen>_<position>` | ❌ needs `syn20827192` |
| 7 | brain | LFQ LC-MS/MS | MSBB / MSSM | PFC | `syn20801227/.../MSSM_Proteomics_PFC_PROTEINOUTPUT.txt.gz` | 6,242 × 306 (22 pool) | `b<batch>_<position>_<specimen>` | ❌ needs `syn21893059` |
| 8 | brain | LC-SRM | ROSMAP | DLPFC | `syn10468856/.../final_srm_data.txt.gz` + `srm_data_1226.txt.gz` | 126 × 1,364 / 68 × 1,340 | `P<plate>_<well>` | ✅ `subject.id` |
| 9 | plasma | SomaScan v4.1 7k | ROSMAP | — | `syn64957327/..._ANML_log10.csv.gz` | 7,289 × 972 | `projid_visit` | ✅ `projid` |
| 10 | csf | Olink Explore 1536 | PPMI + PDBP | — | `PDRD/proteomics-CSF-PPEA-D03/.../olink_explore_format_*.csv.gz` | 1,463 × 1,320 | `<cohort>-<participant>-<visit>-CSF-PPEA-D<nn>` | ✅ `participant_id` |
| 11 | plasma | Olink Explore 1536 | PPMI + PDBP | — | `PDRD/proteomics-PLA-PPEA-D03/.../olink_explore_format_*.csv.gz` | 1,463 × 1,552 | same convention | ✅ `participant_id` |
| 12 | csf | DIA LC-MS | PPMI + PDBP | — | **[updated] landed 2026-09-03** `PDRD/proteomics-CSF-PDIA/CSF-PDIA_protein_batch-corrected.csv` (long form) | 249 × 2,807 | `<cohort>-<participant>-<visit>-CSF-PDIA` | ✅ `participant_id` |
| 13 | plasma | DIA LC-MS | PPMI + PDBP | — | **[updated] landed 2026-09-03** `PDRD/proteomics-PLA-PDIA/PLA-PDIA_protein_batch-corrected.csv` (long form) | 202 × 1,459 | same convention | ✅ `participant_id` |

**[updated] — three of the `person_id` cells above were wrong when written.** Rows 5, 6 and 7 read
"needs `syn21893059`" / "needs `syn20827192`", but both objects were already in the tree under
their plain filenames (`MSBB_biospecimen_metadata.csv`, `MayoRNAseq_biospecimen_metadata.csv`) and
had been all along. Presence was being detected by Synapse-ID *path*, so anything delivered under a
different accession read as absent — the defect tracked as REQUIREMENTS OPS‑2b, now fixed by
content-hash detection. All three strata resolve every non-pool sample. Row 9's `972` is the
delivered column count; 970 ship, after 3 rows with no `projid` at source are removed.

Excluded: `syn61332963` (AMP RA/SLE urine Kiloplex) is prioritization **3** — also contains two
zero-byte files from a truncated download. And the three uncompressed root-level
`rushtmt_fullmulticonsensus_{PSMs,QuanSpectra,MSMSSpectrumInfo}.txt` files (64 GB) are
PSM/spectrum-level intermediates below the protein layer — recommend archiving off this volume.

### 1.1 Priority-2 coverage

All nine priority-2 rows of the tracking file are accounted for.

| Tracking-file row | Inventory # | Buildable |
|---|---|---|
| AMP AD · TMT LC/MS · AMP-AD_DiverseCohorts | 1, 2 | **[updated]** ✅ built — site confirmed as *superior temporal gyrus* for the temporal stratum (N2) |
| AMP AD · LC-SRM · ROSMAP | 8 | ✅ |
| AMP AD · TMT LC/MS · ROSMAP | 3, 4 | **[updated]** ✅ both built — R2 via the §5.1(b) re-derivation |
| AMP AD · SomaScan · ROSMAP | 9 | ✅ (is plasma, not brain) |
| AMP AD · LC-MS · MayoRNAseq | 6 | **[updated]** ✅ built — 200/200 donors resolved |
| AMP AD · LC-MS · MSBB | 7 | **[updated]** ✅ built — 284/284 donors resolved |
| AMP AD · TMT LC/MS · MSBB | 5 | **[updated]** ✅ built — 186/186 donors resolved (4 GIS channels excluded) |
| AMP PDRD · Olink · Unified Cohorts | 10, 11 | ✅ |
| AMP PDRD · LC-MS · Unified Cohorts | 12, 13 | **[updated]** ✅ built — data landed 2026-09-03 |

### 1.2 Metadata availability

Searched the full tree for every Synapse ID the tracking file names in its metadata columns.
**[updated] — seven have since landed** (`syn51757645` biospecimen; `syn21323404`, `syn21893060`,
`syn22344998`, `syn23474101`, `syn23569441`, `syn53185805` assay), and two more were already in the
tree under plain filenames (`syn20827192`, `syn21893059`). `ASSAY.platform` now resolves for six of
the eight strata that lacked it — ROSMAP SRM does not, its metadata records `platform = NA`.
As originally written:

**None are present**, at any depth, as directory or filename:

`syn51757645` `syn20827192` `syn21893059` (biospecimen) · `syn53185805` `syn23569441`
`syn21323404` `syn65414914` `syn65471938` `syn65471939` `syn65473039` `syn23474101`
`syn22344998` `syn21893060` (assay) · `syn73713766` `syn73713767` (clinical)

Only *Data Location* syn IDs were downloaded. Every piece of sample annotation we have arrived
bundled inside a proteomics folder — ROSMAP R1 specimen metadata, the ROSMAP R2 and
DiverseCohorts traits sheets, SomaScan sample + protein metadata, the PDRD Olink sample
inventories, the SRM feature maps. Six of nine datasets are covered by accident of packaging.

**Mayo TCX, MSBB PFC and MSBB PHG have no sample annotation at all** — their run labels cannot
resolve to a donor. Minimum pull: `syn20827192` and `syn21893059`; `syn51757645` additionally
settles the DiverseCohorts site question.

The absent clinical files were **not on the critical path** — clinical is out of scope per §0.
**[updated] — two of them now are, for one column each.** Decision D‑J reclassified post‑mortem
interval as a biospecimen CDE rather than phenotype, so `syn73713766` (Mayo) and `syn73713767`
(MSBB) move onto the critical path **for `pmi` and `pmiUnits` only**, under an allow‑list that
loads nothing else from them. Neither is in the 13‑object delivery, so MSBB and Mayo PMI ship as
`source_not_acquired`. DiverseCohorts is harder: no individual‑metadata object is registered for
it anywhere and `syn51757645` carries no PMI column, so its source is still unidentified. ROSMAP
is unaffected — its PMI is already on disk, already in hours, and joins on `projid`.

### 1.3 Tracking-file corrections

| Tracking file says | Data says |
|---|---|
| ROSMAP SomaScan — brain DLPFC, n=610 | **plasma**, 970 samples / 887 donors *(was 973/888 before the 3 donor‑less rows were dropped)* |
| DiverseCohorts temporal — superior temporal gyrus | **[updated] the tracking file was right.** `syn51757645` records `superior temporal gyrus` for all 280 temporal TMT specimens. An earlier note in this plan and a hardcoded audit verdict both claimed "TCX = temporal cortex (confirmed)" — both wrong. See REQUIREMENTS VER-1b |
| MSBB TMT PHG, n=800 | n=190 (63 control + 127 AD) |
| DiverseCohorts, n=850 | **1,082** frontal + 278 temporal post‑QC *(1,086 frontal ship from the deposit; 4 carry no donor id in either the traits sheet or the assay metadata, which also states they are not controls, so they are dropped as unresolvable rather than mislabelled as pools)* |
| ROSMAP TMT, n=610 | R1 400 channels + R2 224 channels |
| AMP PDRD Olink, n=375 | 398 CSF + 413 plasma participants |
| AMP PDRD untargeted LC-MS, n=800 | **[updated] built.** The batch‑corrected protein‑level DIA landed 2026‑09‑03: 2,807 CSF and 1,459 plasma samples, 937 participants |
| metadata syn IDs populated | **[updated] 13 of 15 present**, detected by content hash rather than by path. The 2 absent are the Mayo and MSBB individual‑metadata objects, now needed for `pmi` alone |

### 1.4 Is the proteomics data itself adequate, per target study?

Separate question from metadata: does each study have a usable protein-level quantitative
matrix on disk?

| Study | Protein-level matrix? | Adequate as-is | Gap |
|---|---|---|---|
| DiverseCohorts frontal | ✅ 9,180 × 1,086 batch-regressed, + 11,748 × 1,191 pre-correction | ✅ | site label unresolved |
| DiverseCohorts temporal | ✅ 9,734 × 278 + 11,003 × 304 | ✅ | site label unresolved |
| ROSMAP TMT R1 | ✅ TAMPOR-corrected 8,817 × 400 (+ A1/C1 variants) | ✅ | — |
| ROSMAP TMT R2 | ⚠️ raw Proteome Discoverer abundance ratios only | **no** | no delivered batch-corrected matrix; we would have to derive one, which breaks comparability with R1 unless recipe §5.1(b) is used |
| MSBB TMT PHG | ✅ 12,147 × 190 linear reporter intensities | ✅ | linear, needs log2; `person_id` blocked |
| Mayo TCX LFQ | ✅ 6,585 × 230 MaxQuant LFQ | ✅ | 30 pool channels; `person_id` blocked |
| MSBB PFC LFQ | ✅ 6,242 × 306 MaxQuant LFQ | ✅ | 22 pool channels; `person_id` blocked |
| ROSMAP SRM | ✅ two panels, 126 × 1,364 and 68 × 1,340 | ✅ | two panels to reconcile (D5) |
| ROSMAP SomaScan | ✅ 7,289 × 972, ANML + non-ANML | ✅ | — |
| PDRD Olink CSF | ✅ 1,463 × 1,320 across 4 panels, 3 formats | ✅ | — |
| PDRD Olink plasma | ✅ 1,463 × 1,552 | ✅ | plain vs `_retracted` unresolved (D2) |
| PDRD DIA CSF | ❌ SDRF only | **no** | matrix not downloaded |
| PDRD DIA plasma | ❌ SDRF only | **no** | matrix not downloaded |

Ten of thirteen have an adequate protein-level matrix. ROSMAP R2 has data but not in a usable
processed form; the two PDRD DIA rows have none.

### 1.5 Deliverable: derived inventory

The tracking file has drifted from the data, so the build emits a machine-generated inventory —
regenerated every run from the files themselves, never transcribed:

- `inventory/derived_inventory.tsv` — one row per (layer × assay × study): computed feature and
  sample counts, participant counts, timepoint sets, native scale, value range, sample-key
  pattern, `person_id` resolvable (bool), metadata-present flags, build status
- `inventory/tracking_file_diff.tsv` — every field where data and tracking file disagree, so
  §1.3 produces itself and stays current
- `inventory/missing_objects.tsv` — §1.2 as data: each referenced object, present or not, and
  which capability its absence blocks

A metadata pass, not a data pass — cheap enough to re-run after every download.

---

## 2. Layers

**One HDF5 file**, with a keyed layer per (assay × matrix) combination. **[updated] — the
filename is datestamped** (`sysbio_proteomics-<MMDDYYYY>.h5`, decision D‑I) rather than the fixed
`amp_proteomics.h5` this section originally named, so no build can overwrite an audited one.
Consumers resolve the newest or pin deliberately.

**[updated] — all eight layers are built.** Row and protein counts are from
`sysbio_proteomics-09042026.h5`.

| Layer key | Matrix | Assay | Studies | Rows | Proteins | Status |
|---|---|---|---|---:|---:|---|
| `tmt_brain_proteomics` | brain | TMT LC-MS/MS | DiverseCohorts ×2, ROSMAP R1+R2, MSBB PHG | 2,160 | 14,916 | ✅ 5 strata |
| `lfq_brain_proteomics` | brain | LFQ LC-MS/MS | MayoRNAseq TCX, MSBB PFC | 536 | 4,437 | ✅ |
| `srm_brain_proteomics` | brain | LC-SRM | ROSMAP | 2,704 | 185 | ✅ |
| `soma_plasma_proteomics` | plasma | SomaScan v4.1 7k | ROSMAP | 970 | 7,289 | ✅ |
| `olink_plasma_proteomics` | plasma | Olink Explore 1536 | PPMI + PDBP | 1,380 | 1,463 | ✅ 172 withdrawn (D‑D) |
| `olink_csf_proteomics` | csf | Olink Explore 1536 | PPMI + PDBP | 1,320 | 1,463 | ✅ |
| `dia_plasma_proteomics` | plasma | DIA LC-MS | PPMI + PDBP | 1,459 | 202 | ✅ one stratum, not two |
| `dia_csf_proteomics` | csf | DIA LC-MS | PPMI + PDBP | 2,807 | 249 | ✅ |

The `person_id` gaps this table originally flagged for MSBB, Mayo and the LFQ layer are closed —
every non-pool sample in every stratum resolves to a donor.

**Formerly deferred, now built** — kept as a record of what the deferral mechanism was for:

1. `dia_csf_proteomics` / `dia_plasma_proteomics` — the batch-corrected protein-level matrices
   landed 2026‑09‑03. Plasma is **one** stratum, not the two planned: the delivered product carries
   the plain `-PLA-` token on every sample and the native/depleted split survives only in the SDRF
   run manifests, which are out of scope. Splitting on a distinction the data does not carry would
   be an assertion, not a harmonisation.
2. **ROSMAP TMT round 2** — the §5.1(b) re-derivation was applied and `rosmap_tmt_r2` resolves
   210/210 donors, so `tmt_brain_proteomics` ships with **5 strata**, not 4.

One note on the naming example: the SomaScan data here is **plasma**, not CSF — the ROSMAP
SomaScan v4.1 deposit is `..._ROSMAP_plasma_Soma7k_...`. The Oh et al. paper the README cites is
*about* a CSF biomarker, but the deposited matrix is plasma. So the layer is
`soma_plasma_proteomics`; there is no `soma_csf_proteomics` in the priority-2 set.

Each layer holds one or more **study strata** — a stratum is one source pipeline, and values are
always native within a stratum (§5).

---

## 3. The uniform column contract

Every table in every layer carries the same columns, in the same order, with the same names and
types. Columns that do not apply to an assay are present and **null** — never absent, never
renamed. This is what makes brain, plasma and CSF parseable by one reader.

### 3.1 Link keys

| Column | Type | Note |
|---|---|---|
| `person_id` | int64 | **the central link.** Never null in a released row |
| `person_source_value` | str | the ID the source AMP uses — `projid`, `individualID`, `participant_id` |
| `specimen_id` | int64 | resolves ASSAY_TO_SPECIMEN |
| `specimen_source_value` | str | original column header — full traceability to the source file |
| `visit_occurrence_id` | int64 | nullable |
| `visit_name` | str | `baseline` \| `follow-up` \| `postmortem` |
| `visit_month` | int32 | months since baseline; null for postmortem. **[updated]** this is the **canonical cross-layer time axis** — present in every layer, populated wherever the dataset is longitudinal, and the only column a consumer should join time on. Where it is derived rather than stated, `visit_month_evidence` says so. Decision D-H, REQUIREMENTS FR-14 |
| `visit_index` | int32 | source-native visit number; null where absent. **[updated]** retained for traceability, never the join key |
| `assay_id` | int64 | resolves ASSAY |
| `file_id` | int64 | resolves FILES — which source file this value came from. **[updated]** the relation is **one-to-many** wherever a stratum's matrix is several files (Olink ships one per panel), so it is carried as a **`sample × file` bridge table**, never a delimited list — decision D-F, REQUIREMENTS FR-13 |

### 3.2 SysBio ASSAY CDEs (Figure 2, green)

`assay_id` · `assay_source_value` · `assay_type` · `platform` · `suspension_type` ·
`analyte_type` · `analysis_pipeline`

### 3.3 SysBio FILES CDEs (Figure 2, green)

`file_id` · `file_name` · `current_version` · `assay_id` · `file_role` · `study` · `grant` ·
`array_type` · `analysis_type` · `biosample_type` · `tissue` · `cell_type` · `species` ·
`processing_status` · `file_format` · `file_size_bytes` · `created_on` · `modified_on` ·
`drs_id`

Note the FILES CDEs already carry `tissue`, `biosample_type`, `cell_type` and `species` — the
biospecimen descriptors come from the extension, so no OMOP `SPECIMEN` build-out is needed here.
`grant` carries the AMP program (`AMP AD` / `AMP PDRD`), which is what the return path in §7.3
filters on.

### 3.4 Harmonized proteomics columns

Uniformly named across every assay — this is the layer that makes an Olink NPX value and a TMT
log-ratio addressable by the same code:

| Column | Type | Note |
|---|---|---|
| `feature_id` | str | **UniProt primary accession** — the cross-assay join key |
| `feature_source_value` | str | original row identifier |
| `gene_symbol` | str | |
| `protein_group` | str | full accession set, `;`-delimited |
| `protein_group_size` | int16 | |
| `measurement_concept_id` | int64 | OMOP concept where resolvable, else 0 |
| `value` | float32 | the harmonized measurement |
| `value_scale` | str | exact scale token — see §3.6 |
| `value_unit` | str | |
| `is_missing` | bool | |
| `lod` | float32 | null where the assay has no LOD concept |
| `below_lod` | bool | flagged, never censored |
| `panel` | str | Olink panel; null elsewhere |
| `panel_lot_nr` | str | Olink; null elsewhere |
| `dilution` | str | SomaScan dilution group; null elsewhere |
| `seq_id` | str | SomaScan; null elsewhere |
| `peptide_sequence` | str | SRM; null elsewhere |

### 3.5 Technical

`batch` · `plate_id` · `channel` · `is_pool` · `qc_status`

`is_pool` marks GIS pools, bridge samples and plate controls — flagged, not deleted.

### 3.6 How each assay fills the contract

Proof that the contract is actually uniform:

| Column | brain/TMT | brain/LFQ | brain/SRM | plasma/SomaScan | plasma+csf/Olink |
|---|---|---|---|---|---|
| `assay_type` | TMT LC-MS/MS | LFQ LC-MS/MS | LC-SRM | SomaScan v4.1 7k | Olink Explore 1536 |
| `platform` | Orbitrap | Orbitrap | triple quadrupole | SOMAmer array | PEA + NGS readout |
| `analyte_type` | protein | protein | protein | protein | protein |
| `suspension_type` | bulk tissue | bulk tissue | bulk tissue | not applicable | not applicable |
| `analysis_pipeline` | Proteome Discoverer + TAMPOR | MaxQuant LFQ | Skyline | SomaLogic ANML | Olink NPX + bridging |
| `tissue` | brain | brain | brain | plasma | plasma / CSF |
| `biosample_type` | postmortem tissue | postmortem tissue | postmortem tissue | fluid | fluid |
| `species` | Homo sapiens | Homo sapiens | Homo sapiens | Homo sapiens | Homo sapiens |
| `value_scale` | `log2_rel_batch_median` | `log2_lfq_intensity` | `log2_ratio` | `log10_rfu_anml` | `npx_log2` |
| `visit_name` | postmortem | postmortem | postmortem | follow-up | baseline / follow-up |
| `lod` | null | null | null | null | populated |
| `panel` | null | null | null | null | populated |
| `dilution` | null | null | null | populated | null |

### 3.7 CDE coverage — what we can actually populate

Graded by evidence, because the difference matters:

- **A — derivable from the data files.** Machine-readable, no assertion.
- **B — documented in a bundled README or script.** True, but transcribed by hand and cited.
- **C — not available.** Needs the assay-metadata Synapse objects (§1.2), a publication, or is
  genuinely not applicable.

**ASSAY CDEs**

| CDE | tmt_brain | lfq_brain | srm_brain | soma_plasma | olink_plasma / csf |
|---|---|---|---|---|---|
| `assay_id` | A mint | A | A | A | A |
| `assay_source_value` | A | A | A | A | A |
| `assay_type` | A | A | A | B | **A** — `assay_type` is a column in the NPX files |
| `platform` | **mixed** — see below | C | partial | B | B |
| `suspension_type` | A | A | A | A | A |
| `analyte_type` | A | A | A | A | **A** — `omics_type` column |
| `analysis_pipeline` | **mixed** — see below | A (MaxQuant; version C) | C | B (ANML) | B (NPX + bridging) |

`platform` and `analysis_pipeline` for brain TMT vary by stratum, and the evidence is uneven:

| Stratum | `platform` | `analysis_pipeline` | Evidence |
|---|---|---|---|
| ROSMAP R2 | **Q Exactive HF-X Orbitrap** (A) | Proteome Discoverer (A) | `..._InputFiles.txt` has `Instrument Name` populated, tune `2.9.4.2954` |
| ROSMAP R1 | Orbitrap FTMS only (partial) | **Proteome Discoverer 2.3.0.522** (A) | `InputFiles` gives the PD version but `Instrument Name` is blank — searched from `.msf`, so the model is lost. `SpecializedTraces` gives scan filter `FTMS + p NSI Full ms [350–1500]`, which fixes the class but not the model |
| DiverseCohorts | C | **FragPipe, TMT 18-plex** (A) | `loader.R` reads `FRAGPIPE_OUTPUT_FOLDER`, `tmt_channels <- 18`, `protein.tsv` |
| MSBB PHG | C | C | xlsx header says only "TMT-LC/LC-MS/MS", 20 batches, BM-36 |
| Mayo / MSSM | C | **MaxQuant** (A) | `proteinGroups.txt` signature: `Mol. weight [kDa]`, `Q-value`, `Only identified by site`, `Reverse`, `Potential contaminant`, `LFQ intensity`, `iBAQ`. Version not recorded |
| ROSMAP SRM | instrument nickname `Smeagol` only (partial) | C | 916 replicate names carry it; a lab alias, not a model |

So `platform` is the weakest CDE in the set — populatable for exactly one of eight layers from
data alone. Pulling the assay-metadata objects in §1.2 is what fixes it.

**FILES CDEs**

| CDE | Coverage | Note |
|---|---|---|
| `file_id`, `assay_id` | A | minted |
| `file_name`, `file_format`, `file_size_bytes` | A | filesystem |
| `file_role`, `analysis_type`, `processing_status` | A | we assign from a controlled list |
| `study`, `grant` | A | from the inventory |
| `tissue`, `biosample_type`, `species` | A | one caveat: DiverseCohorts site is unresolved (§9 blocker 3) |
| `cell_type` | A | `not applicable` — all bulk |
| `array_type` | mixed | Olink: A (panel). SomaScan: B (SOMAmer v4.1). MS assays: not applicable |
| `current_version` | **C** | Synapse object version was not downloaded |
| `created_on`, `modified_on` | **C** | ⚠️ every mtime on disk is `2026-08-21`, i.e. **download** time, not file creation. Using mtime here would be quietly wrong. ROSMAP is the one exception — its PD `InputFiles` carries genuine acquisition dates (`8/19/2019`, `08/03/2021`) |
| `drs_id` | **C** | needs Synapse/GCS DRS resolution. PDRD SDRF carries `gs://` URIs, which is the closest thing present |

**Summary.** Of 7 ASSAY CDEs, 5 are fully populatable across all layers, 1 (`analysis_pipeline`)
for 5 of 6, and 1 (`platform`) for 1 of 6. Of 19 FILES CDEs, 15 are populatable and 4
(`current_version`, `created_on`, `modified_on`, `drs_id`) are not, all four for the same
reason: we have the data files but not the repository objects that describe them.

**Proposed CDE addition.** The five assays differ in detection modality — mass spectrometry,
aptamer, antibody — and that difference matters for interpreting comparability, but Figure 2's
ASSAY has no field for it. Recommend adding `detection_modality` ∈ {mass_spectrometry,
aptamer_affinity, antibody_affinity} to the ASSAY extension. Section 6 of the white paper
explicitly anticipates assay-class-specific fields being added to ASSAY, so this is the
sanctioned mechanism rather than a local invention. **Flagged for the Task Force, not assumed.**

---

## 4. Feature harmonization

`feature_id` is a UniProt primary accession in every layer. Source identifiers vary:

| Source | Form | Resolution |
|---|---|---|
| DiverseCohorts, ROSMAP TMT | `GENE\|UNIPROT` | split; gene may be empty |
| MSBB TMT | `sp\|P04217\|A1BG_HUMAN` | parse middle field |
| Mayo, MSSM (MaxQuant) | `;`-delimited protein groups | leading/majority accession as `feature_id`, full set retained in `protein_group` |
| SomaScan | `SeqId` (e.g. `10000-28`) | via `protein_metadata.csv`; many SeqIds → one UniProt, so `seq_id` stays a secondary key and is **not** silently collapsed |
| Olink | UniProt directly | already primary, already de-duplicated across panels by AMP-PD |
| SRM | `GENE_n` peptide features | needs an explicit gene→UniProt table; phospho-specific features (`tau_PHF1_s404`, `tau_AT8_s202`) stay distinct and are **not** merged into the parent protein |

### 4.1 Olink panels

Panels stay distinct — separate assay runs with their own plates, lots, LODs and QC verdicts.

Measured: cardiometabolic 369, inflammation 365, neurology 364, oncology 365. **Zero UniProt
overlap** — 369+365+364+365 = 1,463 = the exact union. So `panel` is a clean 1:1 attribute of
each feature and no de-duplication is needed on our side.

**CSF and plasma share an identical 1,463-protein feature space** (0 CSF-only, 0 plasma-only),
so the two layers join feature-for-feature with no reconciliation.

QC is per (sample × panel), not per sample — carried on the measurement row via `panel` +
`qc_status` rather than flattened to one per-sample flag. `lod` is per (assay × plate), so it
rides at full measurement resolution.

Any per-feature statistic is computed **within panel** (and within SomaScan dilution group) —
these differ in dynamic range and LOD behaviour.

---

## 5. Scale policy

**Retain fidelity — do not z-score where a stratum is already on a common scale.**

A stratum is one (assay × study). Within a stratum, values are native. Only monotone,
invertible transforms are applied to reach `value` (log2 where a source is linear), and the
chain is recorded in `value_scale` plus a `transform_chain` attribute.

| Layer / assay | `value_scale` | z? |
|---|---|---|
| brain / TMT | per stratum | ✅ only place it is needed — four incompatible recipes |
| brain / LFQ | `log2_lfq_intensity` | ⚠️ D3 |
| brain / SRM | `log2_ratio` | ❌ |
| plasma / SomaScan | `log10_rfu_anml` | ❌ |
| plasma / Olink | `npx_log2` | ❌ |
| csf / Olink | `npx_log2` | ❌ |

Where z is written it is an **additional** column, never a replacement, and `z_center` /
`z_scale` are stored per feature so it is fully reversible. Nothing is imputed — missingness is
explicit in `is_missing`.

### 5.1 Brain TMT — the one real scale problem

Five sources, four recipes: batch-regressed log2 residual (DiverseCohorts), TAMPOR median-polish
log2 ratio (ROSMAP R1), raw Proteome Discoverer abundance ratios (ROSMAP R2), linear TMT
reporter intensities (MSBB). Two ways to satisfy one-scale-per-assay-group, and I recommend
shipping both:

**(a) Per-study strata at native scale — the returnable.** Five strata sharing the §3 contract,
nothing transformed to force agreement. This is also the honest thing to hand an AMP: its own
pipeline's output, unaltered.

**(b) One pooled brain-TMT layer via a defined re-derivation — labelled derived.** All five
reduce to per-batch protein abundance, so one recipe reaches a common scale without any
distribution rescaling:

1. per-channel loading normalization (divide by column total, × median column total)
2. log2
3. per-batch, per-protein centering on the within-batch median

That yields *log2 abundance relative to batch median* — exactly how ROSMAP R1's `C2` is already
constructed, so C2 passes through untouched and the other four are brought to it. A defined
normalization, not a distribution rescale, so it preserves fidelity in a way z does not.

This interacts with **D4**: recipe (b) must start from the *pre-correction*
`normAbundances_post-FP` matrices (11,748 proteins), not `n1086_residual_log2_batch` — regressing
out batch and centering on batch median are different operations and cannot be mixed.

---

## 6. Timepoints

All samples and all timepoints. Nothing collapsed to baseline.

| Layer / assay | Axis |
|---|---|
| plasma / Olink | `visit_month` ∈ {0,3,6,9,12,18,24,30,36,42,48,54,60,72,84,96} — 16 values |
| csf / Olink | `visit_month` ∈ {0,3,6,12,18,24,30,36,48,54,60,84,96} — 13 values |
| plasma / SomaScan | `visit_index` 0–25; `visit_month` null — see D10 |
| brain / all | `visit_name = postmortem`, `visit_month` null |

PDRD longitudinal structure, measured: 398 CSF and 413 plasma participants, and **every CSF
participant also has plasma** (0 CSF-only, 15 plasma-only) — so paired CSF↔plasma analysis is
available for all 398 via `person_id`. Most participants carry 3 timepoints; the distribution
runs 1 → 7. Cohort prefixes `PD-` (PDBP) and `PP-` (PPMI) are parsed out rather than left
embedded in the ID.

On ROSMAP `visit_month`: the only route to a months axis is differencing `age_at_visit` within
donor, and that column is dropped under §0. An elapsed interval is arguably a VISIT_OCCURRENCE
property rather than a phenotype, but it derives from a dropped field — so the conservative
default is `visit_index` only. **D10.**

---

## 7. Outputs

### 7.1 The HDF5

One file, one key per layer:

```
amp_proteomics.h5
├── tmt_brain_proteomics/
│   ├── assay/          ASSAY CDEs (§3.2), one row per stratum
│   ├── files/          FILES CDEs (§3.3), one row per source file
│   ├── sample/         link keys (§3.1) + technical (§3.5), len = n_samples
│   ├── feature/        feature identity (§3.4), len = n_features
│   ├── measurement/    value, is_missing, lod, below_lod — (n_features, n_samples)
│   │                   + value_z, z_center, z_scale where §5 applies
│   └── attrs           value_scale, transform_chain, layer_key, n_strata
├── lfq_brain_proteomics/     … same structure
├── srm_brain_proteomics/     … same structure
├── soma_plasma_proteomics/   … same structure
├── olink_plasma_proteomics/  … same structure
├── olink_csf_proteomics/     … same structure
└── provenance/        source file, syn ID, sha256, transform_log — file-wide
```

Every layer has an identical internal structure, so one reader opens any key. Features ×
samples, chunked, gzip-4, shuffle. Root attributes carry `cdm_version`, `code_sha`,
`build_utc`, and the layer key list.

Features × samples, chunked, gzip-4, shuffle. Attributes carry `value_scale`, `transform_chain`,
`cdm_version`, `code_sha`, `build_utc`.

### 7.2 Tabs (the returnable)

**[updated] — built, and wide only.** This section described the bundles from the first draft and
they had never been written; `build/` held the HDF5, the report and the inventory files and
nothing else, so the section read as delivered when it was not. It is now produced on every build.

One bundle per layer key, three tables, always the same three:

```
tabs/<layer_key>/
├── <layer_key>_matrix.parquet     feature_id + one column per specimen_id
├── <layer_key>_samples.parquet    §3.1 + §3.2 + §3.3 + §3.5, one row per sample
├── <layer_key>_features.parquet   §3.4 identity half, one row per feature
└── <layer_key>_bundle.json        which value surface the matrix carries, and how to join
```

Mirrored as `.tsv.gz`. Parquet is primary — it preserves dtypes and the null distinction that
makes three-state `below_lod` meaningful.

**The long form is dropped** (decision D‑G). It was the natural place to realise §3's uniform
contract literally, but the sizes are the problem this section already anticipated — SomaScan
alone is 7.1 M triples and pooled brain TMT ~25 M — and the contract is realised anyway, because
the runtime column union means every layer's `samples` and `features` carry the same columns in
the wide form too. What the long form would have added is one row shape, at ~35 M rows.

Each bundle records its own `value_surface`, because which surface ships is layer-dependent:
`tmt_brain_proteomics` carries `values_z`, every other layer `values`. Leaving that implicit
would make the matrices look comparable when they are not.

The matrix is written in the same pass that writes the HDF5's `values`, from the same in-memory
array, and the tables in the same pass that writes the HDF5's tables. §8 S8 asks for exactly that:
a standalone HDF5→parquet converter would be simpler and would give up the one property worth
having — **the two cannot drift**.

**SEC‑5 applies to these exactly as to the HDF5.** They carry the same donor-level identifiers,
so they are controlled-access artifacts and not a publishable by-product; `.gitignore` denies
them by name as well as by directory.

### 7.3 Return path

**[updated] — implemented, and demonstrated on every build** rather than described. Output is
`build/returns/<grant>/`.

Every row carries `grant` (the AMP program), `study` and `person_source_value`, so returning
data to an AMP is a filter, not a reconstruction: slice on `grant`, take the `specimen_id` set,
project matrix columns, drop features that go entirely missing. `person_source_value` and
`specimen_source_value` let the receiving AMP re-key to its own identifiers without us shipping
a crosswalk.

The sidecar manifest **demonstrates** the boundary instead of asserting it: it names the grants
present in the build, the grant being returned, and the count of rows from any other grant that
reached the slice — which must be zero, and is checked as an audit row. The current build returns
6,370 AMP AD specimens and 6,966 AMP PDRD specimens, with zero cross-grant rows in either.

It runs on every build rather than at release time because returning one AMP programme's donors
to another is not recoverable after the fact, and a check that only runs when someone remembers
to run it is not a control.

---

## 8. Pipeline

**S0 — Environment.** **[updated]** `h5py` 3.16.0 is installed and is what writes the artifact;
`pyarrow` writes the §7.2 bundles. Pins are in [requirements.txt](requirements.txt). `pytables` is
still not installed and is not needed. As originally written: `h5py` and `pytables` are **not
installed** (pandas 2.3.3, numpy 2.4.4,
scipy 1.17.1, pyarrow 23.0.1, openpyxl 3.1.5 are). 31 GB RAM / 16 cores / 1.5 TB free is ample —
largest single read is 12,147 × 190.

**S1 — Manifest + derived inventory.** One metadata pass producing the `files/` groups, the
audit trail, and the three §1.5 inventory outputs.

**S2 — Acquire.** `syn20827192`, `syn21893059`, `syn51757645`; PDRD DIA matrices from GCS.

**S3 — Readers.** One per inventory row, returning
`(feature_source_value[], specimen_source_value[], values[,], native_scale, transform_chain)`.
Awkward cases: MSBB xlsx publication header (data starts row 5, sample IDs row 4); MaxQuant
`LFQ intensity <sample>` column selection; Proteome Discoverer
`Abundance Ratio: (Fn, 127N) / (Fn, 126)` parsing; Olink long→wide across 4 panels.

For Olink, read the **`olink_explore_format_*`** files, not `matrix_*` — only they carry `LOD`,
`PlateID`, `Panel_Lot_Nr`, `QC_Warning`, `Distribution_QC`, `Outliers_QC`. Tag every feature with
its panel *before* concatenation.

**S4 — Feature harmonization** (§4).

**S5 — Sample harmonization.** `specimen_source_value` → `specimen_id` / `person_id` /
`visit_occurrence_id`. Flag pools rather than dropping them. DiverseCohorts frontal has 107
duplicated `individualID` across 729 mapped channels — replicate specimens per donor, so they
stay distinct `specimen_id` under one `person_id`.

**S6 — Scale** (§5). **S7 — Z layer** where §5 applies.

**S8 — Write + validate.** HDF5 and tabs emitted from the **same in-memory object in one pass**,
so they cannot drift. Validation: dimensions, feature overlap across strata, `person_id`
coverage, missingness per stratum, scale summaries, round-trip to float32 tolerance, and the §0
scope gate — which fails the build on a hit rather than warning.

---

## 9. Blockers and decisions

**[updated] — all resolved or deliberately closed.** Outcomes below; the work list is
[AUDIT_RESPONSE_PLAN.md](AUDIT_RESPONSE_PLAN.md) §4, the tracked non-conformances are
[REQUIREMENTS.md](REQUIREMENTS.md) §8. Original text retained beneath each item where it still
explains the reasoning.

### Blockers

| # | Blocker | Outcome |
|---|---|---|
| 1 | No portal metadata; Mayo TCX, MSBB PFC, MSBB PHG have no route to `person_id` | ✅ **resolved.** All three resolve every non-pool sample — 200/200, 284/284, 190/190. The two biospecimen files were in the tree under plain filenames all along |
| 2 | PDRD untargeted DIA — SDRF manifests only | ✅ **resolved.** The batch-corrected protein-level data landed 2026‑09‑03 and both layers are built: 2,807 CSF + 1,459 plasma samples, 937 participants. Two intake findings changed the design — plasma is one stratum rather than two, and the wide `*_matrix.csv` siblings encode non-detection as literal `0`, so the **long form** is ingested and the matrix treated as a derived view. The SDRFs stay out of scope, so `platform` is declared at evidence `B` rather than read per sample |
| 3 | DiverseCohorts anatomic site unresolved | ✅ **resolved by `syn51757645` — and it reverses the earlier reading.** The site is *superior temporal gyrus*, as the tracking file said |
| 4 | ROSMAP R2 has no batch-corrected matrix locally | ✅ **resolved.** Built via the §5.1(b) re-derivation; 210/210 donors |

### Decisions

| # | Subject | Outcome |
|---|---|---|
| **D2** | Olink plasma variant, plain vs `_retracted` | **Neither — drop the 172 disputed samples** and build from plain. The two variants differ only on `BIOREP_pl1`/`pl2`, and the id set there is identical in both, so the removal does not depend on which is corrected. 11.1% of plasma samples; 39 of 413 participants lose every plasma sample. CSF unaffected |
| **D3** | brain/LFQ z layer | **Native only.** `lfq_brain_proteomics` ships `values`, no `values_z` |
| **D4** | DiverseCohorts input | **Batch-regressed `n1086_residual_log2_batch`** as primary. The pre-correction matrix is not read, so §5.1(b) is unavailable for this stratum |
| **D5** | SRM scope | **Two strata in one layer**; no forced feature union |
| **D6** | 64 GB PSM-level files | **Closed, no action.** Nothing reads them; only acquisition descriptors were sampled. They stay on the volume |
| **D7** | Pooled brain TMT, (a) or (a)+(b) | **Both** — R2 came in via (b) |
| **D8** | Pools and QC-warn samples | **Shipped flagged**, consumer filters. Standing policy: flag, never silently delete. Removals are identity failures only — D2's 172 disputed Olink plasma samples; 3 SomaScan rows with no `projid` at source; and 4 DiverseCohorts DLPFC rows with no donor id in either the traits sheet or the assay metadata, which the build had previously mislabelled as reference pools. The **one** quality-based removal is the 12 DIA plasma runs below the completeness floor, and it exists only because that deposit ships no sample-level QC of its own |
| **D9** | Long-form tabs | **Closed — not wanted.** Wide analysis-ready form only (§7.2) |
| **D10** | ROSMAP time axis | **Derive `visit_month` as `Visit × 12`**, labelled an approximation, with `visit_index` retained. ROSMAP `Visit` is the follow-up year, so this is a structural conversion of an ordinal index and needs no clinical field — SEC-1 untouched. `visit_month` becomes the canonical cross-layer axis (§3.1) |
| **D11** | Post-mortem interval | **PMI is a CDE.** Carried as a biospecimen technical attribute, not phenotype; removed from the clinical deny-list. Deliberate narrowing of §0 |
| **D12** | `detection_modality` CDE | **Closed for now, noted.** The proposal stands — five assays differ in detection modality and ASSAY has no field for it — but it is a Task Force decision. Modality stays implicit in `platform`; not implemented, not assumed |

### Original text


### Blockers

1. **No portal metadata downloaded** (§1.2). Mayo TCX, MSBB PFC, MSBB PHG have no route to
   `person_id` — and `person_id` is the central link, so those three cannot be released at all
   until `syn20827192` and `syn21893059` land. Both are small.
2. **PDRD untargeted DIA** — two priority-2 rows, SDRF manifests only. Data is at
   `gs://amp-pd-proteomics/datasets/proteomics-{CSF,PLA}-PDIA`, Terra-gated.
3. **DiverseCohorts anatomic site** — tracking file says DLPFC + superior temporal gyrus, but the
   file is named `.TCX.` and its traits sheet carries `ProjID` values like `138_TCX` with
   `Dataset ∈ {Mayo, Emory}`. Needs `syn51757645`.
4. **ROSMAP R2 has no batch-corrected matrix locally** — only raw PD abundance ratios.

### Decisions

- **D2 — Olink plasma variant.** `_retracted` vs plain. Measured: identical 1,552 samples,
  identical row counts, both match the sample inventory exactly; mtimes differ by 0.07 s
  (download order, not a release date). They differ in **178 samples (11.5%) being re-keyed** —
  ≈ two plates of 89, matching the README's `BIOREP_pl1`/`pl2` inversion. So `_retracted` is the
  plate-swap-corrected release. But it also clears **every** QC warning in all four panels
  (933/636/861/1,151 → 0), which the README does not describe. Recommend: identifiers from
  `_retracted`, QC flags from plain, confirm with AMP-PD before final. CSF is unaffected.
- **D3 — brain/LFQ z layer.** Mayo and MSSM are both MaxQuant LFQ but separate experiments, so
  absolute intensities are not comparable even though the scale definition matches. Write `value_z`
  or leave native only?
- **D4 — DiverseCohorts input.** Batch-regressed `n1086_residual_log2_batch` (9,180 proteins) as
  primary, or pre-correction `normAbundances_post-FP` (11,748)? Recipe (b) in §5.1 requires the
  latter. Recommend pre-correction as the re-derivation input, residual retained as a secondary
  layer.
- **D5 — SRM scope.** Two files, different panels (126 vs 68 features, ~1,360 samples each,
  partially overlapping subjects). One stratum with a feature union, or two strata?
- **D6 — 64 GB PSM-level files.** Archive off this volume?
- **D7 — Pooled brain TMT.** Ship §5.1 (a) only, or (a) + (b)? Plan assumes both.
- **D8 — Pools and QC-warn samples.** Included with `is_pool` / `qc_status` set, so the consumer
  filters. Confirm.
- **D9 — Long-form tabs.** Every assay, or on request? Pooled brain TMT long form is ~25 M rows.
- **D10 — ROSMAP time axis.** `visit_index` only (default), or derive elapsed months from the
  dropped `age_at_visit`?
- **D11 — Post-mortem interval.** Present in ROSMAP R2 traits; a standard technical covariate for
  brain proteomics, arguably a biospecimen attribute rather than phenotype. Currently dropped
  under §0. Carry it as a specimen attribute instead?
- **D12 — `detection_modality` CDE.** Propose the ASSAY extension field of §3.6 to the Task
  Force, or leave modality implicit in `platform`?
