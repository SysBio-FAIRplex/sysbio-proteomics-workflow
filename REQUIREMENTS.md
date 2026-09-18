# Requirements — AMP Proteomics harmonisation build (`ampprot`)

Companion to [HARMONIZATION_PLAN.md](HARMONIZATION_PLAN.md). The plan says *what will be
built and why*; this document states *what the system must satisfy* to be considered
correct, safe to run, and safe to publish, and records how the current build measures
against each requirement.

Requirements are labelled `ENV`, `IN`, `SEC`, `FR`, `OUT`, `VER`, `OPS`. Conformance as of
the 2026‑08‑22 build is in §8.

---

## 1. Environment (`ENV`)

| id | requirement |
|---|---|
| ENV‑1 | CPython **≥ 3.12**. The code uses PEP 604 unions and `from __future__ import annotations` throughout. Built and verified on 3.12.3. |
| ENV‑2 | Dependencies exactly as pinned in [requirements.txt](requirements.txt): numpy 2.4.4, pandas 2.3.3, h5py 3.16.0, scipy 1.17.1, openpyxl 3.1.5. |
| ENV‑3 | **Linux or WSL2.** [build.py:31-40](src/ampprot/build.py#L31-L40) reads `MemAvailable` from `/proc/meminfo` to size its memory budget. On a platform without `/proc`, it silently falls back to a 4 GB assumption, which will over‑subscribe RAM on a large layer. |
| ENV‑4 | **≥ 16 GB RAM available at run time.** The scheduler targets 45 % of measured available memory ([build.py:271](src/ampprot/build.py#L271)) and bin‑packs strata into waves against per‑reader expansion factors ([build.py:24-27](src/ampprot/build.py#L24-L27)). The MSBB TMT xlsx (×90) and the Olink long‑form pivots (×45) dominate. |
| ENV‑5 | **≥ 85 GB free disk** for the source tree, plus ~250 MB for the output. The tree measured 83.6 GB across 707 files. |
| ENV‑6 | Multi‑core CPU. Concurrency is `min(16, os.cpu_count())` ([build.py:19](src/ampprot/build.py#L19)), using processes rather than threads because openpyxl and much of the pandas work hold the GIL. |
| ENV‑7 | Run as a module from the repository root: `PYTHONPATH=src python -m ampprot.build`. `ROOT` is derived as two parents up from `config.py` ([config.py:13](src/ampprot/config.py#L13)), so the source tree must sit at the repository root. |

## 2. Input data (`IN`)

| id | requirement |
|---|---|
| IN‑1 | All 12 stratum matrix files declared in `STRATA` ([config.py:172-401](src/ampprot/config.py#L172-L401)) must be present at their declared paths relative to `ROOT`. A missing matrix raises during read; there is no skip‑and‑continue path. |
| IN‑2 | Each stratum's annotation file must be present where `annot_path` is set. Person resolution depends on it; absence yields unresolved donors, which VER‑4 fails on. |
| IN‑3 | The UniProt↔symbol reference ([config.py:435](src/ampprot/config.py#L435)) must be present. It is the only route from SRM gene‑symbol features to accessions. |
| IN‑7 | **`syn73713766` (Mayo) and `syn73713767` (MSBB) individual metadata are on the critical path for `pmi` / `pmiUnits` only.** Previously classified out of scope as clinical; decision D‑J changes that for post‑mortem interval alone. They must be read under a strict `usecols` allow‑list — SEC‑1 narrows for PMI, SEC‑2 does not change. Without them MSBB and Mayo carry `pmi_hours` null as `source_not_acquired`. DiverseCohorts has **no** individual‑metadata object registered at all, so its PMI source needs identifying before that stratum can be populated. |
| IN‑4 | `ROSMAP_clinical.csv` (syn3191087) must be present for the ROSMAP `projid`↔`individualID` crosswalk. Its absence is **non‑fatal** — the build logs and continues without the crosswalk ([build.py:332-343](src/ampprot/build.py#L332-L343)). |
| IN‑5 | Source files are read **read‑only**. The build must never modify, move, or delete anything in the source tree. |
| IN‑6 | Filesystem `mtime` is the *download* time, not file creation, and must never be published as a creation date ([manifest.py:122](src/ampprot/manifest.py#L122)). |

## 3. Data handling and disclosure (`SEC`)

This is the controlling section. The inputs are **controlled-access AMP-AD and AMP-PD
data** used under a Data Use Certificate; the derived HDF5 inherits that status.

| id | requirement |
|---|---|
| SEC‑1 | **No clinical or demographic attribute may enter any emitted artifact.** Out of scope: diagnosis, cognition, Braak/CERAD, APOE, education, age, sex, race, ethnicity (plan §0). **Amended — post‑mortem interval is in scope as a CDE.** PMI describes the *specimen*, not the donor, and is a standard technical covariate for brain proteomics, so it is carried as a biospecimen attribute rather than treated as phenotype (decision D‑J). This is a deliberate narrowing of the original boundary and is recorded as such: `pmi` / `postmortem_interval` leave `CLINICAL_DENY`, and the code comment at that constant must say why. Everything else in the list is unchanged. |
| SEC‑2 | Clinical payload must be discarded **at read time**, so it never exists in an in‑memory harmonised object. Every annotation reader must select an explicit allow‑list of columns rather than dropping a deny‑list. |
| SEC‑3 | A **build gate** must re‑check every emitted column name and *fail the build* on a clinical match — not warn ([build.py:299-302](src/ampprot/build.py#L299-L302), [writer.py:223-234](src/ampprot/writer.py#L223-L234)). |
| SEC‑4 | Participant identifiers are carried **as the source defines them**, never re-keyed or pseudonymised ([harmonize.py:411-429](src/ampprot/harmonize.py#L411-L429)). The output is therefore re‑identifiable against the source studies and is **not** de‑identified. |
| SEC‑5 | **The HDF5 must never be committed, pushed, or shared outside the DUC.** It contains donor‑level identifiers for every row, plus a ROSMAP identifier crosswalk. **This extends unchanged to `build/tabs/` and `build/returns/`** — the analysis‑ready bundles and the per‑grant return path are projections of the same rows and carry the same identifiers, so they are controlled artifacts and not a publishable by‑product. Both are denied in [.gitignore](.gitignore) by name as well as by directory. |
| SEC‑6 | **Source data must never be committed.** Enforced deny‑by‑default in [.gitignore](.gitignore). |
| SEC‑7 | Publishable artifacts are limited to: `src/`, `HARMONIZATION_PLAN.md`, `REQUIREMENTS.md`, [AUDIT_RESPONSE_PLAN.md](AUDIT_RESPONSE_PLAN.md), [AUDIT_REPLY.md](AUDIT_REPLY.md), `requirements.txt`, `build/README.md`, `build/build_report.json`, `build/inventory/*.tsv`, `build/inventory/*.json`. Each must be free of participant identifiers. **Every addition to this list requires its own disclosure scan before first release** — see §8. `build/tabs/` and `build/returns/` are **not** on this list and never will be: they carry per‑sample donor identifiers (SEC‑5). |
| SEC‑8 | Inventory and report artifacts must describe sample‑key *structure*, never a real key. `SAMPLE_KEY_PATTERN` ([build.py:73-81](src/ampprot/build.py#L73-L81)) publishes the format so the shape is documented without carrying an identifier. |
| SEC‑9 | `ROSMAP_clinical.csv` (syn3191087) is read twice, each time through an explicit `usecols` allow‑list, and never otherwise: once for the crosswalk (`projid`, `individualID`, `Study`) and once for the specimen CDE (`projid`, `pmi`). `pmi` is the **single** column re‑admitted by D‑J and SEC‑1; the rest of the file — `msex`, `educ`, `race`, `apoe_genotype`, `age_death`, `braaksc`, `ceradsc`, `cogdx`, `mmse` — is never loaded into memory. |

## 4. Functional (`FR`)

| id | requirement |
|---|---|
| FR‑1 | One keyed HDF5 group per (assay × biospecimen matrix) layer; layers declared but not populated are still created, with the reason in attrs ([writer.py:195-203](src/ampprot/writer.py#L195-L203)). |
| FR‑2 | Every feature resolves to a **UniProt primary accession**; unresolvable features are dropped and counted ([harmonize.py:138-237](src/ampprot/harmonize.py#L138-L237)). |
| FR‑3 | Missingness filter — keep features present in ≥ 50 % of samples — applied **per stratum, never after merge** ([harmonize.py:244-256](src/ampprot/harmonize.py#L244-L256)). |
| FR‑4 | Duplicate accessions collapse by mean **except** for SRM (peptides) and SomaScan (aptamers), where the measured entity is finer than the protein ([harmonize.py:35](src/ampprot/harmonize.py#L35)). |
| FR‑5 | Transforms must be **monotone and declared**. log2 only where the source ships linear; rank inverse‑normal only where measured skew exceeds tolerance *and* the layer is being Z‑scored ([harmonize.py:608-637](src/ampprot/harmonize.py#L608-L637)). |
| FR‑6 | Z‑scaling is decided on **measured distributions, not scale labels**, per feature, within stratum, before merge. A single‑stratum layer keeps native scale ([harmonize.py:86-120](src/ampprot/harmonize.py#L86-L120)). |
| FR‑7 | Any disagreement between a declared scale label and the measured distribution must be surfaced, not silently resolved ([harmonize.py:96-102](src/ampprot/harmonize.py#L96-L102)). |
| FR‑8 | `person_id` **is** the source identifier. Prefixes are never stripped. Namespace is recorded per row, and `person_global_key` (`namespace:id`) is the cross‑study join key ([harmonize.py:432-454](src/ampprot/harmonize.py#L432-L454)). |
| FR‑9 | `row_key` must be unique within a layer, composed from real discriminators (person → site → visit → stratum) before any anonymous counter ([writer.py:70-118](src/ampprot/writer.py#L70-L118)). |
| FR‑10 | `values`/`values_z` carry **exactly one null type: float NaN**. Missing *reason* lives in a separate int8 companion array ([harmonize.py:588-605](src/ampprot/harmonize.py#L588-L605)). |
| FR‑11 | Every layer exposes the same CDE and link column surface under the same names; layer‑specific proteomics columns appear only where the assay has them ([config.py:61-121](src/ampprot/config.py#L61-L121)). |
| FR‑12 | The README is **generated from the produced artifact**, never hand‑written prose about intent ([readme.py](src/ampprot/readme.py)). |
| FR‑13 | `samples.file_id` is **one‑to‑many** wherever a stratum's matrix is several files — Olink CSF ships 4 panel files, plasma 4 after WP‑11. The relation is carried as a **`sample × file` bridge table**, never a delimited list column, so `samples` stays one row per sample and the CDM's `ASSAY_INPUT_FILE` maps onto it directly (decision D‑F). `file_role` travels on the bridge row, so selecting matrix files is a filter rather than a filename match. |
| FR‑15 | **Post‑mortem interval is reported on one scale: `pmi_hours`, float, hours with decimals.** Units differ between studies — ADKP ships `pmiUnits` beside `pmi` because MSBB records minutes by convention and ROSMAP hours — so a bare `pmi` column is never trusted. Each value reaches hours either by **conversion** from a stated unit or by **inference** from the observed distribution (PMI is death‑to‑processing, so the plausible band is hours to a few days); which one applies is recorded per row in `pmi_evidence`, alongside `pmi_source_unit`. Ambiguous distributions stay null and are reported undetermined, never guessed. Null carries meaning and must be distinguishable: `not_applicable_antemortem` for fluid layers from living donors is a fact about the specimen, `source_not_acquired` is a gap in what we hold. **Never zero.** A per‑stratum conversion ledger appears in the generated README, and an audit row checks each post‑conversion median against the plausible band — reporting, not auto‑correcting (FR‑7). |
| FR‑14 | **`visit_month` — months since baseline — is the canonical cross‑layer time axis**, present in every layer, populated where the dataset is longitudinal and null for postmortem (decision D‑H). `visit_index` is retained as the source‑native value and is never the join key; `visit_occurrence_id` remains the OMOP link. Where `visit_month` is derived rather than stated, `visit_month_evidence` records which, and the derivation must be declared and reversible. Adding a longitudinal dataset means populating this column, not introducing another axis. |

## 5. Outputs (`OUT`)

| id | requirement |
|---|---|
| OUT‑1 | `build/sysbio_proteomics-<date>.h5` — the harmonised dataset. Controlled access (SEC‑5). **The date is computed at build time, so every build produces its own artifact and no build can overwrite an audited one** (decision D‑I). *Implemented* — `H5_NAME` is derived from the build date ([config.py](src/ampprot/config.py)) and the destructive unlink that the hardcoded literal made necessary is gone, which retires OPS‑1 rather than patching it; a write‑to‑temp‑then‑rename remains so a same‑day rerun cannot destroy a good artifact mid‑write. **The filename is the version** — consumers resolve the latest or pin deliberately, and must not assume a stable name. Each report and generated README names the artifact it describes, and the report also names the preceding build so lineage is explicit. |
| OUT‑2 | `build/build_report.json` — machine‑readable build summary. |
| OUT‑3 | `build/README.md` — generated documentation. |
| OUT‑4 | `build/inventory/{derived_inventory,tracking_file_diff,missing_objects,file_manifest,assertion_audit}.tsv` and `inventory_summary.json`. |
| OUT‑5 | Matrices stored wide (rows = samples), float32, gzip‑4 with shuffle and chunking ([writer.py:64-67](src/ampprot/writer.py#L64-L67)). |
| OUT‑6 | FILES CDEs that cannot be honestly populated (`current_version`, `created_on`, `modified_on`, `drs_id`) are left **null**, never guessed. |

## 6. Verification (`VER`)

| id | requirement |
|---|---|
| VER‑1 | Every claim — tracking file, bundled README, or encoded in `config.py` — is re‑checked against the data; disagreements are reported, not resolved ([verify.py](src/ampprot/verify.py)). |
| VER‑2 | A claim that cannot be checked from data is `UNVERIFIABLE`, never a silent `PASS`. |
| VER‑3 | Declared scale must match observed p1/p99 range ([verify.py:26-34](src/ampprot/verify.py#L26-L34)). |
| VER‑4 | Every non‑pool sample must resolve to a person, or the check fails ([verify.py:71-79](src/ampprot/verify.py#L71-L79)). |
| VER‑5 | Post‑condition, not assumption: all strata in a layer occupy one numeric space **on the surface that actually ships** ([harmonize.py:639-676](src/ampprot/harmonize.py#L639-L676)). |
| VER‑6 | Matrix shape must match the feature and sample table lengths. |
| VER‑7 | The audit is written to `assertion_audit.tsv` and summarised in the report and README. |

## 7. Operations (`OPS`)

| id | requirement |
|---|---|
| OPS‑1 | The build is **deterministic** given the same source tree and pins. |
| OPS‑2 | The metadata pass (manifest, missing objects) is cheap and re‑runnable after any download, so "what is still blocked" stays current ([manifest.py:1-5](src/ampprot/manifest.py#L1-L5)). |
| OPS‑3 | Files above 64 MB are **prefix**‑hashed and the manifest labels them as such, rather than implying a full read ([manifest.py:66-76](src/ampprot/manifest.py#L66-L76)). Row counts are omitted, not estimated, above 32 MB compressed. |
| OPS‑4 | A scope‑gate hit must abort the build before anything is written for that layer. |

---

## 8. Conformance status — 2026‑09‑04 build

**Produced:** 255.2 MB HDF5 `sysbio_proteomics-09042026.h5`, **8 of 8 layers populated**, 14 strata
built, none deferred. Plus 8 `tabs/` bundles and 2 per‑grant returns.
**Audit:** 313 checks — 242 PASS, 18 FAIL, 53 UNVERIFIABLE.

> **The FAIL count rose from 10 to 18, and that is the point.** The previous figure was
> understated: 12 of its 154 rows were hardcoded literals rather than measurements (VER‑1), one of
> which published a verdict the metadata later disproved (VER‑1b), and the six `column_naming` rows
> compared in the wrong direction — certifying a guarantee that did not hold while flagging the one
> layer that did (VER‑7, FR‑11b). Every verdict is now computed. Of the 18 FAILs, **11 are
> tracking‑file disagreements** where the data is right and the hand‑maintained file is stale, and
> the other 7 are recorded findings the build reports rather than resolves: two upstream QC claims
> that do not hold on the delivered files, the ROSMAP R1 platform conflict, three identity‑based
> row removals, and the DIA completeness drop. None is an unexplained failure.

| layer | rows | proteins | strata | surface |
|---|---|---|---|---|
| `tmt_brain_proteomics` | 2,160 | 14,916 | 5 | `values_z` |
| `lfq_brain_proteomics` | 536 | 4,437 | 2 | `values` |
| `srm_brain_proteomics` | 2,704 | 185 | 2 | `values` |
| `soma_plasma_proteomics` | 970 | 7,289 | 1 | `values` |
| `olink_plasma_proteomics` | 1,380 | 1,463 | 1 | `values` |
| `olink_csf_proteomics` | 1,320 | 1,463 | 1 | `values` |
| `dia_plasma_proteomics` | 1,459 | 202 | 1 | `values` |
| `dia_csf_proteomics` | 2,807 | 249 | 1 | `values` |

### What changed against the 2026‑08‑22 build

Verified by direct comparison of the two artifacts, which is possible **because** D‑I stopped a
rebuild overwriting its predecessor:

- **Seven of eight layers are bit‑identical** on every shared row and feature — max |difference| 0,
  no NaN‑pattern changes.
- `tmt_brain_proteomics` native `values` are **also bit‑identical** on all 2,156 shared rows. Its
  `values_z` shifts by at most 0.277, confined entirely to `diversecohorts_dlpfc` — the one stratum
  that lost samples, so its per‑feature mean and SD are recomputed over 1,082 rather than 1,086
  samples. Every other stratum's Z is unchanged to the bit.
- **Row counts:** `soma_plasma` 973 → 970, `tmt_brain` 2,164 → 2,160. Both are the identity‑based
  removals; no other layer changed.
- **Columns:** 7 added to `samples` (`pmi_hours`, `pmi_evidence`, `pmi_source_unit`,
  `brodmann_area_evidence`, `qc_detail`, `qc_warn_count`, `exclude_reason`) and 2 to `features`
  (`colcheck`, `missing_freq`). **None removed.**
- **8 `row_key` values changed**, and all 8 are the two defects being corrected: four
  `Unknown|…|msbb_tmt_phg` keys became `UNRESOLVED|…` when the sentinel stopped being read as a
  donor, and four `UNRESOLVED|…|diversecohorts_dlpfc` keys are gone with the rows they named.

### Disclosure review of the publishable set (SEC‑7)

**[corrected] — the 2026‑08‑22 review recorded "no participant identifiers were found", and that
was wrong.** Three publishable artifacts carried AMP‑PD participant identifiers via the DIA
completeness audit row (DISC‑1). The scan that cleared them did not include an AMP‑PD participant
pattern, so it passed by not looking. That is the failure mode a disclosure scan has to be written
against: an absent rule reads exactly like a clean result.

The scan is now run over the **full** publishable set — 23 files — plus the artifact itself, and
covers four separate things because they fail differently:

1. **identifier‑shaped tokens** in publishable text — `PD-*` / `PP-*`, `R<7+ digits>`,
   `AMPAD_MSSM_*`, and bare 8‑digit projids, with `syn` accessions and hex digests stripped first
   so a real projid cannot hide behind one;
2. **clinical column names in the artifact**, checked independently of the build's own scope gate
   rather than by trusting it;
3. **clinical values in the artifact** — a column can be renamed innocuously and still carry
   phenotype, which a name‑based gate cannot see;
4. **operational leakage** — absolute local paths, email addresses, credentials.

Result on the 2026‑09‑04 build: **PASS on all four.** Separately, the scope gate was re‑verified
against the 19 real clinical column names SEC‑3 was raised about — including the 10 the old
boundary rule let through — and catches all 19, while passing all 11 newly added columns, none of
which are undeclared.

**A `.gitignore` inconsistency was found alongside it:** `AUDIT_RESPONSE_PLAN.md` was named in
SEC‑7 as publishable but was never re‑included after the deny‑by‑default rule, so it would have
been silently excluded from any push. SEC‑7 and `.gitignore` are now consistent, and both name
`AUDIT_REPLY.md`.

`file_manifest.tsv` does publish the full 707‑file listing with Synapse IDs and sha256
prefixes. Synapse accessions are public metadata, so this is judged acceptable, but it does
disclose exactly which controlled objects were obtained — confirm that is acceptable before
the first push.

### Non‑conformances

**Resolution status as of the 2026‑09‑04 build.** Every open item below is closed. The findings
themselves are left as written — they are the record of what was wrong and why — with the outcome
stated here rather than by editing the finding to look as though it never happened.

| ref | was | now | how it was closed |
|---|---|---|---|
| **SEC‑3** | open — high | ✅ resolved | Scope gate matches on substring after stripping non‑alphanumerics; `reagan`, `dcfdx`, `mmse30`, `spanish`, `samplingage` added. Introduced in `report_only` mode over the full emitted column inventory first, and `DENY_EXEMPT` extended for the one legitimate column it newly caught (`reagent_lot` — "re‑AGE‑nt") |
| **SEC‑3** | open — low | ✅ resolved | The allow‑list half is implemented **and now runs**. `undeclared_columns` existed as a function that nothing called, so no allow‑list report was ever produced; it is wired to an audit row that reports — never enforces, per FR‑7 |
| **SEC‑2** | open — low | ⏸ unchanged, verified | `annot_srm_features` still returns the whole file. It reaches `feature_annot` only and `resolve_features` rebuilds an explicit frame, so nothing propagates. Re‑checked: peptide‑level, no clinical columns |
| **VER‑1** | open — high | ✅ resolved | `verify_tracking_file` and `verify_readme_claims` take the built objects and compute. **No verdict in the audit is a literal.** The tracking rows are generated from `tracking_file_diff.tsv` |
| **VER‑1b** | open — high | ✅ resolved | The tracking file was right. `anatomic_site` for the temporal stratum is `superior temporal gyrus`, the contradicting note is deleted, the raw value is kept in `extra` so the correction is reversible, and the hardcoded FAIL is gone with the rest of VER‑1 |
| **VER‑7** | open — high | ✅ resolved | `column_naming` compares against the **union** across layers, so a missing column is a FAIL; `feature_columns` performs a real comparison, one row per layer |
| **FR‑10** | open — high | ✅ resolved | `below_lod` is NaN where `values` is NaN, giving three separable states, documented in the generated README with the reason the old encoding changed results |
| **FR‑11** | open — medium | ✅ resolved | `samples.file_id` resolves into FILES; where a stratum has several matrix files it is null by design and the `sample_files` bridge is authoritative, so a populated `file_id` always means a genuine 1:1 |
| **FR‑11b** | open — high | ✅ resolved | Runtime column union across every layer, computed in the writer rather than from config. Union: `samples` 63, `features` 16 — every layer carries all of them |
| **DATA‑1** | open — high | ✅ closed as UNVERIFIABLE | The two SomaScan objects share md5 `f8a3d03ecd1b57cbb409bc6b054a7890`, so the ANML label has no evidence either way. Values unaffected; the label is reported UNVERIFIABLE with the md5 recorded. **Raised with ADKP — answer outstanding** |
| **DATA‑2** | open — high | ✅ resolved | All 172 disputed samples withdrawn, by **plate** rather than by id list, so the rule survives a re‑release |
| **DATA‑3** | open — medium | ✅ resolved | The QC‑recovery merge is deleted. Removing the 172 removed its reason to exist, and the false premise it published is gone with it |
| **CDE‑1** | open — medium | ✅ resolved | The five assay‑metadata objects are read. `platform` is per sample at evidence **A** for 6 strata, and a computed `cde_conflict` row reports config‑vs‑metadata disagreement without resolving it |
| **CDE‑2** | resolved by data | ✅ built | Now actually spent rather than merely available: `platform` closes at evidence **A** for **6 of the 8** strata that lacked it. ROSMAP SRM's two do not — `syn23569441` carries `platform = NA` on all 1,212 rows, reported as a measured silence in the deposit rather than a gap in the build. **One trap worth publishing:** the assay‑metadata objects and the biospecimen files do not share an assay vocabulary. MSBB's assay object says `LC-MSMS` where its biospecimen file says `label free mass spectrometry`, so filtering the assay object by the biospecimen file's term matches **0 of 306** rows and costs that stratum its platform. A filter that removes everything is indistinguishable from a source that contains nothing, which is why the values found in each object are now recorded rather than filtered against |
| **OPS‑1** | open — low | ✅ retired | The destructive unlink is **gone**, not guarded. `H5_NAME` is datestamped (D‑I) so a build cannot target an audited artifact's name; a write‑to‑temp‑then‑rename covers a same‑day rerun |
| **OPS‑2b** | open — medium | ✅ resolved | Presence is detected by content hash against a known‑md5 table, not by path. 13 of 15 objects present; the 2 absent are the Mayo and MSBB individual‑metadata objects, now needed for `pmi` alone |
| **VER‑9** | open — medium | ✅ resolved | Every `SOURCE_QC` entry is checked. An `applied_upstream=True` claim is verified against the delivered file where the data can show it; an `applied_upstream=False` step FAILs unless the build left a trace. This immediately caught two upstream claims that do not hold — 2 features in DiverseCohorts frontal and **398** in temporal are below the 50% threshold the QC log says was applied |
| **INTEROP‑1** | open — medium | ✅ resolved | `visit_month` is present in every layer with `visit_month_evidence`; `n_timepoints`, `timepoint_values` and `is_longitudinal` are populated in `derived_inventory.tsv` |
| **QC‑1** | open — low | ✅ resolved | `exclude` / `excludeReason` are wired to `qc_status` and an audit row. Still 0 flagged in both files, so nothing changes today — but the flag is live rather than discarded |
| **BIO‑1** | open — low | ✅ resolved | Cross‑assay propagation recovers 123 BA9 + 67 BA10; the rest is convention fill, and `brodmann_area_evidence` keeps inference separable from measurement. The 67 source‑stated BA10 rows were **not** overwritten |
| **FR‑8b** | resolved — drop | ✅ built | The 3 SomaScan rows are dropped with the rule and count recorded. Generalised: **any** non‑pool sample with no donor is removed as an identity failure, which also caught 4 DiverseCohorts rows the build had been shipping mislabelled as reference pools |
| **ENV‑3** | open — low | ⏸ unchanged | `_available_bytes` still falls back to 4 GB without `/proc/meminfo`. Not hit on this platform |
| **ENV‑2**, **SEC‑6** | resolved | ✅ unchanged | Pins in `requirements.txt` (now including `pyarrow` for the bundles); deny‑by‑default `.gitignore`, extended to name `build/tabs/` and `build/returns/` |

**Two findings this round that were not in the audit**, both pre‑existing and both in the
2026‑08‑22 artifact as well:

| ref | severity | finding |
|---|---|---|
| **ID‑1** | medium | **An absent‑value sentinel was travelling as a donor identifier.** MSBB's biospecimen file spells an absent donor `Unknown`; `as_source_id` filtered `""`, `nan`, `NA`, `None` and `<NA>` but not that, so 4 rows of `msbb_tmt_phg` shipped `person_id = "Unknown"` and `person_global_key = "AMP-AD.individualID:Unknown"` — in this build's predecessor too. Contained, because all four are GIS reference channels rather than donors, but a sentinel in an identifier column merges unrelated specimens into one fake person the moment two strata carry it. **Resolved:** an explicit sentinel set, matched case‑insensitively. |
| **DISC‑1** | **high** | **Two publishable artifacts were disclosing participant identifiers.** The DIA sample‑completeness audit row named the samples it dropped; a DIA sample id is shaped `PP-<participant>-<visit>-PLA-PDIA`, so `assertion_audit.tsv`, `build_report.json` and the generated `README.md` each carried AMP‑PD participant identifiers — **in the 2026‑08‑22 build as well**. A newly added donor‑unresolvable row would have done the same with ROSMAP `projid_visit` values. SEC‑8 already required the inventory surface to publish sample‑key *structure* and never a real key, and `SAMPLE_KEY_PATTERN` implements that for the rest of it; these rows bypassed it. **Resolved:** both rows publish the rule and the counts only. The DIA row keeps its diagnostic force by reporting the detected‑protein counts (`[10, 70, 72, …] of 202`) instead of the sample ids. |
| **VER‑10** | low | **`person_resolution` could report a negative count.** It computed unresolved as `n_rows − n_pool − n_resolved`, which assumes a pool never carries a donor id. That stopped being true when `is_pool` became evidence‑based: the MSBB TMT assay metadata flags 4 GIS channels whose biospecimen row still names an individual, and the check reported "‑4 unresolved". **Resolved:** counted on the non‑pool rows directly. |

| ref | severity | finding |
|---|---|---|
| **SEC‑6** | *resolved* | No `.gitignore` and no git repository existed. `git add .` would have staged 83.6 GB of controlled data including `ROSMAP_clinical.csv` (msex, educ, race, apoe_genotype, age_death, braaksc, ceradsc, cogdx, pmi). A deny‑by‑default [.gitignore](.gitignore) is now in place and was verified against a mirror of the tree. |
| **SEC‑3** | **open — high** | The scope gate misses real clinical column names. Its boundary rule ([writer.py:230-231](src/ampprot/writer.py#L230-L231)) requires a deny token to be the whole name or `_`‑delimited, so it passes `braaksc`, `ceradsc`, `niareagansc`, `dcfdx_lv`, `cts_mmse30_lv`, `ad_reagan`, `spanish`, `AgeAtDeath` and `samplingAge` — 10 of 33 real column names probed, including the literal ROSMAP and AD Knowledge Portal spellings. Nothing leaks today only because SEC‑2 holds independently. The gate is documented as the backstop and does not currently hold. Suggested fix: match on substring after stripping non‑alphanumerics, and add `reagan`, `dcfdx`, `mmse30`, `spanish`, `sampling_age`. |
| **SEC‑3** | open — low | Plan §0 specifies "an allow‑list **and** a clinical deny‑list". Only the deny‑list is implemented; `DENY_EXEMPT` ([config.py:32-36](src/ampprot/config.py#L32-L36)) is an exemption set, not an allow‑list. |
| **SEC‑2** | open — low | `annot_srm_features` ([readers.py:569-572](src/ampprot/readers.py#L569-L572)) is the one reader that returns the whole file rather than an allow‑list. It reaches `feature_annot` only, and `resolve_features` rebuilds an explicit frame, so nothing propagates — but it is the single deviation from SEC‑2. The file was checked and is peptide‑level, with no clinical columns. |
| **VER‑1** | **open — high** | `verify_tracking_file` and `verify_readme_claims` ([verify.py:212-249](src/ampprot/verify.py#L212-L249)) take no arguments and return **hardcoded verdicts and observed values**. 12 of 154 audit rows — and **8 of the 10 FAILs** — are literals, not measurements. The module docstring and generated README both state every claim is re‑checked against the data. These entries will not change if the data does, so they can go stale silently. Only 1 computed FAIL exists (`person_resolution` on `rosmap_soma_plasma`). |
| **FR‑11** | open — medium | `samples.file_id` is built as `<stratum_key>:protein_abundance_matrix` ([build.py:202](src/ampprot/build.py#L202)) while `files.file_id` is built as `<stratum_key>:<filename>` ([build.py:152](src/ampprot/build.py#L152)). The two never match, so `file_id` is a dangling foreign key from the sample table into the FILES table. |
| **FR‑10** | open — medium | `below_lod` is written as `(native < lod)` cast to float32 ([writer.py:143](src/ampprot/writer.py#L143)). `NaN < NaN` is `False`, so unmeasured cells are stored as `0.0` — indistinguishable from a genuine "not below LOD". It should be NaN where `values` is NaN. |
| **VER‑7** | open — low | `verify_column_consistency` ([verify.py:188](src/ampprot/verify.py#L188)) reads `PASS if not extra else PASS` — a dead conditional that can never report a difference. If the intent is "report but do not fail", the branch should still be visible in the verdict. |
| **OPS‑1** | open — low | `C.H5_PATH.unlink()` deletes the previous output before the new build starts ([build.py:251-252](src/ampprot/build.py#L251-L252)). A build that fails partway destroys the prior artifact with no backup. |
| **ENV‑3** | open — low | On a platform without `/proc/meminfo`, `_available_bytes` silently returns 4 GB ([build.py:38-40](src/ampprot/build.py#L38-L40)). A wrong‑but‑plausible budget is harder to diagnose than a refusal to run. |
| **ENV‑2** | *resolved* | No `requirements.txt`, `pyproject.toml` or lockfile existed. Pins are now recorded in [requirements.txt](requirements.txt). |

### Findings from the 2026‑09‑02 external audit and its verification

Raised in `harmonization_audit_alan/`, plus what verifying it turned up. Work plan in
[AUDIT_RESPONSE_PLAN.md](AUDIT_RESPONSE_PLAN.md).

| ref | severity | finding |
|---|---|---|
| **FR‑11b** | **open — high** | **The uniform column contract does not hold.** Five of six layers omit columns plan §3 guarantees are present‑and‑null: `samples` spans 53 columns of which `lfq_brain` has 50 and `olink_csf`/`olink_plasma`/`srm_brain` have 47; `features` spans 14 and every layer is missing three to seven. [writer.py `_write_table`](src/ampprot/writer.py#L56) writes whatever columns the frame happens to hold. A reader written against the documented contract raises `KeyError` where it was promised a null. `tmt_brain_proteomics` is the only conforming layer. Repair needs `frac_below_lod` and `olink_lod_flag` added to `FEATURE_COLS_BY_LAYER` — they are emitted but undeclared — then a runtime column union rather than a config‑driven reindex. |
| **VER‑7** | **open — high** *(was low)* | Re‑graded. The dead ternary is the lesser half: [verify.py:186](src/ampprot/verify.py#L186) also compares in the **excess** direction (`set(layer) - common`), so columns a layer is *missing* are never examined. The consequence is that `assertion_audit.tsv` records the inverse of the truth — the three layers missing six columns each are certified "identical to all layers", and `tmt_brain_proteomics`, the only conforming layer, is the one reported as deviating. `feature_columns` is a hardcoded `PASS` performing no comparison at all. |
| **FR‑10** | **open — high** *(was medium)* | Re‑graded: this one propagates into analysis results. `below_lod == 0.0` is written both for a cell measured above the limit and for a cell never measured, so the conventional filter silently retains never‑measured cells. The correct predicate is `below_lod == 0.0 AND missing_reason == 0`, which is not discoverable from the surface. Fix is NaN where `values` is NaN, giving three separable states. |
| **DATA‑1** | **open — high** | **The SomaScan ANML and non‑ANML matrices are byte‑identical.** `OhNM2025_ROSMAP_plasma_Soma7k_protein_level_ANML_log10.csv.gz` and `..._log10.csv.gz` both decompress to 127,604,657 bytes, md5 `f8a3d03ecd1b57cbb409bc6b054a7890`; the corresponding Synapse objects `syn65471938` and `syn65471939` carry that same md5, so this is the deposit and not our download. `value_scale="log10_rfu_anml"` and `transform_chain="source:SomaLogic ANML-normalised RFU\|log10"` ([config.py:353-354](src/ampprot/config.py#L353-L354)) therefore rest on evidence that does not exist — one of the two objects is mislabelled and the data cannot say which. Values are unaffected; the scale *label* must drop to UNVERIFIABLE and the discrepancy be raised with ADKP. |
| **VER‑1b** | **open — high** | **A hardcoded verdict published a wrong result.** `syn51757645` records `tissue = superior temporal gyrus` for all 280 DiverseCohorts `TMT quantitation` temporal specimens (matrix has 278 columns). [config.py:205-206](src/ampprot/config.py#L205-L206) asserts "TCX = temporal cortex (confirmed). The tracking file's 'superior temporal gyrus' is incorrect for this stratum", and [verify.py:222-223](src/ampprot/verify.py#L222-L223) publishes that as a hardcoded FAIL against the tracking file. The tracking file was right. This is VER‑1's predicted failure mode realised, and it moves the OMOP anatomic‑site concept for 508 rows. |
| **CDE‑1** | open — medium | **`platform` is per‑batch, not per‑stratum, and two declared values conflict with the portal metadata.** `syn53185805` records three instruments across DiverseCohorts TMT — Lumos (848 channels), Exploris 240 (423, `mayo_*`), Eclipse (224, `mssm_*`) — and two kit sizes (TMT16 1,072 / TMT18 423), so `analysis_pipeline="FragPipe (TMT 18-plex)"` is wrong for 1,072 of 1,495 channels. `syn21323404` gives ROSMAP TMT R1 as OrbiTrap Fusion / SPS‑MS3 / CID followed by HCD / 60000, against config's "Orbitrap FTMS; HCD MS2" derived from the spectrum file. FR‑7's rule — surface the disagreement, do not resolve it — applies to CDEs and is not currently implemented for them. |
| **OPS‑2b** | open — medium | **Missing‑object detection is path‑based and reported six present files as absent.** Of the 13 objects the audit transferred, six were already in the tree byte‑identical under different accessions: `syn20827192` = `MayoRNAseq_biospecimen_metadata.csv`, `syn21893059` = `MSBB_biospecimen_metadata.csv`, and `syn65414914`/`syn65471938`/`syn65471939`/`syn65473039` under `syn64957327/`. [manifest.py `build_missing_objects`](src/ampprot/manifest.py#L151) looks for a Synapse‑ID path, so anything that arrives bundled under another accession reads as missing. Detection must be content‑hash based. Genuinely new: seven files. |
| **CDE‑2** | *resolved by data* | `ASSAY.platform` — the weakest CDE in the build, populatable for one of six layers from data alone — now closes for six of the eight strata that lacked it (`syn53185805`, `syn21323404`, `syn21893060` Q Exactive HF, `syn22344998` "Q Extrative Plus" *(sic)*, `syn23474101` Q Exactive Plus). ROSMAP SRM does **not** close: `syn23569441` carries `platform = NA` on all 1,212 rows. |
| **DATA‑2** | **open — high** | **Olink plasma: 172 samples have disputed donor identity, and the current build ships them.** The `PLA-PPEA-D03` release ships plain and `_retracted` variants that differ *only* on `BIOREP_pl1_Samplesheet` and `BIOREP_pl2_Samplesheet` — 0 differing rows elsewhere, in all four panels — where 172 samples carry different `participant_id` / `sample_id` / `visit_month`. NPX is untouched. Nothing in the delivered files states which variant is corrected; the `_retracted` choice at [config.py:381](src/ampprot/config.py#L381) is inferred from plate counts, not stated. **Resolution: withhold all 172** (11.1% of plasma samples; 39 of 413 participants lose every plasma sample) and build from plain. The set of ids is identical in both variants, so the removal does not depend on AMP-PD's answer, and the remaining 1,380 samples are variant‑independent. |
| **DATA‑3** | open — medium | **The Olink QC‑recovery merge misattributes QC for exactly those 172 samples.** [readers.py:297-299](src/ampprot/readers.py#L297-L299) states `_retracted` "has blanked [QC] to all‑PASS" and re‑joins QC from the plain files on `(sample_id, UniProt, panel)`. Measured: the QC blanking is in the **common‑format** files (`Cumulative_QC`, 1,151/636/933/861 `WARN` → 0), **not** in the `olink_explore_format` files the build reads, whose QC columns are identical across variants. Because the swap permuted sample ids between the two plates, the merge hands each of the 172 samples the verdict belonging to the other plate's well — 3,006 misattributed verdicts — and is a no‑op for the other 1,380. [verify.py:239](src/ampprot/verify.py#L239) publishes the false premise as a measured FAIL. Removing the 172 (DATA‑2) removes the merge's reason to exist. |
| **VER‑9** | open — medium | **`SOURCE_QC` is declared but never enforced.** The ~30 per‑stratum QC entries at [harmonize.py:262](src/ampprot/harmonize.py#L262) reach `notes["source_qc"]` and the generated README and stop there. Nothing verifies an `applied_upstream=True` claim against the data, and nothing confirms that an `applied_upstream=False` item — the ones the producer left to us — was actually applied. Three are currently unapplied and unflagged: SomaScan `ColCheck` per‑analyte QC, SRM control‑well flagging, and the Olink below‑LOD assay flags (the last is in fact computed, which is exactly why the others being absent is invisible). Same class of defect as VER‑1: a claim published without a measurement behind it. |
| **INTEROP‑1** | open — medium | **No comparable longitudinal axis across layers, and the inventory cannot report one.** `olink_plasma` / `olink_csf` carry `visit_month` (months since baseline); `soma_plasma` carries only `visit_index`, an ordinal follow‑up number; brain layers are postmortem. A consumer joining ROSMAP plasma to PDRD Olink on `person_id` has no shared time column, so nothing longitudinal crosses cohorts. This is the interoperability consequence of plan decision D10, which was left at its conservative default. **It does not require re‑admitting a clinical field:** ROSMAP `Visit` is the follow‑up year (observed 0–14+ across 973 samples), so `visit_month = Visit × 12` is a structural conversion of an ordinal index, not a demographic derivation, and never touches `age_at_visit` or the scope gate. Must be labelled an approximation — ROSMAP follow‑ups are nominally but not exactly annual — with `visit_index` retained. Compounding it: `n_timepoints`, `timepoint_values` and `is_longitudinal` are **blank for every stratum** in the shipped `derived_inventory.tsv`, because [manifest.py:198-200](src/ampprot/manifest.py#L198-L200) reads them from `probes` and [build.py:311-325](src/ampprot/build.py#L311-L325) never writes them — so HARMONIZATION_PLAN §1.5's promised "timepoint sets" are absent and the gap is invisible from the published inventory. |
| **QC‑1** | open — low | **ADKP `exclude` / `excludeReason` are read and discarded.** [readers.py:498](src/ampprot/readers.py#L498) keeps both columns from the biospecimen metadata and nothing downstream consumes them. Verified today as harmless — 0 of 474 MSBB and 0 of 229 Mayo proteomics specimens carry `exclude=TRUE` — but MSBB's `excludeReason` vocabulary includes `sample swap`, so an unwired exclusion flag is a live risk at the next metadata release rather than a cosmetic gap. |
| **FR‑8b** | **resolved — drop** | **3 SomaScan samples cannot resolve to a donor.** `rosmap_soma_plasma` leaves 3 of 973 non‑pool samples unresolved — the only genuine (non‑hardcoded) FAIL in the 2026‑08‑22 audit. Traced to 3 rows whose `projid` is empty in the source sample metadata, so this is a gap in the deposit, not a join defect. **Resolution: drop them**, with the rule and the count recorded, so they do not fall out silently as a join miss. FR‑8/VER‑4 cannot hold for a row with no donor. |
| **BIO‑1** | open — low | **Brodmann area is absent from source for every DiverseCohorts proteomics specimen** — `NA` on all 1,385 `TMT quantitation` rows; the column is populated only on `rnaSeq`/`wholeGenomeSeq`/`10x multiome`. Cross‑assay propagation on donor+tissue recovers 123 → BA9 and 67 → BA10 for DLPFC and nothing for temporal, a ceiling near 12% rather than the 90% requested. Mayo temporal cortex is `NA` throughout; ROSMAP has no biospecimen file on disk. Coverage to ~90% is reached by asserting `DLPFC → BA9` **by convention** (decision D‑B), which requires per‑row `brodmann_area_evidence` so the inference stays separable from measurement — the 67 source‑stated BA10 rows show the convention is not universally true. |

### Scope items closed without action

Two of HARMONIZATION_PLAN §9's decisions are deliberately **not** being taken up in this round, and
are recorded here so the silence is a choice rather than an oversight:

- **D6 — the 64 GB PSM/spectrum‑level files.** Nothing reads them; only acquisition descriptors were ever sampled from them, and those are already recorded in `config.py`. They stay on the volume. Revisit only if a stratum needs PSM‑level input. ENV‑5's disk figure includes them.
- **D12 — the proposed `detection_modality` ASSAY CDE.** The five assays differ in detection modality (mass spectrometry / aptamer / antibody) and Figure 2's ASSAY has no field for it, so the proposal stands — but it is a Task Force decision, not a local one. Modality remains implicit in `platform` until they rule. Not implemented, not assumed.

### Security review

No `eval`, `exec`, `subprocess`, network I/O, or deserialisation of untrusted input. The
only pickling is internal to `ProcessPoolExecutor`. Source files are opened read‑only.
`pandas.read_excel` via openpyxl parses vendor‑supplied `.xlsx.gz` from Synapse — acceptable
given the trusted provenance. No credentials, tokens, or absolute local paths appear in any
publishable artifact.
