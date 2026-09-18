"""Layer definitions, source registry, and the SysBio CDE vocabulary.

Everything the build needs to know about *what* exists lives here. Facts that can be
read off the data files are marked evidence="A"; facts transcribed from a bundled
README or script are "B"; unavailable is None. See HARMONIZATION_PLAN.md section 3.7.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build"
TABS = BUILD / "tabs"          # WP-13 / D-G: the per-layer analysis-ready bundles
RETURNS = BUILD / "returns"    # WP-13 step 4: the per-grant return path

# D-I. **The filename is the version.** Computed at build time rather than hardcoded, so
# a rebuild produces its own artifact and can never overwrite an audited one -- which is
# what retires OPS-1 instead of patching it. Consumers must not pin a filename: resolve
# the newest, or pin deliberately for reproducibility. `build_report.json` and the
# generated README both name the artifact they describe, so a bundle is self-identifying.
H5_STEM = "sysbio_proteomics"
BUILD_DATE = date.today()
H5_NAME = f"{H5_STEM}-{BUILD_DATE:%m%d%Y}.h5"
H5_PATH = BUILD / H5_NAME


def previous_build() -> str | None:
    """The most recent artifact from an earlier build, for the lineage record.

    Published in `build_report.json` so the chain is explicit rather than inferred from
    a directory listing. Only the name is taken -- the file itself is never touched.
    """
    others = sorted((p for p in BUILD.glob(f"{H5_STEM}-*.h5") if p.name != H5_NAME),
                    key=lambda p: p.stat().st_mtime, reverse=True)
    return others[0].name if others else None

# --------------------------------------------------------------------------------------
# Scope gate (plan section 0)
# --------------------------------------------------------------------------------------

# Substrings that must never appear in an emitted column name. Matched after stripping
# non-alphanumerics (SEC-3), so `age_at_visit`, `ageAtVisit` and `age.at.visit` are all
# caught by one token -- the previous word-boundary rule let a separator style through.
#
# `pmi` and `postmortem_interval` were REMOVED from this tuple by decision D-J. Post-
# mortem interval is a biospecimen technical attribute -- it describes the specimen, not
# the donor -- and is a standard covariate for brain proteomics, so it is carried as a
# CDE (`pmi_hours`, always hours with decimals) rather than denied as phenotype. That is
# a deliberate narrowing of plan section 0 and REQUIREMENTS SEC-1, approved 2026-09-03,
# and is recorded as a scope change in both documents and in the generated README. It is
# never to be re-added silently.
CLINICAL_DENY = (
    "dx", "diagnosis", "cogdx", "cogn", "braak", "cerad", "apoe", "educ",
    "age", "sex", "msex", "race", "ethnicity", "mmse", "updrs", "moca",
    "dementia", "reagan", "dcfdx", "mmse30", "spanish", "samplingage",
)

# Columns allowed even though a deny token appears as a substring. Matched on the whole
# normalised column name, so this exempts a specific column and never a pattern.
#
# `reagent_lot` is the only column in the emitted inventory that the substring rule newly
# caught -- "re-AGE-nt" -- when it was run in report_only mode over all 79 emitted
# columns before being switched live.
DENY_EXEMPT = frozenset({
    "assay_id", "assay_type", "assay_source_value", "analysis_type", "analysis_pipeline",
    "array_type", "analyte_type", "storage", "specimen_id", "specimen_source_value",
    "processing_status", "file_format", "message", "reagent_lot",
    # Not currently emitted, exempted ahead of the substring rule catching them.
    "percentage", "coverage", "average", "linkage", "trace", "pct_coverage",
})

# Columns the build legitimately emits that are not part of a layer's declared contract
# surface. Used only by the allow-list report (SEC-3), never by the deny gate.
EXTRA_DECLARED_COLS = frozenset({
    "row_key", "stratum_key", "replicate_index", "relative_path", "synapse_id", "sha256",
    "person_global_key",
    "person_source_namespace", "person_source_namespace_row", "source_dataset",
    "specimen_native_id", "anatomic_site", "anatomic_site_source", "brodmann_area",
    "brodmann_area_evidence", "reagent_lot", "feature_key", "frac_below_lod",
    "olink_lod_flag", "platform_evidence", "pipeline_evidence", "transform_applied",
    "zscored", "value_scale", "value_unit", "transform_chain", "mass_analyzer",
    "activation_type", "ms_order", "scan_filter", "study", "grant", "cohort",
    "visit_month_evidence", "visit_index", "qc_detail", "qc_warn_count", "colcheck",
    "pmi_hours", "pmi_evidence", "pmi_source_unit", "detection_rate", "n_detected",
    "missing_freq", "exclude_reason",
})

# WP-9 / D-B. Declared once, here, so the audit can check it and a consumer can undo it.
# Source-stated Brodmann areas always take precedence; this only fills where the source
# is silent, and every filled row is marked `convention` in `brodmann_area_evidence`.
# WP-11 / D-D. The AMP-PD plate swap put the donor identity of every sample on these two
# plates in dispute, and only AMP-PD can settle which assignment is right. They are
# withdrawn by PLATE -- the physical fact -- rather than by a sample-id list, so the rule
# survives a re-release. This is a deletion, not a QC flag, because what is wrong is the
# row's identity rather than its quality; the standing policy of flagging instead of
# dropping applies to quality, and is unchanged.
OLINK_WITHDRAWN_PLATES = frozenset({
    "BIOREP_pl1_Samplesheet", "BIOREP_pl2_Samplesheet",
})

BRODMANN_CONVENTION = {
    "dorsolateral prefrontal cortex": "BA9",
    "superior temporal gyrus": "BA22",
}

# --------------------------------------------------------------------------------------
# Post-mortem interval (WP-7 / D-J)
# --------------------------------------------------------------------------------------
#
# ONE SCALE: `pmi_hours`, float, hours with decimals -- always. Sub-hour precision is
# real (ROSMAP's median is 6.783 h) and rounding would discard resolution on a covariate
# whose whole use is fine-grained. A bare `pmi` column is never trusted, because units
# differ between studies by convention: MSBB records minutes, ROSMAP hours.
#
# Every stratum reaches hours one of two ways and which one is always recorded per row in
# `pmi_evidence`. NEVER zero -- a zero interval is a claim, not a null.
PMI_HOURS_PER_UNIT = {
    "hour": 1.0, "hours": 1.0, "h": 1.0, "hr": 1.0, "hrs": 1.0,
    "minute": 1 / 60, "minutes": 1 / 60, "min": 1 / 60, "mins": 1 / 60,
    "day": 24.0, "days": 24.0, "d": 24.0,
}

# `converted_from_stated_unit` -- the source names its unit; the conversion is exact and
#     reversible, and the unit as found is kept in `pmi_source_unit`.
# `inferred_from_distribution` -- it does not; the unit is read off the observed median
#     against the plausible band below. Ambiguous stays null, never guessed.
# `not_applicable_antemortem` -- a fact about the specimen: living donor, fluid draw.
# `source_not_acquired`       -- a gap in what we hold. A different null entirely, which
#     is why the two are distinguished per row rather than collapsed.
PMI_EVIDENCE = ("converted_from_stated_unit", "inferred_from_distribution",
                "not_applicable_antemortem", "source_not_acquired")

# Death to processing at a brain bank: hours to a few days. The same durations recorded
# in minutes land near 400, so the bands do not overlap and inference is safe. Used both
# to infer an unstated unit and, post-conversion, as the audit's plausibility check --
# which REPORTS out-of-band, never auto-corrects (FR-7).
PMI_PLAUSIBLE_HOURS = (0.5, 120.0)

# Where PMI comes from, per stratum. `None` means antemortem -- not a gap.
PMI_SOURCE: dict[str, dict | None] = {
    # On disk and already in hours; factor 1, no conversion. Keyed by ROSMAP projid,
    # which is exactly what these four strata carry as `person_source_value`. Read
    # through the same explicit usecols allow-list as every other annotation (SEC-2):
    # the file is overwhelmingly clinical and only `projid` + `pmi` are loaded.
    "rosmap_tmt_r1": {"path": "ROSMAP_clinical.csv", "key": "projid",
                      "column": "pmi", "unit": "hours"},
    "rosmap_tmt_r2": {"path": "ROSMAP_clinical.csv", "key": "projid",
                      "column": "pmi", "unit": "hours"},
    "rosmap_srm_panel1": {"path": "ROSMAP_clinical.csv", "key": "projid",
                          "column": "pmi", "unit": "hours"},
    "rosmap_srm_panel2": {"path": "ROSMAP_clinical.csv", "key": "projid",
                          "column": "pmi", "unit": "hours"},
    # Registered but NOT in the 13-object delivery. These two moved onto the critical
    # path for `pmi` + `pmiUnits` only when D-J reclassified PMI; nothing else in them
    # would ever load. Until they arrive the rows are `source_not_acquired`. MSBB is the
    # study that conventionally records minutes, so it is where conversion stops being
    # theoretical -- the code path is built and exercised the day the file lands.
    "msbb_tmt_phg": {"not_acquired": "syn73713767 (MSBB individual metadata)"},
    "msbb_lfq_pfc": {"not_acquired": "syn73713767 (MSBB individual metadata)"},
    "mayo_lfq_tcx": {"not_acquired": "syn73713766 (MayoRNAseq individual metadata)"},
    # No individual-metadata object is registered anywhere for DiverseCohorts and
    # syn51757645 carries no PMI column, so the source itself is still unidentified.
    "diversecohorts_dlpfc": {"not_acquired": "no DiverseCohorts individual-metadata "
                                             "object is registered; source unidentified"},
    "diversecohorts_temporal": {"not_acquired": "no DiverseCohorts individual-metadata "
                                                "object is registered; source unidentified"},
    # Antemortem fluid from living donors. Verified: no PMI-, postmortem-, death- or
    # autopsy-shaped column exists in any PDRD file or in the column glossary.
    "rosmap_soma_plasma": None,
    "pdrd_olink_plasma": None,
    "pdrd_olink_csf": None,
    "pdrd_dia_plasma": None,
    "pdrd_dia_csf": None,
}

# WP-10 / WP-12. A DIA sample whose detected-protein count falls below this fraction of
# the layer's feature set is dropped: the delivered product carries no sample-level QC of
# its own, and the failures are unambiguous -- one plasma run detected 10 of 202 proteins
# (5%), then a failed cluster at 35-42%, then a clean gap before the tail becomes
# continuous. Above this line the variation is ordinary, not failure. The same 50%
# convention `SOURCE_QC` already uses per feature. CSF loses nothing at this threshold;
# its worst sample is 90% complete.
DIA_MIN_SAMPLE_COMPLETENESS = 0.50

# --------------------------------------------------------------------------------------
# Uniform column contract (plan section 3)
# --------------------------------------------------------------------------------------

LINK_KEYS = [
    "person_id", "person_source_value",
    "specimen_id", "specimen_source_value",
    "visit_occurrence_id", "visit_name", "visit_month", "visit_index",
    "assay_id", "file_id",
]

ASSAY_CDES = [
    "assay_id", "assay_source_value", "assay_type", "platform",
    "suspension_type", "analyte_type", "analysis_pipeline",
]

FILES_CDES = [
    "file_id", "file_name", "current_version", "assay_id", "file_role", "study", "grant",
    "array_type", "analysis_type", "biosample_type", "tissue", "cell_type", "species",
    "processing_status", "file_format", "file_size_bytes", "created_on", "modified_on",
    "drs_id",
]

# Proteomics columns are LAYER-SPECIFIC: a column that does not apply to a layer is
# omitted entirely rather than carried as all-null. The CDEs above are the opposite --
# always present, null where unpopulated, so the CDM surface stays uniform.
FEATURE_COLS_CORE = [
    "feature_id", "feature_source_value", "gene_symbol", "protein_group",
    "protein_group_size", "measurement_concept_id",
]

# `frac_below_lod` and `olink_lod_flag` are declared here because the writer emits them;
# they were emitted while undeclared, which is why a config-driven reindex alone would
# have left the contract broken while the audit reported it fixed. Config is no longer
# the source of truth for the written column set -- the runtime union in build.py is
# (WP-3) -- but leaving these undeclared would keep config lying about what ships.
#
# `missing_freq` (Olink's own per-assay missingness) and `colcheck` (SomaLogic's
# per-analyte column QC) are WP-12: both were read and discarded while `SOURCE_QC`
# declared one of them ours to apply. They are carried as flags, never used to delete a
# feature -- flag, not drop (D8).
FEATURE_COLS_BY_LAYER = {
    "olink_plasma_proteomics": ["panel", "panel_lot_nr", "frac_below_lod",
                                "olink_lod_flag", "missing_freq"],
    "olink_csf_proteomics": ["panel", "panel_lot_nr", "frac_below_lod",
                             "olink_lod_flag", "missing_freq"],
    "soma_plasma_proteomics": ["seq_id", "dilution", "colcheck"],
    "srm_brain_proteomics": ["peptide_sequence"],
}

MEASUREMENT_COLS_CORE = ["value", "value_scale", "value_unit", "is_missing"]

MEASUREMENT_COLS_BY_LAYER = {
    "olink_plasma_proteomics": ["lod", "below_lod"],
    "olink_csf_proteomics": ["lod", "below_lod"],
}

TECHNICAL_COLS_CORE = ["batch", "is_pool", "qc_status"]

TECHNICAL_COLS_BY_LAYER = {
    "tmt_brain_proteomics": ["channel"],
    "olink_plasma_proteomics": ["plate_id"],
    "olink_csf_proteomics": ["plate_id"],
    "srm_brain_proteomics": ["plate_id"],
}


def _dedup(*groups):
    seen, out = set(), []
    for g in groups:
        for c in g:
            if c not in seen:
                seen.add(c)
                out.append(c)
    return out


def feature_cols(layer: str) -> list[str]:
    return _dedup(FEATURE_COLS_CORE, FEATURE_COLS_BY_LAYER.get(layer, []))


def measurement_cols(layer: str) -> list[str]:
    return _dedup(MEASUREMENT_COLS_CORE, MEASUREMENT_COLS_BY_LAYER.get(layer, []))


def technical_cols(layer: str) -> list[str]:
    return _dedup(TECHNICAL_COLS_CORE, TECHNICAL_COLS_BY_LAYER.get(layer, []))


def sample_cols(layer: str) -> list[str]:
    """Link keys + the full CDE surface (always) + layer-specific technical columns."""
    return _dedup(LINK_KEYS, ASSAY_CDES, FILES_CDES, technical_cols(layer))


def long_cols(layer: str) -> list[str]:
    return _dedup(sample_cols(layer), feature_cols(layer), measurement_cols(layer))


# --------------------------------------------------------------------------------------
# Source registry
# --------------------------------------------------------------------------------------

@dataclass
class Stratum:
    """One source pipeline: a study's proteomics matrix plus its annotation."""

    key: str                    # stratum identifier, unique within a layer
    layer: str                  # layer key
    study: str
    grant: str                  # AMP program
    cohort: str
    matrix_path: str            # relative to ROOT
    reader: str                 # dispatch name in readers.py
    tissue: str
    biosample_type: str
    anatomic_site: str
    value_scale: str            # scale AFTER our transform chain
    value_unit: str
    transform_chain: str
    visit_name: str
    # SysBio ASSAY CDEs
    assay_type: str
    platform: str | None
    platform_evidence: str | None
    suspension_type: str
    analyte_type: str
    analysis_pipeline: str | None
    pipeline_evidence: str | None
    # annotation
    annot_path: str | None = None
    annot_reader: str | None = None
    person_resolvable: bool = True
    blocked_reason: str | None = None
    array_type: str | None = None
    extra: dict = field(default_factory=dict)


HUMAN = "Homo sapiens"
BULK = "bulk tissue"
NA = "not applicable"

_DC = "syn53185479/syn59611693"
_R1 = "syn17015098/syn32539359/syn21275338/syn21261728"
_PDRD_C = "PDRD/proteomics-CSF-PPEA-D03/olink-explore/protein-expression"
_PDRD_P = "PDRD/proteomics-PLA-PPEA-D03/olink-explore/protein-expression"
_PDRD_DIA_C = "PDRD/proteomics-CSF-PDIA"
_PDRD_DIA_P = "PDRD/proteomics-PLA-PDIA"

# Assay- and biospecimen-metadata objects landed by WP-1, stored under their own
# Synapse-ID folders. Referenced on the strata they serve so the readers can spend them
# (platform per sample, reagent lot, pool flags -- WP-8; anatomic site and Brodmann
# area -- WP-9) instead of leaving the CDEs at config-declared evidence.
_META_DC_ASSAY = "syn53185805/AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata_260622.csv.gz"
_META_DC_BIOSPEC = "syn51757645/AMP-AD_DiverseCohorts_biospecimen_metadata.csv.gz"
_META_R1_ASSAY = "syn21323404/ROSMAP_assay_proteomics_TMTquantitation_metadata.csv.gz"
_META_SRM_ASSAY = "syn23569441/ROSMAP_assay_proteomics_metadata.csv.gz"
_META_MSBB_TMT_ASSAY = "syn21893060/MSBB_assay_TMT_metadata.csv.gz"
_META_MSBB_LFQ_ASSAY = "syn22344998/MSBB_assay_proteomics_metadata.csv.gz"
_META_MAYO_ASSAY = "syn23474101/MayoRNAseq_assay_proteomics_metadata.csv.gz"
_META_SOMA_FEATURES = "syn64957327/OhNM2025_ROSMAP_plasma_Soma7k_protein_metadata.csv.gz"

STRATA: list[Stratum] = [
    # ---------------- tmt_brain_proteomics ----------------
    Stratum(
        key="diversecohorts_dlpfc", layer="tmt_brain_proteomics",
        study="AMP-AD_DiverseCohorts", grant="AMP AD", cohort="multi-site",
        matrix_path=f"{_DC}/syn55249982/n1086_residual_log2_batch.csv.gz",
        reader="csv_matrix",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="dorsolateral prefrontal cortex",
        value_scale="log2_ratio_total_batch_residual", value_unit="log2 ratio",
        transform_chain="source:log2(protein/total)|batch_regressed",
        visit_name="postmortem",
        assay_type="TMT LC-MS/MS", platform=None, platform_evidence=None,
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline="FragPipe (TMT 18-plex)", pipeline_evidence="A",
        annot_path=f"{_DC}/syn55249982/final_traits_to_be_published_frontal.xlsx.gz",
        annot_reader="dc_traits",
        extra={"assay_meta": _META_DC_ASSAY, "assay_meta_key": ("specimenID",
                                                               "specimen_source_value"),
               "biospecimen_meta": _META_DC_BIOSPEC},
    ),
    Stratum(
        key="diversecohorts_temporal", layer="tmt_brain_proteomics",
        study="AMP-AD_DiverseCohorts", grant="AMP AD", cohort="Mayo + Emory",
        matrix_path=f"{_DC}/syn55249982/n278_residual_log2_batch.TCX.csv.gz",
        reader="csv_matrix",
        tissue="brain", biosample_type="postmortem tissue",
        # N2 (WP-9): syn51757645 gives `tissue = superior temporal gyrus` for all 280
        # TMT-quantitation specimens in this set. The build previously asserted
        # "temporal cortex" and carried a note calling the tracking file wrong; the
        # tracking file was right. The source label is kept in `extra` so the
        # correction is reversible, and the OMOP concept for these 508 rows moves off
        # 4193043 Region of temporal cortex -- flagged to the CDM side.
        anatomic_site="superior temporal gyrus",
        value_scale="log2_ratio_total_batch_residual", value_unit="log2 ratio",
        transform_chain="source:log2(protein/total)|batch_regressed",
        visit_name="postmortem",
        assay_type="TMT LC-MS/MS", platform=None, platform_evidence=None,
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline="FragPipe (TMT 18-plex)", pipeline_evidence="A",
        annot_path=f"{_DC}/syn55249982/final_traits_to_be_published_temporal.xlsx.gz",
        annot_reader="dc_traits",
        extra={"assay_meta": _META_DC_ASSAY, "assay_meta_key": ("specimenID",
                                                               "specimen_source_value"),
               "biospecimen_meta": _META_DC_BIOSPEC,
               "anatomic_site_matrix_code": "TCX",
               "anatomic_site_corrected_from": "temporal cortex"},
    ),
    Stratum(
        key="rosmap_tmt_r1", layer="tmt_brain_proteomics",
        study="ROSMAP", grant="AMP AD", cohort="ROS/MAP",
        matrix_path=(f"{_R1}/syn21266446/"
                     "C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv.gz"),
        reader="csv_matrix",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="dorsolateral prefrontal cortex",
        value_scale="log2_rel_batch_median", value_unit="log2 ratio",
        transform_chain="source:TAMPOR median-polish log2 ratio centred on batch medians",
        visit_name="postmortem",
        assay_type="TMT LC-MS/MS",
        platform="Orbitrap FTMS; HCD MS2 (model not recorded)", platform_evidence="partial",
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline="Proteome Discoverer 2.3.0.522 + TAMPOR", pipeline_evidence="A",
        extra={
            # From rushtmt_fullmulticonsensus_MSMSSpectrumInfo.txt (header region sampled,
            # 400k of ~74M spectra; uniform across the sample). The spectrum-level file
            # itself is NOT ingested -- only these acquisition descriptors.
            "mass_analyzer": "FourierTransform",
            "activation_type": "HCD",
            "ms_order": "MS2",
            "scan_filter": "FTMS + p NSI Full ms [350.0000-1500.0000]",
            "acquisition_evidence": "A (sampled from MSMSSpectrumInfo + SpecializedTraces)",
            # D-C: syn21323404 is trusted over the spectrum-file descriptors above, and
            # the disagreement is reported as a computed `cde_conflict` row (WP-8) rather
            # than silently dropped.
            #
            # Joined on `specimenID` (shaped `ROSMAP.DLPFC.b<batch>.<channel>.<donor>`),
            # NOT on `batchChannel`: the file covers both ROSMAP TMT rounds, and 126 of its
            # batchChannel values carry two different instruments, so that key would
            # quietly import round 2's Q Exactive HF-X into round 1. On specimenID the
            # join is 400/400 and unique.
            #
            # Reading it per sample also narrows D-C itself. D-C recorded the file as
            # saying SPS-MS3 / CID-then-HCD / 60000 for round 1; that is true of 40 of
            # the 400 samples. The other 360 read HCD / MS2 / 30000 -- which AGREES with
            # the spectrum-file descriptors above. So the conflict is real but local, not
            # stratum-wide, and the audit row states the measured split rather than the
            # blanket claim.
            "assay_meta": _META_R1_ASSAY,
            "assay_meta_key": ("specimenID", "specimen_native_id"),
        },
        annot_path=f"{_R1}/syn21266445/rosmap_50batch_specimen_metadata_for_batch_correction.csv.gz",
        annot_reader="rosmap_r1_specimen",
    ),
    Stratum(
        key="msbb_tmt_phg", layer="tmt_brain_proteomics",
        study="MSBB", grant="AMP AD", cohort="Mount Sinai Brain Bank",
        matrix_path="syn21347564/syn24995077/msbb_19batch_normalized_tmt_matrix.xlsx.gz",
        reader="msbb_tmt_xlsx",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="parahippocampal gyrus",
        value_scale="log2_reporter_intensity", value_unit="log2 intensity",
        transform_chain="source:summed TMT reporter ion intensity|log2",
        visit_name="postmortem",
        assay_type="TMT LC-MS/MS", platform=None, platform_evidence=None,
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline=None, pipeline_evidence=None,
        annot_path="MSBB_biospecimen_metadata.csv", annot_reader="adkp_biospecimen",
        person_resolvable=True,
        extra={"assay_filter": "TMT quantitation", "assay_meta": _META_MSBB_TMT_ASSAY,
               "assay_meta_key": ("specimenID", "specimen_source_value")},
    ),

    Stratum(
        key="rosmap_tmt_r2", layer="tmt_brain_proteomics",
        study="ROSMAP", grant="AMP AD", cohort="ROS/MAP",
        matrix_path="syn17015098/syn30390636/syn26051783/"
                    "rosmaptmt_r2_14batch_multiconcensus_Proteins.txt.gz",
        reader="pd_ratio",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="dorsolateral prefrontal cortex",
        value_scale="log2_rel_batch_median", value_unit="log2 ratio",
        transform_chain="source:PD abundance ratio vs 126 GIS|log2|"
                        "centre on per-protein median of per-plex medians",
        visit_name="postmortem",
        assay_type="TMT LC-MS/MS",
        platform="Q Exactive HF-X Orbitrap", platform_evidence="A",
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline="Proteome Discoverer", pipeline_evidence="A",
        annot_path="syn17015098/syn30390636/ROSMAP_Round2_Traits_FINAL.xlsx.gz",
        annot_reader="rosmap_r2_traits",
    ),

    # ---------------- lfq_brain_proteomics ----------------
    Stratum(
        key="mayo_lfq_tcx", layer="lfq_brain_proteomics",
        study="MayoRNAseq", grant="AMP AD", cohort="Mayo Clinic Brain Bank",
        matrix_path="syn7431760/syn7431988/Mayo_Proteomics_TC_proteinoutput.txt.gz",
        reader="maxquant",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="temporal cortex",
        value_scale="log2_lfq_intensity", value_unit="log2 intensity",
        transform_chain="source:MaxQuant LFQ intensity|log2",
        visit_name="postmortem",
        assay_type="LFQ LC-MS/MS", platform=None, platform_evidence=None,
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline="MaxQuant (version not recorded)", pipeline_evidence="A",
        annot_path="MayoRNAseq_biospecimen_metadata.csv", annot_reader="adkp_biospecimen",
        person_resolvable=True,
        extra={"assay_filter": "label free mass spectrometry",
               "assay_meta": _META_MAYO_ASSAY,
               "assay_meta_key": ("specimenID", "specimen_native_id")},
    ),
    Stratum(
        key="msbb_lfq_pfc", layer="lfq_brain_proteomics",
        study="MSBB", grant="AMP AD", cohort="Mount Sinai Brain Bank",
        matrix_path="syn20801227/syn6100410/MSSM_Proteomics_PFC_PROTEINOUTPUT.txt.gz",
        reader="maxquant",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="prefrontal cortex",
        value_scale="log2_lfq_intensity", value_unit="log2 intensity",
        transform_chain="source:MaxQuant LFQ intensity|log2",
        visit_name="postmortem",
        assay_type="LFQ LC-MS/MS", platform=None, platform_evidence=None,
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline="MaxQuant (version not recorded)", pipeline_evidence="A",
        annot_path="MSBB_biospecimen_metadata.csv", annot_reader="adkp_biospecimen",
        person_resolvable=True,
        extra={"assay_filter": "label free mass spectrometry",
               "assay_meta": _META_MSBB_LFQ_ASSAY,
               "assay_meta_key": ("specimenID", "specimen_native_id")},
    ),

    # ---------------- srm_brain_proteomics ----------------
    Stratum(
        key="rosmap_srm_panel1", layer="srm_brain_proteomics",
        study="ROSMAP", grant="AMP AD", cohort="ROS/MAP",
        matrix_path="syn10468856/syn21448467/syn16779052/final_srm_data.txt.gz",
        reader="srm",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="dorsolateral prefrontal cortex",
        value_scale="log2_ratio", value_unit="log2 ratio",
        transform_chain="source:peak-area log ratio, centred",
        visit_name="postmortem",
        assay_type="LC-SRM",
        platform="triple quadrupole (instrument alias 'Smeagol')", platform_evidence="partial",
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline=None, pipeline_evidence=None,
        annot_path="syn10468856/syn21448467/syn16779052/final_srm_data_features.txt.gz",
        annot_reader="srm_features",
        extra={"assay_meta": _META_SRM_ASSAY},
    ),
    Stratum(
        key="rosmap_srm_panel2", layer="srm_brain_proteomics",
        study="ROSMAP", grant="AMP AD", cohort="ROS/MAP",
        matrix_path="syn10468856/syn21448467/syn17057851/srm_data_1226.txt.gz",
        reader="srm",
        tissue="brain", biosample_type="postmortem tissue",
        anatomic_site="dorsolateral prefrontal cortex",
        value_scale="log2_ratio", value_unit="log2 ratio",
        transform_chain="source:peak-area log ratio, centred",
        visit_name="postmortem",
        assay_type="LC-SRM",
        platform="triple quadrupole (instrument alias 'Smeagol')", platform_evidence="partial",
        suspension_type=BULK, analyte_type="protein",
        analysis_pipeline=None, pipeline_evidence=None,
        annot_path="syn10468856/syn21448467/syn17057851/final_srm_data_species.txt.gz",
        annot_reader="srm_features",
        extra={"assay_meta": _META_SRM_ASSAY},
    ),

    # ---------------- soma_plasma_proteomics ----------------
    Stratum(
        key="rosmap_soma_plasma", layer="soma_plasma_proteomics",
        study="ROSMAP", grant="AMP AD", cohort="ROS/MAP",
        matrix_path="syn64957327/OhNM2025_ROSMAP_plasma_Soma7k_protein_level_ANML_log10.csv.gz",
        reader="somascan",
        tissue="plasma", biosample_type="fluid", anatomic_site="blood",
        value_scale="log10_rfu_anml", value_unit="log10 RFU",
        transform_chain="source:SomaLogic ANML-normalised RFU|log10",
        visit_name="follow-up",
        assay_type="SomaScan v4.1 7k", platform="SOMAmer aptamer array", platform_evidence="B",
        suspension_type=NA, analyte_type="protein",
        analysis_pipeline="SomaLogic ANML normalisation", pipeline_evidence="B",
        annot_path="syn64957327/OhNM2025_ROSMAP_plasma_Soma7k_sample_metadata.csv.gz",
        annot_reader="somascan_samples",
        array_type="SOMAmer v4.1 (7k)",
        extra={"feature_meta": _META_SOMA_FEATURES},
    ),

    # ---------------- olink layers ----------------
    Stratum(
        key="pdrd_olink_plasma", layer="olink_plasma_proteomics",
        study="PPMI + PDBP", grant="AMP PDRD", cohort="PPMI/PDBP",
        matrix_path=f"{_PDRD_P}/olink-explore-format",
        reader="olink",
        tissue="plasma", biosample_type="fluid", anatomic_site="blood",
        value_scale="npx_log2", value_unit="NPX (log2)",
        transform_chain="source:Olink NPX with D01<->D02 bridging",
        visit_name="baseline/follow-up",
        assay_type="Olink Explore 1536", platform="PEA with NGS readout", platform_evidence="B",
        suspension_type=NA, analyte_type="protein",
        analysis_pipeline="Olink NPX + Olink Analyze bridging", pipeline_evidence="B",
        annot_path="PDRD/proteomics-PLA-PPEA-D03/proteomics_PLA-PPEA-D03_samples.csv.gz",
        annot_reader="olink_samples",
        array_type="Explore 1536 (4 panels)",
        # WP-11 / D-D: plain, not `_retracted`. With the two disputed plates withdrawn
        # the two variants carry identical identifiers and NPX for every remaining
        # sample, and plain is the release AMP-PD has not superseded.
        extra={"variant": "", "matrix_code": "PLA-PPEA-D03"},
    ),
    Stratum(
        key="pdrd_olink_csf", layer="olink_csf_proteomics",
        study="PPMI + PDBP", grant="AMP PDRD", cohort="PPMI/PDBP",
        matrix_path=f"{_PDRD_C}/olink-explore-format",
        reader="olink",
        tissue="cerebrospinal fluid", biosample_type="fluid",
        anatomic_site="cerebrospinal fluid",
        value_scale="npx_log2", value_unit="NPX (log2)",
        transform_chain="source:Olink NPX with D01<->D02 bridging",
        visit_name="baseline/follow-up",
        assay_type="Olink Explore 1536", platform="PEA with NGS readout", platform_evidence="B",
        suspension_type=NA, analyte_type="protein",
        analysis_pipeline="Olink NPX + Olink Analyze bridging", pipeline_evidence="B",
        annot_path="PDRD/proteomics-CSF-PPEA-D03/proteomics_CSF-PPEA-D03_samples.csv.gz",
        annot_reader="olink_samples",
        array_type="Explore 1536 (4 panels)",
        extra={"variant": "", "matrix_code": "CSF-PPEA-D03"},
    ),

    # ---------------- dia layers (WP-10) ----------------
    # Plasma is ONE stratum, not the two the plan anticipated. The delivered
    # batch-corrected product carries the plain `-PLA-` token on every sample; the
    # native/depleted (PLAN/PLAD) split exists only in the SDRF run manifests, which are
    # out of scope for this build. Splitting a stratum on a distinction the data does not
    # carry would be an assertion, not a harmonisation.
    #
    # `platform` is declared here at evidence B rather than read per sample at evidence
    # A: the SDRFs that state the instrument are not read, so this is transcription from
    # the deposit README and is graded as such. For plasma the instrument is genuinely
    # mixed across the two fractions the product merges, and that is recorded rather than
    # resolved.
    Stratum(
        key="pdrd_dia_csf", layer="dia_csf_proteomics",
        study="PPMI + PDBP", grant="AMP PDRD", cohort="PPMI/PDBP",
        matrix_path=f"{_PDRD_DIA_C}/CSF-PDIA_protein_batch-corrected.csv",
        reader="dia",
        tissue="cerebrospinal fluid", biosample_type="fluid",
        anatomic_site="cerebrospinal fluid",
        value_scale="log2_intensity_batch_corrected", value_unit="log2 intensity",
        transform_chain="source:DIA protein abundance, batch-corrected upstream|log2",
        visit_name="baseline/follow-up",
        assay_type="DIA LC-MS/MS",
        platform="Orbitrap Exploris 480", platform_evidence="B",
        suspension_type=NA, analyte_type="protein",
        analysis_pipeline="fragment -> peptide -> protein roll-up; batch correction "
                          "applied upstream by AMP-PD",
        pipeline_evidence="B",
        extra={"matrix_code": "CSF-PDIA",
               "matrix_sibling": f"{_PDRD_DIA_C}/CSF-PDIA_protein_batch-corrected_matrix.csv"},
    ),
    Stratum(
        key="pdrd_dia_plasma", layer="dia_plasma_proteomics",
        study="PPMI + PDBP", grant="AMP PDRD", cohort="PPMI/PDBP",
        matrix_path=f"{_PDRD_DIA_P}/PLA-PDIA_protein_batch-corrected.csv",
        reader="dia",
        tissue="plasma", biosample_type="fluid", anatomic_site="blood",
        value_scale="log2_intensity_batch_corrected", value_unit="log2 intensity",
        transform_chain="source:DIA protein abundance, batch-corrected upstream|log2",
        visit_name="baseline/follow-up",
        assay_type="DIA LC-MS/MS",
        platform="Orbitrap Exploris 480 and TripleTOF 6600 (fractions merged upstream)",
        platform_evidence="C",
        suspension_type=NA, analyte_type="protein",
        analysis_pipeline="fragment -> peptide -> protein roll-up; batch correction "
                          "applied upstream by AMP-PD",
        pipeline_evidence="B",
        extra={"matrix_code": "PLA-PDIA",
               "matrix_sibling": f"{_PDRD_DIA_P}/PLA-PDIA_protein_batch-corrected_matrix.csv",
               "fraction_note": "the delivered product merges native and depleted plasma "
                                "under one `-PLA-` sample token; the fraction a given "
                                "sample came from is not recoverable from it"},
    ),
]

# Layers declared but not populated in this pass. Keys are created so consumers can
# discover them and so a later pass slots in without a rebuild.
#
# Both DIA layers were here until 2026-09-03, blocked on the quant matrices. The
# batch-corrected protein-level data has since landed and both are built, so the mapping
# is empty -- kept rather than deleted because the mechanism is still the right one for
# the next layer that arrives ahead of its data.
DEFERRED_LAYERS: dict = {}

DEFERRED_STRATA: dict = {}

LAYER_ORDER = [
    "tmt_brain_proteomics",
    "lfq_brain_proteomics",
    "srm_brain_proteomics",
    "soma_plasma_proteomics",
    "olink_plasma_proteomics",
    "olink_csf_proteomics",
    "dia_plasma_proteomics",
    "dia_csf_proteomics",
]

OLINK_PANELS = ("cardiometabolic", "inflammation", "neurology", "oncology")

# UniProt <-> gene symbol reference, found in-tree (90,412 rows). Resolves the SRM
# peptide features, which are the only source keyed by gene symbol alone.
UNIPROT2SYMBOL = (f"{_R1}/syn21266445/Uniprot2Symbol_Human_v03_SimpleforR.csv.gz")

# SRM features the reference cannot resolve: beta-amyloid species, which are APP products.
SRM_MANUAL_UNIPROT = {"bA": "P05067", "bA38": "P05067"}

# syn3191087 - ROSMAP projid <-> individualID bridge. Read for IDs only (see scope gate).
ROSMAP_CROSSWALK = "ROSMAP_clinical.csv"


def strata_for(layer: str) -> list[Stratum]:
    return [s for s in STRATA if s.layer == layer]
