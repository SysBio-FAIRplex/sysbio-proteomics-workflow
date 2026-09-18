"""S3 - per-source readers.

Every reader returns a Slab: a features x samples float32 matrix with source-native
identifiers on both axes, plus whatever per-sample and per-feature annotation the
source carries. Nothing is renamed or rescaled here beyond the monotone transform
declared in the stratum's transform_chain.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import config as C


@dataclass
class Slab:
    values: np.ndarray                      # (n_features, n_samples) float32
    feature_source: list[str]
    sample_source: list[str]
    sample_annot: pd.DataFrame = field(default_factory=pd.DataFrame)
    feature_annot: pd.DataFrame = field(default_factory=pd.DataFrame)
    lod: np.ndarray | None = None           # (n_features, n_samples) float32
    notes: dict = field(default_factory=dict)

    def __post_init__(self):
        nf, ns = self.values.shape
        assert nf == len(self.feature_source), (nf, len(self.feature_source))
        assert ns == len(self.sample_source), (ns, len(self.sample_source))


def _open(path):
    p = C.ROOT / path
    return gzip.open(p, "rb") if p.name.endswith(".gz") else p.open("rb")


def _log2(vals: np.ndarray) -> np.ndarray:
    """log2 of strictly-positive entries; everything else (0, negative, NaN) -> NaN."""
    out = np.full(vals.shape, np.nan, dtype=np.float32)
    pos = np.isfinite(vals) & (vals > 0)
    out[pos] = np.log2(vals[pos])
    return out


def _read_table(path, **kw) -> pd.DataFrame:
    with _open(path) as fh:
        return pd.read_csv(io.BytesIO(fh.read()), **kw)


def _read_excel(path, **kw) -> pd.DataFrame:
    with _open(path) as fh:
        return pd.read_excel(io.BytesIO(fh.read()), **kw)


# --------------------------------------------------------------------------------------
# Matrix readers
# --------------------------------------------------------------------------------------

def read_csv_matrix(s: C.Stratum) -> Slab:
    """Plain features x samples CSV: DiverseCohorts and ROSMAP TMT R1."""
    df = _read_table(s.matrix_path, index_col=0, low_memory=False)
    df.index = df.index.astype(str)
    return Slab(
        values=df.to_numpy(dtype=np.float32),
        feature_source=df.index.tolist(),
        sample_source=[str(c) for c in df.columns],
    )


def read_msbb_tmt_xlsx(s: C.Stratum) -> Slab:
    """MSBB 19-batch xlsx: publication header, sample IDs on row 4, data from row 5.

    Values are summed TMT reporter ion intensities (linear) -> log2 here.
    """
    raw = _read_excel(s.matrix_path, header=None)
    sample_ids = [str(v) for v in raw.iloc[4].tolist()[4:] if pd.notna(v)]
    body = raw.iloc[5:].reset_index(drop=True)
    gene = body.iloc[:, 1].astype(str)
    accession = body.iloc[:, 2].astype(str)
    vals = body.iloc[:, 4:4 + len(sample_ids)].apply(
        pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)

    keep = ~(gene.isin(["nan", "None"]) & accession.isin(["nan", "None"]))
    vals, gene, accession = vals[keep.values], gene[keep].values, accession[keep].values
    vals = _log2(vals)

    return Slab(
        values=vals,
        feature_source=[f"{g}|{a}" for g, a in zip(gene, accession)],
        sample_source=sample_ids,
        feature_annot=pd.DataFrame({"gene_raw": gene, "accession_raw": accession}),
        notes={"transform": "log2 of linear reporter intensity"},
    )


def read_maxquant(s: C.Stratum) -> Slab:
    """MaxQuant proteinGroups.txt: take `LFQ intensity <sample>` columns, log2."""
    df = _read_table(s.matrix_path, sep="\t", low_memory=False)
    lfq = [c for c in df.columns if c.startswith("LFQ intensity ")]
    samples = [c[len("LFQ intensity "):] for c in lfq]

    # MaxQuant decoy / contaminant rows must go before anything else.
    mask = pd.Series(True, index=df.index)
    for col in ("Reverse", "Potential contaminant", "Only identified by site"):
        if col in df.columns:
            mask &= df[col].isna() | (df[col].astype(str).str.strip() == "")
    df = df[mask]

    vals = _log2(df[lfq].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32))

    groups = df["Majority protein IDs"].fillna(df["Protein IDs"]).astype(str)
    is_pool = [bool(re.search(r"gis|pool", x, re.I)) for x in samples]

    return Slab(
        values=vals,
        feature_source=groups.tolist(),
        sample_source=samples,
        sample_annot=pd.DataFrame({
            "specimen_source_value": samples,
            "is_pool": is_pool,
            "batch": [x.split("_")[0] if "_" in x else None for x in samples],
        }),
        feature_annot=pd.DataFrame({"protein_group": groups.values}),
        notes={"n_decoy_contaminant_dropped": int((~mask).sum()),
               "transform": "log2 of MaxQuant LFQ intensity"},
    )


_R2_RATIO = re.compile(r"^Abundance Ratio: \(F(\d+), ([0-9]+[NC]?)\) / \(F\1, 126\)$")


def read_pd_ratio(s: C.Stratum) -> Slab:
    """Proteome Discoverer TMT ratio-to-reference export (ROSMAP round 2).

    Standard TMT practice for a design with an internal reference channel:
      1. take Abundance Ratio (channel / 126 GIS) -- already reference-normalised
      2. log2
      3. centre each protein on the median of its per-plex medians

    Step 3 mirrors how the round-1 product was constructed
    ("log2 abundanceRatio centered on median of batch medians per protein"), so the two
    ROSMAP TMT rounds end up on a comparable footing rather than an ad hoc one.
    """
    df = _read_table(s.matrix_path, sep="\t", low_memory=False)
    ratio_cols, plexes, channels = [], [], []
    for c in df.columns:
        m = _R2_RATIO.match(c)
        if m:
            ratio_cols.append(c)
            plexes.append(int(m.group(1)))
            channels.append(m.group(2))

    acc = df["Accession"].astype(str)
    gene = df["Gene Symbol"].astype(str) if "Gene Symbol" in df.columns else None
    keep = acc.notna() & (acc != "nan")
    df, acc = df[keep], acc[keep]
    if gene is not None:
        gene = gene[keep]

    vals = _log2(df[ratio_cols].apply(pd.to_numeric, errors="coerce")
                 .to_numpy(dtype=np.float32))

    # centre each protein on the median of its per-plex medians
    plex_arr = np.asarray(plexes)
    pres = np.isfinite(vals)
    per_plex, plex_present = [], []
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN rows are expected
        for p in sorted(set(plexes)):
            sel = plex_arr == p
            per_plex.append(np.nanmedian(vals[:, sel], axis=1))
            plex_present.append(pres[:, sel].any(axis=1))
        centre = np.nanmedian(np.vstack(per_plex), axis=0)
    vals = (vals - centre[:, None]).astype(np.float32)

    # Missingness here is per-PLEX identification, not per-sample detection: a protein
    # identified in one plex of the multiconsensus but not another produces an
    # all-or-nothing block. Measured so the README can state it rather than assume it.
    pp = np.vstack(plex_present)                       # (n_plex, n_features)
    n_plex = pp.shape[0]
    plex_hits = pp.sum(axis=0)
    partial = 0
    for i, p in enumerate(sorted(set(plexes))):
        sel = plex_arr == p
        blk = pres[:, sel]
        partial += int((blk.any(axis=1) & ~blk.all(axis=1)).sum())
    structured = {
        "n_plexes": n_plex,
        "mean_plexes_present_per_feature": round(float(plex_hits.mean()), 2),
        "features_partially_present_within_a_plex": partial,
        "missingness_is_plex_structured": bool(partial == 0),
        "n_features_entirely_empty": int((~pres.any(axis=1)).sum()),
    }

    sample_ids = [f"F{p}.{ch}" for p, ch in zip(plexes, channels)]
    return Slab(
        values=vals,
        feature_source=([f"{g}|{a}" for g, a in zip(gene, acc)] if gene is not None
                        else acc.tolist()),
        sample_source=sample_ids,
        sample_annot=pd.DataFrame({
            "specimen_source_value": sample_ids,
            "batch": [f"plex{p}" for p in plexes],
            "channel": channels,
            "is_pool": False,
        }),
        notes={"transform": "log2 of PD abundance ratio vs 126 reference channel",
               "centring": "per-protein median of per-plex medians",
               "n_plexes": len(set(plexes)), "n_channels": len(ratio_cols),
               "structured_missingness": structured},
    )


def read_srm(s: C.Stratum) -> Slab:
    """LC-SRM: samples are rows, features are columns. Transposed to the convention."""
    df = _read_table(s.matrix_path, sep="\t", low_memory=False)
    meta_cols = [c for c in df.columns
                 if c in {"sample_id", "sample_num", "plate", "plate_row", "plate_col",
                          "plate_well", "subject.id", "isControl", "sample.type",
                          "Replicate.Name"} or c.startswith("Unnamed")]
    feat_cols = [c for c in df.columns if c not in meta_cols]

    sample_key = "sample_id" if "sample_id" in df.columns else df.columns[0]
    samples = df[sample_key].astype(str).tolist()

    is_control = (df["isControl"].astype(str).str.upper() == "TRUE"
                  if "isControl" in df.columns else pd.Series(False, index=df.index))
    stype = df["sample.type"].astype(str) if "sample.type" in df.columns else None
    pool = is_control.values | (stype.ne("subjects").values if stype is not None
                                else np.zeros(len(df), bool))

    # subject.id carries non-numeric sentinels ('control', ...) on reference wells.
    # Kept as the source string; only the sentinels are nulled.
    if "subject.id" in df.columns:
        raw = df["subject.id"].astype(str).str.strip()
        numeric_like = raw.str.fullmatch(r"\d+(\.0+)?")
        person = raw.where(numeric_like, None).str.replace(r"\.0+$", "", regex=True)
    else:
        person = None

    annot = pd.DataFrame({
        "specimen_source_value": samples,
        "person_source_value": person,
        "is_pool": pool,
        "plate_id": df["plate"].astype(str) if "plate" in df.columns else None,
    })

    vals = df[feat_cols].apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float32).T
    # WP-12 item 6. `SOURCE_QC` declares control-well flagging as ours to apply; this
    # records that it actually happened and on what evidence, so the claim is measurable
    # instead of merely declared.
    return Slab(values=vals, feature_source=list(feat_cols), sample_source=samples,
                sample_annot=annot,
                notes={"control_wells_flagged": {
                    "n_flagged": int(pool.sum()),
                    "n_samples": len(df),
                    "from_isControl": int(is_control.sum()),
                    "from_sample_type_not_subjects": (
                        int(stype.ne("subjects").sum()) if stype is not None else 0),
                    "columns_used": [c for c in ("isControl", "sample.type")
                                     if c in df.columns]}})


def read_somascan(s: C.Stratum) -> Slab:
    """SomaScan: samples are rows keyed by projid_visit, SeqIds are columns."""
    df = _read_table(s.matrix_path, index_col=0, low_memory=False)
    samples = [str(i) for i in df.index]
    seq_ids = [str(c) for c in df.columns]

    meta = _read_table(s.extra["feature_meta"], low_memory=False)
    meta["SeqId"] = meta["SeqId"].astype(str)
    meta = meta.drop_duplicates("SeqId").set_index("SeqId")
    aligned = meta.reindex(seq_ids)

    feature_annot = pd.DataFrame({
        "seq_id": seq_ids,
        "uniprot_raw": aligned["UniProt"].values,
        "gene_raw": aligned["EntrezGeneSymbol"].values,
        "dilution": aligned["Dilution"].astype(str).values,
        "organism": aligned["Organism"].values,
        "somamer_type": aligned["Type"].values,
        # WP-12 item 2. `ColCheck` is SomaLogic's per-analyte column QC and `SOURCE_QC`
        # has always declared it OURS to apply -- it was read into this frame and then
        # dropped, so the declaration was never true. Carried as a flag beside
        # `frac_below_lod` rather than used to delete analytes: the standing policy is
        # flag, not drop (D8), and SomaLogic leaves the threshold to the study.
        "colcheck": (aligned["ColCheck"].astype(str).str.strip().values
                     if "ColCheck" in aligned.columns else None),
    })
    keep = (feature_annot["organism"].astype(str).str.strip() == "Human") & \
           (feature_annot["somamer_type"].astype(str).str.strip() == "Protein")

    cc = feature_annot.loc[keep, "colcheck"]
    vals = df.to_numpy(dtype=np.float32).T[keep.values]
    return Slab(
        values=vals,
        feature_source=[seq_ids[i] for i in np.flatnonzero(keep.values)],
        sample_source=samples,
        feature_annot=feature_annot[keep].reset_index(drop=True),
        notes={"n_non_human_or_control_dropped": int((~keep).sum()),
               "colcheck_counts": ({str(k): int(v) for k, v in
                                    cc.value_counts(dropna=False).items()}
                                   if cc.notna().any() else None),
               # Recorded as a measured absence, not left to look like a pass: the
               # deposit's sample metadata is clinical only (projid, Visit, msex, age,
               # Diagnosis, cogn, apoe, educ) and carries no QC column at all.
               "sample_level_qc_in_deposit": None},
    )


def read_olink(s: C.Stratum) -> Slab:
    """Olink Explore: read the long `olink_explore_format` files (they carry LOD and QC),
    tag each feature with its panel BEFORE concatenating, then pivot.

    WP-11 / D-D. The build reads the PLAIN files and withdraws every sample on the two
    disputed plates. Two measured facts make that clean rather than a judgement call:
    the sample-id set on those plates is identical in the plain and `_retracted`
    variants, so which rows to drop does not depend on resolving the swap; and off those
    plates the two variants are identical row-for-row, so once the disputed samples are
    gone the variant choice is moot for everything that ships. Plain is the release
    AMP-PD has not superseded.

    The QC-recovery merge that used to run here has been deleted. It joined the plain
    files' QC onto `_retracted` rows on (sample_id, UniProt, panel); because the swap
    permuted sample ids between the two plates, that handed each disputed sample the
    verdict belonging to the other plate's well -- 3,006 misattributed verdicts, all of
    them on the samples now withdrawn, and a no-op for every other sample (N4).
    """
    code = s.extra["matrix_code"]
    variant = s.extra.get("variant", "")
    frames = []

    for panel in C.OLINK_PANELS:
        path = f"{s.matrix_path}/{code}_olink_explore_format_{panel}{variant}.csv.gz"
        d = _read_table(path, low_memory=False)
        d["panel"] = panel
        frames.append(d)

    long = pd.concat(frames, ignore_index=True)
    long["UniProt"] = long["UniProt"].astype(str)
    long["sample_id"] = long["sample_id"].astype(str)

    withdrawn: dict = {}
    if "PlateID" in long.columns:
        on_plate = long["PlateID"].astype(str).isin(C.OLINK_WITHDRAWN_PLATES)
        if on_plate.any():
            ids = sorted(long.loc[on_plate, "sample_id"].unique())
            donors = sorted(long.loc[on_plate, "participant_id"].astype(str).unique())
            kept_donors = set(long.loc[~on_plate, "participant_id"].astype(str))
            withdrawn = {
                "rule": "PlateID in " + ", ".join(sorted(C.OLINK_WITHDRAWN_PLATES)),
                "reason": "donor identity in dispute pending AMP-PD (plate swap, D-D)",
                "n_samples": len(ids),
                "n_samples_total": int(long["sample_id"].nunique()),
                "n_participants_losing_every_sample":
                    len([d for d in donors if d not in kept_donors]),
            }
            long = long.loc[~on_plate].reset_index(drop=True)

    npx = long.pivot_table(index="UniProt", columns="sample_id", values="NPX",
                           aggfunc="mean")
    lod = long.pivot_table(index="UniProt", columns="sample_id", values="LOD",
                           aggfunc="mean").reindex(index=npx.index, columns=npx.columns)

    # WP-12 item 4: the vendor's own per-assay missingness, carried rather than
    # discarded, so it can be cross-checked against the `frac_below_lod` we compute
    # ourselves. Reported side by side; neither is picked over the other.
    agg = {"panel": ("panel", "first"), "panel_lot_nr": ("Panel_Lot_Nr", "first")}
    if "MissingFreq" in long.columns:
        agg["missing_freq"] = ("MissingFreq", "max")
    feat = (long.groupby("UniProt").agg(**agg).reindex(npx.index).reset_index())

    def worst(g):
        vals = set(str(x).upper() for x in g if pd.notna(x))
        return "WARN" if vals - {"PASS"} else "PASS"

    samp = (long.groupby("sample_id")
            .agg(participant_id=("participant_id", "first"),
                 visit_month=("visit_month", "first"),
                 plate_id=("PlateID", "first"),
                 qc_warning=("QC_Warning", worst),
                 distribution_qc=("Distribution_QC", worst),
                 outliers_qc=("Outliers_QC", worst))
            .reindex(npx.columns).reset_index())

    qc_status = np.where(
        (samp.qc_warning.eq("PASS") & samp.distribution_qc.eq("PASS")
         & samp.outliers_qc.eq("PASS")), "PASS", "WARN")

    # WP-12 item 5. `qc_status` stays the one-word summary, but folding three independent
    # vendor checks into it means a consumer inherits our threshold instead of setting
    # their own -- a sample that trips only `Outliers_QC` is not the same as one that
    # trips all three. The three verdicts and their count ship alongside.
    checks = ["qc_warning", "distribution_qc", "outliers_qc"]
    warn_count = sum(samp[c].ne("PASS").astype(int) for c in checks)
    qc_detail = (samp.qc_warning.astype(str) + "|" + samp.distribution_qc.astype(str)
                 + "|" + samp.outliers_qc.astype(str))

    annot = pd.DataFrame({
        "specimen_source_value": samp.sample_id.astype(str),
        "person_source_value": samp.participant_id.astype(str),
        "visit_month": pd.to_numeric(samp.visit_month, errors="coerce").astype("Int32"),
        "plate_id": samp.plate_id.astype(str),
        "qc_status": qc_status,
        "qc_warn_count": warn_count.astype("int32"),
        "qc_detail": "QC_Warning|Distribution_QC|Outliers_QC=" + qc_detail,
        "is_pool": False,
    })

    return Slab(
        values=npx.to_numpy(dtype=np.float32),
        feature_source=[str(i) for i in npx.index],
        sample_source=[str(c) for c in npx.columns],
        sample_annot=annot,
        feature_annot=feat.rename(columns={"UniProt": "uniprot_raw"}),
        lod=lod.to_numpy(dtype=np.float32),
        notes={"variant_used": variant or "plain",
               "n_panels": len(C.OLINK_PANELS),
               "withdrawn_samples": withdrawn or None,
               "qc_verdicts_kept_separately": checks,
               "n_samples_with_any_qc_warn": int((warn_count > 0).sum()),
               "n_samples_by_warn_count": {str(k): int(v) for k, v in
                                           warn_count.value_counts().sort_index().items()},
               "vendor_missing_freq_present": "MissingFreq" in long.columns},
    )


def read_dia(s: C.Stratum) -> Slab:
    """PDRD DIA, batch-corrected protein level (WP-10).

    Reads the LONG form, never the `*_matrix.csv` sibling. The two disagree, and the
    disagreement is not cosmetic: the matrix is dense because it encodes non-detection as
    a literal `0`, while the long form simply omits those rows. Proven on the delivered
    files -- the long forms contain no exact zero among 959,896 observed values, the
    matrix's zero set is cell-for-cell identical to the set the long form omits (36,002
    cells in plasma, 12.1%; 187 in CSF), and observed intensities stop dead at 16.70
    (plasma) and 30.21 (CSF) with nothing in between. Reading the matrix would inject
    36k false zeros that the log2 below turns into NaN or, worse, a consumer reads as a
    real measurement of nothing.

    Sample ids follow the same `<cohort>-<participant>-<visit>-<tissue>-<assay>` grammar
    the Olink reader parses. `visit_month` is carried natively and agrees with the visit
    token in every delivered sample, so this layer lands on the cross-layer time axis at
    evidence `source` rather than derived (WP-14).
    """
    long = _read_table(s.matrix_path, low_memory=False)
    long["UniProt"] = long["UniProt"].astype(str)
    long["sample_id"] = long["sample_id"].astype(str)

    n_rows_in = len(long)
    zeros = int((long["protein_abundance"] == 0).sum())
    negatives = int((long["protein_abundance"] < 0).sum())

    wide = long.pivot_table(index="UniProt", columns="sample_id",
                            values="protein_abundance", aggfunc="mean")

    # Sample-level QC. The deposit ships none of its own -- there is no QC column in the
    # delivery at all -- so completeness is the only sample-level quality signal
    # available, and it separates cleanly (see DIA_MIN_SAMPLE_COMPLETENESS).
    n_feat = wide.shape[0]
    detected = wide.notna().sum(axis=0)
    completeness = detected / n_feat
    keep = completeness >= C.DIA_MIN_SAMPLE_COMPLETENESS
    dropped = sorted(completeness.index[~keep])
    # Detected-protein counts only, never the sample ids. A DIA sample id embeds the
    # AMP-PD participant (`PP-<participant>-<visit>-PLA-PDIA`), so listing the dropped
    # samples publishes participant identifiers into `build_report.json`, the audit TSV
    # and the generated README -- all in the publishable set (SEC-7/SEC-8). The counts
    # carry the whole diagnostic point: they show the failures are unambiguous.
    dropped_counts = sorted(int(detected[sid]) for sid in dropped)
    wide = wide.loc[:, keep.values]

    # A feature that survives in no retained sample carries no information.
    empty_feats = sorted(wide.index[wide.notna().sum(axis=1) == 0])
    if empty_feats:
        wide = wide.drop(index=empty_feats)

    meta = (long.drop_duplicates("sample_id").set_index("sample_id")
            .reindex(wide.columns))
    detected_kept = wide.notna().sum(axis=0)

    annot = pd.DataFrame({
        "specimen_source_value": [str(c) for c in wide.columns],
        "person_source_value": meta["participant_id"].astype(str).values,
        "visit_month": pd.to_numeric(meta["visit_month"], errors="coerce")
                         .astype("Int32").values,
        "visit_month_evidence": "source",
        "qc_status": "PASS",
        "is_pool": False,
        "n_detected": detected_kept.astype("int32").values,
        "detection_rate": (detected_kept / len(wide)).astype("float32").values,
    })

    return Slab(
        # log2 so the layer sits on the same kind of axis as the other intensity layers;
        # the source is a linear batch-corrected abundance.
        values=_log2(wide.to_numpy(dtype=np.float32)),
        feature_source=[str(i) for i in wide.index],
        sample_source=[str(c) for c in wide.columns],
        sample_annot=annot,
        feature_annot=pd.DataFrame({"uniprot_raw": [str(i) for i in wide.index]}),
        notes={
            "source_form": "long (one row per observed sample x protein)",
            "matrix_sibling_not_read": (
                "the *_matrix.csv sibling encodes non-detection as literal 0; reading it "
                "would inject false zeros. Values agree exactly where both are present."),
            "n_long_rows": n_rows_in,
            "n_source_zeros": zeros,
            "n_source_negatives": negatives,
            "qc_min_sample_completeness": C.DIA_MIN_SAMPLE_COMPLETENESS,
            "n_samples_dropped_low_completeness": len(dropped),
            "detected_counts_of_dropped_samples": dropped_counts,
            "n_features_total": int(n_feat),
            "n_features_dropped_empty": len(empty_feats),
            "features_dropped_empty": empty_feats,
        },
    )


MATRIX_READERS = {
    "csv_matrix": read_csv_matrix,
    "msbb_tmt_xlsx": read_msbb_tmt_xlsx,
    "maxquant": read_maxquant,
    "srm": read_srm,
    "pd_ratio": read_pd_ratio,
    "somascan": read_somascan,
    "olink": read_olink,
    "dia": read_dia,
}


# --------------------------------------------------------------------------------------
# Annotation readers -> DataFrame with specimen_source_value + person_source_value
# --------------------------------------------------------------------------------------

def annot_dc_traits(s: C.Stratum) -> pd.DataFrame:
    """DiverseCohorts traits sheet.

    Note on is_pool: the delivered residual matrices are already post-GIS-removal
    (QC_process_*.txt: 1191 -> 1105 non-GIS -> 1086 after PCA outlier removal), so no
    pool channels survive into what we read. A missing individualID here therefore means
    *unresolved donor*, not a reference pool -- the `sd` cohort has 400 rows and zero
    individualID values. Conflating the two would mislabel 371 real samples.
    """
    d = _read_excel(s.annot_path, dtype=str)
    ind = "individualID" if "individualID" in d.columns else "IndividualID"

    def _clean(col):
        if col not in d.columns:
            return pd.Series([None] * len(d), index=d.index, dtype=object)
        return d[col].astype(str).str.strip().replace(
            {"nan": None, "NA": None, "None": None, "": None, "GIS": None})

    individual = _clean(ind)
    proj = _clean("projID")

    # The `sd` contributing dataset carries NO individualID at all -- it is identified by
    # ROSMAP projID instead. Falling back to projID recovers 371 samples that would
    # otherwise have no donor. The two columns are different identifier systems, so the
    # namespace is recorded per row rather than assumed for the stratum.
    person = individual.where(individual.notna(), proj)
    namespace = np.where(individual.notna(), "AMP-AD.individualID",
                         np.where(proj.notna(), "ROSMAP.projid", None))
    is_gis = person.isna()

    return pd.DataFrame({
        "specimen_source_value": d["ID"].astype(str),
        "person_source_value": person,
        "person_source_namespace_row": namespace,
        "batch": (d["Batch"].astype(str) if "Batch" in d.columns
                  else d["ID"].astype(str).str.extract(r"(b\d+)")[0]),
        "channel": (d["Channel"].astype(str) if "Channel" in d.columns
                    else d["ID"].astype(str).str.split(".").str[-1]),
        "source_dataset": d["Dataset"].astype(str) if "Dataset" in d.columns else None,
        "is_pool": is_gis.values,
    })


def load_rosmap_crosswalk() -> pd.DataFrame:
    """ROSMAP projid <-> individualID bridge from syn3191087.

    Scope gate: this file is overwhelmingly clinical (msex, educ, race, apoe_genotype,
    age_death, braaksc, ceradsc, cogdx, pmi, mmse). ONLY the two identifier columns and
    the study label are read; nothing else is loaded into memory.
    """
    d = pd.read_csv(C.ROOT / C.ROSMAP_CROSSWALK, dtype=str,
                    usecols=["projid", "individualID", "Study"])
    d = d.dropna(subset=["projid", "individualID"])
    return pd.DataFrame({
        "rosmap_projid": d["projid"].str.strip(),
        "rosmap_individual_id": d["individualID"].str.strip(),
        "rosmap_study": d["Study"].str.strip(),
    }).drop_duplicates()


def annot_rosmap_r1_specimen(s: C.Stratum) -> pd.DataFrame:
    d = _read_table(s.annot_path, low_memory=False, dtype=str)
    out = pd.DataFrame({
        "specimen_source_value": d["SampleID"].astype(str),
        "person_source_value": d["projid"].astype(str).replace({"NA": None, "nan": None}),
        "batch": d["Batch"].astype(str),
        "channel": d["batch.channel"].astype(str).str.split(".").str[-1],
        "specimen_native_id": d["SpecimenID"].astype(str),
    })
    out["is_pool"] = (d["SpecimenID"].astype(str).str.contains("GIS", case=False, na=False)
                      | out.person_source_value.isna()).values
    return out


def annot_somascan_samples(s: C.Stratum) -> pd.DataFrame:
    # Scope gate: take identifiers only. msex/age/Diagnosis/cogn/apoe/educ are dropped here.
    d = _read_table(s.annot_path, usecols=["projid_visit", "projid", "Visit", "study"],
                    low_memory=False, dtype=str)
    visit_index = pd.to_numeric(d["Visit"], errors="coerce").astype("Int32")

    # WP-14 / D-H. `Visit` is ROSMAP's follow-up year, so `visit_month = Visit * 12` is a
    # structural conversion of an ordinal follow-up number -- it never touches
    # age_at_visit and never trips the scope gate. Without it this layer has no elapsed-
    # time axis at all and nothing longitudinal can cross cohorts: Olink and DIA are in
    # months since baseline, SomaScan had only an ordinal.
    #
    # ROSMAP follow-ups are nominally annual but not exactly, so this is an approximation
    # and is labelled as one. `visit_index` is kept beside it as the source-native value,
    # so the conversion is reversible and is never itself the join key.
    return pd.DataFrame({
        "specimen_source_value": d["projid_visit"].astype(str),
        "person_source_value": d["projid"].astype(str),
        "visit_index": visit_index,
        "visit_month": (visit_index * 12).astype("Int32"),
        "visit_month_evidence": np.where(visit_index.notna(),
                                         "derived_from_visit_index", None),
        "source_dataset": d["study"].astype(str).str.strip(),
        "is_pool": False,
    })


def annot_olink_samples(s: C.Stratum) -> pd.DataFrame:
    d = _read_table(s.annot_path, usecols=["participant_id", "sample_id", "visit_month"],
                    low_memory=False, dtype=str)
    return pd.DataFrame({
        "specimen_source_value": d["sample_id"].astype(str),
        "person_source_value": d["participant_id"].astype(str),
        "visit_month": pd.to_numeric(d["visit_month"], errors="coerce").astype("Int32"),
        "visit_month_evidence": "source",
        "is_pool": False,
    }).drop_duplicates("specimen_source_value")


def _adkp_biospecimen(path: str) -> pd.DataFrame:
    """AD Knowledge Portal biospecimen metadata, scope-gated.

    Drops samplingAge / samplingAgeUnits (demographic, plan section 0) and keeps only
    identifiers plus genuine biospecimen descriptors.
    """
    keep = ["individualID", "specimenID", "specimenIdSource", "organ", "tissue",
            "BrodmannArea", "sampleStatus", "nucleicAcidSource", "cellType",
            "isPostMortem", "assay", "exclude", "excludeReason"]
    d = pd.read_csv(C.ROOT / path, low_memory=False, dtype=str)
    return d[[c for c in keep if c in d.columns]]


def _msbb_specimen_key(x: str) -> str:
    """Matrix `b1_03_1601` -> metadata `b1_1601_03`; `b7r2_03_1752` -> `r2b7_1752_03`."""
    p = str(x).split("_")
    if len(p) != 3:
        return str(x)
    batch = re.sub(r"^b(\d+)r(\d+)$", r"r\2b\1", p[0])
    return f"{batch}_{p[2]}_{p[1]}"


def _mayo_specimen_key(x: str) -> str:
    """Matrix `mayo_b1_001_11` -> metadata `b1_001`."""
    p = str(x).replace("mayo_", "").split("_")
    return f"{p[0]}_{p[1]}" if len(p) >= 2 else str(x)


_KEY_FN = {"msbb_lfq_pfc": _msbb_specimen_key, "mayo_lfq_tcx": _mayo_specimen_key}


def annot_adkp_biospecimen(s: C.Stratum) -> pd.DataFrame:
    """Resolve run labels to individualID via the AD Knowledge Portal biospecimen file."""
    meta = _adkp_biospecimen(s.annot_path)
    if "assay" in meta.columns and s.extra.get("assay_filter"):
        meta = meta[meta.assay == s.extra["assay_filter"]]
    meta = meta.drop_duplicates("specimenID")
    idx = meta.set_index(meta.specimenID.astype(str))

    slab_samples = s.extra["_samples"]
    key_fn = _KEY_FN.get(s.key, lambda x: str(x))
    keys = [key_fn(x) for x in slab_samples]
    al = idx.reindex(keys)

    pool = pd.Series(al.individualID.isna().values) & pd.Series(
        [bool(re.search(r"gis|pool", str(x), re.I)) for x in slab_samples])

    # WP-12 item 3. `exclude` / `excludeReason` were kept in the frame above and consumed
    # by nothing, so an exclusion the portal states was silently ignored. Wired to
    # `qc_status` -- flagged, never dropped (D8), because an exclusion is a quality
    # verdict and the consumer sets their own threshold. No proteomics specimen is
    # currently flagged in either file, so this changes nothing today; MSBB's
    # `excludeReason` vocabulary includes `sample swap`, which is why an unwired flag is a
    # live risk at the next metadata release rather than a cosmetic gap.
    exc = (al["exclude"].astype(str).str.strip().str.upper().isin(["TRUE", "YES", "1"])
           if "exclude" in al.columns else pd.Series(False, index=al.index))
    reason = (al["excludeReason"].values if "excludeReason" in al.columns else None)
    out = pd.DataFrame({
        "specimen_source_value": [str(x) for x in slab_samples],
        "person_source_value": al.individualID.astype(str).where(
            al.individualID.notna(), None).values,
        "specimen_native_id": al.specimenID.values,
        "batch": [str(x).split("_")[0] for x in slab_samples],
        "anatomic_site_source": al.tissue.values if "tissue" in al else None,
        "brodmann_area": al.BrodmannArea.values if "BrodmannArea" in al else None,
        "is_pool": pool.values,
        "qc_status": np.where(exc.values, "EXCLUDED_BY_SOURCE", "not_assessed"),
        "exclude_reason": reason,
    })
    s.extra["_adkp_exclusions"] = {
        "n_flagged": int(exc.sum()),
        "n_specimens": len(out),
        "source": s.annot_path,
        "reasons": ({str(k): int(v) for k, v in
                     pd.Series(reason)[exc.values].value_counts().items()}
                    if reason is not None and exc.any() else {}),
    }
    return out


def annot_rosmap_r2_traits(s: C.Stratum) -> pd.DataFrame:
    """ROSMAP round-2 TMT labelling sheet.

    Scope gate: only Plex / Channel / ProjID / Lot are taken. cogdx_3lev, msex, pmi,
    braaksc and ceradsc are clinical and are dropped at read.
    """
    d = _read_excel(s.annot_path, dtype=str)
    chan = d["Channel"].astype(str).str.replace("_", "", regex=False).str.strip()
    plex = d["Plex"].astype(str).str.strip()
    proj = d["ProjID"].astype(str).str.strip()
    is_gis = proj.str.upper().eq("GIS") | chan.eq("126")
    return pd.DataFrame({
        "specimen_source_value": "F" + plex + "." + chan,
        "person_source_value": proj.where(~is_gis, None),
        "batch": "plex" + plex,
        "channel": chan,
        "reagent_lot": d["Lot #"].astype(str) if "Lot #" in d.columns else None,
        "is_pool": is_gis.values,
    })


# --------------------------------------------------------------------------------------
# Assay metadata (WP-8) -- platform, reagent lot and pool flags, per sample
# --------------------------------------------------------------------------------------

# Columns taken from an AD Knowledge Portal assay-metadata object. Explicit allow-list
# (SEC-2): these files are technical, but the rule is the same for every source.
_ASSAY_META_COLS = ("platform", "isAssayControl", "controlType", "lotNumber",
                    "fragmentation", "resolution", "TMTType", "TMTtype")


def attach_assay_metadata(s: C.Stratum, slab: Slab) -> None:
    """Spend the assay-metadata objects WP-1 landed (WP-8).

    Until this existed the five objects sat on disk, referenced from `config.py` and read
    by nothing, so `platform` shipped null for six strata while the audit reported
    "no evidence available on disk" -- true of the build, not of the deposit.

    What lands per sample: `platform` (evidence A, replacing a config declaration or a
    null), `reagent_lot`, and an evidence-based `is_pool` from `isAssayControl` /
    `controlType = GIS` rather than inferred from a missing donor id.

    The acquisition descriptors (`fragmentation`, `resolution`, TMT type) are NOT emitted
    as columns -- they are summarised into `notes` so the `cde_conflict` check can compare
    them against what `config.py` declares, which is where D-C/N3 surfaces.
    """
    meta_path = s.extra.get("assay_meta")
    key = s.extra.get("assay_meta_key")
    if not meta_path:
        return
    d = _read_table(meta_path, low_memory=False, dtype=str)

    if not key:
        # SRM: the object exists and is read, but carries no platform for any specimen
        # and is keyed by individualID where the matrix is keyed by plate well. Recorded
        # as a measured absence so the audit can say the source is silent rather than
        # that we did not look.
        s.extra["_assay_meta_note"] = {
            "joined": False,
            "reason": "no join key declared; matrix and metadata use different "
                      "identifier systems",
            "n_rows": len(d),
            "platform_populated_in_source": int(d["platform"].notna().sum())
            if "platform" in d.columns else 0,
        }
        slab.notes["assay_meta"] = s.extra["_assay_meta_note"]
        return

    meta_col, sample_col = key
    if meta_col not in d.columns:
        slab.notes["assay_meta"] = {"joined": False,
                                    "reason": f"{meta_col} absent from {meta_path}"}
        return

    # `assay_filter` is NOT applied here. It exists for the biospecimen files, which mix
    # assays in one table, and the two sources do not share a vocabulary: MSBB's
    # biospecimen file says `label free mass spectrometry` where its assay-metadata object
    # says `LC-MSMS`. Applying it to this file silently matched 0 of 306 rows and cost
    # `msbb_lfq_pfc` its platform entirely -- a filter that removes everything looks
    # exactly like a source with nothing in it. Each assay-metadata object is
    # single-assay, so there is nothing to filter; the values are recorded instead, so a
    # future multi-assay object is visible rather than silently over-matched.
    assays_in_file = (sorted({str(x) for x in d["assay"].dropna().unique()})
                      if "assay" in d.columns else [])

    keys = (slab.sample_annot[sample_col] if sample_col in slab.sample_annot.columns
            else pd.Series(slab.sample_source))
    keys = keys.astype(str).values

    idx = d.drop_duplicates(meta_col).set_index(d[meta_col].astype(str))
    al = idx.reindex(keys)
    n_hit = int(al["platform"].notna().sum()) if "platform" in al.columns else 0

    def take(col):
        return al[col].values if col in al.columns else None

    plat = take("platform")
    if plat is not None:
        slab.sample_annot["platform"] = plat
    lot = take("lotNumber")
    if lot is not None:
        # Only overwrite where the metadata actually has a lot; `rosmap_r2_traits`
        # already supplies one from a different source.
        cur = slab.sample_annot.get("reagent_lot")
        slab.sample_annot["reagent_lot"] = (
            pd.Series(lot).where(pd.Series(lot).notna(), cur.values if cur is not None
                                 else None).values)

    # Pool flags from evidence rather than inference. A missing donor id means several
    # different things (see annot_dc_traits), so where the deposit states control status
    # it is authoritative and replaces the inferred flag.
    ctrl, ctype = take("isAssayControl"), take("controlType")
    stated = None
    if ctrl is not None:
        stated = pd.Series(ctrl).astype(str).str.upper().eq("TRUE")
        stated = stated.where(pd.Series(ctrl).notna(), None)
    if ctype is not None:
        gis = pd.Series(ctype).astype(str).str.upper().eq("GIS")
        gis = gis.where(pd.Series(ctype).notna(), None)
        stated = gis if stated is None else stated.astype("boolean").fillna(
            gis.astype("boolean"))
    n_pool_stated = 0
    if stated is not None and stated.notna().any():
        st = stated.astype("boolean")
        cur = slab.sample_annot.get("is_pool")
        base = (pd.Series(cur).astype("boolean") if cur is not None
                else pd.Series([pd.NA] * len(keys), dtype="boolean"))
        merged = st.fillna(base).fillna(False)
        n_pool_stated = int(st.fillna(False).sum())
        slab.sample_annot["is_pool"] = merged.astype(bool).values

    def summarise(col):
        v = take(col)
        if v is None:
            return None
        c = pd.Series(v).dropna().value_counts()
        return {str(k): int(n) for k, n in c.items()} or None

    slab.notes["assay_meta"] = {
        "joined": True,
        "source": meta_path,
        "on": f"{meta_col} -> {sample_col}",
        "n_samples_matched": n_hit,
        "n_samples": len(keys),
        "assay_values_in_object": assays_in_file,
        "platform_observed": summarise("platform"),
        "reagent_lot_observed": summarise("lotNumber"),
        "n_pool_flagged_by_source": n_pool_stated,
        # Acquisition descriptors, kept for the conflict check rather than emitted.
        "acquisition_observed": {
            k: v for k, v in (("fragmentation", summarise("fragmentation")),
                              ("resolution", summarise("resolution")),
                              ("tmt_type", summarise("TMTType") or summarise("TMTtype")))
            if v},
    }


# --------------------------------------------------------------------------------------
# Biospecimen metadata -- anatomic site and Brodmann area (WP-9)
# --------------------------------------------------------------------------------------

def attach_biospecimen_metadata(s: C.Stratum, slab: Slab) -> None:
    """DiverseCohorts anatomic site and Brodmann area from syn51757645 (WP-9 step 2).

    Two joins, in precedence order:

    1. **Direct**, on `specimenID` -- the biospecimen file keys its `TMT quantitation`
       rows exactly as the matrix does (`emdp_b01.127C`), so `tissue` lands for every
       sample. This is what confirms N2: the file says `superior temporal gyrus` for all
       280 temporal specimens, so the tracking file was right and the build's old
       "temporal cortex" claim was wrong.
    2. **Cross-assay propagation**, on donor + tissue -- no `TMT quantitation` row
       carries a `BrodmannArea` at all, but the same donor's rnaSeq / WGS / multiome
       specimens from the same tissue do. Propagating those recovers 123 BA9 and 67 BA10.
       Marked `propagated`, never `source`, because it is another specimen's annotation.

    Anything still empty is left for the convention fill in `harmonize.apply_brodmann`,
    which never overwrites either of the above.
    """
    path = s.extra.get("biospecimen_meta")
    if not path:
        return
    b = _read_table(path, low_memory=False, dtype=str)
    tmt = b[b["assay"] == "TMT quantitation"] if "assay" in b.columns else b

    keys = slab.sample_annot["specimen_source_value"].astype(str).values
    idx = tmt.drop_duplicates("specimenID").set_index(tmt["specimenID"].astype(str))
    al = idx.reindex(keys)

    site = al["tissue"].values if "tissue" in al.columns else None
    ba = al["BrodmannArea"].values if "BrodmannArea" in al.columns else None
    ba = pd.Series(ba if ba is not None else [None] * len(keys))
    evidence = pd.Series([None] * len(keys), dtype=object)
    evidence[ba.notna().values] = "source"

    # Cross-assay propagation on donor + tissue.
    n_prop = 0
    if "individualID" in al.columns and site is not None:
        other = b[b["assay"] != "TMT quantitation"] if "assay" in b.columns else b.iloc[0:0]
        other = other.dropna(subset=["BrodmannArea"])
        if len(other):
            lut = (other.drop_duplicates(["individualID", "tissue"])
                   .set_index(["individualID", "tissue"])["BrodmannArea"].to_dict())
            want = ba.isna().values
            filled = [lut.get((i, t)) if w else None
                      for i, t, w in zip(al["individualID"].values, site, want)]
            got = pd.Series(filled).notna().values
            ba = ba.where(~got, pd.Series(filled))
            evidence[got] = "propagated"
            n_prop = int(got.sum())

    if site is not None:
        slab.sample_annot["anatomic_site_source"] = site
    slab.sample_annot["brodmann_area"] = ba.values
    slab.sample_annot["brodmann_area_evidence"] = evidence.values
    slab.notes["biospecimen_meta"] = {
        "source": path,
        "n_samples_matched": int(pd.Series(site).notna().sum()) if site is not None else 0,
        "n_samples": len(keys),
        "site_values": (sorted({str(x) for x in pd.Series(site).dropna().unique()})
                        if site is not None else []),
        "n_brodmann_source": int((evidence == "source").sum()),
        "n_brodmann_propagated": n_prop,
        "propagation_rule": "same donor + same tissue, from a non-proteomics specimen",
    }


# --------------------------------------------------------------------------------------
# Post-mortem interval (WP-7 / D-J)
# --------------------------------------------------------------------------------------

def attach_pmi(s: C.Stratum, slab: Slab) -> None:
    """Carry post-mortem interval as a biospecimen technical CDE, always in hours.

    D-J reclassified PMI from denied phenotype to an in-scope specimen attribute -- it
    describes the specimen, not the donor. That is a deliberate narrowing of plan section
    0 and REQUIREMENTS SEC-1 and is recorded as one wherever it shows: here, at
    `CLINICAL_DENY`, in both markdowns, and as a remark in the generated README.

    Null here is two different statements and is encoded as two, per row, because a
    source can mix and because the difference matters to a consumer:
      `not_applicable_antemortem`  a fact about the specimen -- a living donor's plasma
                                   or CSF has no post-mortem interval
      `source_not_acquired`        a gap in what we hold -- the object exists, we do not
    Never zero. A zero interval would be a claim.
    """
    n = len(slab.sample_source)
    spec = C.PMI_SOURCE.get(s.key, {"not_acquired": "no PMI source registered"})

    def constant(evidence, note=None):
        slab.sample_annot["pmi_hours"] = np.full(n, np.nan, dtype=np.float32)
        slab.sample_annot["pmi_evidence"] = evidence
        slab.sample_annot["pmi_source_unit"] = None
        slab.notes["pmi"] = {"evidence": evidence, "n_populated": 0, "n_samples": n,
                             "note": note}

    if spec is None:
        return constant("not_applicable_antemortem",
                        "antemortem fluid from a living donor; no PMI exists to record")
    if "not_acquired" in spec:
        return constant("source_not_acquired", spec["not_acquired"])

    # Read through an explicit usecols allow-list: ROSMAP_clinical.csv is overwhelmingly
    # clinical and only the join key and `pmi` are ever loaded into memory (SEC-2).
    d = pd.read_csv(C.ROOT / spec["path"], dtype=str,
                    usecols=[spec["key"], spec["column"]])
    d = d.assign(_k=d[spec["key"]].astype(str).str.strip(),
                 _v=pd.to_numeric(d[spec["column"]], errors="coerce"))
    d = d[d["_k"].notna() & ~d["_k"].isin(["", "nan", "NA"]) & d["_v"].notna()]
    lut = d.drop_duplicates("_k").set_index("_k")["_v"]

    person = slab.sample_annot.get("person_source_value")
    ids = (person.map(lambda v: None if pd.isna(v) else str(v).strip())
           if person is not None else pd.Series([None] * n))
    raw = ids.map(lambda v: lut.get(v) if v is not None else None).astype(float)

    stated = str(spec.get("unit") or "").strip().lower() or None
    if stated:
        factor = C.PMI_HOURS_PER_UNIT.get(stated)
        evidence = "converted_from_stated_unit"
        unit_used = stated
    else:
        # Infer from the distribution. PMI is death-to-processing, so the plausible band
        # is hours to a few days; the same durations in minutes land near 400. Where the
        # median sits in neither band the values stay null and the stratum is reported
        # undetermined -- never guessed.
        med = float(np.nanmedian(raw)) if raw.notna().any() else np.nan
        lo, hi = C.PMI_PLAUSIBLE_HOURS
        if np.isnan(med):
            factor, evidence, unit_used = None, "source_not_acquired", None
        elif lo <= med <= hi:
            factor, evidence, unit_used = 1.0, "inferred_from_distribution", "hours"
        elif lo * 60 <= med <= hi * 60:
            factor, evidence, unit_used = 1 / 60, "inferred_from_distribution", "minutes"
        else:
            factor, evidence, unit_used = None, "inferred_from_distribution", None

    if factor is None:
        hours = pd.Series([np.nan] * n)
        note = (f"unit could not be determined from the distribution "
                f"(median {np.nanmedian(raw) if raw.notna().any() else float('nan'):.3g}); "
                "values withheld rather than guessed")
    else:
        hours = raw * factor
        note = None

    lo, hi = C.PMI_PLAUSIBLE_HOURS
    got = hours.notna()
    slab.sample_annot["pmi_hours"] = hours.astype("float32").values
    slab.sample_annot["pmi_evidence"] = np.where(got, evidence, "source_not_acquired")
    slab.sample_annot["pmi_source_unit"] = np.where(got, unit_used, None)
    slab.notes["pmi"] = {
        "evidence": evidence,
        "source": spec["path"],
        "source_column": spec["column"],
        "joined_on": spec["key"],
        "source_unit": unit_used,
        "unit_stated_by_source": bool(stated),
        "factor_to_hours": factor,
        "n_populated": int(got.sum()),
        "n_samples": n,
        "median_hours": round(float(np.nanmedian(hours)), 3) if got.any() else None,
        "min_hours": round(float(np.nanmin(hours)), 3) if got.any() else None,
        "max_hours": round(float(np.nanmax(hours)), 3) if got.any() else None,
        "n_outside_plausible_band": int(((hours < lo) | (hours > hi)).sum()),
        "plausible_band_hours": [lo, hi],
        "note": note,
    }


def annot_srm_features(s: C.Stratum) -> pd.DataFrame:
    d = _read_table(s.annot_path, sep="\t", low_memory=False)
    d = d.rename(columns={d.columns[0]: "feature_source_value"})
    return d


ANNOT_READERS = {
    "dc_traits": annot_dc_traits,
    "rosmap_r1_specimen": annot_rosmap_r1_specimen,
    "somascan_samples": annot_somascan_samples,
    "olink_samples": annot_olink_samples,
    "srm_features": annot_srm_features,
    "adkp_biospecimen": annot_adkp_biospecimen,
    "rosmap_r2_traits": annot_rosmap_r2_traits,
}

# Annotation readers that need the matrix sample list before they can build a mapping.
NEEDS_SAMPLES = {"adkp_biospecimen"}


def drop_unresolvable_donors(s: C.Stratum, slab: Slab) -> None:
    """Drop non-pool samples the source cannot resolve to a donor.

    Deletion is reserved for rows whose IDENTITY is unusable, never their quality -- the
    standing policy of flagging rather than dropping applies to quality and is unchanged
    (D8). A sample with no donor cannot be joined to anything, cannot be de-duplicated
    against a repeat visit, and cannot be withdrawn if that donor later revokes consent.

    In practice this is the 3 ROSMAP SomaScan rows whose `projid` is empty in the source
    sample metadata. It is a source gap, not a join defect: the projid is genuinely
    absent, and it is NOT recoverable from the `projid_visit` specimen id -- that shortcut
    reproduces the real projid for only 918 of the 970 resolved rows, so using it here
    would invent three donor assignments that are wrong about 5% of the time.

    The specimen ids are deliberately NOT named here or in the audit row. They are
    `projid_visit` values, so naming them would publish three ROSMAP participant
    identifiers into `src/` and `assertion_audit.tsv`, both of which are in the
    publishable set (SEC-7, SEC-8: publish the rule and the count, never a real key). The
    rule is fully reproducible without them, and anyone inside the DUC can recover the
    rows from the source file in one line.
    """
    if not s.person_resolvable:
        return
    ann = slab.sample_annot
    if "person_source_value" not in ann.columns:
        return
    person = ann["person_source_value"]
    unresolved = person.isna() | person.astype(str).str.strip().isin(
        ["", "nan", "NA", "None", "<NA>"])
    pool = (ann["is_pool"].fillna(False).astype(bool) if "is_pool" in ann.columns
            else pd.Series(False, index=ann.index))
    drop = (unresolved & ~pool).to_numpy()
    if not drop.any():
        return

    keep = ~drop
    n_dropped = int(drop.sum())
    slab.values = slab.values[:, keep]
    if slab.lod is not None:
        slab.lod = slab.lod[:, keep]
    slab.sample_source = [x for x, k in zip(slab.sample_source, keep) if k]
    slab.sample_annot = ann.loc[keep].reset_index(drop=True)
    slab.notes["unresolved_donors_dropped"] = {
        "n_dropped": n_dropped,
        # Count and rule only -- never the keys. These notes reach `build_report.json`
        # and the audit TSV, which are publishable (SEC-7/SEC-8).
        "reason": "no donor identifier in the source metadata; identity unusable",
        "rule": "non-pool row whose person_source_value is null or an absent-value "
                "sentinel, after the source annotation has been joined",
    }


def read_stratum(s: C.Stratum) -> Slab:
    slab = MATRIX_READERS[s.reader](s)
    if s.annot_reader in NEEDS_SAMPLES:
        s.extra["_samples"] = slab.sample_source
    if s.annot_reader and s.annot_reader != "srm_features":
        annot = ANNOT_READERS[s.annot_reader](s)
        idx = annot.set_index("specimen_source_value")
        aligned = idx.reindex(slab.sample_source).reset_index()
        aligned = aligned.rename(columns={"index": "specimen_source_value"})
        aligned["specimen_source_value"] = slab.sample_source
        for col in aligned.columns:
            if col == "specimen_source_value":
                continue
            if col in slab.sample_annot.columns:
                slab.sample_annot[col] = slab.sample_annot[col].where(
                    slab.sample_annot[col].notna(), aligned[col])
            else:
                slab.sample_annot[col] = aligned[col].values
        if "specimen_source_value" not in slab.sample_annot.columns:
            slab.sample_annot.insert(0, "specimen_source_value", slab.sample_source)
        if "_adkp_exclusions" in s.extra:
            slab.notes["adkp_exclusions"] = s.extra.pop("_adkp_exclusions")
    elif s.annot_reader == "srm_features":
        slab.feature_annot = ANNOT_READERS["srm_features"](s)
    if slab.sample_annot.empty:
        slab.sample_annot = pd.DataFrame({"specimen_source_value": slab.sample_source})

    # The metadata objects WP-1 landed, spent (WP-8) rather than left on disk. These run
    # after the annotation merge because each of them needs what it produced: the assay
    # join keys off `specimen_native_id`, PMI off `person_source_value`.
    attach_assay_metadata(s, slab)
    attach_biospecimen_metadata(s, slab)
    # Before PMI, so the conversion ledger describes the rows that actually ship rather
    # than counting three specimens that are about to be removed.
    drop_unresolvable_donors(s, slab)
    attach_pmi(s, slab)
    return slab
