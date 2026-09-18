"""S4-S7 - feature resolution, missingness filtering, per-stratum Z, and merge.

Scale policy (plan section 5, as directed):
  1. best transform per subset  -- monotone and invertible, declared in transform_chain
  2. Z-scale each subset SEPARATELY, per feature, on that subset's own samples
  3. only then merge subsets into a layer

Rationale: strata within a layer come from different pipelines on different scales
(brain TMT spans log2 ratios near zero and log2 reporter intensities near 20). Merging
before standardising would let scale dominate biology. Z after transform, computed
within subset, puts every subset on a common dimensionless footing while keeping the
per-feature centre and spread stored so the operation is reversible.

Single-stratum layers are NOT Z-scaled -- there is nothing to merge, so native scale is
retained for fidelity. Both representations are emitted where they differ.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import config as C, readers as R

# Docs precedent (QC_process_Frontal.txt / QC_process_Temporal.txt):
#   "only keep proteins that have missing data in <50% of the samples"
MAX_MISSING_FRACTION = 0.50

# Layers whose features must NOT be collapsed to one row per accession. SRM measures
# peptides (including phospho-specific tau species that are biologically distinct from
# the parent protein); SomaScan measures aptamers, many of which target the same protein
# with different affinity. Collapsing either would destroy the measured entity.
NO_COLLAPSE_LAYERS = frozenset({"srm_brain_proteomics", "soma_plasma_proteomics"})


# A layer is left alone -- no transform, no Z -- when its strata are on the same declared
# scale AND empirically similar in distribution. These are the tolerances for "similar":
# the strata's typical centres must sit within LOC_TOL pooled SDs of each other, and their
# typical spreads must be within a factor of SPREAD_TOL.
LOC_TOL = 0.5
SPREAD_TOL = 2.0


def distributions_similar(parts: list["Harmonised"]) -> dict:
    """Compare strata empirically rather than trusting the declared scale label."""
    if len(parts) < 2:
        return {"n_strata": len(parts), "similar": True,
                "reason": "single stratum - nothing to reconcile"}

    stats = []
    for p in parts:
        v = p.native
        with np.errstate(invalid="ignore"):
            centres = np.nanmedian(v, axis=1)
            spreads = np.nanstd(v, axis=1, ddof=1)
        stats.append({
            "stratum": p.stratum.key,
            "centre": float(np.nanmedian(centres)),
            "spread": float(np.nanmedian(spreads)),
        })

    centres = np.array([s["centre"] for s in stats])
    spreads = np.array([s["spread"] for s in stats])
    pooled = float(np.nanmedian(spreads))
    loc_spread = float(np.nanmax(centres) - np.nanmin(centres))
    loc_ratio = loc_spread / pooled if pooled > 0 else np.inf
    pos = spreads[spreads > 0]
    spread_ratio = float(np.nanmax(pos) / np.nanmin(pos)) if len(pos) else np.inf

    ok = (loc_ratio <= LOC_TOL) and (spread_ratio <= SPREAD_TOL)
    return {
        "n_strata": len(parts),
        "per_stratum": stats,
        "location_range_in_pooled_sd": round(loc_ratio, 3),
        "spread_ratio_max_over_min": round(spread_ratio, 3),
        "similar": bool(ok),
        "reason": ("centres and spreads comparable across strata"
                   if ok else
                   f"strata differ materially (location {loc_ratio:.2f} pooled SD, "
                   f"spread ratio {spread_ratio:.2f})"),
    }


def decide_layer_policy(layer: str, parts: list["Harmonised"]) -> dict:
    """Decide whether a layer needs Z, using declared scale AND measured distribution."""
    scales = sorted({p.stratum.value_scale for p in parts})
    sim = distributions_similar(parts)

    # The decision is made on the measured distributions, not on the scale labels.
    # Labels are carried for documentation, and any disagreement between what the
    # label claims and what the data shows is surfaced rather than silently resolved.
    leave_alone = sim["similar"]
    label_says_same = len(scales) == 1
    disagreement = None
    if label_says_same and not sim["similar"]:
        disagreement = ("declared scales are identical but the measured distributions "
                        "are not - trusting the data")
    elif not label_says_same and sim["similar"]:
        disagreement = ("declared scales differ but the measured distributions are "
                        "comparable - trusting the data, no Z applied")

    return {
        "layer": layer,
        "declared_scales": scales,
        "same_declared_scale": label_says_same,
        "label_vs_data": disagreement,
        "distribution_check": sim,
        "z_applied": not leave_alone,
        "value_semantics": ("native scale, unmodified"
                            if leave_alone else
                            "per-stratum Z (native also retained)"),
        "rationale": (
            "measured centres and spreads are comparable across strata, so the native "
            "values are kept exactly as delivered"
            if leave_alone else
            f"measured distributions differ across strata ({sim['reason']}), so each "
            "stratum is Z-scored on its own samples before merging"),
    }


# --------------------------------------------------------------------------------------
# Feature resolution -> UniProt primary accession
# --------------------------------------------------------------------------------------

_ACC = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"


def _load_uniprot2symbol() -> dict[str, str]:
    d = R._read_table(C.UNIPROT2SYMBOL)
    d.columns = ["uniprot", "symbol"]
    d = d.dropna().astype(str)
    # Prefer the first accession listed for a symbol (reference is symbol-sorted).
    return d.drop_duplicates("symbol").set_index("symbol")["uniprot"].to_dict()


def resolve_features(stratum: C.Stratum, slab: R.Slab,
                     sym2uni: dict[str, str]) -> pd.DataFrame:
    """One row per source feature with feature_id (UniProt) and provenance."""
    src = [str(x) for x in slab.feature_source]
    fa = slab.feature_annot
    n = len(src)
    fid = [None] * n
    gene = [None] * n
    group = [None] * n
    gsize = [1] * n
    peptide = [None] * n

    reader = stratum.reader

    if reader in ("csv_matrix", "pd_ratio"):
        # "GENE|UNIPROT" (gene may be empty, or literal "NA")
        for i, s in enumerate(src):
            g, _, u = s.partition("|")
            g = g.strip().strip('"')
            fid[i] = u.strip() or None
            gene[i] = None if g in ("", "NA", "nan") else g
            group[i] = u.strip() or None

    elif reader == "msbb_tmt_xlsx":
        # "GENE|sp|ACC|NAME_HUMAN"
        for i, s in enumerate(src):
            g, _, rest = s.partition("|")
            m = re.search(_ACC, rest)
            fid[i] = m.group(0) if m else None
            gene[i] = g if g not in ("", "nan", "None") else None
            group[i] = rest

    elif reader == "maxquant":
        # Semicolon-delimited protein group; leading accession is feature_id (D6/#6).
        for i, s in enumerate(src):
            members = [m for m in s.split(";") if m]
            accs = []
            for m in members:
                hit = re.search(_ACC, m)
                if hit:
                    accs.append(hit.group(0))
            fid[i] = accs[0] if accs else None
            group[i] = ";".join(accs) if accs else s
            gsize[i] = max(len(accs), 1)

    elif reader == "srm":
        # Peptide features keyed by gene symbol with a suffix: BIN1_3, tau_PHF1_s404.
        pep = (fa.set_index(fa.columns[0])["Peptide.Sequence"].to_dict()
               if not fa.empty and "Peptide.Sequence" in fa.columns else {})
        for i, s in enumerate(src):
            base = s
            if base.startswith("tau"):
                sym = "MAPT"
            else:
                sym = re.sub(r"_\d+$", "", base)
            fid[i] = C.SRM_MANUAL_UNIPROT.get(sym) or sym2uni.get(sym)
            gene[i] = sym if sym != "MAPT" else "MAPT"
            group[i] = fid[i]
            peptide[i] = pep.get(s)

    elif reader == "somascan":
        u = fa["uniprot_raw"].astype(str).tolist()
        g = fa["gene_raw"].astype(str).tolist()
        for i in range(n):
            accs = re.findall(_ACC, u[i])
            fid[i] = accs[0] if accs else None
            group[i] = ";".join(accs) if accs else None
            gsize[i] = max(len(accs), 1)
            gene[i] = None if g[i] in ("nan", "None", "") else g[i]

    elif reader in ("olink", "dia"):
        # Both deliver a bare UniProt accession per feature, so FR-2 is satisfied without
        # a symbol lookup -- confirmed on the DIA delivery rather than assumed.
        u = fa["uniprot_raw"].astype(str).tolist()
        for i in range(n):
            m = re.search(_ACC, u[i])
            fid[i] = m.group(0) if m else u[i]
            group[i] = fid[i]

    out = pd.DataFrame({
        "feature_id": fid,
        "feature_source_value": src,
        "gene_symbol": gene,
        "protein_group": group,
        "protein_group_size": np.asarray(gsize, dtype=np.int16),
        "measurement_concept_id": np.zeros(n, dtype=np.int64),
    })
    # backfill gene symbol from the reference where the source did not carry one
    uni2sym = {v: k for k, v in sym2uni.items()}
    need = out.gene_symbol.isna() & out.feature_id.notna()
    out.loc[need, "gene_symbol"] = out.loc[need, "feature_id"].map(uni2sym)

    layer_extra = C.FEATURE_COLS_BY_LAYER.get(stratum.layer, [])
    if "panel" in layer_extra:
        out["panel"] = fa["panel"].values if "panel" in fa else None
        out["panel_lot_nr"] = fa["panel_lot_nr"].values if "panel_lot_nr" in fa else None
        # WP-12 item 4: the vendor's own per-assay missingness, carried beside the
        # `frac_below_lod` this build computes so the two can be compared rather than
        # one being quietly preferred.
        out["missing_freq"] = fa["missing_freq"].values if "missing_freq" in fa else None
    if "seq_id" in layer_extra:
        out["seq_id"] = fa["seq_id"].values if "seq_id" in fa else None
        out["dilution"] = fa["dilution"].values if "dilution" in fa else None
        # WP-12 item 2: SomaLogic's per-analyte column QC, declared ours to apply by
        # SOURCE_QC and previously dropped at read.
        out["colcheck"] = fa["colcheck"].values if "colcheck" in fa else None
    if "peptide_sequence" in layer_extra:
        out["peptide_sequence"] = peptide
    return out


# --------------------------------------------------------------------------------------
# Missingness filter and Z
# --------------------------------------------------------------------------------------

def missingness_filter(values: np.ndarray) -> np.ndarray:
    """Keep features present in >= 50% of samples.

    Applied PER STRATUM, never after merging. Each stratum has its own sample set, so
    presence must be judged within it -- a protein measured well by one assay would look
    mostly-missing across a merged layer and be dropped for the wrong reason.

    Threshold follows the source documentation (QC_process_Frontal.txt /
    QC_process_Temporal.txt): "only keep proteins that have missing data in <50% of the
    samples".
    """
    present = np.isfinite(values).mean(axis=1)
    return present >= (1.0 - MAX_MISSING_FRACTION)


# QC each source's own documentation specifies. Anything marked upstream=True was already
# applied by the data producer before deposit; we verify rather than re-apply. Anything
# upstream=False we apply here, per stratum, BEFORE any transform, Z or merge.
SOURCE_QC = {
    "diversecohorts_dlpfc": [
        ("GIS/non-GIS separation", True, "QC_process_Frontal.txt step 1"),
        ("features <50% missing", True, "QC_process_Frontal.txt step 2"),
        ("ratio to total abundance, log2", True, "QC_process_Frontal.txt step 2"),
        ("iterative PCA outlier removal (19 removed)", True, "QC_process_Frontal.txt step 3"),
        ("batch effect regressed out", True, "QC_process_Frontal.txt step 4"),
    ],
    "diversecohorts_temporal": [
        ("GIS/non-GIS separation", True, "QC_process_Temporal.txt step 1"),
        ("features <50% missing", True, "QC_process_Temporal.txt step 2"),
        ("ratio to total abundance, log2", True, "QC_process_Temporal.txt step 2"),
        ("iterative PCA outlier removal (2 removed)", True, "QC_process_Temporal.txt step 3"),
        ("batch effect regressed out", True, "QC_process_Temporal.txt step 4"),
    ],
    "rosmap_tmt_r1": [
        ("GIS channels excluded", True, "delivered C2 matrix is 8 channels/batch"),
        ("TAMPOR median-polish batch correction", True, "TAMPOR.R + batch correction PDF"),
    ],
    "msbb_tmt_phg": [
        ("8 low-quality samples removed (198 -> 190)", True, "xlsx header note"),
    ],
    "mayo_lfq_tcx": [
        ("reverse / contaminant / site-only rows removed", False, "MaxQuant convention"),
    ],
    "msbb_lfq_pfc": [
        ("reverse / contaminant / site-only rows removed", False, "MaxQuant convention"),
    ],
    "rosmap_soma_plasma": [
        ("ANML normalisation", True, "deposit README (SomaLogic recommendation)"),
        ("ColCheck per-analyte QC flagged", False, "SomaScan protein metadata"),
    ],
    "pdrd_olink_plasma": [
        ("D01<->D02 bridging", True, "AMP-PD README"),
        ("QC_Warning / Distribution_QC / Outliers_QC flagged", False, "AMP-PD README"),
        ("below-LOD fraction computed per assay", False,
         "Olink recommendation: consider excluding assays with 25-50%+ below LOD"),
    ],
    "pdrd_olink_csf": [
        ("D01<->D02 bridging", True, "AMP-PD README"),
        ("QC_Warning / Distribution_QC / Outliers_QC flagged", False, "AMP-PD README"),
        ("below-LOD fraction computed per assay", False,
         "Olink recommendation: consider excluding assays with 25-50%+ below LOD"),
    ],
    "rosmap_tmt_r2": [
        ("126 GIS reference channel excluded", True,
         "PD export reports ratios vs 126, so the reference channel is not a sample"),
        ("ratio to internal reference standard", True, "PD Abundance Ratio (ch / 126)"),
        ("log2 of the reference ratio", False, "standard TMT practice"),
        ("per-protein centring on median of per-plex medians", False,
         "matches the round-1 product construction, so both ROSMAP rounds are comparable"),
        ("features <50% missing", False, "same threshold as every other stratum"),
    ],
    "rosmap_srm_panel1": [("control wells flagged", False, "isControl / sample.type")],
    "rosmap_srm_panel2": [("control wells flagged", False, "isControl / sample.type")],
}


# Transform policy: a monotone transform is applied ONLY where the source scale needs it
# to be approximately normal. Sources already delivered on a log scale are Z-scored
# directly, with no transform at all. The decision is verified here rather than asserted:
# we measure skew on what we actually read and record it.
SKEW_TOLERANCE = 1.0


def scale_diagnostics(values: np.ndarray, transform_applied: str) -> dict:
    """Measure whether the (possibly transformed) values are near-normal per feature."""
    finite = np.isfinite(values)
    n = finite.sum(axis=1)
    usable = n >= 8
    if not usable.any():
        return {"transform_applied": transform_applied, "median_skew": None,
                "pct_features_skew_gt_1": None, "normality_ok": None}
    v = np.where(finite, values, np.nan)[usable]
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nanmean(v, axis=1)
        sd = np.nanstd(v, axis=1, ddof=1)
        sd = np.where(sd > 0, sd, np.nan)
        skew = np.nanmean(((v - mu[:, None]) / sd[:, None]) ** 3, axis=1)
    med = float(np.nanmedian(np.abs(skew)))
    frac = float(np.nanmean(np.abs(skew) > SKEW_TOLERANCE))
    return {
        "transform_applied": transform_applied,
        "median_abs_skew": round(med, 3),
        "pct_features_skew_gt_1": round(100 * frac, 1),
        "normality_ok": bool(med <= SKEW_TOLERANCE),
        "transform_rationale": (
            "source already on a log scale; Z-scored directly with no transform"
            if transform_applied == "none" else
            f"source was linear; {transform_applied} applied to reach approximate "
            "normality before Z"),
    }


def rank_inverse_normal(values: np.ndarray) -> np.ndarray:
    """Per-feature rank-based inverse normal transform (Blom).

    Monotone and rank-preserving, so ordering and non-parametric structure survive
    exactly. Used only where the source scale is still materially skewed after any
    log step -- see the normality gate in harmonise_stratum.
    """
    from scipy.stats import norm, rankdata

    out = np.full(values.shape, np.nan, dtype=np.float32)
    for i in range(values.shape[0]):
        row = values[i]
        m = np.isfinite(row)
        k = int(m.sum())
        if k > 3:
            out[i, m] = norm.ppf((rankdata(row[m]) - 0.375) / (k + 0.25))
    return out


def zscore(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-feature Z over this subset's own samples, NaN-aware and reversible."""
    with np.errstate(invalid="ignore"):
        centre = np.nanmean(values, axis=1)
        spread = np.nanstd(values, axis=1, ddof=1)
    safe = np.where(np.isfinite(spread) & (spread > 0), spread, np.nan)
    z = (values - centre[:, None]) / safe[:, None]
    return z.astype(np.float32), centre.astype(np.float32), safe.astype(np.float32)


def observed_platform(p) -> tuple[str | None, str | None]:
    """The stratum's platform as the assay metadata actually reports it (WP-8).

    Returns `(value, evidence)`, falling back to the config declaration where no metadata
    was joined. A stratum whose samples ran on several instruments says so -- the
    DiverseCohorts strata really do span Lumos, Exploris 240 and Eclipse -- because
    collapsing that to one name would be an assertion. The per-sample `platform` column
    is authoritative per row; this is the stratum-level summary of it.
    """
    meta = (getattr(p, "notes", None) or {}).get("assay_meta") or {}
    obs = meta.get("platform_observed") if meta.get("joined") else None
    if not obs:
        return p.stratum.platform, (p.stratum.platform_evidence or "C")
    names = sorted(obs)
    if len(names) == 1:
        return names[0], "A"
    total = sum(obs.values())
    return ("mixed per sample: "
            + "; ".join(f"{k} ({obs[k]}/{total})" for k in names)), "A"


# --------------------------------------------------------------------------------------
# Brodmann area (WP-9 / D-B)
# --------------------------------------------------------------------------------------

def apply_brodmann(s: C.Stratum, samples: pd.DataFrame) -> dict:
    """Fill `brodmann_area` by convention where, and only where, the source is silent.

    Three evidence states, per row, in `brodmann_area_evidence`:
      `source`      the specimen's own annotation stated it
      `propagated`  another specimen from the same donor and tissue stated it (WP-9)
      `convention`  neither did, and `BRODMANN_CONVENTION` supplied it from the site

    **A source-stated value always wins.** The convention never overwrites one -- which
    is what keeps the 67 DiverseCohorts BA10 rows out of the BA9 bucket the convention
    would otherwise sweep them into. Rows where a stated value disagrees with the
    convention are counted and reported rather than reconciled.

    This matters downstream and not only here: PROTEOMICS_REVIEW.md has `brodmann_area`
    taking precedence over the free-text site label, so a CDM loader that reads this
    column without reading its evidence would move ~4,400 DLPFC rows from
    `4195656 Region of frontal cortex` onto a BA9 concept on our inference alone. The
    evidence column is what makes that inference visible and reversible.
    """
    n = len(samples)
    if "brodmann_area" not in samples.columns:
        samples["brodmann_area"] = None
    if "brodmann_area_evidence" not in samples.columns:
        samples["brodmann_area_evidence"] = None

    ba = samples["brodmann_area"].where(
        samples["brodmann_area"].notna()
        & ~samples["brodmann_area"].astype(str).str.strip().isin(["", "nan", "NA", "None"]))
    ev = samples["brodmann_area_evidence"].astype(object)

    # A reader that supplied the value but no evidence (MSBB / Mayo via the ADKP
    # biospecimen file) is stating the specimen's own annotation.
    ev = ev.where(ev.notna() | ba.isna(), "source")

    site = str(s.anatomic_site or "").strip().lower()
    conv = C.BRODMANN_CONVENTION.get(site)
    n_conv = 0
    n_disagree = 0
    if conv is not None:
        disagree = ba.notna() & (ba.astype(str).str.strip().str.upper() != conv.upper())
        n_disagree = int(disagree.sum())
        need = ba.isna()
        n_conv = int(need.sum())
        ba = ba.where(~need, conv)
        ev = ev.where(~need, "convention")

    samples["brodmann_area"] = ba.values
    samples["brodmann_area_evidence"] = ev.values
    return {
        "anatomic_site": s.anatomic_site,
        "convention": conv,
        "n_samples": n,
        "n_source": int((ev == "source").sum()),
        "n_propagated": int((ev == "propagated").sum()),
        "n_convention": n_conv,
        "n_unfilled": int(ba.isna().sum()),
        "n_source_disagrees_with_convention": n_disagree,
        "disagreeing_values": (sorted({str(x) for x in
                                       ba[ba.notna()].astype(str).str.strip().unique()
                                       if conv and str(x).upper() != conv.upper()})
                               if conv else []),
    }


# --------------------------------------------------------------------------------------
# Person identity
# --------------------------------------------------------------------------------------

# Which identifier system each stratum's person_source_value belongs to. IDs are NEVER
# stripped of prefixes -- `PD-`, `PP-`, `R`, `AMPAD_MSSM_` are part of the identifier.
# The namespace is recorded alongside so that two studies reusing the same digits are
# not silently merged into one person.
PERSON_NAMESPACE = {
    "diversecohorts_dlpfc": "AMP-AD.individualID",
    "diversecohorts_temporal": "AMP-AD.individualID",
    "rosmap_tmt_r1": "ROSMAP.projid",
    "rosmap_tmt_r2": "ROSMAP.projid",
    "msbb_tmt_phg": "AMP-AD.individualID",
    "mayo_lfq_tcx": "AMP-AD.individualID",
    "msbb_lfq_pfc": "AMP-AD.individualID",
    "rosmap_srm_panel1": "ROSMAP.projid",
    "rosmap_srm_panel2": "ROSMAP.projid",
    "rosmap_soma_plasma": "ROSMAP.individualID",
    "pdrd_olink_plasma": "AMP-PD.participant_id",
    "pdrd_olink_csf": "AMP-PD.participant_id",
}

_FLOAT_ARTIFACT = re.compile(r"^(\d+)\.0+$")

# Strings a source writes to mean "no identifier". They must become NULL, never travel as
# an identifier -- MSBB's biospecimen file spells its absent donor `Unknown`, and taking
# that literally produced a person called `AMP-AD.individualID:Unknown` shared by every
# row carrying it. It reached the artifact: 4 rows of `msbb_tmt_phg`, all of them GIS
# reference channels, in this build and in its predecessor. Contained, because those rows
# are pools rather than donors, but a sentinel that survives into an identifier column
# will merge unrelated specimens into one fake person as soon as two strata use it.
_ID_SENTINELS = frozenset({
    "", "nan", "na", "n/a", "none", "null", "<na>", "unknown", "unk", "missing",
    "not applicable", "not collected", "not reported",
})


def as_source_id(value) -> str | None:
    """Return the identifier exactly as the source defines it, as a string.

    IDs are strings throughout -- never numeric -- so they survive round-trips and join
    to future metadata drops unchanged. Prefixes (PD-, PP-, R, AMPAD_MSSM_) are part of
    the identifier and are never stripped.

    Two repairs, both restoring the source value rather than altering it: a trailing `.0`,
    which is damage introduced by reading a column containing blanks through pandas
    (1005 -> 1005.0); and an absent-value sentinel, which is the source saying "no donor"
    in words and must become a null rather than a donor named `Unknown`.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    s = str(value).strip()
    if s.lower() in _ID_SENTINELS:
        return None
    m = _FLOAT_ARTIFACT.match(s)
    return m.group(1) if m else s


def attach_person_identity(samples: pd.DataFrame, stratum_key: str) -> pd.DataFrame:
    """person_id IS the source identifier, as a string. No surrogate is invented.

    Because two studies can legitimately reuse the same digits, the identifier system is
    recorded in person_source_namespace and a namespace-qualified `person_global_key` is
    provided for unambiguous cross-study joins. The bare source ID is preserved for
    joining to that study's own future metadata.
    """
    ns_default = PERSON_NAMESPACE.get(stratum_key, "unknown")
    per_row = samples.get("person_source_namespace_row")
    ns = (per_row.fillna(ns_default) if per_row is not None
          else pd.Series([ns_default] * len(samples), index=samples.index))
    raw = samples.get("person_source_value")
    ids = (raw.map(as_source_id) if raw is not None
           else pd.Series([None] * len(samples), index=samples.index))
    samples = samples.copy()
    samples["person_source_value"] = ids.astype("string")
    samples["person_id"] = ids.astype("string")          # string, source-exact
    samples["person_source_namespace"] = pd.Series(ns, index=samples.index).astype("string")
    samples["person_global_key"] = pd.Series(
        np.where(ids.notna(), pd.Series(ns, index=samples.index).astype(str)
                 + ":" + ids.astype(str), None), index=samples.index).astype("string")
    return samples


# --------------------------------------------------------------------------------------
# Per-stratum harmonisation
# --------------------------------------------------------------------------------------

class Harmonised:
    def __init__(self, stratum, features, samples, native, z, centre, spread, lod, notes):
        self.stratum = stratum
        self.features = features          # DataFrame, len == native.shape[0]
        self.samples = samples            # DataFrame, len == native.shape[1]
        self.native = native
        self.z = z
        self.z_center = centre
        self.z_scale = spread
        self.lod = lod
        self.notes = notes
        self.assayed_keys: set[str] = set()
        self.qc_dropped_keys: set[str] = set()


def harmonise_stratum(s: C.Stratum, sym2uni: dict[str, str]) -> Harmonised:
    slab = R.read_stratum(s)
    feats = resolve_features(s, slab, sym2uni)
    values = slab.values
    lod = slab.lod
    notes = dict(slab.notes)

    # 1. drop features we could not resolve to an accession
    ok = feats.feature_id.notna().values
    notes["n_features_unresolved_dropped"] = int((~ok).sum())

    # feature_key is computed BEFORE filtering so we can tell, at merge time, whether a
    # feature was absent from this stratum's assay entirely or was dropped by its QC.
    if feats.feature_id.duplicated().any():
        _dup = feats.feature_id.duplicated(keep=False)
        feats["feature_key"] = np.where(
            _dup, feats.feature_id.astype(str) + "__" + feats.feature_source_value.astype(str),
            feats.feature_id.astype(str))
    else:
        feats["feature_key"] = feats.feature_id.astype(str)

    # 2. docs-precedent missingness filter
    keep_missing = missingness_filter(values)
    notes["n_features_below_50pct_present_dropped"] = int((~keep_missing & ok).sum())

    keep = ok & keep_missing
    assayed_keys = set(feats.feature_key[ok].astype(str))
    qc_dropped_keys = set(feats.feature_key[ok & ~keep_missing].astype(str))

    values = values[keep]
    feats = feats[keep].reset_index(drop=True)
    if lod is not None:
        lod = lod[keep]

    # 3. collapse duplicate accessions within a stratum (mean), keeping provenance.
    #    Skipped for layers where the measured entity is finer than the protein.
    if s.layer not in NO_COLLAPSE_LAYERS and feats.feature_id.duplicated().any():
        order = feats.feature_id.values
        uniq, inv = np.unique(order, return_inverse=True)
        acc = np.zeros((len(uniq), values.shape[1]), dtype=np.float64)
        cnt = np.zeros((len(uniq), values.shape[1]), dtype=np.int32)
        fin = np.isfinite(values)
        np.add.at(acc, inv, np.where(fin, values, 0.0))
        np.add.at(cnt, inv, fin.astype(np.int32))
        with np.errstate(invalid="ignore", divide="ignore"):
            merged = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan).astype(np.float32)
        agg = (feats.assign(_i=inv).groupby("_i")
               .agg({c: "first" for c in feats.columns if c != "feature_id"}))
        agg.insert(0, "feature_id", uniq)
        notes["n_duplicate_accessions_collapsed"] = int(len(feats) - len(uniq))
        if lod is not None:
            lacc = np.zeros((len(uniq), lod.shape[1]))
            lcnt = np.zeros((len(uniq), lod.shape[1]), dtype=np.int32)
            lf = np.isfinite(lod)
            np.add.at(lacc, inv, np.where(lf, lod, 0.0))
            np.add.at(lcnt, inv, lf.astype(np.int32))
            with np.errstate(invalid="ignore", divide="ignore"):
                lod = np.where(lcnt > 0, lacc / np.maximum(lcnt, 1), np.nan).astype(np.float32)
        values, feats = merged, agg.reset_index(drop=True)

    notes["n_features_key_disambiguated"] = int(
        feats.feature_key.astype(str).str.contains("__").sum())

    # Transform decision (measured, not asserted): readers apply log2 only to sources
    # delivered on a linear scale. Everything else reaches this point untransformed.
    applied = "log2" if "log2" in str(notes.get("transform", "")) else "none"
    notes["transform_applied"] = applied
    notes["scale_check"] = scale_diagnostics(values, applied)
    notes["zscored"] = None          # decided at layer level, pass 2

    samples = slab.sample_annot.copy()
    if "specimen_source_value" not in samples.columns:
        samples.insert(0, "specimen_source_value", slab.sample_source)
    samples = samples.reset_index(drop=True)
    samples["stratum_key"] = s.key
    notes["brodmann"] = apply_brodmann(s, samples)

    # Source-specified QC: verify what the producer already applied, apply what they
    # specified but left to the consumer. All of this happens before transform/Z/merge.
    notes["source_qc"] = [
        {"check": name, "applied_upstream": up, "reference": ref}
        for name, up, ref in SOURCE_QC.get(s.key, [])
    ]
    if lod is not None:
        below = np.isfinite(values) & np.isfinite(lod) & (values < lod)
        frac = below.sum(axis=1) / np.maximum(np.isfinite(values).sum(axis=1), 1)
        feats["frac_below_lod"] = frac.astype(np.float32)
        feats["olink_lod_flag"] = np.where(
            frac >= 0.50, "gt50pct_below_lod",
            np.where(frac >= 0.25, "gt25pct_below_lod", "ok"))
        notes["n_assays_gt25pct_below_lod"] = int((frac >= 0.25).sum())
        notes["n_assays_gt50pct_below_lod"] = int((frac >= 0.50).sum())
        notes["lod_policy"] = ("assays are flagged, not dropped - Olink leaves the "
                               "exclusion threshold to the study")

    notes["n_features_final"] = int(values.shape[0])
    notes["n_samples"] = int(values.shape[1])

    # A stratum with no surviving features is always a wiring fault, not a data fact --
    # most often a reader with no branch in `resolve_features`, which leaves every
    # accession unresolved and silently drops the lot. Say so here rather than failing
    # several lines later inside a numpy reduction over an empty array.
    if values.shape[0] == 0:
        raise ValueError(
            f"{s.key}: no features survived resolution "
            f"({notes['n_features_unresolved_dropped']} unresolved, "
            f"{notes['n_features_below_50pct_present_dropped']} below the missingness "
            f"floor). Does resolve_features handle reader {s.reader!r}?")

    # Sparsity of what actually ships, after filtering.
    fin = np.isfinite(values)
    per_feat = fin.mean(axis=1)
    notes["sparsity"] = {
        "pct_missing_overall": round(float(100 * (1 - fin.mean())), 2),
        "pct_features_fully_observed": round(float(100 * (per_feat == 1).mean()), 1),
        "median_feature_presence_pct": round(float(100 * np.median(per_feat)), 1),
        "min_feature_presence_pct": round(float(100 * per_feat.min()), 1),
    }
    h = Harmonised(s, feats, samples, values, None, None, None, lod, notes)
    h.assayed_keys = assayed_keys
    h.qc_dropped_keys = qc_dropped_keys
    return h


# Missing-reason codes. These live in a SEPARATE int8 matrix, never inside `values`.
# `values` stays pure float32 with an ordinary NaN in every missing cell, so numpy,
# pandas and sklearn all behave normally; anyone who does not care about the reason can
# ignore this array entirely. Encoding two null kinds inside the float matrix would break
# ordinary analysis -- Python has only one usable float null, NaN payloads are not
# preserved by arithmetic or dtype changes, and a numeric sentinel would silently poison
# means, correlations and PCA.
MISSING_OBSERVED = 0        # value present
MISSING_NOT_DETECTED = 1    # assayed in this stratum, not quantified in this sample (NA)
MISSING_NOT_ASSAYED = 2     # column exists only because another stratum measures it
MISSING_FAILED_QC = 3       # present in this stratum's source but dropped by its QC

MISSING_LEGEND = {
    MISSING_OBSERVED: "observed",
    MISSING_NOT_DETECTED: "not_detected_in_sample (assayed in this stratum)",
    MISSING_NOT_ASSAYED: "not_assayed_in_this_stratum (column comes from another stratum)",
    MISSING_FAILED_QC: "failed_qc_in_this_stratum (present in source, dropped by QC)",
}


def apply_layer_policy(parts: list[Harmonised], policy: dict) -> None:
    """Pass 2. Given the layer decision, apply the normality gate and Z where needed.

    Nothing here touches a layer that keeps native values -- those are emitted exactly
    as delivered.
    """
    if not policy["z_applied"]:
        for p in parts:
            p.notes["zscored"] = False
            p.notes["value_semantics"] = "native scale, unmodified"
        return

    for p in parts:
        # Normality gate: Z on a badly skewed feature is not interpretable, so a
        # rank-based inverse normal transform is applied first where the measurement
        # says it is needed.
        if p.notes["scale_check"]["normality_ok"] is False:
            p.native = rank_inverse_normal(p.native)
            applied = p.notes["transform_applied"]
            applied = "rank_int" if applied == "none" else f"{applied}+rank_int"
            p.notes["transform_applied"] = applied
            p.notes["scale_check_after_transform"] = scale_diagnostics(p.native, applied)
            p.notes["rank_int_reason"] = (
                "median |skew| exceeded tolerance and this stratum is Z-scored; a "
                "monotone rank-based inverse normal transform was applied first so the "
                "Z score is interpretable")
        p.z, p.z_center, p.z_scale = zscore(p.native)
        p.notes["zscored"] = True
        p.notes["value_semantics"] = "per-stratum Z (native retained alongside)"


def verify_common_numeric_space(parts: list[Harmonised], policy: dict) -> dict:
    """Post-condition: every stratum in a layer must share one numeric space at merge.

    Checked on the values that are actually merged -- Z where Z was applied, native
    otherwise -- so the guarantee is about what ships, not about intent.
    """
    if len(parts) < 2:
        return {"verified": True, "n_strata": len(parts),
                "reason": "single stratum, trivially one space",
                "merged_surface": "values_z" if policy["z_applied"] else "values"}

    surface, stats = [], []
    for p in parts:
        v = p.z if (policy["z_applied"] and p.z is not None) else p.native
        surface.append(v)
        with np.errstate(invalid="ignore"):
            stats.append({
                "stratum": p.stratum.key,
                "centre": float(np.nanmedian(np.nanmedian(v, axis=1))),
                "spread": float(np.nanmedian(np.nanstd(v, axis=1, ddof=1))),
            })
    centres = np.array([s["centre"] for s in stats])
    spreads = np.array([s["spread"] for s in stats])
    pooled = float(np.nanmedian(spreads))
    loc = float(np.nanmax(centres) - np.nanmin(centres)) / (pooled if pooled > 0 else np.inf)
    pos = spreads[spreads > 0]
    sr = float(np.nanmax(pos) / np.nanmin(pos)) if len(pos) else np.inf
    ok = (loc <= LOC_TOL) and (sr <= SPREAD_TOL)
    return {
        "verified": bool(ok),
        "n_strata": len(parts),
        "merged_surface": "values_z" if policy["z_applied"] else "values",
        "per_stratum": stats,
        "location_range_in_pooled_sd": round(loc, 3),
        "spread_ratio_max_over_min": round(sr, 3),
        "reason": ("all strata occupy one numeric space at merge" if ok else
                   "strata still differ after the policy was applied"),
    }


def merge_layer(parts: list[Harmonised], policy: dict) -> dict:
    """Union features across strata (on feature_key), concatenate samples."""
    keys = sorted({k for p in parts for k in p.features.feature_key})
    pos = {k: i for i, k in enumerate(keys)}
    n_f, n_s = len(keys), sum(p.native.shape[1] for p in parts)

    native = np.full((n_f, n_s), np.nan, dtype=np.float32)
    zed = (np.full((n_f, n_s), np.nan, dtype=np.float32)
           if policy["z_applied"] else None)
    lod = (np.full((n_f, n_s), np.nan, dtype=np.float32)
           if any(p.lod is not None for p in parts) else None)

    # Companion reason matrix: same shape, int8, never mixed into `values`.
    reason = np.full((n_f, n_s), MISSING_NOT_ASSAYED, dtype=np.int8)

    col = 0
    samp_frames, feat_frames = [], []
    for p in parts:
        rows = np.array([pos[k] for k in p.features.feature_key])
        w = p.native.shape[1]
        native[rows, col:col + w] = p.native
        if zed is not None and p.z is not None:
            zed[rows, col:col + w] = p.z
        if lod is not None and p.lod is not None:
            lod[rows, col:col + w] = p.lod

        block = np.full((n_f, w), MISSING_NOT_ASSAYED, dtype=np.int8)
        # features this stratum kept: observed vs not detected in that sample
        block[rows] = np.where(np.isfinite(p.native),
                               MISSING_OBSERVED, MISSING_NOT_DETECTED).astype(np.int8)
        # features this stratum had in source but dropped at QC
        qc_rows = [pos[k] for k in getattr(p, "qc_dropped_keys", set()) if k in pos]
        if qc_rows:
            block[np.array(qc_rows)] = MISSING_FAILED_QC
        reason[:, col:col + w] = block

        sf = p.samples.copy()
        sf["_col"] = np.arange(col, col + w)
        samp_frames.append(sf)
        feat_frames.append(p.features.assign(_row=rows))
        col += w

    features = (pd.concat(feat_frames, ignore_index=True)
                .drop_duplicates("_row").set_index("_row")
                .reindex(range(n_f)).reset_index(drop=True))
    features["feature_key"] = keys
    samples = (pd.concat(samp_frames, ignore_index=True)
               .sort_values("_col").drop(columns=["_col"]).reset_index(drop=True))
    return {"features": features, "samples": samples,
            "native": native, "z": zed, "lod": lod, "reason": reason,
            "policy": policy,
            "numeric_space": verify_common_numeric_space(parts, policy)}
