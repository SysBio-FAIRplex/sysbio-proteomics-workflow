"""Assertion audit - "trust no one", including this codebase.

Every claim encoded in config.py, transcribed from a README, or asserted by the
tracking file is re-checked against the data itself. Disagreements are reported, not
silently resolved. A claim that cannot be checked from the data is reported as
UNVERIFIABLE rather than quietly passing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C

PASS, FAIL, UNVERIFIABLE = "PASS", "FAIL", "UNVERIFIABLE"


def _r(check, subject, claim, source, verdict, observed=""):
    return {"check": check, "subject": subject, "claim": claim,
            "claim_source": source, "verdict": verdict, "observed": observed}


# Plausible value ranges for each declared scale, used to test the scale label against
# the numbers actually present.
SCALE_EXPECT = {
    "log2_ratio_total_batch_residual": (-6, 6, "log2 ratio centred near 0"),
    "log2_rel_batch_median": (-8, 8, "log2 ratio centred near 0"),
    "log2_reporter_intensity": (5, 40, "log2 of a linear intensity"),
    "log2_lfq_intensity": (5, 45, "log2 of a linear intensity"),
    "log2_ratio": (-8, 8, "log2 ratio centred near 0"),
    "log10_rfu_anml": (0, 8, "log10 RFU"),
    "npx_log2": (-15, 20, "Olink NPX, log2"),
    "log2_intensity_batch_corrected": (2, 32, "log2 of a batch-corrected linear intensity"),
}


def verify_stratum(h) -> list[dict]:
    """Audit one harmonised stratum against what config.py claims about it."""
    s, v, out = h.stratum, h.native, []

    # 1. declared scale vs observed range
    lo, hi, desc = SCALE_EXPECT.get(s.value_scale, (None, None, None))
    p1, p99 = np.nanpercentile(v, [1, 99])
    if lo is None:
        out.append(_r("scale_range", s.key, s.value_scale, "config", UNVERIFIABLE,
                      "no expectation registered for this scale"))
    else:
        ok = (p1 >= lo) and (p99 <= hi)
        out.append(_r("scale_range", s.key, f"{s.value_scale} ({desc})", "config",
                      PASS if ok else FAIL, f"p1={p1:.2f} p99={p99:.2f}, expected [{lo},{hi}]"))

    # 2. centred-at-zero claim for anything calling itself a ratio
    if "ratio" in s.value_scale and "intensity" not in s.value_scale:
        med = float(np.nanmedian(v))
        ok = abs(med) < 0.5
        out.append(_r("ratio_centred", s.key, "ratio scale should centre near 0",
                      "config", PASS if ok else FAIL, f"median={med:.4f}"))

    # 3. tissue claim vs whatever the biospecimen metadata says
    site_obs = h.samples.get("anatomic_site_source")
    if site_obs is not None and site_obs.notna().any():
        obs = set(str(x).lower() for x in site_obs.dropna().unique())
        ok = any(s.anatomic_site.lower() in o or o in s.anatomic_site.lower() for o in obs)
        out.append(_r("anatomic_site", s.key, s.anatomic_site, "config",
                      PASS if ok else FAIL, f"metadata says {sorted(obs)}"))
    else:
        out.append(_r("anatomic_site", s.key, s.anatomic_site, "config", UNVERIFIABLE,
                      "no anatomic site in the source annotation"))

    # 4. person resolution actually achieved.
    #
    # Counted on the NON-POOL rows specifically, rather than subtracting the pool count
    # from a whole-stratum resolved total. That shortcut assumed a pool never carries a
    # donor id, which stopped being true once `is_pool` became evidence-based (WP-8): the
    # MSBB TMT assay metadata flags 4 GIS channels whose biospecimen row still names an
    # individual, and the old arithmetic reported "-4 unresolved" for them.
    per = h.samples.get("person_source_value")
    pool = h.samples.get("is_pool")
    idx = h.samples.index
    pool_mask = (pool.fillna(False).astype(bool) if pool is not None
                 else pd.Series(False, index=idx))
    res_mask = (per.notna() if per is not None else pd.Series(False, index=idx))
    n_nonpool = int((~pool_mask).sum())
    nres = int((res_mask & ~pool_mask).sum())
    unres = n_nonpool - nres
    out.append(_r("person_resolution", s.key, "person_resolvable=%s" % s.person_resolvable,
                  "config", PASS if unres == 0 else FAIL,
                  f"{nres}/{n_nonpool} non-pool samples resolved, {unres} unresolved "
                  f"({int(pool_mask.sum())} pool/control rows excluded)"))

    # 5. pipeline / platform evidence grade. `platform` is read per sample from the assay
    #    metadata where one exists (WP-8), so its claim source is the deposit rather than
    #    config for six strata -- and where the object was read but says nothing, the
    #    audit reports a measured silence instead of "no evidence available on disk",
    #    which used to be a statement about the build rather than about the deposit.
    from .harmonize import observed_platform
    plat, plat_ev = observed_platform(h)
    meta = (h.notes or {}).get("assay_meta") or {}
    for field, val, ev, src in [
            ("analysis_pipeline", s.analysis_pipeline, s.pipeline_evidence, "config"),
            ("platform", plat, plat_ev, "assay metadata" if meta.get("joined") else "config")]:
        if val is None:
            observed = "no evidence available on disk"
            if field == "platform" and meta and not meta.get("joined"):
                observed = (f"{meta.get('reason')}; the object holds "
                            f"{meta.get('platform_populated_in_source', 0)} populated "
                            f"platform values across {meta.get('n_rows', 0)} rows")
            out.append(_r(f"cde_{field}", s.key, "not populated", "config",
                          UNVERIFIABLE, observed))
        elif ev == "B":
            out.append(_r(f"cde_{field}", s.key, val, "bundled README", UNVERIFIABLE,
                          "transcribed from documentation, not derivable from data"))
        else:
            detail = f"evidence={ev}"
            if field == "platform" and meta.get("joined"):
                detail = (f"evidence={ev}, read per sample from {meta['source']} "
                          f"on {meta['on']} ({meta['n_samples_matched']}/"
                          f"{meta['n_samples']} matched)")
            out.append(_r(f"cde_{field}", s.key, val, src, PASS, detail))

    out += verify_cde_conflict(h)
    out += verify_brodmann(h)
    out += verify_pmi(h)

    # 6. feature identity
    nf = len(h.features)
    nuniq = h.features.feature_id.nunique()
    out.append(_r("feature_ids", s.key, "every feature resolves to a UniProt accession",
                  "contract", PASS if h.features.feature_id.notna().all() else FAIL,
                  f"{nf} features, {nuniq} distinct accessions"))
    return out


def verify_cde_conflict(h) -> list[dict]:
    """Config's declared CDEs against what the assay metadata says (WP-8 step 2).

    Reported, never resolved (FR-7 applied to CDEs). This is where N3 / D-C surfaces, and
    reading the metadata per sample makes the finding sharper than the decision that
    prompted it: D-C recorded syn21323404 as saying SPS-MS3 / CID-then-HCD / 60000 for
    ROSMAP round 1, and 40 of its 400 samples do -- the other 360 read HCD / MS2 / 30000,
    which agrees with the acquisition descriptors sampled from the spectrum file. So the
    conflict is real but local, and the row says so rather than repeating the blanket
    claim.
    """
    s, meta, out = h.stratum, (h.notes or {}).get("assay_meta") or {}, []
    if not meta.get("joined"):
        return out

    obs = meta.get("platform_observed") or {}
    if obs and s.platform:
        declared = str(s.platform).strip().lower()
        agree = any(declared in k.lower() or k.lower() in declared for k in obs)
        out.append(_r("cde_conflict", s.key,
                      f"config declares platform = {s.platform!r}", "config",
                      PASS if agree else FAIL,
                      f"assay metadata says {obs}"
                      + ("" if agree else " - reported, not resolved")))

    # Acquisition descriptors sampled from the spectrum files, against the deposit's own.
    acq = meta.get("acquisition_observed") or {}
    declared_acq = {k: s.extra.get(k) for k in ("activation_type", "ms_order")
                    if s.extra.get(k)}
    if acq and declared_acq:
        detail = "; ".join(f"{k}={v}" for k, v in sorted(acq.items()))
        stated = "; ".join(f"{k}={v}" for k, v in sorted(declared_acq.items()))
        out.append(_r("cde_conflict", s.key,
                      f"spectrum-file descriptors say {stated}",
                      "sampled spectrum file", UNVERIFIABLE,
                      f"assay metadata says {detail} - the two disagree for part of the "
                      "stratum; the metadata is trusted per D-C and the disagreement is "
                      "recorded rather than dropped"))
    return out


def verify_brodmann(h) -> list[dict]:
    """The convention fill, reported as an assertion rather than a measurement (D-B).

    UNVERIFIABLE by construction where the fill is by convention: nothing in the data can
    confirm it. The row exists so a consumer -- particularly the CDM loader, which takes
    `brodmann_area` over the free-text site label -- can see how many rows are inference
    and refuse to move them onto a BA concept on our say-so.
    """
    b = (h.notes or {}).get("brodmann") or {}
    if not b:
        return []
    n_inferred = b.get("n_convention", 0)
    detail = (f"{b.get('n_source', 0)} source, {b.get('n_propagated', 0)} propagated, "
              f"{n_inferred} by convention"
              + (f" ({b['convention']} from {b['anatomic_site']!r})"
                 if b.get("convention") else "")
              + f", {b.get('n_unfilled', 0)} unfilled")
    if b.get("n_source_disagrees_with_convention"):
        detail += (f"; {b['n_source_disagrees_with_convention']} rows carry a stated or "
                   f"propagated value that disagrees with the convention and were NOT "
                   f"overwritten ({b.get('disagreeing_values')})")
    return [_r("brodmann_area", h.stratum.key,
               "convention fills only where the source is silent", "config (D-B)",
               PASS if n_inferred == 0 else UNVERIFIABLE, detail)]


def verify_pmi(h) -> list[dict]:
    """Post-mortem interval, per stratum (WP-7 / D-J).

    Two checks in one row set: that the evidence state is one the contract declares, and
    that the post-conversion median lands inside the plausible band. Out-of-band REPORTS
    -- it never auto-corrects, which is FR-7's rule applied to units. A wrong unit is the
    failure mode this exists to catch: MSBB records minutes by convention, so a stratum
    silently read as hours would be off by 60x on a covariate people model with.
    """
    p = (h.notes or {}).get("pmi") or {}
    if not p:
        return []
    ev, key = p.get("evidence"), h.stratum.key
    if ev == "not_applicable_antemortem":
        return [_r("pmi", key, "no post-mortem interval exists for this specimen type",
                   "config (D-J)", PASS, p.get("note") or "antemortem fluid")]
    if ev == "source_not_acquired" or not p.get("n_populated"):
        return [_r("pmi", key, "post-mortem interval carried where the source has it",
                   "config (D-J)", UNVERIFIABLE,
                   f"not acquired: {p.get('note') or 'source not held'}")]

    lo, hi = p.get("plausible_band_hours", C.PMI_PLAUSIBLE_HOURS)
    med, n_out = p.get("median_hours"), p.get("n_outside_plausible_band", 0)
    in_band = med is not None and lo <= med <= hi
    return [_r("pmi", key,
               f"pmi_hours in hours via {ev} (factor {p.get('factor_to_hours')})",
               "config (D-J)", PASS if in_band else FAIL,
               f"{p['n_populated']}/{p['n_samples']} populated from "
               f"{p.get('source')}:{p.get('source_column')} as {p.get('source_unit')}; "
               f"median {med} h, range {p.get('min_hours')}-{p.get('max_hours')} h; "
               f"{n_out} value(s) outside the plausible band [{lo}, {hi}] h "
               f"- reported, not corrected")]


def verify_undeclared_columns(undeclared: list[str], n_checked: int) -> list[dict]:
    """SEC-3's other half: the allow-list, not merely the exemption set (WP-7 item 4).

    `DENY_EXEMPT` says which columns survive the deny rule; it does not say which columns
    are *supposed* to exist. Without this, a new column can appear in the artifact and in
    every downstream consumer without anyone having declared it.

    **Reported, never enforced.** An undeclared column is a documentation gap, not a scope
    breach, and stopping a build over one would punish the wrong thing (FR-7's rule). The
    deny gate is what fails a build; this is what keeps the contract honest.
    """
    return [_r("undeclared_columns", "column union",
               "every emitted column is declared in the contract surface", "SEC-3",
               PASS if not undeclared else UNVERIFIABLE,
               f"{n_checked} columns checked, all declared" if not undeclared else
               f"{len(undeclared)} of {n_checked} emitted columns are outside the "
               f"declared surface: {undeclared} - reported, not enforced; declare them "
               "in `long_cols` or `EXTRA_DECLARED_COLS`")]


def verify_bundles(bundles: dict[str, dict], layer_reports: list[dict]) -> list[dict]:
    """WP-13: the bundle must exist and must agree with the HDF5 it was cut from.

    Checked as shape rather than as a promise, because "the bundle is a projection of the
    same in-memory object" is exactly the kind of claim that stays true in the docstring
    long after it stops being true in the code.
    """
    out = []
    by_layer = {r["layer"]: r for r in layer_reports}
    for layer, b in sorted(bundles.items()):
        rep, m = by_layer.get(layer, {}), b.get("matrix", {})
        ok = (m.get("n_features") == rep.get("n_proteins")
              and m.get("n_specimens") == rep.get("n_rows")
              and b.get("samples", {}).get("n_rows") == rep.get("n_rows")
              and b.get("features", {}).get("n_rows") == rep.get("n_proteins"))
        out.append(_r("tabs_bundle", layer,
                      "wide bundle matches the HDF5 layer it was cut from", "contract",
                      PASS if ok else FAIL,
                      f"matrix {m.get('n_features')}f x {m.get('n_specimens')}s on "
                      f"{m.get('value_surface')} keyed by {m.get('column_basis')}; "
                      f"samples {b.get('samples', {}).get('n_rows')}, features "
                      f"{b.get('features', {}).get('n_rows')}; HDF5 has "
                      f"{rep.get('n_proteins')}f x {rep.get('n_rows')}s"))
    return out


def verify_returns(returns: list[dict]) -> list[dict]:
    """WP-13 step 4: the return path, demonstrated rather than designed.

    The boundary is the whole point of section 7.3 -- returning one AMP programme's
    donors to another is not recoverable after the fact -- so it is computed on every
    build, not at release time.
    """
    return [_r("return_boundary", r["grant"],
               "a grant's return contains no other grant's specimens", "computed",
               PASS if r["n_rows_from_another_grant"] == 0 else FAIL,
               f"{r['n_specimens_returned']} specimens across "
               f"{len(r['layers'])} layer(s); {r['n_rows_from_another_grant']} rows from "
               f"another grant (grants in build: {r['grants_present_in_build']})")
            for r in returns]


def verify_layer(layer: str, merged: dict, parts) -> list[dict]:
    out = []
    pol = merged["policy"]

    # Z decision must follow the measurement, not the label
    dc = pol["distribution_check"]
    out.append(_r("z_decision", layer,
                  f"z_applied={pol['z_applied']}", "measured distributions", PASS,
                  pol["rationale"]))
    if pol.get("label_vs_data"):
        out.append(_r("label_vs_data", layer, "declared scale agrees with the data",
                      "config", FAIL, pol["label_vs_data"]))

    # If Z was applied, every stratum should now be centred at 0 with unit spread
    if pol["z_applied"]:
        for p in parts:
            if p.z is None:
                continue
            med = float(np.nanmedian(p.z))
            sd = float(np.nanmedian(np.nanstd(p.z, axis=1, ddof=1)))
            ok = abs(med) < 0.2 and abs(sd - 1) < 0.2
            out.append(_r("z_result", p.stratum.key, "Z gives centre~0 spread~1",
                          "computed", PASS if ok else FAIL,
                          f"median={med:.3f} spread={sd:.3f}"))
        # And a normality gate should have left nothing badly skewed
        for p in parts:
            sc = p.notes.get("scale_check_after_transform") or p.notes["scale_check"]
            ok = sc["normality_ok"] is not False
            out.append(_r("normality_before_z", p.stratum.key,
                          "approximately normal before Z", "computed",
                          PASS if ok else FAIL,
                          f"median|skew|={sc['median_abs_skew']} "
                          f"transform={sc['transform_applied']}"))
    else:
        for p in parts:
            out.append(_r("native_preserved", p.stratum.key,
                          "no transform or Z beyond what the source needed", "policy",
                          PASS, f"transform={p.notes['transform_applied']}, z=False"))

    # Post-condition, not an assumption: the strata that are about to be merged must
    # actually occupy one numeric space on the surface that ships.
    space = merged.get("numeric_space", {})
    out.append(_r("common_numeric_space", layer,
                  "all strata share one numeric space at merge", "post-condition",
                  PASS if space.get("verified") else FAIL,
                  f"{space.get('reason')} (surface={space.get('merged_surface')}, "
                  f"loc={space.get('location_range_in_pooled_sd')}, "
                  f"spread_ratio={space.get('spread_ratio_max_over_min')})"))

    # Row key must be unique
    from .writer import build_row_key
    rk = build_row_key(merged["samples"])
    out.append(_r("row_key_unique", layer, "one row per person+visit", "contract",
                  PASS if rk.is_unique else FAIL,
                  f"{len(rk)} rows, {rk.nunique()} distinct"))

    # Wide table shape must match the feature and sample tables
    ok = (merged["native"].shape[0] == len(merged["features"])
          and merged["native"].shape[1] == len(merged["samples"]))
    out.append(_r("shape_consistency", layer, "matrix matches feature/sample tables",
                  "contract", PASS if ok else FAIL,
                  f"{merged['native'].shape} vs {len(merged['features'])}f "
                  f"x {len(merged['samples'])}s"))
    return out


def verify_file_ids(layer: str, samples: pd.DataFrame, files_df: pd.DataFrame) -> list[dict]:
    """FR-11: every `samples.file_id` must name a row that exists in the FILES table.

    Previously the column held a composed role placeholder (`<stratum>:protein_abundance
    _matrix`) that matched nothing, so the link was unresolvable from the artifact alone.
    A null is legitimate here -- it means the stratum has several matrix files and the
    `sample_files` bridge carries the relation -- so nulls are reported, not failed.
    """
    known = set(files_df["file_id"].dropna().astype(str)) if len(files_df) else set()
    fid = samples.get("file_id")
    if fid is None:
        return [_r("file_id_resolves", layer, "samples.file_id resolves into FILES",
                   "contract", FAIL, "column absent")]
    present = fid.dropna().astype(str)
    unresolved = sorted(set(present) - known)
    n_null = int(fid.isna().sum())
    return [_r("file_id_resolves", layer, "samples.file_id resolves into FILES",
               "contract", PASS if not unresolved else FAIL,
               f"{len(present)}/{len(fid)} populated and all resolving"
               if not unresolved else
               f"{len(unresolved)} unresolved id(s): {unresolved[:3]}"),
            _r("file_id_coverage", layer,
               "null file_id means a one-to-many stratum covered by sample_files",
               "contract", PASS if n_null == 0 else UNVERIFIABLE,
               f"{n_null} rows null (multi-file strata); bridge carries the relation"
               if n_null else "every row resolves 1:1")]


def verify_column_consistency(layer_samples: dict[str, list[str]],
                              layer_features: dict[str, list[str]]) -> list[dict]:
    """Every layer must expose the same CDE surface under the same names."""
    from . import config as C
    out = []
    required = set(C.LINK_KEYS) | set(C.ASSAY_CDES) | set(C.FILES_CDES)

    for layer, cols in sorted(layer_samples.items()):
        missing = sorted(required - set(cols))
        out.append(_r("cde_columns_present", layer,
                      f"all {len(required)} CDE/link columns present", "contract",
                      PASS if not missing else FAIL,
                      "complete" if not missing else f"missing {missing}"))

    # Names must not drift. The comparison runs against the UNION across layers, not the
    # intersection: a column a layer is *missing* is the contract violation, and the old
    # `set(layer) - common` direction could only ever see columns held in excess. Checked
    # against the shipped artifact, that inversion reported PASS on 11 of 12 rows it got
    # wrong -- certifying the three layers missing six columns each as "identical to all
    # layers" while reporting the one conforming layer as the deviant one (VER-7).
    all_layers = sorted(layer_samples)
    common = set.intersection(*(set(layer_samples[l]) for l in all_layers))
    union = set().union(*(set(layer_samples[l]) for l in all_layers))
    for layer in all_layers:
        absent = sorted(union - set(layer_samples[layer]))
        extra = sorted(set(layer_samples[layer]) - common)
        detail = []
        if absent:
            detail.append(f"missing {len(absent)}: {absent}")
        if extra:
            detail.append(f"{len(extra)} not in every layer: {extra}")
        out.append(_r("column_naming", layer, "column set matches the other layers",
                      "contract", PASS if not absent else FAIL,
                      "identical to all layers" if not detail else "; ".join(detail)))

    # Same treatment for the feature tables, which previously asserted PASS without
    # performing any comparison at all. One row per layer rather than one global row, so
    # the audit names which layer is short and by what.
    funion = set().union(*(set(layer_features[l]) for l in all_layers))
    for layer in all_layers:
        absent = sorted(funion - set(layer_features[layer]))
        out.append(_r("feature_columns", layer,
                      f"feature table carries all {len(funion)} columns", "contract",
                      PASS if not absent else FAIL,
                      "complete" if not absent else f"missing {len(absent)}: {absent}"))
    return out


# Traces a `SOURCE_QC` step leaves in `notes` when it actually runs. A step declared
# ours to apply (`applied_upstream=False`) that leaves no trace has not been applied, and
# is reported as a FAIL rather than being taken on the declaration's word -- VER-1's rule
# moved from claims about data to claims about processing.
_QC_TRACE = {
    "reverse / contaminant": "n_decoy_contaminant_dropped",
    "colcheck": "colcheck_counts",
    "qc_warning": "n_samples_with_any_qc_warn",
    "below-lod fraction": "n_assays_gt25pct_below_lod",
    "log2 of the reference ratio": "transform_applied",
    "per-protein centring": "centring",
    "features <50% missing": "n_features_below_50pct_present_dropped",
    "control wells flagged": "control_wells_flagged",
}


def verify_source_qc(h) -> list[dict]:
    """Every `SOURCE_QC` entry, checked instead of merely recited (WP-12 item 1, VER-9).

    The ~30 entries used to reach `notes["source_qc"]` and the generated README and stop
    there: nothing confirmed that a producer's `applied_upstream=True` claim held in the
    delivered file, and nothing confirmed that an `applied_upstream=False` step -- the
    ones the producer explicitly left to the consumer -- had actually been performed.
    Three were in fact unapplied and unflagged, which was invisible precisely because one
    of them (the Olink below-LOD flags) *was* computed and made the section look live.
    """
    s, notes, out = h.stratum, h.notes, []
    for q in notes.get("source_qc", []):
        name, upstream, ref = q["check"], q["applied_upstream"], q["reference"]
        low = name.lower()

        if upstream:
            # Verified where the delivered data can actually show it.
            if "50%" in low or "50 %" in low:
                n = notes.get("n_features_below_50pct_present_dropped", 0)
                out.append(_r("source_qc", s.key, f"upstream: {name}", ref,
                              PASS if not n else FAIL,
                              "delivered file satisfies the threshold" if not n else
                              f"{n} features in the delivered file are below 50% present, "
                              "so the claimed upstream filter does not hold on what was "
                              "shipped; re-applied here"))
            elif "gis" in low or "reference channel" in low:
                pool = h.samples.get("is_pool")
                n = int(pool.fillna(False).astype(bool).sum()) if pool is not None else 0
                out.append(_r("source_qc", s.key, f"upstream: {name}", ref,
                              PASS if n == 0 else FAIL,
                              f"{n} pool/control channel(s) present in the delivered "
                              f"matrix" if n else "no pool channels survive into the "
                              "delivered matrix"))
            else:
                out.append(_r("source_qc", s.key, f"upstream: {name}", ref, UNVERIFIABLE,
                              "the producer states this was applied before deposit and "
                              "the delivered file carries no trace either way"))
            continue

        key = next((v for k, v in _QC_TRACE.items() if k in low), None)
        trace = notes.get(key) if key else None
        done = trace is not None and trace != {} and trace is not False
        out.append(_r("source_qc", s.key, f"ours to apply: {name}", ref,
                      PASS if done else FAIL,
                      f"applied; trace `{key}` = "
                      f"{str(trace)[:90]}" if done else
                      "declared ours to apply but the build leaves no trace of having "
                      "done it"))
    return out


def verify_qc_metadata(h) -> list[dict]:
    """The QC signals the deposits ship that this build previously read and discarded.

    Each is reported, never used to delete anything: flagging is the standing policy and
    the exclusion threshold belongs to the study, not to the harmoniser (D8).
    """
    s, notes, out = h.stratum, h.notes, []

    cc = notes.get("colcheck_counts")
    if cc is not None:
        flagged = sum(v for k, v in cc.items() if str(k).upper() != "PASS")
        out.append(_r("colcheck", s.key,
                      "SomaLogic per-analyte column QC is carried, not applied",
                      "SomaScan protein metadata", PASS,
                      f"{cc}; {flagged} analyte(s) flagged and shipped with the flag"))
    if "sample_level_qc_in_deposit" in notes and notes["sample_level_qc_in_deposit"] is None:
        out.append(_r("sample_qc_available", s.key,
                      "the deposit provides sample-level QC", "deposit", UNVERIFIABLE,
                      "it does not - the sample metadata is clinical only and carries no "
                      "QC column, so no sample-level verdict can be derived from it"))

    ex = notes.get("adkp_exclusions")
    if ex is not None:
        out.append(_r("source_exclusion", s.key,
                      "portal `exclude` flags are wired to qc_status", ex["source"],
                      PASS,
                      f"{ex['n_flagged']} of {ex['n_specimens']} specimens flagged"
                      + (f": {ex['reasons']}" if ex["reasons"] else
                         " - none today, but the flag is now live rather than discarded")))

    cw = notes.get("control_wells_flagged")
    if cw is not None:
        out.append(_r("control_wells", s.key,
                      "control wells are flagged from the source columns",
                      "isControl / sample.type", PASS if cw["n_flagged"] else FAIL,
                      f"{cw['n_flagged']} of {cw['n_samples']} flagged from "
                      f"{cw['columns_used']} (isControl {cw['from_isControl']}, "
                      f"sample.type != subjects {cw['from_sample_type_not_subjects']})"))

    wc = notes.get("n_samples_by_warn_count")
    if wc is not None:
        out.append(_r("qc_verdicts_separable", s.key,
                      "the three vendor QC verdicts ship separately, not folded into one",
                      "contract", PASS,
                      f"samples by number of failing checks: {wc}; "
                      f"{notes.get('n_samples_with_any_qc_warn', 0)} sample(s) trip at "
                      "least one. `qc_status` remains the summary; `qc_detail` and "
                      "`qc_warn_count` carry the detail"))

    # Vendor missingness against ours, reported side by side rather than reconciled.
    if "missing_freq" in getattr(h, "features", pd.DataFrame()).columns \
            and "frac_below_lod" in h.features.columns:
        mf = pd.to_numeric(h.features["missing_freq"], errors="coerce")
        fb = pd.to_numeric(h.features["frac_below_lod"], errors="coerce")
        both = mf.notna() & fb.notna()
        if both.any():
            d = (mf[both] - fb[both]).abs()
            corr = float(mf[both].corr(fb[both])) if both.sum() > 2 else float("nan")
            out.append(_r("missingness_cross_check", s.key,
                          "vendor MissingFreq agrees with the computed frac_below_lod",
                          "computed", PASS if d.median() <= 0.05 else FAIL,
                          f"n={int(both.sum())}, median |difference| {d.median():.4f}, "
                          f"max {d.max():.4f}, correlation {corr:.3f} - both ship; "
                          "neither is preferred over the other"))
    return out


def verify_sparsity(h) -> list[dict]:
    """No stratum should ship extreme sparsity after filtering."""
    sp = h.notes.get("sparsity", {})
    pct = sp.get("pct_missing_overall")
    if pct is None:
        return []
    return [_r("sparsity", h.stratum.key, "post-filter missingness stays moderate",
               "computed", PASS if pct <= 35 else FAIL,
               f"{pct}% missing overall, median feature present "
               f"{sp.get('median_feature_presence_pct')}%, "
               f"min {sp.get('min_feature_presence_pct')}%")]


def verify_exclusions(h) -> list[dict]:
    """Rows removed from a stratum, reported as data rather than falling out silently.

    Deletion is reserved for the two cases where a row's IDENTITY is unusable, not its
    quality: the Olink plate swap (D-D) and a sample the source cannot resolve to a
    donor. The DIA completeness drop is the one quality-based exclusion, and it exists
    because that deposit ships no sample-level QC of its own. Everything else is flagged
    and shipped for the consumer to filter (D8).
    """
    out, n = [], h.stratum.key

    w = h.notes.get("withdrawn_samples")
    if w:
        out.append(_r("samples_withdrawn", n,
                      "samples withdrawn for disputed identity", "computed", FAIL,
                      f"{w['n_samples']} of {w['n_samples_total']} samples "
                      f"({w['n_samples']/max(w['n_samples_total'],1):.1%}) withdrawn; "
                      f"{w['n_participants_losing_every_sample']} participants lose every "
                      f"sample in this stratum. Rule: {w['rule']}. {w['reason']}"))

    # SEC-8 throughout this function: publish the rule and the counts, never a real key.
    # Both of these exclusions are keyed by identifiers that embed a participant -- a
    # SomaScan `projid_visit`, a DIA `PP-<participant>-<visit>-...` -- and these rows go
    # straight into `assertion_audit.tsv` and the generated README, which are publishable.
    u = h.notes.get("unresolved_donors_dropped")
    if u:
        out.append(_r("unresolved_donor_dropped", n,
                      "a sample the source cannot resolve to a donor is removed",
                      "computed", FAIL,
                      f"{u['n_dropped']} sample(s) dropped. {u['reason']}. Rule: "
                      f"{u['rule']}. Specimen ids are withheld from this artifact "
                      "(SEC-8) - they embed a participant identifier, and the rule "
                      "reproduces the set exactly"))

    n_drop = h.notes.get("n_samples_dropped_low_completeness")
    if n_drop is not None:
        thr = h.notes.get("qc_min_sample_completeness")
        counts = h.notes.get("detected_counts_of_dropped_samples") or []
        n_feat = h.notes.get("n_features_total")
        out.append(_r("qc_sample_completeness", n,
                      f"samples detecting < {thr:.0%} of the layer's features are dropped",
                      "config", PASS if n_drop == 0 else FAIL,
                      f"{n_drop} sample(s) dropped"
                      + (f"; they detected {counts[:5]} of {n_feat} features "
                         "(lowest first) - sample ids withheld, they embed a participant "
                         "identifier (SEC-8)" if counts else
                         "; none below threshold")))

    n_feat = h.notes.get("n_features_dropped_empty")
    if n_feat:
        out.append(_r("features_dropped_empty", n,
                      "a feature detected in no retained sample carries no information",
                      "computed", FAIL,
                      f"{n_feat} feature(s) dropped: "
                      f"{h.notes.get('features_dropped_empty', [])[:5]}"))

    # The zero-sentinel finding, recorded so the artifact carries it (WP-10).
    if h.notes.get("matrix_sibling_not_read"):
        out.append(_r("source_form", n, "read the long form, not the dense matrix",
                      "computed", PASS, h.notes["matrix_sibling_not_read"]))
    return out


def verify_tracking_file(diff: pd.DataFrame) -> list[dict]:
    """The tracking file's own assertions, re-checked -- generated from the computed
    diff rather than restated beside it.

    Every row here used to be a literal. One of them ("DiverseCohorts temporal ->
    TCX = temporal cortex (confirmed)") published a wrong FAIL for three weeks because a
    hardcoded verdict cannot change when the data does; the tracking file had been right
    all along (N2). Deriving these from `tracking_file_diff` is the fix for the class,
    not just for that row.
    """
    out = []
    for r in diff.itertuples():
        if r.agrees is None or (isinstance(r.agrees, float) and pd.isna(r.agrees)):
            verdict, observed = UNVERIFIABLE, (r.note or "not computable from the build")
        elif r.agrees:
            verdict, observed = PASS, f"derived {r.derived_value}"
        else:
            verdict, observed = FAIL, f"derived {r.derived_value} ({r.note})"
        out.append(_r("tracking", r.tracking_row, f"{r.field} = {r.tracking_value}",
                      "tracking file", verdict, observed))
    return out


def verify_readme_claims(layer_reports: list[dict], probes: dict) -> list[dict]:
    """Claims taken from bundled READMEs, which are documentation and not data.

    A README claim is only ever PASS here when the artifact itself can show it. Anything
    that needs a fact we do not hold is UNVERIFIABLE with the reason (VER-2) -- never a
    quiet pass, and never a verdict typed in by hand.
    """
    out = []
    by_layer = {r["layer"]: r for r in layer_reports}

    # Olink panels partition the UniProt space -- checkable from the built feature table.
    for layer in ("olink_plasma_proteomics", "olink_csf_proteomics"):
        rep = by_layer.get(layer)
        if rep is None:
            continue
        out.append(_r("readme", layer,
                      "duplicate UniProts across panels were averaged", "AMP-PD README",
                      UNVERIFIABLE,
                      f"{rep['n_proteins']} distinct accessions after merge; whether any "
                      "were averaged cannot be recovered from the delivered matrices"))

    # The SomaScan scale label. N1: the ANML and non-ANML objects are byte-identical, so
    # the normalisation label has no evidence behind it either way.
    if "soma_plasma_proteomics" in by_layer:
        out.append(_r("readme", "soma_plasma_proteomics", "ANML normalised, log10",
                      "deposit README", UNVERIFIABLE,
                      "the ANML and non-ANML deposits share md5 "
                      "f8a3d03ecd1b57cbb409bc6b054a7890, so one of the two is "
                      "mislabelled and the data cannot say which (N1, D-E); the log10 "
                      "half of the claim is consistent with the observed range"))

    # CSF has no _retracted variant, so nothing on disk can confirm it was unaffected.
    if "olink_csf_proteomics" in by_layer:
        out.append(_r("readme", "olink_csf_proteomics", "no correction issued for CSF",
                      "inferred", UNVERIFIABLE,
                      "CSF ships no _retracted variant; absence of a correction file is "
                      "not evidence that CSF was unaffected"))
    return out


def summarise(rows: list[dict]) -> dict:
    df = pd.DataFrame(rows)
    return {
        "n_checks": len(df),
        "n_pass": int((df.verdict == PASS).sum()),
        "n_fail": int((df.verdict == FAIL).sum()),
        "n_unverifiable": int((df.verdict == UNVERIFIABLE).sum()),
    }
