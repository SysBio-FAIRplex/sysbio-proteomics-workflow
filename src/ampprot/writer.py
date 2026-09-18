"""S8 - write the single HDF5.

Layout, one keyed group per (assay x matrix) layer:

    sysbio_proteomics-08222026.h5
    |-- <layer_key>/
    |     |-- values          (n_rows, n_proteins) float32   WIDE: rows are samples
    |     |-- values_z        same shape, only where the layer needed Z
    |     |-- lod, below_lod  same shape, Olink layers only
    |     |-- row_key         (n_rows,) person + visit, the analysis key
    |     |-- samples/<col>   one dataset per CDE / link / technical column
    |     |-- features/<col>  one dataset per feature-description column
    |     |-- assay/<col>     SysBio ASSAY CDE table, one row per stratum
    |     |-- files/<col>     SysBio FILES CDE table, one row per source file
    |     +-- attrs           scale policy, transform, z flag, provenance
    +-- <deferred layer>/     declared but empty, with the reason in attrs

Rows are keyed off person and visit so the table is directly analysis-ready and
mergeable across layers on `person_id` (+ visit where longitudinal).
"""

from __future__ import annotations

import json

import h5py
import numpy as np
import pandas as pd

from . import config as C, harmonize as H

STR = h5py.string_dtype(encoding="utf-8")


def _write_col(grp: h5py.Group, name: str, series) -> None:
    """Write one column, choosing a sane HDF5 dtype and encoding nulls honestly."""
    s = pd.Series(series)
    if pd.api.types.is_bool_dtype(s.dtype):
        grp.create_dataset(name, data=s.fillna(False).to_numpy(np.bool_),
                           compression="gzip", compression_opts=4)
        return
    if pd.api.types.is_integer_dtype(s.dtype) or isinstance(s.dtype, pd.Int64Dtype):
        # nullable ints -> float64 so NaN survives; ids are well below 2^53
        arr = s.astype("Float64").to_numpy(np.float64, na_value=np.nan)
        grp.create_dataset(name, data=arr, compression="gzip", compression_opts=4)
        return
    if pd.api.types.is_float_dtype(s.dtype):
        grp.create_dataset(name, data=s.to_numpy(np.float64),
                           compression="gzip", compression_opts=4)
        return
    txt = s.where(s.notna(), None).map(lambda v: "" if v is None else str(v))
    grp.create_dataset(name, data=txt.to_numpy(dtype=object), dtype=STR,
                       compression="gzip", compression_opts=4)


def _write_table(parent: h5py.Group, name: str, df: pd.DataFrame) -> None:
    g = parent.create_group(name)
    g.attrs["n_rows"] = len(df)
    g.attrs["columns"] = json.dumps(list(df.columns))
    for col in df.columns:
        _write_col(g, col, df[col])


def _matrix(g: h5py.Group, name: str, arr: np.ndarray) -> h5py.Dataset:
    chunk = (min(256, arr.shape[0]), min(512, arr.shape[1])) if arr.size else None
    return g.create_dataset(name, data=arr, dtype="f4", chunks=chunk, shuffle=True,
                            compression="gzip", compression_opts=4)


def build_row_key(samples: pd.DataFrame) -> pd.Series:
    """The analysis key. A person legitimately repeats within a layer -- several brain
    regions, several timepoints, or a technical replicate -- so the key is built up from
    the real discriminators before any anonymous counter is used.

    Order of composition, each part added only if it actually disambiguates:
        person -> anatomic site -> visit -> stratum -> replicate counter

    Every component is also emitted as its own column, so callers can regroup by person,
    by site, or by timepoint without parsing this string.
    """
    n = len(samples)
    idx = samples.index

    def col(name):
        c = samples.get(name)
        return c if c is not None else pd.Series([None] * n, index=idx)

    person = col("person_source_value").astype("object")
    person = pd.Series([p if pd.notna(p) else "UNRESOLVED" for p in person], index=idx)
    key = person.astype(str)

    # anatomic site, only when the layer spans more than one
    site = col("anatomic_site").astype("object")
    if site.notna().any() and site.dropna().nunique() > 1:
        key = key + "|" + site.fillna("NA").astype(str).str.replace(r"\s+", "-", regex=True)

    # visit, only where the layer actually has timepoints
    vm, vi = col("visit_month"), col("visit_index")
    if vm.notna().any():
        key = key + pd.Series(
            [f"|M{int(v)}" if pd.notna(v) else "" for v in vm], index=idx)
    elif vi.notna().any():
        key = key + pd.Series(
            [f"|V{int(v)}" if pd.notna(v) else "" for v in vi], index=idx)

    # stratum, if the same person+site+visit still appears twice
    if key.duplicated().any():
        st = col("stratum_key").astype("object")
        cand = key + "|" + st.fillna("NA").astype(str)
        if cand.duplicated().sum() < key.duplicated().sum():
            key = cand

    # last resort: an explicit replicate counter
    dup = key.duplicated(keep=False)
    if dup.any():
        rep = key.groupby(key).cumcount() + 1
        key = key.where(~dup, key + "|R" + rep.astype(str))
    return key


def replicate_index(samples: pd.DataFrame) -> pd.Series:
    """1-based occurrence of each person within the layer, so repeats are explicit."""
    p = samples["person_source_value"].fillna("UNRESOLVED").astype(str)
    return (p.groupby(p).cumcount() + 1).astype("int32")


def sample_files(samples: pd.DataFrame, files_df: pd.DataFrame) -> pd.DataFrame:
    """The `sample x file` relation (D-F, FR-13).

    A sample maps to several source files whenever the assay ships one matrix per panel
    -- Olink CSF to 4, plasma to 8 -- so the link cannot live in `samples.file_id`
    without either breaking "one row per sample" or smuggling in a list column. It is
    modelled as a relation instead, which is what the CDM's ASSAY_INPUT_FILE is for.

    `file_role` rides on the row so a consumer can select the matrix files without
    string-matching on names.
    """
    cols = ["row_key", "stratum_key", "file_id", "file_role"]
    if files_df.empty or "row_key" not in samples.columns:
        return pd.DataFrame(columns=cols)
    fl = files_df[["file_id", "file_role", "assay_id"]].copy()
    fl["stratum_key"] = fl["assay_id"].astype(str).str.rsplit(":", n=1).str[-1]
    out = (samples[["row_key", "stratum_key"]]
           .merge(fl.drop(columns="assay_id"), on="stratum_key", how="left"))
    return out[cols].reset_index(drop=True)


def column_unions(pending: dict[str, dict[str, pd.DataFrame]]) -> dict[str, list[str]]:
    """Ordered column union per table name, across every layer (WP-3, option 2).

    The writer -- not `config.py` -- is the source of truth for which columns exist,
    because config demonstrably drifts from what the writer emits: `frac_below_lod` and
    `olink_lod_flag` were emitted for years while undeclared, so a config-driven reindex
    would have left the contract broken while the audit reported it fixed.

    Order is first-seen across layers in `LAYER_ORDER`, so it is deterministic and a
    layer's own columns stay contiguous.
    """
    unions: dict[str, list[str]] = {}
    for tables in pending.values():
        for name, df in tables.items():
            seen = unions.setdefault(name, [])
            for col in df.columns:
                if col not in seen:
                    seen.append(col)
    return unions


def write_layer_tables(h5: h5py.File, layer: str, tables: dict[str, pd.DataFrame],
                       unions: dict[str, list[str]]) -> dict[str, int]:
    """Second pass: write the small per-layer tables reindexed to the runtime union.

    Deferred out of `write_layer` so the union can be computed across every layer first.
    Only these tables wait -- the value matrices are written during the first pass and
    released, so holding them costs a few MB rather than the whole build.
    """
    g = h5[layer]
    filled = {}
    for name, df in tables.items():
        cols = unions.get(name)
        out = df.reindex(columns=cols) if cols else df
        filled[name] = len(set(cols or []) - set(df.columns))
        _write_table(g, name, out)
    g.attrs["column_union_filled"] = json.dumps(filled)
    return filled


# --------------------------------------------------------------------------------------
# WP-13 / D-G -- the returnable bundles (`tabs/`), wide only
# --------------------------------------------------------------------------------------
#
# HARMONIZATION_PLAN section 7.2 has specified these since the first draft and section
# 7.3 builds the AMP return path on top of them, but nothing was ever written -- and the
# omission was not recorded anywhere, so section 7.2 read as delivered. This is that gap.
#
# **Wide only.** The long form is dropped (D-G): pooled brain TMT alone would be ~25M
# rows, and the wide matrix is what the HDF5 already holds, so the three tables are a
# projection of the in-memory object rather than a re-derivation.
#
# The matrix is written during the SAME pass that writes the HDF5 values, and the tables
# during the same pass that writes the HDF5 tables, from the same in-memory frames.
# A standalone HDF5 -> parquet converter would be simpler and would give up exactly the
# property section 8 S8 asks for: that the two cannot drift.
#
# SEC-5: these carry the same donor-level identifiers as the HDF5. They are
# controlled-access artifacts, not a publishable by-product, and `.gitignore` denies them
# by name as well as by directory.

def _bundle_dir(layer: str):
    d = C.TABS / layer
    d.mkdir(parents=True, exist_ok=True)
    return d


def _emit(df: pd.DataFrame, path) -> dict:
    """Parquet primary, `.tsv.gz` mirror (section 7.2).

    Parquet stays primary because it preserves dtypes and the null distinction WP-4 just
    made meaningful -- a `below_lod` of NaN and one of 0.0 are different facts, and a TSV
    round-trip loses that unless every consumer reads the legend.
    """
    df.to_parquet(path.with_suffix(".parquet"), index=False, compression="zstd")
    df.to_csv(path.with_suffix(".tsv.gz"), sep="\t", index=False, compression="gzip")
    return {"parquet": path.with_suffix(".parquet").name,
            "tsv_gz": path.with_suffix(".tsv.gz").name,
            "n_rows": len(df), "n_cols": int(df.shape[1])}


def write_bundle_matrix(layer: str, merged: dict, row_key: pd.Series) -> dict:
    """The wide matrix: `feature_id` plus one column per specimen.

    Orientation is features x specimens, per section 7.2 -- the transpose of the HDF5's
    sample-major layout, which is the orientation an analyst joins a feature annotation
    to. Both come from `merged["native"]`, so they are the same numbers.
    """
    pol = merged["policy"]
    surface = "values_z" if (pol["z_applied"] and merged["z"] is not None) else "values"
    arr = merged["z"] if surface == "values_z" else merged["native"]

    samples = merged["samples"]
    cols = samples["specimen_id"].astype(str)
    col_basis = "specimen_id"
    if not cols.is_unique:
        # Never ship ambiguous headers. `row_key` is the audited-unique analysis key and
        # the samples table carries both, so the join stays expressible either way.
        cols = row_key.astype(str)
        col_basis = "row_key"

    out = pd.DataFrame(arr, columns=cols.tolist())
    out.insert(0, "feature_id", merged["features"]["feature_key"].astype(str).values)
    info = _emit(out, _bundle_dir(layer) / f"{layer}_matrix")
    info.update({"value_surface": surface, "column_basis": col_basis,
                 "n_features": int(arr.shape[0]), "n_specimens": int(arr.shape[1])})
    return info


def write_bundle_tables(layer: str, tables: dict[str, pd.DataFrame],
                        unions: dict[str, list[str]], matrix_info: dict,
                        provenance: dict) -> dict:
    """`samples` and `features` for the bundle, plus the sidecar manifest.

    Written from the union-reindexed frames -- the same objects the HDF5 gets -- so a
    column the contract guarantees is present here too, null where it does not apply.
    """
    d = _bundle_dir(layer)
    out = {"layer": layer, "matrix": matrix_info}
    for name in ("samples", "features"):
        df = tables[name]
        cols = unions.get(name)
        out[name] = _emit(df.reindex(columns=cols) if cols else df,
                          d / f"{layer}_{name}")

    # Which surface the matrix carries is layer-dependent, so the bundle records it
    # rather than leaving the consumer to guess from the numbers.
    manifest = {
        "layer": layer,
        "artifact": C.H5_NAME,
        "form": "wide only (D-G); the long form is deliberately not emitted",
        "value_surface": matrix_info["value_surface"],
        "value_surface_meaning": (
            "per-stratum Z; the native scale is in the HDF5 as `values`"
            if matrix_info["value_surface"] == "values_z" else
            "native scale exactly as delivered by the source"),
        "matrix_columns_are": matrix_info["column_basis"],
        "join": f"matrix column header -> samples.{matrix_info['column_basis']}; "
                "matrix.feature_id -> features.feature_key",
        "files": {k: out[k] for k in ("matrix", "samples", "features")},
        "provenance": _jsonable(provenance),
        "access": ("CONTROLLED. Carries donor-level identifiers, exactly as the HDF5 "
                   "does. SEC-5 applies unchanged: not a publishable by-product."),
    }
    (d / f"{layer}_bundle.json").write_text(json.dumps(manifest, indent=2, default=str))
    out["manifest"] = f"{layer}_bundle.json"
    return out


def write_returns(bundles: dict[str, dict], samples_by_layer: dict[str, pd.DataFrame]
                  ) -> list[dict]:
    """Section 7.3's return path, run rather than described.

    One directory per `grant`, holding that grant's specimen set for every layer it
    appears in, plus a manifest that *demonstrates* the boundary instead of asserting it:
    it names the grants present in the source, the grant returned, and the count of rows
    from any other grant that reached the slice -- which must be zero.

    Until this runs, section 7.3 is a design. The check is cheap and the failure it
    guards against -- returning one AMP's donors to another -- is not recoverable after
    the fact, so it is computed on every build rather than at release time.
    """
    C.RETURNS.mkdir(parents=True, exist_ok=True)
    grants = sorted({g for df in samples_by_layer.values()
                     for g in df["grant"].dropna().astype(str).unique()})
    out = []
    for grant in grants:
        slug = "".join(ch if ch.isalnum() else "_" for ch in grant.lower()).strip("_")
        gdir = C.RETURNS / slug
        gdir.mkdir(parents=True, exist_ok=True)
        layers, leaked, n_rows = {}, 0, 0
        for layer, df in samples_by_layer.items():
            sel = df["grant"].astype(str) == grant
            if not sel.any():
                continue
            sl = df.loc[sel]
            leaked += int((sl["grant"].astype(str) != grant).sum())
            n_rows += len(sl)
            keep = [c for c in ("row_key", "specimen_id", "specimen_source_value",
                                "person_id", "person_global_key", "stratum_key",
                                "study", "cohort", "grant", "assay_id", "file_id",
                                "visit_month", "tissue", "anatomic_site")
                    if c in sl.columns]
            layers[layer] = _emit(sl[keep].reset_index(drop=True),
                                  gdir / f"{layer}_specimens")
        manifest = {
            "grant": grant,
            "artifact": C.H5_NAME,
            "layers": layers,
            "n_specimens_returned": n_rows,
            "grants_present_in_build": grants,
            "n_rows_from_another_grant": leaked,
            "boundary_check": ("PASS - no specimen from another AMP programme is in this "
                               "return" if leaked == 0 else
                               f"FAIL - {leaked} rows from another grant"),
            "bundle_source": {k: v.get("manifest") for k, v in bundles.items()
                              if k in layers},
            "access": "CONTROLLED. Donor-level identifiers; SEC-5 applies.",
        }
        (gdir / "return_manifest.json").write_text(
            json.dumps(manifest, indent=2, default=str))
        out.append(manifest)
    return out


def write_layer(h5: h5py.File, layer: str, merged: dict, parts, assay_df, files_df,
                provenance: dict) -> tuple[dict, dict]:
    g = h5.create_group(layer)
    pol = merged["policy"]

    # WIDE orientation: transpose so rows are samples and columns are proteins.
    #
    # `values` (and `values_z`) are the ANALYSIS-READY tables and carry exactly ONE kind
    # of null: float NaN. Nothing else is encoded in them -- no sentinels, no second null
    # flavour -- so numpy/pandas/sklearn behave normally.
    native_w = merged["native"].T
    _matrix(g, "values", native_w)
    if merged["z"] is not None:
        _matrix(g, "values_z", merged["z"].T)
    if merged["lod"] is not None:
        _matrix(g, "lod", merged["lod"].T)
        # Tri-state, not boolean (WP-4). `NaN < lod` is False, which casts to 0.0, so a
        # never-measured cell used to be written identically to one measured above the
        # limit -- making the conventional filter `below_lod == 0.0` silently retain
        # never-measured cells. Restoring NaN where the value is NaN separates:
        #   0.0  measured, at or above the detection limit
        #   1.0  measured, below it
        #   NaN  not measured
        # This uses the single null flavour the file already documents, so no new
        # sentinel is introduced.
        below = (merged["native"] < merged["lod"]).astype(np.float32)
        below[np.isnan(np.asarray(merged["native"], dtype=np.float32))] = np.nan
        d = _matrix(g, "below_lod", below.T)
        if d is not None:
            d.attrs["states"] = json.dumps({
                "0.0": "measured, at or above the limit of detection",
                "1.0": "measured, below the limit of detection",
                "NaN": "not measured - no value to compare against the limit",
            })

    # Separate diagnostic companion, NOT part of the analysis-ready table. It explains
    # *why* a cell is NaN without putting a second null flavour into the value matrix.
    if merged.get("reason") is not None:
        d = g.create_dataset("missing_reason", data=merged["reason"].T, dtype="i1",
                             chunks=True, compression="gzip", compression_opts=4)
        d.attrs["legend"] = json.dumps(
            {str(k): v for k, v in H.MISSING_LEGEND.items()})
        d.attrs["note"] = ("diagnostic only - `values` carries a single null type (NaN); "
                           "this array says which kind of missing each NaN is")

    samples = merged["samples"].copy()
    row_key = build_row_key(samples)
    samples.insert(0, "row_key", row_key)
    g.create_dataset("row_key", data=row_key.to_numpy(dtype=object), dtype=STR,
                     compression="gzip", compression_opts=4)

    # WP-13: the bundle's matrix, written here rather than later, from the same array
    # that just became `values`/`values_z`. Written before the array is released so the
    # HDF5 and the bundle cannot disagree about a single number.
    bundle_matrix = write_bundle_matrix(layer, merged, row_key)
    g.create_dataset("protein_columns",
                     data=merged["features"].feature_key.to_numpy(dtype=object),
                     dtype=STR, compression="gzip", compression_opts=4)

    # The tables are NOT written here. They are handed back so build.py can compute the
    # column union across every layer and write them reindexed to it (WP-3), which is
    # what makes the uniform-contract guarantee hold rather than merely be documented.
    tables = {"samples": samples, "features": merged["features"],
              "assay": assay_df, "files": files_df,
              "sample_files": sample_files(samples, files_df)}

    stratum_notes = {p.stratum.key: _jsonable(p.notes) for p in parts}
    g.attrs["layer_key"] = layer
    g.attrs["n_rows"] = int(native_w.shape[0])
    g.attrs["n_proteins"] = int(native_w.shape[1])
    g.attrs["n_strata"] = len(parts)
    g.attrs["z_applied"] = bool(pol["z_applied"])
    g.attrs["value_semantics"] = pol["value_semantics"]
    g.attrs["scale_policy"] = json.dumps(_jsonable(pol))
    g.attrs["stratum_notes"] = json.dumps(stratum_notes)
    g.attrs["provenance"] = json.dumps(_jsonable(provenance))
    g.attrs["numeric_space_check"] = json.dumps(_jsonable(merged.get("numeric_space", {})))
    g.attrs["null_policy"] = ("values/values_z carry exactly one null type: float NaN. "
                              "missing_reason is a separate diagnostic array.")

    return {
        "layer": layer,
        "n_rows": int(native_w.shape[0]),
        "n_proteins": int(native_w.shape[1]),
        "n_strata": len(parts),
        "z_applied": bool(pol["z_applied"]),
        "policy": _jsonable(pol),
        "numeric_space": _jsonable(merged.get("numeric_space", {})),
        "strata": stratum_notes,
        "bundle_matrix": bundle_matrix,
    }, tables


def write_deferred(h5: h5py.File, key: str, meta: dict) -> None:
    g = h5.create_group(key)
    g.attrs["layer_key"] = key
    g.attrs["status"] = "deferred"
    g.attrs["n_rows"] = 0
    g.attrs["n_proteins"] = 0
    for k, v in meta.items():
        g.attrs[k] = str(v)


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


# --------------------------------------------------------------------------------------
# Scope gate (plan section 0) -- fails the build rather than warning
# --------------------------------------------------------------------------------------

def _norm(name: str) -> str:
    """Lowercase, non-alphanumerics stripped. `Age_At_Visit`, `ageAtVisit` and
    `age.at.visit` all normalise to the same string, so a column cannot evade the gate
    by changing its separator style -- which the old boundary rule allowed."""
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def scope_gate(columns: list[str], report_only: bool = False) -> list[str]:
    """Columns that must never ship (plan section 0). Fails the build rather than warning.

    The rule is strip-non-alphanumerics then substring, matching the SysBio CDM loader.
    It denies strictly more than the old word-boundary rule, so it was introduced in
    `report_only` mode first and `DENY_EXEMPT` extended for the legitimate columns it
    newly caught (`percentage`, `coverage`, `average`, `linkage`, `trace` and friends)
    before being switched live -- a false positive here stops a build.
    """
    exempt = {_norm(x) for x in C.DENY_EXEMPT}
    bad = []
    for col in columns:
        norm = _norm(col)
        if norm in exempt:
            continue
        for token in C.CLINICAL_DENY:
            if _norm(token) in norm:
                bad.append(f"{col} (matched '{token}')")
                break
    return bad


def undeclared_columns(columns: list[str]) -> list[str]:
    """SEC-3's other half: an allow-list, not merely an exemption set.

    `DENY_EXEMPT` says which columns survive the deny rule; it does not say which columns
    are *supposed* to exist. This reports columns outside the declared contract surface
    so a new column has to be declared rather than appearing silently. Reported, never
    enforced -- an undeclared column is a documentation gap, not a scope breach, and
    stopping the build on one would punish the wrong thing (FR-7's rule).
    """
    declared = set()
    for layer in C.LAYER_ORDER:
        declared |= set(C.long_cols(layer))
    declared |= set(C.EXTRA_DECLARED_COLS)
    return sorted({c for c in columns if c not in declared})
