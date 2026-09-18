"""Orchestrator: read -> harmonise -> decide policy -> verify -> write one HDF5."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd
import gc
import os
import pathlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

# Memory is the real constraint, not core count: a stratum peaks around 1-2 GB against
# ~29 GB available, so let the budget bind and keep the worker ceiling generous.
WORKERS = min(16, os.cpu_count() or 8)

# Rough peak resident cost of harmonising one stratum, as a multiple of its compressed
# source size. Decompression plus pandas object overhead dominates; xlsx via openpyxl and
# the Olink long-form pivots are the expensive ones.
EXPANSION = {
    "csv_matrix": 12, "msbb_tmt_xlsx": 90, "maxquant": 30,
    "srm": 15, "somascan": 12, "olink": 45,
    # DIA reads an uncompressed long-form CSV and pivots it; the pivot is small (a few
    # hundred features) but the parsed long frame is several times the file on disk.
    "dia": 8,
}
MIN_FOOTPRINT = 128 << 20     # never assume a stratum is free


def _available_bytes() -> int:
    """Available RAM, from the kernel rather than assumed."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return 4 << 30


def _footprint(st: C.Stratum) -> int:
    p = C.ROOT / st.matrix_path
    if p.is_dir():                      # Olink: four panel files
        size = sum(f.stat().st_size for f in p.glob("*.csv.gz"))
    elif p.exists():
        size = p.stat().st_size
    else:
        size = MIN_FOOTPRINT
    return max(int(size * EXPANSION.get(st.reader, 20)), MIN_FOOTPRINT)


def _schedule(strata: list[C.Stratum], budget: int) -> list[list[C.Stratum]]:
    """Greedy bin-packing so concurrent strata stay inside the memory budget."""
    waves, wave, used = [], [], 0
    for st in sorted(strata, key=_footprint, reverse=True):
        need = _footprint(st)
        if wave and (used + need > budget or len(wave) >= WORKERS):
            waves.append(wave)
            wave, used = [], 0
        wave.append(st)
        used += need
    if wave:
        waves.append(wave)
    return waves

from . import config as C, harmonize as H, manifest as M, readers as R, verify as V, writer as W


# Structural shape of each source's sample identifier. Published instead of a real key
# so the inventory documents the format without carrying a participant or specimen ID.
SAMPLE_KEY_PATTERN = {
    "csv_matrix": "[<dataset>_]b<batch>.<tmt_channel>",
    "pd_ratio": "F<plex>.<tmt_channel>",
    "msbb_tmt_xlsx": "TMT<batch>_<group>_<specimen>",
    "maxquant": "[<study>_]b<batch>_<position>_<specimen>",
    "srm": "P<plate>_<well>",
    "somascan": "<individual_id>_<visit>",
    "olink": "<cohort>-<participant>-<visit><timepoint>-<tissue>-<assay>",
}


def _harmonise_one(stratum_key: str):
    """Module-level so a process pool can pickle it. Returns a Harmonised object."""
    st = next(s for s in C.STRATA if s.key == stratum_key)
    h = H.harmonise_stratum(st, _SYM2UNI)
    h.samples = H.attach_person_identity(h.samples, st.key)
    return h


_SYM2UNI: dict = {}


def _init_worker(sym2uni):
    global _SYM2UNI
    _SYM2UNI = sym2uni


def _assay_table(parts) -> pd.DataFrame:
    rows = []
    for i, p in enumerate(parts, start=1):
        s = p.stratum
        plat, plat_ev = H.observed_platform(p)
        rows.append({
            "assay_id": f"{s.layer}:{s.key}",
            "assay_source_value": f"{s.study} {s.assay_type} {s.anatomic_site}".strip(),
            "assay_type": s.assay_type,
            "platform": plat,
            "platform_evidence": plat_ev,
            "suspension_type": s.suspension_type,
            "analyte_type": s.analyte_type,
            "analysis_pipeline": s.analysis_pipeline,
            "pipeline_evidence": s.pipeline_evidence or "C",
            "stratum_key": s.key,
            "study": s.study,
            "grant": s.grant,
            "cohort": s.cohort,
            "transform_applied": p.notes.get("transform_applied"),
            "zscored": p.notes.get("zscored"),
            "value_scale": s.value_scale,
            "value_unit": s.value_unit,
            "transform_chain": s.transform_chain,
            "mass_analyzer": s.extra.get("mass_analyzer"),
            "activation_type": s.extra.get("activation_type"),
            "ms_order": s.extra.get("ms_order"),
            "scan_filter": s.extra.get("scan_filter"),
            "array_type": s.array_type,
        })
    return pd.DataFrame(rows)


def _files_table(parts, file_manifest: pd.DataFrame) -> pd.DataFrame:
    idx = file_manifest.set_index("relative_path")
    rows = []
    for p in parts:
        s = p.stratum
        # A stratum's "matrix" may be several files (Olink ships one per panel), so emit
        # one FILES row per actual file rather than one per stratum.
        entries = []
        mp = C.ROOT / s.matrix_path
        if mp.is_dir():
            # Only the files actually read. Globbing every *.csv.gz in the directory made
            # Olink plasma emit 8 FILES rows for the 4 files the reader opens, listing
            # the `_retracted` siblings as inputs to a build that does not read them
            # (WP-11 step 6). They stay on disk untouched (IN-5); they are just not
            # claimed as provenance.
            code, variant = s.extra.get("matrix_code"), s.extra.get("variant", "")
            names = [f"{code}_olink_explore_format_{panel}{variant}.csv.gz"
                     for panel in C.OLINK_PANELS] if code else None
            files_used = ([mp / n for n in names if (mp / n).exists()] if names
                          else sorted(mp.glob("*.csv.gz")))
            for fp in files_used:
                entries.append(("protein abundance matrix",
                                fp.relative_to(C.ROOT).as_posix()))
        else:
            entries.append(("protein abundance matrix", s.matrix_path))
        if s.annot_path:
            entries.append(("sample annotation", s.annot_path))

        for role, path in entries:
            meta = idx.loc[path].to_dict() if path in idx.index else {}
            rows.append({
                "file_id": f"{s.key}:{path.rsplit('/', 1)[-1]}",
                "file_name": meta.get("file_name", path.rsplit("/", 1)[-1]),
                "current_version": None,          # Synapse version not downloaded
                "assay_id": f"{s.layer}:{s.key}",
                "file_role": role,
                "study": s.study,
                "grant": s.grant,
                "array_type": s.array_type,
                "analysis_type": "protein quantification",
                "biosample_type": s.biosample_type,
                "tissue": s.tissue,
                "cell_type": C.NA,
                "species": C.HUMAN,
                "processing_status": "harmonised",
                "file_format": meta.get("file_format"),
                "file_size_bytes": meta.get("file_size_bytes"),
                "created_on": None,               # mtime is download time, not creation
                "modified_on": None,
                "drs_id": None,
                "relative_path": path,
                "synapse_id": meta.get("synapse_id"),
                "sha256": meta.get("sha256"),
            })
    return pd.DataFrame(rows)


def _layer_samples(merged, layer, files_df: pd.DataFrame) -> pd.DataFrame:
    """Assemble the sample table to the uniform contract."""
    s = merged["samples"].copy()
    strata = {x.key: x for x in C.strata_for(layer)}

    # FR-11: resolve the real file identifier instead of composing a role placeholder.
    # `_files_table` already knows the ids, so it is built first and passed in. The
    # relation is one-to-many for Olink (4 CSF panel files, 8 plasma), and a list column
    # would break "one row per sample" -- so where a stratum has several matrix files
    # this column stays null and the `sample_files` bridge (D-F) is authoritative. A
    # populated `file_id` therefore always means a genuine 1:1, never a guess.
    matrix_files: dict[str, list[dict]] = {}
    for _, r in files_df.iterrows():
        if r.get("file_role") != "protein abundance matrix":
            continue
        matrix_files.setdefault(str(r["assay_id"]).rsplit(":", 1)[-1], []).append(r.to_dict())

    def _sole(key, field):
        hits = matrix_files.get(key, [])
        return hits[0].get(field) if len(hits) == 1 else None

    # `platform` is the one CDE a source can state PER SAMPLE, and six strata now do
    # (WP-8). The config declaration is a fallback for the rows the metadata does not
    # reach, never an overwrite of a value read off the deposit -- filling first and
    # declaring second is what makes `platform_evidence=A` mean something.
    declared_platform = pd.Series([strata[k].platform for k in s.stratum_key],
                                  index=s.index, dtype=object)
    if "platform" in s.columns:
        s["platform"] = s["platform"].where(s["platform"].notna(), declared_platform)
    else:
        s["platform"] = declared_platform

    # stratum-derived CDE columns
    for col, get in [
        ("assay_id", lambda st: f"{layer}:{st.key}"),
        ("assay_source_value", lambda st: f"{st.study} {st.assay_type} {st.anatomic_site}"),
        ("assay_type", lambda st: st.assay_type),
        ("suspension_type", lambda st: st.suspension_type),
        ("analyte_type", lambda st: st.analyte_type),
        ("analysis_pipeline", lambda st: st.analysis_pipeline),
        ("study", lambda st: st.study),
        ("grant", lambda st: st.grant),
        ("cohort", lambda st: st.cohort),
        ("tissue", lambda st: st.tissue),
        ("biosample_type", lambda st: st.biosample_type),
        ("anatomic_site", lambda st: st.anatomic_site),
        ("visit_name", lambda st: st.visit_name),
        ("value_scale", lambda st: st.value_scale),
        ("value_unit", lambda st: st.value_unit),
        ("array_type", lambda st: st.array_type),
    ]:
        s[col] = [get(strata[k]) for k in s.stratum_key]

    for col, field in [("file_id", "file_id"), ("file_name", "file_name"),
                       ("file_format", "file_format"),
                       ("file_size_bytes", "file_size_bytes")]:
        s[col] = [_sole(k, field) for k in s.stratum_key]

    s["cell_type"] = C.NA
    s["species"] = C.HUMAN
    s["analysis_type"] = "protein quantification"
    s["processing_status"] = "harmonised"
    s["current_version"] = None
    s["created_on"] = None
    s["modified_on"] = None
    s["drs_id"] = None
    s["file_role"] = "protein abundance matrix"

    for col in ("visit_month", "visit_index", "batch", "channel", "plate_id",
                "qc_status", "is_pool"):
        if col not in s.columns:
            s[col] = None
    s["is_pool"] = s["is_pool"].fillna(False).astype(bool)
    s["qc_status"] = s["qc_status"].fillna("not_assessed")
    s["specimen_id"] = (s.stratum_key.astype(str) + ":"
                        + s.specimen_source_value.astype(str))
    s["visit_occurrence_id"] = np.where(
        s.visit_month.notna(),
        s.specimen_id.astype(str) + "|M" + s.visit_month.astype(str),
        np.where(s.visit_index.notna(),
                 s.specimen_id.astype(str) + "|V" + s.visit_index.astype(str), None))
    s["replicate_index"] = W.replicate_index(s)
    return s


def run(verbose=True) -> dict:
    t0 = datetime.now(timezone.utc)
    C.BUILD.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("S1  manifest ...", flush=True)
    files = M.build_file_manifest()
    missing = M.build_missing_objects()

    sym2uni = H._load_uniprot2symbol()
    audit, layer_reports, probes = [], [], {}
    _layer_cols, _feat_cols = {}, {}
    # Small per-layer tables, held back so the column union can be computed across every
    # layer before any of them is written (WP-3). The value matrices are not held.
    pending_tables: dict[str, dict] = {}
    # WP-13: bundle descriptors and the union-reindexed sample tables, held so the
    # per-grant return path can be built and its boundary demonstrated once every layer
    # is final. Small -- descriptors and sample rows, never a value matrix.
    bundles: dict[str, dict] = {}
    samples_by_layer: dict[str, pd.DataFrame] = {}
    # The tracking-file and README audits are NOT emitted here any more. They are
    # computed from the built objects after the build (WP-6), because a verdict written
    # before the data is read is a transcription, not a check.

    # D-I. The unconditional `H5_PATH.unlink()` that used to stand here existed only
    # because the filename was a fixed literal, so every build had to destroy the last
    # one before it could start. With the name datestamped that is gone -- which retires
    # OPS-1 rather than patching it. The temp-then-rename remains, because a SAME-DAY
    # rerun still resolves to the same name and a crash mid-write would otherwise leave a
    # truncated file where a good artifact was.
    previous = C.previous_build()
    tmp_path = C.H5_PATH.with_suffix(".h5.partial")
    if tmp_path.exists():
        tmp_path.unlink()

    with h5py.File(tmp_path, "w") as h5:
        for layer in C.LAYER_ORDER:
            strata = C.strata_for(layer)
            if not strata:
                W.write_deferred(h5, layer, C.DEFERRED_LAYERS.get(layer, {}))
                if verbose:
                    print(f"  {layer:26s} deferred", flush=True)
                continue

            # Strata are independent and the work is dominated by decompression and
            # parsing, so a thread pool overlaps the I/O. Concurrency is bounded by a
            # measured memory budget rather than a fixed worker count, so a layer with
            # several expensive strata degrades to fewer parallel readers instead of
            # exhausting RAM.
            # Processes, not threads: openpyxl is pure Python and much of the pandas
            # work holds the GIL, so threads serialise. Separate interpreters give real
            # parallelism; results are 10-40 MB per stratum, cheap to pickle back.
            budget = int(_available_bytes() * 0.45)
            waves = _schedule(strata, budget)
            by_key = {}
            for wave in waves:
                if len(wave) > 1:
                    with ProcessPoolExecutor(max_workers=len(wave),
                                             initializer=_init_worker,
                                             initargs=(sym2uni,)) as pool:
                        for h in pool.map(_harmonise_one, [s.key for s in wave]):
                            by_key[h.stratum.key] = h
                else:
                    _init_worker(sym2uni)
                    h = _harmonise_one(wave[0].key)
                    by_key[h.stratum.key] = h
                gc.collect()
            parts = [by_key[st.key] for st in strata]     # restore declared order
            if verbose and len(waves) > 1:
                print(f"      ({len(strata)} strata in {len(waves)} memory-bounded waves, "
                      f"budget {budget/1e9:.1f} GB)", flush=True)
            for h in parts:
                audit += (V.verify_stratum(h) + V.verify_sparsity(h)
                          + V.verify_exclusions(h) + V.verify_source_qc(h)
                          + V.verify_qc_metadata(h))

            policy = H.decide_layer_policy(layer, parts)
            H.apply_layer_policy(parts, policy)
            merged = H.merge_layer(parts, policy)
            # Files first: `samples.file_id` resolves against it (FR-11, WP-5).
            files_df = _files_table(parts, files)
            merged["samples"] = _layer_samples(merged, layer, files_df)
            audit += V.verify_layer(layer, merged, parts)
            audit += V.verify_file_ids(layer, merged["samples"], files_df)

            bad = W.scope_gate(list(merged["samples"].columns)
                               + list(merged["features"].columns))
            if bad:
                raise SystemExit(f"SCOPE GATE FAILED for {layer}: {bad}")

            rep, tables = W.write_layer(
                h5, layer, merged, parts, _assay_table(parts), files_df,
                {"source_files": [p.stratum.matrix_path for p in parts]})
            layer_reports.append(rep)
            pending_tables[layer] = tables
            for p in parts:
                # WP-14 step 4: `build_derived_inventory` reads these three keys from
                # `probes`, but nothing ever wrote them, so all three shipped blank for
                # every stratum while HARMONIZATION_PLAN section 1.5 promised timepoint
                # sets there.
                vm = p.samples.get("visit_month")
                tp = sorted(int(v) for v in vm.dropna().unique()) if vm is not None else []
                probes[p.stratum.key] = {
                    "n_features": int(p.native.shape[0]),
                    "n_samples": int(p.native.shape[1]),
                    "n_participants": int(p.samples.person_id.dropna().nunique()),
                    "n_pool": int(p.samples.get("is_pool",
                                                pd.Series(dtype=bool)).fillna(False).sum()),
                    "n_timepoints": len(tp),
                    "timepoint_values": ",".join(str(t) for t in tp) or None,
                    "is_longitudinal": len(tp) > 1,
                    "value_p1": float(np.nanpercentile(p.native, 1)),
                    "value_p99": float(np.nanpercentile(p.native, 99)),
                    "pct_missing": round(float(100 * np.isnan(p.native).mean()), 2),
                    "sample_key_pattern": SAMPLE_KEY_PATTERN.get(
                        p.stratum.reader, "<opaque source identifier>"),
                    "build_status": "built",
                    "transform_applied": p.notes.get("transform_applied"),
                    "zscored": p.notes.get("zscored"),
                }
            if verbose:
                print(f"  {layer:26s} {rep['n_rows']:5d} rows x {rep['n_proteins']:6d} "
                      f"proteins  z={rep['z_applied']}", flush=True)

        # Second write pass (WP-3). Every layer's samples/features/assay/files table is
        # reindexed to the union of columns actually emitted across all layers, so a
        # column that does not apply to a layer ships present-and-null rather than
        # absent. This is what makes the uniform-contract guarantee hold; before it,
        # five of six layers were short between three and six sample columns and all six
        # were short feature columns, while the audit reported the contract satisfied.
        if pending_tables:
            unions = W.column_unions(pending_tables)
            all_cols = unions.get("samples", []) + unions.get("features", [])
            bad = W.scope_gate(all_cols)
            if bad:
                raise SystemExit(f"SCOPE GATE FAILED on the column union: {bad}")
            # SEC-3's allow-list half (WP-7 item 4). The function existed and nothing
            # called it, so no allow-list report was ever produced -- the same
            # defined-but-never-run failure as the metadata objects WP-8 landed.
            audit += V.verify_undeclared_columns(W.undeclared_columns(all_cols),
                                                 len(all_cols))
            by_layer = {r["layer"]: r for r in layer_reports}
            for layer, tables in pending_tables.items():
                filled = W.write_layer_tables(h5, layer, tables, unions)
                # WP-13: the bundle's samples/features, from the same union-reindexed
                # frames the HDF5 just received.
                bundles[layer] = W.write_bundle_tables(
                    layer, tables, unions, by_layer[layer]["bundle_matrix"],
                    {"source_files": by_layer[layer].get("strata", {})})
                samples_by_layer[layer] = tables["samples"]
                _layer_cols[layer] = list(unions["samples"])
                _feat_cols[layer] = list(unions["features"])
                # Record the repair as data, not as a claim: how many columns each layer
                # gained by the union, so the fix is measurable in the audit itself.
                audit.append({
                    "check": "column_union_fill", "subject": layer,
                    "claim": "layer reindexed to the cross-layer column union",
                    "claim_source": "computed",
                    "verdict": V.PASS,
                    "observed": "; ".join(f"{n}: +{c}" for n, c in sorted(filled.items())
                                          if c) or "already complete",
                })
            if verbose:
                print(f"  column union: "
                      + ", ".join(f"{n}={len(c)}" for n, c in sorted(unions.items())),
                      flush=True)

        # ROSMAP projid <-> individualID bridge (syn3191087), identifiers only.
        # Lets ROSMAP plasma (individualID namespace) join to ROSMAP brain (projid).
        try:
            cw = R.load_rosmap_crosswalk()
            W._write_table(h5, "person_crosswalk_rosmap", cw)
            h5["person_crosswalk_rosmap"].attrs["source"] = "syn3191087 (ROSMAP_clinical.csv)"
            h5["person_crosswalk_rosmap"].attrs["note"] = (
                "identifier columns only - the source file is predominantly clinical and "
                "the rest is dropped at read per the scope gate")
            if verbose:
                print(f"  person_crosswalk_rosmap    {len(cw):5d} projid<->individualID pairs",
                      flush=True)
        except Exception as e:
            print(f"  crosswalk unavailable: {e}", flush=True)

        h5.attrs["build_utc"] = t0.isoformat()
        h5.attrs["cdm"] = "SysBio CDM (OMOP v5.4 + ASSAY/FILES extension)"
        h5.attrs["layers"] = json.dumps(C.LAYER_ORDER)
        h5.attrs["id_policy"] = ("all identifiers are strings, exactly as defined by the "
                                 "source; prefixes are never stripped")
        h5.attrs["artifact_name"] = C.H5_NAME
        h5.attrs["previous_build"] = previous or ""
        h5.attrs["versioning"] = ("the filename is the version (D-I) - resolve the newest "
                                  "rather than pinning a name, or pin deliberately for "
                                  "reproducibility")

    # Only now does the artifact take its real name: a reader can never observe a
    # half-written file under it.
    tmp_path.replace(C.H5_PATH)

    if _layer_cols:
        audit += V.verify_column_consistency(_layer_cols, _feat_cols)

    # WP-13 step 4: the return path, demonstrated once per grant on every build.
    returns = W.write_returns(bundles, samples_by_layer) if samples_by_layer else []
    if verbose and returns:
        print("  returns: " + ", ".join(
            f"{r['grant']}={r['n_specimens_returned']}" for r in returns), flush=True)

    inv = M.build_derived_inventory(probes)
    diff = M.build_tracking_diff(inv, missing)
    M.write_inventory(inv, diff, missing, files, C.BUILD / "inventory")

    # Computed, not transcribed (WP-6, VER-1).
    audit += V.verify_tracking_file(diff)
    audit += V.verify_readme_claims(layer_reports, probes)
    audit += V.verify_bundles(bundles, layer_reports)
    audit += V.verify_returns(returns)

    audit_df = pd.DataFrame(audit)
    audit_df.to_csv(C.BUILD / "inventory" / "assertion_audit.tsv", sep="\t", index=False)

    report = {
        "build_utc": t0.isoformat(),
        "h5": str(C.H5_PATH.relative_to(C.ROOT)),
        "h5_bytes": C.H5_PATH.stat().st_size,
        # D-I: the lineage is published rather than inferred from a directory listing.
        "artifact_name": C.H5_NAME,
        "previous_build": previous,
        "versioning": ("the filename is the version - resolve the newest rather than "
                       "pinning a name, or pin deliberately for reproducibility"),
        "layers": layer_reports,
        "tabs": bundles,
        "returns": returns,
        "audit": V.summarise(audit),
    }
    (C.BUILD / "build_report.json").write_text(json.dumps(report, indent=2, default=str))

    from . import readme as RM
    report["readme"] = str(pathlib.Path(
        RM.write(report, audit_df, layer_reports, inv)).relative_to(C.ROOT))
    (C.BUILD / "build_report.json").write_text(json.dumps(report, indent=2, default=str))
    return report, audit_df, layer_reports


if __name__ == "__main__":
    rep, audit, layers = run()
    print("\n" + json.dumps(rep["audit"], indent=2))
    print(f"\nwrote {rep['h5']}  ({rep['h5_bytes']/1e6:.1f} MB)")
    fails = audit[audit.verdict == "FAIL"]
    if len(fails):
        print(f"\n{len(fails)} assertion FAILs:")
        for _, r in fails.iterrows():
            print(f"  [{r['check']}] {r['subject']}: claim={r['claim']!r} "
                  f"({r['claim_source']}) -> {r['observed']}")
    sys.exit(0)
