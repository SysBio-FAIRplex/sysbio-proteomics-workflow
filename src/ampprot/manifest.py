"""S1 - manifest, FILES CDE rows, and the data-derived inventory.

A metadata pass, not a data pass. Cheap enough to re-run after every download, so
"what is still blocked" stays a current fact rather than a stale note.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import config as C

# Synapse IDs the tracking file references in its metadata columns.
REFERENCED_OBJECTS = {
    "syn51757645": ("biospecimen", "DiverseCohorts biospecimen metadata",
                    "anatomic site assignment for diversecohorts_* strata"),
    "syn20827192": ("biospecimen", "Mayo biospecimen metadata",
                    "person_id for mayo_lfq_tcx"),
    "syn21893059": ("biospecimen", "MSBB biospecimen metadata",
                    "person_id for msbb_tmt_phg and msbb_lfq_pfc"),
    "syn53185805": ("assay", "DiverseCohorts assay metadata", "ASSAY.platform"),
    "syn23569441": ("assay", "ROSMAP SRM assay metadata", "ASSAY.platform"),
    "syn21323404": ("assay", "ROSMAP TMT assay metadata", "ASSAY.platform"),
    "syn65414914": ("assay", "ROSMAP SomaScan assay metadata", "ASSAY.platform"),
    "syn65471938": ("assay", "ROSMAP SomaScan assay metadata", "ASSAY.platform"),
    "syn65471939": ("assay", "ROSMAP SomaScan assay metadata", "ASSAY.platform"),
    "syn65473039": ("assay", "ROSMAP SomaScan assay metadata", "ASSAY.platform"),
    "syn23474101": ("assay", "Mayo assay metadata", "ASSAY.platform"),
    "syn22344998": ("assay", "MSBB LFQ assay metadata", "ASSAY.platform"),
    "syn21893060": ("assay", "MSBB TMT assay metadata", "ASSAY.platform"),
    "syn73713766": ("clinical", "Mayo individual metadata",
                    "pmi + pmiUnits for mayo_lfq_tcx (D-J); nothing else read"),
    "syn73713767": ("clinical", "MSBB individual metadata",
                    "pmi + pmiUnits for msbb_tmt_phg and msbb_lfq_pfc (D-J); nothing else read"),
}

# md5 of each referenced object's canonical (DECOMPRESSED) content, with that content's
# byte length. Transcribed from MANIFEST.tsv in the 13-object delivery.
#
# Presence is detected by content hash, never by path. Six of these objects were
# reported absent by the old path-token check while sitting in the tree the whole time,
# because a deposit stores them under its own Synapse-ID folder rather than under theirs
# -- the four SomaScan objects live under syn64957327/, and the two ADKP biospecimen
# tables at the repository root. Path detection would fail the same way on the next
# transfer; a content hash cannot.
OBJECT_CONTENT_MD5: dict[str, tuple[str, int]] = {
    "syn21323404": ("9ae260468fbc0ab6390cb64ffceb5a6f", 72105),
    "syn21893060": ("6b566d1789c428d2d4b053d816356bec", 14300),
    "syn22344998": ("521f3dfdaa961dfd0f776a9122d4fbd4", 16959),
    "syn23474101": ("01ae27cdb09ff5f8d2a70aa1f26f49ea", 14768),
    "syn23569441": ("664b1d28505a850b0273cfe6b0a1bd1b", 35199),
    "syn53185805": ("02c68cc4212a13594515649aa3df0516", 140859),
    "syn51757645": ("8d5c3cf33864c24bdc38586e63d67a08", 651374),
    "syn65414914": ("e40e9fbbe02147c3aa02eb3566b179c1", 2512765),
    "syn65471938": ("f8a3d03ecd1b57cbb409bc6b054a7890", 127604657),
    "syn65471939": ("f8a3d03ecd1b57cbb409bc6b054a7890", 127604657),
    "syn65473039": ("9d77130afa94275672767ace72362edb", 91548),
    "syn20827192": ("676054d647890f68a1cbd34cf5b7b901", 171659),
    "syn21893059": ("4d9e0fcb284d6f3614b2fdc601325b78", 444996),
    # syn73713766 / syn73713767 were not in the delivery -- no hash to check against.
}

# What the hand-maintained tracking file asserts, one entry per row of it, plus the
# stratum-key prefix that selects the strata each row is about. The prefix is the only
# thing added to the transcription -- it turns the assertions into something the build
# can CHECK rather than something a human restates beside the data. Every derived value
# in the diff is computed from the inventory; nothing here is a verdict (VER-1).
TRACKING_ASSERTIONS = [
    dict(grant="AMP AD", assay="TMT LC/MS", study="AMP-AD_DiverseCohorts",
         n_individuals="850?",
         tissue="Brain - dorsolateral prefrontal cortex, superior temporal gyrus",
         prefix="diversecohorts"),
    dict(grant="AMP AD", assay="LC-SRM", study="ROSMAP", n_individuals="?",
         tissue="Brain - dorsolateral prefrontal cortex", prefix="rosmap_srm"),
    dict(grant="AMP AD", assay="TMT LC/MS", study="ROSMAP", n_individuals="610?",
         tissue="Brain - dorsolateral prefrontal cortex", prefix="rosmap_tmt"),
    dict(grant="AMP AD", assay="SomaScan", study="ROSMAP", n_individuals="610?",
         tissue="Brain - dorsolateral prefrontal cortex", prefix="rosmap_soma"),
    dict(grant="AMP AD", assay="LC-MS", study="MayoRNAseq", n_individuals="230",
         tissue="Brain - temporal cortex", prefix="mayo_lfq"),
    dict(grant="AMP AD", assay="LC-MS", study="MSBB", n_individuals="308",
         tissue="Brain - prefrontal cortex", prefix="msbb_lfq"),
    dict(grant="AMP AD", assay="TMT LC/MS", study="MSBB", n_individuals="800",
         tissue="Brain - parahippocampal gyrus", prefix="msbb_tmt"),
    dict(grant="AMP PDRD", assay="Olink", study="Unified Cohorts", n_individuals="375?",
         tissue="Plasma, CSF", prefix="pdrd_olink"),
    dict(grant="AMP PDRD", assay="LC-MS", study="Unified Cohorts", n_individuals="800?",
         tissue="Plasma, CSF", prefix="pdrd_dia"),
]


# Hashing and row-counting are I/O bound on a Windows mount, so both are capped and
# run on a thread pool. Anything above these thresholds is sampled, and the manifest
# records that it was sampled rather than implying a full read.
HASH_PREFIX_BYTES = 64 << 20          # 64 MB
ROWCOUNT_MAX_COMPRESSED = 32 << 20    # only count rows for files below this
MANIFEST_WORKERS = 12


def sha256(path: Path, limit_bytes: int | None = HASH_PREFIX_BYTES) -> str:
    """Full-file sha256, or a prefix hash for large files (labelled as a prefix)."""
    h = hashlib.sha256()
    read = 0
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
            read += len(chunk)
            if limit_bytes is not None and read >= limit_bytes:
                return f"sha256:prefix{limit_bytes}:{h.hexdigest()}"
    return f"sha256:{h.hexdigest()}"


def _peek_shape(path: Path, size: int) -> tuple[int | None, int | None, str | None]:
    """(n_rows, n_cols, delimiter). Row count only for files small enough to scan."""
    inner = path.name[:-3] if path.name.endswith(".gz") else path.name
    if inner.endswith((".xlsx", ".pdf", ".R")):
        return None, None, None
    opener = gzip.open if path.name.endswith(".gz") else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            first = fh.readline()
            if not first:
                return 0, 0, None
            delim = "\t" if first.count("\t") > first.count(",") else ","
            ncol = len(first.split(delim))
            if size > ROWCOUNT_MAX_COMPRESSED:
                return None, ncol, delim       # too big to count honestly, so we do not
            nrow = 1 + sum(1 for _ in fh)
        return nrow, ncol, delim
    except Exception:
        return None, None, None


def _manifest_row(args):
    p, root = args
    rel = p.relative_to(root).as_posix()
    st = p.stat()
    nrow, ncol, delim = _peek_shape(p, st.st_size)
    syn = next((part for part in p.parts if part.startswith("syn")
                and part[3:].isdigit()), None)
    return {
        "file_name": p.name,
        "file_format": _fmt(p.name),
        "file_size_bytes": st.st_size,
        "species": C.HUMAN,
        "cell_type": C.NA,
        "relative_path": rel,
        "synapse_id": syn,
        "sha256": sha256(p),
        "sha256_is_prefix": st.st_size > HASH_PREFIX_BYTES,
        "n_rows": nrow,
        "n_rows_counted": nrow is not None,
        "n_cols": ncol,
        "delimiter": delim,
        "current_version": None,
        "created_on": None,   # mtime is DOWNLOAD time, not creation - never use it here
        "modified_on": None,
        "drs_id": None,
        "fs_mtime_download": datetime.fromtimestamp(
            st.st_mtime, tz=timezone.utc).isoformat(),
    }


def build_file_manifest(root: Path = C.ROOT, workers: int = MANIFEST_WORKERS
                        ) -> pd.DataFrame:
    """Every source file on disk, with the FILES CDE fields we can honestly populate."""
    from concurrent.futures import ThreadPoolExecutor

    skip_dirs = {".venv", "build", "src", ".git", "__pycache__"}
    paths = [p for p in root.rglob("*")
             if p.is_file()
             and not any(part in skip_dirs for part in p.relative_to(root).parts)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(_manifest_row, ((p, root) for p in paths)))
    return pd.DataFrame(rows).sort_values("relative_path").reset_index(drop=True)


def _fmt(name: str) -> str:
    n = name[:-3] if name.endswith(".gz") else name
    ext = n.rsplit(".", 1)[-1].lower() if "." in n else "unknown"
    return {"txt": "tsv", "csv": "csv", "xlsx": "xlsx", "tsv": "tsv",
            "pdf": "pdf", "r": "R script"}.get(ext, ext)


def _uncompressed_size(path: Path) -> int:
    """Byte length of a file's content, reading gzip's ISIZE trailer rather than
    inflating. ISIZE is modulo 2**32, so it is a pre-filter for the hash, not proof."""
    if path.suffix != ".gz":
        return path.stat().st_size
    try:
        with path.open("rb") as fh:
            fh.seek(-4, 2)
            return int.from_bytes(fh.read(4), "little")
    except OSError:
        return -1


def content_md5(path: Path) -> str:
    """md5 of a file's DECOMPRESSED content, so a gzipped copy of an object in the tree
    hashes identically to the plain file the depositor published."""
    opener = gzip.open if path.suffix == ".gz" else open
    h = hashlib.md5()
    with opener(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _content_index(root: Path) -> dict[str, list[str]]:
    """md5 -> relative paths, over files whose decompressed size matches an object we
    are looking for. The size pre-filter keeps this off the ~64 GB spectrum-level files:
    without it, indexing the tree by content would cost more than the build."""
    wanted = {size for _, size in OBJECT_CONTENT_MD5.values()}
    index: dict[str, list[str]] = {}
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        parts = p.relative_to(root).parts
        if any(part in {".venv", "build", "src", ".git", "__pycache__"} for part in parts):
            continue
        if _uncompressed_size(p) not in wanted:
            continue
        try:
            index.setdefault(content_md5(p), []).append(p.relative_to(root).as_posix())
        except OSError:
            continue
    return index


def build_missing_objects(root: Path = C.ROOT) -> pd.DataFrame:
    """The section 1.2 audit, as data. Presence is a content fact (VER-1)."""
    index = _content_index(root)
    rows = []
    for syn, (kind, desc, blocks) in REFERENCED_OBJECTS.items():
        want = OBJECT_CONTENT_MD5.get(syn)
        found = index.get(want[0], []) if want else []
        rows.append({
            "synapse_id": syn,
            "metadata_kind": kind,
            "description": desc,
            "present_on_disk": bool(found),
            "found_at": "; ".join(found) or None,
            "content_md5": want[0] if want else None,
            "detection": "content_hash" if want else "no_reference_hash",
            "blocks": blocks,
            "on_critical_path": kind != "clinical",
        })
    return pd.DataFrame(rows).sort_values(
        ["on_critical_path", "metadata_kind", "synapse_id"], ascending=[False, True, True])


def build_derived_inventory(probes: dict[str, dict]) -> pd.DataFrame:
    """One row per stratum, every field computed from the data by the readers."""
    rows = []
    for s in C.STRATA:
        pr = probes.get(s.key, {})
        rows.append({
            "layer_key": s.layer,
            "stratum_key": s.key,
            "grant": s.grant,
            "study": s.study,
            "cohort": s.cohort,
            "tissue": s.tissue,
            "anatomic_site": s.anatomic_site,
            "assay_type": s.assay_type,
            "platform": s.platform,
            "platform_evidence": s.platform_evidence or "C",
            "analysis_pipeline": s.analysis_pipeline,
            "pipeline_evidence": s.pipeline_evidence or "C",
            "primary_file": s.matrix_path,
            "n_features": pr.get("n_features"),
            "n_samples": pr.get("n_samples"),
            "n_participants": pr.get("n_participants"),
            "n_pool_samples": pr.get("n_pool"),
            "n_timepoints": pr.get("n_timepoints"),
            "timepoint_values": pr.get("timepoint_values"),
            "is_longitudinal": pr.get("is_longitudinal"),
            "native_scale": s.value_scale,
            "transform_chain": s.transform_chain,
            "value_p1": pr.get("value_p1"),
            "value_p99": pr.get("value_p99"),
            "pct_missing": pr.get("pct_missing"),
            "sample_key_pattern": pr.get("sample_key_pattern"),
            "person_key_source": s.annot_reader,
            "person_id_resolvable": s.person_resolvable,
            "biospecimen_metadata_present": s.annot_path is not None,
            "clinical_metadata_present": False,   # by design, plan section 0
            "build_status": pr.get("build_status", "not_built"),
            "blocked_reason": s.blocked_reason,
        })
    for key, meta in C.DEFERRED_STRATA.items():
        rows.append({
            "layer_key": meta["layer"], "stratum_key": key, "grant": meta["grant"],
            "study": meta["study"], "platform": meta.get("platform"),
            "platform_evidence": "A",
            "analysis_pipeline": meta.get("analysis_pipeline"), "pipeline_evidence": "A",
            "build_status": "deferred", "blocked_reason": meta["reason"],
            "person_id_resolvable": True,
        })
    for key, meta in C.DEFERRED_LAYERS.items():
        rows.append({
            "layer_key": key, "stratum_key": None, "grant": meta["grant"],
            "tissue": meta["tissue"], "assay_type": meta["assay_type"],
            "build_status": "deferred", "blocked_reason": meta["reason"],
        })
    return pd.DataFrame(rows)


def _stated_int(v: str) -> int | None:
    """The tracking file writes approximations as '610?' and unknowns as '?'."""
    digits = "".join(ch for ch in str(v) if ch.isdigit())
    return int(digits) if digits else None


def build_tracking_diff(inv: pd.DataFrame, missing: pd.DataFrame | None = None
                        ) -> pd.DataFrame:
    """Where the derived inventory disagrees with the hand-maintained tracking file.

    Every derived value is computed from `inv`. Nothing is transcribed: a hardcoded row
    here published a wrong verdict about the DiverseCohorts temporal site for three weeks
    (N2), which is the precise failure VER-1 predicts.
    """
    built = inv[inv.stratum_key.notna() & (inv.build_status == "built")]
    rows = []

    for a in TRACKING_ASSERTIONS:
        sel = built[built.stratum_key.str.startswith(a["prefix"], na=False)]
        label = f"{a['grant']} / {a['assay']} / {a['study']}"
        stated = _stated_int(a["n_individuals"])

        if sel.empty:
            rows.append(dict(
                tracking_row=label, field="Number Individuals",
                tracking_value=a["n_individuals"], derived_value=None, agrees=None,
                note=f"no built stratum matches prefix '{a['prefix']}'"))
            continue

        derived = sel["n_participants"].dropna()
        derived = int(derived.sum()) if len(derived) else None
        rows.append(dict(
            tracking_row=label, field="Number Individuals",
            tracking_value=a["n_individuals"], derived_value=derived,
            agrees=None if (stated is None or derived is None) else stated == derived,
            note="; ".join(f"{r.stratum_key}={r.n_participants:g}"
                           for r in sel.itertuples() if pd.notna(r.n_participants))))

        sites = sorted({str(x) for x in sel["anatomic_site"].dropna().unique()})
        stated_site = a["tissue"].lower()
        agrees_site = all(any(tok in stated_site for tok in s.lower().split())
                          for s in sites) if sites else None
        rows.append(dict(
            tracking_row=label, field="Biospecimen source",
            tracking_value=a["tissue"], derived_value=", ".join(sites) or None,
            agrees=agrees_site,
            note="derived from the biospecimen metadata, not from the matrix filename"))

    if missing is not None and len(missing):
        n_present = int(missing.present_on_disk.sum())
        rows.append(dict(
            tracking_row="all rows", field="metadata syn IDs",
            tracking_value=f"{len(missing)} objects referenced",
            derived_value=f"{n_present} present",
            agrees=n_present == len(missing),
            note="presence detected by content hash, not by path"))
    return pd.DataFrame(rows)


def write_inventory(inv, diff, missing, files, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for name, df in [("derived_inventory", inv), ("tracking_file_diff", diff),
                     ("missing_objects", missing), ("file_manifest", files)]:
        df.to_csv(outdir / f"{name}.tsv", sep="\t", index=False)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_source_files": int(len(files)),
        "total_bytes": int(files.file_size_bytes.sum()),
        "n_strata_built": int((inv.build_status == "built").sum()),
        "n_strata_deferred": int((inv.build_status == "deferred").sum()),
        "n_referenced_objects": int(len(missing)),
        "n_objects_present": int(missing.present_on_disk.sum()),
        # `agrees` is tri-state: True, False, or None where the build cannot compute the
        # comparison. Only an explicit False is a disagreement -- None is an
        # UNVERIFIABLE, and counting it as a disagreement would overstate the diff.
        "tracking_file_disagreements": int((diff.agrees == False).sum()),  # noqa: E712
        "tracking_file_unverifiable": int(diff.agrees.isna().sum()),
    }
    (outdir / "inventory_summary.json").write_text(json.dumps(summary, indent=2))
