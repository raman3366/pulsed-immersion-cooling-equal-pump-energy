#!/usr/bin/env python3
"""Readers for the K0 housing case monitors (WP 2.5c).

The monitors are rendered by scripts/model/housing_physics.monitors() into
system/monitors: every function object writes
postProcessing/<region>/<name>/<startTime>/<file>.dat (v2606 region
function objects). A restarted run starts a new <startTime> segment; the
stitched series truncates each segment at the start of the next one (the
NEW-25 rule used by gates.continuity_gate, which judges the splice itself).

- case_meta(case)                      system/caseMeta.yaml
- segments(case, region, name)         [(label, FOData)] sorted by start
- series(case, region, name)           stitched (t, {column: array}); an
                                       INIT_FROM member is prefixed with its
                                       shared state's samples up to t_pre
- column(case, region, name, pattern)  one stitched column by name match
- minmax(case, region, name)            stitched fieldMinMax (value + location)
- solver_logs(case)                    log.chtMultiRegionFoam[.heat][.zst]
                                       in run order
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import v1_io  # noqa: E402  (header-aware FO reader, log reader)


def case_meta(case) -> dict:
    return yaml.safe_load((Path(case) / "system/caseMeta.yaml").read_text())


def claim_status(meta: dict) -> tuple[bool, list]:
    """(claim_eligible, notes) of a case from its caseMeta (written by
    housing_physics.case_meta: diagnostics, sensitivities and mesh options
    pending their D-row, e.g. K0_ZLEVER=module, are not claim-eligible). A
    caseMeta without the record is NOT eligible. scripts/post/claims.py
    refuses a filled claim that cites a non-eligible run."""
    ok = meta.get("claim_eligible") is True
    notes = list(meta.get("claim_notes") or [])
    if "claim_eligible" not in meta:
        notes.append("caseMeta has no claim_eligible record")
    return ok, notes


def _dirs(case, region, name):
    base = Path(case) / "postProcessing"
    for cand in (base / region / name, base / name / region, base / name):
        if cand.is_dir():
            return cand
    return None


def _segment_file(sd, name):
    """The one .dat file of a monitor segment directory. v2606 writeFile
    (writeFile.C newFileAtTime) opens <file>_<startTime>.dat beside an
    existing <file>.dat when a second solver run starts at the same time
    (a resume from the same written time after a crash before the next
    write): two runs in one segment. Stitching those is not implemented,
    and reading either file alone would be wrong, so it is an input error
    (None for a directory without .dat)."""
    files = sorted(sd.glob("*.dat"))
    if len(files) > 1:
        raise ValueError(
            f"{name}: segment {sd} holds {len(files)} files "
            f"({', '.join(f.name for f in files)}): two solver runs started "
            f"at t = {sd.name} (a second resume from the same written time); "
            "same-start segments are not stitched - keep the later run's "
            "file or rerun from an earlier write")
    return files[0] if files else None


def segments(case, region, name):
    d = _dirs(case, region, name)
    if d is None:
        raise FileNotFoundError(f"{case}: no postProcessing for {region}/{name}")
    segs = []
    for sd in d.iterdir():
        try:
            label = float(sd.name)
        except ValueError:
            continue
        f = _segment_file(sd, name)
        if f is None:
            continue
        segs.append((label, v1_io.read_fo(f)))
    segs.sort(key=lambda s: s[0])
    if not segs:
        raise FileNotFoundError(f"{case}: {region}/{name} has no .dat segment")
    return segs


def init_source(case):
    """(source case, t_pre) of an INIT_FROM member (init_from.json), else
    (None, None)."""
    j = Path(case) / "init_from.json"
    if not j.exists():
        return None, None
    import json
    rec = json.loads(j.read_text())
    repo = next((q for q in Path(case).resolve().parents
                 if (q / "data/reference/properties.yaml").exists()),
                Path(os.environ["PULSED_REPO"]) if os.environ.get("PULSED_REPO") else None)
    cands = [Path(rec["source"])]
    if rec.get("source_rel") and repo is not None:
        cands.append(repo / rec["source_rel"])       # pulled to another host
    cands.append(Path(case).resolve().parent / Path(rec["source"]).name)
    src = next((c for c in cands if c.is_dir()), cands[0])
    return src, float(rec["t_pre"])


def series(case, region, name, with_source=True):
    """Stitched (t, {column: values}); segment k is truncated at the first
    time of segment k+1. An INIT_FROM member's monitors start after t_pre
    (function objects do not write at a restart's start time), so its
    series is prefixed with the shared state's own samples up to t_pre
    (the same state: the member copied it; a tilt seed is accounted
    separately, tilt.json)."""
    t, cols = _series(case, region, name)
    src, t_pre = init_source(case) if with_source else (None, None)
    if src is None or not src.is_dir():
        return t, cols
    ts, cs = _series(src, region, name)
    keep = ts <= t_pre + 1e-9 * max(1.0, t_pre)
    mine = t > t_pre + 1e-9 * max(1.0, t_pre)
    return (np.concatenate([ts[keep], t[mine]]),
            {c: np.concatenate([cs[c][keep], v[mine]]) for c, v in cols.items()})


def _series(case, region, name):
    segs = segments(case, region, name)
    cols = segs[0][1].columns
    parts_t, parts = [], {c: [] for c in cols[1:]}
    for k, (_, fo) in enumerate(segs):
        if fo.columns != cols:
            raise ValueError(f"{name}: columns changed between segments "
                             f"({cols} vs {fo.columns})")
        t = fo.t
        if t.size == 0:
            continue
        keep = np.ones(t.size, bool)
        if k + 1 < len(segs) and segs[k + 1][1].t.size:
            keep = t < segs[k + 1][1].t[0] - 1e-12
        parts_t.append(t[keep])
        for c in cols[1:]:
            parts[c].append(np.asarray(fo.data[c])[keep])
    t = np.concatenate(parts_t) if parts_t else np.array([])
    return t, {c: np.concatenate(v) for c, v in parts.items()}


def column(case, region, name, pattern=None):
    """(t, values) of the single column whose name matches `pattern`
    (regex search; default: the only data column)."""
    t, cols = series(case, region, name)
    keys = [c for c in cols if pattern is None or re.search(pattern, c)]
    if len(keys) != 1:
        raise KeyError(f"{name}: column {pattern!r} matches {keys} of "
                       f"{list(cols)}")
    return t, cols[keys[0]]


# the solver logs a housing case writes, in run order: stage 1 or the single
# run (resumes append to it), then the heat stage of the two-stage path
LOG_NAMES = ("log.chtMultiRegionFoam", "log.chtMultiRegionFoam.heat")


def solver_logs(case):
    """Solver logs in run order (LOG_NAMES). Each is the plain file or, if a
    pull carried only the compressed copy, its '.zst' (read through zstd).
    Other files matching log.chtMultiRegionFoam* (a renamed rejected run,
    a .zst beside its plain log) are not solver logs of the case."""
    case = Path(case)
    out = []
    for name in LOG_NAMES:
        for cand in (case / name, case / f"{name}.zst"):
            if cand.is_file():
                out.append(cand)
                break
    return out


def read_logs(case):
    return [v1_io.read_log(p) for p in solver_logs(case)]


open_text = v1_io.open_text


def minmax(case, region, name):
    """Stitched fieldMinMax output {field: v1_io.MinMax} (location true):
    segment k is truncated at the first time of segment k+1."""
    d = _dirs(case, region, name)
    if d is None:
        raise FileNotFoundError(f"{case}: no postProcessing for {region}/{name}")
    segs = []
    for sd in d.iterdir():
        try:
            label = float(sd.name)
        except ValueError:
            continue
        f = _segment_file(sd, name)
        if f is not None:
            segs.append((label, v1_io.read_field_minmax(f)))
    segs.sort(key=lambda x: x[0])
    if not segs:
        raise FileNotFoundError(f"{case}: {region}/{name} has no .dat segment")
    out = {}
    for fld in dict.fromkeys(k for _, s in segs for k in s):
        parts = []
        for k, (_, s) in enumerate(segs):
            if fld not in s or s[fld].t.size == 0:
                continue
            m = s[fld]
            keep = np.ones(m.t.size, bool)
            nxt = next((segs[j][1][fld].t[0] for j in range(k + 1, len(segs))
                        if fld in segs[j][1] and segs[j][1][fld].t.size), None)
            if nxt is not None:
                keep = m.t < nxt - 1e-12
            parts.append(m if keep.all() else v1_io.MinMax(
                m.t[keep], m.min[keep], m.min_loc[keep], m.max[keep],
                m.max_loc[keep]))
        out[fld] = v1_io.MinMax(*(np.concatenate([getattr(p, a) for p in parts])
                                  for a in ("t", "min", "min_loc", "max",
                                            "max_loc")))
    return out
