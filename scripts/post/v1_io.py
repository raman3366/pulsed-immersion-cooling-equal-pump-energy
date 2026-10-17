#!/usr/bin/env python3
"""Readers for the v1 (AWS, 2026-09) raw OpenFOAM outputs (WP 2.9).

Used by v1_rescore.py and v1_fields.py; unit-tested in tests/post/test_v1_io.py.

Everything here is header-aware, so a vector column '(x y z)' is ONE column
(the v1 extract_metrics.py split it into three tokens, which is how
Tmax_fluid / magUmax_fluid ended up holding the y coordinate of the minimum
location - erratum item E-3).

- read_fo(path)              any function-object .dat / probes file
- read_field_minmax(path)    fieldMinMax.dat -> per field: min, max, both locations
- read_log(path)             solver log (full or the v1 'slim' log): step times,
                             deltaT, Courant max, run start times, header date
- list_segments(case_dir)    every postProcessing/<fo>/<startTime>/<file>
- run_start_index()          OpenFOAM time index at the start of a run (restarts)
- expected_sample_times()    step times at which an every-N-steps FO writes
"""

from __future__ import annotations

import contextlib
import math
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_TOKEN = re.compile(r"\([^()]*\)|\S+")
_NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
_META = re.compile(r"^#\s*([A-Za-z][\w ()]*?)\s*:\s*(.*?)\s*$")


@dataclass
class FOData:
    path: Path
    meta: dict                  # '# Key : value' header lines (Cells, Faces, Volume, Area, ...)
    columns: list               # column names (duplicates made unique by position)
    data: dict                  # name -> np.ndarray (n,) float, (n, 3) vector, or list[str]

    @property
    def t(self) -> np.ndarray:
        return self.data[self.columns[0]]

    def __len__(self):
        return len(self.t)


def _header_columns(line: str) -> list:
    """Column names of a '#' header line. A repeated name (fieldMinMax writes
    'processor' after each location) takes the qualifier of the column before
    it: 'location(min) processor' -> 'processor(min)'."""
    names = line.lstrip("#").split()
    out = []
    for i, n in enumerate(names):
        if names.count(n) > 1:
            m = re.match(r".*\((.*)\)$", names[i - 1]) if i else None
            n = f"{n}({m.group(1)})" if m else f"{n}#{i}"
        out.append(n)
    return out


def _parse_token(tok: str):
    if tok.startswith("("):
        return [float(x) for x in tok[1:-1].split()]
    try:
        return float(tok)
    except ValueError:
        return tok


def read_fo(path) -> FOData:
    """Read a function-object output file.

    The column names come from the last '#' line before the data; a data
    token is either '( ... )' (one vector column) or a bare word. Scalar
    columns become float arrays, vector columns (n, 3) arrays, word columns
    lists of str. A row whose token count differs from the header is an
    error (it is never silently re-aligned)."""
    path = Path(path)
    meta, header, rows = {}, None, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = _META.match(s)
                if m and rows == []:
                    meta[m.group(1).strip()] = m.group(2)
                header = s
                continue
            rows.append(s)
    if header is None:
        raise ValueError(f"{path}: no header line")
    cols = _header_columns(header)
    if not rows:
        return FOData(path, meta, cols, {c: np.array([]) for c in cols})
    simple = "(" not in rows[0] and all(_NUM.match(t) for t in rows[0].split())
    if simple:
        flat = " ".join(rows).split()
        if len(flat) != len(rows) * len(cols):
            # fall back to the checked per-row path to report the bad row
            simple = False
        else:
            arr = np.array(flat, dtype=float).reshape(len(rows), len(cols))
            return FOData(path, meta, cols, {c: arr[:, j] for j, c in enumerate(cols)})
    parsed = []
    for r in rows:
        toks = _TOKEN.findall(r)
        if len(toks) != len(cols):
            raise ValueError(f"{path}: row has {len(toks)} tokens, header "
                             f"{len(cols)} ({cols}): {r!r}")
        parsed.append([_parse_token(t) for t in toks])
    data = {}
    for j, c in enumerate(cols):
        vals = [p[j] for p in parsed]
        if all(isinstance(v, float) for v in vals):
            data[c] = np.array(vals, dtype=float)
        elif all(isinstance(v, list) for v in vals):
            data[c] = np.array(vals, dtype=float)
        else:
            data[c] = [v if isinstance(v, str) else repr(v) for v in vals]
    return FOData(path, meta, cols, data)


@dataclass
class MinMax:
    t: np.ndarray
    min: np.ndarray
    min_loc: np.ndarray        # (n, 3)
    max: np.ndarray
    max_loc: np.ndarray        # (n, 3)


def read_field_minmax(path) -> dict:
    """fieldMinMax.dat -> {field name: MinMax}. Header of v2606:
    Time field min location(min) [processor] max location(max) [processor]."""
    fo = read_fo(path)
    need = ("field", "min", "location(min)", "max", "location(max)")
    miss = [c for c in need if c not in fo.data]
    if miss:
        raise ValueError(f"{path}: not a fieldMinMax file (missing {miss}; "
                         f"columns {fo.columns})")
    names = fo.data["field"]
    out = {}
    for fld in dict.fromkeys(names):
        m = np.array([n == fld for n in names])
        out[fld] = MinMax(fo.t[m], fo.data["min"][m], fo.data["location(min)"][m],
                          fo.data["max"][m], fo.data["location(max)"][m])
    return out


# ------------------------------------------------------------------ logs

_TIME_LINE = re.compile(r"^Time = ([-+0-9.eE]+)\s*$")
_DT_LINE = re.compile(r"^deltaT = ([-+0-9.eE]+)\s*$")
_CO_LINE = re.compile(r"Courant Number mean:\s*\S+\s+max:\s*([-+0-9.eE]+)")
_DATE = re.compile(r"^Date\s*:\s*(.+?)\s*$")


@dataclass
class LogInfo:
    path: Path
    step_times: np.ndarray          # every 'Time = ' value, in log order
    step_dt: np.ndarray             # deltaT of each step (last 'deltaT =' before it)
    courant_max: np.ndarray         # every 'Courant Number ... max:' value
    run_starts: list                # start time of every run in the log
    run_first_index: list           # index into step_times of each run's first step
    date: str | None = None
    ended: bool = False

    @property
    def end_time(self) -> float:
        return float(self.step_times[-1])


@contextlib.contextmanager
def open_text(path):
    """Text lines of a (possibly zstd-compressed) log: a '.zst' file is
    streamed through 'zstd -dc' (the run driver compresses the solver
    logs, and a pull may carry only the .zst)."""
    path = Path(path)
    if path.suffix != ".zst":
        with open(path, errors="replace") as f:
            yield f
        return
    proc = subprocess.Popen(["zstd", "-q", "-dc", "--", str(path)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, errors="replace")
    try:
        yield proc.stdout
        while proc.stdout.read(1 << 20):    # drain before wait (no deadlock)
            pass
    finally:
        proc.stdout.close()
        err = proc.stderr.read()
        proc.stderr.close()
        rc = proc.wait()
    if rc != 0:
        raise ValueError(f"{path}: zstd -dc failed ({rc}): {err.strip()}")


def read_log(path) -> LogInfo:
    """Solver log (a full chtMultiRegionFoam log or the v1 'slim' log that
    keeps only the Courant / deltaT / Time lines; '.zst' is decompressed on
    the fly). A new run starts where a time does not increase (a restart
    from an earlier written time appended to the same log) or where a
    solver banner's 'Create time' follows earlier steps (a resume appended
    after a clean stop at a written time, whose times keep increasing)."""
    path = Path(path)
    times, dts, co = [], [], []
    creates = []
    last_dt = math.nan
    date = None
    ended = False
    with open_text(path) as f:
        for line in f:
            if line.startswith("Create time"):
                creates.append(len(times))
                continue
            if line.startswith("Time = "):
                m = _TIME_LINE.match(line)
                if m:
                    times.append(float(m.group(1)))
                    dts.append(last_dt)
                continue
            if line.startswith("deltaT = "):
                m = _DT_LINE.match(line)
                if m:
                    last_dt = float(m.group(1))
                continue
            if "Courant Number" in line:
                m = _CO_LINE.search(line)
                if m:
                    co.append(float(m.group(1)))
                continue
            if date is None and line.startswith("Date"):
                m = _DATE.match(line)
                if m:
                    date = m.group(1)
            elif line.strip() == "End":
                ended = True
    t = np.array(times, dtype=float)
    d = np.array(dts, dtype=float)
    if t.size == 0:
        raise ValueError(f"{path}: no 'Time =' lines")
    firsts = sorted({0, *(int(i) + 1 for i in np.nonzero(np.diff(t) <= 0)[0]),
                     *(c for c in creates if 0 < c < t.size)})
    starts = []
    for i in firsts:
        dt0 = d[i] if math.isfinite(d[i]) else (t[i + 1] - t[i] if i + 1 < t.size else 0.0)
        starts.append(float(t[i] - dt0))
    return LogInfo(path, t, d, np.array(co, dtype=float), starts, firsts, date, ended)


def run_start_index(log: LogInfo, run: int, rtol: float = 1e-7) -> int | None:
    """OpenFOAM time index at the start of one run of the log (timeIndex;
    a restart restores it from <time>/uniform/time 'index', Time.C).

    A run from t = 0 starts at index 0. A later run restarts from a time
    an earlier run of the same log wrote: its start index is that run's
    start index plus the number of steps it took to reach the restart time.
    None if neither holds (a log that begins with a restart, or a restart
    time that is not a step of an earlier run in this log): the index is
    then unknown."""
    k = run % len(log.run_starts)
    if abs(log.run_starts[k]) <= 1e-9:
        return 0
    if k == 0:
        return None
    t_r = log.run_starts[k]
    for j in range(k - 1, -1, -1):
        steps = log.step_times[log.run_first_index[j]:log.run_first_index[j + 1]]
        hit = np.nonzero(np.abs(steps - t_r) <= rtol * max(1.0, abs(t_r)) + 1e-12)[0]
        if hit.size:
            base = run_start_index(log, j, rtol)
            return None if base is None else base + int(hit[0]) + 1
    return None


def expected_sample_times(log: LogInfo, every: int, run: int = -1,
                          start_index: int | None = None) -> np.ndarray:
    """Step times at which an FO with writeControl timeStep / writeInterval
    `every` writes during one run of the log.

    OpenFOAM v2606 fires such an FO when the ABSOLUTE time index is a
    multiple of `every` (timeControl.C: !(time_.timeIndex() % intInterval_)),
    and a restart restores the index of the restart time (Time.C reads
    'index' from <time>/uniform/time). A restart at an index that is not a
    multiple of `every` therefore writes first after fewer than `every`
    steps (both v1 field re-runs restarted at 870 s write first at 870.025 s,
    their 2nd step). start_index: the time index at the start of the run;
    default run_start_index() (0 for the first run of the log). ValueError
    if that index is unknown."""
    k = run % len(log.run_starts)
    i0 = log.run_first_index[k]
    i1 = (log.run_first_index[k + 1] if k + 1 < len(log.run_first_index)
          else log.step_times.size)
    steps = log.step_times[i0:i1]
    if start_index is None:
        start_index = run_start_index(log, k)
        if start_index is None:
            raise ValueError(f"{log.path}: time index at the start of run {k} (t = "
                             f"{log.run_starts[k]:g} s) is unknown; pass start_index")
    idx = start_index + np.arange(1, steps.size + 1)
    return steps[idx % every == 0]


def match_times(sample_t, ref_t, rtol: float = 1e-7) -> np.ndarray:
    """Boolean mask: which sample times appear in ref_t (to print precision)."""
    s = np.asarray(sample_t, dtype=float)
    r = np.sort(np.asarray(ref_t, dtype=float))
    if r.size == 0:
        return np.zeros(s.size, bool)
    j = np.clip(np.searchsorted(r, s), 1, r.size - 1)
    near = np.minimum(np.abs(r[j] - s), np.abs(r[j - 1] - s))
    return near <= rtol * np.maximum(1.0, np.abs(s)) + 1e-12


# ------------------------------------------------------------------ segments

def _is_number(name: str) -> bool:
    return bool(_NUM.match(name))


@dataclass
class Segment:
    fo: str                     # e.g. 'fluid/TvolAverage_fluid', 'T2probe/cell1'
    label: float                # the <startTime> directory
    label_str: str
    path: Path
    mtime: float
    kind: str = "monitor"       # 'monitor' | 'sampler' (surfaces etc., not a series)
    extra: dict = field(default_factory=dict)


def list_segments(case_dir) -> list:
    """Every postProcessing/<fo...>/<startTime>/<file>. Dot files (macOS
    .DS_Store, AppleDouble) are ignored; files not at <numeric>/<file> are
    reported with kind 'misplaced'."""
    pp = Path(case_dir) / "postProcessing"
    out = []
    if not pp.is_dir():
        return out
    for dirpath, dirnames, filenames in os.walk(pp):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for fn in sorted(filenames):
            if fn.startswith("."):
                continue
            p = Path(dirpath) / fn
            rel = p.relative_to(pp)
            parts = rel.parts
            if len(parts) >= 3 and _is_number(parts[-2]):
                fo = "/".join(parts[:-2])
                kind = "monitor" if (fn.endswith(".dat") or fn in ("T", "U", "p")) else "sampler"
                out.append(Segment(fo, float(parts[-2]), parts[-2], p, p.stat().st_mtime, kind))
            elif len(parts) >= 2 and _is_number(parts[0]):
                # postProcessing/<time>/<surfaces...> (sampled surfaces)
                out.append(Segment(parts[0] if len(parts) == 2 else "/".join(parts[1:-1]),
                                   float(parts[0]), parts[0], p, p.stat().st_mtime, "sampler"))
            else:
                out.append(Segment("/".join(parts[:-1]), math.nan, "", p,
                                   p.stat().st_mtime, "misplaced"))
    return out


def fo_every(fo_path: str) -> int | None:
    """writeInterval (in steps) of the v1 monitors (system/monitors of the
    v1 template, generated by makePhysics.py): volFieldValue / fieldMinMax /
    probes / wallLossInt every 10 steps, inlet/outlet surfaceFieldValue every
    5, the wallHeatFlux FO writes its .dat every step. None = unknown."""
    name = fo_path.split("/")[-1] if not fo_path.startswith("T2probe") else "T2probe"
    if name in ("sum_inlet", "sum_outlet", "areaAverage_inlet", "areaAverage_outlet",
                "TphiSum_inlet", "TphiSum_outlet"):
        return 5
    if name == "wallLoss_fluid":
        return 1
    if (name.startswith(("Tmax_cell", "Tmin_cell", "TvolAverage_", "rhoh_"))
            or name in ("TfluidMinMax", "wallLossInt_fluid", "T2probe")):
        return 10
    return None
