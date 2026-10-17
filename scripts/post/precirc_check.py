#!/usr/bin/env python3
"""NEW-14 stationarity of a housing pre-circulation (WP 2.5c).

Rule (plan NEW-14 as amended by PROJECT.md §2.5 NEW-14.1, 2026-09-28): the
isothermal pre-circulated flow is stationary when the port pressure-
difference amplitude (steady: its mean) and the mean fluid kinetic energy
change by < 1 % between consecutive windows over the last 3 changes (4
complete windows). A window is the smallest whole number of periods that
lasts >= 10 s (steady: 10 s), so a short period cannot hide a slow drift.
The default pre-run is 120 s (pulsed: whole periods >= max(120 s, ramp + 5
periods)); a claim-bearing state that fails is re-run once with t_pre
doubled, and a second failure is a hard failure (the class is reported).

Implementation, on the monitors of [0, t_pre] (housing_physics.monitors):
  dp*(t) = p*_in - p*_out, p* = p - rho_in gh (D-P), from the every-step
           areaAverage(p) and areaAverage(gh) of each port (exact: p* is
           linear in both)
  KE(t)  = volIntegrate(rho |U|^2 / 2) over the fluid (KE_fluid)
Windows: the smallest whole number of periods >= 10 s (pulsed) or 10 s
(steady). Per window: dp* amplitude (max - min; steady: the mean) and the
mean KE. The last 3
window-to-window relative changes of both must be < 1 % (4 complete
windows needed). A metric below its no-flow floor passes with a note.

Exit 0 = stationary, 1 = not stationary (or too short), 2 = input error;
--report-only always exits 0 after printing (smokes).

Usage: precirc_check.py CASE [--report-only] [--json OUT]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import housing_io as HIO  # noqa: E402

REL_TOL = 0.01          # NEW-14: < 1 % change
N_CHANGES = 3           # over the last 3 consecutive window-to-window changes
MIN_WINDOW_S = 10.0     # NEW-14.1: windows last >= 10 s (steady: exactly 10 s)
STEADY_WINDOW_S = MIN_WINDOW_S
DOUBLE_ONCE = 2.0       # NEW-14.1: a failed state is re-run once at 2 x t_pre


def window_length(period):
    """NEW-14.1 window: the smallest whole number of periods lasting >=
    MIN_WINDOW_S (steady: MIN_WINDOW_S)."""
    if not period:
        return MIN_WINDOW_S
    return math.ceil(MIN_WINDOW_S / period - 1e-9) * period
FLOOR_DP = 1e-9         # Pa: no-flow floor of the dp* metric
FLOOR_KE = 1e-15        # J: no-flow floor of the KE metric


def port_pstar(case, meta, port):
    t, cols = HIO.series(case, "fluid", f"areaAverage_{port}")
    p = next(v for k, v in cols.items() if k == "areaAverage(p)")
    gh = next(v for k, v in cols.items() if k == "areaAverage(gh)")
    return t, p - meta["rho_inlet"] * gh


def window_metrics(t, y, t0, t1, win, kind):
    """[(t_end, metric)] over complete windows in [t0, t1]."""
    out = []
    n = int(math.floor((t1 - t0) / win + 1e-9))
    for k in range(n):
        a, b = t0 + k * win, t0 + (k + 1) * win
        m = (t > a + 1e-12) & (t <= b + 1e-12)
        if m.sum() < 2:
            continue
        yy = y[m]
        out.append((b, float(np.max(yy) - np.min(yy)) if kind == "amp"
                    else float(np.mean(yy))))
    return out


def rel_changes(ws, floor):
    ch = []
    for (_, a), (tb, b) in zip(ws[:-1], ws[1:]):
        ref = max(abs(b), abs(a))
        ch.append((tb, 0.0 if ref < floor else abs(b - a) / ref))
    return ch


def check(case):
    meta = HIO.case_meta(case)
    t_pre = float(meta["t_pre"])
    P = meta.get("period_s")
    win = window_length(P)
    ti, pi = port_pstar(case, meta, "inlet")
    to, po = port_pstar(case, meta, "outlet")
    n = min(ti.size, to.size)
    if n == 0 or np.any(np.abs(ti[:n] - to[:n]) > 1e-9 * max(1.0, t_pre)):
        raise ValueError("port monitors are not on one time base")
    t, dp = ti[:n], pi[:n] - po[:n]
    tk, ke = HIO.column(case, "fluid", "KE_fluid")
    sel, selk = t <= t_pre + 1e-9, tk <= t_pre + 1e-9
    kind = "amp" if P else "mean"
    wdp = window_metrics(t[sel], dp[sel], 0.0, t_pre, win, kind)
    wke = window_metrics(tk[selk], ke[selk], 0.0, t_pre, win, "mean")
    cdp, cke = rel_changes(wdp, FLOOR_DP), rel_changes(wke, FLOOR_KE)
    notes = []
    enough = len(cdp) >= N_CHANGES and len(cke) >= N_CHANGES
    if not enough:
        notes.append(f"only {len(wdp)} complete {win:g} s window(s) in the "
                     f"{t_pre:g} s pre-run; NEW-14 needs {N_CHANGES + 1}")
    last_dp = [c for _, c in cdp[-N_CHANGES:]]
    last_ke = [c for _, c in cke[-N_CHANGES:]]
    ok = enough and max(last_dp + last_ke, default=math.inf) < REL_TOL
    if wdp and max(abs(v) for _, v in wdp) < FLOOR_DP:
        notes.append("dp* below the no-flow floor")
    if not ok:
        rerun = DOUBLE_ONCE * t_pre
        if P:
            rerun = math.ceil(rerun / P - 1e-9) * P
        notes.append(f"NOT stationary: re-run the pre-circulation once at "
                     f"t_pre = {rerun:g} s (NEW-14.1); a second failure is a "
                     "hard failure, reported")
    return dict(case=str(case), t_pre=t_pre, window_s=win, metric_dp=kind,
                windows_dp=wdp, windows_ke=wke,
                last_changes_dp=last_dp, last_changes_ke=last_ke,
                rel_tol=REL_TOL, stationary=bool(ok), notes=notes)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case", type=Path)
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args(argv)
    try:
        r = check(a.case)
    except (OSError, KeyError, ValueError, StopIteration) as e:
        print(f"precirc_check: input error: {e!r}", file=sys.stderr)
        return 0 if a.report_only else 2
    if a.json:
        a.json.write_text(json.dumps(r, indent=1) + "\n")
    tag = "STATIONARY" if r["stationary"] else "NOT STATIONARY"
    print(f"precirc_check {tag}: t_pre {r['t_pre']:g} s, windows "
          f"{r['window_s']:g} s; last dp* changes "
          f"{[round(100 * c, 3) for c in r['last_changes_dp']]} %, KE "
          f"{[round(100 * c, 3) for c in r['last_changes_ke']]} % "
          f"(tol {100 * REL_TOL:g} %)")
    for n in r["notes"]:
        print("  note:", n)
    if a.report_only:
        return 0
    return 0 if r["stationary"] else 1


if __name__ == "__main__":
    sys.exit(main())
