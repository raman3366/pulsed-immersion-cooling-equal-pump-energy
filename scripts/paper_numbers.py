#!/usr/bin/env python3
"""Derived numbers of the paper from a reproduction: the pumping-energy
ratios, the numerical checks, the matched-energy and mechanism comparisons and
the equal-energy histories (the paper's Tables I-III and Figs. 3-5).

    python3 scripts/paper_numbers.py [RESULTS_DIR]

RESULTS_DIR is the AISEE_OUT directory of a scorer run (default work/results;
it holds aisee_results.csv and the gate records aisee_cache/<case>.json), and
AISEE_CASE_ROOT the extracted case directories (default work/cases), as for
scripts/post/aisee_results.py. Prints:

1. per case, from the gate records: enthalpy-ledger residual, mean-flow
   error, Courant maxima by region (claim-region bulk, thin layer, outflow
   column; reported, not gated), liquid minimum below the inlet temperature,
   isothermal port-work budget residual;
2. energy ratios (Tables I and III): r_QS (viscous-limit pumping energy relative to steady
   90 mL/min, (<Q^2>/Qbar^2)(Qbar/90)^2) and r (measured E_pump relative to
   steady 90 mL/min);
3. matched pumping energy: each pulsed case against the steady curve
   interpolated linearly in log E_pump (peak temperature, end spread, and the
   equivalent steady flow by log-log interpolation of Qbar against E_pump);
4. the mechanism comparison at 90 mL/min (pulsed minus steady), the intervals
   in which each pulsed case (and C1 at dt/2) is less uniform than steady flow
   (1 s grid after heat-on; largest excess) and the row means at the end
   (rows are the cell pairs (0,3), (1,4), (2,5), upstream first);
5. numerical checks (Table II): the C1 time-step pair at t = 330 s (aisee_results.at_time),
   the 360 mL/min time-step and limited-gradient pairs, and the three meshes at
   t = 300 s with the grid convergence index of the spread (Celik et al.,
   J. Fluids Eng. 130 (2008) 078001, as implemented in scripts/post/gci2.py;
   cell size (1/N)^(1/3) from the cell counts in each case's k0_mesh.json; the
   GCI of the spread on the fine and the medium mesh, and the observed order of
   T_max); and the largest change of the C1 spread at dt/2 over 100-250 s after
   heat-on (1 s grid);
6. the equal-energy histories (Fig. 4): each pulsed case against the steady
   histories interpolated linearly in log E_pump between the two bracketing
   steady flows at each time (aisee_figures.equal_energy_steady): the largest
   T_max difference over the first 200 s after heat-on, its value at the end
   and whether it grows monotonically after 200 s (10 s grid); the range of the
   relative change of the spread from 60 s; and the intervals in which the
   pulsed case is less uniform, with the largest excess.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts" / "post"))
import aisee_results as AR      # noqa: E402
import gci2                     # noqa: E402

T_IN = 298.15
TAU_V = 0.05                     # valve edge (s)
QS = {"steady": lambda: 1.0,     # <Q^2>/Qbar^2, paper Eq. (1), "Flow Schedules and Cases"
      "C3": lambda: 1 + 0.6 ** 2 / 2,          # sine, A 0.6
      "C1": lambda: 1 / 0.25 - TAU_V * 1.0 / (3 * 0.25 ** 2),
      "C2p": lambda: 1 / 0.5 - TAU_V * 0.5 / (3 * 0.5 ** 2)}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    res = Path(argv[0]) if argv else AR.OUT
    rows = {r["key"]: r for r in csv.DictReader((res / "aisee_results.csv").open())}
    cs = {c["key"]: c for c in AR.cases()}
    f = lambda k, c: float(rows[k][c])    # noqa: E731

    print("== 1. gate records (Courant: reported, not gated)")
    for k in rows:
        d = {x["name"]: x["detail"] for x in json.loads(
            (res / "aisee_cache" / f"{cs[k]['path'].name}.json").read_text())["gates"]}
        led = re.search(r"= ([+-][\d.]+) % of E_inj", d.get("ledger", ""))
        mf = re.search(r"\(([+-][\d.]+) %", d.get("mean_flow", ""))
        co = re.search(r"bulk max ([\d.]+).*thin layer max ([\d.]+).*outflow column max ([\d.]+)",
                       d.get("courant", ""))
        cog = re.search(r"Courant max ([\d.]+)", d.get("courant", ""))
        mp = re.search(r"in \[([\d.]+), ([\d.]+)\]", d.get("max_principle", ""))
        cl = re.search(r"= ([+-][\d.]+) % of", d.get("hydraulic_closure[pre]", ""))
        cot = (f"bulk {co.group(1)} thin {co.group(2)} out {co.group(3)}" if co
               else f"global {cog.group(1)}" if cog else "-")
        print(f"{k:24s} ledger {led.group(1) if led else '-':>7} %  mean flow "
              f"{mf.group(1) if mf else '-':>7} %  Co {cot:30s} liquid min - T_in "
              f"{float(mp.group(1)) - T_IN if mp else float('nan'):+.3f} K  budget "
              f"{cl.group(1) + ' %' if cl else '-'}")

    print("\n== 2. energy ratios (r_QS: viscous-limit estimate; r: measured)")
    e90 = f("steady_Q90", "E_pump_mJ")
    for k, r in rows.items():
        if r["role"] not in ("curve", "mean-flow", "matched-energy"):
            continue
        q = float(r["Qbar_mL_min"])
        print(f"{k:12s} r_QS {QS[r['arm']]() * (q / 90) ** 2:6.2f}   r {f(k, 'E_pump_mJ') / e90:6.2f}"
              f"   E_pump {f(k, 'E_pump_mJ'):.4g} mJ")

    print("\n== 3. matched pumping energy (steady curve interpolated in log E_pump)")
    for q360 in ("steady_Q360", "steady_Q360_cand2"):
        curve = ["steady_Q90", "steady_Q135", "steady_Q180", q360]
        le = np.log([f(c, "E_pump_mJ") for c in curve])
        print(f"-- steady curve with {q360}")
        for k in ("C3_Q82.9", "C3_Q90", "C2p_Q64.2", "C1_Q46.6", "C1_Q90"):
            e = f(k, "E_pump_mJ")
            tm = float(np.interp(np.log(e), le, [f(c, "Tmax_peak_C") for c in curve]))
            dt = float(np.interp(np.log(e), le, [f(c, "dTcell_end_K") for c in curve]))
            qe = float(np.exp(np.interp(np.log(e), le, np.log([90, 135, 180, 360]))))
            T, D = f(k, "Tmax_peak_C"), f(k, "dTcell_end_K")
            print(f"{k:10s} E {e:.3f} mJ ~ steady {qe:.0f} mL/min: T_max,peak {T:.3f} vs {tm:.3f} C "
                  f"({T - tm:+.2f} K); dT_cell,end {1e3 * D:.1f} vs {1e3 * dt:.1f} mK "
                  f"({100 * (D / dt - 1):+.1f} %)")

    print("\n== 4. mechanism at 90 mL/min (pulsed minus steady)")
    for k in ("C1_Q90", "C3_Q90"):
        print(f"{k}: T_max,peak {f(k, 'Tmax_peak_C') - f('steady_Q90', 'Tmax_peak_C'):+.3f} K, "
              f"T_module,end {f(k, 'Tmodule_end_C') - f('steady_Q90', 'Tmodule_end_C'):+.3f} K, "
              f"dT_cell,end {100 * (f(k, 'dTcell_end_K') / f('steady_Q90', 'dTcell_end_K') - 1):+.1f} %, "
              f"time-mean dT_cell {100 * (f(k, 'dTcell_mean_K') / f('steady_Q90', 'dTcell_mean_K') - 1):+.1f} %, "
              f"E_pump x{f(k, 'E_pump_mJ') / e90:.2f}")
    grid = np.arange(0.0, 900.0 + 0.5, 1.0)          # s after heat-on

    def on_grid(k):
        t, smax, dtc, _ = AR.series(cs[k]["path"])
        return np.interp(grid + 30.0, t, dtc)

    def spans(excess):
        out, i = [], 0
        while i < excess.size:
            if excess[i] > 1e-6:                     # above 1 uK (round-off at heat-on)
                j = i
                while j + 1 < excess.size and excess[j + 1] > 1e-6:
                    j += 1
                out.append(f"{grid[i]:.0f}-{grid[j]:.0f} s ({1e3 * excess[i:j + 1].max():+.2f} mK)")
                i = j + 1
            else:
                i += 1
        return ", ".join(out) or "never"

    d90 = on_grid("steady_Q90")
    for k in ("C1_Q90", "C3_Q90", "C1_Q90_dt6.25_h300"):
        d = on_grid(k)
        n = int(np.searchsorted(grid, 300.0, "right")) if k.endswith("h300") else grid.size
        print(f"{k} less uniform than steady 90: {spans((d - d90)[:n])}")
    for k in ("steady_Q90", "C1_Q90", "C3_Q90"):
        t, _, _, _ = AR.series(cs[k]["path"])
        tav = [AR.HIO.column(cs[k]["path"], f"cell{i}", f"TvolAverage_cell{i}")[1][-1] for i in range(6)]
        rows3 = [(tav[i] + tav[i + 3]) / 2 for i in range(3)]
        print(f"{k} end row means above the upstream row: middle {1e3 * (rows3[1] - rows3[0]):.1f} mK, "
              f"downstream {1e3 * (rows3[2] - rows3[0]):.1f} mK")

    print("\n== 5. numerical checks")
    a = AR.at_time(cs["C1_Q90"]["path"], 330.0)
    b = AR.at_time(cs["C1_Q90_dt6.25_h300"]["path"], 330.0)
    print(f"C1 90, dt 12.5 -> 6.25 ms at t = {b['t']:g} s: T_max {a['Tmax_C']:.4f} -> {b['Tmax_C']:.4f} C "
          f"({1e3 * (b['Tmax_C'] - a['Tmax_C']):+.1f} mK), dT_cell {1e3 * a['dTcell_K']:.2f} -> "
          f"{1e3 * b['dTcell_K']:.2f} mK ({100 * (b['dTcell_K'] / a['dTcell_K'] - 1):+.1f} %), E_pump "
          f"{f('C1_Q90', 'E_pump_mJ'):.4f} -> {f('C1_Q90_dt6.25_h300', 'E_pump_mJ'):.4f} mJ "
          f"({100 * (f('C1_Q90_dt6.25_h300', 'E_pump_mJ') / f('C1_Q90', 'E_pump_mJ') - 1):+.1f} %)")
    for k, what in (("steady_Q360_dt6.25", "dt 6.25 ms"), ("steady_Q360_cand2", "limited gradients")):
        print(f"steady 360, {what}: T_max,peak {1e3 * (f(k, 'Tmax_peak_C') - f('steady_Q360', 'Tmax_peak_C')):+.1f} mK, "
              f"dT_cell,end {100 * (f(k, 'dTcell_end_K') / f('steady_Q360', 'dTcell_end_K') - 1):+.1f} %")
    keys = ["steady_Q90_coarse_h270", "steady_Q90", "steady_Q90_fine_h270"]
    vals = [AR.at_time(cs[k]["path"], 300.0) for k in keys]
    n = [json.loads((cs[k]["path"] / "k0_mesh.json").read_text())["cells_total"] for k in keys]
    print("three meshes, steady 90 at t = 300 s (coarse / medium / fine; "
          f"{n[0]:,} / {n[1]:,} / {n[2]:,} cells):")
    for lab, key, sc in (("T_max (C)", "Tmax_C", 1), ("T_module (C)", "Tmodule_C", 1),
                         ("dT_cell (mK)", "dTcell_K", 1e3)):
        print(f"  {lab:13s} " + " / ".join(f"{sc * v[key]:.4f}" for v in vals))
    d = [1e3 * v["dTcell_K"] for v in vals]           # mK, coarse / medium / fine
    g = gci2.celik_gci([d[2], d[1], d[0]], [gci2.cell_size(x) for x in (n[2], n[1], n[0])],
                       rise=d[2])
    print(f"  r21 {g.r21:.3f}, r32 {g.r32:.3f}; dT_cell ({g.convergence}): observed order "
          f"{g.p_observed:.2f}, extrapolated {g.phi_ext:.1f} mK, GCI_fine "
          f"{100 * g.U_fine / d[2]:.1f} %, GCI_medium {100 * g.U_medium / d[1]:.1f} % (Fs {g.fs:g})")
    tm = [v["Tmax_C"] for v in vals]
    gt = gci2.celik_gci([tm[2], tm[1], tm[0]], [gci2.cell_size(x) for x in (n[2], n[1], n[0])],
                        rise=tm[2] - (T_IN - 273.15))
    # Celik's fixed-point order on T_max regardless of the class gci2 assigns
    # (e21 and e32 nearly equal: no asymptotic range, so no extrapolation)
    e21, e32 = tm[1] - tm[2], tm[0] - tm[1]
    sgn = np.sign(e32 / e21)
    pp = 1.0
    for _ in range(200):
        q = np.log((gt.r21 ** pp - sgn) / (gt.r32 ** pp - sgn))
        pp = abs(np.log(abs(e32 / e21)) + q) / np.log(gt.r21)
    print(f"  T_max: changes per level {1e3 * e32:.1f} / {1e3 * e21:.1f} mK (gci2: {gt.convergence}); "
          f"Celik fixed-point order {pp:.2f}: no asymptotic range")
    g1 = on_grid("C1_Q90")
    g2 = on_grid("C1_Q90_dt6.25_h300")
    w = (grid >= 100.0) & (grid <= 250.0)
    dd = g2[w] - g1[w]
    i = int(np.argmax(np.abs(dd)))
    print(f"  C1 90 spread at dt/2 over 100-250 s: largest change {1e3 * dd[i]:+.1f} mK "
          f"({100 * dd[i] / g1[w][i]:+.1f} %) at {grid[w][i]:.0f} s")
    print(f"  T_max change per level: {1e3 * (vals[1]['Tmax_C'] - vals[0]['Tmax_C']):+.1f}, "
          f"{1e3 * (vals[2]['Tmax_C'] - vals[1]['Tmax_C']):+.1f} mK; E_pump "
          + " / ".join(f"{f(k, 'E_pump_mJ'):.4f}" for k in keys) + " mJ")

    print("\n== 6. equal-energy histories (t: time after heat-on)")
    import aisee_figures as AF
    frows = [dict(r, E_pump_mJ=float(r["E_pump_mJ"])) for r in rows.values()]
    cache = {}

    def series_of(k):
        if k not in cache:
            t, smax, dtc, _ = AR.series(cs[k]["path"])
            cache[k] = (t, smax, dtc)
        return cache[k]

    for k in ("C1_Q46.6", "C2p_Q64.2", "C3_Q82.9", "C1_Q90", "C3_Q90"):
        t, smax, dtc = series_of(k)
        te, ae, be, (lo, hi, w) = AF.equal_energy_steady(frows, series_of, f(k, "E_pump_mJ"))
        h = t > 30.0
        th = t[h] - 30.0
        dT = smax[h] - np.interp(t[h], te, ae)
        rel = 100 * (dtc[h] / np.interp(t[h], te, be) - 1)
        grid = np.arange(200.0, 901.0, 10.0)
        mono = bool(np.all(np.diff(np.interp(grid, th, dT)) > 0))
        m60 = th >= 60.0
        less = rel[m60] > 0
        edges = np.flatnonzero(np.diff(less.astype(int)))
        spans = []
        for i0, i1 in zip(np.r_[0, edges + 1], np.r_[edges + 1, less.size]):
            if less[i0]:
                seg = rel[m60][i0:i1]
                spans.append(f"{th[m60][i0]:.0f}-{th[m60][i1 - 1]:.0f} s (max {seg.max():+.1f} %)")
        print(f"{k:10s} steady {lo}/{hi} w {w:.3f}: T_max diff |max| over 0-200 s "
              f"{np.abs(dT[th <= 200]).max():.4f} K, at 900 s {dT[-1]:+.2f} K, increasing after "
              f"200 s: {mono}; spread change from 60 s {rel[m60].min():+.1f} to "
              f"{rel[m60].max():+.1f} %; less uniform: {', '.join(spans) or 'never'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
