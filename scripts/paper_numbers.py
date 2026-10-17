#!/usr/bin/env python3
"""Derived numbers of the paper from a reproduction: the ratios of Table I,
the numerical checks of Table II and the comparisons of Section III.

    python3 scripts/paper_numbers.py [RESULTS_DIR]

RESULTS_DIR is the AISEE_OUT directory of a scorer run (default work/results;
it holds aisee_results.csv and the gate records aisee_cache/<case>.json), and
AISEE_CASE_ROOT the extracted case directories (default work/cases), as for
scripts/post/aisee_results.py. Prints:

1. per case, from the gate records: enthalpy-ledger residual, mean-flow
   error, Courant maxima by region (claim-region bulk, thin layer, outflow
   column; reported, not gated), liquid minimum below the inlet temperature,
   isothermal port-work budget residual;
2. Table I ratios: r_QS (viscous-limit pumping energy relative to steady
   90 mL/min, (<Q^2>/Qbar^2)(Qbar/90)^2) and r (measured E_pump relative to
   steady 90 mL/min);
3. matched pumping energy: each pulsed case against the steady curve
   interpolated linearly in log E_pump (peak temperature, end spread, and the
   equivalent steady flow by log-log interpolation of Qbar against E_pump);
4. the mechanism comparison at 90 mL/min (pulsed minus steady);
5. Table II checks: the C1 time-step pair at t = 330 s (aisee_results.at_time),
   the 360 mL/min time-step and limited-gradient pairs, and the three meshes at
   t = 300 s with the grid convergence index of the spread (Celik et al.,
   J. Fluids Eng. 130 (2008) 078001, as implemented in scripts/post/gci2.py;
   cell size (1/N)^(1/3) from the cell counts in each case's k0_mesh.json).
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
QS = {"steady": lambda: 1.0,     # <Q^2>/Qbar^2, paper Eq. (1)
      "C3": lambda: 1 + 0.6 ** 2 / 2,
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

    print("\n== 2. Table I ratios (r_QS: viscous-limit estimate; r: measured)")
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

    print("\n== 5. Table II checks")
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
          f"{100 * g.U_fine / d[2]:.1f} % (Fs {g.fs:g})")
    print(f"  T_max change per level: {1e3 * (vals[1]['Tmax_C'] - vals[0]['Tmax_C']):+.1f}, "
          f"{1e3 * (vals[2]['Tmax_C'] - vals[1]['Tmax_C']):+.1f} mK; E_pump "
          + " / ".join(f"{f(k, 'E_pump_mJ'):.4f}" for k in keys) + " mJ")
    return 0


if __name__ == "__main__":
    sys.exit(main())
