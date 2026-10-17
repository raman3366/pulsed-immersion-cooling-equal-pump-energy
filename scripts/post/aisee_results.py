#!/usr/bin/env python3
"""AISEE 2027 conference study (D-AIS.1): results table and figures.

Scores every case of the study with the registered housing scorer
(housing_metrics.evaluate: the rule-4 hard gates), then derives the
conference quantities:

- T_max,peak: the largest solid temperature of the six cells over the run
  (max over time of max_i Tmax_cell_i), and its rise over the inlet;
- dT_cell,end: the spread of the six cells' volume-averaged temperatures at
  the end of the discharge; T_module,end: their volume-weighted mean;
- E_pump: the forced-flow hydraulic work of the waveform, measured on the
  case's own isothermal pre-circulation window [ramp, t_pre] (the registered
  `hydraulic_pre` record: port work less the kinetic-energy change, over
  whole periods), as a mean power times the 900 s discharge. In the heated
  window buoyancy does more work than the pump (GDT1: E_port < 0, W_b > 0),
  so the heated port work is not a pump energy (housing_metrics' note).

Every case must pass all of the scorer's gates (energy ledger <= 2 %, mean
flow <= 0.5 %, fixed dt, lattice, Courant as D-F.1 reports it, continuity);
a failing case stops the script (rule 4: hard gates) unless it is named in
--allow-fail, and is then flagged in the table, never silently plotted.

    python3 scripts/post/aisee_results.py            # table + figures
    python3 scripts/post/aisee_results.py --no-figs  # table only

Paths (public data package): AISEE_CASE_ROOT is the directory that holds the
extracted case directories (data/cases/*.tar.xz; default work/cases), and
AISEE_OUT the output directory for aisee_results.csv, the gate cache and the
figures (default work/results). scripts/reproduce.sh sets both.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "post"))
import housing_io as HIO          # noqa: E402
import housing_metrics as HM      # noqa: E402

ARCHIVE = Path(os.environ.get("AISEE_CASE_ROOT", REPO / "work/cases"))
OUT = Path(os.environ.get("AISEE_OUT", REPO / "work/results"))
OUT_CSV = OUT / "aisee_results.csv"
CACHE = OUT / "aisee_cache"
FIGDIR = OUT / "figures"
HEAT_S = 900.0
N_CELLS = 6
NOISE_TMAX_K = 0.004      # GDT1: dt pair 0.0038 K, seed pair 0.0014 K (volume T_max,peak)
NOISE_DTCELL_K = 0.0004   # GDT1: dt and seed pairs, dT_cell,end

# (key, role, arm, Qbar mL/min, grid, dt ms, case directory candidates)
REUSED = [
    ("steady_Q90", "curve", "steady", 90, "medium", 12.5,
     [ARCHIVE / "h_dtladder_steady_Q90_medium_dt12.5"]),
    ("C1_Q90", "mean-flow", "C1", 90, "medium", 12.5,
     [ARCHIVE / "h_dtladder_sq1D0.25_Q90_medium_dt12.5"]),
    ("C3_Q90", "mean-flow", "C3", 90, "medium", 12.5,
     [ARCHIVE / "h_dtladder_sin2A0.6_Q90_medium_dt12.5"]),
]
NEW = [
    ("steady_Q135", "curve", "steady", 135, "medium", 12.5, "ais2",
     "h_aisee_steady_Q135_medium_dt12.5"),
    ("steady_Q180", "curve", "steady", 180, "medium", 12.5, "ais2",
     "h_aisee_steady_Q180_medium_dt12.5"),
    ("steady_Q360", "curve", "steady", 360, "medium", 12.5, "ais1b",
     "h_aisee_steady_Q360_medium_dt12.5"),
    ("C1_Q46.6", "matched-energy", "C1", 46.6, "medium", 12.5, "ais2",
     "h_aisee_sq1D0.25_Q46.6_medium_dt12.5"),
    ("C2p_Q64.2", "matched-energy", "C2p", 64.2, "medium", 12.5, "ais2",
     "h_aisee_sq0.5D0.5_Q64.2_medium_dt12.5"),
    ("C3_Q82.9", "matched-energy", "C3", 82.9, "medium", 12.5, "ais2",
     "h_aisee_sin2A0.6_Q82.9_medium_dt12.5"),
    ("steady_Q90_coarse_h270", "grid", "steady", 90, "coarse", 12.5, "ais2",
     "h_aisee_steady_Q90_coarse_dt12.5_h270"),
    ("steady_Q90_fine_h270", "grid", "steady", 90, "fine", 12.5, "ais2",
     "h_aisee_steady_Q90_fine_dt12.5_h270"),
    ("C1_Q90_dt6.25_h300", "dt", "C1", 90, "medium", 6.25, "ais2",
     "h_aisee_sq1D0.25_Q90_medium_dt6.25_h300"),
    ("steady_Q360_dt6.25", "dt", "steady", 360, "medium", 6.25, "ais3",
     "h_aisee_steady_Q360_medium_dt6.25"),
    ("steady_Q360_cand2", "scheme", "steady", 360, "medium", 12.5, "ais4",
     "h_aisee_steady_Q360_medium_dt12.5_cand2"),
]
ARM_LABEL = {"steady": "steady", "C1": "square 1 Hz, D 0.25",
             "C2p": "square 0.5 Hz, D 0.5", "C3": "sine 2 Hz, A 0.6"}


def cases():
    out = []
    for key, role, arm, q, grid, dt, cands in REUSED:
        out.append(dict(key=key, role=role, arm=arm, q=q, grid=grid, dt=dt,
                        path=next((c for c in cands if c.exists()), None),
                        source="reused"))
    for key, role, arm, q, grid, dt, part, name in NEW:
        cands = [ARCHIVE / name]
        out.append(dict(key=key, role=role, arm=arm, q=q, grid=grid, dt=dt,
                        path=next((c for c in cands if (c / "case.env").exists()), None),
                        source=part))
    return out


def _done(case: Path) -> bool:
    """A finished run: its solver log reached the case's end time."""
    try:
        meta = HIO.case_meta(case)
        t, _ = HIO.column(case, "fluid", "Tmax_fluid")
        return t.size > 0 and t[-1] >= float(meta["t_end"]) - 1e-6
    except (OSError, KeyError, ValueError):
        return False


def score(case: Path, refresh=False) -> dict:
    """The scorer's gates and info, cached by the log's size and mtime."""
    CACHE.mkdir(parents=True, exist_ok=True)
    logs = sorted(case.glob("log.chtMultiRegionFoam*"))
    stamp = [(p.name, p.stat().st_size, int(p.stat().st_mtime)) for p in logs]
    cf = CACHE / f"{case.name}.json"
    if cf.exists() and not refresh:
        c = json.loads(cf.read_text())
        if c.get("stamp") == [list(s) for s in stamp]:
            return c
    res, info = HM.evaluate(case)
    gates = [{"name": n, "ok": bool(r.ok), "detail": str(r.detail)} for n, r in res]
    out = {"stamp": stamp, "gates": gates, "info": json.loads(
        json.dumps(info, default=lambda o: None if isinstance(o, float) and not np.isfinite(o) else str(o)))}
    cf.write_text(json.dumps(out, indent=1, default=str))
    return out


def series(case: Path):
    """t, T_max(t) over the six cells, dT_cell(t), module mean T(t)."""
    t, smax, tav = None, None, []
    for i in range(N_CELLS):
        ts, v = HIO.column(case, f"cell{i}", f"Tmax_cell{i}")
        t = ts if t is None else t
        smax = v if smax is None else np.maximum(smax, v)
        tav.append(HIO.column(case, f"cell{i}", f"TvolAverage_cell{i}")[1])
    vols = [float(HIO.segments(case, f"cell{i}", f"TvolAverage_cell{i}")[0][1]
                  .meta["Volume"]) for i in range(N_CELLS)]
    tav = np.vstack(tav)
    mod = (np.asarray(vols)[:, None] * tav).sum(0) / sum(vols)
    return t, smax, tav.max(0) - tav.min(0), mod


def quantities(c: dict, sc: dict) -> dict:
    case = c["path"]
    meta = HIO.case_meta(case)
    t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
    t_in = float(meta["T_inlet"])
    t, smax, dtc, mod = series(case)
    info = sc["info"]
    hp = info.get("hydraulic_pre") or {}
    a, b = hp.get("window", [np.nan, np.nan])
    # square waves: less the registered BDF2 valve-edge loss (E_edge_model,
    # D-P.7 (5)), a time-discretisation dissipation at the 4-step edges, not
    # a pump cost (D-AIS.1 (7); checked by the 6.25 ms C1 run)
    e_edge = float(hp.get("E_edge_model", np.nan))
    e_edge = e_edge if np.isfinite(e_edge) else 0.0
    e_raw = float(hp.get("E_port", np.nan)) - float(hp.get("dKE", 0.0))
    e_forced = e_raw - e_edge
    span = float(b) - float(a) if b == b else np.nan
    p_pump = e_forced / span
    heated = t > t_pre
    failed = [g["name"] for g in sc["gates"] if not g["ok"]]
    return dict(
        key=c["key"], role=c["role"], arm=c["arm"], waveform=ARM_LABEL[c["arm"]],
        Qbar_mL_min=c["q"], grid=c["grid"], dt_ms=c["dt"], source=c["source"],
        t_heat_s=round(t_end - t_pre, 6),
        Tmax_peak_C=float(smax[heated].max()) - 273.15,
        Tmax_rise_K=float(smax[heated].max()) - t_in,
        dTcell_end_K=float(dtc[-1]),
        dTcell_mean_K=float(np.trapz(dtc[heated], t[heated]) /
                            (t[heated][-1] - t[heated][0])),
        Tmodule_end_C=float(mod[-1]) - 273.15,
        P_pump_uW=p_pump * 1e6,
        E_pump_mJ=p_pump * HEAT_S * 1e3,
        E_pump_raw_mJ=e_raw / span * HEAT_S * 1e3,
        E_edge_removed_mJ=e_edge / span * HEAT_S * 1e3,
        pump_window_s=f"[{float(a):g}, {float(b):g}]",
        ledger_rel_pct=100 * float((info.get("ledger") or {}).get("residual_rel", np.nan))
        if isinstance(info.get("ledger"), dict) else np.nan,
        gates_passed=f"{len(sc['gates']) - len(failed)}/{len(sc['gates'])}",
        gates_failed=";".join(failed),
    )


def at_time(case: Path, t_s: float):
    """T_max(t), dT_cell and module mean at the sample nearest t_s."""
    t, smax, dtc, mod = series(case)
    i = int(np.argmin(np.abs(t - t_s)))
    return dict(t=float(t[i]), Tmax_C=float(smax[i]) - 273.15,
                dTcell_K=float(dtc[i]), Tmodule_C=float(mod[i]) - 273.15)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-figs", action="store_true")
    ap.add_argument("--refresh", action="store_true", help="re-score cached cases")
    ap.add_argument("--allow-fail", nargs="*", default=[],
                    help="case keys whose failed gates are flagged, not fatal")
    a = ap.parse_args(argv)
    rows, missing = [], []
    for c in cases():
        if c["path"] is None or not _done(c["path"]):
            missing.append(c["key"])
            continue
        sc = score(c["path"], a.refresh)
        r = quantities(c, sc)
        if r["gates_failed"] and c["key"] not in a.allow_fail:
            sys.exit(f"aisee_results: {c['key']} fails gates {r['gates_failed']} "
                     f"(rule 4); fix it or name it in --allow-fail to flag it")
        rows.append(r)
        print(f"{r['key']:<24} Tmax {r['Tmax_peak_C']:.3f} C  dTcell {r['dTcell_end_K']:.4f} K  "
              f"E_pump {r['E_pump_mJ']:.4g} mJ  gates {r['gates_passed']}"
              + (f"  FAILED {r['gates_failed']}" if r["gates_failed"] else ""))
    if missing:
        print("not yet available: " + ", ".join(missing))
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {OUT_CSV} ({len(rows)} cases)")
    if not a.no_figs:
        import aisee_figures
        aisee_figures.make_all(rows, {c["key"]: c["path"] for c in cases()
                                      if c["path"] is not None and c["key"] not in missing},
                               FIGDIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
