#!/usr/bin/env python3
"""Hard gates for the resubmission campaign (WP 2.4; rule 4,
PROJECT.md §9, plan D-F / D-T / NEW-25 / NEW-26, register R-01).

Every gate returns ``GateResult(ok, detail)`` — a 2-tuple, so
``ok, detail = gate(...)`` works. There is NO waiver argument anywhere: D-F
fixes zero waivers and a failure means a rerun. Non-finite or missing input
FAILS a gate (it never passes silently).

Gates
- ledger_gate / ledger_series_gate: whole-run energy ledger closure
  |E_inj - sum(terms)| <= 2 % of the energy the solver injected (D-F, R-01(d));
  the series form optionally checks the per-interval state continuity of
  R-01(b).
- source_power_gate: power actually injected (meshed volume x q''') equals
  the intended power (the v1 0.985 x Eq. 5 source-volume shortfall).
- lattice_gate: fixed dt from 25*2^-k ms; period, 1/4 and 1/8 period, ON
  time, valve edge and sampling interval are integer step counts; a square
  wave's valve edge is exactly tau_v = min(50 ms, ON/5) and >= 2 steps;
  >= N_min steps per period (N_min in {40, 80, 160}); phase samplers
  n_phi >= 8, or 16 for mechanism cases (D-T, NEW-2).
- fixed_dt_gate: >= 99.9 % of solver steps at the fixed dt and a step count
  consistent with the elapsed time (D-F; NEW-26 true fixed-dt mode).
- max_principle_gate: T_in - 0.05 K <= T_fluid <= max T_solid (D-F, D-J).
- continuity_gate: monitor segments — no backward time, no duplicate time
  stamps or segments, no gaps > 2 output intervals, no cadence change,
  segment starts are genuine (re)start times of the case's own log and the
  series ends where the log ends (both REQUIRED), and an overlapping restart
  segment agrees with its predecessor and splices without a jump (the v1
  stale-segment trap, R-01(a), NEW-25).
- courant_gate: max Courant <= the guard (D-F, §9); a DIAGNOSTIC run (ladder
  or coupling rung, never claim-eligible) is reported, not gated (D-F.1).
- thin_layer_courant_gate: the Courant rule of the K0 housing family by
  region (PROJECT.md §2.5 D-F.1, per class since D-F.5, three regions since
  §2.7 D-F.9, 2026-10-02): Co <= B (COURANT_THIN_B_CLASS; provisionally
  COURANT_THIN_B = COURANT_THIN_MARGIN x COURANT_THIN_TESTED for every
  class) in the thin layer T (the cap-plane band of housing_physics.
  THIN_LAYER and the mesh-defined headspace set T_head, THIN_HEAD), <= B_out
  (COURANT_OUT_B_CLASS: the D-F.5 (4) rule applied to the outflow column
  R_out, provisionally 1 for every class) in R_out, and <= 1 (COURANT_LIMIT,
  COURANT_BULK_B_CLASS = 1 for every class: D-F.5 (4)'s rule-updated B_bulk
  is withdrawn) in the claim-region bulk (fluid - T - R_out); T is scaled
  by the cap-cell level ratio and R_out and the bulk by r_h^k off the
  Delta-t ladder's level (never below 1), and the log's every-step global
  maximum must stay <= max(B, B_out, 1), each at its region's scaling. The
  class is the waveform: steady flows take the steady ladder's bounds,
  square waves C1's, sine waves C3's; a spec (no class) is judged at the
  smallest bounds (strict).
- mean_flow_gate: Qbar over whole periods within 0.5 % of intent (§9, D-T).
- hydraulic_closure_gate: |E_port - (Phi - W_b + dKE)| <= TOL_level x
  max(|Phi|, |W_b|) (D-P closure as re-specified in WP 2.5c and amended by
  its review, per waveform class since D-P.7 (2026-10-02): judged on the
  developed ISOTHERMAL windows of STEADY flow - isothermal runs and the
  pre-circulation of every heated run; sine windows report it (GDT2 C3
  +25.7 % on medium) and gate its necessary checks below; heated windows
  and square waves are reported, not gated, by the caller,
  housing_metrics / pump_energy.measured; the allowance covers the
  cell-centred Phi's discretisation difference, P00b medium +7.1...+7.4 %
  at 12.5 ms and coarse 18.7 %). One tolerance per grid level
  (HYDRAULIC_TOL_LEVEL, PROJECT.md §2.5 D-P.4: the residual is a spatial
  discretisation difference that falls with refinement, so one value set
  on medium would fail every coarse run); HYDRAULIC_TOL is the production
  (medium) value. PROVISIONAL: D-P.2 / D-P.4 fix them by a calibration
  rule before the first claim-bearing run; a level whose value is None is
  reported there, not gated. A spec (no level) is judged at the smallest
  value (strict).
- closure_phase_sign_gate / closure_stationarity_gate /
  closure_dt_invariance_gate: D-P.7 (3) necessary checks of a sine
  window's closure: every phase bin's per-period residual >= -0.5 % of the
  per-period Phi; each >= 4 s whole-period sub-window's r/Phi within +-2
  points of the window's; the G-dt pair's r/Phi within +-3 points. They
  catch temporal, phase-offset and start-up errors and cannot bound a
  stationary spatial bias in E_port (the face-based-Phi re-examination,
  D-P.7 (4), does). Constants locked and in POLICY (since the
  pre-registration v1.2 re-render, 2026-10-02).
- max_principle_gate (D-F.10): the bound unchanged for every class and
  level; the detail names the worst sample's time when ``times`` is given.

Pre-registered tolerances are module constants (LEDGER_TOL, COURANT_LIMIT,
...). They are LOCKED: a spec file may give data (arrays, times, dt, the
case's timing, write precision) but never a tolerance nor a switch that
exempts or rescales a check (diagnostic, edge_sensitivity, the grid level),
and the CLI has no flag for a tolerance (LOCKED lists them per gate; the
CLI's --diagnostic is a manual report of a DIAGNOSTIC rung, never a
housing_metrics score). An unknown, locked or missing
spec argument is an input error: that gate is reported as INPUT-ERROR, the
other gates still run and are reported, and the exit code is 2.

CLI (exit 0 = all pass, 1 = a gate failed, 2 = input error):
    gates.py spec case_gates.json [--out results.json]
    gates.py courant LOG [--diagnostic]
    gates.py fixed-dt LOG --dt 0.00625 [--t-start T]
    gates.py lattice --dt 0.00625 [--period 1] [--on-time 0.25] [--sample 0.0625]
                     [--valve-edge 0.05] [--n-min 40] [--n-phase-min 8]
                     [--diagnostic] [--edge-sensitivity]
    gates.py continuity --glob 'postProcessing/x/*/volFieldValue.dat'
                     --log log.1 [--log log.2] [--column 1] [--dt 0.00625]
                     [--write-precision 8] [--scan-jumps]
    gates.py ledger --e-injected J --term dH_cells=J --term out=J ...
    gates.py max-principle --csv ts.csv --t-in 298.15 --fluid-min C
                     --fluid-max C --solid-max C
A spec file is {"case": ..., "gates": [{"gate": "<name>", ...kwargs}]};
any array argument may be {"csv": path, "column": name} instead of a list.

WP 2.4 status and hand-offs (review of prep/2.4, 2026-09-25)
- WP 2.4 is NOT complete against plan P-06. Still missing: the per-step
  u_iter extractor (D-F: change of interface flux and T between the last
  two outer correctors, accumulated over the run), the claims classifier,
  the manifest reader, a fieldMinMax parser (for max_principle inputs) and
  per-row gate status (one gate record per manifest row). The tracker must
  not show 2.4 as done.
- WP 2.5 must set, explicitly, in every campaign controlDict:
  * timePrecision >= min_time_precision(t_end, dt) (v2606 raises it
    automatically only at write times, Time.C:1191, so 'Time =' lines and
    time names can be rounded in between and the fixed-dt, cadence and
    splice checks would fail spuriously);
  * writePrecision (recommended >= 8: 6 digits is a 1 mK quantum at
    ~300 K, the size of the effects and of the 5 mK noise floor), passed to
    the continuity gate as write_precision - the splice/jump tolerance floor
    is the print quantum it implies (DEFAULT_WRITE_PRECISION 6 until then).
- The ledger interval factor (ledger_series_gate interval_factor) is
  RETIRED for the housing family (PROJECT.md §2.5 D-F.2, 2026-09-28):
  LEDGER_INTERVAL_FACTOR stays None; state continuity (R-01(b)) is the
  continuity gate's and the discrete whole-run ledger's (housing_metrics);
  a spec cannot set it.
"""

from __future__ import annotations

import argparse
import glob as _glob
import inspect
import json
import math
import re
import sys
from pathlib import Path
from typing import NamedTuple, Sequence

import numpy as np


class GateResult(NamedTuple):
    ok: bool
    detail: str


class InputError(GateResult):
    """A gate that could not be evaluated because its spec was invalid
    (unknown/locked/missing argument, unreadable data). Counts as a failure
    and makes the run exit 2."""
    __slots__ = ()


def _fail(msg: str) -> GateResult:
    return GateResult(False, msg)


# --------------------------------------------------------------------------
# pre-registered constants (LOCKED: not settable from a spec file)
# --------------------------------------------------------------------------

LEDGER_TOL = 0.02                 # D-F: |residual| <= 2 % of injected energy
LEDGER_INTERVAL_FACTOR = None     # R-01(b): retired for the housing family (D-F.2): the
                                  # naive per-interval form carries half a step of
                                  # storage at every heat-on/restart; continuity gate
                                  # + discrete ledger cover state continuity
SOURCE_POWER_RTOL = 1e-4          # meshed source volume x q_vol vs intended power
HYDRAULIC_TOL = 0.10              # D-P closure (WP 2.5c re-specification;
                                  # was 0.05: P00b medium isothermal steady
                                  # closes at +7.1...+7.4 % at 12.5 ms).
                                  # PROVISIONAL: one grid at one dt; no grid
                                  # order is established (the coarse/medium
                                  # ramp pair mixes h and dt). To be fixed
                                  # from a same-dt coarse/medium/fine
                                  # isothermal check and the impedance block
                                  # before pre-registration; sine C3 at
                                  # 12.5 ms is estimated at ~9.9 %.
                                  # D-P.2 (2026-09-28): set from the
                                  # closurecal block by the pre-registered
                                  # rule (smallest of 5/7.5/10/12.5/15 %
                                  # >= 1.25 x the production-grid residual,
                                  # falling with refinement); GDT1: +7.32 %
                                  # at 12.5 ms, +7.36 % at 6.25 ms.
# D-P.4 (2026-09-29): one tolerance per grid level (the closure residual is the
# cell-centred Phi's discretisation difference and falls with refinement);
# PROVISIONAL until the closurecal block (recipe_rules.hydraulic_tol): medium
# as above; coarse from the P00b coarse ramp window (18.7 %, mixed 25/12.5 ms):
# the smallest candidate >= 1.25 x 18.7 % = 23.4 %; fine and xfine take
# medium's until measured (the residual falls with refinement)
HYDRAULIC_TOL_LEVEL = {"coarse": 0.25, "medium": HYDRAULIC_TOL, "fine": 0.10, "xfine": 0.10}
# D-P.7 (PROJECT.md §2.7, 2026-10-02, after GDT2 and before the OSF deposit):
# the closure tolerance above is the STEADY class's; a sine window's closure
# against the cell-centred Phi is reported, and these necessary checks of it
# are gated on isothermal sine windows (closure_phase_sign_gate,
# closure_stationarity_gate) and on the G-dt pair (closure_dt_invariance_gate).
# Locked and in POLICY (quoted in prereg.yaml and parts/campaign.yaml
# gates.policy since the pre-registration v1.2, 2026-10-02).
CLOSURE_PHASE_BINS = 8            # (a) phase bins of whole periods
CLOSURE_SIGN_MIN_REL = -0.005     # (a) every bin's per-period r >= -0.5 % of the per-period Phi
CLOSURE_SUB_WINDOW_S = 4.0        # (b) sub-windows: the smallest whole number of periods >= 4 s
CLOSURE_STATIONARITY_PP = 0.02    # (b) each sub-window's r/Phi within +-2 points of the window's
CLOSURE_DT_INVARIANCE_PP = 0.03   # (c) r/Phi at dt and dt/2 within +-3 points (G-dt pair)
FIXED_DT_RTOL = 1e-3              # step equal to dt within 0.1 %
FIXED_DT_MIN_FRACTION = 0.999     # NEW-26
UNDERSHOOT_K = 0.05               # D-F / D-J: T_fluid >= T_in - 0.05 K
OVERSHOOT_K = 0.0                 # T_fluid <= running max T_solid
COURANT_LIMIT = 1.0               # §9 / NEW-2 Courant guard (the claim-region bulk limit,
                                  # D-F.1 (2) as amended by D-F.9 (3))
COURANT_THIN_TESTED = 4.48        # D-F.1: the largest thin-layer Co_max of a medium ladder
                                  # or coupling rung that passed against its half-dt
                                  # partner (GDT1 steady Q90 12.5 ms rung; no rung has
                                  # failed), PROVISIONAL: re-evaluated by the same rule at
                                  # G-dt and FREEZE, changed only by a dated D-row
COURANT_THIN_MARGIN = 1.25        # D-F.1: guard band for the realisation spread of Co_max
                                  # (4.33 / 4.48 on two steady Q90 realisations at 12.5 ms);
                                  # admissibility of a dt still needs the class's own
                                  # NEW-2 ladder verdict or its NEW-2 bounding rule
COURANT_THIN_B = round(COURANT_THIN_MARGIN * COURANT_THIN_TESTED, 2)   # 5.60
# D-F.5 (2026-09-29): the bounds per class (steady ladder -> steady flows, C1 ->
# square waves, C3 -> sine waves), each by the D-F.1 (3) rule on its own
# admissible medium rungs (recipe_rules.courant_bounds). PROVISIONAL: B is the
# steady-tested 5.60 for every class until the pulsed classes' G-dt rungs
COURANT_THIN_B_CLASS = {"steady": COURANT_THIN_B, "square": COURANT_THIN_B,
                        "sine": COURANT_THIN_B}
# D-F.9 (3) (2026-10-02): the claim-region bulk (fluid - T - R_out) is bounded
# by COURANT_LIMIT for every class; D-F.5 (4)'s rule-updated B_bulk is
# withdrawn (its rule now sets B_out). Kept per class for the class-keyed
# readers (housing_metrics, manuscript_params); never rule-updated
COURANT_BULK_B_CLASS = {"steady": COURANT_LIMIT, "square": COURANT_LIMIT,
                        "sine": COURANT_LIMIT}
# D-F.9 (2): the outflow column R_out = {x >= x_max - x_port} U {z > z_top}
# (the D-G outlet slot column and the riser, downstream of every claim
# quantity) is bounded by its OWN tested value: B_out,c = max(1, 1.25 x the
# R_out Co_max of an admissible medium ladder rung of class c), to 2 decimals,
# capped one step below the smallest R_out Co_max of a failing rung
# (recipe_rules.outflow_b), level-scaled as the bulk. PROVISIONAL 1 for every
# class until its pair passes (GDT2's riser-entrance 3.46 sine / 2.78 square
# give 4.33 / 3.47 only after their G-dt pairs pass); changed only by a dated
# D-row applying the rule
COURANT_OUT_B_CLASS = {"steady": COURANT_LIMIT, "square": COURANT_LIMIT,
                       "sine": COURANT_LIMIT}
# D-F.1 level scaling off the dt-ladder level (the D-C family runs at one dt):
# bulk = r_h^k (r_h 1.35, D-C), also the outflow column's (D-F.9 (2)) and the
# T_head threshold's (D-F.9 (1)); thin layer = the cap-plane first-cell ratio
# of the K0 design (0.050 / 0.0355 / 0.0264 mm at medium / fine / xfine); a
# level coarser than the ladder level keeps the ladder level's limits (ratio
# >= 1)
COURANT_LEVEL_BULK = {"coarse": 1.0 / 1.35, "medium": 1.0, "fine": 1.35, "xfine": 1.8225}
COURANT_LEVEL_THIN = {"coarse": 0.728, "medium": 1.0, "fine": 1.408, "xfine": 1.894}
MEAN_FLOW_RTOL = 0.005            # §9: Qbar within 0.5 % of intent
MEAN_FLOW_ZERO_ATOL = 1e-12       # m^3/s, |Qbar| bound for a zero-flow case
DT_BASE = 0.025                   # s; D-T ladder 25 * 2^-k ms
N_MIN_ALLOWED = (40, 80, 160)     # NEW-2 N_min candidates
N_PHASE_ALLOWED = (8, 16)         # D-T: >= 8; 16 for NEW-20 mechanism cases
VALVE_EDGE_MAX = 0.050            # D-T: tau_v = min(50 ms, ON/5)
VALVE_EDGE_ON_FRACTION = 0.2
VALVE_EDGE_MIN_STEPS = 2
MAX_GAP_FACTOR = 2.0              # continuity
CADENCE_RTOL = 0.05
JUMP_FACTOR = 5.0
JUMP_WINDOW = 5
DEFAULT_WRITE_PRECISION = 6       # OpenFOAM default; WP 2.5 must set >= 8 explicitly
MIN_WRITE_PRECISION = 6


def _arr(x, name: str) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    if a.size == 0:
        raise ValueError(f"{name}: empty")
    return a


def _finite(*arrays) -> bool:
    return all(np.all(np.isfinite(np.asarray(a, dtype=float))) for a in arrays)


# --------------------------------------------------------------------------
# energy
# --------------------------------------------------------------------------

def ledger_gate(e_injected: float, terms: dict, tol: float = LEDGER_TOL) -> GateResult:
    """Whole-run ledger: E_inj = sum(terms) within tol * E_inj.

    ``terms`` are the signed sinks/stores of the same energy in J (e.g.
    dH_cells, dH_fluid, net advected enthalpy out, inlet/outlet conduction,
    wall loss), each positive when it takes up injected energy.
    """
    if not terms:
        return _fail("ledger: no terms")
    vals = {k: float(v) for k, v in terms.items()}
    if not _finite([e_injected], list(vals.values())):
        return _fail(f"ledger: non-finite input (E_inj={e_injected}, terms={vals})")
    if not e_injected > 0:
        return _fail(f"ledger: injected energy must be > 0 (got {e_injected})")
    total = sum(vals.values())
    res = e_injected - total
    rel = res / e_injected
    txt = ", ".join(f"{k}={v:.6g}" for k, v in vals.items())
    detail = (f"ledger residual {res:+.6g} J = {100 * rel:+.3f} % of E_inj "
              f"{e_injected:.6g} J (tol {100 * tol:.2f} %); {txt}")
    return GateResult(abs(rel) <= tol, detail)


def _trapz(y, t):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def ledger_series_gate(t, p_inj, stored, fluxes_out: dict | None = None,
                       tol: float = LEDGER_TOL,
                       interval_factor: float | None = LEDGER_INTERVAL_FACTOR) -> GateResult:
    """Ledger from time series on one time base.

    t [s], p_inj [W] injected power, stored [J] total stored enthalpy (solid +
    fluid, one reference), fluxes_out {name: W} signed outflows (advected
    enthalpy, conduction, walls). Whole run: int P dt = dH + sum int F dt.
    If ``interval_factor`` is set, every sampling interval must also close:
    |dH_i - int_i (P - sum F) dt| <= interval_factor * max(P) * dt_i
    (state continuity, R-01(b); the factor is calibrated before use).
    """
    t = _arr(t, "t")
    p = _arr(p_inj, "p_inj")
    H = _arr(stored, "stored")
    fl = {k: _arr(v, k) for k, v in (fluxes_out or {}).items()}
    n = t.size
    if n < 2 or p.size != n or H.size != n or any(v.size != n for v in fl.values()):
        return _fail("ledger series: arrays must share one time base (>= 2 samples)")
    if not _finite(t, p, H, *fl.values()):
        return _fail("ledger series: non-finite samples")
    if np.any(np.diff(t) <= 0):
        return _fail("ledger series: time not strictly increasing")
    e_inj = _trapz(p, t)
    terms = {"dH_stored": float(H[-1] - H[0])}
    terms.update({f"out_{k}": _trapz(v, t) for k, v in fl.items()})
    ok, detail = ledger_gate(e_inj, terms, tol)
    if interval_factor is not None:
        net = p - sum(fl.values()) if fl else p
        dt = np.diff(t)
        r = np.diff(H) - 0.5 * (net[1:] + net[:-1]) * dt
        lim = interval_factor * float(np.max(np.abs(p))) * dt
        bad = np.nonzero(np.abs(r) > lim)[0]
        if bad.size:
            i = int(bad[0])
            ok = False
            detail += (f"; state continuity broken at {bad.size} interval(s), first "
                       f"[{t[i]:.6g}, {t[i + 1]:.6g}] s: residual {r[i]:+.4g} J "
                       f"> {lim[i]:.4g} J")
        else:
            detail += f"; per-interval residual max {np.max(np.abs(r)):.4g} J within limit"
    else:
        detail += "; per-interval state continuity NOT checked (factor not yet calibrated)"
    return GateResult(ok, detail)


def source_power_gate(p_injected: float, p_intended: float,
                      rtol: float = SOURCE_POWER_RTOL) -> GateResult:
    """Injected power (meshed source volume x q''') equals the intended power."""
    if not _finite([p_injected, p_intended]) or not p_intended > 0:
        return _fail(f"source power: invalid input ({p_injected}, {p_intended})")
    rel = p_injected / p_intended - 1.0
    return GateResult(abs(rel) <= rtol,
                      f"source power {p_injected:.8g} W vs intended {p_intended:.8g} W "
                      f"({rel:+.3e}, tol {rtol:.1e})")


def hydraulic_closure_gate(e_pump: float, phi_diss: float, w_buoy: float,
                           d_ke: float, level: str | None = None,
                           tol: float | None = None) -> GateResult:
    """D-P: |E_port - (Phi - W_b + dKE)| <= tol * max(|Phi|, |W_b|).

    e_pump is the measured port work E_port over the closure window. The
    gate is a check of the port-work instrument at the scale of Phi and is
    applied only to developed isothermal windows of STEADY flow (D-P.7;
    the caller's scope: housing_metrics / pump_energy.measured; sine
    windows report the closure and gate its necessary checks). A square
    wave's residual is dominated by the backward scheme's dissipation of
    the kinetic energy its valve edges exchange (P00b valve case: 2.1 Phi),
    which a reference widened by that energy would let an E_port error of
    several Phi through; square waves are reported, not gated.

    tol None: the grid level's HYDRAULIC_TOL_LEVEL (D-P.4; housing_metrics
    passes the case's grid); no level: the smallest level value (strict). A
    level whose tolerance is None (D-P.4: calibrated as reported-only) is
    reported and passes."""
    if not _finite([e_pump, phi_diss, w_buoy, d_ke]):
        return _fail("hydraulic closure: non-finite input")
    ref = max(abs(phi_diss), abs(w_buoy))
    if ref == 0.0:
        return _fail("hydraulic closure: Phi and W_b both zero")
    res = e_pump - (phi_diss - w_buoy + d_ke)
    where = "strict: smallest level value" if level is None else f"level {level}"
    if tol is None:
        if level is None:
            tol = min(v for v in HYDRAULIC_TOL_LEVEL.values() if v is not None)
        elif level not in HYDRAULIC_TOL_LEVEL:
            return _fail(f"hydraulic closure: unknown grid level {level!r}")
        else:
            tol = HYDRAULIC_TOL_LEVEL[level]
    det = (f"hydraulic closure residual {res:+.4g} J = "
           f"{100 * res / ref:+.2f} % of max(|Phi|,|W_b|) {ref:.4g} J")
    if tol is None:
        return GateResult(True, det + f" ({where}: reported, not gated, D-P.4)")
    return GateResult(abs(res) <= tol * ref, det + f" (tol {100 * tol:.1f} %, {where}, provisional)")


def closure_phase_sign_gate(r_bins_per_period, phi_per_period: float,
                            min_rel: float = CLOSURE_SIGN_MIN_REL) -> GateResult:
    """D-P.7 (3a), sine windows: the per-period closure residual summed over
    each phase bin of whole periods is >= min_rel x the per-period Phi in
    every bin (dissipation is one-signed; a phase or time-level error of
    E_port need not be)."""
    r = _arr(r_bins_per_period, "r_bins_per_period")
    if not _finite(r, [phi_per_period]) or not phi_per_period > 0:
        return _fail(f"closure phase sign: invalid input (Phi per period {phi_per_period})")
    rel = r / phi_per_period
    j = int(np.argmin(rel))
    det = (f"closure residual by phase bin ({r.size} bins): min {100 * rel[j]:+.3f} % of "
           f"the per-period Phi (bin {j}), max {100 * rel.max():+.3f} % "
           f"(limit {100 * min_rel:+.2f} %)")
    return GateResult(bool(rel.min() >= min_rel), det)


def closure_stationarity_gate(sub_r_over_phi, window_r_over_phi: float,
                              tol: float = CLOSURE_STATIONARITY_PP) -> GateResult:
    """D-P.7 (3b), sine windows: r/Phi of every whole-period sub-window
    (>= CLOSURE_SUB_WINDOW_S) lies within +-tol (absolute, fraction of Phi)
    of the window's r/Phi; at least two sub-windows."""
    s = _arr(sub_r_over_phi, "sub_r_over_phi")
    if s.size < 2:
        return _fail(f"closure stationarity: {s.size} sub-window(s); two are needed")
    if not _finite(s, [window_r_over_phi]):
        return _fail("closure stationarity: non-finite input")
    dev = s - window_r_over_phi
    j = int(np.argmax(np.abs(dev)))
    det = (f"closure r/Phi of {s.size} sub-windows {100 * s.min():+.2f} ... "
           f"{100 * s.max():+.2f} % against the window's {100 * window_r_over_phi:+.2f} %: "
           f"largest departure {100 * dev[j]:+.2f} points (tol +-{100 * tol:.1f})")
    return GateResult(bool(np.abs(dev).max() <= tol), det)


def closure_dt_invariance_gate(r_over_phi_dt: float, r_over_phi_half: float,
                               tol: float = CLOSURE_DT_INVARIANCE_PP) -> GateResult:
    """D-P.7 (3c), the sine class at G-dt: the isothermal closure r/Phi of
    the dt rung and its dt/2 partner agree within +-tol (absolute, fraction
    of Phi). Outside it the excess is temporal (backward's damping, a
    flux-pressure offset tied to the step)."""
    if not _finite([r_over_phi_dt, r_over_phi_half]):
        return _fail("closure dt invariance: non-finite input")
    d = r_over_phi_dt - r_over_phi_half
    return GateResult(bool(abs(d) <= tol),
                      f"closure r/Phi {100 * r_over_phi_dt:+.2f} % at dt against "
                      f"{100 * r_over_phi_half:+.2f} % at dt/2: {100 * d:+.2f} points "
                      f"(tol +-{100 * tol:.1f})")


# --------------------------------------------------------------------------
# time-step lattice and fixed dt
# --------------------------------------------------------------------------

def _steps(duration: float, dt: float, rtol: float = 1e-9):
    x = duration / dt
    n = round(x)
    return n, abs(x - n) <= rtol * max(1.0, abs(x))


def phase_samplers(steps_per_period: int) -> int:
    """D-T: n_phi = largest power of two <= 16 dividing the steps per period."""
    for n in (16, 8, 4, 2, 1):
        if steps_per_period % n == 0:
            return n
    return 1


def valve_edge_required(on_time: float) -> float:
    """D-T: tau_v = min(50 ms, ON/5), a physical time."""
    return min(VALVE_EDGE_MAX, VALVE_EDGE_ON_FRACTION * on_time)


def lattice_gate(dt: float, period: float | None = None,
                 on_time: float | None = None, sample: float | None = None,
                 valve_edge: float | None = None, n_min: int = 40,
                 n_phase_min: int = 8, dt_ladder: bool = True,
                 diagnostic: bool = False,
                 edge_sensitivity: bool = False) -> GateResult:
    """D-T lattice: every timing quantity an integer number of fixed steps.

    period=None means a steady case (period, ON time and phase checks skip).
    on_time given means a square wave: its valve_edge must be given and equal
    tau_v = min(50 ms, ON/5) exactly, with >= 2 steps. ``edge_sensitivity``
    (the ideal/25/200 ms edge runs) lifts only the tau_v value (an ideal
    edge is valve_edge=None); ``diagnostic`` (ladder/coupling runs, NEW-2)
    lifts only N_min. Both are written into the detail: such runs are not
    claim-bearing.
    """
    if not (_finite([dt]) and dt > 0):
        return _fail(f"lattice: invalid dt {dt}")
    fails, info = [], []
    if n_min not in N_MIN_ALLOWED:
        fails.append(f"n_min {n_min} is not a pre-registered N_min {N_MIN_ALLOWED}")
    if n_phase_min not in N_PHASE_ALLOWED:
        fails.append(f"n_phase_min {n_phase_min} not in {N_PHASE_ALLOWED}")
    if dt_ladder:
        k = math.log2(DT_BASE / dt)
        if abs(k - round(k)) > 1e-9 or round(k) < 0:
            fails.append(f"dt {dt:g} s is not 25*2^-k ms")
    checks = {"sample": sample, "valve_edge": valve_edge}
    if period is not None:
        checks.update(period=period, quarter_period=period / 4,
                      eighth_period=period / 8, on_time=on_time)
    elif on_time is not None:
        fails.append("on_time given for a steady case")
    for name, dur in checks.items():
        if dur is None:
            continue
        if not (_finite([dur]) and dur > 0):
            fails.append(f"{name} invalid ({dur})")
            continue
        n, ok = _steps(dur, dt)
        if not ok:
            fails.append(f"{name} {dur:g} s = {dur / dt:.6g} steps (not an integer)")
        else:
            info.append(f"{name}={n}")
    if valve_edge is not None and _finite([valve_edge]):
        n, ok = _steps(valve_edge, dt)
        if ok and n < VALVE_EDGE_MIN_STEPS:
            fails.append(f"valve edge {n} step(s) < {VALVE_EDGE_MIN_STEPS}")
    if on_time is not None and period is not None and _finite([on_time]) and on_time > 0:
        tau = valve_edge_required(on_time)
        if edge_sensitivity:
            info.append(f"EDGE-SENSITIVITY run (tau_v {valve_edge} s; D-T value "
                        f"{tau:g} s not enforced; not claim-bearing)")
        elif valve_edge is None:
            fails.append(f"square wave without valve_edge: D-T tau_v = "
                         f"min(50 ms, ON/5) = {tau:g} s required")
        elif _finite([valve_edge]) and abs(valve_edge - tau) > 1e-9 * max(1.0, tau):
            fails.append(f"valve edge {valve_edge:g} s != D-T tau_v = "
                         f"min(50 ms, ON/5) = {tau:g} s")
    elif valve_edge is not None:
        fails.append("valve_edge given for a case without an ON time")
    if period is not None and _finite([period]) and period > 0:
        n, ok = _steps(period, dt)
        if ok:
            if n < n_min:
                if diagnostic:
                    info.append(f"DIAGNOSTIC run: {n} steps/period < N_min {n_min} "
                                f"allowed (not claim-bearing)")
                else:
                    fails.append(f"{n} steps/period < N_min {n_min}")
            nphi = phase_samplers(n)
            info.append(f"n_phi={nphi}")
            if nphi < n_phase_min:
                fails.append(f"phase samplers n_phi {nphi} < {n_phase_min}")
    elif diagnostic:
        info.append("DIAGNOSTIC run")
    if fails:
        return _fail("lattice: " + "; ".join(fails))
    return GateResult(True, f"lattice dt={dt:g} s: " + ", ".join(info or ["steady"]))


def _decimals(x: float, max_dec: int = 15) -> int:
    """Decimal places needed to write x exactly (x a finite decimal)."""
    for k in range(max_dec + 1):
        v = x * 10 ** k
        if abs(v - round(v)) <= 1e-9 * max(1.0, abs(v)):
            return k
    raise ValueError(f"{x!r} is not a short decimal")


def min_time_precision(t_end: float, dt: float) -> int:
    """Significant digits timePrecision must have so that every time
    t = n dt <= t_end prints exactly (hand-off to WP 2.5): the integer
    digits of t_end plus the decimals of dt (8 for 900 s at 6.25 ms,
    10 at 1.5625 ms)."""
    if not (t_end > 0 and dt > 0):
        raise ValueError("t_end and dt must be > 0")
    return int(math.floor(math.log10(t_end))) + 1 + _decimals(dt)


def fixed_dt_gate(step_times, dt: float, t_start: float | None = None,
                  rtol: float = FIXED_DT_RTOL,
                  min_fraction: float = FIXED_DT_MIN_FRACTION) -> GateResult:
    """Solver step times (the log's 'Time =' values): >= min_fraction of steps
    equal dt within rtol, and the step count matches (t_end - t_0)/dt.

    v2606 prints no 'deltaT =' line with adjustTimeStep off
    (setMultiRegionDeltaT.H), so steps come from 'Time =' lines; the case
    must write times with enough digits: timePrecision >=
    min_time_precision(t_end, dt) (too few digits make steps uneven and the
    gate fails - it cannot pass spuriously).
    """
    t = _arr(step_times, "step_times")
    if t_start is not None:
        t = np.concatenate([[float(t_start)], t])
    if t.size < 2 or not _finite(t) or not (_finite([dt]) and dt > 0):
        return _fail("fixed dt: need >= 2 finite step times and dt > 0")
    d = np.diff(t)
    if np.any(d <= 0):
        return _fail(f"fixed dt: time not strictly increasing ({int(np.sum(d <= 0))} steps)")
    on = np.abs(d / dt - 1.0) <= rtol
    frac = float(on.mean())
    n_expect = (t[-1] - t[0]) / dt
    fails = []
    if frac < min_fraction:
        fails.append(f"only {100 * frac:.3f} % of steps at dt (need {100 * min_fraction:.1f} %)")
    if abs(n_expect - d.size) > 0.5:
        fails.append(f"{d.size} steps over {t[-1] - t[0]:.6g} s but {n_expect:.3f} expected")
    detail = (f"fixed dt {dt:g} s: {100 * frac:.3f} % of {d.size} steps on dt, "
              f"step range [{d.min():.6g}, {d.max():.6g}] s")
    return GateResult(not fails, detail + ("; " + "; ".join(fails) if fails else ""))


# --------------------------------------------------------------------------
# physics
# --------------------------------------------------------------------------

def max_principle_gate(t_fluid_min, t_fluid_max, t_in: float, t_solid_max,
                       undershoot: float = UNDERSHOOT_K, overshoot: float = OVERSHOOT_K,
                       running_max: bool = True, times=None) -> GateResult:
    """T_in - undershoot <= T_fluid_min(t) and T_fluid_max(t) <= bound(t).

    bound(t) = max(T_in, running max of T_solid_max up to t) + overshoot:
    heat reaches the fluid only through the solids and the inlet, so the
    fluid can never exceed the hottest solid state seen so far
    (``running_max=False`` compares sample by sample). ``times`` (the
    samples' times, optional data) makes the detail name the worst
    sample's TIME and the span of the failing samples instead of an index
    into the series (D-F.10 (4): GDT2's first reading took the heated-sample
    index 825 for t = 52 s; the sample was at t = 81.625 s)."""
    fmin, fmax, smax = (_arr(a, n) for a, n in ((t_fluid_min, "T_fluid_min"),
                                                  (t_fluid_max, "T_fluid_max"),
                                                  (t_solid_max, "T_solid_max")))
    if not (fmin.size == fmax.size == smax.size):
        return _fail("max principle: series lengths differ")
    if not _finite(fmin, fmax, smax, [t_in]):
        return _fail("max principle: non-finite input")
    tt = None
    if times is not None:
        tt = _arr(times, "times")
        if tt.size != fmin.size or not _finite(tt):
            return _fail("max principle: times do not match the series")
    lo = t_in - undershoot
    sbound = np.maximum.accumulate(smax) if running_max else smax
    hi = np.maximum(sbound, t_in) + overshoot
    under = fmin < lo
    over = fmax > hi

    def at(i, mask):
        if tt is None:
            return f"sample {i}"
        k = np.flatnonzero(mask)
        return (f"t = {tt[i]:.10g} s (failing samples t = {tt[k].min():.10g} ... "
                f"{tt[k].max():.10g} s)")
    detail = (f"fluid T in [{fmin.min():.6g}, {fmax.max():.6g}] K; "
              f"lower bound {lo:.6g} K (T_in - {undershoot:g}), "
              f"upper bound max T_solid")
    fails = []
    if under.any():
        i = int(np.argmin(fmin))
        fails.append(f"undershoot {fmin[i] - t_in:+.4g} K at {at(i, under)} "
                     f"({int(under.sum())} samples)")
    if over.any():
        i = int(np.argmax(fmax - hi))
        fails.append(f"fluid max exceeds solid bound by {fmax[i] - hi[i]:+.4g} K at "
                     f"{at(i, over)} ({int(over.sum())} samples)")
    return GateResult(not fails, detail + ("; " + "; ".join(fails) if fails else ""))


def courant_gate(co_max_series, limit: float = COURANT_LIMIT,
                 diagnostic: bool = False) -> GateResult:
    """max Courant over the run <= limit (the pre-registered guard).

    ``diagnostic`` (a DIAGNOSTIC=1 ladder or coupling rung: Delta t is its
    independent variable, and it is never claim-eligible) reports the value
    and passes (D-F.1 (5)); the detail says so."""
    co = _arr(co_max_series, "Courant")
    if not _finite(co) or not (_finite([limit]) and limit > 0):
        return _fail("Courant: non-finite input")
    m = float(co.max())
    detail = f"Courant max {m:.4g} over {co.size} samples (limit {limit:g})"
    if diagnostic:
        return GateResult(True, detail + "; DIAGNOSTIC run: reported, not gated "
                                          "(D-F.1), not claim-eligible")
    return GateResult(m <= limit, detail)


def courant_level_ratios(level: str, ladder_level: str = "medium"):
    """(bulk, thin) limit multipliers of a grid level relative to the level
    of the Delta-t ladder (D-F.1 (2)); never below 1."""
    for lv in (level, ladder_level):
        if lv not in COURANT_LEVEL_BULK:
            raise ValueError(f"unknown grid level {lv!r} (one of "
                             f"{sorted(COURANT_LEVEL_BULK)})")
    rb = COURANT_LEVEL_BULK[level] / COURANT_LEVEL_BULK[ladder_level]
    rt = COURANT_LEVEL_THIN[level] / COURANT_LEVEL_THIN[ladder_level]
    return max(1.0, rb), max(1.0, rt)


def courant_class_bounds(waveform: str | None = None):
    """D-F.5 / D-F.9: (claim-region bulk bound, thin-layer B) of a waveform
    class (steady, square, sine); None: the strict pair (the smallest bulk
    bound, and the smallest B, or None when a class has no B). The bulk
    bound is COURANT_LIMIT for every class since D-F.9 (3); the outflow
    column's bound is courant_out_bound."""
    if waveform is None:
        bs = list(COURANT_THIN_B_CLASS.values())
        return (min(COURANT_BULK_B_CLASS.values()),
                None if any(b is None for b in bs) else min(bs))
    if waveform not in COURANT_THIN_B_CLASS:
        raise ValueError(f"unknown waveform class {waveform!r} (one of "
                         f"{sorted(COURANT_THIN_B_CLASS)})")
    return COURANT_BULK_B_CLASS[waveform], COURANT_THIN_B_CLASS[waveform]


def courant_out_bound(waveform: str | None = None) -> float:
    """D-F.9 (2): B_out of a waveform class; None: the strict (smallest)."""
    if waveform is None:
        return min(COURANT_OUT_B_CLASS.values())
    if waveform not in COURANT_OUT_B_CLASS:
        raise ValueError(f"unknown waveform class {waveform!r} (one of "
                         f"{sorted(COURANT_OUT_B_CLASS)})")
    return COURANT_OUT_B_CLASS[waveform]


_STRICT = object()


def thin_layer_courant_gate(co_bulk_series, co_thin_series, co_out_series,
                            co_log_series=None, level: str = "medium",
                            ladder_level: str = "medium", thin_limit=_STRICT,
                            out_limit=_STRICT,
                            diagnostic: bool = False) -> GateResult:
    """The Courant rule of the K0 housing family by region (D-F.1, per class
    D-F.5, regions D-F.9).

    co_bulk_series / co_thin_series / co_out_series: the maximum Courant
    number over the claim-region bulk (fluid - T - R_out), the thin layer T
    (cap-plane band + T_head) and the outflow column R_out per monitor
    sample (CoMinMax_fluid's CoBulk / CoThin / CoOut since D-F.9);
    co_log_series: the solver log's every-step global maximum (optional).
    Passes iff max(bulk) <= COURANT_LIMIT * r_bulk, max(thin) <= B * r_thin,
    max(out) <= B_out * r_bulk and max(log) <= max(B * r_thin, B_out * r_bulk,
    COURANT_LIMIT * r_bulk), r = courant_level_ratios(level, ladder_level).
    B and B_out default to the strict class values (courant_class_bounds(
    None), courant_out_bound(None)); housing_metrics passes the case's class.
    B None means no bound is fixed: the thin layer then takes the bulk limit
    (the strict gate). The claim-region bulk limit is not an argument: it is
    1 for every class (D-F.9 (3)). ``diagnostic`` reports and passes, as
    courant_gate."""
    bulk_limit = COURANT_LIMIT
    thin_limit = courant_class_bounds(None)[1] if thin_limit is _STRICT \
        else thin_limit
    out_limit = courant_out_bound(None) if out_limit is _STRICT else out_limit
    cb = _arr(co_bulk_series, "Courant bulk")
    ct = _arr(co_thin_series, "Courant thin layer")
    co = _arr(co_out_series, "Courant outflow column")
    cl = None if co_log_series is None else _arr(co_log_series, "Courant log")
    if not _finite(cb, ct, co) or (cl is not None and not _finite(cl)):
        return _fail("thin-layer Courant: non-finite input")
    b = bulk_limit if thin_limit is None else thin_limit
    if not (_finite([b, out_limit]) and b > 0 and out_limit >= bulk_limit):
        return _fail("thin-layer Courant: invalid limits")
    try:
        rb, rt = courant_level_ratios(level, ladder_level)
    except ValueError as e:
        return _fail(f"thin-layer Courant: {e}")
    lim_b, lim_t, lim_o = bulk_limit * rb, b * rt, out_limit * rb
    lim_l = max(lim_b, lim_t, lim_o)
    mb, mt, mo = float(cb.max()), float(ct.max()), float(co.max())
    ok = mb <= lim_b and mt <= lim_t and mo <= lim_o
    detail = (f"Courant claim-region bulk max {mb:.4g} (limit {lim_b:.4g}), "
              f"thin layer max {mt:.4g} (limit {lim_t:.4g}"
              f"{'' if thin_limit is not None else ', B not fixed: strict'}), "
              f"outflow column max {mo:.4g} (limit {lim_o:.4g}) over {cb.size} "
              f"samples; level {level} vs ladder {ladder_level}")
    if cl is not None:
        ml = float(cl.max())
        ok = ok and ml <= lim_l
        detail += f"; log every-step max {ml:.4g} (limit {lim_l:.4g})"
    if diagnostic:
        return GateResult(True, detail + "; DIAGNOSTIC run: reported, not gated "
                                          "(D-F.1), not claim-eligible")
    return GateResult(ok, detail)


def mean_flow_gate(t, q, q_target: float, period: float | None = None,
                   t_start: float | None = None, rtol: float = MEAN_FLOW_RTOL,
                   atol: float | None = None) -> GateResult:
    """Qbar within rtol of intent. Pulsed: averaged over the largest whole
    number of periods ending at the last sample (after t_start); steady:
    over [t_start, t_end]. atol defaults to 0, or to MEAN_FLOW_ZERO_ATOL
    (m^3/s) for a zero-flow target."""
    t = _arr(t, "t")
    q = _arr(q, "Q")
    if t.size != q.size or t.size < 2 or not _finite(t, q):
        return _fail("mean flow: need matching finite t and Q")
    if np.any(np.diff(t) <= 0):
        return _fail("mean flow: time not strictly increasing")
    if atol is None:
        atol = MEAN_FLOW_ZERO_ATOL if q_target == 0 else 0.0
    if q_target == 0 and atol <= 0:
        return _fail("mean flow: target 0 needs an absolute tolerance")
    t0 = t[0] if t_start is None else max(float(t_start), t[0])
    t1 = t[-1]
    if period is not None:
        nper = math.floor((t1 - t0) / period + 1e-9)
        if nper < 1:
            return _fail(f"mean flow: window shorter than one period ({t1 - t0:g} s)")
        t0 = t1 - nper * period
    grid = np.concatenate([[t0], t[(t > t0) & (t < t1)], [t1]])
    qq = np.interp(grid, t, q)
    qbar = _trapz(qq, grid) / (t1 - t0)
    err = abs(qbar - q_target)
    ok = err <= rtol * abs(q_target) + atol
    rel = f"{100 * (qbar / q_target - 1):+.3f} %" if q_target else f"{qbar:+.4g} abs"
    win = f"{nper} periods" if period is not None else f"{t1 - t0:g} s"
    return GateResult(ok, f"Qbar {qbar:.6g} vs intent {q_target:.6g} ({rel}; "
                          f"tol {100 * rtol:.2f} %) over {win}")


# --------------------------------------------------------------------------
# time-series continuity / segment provenance
# --------------------------------------------------------------------------

def print_quantum(values, write_precision: int) -> float:
    """Largest print quantum of values written with ``write_precision``
    significant digits (OpenFOAM general format): 10^(floor(log10 max|v|)
    + 1 - P); 1e-3 K for 6 digits at ~300 K, 1e-5 K for 8."""
    v = np.abs(np.asarray(values, dtype=float))
    vmax = float(v.max()) if v.size else 0.0
    if vmax == 0.0 or not math.isfinite(vmax):
        return 0.0
    return 10.0 ** (math.floor(math.log10(vmax)) + 1 - int(write_precision))


def continuity_gate(segments: Sequence, log_start_times: Sequence[float] | None = None,
                    log_end_time: float | None = None,
                    max_gap_factor: float = MAX_GAP_FACTOR,
                    cadence_rtol: float | None = CADENCE_RTOL,
                    expected_interval: float | None = None,
                    time_tol: float | None = None, jump_factor: float = JUMP_FACTOR,
                    overlap_atol: float | None = None, window: int = JUMP_WINDOW,
                    scan_jumps: bool = False, jump_atol: float | None = None,
                    write_precision: int = DEFAULT_WRITE_PRECISION,
                    require_provenance: bool = True) -> GateResult:
    """Monitor-segment continuity and provenance.

    segments: [(label, times[, values]), ...] — one entry per postProcessing
    <startTime>/ directory; label is that start time. Checks:
      * no empty segment, no duplicate labels, no segment contained in another
        (duplicate segments);
      * times strictly increasing inside every segment;
      * stitched series (segment k truncated at the start of k+1, NEW-25)
        strictly increasing, no gap > max_gap_factor x the output interval,
        and (cadence_rtol) every interval within rtol of the nominal one;
      * provenance (REQUIRED unless require_provenance=False, which is for
        unit use only and is written into the detail): every label equals a
        start/restart time of the case's own log (log_start_times), and the
        series ends at log_end_time;
      * with values, at every splice k -> k+1, with the tolerance scale taken
        from the EARLIER segment only (its last ``window`` sample-to-sample
        changes before the splice and on its overlapped part): the first
        ``window`` overlapping samples must agree (overlap_atol, default
        jump_factor x scale; later overlap samples are not compared because
        a restart from written fields diverges slowly from the original),
        the first new sample must continue the linear extrapolation of k,
        and the first ``window`` changes of k+1 must stay within
        jump_factor x scale — a step right after an unchanged first sample
        is caught;
      * scan_jumps: also flag isolated jumps anywhere in the stitched values
        (|dv_i| > jump_factor x the largest neighbouring change within
        ``window`` samples) — for series already stitched upstream.
    ``jump_atol`` (default: the print quantum implied by ``write_precision``
    significant digits, see print_quantum) is an absolute floor for every
    value tolerance, so a one-digit rounding difference on a flat monitor
    does not fail a genuine restart.
    """
    if not segments:
        return _fail("continuity: no segments")
    if int(write_precision) < MIN_WRITE_PRECISION:
        return _fail(f"continuity: write_precision {write_precision} < "
                     f"{MIN_WRITE_PRECISION}")
    segs = []
    for s in segments:
        label, tt = float(s[0]), np.asarray(s[1], dtype=float).ravel()
        vv = np.asarray(s[2], dtype=float).ravel() if len(s) > 2 and s[2] is not None else None
        if tt.size == 0:
            return _fail(f"continuity: segment {label:g} is empty")
        if vv is not None and vv.size != tt.size:
            return _fail(f"continuity: segment {label:g} times/values length differ")
        if not _finite(tt) or (vv is not None and not _finite(vv)):
            return _fail(f"continuity: segment {label:g} has non-finite samples")
        segs.append((label, tt, vv))
    fails = []
    notes = []
    labels = [s[0] for s in segs]
    if len(set(labels)) != len(labels):
        fails.append("duplicate segment labels " + str(sorted(
            {l for l in labels if labels.count(l) > 1})))
    segs.sort(key=lambda s: (s[0], s[1][0]))
    have_v = all(s[2] is not None for s in segs)
    if jump_atol is None:
        jump_atol = (print_quantum(np.concatenate([s[2] for s in segs]), write_precision)
                     if have_v else 0.0)

    for lab, tt, _ in segs:
        d = np.diff(tt)
        if np.any(d <= 0):
            i = int(np.argmax(d <= 0))
            fails.append(f"segment {lab:g}: time not increasing at {tt[i]:.6g} -> {tt[i + 1]:.6g}")
    all_d = np.concatenate([np.diff(s[1]) for s in segs if s[1].size > 1] or [np.array([])])
    all_d = all_d[all_d > 0]
    dt_out = expected_interval or (float(np.median(all_d)) if all_d.size else None)
    if dt_out is None:
        return _fail("continuity: cannot infer the output interval (single samples)")
    tol_t = time_tol if time_tol is not None else 0.5 * dt_out

    for i, (li, ti, _) in enumerate(segs):
        for j, (lj, tj, _) in enumerate(segs):
            if i != j and tj[0] >= ti[0] - tol_t and tj[-1] <= ti[-1] + tol_t and \
                    (tj.size < ti.size or (tj.size == ti.size and j > i)):
                fails.append(f"segment {lj:g} duplicates part of segment {li:g}")
    prov_missing = [n for n, v in (("log_start_times", log_start_times),
                                   ("log_end_time", log_end_time)) if v is None]
    if prov_missing:
        msg = "provenance NOT checked (" + ", ".join(prov_missing) + " missing)"
        (fails if require_provenance else notes).append(msg)
    if log_start_times is not None:
        starts = np.asarray(log_start_times, dtype=float)
        for lab, _, _ in segs:
            if starts.size == 0 or np.min(np.abs(starts - lab)) > tol_t:
                fails.append(f"segment {lab:g} is not a start/restart time of the case "
                             f"log (stale segment)")

    # stitch: segment k is truncated where segment k+1 starts (NEW-25)
    t_parts, v_parts = [], []
    for k, (lab, tt, vv) in enumerate(segs):
        last = k + 1 == len(segs)
        nxt = math.inf if last else float(segs[k + 1][1][0])
        keep = np.ones(tt.size, bool) if last else tt < nxt - tol_t
        if have_v and not last:
            nt, nv = segs[k + 1][1], segs[k + 1][2]
            nk = int(keep.sum())
            # tolerance scale from the earlier segment only, around the splice
            ctx = vv[max(0, nk - window):min(tt.size, nk + window)]
            loc = np.abs(np.diff(ctx))
            ref = abs(float(vv[nk - 1])) if nk else abs(float(nv[0]))
            if loc.size == 0:
                fails.append(f"splice at t={nt[0]:.6g}: segment {lab:g} has < 2 samples "
                             f"near the splice (cannot verify it)")
            scale = max([float(loc.max()) if loc.size else 0.0,
                         1e-9 * max(1.0, ref), jump_atol])
            tol_v = overlap_atol if overlap_atol is not None else jump_factor * scale
            # overlap agreement at common time stamps (first `window` of them)
            ov_t, ov_v = tt[~keep], vv[~keep]
            if ov_t.size:
                hi_i = np.clip(np.searchsorted(nt, ov_t), 0, nt.size - 1)
                lo_i = np.clip(hi_i - 1, 0, nt.size - 1)
                idx = np.where(np.abs(nt[lo_i] - ov_t) < np.abs(nt[hi_i] - ov_t), lo_i, hi_i)
                near = np.nonzero(np.abs(nt[idx] - ov_t) <= 1e-6 * max(1.0, abs(nxt)))[0]
                near = near[:window]
                if near.size:
                    diff = np.abs(ov_v[near] - nv[idx[near]])
                    if diff.max() > tol_v:
                        fails.append(f"segment {segs[k + 1][0]:g} disagrees with the "
                                     f"overlapped segment {lab:g} by {diff.max():.4g} "
                                     f"(tol {tol_v:.3g}): not the same run")
            lim = jump_factor * scale + (overlap_atol or 0.0) + jump_atol
            # splice jump: first kept value of k+1 vs linear extrapolation of k
            if nk >= 2:
                tk, vk = tt[keep][-2:], vv[keep][-2:]
                slope = (vk[1] - vk[0]) / (tk[1] - tk[0])
                pred = vk[1] + slope * (nt[0] - tk[1])
                jump = abs(nv[0] - pred)
                if jump > lim:
                    fails.append(f"jump {nv[0] - pred:+.4g} at splice t={nt[0]:.6g} "
                                 f"(local change {scale:.3g})")
            # the first changes after the splice (a step after an unchanged
            # first sample)
            post = np.abs(np.diff(nv[:window + 1]))
            if post.size and post.max() > lim:
                i = int(np.argmax(post))
                fails.append(f"jump {nv[i + 1] - nv[i]:+.4g} at t={nt[i + 1]:.6g}, "
                             f"{i + 1} sample(s) after the splice at t={nt[0]:.6g} "
                             f"(local change before the splice {scale:.3g})")
        t_parts.append(tt[keep])
        if have_v:
            v_parts.append(vv[keep])
    ts = np.concatenate(t_parts)
    d = np.diff(ts)
    if np.any(d <= 0):
        fails.append(f"stitched series goes backward/duplicates at "
                     f"{int(np.sum(d <= 0))} point(s)")
    big = d > max_gap_factor * dt_out
    if big.any():
        i = int(np.argmax(big))
        fails.append(f"gap {d[i]:.6g} s at t={ts[i]:.6g} > {max_gap_factor:g} x "
                     f"interval {dt_out:.6g} s")
    if cadence_rtol is not None:
        off = (d > 0) & ~big & (np.abs(d / dt_out - 1.0) > cadence_rtol)
        if off.any():
            i = int(np.argmax(off))
            fails.append(f"cadence change: interval {d[i]:.6g} s at t={ts[i]:.6g} vs "
                         f"nominal {dt_out:.6g} s ({int(off.sum())} intervals)")
    if scan_jumps and have_v and ts.size > 2:
        vs = np.concatenate(v_parts)
        dv = np.abs(np.diff(vs))
        floor = max(1e-9 * max(1.0, float(np.max(np.abs(vs)))), jump_atol)
        pad = np.concatenate([np.zeros(window), dv, np.zeros(window)])
        win = np.lib.stride_tricks.sliding_window_view(pad, window)
        left = win[:dv.size].max(axis=1)                 # dv[i-window:i]
        right = win[window + 1:window + 1 + dv.size].max(axis=1)   # dv[i+1:i+1+window]
        ref = np.maximum(np.maximum(left, right), floor)
        hit = np.nonzero((dv > floor) & (dv > jump_factor * ref))[0]
        if hit.size:
            i = int(hit[0])
            fails.append(f"isolated jump {vs[i + 1] - vs[i]:+.4g} at t={ts[i + 1]:.6g} "
                         f"(neighbouring change {ref[i]:.3g}; {hit.size} such)")
    if log_end_time is not None and abs(ts[-1] - log_end_time) > tol_t:
        fails.append(f"monitor ends at {ts[-1]:.6g} s, log at {log_end_time:.6g} s")
    detail = (f"{len(segs)} segment(s), {ts.size} samples over "
              f"[{ts[0]:.6g}, {ts[-1]:.6g}] s, interval {dt_out:.6g} s")
    if have_v:
        detail += f", value floor {jump_atol:.3g}"
    else:
        detail += " (times only: overlap/jump checks skipped)"
    if notes:
        detail += "; " + "; ".join(notes)
    return GateResult(not fails, detail + ("; " + "; ".join(fails) if fails else ""))


# --------------------------------------------------------------------------
# readers for the CLI
# --------------------------------------------------------------------------

_TIME_RE = re.compile(r"^Time = ([-+0-9.eE]+)\s*$", re.M)
_CO_RE = re.compile(r"Courant Number mean:\s*\S+\s+max:\s*([-+0-9.eE]+)")


def log_times(text: str) -> np.ndarray:
    """Solver step times from a log ('Time = ...' lines, not ExecutionTime)."""
    return np.array([float(m) for m in _TIME_RE.findall(text)])


def log_courant(text: str) -> np.ndarray:
    """Every 'Courant Number mean: .. max: ..' value (v2606 prints it per
    fluid region, prefixed 'Region: <name>')."""
    return np.array([float(m) for m in _CO_RE.findall(text)])


def log_start_time(text: str, dt: float | None = None) -> float:
    """Start time of the run a log belongs to: first step time minus dt."""
    t = log_times(text)
    if t.size == 0:
        raise ValueError("log has no 'Time =' lines")
    if dt is None:
        if t.size < 2:
            raise ValueError("need dt or >= 2 steps")
        dt = t[1] - t[0]
    return float(t[0] - dt)


def read_dat(path, column: int = 1):
    """(times, values) of an OpenFOAM function-object .dat file; parentheses
    are treated as separators so vector columns count as separate tokens."""
    tt, vv = [], []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        tok = s.replace("(", " ").replace(")", " ").split()
        try:
            tt.append(float(tok[0]))
            vv.append(float(tok[column]))
        except (ValueError, IndexError) as e:
            raise ValueError(f"{path}: cannot parse line {s!r}") from e
    return np.array(tt), np.array(vv)


def read_segments(pattern: str, column: int = 1):
    """[(label, times, values)] from postProcessing/<fo>/<startTime>/<file>."""
    hits = sorted(_glob.glob(pattern))
    if not hits:
        raise ValueError(f"no files match {pattern}")
    out = []
    for h in hits:
        label = float(Path(h).parent.name)
        t, v = read_dat(h, column)
        out.append((label, t, v))
    return out


def _load_array(v, base: Path):
    if isinstance(v, dict) and "csv" in v:
        import pandas as pd
        p = Path(v["csv"])
        df = pd.read_csv(p if p.is_absolute() else base / p)
        return df[v["column"]].to_numpy(dtype=float)
    return v


GATES = {
    "ledger": ledger_gate, "ledger_series": ledger_series_gate,
    "source_power": source_power_gate, "hydraulic_closure": hydraulic_closure_gate,
    "lattice": lattice_gate, "fixed_dt": fixed_dt_gate,
    "max_principle": max_principle_gate, "courant": courant_gate,
    "thin_layer_courant": thin_layer_courant_gate,
    "mean_flow": mean_flow_gate, "continuity": continuity_gate,
    "closure_phase_sign": closure_phase_sign_gate,
    "closure_stationarity": closure_stationarity_gate,
    "closure_dt_invariance": closure_dt_invariance_gate,
}

# pre-registered policy arguments a spec file may NOT set (review of prep/2.4:
# a spec that loosened tol or switched a check off was a waiver by another
# route). The exemption and scaling switches are locked too (review of
# prep/2.10, NEW-10.3): `diagnostic` passes any Courant number and lifts
# N_min, `edge_sensitivity` lifts the valve-edge rule, and `level` scales the
# Courant limits by up to 1.9x and selects the closure tolerance (D-P.4); a
# spec is judged strictly (medium Courant limits of the strictest class, the
# smallest closure tolerance, not diagnostic). housing_metrics sets them from
# the case's own caseMeta.
LOCKED = {
    "ledger": {"tol"},
    "ledger_series": {"tol", "interval_factor"},
    "source_power": {"rtol"},
    "hydraulic_closure": {"tol", "level"},
    "lattice": {"dt_ladder", "diagnostic", "edge_sensitivity"},
    "fixed_dt": {"rtol", "min_fraction"},
    "max_principle": {"undershoot", "overshoot", "running_max"},
    "courant": {"limit", "diagnostic"},
    "thin_layer_courant": {"thin_limit", "out_limit", "ladder_level", "level",
                           "diagnostic"},
    "mean_flow": {"rtol", "atol"},
    "continuity": {"max_gap_factor", "cadence_rtol", "time_tol", "jump_factor",
                   "overlap_atol", "window", "jump_atol", "require_provenance"},
    "closure_phase_sign": {"min_rel"},
    "closure_stationarity": {"tol"},
    "closure_dt_invariance": {"tol"},
}

POLICY = {
    "LEDGER_TOL": LEDGER_TOL, "LEDGER_INTERVAL_FACTOR": LEDGER_INTERVAL_FACTOR,
    "SOURCE_POWER_RTOL": SOURCE_POWER_RTOL, "HYDRAULIC_TOL": HYDRAULIC_TOL,
    "HYDRAULIC_TOL_LEVEL": HYDRAULIC_TOL_LEVEL,
    "CLOSURE_PHASE_BINS": CLOSURE_PHASE_BINS, "CLOSURE_SIGN_MIN_REL": CLOSURE_SIGN_MIN_REL,
    "CLOSURE_SUB_WINDOW_S": CLOSURE_SUB_WINDOW_S,
    "CLOSURE_STATIONARITY_PP": CLOSURE_STATIONARITY_PP,
    "CLOSURE_DT_INVARIANCE_PP": CLOSURE_DT_INVARIANCE_PP,
    "FIXED_DT_RTOL": FIXED_DT_RTOL, "FIXED_DT_MIN_FRACTION": FIXED_DT_MIN_FRACTION,
    "UNDERSHOOT_K": UNDERSHOOT_K, "OVERSHOOT_K": OVERSHOOT_K,
    "COURANT_LIMIT": COURANT_LIMIT, "COURANT_THIN_TESTED": COURANT_THIN_TESTED,
    "COURANT_THIN_MARGIN": COURANT_THIN_MARGIN, "COURANT_THIN_B": COURANT_THIN_B,
    "COURANT_THIN_B_CLASS": COURANT_THIN_B_CLASS, "COURANT_BULK_B_CLASS": COURANT_BULK_B_CLASS,
    "COURANT_OUT_B_CLASS": COURANT_OUT_B_CLASS,
    "COURANT_LEVEL_BULK": COURANT_LEVEL_BULK, "COURANT_LEVEL_THIN": COURANT_LEVEL_THIN,
    "MEAN_FLOW_RTOL": MEAN_FLOW_RTOL,
    "MEAN_FLOW_ZERO_ATOL": MEAN_FLOW_ZERO_ATOL, "N_MIN_ALLOWED": N_MIN_ALLOWED,
    "N_PHASE_ALLOWED": N_PHASE_ALLOWED, "VALVE_EDGE_MAX": VALVE_EDGE_MAX,
    "MAX_GAP_FACTOR": MAX_GAP_FACTOR, "CADENCE_RTOL": CADENCE_RTOL,
    "JUMP_FACTOR": JUMP_FACTOR, "JUMP_WINDOW": JUMP_WINDOW,
    "MIN_WRITE_PRECISION": MIN_WRITE_PRECISION,
}


def _check_args(name: str, kw: dict) -> list:
    """Problems with a spec entry's arguments (empty list = fine)."""
    sig = inspect.signature(GATES[name])
    params = sig.parameters
    probs = []
    unknown = sorted(set(kw) - set(params))
    if unknown:
        probs.append(f"unknown argument(s) {unknown}")
    locked = sorted(set(kw) & LOCKED[name])
    if locked:
        probs.append(f"pre-registered argument(s) {locked} cannot be set from a spec")
    missing = sorted(n for n, prm in params.items()
                     if prm.default is inspect.Parameter.empty and n not in kw)
    if missing:
        probs.append(f"missing required argument(s) {missing}")
    return probs


def run_spec(spec: dict, base: Path = Path(".")) -> list:
    """Run every gate of a spec dict; returns [(label, GateResult)]. A gate
    whose entry is invalid gets an InputError result; the others still run."""
    results = []
    entries = spec.get("gates", [])
    if not isinstance(entries, list):
        return [("spec", InputError(False, "input error: 'gates' must be a list"))]
    for i, g in enumerate(entries):
        if not isinstance(g, dict) or "gate" not in g:
            results.append((f"gates[{i}]", InputError(
                False, "input error: entry needs a 'gate' name")))
            continue
        g = dict(g)
        name = g.pop("gate")
        label = g.pop("label", name)
        if name not in GATES:
            results.append((label, InputError(False, f"input error: unknown gate {name!r}")))
            continue
        probs = _check_args(name, g)
        if probs:
            results.append((label, InputError(False, f"input error ({name}): "
                                              + "; ".join(probs))))
            continue
        try:
            kw = {}
            for k, v in g.items():
                if k == "segments":
                    v = [(s["label"], _load_array(s["times"], base),
                          _load_array(s.get("values"), base)) for s in v]
                elif k == "fluxes_out":
                    v = {kk: _load_array(vv, base) for kk, vv in v.items()}
                else:
                    v = _load_array(v, base)
                kw[k] = v
        except (OSError, KeyError, TypeError, ValueError) as e:
            results.append((label, InputError(False, f"input error ({name}): "
                                              f"cannot load data: {e!r}")))
            continue
        try:
            results.append((label, GATES[name](**kw)))
        except TypeError as e:
            results.append((label, InputError(False, f"input error ({name}): {e}")))
        except ValueError as e:
            results.append((label, _fail(f"{name}: invalid input: {e}")))
    return results


def _report(results) -> int:
    for name, r in results:
        tag = "INPUT-ERROR" if isinstance(r, InputError) else ("PASS" if r.ok else "FAIL")
        print(f"{tag} {name}: {r.detail}")
    n_fail = sum(not r.ok for _, r in results)
    n_input = sum(isinstance(r, InputError) for _, r in results)
    print(f"gates: {len(results) - n_fail}/{len(results)} passed"
          + (f", {n_input} input error(s)" if n_input else ""))
    if n_input or not results:
        return 2
    return 1 if n_fail else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Campaign hard gates (zero waivers).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("spec")
    s.add_argument("file")
    s.add_argument("--out")
    s = sub.add_parser("courant")
    s.add_argument("log")
    s.add_argument("--diagnostic", action="store_true",
                   help="a DIAGNOSTIC=1 ladder/coupling rung: reported, not gated (D-F.1)")
    s = sub.add_parser("fixed-dt")
    s.add_argument("log")
    s.add_argument("--dt", type=float, required=True)
    s.add_argument("--t-start", type=float)
    s = sub.add_parser("lattice")
    s.add_argument("--dt", type=float, required=True)
    for k in ("period", "on-time", "sample", "valve-edge"):
        s.add_argument(f"--{k}", type=float)
    s.add_argument("--n-min", type=int, default=40, choices=N_MIN_ALLOWED)
    s.add_argument("--n-phase-min", type=int, default=8, choices=N_PHASE_ALLOWED)
    s.add_argument("--diagnostic", action="store_true")
    s.add_argument("--edge-sensitivity", action="store_true")
    s = sub.add_parser("continuity")
    s.add_argument("--glob", required=True)
    s.add_argument("--column", type=int, default=1)
    s.add_argument("--log", action="append", required=True,
                   help="every solver log of the case (provenance is required)")
    s.add_argument("--dt", type=float)
    s.add_argument("--write-precision", type=int, default=DEFAULT_WRITE_PRECISION)
    s.add_argument("--scan-jumps", action="store_true")
    s = sub.add_parser("ledger")
    s.add_argument("--e-injected", type=float, required=True)
    s.add_argument("--term", action="append", required=True, help="name=J")
    s = sub.add_parser("max-principle")
    s.add_argument("--csv", required=True)
    s.add_argument("--t-in", type=float, required=True)
    s.add_argument("--fluid-min", required=True)
    s.add_argument("--fluid-max", required=True)
    s.add_argument("--solid-max", required=True)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "spec":
            p = Path(a.file)
            spec = json.loads(p.read_text())
            if not isinstance(spec, dict):
                raise ValueError("spec must be a JSON object")
            res = run_spec(spec, p.parent)
            if a.out:
                Path(a.out).write_text(json.dumps(
                    {"case": spec.get("case"), "policy": POLICY,
                     "gates": [{"gate": n, "ok": bool(r.ok),
                                "input_error": isinstance(r, InputError),
                                "detail": r.detail} for n, r in res]},
                    indent=2) + "\n")
        elif a.cmd == "courant":
            res = [("courant", courant_gate(log_courant(Path(a.log).read_text()),
                                            diagnostic=a.diagnostic))]
        elif a.cmd == "fixed-dt":
            res = [("fixed_dt", fixed_dt_gate(log_times(Path(a.log).read_text()),
                                              a.dt, a.t_start))]
        elif a.cmd == "lattice":
            res = [("lattice", lattice_gate(a.dt, a.period, a.on_time, a.sample,
                                            a.valve_edge, a.n_min, a.n_phase_min,
                                            diagnostic=a.diagnostic,
                                            edge_sensitivity=a.edge_sensitivity))]
        elif a.cmd == "continuity":
            segs = read_segments(a.glob, a.column)
            texts = [Path(l).read_text() for l in a.log]
            starts = [log_start_time(t, a.dt) for t in texts]
            end = float(max(log_times(t)[-1] for t in texts))
            res = [("continuity", continuity_gate(segs, starts, end,
                                                  scan_jumps=a.scan_jumps,
                                                  write_precision=a.write_precision))]
        elif a.cmd == "ledger":
            terms = {}
            for t in a.term:
                k, sep, v = t.partition("=")
                if not sep:
                    raise ValueError(f"--term {t!r} is not name=J")
                terms[k] = float(v)
            res = [("ledger", ledger_gate(a.e_injected, terms))]
        else:  # max-principle
            import pandas as pd
            df = pd.read_csv(a.csv)
            res = [("max_principle", max_principle_gate(
                df[a.fluid_min], df[a.fluid_max], a.t_in, df[a.solid_max]))]
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as e:
        print(f"gates: input error: {e}", file=sys.stderr)
        return 2
    return _report(res)


if __name__ == "__main__":
    sys.exit(main())
