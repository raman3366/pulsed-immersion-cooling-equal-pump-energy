#!/usr/bin/env python3
"""Hard gates of a K0 housing case from its monitors (WP 2.5c; rule 4,
PROJECT.md §9, plan D-F, D-J, D-P, D-T, NEW-15, NEW-25).

Reads the housing monitors (scripts/model/housing_physics.monitors), the
solver logs and system/caseMeta.yaml, builds the arguments of the
scripts/post/gates.py gates and runs them. No tolerance is set here: the
gates' pre-registered constants apply (gates.LOCKED).

Gates evaluated
- lattice        caseMeta lattice (dt, period, ON time, valve edge, sampler
                 interval, N_min, n_phi) - the generator already asserted it
- fixed_dt       every solver log (one per stage)
- courant        the Courant rule by region (PROJECT.md §2.5 D-F.1, per
                 waveform class D-F.5, regions §2.7 D-F.9): CoMinMax_fluid's
                 CoThin <= the class's B (thin layer: cap-plane band +
                 T_head), CoOut <= its B_out (outflow column) and CoBulk
                 (claim region) <= 1, at the level-scaled limits of
                 gates.thin_layer_courant_gate, and every log line <=
                 max(B, B_out, 1); a case rendered before the D-F.9 monitors
                 existed (GDT1, GDT2: CoBulk held R_out and the sleeve edges
                 together) takes the strict global gate (every log line <=
                 1), and its D-F.9 regions are reported as bounds from the
                 location of each sample's bulk maximum (legacy_regions,
                 with the D-F.9 thin layer T from the located CoThin band
                 maximum, thin_T). ladder_courant gives a ladder rung's
                 {co_thin, co_out} for recipe_rules.courant_bounds (the
                 unlocated every-step log maximum never enters a bound).
                 CoBulkSleeves (the claim region under D-F.1's sleeves) is
                 the D-F.9 (1) sensitivity column, reported, never gated.
                 DIAGNOSTIC runs (ladder, coupling) are reported, not gated.
                 Where the maxima sit is reported
- max_principle  T_in - 0.05 K <= T_fluid(t) <= running max T_solid(t)
                 (volFieldValue min/max of the fluid; max over the six cells)
                 on the heated samples t > t_pre (not for isothermal runs);
                 D-F.10: the bound unchanged for every class and level, no
                 class tolerance, no exclusion window; reported beside it
                 (info max_principle_detail): the worst sample's time (never
                 its index), location (TfluidMinMax), waveform phase,
                 interstice flag, the span below -20 mK, the like-for-like
                 minima per window (all / interstice core) and the overlap
                 with every D-F.4 candidate of the primary window;
                 undershoot_attribution applies D-F.10 (3) to a rung and
                 its half-Dt partner (Dt-controlled or h transport)
- mean_flow      Qbar = -sum(phi/rho) at the inlet over whole periods of
                 [max(ramp, t_pre), t_end] ([ramp, t_pre] for a
                 pre-circulation-only run); a run shorter than that (smokes)
                 is judged as delivery against the ramped intent over its
                 whole periods
- ledger         heated stage [t_pre, t_end] in the solver's discrete form
                 (heated runs only):
                   E_inj = sum_cells sum_steps P(t_n) dt   (the solver
                           evaluates the source at each step's new time;
                           P from caseMeta = the fvOptions Function1)
                 stores (backward ddt: change of the invariant H + dH/2,
                 dH = one step's increment; dH = 0 at a restart,
                 scheme_storage):
                   dH_cells   sum_i rho_cell (volIntegrate h)_i
                   dH_fluid   weightedVolIntegrate(rho h)
                   dKE        volIntegrate(rho |U|^2/2)
                   restart_loss  dH^r/2 of every store dropped by a restart
                              inside the window (resumes)
                 step sums sum_n x^n dt (new time level, step_sum):
                   out_h      sum_ports sum(phi h)   (advection out)
                   out_K      sum_ports sum(phi |U|^2/2)
                   out_cond   k dT/dx over the inlet (conduction out
                              through the fixedValue inlet; cell-gradient
                              estimate)
                   out_walls  -wallHeatFlux on the housing walls
                              (adiabatic: 0)
                   src_rhoUg  -int rho U.g dV   (energy-equation source of
                              v2606 chtMultiRegionFoam EEqn.H)
                   src_dpdt   -int dpdt dV
                 tilt       -(enthalpy added by the tilt seed, tilt.json)
                 The naive form (Delta H, trapezoids) is reported beside it.
- source_power   D-J (i): per cell, the power the solver injects
                 (value / VDash x meshed volume; VDash from the log's
                 "- selected ... cell(s) with volume", the meshed volume
                 from the monitor header) equals the model power
- hydraulic_closure  D-P validation column (pump_energy.measured):
                 measured port work E_port = int -sum_ports sum phi
                 (p*/rho + |U|^2/2) dt, p* = p - rho_in gh, against
                 Phi - W_b + dKE (Phi = viscous dissipation, W_b = buoyancy
                 (chimney) work int (rho - rho_in) U.g dV) over the
                 developed window, whole periods, on max(|Phi|, |W_b|),
                 at the tolerance of the case's grid level (D-P.4).
                 Per waveform class (D-P.7, 2026-10-02): GATED on isothermal
                 windows of STEADY flow only - the isothermal run's window
                 ('hydraulic_closure') and the pre-circulation [ramp, t_pre]
                 of every heated run ('hydraulic_closure[pre]'). SINE
                 isothermal windows report the closure and gate its
                 necessary checks: 'closure_sign' (every phase bin's
                 per-period residual >= -0.5 % of the per-period Phi) and
                 'closure_stationarity' (>= 4 s sub-windows within +-2
                 points of the window's r/Phi; reported when the window
                 holds fewer than two); the dt invariance is the G-dt pair
                 check (pump_energy sine-dt). SQUARE waves report it split
                 by valve segment. Heated windows and runs without a
                 developed window report it. D-P.9 (2026-10-06): a case
                 with the face-based Phi monitor (PHI_FACE=1) also reports
                 the closure against the face-based Phi (pump_energy.
                 face_closure), never gated: the cell-centred Phi stays
                 the gate's instrument. A face-based monitor that fails
                 its own checks is printed as such after the registered
                 report, which it never suppresses, and the exit code is
                 then at least 2 (an input error of that instrument)
- continuity     key monitors: segments vs the solver logs (NEW-25)
- precirc        NEW-14.1 (c) / NEW-14.2: the stationarity rule of
                 precirc_check.py on the case's pre-circulation [0, t_pre]
                 (its own monitors, or the shared state's through
                 INIT_FROM), recomputed here; a GATE on every claim-eligible
                 run with a pre-circulation (claim eligibility requires
                 PRECIRC_CHECK=strict), reported for the others (smokes,
                 DIAGNOSTIC rungs, closurecal), and not applicable to zero
                 flow. Allrun.heat's precirc_check.json is reported beside it

Claim eligibility (caseMeta claim_eligible / claim_notes, housing_io.
claim_status) is reported and written to the JSON top level; it is not a
gate here (diagnostic runs are scored too), but scripts/post/claims.py
refuses a filled claim that cites a run that is not claim-eligible.

Reported beside the gates: per-cell interface balance (E_inj,i against
dH_i + int qInt_i dt, fluid-side flux), module-mean and surface T, fluid
T range, port backflow fraction, the pre-circulation check (if present),
and (D-G.1 (4)) the outlet backflow per window ([0, t_pre], [t_pre,
t_pre + 70 s], [t_pre + 70 s, t_end], [t_pre, t_end]) with the riser
exchange-flow bursts (|U|max > 60 mm/s at z > 200 mm: time, peak |U|,
peak KE against the developed median).

Exit 0 = every gate passes, 1 = a gate fails, 2 = input error.

Usage: housing_metrics.py CASE [--json OUT] [--steady-ref STEADY_CASE]
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "model"))
import gates as G  # noqa: E402
import housing_io as HIO  # noqa: E402
import housing_physics as HP  # noqa: E402
import precirc_check as PC  # noqa: E402
import pump_energy as PE  # noqa: E402

N_CELLS = HP.N_CELLS
CONTINUITY_MONITORS = [("fluid", "rhoh_fluid"), ("fluid", "Tmax_fluid"),
                       ("fluid", "hphi_outlet"), ("fluid", "qInt_cell0"),
                       ("cell0", "rhoh_cell0"), ("cell0", "Tsurfmax_cell0")]


def _trapz(y, t):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def _window(t, y, t0, t1, what):
    """Samples of y on [t0, t1]; both ends must be samples (to 1e-9 s)."""
    tol = 1e-9 * max(1.0, abs(t1))
    m = (t >= t0 - tol) & (t <= t1 + tol)
    tt, yy = t[m], np.asarray(y)[m]
    if tt.size < 2 or abs(tt[0] - t0) > tol or abs(tt[-1] - t1) > tol:
        raise ValueError(f"{what}: no samples at both ends of [{t0:g}, "
                         f"{t1:g}] s (have {tt[:1]}..{tt[-1:]})")
    return tt, yy


def _at(t, y, t0, what):
    i = int(np.argmin(np.abs(t - t0)))
    if abs(t[i] - t0) > 1e-9 * max(1.0, abs(t0)):
        raise ValueError(f"{what}: no sample at t = {t0:g} s")
    return float(np.asarray(y)[i])


def _col(case, region, name, pattern=None):
    return HIO.column(case, region, name, pattern)


def solver_steps(case):
    """(logs, runs): runs = [(start, step times, log)], one per solver run,
    in execution order (the logs in run order, housing_io.LOG_NAMES, and
    each log's runs in the order they were appended). The step times are
    every step the run took, including steps a later run replays
    (effective_steps removes those)."""
    logs = HIO.read_logs(case)
    if not logs:
        raise FileNotFoundError(
            f"{case}: no solver log (expected one of "
            f"{', '.join(n + '[.zst]' for n in HIO.LOG_NAMES)}); a case "
            "directory with only case.env and system/ has not run")
    runs = []
    for lg in logs:
        for k, i0 in enumerate(lg.run_first_index):
            i1 = (lg.run_first_index[k + 1] if k + 1 < len(lg.run_first_index)
                  else lg.step_times.size)
            runs.append((lg.run_starts[k], lg.step_times[i0:i1], lg))
    return logs, runs


def effective_steps(runs, dt):
    """The steps that made the stitched solution, as [(start, step times,
    log)] in time order (review of WP 2.5c, MUST-FIX 1).

    A run killed after its last write and resumed from latestTime leaves
    in the log its own steps past the write time and then the same steps
    again from the resumed run (one to two write intervals: WRITE_S 30 s,
    and the run driver renames an incomplete latest time). The stitched
    monitors count those steps once (housing_io truncates a segment at the
    start of the next), so the step list must too: each run keeps its steps
    t <= (start of the next run) + dt/2, the next run's start being the
    written time it resumed from. A run whose steps a later run replays
    entirely is dropped (it contributes no step and no restart).

    Runs must start in non-decreasing time order in execution order: a run
    that restarts before an earlier run's start would supersede that run's
    later steps and monitor segments, which neither this cut nor the
    segment stitching (sorted by start time) handles - an input error."""
    starts = [float(r[0]) for r in runs]
    for k in range(1, len(starts)):
        if starts[k] < starts[k - 1] - 0.5 * dt:
            raise ValueError(
                f"solver run {k} starts at {starts[k]:.10g} s, before run "
                f"{k - 1} ({starts[k - 1]:.10g} s): a restart from an earlier "
                "time than a previous run's start is not a resume; the "
                "stitched monitors and the step list would disagree")
    out = []
    for k, (s, steps, lg) in enumerate(runs):
        steps = np.asarray(steps, float)
        if k + 1 < len(runs):
            steps = steps[steps <= starts[k + 1] + 0.5 * dt]
        if steps.size:
            out.append((s, steps, lg))
    t = np.concatenate([r[1] for r in out]) if out else np.array([])
    if t.size > 1 and np.any(np.diff(t) <= 0.5 * dt):
        i = int(np.argmax(np.diff(t) <= 0.5 * dt))
        raise ValueError(f"solver steps not increasing after the resume cut "
                         f"at {t[i]:.10g} -> {t[i + 1]:.10g} s")
    return out


def injected_energy(meta, step_times):
    """Discrete injected energy per cell: sum over the steps of the heated
    stage of P(t_n - t_pre) dt (the source window of caseMeta.heat)."""
    h = meta["heat"]
    hm = HP.HeatModel(h["kind"], coeffs=h.get("coeffs"), table=h.get("table"),
                      scale=h.get("scale", 1.0), const_w=h.get("const_w", 0.0))
    dt = meta["lattice"]["dt"]
    t = np.asarray(step_times, float)
    t0, dur = h["t_on"], h["duration"]
    on = (t >= t0) & (t <= t0 + dur)
    p = np.where(on, hm.power(t - meta["t_pre"]), 0.0)
    return float(np.sum(p) * dt), t, p


def selected_volumes(logs):
    """VDash per heat option ('- selected N cell(s) with volume V' of
    cellSetOption::setVol), one list per solver run that printed any (a
    run starts at its banner's 'Create time': the stage-1 log, the heat
    stage's log and every resume appended to a log each print them)."""
    runs, cur = [], []
    rx = re.compile(r"- selected (\d+) cell\(s\) with volume ([-+0-9.eE]+)")
    for lg in logs:
        with HIO.open_text(lg.path) as f:
            for line in f:
                if line.startswith("Create time"):
                    if cur:
                        runs.append(cur)
                    cur = []
                    continue
                m = rx.search(line)
                if m:
                    cur.append((int(m.group(1)), float(m.group(2))))
        if cur:
            runs.append(cur)
        cur = []
    return runs


# --------------------------------------------------------------------------
# discrete (scheme-consistent) ledger bookkeeping
# --------------------------------------------------------------------------

def _sample_index(t, x):
    i = np.flatnonzero(np.abs(t - x) <= 1e-9 * max(1.0, abs(x)))
    return int(i[0]) if i.size else None


def step_increment(t, y, t_at, dt):
    """Per-step increment of a monitored store at the step ending at t_at:
    the change over the last monitor interval ending at t_at divided by its
    number of steps (exact for every-step samples; the mean increment of
    the last M steps otherwise). nan without a sample at t_at and one
    before it."""
    i = _sample_index(t, t_at)
    if i is None or i == 0:
        return math.nan
    n = (t[i] - t[i - 1]) / dt
    return float((y[i] - y[i - 1]) / n)


def scheme_storage(t, y, t0, t1, dt, starts, what):
    """Storage of a monitored store H over the steps ending in (t0, t1], as
    the solver's backward (BDF2) ddt at fixed dt accumulates it:

        sum_n (3 H^n - 4 H^(n-1) + H^(n-2)) / 2
            = H^N - H^0 + (dH^N - dH^0) / 2        (dH^n = H^n - H^(n-1)),

    i.e. the change of the discrete invariant H + dH/2. dH^0 = 0 when a
    solver run starts at t0: a restart (Allrun.heat, INIT_FROM, a resume)
    starts the backward scheme with the old-old h equal to the old h, h
    having no written old time (NEW-25; measured: the 2-step storage ratio
    of the P00b steady smoke is 0.7785 against the recursion's 14/18). A
    restart at t_r inside the window drops dH^r/2 of the invariant, which
    the solver never stores: ``restart_loss``. Returns a dict."""
    tt, yy = _window(t, y, t0, t1, what)
    d_end = step_increment(t, y, t1, dt)
    restart0 = any(abs(s_ - t0) <= 0.5 * dt for s_ in starts)
    d_start = 0.0 if restart0 else step_increment(t, y, t0, dt)
    inner = [s_ for s_ in starts if t0 + 0.5 * dt < s_ < t1 - 0.5 * dt]
    d_r = []
    for s_ in inner:
        j = np.flatnonzero(t <= s_ + 1e-9 * max(1.0, s_))
        d_r.append(step_increment(t, y, float(t[j[-1]]), dt) if j.size
                   else math.nan)
    dH = float(yy[-1] - yy[0])
    return dict(dH=dH, d_start=d_start, d_end=d_end, restart_at_t0=restart0,
                restarts=inner, d_restart=d_r,
                storage=dH + 0.5 * (d_end - d_start),
                restart_loss=0.5 * float(np.sum(d_r)) if d_r else 0.0)


def step_sum(t, y, t0, t1, dt, what):
    """The solver's step sum sum_n y^n dt over the steps ending in (t0, t1]
    (fluxes and sources at each step's new time level: backward and the
    implicit coupling evaluate them there). Every-step samples give it
    exactly as trapezoid + dt/2 (y(t1) - y(t0)); coarser samples (every M
    steps) to O(M dt^2 y'')."""
    tt, yy = _window(t, y, t0, t1, what)
    return _trapz(yy, tt) + 0.5 * dt * float(yy[-1] - yy[0])


def ledger(case, meta, runs, all_steps):
    """The D-F energy ledger of the heated stage [t_pre, t_end] in the
    solver's own discrete form (WP 2.5c completion; D-F ledger row of
    docs/model/housing-physics.md section 11):

        E_inj = sum_n sum_cells P(t_n - t_pre) dt
              = S(cells) + S(fluid rho h) + S(KE) + restart losses
                + sum_n [ sum_ports (phi h + phi K) + Q_cond,in - Q_walls
                          - int (rho U.g + dpdt) dV ]^n dt  - E_tilt

    S = scheme_storage (backward: change of H + dH/2), the bracket =
    step_sum. The naive form (Delta H, trapezoids) leaves about half of one
    step's stored energy (the invariant's dH^N/2 at the window end: +2.4 %
    of the 20-step P00b smokes; 0.5 f dt / t_dis, about 3e-6, of a 900 s
    discharge at 6.25 ms with a stored fraction f of about 0.9); it is
    reported beside the gate. `runs` are the effective runs (effective_steps)
    and `all_steps` their steps: a crash-resume's replayed steps once.
    Returns (GateResult, info dict)."""
    dt = meta["lattice"]["dt"]
    t0, t_end = float(meta["t_pre"]), float(meta["t_end"])
    starts = [float(r[0]) for r in runs]
    heat_steps = all_steps[(all_steps > t0 + 0.5 * dt)
                           & (all_steps <= t_end + 0.5 * dt)]
    e_cell, _, _ = injected_energy(meta, heat_steps)
    e_inj = N_CELLS * e_cell
    import yaml
    cell = yaml.safe_load((HP.REPO / "data/reference/properties.yaml")
                          .read_text())["cell"]
    rho_c = float(cell["density"])
    tilt = case / "tilt.json"
    tilt_rec = json.loads(tilt.read_text()) if tilt.exists() else None
    terms, naive, per_cell, stores = {}, {}, {}, {}
    dh_cells = dh_cells_naive = loss = 0.0
    for i in range(N_CELLS):
        th, hh = _col(case, f"cell{i}", f"rhoh_cell{i}")
        st = scheme_storage(th, rho_c * hh, t0, t_end, dt, starts,
                            f"rhoh_cell{i}")
        tq_, qi = _col(case, "fluid", f"qInt_cell{i}")
        eq = step_sum(tq_, qi, t0, t_end, dt, f"qInt_cell{i}")
        e_tilt = (float(tilt_rec["regions"].get(f"cell{i}", {}).get("dH_J", 0))
                  if tilt_rec else 0.0)
        r = e_cell + e_tilt - st["storage"] - st["restart_loss"] - eq
        per_cell[f"cell{i}"] = dict(
            E_inj=e_cell, E_tilt=e_tilt, dH=st["dH"], storage=st["storage"],
            restart_loss=st["restart_loss"], E_interface=eq,
            residual_rel=r / e_cell if e_cell else math.nan)
        dh_cells += st["storage"]
        dh_cells_naive += st["dH"]
        loss += st["restart_loss"]
        stores[f"cell{i}"] = st
    terms["dH_cells"] = dh_cells
    naive["dH_cells"] = dh_cells_naive
    for key, name in (("dH_fluid", "rhoh_fluid"), ("dKE", "KE_fluid")):
        ts_, v = _col(case, "fluid", name)
        st = scheme_storage(ts_, v, t0, t_end, dt, starts, name)
        terms[key] = st["storage"]
        naive[key] = st["dH"]
        loss += st["restart_loss"]
        stores[name] = st
    for key, name in (("out_h", "hphi"), ("out_K", "kPhi")):
        tot = tot_n = 0.0
        for port in ("inlet", "outlet"):
            tp, vp = _col(case, "fluid", f"{name}_{port}")
            tot += step_sum(tp, vp, t0, t_end, dt, f"{name}_{port}")
            tot_n += _trapz(*reversed(_window(tp, vp, t0, t_end, name)))
        terms[key], naive[key] = tot, tot_n
    for key, name, sgn in (("out_cond", "qCond_inlet", 1.0),
                           ("out_walls", "qWall_fluid", -1.0),
                           ("src_rhoUg", "rhoUg_fluid", -1.0),
                           ("src_dpdt", "dpdt_fluid", -1.0)):
        ts_, v = _col(case, "fluid", name)
        terms[key] = sgn * step_sum(ts_, v, t0, t_end, dt, name)
        naive[key] = sgn * _trapz(*reversed(_window(ts_, v, t0, t_end, name)))
    if loss:
        terms["restart_loss"] = loss
    if tilt_rec:
        terms["tilt"] = naive["tilt"] = -float(tilt_rec["dH_total_J"])
    if e_inj > 0:
        gate = G.ledger_gate(e_inj, terms)
    else:
        gate = G.InputError(False, "ledger: no injected energy in the window "
                                   "(unheated)")
    res_ = e_inj - sum(terms.values())
    res_n = e_inj - sum(naive.values())
    return gate, dict(
        window=[t0, t_end], form="backward (BDF2) discrete",
        E_inj=e_inj, terms=terms, residual=res_,
        residual_rel=res_ / e_inj if e_inj else math.nan,
        naive_terms=naive, naive_residual=res_n,
        naive_residual_rel=res_n / e_inj if e_inj else math.nan,
        restart_at_t_pre=any(abs(s_ - t0) <= 0.5 * dt for s_ in starts),
        restarts_inside=[s_ for s_ in starts
                         if t0 + 0.5 * dt < s_ < t_end - 0.5 * dt],
        store_increments={k: dict(d_start=v["d_start"], d_end=v["d_end"],
                                  d_restart=v["d_restart"])
                          for k, v in stores.items()},
        per_cell=per_cell)


# --------------------------------------------------------------------------
# D-F.10 (4): the maximum-principle detail (moved here from the GDT2
# deep-dive, scripts/gdt/gdt2_analysis.py, which imports it back)
# --------------------------------------------------------------------------

UNDERSHOOT_REPORT_K = 0.020       # D-F.10 (4): the span below -20 mK is reported
MP_EARLY_S = 100.0                # D-F.10 (3): the early heated window of the
                                  # like-for-like minima, (t_pre, t_pre + 100 s]
                                  # (GDT2's C1 interstice episode: 30-130 s)
MP_ATTRIBUTION_RATIO = 0.5        # D-F.10 (3): Dt-controlled iff the Dt/2 rung's
                                  # interstice undershoot in that window is within
                                  # the bound and at most half the Dt rung's
HEAT_EARLY_S = 70.0               # D-G.1 (4) / D-F.10: the early heated window
                                  # [t_pre, t_pre + 70 s] (= [30, 100] s on the 30 s
                                  # ladder protocol of GDT1/GDT2)
BURST_U = 0.06                    # D-G.1 (4): a riser burst = |U|max > 60 mm/s ...
BURST_Z = 0.200                   # ... at z > 200 mm (the riser above the module)
BURST_GAP_SAMPLES = 40            # samples apart that still belong to one burst
KE_REF_DELAY_S = 270.0            # the developed KE reference starts 270 s after
                                  # heat-on ([300, 930] s on the 30 s protocol)


def module_geometry_mm(g=None):
    """(cell axes [(x, y)], pitch, cell radius) in mm from the geometry YAML."""
    g = HP.KB.load_geometry() if g is None else g
    axes = [(1e3 * x, 1e3 * y) for x, y in HP.cell_axes(g)]
    return axes, float(g["D"] + g["gap"]), 0.5 * float(g["D"])


def interstice_centres(axes=None, pitch=None):
    """Centres of the four-cell interstices [mm]: one per unit square of
    axes (the K0 module: (-10, 0) and (10, 0))."""
    if axes is None or pitch is None:
        a, p, _ = module_geometry_mm()
        axes = a if axes is None else axes
        pitch = p if pitch is None else pitch
    s = {(round(x, 6), round(y, 6)) for x, y in axes}
    return [(x + 0.5 * pitch, y + 0.5 * pitch) for x, y in sorted(s)
            if {(round(x + pitch, 6), y), (x, round(y + pitch, 6)),
                (round(x + pitch, 6), round(y + pitch, 6))} <= s]


def _interstice_r_mm():
    _, pitch, r_cell = module_geometry_mm()
    return math.sqrt(2.0) * 0.5 * pitch - r_cell


# the interstice core: horizontally within the centre-to-surface distance
# (sqrt(2) pitch / 2 - R = 5.14 mm) of a four-cell interstice centre
INTERSTICE_R_MM = _interstice_r_mm()


def in_interstice(loc_mm, centres=None, r_mm=INTERSTICE_R_MM):
    """True for points (n, 3) [mm] inside a four-cell interstice core."""
    c = np.asarray(interstice_centres() if centres is None else centres, float)
    p = np.atleast_2d(np.asarray(loc_mm, float))
    d = np.sqrt((p[:, 0, None] - c[:, 0]) ** 2 + (p[:, 1, None] - c[:, 1]) ** 2)
    return d.min(axis=1) < r_mm


def undershoot_minima(t, dev_mk, loc_mm, windows):
    """D-F.10 (3): worst fluid undershoot (T_min - T_in, mK) per time window,
    over all samples and over the samples whose minimum lies in a four-cell
    interstice core, with time and location: the like-for-like comparison of
    two runs (same place, same window). A window or subset without samples
    gives None."""
    t = np.asarray(t, float)
    dev = np.asarray(dev_mk, float)
    loc = np.atleast_2d(np.asarray(loc_mm, float))
    inter = in_interstice(loc)
    out = {}
    for a, b in windows:
        w = (t > a) & (t <= b)
        row = {}
        for lab, k in (("all", w), ("interstice", w & inter)):
            if not k.any():
                row[lab] = None
                continue
            j = int(np.flatnonzero(k)[np.argmin(dev[k])])
            row[lab] = dict(worst_mK=float(dev[j]), t=float(t[j]),
                            location_mm=[round(float(v), 2) for v in loc[j]],
                            n=int(k.sum()))
        out[f"{a:g}-{b:g}"] = row
    return out


def primary_windows(meta):
    """D-F.4 candidates of the primary uniformity window, whole periods
    ending at t_end: {'last 20 %': (t0, t_end), ...} (campaign_stats.
    WINDOWS_PCT; the window itself is chosen at FREEZE, so every candidate
    is reported)."""
    import campaign_stats as CS
    t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
    P = meta.get("period_s")
    dur = t_end - t_pre
    out = {}
    for w in CS.WINDOWS_PCT:
        span = w / 100.0 * dur
        if P:
            span = math.floor(span / P + 1e-9) * P
        out[f"last {w:g} %"] = (t_end - span, t_end)
    return out


def max_principle_detail(case, meta, tf, fmin, heated):
    """D-F.10 (4): the worst heated sample of the fluid undershoot below T_in
    - its time (never an index), location (TfluidMinMax), waveform phase,
    whether it lies in a four-cell interstice core - the span of heated
    samples below -20 mK, the like-for-like minima (all / interstice) per
    window and the overlap with every D-F.4 candidate of the primary window.
    The bound itself is gated by gates.max_principle_gate, unchanged."""
    t_in = float(meta["T_inlet"])
    t, dev = tf[heated], (fmin[heated] - t_in)
    if t.size == 0:
        return dict(note="no heated samples")
    i = int(np.argmin(dev))
    P = meta.get("period_s")
    dt = float(meta["lattice"]["dt"])
    rep = dict(worst_K=float(dev[i]), t=float(t[i]),
               t_after_heat_on=float(t[i] - float(meta["t_pre"])))
    if P:
        spp = int(round(P / dt))
        rep["phase"] = float((PE.step_index(t[i], dt) % spp) / spp)
    below = dev < -UNDERSHOOT_REPORT_K
    rep["below_20mK"] = dict(
        n=int(below.sum()),
        t_first=float(t[below].min()) if below.any() else None,
        t_last=float(t[below].max()) if below.any() else None)
    loc = None
    try:
        mm = HIO.minmax(case, "fluid", "TfluidMinMax")["T"]
        if mm.t.size == tf.size and np.allclose(mm.t, tf):
            loc = 1e3 * np.asarray(mm.min_loc, float)[heated]
        else:
            rep["location_note"] = "TfluidMinMax and Tmin_fluid not on one time base"
    except (FileNotFoundError, KeyError) as e:
        rep["location_note"] = f"no TfluidMinMax monitor ({e!s})"
    if loc is not None:
        rep["location_mm"] = [round(float(v), 3) for v in loc[i]]
        rep["in_interstice"] = bool(in_interstice(loc[i])[0])
        t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
        early = min(t_end, t_pre + MP_EARLY_S)              # GDT2: (30, 130] s
        wins = [(t_pre, early), (early, t_end), (t_pre, t_end)]
        rep["minima_mK"] = undershoot_minima(t, 1e3 * dev, loc,
                                             [w for w in wins if w[1] > w[0]])
        # the D-F.10 (3) episode window (undershoot_attribution reads it)
        rep["episode_window"] = f"{t_pre:g}-{early:g}"
    ov = {}
    for name, (a, b) in primary_windows(meta).items():
        w = (t > a - 1e-9) & (t <= b + 1e-9)
        ov[name] = dict(window=[a, b],
                        worst_K=float(dev[w].min()) if w.any() else None,
                        holds_worst=bool(a - 1e-9 < t[i] <= b + 1e-9),
                        n_below_20mK=int((w & below).sum()))
    rep["primary_window_overlap"] = ov
    return rep


def undershoot_attribution(detail_dt, detail_half, bound_K=None,
                           ratio=MP_ATTRIBUTION_RATIO):
    """D-F.10 (3): attribution of a rung's interstice undershoot by its
    class's half-Dt partner, like for like (the C1 12.5 / 6.25 ms pair at
    G-Dt). ``detail_dt`` / ``detail_half``: max_principle_detail of the rung
    at Dt and of its Dt/2 partner. Both are read over the same episode window
    (t_pre, t_pre + MP_EARLY_S] ("episode_window"; the 30 s ladder protocol
    of both rungs gives (30, 130] s) and in the four-cell interstice cores
    (minima_mK[...]["interstice"]); u = max(0, -worst) [mK].

    - The Dt rung's interstice undershoot within the bound (u_dt <= 1e3 x
      gates.UNDERSHOOT_K): verdict "within_bound", nothing to attribute.
    - Else "dt_controlled" iff u_half <= the bound and u_half <= ratio x u_dt
      (ratio 0.5): the class's Dt* must then satisfy the bound too (NEW-2.1
      (c), NEW-2.4); otherwise "h_transport": a property of the energy
      transport, which only the D-J (ii) scheme gate can clear.
    - The Dt/2 rung's interstice value is unobserved when no sample's
      minimum lies in an interstice core in that window (the monitor locates
      the global minimum only); it is then bounded by the window's global
      minimum (an interstice undershoot is never larger), and the verdict
      uses that bound ("basis_half": "all"). A Dt rung without an interstice
      minimum in the window, or reports over different windows, are input
      errors (ValueError)."""
    bound = 1e3 * float(G.UNDERSHOOT_K if bound_K is None else bound_K)

    def window(d, lab):
        w = d.get("episode_window")
        if w is None or "minima_mK" not in d:
            raise ValueError(f"undershoot attribution: the {lab} report has no "
                             "located minima (TfluidMinMax) or no episode window")
        return w
    w_dt, w_half = window(detail_dt, "Dt"), window(detail_half, "Dt/2")
    if w_dt != w_half:
        raise ValueError(f"undershoot attribution: episode windows differ "
                         f"({w_dt} against {w_half}); like for like needs one window")
    row_dt, row_half = detail_dt["minima_mK"][w_dt], detail_half["minima_mK"][w_half]
    if row_dt.get("interstice") is None:
        raise ValueError(f"undershoot attribution: the Dt rung has no interstice "
                         f"minimum in {w_dt} s")
    u_dt = max(0.0, -float(row_dt["interstice"]["worst_mK"]))
    basis = "interstice"
    half = row_half.get("interstice")
    if half is None:
        half, basis = row_half.get("all"), "all"
        if half is None:
            raise ValueError(f"undershoot attribution: the Dt/2 rung has no sample "
                             f"in {w_half} s")
    u_half = max(0.0, -float(half["worst_mK"]))
    out = dict(window_s=w_dt, bound_mK=bound, ratio=ratio, u_dt_mK=u_dt,
               u_half_mK=u_half, basis_half=basis,
               at_dt=row_dt["interstice"], at_half=half)
    if u_dt <= bound:
        out.update(verdict="within_bound",
                   note="the Dt rung's interstice undershoot is within the bound")
    elif u_half <= bound and u_half <= ratio * u_dt:
        out.update(verdict="dt_controlled",
                   note="Dt-controlled: the class's Dt* must satisfy the bound "
                        "(NEW-2.1 (c), NEW-2.4)")
    else:
        out.update(verdict="h_transport",
                   note="a property of the h transport: only the D-J (ii) scheme "
                        "gate can clear it")
    return out


# --------------------------------------------------------------------------
# D-G.1 (4): outlet backflow per window and the riser exchange-flow bursts
# --------------------------------------------------------------------------

def group_runs(idx, gap):
    """Split sorted indices into runs whose consecutive members are at most
    `gap` apart."""
    idx = np.asarray(idx, int)
    if idx.size == 0:
        return []
    return np.split(idx, np.flatnonzero(np.diff(idx) > gap) + 1)


def outlet_windows(meta):
    """D-G.1 (4): [0, t_pre], [t_pre, t_pre + 70 s], [t_pre + 70 s, t_end]
    and the heated stage [t_pre, t_end] (an isothermal or pre-circulation
    run: the whole run)."""
    t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
    if t_end <= t_pre + 1e-9:
        return [(0.0, t_end)]
    mid = min(t_end, t_pre + HEAT_EARLY_S)
    w = [(0.0, t_pre), (t_pre, mid), (mid, t_end), (t_pre, t_end)]
    return [x for x in w if x[1] > x[0] + 1e-9]


def outlet_flow(case, meta, windows=None, burst_u=BURST_U, burst_z=BURST_Z):
    """D-G.1 (4): the outlet backflow fraction (sumMag - sum)/2 over sumMag
    of the outlet mass flux per window, and the riser exchange-flow bursts
    (|U|max > burst_u at z > burst_z, TfluidMinMax), each with its time,
    peak |U|, peak fluid KE and that peak's ratio to the median KE of the
    developed heated window [t_pre + 270 s, t_end] (the heated stage if
    shorter than 30 s; the whole run if isothermal)."""
    windows = outlet_windows(meta) if windows is None else windows
    ts, s = HIO.column(case, "fluid", "sum_outlet")
    _, sm = HIO.column(case, "fluid", "sumMag_outlet")
    rep = dict(backflow={}, bursts=None, burst_u_mm_s=1e3 * burst_u,
               burst_z_mm=1e3 * burst_z)
    for a, b in windows:
        m = (ts >= a - 1e-9) & (ts <= b + 1e-9)
        den = _trapz(sm[m], ts[m]) if m.sum() > 1 else 0.0
        rep["backflow"][f"{a:g}-{b:g}"] = (0.5 * _trapz((sm - s)[m], ts[m]) / den
                                           if den > 0 else None)
    try:
        U = HIO.minmax(case, "fluid", "TfluidMinMax")["mag(U)"]
    except (FileNotFoundError, KeyError) as e:
        rep["burst_note"] = f"no TfluidMinMax mag(U) monitor ({e!s})"
        return rep
    tk, ke = HIO.column(case, "fluid", "KE_fluid")
    t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
    a = t_pre + KE_REF_DELAY_S
    if meta["heat"]["kind"] == "none" or t_end <= t_pre + 1e-9:
        ref_w = (float(meta["ramp_s"]), t_end)
    elif t_end - a >= 30.0:
        ref_w = (a, t_end)
    else:
        ref_w = (t_pre, t_end)
    mref = (tk > ref_w[0]) & (tk <= ref_w[1])
    ke_ref = float(np.median(ke[mref])) if mref.any() else None
    rep.update(ke_ref_median=ke_ref, ke_ref_window=list(ref_w), bursts=[])
    loc = np.asarray(U.max_loc, float)
    k = np.flatnonzero((U.max > burst_u) & (loc[:, 2] > burst_z))
    for g in group_runs(k, BURST_GAP_SAMPLES):
        j = g[np.argmax(U.max[g])]
        mk = (tk >= U.t[g[0]] - 1e-9) & (tk <= U.t[g[-1]] + 1e-9)
        kmax = float(ke[mk].max()) if mk.any() else None
        rep["bursts"].append(dict(
            t0=float(U.t[g[0]]), t1=float(U.t[g[-1]]),
            u_max_mm_s=1e3 * float(U.max[j]), t_umax=float(U.t[j]),
            location_mm=[round(1e3 * float(v), 2) for v in loc[j]],
            ke_max=kmax,
            ke_ratio=(kmax / ke_ref if kmax is not None and ke_ref else None)))
    return rep


def closure_gates(hy, level, tag=""):
    """D-P.7: the gates of one closure window by waveform class: steady
    isothermal -> hydraulic_closure (the level's tolerance, D-P.4); sine
    isothermal -> the necessary checks closure_sign and, with >= 2
    sub-windows, closure_stationarity; square and heated windows are
    reported only."""
    out = []
    if hy["gated"]:
        out.append((f"hydraulic_closure{tag}", G.hydraulic_closure_gate(
            hy["E_port"], hy["Phi"], hy["W_b"], hy["dKE"], level=level)))
    sc = hy.get("sine_checks")
    if sc and sc["gated"]:
        out.append((f"closure_sign{tag}", G.closure_phase_sign_gate(
            sc["sign"]["r_bins_per_period"], sc["phi_per_period"])))
        st = sc["stationarity"]
        if st["ok"] is not None:
            out.append((f"closure_stationarity{tag}", G.closure_stationarity_gate(
                [s_["r_over_phi"] for s_ in st["sub_windows"]],
                st["window_r_over_phi"])))
    return out


# --------------------------------------------------------------------------
# D-F.9: the Courant regions (thin layer T, outflow column, claim-region bulk)
# --------------------------------------------------------------------------

REGION_FIELDS =("CoThin", "CoOut", "CoBulk")     # D-F.9 (4): the gated regions


def _in_outflow(loc_m, oc=None):
    """D-F.9 (2) R_out membership of fieldMinMax locations (n, 3) [m]."""
    oc = HP.outflow_column() if oc is None else oc
    p = np.atleast_2d(np.asarray(loc_m, float))
    return (p[:, 0] >= oc["x_from"]) | (p[:, 2] > oc["z_above"])


def _in_head_envelope(loc_m):
    """Inside the envelope that holds T_head (layer B within pitch/sqrt 2 of
    a cell axis, housing_physics.THIN_HEAD): a point outside it is never in
    T_head, whatever the mesh; inside it, membership needs the cell width."""
    th = HP.THIN_HEAD
    p = np.atleast_2d(np.asarray(loc_m, float))
    ax = np.asarray(HP.cell_axes(), float)
    r = np.sqrt((p[:, 0, None] - ax[:, 0]) ** 2
                + (p[:, 1, None] - ax[:, 1]) ** 2).min(axis=1)
    return (p[:, 2] >= th["z"][0]) & (p[:, 2] <= th["z"][1]) & (r <= th["r"])


def _in_band(loc_m, tl=HP.THIN_LAYER):
    p = np.atleast_2d(np.asarray(loc_m, float))
    return (np.abs(p[:, 0]) <= tl["band_x_half"]) & (p[:, 2] >= tl["band_z"][0]) \
        & (p[:, 2] <= tl["band_z"][1])


def _peak(m, sel=None):
    """{max, t, location_m} of a fieldMinMax series over the samples sel."""
    idx = np.arange(m.t.size) if sel is None else np.flatnonzero(sel)
    if not idx.size:
        return None
    j = int(idx[np.argmax(m.max[idx])])
    return dict(max=float(m.max[j]), t=float(m.t[j]),
                location_m=[float(x) for x in m.max_loc[j]])


def legacy_regions(bulk, rb=1.0, thin=None):
    """D-F.9 regions bounded from a pre-D-F.9 CoMinMax_fluid (GDT1/GDT2 era:
    CoBulk = the whole fluid outside the D-F.1 thin layer, band + sleeves,
    with the location of each sample's maximum). Per sample, the region
    holding the bulk maximum carries exactly that value and every other
    bulk region at most that value, so each D-F.9 region's run maximum lies
    in [lower, upper]: lower = the largest sample maximum located in the
    region, upper = the largest bulk maximum over all samples; it is exact
    when the run's bulk maximum lies in the region. The D-F.9 claim-region
    bulk excludes T_head, whose cells the monitor does not identify: a
    sample located outside R_out, the band and the T_head envelope is
    claim-region bulk on any mesh; one inside the envelope is unresolved
    without the case's T_head set. Reported, never gated: the gate of a
    case without the D-F.9 monitors is the strict global one.

    ``thin`` (the legacy CoThin series, T = band + sleeves) adds thin_T, the
    D-F.9 thin layer T = band + T_head: D-F.9's T lies inside the legacy
    thin layer and bulk together, so its run maximum lies in [lower, upper],
    lower = the largest CoThin sample located in the cap-plane band (band
    cells are in T on any mesh), upper = the largest of the CoThin and
    CoBulk maxima; exact when the run's monitored maximum is a band sample
    (GDT2: C3 4.586, C1 5.105, both at z 65.02 mm). This lower value is the
    rung's co_thin in D-F.1 (3) (ladder_courant; conservative when not
    exact). The solver log's every-step maximum has no location and never
    enters it (GDT2 C1: 5.675 at field steps 2-4 of the period, which the
    5-step monitor does not sample)."""
    upper = float(bulk.max.max())
    out_m = _in_outflow(bulk.max_loc)
    env = _in_head_envelope(bulk.max_loc) & ~out_m
    claim_any = ~out_m & ~env & ~_in_band(bulk.max_loc)

    def span(sel, extra=None):
        pk = _peak(bulk, sel)
        d = dict(lower=None if pk is None else pk["max"], upper=upper,
                 exact=pk is not None and pk["max"] == upper, at_lower=pk)
        if extra:
            d.update(extra)
        return d
    out = dict(
        note="pre-D-F.9 monitors (CoThin, CoBulk with T = the D-F.1 band and "
             "sleeves): regions bounded from the location of each sample's "
             "bulk maximum; the gate is the strict global one",
        outflow_column=span(out_m),
        claim_bulk_sleeves=span(~out_m),
        claim_bulk=span(claim_any, dict(
            unresolved_samples_gt_limit=int(np.sum(env & (bulk.max > rb))),
            unresolved_max=(float(bulk.max[env].max()) if env.any() else None),
            unresolved_note="samples in the T_head envelope (layer B within "
                            "pitch/sqrt 2 of a cell axis): T_head or claim-"
                            "region bulk by cell width, which needs the "
                            "case's mesh (k0_measure.py --thin-head)")))
    if thin is not None:
        up_t = max(float(thin.max.max()), upper)
        pk = _peak(thin, _in_band(thin.max_loc))
        out["thin_T"] = dict(
            lower=None if pk is None else pk["max"], upper=up_t,
            exact=pk is not None and pk["max"] == up_t, at_lower=pk,
            note="D-F.9 thin layer (band + T_head) from the legacy CoThin "
                 "(band + sleeves): lower = the largest CoThin sample located "
                 "in the cap-plane band, the rung's co_thin in D-F.1 (3); the "
                 "unlocated every-step log maximum never enters it")
    return out


def ladder_courant(co):
    """D-F.1 (3), D-F.9 (2): {co_thin, co_out} of a medium ladder rung for
    recipe_rules.courant_bounds, from courant_location's report. With the
    D-F.9 monitors: the CoThin and CoOut maxima. With the legacy monitors
    (GDT2's 12.5 ms rungs): co_thin = legacy_regions' thin_T lower value
    (the located cap-plane band maximum; exact when it is the run's
    monitored maximum, as on both GDT2 rungs) and co_out = the outflow
    column's maximum where it is exact (the run's bulk maximum lies in
    R_out), else None (the rung does not enter B_out; recipe_rules.
    outflow_b). With the global Co monitor only (GDT1): the run's maximum
    is exact in the region that holds its location, so co_thin = that
    maximum when it lies in the cap-plane band (GDT1: 4.482, 2.229 at
    z 64.97-65.02 mm, the rungs of D-F.1's B 5.60) and co_out None (no
    R_out reading, D-F.9 (5)); a global maximum outside the band bounds the
    thin layer from above only, an input error. The solver log's every-step
    maximum is unlocated and never enters. Each value carries its exact
    flag. A report without a located monitor is an input error."""
    if co.get("regions") == "D-F.9":
        return dict(co_thin=float(co["CoThin"]["max"]), co_out=float(co["CoOut"]["max"]),
                    exact_thin=True, exact_out=True, source="D-F.9 monitors")
    lg = co.get("legacy_D-F.9")
    if not lg and co.get("location_m") is not None and "regions" not in co:
        loc = co["location_m"]
        if not bool(_in_band(loc)[0]):
            raise ValueError("ladder_courant: the global Co maximum is not located in "
                             "the cap-plane band; it bounds the thin layer from above only")
        # the band lies outside R_out: a global-only run has no R_out reading
        return dict(co_thin=float(co["Co_max_monitor"]), exact_thin=True,
                    co_out=None, exact_out=False,
                    source="global Co monitor (maximum located in the band)")
    if not lg or "thin_T" not in lg:
        raise ValueError("ladder_courant: no located Courant monitor (CoThin or the "
                         "global Co with location); a value without location bounds "
                         "no region")
    tt, oc = lg["thin_T"], lg["outflow_column"]
    if tt["lower"] is None:
        raise ValueError("ladder_courant: no CoThin sample located in the cap-plane "
                         "band; the D-F.9 thin-layer maximum is not bounded below")
    return dict(co_thin=float(tt["lower"]), exact_thin=bool(tt["exact"]),
                co_out=float(oc["lower"]) if oc["exact"] else None,
                exact_out=bool(oc["exact"]), source="legacy monitors (legacy_regions)")


def courant_location(case, meta=None):
    """Where the fluid Courant number peaks: the CoMinMax_fluid monitor
    (fieldMinMax of the CourantNo field, location true) and its regions.
    D-F.9 monitors (CoThin, CoOut, CoBulk, CoBulkSleeves): each region's
    maximum with time and location, the case's class bounds, and the
    claim-region bulk under D-F.1's sleeves as the D-F.9 (1) sensitivity
    column (reported against the bulk limit, never gated). Pre-D-F.9
    monitors (CoThin, CoBulk; GDT1/GDT2): their maxima and the D-F.9
    regions bounded by legacy_regions. The P00b hot layer: the z = 65 mm
    cap-plane cells."""
    try:
        mmd = HIO.minmax(case, "fluid", "CoMinMax_fluid")
        mm = mmd["Co"]
    except (FileNotFoundError, KeyError) as e:
        return dict(note=f"no CoMinMax_fluid monitor ({e!s}); cases rendered "
                         "before WP 2.5c completion lack it")
    k = int(np.argmax(mm.max))
    out = dict(Co_max_monitor=float(mm.max[k]), t=float(mm.t[k]),
               location_m=[float(x) for x in mm.max_loc[k]],
               samples=int(mm.t.size),
               Co_max_last=float(mm.max[-1]),
               location_last_m=[float(x) for x in mm.max_loc[-1]])
    for fld in HP.COURANT_FIELDS[1:]:
        if fld in mmd:
            out[fld] = _peak(mmd[fld])
    meta = {} if meta is None else meta
    try:
        rb = G.courant_level_ratios(meta.get("grid", "medium"))[0]
    except ValueError:
        rb = None
    if all(f in mmd for f in REGION_FIELDS):
        wf = meta.get("waveform")
        out["regions"] = "D-F.9"
        if wf in G.COURANT_THIN_B_CLASS:
            out["bounds_medium"] = dict(
                thin_B=G.courant_class_bounds(wf)[1],
                outflow_B_out=G.courant_out_bound(wf),
                claim_bulk=G.courant_class_bounds(wf)[0])
        if "CoBulkSleeves" in mmd and rb is not None:
            s = out["CoBulkSleeves"]
            lim = G.COURANT_LIMIT * rb
            out["CoBulkSleeves"] = dict(
                s, limit=lim, within=s["max"] <= lim,
                note="D-F.9 (1) sensitivity: the claim-region bulk with D-F.1's "
                     "registered sleeves (8 <= r <= 10 mm) in place of T_head; "
                     "reported, not gated")
    elif "CoBulk" in mmd:
        out["regions"] = "D-F.1 (legacy)"
        out["legacy_D-F.9"] = legacy_regions(mmd["CoBulk"],
                                             rb if rb is not None else 1.0,
                                             thin=mmd.get("CoThin"))
    return out


def courant_gate(case, meta, logs):
    """The Courant rule by region (D-F.1, per class D-F.5, regions D-F.9)
    where the case has the D-F.9 monitors (CoThin, CoOut, CoBulk): the
    bounds of the case's waveform class (B, B_out, claim-region bulk 1),
    level-scaled, and the log's every-step maximum <= max(B, B_out, 1).
    Without them (cases rendered before D-F.9, e.g. GDT1/GDT2, whose CoBulk
    held R_out and the sleeve edges together; or no monitor at all) the
    strict global gate on every log line (Co <= 1 everywhere): the regions
    such a case can only bound are reported by courant_location."""
    co_log = np.concatenate([lg.courant_max for lg in logs])
    diag = bool(meta["diagnostic"])
    try:
        mmd = HIO.minmax(case, "fluid", "CoMinMax_fluid")
    except FileNotFoundError:
        mmd = {}
    if all(f in mmd for f in REGION_FIELDS):
        wf = meta.get("waveform", "steady")
        _, b_thin = G.courant_class_bounds(wf)
        return G.thin_layer_courant_gate(
            mmd["CoBulk"].max, mmd["CoThin"].max, mmd["CoOut"].max, co_log,
            level=meta.get("grid", "medium"), thin_limit=b_thin,
            out_limit=G.courant_out_bound(wf), diagnostic=diag)
    r = G.courant_gate(co_log, diagnostic=diag)
    why = ("pre-D-F.9 region monitors (no CoOut)" if "CoBulk" in mmd
           else "no Courant region monitors")
    return G.GateResult(r.ok, f"{r.detail}; {why}: the strict global gate "
                              "(D-F.9)")


def evaluate(case):
    case = Path(case)
    meta = HIO.case_meta(case)
    lat = meta["lattice"]
    dt = lat["dt"]
    t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
    res, info = [], {}
    info["claim_eligible"], info["claim_notes"] = HIO.claim_status(meta)
    logs, runs = solver_steps(case)
    # the steps that made the stitched solution: a crash-resume's replayed
    # steps count once (ledger E_inj, completeness, continuity)
    eff = effective_steps(runs, dt)
    all_steps = np.concatenate([r[1] for r in eff])
    n_replayed = sum(r[1].size for r in runs) - all_steps.size
    if n_replayed:
        info["replayed_steps"] = (f"{n_replayed} step(s) of a run superseded "
                                  "by a resume from an earlier written time "
                                  "counted once")
    if abs(all_steps[-1] - t_end) > 0.5 * dt:
        info["incomplete"] = f"last step {all_steps[-1]:g} s, t_end {t_end:g} s"

    # ---- lattice -----------------------------------------------------------
    P = meta.get("period_s")
    on = meta["duty"] * P if meta["waveform"] == "square" else None
    edge = lat.get("valve_edge_s")
    res.append(("lattice", G.lattice_gate(
        dt, period=P, on_time=on, sample=lat["sampler_steps"] * dt,
        valve_edge=edge, n_min=meta["n_min"],
        n_phase_min=meta["n_phi_min"] if P else 8,
        diagnostic=meta["diagnostic"],
        edge_sensitivity=(meta["waveform"] == "square" and meta["edge"] != "dt"))))

    # ---- fixed dt, Courant --------------------------------------------------
    # every step every run took (replayed steps included: each was a solver
    # step and must be at the fixed dt)
    for k, (start, steps, lg) in enumerate(runs):
        res.append((f"fixed_dt[{Path(lg.path).name.removesuffix('.zst')}#{k}]",
                     G.fixed_dt_gate(steps, dt, t_start=start)))
    res.append(("courant", courant_gate(case, meta, logs)))
    info["courant"] = courant_location(case, meta)

    # ---- max principle ------------------------------------------------------
    tf, fmin = _col(case, "fluid", "Tmin_fluid")
    _, fmax = _col(case, "fluid", "Tmax_fluid")
    smax = None
    for i in range(N_CELLS):
        ts, v = _col(case, f"cell{i}", f"Tmax_cell{i}")
        if ts.size != tf.size or np.any(np.abs(ts - tf) > 1e-9):
            raise ValueError("cell and fluid monitors not on one time base")
        smax = v if smax is None else np.maximum(smax, v)
    # judged on the heated samples (t > t_pre): before heat-on the solids sit
    # at T_in and the fluid departs from it only by the energy equation's
    # pressure-work and potential-energy terms (dpdt, rho U.g: micro- to
    # milli-kelvin, both signs), which the D-F bound is not about
    heated = tf > t_pre + 0.5 * dt
    if meta["heat"]["kind"] != "none" and heated.sum():
        # D-F.10: the bound unchanged for every class and level; the detail
        # names the worst sample's time, place, phase and span
        res.append(("max_principle", G.max_principle_gate(
            fmin[heated], fmax[heated], meta["T_inlet"], smax[heated],
            times=tf[heated])))
        info["max_principle_detail"] = max_principle_detail(
            case, meta, tf, fmin, heated)
    else:
        info["max_principle"] = ("not evaluated: no heated samples "
                                 "(isothermal run or pre-circulation only)")
    info["T_fluid_range_K"] = [float(fmin.min()), float(fmax.max())]
    info["T_fluid_range_pre_K"] = [float(fmin[~heated].min()),
                                   float(fmax[~heated].max())] \
        if (~heated).sum() else None
    info["T_solid_max_K"] = float(smax.max())

    # ---- mean flow ----------------------------------------------------------
    # Qbar over whole periods of [max(ramp, t_pre), t_end]; a run whose
    # heated stage ends before one whole period after the start-up ramp
    # (smokes) is judged as DELIVERY AGAINST THE RAMPED INTENT instead:
    # Qbar + mean(Q - Q_intent(t)) over the whole periods of the run, with
    # Q_intent the same waveform x ramp the inlet Function1 evaluates
    tq, vq = _col(case, "fluid", "volFlux_inlet")
    q = -vq
    # a pre-circulation-only run (t_end = t_pre: a shared NEW-14 state, the
    # P00b valve check) is judged after its ramp
    t_a = (max(meta["ramp_s"], t_pre) if t_end > t_pre + 0.5 * dt
           else meta["ramp_s"])
    if (P and t_end - t_a >= P - 1e-9) or (not P and t_end - t_a > 1e-9):
        m = tq >= t_a - 1e-9
        res.append(("mean_flow", G.mean_flow_gate(
            tq[m], q[m], meta["qbar_m3s"], period=P, t_start=t_a)))
    else:
        qf, _, _, _ = PE.waveform_from_meta(meta)
        intent = qf(tq) * np.clip(tq / meta["ramp_s"], 0.0, 1.0)
        ok, det = G.mean_flow_gate(tq, q - intent + meta["qbar_m3s"],
                                   meta["qbar_m3s"], period=P)
        res.append(("mean_flow", G.GateResult(
            ok, "delivery against the ramped intent (run shorter than a "
                "period after the ramp): " + det)))

    # ---- ledger (heated stage) ----------------------------------------------
    if meta["heat"]["kind"] == "none":
        info["ledger"] = "not evaluated: isothermal run (no injected energy)"
    else:
        ok_led, led = ledger(case, meta, eff, all_steps)
        res.append(("ledger", ok_led))
        info["ledger"] = led

    # ---- source power (D-J (i)) ----------------------------------------------
    sel_runs = selected_volumes(logs)
    sel = sel_runs[0] if sel_runs else []
    if meta["heat"]["kind"] != "none":
        if len(sel) != N_CELLS:
            res.append(("source_power", G.InputError(
                False, f"source power: {len(sel)} heat selections in the "
                       f"first solver run, {N_CELLS} expected")))
        elif any(r != sel for r in sel_runs[1:]):
            res.append(("source_power", G.GateResult(
                False, f"source power: the heat selections differ between "
                       f"solver runs ({len(sel_runs)} runs): {sel_runs}")))
        else:
            vmesh = []
            for i in range(N_CELLS):
                segs = HIO.segments(case, f"cell{i}", f"Tmax_cell{i}")
                vmesh.append(float(segs[0][1].meta["Volume"]))
            # regions are constructed in regionProperties order cell0..5
            p_model = 1.0
            p_inj = [p_model * vm / vd for vm, (_, vd) in zip(vmesh, sel)]
            worst = max(p_inj, key=lambda x: abs(x - 1))
            ok, det = G.source_power_gate(worst, p_model)
            res.append(("source_power", G.GateResult(
                ok, det + f"; VDash {[v for _, v in sel]} m3, meshed "
                          f"{vmesh} m3 (volumeMode absolute)")))

    # ---- hydraulic work (D-P validation column; per class, D-P.7) ---------------
    level = meta.get("grid", "medium")
    hy = PE.measured(case)
    res.extend(closure_gates(hy, level))
    hy["note"] = ("E_port is the D-P validation column; the primary energy "
                  "basis is E_diss, the dissipative work of the isothermal "
                  "impedance model (D-P.1, D-P.3: scripts/post/pump_energy.py "
                  "primary). A negative E_port is chimney-driven flow, never "
                  "a pump saving")
    info["hydraulic"] = hy
    # the port-work instrument of a HEATED case, checked on its isothermal
    # pre-circulation [ramp, t_pre] (review SHOULD-FIX 2; the closure gated on
    # steady, the necessary checks on sine, D-P.7; an INIT_FROM member's
    # pre-circulation is its shared state's)
    if meta["heat"]["kind"] != "none" and t_pre > 0:
        hp = PE.measured(case, stage="pre")
        res.extend(closure_gates(hp, level, "[pre]"))
        info["hydraulic_pre"] = hp

    # ---- continuity ------------------------------------------------------------
    # A heat-on after a pre-circulation run (Allrun.heat, INIT_FROM) is a NEW
    # run from a written state (NEW-25 INIT_FROM semantics, never a resume):
    # the state may carry the tilt seed, so each stage is judged on its own
    # segments against its own logs.
    two_stage = meta["run_stage"] == "heat"
    stages = ([("pre", lambda x: x < t_pre - 0.5 * dt),
               ("heat", lambda x: x >= t_pre - 0.5 * dt)] if two_stage
              else [("run", lambda x: True)])
    for region, name in CONTINUITY_MONITORS:
        try:
            segs_all = [(lab, fo.t, np.asarray(fo.data[fo.columns[1]], float))
                        for lab, fo in HIO.segments(case, region, name)]
        except (FileNotFoundError, KeyError) as e:
            res.append((f"continuity[{name}]", G.InputError(False, repr(e))))
            continue
        for tag, inside in stages:
            segs = [sg for sg in segs_all if inside(sg[0])]
            # the runs in time order after the resume cut: the series ends
            # where the last effective run ends (a crash-resume's replayed
            # steps are not a later end)
            st_runs = [r for r in eff if inside(r[0])]
            if not segs and not st_runs:
                continue
            if not segs or not st_runs:
                res.append((f"continuity[{name}:{tag}]", G.InputError(
                    False, f"stage {tag}: {len(segs)} segment(s), "
                           f"{len(st_runs)} run(s)")))
                continue
            res.append((f"continuity[{name}:{tag}]", G.continuity_gate(
                segs, log_start_times=[r[0] for r in st_runs],
                log_end_time=float(st_runs[-1][1][-1]), write_precision=17)))

    # ---- reported quantities ------------------------------------------------------
    tsm = [_col(case, f"cell{i}", f"Tsurfmax_cell{i}") for i in range(N_CELLS)]
    info["T_surface_max_K"] = float(max(v.max() for _, v in tsm))
    vols = [float(HIO.segments(case, f"cell{i}", f"TvolAverage_cell{i}")[0][1]
                  .meta["Volume"]) for i in range(N_CELLS)]
    tav = [_col(case, f"cell{i}", f"TvolAverage_cell{i}")[1]
           for i in range(N_CELLS)]
    mod = sum(v * a for v, a in zip(vols, tav)) / sum(vols)
    info["T_module_mean_end_K"] = float(mod[-1])
    info["dT_cell_end_K"] = float(max(a[-1] for a in tav) -
                                  min(a[-1] for a in tav))
    ts_, s_ = _col(case, "fluid", "sum_outlet")
    _, sm_ = _col(case, "fluid", "sumMag_outlet")
    back = 0.5 * (sm_ - s_)
    info["outlet_backflow_fraction"] = (_trapz(back, ts_) /
                                        max(_trapz(sm_, ts_), 1e-300))
    # D-G.1 (4): backflow per window and the riser bursts, beside the gates
    info["outlet"] = outlet_flow(case, meta)
    pc = case / "precirc_check.json"
    if pc.exists():
        info["precirc_check"] = json.loads(pc.read_text())
    gate, rep = precirc_gate(case, meta, info["claim_eligible"])
    if gate is not None:
        res.append(("precirc", gate))
    info["precirc"] = rep
    return res, info


def precirc_gate(case, meta, claim_eligible):
    """NEW-14.1 (c) / NEW-14.2: (GateResult or None, report). The rule of
    precirc_check.check on [0, t_pre] (the case's own monitors, or the
    shared state's through INIT_FROM), recomputed by this registered scorer.
    A gate for a claim-eligible run (input problems are an input error);
    otherwise reported only. No pre-circulation (zero flow, t_pre 0): not
    applicable (an impulsive start of a flowing case is not claim-eligible,
    housing_physics.Config.claim_eligible)."""
    t_pre = float(meta["t_pre"])
    if t_pre <= 0 or float(meta["qbar_m3s"]) == 0.0:
        return None, dict(applicable=False, note="no pre-circulation (zero flow "
                                                 "or t_pre = 0)")
    try:
        r = PC.check(case)
    except (OSError, KeyError, ValueError, StopIteration) as e:
        if claim_eligible:
            return G.InputError(False, f"pre-circulation check: {e!r}"), dict(
                applicable=True, error=repr(e))
        return None, dict(applicable=True, error=repr(e), gated=False)
    det = (f"pre-circulation t_pre {t_pre:g} s, windows {r['window_s']:g} s: last "
           f"dp* changes {[round(100 * c, 3) for c in r['last_changes_dp']]} %, KE "
           f"{[round(100 * c, 3) for c in r['last_changes_ke']]} % (tol "
           f"{100 * PC.REL_TOL:g} %); "
           + ("stationary" if r["stationary"] else "NOT stationary: "
              + "; ".join(r["notes"])))
    rep = dict(applicable=True, gated=bool(claim_eligible), **r)
    if claim_eligible:
        return G.GateResult(bool(r["stationary"]), det), rep
    return None, rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case", type=Path)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--steady-ref", type=Path,
                    help="D-P.7 (6): a steady case of the same protocol; the "
                         "per-period closure excess of this (pulsed) case "
                         "over it is reported for each window")
    a = ap.parse_args(argv)
    try:
        res, info = evaluate(a.case)
        if a.steady_ref is not None:
            for key, stage in (("hydraulic", "run"), ("hydraulic_pre", "pre")):
                hy = info.get(key)
                if hy is None or not hy.get("period_s"):
                    continue
                ref = PE.measured(a.steady_ref, stage)
                hy["excess_over_steady"] = dict(
                    PE.excess_over_steady(hy, ref), steady_case=str(a.steady_ref),
                    steady_window=ref["window"])
    except (OSError, KeyError, ValueError, StopIteration) as e:
        print(f"housing_metrics: input error: {e!r}", file=sys.stderr)
        return 2
    code = G._report(res)
    face_failed = False
    led = info["ledger"]
    if isinstance(led, dict):
        print(f"ledger ({led['form']}): E_inj {led['E_inj']:.6g} J, residual "
              f"{led['residual']:+.4g} J ({100 * led['residual_rel']:+.3f} %); "
              + ", ".join(f"{k} {v:+.4g}" for k, v in led["terms"].items()))
        print(f"ledger, naive form (Delta H, trapezoids; not the gate): "
              f"residual {led['naive_residual']:+.4g} J "
              f"({100 * led['naive_residual_rel']:+.3f} %)")
    else:
        print(f"ledger: {led}")
    for key in ("hydraulic", "hydraulic_pre"):
        hy = info.get(key)
        if hy is None:
            continue
        edge = (f", edge-loss estimate {hy['E_edge_model']:.4g} J (residual "
                f"less it {100 * hy['closure_less_edge_model_rel_phi']:+.1f} "
                f"% of Phi)" if math.isfinite(hy["E_edge_model"]) else "")
        print(f"{key} [{hy['window'][0]:g}, {hy['window'][1]:g}] s "
              f"({hy['window_note']}; "
              f"{'GATED' if hy['gated'] else hy['gate_note']}): E_port "
              f"{hy['E_port']:+.4g} J, Phi {hy['Phi']:.4g} J, W_b "
              f"{hy['W_b']:+.4g} J, dKE {hy['dKE']:+.4g} J, KE exchanged "
              f"{hy['KE_exchanged']:.4g} J, residual "
              f"{hy['closure_residual']:+.4g} J "
              f"({100 * hy['closure_rel']:+.2f} % of max(|Phi|,|W_b|))"
              + edge)
        pf = hy.get("phi_face")
        if pf and pf.get("error"):
            face_failed = True
            print(f"  face-based Phi (D-P.9, reported, never gated): HARD FAILURE of the "
                  f"instrument, no face-based number reported: {pf['error']}")
        elif pf:
            print(f"  face-based Phi (D-P.9, reported, never gated): Phi_face "
                  f"{pf['Phi_face']:.4g} J ({pf['face_over_cell']:.4f} x the cell-centred "
                  f"Phi; wall half-cells {100 * pf['wall_share']:.1f} %, port ducts "
                  f"{100 * pf['duct_share_face']:.1f} % against {100 * pf['duct_share_cell']:.1f} "
                  f"% cell-centred), residual {pf['closure_residual_face']:+.4g} J: r/Phi face "
                  f"{100 * pf['r_over_phi_face']:+.2f} % against cell-centred "
                  f"{100 * pf['r_over_phi_cell']:+.2f} %; W_bound (viscous work at the ports) "
                  f"{pf['W_bound']:+.4g} J ({100 * pf['W_bound_over_phi_face']:+.2f} % of "
                  f"Phi_face); identity {pf['identity_max_rel']:.1e}, cell column vs Phi_fluid "
                  f"{pf['cell_vs_registered_rel']:+.1e}")
        sc = hy.get("sine_checks")
        if sc:
            st = sc["stationarity"]
            print(f"  sine (D-P.7 (3), {'gated' if sc['gated'] else 'reported'}): r/Phi "
                  f"{100 * sc['r_over_phi']:+.2f} %, {sc['r_per_period']:.4g} J per "
                  f"period; sign: {sc['sign']['detail']}; stationarity: {st['detail']}")
        sq = hy.get("square_split")
        if sq and "segments" in sq:
            print(f"  square by valve segment (D-P.7 (5), per period): edges r "
                  f"{sq['edges_r']:.4g} J ({100 * sq['edges_share_of_r']:.1f} % of r), "
                  f"less plateau-like {sq['edges_less_plateau_like']:.4g} J against "
                  f"KE_exchanged/n^2 {sq['edge_model_per_period']:.4g} J "
                  f"({100 * sq['edges_vs_model_rel']:+.1f} %); on r/Phi "
                  f"{100 * sq['on_r_over_phi']:.1f} %, off "
                  f"{100 * sq['off_r_over_phi']:.1f} %")
        ex = hy.get("excess_over_steady")
        if ex:
            print(f"  excess over steady (D-P.7 (6), per period): "
                  f"{ex['excess_per_time']:.4g} J (steady per unit time), "
                  f"{ex['excess_per_phi']:.4g} J (in proportion to Phi)")
    mp = info.get("max_principle_detail")
    if isinstance(mp, dict) and "worst_K" in mp:
        where = (f" at ({', '.join(f'{v:.2f}' for v in mp['location_mm'])}) mm"
                 f"{' (interstice)' if mp.get('in_interstice') else ''}"
                 if "location_mm" in mp else "")
        ph = f", phase {mp['phase']:.3f}" if "phase" in mp else ""
        b = mp["below_20mK"]
        span = (f"; below -20 mK: {b['n']} samples, t {b['t_first']:g}...{b['t_last']:g} s"
                if b["n"] else "; none below -20 mK")
        print(f"max principle detail (D-F.10): worst {1e3 * mp['worst_K']:+.2f} mK at t "
              f"{mp['t']:g} s ({mp['t_after_heat_on']:g} s after heat-on){ph}{where}{span}")
    ou = info.get("outlet")
    if ou:
        bf = ", ".join(f"[{k}] s {100 * v:.1f} %" if v is not None else f"[{k}] s n/a"
                       for k, v in ou["backflow"].items())
        nb = (f"; {len(ou['bursts'])} riser burst(s) (|U|max > "
              f"{ou['burst_u_mm_s']:g} mm/s at z > {ou['burst_z_mm']:g} mm)"
              + (", largest " + f"{max(b_['u_max_mm_s'] for b_ in ou['bursts']):.0f} mm/s"
                 if ou["bursts"] else "")
              if ou.get("bursts") is not None else f"; {ou.get('burst_note', '')}")
        print(f"outlet backflow (D-G.1 (4)): {bf}{nb}")
    if info["claim_eligible"]:
        print("claim-eligible: yes")
    else:
        print("NOT claim-eligible (claims.py refuses rows citing this run): "
              + "; ".join(info["claim_notes"]))
    co = info.get("courant", {})
    if "location_m" in co:
        x, y, z = co["location_m"]
        print(f"Courant monitor: max {co['Co_max_monitor']:.4g} at t "
              f"{co['t']:g} s, (x, y, z) = ({1e3 * x:.2f}, {1e3 * y:.2f}, "
              f"{1e3 * z:.2f}) mm")
    if co.get("regions") == "D-F.9":
        print("Courant regions (D-F.9): "
              + ", ".join(f"{f} {co[f]['max']:.4g}" for f in HP.COURANT_FIELDS[1:]
                          if co.get(f))
              + " (CoBulkSleeves: sensitivity, not gated)")
    elif "legacy_D-F.9" in co:
        lg = co["legacy_D-F.9"]

        def rng(d):
            lo = "n/a" if d["lower"] is None else f"{d['lower']:.4g}"
            return (f"{lo}" if d["exact"] else f"[{lo}, {d['upper']:.4g}]")
        tt = (f"thin layer T (band + T_head) {rng(lg['thin_T'])}, "
              if "thin_T" in lg else "")
        print("Courant regions from pre-D-F.9 monitors (bounds; the gate is "
              f"strict): {tt}outflow column {rng(lg['outflow_column'])}, claim-region "
              f"bulk under the D-F.1 sleeves {rng(lg['claim_bulk_sleeves'])}, "
              f"claim-region bulk {rng(lg['claim_bulk'])} with "
              f"{lg['claim_bulk']['unresolved_samples_gt_limit']} sample(s) "
              "above the limit unresolved (T_head envelope)")
    if a.json:
        a.json.write_text(json.dumps(dict(
            case=str(a.case), claim_eligible=info["claim_eligible"],
            claim_notes=info["claim_notes"], gates=[dict(gate=n, ok=bool(r.ok),
                                          input_error=isinstance(r, G.InputError),
                                          detail=r.detail) for n, r in res],
            info=info, policy=G.POLICY),
            indent=1, default=float) + "\n")
    return max(code, 2) if face_failed else code


if __name__ == "__main__":
    sys.exit(main())
