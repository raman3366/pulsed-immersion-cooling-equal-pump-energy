#!/usr/bin/env python3
"""Pump energy of the housing (plan D-P as decided for WP 2.5c, amended by
PROJECT.md §2.5 D-P.1 and D-P.3 on 2026-09-28; D-H, D-H.1).

PRIMARY basis (D-P.1) - the DISSIPATIVE hydraulic work of an ISOTHERMAL
IMPEDANCE MODEL, E_primary = E_diss = int Q dp dt over whole periods (the
in-phase work: R0, K Q|Q| and the frequency-dependent R(w); the reactive
part cancels). It sets the F-E1 comparators (Q_eq, D-H.1), the fronts and
PEC. The non-regenerative drive work E_nr = int max(0, Q dp) dt (the key
E_pump below, which adds the inertial work of the MODELLED DOMAIN, valve
kinetic energy included) is reported as a column labelled domain-inertance
dependent: most of the housing's inertance is its inlet channel and riser,
whose length the D-G port gate sets for boundary independence, and with it
E_nr/E_diss is 10-42 and Q_eq 315-1134 mL/min for the four arms (D-P.1
evidence), so it cannot set a comparator. The system basis (a documented
reference loop) is exploratory (NEW-13.2, D-P.3) and not implemented here.

D-P.5 (2026-09-29): the in-phase excess Re dZ is clamped at >= 0 (a fitted
R(f) below the quasi-steady r_qs would make a pulsed arm cheaper than its
<Q^2> R0 work, which the energy-honesty rule excludes); Impedance.below_qs
reports the clamped fit frequencies.

D-P.3 (2026-09-28) - the registered use of the model:
  * outside the fitted frequencies [f_min, f_max] the in-phase excess
    Re dZ(f) = R(f) - r_qs is evaluated twice: held (clamped) and scaled as
    sqrt(f / f_max) above / sqrt(f / f_min) below (the Stokes-layer growth
    of the wall resistance, and its quasi-steady limit); the arm's E_diss
    is the LARGER of the two (primary_energy_conservative): a pulsed arm is
    never made cheaper by the extrapolation. The dissipative share of the
    harmonics outside the fitted range is reported;
  * u(E_diss) propagates the fit covariance (R0, K from the steady fit
    residuals; R(f_i) from the period-to-period scatter of the fundamental
    fit) through E_diss by central differences (ediss_uncertainty); the
    energy axis of the H_E front test is resolved beyond
    max(2 u(E_S - E_P), 2 % of the larger energy) (he_energy_u95);
  * the isothermal square-wave runs of the impedance block (C1, C2') check
    the model: it passes iff the model's conservative E_diss is >= 98 % of
    the measured isothermal port work less the valve-edge temporal
    dissipation estimate over the same whole periods (model_check); on a
    failure the arm's E_diss is that measured value (the larger) and Q_eq is
    re-derived from it.

The model (D-P as decided 2026-09-25):
The housing's measured port work is dominated by the chimney: with the
outlet riser warmer than the inlet, the buoyancy work W_b exceeds the
viscous loss by ~600x at Qbar 90 mL/min (K0-design §10.10 item 9), so the
measured module work E_port = Phi - W_b + dKE is negative and says nothing
about what a pump spends. The pump energy is therefore modelled on the
housing's isothermal impedance, identified from isothermal (no-heat) runs
of the same housing; its in-phase (dissipative) work is the primary basis
(D-P.1), the non-regenerative drive work E_nr a reported column:

    dp(t) = R0 Q + K Q|Q| + sum_k Re[ dZ(k w) Qhat_k e^{i k w t} ]
    dZ(w) = Z(w) - (R0 + 2 K Qbar_fit)   (dynamic excess; -> 0 as w -> 0)
    Z(w)  = R(w) + i w L(w)              (fitted per frequency)
    E_pump = int max(0, Q dp) dt         (E_nr, non-regenerative drive;
                                          always >= 0; independent of
                                          buoyancy by construction; reported)
    E_visc     = int Q dp dt             (in-phase work; the reactive part
                                          integrates to zero over periods)
                                          = E_diss = E_primary (D-P.1)
    E_inertial = E_pump - E_visc >= 0    (inertial work a non-regenerative
                                          drive does not recover; for a
                                          square wave ~ the valve kinetic
                                          energy n_open L Q_on^2 / 2, reported
                                          as E_valve_KE for comparison)
R0, K come from steady isothermal runs (dp = P_port / Q), R(w), L(w) from
the fundamental of oscillatory isothermal runs (sine, A 0.6, 0.5/1/2 Hz;
generate_cases.py --family housing --blocks impedance). Between fitted
frequencies R and L are interpolated in log f; beyond them L is held and
the in-phase excess is held (extrap "clamp") or scaled with sqrt(f)
(extrap "sqrt", D-P.3); the conservative E_diss is the larger of the two.
Reported: the |dZ|-weighted share of clamped harmonics
(extrapolated_fraction, which includes the reactive part) and the share of
the dynamic DISSIPATIVE work carried by harmonics outside the fitted range
(dissipative_extrapolated_share, the one that bounds the E_diss error).

VALIDATION column - measured port work (every-step monitors):
    E_port = int -sum_ports sum_f phi_f (p*_f / rho_f + |U_f|^2 / 2) dt,
    p* = p - rho(T_in) gh (the hydrostatic column at the inlet density,
    identical at both ports), reported with the buoyancy (chimney) work
    W_b = int int (rho - rho_in) U.g dV dt as its own term, the viscous
    dissipation Phi and dKE (closure Phi - W_b + dKE, gates.py
    hydraulic_closure_gate), and the forced part E_port + W_b. A negative
    E_port is chimney-driven flow and is never reported as a pump saving.
    The closure (re-specified in WP 2.5c on the P00b data, amended by its
    review) runs over a developed window, whole periods (closure_window),
    on max(|Phi|, |W_b|), and is a GATE only where it checks the port-work
    instrument at the scale of Phi: ISOTHERMAL windows of STEADY flow - the
    isothermal runs' measurement window and the pre-circulation [ramp,
    t_pre] of every heated run (stage "pre"). There the residual is the
    cell-centred Phi's discretisation difference (+7 % on the medium
    housing at 12.5 ms; one tolerance per grid level, provisional,
    gates.HYDRAULIC_TOL_LEVEL, D-P.4).
    D-P.7 (2026-10-02, after GDT2): per waveform class. SINE isothermal
    windows report the closure (GDT2 C3: +25.7 % of Phi on medium, above
    D-P.2's 0.15 ceiling when pooled with steady) and gate its necessary
    checks (phase-bin sign, sub-window stationarity; dt invariance at G-dt
    as a pair check, sine_dt_invariance), which cannot bound a stationary
    spatial bias in E_port: the face-based-Phi re-examination of D-P.2 /
    D-P.4 runs before Stage B, and no Stage-B fit sets a Q_eq before it
    reports (q_eq_registered). SQUARE windows are split by valve segment
    (square_closure_split).
    D-P.8 (2026-10-02): the steady fit is sign-guarded (fit_steady_guarded:
    R0 > 0; K < 0 refused beyond 2 SE, clamped to 0 within 2 SE) and an
    arm's E_diss is checked against steady flow at its mean (jensen_floor).
    REPORTED, never gated: heated windows (the residual measures the
    discretisation of the buoyant flow: numerical dissipation in the
    plumes, face-based buoyancy of the p_rgh formulation against the
    cell-centred W_b, not separated) and square waves (the backward
    scheme's dissipation of the kinetic energy the valve edges exchange
    dominates the residual - 2.1 Phi on the P00b valve case - so the check
    would not resolve the instrument at the scale of Phi; the edge-loss
    estimate KE_exchanged / n_edge^2 is reported beside it).

Energy honesty (study rule): at matched mean flow a pulsed arm always costs
at least the steady flow's viscous work (<Q^2> >= <Q>^2); the primary
basis adds inertial work on top.

Usage:
    pump_energy.py fit --steady CASE [CASE ...] --osc CASE [...] --out imp.json
    pump_energy.py primary --impedance imp.json --case CASE
    pump_energy.py primary --impedance imp.json --waveform sine --qbar 90
                   --freq 2 --amp 0.6 [--duty D] [--edge-ms 50] [--t 900]
    pump_energy.py qeq --impedance imp.json --energy J [--t 900]
                   (J = the arm's E_primary: D-H.1 matches on E_diss)
    pump_energy.py measured CASE [--stage pre]
    pump_energy.py sine-dt CASE_DT CASE_DT_HALF [--stage pre|run]
                   (D-P.7 (3c): exit 1 when the dt invariance fails)
    pump_energy.py closurediag CASE [CASE ...] [--stage run|pre] [--json F]
                   (D-P.9: the closure against the cell-centred and the
                   face-based Phi, and a sine case's excess over a steady
                   one in both forms; reported, never gated)

D-P.9 (2026-10-06): the face-based Phi (housing_physics.phi_face_monitor,
PHI_FACE=1) is the dissipation the discrete viscous operator does on the
current U, sum_f t_f . dU_f over every face (near-wall half-cells
included); measured() reports the closure against it beside the registered
cell-centred one (``phi_face``, face_closure), never gated, with W_bound,
the viscous traction's work at the ports, beside r_face. A face-based
monitor that fails its own checks is reported as such (``phi_face.error``)
and never hides the registered closure (internal review, 6 Oct).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "model"))
import housing_io as HIO  # noqa: E402

N_SYNTH = 4096          # samples per period for the harmonic synthesis
ML_MIN = 1.0 / 6e7      # m3/s per mL/min


# --------------------------------------------------------------------------
# impedance model
# --------------------------------------------------------------------------

EXTRAPOLATIONS = ("clamp", "sqrt")      # D-P.3: the two evaluations outside the fitted range
ENERGY_FLOOR_REL = 0.02                 # D-P.3 / NEW-13.2: the energy axis is never resolved
                                        # finer than the D-H 2 % matching tolerance
ENERGY_K = 2.0                          # coverage factor of U95(E_diss) (Type B, k = 2)
MODEL_CHECK_MIN_RATIO = 0.98            # D-P.3: model E_diss >= 98 % of the measured edged work
FD_REL_STEP = 1e-4                      # central-difference step of the E_diss Jacobian
FIT_K_NEG_SE = 2.0                      # D-P.8: K < 0 beyond 2 SE refused, within 2 SE clamped to 0
FIT_K_ROUNDOFF_REL = 1e-9               # D-P.8: |K| q_max <= 1e-9 R0 is round-off of a linear law
JENSEN_RTOL = 1e-9                      # D-P.8: E_diss(arm) >= E_diss(steady at its Qbar) to round-off
# D-P.9: the face-based Phi monitor's own checks (hard failures of the reader:
# a face-based number that fails them is never reported). The identity is the
# same dissipation from the solver's fvc operators (round-off: ~1e-14 on the
# Pi validation); the cell-centred column repeats Phi_fluid's integrand with
# thermo:mu (~1e-15)
PHI_FACE_IDENTITY_RTOL = 1e-8
PHI_FACE_CELL_RTOL = 1e-6
PHI_FACE_KEYS = ("Phi_face", "Phi_face_lap", "Phi_face_T", "Phi_face_wall",
                 "Phi_face_duct", "Phi_cell", "Phi_cell_lap", "Phi_cell_duct",
                 "Phi_identity", "W_bound")
# what face_closure's data and checks can raise; measured() turns it into
# phi_face.error so that the registered closure is still reported
PHI_FACE_ERRORS = (OSError, KeyError, ValueError, IndexError)


@dataclass
class Impedance:
    R0: float                       # Pa s/m3
    K: float = 0.0                  # Pa s2/m6
    qbar_fit: float = 0.0           # m3/s, mean flow of the oscillatory fits
    freqs: list = field(default_factory=list)   # Hz
    R: list = field(default_factory=list)       # Pa s/m3 (in-phase, at f)
    L: list = field(default_factory=list)       # Pa s2/m3 (inertance, at f)
    source: dict = field(default_factory=dict)
    # D-P.3 fit uncertainty: covariance of (R0, K) from the steady fit's
    # residuals (2 x 2; [] = not available) and the standard error of each
    # R(f_i) from the period-to-period scatter of the fundamental fit
    cov_steady: list = field(default_factory=list)
    R_se: list = field(default_factory=list)

    def dz(self, f, extrap: str = "clamp"):
        """Dynamic excess dZ(f) (complex) for f > 0 (array). Inside the
        fitted range R and L are interpolated in log f. Outside it L is held
        and the in-phase excess Re dZ = R - r_qs is held (extrap "clamp") or
        scaled with sqrt(f / f_max) above and sqrt(f / f_min) below (extrap
        "sqrt", D-P.3: Stokes-layer growth, quasi-steady limit). Re dZ is
        clamped at >= 0 (D-P.5). Returns (dZ, mask of the harmonics outside
        the fitted range)."""
        if extrap not in EXTRAPOLATIONS:
            raise ValueError(f"extrap {extrap!r}: one of {EXTRAPOLATIONS}")
        f = np.atleast_1d(np.asarray(f, float))
        if not self.freqs:
            return np.zeros(f.shape, complex), np.zeros(f.shape, bool)
        fr = np.asarray(self.freqs, float)
        o = np.argsort(fr)
        lf = np.log(fr[o])
        R = np.interp(np.log(f), lf, np.asarray(self.R, float)[o])
        L = np.interp(np.log(f), lf, np.asarray(self.L, float)[o])
        r_qs = self.R0 + 2.0 * self.K * self.qbar_fit
        lo = f < fr.min() * (1 - 1e-12)
        hi = f > fr.max() * (1 + 1e-12)
        # D-P.5: the in-phase excess is never negative (a fitted R(f) below
        # the quasi-steady r_qs is fit noise: an oscillating flow cannot
        # dissipate less than its quasi-steady limit); clamped at 0, so a
        # pulsed waveform is never cheaper than <Q^2> at R0 (study rule)
        ex = np.maximum(R - r_qs, 0.0)
        if extrap == "sqrt":
            ex = np.where(hi, ex * np.sqrt(f / fr.max()), ex)
            ex = np.where(lo, ex * np.sqrt(f / fr.min()), ex)
        return ex + 1j * 2 * np.pi * f * L, (lo | hi)

    def below_qs(self):
        """D-P.5: the fitted frequencies whose R(f) lies below r_qs (their
        in-phase excess is clamped at 0 by dz; reported with the fit)."""
        r_qs = self.R0 + 2.0 * self.K * self.qbar_fit
        return [float(f) for f, r in zip(self.freqs, self.R) if float(r) < r_qs]

    def to_json(self):
        return asdict(self)

    @classmethod
    def from_json(cls, d):
        return cls(**d)

    # ---- parameter vector of E_diss (D-P.3 uncertainty) --------------------
    def theta(self):
        """(R0, K, R(f_1), ..., R(f_m)): the parameters E_diss depends on."""
        return np.array([self.R0, self.K] + [float(r) for r in self.R], float)

    def with_theta(self, th):
        th = np.asarray(th, float)
        return Impedance(R0=float(th[0]), K=float(th[1]), qbar_fit=self.qbar_fit,
                         freqs=list(self.freqs), R=[float(x) for x in th[2:]],
                         L=list(self.L), source=self.source,
                         cov_steady=self.cov_steady, R_se=self.R_se)

    def theta_cov(self):
        """Covariance of theta: blockdiag(cov(R0, K), diag(R_se^2)) (the
        steady and oscillatory fits use different runs: independent)."""
        m = len(self.R)
        if not self.cov_steady or len(self.R_se) != m:
            raise ValueError("impedance fit carries no uncertainty (cov_steady, R_se): "
                             "refit with fit_cases (>= 3 steady flows, >= 2 whole periods)")
        c = np.zeros((2 + m, 2 + m))
        c[:2, :2] = np.asarray(self.cov_steady, float)
        c[2:, 2:] = np.diag(np.asarray(self.R_se, float) ** 2)
        return c


def fit_steady(q, dp):
    """Least squares dp = R0 q + K q|q| (K = 0 with one flow), through the
    D-P.8 sign guard (fit_steady_guarded): R0 <= 0 and K < 0 beyond 2 SE
    are refused, K < 0 within 2 SE is clamped to 0 and R0 refitted."""
    g = fit_steady_guarded(q, dp)
    return g["R0"], g["K"]


def fit_steady_guarded(q, dp, k_se_mult: float | None = None) -> dict:
    """D-P.8 (PROJECT.md §2.5a, 2026-10-02): the sign guard of the steady
    impedance fit dp = R0 q + K q|q|. The energy bound of the primary basis
    (E_diss of a pulsed arm >= E_diss of steady flow at its mean, Jensen)
    rests on R0 > 0 and K >= 0:
      * R0 <= 0 is refused (ValueError), before and after a clamp;
      * K < 0 beyond FIT_K_NEG_SE standard errors (K + 2 SE(K) < 0) is
        refused, and so is any K < 0 whose SE is unknown (two flows: no
        residual dof);
      * K < 0 within 2 SE is clamped to 0 and R0 refitted by least squares
        on dp = R0 q (K_clamped True; the unclamped K and its SE are
        reported). The covariance then carries var(R0) of the refit and
        the unclamped var(K): the clamp does not remove K's uncertainty;
      * a K < 0 whose quadratic term is below FIT_K_ROUNDOFF_REL (1e-9) of
        the linear one at the largest fitted flow is the round-off of an
        exactly linear law (its SE is round-off too) and is clamped the same
        way (K_roundoff True).
    Returns dict(R0, K, cov (2 x 2 or None), K_fit, K_se, K_clamped,
    K_roundoff, n)."""
    mult = FIT_K_NEG_SE if k_se_mult is None else float(k_se_mult)
    q, dp = np.asarray(q, float), np.asarray(dp, float)
    if q.size == 0 or q.size != dp.size:
        raise ValueError("impedance fit: q and dp must be non-empty and of equal length")
    if not (np.all(np.isfinite(q)) and np.all(np.isfinite(dp))):
        raise ValueError("impedance fit: non-finite input")
    r0, k, cov = fit_steady_cov(q, dp)
    k_se = math.sqrt(cov[1][1]) if cov is not None and cov[1][1] >= 0 else None
    out = dict(R0=r0, K=k, cov=cov, K_fit=k, K_se=k_se, K_clamped=False,
               K_roundoff=False, n=int(q.size))
    # a quadratic term below FIT_K_ROUNDOFF_REL of the linear one at every
    # fitted flow is round-off of an exactly linear law (its SE is round-off
    # too): clamped like a K within 2 SE, never refused
    roundoff = (k < 0 and abs(k) * float(np.max(np.abs(q)))
                <= FIT_K_ROUNDOFF_REL * abs(r0))
    if k < 0:
        if not roundoff and k_se is None:
            raise ValueError(
                f"impedance fit (D-P.8): K = {k:.4g} Pa s2/m6 < 0 with no standard "
                f"error ({q.size} steady flows; >= 3 are needed to bound its sign): "
                "refused")
        if not roundoff and k + mult * k_se < 0:
            raise ValueError(
                f"impedance fit (D-P.8): K = {k:.4g} Pa s2/m6 < 0 beyond {mult:g} SE "
                f"(SE {k_se:.4g}): a negative quadratic resistance is not a "
                "physical housing; refused")
        # within 2 SE (or round-off): clamp K to 0 and refit R0 on dp = R0 q
        qq = float(q @ q)
        r0c = float(q @ dp) / qq
        res = dp - r0c * q
        var_r0 = (float(res @ res) / (q.size - 1)) / qq if q.size > 1 else 0.0
        out.update(R0=r0c, K=0.0, K_clamped=True, K_roundoff=bool(roundoff),
                   cov=([[var_r0, 0.0], [0.0, float(k_se) ** 2]]
                        if k_se is not None else None))
    if not out["R0"] > 0:
        raise ValueError(f"impedance fit (D-P.8): R0 = {out['R0']:.4g} Pa s/m3 <= 0: "
                         "refused")
    return out


def fit_steady_cov(q, dp):
    """The UNGUARDED least squares of fit_steady plus the covariance of
    (R0, K): s^2 (A^T A)^-1 with s^2 = RSS / (n - 2) (None with fewer than 3
    flows: no residual dof). The registered fit is fit_steady_guarded
    (D-P.8), which calls this."""
    q, dp = np.asarray(q, float), np.asarray(dp, float)
    if q.size == 1:
        return float(dp[0] / q[0]), 0.0, None
    A = np.column_stack([q, q * np.abs(q)])
    (r0, k), *_ = np.linalg.lstsq(A, dp, rcond=None)
    cov = None
    if q.size > 2:
        res = dp - A @ np.array([r0, k])
        s2 = float(res @ res) / (q.size - 2)
        cov = (s2 * np.linalg.inv(A.T @ A)).tolist()
    return float(r0), float(k), cov


def whole_periods(t, period, t_end=None):
    """Mask of the largest whole number of periods ending at t_end."""
    t = np.asarray(t, float)
    t_end = t[-1] if t_end is None else t_end
    n = math.floor((t_end - t[0]) / period + 1e-9)
    if n < 1:
        raise ValueError("less than one whole period of samples")
    t0 = t_end - n * period
    return (t > t0 + 1e-12) & (t <= t_end + 1e-12), n


def _fundamental(t, q, dp, w):
    A = np.column_stack([np.ones_like(t), np.sin(w * t), np.cos(w * t)])
    cq, *_ = np.linalg.lstsq(A, q, rcond=None)
    cp, *_ = np.linalg.lstsq(A, dp, rcond=None)
    # x = x0 + xs sin(wt) + xc cos(wt) = x0 + Re[(xc - i xs) e^{iwt}]
    qhat = cq[2] - 1j * cq[1]
    phat = cp[2] - 1j * cp[1]
    return complex(phat / qhat), float(cq[0]), float(cp[0])


def fit_oscillation(t, q, dp, freq):
    """Fundamental of dp against the oscillating flow over whole periods:
    dp ~ a0 + a sin(wt) + b cos(wt) with q ~ q0 + qs sin + qc cos.
    Returns (Z = R + i w L, qbar, mean dp)."""
    t, q, dp = (np.asarray(x, float) for x in (t, q, dp))
    m, _ = whole_periods(t, 1.0 / freq)
    return _fundamental(t[m], q[m], dp[m], 2 * np.pi * freq)


def oscillation_se(t, q, dp, freq):
    """D-P.3: standard error of R = Re Z from the period-to-period scatter
    of the fundamental fit (Z of each whole period; SD / sqrt(n)); None
    with fewer than 2 whole periods."""
    t, q, dp = (np.asarray(x, float) for x in (t, q, dp))
    P = 1.0 / freq
    m, n = whole_periods(t, P)
    if n < 2:
        return None
    t, q, dp = t[m], q[m], dp[m]
    t0 = t[0] - (t[1] - t[0])
    zs = []
    for j in range(n):
        s = (t > t0 + j * P + 1e-12) & (t <= t0 + (j + 1) * P + 1e-12)
        if s.sum() >= 4:
            zs.append(_fundamental(t[s], q[s], dp[s], 2 * np.pi * freq)[0].real)
    if len(zs) < 2:
        return None
    return float(np.std(zs, ddof=1) / math.sqrt(len(zs)))


def model_dp(imp, q_period, period, extrap: str = "clamp"):
    """Pressure difference of the impedance model for one period of q.
    Returns (dp, dp_res = its in-phase part, the |Qhat|^2 |dZ|-weighted
    share of the harmonics outside the fitted range (reactive part
    included), the share of the dynamic DISSIPATIVE work carried by them,
    sum |Re dZ_k| |Qhat_k|^2 outside / over all harmonics)."""
    q = np.asarray(q_period, float)
    n = q.size
    Q = np.fft.rfft(q)
    k = np.arange(Q.size)
    dyn = np.zeros(Q.size, complex)
    dyn_r = np.zeros(Q.size, complex)
    xfrac = dfrac = 0.0
    if period is not None and Q.size > 1:
        dz, out = imp.dz(k[1:] / period, extrap)
        dyn[1:] = dz * Q[1:]
        dyn_r[1:] = dz.real * Q[1:]
        a2 = np.abs(Q[1:]) ** 2
        wk = a2 * np.abs(dz)
        tot = float(wk.sum())
        xfrac = float(wk[out].sum()) / tot if tot > 0 else 0.0
        wd = a2 * np.abs(dz.real)
        totd = float(wd.sum())
        dfrac = float(wd[out].sum()) / totd if totd > 0 else 0.0
    dp_dyn = np.fft.irfft(dyn, n=n)
    dp = imp.R0 * q + imp.K * q * np.abs(q) + dp_dyn
    dp_res = imp.R0 * q + imp.K * q * np.abs(q) + np.fft.irfft(dyn_r, n=n)
    return dp, dp_res, xfrac, dfrac


def primary_energy(imp, q_of_t, period, t_total, n=N_SYNTH, q_on=None,
                   n_open=None, extrap: str = "clamp"):
    """Energies of the impedance model over t_total for a periodic flow
    q_of_t (steady: period None): E_primary = E_diss = E_visc (D-P.1, the
    primary basis), E_pump = E_nr = int max(0, Q dp) dt (non-regenerative
    drive, domain-inertance dependent, reported), E_inertial = E_nr -
    E_diss. ``extrap``: the evaluation outside the fitted range (D-P.3;
    primary_energy_conservative takes the larger of both)."""
    if period is None:
        q = float(np.asarray(q_of_t(np.array([0.0])))[0])
        dp = imp.R0 * q + imp.K * q * abs(q)
        e = max(0.0, q * dp) * t_total
        ed = q * dp * t_total
        return dict(E_primary=ed, E_diss=ed, E_pump=e, E_nr=e, E_visc=ed,
                    E_inertial=0.0, E_valve_KE=0.0, periods=None,
                    extrapolated_fraction=0.0,
                    dissipative_extrapolated_share=0.0, extrap=extrap,
                    qbar=q, q2bar=q * q)
    t = np.arange(n) * (period / n)
    q = np.asarray(q_of_t(t), float)
    dp, dp_res, xfrac, dfrac = model_dp(imp, q, period, extrap)
    p = q * dp
    e_per = float(np.sum(np.maximum(p, 0.0)) * period / n)
    e_vis = float(np.sum(q * dp_res) * period / n)
    nper = t_total / period
    out = dict(E_primary=e_vis * nper, E_diss=e_vis * nper,
               E_pump=e_per * nper, E_nr=e_per * nper, E_visc=e_vis * nper,
               E_inertial=(e_per - e_vis) * nper, periods=nper,
               extrapolated_fraction=xfrac,
               dissipative_extrapolated_share=dfrac, extrap=extrap,
               qbar=float(q.mean()), q2bar=float(np.mean(q * q)))
    if q_on is not None and n_open is not None and imp.freqs:
        fr = np.asarray(imp.freqs)
        L1 = float(np.interp(np.log(1.0 / period), np.log(np.sort(fr)),
                             np.asarray(imp.L)[np.argsort(fr)]))
        out["E_valve_KE"] = n_open * 0.5 * L1 * q_on ** 2 * nper
    else:
        out["E_valve_KE"] = 0.0
    return out


def primary_energy_conservative(imp, q_of_t, period, t_total, n=N_SYNTH,
                                q_on=None, n_open=None):
    """D-P.3: the registered E_diss of a waveform = the larger of the clamp
    and sqrt evaluations (a pulsed arm is never made cheaper by the model's
    extrapolation). Returns the chosen evaluation's dict with E_diss_clamp,
    E_diss_sqrt and model_form_rel = |sqrt - clamp| / max."""
    ev = {x: primary_energy(imp, q_of_t, period, t_total, n, q_on, n_open, x)
          for x in EXTRAPOLATIONS}
    pick = max(EXTRAPOLATIONS, key=lambda x: ev[x]["E_diss"])
    out = dict(ev[pick])
    a, b = ev["clamp"]["E_diss"], ev["sqrt"]["E_diss"]
    out.update(E_diss_clamp=a, E_diss_sqrt=b,
               model_form_rel=(abs(a - b) / max(abs(a), abs(b))
                               if max(abs(a), abs(b)) > 0 else 0.0))
    return out


def _ediss(imp, wave, t_total, extrap):
    q_of_t, period = wave[0], wave[1]
    return primary_energy(imp, q_of_t, period, t_total, extrap=extrap)["E_diss"]


def ediss_uncertainty(imp, waves, t_total):
    """D-P.3: E_diss (conservative evaluation) of each waveform and the
    Jacobian-propagated covariance of the fit parameters theta = (R0, K,
    R(f_i)). ``waves`` maps a name to (q_of_t, period, ...) (waveform()'s
    tuple). Returns dict(E={name: E}, J={name: dE/dtheta}, cov=theta cov,
    u={name: standard u(E)}). The derivative is taken on the evaluation
    (clamp or sqrt) chosen at the nominal parameters."""
    th = imp.theta()
    cov = imp.theta_cov()
    E, J, u, pick = {}, {}, {}, {}
    for name, w in waves.items():
        ev = {x: _ediss(imp, w, t_total, x) for x in EXTRAPOLATIONS}
        pick[name] = max(EXTRAPOLATIONS, key=lambda x: ev[x])
        E[name] = ev[pick[name]]
        g = np.zeros(th.size)
        for j in range(th.size):
            # E_diss is linear in theta (dp_res is linear in R0, K and Re dZ),
            # so the central difference is exact for any step (piecewise: at
            # an R(f_i) clamped at r_qs, D-P.5, the one-sided slopes differ)
            h = FD_REL_STEP * max(abs(th[j]), 1.0)
            tp, tm = th.copy(), th.copy()
            tp[j] += h
            tm[j] -= h
            g[j] = (_ediss(imp.with_theta(tp), w, t_total, pick[name])
                    - _ediss(imp.with_theta(tm), w, t_total, pick[name])) / (2 * h)
        J[name] = g
        u[name] = float(math.sqrt(max(0.0, g @ cov @ g)))
    return dict(E=E, J=J, cov=cov, u=u, extrap=pick)


def he_energy_u95(imp, wave_s, wave_p, t_total):
    """NEW-13.2 / D-P.3: U95 of E_diss(S) - E_diss(P) for the H_E front test:
    max(k u(E_S - E_P), 2 % of the larger energy), u from the joint
    Jacobian (the shared R0, K correlate the two energies), k = 2.
    Returns dict(E_S, E_P, diff, u, U95, resolved = diff > U95)."""
    r = ediss_uncertainty(imp, {"S": wave_s, "P": wave_p}, t_total)
    g = r["J"]["S"] - r["J"]["P"]
    u = float(math.sqrt(max(0.0, g @ r["cov"] @ g)))
    es, ep = r["E"]["S"], r["E"]["P"]
    U = max(ENERGY_K * u, ENERGY_FLOOR_REL * max(abs(es), abs(ep)))
    return dict(E_S=es, E_P=ep, diff=es - ep, u=u, U95=U, resolved=bool(es - ep > U))


def model_check(e_model: float, e_port: float, e_edge: float) -> dict:
    """D-P.3 model check on an isothermal square-wave run of the impedance
    block (C1, C2'): over the same whole periods, the model's conservative
    E_diss against the measured port work less the valve-edge temporal
    dissipation estimate (KE_exchanged / n_edge^2, pump_energy.measured).
    Passes iff e_model >= 0.98 x (e_port - e_edge); the arm's E_diss is then
    the model's, else the measured value (the larger: the check can only
    make the arm costlier)."""
    vals = (e_model, e_port, e_edge)
    if not all(math.isfinite(float(v)) for v in vals):
        raise ValueError("model check: non-finite input")
    meas = float(e_port) - float(e_edge)
    if not meas > 0:
        raise ValueError(f"model check: measured edged work {meas:g} J <= 0")
    ok = float(e_model) >= MODEL_CHECK_MIN_RATIO * meas
    return dict(passed=bool(ok), e_model=float(e_model), e_measured=meas,
                ratio=float(e_model) / meas,
                e_used=float(e_model) if ok else max(float(e_model), meas))


def jensen_floor(imp, q_of_t, period, t_total, e_diss=None) -> dict:
    """D-P.8: the energy bound of the primary basis, checked after the fit
    guard: a waveform's E_diss is >= the E_diss of steady flow at its own
    mean flow (<Q^2> >= Qbar^2, <Q^3> >= Qbar^3 for Q >= 0, Re dZ >= 0 by
    D-P.5). Returns dict(E_diss, E_steady, ratio); raises ValueError when the
    bound fails beyond round-off (a pulsed arm cheaper than steady flow at
    the same mean is a fit or model error, never a result)."""
    if e_diss is None:
        e_diss = primary_energy_conservative(imp, q_of_t, period, t_total)["E_diss"]
    if period is None:
        qbar = float(np.asarray(q_of_t(np.array([0.0])))[0])
    else:
        t = np.arange(N_SYNTH) * (period / N_SYNTH)
        qbar = float(np.mean(np.asarray(q_of_t(t), float)))
    e_st = (imp.R0 * qbar + imp.K * qbar * abs(qbar)) * qbar * t_total
    if e_diss < e_st * (1 - JENSEN_RTOL) - 1e-300:
        raise ValueError(f"D-P.8 energy bound violated: E_diss {e_diss:.6g} J < "
                         f"E_diss of steady flow at the same mean {e_st:.6g} J "
                         f"(R0 {imp.R0:.4g}, K {imp.K:.4g}): the fit or the "
                         "model is wrong")
    return dict(E_diss=float(e_diss), E_steady=float(e_st),
                ratio=float(e_diss / e_st) if e_st > 0 else math.inf)


def q_eq_registered(imp, energy, t_total):
    """D-P.7 (4): the Q_eq of a Stage-B impedance fit (one with oscillatory
    runs, fit_cases) is refused until the sine closure re-examination
    (face-based Phi; Stokes-layer vs limitedLinear separation) has reported
    and is recorded in the fit's source as ``dp7_reexamination`` (the report
    and its dated row). A model without oscillatory runs (the provisional
    lumped impedance) is not a Stage-B fit and passes through."""
    src = imp.source or {}
    if (imp.freqs or src.get("osc")) and not src.get("dp7_reexamination"):
        raise ValueError("Q_eq refused (D-P.7 (4)): this impedance fit has oscillatory "
                         "runs and no sine closure re-examination record "
                         "(source.dp7_reexamination); no Stage-B fit sets a Q_eq "
                         "before the re-examination reports")
    return q_eq(imp, energy, t_total)


def q_eq(imp, energy, t_total):
    """Steady flow with the same primary energy (D-H.1 comparators: pass the
    arm's E_primary = E_diss, D-P.1): (R0 Q + K Q^2) Q t = E. The registered
    entry for a Stage-B fit is q_eq_registered (D-P.7 (4))."""
    if energy <= 0:
        return 0.0
    c = np.roots([imp.K * t_total, imp.R0 * t_total, 0.0, -energy]) \
        if imp.K else np.array([math.sqrt(energy / (imp.R0 * t_total))])
    r = [x.real for x in np.atleast_1d(c) if abs(x.imag) < 1e-9 * abs(x) and x.real > 0]
    if not r:
        raise ValueError("no positive steady flow for this energy")
    return float(min(r))


# --------------------------------------------------------------------------
# waveforms (same conventions as housing_physics / qs_model)
# --------------------------------------------------------------------------

def waveform(kind, qbar, freq=None, amp=None, duty=None, edge=None):
    """(q_of_t, period, q_on, n_open_per_period) in SI."""
    if kind == "steady":
        return (lambda t: np.full_like(np.asarray(t, float), qbar)), None, None, None
    P = 1.0 / freq
    if kind == "sine":
        return (lambda t: qbar * (1 + amp * np.sin(2 * np.pi * freq * np.asarray(t, float))),
                P, None, None)
    q_on = qbar / duty
    tau = edge if edge else 0.0
    h = 0.5 * tau
    tt = [0.0, h, duty * P - h, duty * P + h, P - h, P]
    qq = [0.5 * q_on, q_on, q_on, 0.0, 0.0, 0.5 * q_on]
    if tau == 0:
        return (lambda t: np.where(np.mod(np.asarray(t, float), P) < duty * P, q_on, 0.0),
                P, q_on, 1)
    return (lambda t: np.interp(np.mod(np.asarray(t, float), P), tt, qq)), P, q_on, 1


def waveform_from_meta(meta):
    kind = meta["waveform"]
    edge = meta["lattice"].get("valve_edge_s")
    return waveform(kind, meta["qbar_m3s"], meta.get("freq_hz") or None,
                    meta.get("amp"), meta.get("duty"), edge)


# --------------------------------------------------------------------------
# the closure split per waveform class (D-P.7; moved here from
# scripts/gdt/gdt2_analysis.py, the GDT2 deep-dive, which imports them back)
# --------------------------------------------------------------------------

def step_index(t, dt):
    """Integer step number of times on the fixed-dt lattice."""
    return np.rint(np.asarray(t, float) / dt).astype(np.int64)


def phase_bin(t, dt, steps_per_period, nbins):
    """Phase bin (0..nbins-1) of each time by integer step-in-period
    arithmetic: bin k holds steps [k, k+1) * steps_per_period / nbins of the
    period (no floating-point phase rounding)."""
    s = step_index(t, dt) % steps_per_period
    return (s * nbins) // steps_per_period


def interval_closure(t_port, p_port, t_h, ke, phi, wb):
    """Per-interval hydraulic closure on the hydraulic-monitor grid t_h:
    E_port over [t_k, t_k+1] (trapezoid on the every-step port series),
    Phi and W_b by trapezoid of their two end samples, dKE the change. The
    sums over any run of intervals equal measured()'s window integrals.
    Returns a dict of arrays (t0, t1, ep, phi, wb, dke, r)."""
    t_port = np.asarray(t_port, float)
    p_port = np.asarray(p_port, float)
    t_h = np.asarray(t_h, float)
    phi = np.asarray(phi, float)
    wb = np.asarray(wb, float)
    ke = np.asarray(ke, float)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (p_port[1:] + p_port[:-1])
                                           * np.diff(t_port))])
    idx = np.searchsorted(t_port, t_h - 1e-9)
    ok = idx < t_port.size
    ok[ok] &= np.abs(t_port[idx[ok]] - t_h[ok]) <= 1e-9 * np.maximum(1.0, t_h[ok])
    if not ok.all():
        raise ValueError(f"{int((~ok).sum())} hydraulic samples are not port "
                         "samples (the port series must be every step)")
    dtk = np.diff(t_h)
    ep = np.diff(cum[idx])
    phik = 0.5 * (phi[1:] + phi[:-1]) * dtk
    wbk = 0.5 * (wb[1:] + wb[:-1]) * dtk
    dke = np.diff(ke)
    return dict(t0=t_h[:-1], t1=t_h[1:], ep=ep, phi=phik, wb=wbk, dke=dke,
                r=ep - (phik - wbk + dke))


def window_sums(iv, a, b, keys=("ep", "phi", "wb", "dke", "r")):
    """Sums of the interval terms over the intervals inside [a, b]."""
    m = (iv["t0"] >= a - 1e-9) & (iv["t1"] <= b + 1e-9)
    return {k: float(iv[k][m].sum()) for k in keys}


def phase_sums(iv, a, b, dt, steps_per_period, nbins,
               keys=("ep", "phi", "wb", "dke", "r")):
    """Per-period sums by phase bin of the interval start over whole periods
    of [a, b] (b - a a whole number of periods)."""
    m = (iv["t0"] >= a - 1e-9) & (iv["t1"] <= b + 1e-9)
    nper = (b - a) / (steps_per_period * dt)
    if abs(nper - round(nper)) > 1e-6:
        raise ValueError(f"[{a}, {b}] is not a whole number of periods")
    pb = phase_bin(iv["t0"][m], dt, steps_per_period, nbins)
    return {k: [float(iv[k][m][pb == j].sum() / nper) for j in range(nbins)]
            for k in keys}


def square_segments(steps_per_period, edge_steps, on_steps):
    """Segment label per step-in-period of a trapezoid square wave with the
    waveform() convention (opening edge centred on t = 0, closing edge
    centred on t = D P): a step is labelled by the interval it ENDS, i.e. the
    interval (s - 1, s] in steps (the state after step s)."""
    h = edge_steps // 2
    lab = []
    for s in range(steps_per_period):
        e = (s - 1) % steps_per_period          # interval [e, e + 1]
        if e >= steps_per_period - h or e < h:
            lab.append("opening_edge")
        elif e < on_steps - h:
            lab.append("on")
        elif e < on_steps + h:
            lab.append("closing_edge")
        else:
            lab.append("off")
    return lab


SQUARE_SEGMENTS = ("opening_edge", "on", "closing_edge", "off")


def _rel(a, b):
    return float(a) / float(b) if b else math.nan


def sub_windows(t0, t1, period, min_s):
    """D-P.7 (3b): consecutive sub-windows of [t0, t1] of the smallest whole
    number of periods lasting >= min_s (steady: min_s), the remainder merged
    into the last one, so every sub-window lasts >= min_s and they cover the
    window. [] if the window is shorter than one sub-window."""
    span = (math.ceil(min_s / period - 1e-9) * period) if period else min_s
    n = math.floor((t1 - t0) / span + 1e-9)
    if n < 1:
        return []
    edges = [t0 + k * span for k in range(n)] + [t1]
    return list(zip(edges[:-1], edges[1:]))


def sine_closure_checks(iv, t0, t1, meta, isothermal=True):
    """D-P.7 (3): the necessary checks of a sine window's closure against the
    cell-centred Phi (the closure itself is REPORTED, not gated):
    (a) sign: the per-period residual summed over each of
        gates.CLOSURE_PHASE_BINS phase bins of whole periods is
        >= gates.CLOSURE_SIGN_MIN_REL x the per-period Phi in every bin;
    (b) stationarity: r/Phi of every whole-period sub-window of >=
        gates.CLOSURE_SUB_WINDOW_S lies within gates.CLOSURE_STATIONARITY_PP
        of the window's r/Phi (not evaluable with fewer than two
        sub-windows: reported, as an undeveloped window reports the closure);
    (c) dt invariance at G-dt is a pair check (sine_dt_invariance).
    They catch temporal, phase-offset and start-up errors; a stationary
    spatial excess passes all three, so they cannot bound a bias in E_port
    (D-P.7 (4)). Gated (by housing_metrics) on isothermal windows only."""
    import gates as G
    dt = float(meta["lattice"]["dt"])
    P = float(meta["period_s"])
    spp = int(round(P / dt))
    nbins = G.CLOSURE_PHASE_BINS
    w = window_sums(iv, t0, t1)
    nper = (t1 - t0) / P
    ps = phase_sums(iv, t0, t1, dt, spp, nbins)
    phi_per = w["phi"] / nper
    sign = G.closure_phase_sign_gate(ps["r"], phi_per)
    subs = []
    for a, b in sub_windows(t0, t1, P, G.CLOSURE_SUB_WINDOW_S):
        s = window_sums(iv, a, b)
        subs.append(dict(t0=float(a), t1=float(b), r=s["r"], phi=s["phi"],
                         r_over_phi=_rel(s["r"], s["phi"])))
    win_rel = _rel(w["r"], w["phi"])
    if len(subs) >= 2:
        stat = G.closure_stationarity_gate([s["r_over_phi"] for s in subs], win_rel)
        stat_ok, stat_det = bool(stat.ok), stat.detail
    else:
        stat_ok, stat_det = None, (f"not evaluable: {len(subs)} sub-window(s) of >= "
                                   f"{G.CLOSURE_SUB_WINDOW_S:g} s in [{t0:g}, {t1:g}] s "
                                   "(two are needed)")
    return dict(
        gated=bool(isothermal), r_over_phi=win_rel, periods=nper,
        r_per_period=w["r"] / nper, phi_per_period=phi_per,
        sign=dict(ok=bool(sign.ok), detail=sign.detail, bins=nbins,
                  r_bins_per_period=ps["r"], phi_bins_per_period=ps["phi"],
                  min_bin_rel=(min(ps["r"]) / phi_per if phi_per else math.nan),
                  threshold_rel=G.CLOSURE_SIGN_MIN_REL),
        stationarity=dict(ok=stat_ok, detail=stat_det, sub_windows=subs,
                          window_r_over_phi=win_rel,
                          max_abs_dev=(max(abs(s["r_over_phi"] - win_rel) for s in subs)
                                       if subs else None),
                          tol=G.CLOSURE_STATIONARITY_PP,
                          min_s=G.CLOSURE_SUB_WINDOW_S))


def square_closure_split(iv, t0, t1, meta, ke_exchanged):
    """D-P.7 (5): a square wave's closure (reported) split by valve segment,
    per period over whole periods of [t0, t1] (hydraulic monitors every
    step): opening edge, on-plateau, closing edge, off-plateau. The edge
    segments are set against the backward edge model KE_exchanged / n^2
    after removing the plateau-like excess (the edges' Phi x the plateaus'
    pooled r/Phi); the on- and off-plateau remainders are reported as
    fractions of their Phi."""
    lat = meta["lattice"]
    dt = float(lat["dt"])
    P = float(meta["period_s"])
    spp = int(round(P / dt))
    n_edge, on_steps = lat.get("valve_edge_steps"), lat.get("on_steps")
    if not n_edge or not on_steps:
        return dict(note="no valve-edge lattice (ideal edges): segments not split")
    m = (iv["t0"] >= t0 - 1e-9) & (iv["t1"] <= t1 + 1e-9)
    if np.any(np.abs(np.diff(iv["t0"][m]) - dt) > 1e-9 * max(1.0, dt)):
        return dict(note="segments need every-step hydraulic monitors "
                         "(hydraulic_steps 1): not split")
    lab = np.array(square_segments(spp, n_edge, on_steps))
    nper = (t1 - t0) / P
    seg = lab[step_index(iv["t1"][m], dt) % spp]
    out = {}
    for name in SQUARE_SEGMENTS:
        k = seg == name
        out[name] = {q: float(iv[q][m][k].sum() / nper)
                     for q in ("ep", "phi", "wb", "dke", "r")}
        out[name]["steps"] = int(round(k.sum() / nper))
    on, off = out["on"], out["off"]
    plat_rel = _rel(on["r"] + off["r"], on["phi"] + off["phi"])
    e_r = out["opening_edge"]["r"] + out["closing_edge"]["r"]
    e_phi = out["opening_edge"]["phi"] + out["closing_edge"]["phi"]
    model = (float(ke_exchanged) / nper) / n_edge ** 2
    less = e_r - e_phi * plat_rel
    return dict(segments=out, periods=nper, valve_edge_steps=int(n_edge),
                r_per_period=sum(out[s]["r"] for s in SQUARE_SEGMENTS),
                on_r_over_phi=_rel(on["r"], on["phi"]),
                off_r_over_phi=_rel(off["r"], off["phi"]),
                plateau_r_over_phi=plat_rel,
                edges_r=e_r, edges_phi=e_phi,
                edges_plateau_like=e_phi * plat_rel,
                edges_less_plateau_like=less,
                edge_model_per_period=model,
                edges_vs_model_rel=(less / model - 1.0) if model else math.nan,
                edges_share_of_r=_rel(e_r, e_r + on["r"] + off["r"]))


def excess_over_steady(m_pulsed, m_steady):
    """D-P.7 (6): per-period closure excess of a pulsed window over a steady
    reference window of the same kind (measured() dicts): the steady
    residual subtracted per unit time and in proportion to Phi (the two
    conventions differ by ~13 % on GDT2's C3)."""
    w0, w1 = m_pulsed["window"]
    s0, s1 = m_steady["window"]
    P = m_pulsed.get("period_s")
    if not P:
        raise ValueError("excess over steady: the pulsed window has no period")
    nper = (w1 - w0) / P
    r_rate = m_steady["closure_residual"] / (s1 - s0)
    r_frac = _rel(m_steady["closure_residual"], m_steady["Phi"])
    r_per = m_pulsed["closure_residual"] / nper
    return dict(per_period=r_per,
                excess_per_time=r_per - r_rate * P,
                excess_per_phi=r_per - r_frac * m_pulsed["Phi"] / nper,
                steady_rel_phi=r_frac, steady_rate_W=r_rate)


def sine_dt_invariance(case_dt, case_half, stage="pre"):
    """D-P.7 (3c): the dt invariance of a sine class's isothermal closure at
    G-dt: the r/Phi of the [ramp, t_pre] window (stage 'pre'; 'run' for an
    isothermal run's window) of the dt rung against its dt/2 partner,
    gates.closure_dt_invariance_gate. Returns (GateResult, report)."""
    import gates as G
    a, b = measured(case_dt, stage), measured(case_half, stage)
    for m, c in ((a, case_dt), (b, case_half)):
        if m["waveform"] != "sine" or not m["developed"]:
            raise ValueError(f"{c}: not a developed sine closure window ({m['gate_note']})")
    ra, rb = _rel(a["closure_residual"], a["Phi"]), _rel(b["closure_residual"], b["Phi"])
    gate = G.closure_dt_invariance_gate(ra, rb)
    return gate, dict(r_over_phi_dt=ra, r_over_phi_half=rb, diff=ra - rb,
                      windows=[a["window"], b["window"]])


def sine_closure_dt(r_over_phi_by_dt):
    """D-P.7 (3c): the largest ladder dt whose sine closure passes the dt
    invariance against its half (dt -> r/Phi of the class's isothermal
    window). Outside it the excess is temporal and the sine class's
    port-work quantities (Stage-B oscillatory impedance runs, D-H.1 R(w))
    run at that dt. Returns dict(dt, pairs) (dt None if no pair passes)."""
    import gates as G
    d = {float(k): float(v) for k, v in r_over_phi_by_dt.items()}
    pairs = {}
    for dt in sorted(d, reverse=True):
        half = next((h for h in d if abs(h - 0.5 * dt) <= 1e-9 * dt), None)
        if half is not None:
            pairs[dt] = bool(G.closure_dt_invariance_gate(d[dt], d[half]).ok)
    best = max((dt for dt, ok in pairs.items() if ok), default=None)
    return dict(dt=best, pairs=pairs)


# --------------------------------------------------------------------------
# measurements from housing cases
# --------------------------------------------------------------------------

def port_power(case, meta=None):
    """(t, P_port [W], Q_in [m3/s]) every step: P_port = -sum_ports sum
    phi (p*/rho + |U|^2/2), p* = p - rho_in gh."""
    meta = meta or HIO.case_meta(case)
    rho_in = meta["rho_inlet"]
    parts = []
    for port in ("inlet", "outlet"):
        t, cols = HIO.series(case, "fluid", f"pPhi_{port}")
        tk, k = HIO.column(case, "fluid", f"kPhi_{port}")
        parts.append((t, cols["weightedSum(p)"] - rho_in * cols["weightedSum(gh)"],
                      tk, k))
    tq, vq = HIO.column(case, "fluid", "volFlux_inlet")
    n = min(min(p[0].size, p[2].size) for p in parts)
    n = min(n, tq.size)
    t0 = parts[0][0][:n]
    for t, _, tk, _ in parts:
        if np.any(np.abs(t[:n] - t0) > 1e-9) or np.any(np.abs(tk[:n] - t0) > 1e-9):
            raise ValueError(f"{case}: port monitors are not on one time base")
    if np.any(np.abs(tq[:n] - t0) > 1e-9):
        raise ValueError(f"{case}: inlet flux monitor not on the port time base")
    tot = -sum(pw[:n] + k[:n] for _, pw, _, k in parts)
    return t0, tot, -vq[:n]


def closure_window(meta, t_ports, t_ke, stage="run"):
    """Window of the D-P hydraulic closure (re-specified in WP 2.5c, D-P
    closure row of docs/model/housing-physics.md section 11).

    Developed flow only, after the start-up ramp:
    - stage "run": after the NEW-14 pre-circulation for a run with a
      measurement or heated stage, [max(t_pre, ramp), t_end], or
      [ramp, t_pre] for a pre-circulation-only run (a shared NEW-14 state,
      the P00b valve check);
    - stage "pre": the isothermal pre-circulation [ramp, t_pre] of a run
      that goes on to heat (the port-work instrument check of a heated
      case, review SHOULD-FIX 2);
    whole periods ending at the last common port/KE sample inside the
    stage (steady: the whole window). Returns (t0, t1, developed, why)."""
    tol = 1e-9
    t_pre, t_end = float(meta["t_pre"]), float(meta["t_end"])
    ramp = float(meta["ramp_s"])
    t_last = min(float(t_ports[-1]), float(t_ke[-1]))
    if stage == "pre":
        t_last = min(t_last, t_pre)
    elif stage != "run":
        raise ValueError(f"closure_window: stage {stage!r} (run or pre)")
    inside = t_ke[t_ke <= t_last + tol]
    if inside.size == 0:
        return (float(t_ke[0]), float(t_ke[0]), False,
                f"no KE sample by {t_last:g} s (stage {stage})")
    t1 = float(inside[-1])                              # a KE sample
    t_a = (ramp if stage == "pre" or t_end <= t_pre + tol
           else max(t_pre, ramp))
    P = meta.get("period_s")
    # without a developed window the closure is reported (never gated) over
    # the stage after t_pre, or the whole run for a pre-circulation only
    # (or the pre-circulation stage)
    t_rep = (t_pre if stage == "run" and t1 > t_pre + tol
             else float(t_ke[0]))
    if t1 <= t_a + tol:
        return t_rep, t1, False, (f"no developed window: the run ends at "
                                  f"{t1:g} s, before {t_a:g} s (start-up "
                                  f"ramp {ramp:g} s, pre-circulation "
                                  f"{t_pre:g} s)")
    if P:
        n = math.floor((t1 - t_a) / P + 1e-9)
        if n < 1:
            return t_rep, t1, False, (f"no whole period in the developed "
                                      f"window [{t_a:g}, {t1:g}] s")
        return t1 - n * P, t1, True, f"{n} whole periods"
    return t_a, t1, True, f"{t1 - t_a:g} s steady"


def phi_face_series(case, meta=None):
    """(t, {column: array}) of the D-P.9 face-based Phi monitor
    (housing_physics.phi_face_monitor, PHI_FACE=1), or None when the case
    did not render it (caseMeta 'phi_face' absent or null). A case that
    rendered it but has no output, or lacks a column, is an input error."""
    meta = meta or HIO.case_meta(case)
    pf = meta.get("phi_face")
    if not pf:
        return None
    t, cols = HIO.series(case, "fluid", pf.get("monitor", "PhiFace_fluid"))
    missing = [k for k in PHI_FACE_KEYS if k not in cols]
    if missing:
        raise ValueError(f"{case}: the face-based Phi monitor lacks {missing}")
    return t, cols


def face_closure(case, meta, t_port, p_port, tk, ke, t0, t1, e_port, w_b,
                 dke, phi_cell, closure_cell, developed):
    """D-P.9 (D-P.7 (4)): the hydraulic closure against the face-based Phi
    (the dissipation the discrete viscous operator does, near-wall
    half-cells included; housing_physics.phi_face_monitor) beside the
    registered cell-centred one, over the same window. REPORTED, never
    gated: the cell-centred Phi stays the gate's instrument (D-P closure
    row, D-P.2/D-P.4/D-P.7). Hard failures (ValueError): the monitor off
    the KE time base, its identity column (the same dissipation from the
    solver's fvc operators) off Phi_face by more than PHI_FACE_IDENTITY_RTOL,
    or its cell-centred column off the registered Phi_fluid by more than
    PHI_FACE_CELL_RTOL over the window. Window integrals in J; the sine
    split (per period, sub-windows of >= gates.CLOSURE_SUB_WINDOW_S, phase
    bins) is reported for developed periodic windows. None without the
    monitor.

    W_bound (internal review, 6 Oct): the viscous traction's work
    at the ports, sum over the non-coupled boundary faces of t_b . U_b (zero
    on no-slip walls). The viscous term's work on the resolved KE is
    W_bound - Phi_face, so r_face + W_bound = E_port + W_b + W_bound -
    Phi_face - dKE is the remainder against the full discrete viscous work
    (``closure_residual_face_with_W_bound``). r_face itself is D-P.9's
    definition, unchanged; W_bound is reported beside it."""
    s = phi_face_series(case, meta)
    if s is None:
        return None
    tf, cols = s
    if not (tf.size == tk.size and np.allclose(tf, tk, rtol=0, atol=1e-9)):
        raise ValueError(f"{case}: the face-based Phi monitor is not on the KE "
                         "time base (hydraulic cadence)")
    m = (tf >= t0 - 1e-9) & (tf <= t1 + 1e-9)
    if m.sum() < 2:
        raise ValueError(f"{case}: fewer than two face-based Phi samples in "
                         f"[{t0:g}, {t1:g}] s")
    trap = lambda y, x: float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))  # noqa: E731
    w = {k: trap(np.asarray(cols[k], float)[m], tf[m]) for k in PHI_FACE_KEYS}
    face, ident = np.asarray(cols["Phi_face"], float), np.asarray(cols["Phi_identity"], float)
    scale = float(np.max(np.abs(face))) if face.size else 0.0
    id_dev = float(np.max(np.abs(ident - face)) / scale) if scale > 0 else 0.0
    if id_dev > PHI_FACE_IDENTITY_RTOL:
        raise ValueError(f"{case}: face-based Phi fails its identity check "
                         f"(max |Phi_identity - Phi_face| / max |Phi_face| = "
                         f"{id_dev:.3g} > {PHI_FACE_IDENTITY_RTOL:g})")
    if phi_cell:
        cell_dev = w["Phi_cell"] / phi_cell - 1.0
    else:                               # no flow: both zero, or a mismatch
        cell_dev = 0.0 if w["Phi_cell"] == 0 else math.inf
    if not (math.isfinite(cell_dev) and abs(cell_dev) <= PHI_FACE_CELL_RTOL):
        raise ValueError(f"{case}: the face-based monitor's cell-centred column "
                         f"differs from Phi_fluid by {cell_dev:.3g} over the window "
                         f"(> {PHI_FACE_CELL_RTOL:g}): not the same instrument base")
    pf = w["Phi_face"]
    r_face = e_port - (pf - w_b + dke)
    ref = max(abs(pf), abs(w_b))
    out = dict(
        error=None,
        window=[float(t0), float(t1)], **{k: w[k] for k in PHI_FACE_KEYS},
        Phi_cell_registered=float(phi_cell),
        identity_max_rel=id_dev, cell_vs_registered_rel=float(cell_dev),
        face_over_cell=_rel(pf, phi_cell),
        lap_face_over_cell=_rel(w["Phi_face_lap"], w["Phi_cell_lap"]),
        wall_share=_rel(w["Phi_face_wall"], pf),
        duct_share_face=_rel(w["Phi_face_duct"], pf),
        duct_share_cell=_rel(w["Phi_cell_duct"], phi_cell),
        closure_residual_face=r_face,
        closure_rel_face=r_face / ref if ref else math.nan,
        r_over_phi_face=_rel(r_face, pf),
        r_over_phi_cell=_rel(closure_cell, phi_cell),
        W_bound_over_phi_face=_rel(w["W_bound"], pf),
        closure_residual_face_with_W_bound=r_face + w["W_bound"],
        r_with_W_bound_over_phi_face=_rel(r_face + w["W_bound"], pf),
        sine=None, note=("reported, never gated (D-P.9): the cell-centred Phi "
                         "stays the gate's instrument"))
    P = meta.get("period_s")
    if developed and P and meta["waveform"] == "sine":
        import gates as G
        twb, wbs = HIO.column(case, "fluid", "Wb_fluid")
        if not (twb.size == tk.size and np.allclose(twb, tk, rtol=0, atol=1e-9)):
            raise ValueError(f"{case}: W_b is not on the KE time base")
        iv = interval_closure(t_port, p_port, tk, ke, cols["Phi_face"], wbs)
        dt = float(meta["lattice"]["dt"])
        spp = int(round(P / dt))
        nper = (t1 - t0) / P
        ws = window_sums(iv, t0, t1)
        subs = []
        for a, b in sub_windows(t0, t1, P, G.CLOSURE_SUB_WINDOW_S):
            sw = window_sums(iv, a, b)
            subs.append(dict(t0=float(a), t1=float(b),
                             r_over_phi=_rel(sw["r"], sw["phi"])))
        ps = phase_sums(iv, t0, t1, dt, spp, G.CLOSURE_PHASE_BINS)
        out["sine"] = dict(periods=nper, r_per_period=ws["r"] / nper,
                           phi_per_period=ws["phi"] / nper,
                           sub_windows=subs,
                           r_bins_per_period=ps["r"],
                           phi_bins_per_period=ps["phi"])
    return out


def measured(case, stage="run"):
    """Measured port work over the closure window (closure_window: developed
    flow, whole periods; stage "run" or "pre"), with the chimney term
    separated and the kinetic energy the flow exchanges through the ports
    (KE_exchanged), and the closure split per waveform class (D-P.7,
    PROJECT.md §2.5a, 2026-10-02):
    - STEADY: ``gated`` True for a developed isothermal window (an
      isothermal run, or the pre-circulation of any run): the closure checks
      the port-work instrument at the scale of Phi (D-P closure row, D-P.2 /
      D-P.4 tolerance per level);
    - SINE: the closure against the cell-centred Phi is REPORTED, not gated,
      on every level (GDT2: +25.7 % on medium at 12.5 ms, stationary and
      positive in every phase bin), with the D-P.7 (3) necessary checks
      ``sine_checks`` (phase-bin sign, sub-window stationarity), which
      housing_metrics gates on isothermal windows; the dt invariance is the
      pair check sine_dt_invariance;
    - SQUARE: reported, split by valve segment (``square_split``: the
      edges against KE_exchanged / n^2 after the plateau-like excess, the
      plateaus as fractions of their Phi).
    Heated windows are reported for every class (the residual measures the
    discretisation of the buoyant flow, P00b: +9 % -> +68 % of Phi as the
    plumes develop), pulsed ones with their per-period residual.
    D-P.9 (2026-10-06): a case with the face-based Phi monitor (PHI_FACE=1)
    also reports ``phi_face`` (face_closure): the same closure against the
    face-based Phi, never gated; None otherwise. A monitor whose data or
    own checks fail (face_closure's hard failures) gives ``phi_face`` =
    {error: ..., ...} and no face-based number, and the registered closure
    above is returned unchanged: an unregistered instrument never hides the
    registered one (internal review, 6 Oct)."""
    meta = HIO.case_meta(case)
    t, pw, q = port_power(case, meta)
    tk, ke = HIO.column(case, "fluid", "KE_fluid")
    t0, t1, developed, why = closure_window(meta, t, tk, stage)
    iso = meta["heat"]["kind"] == "none" or stage == "pre"
    square = meta["waveform"] == "square"
    m = (t >= t0 - 1e-9) & (t <= t1 + 1e-9)
    tt, pp = t[m], pw[m]
    trap = lambda y, x: float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))  # noqa: E731
    e_port = trap(pp, tt) if tt.size >= 2 else math.nan

    def vol_int(name):
        ts, v = HIO.column(case, "fluid", name)
        mm = (ts >= t0 - 1e-9) & (ts <= t1 + 1e-9)
        return trap(v[mm], ts[mm]) if mm.sum() >= 2 else math.nan

    w_b = vol_int("Wb_fluid")
    phi = vol_int("Phi_fluid")
    mk = (tk >= t0 - 1e-9) & (tk <= t1 + 1e-9)
    kw = ke[mk]
    dke = float(kw[-1] - kw[0]) if kw.size >= 2 else math.nan
    ke_tv = float(np.sum(np.abs(np.diff(kw)))) if kw.size >= 2 else math.nan
    # kinetic energy the flow exchanges through the ports over the window
    # (total variation of the fluid KE: uptake while the flow accelerates,
    # return while it decelerates). It is "returned" only in the fluid's own
    # mechanical-energy balance: a NON-REGENERATIVE pump (the D-P primary
    # basis) pays the opening edge's KE and recovers nothing at the closing
    # edge (E_inertial >= 0). Backward at fixed dt dissipates part of it:
    # for a lumped inertance driven through an n-step trapezoid edge the
    # loss per period is (2/n^2) x the KE jump, i.e. KE_exchanged / n^2
    # over whole periods (docs/model/housing-physics.md section 4.1). This
    # temporal dissipation is in the square-wave E_port, which therefore
    # errs on the costly side (never toward a saving). Smooth sines: <= 3 %
    # of Phi at 12.5 ms and 2 Hz by the same estimate.
    n_edge = meta["lattice"].get("valve_edge_steps")
    e_edge = (ke_tv / n_edge ** 2 if square and n_edge else math.nan)
    # port-based chimney estimate: (rho_in - rho_out) g dz Q, rho_out the
    # flux-weighted outlet density (mass / volume flux)
    ts, mo = HIO.column(case, "fluid", "sum_outlet")
    tv, vo = HIO.column(case, "fluid", "volFlux_outlet")
    n = min(ts.size, tv.size)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_out = np.where(np.abs(vo[:n]) > 0, mo[:n] / vo[:n], np.nan)
    ports = meta["ports"]
    dz = ports["z_outlet_m"] - ports["z_inlet_centroid_m"]
    chim = (meta["rho_inlet"] - rho_out) * 9.81 * dz * vo[:n]
    mc = ((ts[:n] >= t0 - 1e-9) & (ts[:n] <= t1 + 1e-9) & np.isfinite(chim))
    w_port_chimney = trap(chim[mc], ts[:n][mc]) if mc.sum() >= 2 else math.nan
    closure = e_port - (phi - w_b + dke)
    ref = max(abs(phi), abs(w_b))
    wave = meta["waveform"]
    sine = wave == "sine"
    # D-P.7: the closure is GATED on developed isothermal STEADY windows only
    gated = bool(iso and developed and wave == "steady")
    if gated:
        gate_note = f"gated: isothermal steady window ({stage}), developed"
    elif not developed:
        gate_note = "reported, not gated: no developed window"
    elif not iso:
        gate_note = ("reported, not gated: heated run (the residual "
                     "measures the discretisation of the buoyant flow, not "
                     "the port-work instrument)")
    elif sine:
        gate_note = ("reported, not gated: sine wave (D-P.7 (3): the closure "
                     "against the cell-centred Phi is reported until the "
                     "registered re-examination; its necessary checks - "
                     "phase-bin sign, stationarity - are gated)")
    else:
        gate_note = ("reported, not gated: square wave (the valve edges' "
                     "temporal dissipation, ~KE_exchanged / n_edge^2, "
                     "dominates the residual; D-P.7 (5) splits it by valve "
                     "segment)")
    # D-P.7: the per-interval split (phase bins, sub-windows, valve segments)
    split = dict(sine_checks=None, square_split=None, r_per_period=None,
                 split_note="")
    P = meta.get("period_s")
    if developed and P:
        try:
            tph, phis = HIO.column(case, "fluid", "Phi_fluid")
            twb, wbs = HIO.column(case, "fluid", "Wb_fluid")
            if not (tph.size == tk.size == twb.size and np.allclose(tph, tk)
                    and np.allclose(twb, tk)):
                raise ValueError("KE, Phi and W_b are not on one time base")
            iv = interval_closure(t, pw, tk, ke, phis, wbs)
            split["r_per_period"] = closure / ((t1 - t0) / P)
            if sine:
                split["sine_checks"] = sine_closure_checks(iv, t0, t1, meta,
                                                           isothermal=iso)
            elif square:
                split["square_split"] = square_closure_split(iv, t0, t1, meta,
                                                             ke_tv)
        except ValueError as e:
            if sine and iso:
                # the necessary checks of a gated class are never skipped
                raise
            split["split_note"] = f"closure split not available: {e}"
    edge_rest = closure - e_edge if np.isfinite(e_edge) else math.nan
    # D-P.9: the closure against the face-based Phi, reported beside the
    # registered cell-centred one (None without the PhiFace monitor)
    try:
        face = face_closure(case, meta, t, pw, tk, ke, t0, t1, e_port, w_b, dke,
                            phi, closure, developed)
    except PHI_FACE_ERRORS as e:
        face = dict(error=f"{type(e).__name__}: {e}", window=[float(t0), float(t1)],
                    note=("the face-based monitor failed its data or its own checks "
                          "(D-P.9 (3)): no face-based number is reported; the "
                          "registered cell-centred closure is unaffected"))
    return dict(window=[float(t0), float(t1)], stage=stage,
                waveform=wave, period_s=P,
                developed=bool(developed), window_note=why, gated=gated,
                gate_note=gate_note, **split,
                E_port=e_port, W_b=w_b,
                W_b_port_estimate=w_port_chimney, Phi=phi, dKE=dke,
                KE_exchanged=ke_tv, valve_edge_steps=n_edge,
                E_edge_model=e_edge,
                closure_less_edge_model=edge_rest,
                closure_less_edge_model_rel_phi=(edge_rest / abs(phi)
                                                 if phi else math.nan),
                E_forced_measured=e_port + w_b, closure_residual=closure,
                closure_rel=closure / ref if ref else math.nan,
                phi_face=face,
                negative_E_port=bool(e_port < 0),
                note=("negative: chimney-driven flow, never a pump saving"
                      if e_port < 0 else ""))


def closurediag_report(cases, stage="run"):
    """D-P.9: the closure of each case against both Phi forms over its
    closure window (measured(); each case needs the face-based monitor), with
    W_bound (the viscous traction's work at the ports) beside r_face, and for
    every sine case the excess of its r/Phi over a steady case's in the
    list, in proportion to Phi (D-P.7 (6)'s second convention), in both
    forms. Reported, never gated. Returns a list of dict rows; a case whose
    face-based monitor fails its data or own checks gives a row with `error`
    and its registered cell-centred numbers only (the CLI then exits 1)."""
    rows = []
    for c in cases:
        m = measured(c, stage)
        f = m.get("phi_face")
        if f is None:
            raise ValueError(f"{c}: no face-based Phi monitor (PHI_FACE=1)")
        meta = HIO.case_meta(c)
        row = dict(case=str(c), waveform=m["waveform"],
                   div_u=meta.get("div_u", "scheme"),
                   scheme=meta.get("scheme"), window=m["window"],
                   developed=m["developed"], E_port=m["E_port"],
                   W_b=m["W_b"], dKE=m["dKE"], Phi_cell=m["Phi"],
                   r_over_phi_cell=_rel(m["closure_residual"], m["Phi"]),
                   error=f.get("error"))
        if not row["error"]:
            row.update({k: f[k] for k in (
                "Phi_face", "face_over_cell", "wall_share", "duct_share_face",
                "duct_share_cell", "r_over_phi_face", "identity_max_rel", "W_bound",
                "W_bound_over_phi_face", "r_with_W_bound_over_phi_face")})
        rows.append(row)
    steady = [r for r in rows if r["waveform"] == "steady"]
    for r in rows:
        if r["waveform"] == "sine" and steady:
            s0 = steady[0]
            r["excess_over_steady_cell"] = r["r_over_phi_cell"] - s0["r_over_phi_cell"]
            if not (r["error"] or s0["error"]):
                r["excess_over_steady_face"] = r["r_over_phi_face"] - s0["r_over_phi_face"]
            r["steady_case"] = s0["case"]
    return rows


def isothermal_dp(case):
    """(t, Q, dp_eff = P_port / Q) of an isothermal impedance run over its
    measurement window [t_pre, t_end] (whole run if no pre-circulation)."""
    meta = HIO.case_meta(case)
    if meta["heat"]["kind"] != "none":
        raise ValueError(f"{case}: not an isothermal (HEAT_MODEL none) run")
    t, pw, q = port_power(case, meta)
    t0 = max(meta["t_pre"], meta["ramp_s"])
    m = t >= t0 - 1e-9
    if np.any(q[m] <= 0):
        raise ValueError(f"{case}: Q <= 0 in the window; dp = P/Q undefined")
    return t[m], q[m], pw[m] / q[m], meta


def fit_cases(steady_cases, osc_cases):
    """The impedance fit of the Stage-B impedance block (D-H.1), with the
    D-P.3 fit uncertainty (cov_steady from >= 3 steady flows, R_se from the
    period-to-period scatter of each oscillatory run)."""
    qs, dps = [], []
    for c in steady_cases:
        t, q, dp, _ = isothermal_dp(c)
        qs.append(float(np.mean(q)))
        dps.append(float(np.mean(dp)))
    # D-P.8: the sign guard (R0 > 0; K < 0 refused beyond 2 SE, clamped
    # to 0 within 2 SE with R0 refitted; the unclamped K reported)
    g = fit_steady_guarded(qs, dps)
    imp = Impedance(R0=g["R0"], K=g["K"],
                    source=dict(steady=[str(c) for c in steady_cases],
                                steady_q=qs, steady_dp=dps,
                                osc=[str(c) for c in osc_cases],
                                K_fit=g["K_fit"], K_se=g["K_se"],
                                K_clamped=g["K_clamped"],
                                K_roundoff=g["K_roundoff"],
                                # D-P.7 (4): no Q_eq from a fit with oscillatory
                                # runs before the sine closure re-examination
                                # reports (q_eq_registered); filled with that
                                # report and its dated row
                                dp7_reexamination=None),
                    cov_steady=g["cov"] or [])
    fr, R, L, qb, se = [], [], [], [], []
    for c in osc_cases:
        t, q, dp, meta = isothermal_dp(c)
        f = float(meta["freq_hz"])
        if any(abs(f - x) < 1e-12 for x in fr):
            raise ValueError(f"impedance fit: two runs at {f:g} Hz; the 2 Hz run at "
                             "dt*/2 is the D-P.5 damping sensitivity (point_sensitivity), "
                             "not a fit point")
        Z, qbar, _ = fit_oscillation(t, q, dp, f)
        fr.append(f)
        R.append(Z.real)
        L.append(Z.imag / (2 * np.pi * f))
        qb.append(qbar)
        se.append(oscillation_se(t, q, dp, f))
    imp.freqs, imp.R, imp.L = fr, R, L
    imp.R_se = [] if any(x is None for x in se) else se
    imp.qbar_fit = float(np.mean(qb)) if qb else 0.0
    imp.source["R_below_qs_Hz"] = imp.below_qs()          # D-P.5: clamped points
    return imp


def point_sensitivity(imp, case):
    """D-P.5: the impedance with one fitted frequency's (R, L) replaced by
    those of another isothermal oscillatory run at the same frequency (the
    2 Hz sine at dt*/2): the backward scheme's own damping adds an in-phase
    resistance of order L w theta^3 / 4 at the arms' dt*, which the fit at
    dt* keeps (conservative: it can only raise a pulsed arm's E_diss); the
    E_diss and Q_eq recomputed with this impedance are reported beside the
    registered ones. Returns (impedance, {f, R_fit, R_case, dR_rel})."""
    t, q, dp, meta = isothermal_dp(case)
    f = float(meta["freq_hz"])
    idx = [i for i, x in enumerate(imp.freqs) if abs(float(x) - f) < 1e-12]
    if not idx:
        raise ValueError(f"point_sensitivity: {f:g} Hz is not a fitted frequency")
    Z, _, _ = fit_oscillation(t, q, dp, f)
    i = idx[0]
    R, L = list(imp.R), list(imp.L)
    r_fit = float(R[i])
    R[i], L[i] = float(Z.real), float(Z.imag / (2 * np.pi * f))
    alt = Impedance(R0=imp.R0, K=imp.K, qbar_fit=imp.qbar_fit, freqs=list(imp.freqs), R=R, L=L,
                    source=dict(imp.source, point_sensitivity=str(case)),
                    cov_steady=imp.cov_steady, R_se=imp.R_se)
    return alt, dict(f=f, R_fit=r_fit, R_case=R[i], dR_rel=(R[i] - r_fit) / imp.R0)


def arm_energy(imp, meta, t_total=None, check_case=None):
    """D-H.1 / D-P.3: the registered E_diss of an arm as run (its caseMeta:
    the waveform with its valve edges) over t_total (default the case's
    discharge): the conservative evaluation, replaced by the measured edged
    work of its isothermal model-check run (``check_case``, the impedance
    block's C1 / C2' run) when the model check fails. Returns the
    conservative dict plus E_arm (the value Q_eq inverts) and model_check."""
    qf, P, q_on, n_open = waveform_from_meta(meta)
    t_total = float(meta["t_dis"]) if t_total is None else float(t_total)
    r = primary_energy_conservative(imp, qf, P, t_total, q_on=q_on, n_open=n_open)
    # D-P.8: the arm's E_diss is never below steady flow's at its own mean
    r["jensen"] = jensen_floor(imp, qf, P, t_total, e_diss=r["E_diss"])
    r["E_arm"], r["model_check"] = r["E_diss"], None
    if check_case is not None:
        m = measured(check_case)
        cmeta = HIO.case_meta(check_case)
        cqf, cP, _, _ = waveform_from_meta(cmeta)
        w0, w1 = m["window"]
        e_mod = primary_energy_conservative(imp, cqf, cP, w1 - w0)["E_diss"]
        e_edge = m["E_edge_model"] if math.isfinite(m["E_edge_model"]) else 0.0
        mc = model_check(e_mod, m["E_port"], e_edge)
        r["model_check"] = mc
        if not mc["passed"]:
            # the measured edged work per unit time, over the discharge
            r["E_arm"] = max(r["E_diss"], mc["e_measured"] * t_total / (w1 - w0))
    return r


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--steady", nargs="+", type=Path, required=True)
    f.add_argument("--osc", nargs="*", type=Path, default=[])
    f.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("primary")
    p.add_argument("--impedance", type=Path, required=True)
    p.add_argument("--case", type=Path)
    p.add_argument("--waveform", choices=("steady", "sine", "square"))
    p.add_argument("--qbar", type=float, help="mL/min")
    p.add_argument("--freq", type=float)
    p.add_argument("--amp", type=float)
    p.add_argument("--duty", type=float)
    p.add_argument("--edge-ms", type=float, default=None)
    p.add_argument("--t", type=float, default=900.0)
    qe = sub.add_parser("qeq")
    qe.add_argument("--impedance", type=Path, required=True)
    qe.add_argument("--energy", type=float, required=True)
    qe.add_argument("--t", type=float, default=900.0)
    m = sub.add_parser("measured")
    m.add_argument("case", type=Path)
    m.add_argument("--stage", choices=("run", "pre"), default="run",
                   help="pre: the isothermal pre-circulation [ramp, t_pre]")
    cd = sub.add_parser("closurediag", help="D-P.9: the closure against the "
                        "cell-centred and the face-based Phi (reported)")
    cd.add_argument("cases", type=Path, nargs="+")
    cd.add_argument("--stage", choices=("run", "pre"), default="run")
    cd.add_argument("--json", type=Path)
    sd = sub.add_parser("sine-dt", help="D-P.7 (3c): dt invariance of a sine "
                        "class's isothermal closure, dt rung vs its dt/2 partner")
    sd.add_argument("case_dt", type=Path)
    sd.add_argument("case_half", type=Path)
    sd.add_argument("--stage", choices=("run", "pre"), default="pre")
    a = ap.parse_args(argv)
    if a.cmd == "fit":
        imp = fit_cases(a.steady, a.osc)
        a.out.write_text(json.dumps(imp.to_json(), indent=1) + "\n")
        print(json.dumps(imp.to_json(), indent=1))
        return 0
    if a.cmd == "measured":
        print(json.dumps(measured(a.case, a.stage), indent=1, default=float))
        return 0
    if a.cmd == "closurediag":
        rows = closurediag_report(a.cases, a.stage)
        for r in rows:
            head = (f"{Path(r['case']).name} [{r['window'][0]:g}, {r['window'][1]:g}] s "
                    f"({'developed' if r['developed'] else 'NOT developed'}; div(phi,U) "
                    f"{r['div_u']}): E_port {r['E_port']:.4g} J, ")
            ex = ""
            if "excess_over_steady_cell" in r:
                ex = (f"; excess over steady (per Phi) cell "
                      f"{100 * r['excess_over_steady_cell']:+.2f}"
                      + (f" / face {100 * r['excess_over_steady_face']:+.2f}"
                         if "excess_over_steady_face" in r else "") + " points")
            if r["error"]:
                print(head + f"Phi cell {r['Phi_cell']:.4g} J, r/Phi cell "
                      f"{100 * r['r_over_phi_cell']:+.2f} %{ex}; face-based Phi: HARD FAILURE "
                      f"of the instrument, no face-based number: {r['error']}")
                continue
            print(head + f"Phi cell {r['Phi_cell']:.4g} J / face {r['Phi_face']:.4g} J "
                  f"(face/cell {r['face_over_cell']:.4f}, wall share "
                  f"{100 * r['wall_share']:.1f} %), r/Phi cell "
                  f"{100 * r['r_over_phi_cell']:+.2f} % / face "
                  f"{100 * r['r_over_phi_face']:+.2f} %; W_bound {r['W_bound']:+.4g} J "
                  f"({100 * r['W_bound_over_phi_face']:+.2f} % of Phi_face; r_face + W_bound "
                  f"{100 * r['r_with_W_bound_over_phi_face']:+.2f} % of it){ex}")
        if a.json:
            a.json.write_text(json.dumps(rows, indent=1, default=float) + "\n")
        return 1 if any(r["error"] for r in rows) else 0
    if a.cmd == "sine-dt":
        gate, rep = sine_dt_invariance(a.case_dt, a.case_half, a.stage)
        print(("PASS" if gate.ok else "FAIL") + f" closure_dt_invariance: {gate.detail}")
        print(json.dumps(rep, indent=1, default=float))
        return 0 if gate.ok else 1
    imp = Impedance.from_json(json.loads(a.impedance.read_text()))
    if a.cmd == "qeq":
        # D-P.7 (4): a Stage-B fit sets no Q_eq before the re-examination
        print(f"Q_eq = {q_eq_registered(imp, a.energy, a.t) / ML_MIN:.6g} mL/min")
        return 0
    if a.case:
        meta = HIO.case_meta(a.case)
        qf, P, q_on, n_open = waveform_from_meta(meta)
        t_total = meta["t_dis"]
    else:
        edge = a.edge_ms * 1e-3 if a.edge_ms is not None else (
            min(0.05, a.duty / a.freq / 5) if a.waveform == "square" else None)
        qf, P, q_on, n_open = waveform(a.waveform, a.qbar * ML_MIN, a.freq,
                                       a.amp, a.duty, edge)
        t_total = a.t
    # the registered evaluation (D-P.3): the larger of clamp and sqrt
    r = primary_energy_conservative(imp, qf, P, t_total, q_on=q_on, n_open=n_open)
    print(json.dumps(r, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
