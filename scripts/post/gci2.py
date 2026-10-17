#!/usr/bin/env python3
"""Grid / time-step convergence v2 (WP 2.4; plan D-C, NEW-2; register R-01(j), R-23).

Replaces gci.py (3 grids, absolute-degC normalisation, no exit code) for the
new campaign. Works for any refinement parameter h: cell size (V/N)^(1/3)
for grid families, the time step for Delta-t ladders.

Methods
- Celik et al. (2008, ASME J. Fluids Eng. 130:078001) on a grid triplet, with
  the convergence class from R = eps21/eps32 (Eca & Hoekstra 2014 Eq. 2):
  monotonic 0<R<1 -> GCI; oscillatory -1<R<0 -> Stern et al. (2001) band
  U = (S_U - S_L)/2 over the triplet; divergent (R >= 1 or R <= -1) -> U = inf.
  Pre-registered order/safety-factor rule for the triplet, aligned with the
  Eca-Hoekstra rule below so that 3- and 4-level studies cannot disagree by
  an order of magnitude: the order used is min(p_obs, p_formal) and
  Fs = 1.25 iff 0.5 <= p_obs < p_formal + 0.1 (2.1 at p_formal = 2), else 3.
  eps21 = 0 exactly (fine = medium) with no noise floor and no print
  resolution is "degenerate": U = inf (the GCI would otherwise be 0 with
  p = inf). Give ``se`` or ``resolution`` to get a finite band.
- Eca & Hoekstra (2014, J. Comput. Phys. 262:104-130, Appendix A) least
  squares for >= 4 levels: delta_RE = alpha h^p (weighted and unweighted);
  0.5 <= p <= 2 -> delta_RE; p > 2 -> best of delta_1, delta_2; p < 0.5 or
  no p > 0 -> best of delta_1, delta_2, delta_12 (all weighted and
  unweighted; "best" = smallest standard deviation). Data range parameter
  Dphi = (phi_max - phi_min)/(n_g - 1); Fs = 1.25 iff 0.5 <= p < 2.1 and
  sigma < Dphi, else 3; U = Fs eps + sigma + |phi_i - phi_fit| if
  sigma < Dphi, else U = 3 (sigma/Dphi)(eps + sigma + |phi_i - phi_fit|).

Project rules layered on top
- Every relative error is normalised by the supplied temperature RISE, never
  by absolute degC (the v1 gci.py:39-40 defect).
- Noise floor (D-C, R-01(j)): if a level is an ensemble mean with standard
  error se_i, the error used for level i is max(|eps_i|, 2 SE_pair) with
  SE_pair = sqrt(se_i^2 + se_partner^2) (partner = next finer level; the
  finest level pairs with the second finest). ``resolution`` (the print
  quantum of phi, e.g. 1e-3 K for 6 significant digits at ~300 K) is a
  second floor: a level difference <= resolution is indistinguishable from
  zero. In the Celik branch, a triplet whose |eps21| or |eps32| is below its
  floor is "noise": the observed order is not used, p = p_formal and U uses
  fs_noise (default 1.25, the plan's tost_a7 u_num model: Fs * F * |level
  difference| with F = r^p/(r^p - 1)).
- Noise-aware gate threshold (D-C for grids with F = 2.2, NEW-2 for dt
  ladders with F = 4/3): U <= max(5 mK, 2 F Fs SE_pair) — the CLI's
  --noise-aware mode, for Delta-T_mid and pair-averaged Delta-T.
- Standard uncertainty u = U / k_num. ASME V&V 20-2009 recommends k = 2 for a
  Gaussian error distribution centred on the fine-grid value and 1.1-1.15 for
  a distribution centred on the extrapolated value (shifted Gaussian); the
  default here is the conservative 1.15, which the pre-registration fixes
  (PROJECT.md §2.5 D-F.3, docs/prereg/prereg.yaml; tests/prereg).

CLI (exit 0 = pass, 1 = D-C gate failed, 2 = input error):
    gci2.py --phi 31.81 31.28 31.23 31.21 --cells 218e3 468e3 1048e3 2.3e6 \
            --rise 6.5 [--se ...] [--resolution 1e-3] --level 1 \
            --max-U-K 0.07 --max-U-pct-rise 1.0
    gci2.py --phi ... --h ... --rise 6.5 --se ... --level 1 --noise-aware \
            [--F 2.2] [--fs-noise 1.25] [--noise-floor-K 0.005]
The level index refers to the input order. Output: JSON (stdout or --out).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass

import numpy as np
from scipy import optimize

K_NUM = 1.15            # V&V 20 conversion of a GCI-type band to a standard uncertainty
FS_GOOD, FS_BAD = 1.25, 3.0
P_LOW = 0.5             # below: Fs = 3 (Eca & Hoekstra 2014)
P_FS_MARGIN = 0.1       # Fs = 1.25 up to p_formal + 0.1 (2.1 at p_formal = 2)
NOISE_FLOOR_K = 0.005   # D-C / NEW-2 absolute floor of the noise-aware threshold
P_SEARCH = (-8.0, 8.0)  # observed-order search interval for delta_RE


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def cell_size(n_cells, volume: float = 1.0):
    """Representative cell size h = (V / N)^(1/3) (Celik 2008 Eq. 1)."""
    n = np.asarray(n_cells, dtype=float)
    if np.any(n <= 0) or volume <= 0:
        raise ValueError("cell counts and volume must be > 0")
    return (volume / n) ** (1.0 / 3.0)


def richardson_factor(r: float, p: float = 2.0) -> float:
    """F = r^p / (r^p - 1): error of the COARSER level of a pair in units of
    the level difference (2.216 for r = 1.35, p = 2; 4/3 for r = 2, p = 2)."""
    rp = r ** p
    if rp <= 1.0:
        raise ValueError("need r^p > 1")
    return rp / (rp - 1.0)


def noise_threshold(se_diff: float, F: float, fs: float = FS_GOOD,
                    floor: float = NOISE_FLOOR_K) -> float:
    """Noise-aware gate threshold of D-C / NEW-2: max(floor, 2 F Fs SE)."""
    return max(floor, 2.0 * F * fs * se_diff)


def _prep(phi, h, rise, se):
    phi = np.asarray(phi, dtype=float).ravel()
    h = np.asarray(h, dtype=float).ravel()
    if phi.size != h.size:
        raise ValueError("phi and h must have the same length")
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(h))):
        raise ValueError("phi and h must be finite")
    if np.any(h <= 0):
        raise ValueError("h must be > 0")
    if np.unique(h).size != h.size:
        raise ValueError("h values must be distinct")
    if not (rise is not None and math.isfinite(rise) and rise > 0):
        raise ValueError("a positive temperature rise is required for normalisation")
    if se is None:
        se = np.zeros_like(phi)
    else:
        se = np.asarray(se, dtype=float).ravel()
        if se.size != phi.size or np.any(~np.isfinite(se)) or np.any(se < 0):
            raise ValueError("se must be finite, >= 0 and one per level")
    order = np.argsort(h)                     # finest first
    return phi[order], h[order], se[order], order


def _partner(i: int) -> int:
    """Sorted index of the level a level is differenced with (next finer;
    the finest pairs with the second finest)."""
    return i - 1 if i > 0 else 1


def _pair_floors(se_sorted, resolution: float = 0.0):
    """max(2 * SE of each level's difference with its partner, resolution)."""
    n = se_sorted.size
    out = np.empty(n)
    for i in range(n):
        out[i] = max(2.0 * math.hypot(se_sorted[i], se_sorted[_partner(i)]),
                     resolution)
    return out


def _check_resolution(resolution) -> float:
    q = float(resolution)
    if not (math.isfinite(q) and q >= 0):
        raise ValueError("resolution must be finite and >= 0")
    return q


# --------------------------------------------------------------------------
# Celik triplet
# --------------------------------------------------------------------------

@dataclass
class CelikResult:
    h: tuple                 # fine, medium, coarse
    phi: tuple
    rise: float
    r21: float
    r32: float
    eps21: float             # phi_medium - phi_fine
    eps32: float             # phi_coarse - phi_medium
    R: float                 # eps21 / eps32
    convergence: str         # monotonic | oscillatory | monotonic-divergent |
                             # oscillatory-divergent | noise | degenerate
    p: float                 # order used for the band
    p_observed: float        # Celik fixed-point order (nan if not defined)
    phi_ext: float
    floor21: float           # max(2 SE_pair, resolution)
    floor32: float
    resolution: float
    fs: float
    U_fine: float            # expanded band, K
    U_medium: float
    U_fine_pct_rise: float
    U_medium_pct_rise: float
    u_fine: float            # standard, K
    u_medium: float
    ea21_pct_rise: float     # |eps21| / rise
    eext21_pct_rise: float   # |phi_ext - phi_fine| / rise
    k_num: float
    notes: tuple


def _celik_order(e21, e32, r21, r32, s):
    a = abs(e32 / e21)
    lr = math.log(r21)
    p = max(abs(math.log(a)) / lr, 1e-9)
    for _ in range(2000):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = max(abs(math.log(a) + q) / lr, 1e-9)
        if abs(p_new - p) < 1e-12:
            return p_new, True
        p = 0.5 * (p + p_new)                  # damped fixed point
    return p, False


def celik_gci(phi, h, rise: float, se=None, fs: float = FS_GOOD,
              p_formal: float = 2.0, fs_noise: float = FS_GOOD,
              k_num: float = K_NUM, resolution: float = 0.0) -> CelikResult:
    """Celik (2008) GCI on exactly three levels (any input order).

    ``resolution`` is the print quantum of phi (K); level differences at or
    below it are treated as zero (noise branch). See the module docstring for
    the order / safety-factor rule and the degenerate eps21 = 0 case.
    """
    phi_s, h_s, se_s, _ = _prep(phi, h, rise, se)
    if phi_s.size != 3:
        raise ValueError("celik_gci needs exactly three levels")
    q = _check_resolution(resolution)
    p1, p2, p3 = phi_s
    h1, h2, h3 = h_s
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = p2 - p1, p3 - p2
    s21 = 2.0 * math.hypot(se_s[0], se_s[1])
    s32 = 2.0 * math.hypot(se_s[1], se_s[2])
    f21, f32 = max(s21, q), max(s32, q)
    notes = []
    R = e21 / e32 if e32 != 0 else (math.inf if e21 != 0 else math.nan)
    p_obs = math.nan
    phi_ext = math.nan
    fs_used = fs
    q_tol = q + 1e-9 * max(1.0, float(np.max(np.abs(phi_s))))   # float slack of phi_i - phi_j
    below = lambda e, s_floor: abs(e) < s_floor or (q > 0 and abs(e) <= q_tol)
    noisy = below(e21, s21) or below(e32, s32)

    if noisy:
        conv, p = "noise", p_formal
        fs_used = fs_noise
        e = max(abs(e21), f21)
        rp = r21 ** p
        U_f, U_m = fs_used * e / (rp - 1.0), fs_used * e * rp / (rp - 1.0)
        phi_ext = p1 - e21 / (rp - 1.0)
        notes.append("increment at or below the noise/resolution floor: formal order used")
        if below(e21, s21) and not below(e32, s32):
            notes.append("coarse increment resolved but fine increment at the floor: "
                         "check the asymptotic range (fourth level or more members)")
    elif abs(e21) <= q_tol:
        # fine = medium (to float slack) and no floor: Celik would give p = inf, U = 0
        conv, p, U_f, U_m = "degenerate", math.nan, math.inf, math.inf
        notes.append("eps21 = 0 with no se and no resolution: GCI undefined; "
                     "supply se (ensemble levels) or the print resolution")
    elif not math.isfinite(R) or R >= 1.0:
        conv, p, U_f, U_m = "monotonic-divergent", math.nan, math.inf, math.inf
    elif R <= -1.0:
        conv, p, U_f, U_m = "oscillatory-divergent", math.nan, math.inf, math.inf
    elif R < 0.0:
        conv = "oscillatory"
        p_obs, ok = _celik_order(e21, e32, r21, r32, -1.0)
        p = p_obs
        U_f = U_m = 0.5 * (max(phi_s) - min(phi_s))
        notes.append("oscillatory convergence: Stern band (S_U - S_L)/2")
    else:
        conv = "monotonic"
        p_obs, ok = _celik_order(e21, e32, r21, r32, 1.0)
        if not ok:
            notes.append("Celik fixed point did not converge")
        if not (P_LOW <= p_obs < p_formal + P_FS_MARGIN):
            fs_used = FS_BAD
            notes.append(f"observed order {p_obs:.3g} outside [{P_LOW}, "
                         f"{p_formal + P_FS_MARGIN:.3g}): Fs = {FS_BAD:g}")
        p = min(p_obs, p_formal)
        if p_obs > p_formal:
            notes.append(f"observed order {p_obs:.3g} > formal: formal order "
                         f"{p_formal:g} used (Eca-Hoekstra delta_{p_formal:g} analogue)")
        rp = r21 ** p
        phi_ext = (rp * p1 - p2) / (rp - 1.0)
        U_f = fs_used * abs(e21) / (rp - 1.0)
        U_m = fs_used * abs(e21) * rp / (rp - 1.0)
    if q == 0.0 and not np.any(se_s > 0):
        notes.append("no se and no resolution supplied: no noise floor applied")
    pct = lambda v: 100.0 * v / rise
    return CelikResult(
        tuple(h_s), tuple(phi_s), rise, r21, r32, e21, e32, R, conv, p, p_obs,
        phi_ext, f21, f32, q, fs_used, U_f, U_m, pct(U_f), pct(U_m),
        U_f / k_num, U_m / k_num, pct(abs(e21)),
        pct(abs(phi_ext - p1)) if math.isfinite(phi_ext) else math.nan,
        k_num, tuple(notes))


# --------------------------------------------------------------------------
# Eca & Hoekstra (2014) least squares
# --------------------------------------------------------------------------

@dataclass
class Fit:
    name: str                # RE | 1 | 2 | 12
    weighted: bool
    phi0: float
    coef: tuple              # alpha (RE/1/2) or (alpha1, alpha2) (12)
    p: float                 # RE order; 1/2 fixed; nan for 12
    sigma: float
    err: np.ndarray          # signed fitted error phi_fit(h_i) - phi0
    resid: np.ndarray        # phi_i - phi_fit(h_i)


def _weights(h, weighted):
    if not weighted:
        return np.ones_like(h), np.ones_like(h)
    w = (1.0 / h) / np.sum(1.0 / h)
    return w, h.size * w


def _linear_fit(hn, phi, basis, weighted, npar, name, p):
    w, nw = _weights(hn, weighted)
    X = np.column_stack([np.ones_like(hn)] + basis)
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], phi * sw, rcond=None)
    fit = X @ beta
    resid = phi - fit
    sigma = math.sqrt(float(np.sum(nw * resid ** 2)) / (hn.size - npar))
    return Fit(name, weighted, float(beta[0]), tuple(float(b) for b in beta[1:]),
               p, sigma, fit - beta[0], resid)


def _re_fit(hn, phi, weighted):
    w, _ = _weights(hn, weighted)
    sw = np.sqrt(w)

    def ssr(p):
        X = np.column_stack([np.ones_like(hn), hn ** p])
        beta, *_ = np.linalg.lstsq(X * sw[:, None], phi * sw, rcond=None)
        return float(np.sum(w * (phi - X @ beta) ** 2))

    grid = np.linspace(*P_SEARCH, 3201)
    grid = grid[np.abs(grid) > 1e-3]
    vals = np.array([ssr(p) for p in grid])
    i = int(np.argmin(vals))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, grid.size - 1)]
    res = optimize.minimize_scalar(ssr, bounds=(lo, hi), method="bounded",
                                   options={"xatol": 1e-10})
    p = float(res.x) if res.fun <= vals[i] else float(grid[i])
    fit = _linear_fit(hn, phi, [hn ** p], weighted, 3, "RE", p)
    return fit


@dataclass
class EcaHoekstraResult:
    h: tuple                 # input order
    phi: tuple
    rise: float
    p: float                 # observed order of the selected RE fit (nan: anomalous)
    estimator: str           # RE | 1 | 2 | 12
    weighted: bool
    phi0: float
    sigma: float
    delta_phi: float         # data range parameter
    good_fit: bool           # sigma < delta_phi
    fs: float
    eps: tuple               # |fitted error| per level, K (before the floor)
    floor: tuple             # max(2 * SE_pair, resolution) per level
    noise_levels: tuple      # input indices where the floor binds
    U: tuple                 # expanded band per level, K
    U_pct_rise: tuple
    u: tuple                 # standard uncertainty per level, K
    k_num: float
    fits: tuple              # (name, weighted, p, sigma) of every fit tried
    notes: tuple


def eca_hoekstra(phi, h, rise: float, se=None, k_num: float = K_NUM,
                 resolution: float = 0.0) -> EcaHoekstraResult:
    """Eca & Hoekstra (2014) uncertainty for every level of a >= 4-level study."""
    phi_s, h_s, se_s, order = _prep(phi, h, rise, se)
    q = _check_resolution(resolution)
    ng = phi_s.size
    if ng < 4:
        raise ValueError("Eca-Hoekstra least squares needs >= 4 levels")
    hn = h_s / h_s[0]
    notes = []
    tried = []

    re_fits = [_re_fit(hn, phi_s, wgt) for wgt in (False, True)]
    tried += re_fits
    positive = [f for f in re_fits if f.p > 0]
    in_range = [f for f in positive if 0.5 <= f.p <= 2.0]
    fixed = lambda names: [_linear_fit(hn, phi_s, b, wgt, npar, nm, pp)
                           for wgt in (False, True)
                           for nm, b, npar, pp in (
                               ("1", [hn], 2, 1.0), ("2", [hn ** 2], 2, 2.0),
                               ("12", [hn, hn ** 2], 3, math.nan))
                           if nm in names]
    if in_range:
        sel = min(in_range, key=lambda f: f.sigma)
        p_obs = sel.p
    elif positive:
        best = min(positive, key=lambda f: f.sigma)
        p_obs = best.p
        alt = fixed(("1", "2")) if p_obs > 2.0 else fixed(("1", "2", "12"))
        tried += alt
        sel = min(alt, key=lambda f: f.sigma)
        notes.append(f"observed order {p_obs:.3g} outside [0.5, 2]: "
                     f"fixed-exponent estimator delta_{sel.name}")
    else:
        p_obs = math.nan
        alt = fixed(("1", "2", "12"))
        tried += alt
        sel = min(alt, key=lambda f: f.sigma)
        notes.append("anomalous behaviour (no fit with p > 0)")
    if any(abs(f.p) >= P_SEARCH[1] - 1e-6 for f in re_fits):
        notes.append("an RE fit reached the order search bound")

    delta_phi = (phi_s.max() - phi_s.min()) / (ng - 1)
    eps = np.abs(sel.err)
    floor = _pair_floors(se_s, q)
    eps_eff = np.maximum(eps, floor)
    absres = np.abs(sel.resid)
    if delta_phi == 0.0:
        good, fs = True, FS_BAD
        if np.all(floor == 0.0):
            U = np.full(ng, math.inf)
            notes.append("degenerate: all levels identical and no se/resolution "
                         "floor; supply se or the print resolution")
        else:
            U = FS_BAD * floor
            notes.append("flat data (all levels identical)")
    else:
        good = sel.sigma < delta_phi
        fs = FS_GOOD if (math.isfinite(p_obs) and 0.5 <= p_obs < 2.1 and good) else FS_BAD
        if good:
            U = fs * eps_eff + sel.sigma + absres
        else:
            U = 3.0 * (sel.sigma / delta_phi) * (eps_eff + sel.sigma + absres)
            notes.append("sigma >= data range: scattered-data band")

    inv = np.empty(ng, dtype=int)
    inv[order] = np.arange(ng)                   # sorted index of each input level
    back = lambda a: tuple(float(a[inv[i]]) for i in range(ng))
    noise_idx = tuple(int(i) for i in range(ng) if eps[inv[i]] < floor[inv[i]])
    return EcaHoekstraResult(
        back(h_s), back(phi_s), rise, p_obs, sel.name, sel.weighted, sel.phi0,
        sel.sigma, delta_phi, bool(good), fs, back(eps), back(floor), noise_idx,
        back(U), back(100.0 * U / rise), back(U / k_num), k_num,
        tuple((f.name, f.weighted, f.p, f.sigma) for f in tried), tuple(notes))


# --------------------------------------------------------------------------
# dispatcher
# --------------------------------------------------------------------------

def grid_uncertainty(phi, h, rise: float, se=None, method: str = "auto",
                     k_num: float = K_NUM, p_formal: float = 2.0,
                     resolution: float = 0.0) -> dict:
    """Per-level expanded (U) and standard (u) uncertainty, in the input order.

    method 'auto': Eca-Hoekstra for >= 4 levels, Celik for 3. Every
    consecutive triplet is also classified with Celik (oscillation check).
    """
    phi_a = np.asarray(phi, dtype=float).ravel()
    h_a = np.asarray(h, dtype=float).ravel()
    ng = phi_a.size
    if ng < 3:
        raise ValueError("need >= 3 levels")
    if method == "auto":
        method = "eh" if ng >= 4 else "celik"
    phi_s, h_s, se_s, order = _prep(phi_a, h_a, rise, se)
    triplets = []
    for i in range(ng - 2):
        c = celik_gci(phi_s[i:i + 3], h_s[i:i + 3], rise, se_s[i:i + 3],
                      p_formal=p_formal, k_num=k_num, resolution=resolution)
        triplets.append(asdict(c))
    out = {"method": method, "n_levels": ng, "rise": rise, "k_num": k_num,
           "resolution": resolution, "triplets": triplets}
    if method == "eh":
        r = eca_hoekstra(phi_a, h_a, rise, se, k_num, resolution)
        out.update(result=asdict(r), U=list(r.U), u=list(r.u),
                   U_pct_rise=list(r.U_pct_rise))
    elif method == "celik":
        if ng != 3:
            raise ValueError("method 'celik' needs exactly 3 levels")
        c = celik_gci(phi_a, h_a, rise, se, p_formal=p_formal, k_num=k_num,
                      resolution=resolution)
        inv = np.empty(3, dtype=int)
        inv[order] = np.arange(3)
        sorted_U = [c.U_fine, c.U_medium, math.nan]   # the coarse level has no GCI
        U = [sorted_U[inv[i]] for i in range(3)]
        out.update(result=asdict(c), U=U, u=[v / k_num for v in U],
                   U_pct_rise=[100.0 * v / rise for v in U])
    else:
        raise ValueError("method must be auto, eh or celik")
    return out


def level_noise_threshold(h, se, level: int, p_formal: float = 2.0,
                          F: float | None = None, fs: float = FS_GOOD,
                          floor: float = NOISE_FLOOR_K) -> dict:
    """Noise-aware threshold (D-C / NEW-2) for the input-order ``level``.

    SE_pair = sqrt(se_level^2 + se_partner^2) with the partner of
    _pair_floors. F defaults to the factor that maps a level difference to
    the error of the gated level at the formal order: r^p/(r^p - 1) for a
    coarser level of its pair (2.216 at r = 1.35, 4/3 at r = 2) and
    1/(r^p - 1) for the finest level. D-C fixes F = 2.2 for grids and
    NEW-2 F = 4/3 for dt ladders; pass F to use the plan's rounded value.
    """
    h = np.asarray(h, dtype=float).ravel()
    if se is None:
        raise ValueError("the noise-aware threshold needs se for every level")
    se = np.asarray(se, dtype=float).ravel()
    if se.size != h.size or np.any(~np.isfinite(se)) or np.any(se < 0):
        raise ValueError("se must be finite, >= 0 and one per level")
    if not 0 <= level < h.size:
        raise ValueError("level out of range")
    order = np.argsort(h)
    inv = np.empty(h.size, dtype=int)
    inv[order] = np.arange(h.size)
    i = int(inv[level])
    j = _partner(i)
    lvl, par = int(order[i]), int(order[j])
    se_pair = math.hypot(se[lvl], se[par])
    if F is None:
        r = h[order[max(i, j)]] / h[order[min(i, j)]]
        rp = r ** p_formal
        F = rp / (rp - 1.0) if i > 0 else 1.0 / (rp - 1.0)
    if not (math.isfinite(F) and F > 0):
        raise ValueError("F must be > 0")
    thr = noise_threshold(se_pair, F, fs, floor)
    return {"se_pair": se_pair, "F": F, "fs": fs, "floor": floor, "threshold": thr}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else str(v)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--phi", type=float, nargs="+", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--cells", type=float, nargs="+", help="cell counts per level")
    g.add_argument("--h", type=float, nargs="+", help="refinement parameter (e.g. dt)")
    ap.add_argument("--volume", type=float, default=1.0)
    ap.add_argument("--rise", type=float, required=True,
                    help="temperature rise used to normalise every %% (K)")
    ap.add_argument("--se", type=float, nargs="+", default=None,
                    help="standard error of each level (ensemble means)")
    ap.add_argument("--method", choices=("auto", "eh", "celik"), default="auto")
    ap.add_argument("--k-num", type=float, default=K_NUM)
    ap.add_argument("--p-formal", type=float, default=2.0)
    ap.add_argument("--level", type=int, default=None,
                    help="input index of the production level to gate")
    ap.add_argument("--resolution", type=float, default=0.0,
                    help="print quantum of phi in K (e.g. 1e-3 for 6 significant "
                         "digits near 300 K); differences <= it count as zero")
    thr = ap.add_mutually_exclusive_group()
    thr.add_argument("--max-U-K", type=float, default=None,
                     help="fixed threshold (D-C: 0.07 K on surface T_max)")
    thr.add_argument("--noise-aware", action="store_true",
                     help="threshold max(floor, 2 F Fs SE_pair) (D-C / NEW-2, "
                          "for Delta-T_mid and pair-averaged Delta-T); needs --se")
    ap.add_argument("--F", type=float, default=None,
                    help="noise-aware F (default from the level's refinement ratio "
                         "at --p-formal; D-C uses 2.2, NEW-2 4/3)")
    ap.add_argument("--fs-noise", type=float, default=FS_GOOD,
                    help="noise-aware Fs (default 1.25)")
    ap.add_argument("--noise-floor-K", type=float, default=NOISE_FLOOR_K,
                    help="noise-aware absolute floor (default 0.005 K)")
    ap.add_argument("--max-U-pct-rise", type=float, default=None)
    ap.add_argument("--out", default=None, help="write the JSON result here")
    a = ap.parse_args(argv)
    if a.noise_aware and a.level is None:
        print("gci2: input error: --noise-aware needs --level", file=sys.stderr)
        return 2
    try:
        h = cell_size(a.cells, a.volume) if a.cells is not None else a.h
        res = grid_uncertainty(a.phi, h, a.rise, a.se, a.method, a.k_num, a.p_formal,
                               a.resolution)
        noise = None
        if a.noise_aware:
            noise = level_noise_threshold(h, a.se, a.level, a.p_formal, a.F,
                                          a.fs_noise, a.noise_floor_K)
    except ValueError as e:
        print(f"gci2: input error: {e}", file=sys.stderr)
        return 2
    rc = 0
    if a.level is not None:
        if not 0 <= a.level < len(a.phi):
            print("gci2: --level out of range", file=sys.stderr)
            return 2
        U = res["U"][a.level]
        Up = res["U_pct_rise"][a.level]
        fails = []
        if not math.isfinite(U):
            fails.append(f"U undefined ({U})")
        if a.max_U_K is not None and not U <= a.max_U_K:
            fails.append(f"U {U:.4g} K > {a.max_U_K} K")
        if noise is not None and not U <= noise["threshold"]:
            fails.append(f"U {U:.4g} K > noise-aware threshold "
                         f"{noise['threshold']:.4g} K (max({noise['floor']:g}, "
                         f"2 x {noise['F']:.4g} x {noise['fs']:g} x SE_pair "
                         f"{noise['se_pair']:.4g}))")
        if a.max_U_pct_rise is not None and not Up <= a.max_U_pct_rise:
            fails.append(f"U {Up:.3g} % of rise > {a.max_U_pct_rise} %")
        res["gate"] = {"level": a.level, "U_K": U, "U_pct_rise": Up,
                       "threshold_K": noise["threshold"] if noise else a.max_U_K,
                       "noise_aware": noise, "ok": not fails, "fails": fails}
        rc = 1 if fails else 0
    text = json.dumps(_jsonable(res), indent=2)
    if a.out:
        with open(a.out, "w") as f:
            f.write(text + "\n")
    else:
        print(text)
    if a.level is not None:
        g = res["gate"]
        print(f"GCI2 {'PASS' if g['ok'] else 'FAIL'} level {a.level}: "
              f"U = {g['U_K']:.4g} K ({g['U_pct_rise']:.3g} % of rise)"
              + ("" if g["ok"] else " — " + "; ".join(g["fails"])), file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
