"""D-I near-wall targets and property-derived layer thicknesses (WP 2.5a).

One place for the numbers that k0_blockmesh.py (design report) and
k0_measure.py (measurement + --gate) both need, so the two cannot disagree.

Sources:
  * D-I rule text (MASTER-PLAN-full-scope.md §2, D-I "Rule"): >= 5 cells
    inside delta_T (band 0.37-0.43 mm), >= 3 cells inside the Stokes layer at
    the highest frequency run on that mesh ("0.18 mm at 4 Hz; 0.13 mm for
    conditional 8 Hz"), last-layer/core <= 2.0, non-orthogonality <= 65,
    skewness <= 4, checkMesh -allRegions clean; D-C: r_N >= 1.3 per region.
  * Stokes thickness delta_S = sqrt(2 nu / omega) = sqrt(nu / (pi f)) with
    nu = mu / rho of the coolant from data/reference/properties.yaml (the
    25 C anchors today; D-R will lower nu at the film temperature, so pass
    the film-temperature nu explicitly once D-R lands).
  * Frequencies: the highest CORE frequency is 2 Hz (C3, sine 2 Hz A0.6,
    RESUBMISSION-PLAN §2.2); 4 Hz is PLUS-4 only.
"""

from __future__ import annotations

import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROPS = REPO / "data" / "reference" / "properties.yaml"

FILM_BAND_MM = (0.37, 0.43)     # D-I delta_T band; cells are counted inside
FILM_MM = FILM_BAND_MM[0]       # the LOWER bound (the conservative reading)
FILM_MIN_CELLS = 5
STOKES_MIN_CELLS = 3
# D-I rule text quotes rounded thicknesses; the gate uses the thinner of the
# quoted value and sqrt(nu / (pi f)) at the nu in use
DI_STOKES_QUOTED_MM = {4.0: 0.18, 8.0: 0.13}
CORE_MAX_FREQ_HZ = 2.0          # C3; 4 Hz = PLUS-4 only
REPORT_FREQS_HZ = (2.0, 4.0)
NEARWALL_RATIO_MAX = 2.0        # D-I last-layer/adjacent-core analogue
NON_ORTHO_MAX = 65.0
SKEWNESS_MAX = 4.0
RN_MIN = 1.3                    # D-C: achieved r_N per region
# far-field expansion-ratio REGRESSION GUARD, not a quality target: the
# up/downstream y-spacing jump at y = +-20 mm is 4.8 at medium (K0-design
# §7 item 4, fix deferred to WP 2.5b, whose target is <= 2.0)
ER_GUARD_ALL_FLUID = 5.0
ER_MAX_TILES = 1.3              # module tiles (z < cap), every level

# ---- D-G housing (WP 2.5b) -------------------------------------------------
# Port regions (inlet channel, the two slot columns, riser) carry the whole
# through-flow; their expansion ratio is held to the D-I near-wall limit
# (the 2.5a far-field 4.8 jump is removed there by the far-line
# y-distribution). The U2/D2 H-blocks keep the all-fluid guard above: their
# module face still meets the 3 mm strip's 0.26 mm cells (K0-design §10).
ER_MAX_PORTS = 2.0
# all-fluid guard of the housing: the module-end strip lines (strip_Cend)
# cut the U2/D2 module-face jump from 4.8 to 2.35; <= 2.0 (the 2.5a review's
# ask for 2.5b) would need graded perimeter tile sides, i.e. a change of the
# module tiles (not done: far field, no claim quantity). ACCEPTED 2026-09-25
# (delegated decision; WP 2.5c, K0-design §10.10 item 5, D-I addendum row):
# all-fluid 2.33 at medium with the module tiles unchanged; the port regions
# stay held to ER_MAX_PORTS. 2.5 is the regression guard, not a target.
ER_GUARD_ALL_FLUID_HOUSING = 2.5
# Outflow column (D1 slot column + riser): the outflow must not set the
# Courant limit. Judged at the highest CORE flow (Q 500 mL/min, plan NEW-4)
# on the smallest z-cell of the column, with the PEAK velocity: the mean
# slot velocity Q / (10 x 46 mm) times the developed-laminar peak/mean
# ratio of the 10 x 46 mm section (series solution, laminar_peak_to_mean:
# 1.7355). The riser foot is developing and skewed and can exceed it (the
# inlet-elbow |U| max was 1.8x the slot mean in the coarse smoke); the
# binding check is the solver's Co_max (k0_smoke_log.py / gates.py, §9).
#   Co_peak = PEAK * (Q / slot area) * dt / dz_min <= 1
# at the level's Q 500 time step: dt* = 6.25 ms (plan §1.1) on coarse,
# medium and fine; dt*/2 on xfine, where the plan runs Q 500 only in the
# NEW-19 resolution test at dt*/2 (its D-C cases peak at the C1 on-flow,
# 360 mL/min, at dt*: Co_peak 0.87-0.88). xfine at Q 500 and dt* is
# 1.21-1.22 (not allowed); K0-design §10.5 lists the largest dt per level
# and flow (NEW-19's Q 2000 needs dt*/4 on fine and dt*/8 on xfine).
CO_OUTFLOW_Q_M3S = 500e-6 / 60.0
CO_OUTFLOW_DT_S = 6.25e-3
CO_OUTFLOW_DT_BY_LEVEL = {"coarse": 6.25e-3, "medium": 6.25e-3,
                          "fine": 6.25e-3, "xfine": 3.125e-3}
CO_OUTFLOW_MAX = 1.0
# z-lever: the cap-plane (z = H) cells of the port columns must be at least
# this multiple of the module's cap-plane cells (0.05 mm at medium), i.e.
# the 0.05 mm Courant hot layer does not extend into the ports
CAPPLANE_LEVER_MIN = 4.0

# ---- module-interior z-lever (K0_ZLEVER=module; K0-design §10.14) ------------
# P00b (docs/results/P00b-summary.md §3): on the medium housing at 12.5 ms
# the heated plumes cross the 0.05 mm cap-plane layer at Co 4.2 in the
# fluid beside the cells; the summary gives only the box x ~ -12 ... -15 mm,
# y ~ +-1 ... +-16 mm. That box does NOT place the samples in wall
# distance d: over it d runs from -1.19 mm (inside a cell, at (-15, +-16))
# through 0.22-1.0 mm (the y ~ +-16 end) to 1.30-3.04 mm (the y ~ +-1
# end); the per-sample d comes from the run's CoMinMax_fluid locations
# (scripts/mesh/k0_co_lever.py samples). w = Co dz / dt = 4.2 x 0.05 mm /
# 12.5 ms = 16.8 ~ 17 mm/s is the vertical velocity that Co_max implies
# for a purely vertical crossing of the 0.05 mm cells (an identity, not a
# check of the velocity). The design bands below judge the lever by the
# Courant number w dt / dz of that plume on the fluid cells either side
# of the cap plane, per wall-distance band (PLUME_D_BINS_MM) and for the
# design band d >= PLUME_D_MIN_MM, at the three Delta-t branches of the
# re-plan. Which band the P00b maximum sits in, and w nearer the wall, are
# not known from the summary; the solver's Co_max (§9 gate) stays binding.
PLUME_W_M_S = 17e-3
PLUME_DT_S = (12.5e-3, 6.25e-3, 3.125e-3)
PLUME_D_MIN_MM = 1.0
PLUME_D_BINS_MM = (0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 6.0)
# design limit of the module lever: the z-lines of the fluid ring tilt
# between the fine cylinder line (r = R) and the coarser tile-edge line;
# on the gap ray (1 mm, the thinnest ring) the tilt is atan(D / 1 mm), D
# the largest offset between the two lines' points (K0-design §10.14)
ZLEVER_TILT_MAX_DEG = 30.0
# gate item when the module lever is requested (K0_ZLEVER=module): the
# largest cap-plane cell of the module region (the tile-corner lines) must
# be >= CAPPLANE_LEVER_MIN x the cap cell, and the plume band's smallest
# cap-plane cell >= ZLEVER_MODULE_BAND_MIN x the cap cell; when it is not
# requested the module region must carry the uniform cap layer (ratio 1)
ZLEVER_MODULE_BAND_MIN = 1.5
ZLEVERS = ("far", "module")


def co_z(w_m_s: float, dt_s: float, dz_m: float) -> float:
    """Vertical Courant number w dt / dz."""
    return w_m_s * dt_s / dz_m


def dt_co1(w_m_s: float, dz_m: float) -> float:
    """Largest dt [s] with w dt / dz <= 1."""
    return dz_m / w_m_s


def plume_courant(dz_mm: float, w_m_s: float = PLUME_W_M_S) -> dict:
    """Co of the P00b plume at the PLUME_DT_S branches and the largest dt
    for Co <= 1, for a vertical cell height dz [mm]."""
    dz = dz_mm * 1e-3
    return dict(dz_mm=dz_mm, w_mm_s=1e3 * w_m_s,
                co={f"{1e3 * dt:g}ms": co_z(w_m_s, dt, dz)
                    for dt in PLUME_DT_S},
                dt_co1_ms=1e3 * dt_co1(w_m_s, dz))


def laminar_peak_to_mean(aspect: float, terms: int = 200) -> float:
    """u_max / u_mean of developed laminar flow in a rectangular duct of
    side ratio `aspect` = short/long (0 < aspect <= 1): the classical
    Fourier-series solution (White, Viscous Fluid Flow, eq. 3-48); 1.5 for
    parallel plates, 2.096 for a square."""
    if not 0 < aspect <= 1:
        raise ValueError(f"aspect {aspect}: need 0 < short/long <= 1")
    a, b = 1.0, aspect                   # half-sides, |y| < a, |z| < b
    s_q = s_c = 0.0
    for k in range(terms):
        i = 2 * k + 1
        x = i * math.pi * a / (2 * b)
        s_q += math.tanh(x) / i ** 5
        s_c += (-1) ** k * (1 - (1 / math.cosh(x) if x < 700 else 0.0)) \
            / i ** 3
    u_mean = (b ** 2 / 3) * (1 - 192 * b / (math.pi ** 5 * a) * s_q)
    u_max = 16 * b ** 2 / math.pi ** 3 * s_c
    return u_max / u_mean


# the 10 x 46 mm slot / riser section (D-G): 1.7355
CO_OUTFLOW_PEAK_FACTOR = laminar_peak_to_mean(10.0 / 46.0)


def co_outflow_dt(level) -> float:
    """Q 500 reference time step of a K0 level (dt* if unknown)."""
    return CO_OUTFLOW_DT_BY_LEVEL.get(level, CO_OUTFLOW_DT_S)


def coolant_nu(props_path: Path = PROPS) -> float:
    """Kinematic viscosity mu / rho [m2/s] from properties.yaml (hard
    failure if missing: never a silent fallback)."""
    import yaml
    c = yaml.safe_load(Path(props_path).read_text())["coolant"]
    return float(c["dynamic_viscosity"]) / float(c["density"])


def stokes_thickness(f_hz: float, nu: float) -> float:
    """delta_S = sqrt(2 nu / omega) [m]."""
    return math.sqrt(nu / (math.pi * f_hz))


def stokes_gate_thickness(f_hz: float, nu: float) -> float:
    """Thickness the D-I gate counts cells in at f [m]: the thinner of the
    computed layer and the value the D-I text quotes (if it quotes one)."""
    d = stokes_thickness(f_hz, nu)
    q = DI_STOKES_QUOTED_MM.get(float(f_hz))
    return min(d, q * 1e-3) if q is not None else d


def nu_min_for(f_hz: float, depth: float) -> float:
    """Smallest nu for which the Stokes layer at f still reaches `depth`
    [m] (e.g. the end of the 3rd wall cell): nu = pi f depth^2."""
    return math.pi * f_hz * depth ** 2
