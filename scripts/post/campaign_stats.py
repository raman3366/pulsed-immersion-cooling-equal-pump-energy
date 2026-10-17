#!/usr/bin/env python3
"""Statistics for the resubmission campaign (WP 2.4).

Pure functions, no I/O. Implements the rules of the plan of record
(docs/RESUBMISSION-PLAN-2026-09-25.md §3) and the audit register:

- D-F   U95 of a contrast = u_ens (Welch t) (+) u_num, u_num = u_grid (+) u_time
        (+) u_iter in quadrature; power and minimum detectable effect (MDE)
        are stated WITH u_num included.
- D-L   Brown-Forsythe/Levene gate before pooling variances (rejected or
        untestable -> no pooled s); VM-effect one-way ANOVA; sqrt(k)
        comparator allocation; N_arm rule (n_arm); exact Mann-Whitney for
        the small-N amplitude-reversal pair.
- NEW-13  Holm within confirmatory families, Benjamini-Hochberg for the
        exploratory rest; TOST "practically negligible" with the full U95.
- R-23  standard uncertainties combined in quadrature, expanded with a
        Student-t coverage factor at the Welch-Satterthwaite effective dof.

Conventions: every uncertainty argument named ``u_*`` is a STANDARD
uncertainty (k = 1). Expanded bands (GCI, Stern, E&H U) must be converted
first (gci2.py returns both). All tests are two-sided unless stated.

Brown-Forsythe p-values are Monte Carlo p-values under the normal null at the
observed group sizes (BF_NSIM draws, fixed seed BF_SEED, so a result is
reproducible). The classical F reference is reported alongside as p_F but
never gates: with groups of n <= 3 it is badly mis-sized (all-N2 pools always
reject because each N2 group's two median deviations are equal; all-N3 pools
never reject because of the structural zero; six N2 + one N10 rejects
~22 % at nominal 5 %). A pool in which every group has n = 2 is untestable
by any median-deviation test and is refused (no pooled s).

PRE-REGISTERED CONVENTIONS (WP 2.10 part 2, 2026-09-28; PROJECT.md §2.5,
docs/prereg/PREREGISTRATION.md). The three conventions this module used to
leave open are fixed; the registered constants below are quoted in
docs/prereg/prereg.yaml and tests/prereg assert equality.
1. Equivalence and classes (D-F.3). A contrast's U95 = t_{0.975, nu_eff} x
   sqrt(u_ens^2 + u_num^2) (expanded_uncertainty; u_num standard = U / k_num,
   gci2.K_NUM 1.15, infinite dof). "equivalent" (SESOI) and "practically
   negligible" (0.1 K) hold iff |estimate| + U95 < margin: the 95 % band
   inside the margin, i.e. two one-sided tests at TOST_ALPHA = 0.025 each
   on u_c (tost(..., alpha=TOST_ALPHA)). tost_expanded() (the plan's
   1.645 se (+) U_num model, review of prep/2.4) is a planning model only and
   never classifies. The planning pass rate of the S3 TOST at s = 20 mK, N8
   vs N16, A7 members at N4 is 0.65 / 0.46 at a true Delta of 0 / 30 mK under
   the registered rule (0.92 / 0.80 under the plan model; 20,000 draws,
   docs/prereg/PREREGISTRATION.md §2.5 and §7.4). classify_contrast() applies the class rules
   with the registered precedence.
2. N_arm (D-L.1): the plan convention (N_ARM_CONVENTION = "plan": s_up =
   1.25 s_pilot, normal approximation + 1), which gives the filed knife
   edges 20.8 / 22.2 mK; n_arm_registered().
3. All-N2 pools (D-L.2): never tested for homogeneity and never pooled for a
   confirmatory contrast. A level's standard error uses max(s_ref, s_level)
   (level_sd(): s_ref the Brown-Forsythe-pooled stage-1 s, s_level the
   untested pooled within-case s of the level's N2 cases), and variance
   persistence is a two-sided F test of s_level^2 / s_ref^2 at 5 %
   (reported, never a gate).
Also registered here: the primary-window selection (select_window, D-F.4),
the precision-extension trigger (precision_extension, D-L.3) and its group
rule (precision_extension_groups, D-L.4), the per-arm gatekeeping of the
secondary families (holm_gatekept, NEW-13.1), the storyline assignment
(storyline, NEW-13.1), the H_E front test (he_test, pair_contrast,
bridge_offset; NEW-13.2), the single-arm comparator size
(comparator_n_single, D-L.4), the VM-effect term (u_ens_with_vm, D-L.4),
the heat-model label (heat_model_dependent, D-A(4)) and the mechanism and
headspace shares (mechanism_causal, headspace_triggered).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, NamedTuple, Sequence

import numpy as np
from scipy import integrate, optimize
from scipy import stats as _st

INF = math.inf
BF_NSIM = 20000          # Monte Carlo draws for the Brown-Forsythe p-value
BF_SEED = 20260925       # fixed: the p-value is a deterministic function of the data


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _sample(x, name="sample", min_n=1) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    if a.size < min_n:
        raise ValueError(f"{name}: need >= {min_n} values, got {a.size}")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name}: non-finite values")
    return a


def _u_total(u) -> float:
    """Standard uncertainty from a scalar or an iterable of independent
    components (combined in quadrature)."""
    if u is None:
        return 0.0
    if np.ndim(u) == 0:
        return combine_quadrature(float(u))
    return combine_quadrature(*[float(v) for v in u])


def combine_quadrature(*u: float) -> float:
    """sqrt(sum u_i^2) of independent standard uncertainties."""
    vals = np.asarray(u, dtype=float).ravel()
    if not np.all(np.isfinite(vals)) or np.any(vals < 0):
        raise ValueError(f"uncertainties must be finite and >= 0: {vals}")
    return float(math.sqrt(float(np.sum(vals ** 2))))


def welch_satterthwaite(u: Sequence[float], dof: Sequence[float]) -> float:
    """Effective dof of a quadrature sum: (sum u^2)^2 / sum(u^4 / nu).
    Components with dof = inf (Type B, e.g. u_num) contribute no term.
    Returns inf when every nonzero component has infinite dof."""
    u = np.asarray(u, dtype=float)
    dof = np.asarray(dof, dtype=float)
    if u.shape != dof.shape:
        raise ValueError("u and dof must have the same length")
    if np.any(dof <= 0):
        raise ValueError("dof must be > 0")
    num = float(np.sum(u ** 2)) ** 2
    finite = np.isfinite(dof) & (u > 0)
    den = float(np.sum(u[finite] ** 4 / dof[finite]))
    if den == 0.0:
        return INF
    return num / den


def _t_ppf(q: float, dof: float) -> float:
    return float(_st.norm.ppf(q)) if not np.isfinite(dof) else float(_st.t.ppf(q, dof))


def _t_sf(x: float, dof: float) -> float:
    return float(_st.norm.sf(x)) if not np.isfinite(dof) else float(_st.t.sf(x, dof))


def _t_cdf(x: float, dof: float) -> float:
    return float(_st.norm.cdf(x)) if not np.isfinite(dof) else float(_st.t.cdf(x, dof))


# --------------------------------------------------------------------------
# Welch two-sample t
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class WelchResult:
    diff: float          # mean(x) - mean(y)
    se: float            # standard error of diff (= u_ens)
    dof: float           # Welch-Satterthwaite dof
    t: float
    p: float             # two-sided
    ci_low: float        # (1 - alpha) two-sided CI of diff
    ci_high: float
    alpha: float
    n1: int
    n2: int


def welch_from_stats(m1: float, s1: float, n1: int, m2: float, s2: float,
                     n2: int, alpha: float = 0.05) -> WelchResult:
    """Welch test from summary statistics (s = sample SD, ddof = 1)."""
    if n1 < 2 or n2 < 2:
        raise ValueError("Welch test needs n >= 2 in both groups")
    if s1 < 0 or s2 < 0:
        raise ValueError("standard deviations must be >= 0")
    v1, v2 = s1 ** 2 / n1, s2 ** 2 / n2
    se = math.sqrt(v1 + v2)
    if se == 0.0:
        raise ValueError("zero variance in both groups: Welch t undefined")
    dof = welch_satterthwaite([math.sqrt(v1), math.sqrt(v2)], [n1 - 1, n2 - 1])
    diff = m1 - m2
    t = diff / se
    p = min(1.0, 2.0 * _t_sf(abs(t), dof))
    tc = _t_ppf(1.0 - alpha / 2.0, dof)
    return WelchResult(diff, se, dof, t, p, diff - tc * se, diff + tc * se,
                       alpha, int(n1), int(n2))


def welch_t(x, y, alpha: float = 0.05) -> WelchResult:
    """Welch two-sample t on raw samples (two-sided p and CI)."""
    x = _sample(x, "x", 2)
    y = _sample(y, "y", 2)
    return welch_from_stats(float(x.mean()), float(x.std(ddof=1)), x.size,
                            float(y.mean()), float(y.std(ddof=1)), y.size, alpha)


# --------------------------------------------------------------------------
# U95 of a contrast (D-F)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class ExpandedU:
    u_c: float           # combined standard uncertainty
    dof_eff: float       # Welch-Satterthwaite effective dof
    k: float             # coverage factor t_{(1+cov)/2, dof_eff}
    U: float             # expanded uncertainty k * u_c
    coverage: float


def expanded_uncertainty(u_ens: float, dof_ens: float, u_num=0.0,
                         coverage: float = 0.95,
                         u_num_dof: float = INF) -> ExpandedU:
    """U = k * sqrt(u_ens^2 + u_num^2); k from Student t at the effective dof.
    ``u_num`` may be a scalar or the components (u_grid, u_time, u_iter)."""
    un = _u_total(u_num)
    if u_ens < 0 or not math.isfinite(u_ens):
        raise ValueError("u_ens must be finite and >= 0")
    u_c = combine_quadrature(u_ens, un)
    dof_eff = welch_satterthwaite([u_ens, un], [dof_ens, u_num_dof])
    k = _t_ppf(0.5 + coverage / 2.0, dof_eff)
    return ExpandedU(u_c, dof_eff, k, k * u_c, coverage)


@dataclass(frozen=True)
class ContrastResult:
    diff: float
    u_ens: float
    dof_ens: float
    u_num: float
    u_c: float
    dof_eff: float
    k: float
    U95: float
    ci_low: float
    ci_high: float
    resolved: bool       # the U95 band excludes zero (before multiplicity)
    p: float             # two-sided p of diff / u_c at dof_eff


def contrast_u95(x, y, u_num=0.0, coverage: float = 0.95,
                 u_num_dof: float = INF) -> ContrastResult:
    """D-F contrast: Welch u_ens (+) u_num, expanded at the effective dof."""
    w = welch_t(x, y, alpha=1.0 - coverage)
    e = expanded_uncertainty(w.se, w.dof, u_num, coverage, u_num_dof)
    un = _u_total(u_num)
    p = min(1.0, 2.0 * _t_sf(abs(w.diff) / e.u_c, e.dof_eff)) if e.u_c > 0 else 0.0
    return ContrastResult(w.diff, w.se, w.dof, un, e.u_c, e.dof_eff, e.k, e.U,
                          w.diff - e.U, w.diff + e.U, abs(w.diff) > e.U, p)


# --------------------------------------------------------------------------
# variance homogeneity and pooling (D-L)
# --------------------------------------------------------------------------

class BFResult(NamedTuple):
    W: float             # Brown-Forsythe statistic (ANOVA F of |x - median_i|)
    p: float             # Monte Carlo p under the normal null (the gate; nan if untestable)
    p_F: float           # classical F-reference p (scipy levene, center='median'); report only
    testable: bool
    n_sim: int
    note: str


def _bf_stat(ssb: np.ndarray, ssw: np.ndarray, k: int, N: int) -> np.ndarray:
    ssb = np.asarray(ssb, dtype=float)
    ssw = np.asarray(ssw, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        W = (ssb / (k - 1)) / (ssw / (N - k))
    W = np.where(ssw > 0, W, np.where(ssb > 0, INF, 0.0))
    return W


def _bf_parts(zs):
    N = sum(z.size for z in zs)
    grand = float(sum(z.sum() for z in zs)) / N
    ssb = float(sum(z.size * (z.mean() - grand) ** 2 for z in zs))
    ssw = float(sum(((z - z.mean()) ** 2).sum() for z in zs))
    return ssb, ssw


def _bf_null(sizes, n_sim: int, seed: int, chunk: int = 2000) -> np.ndarray:
    """Brown-Forsythe statistics of n_sim normal-null data sets with the
    given group sizes (W is location- and scale-invariant, so N(0, 1) is
    the exact null for normal members with equal variances)."""
    rng = np.random.default_rng(seed)
    sizes = [int(n) for n in sizes]
    k, N = len(sizes), sum(sizes)
    uniq = sorted(set(sizes))
    count = {n: sizes.count(n) for n in uniq}
    out = []
    done = 0
    while done < n_sim:
        m = min(chunk, n_sim - done)
        zsum = np.zeros(m)
        parts = []
        for n in uniq:
            x = rng.standard_normal((m, count[n], n))
            z = np.abs(x - np.median(x, axis=2, keepdims=True))
            zbar = z.mean(axis=2)                                  # (m, count)
            ss = ((z - zbar[..., None]) ** 2).sum(axis=(1, 2))     # (m,)
            parts.append((n, zbar, ss))
            zsum += n * zbar.sum(axis=1)
        grand = zsum / N
        ssb = sum((n * (zbar - grand[:, None]) ** 2).sum(axis=1) for n, zbar, _ in parts)
        ssw = sum(ss for _, _, ss in parts)
        out.append(_bf_stat(ssb, ssw, k, N))
        done += m
    return np.concatenate(out)


def brown_forsythe(*groups, n_sim: int = BF_NSIM, seed: int = BF_SEED) -> BFResult:
    """Brown-Forsythe (median-centred Levene) test of equal variances.

    The p-value that gates is a Monte Carlo p under the normal null at the
    observed group sizes, (1 + #{W_null >= W}) / (1 + n_sim): exact-size for
    any mix of small groups. A pool in which every group has n = 2 is
    untestable (W is infinite for every data set: each N2 group's two
    deviations are equal) and returns testable = False with p = nan.
    """
    gs = [_sample(g, f"group{i}", 2) for i, g in enumerate(groups)]
    if len(gs) < 2:
        raise ValueError("need >= 2 groups")
    if int(n_sim) < 99:
        raise ValueError("n_sim must be >= 99")
    k = len(gs)
    N = sum(g.size for g in gs)
    if all(np.all(g == np.median(g)) for g in gs):
        return BFResult(0.0, 1.0, 1.0, True, 0,
                        "every deviation zero: no evidence of heterogeneity")
    zs = [np.abs(g - np.median(g)) for g in gs]
    ssb, ssw = _bf_parts(zs)
    W = float(_bf_stat(ssb, ssw, k, N))
    if np.isfinite(W):
        p_F = float(_st.f.sf(W, k - 1, N - k))
    else:
        p_F = 0.0
    if max(g.size for g in gs) < 3:
        return BFResult(W, math.nan, p_F, False, 0,
                        "untestable: every group has n = 2 (median deviations "
                        "equal within each group)")
    null = _bf_null([g.size for g in gs], int(n_sim), int(seed))
    p = float((1 + np.count_nonzero(null >= W * (1.0 - 1e-12))) / (1 + null.size))
    note = "Monte Carlo p (normal null, observed group sizes)"
    if min(g.size for g in gs) <= 3:
        note += f"; small groups present - classical F p {p_F:.3g} not used"
    return BFResult(W, p, p_F, True, int(n_sim), note)


@dataclass(frozen=True)
class PooledResult:
    s2: float            # pooled variance (nan unless status == 'pooled')
    s: float             # pooled SD (nan unless status == 'pooled')
    dof: int
    levene_W: float
    levene_p: float      # Monte Carlo Brown-Forsythe p (nan if untestable)
    pooled_ok: bool      # status == 'pooled'
    group_s: tuple = field(default_factory=tuple)
    status: str = "pooled"   # pooled | rejected | untestable
    levene_p_F: float = math.nan
    note: str = ""


def pooled_variance(groups: Iterable, alpha: float = 0.05,
                    n_sim: int = BF_NSIM, seed: int = BF_SEED) -> PooledResult:
    """Pooled within-group variance sum((n_i-1) s_i^2) / sum(n_i-1), released
    only if the Brown-Forsythe test is possible and does not reject
    homogeneity at ``alpha`` (D-L: "Levene before pooling; if rejected, no
    pooled s"). An untestable pool (every group n = 2) is refused as well."""
    gs = [_sample(g, f"group{i}", 2) for i, g in enumerate(groups)]
    bf = brown_forsythe(*gs, n_sim=n_sim, seed=seed)
    dof = int(sum(g.size - 1 for g in gs))
    s2 = float(sum((g.size - 1) * g.var(ddof=1) for g in gs) / dof)
    gsd = tuple(float(g.std(ddof=1)) for g in gs)
    if not bf.testable:
        return PooledResult(math.nan, math.nan, dof, bf.W, bf.p, False, gsd,
                            "untestable", bf.p_F, bf.note)
    if not bf.p >= alpha:
        return PooledResult(math.nan, math.nan, dof, bf.W, bf.p, False, gsd,
                            "rejected", bf.p_F, bf.note)
    return PooledResult(s2, math.sqrt(s2), dof, bf.W, bf.p, True, gsd,
                        "pooled", bf.p_F, bf.note)


@dataclass(frozen=True)
class AnovaResult:
    F: float
    p: float
    df_between: int
    df_within: int
    ms_between: float
    ms_within: float
    var_between: float   # method-of-moments variance component, >= 0


def one_way_anova(*groups) -> AnovaResult:
    """One-way ANOVA (D-L VM-effect test) plus the between-group variance
    component (MSB - MSW)/n0, n0 = (N - sum n_i^2 / N)/(k - 1), floored at 0.
    A significant VM effect adds var_between to u_ens^2 (per member)."""
    gs = [_sample(g, f"group{i}", 1) for i, g in enumerate(groups)]
    k = len(gs)
    n = np.array([g.size for g in gs], dtype=float)
    N = float(n.sum())
    if k < 2 or N - k < 1:
        raise ValueError("ANOVA needs >= 2 groups and N > k")
    grand = float(np.concatenate(gs).mean())
    ssb = float(sum(g.size * (g.mean() - grand) ** 2 for g in gs))
    ssw = float(sum(((g - g.mean()) ** 2).sum() for g in gs))
    dfb, dfw = k - 1, int(N - k)
    msb, msw = ssb / dfb, ssw / dfw
    if msw == 0.0:
        F, p = (INF, 0.0) if msb > 0 else (math.nan, 1.0)
    else:
        F = msb / msw
        p = float(_st.f.sf(F, dfb, dfw))
    n0 = (N - float((n ** 2).sum()) / N) / (k - 1)
    var_b = max(0.0, (msb - msw) / n0)
    return AnovaResult(F, p, dfb, dfw, msb, msw, var_b)


# --------------------------------------------------------------------------
# multiplicity (NEW-13)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class MultipleResult:
    adjusted: np.ndarray
    reject: np.ndarray
    alpha: float


def _pvals(p) -> np.ndarray:
    p = np.asarray(p, dtype=float).ravel()
    if p.size == 0:
        raise ValueError("no p-values")
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must be finite and in [0, 1]")
    return p


def holm(p, alpha: float = 0.05) -> MultipleResult:
    """Holm step-down FWER control; adjusted p in the input order."""
    p = _pvals(p)
    m = p.size
    order = np.argsort(p, kind="stable")
    adj_sorted = np.minimum(1.0, np.maximum.accumulate((m - np.arange(m)) * p[order]))
    adj = np.empty(m)
    adj[order] = adj_sorted
    return MultipleResult(adj, adj <= alpha, alpha)


def benjamini_hochberg(p, alpha: float = 0.05) -> MultipleResult:
    """Benjamini-Hochberg step-up FDR control; adjusted p in the input order."""
    p = _pvals(p)
    m = p.size
    order = np.argsort(p, kind="stable")
    ranked = p[order] * m / np.arange(1, m + 1)
    adj_sorted = np.minimum(1.0, np.minimum.accumulate(ranked[::-1])[::-1])
    adj = np.empty(m)
    adj[order] = adj_sorted
    return MultipleResult(adj, adj <= alpha, alpha)


# --------------------------------------------------------------------------
# TOST equivalence (D-F classes; NEW-13 "practically negligible")
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class TostResult:
    diff: float
    margin: float
    u_c: float           # sqrt(u_ens^2 + u_num^2)
    dof_eff: float
    t_crit: float        # t_{1-alpha, dof_eff}
    ci_low: float        # (1 - 2 alpha) interval
    ci_high: float
    p_lower: float       # H0: diff <= -margin
    p_upper: float       # H0: diff >= +margin
    p: float             # max(p_lower, p_upper)
    equivalent: bool
    alpha: float


def tost(x=None, y=None, *, margin: float, alpha: float = 0.05, u_num=0.0,
         diff: float | None = None, se: float | None = None,
         dof: float | None = None, u_num_dof: float = INF) -> TostResult:
    """Two one-sided tests for |true diff| < margin.

    Give samples ``x, y`` (Welch u_ens and dof) or a summary ``diff, se, dof``.
    ``u_num`` (standard, scalar or components) is added in quadrature and the
    dof become the Welch-Satterthwaite effective dof, so equivalence is only
    declared when the (1 - 2 alpha) interval of diff, widened by u_num, lies
    inside (-margin, +margin) — the "full U95" rule of NEW-13.
    """
    if not (margin > 0 and math.isfinite(margin)):
        raise ValueError("margin must be > 0")
    if x is not None or y is not None:
        if x is None or y is None or diff is not None:
            raise ValueError("give both samples x, y or (diff, se, dof), not a mix")
        w = welch_t(x, y)
        d, s, nu = w.diff, w.se, w.dof
    else:
        if diff is None or se is None or dof is None:
            raise ValueError("summary mode needs diff, se and dof")
        if se < 0 or dof <= 0:
            raise ValueError("se must be >= 0 and dof > 0")
        d, s, nu = float(diff), float(se), float(dof)
    un = _u_total(u_num)
    u_c = combine_quadrature(s, un)
    nu_eff = welch_satterthwaite([s, un], [nu, u_num_dof]) if u_c > 0 else INF
    tc = _t_ppf(1.0 - alpha, nu_eff)
    if u_c == 0.0:
        inside = abs(d) < margin
        pl = 0.0 if d > -margin else 1.0
        pu = 0.0 if d < margin else 1.0
        return TostResult(d, margin, 0.0, nu_eff, tc, d, d, pl, pu, max(pl, pu),
                          inside, alpha)
    pl = _t_sf((d + margin) / u_c, nu_eff)
    pu = _t_cdf((d - margin) / u_c, nu_eff)
    p = max(pl, pu)
    return TostResult(d, margin, u_c, nu_eff, tc, d - tc * u_c, d + tc * u_c,
                      pl, pu, p, bool(p < alpha), alpha)


@dataclass(frozen=True)
class TostExpandedResult:
    diff: float
    margin: float
    se: float
    dof: float
    k_ens: float         # t_{1-alpha, dof} (z when dof = inf)
    U_num: float         # EXPANDED numerical band, used as is
    band: float          # sqrt((k_ens se)^2 + U_num^2)
    equivalent: bool     # |diff| + band <= margin
    alpha: float


def tost_expanded(diff: float, se: float, dof: float, U_num: float,
                  margin: float, alpha: float = 0.05) -> TostExpandedResult:
    """The plan's TOST rule (models/tost_a7.py): equivalent iff
    |d| + sqrt((t_{1-alpha, dof} se)^2 + U_num^2) <= margin, with the
    EXPANDED numerical band U_num (e.g. Fs F |level difference| from gci2)
    added as is. dof = inf reproduces the model's 1.645. Compare tost(),
    which works from the standard u_num = U / k_num inside t * u_c; which
    convention is confirmatory is an open pre-registration item (module
    docstring)."""
    if not (margin > 0 and math.isfinite(margin)):
        raise ValueError("margin must be > 0")
    vals = (diff, se, U_num)
    if not all(math.isfinite(float(v)) for v in vals) or se < 0 or U_num < 0:
        raise ValueError("diff, se and U_num must be finite; se, U_num >= 0")
    if not dof > 0:
        raise ValueError("dof must be > 0")
    k = _t_ppf(1.0 - alpha, dof)
    band = math.hypot(k * se, U_num)
    return TostExpandedResult(float(diff), margin, float(se), float(dof), k,
                              float(U_num), band, bool(abs(diff) + band <= margin),
                              alpha)


# --------------------------------------------------------------------------
# exact Mann-Whitney (D-L: 1 Hz A0.4 vs A0.6 reversal, N4 vs N4)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class MannWhitneyResult:
    U: float             # U of x: #(x > y) + 0.5 #(x == y)
    p: float
    n1: int
    n2: int
    ties: bool
    alternative: str


MW_MAX_N = 60            # int64 subset counts stay exact up to C(60, 30)


def mann_whitney_exact(x, y, alternative: str = "two-sided") -> MannWhitneyResult:
    """Exact permutation Mann-Whitney U test for small samples.

    The null distribution of the rank sum is enumerated exactly (dynamic
    programme over subsets of the pooled midranks), so ties are handled by
    the exact conditional permutation distribution. Two-sided p =
    min(1, 2 min(P(U >= u), P(U <= u))) — identical to scipy's exact
    method when there are no ties. 'greater' tests x stochastically larger.
    """
    x = _sample(x, "x", 1)
    y = _sample(y, "y", 1)
    n1, n2 = x.size, y.size
    n = n1 + n2
    if n > MW_MAX_N:
        raise ValueError(f"exact enumeration limited to n1 + n2 <= {MW_MAX_N}")
    if alternative not in ("two-sided", "greater", "less"):
        raise ValueError("alternative must be two-sided, greater or less")
    pooled = np.concatenate([x, y])
    r2 = np.rint(2.0 * _st.rankdata(pooled)).astype(np.int64)   # doubled midranks
    ties = bool(np.unique(pooled).size < n)
    smax = int(r2.sum())
    ways = np.zeros((n1 + 1, smax + 1), dtype=np.int64)
    ways[0, 0] = 1
    for i, r in enumerate(r2):
        for k in range(min(i + 1, n1), 0, -1):
            ways[k, r:] += ways[k - 1, :smax + 1 - r]
    dist = ways[n1]
    total = float(dist.sum())                     # = C(n, n1)
    s_obs = int(r2[:n1].sum())
    p_ge = float(dist[s_obs:].sum()) / total
    p_le = float(dist[:s_obs + 1].sum()) / total
    U = s_obs / 2.0 - n1 * (n1 + 1) / 2.0
    if alternative == "greater":
        p = p_ge
    elif alternative == "less":
        p = p_le
    else:
        p = min(1.0, 2.0 * min(p_ge, p_le))
    return MannWhitneyResult(U, p, n1, n2, ties, alternative)


def complete_separation(x, y) -> bool:
    """True iff min(x) > max(y) — the D-L amplitude-reversal claim criterion
    ("min(A0.4) > max(A0.6) over all members")."""
    return bool(_sample(x, "x").min() > _sample(y, "y").max())


# --------------------------------------------------------------------------
# allocation, power and minimum detectable effect (D-L, D-F, B1)
# --------------------------------------------------------------------------

def sqrt_k_allocation(n_arm: int, k: int, cap: int | None = None) -> int:
    """Dunnett-style shared-comparator size ceil(sqrt(k) * n_arm), optionally
    capped (D-L: baseline min(40, ceil(sqrt(k) N_arm)); a comparator shared
    by k_c arms gets ceil(sqrt(k_c) N_arm))."""
    if int(n_arm) != n_arm or n_arm < 1 or int(k) != k or k < 1:
        raise ValueError("n_arm and k must be integers >= 1")
    n = int(math.ceil(math.sqrt(k) * n_arm - 1e-9))
    return min(n, int(cap)) if cap is not None else n


def power_two_sample(delta: float, s: float, n1: float, n2: float,
                     alpha: float = 0.05, u_num=0.0, two_sided: bool = True,
                     s2: float | None = None, dof_mode: str = "welch") -> float:
    """Probability that a contrast with true difference ``delta`` is resolved,
    i.e. |est| > t_{nu_eff} * sqrt(se_hat^2 + u_num^2) (one-sided: est > ...).

    Members are N(mu, s^2) (group 2: s2, default s). se_hat^2 is modelled as
    se^2 chi2_nu / nu with nu the Welch-Satterthwaite dof at the true
    variances (exact for equal n and s); ``dof_mode="pooled"`` uses
    n1 + n2 - 2 instead (pooled-variance t after a passed Levene gate).
    u_num is a fixed standard uncertainty (Type B, infinite dof). With
    u_num = 0 this is the exact noncentral-t power. The integral over chi2
    is done numerically.
    """
    s2 = s if s2 is None else s2
    if s <= 0 or s2 <= 0 or n1 < 2 or n2 < 2:
        raise ValueError("need s > 0 and n >= 2 per group")
    v1, v2 = s * s / n1, s2 * s2 / n2
    se = math.sqrt(v1 + v2)
    if dof_mode == "welch":
        nu = welch_satterthwaite([math.sqrt(v1), math.sqrt(v2)], [n1 - 1, n2 - 1])
    elif dof_mode == "pooled":
        nu = float(n1 + n2 - 2)
    else:
        raise ValueError("dof_mode must be 'welch' or 'pooled'")
    un = _u_total(u_num)
    q_tail = 1.0 - (alpha / 2.0 if two_sided else alpha)
    d = abs(delta) if two_sided else delta

    def resolved_prob(u: float) -> float:
        q = float(_st.chi2.ppf(u, nu)) / nu            # se_hat^2 / se^2
        se_hat2 = q * se * se
        if un > 0:
            nu_eff = nu * ((se_hat2 + un * un) / se_hat2) ** 2 if se_hat2 > 0 else INF
        else:
            nu_eff = nu
        c = _t_ppf(q_tail, nu_eff) * math.sqrt(se_hat2 + un * un)
        pr = float(_st.norm.cdf((d - c) / se))
        if two_sided:
            pr += float(_st.norm.cdf((-d - c) / se))
        return pr

    val, _ = integrate.quad(resolved_prob, 0.0, 1.0, limit=200,
                            epsabs=1e-10, epsrel=1e-8)
    return float(min(1.0, max(0.0, val)))


def mde(s: float, n1: float, n2: float, alpha: float = 0.05,
        power: float = 0.80, u_num=0.0, two_sided: bool = True,
        s2: float | None = None, dof_mode: str = "welch") -> float:
    """Minimum detectable effect: smallest |delta| resolved with probability
    ``power`` under power_two_sample (u_num included, D-F/B1)."""
    if not 0.0 < power < 1.0:
        raise ValueError("power must be in (0, 1)")
    f = lambda d: power_two_sample(d, s, n1, n2, alpha, u_num, two_sided, s2,
                                   dof_mode) - power
    s2_ = s if s2 is None else s2
    hi = 4.0 * (math.sqrt(s * s / n1 + s2_ * s2_ / n2) + _u_total(u_num)) + 1e-12
    while f(hi) < 0:
        hi *= 2.0
        if hi > 1e12:
            raise RuntimeError("MDE bracket failed")
    return float(optimize.brentq(f, 0.0, hi, xtol=1e-12, rtol=1e-10))


def n_for_power(delta: float, s: float, alpha: float = 0.05,
                power: float = 0.80, ratio: float = 1.0, u_num=0.0,
                two_sided: bool = True, n_min: int = 2,
                n_max: int = 10000, dof_mode: str = "welch") -> int:
    """Smallest n1 (with n2 = ceil(ratio * n1)) whose power reaches ``power``.
    Raises if u_num alone makes the target unreachable within n_max.
    Power is monotone in n, so the search is galloping + bisection."""
    def ok(n: int) -> bool:
        n2 = max(2, int(math.ceil(ratio * n - 1e-9)))
        return power_two_sample(delta, s, n, n2, alpha, u_num, two_sided,
                                dof_mode=dof_mode) >= power
    lo = max(2, n_min)
    if lo > n_max:
        raise ValueError("n_min > n_max")
    if ok(lo):
        return lo
    hi = lo
    while True:                       # invariant: not ok(lo)
        hi = min(n_max, max(hi + 1, 2 * hi))
        if ok(hi):
            break
        if hi == n_max:
            raise ValueError(f"power {power} not reachable with n <= {n_max}")
        lo = hi
    while hi - lo > 1:                # not ok(lo), ok(hi)
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi


def upper_sd_bound(s: float, dof: float, confidence: float = 0.80) -> float:
    """One-sided upper confidence bound of an SD: s sqrt(nu / chi2_{1-c, nu})
    (D-L sizes N on "the upper-80 % pilot s")."""
    if s < 0 or dof <= 0:
        raise ValueError("s >= 0 and dof > 0 required")
    return float(s * math.sqrt(dof / _st.chi2.ppf(1.0 - confidence, dof)))


N_ARM_DELTA = 0.040      # K, D-L detectable effect
N_ARM_K = 4              # confirmatory arms (NEW-13), Holm alpha / k
N_ARM_FLOOR, N_ARM_CAP = 6, 10
N_ARM_CONVENTIONS = ("plan", "pooled", "welch")
PLAN_S_UP = 1.25         # plan model: upper-80 % pilot s taken as 1.25 s


@dataclass(frozen=True)
class NArmResult:
    n_arm: int           # max(floor, n80) capped
    n80: int             # members for 80 % power before floor / cap
    comparator: int      # shared comparator ceil(sqrt(k) N_arm) (baseline)
    s_pilot: float
    dof_pilot: float
    s_upper: float       # the s the power is computed at
    convention: str
    alpha_test: float    # two-sided alpha per contrast (alpha / k)
    capped: bool


def n_arm(s_pilot: float, dof_pilot: float, *, convention: str,
          k: int = N_ARM_K, delta: float = N_ARM_DELTA, alpha: float = 0.05,
          power: float = 0.80, floor: int = N_ARM_FLOOR, cap: int = N_ARM_CAP,
          confidence: float = 0.80) -> NArmResult:
    """D-L: N_arm = max(6, n80) at |Delta| = 40 mK, two-sided Holm alpha/k,
    80 % power, on the upper-80 % pilot s, variance factor (1 + 1/sqrt(k))
    (arm vs a comparator of sqrt(k) N_arm members), capped at 10.

    ``convention`` (no default here; the pre-registered one is
    N_ARM_CONVENTION = 'plan', D-L.1, applied by n_arm_registered) fixes how
    the upper s and the test's dof are taken:
      'plan'   the plan-of-record model: s_up = 1.25 s (dof_pilot recorded
               only), n = (1 + 1/sqrt(k)) (z_{1-alpha/2k} + z_power)^2
               (s_up/delta)^2, N = ceil(n) + 1 (t correction). Reproduces
               the stated knife edges: 8, 9 from 20.8 mK, 10 from 22.2 mK.
      'pooled' s_up = upper_sd_bound(s, dof_pilot, 0.80); exact power of the
               pooled-variance t (dof n1 + n2 - 2) with n2 = ceil(sqrt(k) n1).
      'welch'  as 'pooled' with the Welch-Satterthwaite dof.
    """
    if convention not in N_ARM_CONVENTIONS:
        raise ValueError(f"convention must be one of {N_ARM_CONVENTIONS}")
    if not (s_pilot > 0 and math.isfinite(s_pilot)) or not dof_pilot > 0:
        raise ValueError("s_pilot > 0 and dof_pilot > 0 required")
    if int(k) != k or k < 1 or not delta > 0:
        raise ValueError("k >= 1 integer and delta > 0 required")
    a = alpha / k
    ratio = math.sqrt(k)
    if convention == "plan":
        s_up = PLAN_S_UP * s_pilot
        z = float(_st.norm.ppf(1.0 - a / 2.0) + _st.norm.ppf(power))
        n = (1.0 + 1.0 / ratio) * z * z * (s_up / delta) ** 2
        n80 = int(math.ceil(n - 1e-12)) + 1
    else:
        s_up = upper_sd_bound(s_pilot, dof_pilot, confidence)
        n80 = n_for_power(delta, s_up, alpha=a, power=power, ratio=ratio,
                          dof_mode=convention)
    na = min(cap, max(floor, n80))
    return NArmResult(na, n80, sqrt_k_allocation(na, k), float(s_pilot),
                      float(dof_pilot), s_up, convention, a, n80 > cap)


# --------------------------------------------------------------------------
# block-average effective sample size (serially correlated windows)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class BlockResult:
    n: int
    n_eff: float          # s^2 / se_opt^2, clipped to [1, n]
    se: float             # standard error of the mean at the optimal block
    block_size: int
    n_blocks: int
    tau_int: float        # integrated autocorrelation time in samples = n / (2 n_eff)
    converged: bool       # optimal-block criterion met before min_blocks
    table: tuple          # ((block_size, n_blocks, se, se_err), ...)


def block_average_neff(y, min_blocks: int = 8) -> BlockResult:
    """Flyvbjerg-Petersen (1989) blocking of a uniformly sampled series.

    Block sizes 1, 2, 4, ...; the standard error of the mean at block size B
    uses the trailing floor(n/B)*B samples. The optimal block is the
    smallest B with B^3 > 2 n (se_B / se_1)^4 (Lee et al. 2011, the
    criterion used by pyblock). If no B with >= min_blocks blocks meets it,
    the largest SE over those blocks is returned and converged = False
    (conservative). N_eff = var(y) / se_opt^2.
    """
    y = _sample(y, "series", 2 * min_blocks)
    n = y.size
    s = float(y.std(ddof=1))
    rows = []
    b = 1
    while n // b >= min_blocks:
        nb = n // b
        blocks = y[n - nb * b:].reshape(nb, b).mean(axis=1)
        se_b = float(blocks.std(ddof=1) / math.sqrt(nb))
        rows.append((b, nb, se_b, se_b / math.sqrt(2.0 * (nb - 1))))
        b *= 2
    se1 = rows[0][2]
    if se1 == 0.0:
        return BlockResult(n, float(n), 0.0, 1, n, 0.5, True, tuple(rows))
    opt = None
    for row in rows:
        if row[0] ** 3 > 2.0 * n * (row[2] / se1) ** 4:
            opt = row
            break
    converged = opt is not None
    if opt is None:
        opt = max(rows, key=lambda r: r[2])
    se = opt[2]
    n_eff = float(min(n, max(1.0, s * s / (se * se)))) if se > 0 else float(n)
    return BlockResult(n, n_eff, se, opt[0], opt[1], n / (2.0 * n_eff),
                       converged, tuple(rows))


# --------------------------------------------------------------------------
# pre-registered conventions (WP 2.10 part 2, 2026-09-28; PROJECT.md §2.5
# D-F.3, D-F.4, D-L.1, D-L.2, D-L.3, NEW-13.1; docs/prereg/prereg.yaml)
# --------------------------------------------------------------------------

FAMILY_ALPHA = 0.05      # FWER per confirmatory family (Holm), FDR (BH) per exploratory metric
COVERAGE = 0.95          # U95 of every contrast
TOST_ALPHA = 0.025       # D-F.3: each one-sided test; = the 95 % band inside the margin
NEGLIGIBLE_K = 0.10      # "practically negligible": |est| + U95 < 0.1 K
RELEVANT_ABS_K = 0.10    # "practically relevant": resolved, |est| >= 0.1 K ...
RELEVANT_REL_PCT = 10.0  # ... and >= 10 % of the reference (baseline or comparator)
# SESOI per metric kind: (absolute K, relative % of the reference or None);
# the bound is max(absolute, relative x |reference|) (D-F)
SESOI = {"tmax": (0.05, None), "dT": (0.01, 5.0), "module": (0.02, 5.0)}
CLASS_ORDER = ("practically_relevant", "resolved", "equivalent",
               "practically_negligible", "inconclusive")
N_ARM_CONVENTION = "plan"            # D-L.1
WINDOWS_PCT = (20, 50, 100)          # D-F.4 candidates, in tie-break order
PRECISION_HALF_WIDTH_K = {"uniformity": 0.010, "tmax": 0.025}   # D-L.3
PRECISION_ENS_SHARE = 0.5            # D-L.3: u_ens^2 >= 50 % of u_c^2
PRECISION_EXTRA_MEMBERS = 2          # D-L.3: +2 per group, one round only
PERSISTENCE_ALPHA = 0.05             # D-L.2: two-sided F test (reported)
CONFIRMATORY_ARMS = ("C1", "C2p", "C3", "C5")   # NEW-13 (k = 4)
S1_ARM = "C5"                        # NEW-13 storyline S1: C5 resolved in F-M1


def n_arm_registered(s_pilot: float, dof_pilot: float, **kw) -> NArmResult:
    """D-L.1: n_arm under the registered convention (plan model)."""
    return n_arm(s_pilot, dof_pilot, convention=N_ARM_CONVENTION, **kw)


def sesoi_bound(kind: str, ref: float | None = None) -> float:
    """SESOI of a metric kind ('tmax', 'dT' = Delta T_cell and sigma_T,
    'module' = module surface max - min) at reference value ``ref`` (the
    steady baseline or comparator value of the same metric, K)."""
    if kind not in SESOI:
        raise ValueError(f"kind must be one of {sorted(SESOI)}")
    a, rel = SESOI[kind]
    if rel is None or ref is None:
        return a
    return max(a, rel / 100.0 * abs(float(ref)))


def classify_contrast(est: float, U95: float, p_adj: float | None, *,
                      kind: str | None = None, ref: float | None = None,
                      alpha: float = FAMILY_ALPHA) -> str:
    """D-F / D-F.3 class of one contrast, registered precedence
    practically_relevant > resolved > equivalent > practically_negligible >
    inconclusive.

    est, U95 in K (U95 from contrast_u95 / expanded_uncertainty: the full
    band incl. u_num); p_adj the family-adjusted p (Holm for a confirmatory
    family, BH for an exploratory metric; None = not tested, e.g. a gated-
    out secondary contrast, which can never be resolved); kind the SESOI
    kind (None: no SESOI, so never 'equivalent'); ref the reference value
    for the relative thresholds (None: no relative test, so never
    'practically_relevant').
      resolved              p_adj <= alpha (and |est| > U95)
      practically_relevant  resolved, |est| >= 0.1 K and >= 10 % of |ref|
      equivalent            |est| + U95 < SESOI(kind, ref)
      practically_negligible |est| + U95 < 0.1 K
    """
    if not (math.isfinite(est) and math.isfinite(U95) and U95 >= 0):
        raise ValueError("est must be finite and U95 finite and >= 0")
    if p_adj is not None and not (0.0 <= p_adj <= 1.0):
        raise ValueError("p_adj must be in [0, 1] or None")
    resolved = p_adj is not None and p_adj <= alpha and abs(est) > U95
    if resolved:
        if (ref is not None and ref != 0 and abs(est) >= RELEVANT_ABS_K
                and 100.0 * abs(est) / abs(ref) >= RELEVANT_REL_PCT):
            return "practically_relevant"
        return "resolved"
    if kind is not None and abs(est) + U95 < sesoi_bound(kind, ref):
        return "equivalent"
    if abs(est) + U95 < NEGLIGIBLE_K:
        return "practically_negligible"
    return "inconclusive"


@dataclass(frozen=True)
class GatekeptResult:
    tested: dict             # key -> bool (gate open)
    p_adj: dict              # key -> Holm-adjusted p over the tested set (None if not tested)
    reject: dict             # key -> bool


def holm_gatekept(p_secondary: dict, gate: dict, alpha: float = FAMILY_ALPHA) -> GatekeptResult:
    """NEW-13.1: a secondary family (F-E2, F-M2) is tested only for the arms
    its primary family rejected. ``p_secondary`` maps (arm, metric) -> raw
    two-sided p; ``gate`` maps arm -> True iff the arm's primary hypothesis
    was rejected. Holm at ``alpha`` over the tested keys only; untested keys
    get p_adj None and no rejection (never 'resolved')."""
    tested = {k: bool(gate.get(k[0], False)) for k in p_secondary}
    keys = [k for k in p_secondary if tested[k]]
    p_adj = {k: None for k in p_secondary}
    reject = {k: False for k in p_secondary}
    if keys:
        h = holm([p_secondary[k] for k in keys], alpha)
        for k, pa, r in zip(keys, h.adjusted, h.reject):
            p_adj[k], reject[k] = float(pa), bool(r)
    return GatekeptResult(tested, p_adj, reject)


def storyline(resolved_fm1: Iterable[str], he_rejected_primary: bool | None = None,
              he_rejected_system: bool | None = None) -> dict:
    """NEW-13 (rescaled to k = 4) / NEW-13.1: S1 iff C5 is resolved in F-M1;
    S2 iff the F-M1 resolved set is non-empty and inside {C1, C2', C3}; S3
    iff it is empty. The S4 overlay holds iff H_E (he_test) is rejected on
    the primary basis and the rejection also holds on the system basis
    (None = not yet evaluated -> S4 False). NEW-13.2: the system basis has
    no registered reference loop, so S4 is an EXPLORATORY overlay
    (s4_exploratory True); the confirmatory front statement is H_E on the
    primary basis (he_rejected_primary) and names its basis."""
    res = set(resolved_fm1)
    bad = res - set(CONFIRMATORY_ARMS)
    if bad:
        raise ValueError(f"not confirmatory arms: {sorted(bad)}")
    story = "S1" if S1_ARM in res else ("S2" if res else "S3")
    s4 = bool(he_rejected_primary) and bool(he_rejected_system)
    return dict(story=story, s4=s4, s4_exploratory=True,
                he_primary=bool(he_rejected_primary), resolved=sorted(res))


def select_window(groups_by_window: dict, n_sim: int = BF_NSIM,
                  seed: int = BF_SEED) -> dict:
    """D-F.4: the primary window of a uniformity metric, chosen once at
    FREEZE from the stage-1 members, blind to every delta. For each
    candidate window w (percent of the discharge at its end, WINDOWS_PCT),
    ``groups_by_window[w]`` is a list of per-class arrays of member window
    means. The window's s is the Brown-Forsythe-pooled within-class SD; if
    pooling is refused (rejected or untestable), the largest class SD. The
    smallest s wins; an exact tie keeps the earlier window of WINDOWS_PCT
    (last 20 % first)."""
    s_by, how = {}, {}
    for w in WINDOWS_PCT:
        if w not in groups_by_window:
            continue
        pr = pooled_variance(groups_by_window[w], n_sim=n_sim, seed=seed)
        if pr.pooled_ok:
            s_by[w], how[w] = pr.s, "pooled"
        else:
            s_by[w], how[w] = max(pr.group_s), f"max class s ({pr.status})"
    unknown = set(groups_by_window) - set(WINDOWS_PCT)
    if unknown or not s_by:
        raise ValueError(f"windows must be among {WINDOWS_PCT} (got {sorted(groups_by_window)})")
    best = None
    for w in WINDOWS_PCT:
        if w in s_by and (best is None or s_by[w] < s_by[best]):
            best = w
    return dict(window=best, s=s_by, how=how)


def precision_extension(U95: float, kind: str, u_ens: float, u_c: float) -> bool:
    """D-L.3: the one precision-extension round (+2 members to each group
    of the contrast) fires iff the 95 % half-width U95 exceeds 0.025 K
    (surface T_max,peak, kind 'tmax') or 0.010 K (uniformity metrics) AND
    the ensemble term is >= 50 % of u_c^2. Variance components only: the
    estimate never enters."""
    if kind not in ("tmax", "dT", "module"):
        raise ValueError("kind must be 'tmax', 'dT' or 'module'")
    if not all(math.isfinite(v) and v >= 0 for v in (U95, u_ens, u_c)):
        raise ValueError("U95, u_ens and u_c must be finite and >= 0")
    thr = PRECISION_HALF_WIDTH_K["tmax" if kind == "tmax" else "uniformity"]
    return bool(U95 > thr and u_ens * u_ens >= PRECISION_ENS_SHARE * u_c * u_c)


@dataclass(frozen=True)
class LevelSD:
    s_level: float       # untested pooled within-case SD of the level's groups
    dof_level: int
    s_ref: float
    dof_ref: float
    s_use: float         # max(s_ref, s_level): the SD a level SE uses
    F: float             # s_level^2 / s_ref^2
    p: float             # two-sided F-test p (variance persistence, reported)
    persists: bool       # p >= PERSISTENCE_ALPHA


def level_sd(groups: Iterable, s_ref: float, dof_ref: float,
             alpha: float = PERSISTENCE_ALPHA) -> LevelSD:
    """D-L.2: the SD of a grid level or Delta-t rung whose cases run at N2
    (an all-N2 pool is untestable by Brown-Forsythe, so it is never pooled
    for a confirmatory contrast). s_level = sqrt(sum (n_i - 1) s_i^2 /
    sum (n_i - 1)) without a homogeneity test; the level's standard errors
    use s_use = max(s_ref, s_level), s_ref being the Brown-Forsythe-pooled
    stage-1 s; variance persistence is the two-sided F test of
    s_level^2 / s_ref^2 (reported, never a gate)."""
    gs = [_sample(g, f"group{i}", 2) for i, g in enumerate(groups)]
    if not gs:
        raise ValueError("need >= 1 group")
    if not (s_ref > 0 and math.isfinite(s_ref) and dof_ref > 0):
        raise ValueError("s_ref > 0 and dof_ref > 0 required")
    dof = int(sum(g.size - 1 for g in gs))
    s2 = float(sum((g.size - 1) * g.var(ddof=1) for g in gs) / dof)
    F = s2 / (s_ref * s_ref)
    p = float(min(1.0, 2.0 * min(_st.f.cdf(F, dof, dof_ref), _st.f.sf(F, dof, dof_ref))))
    s_lv = math.sqrt(s2)
    return LevelSD(s_lv, dof, float(s_ref), float(dof_ref), max(float(s_ref), s_lv),
                   float(F), p, bool(p >= alpha))


# --------------------------------------------------------------------------
# operational definitions fixed before the deposit (2026-09-28; PROJECT.md
# §2.5 D-L.4, NEW-13.2, NEW-30.1, D-A(4), NEW-20; docs/prereg/prereg.yaml)
# --------------------------------------------------------------------------

VM_ANOVA_ALPHA = 0.05                # D-L.4: a VM effect at p < 0.05 adds its component
HE_ENERGY_FLOOR_REL = 0.02           # NEW-13.2: = pump_energy.ENERGY_FLOOR_REL (D-H 2 %)
HEAT_LABEL_ABS_K = 0.005             # D-A(4): "heat-model-dependent" iff the delta moves by
HEAT_LABEL_REL = 0.25                # more than max(0.005 K, 0.25 |delta|, U95)
MECHANISM_SHARE = 0.50               # NEW-20: causal wording iff the attribution explains
                                     # >= 50 % of the delta in every member
HEADSPACE_SHARE = 0.50               # conditional 20 mm headspace runs iff the mechanism
                                     # attributes >= 50 % of a confirmatory delta to the ceiling layer


def comparator_n_single(n_arm_: int, s_pilot: float, *, k: int = N_ARM_K,
                        delta: float = N_ARM_DELTA, alpha: float = FAMILY_ALPHA,
                        power: float = 0.80) -> int:
    """D-L.4: size of a single-arm comparator (Q127, Q174) under the plan
    model (D-L.1): the smallest n_c with 80 % plan-model power for the
    contrast of N_arm members against it, z^2 (s_up/delta)^2 (1/N_arm +
    1/n_c) <= 1, s_up = 1.25 s_pilot, z = z_{1-alpha/2k} + z_power, then +1
    (the plan's t correction); never below the N_arm floor (6) and never
    above the baseline size ceil(sqrt(k) N_arm), which it takes when the arm
    alone exhausts the power budget (N_arm capped). N11 at s = 20 mK,
    N_arm 8 (the D-L central value)."""
    if int(n_arm_) != n_arm_ or n_arm_ < 1:
        raise ValueError("n_arm must be an integer >= 1")
    if not (s_pilot > 0 and math.isfinite(s_pilot)):
        raise ValueError("s_pilot > 0 required")
    cap = sqrt_k_allocation(int(n_arm_), k)
    z = float(_st.norm.ppf(1.0 - alpha / (2.0 * k)) + _st.norm.ppf(power))
    c = z * z * (PLAN_S_UP * s_pilot / delta) ** 2
    need = 1.0 / c - 1.0 / n_arm_
    if need <= 0:
        return cap
    n_c = int(math.ceil(1.0 / need - 1e-12)) + 1
    return int(min(cap, max(N_ARM_FLOOR, n_c)))


def u_ens_with_vm(u_ens: float, dof_ens: float, var_vm: float, df_vm: float,
                  n_a: int, n_b: int) -> tuple[float, float]:
    """D-L.4: a VM effect (one_way_anova on the baseline members, p <
    VM_ANOVA_ALPHA) adds its between-VM variance component to each group
    mean: u_ens^2 + var_vm (1/n_a + 1/n_b), with the Welch-Satterthwaite
    dof of the two terms (the VM term at the ANOVA's between dof).
    Conservative: the within-group SDs already carry part of it."""
    if not all(math.isfinite(v) and v >= 0 for v in (u_ens, var_vm)) or n_a < 1 or n_b < 1:
        raise ValueError("u_ens, var_vm >= 0 finite and n >= 1 required")
    u_vm = math.sqrt(var_vm * (1.0 / n_a + 1.0 / n_b))
    u = combine_quadrature(u_ens, u_vm)
    return u, welch_satterthwaite([u_ens, u_vm], [dof_ens, max(df_vm, 1e-12)])


def precision_extension_groups(fired: dict) -> set:
    """D-L.4 / D-L.3: the one extension round adds PRECISION_EXTRA_MEMBERS
    to every group of every contrast whose precision_extension fired, ONCE
    per group (a shared baseline or comparator gets +2, not +2 per
    contrast). ``fired`` maps contrast id -> (group_a, group_b, fired)."""
    out = set()
    for cid, (ga, gb, f) in fired.items():
        if f:
            out.update((ga, gb))
    return out


def heat_model_dependent(delta_hstar: float, delta_alt: float, U95: float) -> bool:
    """D-A(4): a confirmatory delta is labelled heat-model-dependent
    iff the alternative heat input moves it by more than max(0.005 K,
    0.25 |delta_H*|, U95 of the H* delta)."""
    if not all(math.isfinite(v) for v in (delta_hstar, delta_alt, U95)) or U95 < 0:
        raise ValueError("finite deltas and U95 >= 0 required")
    return abs(delta_alt - delta_hstar) > max(HEAT_LABEL_ABS_K, HEAT_LABEL_REL * abs(delta_hstar), U95)


def mechanism_causal(shares: Iterable[float]) -> bool:
    """NEW-20: causal wording iff the phase-sampler attribution explains >=
    50 % of the delta in EVERY member."""
    sh = [float(x) for x in shares]
    if not sh or not all(math.isfinite(x) for x in sh):
        raise ValueError("need finite shares")
    return all(x >= MECHANISM_SHARE for x in sh)


def headspace_triggered(ceiling_share: float) -> bool:
    """The conditional 20 mm headspace runs fire iff the mechanism attributes
    >= 50 % of a confirmatory delta to the ceiling layer."""
    if not math.isfinite(ceiling_share):
        raise ValueError("finite share required")
    return ceiling_share >= HEADSPACE_SHARE


# ---- H_E front test (NEW-13.2) ---------------------------------------------

@dataclass(frozen=True)
class FrontPoint:
    name: str
    pulsed: bool
    E: float                 # E_diss (J): pump_energy's conservative evaluation (D-P.3)
    tier_ms: float           # the run's Delta t (a Courant tier differs from Delta t*)


@dataclass(frozen=True)
class PairContrast:
    diff: float              # Delta T_cell(S) - Delta T_cell(P), K (bridged across tiers)
    u_c: float
    dof_eff: float
    U95: float
    p: float                 # two-sided p of diff / u_c at dof_eff


def bridge_offset(at_tier_p, at_tier_s) -> tuple[float, float, float]:
    """NEW-13.2: the Courant-tier bridge of a steady point S run at another
    Delta t than the pulsed point P: the bridging steady pair (one steady
    flow at both tiers) gives offset = mean(at P's tier) - mean(at S's
    tier), with its Welch standard error and dof, added to Delta T(S)."""
    w = welch_t(at_tier_p, at_tier_s)
    return float(w.diff), float(w.se), float(w.dof)


def pair_contrast(s_members, p_members, u_num=0.0, bridge=None) -> PairContrast:
    """Delta T_cell(S) - Delta T_cell(P) of one H_E pair: Welch u_ens, u_num
    (standard, of the pair's difference) at infinite dof, and, across Delta-t
    tiers, the bridge (offset, se, dof) of bridge_offset: its offset added
    to the estimate and its se as a component with its own dof."""
    w = welch_t(s_members, p_members)
    un = _u_total(u_num)
    comps, dofs = [w.se, un], [w.dof, INF]
    d = w.diff
    if bridge is not None:
        off, se_b, dof_b = bridge
        d += float(off)
        comps.append(float(se_b))
        dofs.append(float(dof_b))
    u_c = combine_quadrature(*comps)
    if u_c == 0:
        raise ValueError("pair contrast: zero uncertainty")
    nu = welch_satterthwaite(comps, dofs)
    k = _t_ppf(0.975, nu)
    p = min(1.0, 2.0 * _t_sf(abs(d) / u_c, nu))
    return PairContrast(float(d), u_c, nu, k * u_c, p)


def he_test(points: Sequence[FrontPoint], energy_u95: dict, contrasts: dict,
            alpha: float = FAMILY_ALPHA) -> dict:
    """NEW-13.2: H_E = "no pulsed computed point dominates the steady front
    on the primary energy basis", on Delta T_cell (primary window) of
    computed points only (Delta T(Q) never interpolated).

    ``energy_u95[(S, P)]``: U95 of E_diss(S) - E_diss(P)
    (pump_energy.he_energy_u95: max(2 u, 2 % of the larger energy));
    ``contrasts[(S, P)]``: PairContrast of Delta T(S) - Delta T(P) (tier-
    bridged where S and P ran at different Delta t), or None when the pair
    is not comparable (different tiers and no bridging pair): such a pair
    is never a win, and as a rival it keeps H_E. Both mappings need an
    entry for EVERY steady-pulsed pair.

    (1) Candidate pairs: E_diss(P) < E_diss(S) (point estimates); Holm at
        alpha over their two-sided p. P beats S iff E(S) - E(P) >
        U95_E(S, P), diff > 0, the Holm-adjusted p <= alpha and diff >
        U95.
    (2) P dominates the front iff it beats some S AND every steady point S'
        not resolved more expensive than P (E(S') - E(P) <= U95_E(S', P):
        cheaper, equal or unresolved) is resolved worse than P, diff(S', P)
        > U95(S', P) (intersection-union at the unadjusted 95 % level: a
        cheaper steady point at or below P's Delta T_cell within its U95
        keeps H_E). Steady Q0 is always such an S'.
    H_E is rejected iff some P dominates the front."""
    names = [pt.name for pt in points]
    if len(set(names)) != len(names):
        raise ValueError("front points need unique names")
    steady = [pt for pt in points if not pt.pulsed]
    pulsed = [pt for pt in points if pt.pulsed]
    if not steady or not pulsed:
        raise ValueError("need steady and pulsed points")
    for S in steady:
        for P in pulsed:
            if (S.name, P.name) not in energy_u95 or (S.name, P.name) not in contrasts:
                raise ValueError(f"missing pair ({S.name}, {P.name})")
    cand = [(S, P) for P in pulsed for S in steady
            if P.E < S.E and contrasts[(S.name, P.name)] is not None]
    beats = {}
    padj = {}
    if cand:
        h = holm([contrasts[(S.name, P.name)].p for S, P in cand], alpha)
        for (S, P), pa in zip(cand, h.adjusted):
            c = contrasts[(S.name, P.name)]
            padj[(S.name, P.name)] = float(pa)
            if (S.E - P.E > energy_u95[(S.name, P.name)] and c.diff > 0
                    and pa <= alpha and c.diff > c.U95):
                beats.setdefault(P.name, []).append(S.name)
    front = {}
    for P in pulsed:
        if P.name not in beats:
            continue
        rivals = [S for S in steady if S.E - P.E <= energy_u95[(S.name, P.name)]]
        kept_by = [S.name for S in rivals
                   if contrasts[(S.name, P.name)] is None
                   or not contrasts[(S.name, P.name)].diff > contrasts[(S.name, P.name)].U95]
        front[P.name] = dict(beats=sorted(beats[P.name]), rivals=[S.name for S in rivals],
                             kept_by=kept_by, dominates=not kept_by)
    rejected = any(v["dominates"] for v in front.values())
    return dict(rejected=bool(rejected), candidates=[(S.name, P.name) for S, P in cand],
                p_adj=padj, front=front)


# --------------------------------------------------------------------------
# numerical uncertainty of a contrast and the class noise reference
# (PROJECT.md §2.5 D-F.6, D-L.5; review of prep/2.10 part 3, MUST-FIX 2,
# SHOULD-FIX 5, 9; 2026-09-29)
# --------------------------------------------------------------------------

STAGE1_BOUND_BY = {"C2p": "C1", "C5": "C3", "Q174": "steady"}   # D-L.5: the bounding class
S_REF_MIN_MEMBERS = 4        # D-L.5: a stage-1 class with fewer members takes max(own, bounding) s
GRID_LEVELS = ("coarse", "medium", "fine", "xfine")
UNUM_FLOOR_SE_MULT = 2.0     # §7.1 / D-C: the u_num floor is 2 SE of the level (rung) difference


def class_s_ref(stage1: dict, cls: str) -> tuple[float, int, str]:
    """D-L.5 (amends D-L.2's s_ref): the noise reference of a class and one
    quantity = the stage-1 s of the class ITSELF (not a pool over classes:
    steady and pulsed noise differ by an order of magnitude, GDT1 0.6 mK
    against the 20 mK planning value), with its dof. A class with fewer
    than 4 stage-1 members (C2', C5, steady Q174: N2) takes the larger of
    its own s and that of its bounding class (C1, C3, steady Q90), with
    the dof of the one taken. ``stage1`` maps class -> member values
    (classes steady, C1, C3, C2p, C5, Q174). Returns (s, dof, how)."""
    if cls not in stage1:
        raise ValueError(f"class_s_ref: no stage-1 members of {cls!r}")
    own = _sample(stage1[cls], cls, 2)
    s, d, how = float(own.std(ddof=1)), int(own.size - 1), f"own stage-1 s ({cls}, N{own.size})"
    if own.size < S_REF_MIN_MEMBERS:
        b = STAGE1_BOUND_BY.get(cls)
        if b is None or b not in stage1:
            raise ValueError(f"class_s_ref: {cls} has N{own.size} < {S_REF_MIN_MEMBERS} "
                             "and no bounding class")
        bm = _sample(stage1[b], b, 2)
        sb = float(bm.std(ddof=1))
        if sb > s:
            s, d, how = sb, int(bm.size - 1), f"bounding class {b} (N{bm.size}) > own {cls}"
    return s, d, how


def _s_use(members: np.ndarray, s_ref) -> float:
    """max(s_ref, the group's own SD) (D-L.2 / D-L.5); s_ref may be None."""
    own = float(members.std(ddof=1)) if members.size >= 2 else None
    c = [x for x in (None if s_ref is None else float(s_ref[0] if isinstance(s_ref, (tuple, list))
                                                        else s_ref), own) if x is not None]
    if not c:
        raise ValueError("a single member needs a class s_ref")
    return max(c)


def _group_mean_se(members, s_ref) -> tuple[float, float]:
    m = _sample(members, "group", 1)
    return float(m.mean()), _s_use(m, s_ref) / math.sqrt(m.size)


def grid_delta_u(arm: dict, cells: dict, rise: float, s_ref_arm, ref: dict | None = None,
                 s_ref_ref=None, production: str = "medium") -> dict:
    """D-F.6 (1): u_grid of a delta (or, with ``ref`` None, of an absolute
    quantity) from its per-level series. ``arm`` / ``ref`` map grid level ->
    member values; the levels both hold (>= 3; C1, C3 and steady Q90 have
    four) enter in GRID_LEVELS order. Delta_l = mean(arm_l) - mean(ref_l),
    se_l = sqrt(s_use,a^2 / n_a + s_use,r^2 / n_r) with s_use = max(the
    class s_ref, the group's own SD) (D-L.5). The series goes to
    gci2.grid_uncertainty (Eca-Hoekstra on 4 levels, Celik on 3), whose
    floor max(|eps|, 2 SE_pair) keeps noise from making u_grid vanish.
    ``cells`` maps level -> cell count (h = N^(-1/3)); ``rise`` normalises
    the percentages only. Returns dict(U, u, levels, deltas, se, method) at
    the production level."""
    import gci2
    levels = [lv for lv in GRID_LEVELS if lv in arm and (ref is None or lv in ref)]
    if len(levels) < 3:
        raise ValueError(f"grid_delta_u: need >= 3 common levels (got {levels})")
    if production not in levels:
        raise ValueError(f"grid_delta_u: production level {production!r} missing")
    deltas, se = [], []
    for lv in levels:
        ma, sa = _group_mean_se(arm[lv], s_ref_arm)
        if ref is None:
            deltas.append(ma)
            se.append(sa)
        else:
            mr, sr = _group_mean_se(ref[lv], s_ref_ref)
            deltas.append(ma - mr)
            se.append(math.hypot(sa, sr))
    h = gci2.cell_size([float(cells[lv]) for lv in levels])
    r = gci2.grid_uncertainty(deltas, h, rise, se=se)
    i = levels.index(production)
    U = float(r["U"][i])
    return dict(U=U, u=U / gci2.K_NUM, levels=levels, deltas=deltas, se=se, method=r["method"])


def time_delta_u(arm: dict, s_ref_arm, ref: dict | None = None, s_ref_ref=None) -> dict:
    """D-F.6 (2): u_time of a delta (or, with ``ref`` None, of an absolute
    quantity) from the FREEZE ladder: ``arm`` / ``ref`` map 'star' (the
    Delta t* rung: the stage-1 members), 'half' (Delta t*/2, FREEZE N2) and
    optionally 'quarter' (Delta t*/4) -> member values of that class
    (steady Q90 for the reference of an F-M1 delta). d = Delta(half) -
    Delta(star), SE_pair = hypot(se_star, se_half) (se as grid_delta_u).
    U_time = Fs (4/3) max(|d|, 2 SE_pair), Fs = 1.25 if |d| < 2 SE_pair
    (gci2's noise branch: formal order) or a quarter rung gives a monotone
    observed order in [0.5, 2.1) with its increment beyond 2 SE, else 3 -
    the U of recipe_rules.freeze_pair_ok WITH the §7.1 floor, so noise
    never makes u_time vanish. u_time = U_time / 1.15."""
    import gci2
    F = gci2.richardson_factor(2.0, 2.0)

    def delta(rung):
        ma, sa = _group_mean_se(arm[rung], s_ref_arm)
        if ref is None:
            return ma, sa
        mr, sr = _group_mean_se(ref[rung], s_ref_ref)
        return ma - mr, math.hypot(sa, sr)

    for rung in ("star", "half"):
        if rung not in arm or (ref is not None and rung not in ref):
            raise ValueError(f"time_delta_u: rung {rung!r} missing")
    d_s, se_s = delta("star")
    d_h, se_h = delta("half")
    e32, se32 = d_h - d_s, math.hypot(se_s, se_h)
    fs, why = gci2.FS_BAD, "two rungs: no observed order"
    if abs(e32) < UNUM_FLOOR_SE_MULT * se32:
        fs, why = gci2.FS_GOOD, "difference within 2 SE_pair (noise branch, formal order)"
    elif "quarter" in arm and (ref is None or "quarter" in ref):
        d_q, se_q = delta("quarter")
        e21, se21 = d_q - d_h, math.hypot(se_h, se_q)
        if abs(e21) < UNUM_FLOOR_SE_MULT * se21:
            why = "finer increment within the noise: order not observable"
        elif e32 / e21 <= 0:
            why = "oscillatory"
        else:
            p = math.log(e32 / e21) / math.log(2.0)
            ok = gci2.P_LOW <= p < 2.0 + gci2.P_FS_MARGIN
            fs, why = (gci2.FS_GOOD if ok else gci2.FS_BAD), f"observed order {p:.3g}"
    U = fs * F * max(abs(e32), UNUM_FLOOR_SE_MULT * se32)
    return dict(U=U, u=U / gci2.K_NUM, diff=e32, se_pair=se32, Fs=fs, fs_note=why)


def family_flow_u(q: float, table: dict) -> tuple[float, str]:
    """D-F.6 (3): the absolute numerical u of a steady flow without its own
    grid or Delta-t study (a comparator Q_eq off the family flows, a front
    point): the family flows that bracket it, the larger of their u; below
    or above every family flow, the nearest one's. ``table`` maps family
    flow (mL/min) -> u of the quantity."""
    if not table:
        raise ValueError("family_flow_u: empty table")
    fl = sorted(float(x) for x in table)
    u = {float(k): float(v) for k, v in table.items()}
    q = float(q)
    for f in fl:
        if abs(f - q) < 1e-9:
            return u[f], f"family flow {f:g}"
    lo = [f for f in fl if f < q]
    hi = [f for f in fl if f > q]
    if lo and hi:
        a, b = lo[-1], hi[0]
        return max(u[a], u[b]), f"bracketed by {a:g} and {b:g}: the larger"
    f = hi[0] if hi else lo[-1]
    return u[f], f"nearest family flow {f:g}"


def delta_unum(u_grid, u_time, u_iter) -> dict:
    """D-F.6: u_num = sqrt(u_grid^2 + u_time^2 + u_iter^2) of a contrast
    (standard uncertainties; a dict argument contributes its 'u'), the
    input of contrast_u95 / pair_contrast. Which estimator each component
    comes from is fixed by D-F.6:
      F-M1 delta (arm - steady Q90): u_grid = grid_delta_u on the family
        levels both hold (C1, C3: 4; C2', C5: 3); u_time = time_delta_u on
        the FREEZE ladders of the arm and steady Q90 (C2' takes C1's value,
        C5 C3's: the NEW-2 bounding rule).
      F-E1 delta (arm - steady Q_eq): u_grid by grid_delta_u against the
        family's comparator when Q_eq is a family flow, else the arm's
        absolute u_grid (+) family_flow_u(Q_eq) in quadrature; u_time = the
        arm's absolute u_time (+) steady Q90's absolute u_time (the only
        steady Delta-t study) in quadrature - no cancellation is assumed.
      H_E pairs and every other contrast: the absolute rule as F-E1.
      u_iter: recipe_rules.u_iter (one value for every delta, D-J.3/D-J.5)."""
    parts = []
    for name, x in (("u_grid", u_grid), ("u_time", u_time), ("u_iter", u_iter)):
        v = float(x["u"]) if isinstance(x, dict) else float(x)
        if not (math.isfinite(v) and v >= 0):
            raise ValueError(f"delta_unum: {name} must be finite and >= 0")
        parts.append(v)
    return dict(u_grid=parts[0], u_time=parts[1], u_iter=parts[2],
                u_num=combine_quadrature(*parts))
