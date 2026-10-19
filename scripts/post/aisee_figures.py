#!/usr/bin/env python3
"""Figures of the paper (D-AIS.1), called by aisee_results.py.

fig_geometry.pdf   side view (x-z) of the housing, inlet channel and riser,
                   with a top view (x-y) of the 2 x 3 array as an inset;
                   dimensions from scripts/mesh/k0_blockmesh.py (GEOM,
                   GEOM_HOUSING, port_geometry) and properties.yaml
fig_pareto.pdf     T_max,pk (a) and dT_cell at the end of the discharge (b)
                   against the pumping energy: the steady curve and the
                   pulsed arms (matched mean flow: filled; reduced mean flow
                   aimed at steady Q90's energy: open)
fig_histories.pdf  T_max(t) relative to steady Q90 and dT_cell(t) over the
                   discharge, steady Q90 against the pulsed arms at Q90
fig_waveforms.pdf  the imposed inlet flow Q(t)/Qbar of C1, C2 and C3 over 2 s (case model)
fig_matched.pdf    every pulsed case against steady flow of equal pumping
                   energy over the discharge (equal_energy_steady: the steady
                   histories interpolated linearly in log E_pump, the rule of
                   the results table's delta columns): T_max(t) difference and
                   relative change of dT_cell(t); text at 8 pt (IEEE guidance)

Okabe-Ito hues (validated, scripts in the dataviz skill; adjacent CVD pair
in the 6-8 band), so every series also has its own marker and line style
(grayscale print) and the values are in the paper's table. A case that
failed a gate carries the mark the results table gives it (rule 4: flagged, never
plotted silently): dagger = hydraulic_closure[pre], double dagger =
max_principle, asterisk = any other gate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Polygon, Rectangle  # noqa: E402

REPO = Path(__file__).resolve().parents[2]

COLOR = {"steady": "#0072B2", "C1": "#E69F00", "C2p": "#009E73", "C3": "#CC79A7"}
MARKER = {"steady": "o", "C1": "s", "C2p": "D", "C3": "^"}
LS = {"steady": "-", "C1": "--", "C2p": "-.", "C3": ":"}
# case names as in the paper's text and tables (C2p is the paper's C2)
LABEL = {"steady": "Steady", "C1": "C1: square 1 Hz, $D$ = 0.25",
         "C2p": "C2: square 0.5 Hz, $D$ = 0.5", "C3": "C3: sine 2 Hz, $A$ = 0.6"}
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e3e3e3"
NOISE_TMAX_K = 0.004      # GDT1 dt and seed pairs, volume T_max,peak (aisee_results.py)
COL_W = 3.5               # IEEE column width, inches
FIG1_W = 3.45             # fig_geometry width (printed at its natural size, 8 pt labels)

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.2, "pdf.fonttype": 42,
    # Times-like text and math, to match the IEEE body font (ships with matplotlib)
    "font.family": "serif", "font.serif": ["STIXGeneral"], "mathtext.fontset": "stix",
})


def _flags(r) -> str:
    """The results table's marks for a failed gate: dagger the port-work residual
    (hydraulic_closure[pre]), double dagger the liquid bound (max_principle),
    an asterisk for any other gate, so nothing failed is plotted unmarked."""
    failed = [g for g in (r.get("gates_failed") or "").split(";") if g]
    mark = ("†" if "hydraulic_closure[pre]" in failed else "") + \
           ("‡" if "max_principle" in failed else "")
    other = [g for g in failed if g not in ("hydraulic_closure[pre]", "max_principle")]
    return mark + ("*" if other else "")


def _tag(r):
    """Point label: Qbar with the case's marks."""
    return f"{r['Qbar_mL_min']:g}" + _flags(r)


# ------------------------------------------------------------ fig_geometry
def _geometry():
    sys.path.insert(0, str(REPO / "scripts" / "mesh"))
    import k0_blockmesh as K
    g = K.load_geometry()
    g.update(K.port_geometry(g, 1))
    return g


def fig_geometry(path):
    """Side view (x-z, to scale) with the channel and riser cut short, and a
    top view (x-y) inset at reduced scale."""
    g = _geometry()
    D, H, gap = g["D"], g["H"], g["gap"]
    x0, x1, zt, xp, hc = g["x_min"], g["x_max"], g["z_top"], g["x_port"], g["h_ch"]
    yh, L = g["y_half"], g["L"]
    pitch = D + gap
    xs = [(-(g["cols"] - 1) / 2 + c) * pitch for c in range(g["cols"])]
    ys = [(-(g["rows"] - 1) / 2 + r) * pitch for r in range(g["rows"])]
    xl, ztop = x0 - 34.0, zt + 16.0       # cut ends of the channel and the riser
    xlim, ylim = (xl - 3, x1 + 13), (-hc - 2.5, ztop + 2.5)
    w = FIG1_W
    fig = plt.figure(figsize=(w, w * (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_aspect("equal")
    ax.axis("off")
    liquid, wall, cellc = "#dceaf5", INK, "#9a9a9a"
    kw = dict(lw=0.8, ec=wall)
    # liquid regions: housing, inlet channel (with the elbow under the slot), riser
    ax.add_patch(Rectangle((x0, 0), x1 - x0, zt, fc=liquid, **kw))
    ax.add_patch(Rectangle((xl, -hc), x0 + xp - xl, hc, fc=liquid, **kw))
    # the riser ends on its cut line (drawn below), so no wall stub shows above it
    _k = 4.0 / (xp + 4.0)   # slope of the riser cut line
    ax.add_patch(Polygon([(x1 - xp, zt), (x1, zt), (x1, ztop - 2 + _k * (xp + 2)),
                          (x1 - xp, ztop - 2 + _k * 2)], closed=True, fc=liquid, **kw))
    # open the full-width floor and lid slots
    ax.add_patch(Rectangle((x0 + 0.4, -0.6), xp - 0.8, 1.2, fc=liquid, ec="none", zorder=3))
    ax.add_patch(Rectangle((x1 - xp + 0.4, zt - 0.6), xp - 0.8, 1.2, fc=liquid, ec="none",
                           zorder=3))
    # cut marks: the straight lengths L continue beyond
    for xa, za, xb, zb in ((xl - 2, -hc - 2, xl + 2, 2), (x1 - xp - 2, ztop - 2, x1 + 2, ztop + 2)):
        ax.plot([xa, xb], [za, zb], color="white", lw=2.6, zorder=4)
        ax.plot([xa, xb], [za, zb], color=wall, lw=0.6, zorder=5)
    # cells (one row seen from the side)
    for xc in xs:
        ax.add_patch(Rectangle((xc - D / 2, 0), D, H, fc=cellc, ec=INK, lw=0.6, zorder=2))
    fs = 8
    tx = dict(fontsize=fs, color=INK, zorder=6)
    # flow arrows
    ar = dict(arrowstyle="-|>", color=COLOR["steady"], lw=0.9, mutation_scale=6)
    ax.annotate("", (xl + 24, -hc / 2), (xl + 4, -hc / 2), arrowprops=ar, zorder=6)
    ax.annotate("", (x0 + xp / 2, 9), (x0 + xp / 2, -3), arrowprops=ar, zorder=6)
    ax.annotate("", (x1 - xp / 2, ztop - 3), (x1 - xp / 2, zt - 8), arrowprops=ar, zorder=6)
    ax.text(xl + 1, 2, "inlet\nchannel", va="bottom", ha="left", linespacing=1.0, **tx)
    ax.text(x0 + xp + 2, 1.5, "floor slot", va="bottom", ha="left", **tx)
    ax.text(x1 - xp - 4, zt + 1.2, "lid slot", va="bottom", ha="right", **tx)
    ax.text(x1 - xp - 4, zt + 8.0, "riser", va="bottom", ha="right", **tx)
    # headspace above the cells: a dimension arrow beside the last cell, labelled below it
    xh = xs[-1] + D / 2 + 3.0
    ax.annotate("", (xh, H), (xh, zt), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.5,
                                                       mutation_scale=4, shrinkA=0, shrinkB=0))
    ax.text(xh - 1.0, H - 0.8, f"{zt - H:.0f} mm headspace", va="top", ha="left", fontsize=fs,
            color=INK, zorder=6)
    ax.text(xs[0], H / 2, "cell", va="center", ha="center", fontsize=fs, color="white",
            zorder=6, rotation=90)
    # height
    dim = dict(arrowstyle="<->", color=MUTED, lw=0.5, mutation_scale=5, shrinkA=0, shrinkB=0)
    ax.annotate("", (x1 + 5, 0), (x1 + 5, zt), arrowprops=dim)
    ax.text(x1 + 7, zt / 2, f"{zt:.0f} mm", va="center", ha="left", rotation=90, **tx)
    # gravity
    ax.annotate("", (xl + 6, 36), (xl + 6, 60), arrowprops=dict(arrowstyle="-|>", color=INK,
                                                               lw=0.7, mutation_scale=6))
    ax.text(xl + 9, 48, r"$g$", va="center", ha="left", **tx)
    # top view inset (x-y), reduced scale, above the housing
    s, ox, oz = 0.29, -30.0, zt + 9.2
    ax.add_patch(Rectangle((ox + s * x0, oz - s * yh), s * (x1 - x0), 2 * s * yh, fc=liquid,
                           lw=0.6, ec=wall, zorder=3))
    for xc in xs:
        for yc in ys:
            ax.add_patch(Circle((ox + s * xc, oz + s * yc), s * D / 2, fc=cellc, ec=INK,
                                lw=0.4, zorder=4))
    for xa in (x0, x1 - xp):      # the full-width slots
        ax.add_patch(Rectangle((ox + s * xa, oz - s * yh), s * xp, 2 * s * yh, fc="white",
                               ec=wall, lw=0.4, hatch="////", zorder=4))
    ax.annotate("", (ox + s * x1 + 9, oz), (ox + s * x1 + 1, oz),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.6, mutation_scale=5))
    ax.text(ox + s * x1 + 10, oz, r"$x$", va="center", ha="left", **tx)
    ax.text(ox + s * x0 - 2, oz, "top view, reduced\nscale; slots hatched", va="center",
            ha="right", linespacing=0.95, fontsize=fs, color=INK, zorder=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.savefig(path)
    plt.close(fig)


# -------------------------------------------------------------- fig_pareto
# label offsets (points) of the pulsed symbols in panel (a), clear of their neighbours
PULSED_LABEL_AT = {"C3_Q82.9": (-8, 5), "C3_Q90": (5, 3), "C2p_Q64.2": (6, -2)}
# the steady 135 label of panel (b) sits above its own dot, clear of the C3 and C2 symbols
STEADY_LABEL_AT = {("Tmax_peak_C", 90): (3, 3), ("dTcell_end_mK", 360): (-26, 4),
                   ("dTcell_end_mK", 135): (-6, 7)}


def _main_rows(rows):
    return [r for r in rows if r["role"] in ("curve", "mean-flow", "matched-energy")]


def fig_pareto(rows, path):
    rows = [dict(r, dTcell_end_mK=1e3 * r["dTcell_end_K"]) for r in _main_rows(rows)]
    # one column, panels stacked on a shared energy axis
    fig, axs = plt.subplots(2, 1, figsize=(COL_W, 2.8), sharex=True)
    for ax, key, lab in ((axs[0], "Tmax_peak_C", r"$T_{\max,\mathrm{pk}}$ (°C)"),
                         (axs[1], "dTcell_end_mK", r"$\Delta T_{\mathrm{cell}}$ (mK)")):
        st = sorted([r for r in rows if r["arm"] == "steady"], key=lambda r: r["E_pump_mJ"])
        if st:
            x = [r["E_pump_mJ"] for r in st]
            y = [r[key] for r in st]
            ax.plot(x, y, color=COLOR["steady"], ls="-", marker="o", ms=4.5,
                    mfc=COLOR["steady"], mec="white", mew=0.6, label="Steady", zorder=3)
            # the run-to-run noise (NOISE_TMAX_K, 0.4 mK in dT_cell) is smaller than the
            # symbols at this scale: stated in the caption, not drawn
            for r in st:
                at = STEADY_LABEL_AT.get((key, r["Qbar_mL_min"]),
                                         (3, 3) if key == "Tmax_peak_C" or r["Qbar_mL_min"] < 100
                                         else (-18, 5))
                ax.annotate(_tag(r), (r["E_pump_mJ"], r[key]), textcoords="offset points",
                            fontsize=8, color=INK, xytext=at)
        for arm in ("C1", "C2p", "C3"):
            for role, filled in (("mean-flow", True), ("matched-energy", False)):
                pts = [r for r in rows if r["arm"] == arm and r["role"] == role]
                if not pts:
                    continue
                ax.plot([r["E_pump_mJ"] for r in pts], [r[key] for r in pts],
                        ls="none", marker=MARKER[arm], ms=5.5 if key == "Tmax_peak_C" else 4.2,
                        mfc=COLOR[arm] if filled else "white", mec=COLOR[arm], mew=1.1,
                        label=f"{LABEL[arm]}, " + (r"$\bar Q$ = 90 mL/min" if filled else
                                                   r"reduced $\bar Q$"), zorder=4)
                # pulsed points carry their mean flow in panel (a) only (panel (b) is
                # crowded near 0.3-0.4 mJ and the legend identifies each point); a
                # failed check is marked in both panels
                for r in pts:
                    if key == "Tmax_peak_C":
                        ax.annotate(_tag(r), (r["E_pump_mJ"], r[key]), textcoords="offset points",
                                    fontsize=8, color=INK,
                                    xytext=PULSED_LABEL_AT.get(r["key"], (4, -9)))
                    elif _flags(r):
                        ax.annotate(_flags(r), (r["E_pump_mJ"], r[key]), textcoords="offset points",
                                    fontsize=8, color=INK, xytext=(5, -3))
        ax.set_xscale("log")
        from matplotlib.ticker import LogLocator, FuncFormatter
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
        ax.set_ylabel(lab)
    axs[1].set_xlabel(r"Pumping energy over the 900 s discharge, $E_{\mathrm{pump}}$ (mJ)")
    # fixed bottom: the C1 90 symbol and its mark stay clear of the spine
    axs[1].set_ylim(50, 420)
    # panel titles inside the axes, in corners clear of the data
    for ax, t, x, ha in ((axs[0], "(a) Peak cell temperature", 0.99, "right"),
                         (axs[1], "(b) Final spread", 0.5, "left")):
        ax.text(x, 0.96, t, transform=ax.transAxes, ha=ha, va="top", fontsize=8,
                color=INK, bbox=dict(fc="white", ec="none", pad=0.5))
    # one legend entry per arm; filled (mean flow 90 mL/min) and open (reduced mean
    # flow) symbols are explained in the paper's caption
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=COLOR["steady"], marker="o", ms=4.5, mec="white", mew=0.6,
                 label="Steady")]
    # open handles: every arm has an open (reduced mean flow) point, C2 no filled one;
    # filled/open is explained in the caption
    hs += [Line2D([], [], ls="none", marker=MARKER[a], ms=5.5, mfc="white", mec=COLOR[a],
                  mew=1.1, label=LABEL[a]) for a in ("C1", "C2p", "C3")]
    # legend in the empty upper left of panel (b)
    axs[1].legend(handles=hs, loc="upper left", frameon=False, fontsize=8,
                  handlelength=1.3, labelspacing=0.15, borderaxespad=0.2,
                  handletextpad=0.4)
    fig.tight_layout(pad=0.3, h_pad=0.3)
    fig.savefig(path)
    plt.close(fig)


# ----------------------------------------------------------- fig_histories
def fig_histories(rows, paths, path):
    import aisee_results as AR
    # the mechanism figure: the same mean flow, Q90 (matched-energy arms: fig_pareto, fig_matched)
    keys = [r["key"] for r in _main_rows(rows)
            if r["key"] in paths and (r["role"] == "mean-flow" or r["key"] == "steady_Q90")]
    if "steady_Q90" not in keys:
        return
    t_pre = 30.0
    ts, smax_s, _, _ = AR.series(paths["steady_Q90"])
    fig, axs = plt.subplots(2, 1, figsize=(COL_W, 1.8), sharex=True)
    for k in ["steady_Q90"] + [k for k in keys if k != "steady_Q90"]:
        r = next(x for x in rows if x["key"] == k)
        t, smax, dtc, _ = AR.series(paths[k])
        m = t > t_pre
        arm = r["arm"]
        style = dict(color=COLOR[arm], ls=LS[arm], lw=1.2)
        # short legend labels (waveform parameters are in the cases table), clear of the curves
        lab = {"steady": "Steady", "C1": "C1: square 1 Hz", "C3": "C3: sine 2 Hz",
               "C2p": "C2: square 0.5 Hz"}[arm] + _flags(r)
        # top: T_max relative to steady Q90 at the same time since heat-on
        rel = smax[m] - np.interp(t[m], ts, smax_s)
        axs[0].plot(t[m] - t_pre, rel, label=lab, **style)
        axs[1].plot(t[m] - t_pre, 1e3 * dtc[m], **style)
    # change from steady Q90 at the same time (caption); "delta" is the results table's equal-energy symbol
    # top label split so only "(K)" sits beside the axis, clear of the bottom label
    axs[0].set_ylabel(r"$T_{\max}$ change" + "\n(K)", linespacing=1.0)
    axs[1].set_ylabel(r"$\Delta T_{\mathrm{cell}}$ (mK)")
    axs[1].set_xlabel("Time since heat-on (s)")
    axs[1].set_xlim(0, 900)
    h, lbl = axs[0].get_legend_handles_labels()
    # one legend row above the panels, clear of the curves
    fig.legend(h, lbl, frameon=False, loc="upper center", ncol=3, fontsize=8,
               handlelength=2.0, columnspacing=1.2, bbox_to_anchor=(0.55, 1.0),
               borderaxespad=0.1)
    fig.tight_layout(pad=0.3, h_pad=0.4, rect=(0, 0, 1, 0.88))
    fig.align_ylabels(axs)
    fig.savefig(path)
    plt.close(fig)


# ------------------------------------------------------------- fig_waveforms
def fig_waveforms(rows, paths, path):
    """The imposed inlet flow Q(t)/Qbar of each pulsed arm over 2 s, from the
    case's own switches (case.env) through the case model
    (housing_physics.waveform_q: the solver's inlet table, 50 ms valve ramps),
    without the start-up ramp; the legend gives the measured <Q^2>/Qbar^2."""
    import re
    sys.path.insert(0, str(REPO / "scripts" / "model"))
    import housing_physics as HP
    keys = {r["arm"]: r["key"] for r in _main_rows(rows)
            if r["role"] == "matched-energy" and r["key"] in paths}
    if not keys:
        return
    t = np.linspace(0.0, 2.0, 8001)
    fig, ax = plt.subplots(figsize=(COL_W, 1.75))
    ax.axhline(1.0, color=COLOR["steady"], lw=1.2, label=r"Steady (1)")
    for arm in ("C1", "C2p", "C3"):
        if arm not in keys:
            continue
        env = dict(re.findall(r"^(?:export\s+)?(\w+)=(.*)$",
                              (Path(paths[keys[arm]]) / "case.env").read_text(), re.M))
        cfg = HP.Config.from_env(env)
        q = HP.waveform_q(cfg, t, ramp=False)
        qs = np.mean(HP.waveform_q(cfg, np.linspace(0.0, 4.0, 40001)[:-1], ramp=False) ** 2) / \
            (cfg.qbar_ml_min / 6e7) ** 2
        lab = {"C1": "C1: square 1 Hz", "C2p": "C2: square 0.5 Hz", "C3": "C3: sine 2 Hz"}[arm]
        ax.plot(t, q / (cfg.qbar_ml_min / 6e7), color=COLOR[arm], ls=LS[arm], lw=1.2,
                label=f"{lab} ({qs:.2f})")
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.2, 4.4)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"$Q/\bar Q$")
    # legend in two rows above the axes, clear of the pulses
    fig.legend(*ax.get_legend_handles_labels(), loc="upper center", ncol=2, frameon=False,
               fontsize=8, handlelength=2.0, columnspacing=1.0, labelspacing=0.15,
               bbox_to_anchor=(0.55, 1.0), borderaxespad=0.1)
    fig.tight_layout(pad=0.3, rect=(0, 0, 1, 0.83))
    fig.savefig(path)
    plt.close(fig)


# ------------------------------------------------------------- fig_matched
def equal_energy_steady(rows, series_of, e):
    """The steady history at pumping energy e: the computed steady curve's
    histories interpolated linearly in log E_pump between the two bracketing
    steady flows, at each time, the rule of Table III's delta columns.
    series_of(key) -> (t, T_max(t), dT_cell(t)); returns (t, T_max, dT_cell,
    (key_lo, key_hi, weight)) on the lower steady flow's time grid."""
    st = sorted([r for r in _main_rows(rows) if r["arm"] == "steady"],
                key=lambda r: r["E_pump_mJ"])
    lo = max((r for r in st if r["E_pump_mJ"] <= e), key=lambda r: r["E_pump_mJ"])
    hi = min((r for r in st if r["E_pump_mJ"] > e), key=lambda r: r["E_pump_mJ"])
    w = (np.log(e) - np.log(lo["E_pump_mJ"])) / (np.log(hi["E_pump_mJ"]) - np.log(lo["E_pump_mJ"]))
    t, a_lo, b_lo = series_of(lo["key"])
    th, a_hi, b_hi = series_of(hi["key"])
    a = (1 - w) * a_lo + w * np.interp(t, th, a_hi)
    b = (1 - w) * b_lo + w * np.interp(t, th, b_hi)
    return t, a, b, (lo["key"], hi["key"], float(w))


def fig_matched_histories(rows, paths, path):
    """T_max(t) of every pulsed case less that of steady flow of equal pumping
    energy (equal_energy_steady), and the relative change of its spread
    dT_cell(t), over the discharge; left the reduced mean flows, right the
    pulsed runs at Q90. The spread change starts SPREAD_FROM_S after heat-on,
    when the equal-energy steady spreads exceed 14 mK (a ratio of the first
    mK-level spreads would dominate the scale)."""
    import aisee_results as AR
    t_pre = 30.0
    SPREAD_FROM_S = 60.0
    cache = {}

    def series_of(k):
        if k not in cache:
            t, smax, dtc, _ = AR.series(paths[k])
            cache[k] = (t, smax, dtc)
        return cache[k]

    main = _main_rows(rows)
    groups = [("Reduced mean flow", [r for r in main if r["role"] == "matched-energy"]),
              (r"$\bar Q$ = 90 mL/min", [r for r in main if r["role"] == "mean-flow"])]
    if any(r["key"] not in paths for _, g in groups for r in g) or not groups[0][1]:
        return
    fig, axs = plt.subplots(2, 2, figsize=(2 * COL_W + 0.16, 2.3), sharex=True, sharey="row")
    for j, (title, grp) in enumerate(groups):
        for r in sorted(grp, key=lambda r: ("C1", "C2p", "C3").index(r["arm"])):
            t, smax, dtc = series_of(r["key"])
            te, ae, be, _ = equal_energy_steady(rows, series_of, r["E_pump_mJ"])
            m = t > t_pre
            style = dict(color=COLOR[r["arm"]], ls=LS[r["arm"]], lw=1.1)
            lab = {"C1": "C1: square 1 Hz", "C2p": "C2: square 0.5 Hz",
                   "C3": "C3: sine 2 Hz"}[r["arm"]] + f", {r['Qbar_mL_min']:g} mL/min" + _flags(r)
            axs[0, j].plot(t[m] - t_pre, smax[m] - np.interp(t[m], te, ae), label=lab, **style)
            ms = t > t_pre + SPREAD_FROM_S
            axs[1, j].plot(t[ms] - t_pre, 100 * (dtc[ms] / np.interp(t[ms], te, be) - 1), **style)
        axs[0, j].set_title(title, fontsize=8, color=INK, pad=3)
        axs[0, j].legend(loc="upper left", frameon=False, fontsize=8, handlelength=2.2,
                         labelspacing=0.2, borderaxespad=0.2)
        for ax in axs[:, j]:
            ax.axhline(0, color=MUTED, lw=0.6, zorder=1)
        axs[1, j].set_xlabel("Time since heat-on (s)")
        axs[1, j].set_xlim(0, 900)
        axs[1, j].set_xticks(range(0, 901, 150))
    axs[1, 0].set_yticks(range(-100, 51, 50))
    # two-line labels, so the two rows' labels stay clear of each other
    axs[0, 0].set_ylabel(r"$T_{\max}$ difference" + "\n(K)", linespacing=1.0)
    axs[1, 0].set_ylabel(r"$\Delta T_{\mathrm{cell}}$ change" + "\n(%)", linespacing=1.0)
    fig.tight_layout(pad=0.3, h_pad=0.4, w_pad=1.2)
    fig.align_ylabels(axs[:, 0])
    fig.savefig(path)
    plt.close(fig)


def make_all(rows, paths, figdir):
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig_geometry(figdir / "fig_geometry.pdf")
    fig_pareto(rows, figdir / "fig_pareto.pdf")
    fig_waveforms(rows, paths, figdir / "fig_waveforms.pdf")
    fig_histories(rows, paths, figdir / "fig_histories.pdf")
    fig_matched_histories(rows, paths, figdir / "fig_matched.pdf")
    print(f"figures -> {figdir}")
