#!/usr/bin/env python3
"""AISEE 2027 conference figures (D-AIS.1), called by aisee_results.py.

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

Okabe-Ito hues (validated, scripts in the dataviz skill; adjacent CVD pair
in the 6-8 band), so every series also has its own marker and line style
(grayscale print) and the values are in the paper's table. A case that
failed a gate carries the mark Table I gives it (rule 4: flagged, never
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
FIG1_W = 2.95             # Fig. 1 width (printed at its natural size, 7 pt labels)

plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.2, "pdf.fonttype": 42,
    # Times-like text and math, to match the IEEE body font (ships with matplotlib)
    "font.family": "serif", "font.serif": ["STIXGeneral"], "mathtext.fontset": "stix",
})


def _flags(r) -> str:
    """Table I's marks for a failed gate: dagger the port-work residual
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


# --------------------------------------------------------------------- Fig. 1
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
    fs = 7
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
    ax.text(xs[-1] + D / 2 + 1.5, (H + zt) / 2 - 0.3, f"{zt - H:.0f} mm headspace", va="center",
            ha="left", fontsize=6.5, color=INK, zorder=6)
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
            ha="right", linespacing=0.95, fontsize=6.5, color=INK, zorder=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------- Fig. 2
# label offsets (points) of the pulsed symbols in panel (a), clear of their neighbours
PULSED_LABEL_AT = {"C3_Q82.9": (-8, 5), "C3_Q90": (5, 3), "C2p_Q64.2": (6, -2)}
# the steady 135 label of panel (b) sits above its own dot, clear of the C3 and C2 symbols
STEADY_LABEL_AT = {("Tmax_peak_C", 90): (3, 3), ("dTcell_end_mK", 360): (-26, 4),
                   ("dTcell_end_mK", 135): (-6, 7)}


def _main_rows(rows):
    return [r for r in rows if r["role"] in ("curve", "mean-flow", "matched-energy")]


def fig_pareto(rows, path):
    rows = [dict(r, dTcell_end_mK=1e3 * r["dTcell_end_K"]) for r in _main_rows(rows)]
    # one column, panels stacked on a shared energy axis (6-page limit)
    fig, axs = plt.subplots(2, 1, figsize=(COL_W, 2.2), sharex=True)
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
                            fontsize=7, color=INK, xytext=at)
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
                                    fontsize=7, color=INK,
                                    xytext=PULSED_LABEL_AT.get(r["key"], (4, -9)))
                    elif _flags(r):
                        ax.annotate(_flags(r), (r["E_pump_mJ"], r[key]), textcoords="offset points",
                                    fontsize=7, color=INK, xytext=(5, -3))
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
        ax.text(x, 0.96, t, transform=ax.transAxes, ha=ha, va="top", fontsize=7.5,
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
    axs[1].legend(handles=hs, loc="upper left", frameon=False, fontsize=6.8,
                  handlelength=1.3, labelspacing=0.15, borderaxespad=0.2,
                  handletextpad=0.4)
    fig.tight_layout(pad=0.3, h_pad=0.3)
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------- Fig. 3
def fig_histories(rows, paths, path):
    import aisee_results as AR
    # the mechanism figure: the same mean flow, Q90 (matched-energy arms are in Fig. 2)
    keys = [r["key"] for r in _main_rows(rows)
            if r["key"] in paths and (r["role"] == "mean-flow" or r["key"] == "steady_Q90")]
    if "steady_Q90" not in keys:
        return
    t_pre = 30.0
    ts, smax_s, _, _ = AR.series(paths["steady_Q90"])
    fig, axs = plt.subplots(2, 1, figsize=(COL_W, 1.5), sharex=True)
    for k in ["steady_Q90"] + [k for k in keys if k != "steady_Q90"]:
        r = next(x for x in rows if x["key"] == k)
        t, smax, dtc, _ = AR.series(paths[k])
        m = t > t_pre
        arm = r["arm"]
        style = dict(color=COLOR[arm], ls=LS[arm], lw=1.2)
        # short legend labels (waveform parameters are in Table I), clear of the curves
        lab = {"steady": "Steady", "C1": "C1: square 1 Hz", "C3": "C3: sine 2 Hz",
               "C2p": "C2: square 0.5 Hz"}[arm] + _flags(r)
        # top: T_max relative to steady Q90 at the same time since heat-on
        rel = smax[m] - np.interp(t[m], ts, smax_s)
        axs[0].plot(t[m] - t_pre, rel, label=lab, **style)
        axs[1].plot(t[m] - t_pre, 1e3 * dtc[m], **style)
    # change from steady Q90 at the same time (caption); "delta" is Table I's equal-energy symbol
    # top label split so only "(K)" sits beside the axis, clear of the bottom label
    axs[0].set_ylabel(r"$T_{\max}$ change" + "\n(K)", linespacing=1.0)
    axs[1].set_ylabel(r"$\Delta T_{\mathrm{cell}}$ (mK)")
    axs[1].set_xlabel("Time since heat-on (s)")
    axs[1].set_xlim(0, 900)
    h, lbl = axs[0].get_legend_handles_labels()
    # one legend row above the panels, clear of the curves
    fig.legend(h, lbl, frameon=False, loc="upper center", ncol=3, fontsize=7,
               handlelength=2.0, columnspacing=1.2, bbox_to_anchor=(0.55, 1.0),
               borderaxespad=0.1)
    fig.tight_layout(pad=0.3, h_pad=0.4, rect=(0, 0, 1, 0.88))
    fig.align_ylabels(axs)
    fig.savefig(path)
    plt.close(fig)


def make_all(rows, paths, figdir):
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig_geometry(figdir / "fig_geometry.pdf")
    fig_pareto(rows, figdir / "fig_pareto.pdf")
    fig_histories(rows, paths, figdir / "fig_histories.pdf")
    print(f"figures -> {figdir}")
