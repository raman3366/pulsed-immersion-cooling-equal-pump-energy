#!/usr/bin/env python3
"""AISEE 2027 conference figures (D-AIS.1), called by aisee_results.py.

fig_waveforms.pdf  imposed flow waveforms over two seconds
fig_pareto.pdf     pump energy against T_max,peak (a) and dT_cell,end (b):
                   the steady curve and the pulsed arms (matched mean flow:
                   filled; reduced mean flow aimed at steady Q90's energy:
                   open)
fig_histories.pdf  T_max(t) and dT_cell(t) over the discharge, steady Q90
                   against the pulsed arms

Okabe-Ito hues (validated, scripts in the dataviz skill; adjacent CVD pair
in the 6-8 band), so every series also has its own marker and line style
(grayscale print) and the values are in the paper's table.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLOR = {"steady": "#0072B2", "C1": "#E69F00", "C2p": "#009E73", "C3": "#CC79A7"}
MARKER = {"steady": "o", "C1": "s", "C2p": "D", "C3": "^"}
LS = {"steady": "-", "C1": "--", "C2p": "-.", "C3": ":"}
LABEL = {"steady": "steady", "C1": "square 1 Hz, D = 0.25",
         "C2p": "square 0.5 Hz, D = 0.5", "C3": "sine 2 Hz, A = 0.6"}
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e3e3e3"
NOISE_TMAX_K = 0.004      # GDT1 dt and seed pairs, volume T_max,peak (aisee_results.py)
COL_W = 3.5               # IEEE column width, inches

plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.2, "pdf.fonttype": 42, "font.family": "serif",
})


def q_wave(arm, qbar, t, tau_v=0.05):
    """Imposed volumetric flow (mL/min): the generator's waveforms after the
    ramp; square = mean-preserving trapezoid, edges tau_v centred on the
    switching instants (housing_physics.square_table)."""
    if arm == "steady":
        return np.full_like(t, qbar)
    if arm == "C3":
        return qbar * (1 + 0.6 * np.sin(2 * np.pi * 2.0 * t))
    f, d = (1.0, 0.25) if arm == "C1" else (0.5, 0.5)
    p = 1 / f
    q_on = qbar / d
    ph = np.mod(t, p)
    rise = np.clip((ph + tau_v / 2) / tau_v, 0, 1)
    fall = np.clip((d * p + tau_v / 2 - ph) / tau_v, 0, 1)
    wrap = np.clip((ph - (p - tau_v / 2)) / tau_v, 0, 1)
    return q_on * np.maximum(np.minimum(rise, fall), wrap)


def fig_waveforms(path):
    t = np.linspace(0, 2, 4001)
    fig, ax = plt.subplots(figsize=(COL_W, 1.95))
    for arm, q in (("C1", 90), ("C2p", 90), ("C3", 90), ("steady", 90)):
        ax.plot(t, q_wave(arm, q, t), color=COLOR[arm], ls=LS[arm], lw=1.1,
                label=LABEL[arm])
    ax.set_xlabel("time in period window $t$ (s)")
    ax.set_ylabel(r"$Q$ (mL min$^{-1}$)")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 380)
    ax.set_title(r"imposed inlet flow at $\bar Q$ = 90 mL min$^{-1}$", color=INK)
    h, lbl = ax.get_legend_handles_labels()
    fig.legend(h, lbl, loc="lower center", ncol=2, frameon=False, handlelength=2.4,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(pad=0.3, rect=(0, 0.2, 1, 1))
    fig.savefig(path)
    plt.close(fig)


# label offsets (points) of the pulsed symbols in panel (a), clear of their neighbours
PULSED_LABEL_AT = {"C3_Q82.9": (-17, 3)}


def _main_rows(rows):
    return [r for r in rows if r["role"] in ("curve", "mean-flow", "matched-energy")]


def _tag(r):
    """Point label: Qbar, with a dagger when the case failed a gate (--allow-fail;
    rule 4: flagged, never plotted silently; Table I names the gate)."""
    return f"{r['Qbar_mL_min']:g}" + ("†" if r.get("gates_failed") else "")


def fig_pareto(rows, path):
    rows = _main_rows(rows)
    fig, axs = plt.subplots(1, 2, figsize=(COL_W * 2.05, 2.25))
    for ax, key, lab in ((axs[0], "Tmax_peak_C", r"$T_{\max,\mathrm{peak}}$ ($^\circ$C)"),
                         (axs[1], "dTcell_end_K", r"$\Delta T_{\mathrm{cell,end}}$ (K)")):
        st = sorted([r for r in rows if r["arm"] == "steady"], key=lambda r: r["E_pump_mJ"])
        if st:
            x = [r["E_pump_mJ"] for r in st]
            y = [r[key] for r in st]
            ax.plot(x, y, color=COLOR["steady"], ls="-", marker="o", ms=4.5,
                    mfc=COLOR["steady"], mec="white", mew=0.6, label="steady", zorder=3)
            # the run-to-run noise (NOISE_TMAX_K, 0.4 mK in dT_cell) is smaller than the
            # symbols at this scale: stated in the caption, not drawn
            for r in st:
                ax.annotate(_tag(r), (r["E_pump_mJ"], r[key]),
                            textcoords="offset points", fontsize=6, color=MUTED,
                            xytext=(3, 3) if key == "Tmax_peak_C" or r["Qbar_mL_min"] < 100
                            else (-16, 5))
        for arm in ("C1", "C2p", "C3"):
            for role, filled in (("mean-flow", True), ("matched-energy", False)):
                pts = [r for r in rows if r["arm"] == arm and r["role"] == role]
                if not pts:
                    continue
                ax.plot([r["E_pump_mJ"] for r in pts], [r[key] for r in pts],
                        ls="none", marker=MARKER[arm], ms=5.5,
                        mfc=COLOR[arm] if filled else "white", mec=COLOR[arm], mew=1.2,
                        label=f"{LABEL[arm]}, " + (r"$\bar Q$ = 90" if filled else
                                                   r"reduced $\bar Q$"), zorder=4)
                # pulsed points are labelled in panel (a) only: panel (b) is crowded
                # near 0.3-0.5 mJ, and the legend identifies each point (one per
                # arm and role)
                for r in pts if key == "Tmax_peak_C" else []:
                    ax.annotate(_tag(r), (r["E_pump_mJ"], r[key]),
                                textcoords="offset points", fontsize=6, color=MUTED,
                                xytext=PULSED_LABEL_AT.get(r["key"], (4, -8)))
        ax.set_xscale("log")
        from matplotlib.ticker import LogLocator, FuncFormatter
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
        ax.set_xlabel(r"pump energy over the 900 s discharge $E_{\mathrm{pump}}$ (mJ)")
        ax.set_ylabel(lab)
    axs[0].set_title("(a) peak cell temperature", color=INK, loc="left")
    axs[1].set_title("(b) cell-to-cell spread at end of discharge", color=INK, loc="left")
    h, lbl = axs[0].get_legend_handles_labels()
    fig.legend(h, lbl, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02), handlelength=1.8)
    fig.tight_layout(pad=0.3, rect=(0, 0.13, 1, 1))
    fig.savefig(path)
    plt.close(fig)


def fig_histories(rows, paths, path):
    import aisee_results as AR
    # the mechanism figure: the same mean flow, Q90 (matched-energy arms are in Fig. 2)
    keys = [r["key"] for r in _main_rows(rows)
            if r["key"] in paths and (r["role"] == "mean-flow" or r["key"] == "steady_Q90")]
    if not keys:
        return
    fig, axs = plt.subplots(2, 1, figsize=(COL_W, 3.3), sharex=True)
    for k in keys:
        r = next(x for x in rows if x["key"] == k)
        t, smax, dtc, _ = AR.series(paths[k])
        t_pre = 30.0
        m = t > t_pre
        arm = r["arm"]
        style = dict(color=COLOR[arm], ls=LS[arm] if r["role"] != "matched-energy" else ":",
                     lw=1.0 if r["role"] == "matched-energy" else 1.2)
        lab = f"{LABEL[arm]}, {r['Qbar_mL_min']:g}"
        axs[0].plot(t[m] - t_pre, smax[m] - 273.15, label=lab, **style)
        axs[1].plot(t[m] - t_pre, dtc[m], **style)
    axs[0].set_ylabel(r"$T_{\max}$ ($^\circ$C)")
    axs[1].set_ylabel(r"$\Delta T_{\mathrm{cell}}$ (K)")
    axs[1].set_xlabel("time since heat-on (s)")
    axs[0].legend(frameon=False, loc="lower right", fontsize=6)
    fig.tight_layout(pad=0.3)
    fig.savefig(path)
    plt.close(fig)


def make_all(rows, paths, figdir):
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig_waveforms(figdir / "fig_waveforms.pdf")
    fig_pareto(rows, figdir / "fig_pareto.pdf")
    fig_histories(rows, paths, figdir / "fig_histories.pdf")
    print(f"figures -> {figdir}")
