#!/usr/bin/env python3
"""Measure a K0 mesh from its polyMesh (WP 2.5a; plan D-I acceptance numbers).

Reads the undivided blockMesh mesh (constant/polyMesh + cellZones: the
fluid / cell<k> zones that splitMeshRegions -cellZones turns into regions)
and reports, from the mesh itself (not from the generator's design):

  * per-zone cell counts;
  * first fluid cell height on every fluid-solid interface face
    (2 x distance from the face to the fluid cell centre along the face
    normal), split into lateral surface and top cap;
  * wall-normal stacks walked cell by cell from every interface face:
    cells inside the film (D-I band 0.37-0.43 mm, counted at 0.37 mm) and
    inside the Stokes layers sqrt(nu/(pi f)) at 2 and 4 Hz (nu = mu/rho from
    properties.yaml, or --nu for the D-R film temperature), the D-I gate
    thickness (thinner of that and the quoted 0.18 mm at 4 Hz), the wall
    distance where the 5th / 3rd cell ends (the margins), largest successive
    thickness ratio in the first 1 mm, and, for faces whose walk lands on
    another cell's surface within 2.5 mm, the number of cells across that
    gap (the 2 mm gaps);
  * over all internal fluid faces: expansion ratio d_N/d_P (face-normal
    distances to the two cell centres) and non-orthogonality, with the
    location of the worst face;
  * meshed volume per zone (divergence theorem, exact for the planar faces
    of this extruded mesh) and each cell zone against pi r^2 h: the heat
    source is normalised by pi r^2 h until D-J (i), so this ratio is the
    injected-power factor of the level;
  * optionally the checkMesh -allRegions summary parsed from its log.

--gate applies the D-I / D-C mesh thresholds (k0_targets.py) and exits 1 on
any failure: film >= 5 cells at 0.37 mm, Stokes >= 3 cells at the highest
frequency run on the mesh (--max-freq, default 2 Hz = CORE; 4 Hz for a
PLUS-4 mesh), near-wall ratio <= 2.0, checkMesh OK in all 7 regions with
non-orthogonality <= 65 and skewness <= 4, counts equal to a committed K0
level with r_N >= 1.3 per region to its neighbours, tile expansion ratio
< 1.3 and the all-fluid expansion-ratio regression guard.

The configuration is detected from the mesh (WP 2.5b): cells below the
floor mean the D-G housing, whose port length L follows from the inlet
position. For the housing the gate adds: patches exactly inlet / outlet /
walls; the inlet is the full 46 x 10 mm channel entry at x = -90 - L and
the outlet the full 10 x 46 mm riser top at z = 72 + L; the meshed volume
equals housing + channel + riser; expansion ratio <= 2.0 in the port
regions (channel, slot columns, riser); the z-lever is active (port-column
cells at the cap plane >= 4 x the module's); and the highest CORE flow
does not set the Courant limit in the outflow column: the developed-laminar
peak (1.7355 x the mean slot velocity) at Q 500 mL/min gives Co <= 1 at the
level's Q 500 time step (dt* 6.25 ms; dt*/2 on xfine, k0_targets). Counts
are matched against the housing family at the detected port length.

The Courant hot layer of P00b (K0-design §10.14): every fluid face in the
cap plane z = H inside the module region (fluid ring and strips) with its
thinner cell's height and its horizontal distance d to the nearest cell
wall, reported per d band and for the design band d >= 1 mm (whether the
P00b maxima lie in it is not known from its summary: k0_co_lever.py
samples), with the Courant number of the P00b plume (w = 17 mm/s) at
dt 12.5 / 6.25 / 3.125 ms and the largest dt for Co <= 1
(`capplane_plume`). --zlever far|module (K0_ZLEVER, from Allrun.pre) adds
a gate item that the mesh carries the z-lever the case asked for: the
module-interior lever (module), or the uniform cap layer (far, the
committed design).

Cell centres are vertex averages (exact for parallelepipeds; the thin
wall cells are near-parallelepipeds), so non-orthogonality here can differ
from checkMesh by a fraction of a degree; checkMesh remains the authority
for the quality gate.

D-F.9 (1) (PROJECT.md §2.7, 2026-10-02): --thin-head --level LEVEL writes
the headspace part T_head of the Courant rule's thin layer from the split
case's own fluid region mesh (housing_physics.THIN_HEAD: layer-B fluid
cells within pitch/sqrt 2 of a cell axis whose smallest horizontal extent,
V / largest side-face area, is below 0.2 mm / r_h^k) as the cellSet
constant/fluid/polyMesh/sets/thinHeadSet plus thin_head.json; Allrun.pre
turns it into the fluid cellZone thinHead (topoSet) that the CoThin mask of
housing_physics.monitors reads. It is a mandatory Allrun.pre step (not
skipped by K0_MEASURE=0) and refuses an empty set.

Usage: k0_measure.py CASE [--json out.json] [--checkmesh log.checkMesh]
                       [--nu NU] [--gate [--max-freq HZ] [--level LEVEL]
                                         [--zlever far|module]]
       k0_measure.py CASE --thin-head --level LEVEL
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import math

import numpy as np

import k0_targets as T

FILM = T.FILM_MM * 1e-3   # D-I film band lower bound [m]
NEAR = 1.0e-3        # expansion ratios reported within this distance [m]
GAP_MAX = 2.5e-3     # a walk that meets another solid within this = a gap


# ---------------------------------------------------------------------------
# ascii polyMesh reader
# ---------------------------------------------------------------------------

def _body(path):
    txt = Path(path).read_text()
    i = txt.find("FoamFile")
    i = txt.find("}", i) + 1
    txt = re.sub(r"//.*", "", txt[i:])
    m = re.search(r"(\d+)\s*\(", txt)
    n = int(m.group(1))
    return n, txt[m.end():txt.rfind(")")]


def read_points(path):
    n, body = _body(path)
    a = np.array(body.replace("(", " ").replace(")", " ").split(), float)
    return a.reshape(n, 3)


def read_faces(path):
    n, body = _body(path)
    a = np.array(body.replace("(", " ").replace(")", " ").split(), np.int64)
    a = a.reshape(n, 5)
    if not np.all(a[:, 0] == 4):
        raise SystemExit("non-quad faces: this reader handles hex meshes only")
    return a[:, 1:]


def read_labels(path):
    n, body = _body(path)
    return np.array(body.split(), np.int64)[:n]


def read_boundary(path):
    _, body = _body(path)
    out = {}
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", body):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        out[m.group(1)] = (int(d["startFace"]), int(d["nFaces"]))
    return out


def read_cellzones(path, ncells):
    txt = Path(path).read_text()
    i = txt.find("}", txt.find("FoamFile")) + 1
    body = txt[i:]
    zone = np.full(ncells, -1, np.int64)
    names = []
    for m in re.finditer(r"(\w+)\s*\{\s*type\s+cellZone;\s*cellLabels\s+"
                         r"List<label>\s*(\d+)\s*\(([^)]*)\)", body):
        k = len(names)
        names.append(m.group(1))
        lab = np.array(m.group(3).split(), np.int64)
        assert lab.size == int(m.group(2))
        zone[lab] = k
    return names, zone


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------

def face_geometry(P, F):
    p = P[F]                                     # (nF, 4, 3)
    ctr = p.mean(axis=1)
    area = 0.5 * np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
    mag = np.linalg.norm(area, axis=1)
    return ctr, area / mag[:, None], mag


def cell_centres(P, F, own, nei, ncells):
    fp = P[F].sum(axis=1)                        # sum of the 4 face points
    s = np.zeros((ncells, 3))
    c = np.zeros(ncells)
    np.add.at(s, own, fp)
    np.add.at(c, own, 4)
    ni = nei.size
    np.add.at(s, nei, fp[:ni])
    np.add.at(c, nei, 4)
    return s / c[:, None]


def cell_face_table(own, nei, ncells):
    """(ncells, 6) face labels per cell (hex)."""
    nf = own.size
    cells = np.concatenate([own, nei])
    faces = np.concatenate([np.arange(nf), np.arange(nei.size)])
    order = np.lexsort((faces, cells))
    cells, faces = cells[order], faces[order]
    cnt = np.bincount(cells, minlength=ncells)
    if not np.all(cnt == 6):
        raise SystemExit("cells with != 6 faces: hex meshes only")
    return faces.reshape(ncells, 6)


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------

def walk_stacks(start_face, start_cell, normal, F, own, nei, cf, Cf, zone,
                fluid_id, max_steps=80):
    """March from each start face through the fluid, always leaving a hex
    through the face opposite to the entry face (no shared vertex).
    Returns per stack: thickness list (along the start normal) and the zone
    that ended the walk (-2 boundary, k solid zone, -1 still fluid)."""
    n = start_face.size
    ni = nei.size
    entry = start_face.copy()
    cell = start_cell.copy()
    pos = Cf[start_face]
    alive = np.ones(n, bool)
    thick = np.full((n, max_steps), np.nan)
    end_zone = np.full(n, -1, np.int64)
    for step in range(max_steps):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        c = cell[idx]
        cand = cf[c]                                    # (m, 6)
        shared = (F[cand][:, :, :, None] ==
                  F[entry[idx]][:, None, None, :]).any(axis=(2, 3))
        shared |= cand == entry[idx][:, None]
        # exactly one face shares no vertex with the entry face
        opp = cand[np.arange(idx.size), np.argmin(shared, axis=1)]
        d = np.einsum("ij,ij->i", Cf[opp] - pos[idx], normal[idx])
        thick[idx, step] = np.abs(d)
        pos[idx] = Cf[opp]
        entry[idx] = opp
        internal = opp < ni
        nxt = np.where(internal,
                       np.where(own[np.minimum(opp, own.size - 1)] == c,
                                nei[np.minimum(opp, ni - 1)], own[opp]), -1)
        stop_b = ~internal
        end_zone[idx[stop_b]] = -2
        alive[idx[stop_b]] = False
        ok = internal
        zn = np.where(ok, zone[np.maximum(nxt, 0)], -2)
        stop_s = ok & (zn != fluid_id)
        end_zone[idx[stop_s]] = zn[stop_s]
        alive[idx[stop_s]] = False
        go = ok & (zn == fluid_id)
        cell[idx[go]] = nxt[go]
    return thick, end_zone


def stokes_layers(nu, freqs=T.REPORT_FREQS_HZ):
    """{label: thickness [m]}: computed layers at the frequencies and the
    D-I gate thickness at each (k0_targets)."""
    out = {}
    for f in freqs:
        out[f"stokes_{f:g}Hz"] = T.stokes_thickness(f, nu)
        out[f"stokes_{f:g}Hz_gate"] = T.stokes_gate_thickness(f, nu)
    return out


def stack_stats(thick, end_zone, stokes):
    t = np.nan_to_num(thick, nan=0.0)
    cum = np.cumsum(t, axis=1)
    valid = ~np.isnan(thick)
    in_film = ((cum <= FILM * (1 + 1e-9)) & valid).sum(axis=1)
    in_stokes = {k: ((cum <= d * (1 + 1e-9)) & valid).sum(axis=1)
                 for k, d in stokes.items()}
    # wall distance where the N-th cell ends (the margin to a threshold)
    end_film_n = np.where(valid[:, T.FILM_MIN_CELLS - 1],
                          cum[:, T.FILM_MIN_CELLS - 1], np.nan)
    end_stokes_n = np.where(valid[:, T.STOKES_MIN_CELLS - 1],
                            cum[:, T.STOKES_MIN_CELLS - 1], np.nan)
    near = (cum - t < NEAR) & valid
    r = np.full(thick.shape, np.nan)
    r[:, 1:] = np.maximum(thick[:, 1:] / thick[:, :-1],
                          thick[:, :-1] / thick[:, 1:])
    r_near = np.where(near[:, :], r, np.nan)
    qmax = np.nanmax(np.where(np.isnan(r_near), 1.0, r_near), axis=1)
    total = np.where(valid, t, 0).sum(axis=1)
    ncell = valid.sum(axis=1)
    gap = (end_zone >= 0) & (total <= GAP_MAX)
    return dict(h1=thick[:, 0], in_film=in_film, in_stokes=in_stokes,
                end_film_n=end_film_n, end_stokes_n=end_stokes_n,
                q_near=qmax,
                gap=gap, gap_cells=np.where(gap, ncell, 0),
                gap_width=np.where(gap, total, np.nan))


def summarise(x):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if x.size == 0:
        return None
    return dict(n=int(x.size), min=float(x.min()), mean=float(x.mean()),
                max=float(x.max()))


def parse_checkmesh(log):
    """Per-region checkMesh -allRegions summary from its log."""
    txt = Path(log).read_text()
    out = {}
    # checkMesh -allRegions prints one "Mesh stats <region>" block per region
    # (a single-region run prints "Mesh stats" without a name)
    for blk in re.split(r"\nMesh stats", txt)[1:]:
        first, _, body = blk.partition("\n")
        region = first.strip() or "region0"
        blk = body
        g = lambda pat: (re.search(pat, blk) or [None, None])[1]
        out[region] = dict(
            cells=int(g(r"cells:\s+(\d+)") or -1),
            patches=int(g(r"boundary patches:\s+(\d+)") or -1),
            mesh_ok=("Mesh OK." in blk) and ("Failed" not in blk),
            max_non_ortho=float(g(r"non-orthogonality Max:\s*([\d.eE+-]+)")
                                or "nan"),
            max_skewness=float(g(r"Max skewness = ([\d.eE+-]+)") or "nan"),
            max_aspect=float(g(r"Max aspect ratio = ([\d.eE+-]+)") or "nan"),
        )
    return out


def zone_volumes(F, Cf, Nf, Af, own, nei, zone, names, ncells):
    """Per-zone volume, divergence theorem V = 1/3 sum(Cf . Sf) with
    outward face normals (the owner->neighbour normal is outward for the
    owner). Exact for planar faces, which all faces of this 2.5-D
    extruded hex mesh are."""
    cs = np.einsum("ij,ij->i", Cf, Nf) * Af / 3.0
    v = np.zeros(ncells)
    np.add.at(v, own, cs)
    np.add.at(v, nei, -cs[:nei.size])
    return {n: float(v[zone == k].sum()) for k, n in enumerate(names)}, v


def housing_geometry(P):
    """None for the open duct; for the D-G housing (cells below the floor)
    its geometry [m] from k0_blockmesh.GEOM and the port length from the
    inlet position."""
    if P[:, 2].min() > -1e-9:
        return None
    import k0_blockmesh as K
    g = K.load_geometry()
    L = 1e3 * (g["x_min"] * 1e-3 - float(P[:, 0].min()))       # [mm]
    pg = K.port_geometry(g, 1.0)
    W = 2 * g["y_half"]
    vol = (g["x_max"] - g["x_min"]) * W * g["z_top"] \
        + (L + g["x_port"]) * W * g["h_ch"] + g["x_port"] * W * L
    gm = {k: g[k] * 1e-3 for k in ("x_min", "x_max", "x_port", "y_half",
                                   "z_top", "h_ch", "H")}
    return dict(g=gm, L_mm=L, port_mult=L / pg["L1"], L1_mm=pg["L1"],
                Dh_mm=pg["Dh_slot"], volume_m3=vol * 1e-9)


def patch_report(bnd, Cf, Nf, Af):
    """Per patch: faces, area, centre-coordinate ranges, mean unit normal."""
    out = {}
    for name, (start, n) in bnd.items():
        sl = slice(start, start + n)
        a = Af[sl]
        nrm = (Nf[sl] * a[:, None]).sum(axis=0) / a.sum()
        out[name] = dict(
            faces=int(n), area_m2=float(a.sum()),
            centre_min=Cf[sl].min(axis=0).tolist(),
            centre_max=Cf[sl].max(axis=0).tolist(),
            mean_normal=nrm.round(12).tolist())
    return out


def capplane_dz(Cf, Nf, dP, dN, ff, zcap, xm):
    """Heights of the fluid cells on both sides of the cap plane z = zcap
    (the internal fluid faces lying in it), min per x-region: module
    (|x| <= module half-length), bridge (the U2/D2 H-blocks) and ports
    (|x| >= module half-length + 50 mm, the slot columns). [mm]"""
    f = Cf[ff]
    on = (np.abs(Nf[ff][:, 2]) > 0.99) & (np.abs(f[:, 2] - zcap) < 1e-9)
    h = 2 * np.minimum(dP, dN)
    ax = np.abs(f[:, 0])
    out = {}
    for k, m in (("module", ax <= xm + 1e-9),
                 ("bridge", (ax > xm + 1e-9) & (ax < xm + 0.05 - 1e-9)),
                 ("ports", ax >= xm + 0.05 - 1e-9)):
        mm = on & m
        if mm.any():
            out[k] = float(1e3 * h[mm].min())
    return out


def zone_axes(C, zone, names):
    """(x, y) axis of every cell zone [m]: the mean of its cell centres
    (exact for the symmetric O-grid of a K0 cylinder)."""
    return np.array([C[zone == k][:, :2].mean(axis=0)
                     for k, n in enumerate(names) if n != "fluid"])


def capplane_plume(Cf, Nf, dP, dN, ff, zcap, xm, axes, R, h_cap_mm):
    """The Courant hot layer of P00b (§3) under the module z-lever
    (K0-design §10.14): every internal fluid face in the cap plane z = zcap
    inside the module region |x| <= xm (tiles and strips; over a cap the
    plane is the fluid-solid interface, not a fluid face). Per face: the
    height of the thinner of its two cells (2 x the smaller centre
    distance) and the horizontal distance d from the face centre to the
    nearest cell wall (|xy - axis| - R). Reported per d band
    (k0_targets.PLUME_D_BINS_MM) and for the plume band d >= PLUME_D_MIN_MM
    (a design band, not established as the P00b maxima's; see
    k0_co_lever.py) and d >= 0.5 mm: the smallest height, where, its
    ratio to the cap cell, and the P00b plume's Courant number at the
    PLUME_DT_S branches with the largest dt for Co <= 1 (w = PLUME_W_M_S).
    [m in, mm out]"""
    f = Cf[ff]
    on = ((np.abs(Nf[ff][:, 2]) > 0.99) & (np.abs(f[:, 2] - zcap) < 1e-9)
          & (np.abs(f[:, 0]) <= xm + 1e-9))
    h = 1e3 * 2 * np.minimum(dP, dN)[on]
    xy = f[on][:, :2]
    d = 1e3 * (np.min(np.linalg.norm(xy[:, None, :] - axes[None, :, :],
                                     axis=2), axis=1) - R)
    out = dict(w_mm_s=1e3 * T.PLUME_W_M_S, module_cap_mm=h_cap_mm,
               faces=int(on.sum()), max_dz_mm=float(h.max()),
               max_ratio_to_cap=float(h.max() / h_cap_mm),
               d_min_mm=float(d.min()), bins=[])

    def pick(msk):
        if not msk.any():
            return dict(faces=0, min_dz_mm=None, ratio_to_cap=0.0)
        j = int(np.argmin(np.where(msk, h, np.inf)))
        return dict(faces=int(msk.sum()), min_dz_mm=float(h[j]),
                    ratio_to_cap=float(h[j] / h_cap_mm),
                    at_mm=(1e3 * xy[j]).round(3).tolist(),
                    d_at_mm=round(float(d[j]), 4),
                    **{k: v for k, v in T.plume_courant(float(h[j])).items()
                       if k != "dz_mm"})
    edges = list(T.PLUME_D_BINS_MM) + [np.inf]
    for lo, hi in zip(edges[:-1], edges[1:]):
        msk = (d >= lo) & (d < hi)
        if msk.any():
            out["bins"].append(dict(d_mm=[lo, None if hi == np.inf else hi],
                                    **pick(msk)))
    out["band"] = dict(d_min_mm=T.PLUME_D_MIN_MM,
                       **pick(d >= T.PLUME_D_MIN_MM))
    out["band_0p5"] = dict(d_min_mm=0.5, **pick(d >= 0.5))
    return out


def outflow_column(P, F, C, cf, hs):
    """Smallest z-extent of the cells of the outflow column (D1 slot column
    and riser) and the Courant number there of the highest CORE flow
    (k0_targets): of the mean slot velocity at dt*, and the largest dt at
    which the developed-laminar peak (CO_OUTFLOW_PEAK_FACTOR x mean) stays
    at Co 1. The gate scales co_mean_outflow to the peak at the level's
    Q 500 time step."""
    g = hs["g"]
    col = np.nonzero(C[:, 0] > g["x_max"] - g["x_port"])[0]
    z = P[F[cf[col]]][..., 2].reshape(col.size, -1)
    dz = z.max(axis=1) - z.min(axis=1)
    area = g["x_port"] * 2 * g["y_half"]
    w = T.CO_OUTFLOW_Q_M3S / area
    i = int(np.argmin(dz))
    return dict(cells=int(col.size), dz_min_mm=float(1e3 * dz[i]),
                dz_min_at_mm=(1e3 * C[col[i]]).round(3).tolist(),
                w_mean_mm_s=1e3 * w,
                co_mean_outflow=float(w * T.CO_OUTFLOW_DT_S / dz[i]),
                peak_factor=T.CO_OUTFLOW_PEAK_FACTOR,
                dt_max_peak_s=float(T.CO_OUTFLOW_MAX * dz[i]
                                    / (T.CO_OUTFLOW_PEAK_FACTOR * w)),
                q_m3s=T.CO_OUTFLOW_Q_M3S, dt_s=T.CO_OUTFLOW_DT_S)


def measure(case, nu=None, freqs=T.REPORT_FREQS_HZ):
    nu = T.coolant_nu() if nu is None else nu
    freqs = sorted({float(f) for f in freqs})
    stokes = stokes_layers(nu, freqs)
    pm = Path(case) / "constant" / "polyMesh"
    P = read_points(pm / "points")
    F = read_faces(pm / "faces")
    own = read_labels(pm / "owner")
    nei = read_labels(pm / "neighbour")
    ncells = int(own.max()) + 1
    names, zone = read_cellzones(pm / "cellZones", ncells)
    bnd = read_boundary(pm / "boundary")
    fluid_id = names.index("fluid")
    Cf, Nf, Af = face_geometry(P, F)
    C = cell_centres(P, F, own, nei, ncells)
    cf = cell_face_table(own, nei, ncells)
    ni = nei.size
    res = dict(case=str(case), cells={n: int((zone == k).sum())
                                      for k, n in enumerate(names)},
               cells_total=int(ncells), patches={k: v[1] for k, v in
                                                 bnd.items()})

    # interface faces (fluid | solid), oriented into the fluid
    zo, zn = zone[own[:ni]], zone[nei]
    iface = np.nonzero((zo != zn) & ((zo == fluid_id) | (zn == fluid_id)))[0]
    fl_is_own = zo[iface] == fluid_id
    fcell = np.where(fl_is_own, own[iface], nei[iface])
    n_in = Nf[iface] * np.where(fl_is_own, -1.0, 1.0)[:, None]  # into fluid
    # (owner->neighbour normal points out of the owner; into the fluid when
    #  the fluid is the neighbour)
    d = np.einsum("ij,ij->i", C[fcell] - Cf[iface], n_in)
    if np.any(d <= 0):
        raise SystemExit("interface normal orientation error")
    cap = np.abs(n_in[:, 2]) > 0.9
    thick, endz = walk_stacks(iface, fcell, n_in, F, own, nei, cf, Cf, zone,
                              fluid_id)
    st = stack_stats(thick, endz, stokes)
    h1_centre = 2 * d
    res["interface_faces"] = dict(lateral=int((~cap).sum()),
                                  cap=int(cap.sum()))
    res["h1_centre_mm"] = {k: summarise(1e3 * h1_centre[m]) for k, m in
                           (("lateral", ~cap), ("cap", cap))}
    res["h1_walk_mm"] = {k: summarise(1e3 * st["h1"][m]) for k, m in
                         (("lateral", ~cap), ("cap", cap))}
    parts = (("lateral", ~cap), ("cap", cap))
    res["nearwall_targets"] = dict(
        nu_m2s=nu, film_mm=1e3 * FILM, film_band_mm=list(T.FILM_BAND_MM),
        stokes_mm={k: 1e3 * v for k, v in stokes.items()})
    res["cells_in_stokes_layer"] = {
        s: {k: summarise(st["in_stokes"][s][m]) for k, m in parts}
        for s in stokes}
    res["cells_in_film"] = {k: summarise(st["in_film"][m]) for k, m in parts}
    # margins: the largest wall distance (over faces) at which the 5th cell
    # (film) / 3rd cell (Stokes) ends; a threshold below it loses a cell
    res["film_cell5_end_mm"] = {k: summarise(1e3 * st["end_film_n"][m])
                                for k, m in parts}
    res["stokes_cell3_end_mm"] = {k: summarise(1e3 * st["end_stokes_n"][m])
                                  for k, m in parts}
    e3 = float(np.nanmax(st["end_stokes_n"]))
    res["stokes_nu_min_m2s"] = {
        f"{f:g}Hz": T.nu_min_for(f, e3) for f in freqs}
    res["max_ratio_first_1mm"] = {k: summarise(st["q_near"][m])
                                  for k, m in parts}
    g = st["gap"] & ~cap
    res["gap_stacks"] = dict(
        n_faces=int(g.sum()),
        cells_across=summarise(st["gap_cells"][g]),
        width_mm=summarise(1e3 * st["gap_width"][g]))
    # per gap: which two cells each gap stack joins (7 internal gaps in the
    # 2 x 3 module); cells across and narrowest width per pair
    solid_start = np.where(fl_is_own, zn[iface], zo[iface])
    res["gaps"] = {}
    for a_, b_ in sorted({tuple(sorted((int(s), int(e)))) for s, e in
                          zip(solid_start[g], endz[g])}):
        m = g & (((solid_start == a_) & (endz == b_))
                 | ((solid_start == b_) & (endz == a_)))
        res["gaps"][f"{names[a_]}-{names[b_]}"] = dict(
            faces=int(m.sum()),
            cells_across_min=int(st["gap_cells"][m].min()),
            cells_across_max=int(st["gap_cells"][m].max()),
            min_width_mm=float(1e3 * np.nanmin(st["gap_width"][m])))
    # expansion ratio and non-orthogonality on internal fluid-fluid faces
    ff = np.nonzero((zo == fluid_id) & (zn == fluid_id))[0]
    dP = np.abs(np.einsum("ij,ij->i", C[own[ff]] - Cf[ff], Nf[ff]))
    dN = np.abs(np.einsum("ij,ij->i", C[nei[ff]] - Cf[ff], Nf[ff]))
    er = np.maximum(dP, dN) / np.minimum(dP, dN)
    dv = C[nei[ff]] - C[own[ff]]
    cosang = np.einsum("ij,ij->i", dv, Nf[ff]) / np.linalg.norm(dv, axis=1)
    no = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
    ie, io = int(np.argmax(er)), int(np.argmax(no))
    res["fluid_expansion_ratio"] = dict(
        max=float(er.max()), p999=float(np.percentile(er, 99.9)),
        p99=float(np.percentile(er, 99)),
        frac_gt_1p2=float((er > 1.2 + 1e-6).mean()),
        worst_at_mm=(1e3 * Cf[ff[ie]]).round(3).tolist())
    res["fluid_non_ortho_deg"] = dict(
        max=float(no.max()), mean=float(no.mean()),
        p99=float(np.percentile(no, 99)),
        worst_at_mm=(1e3 * Cf[ff[io]]).round(3).tolist())
    # where do the large ratios live? module tiles / side strips / up-down
    # stream H-blocks, and the headspace (z > cell top) separately
    xf, yf, zf = (1e3 * Cf[ff]).T
    xm = 1e3 * float(np.max(np.abs(Cf[iface][:, 0]))) + 1.0 + 1e-6
    ym = 1e3 * float(np.max(np.abs(Cf[iface][:, 1]))) + 1.0 + 1e-6
    ztop = 1e3 * float(Cf[iface][cap][:, 2].max()) if cap.any() else 1e9
    head = zf > ztop + 1e-6
    hs = housing_geometry(P)
    res["config"] = "housing" if hs else "duct"
    if hs:
        res["housing"] = {k: v for k, v in hs.items() if k != "g"}
        g = hs["g"]
        zt = 1e3 * g["z_top"]
        chan = zf < -1e-6
        riser = zf > zt + 1e-6
        portcol = ~chan & ~riser & (np.abs(xf) >= 1e3 * (
            g["x_max"] - g["x_port"]) - 1e-6)
        box = ~chan & ~riser & ~portcol
    else:
        box = np.ones(xf.shape, bool)
    regions = {
        "tiles": (np.abs(xf) <= xm) & (np.abs(yf) <= ym) & ~head & box,
        "strips": (np.abs(xf) <= xm) & (np.abs(yf) > ym) & ~head & box,
        "updown": (np.abs(xf) > xm) & ~head & box,
        "headspace_over_module": head & (np.abs(xf) <= xm) & box,
        "headspace_updown": head & (np.abs(xf) > xm) & box,
    }
    if hs:
        regions.update(port_columns=portcol, channel=chan, riser=riser,
                       ports=portcol | chan | riser)
    res["fluid_expansion_by_region"] = {}
    for k, msk in regions.items():
        if not msk.any():
            continue
        e = er[msk]
        j = int(np.argmax(e))
        res["fluid_expansion_by_region"][k] = dict(
            faces=int(msk.sum()), max=float(e.max()),
            p99=float(np.percentile(e, 99)),
            frac_gt_1p2=float((e > 1.2 + 1e-6).mean()),
            frac_gt_1p5=float((e > 1.5).mean()),
            frac_gt_2=float((e > 2.0).mean()),
            worst_at_mm=(1e3 * Cf[ff[msk]][j]).round(3).tolist(),
            non_ortho_max=float(no[msk].max()))
    # meshed volumes: each cell zone against pi r^2 h (properties.yaml)
    vol, _ = zone_volumes(F, Cf, Nf, Af, own, nei, zone, names, ncells)
    import yaml
    cellp = yaml.safe_load(T.PROPS.read_text())["cell"]
    v_cyl = math.pi * (float(cellp["diameter"]) / 2) ** 2 \
        * float(cellp["height"])
    lo, hi = P.min(axis=0), P.max(axis=0)
    res["volume_m3"] = dict(
        zones=vol, total=float(sum(vol.values())),
        bounding_box=float(np.prod(hi - lo)),
        pi_r2_h=v_cyl,
        cell_over_pi_r2_h={n: v / v_cyl for n, v in vol.items()
                           if n != "fluid"})
    # cap-plane (z = cell top) cell heights by x-region: the Courant hot
    # layer of 2.5a (K0-design §7 item 5) and the housing z-lever
    res["capplane_dz_mm"] = capplane_dz(Cf, Nf, dP, dN, ff, ztop * 1e-3,
                                        xm * 1e-3)
    # the plume-crossing cells beside the cells (module z-lever, §10.14)
    if cap.any():
        res["capplane_plume"] = capplane_plume(
            Cf, Nf, dP, dN, ff, ztop * 1e-3, xm * 1e-3,
            zone_axes(C, zone, names), 0.5 * float(cellp["diameter"]),
            res["h1_walk_mm"]["cap"]["min"])
    if hs:
        res["housing"]["patches"] = patch_report(bnd, Cf, Nf, Af)
        res["volume_m3"]["housing_analytic"] = hs["volume_m3"]
        res["outflow_column"] = outflow_column(P, F, C, cf, hs)
    return res, dict(Cf=Cf, iface=iface, cap=cap, st=st, er=er, ff=ff,
                     no=no)


# ---------------------------------------------------------------------------
# D-F.9 (1): the headspace thin set T_head of the Courant rule (housing)
# ---------------------------------------------------------------------------

def cell_extents(P, F, own, nei, ncells):
    """Per hex cell of a polyMesh: volume V (divergence theorem, exact for
    planar faces), centre (vertex average, as cell_centres), smallest
    horizontal extent h_hor = V / (largest side-face area; side = |n_z| <
    0.5) and vertical extent h_z = V / (largest top/bottom face area). On
    the 2.5-D K0 mesh every face is vertical or horizontal, so h_hor is the
    cell's width across its narrowest horizontal direction."""
    Cf, Nf, Af = face_geometry(P, F)
    xs = np.einsum("ij,ij->i", Cf, Nf) * Af / 3.0
    ni = nei.size
    V = np.zeros(ncells)
    np.add.at(V, own, xs)
    np.add.at(V, nei, -xs[:ni])
    side = np.abs(Nf[:, 2]) < 0.5
    a_side = np.zeros(ncells)
    a_top = np.zeros(ncells)
    for lab, n in ((own, own.size), (nei, ni)):
        a, sd = Af[:n], side[:n]
        np.maximum.at(a_side, lab[sd], a[sd])
        np.maximum.at(a_top, lab[~sd], a[~sd])
    with np.errstate(divide="ignore", invalid="ignore"):
        return dict(V=V, C=cell_centres(P, F, own, nei, ncells),
                    h_hor=V / a_side, h_z=V / a_top)


def _housing_physics():
    model = str(Path(__file__).resolve().parents[1] / "model")
    if model not in sys.path:
        sys.path.insert(0, model)
    import housing_physics as HP
    return HP


def thin_head_mask(C, h_hor, level, fluid=None, axes=None):
    """D-F.9 (1) T_head (housing_physics.THIN_HEAD): the fluid cells whose
    centre lies in layer B (65 <= z <= 72 mm) within the tile half-diagonal
    (pitch / sqrt 2 = 14.14 mm) of a cell axis and whose smallest horizontal
    extent is below the level's threshold (0.2 mm on medium, / r_h^k;
    housing_physics.thin_head_h). C [m], h_hor [m]."""
    HP = _housing_physics()
    th = HP.THIN_HEAD
    axes = np.asarray(HP.cell_axes() if axes is None else axes, float)
    C = np.asarray(C, float)
    r = np.sqrt((C[:, 0, None] - axes[:, 0]) ** 2
                + (C[:, 1, None] - axes[:, 1]) ** 2).min(axis=1)
    m = ((C[:, 2] >= th["z"][0]) & (C[:, 2] <= th["z"][1]) & (r <= th["r"])
         & (np.asarray(h_hor, float) < HP.thin_head_h(level)))
    return m if fluid is None else m & np.asarray(fluid, bool)


def thin_head(polymesh, level):
    """T_head on a polyMesh directory (the fluid region mesh of a split case,
    whose cells are all fluid, or the undivided mesh, whose 'fluid' cellZone
    selects them): (sorted cell labels, summary)."""
    HP = _housing_physics()
    pm = Path(polymesh)
    P = read_points(pm / "points")
    F = read_faces(pm / "faces")
    own = read_labels(pm / "owner")
    nei = read_labels(pm / "neighbour")
    nc = int(own.max()) + 1
    fluid = None
    if (pm / "cellZones").exists():
        names, zone = read_cellzones(pm / "cellZones", nc)
        if "fluid" in names:
            fluid = zone == names.index("fluid")
    ext = cell_extents(P, F, own, nei, nc)
    T = thin_head_mask(ext["C"], ext["h_hor"], level, fluid)
    lab = np.flatnonzero(T)
    C_mm = 1e3 * ext["C"]
    axes = 1e3 * np.asarray(HP.cell_axes(), float)
    r = np.sqrt((C_mm[:, 0, None] - axes[:, 0]) ** 2
                + (C_mm[:, 1, None] - axes[:, 1]) ** 2).min(axis=1)
    tl = HP.THIN_LAYER
    zb = (C_mm[:, 2] >= 1e3 * tl["sleeve_z"][0]) \
        & (C_mm[:, 2] <= 1e3 * tl["sleeve_z"][1])
    sleeve = zb & (r >= 1e3 * tl["sleeve_r"][0]) & (r <= 1e3 * tl["sleeve_r"][1])
    if fluid is not None:
        sleeve &= fluid
    summ = dict(
        rule="D-F.9 (1)", mesh=str(pm), level=level, cells_mesh=nc,
        h_hor_below_mm=1e3 * HP.thin_head_h(level),
        z_mm=[1e3 * v for v in HP.THIN_HEAD["z"]],
        r_max_mm=1e3 * HP.THIN_HEAD["r"], n_cells=int(lab.size),
        r_mm=[float(r[T].min()), float(r[T].max())] if lab.size else None,
        h_hor_mm=([float(1e3 * ext["h_hor"][T].min()),
                   float(1e3 * ext["h_hor"][T].max())] if lab.size else None),
        n_beyond_sleeve_r=int((r[T] > 1e3 * tl["sleeve_r"][1]).sum()),
        sleeve_cells=int(sleeve.sum()),
        sleeve_cells_in_T_head=int((sleeve & T).sum()))
    return lab, summ


def write_cellset(path, labels, name, location):
    """An ascii OpenFOAM cellSet (labelHashSet list) file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    head = ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
            "    arch        \"LSB;label=32;scalar=64\";\n"
            "    class       cellSet;\n"
            f"    location    \"{location}\";\n    object      {name};\n}}\n"
            "// written by scripts/mesh/k0_measure.py --thin-head (D-F.9 (1))\n\n")
    body = f"{len(labels)}\n(\n" + "".join(f"{int(i)}\n" for i in labels) + ")\n"
    path.write_text(head + body)


def write_thin_head(case, level, region="fluid"):
    """k0_measure.py CASE --thin-head --level LEVEL: T_head on the case's own
    region mesh, written as the cellSet constant/<region>/polyMesh/sets/
    <THIN_HEAD set> (Allrun.pre's topoSet makes it the cellZone the CoThin
    mask reads) and summarised in CASE/thin_head.json. An empty set is a
    hard failure: the monitors would silently fall back to the band."""
    HP = _housing_physics()
    case = Path(case)
    pm = case / "constant" / region / "polyMesh"
    lab, summ = thin_head(pm, level)
    if not lab.size:
        raise SystemExit(f"k0_measure: T_head is empty on {pm} at level "
                         f"{level} (D-F.9 (1)): refusing to write the set")
    name = HP.THIN_HEAD["set"]
    write_cellset(pm / "sets" / name, lab, name,
                  f"constant/{region}/polyMesh/sets")
    summ.update(cellset=str((pm / "sets" / name).relative_to(case)),
                zone=HP.THIN_HEAD["zone"])
    (case / "thin_head.json").write_text(json.dumps(summ, indent=1) + "\n")
    return summ


def family_level(cells, config="duct", port_mult=1.0):
    """(level, {level: per-zone counts}) of the committed K0 family whose
    exact block-product counts equal the measured ones (level None if no
    level matches, e.g. a --set exploration variant). The housing family is
    taken at the given port multiple (None if it has no committed counts)."""
    import k0_blockmesh as K
    if config == "housing":
        try:
            topo = K.Topology(K.load_geometry(), K.DESIGN_HOUSING, "housing",
                              port_mult)
        except ValueError:
            return None, {}
    else:
        topo = K.Topology(K.load_geometry(), K.DESIGN)
    fam = {lv: K.render(topo, lv)[1] for lv in K.LEVELS}
    hit = [lv for lv, c in fam.items() if c == cells]
    return (hit[0] if hit else None), fam


def gate(res, max_freq=T.CORE_MAX_FREQ_HZ, level=None, zlever=None):
    """D-I / D-C mesh gate. Returns [(item, value, limit, ok)]. zlever
    ('far' or 'module', K0_ZLEVER): the z-lever the case asked for; the
    gate then also checks that the mesh carries it (zlever_gate)."""
    items = []

    def add(name, value, limit, ok):
        items.append((name, value, limit, bool(ok)))

    def mn(d):
        return min(v["min"] for v in d.values() if v)

    def mx(d):
        return max(v["max"] for v in d.values() if v)

    add(f"film cells within {T.FILM_MM} mm (min, lateral+cap)",
        mn(res["cells_in_film"]), f">= {T.FILM_MIN_CELLS}",
        mn(res["cells_in_film"]) >= T.FILM_MIN_CELLS)
    nu = res["nearwall_targets"]["nu_m2s"]
    d = T.stokes_gate_thickness(max_freq, nu)
    key = f"stokes_{max_freq:g}Hz_gate"
    if key not in res["cells_in_stokes_layer"]:
        raise SystemExit(f"gate at {max_freq:g} Hz: measure() was not run "
                         f"with that frequency")
    n3 = mn(res["cells_in_stokes_layer"][key])
    add(f"Stokes cells within {1e3 * d:.4f} mm ({max_freq:g} Hz; min)",
        n3, f">= {T.STOKES_MIN_CELLS}", n3 >= T.STOKES_MIN_CELLS)
    q = mx(res["max_ratio_first_1mm"])
    add("near-wall ratio, first 1 mm (max)", round(q, 4),
        f"<= {T.NEARWALL_RATIO_MAX}", q <= T.NEARWALL_RATIO_MAX)
    cm = res.get("checkMesh")
    zones = sorted(res["cells"])
    if not cm:
        add("checkMesh -allRegions", "not given", "--checkmesh LOG", False)
    else:
        add("checkMesh regions", sorted(cm), zones, sorted(cm) == zones)
        add("checkMesh Mesh OK (all regions)",
            [r for r, v in cm.items() if not v["mesh_ok"]] or "all",
            "all", all(v["mesh_ok"] for v in cm.values()))
        no = max(v["max_non_ortho"] for v in cm.values())
        add("max non-orthogonality [deg]", no, f"<= {T.NON_ORTHO_MAX}",
            no <= T.NON_ORTHO_MAX)
        sk = max(v["max_skewness"] for v in cm.values())
        add("max skewness", sk, f"<= {T.SKEWNESS_MAX}", sk <= T.SKEWNESS_MAX)
        add("checkMesh counts == measured",
            {r: v["cells"] for r, v in cm.items()}, "equal",
            all(cm.get(z, {}).get("cells") == n
                for z, n in res["cells"].items()))
    hs = res.get("housing")
    pm = round(hs["port_mult"], 6) if hs else 1.0
    found, fam = family_level(res["cells"], res.get("config", "duct"), pm)
    if level is not None and found != level:
        add("K0 family level", found, level, False)
    elif found is None:
        add("K0 family level", "no committed level matches the counts",
            "coarse|medium|fine|xfine", False)
    else:
        lv = list(fam)
        i = lv.index(found)
        rn = {}
        for a_, b_ in ((i - 1, i), (i, i + 1)):
            if a_ < 0 or b_ >= len(lv):
                continue
            for z in fam[found]:
                rn[f"{lv[a_]}->{lv[b_]} {z}"] = (
                    fam[lv[b_]][z] / fam[lv[a_]][z]) ** (1 / 3)
        worst = min(rn, key=rn.get)
        add(f"r_N per region to neighbours of {found} (min: {worst})",
            round(rn[worst], 4), f">= {T.RN_MIN}", rn[worst] >= T.RN_MIN)
    er = res["fluid_expansion_by_region"]
    add("expansion ratio, module tiles (max)", round(er["tiles"]["max"], 4),
        f"< {T.ER_MAX_TILES}", er["tiles"]["max"] < T.ER_MAX_TILES)
    e = res["fluid_expansion_ratio"]["max"]
    guard = T.ER_GUARD_ALL_FLUID_HOUSING if hs else T.ER_GUARD_ALL_FLUID
    add("expansion ratio, all fluid (max; regression guard"
        + (", 2.33 at medium accepted 2026-09-25 instead of the 2.0 "
           "target, K0-design §10.6)" if hs else ")"),
        round(e, 4), f"<= {guard}", e <= guard)
    if hs:
        items += housing_gate(res, level if level is not None else found)
    if zlever is not None:
        items += zlever_gate(res, zlever)
    return items


def zlever_gate(res, zlever):
    """The module-interior z-lever is what the case asked for (K0-design
    §10.14). 'module': the module region's largest cap-plane cell (the
    tile-corner lines) >= CAPPLANE_LEVER_MIN x the cap cell, and the plume
    band's (d >= PLUME_D_MIN_MM) smallest >= ZLEVER_MODULE_BAND_MIN x the
    cap cell. 'far' (the committed design): the module region carries the
    uniform cap layer (largest cap-plane cell = the cap cell, 2 %)."""
    if zlever not in T.ZLEVERS:
        raise SystemExit(f"--zlever {zlever!r}: one of {T.ZLEVERS}")
    items = []
    cp = res.get("capplane_plume")
    if cp is None:
        return [("z-lever of the module region", "no cap-plane faces",
                 zlever, False)]
    mx = cp["max_ratio_to_cap"]
    if zlever == "module":
        items.append(("module z-lever: largest module cap-plane cell / cap "
                      "cell", round(mx, 3), f">= {T.CAPPLANE_LEVER_MIN}",
                      mx >= T.CAPPLANE_LEVER_MIN))
        b = cp["band"]["ratio_to_cap"]
        items.append((f"module z-lever: plume band (d >= "
                      f"{T.PLUME_D_MIN_MM:g} mm) smallest cap-plane cell / "
                      f"cap cell", round(b, 3),
                      f">= {T.ZLEVER_MODULE_BAND_MIN}",
                      b >= T.ZLEVER_MODULE_BAND_MIN))
    else:
        items.append(("module z-lever off: largest module cap-plane cell / "
                      "cap cell (uniform cap layer)", round(mx, 3),
                      "<= 1.02", mx <= 1.02))
    return items


def housing_gate(res, level=None):
    """D-G housing items (WP 2.5b). Returns [(item, value, limit, ok)].
    level (the expected or detected K0 level) sets the Q 500 time step of
    the outflow-column Courant item (k0_targets.co_outflow_dt; dt* if
    unknown)."""
    items = []

    def add(name, value, limit, ok):
        items.append((name, value, limit, bool(ok)))

    hs = res["housing"]
    import k0_blockmesh as K
    g = K.load_geometry()
    L = hs["L_mm"]
    W = 2 * g["y_half"]
    pt = hs["patches"]
    add("housing patches", sorted(pt), ["inlet", "outlet", "walls"],
        sorted(pt) == ["inlet", "outlet", "walls"])
    for name, area, axis, pos, nrm in (
            ("inlet", W * g["h_ch"], 0, g["x_min"] - L, (-1, 0, 0)),
            ("outlet", g["x_port"] * W, 2, g["z_top"] + L, (0, 0, 1))):
        p = pt.get(name)
        ok = p is not None and math.isclose(
            1e6 * p["area_m2"], area, rel_tol=1e-6) and all(
            abs(1e3 * p[k][axis] - pos) < 1e-6
            for k in ("centre_min", "centre_max")) and all(
            abs(a - b) < 1e-9 for a, b in zip(p["mean_normal"], nrm))
        add(f"{name}: full {area:g} mm2 port at "
            f"{'xyz'[axis]} = {pos:.6g} mm, normal {nrm}",
            None if p is None else dict(
                area_mm2=round(1e6 * p["area_m2"], 6),
                at_mm=[round(1e3 * p[k][axis], 6)
                       for k in ("centre_min", "centre_max")],
                normal=p["mean_normal"]), "exact", ok)
    v = res["volume_m3"]
    add("meshed volume == housing + channel + riser",
        dict(meshed=v["total"], analytic=v["housing_analytic"]),
        "rel 1e-9", math.isclose(v["total"], v["housing_analytic"],
                                 rel_tol=1e-9))
    er = res["fluid_expansion_by_region"]["ports"]["max"]
    add("expansion ratio, ports (channel, slot columns, riser; max)",
        round(er, 4), f"<= {T.ER_MAX_PORTS}", er <= T.ER_MAX_PORTS)
    cz = res["capplane_dz_mm"]
    lev = cz["ports"] / cz["module"]
    add("z-lever: port-column / module cap-plane cell height",
        round(lev, 3), f">= {T.CAPPLANE_LEVER_MIN}",
        lev >= T.CAPPLANE_LEVER_MIN)
    dt = T.co_outflow_dt(level)
    co = res["outflow_column"]["co_mean_outflow"] * T.CO_OUTFLOW_PEAK_FACTOR \
        * dt / T.CO_OUTFLOW_DT_S
    add(f"outflow column peak Co ({T.CO_OUTFLOW_PEAK_FACTOR:.4g} x mean, Q "
        f"{T.CO_OUTFLOW_Q_M3S * 6e7:g} mL/min, dt {dt:g} s at {level})",
        round(co, 4), f"<= {T.CO_OUTFLOW_MAX}", co <= T.CO_OUTFLOW_MAX)
    return items


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--checkmesh", type=Path,
                    help="checkMesh -allRegions log to parse")
    ap.add_argument("--nu", type=float,
                    help="coolant kinematic viscosity [m2/s] for the Stokes "
                         "layers (default mu/rho from properties.yaml; D-R: "
                         "pass the film-temperature value)")
    ap.add_argument("--gate", action="store_true",
                    help="apply the D-I/D-C thresholds; exit 1 on failure")
    ap.add_argument("--max-freq", type=float, default=T.CORE_MAX_FREQ_HZ,
                    help="highest pulsation frequency run on this mesh [Hz] "
                         "(default 2 = CORE; 4 for a PLUS-4 mesh)")
    ap.add_argument("--level", help="expected K0 level (gate fails if the "
                                    "counts are another level's)")
    ap.add_argument("--zlever", choices=T.ZLEVERS,
                    help="z-lever the case asked for (K0_ZLEVER; the gate "
                         "checks that the mesh carries it: 'module' = the "
                         "module-interior lever of K0-design §10.14, 'far' "
                         "= the committed uniform cap layer in the module)")
    ap.add_argument("--thin-head", action="store_true",
                    help="D-F.9 (1) only: write the headspace thin set T_head "
                         "of the split case's fluid region mesh as the cellSet "
                         "constant/fluid/polyMesh/sets/thinHeadSet and "
                         "thin_head.json (needs --level: the threshold "
                         "0.2 mm / r_h^k); no other measurement")
    a = ap.parse_args(argv)
    if a.thin_head:
        if not a.level:
            ap.error("--thin-head needs --level (the grid level of the "
                     "threshold)")
        print(json.dumps(write_thin_head(a.case, a.level), indent=1))
        return 0
    res, _ = measure(a.case, a.nu,
                     freqs=set(T.REPORT_FREQS_HZ) | {a.max_freq})
    if a.checkmesh:
        res["checkMesh"] = parse_checkmesh(a.checkmesh)
    rc = 0
    if a.gate:
        items = gate(res, a.max_freq, a.level, a.zlever)
        res["gate"] = dict(max_freq_hz=a.max_freq, passed=all(
            i[3] for i in items), items=[dict(item=n, value=v, limit=l,
                                               ok=ok) for n, v, l, ok in items])
        rc = 0 if res["gate"]["passed"] else 1
    txt = json.dumps(res, indent=1, default=str)
    if a.json:
        a.json.write_text(txt)
    print(txt)
    if a.gate:
        print("\nD-I / D-C mesh gate (k0_targets.py):", file=sys.stderr)
        for n, v, l, ok in items:
            print(f"  {'PASS' if ok else 'FAIL'}  {n}: {v}  [{l}]",
                  file=sys.stderr)
        print(f"GATE {'PASSED' if rc == 0 else 'FAILED'}", file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
