#!/usr/bin/env python3
"""K0 structured blockMesh generator for the 6-cell 18650 module (WP 2.5a).

Plan: D-I (K0 = scripted structured O/H blockMesh, first mesh candidate),
D-C (4-level geometrically similar family, r_h = 1.35), NEW-1. Design and
numbers: docs/mesh/K0-design.md.

Topology (2.5-D: one 2-D multi-block layout extruded through two z-layers)
--------------------------------------------------------------------------
Each cylinder sits in a square tile of side = pitch (20 mm). Per tile:
  * solid O-grid: 2x2 core blocks + 8 solid-ring blocks up to the cylinder
    surface (zone cell<k> in layer A);
  * fluid O-grid ring: 8 blocks from the cylinder surface to the tile edge,
    block edges at 0/45/90...deg, so every gap (tile mid-side, 1 mm from the
    wall) and every interstitial diagonal (tile corner, 5.14 mm from the
    wall) is a block edge with its own wall-normal grading.
H-blocks fill the side strips (tile edge -> side wall) and the upstream /
downstream sections (module -> inlet / outlet), each split into a 10 mm
"port column" at the duct end (x in [-90,-80] and [80,90]) so that WP 2.5b
can attach the D-G inlet channel below the floor and the riser above the lid
to those columns without touching the rest of the topology.
z-layers: A = [0, 65] mm (cells stand on the floor), B = [65, 72] mm
(headspace). Solid blocks exist only in layer A; the blocks above the cells
are fluid, so the top caps are conformal interfaces too. splitMeshRegions
-cellZones therefore produces fluid + cell0..cell5 with 1:1 conformal
mapped interfaces.

Refinement (geometric similarity)
---------------------------------
Every block-direction count belongs to a count class (TH, R, S, ST, U1, U2,
D1, D2, ZA, ZB). The design is fixed at the MEDIUM level; level k (coarse
-1, medium 0, fine +1, xfine +2) uses n = round(n_medium * 1.35**k). Every
graded edge is a fixed piecewise-geometric mapping x(s) (per section: length
fraction f, cell fraction nu, kappa = n ln q) designed at medium and SAMPLED
at the level's count, so each level is the same mapping on a uniformly
refined computational grid (Roache's systematic refinement): first-cell
heights scale ~1/r, expansion ratios as q**(1/r). Integer rounding of the
counts (and of the per-section counts) is the only departure from exact
similarity; it is reported per class and per region (--report).

Configurations (--config)
-------------------------
  duct     (default; WP 2.5a, cases/k0proto) the open production duct of
           D-10: inlet = the x = -90 face, outlet = the x = +90 face.
  housing  (WP 2.5b, cases/k0housing) the closed module housing of D-G:
           closed end walls; a full-width 10 mm inlet slot in the floor at
           x in [-90, -80], fed from below by a horizontal 46 x 10 mm inlet
           channel of straight length L upstream of the slot (inlet patch at
           x = -90 - L, z in [-10, 0]); a full-width 10 mm outlet slot in
           the lid at x in [80, 90] draining into a vertical 10 x 46 mm
           riser of length L (outlet patch at z = 72 + L). L = port_mult x
           10 D_h (D_h = 16.43 mm, the D-G start length L = 164.3 mm);
           --port-mult 1 / 2 / 4 gives the D-G port-gate variants L/2L/4L.
           Housing-only design changes (DESIGN_HOUSING): port columns
           U1/D1 = 8 (they now carry the slot and riser flow and follow the
           family); U2/D2 graded two-sided down to the port cell (split
           chosen for the smoothest ratios over the whole level family);
           the far-field lines (port columns, channel, riser) get their own
           smooth y-distribution (side strip 5 mm wide instead of 3 mm, the
           U2/D2 blocks bridge the two); the module-end strip lines
           x = +-30 mm shrink toward the wall at exactly Q_MAX
           (strip_Cend: the U2/D2 module-face jump 4.8 -> 2.3); and the
           per-edge z-grading lever
           (zlever): every vertical vertex line with |x| > 30 mm carries a
           "far" z-mapping without the 0.05 mm cap-plane layer (the module
           tiles and strips keep it), so the Courant hot layer and the
           aspect-ratio-70 cells leave the H-blocks, port columns and the
           outflow column. Cell counts per layer are global in a conformal
           multi-block mesh, so the lever is count-neutral.
           OPTION (off by default; --zlever module, env K0_ZLEVER=module in
           the housing case; K0-design §10.14): the module-interior
           z-lever. The cylinder lines keep the 0.05 mm cap-plane cells
           (every cap-face column: the D-I cap film and Stokes layer are
           unchanged); the tile-edge lines (gap throats, tile corners) and
           the strip side-wall lines take coarser cap-plane cells
           (LeverMapping: one fixed map of the module mapping per line
           family, so the z-line tilt across the fluid ring is the same at
           every level), and the fluid ring blocks blend the two. Also
           count-neutral.

Usage
-----
    k0_blockmesh.py LEVEL [-o system/blockMeshDict] [--report report.json]
                    [--manifest constant/triSurface/MANIFEST.sha256]
                    [--nu NU]          # film-temperature nu (D-R); default
                                       # mu/rho from properties.yaml
                    [--config housing [--port-mult 1|2|4]
                                      [--zlever far|module]]
    k0_blockmesh.py --table [--config housing]   # counts, all levels
    k0_blockmesh.py --table --minimize --set dz_max=3.0 --set n.TH=7
                                       # explore variants (never for a
                                       # claim-bearing mesh: use DESIGN)
LEVEL in {coarse, medium, fine, xfine}. --manifest cross-checks GEOM against
the make_geometry.py manifest that makePhysics reads (rows, cols, gap, d, h,
z0): the two geometry sources must agree or the mesh is not written.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import k0_targets as T

REPO = Path(__file__).resolve().parents[2]
PROPS = REPO / "data" / "reference" / "properties.yaml"

R_REFINE = 1.35                       # D-C: r_h per direction
LEVELS = {"coarse": -1, "medium": 0, "fine": 1, "xfine": 2}
Q_MAX = 1.2                           # target max expansion ratio (medium)

# ---- geometry [mm] ---------------------------------------------------------
# Cell D/H are replaced by properties.yaml values (load_geometry); array and
# duct as the production duct (D-02 gap, D-10 duct, identical to
# cases/template/makeMesh.py).
GEOM = dict(
    D=18.0, H=65.0, gap=2.0, rows=2, cols=3,
    x_min=-90.0, x_max=90.0, y_half=23.0, z_top=72.0,
    x_port=10.0,          # end "port columns" for WP 2.5b (inlet slot / riser)
    rho_core_mid=0.72,    # core mid-side vertex radius / R   (at 0, 90, ...)
    rho_core_corner=0.78, # core corner vertex radius / R     (at 45, 135, ...)
)

# ---- design at the MEDIUM level (targets; the counts are the family anchor) -
DESIGN = dict(
    h1_wall=0.05,         # first fluid cell on the lateral cell surface [mm]
    h1_cap=0.05,          # first fluid cell above the top cap (z = 65) [mm]
    h_floor=0.50,         # first cell on the floor (adiabatic) [mm]
    h_lid=0.25,           # first cell below the lid (adiabatic) [mm]
    dz_max=2.5,           # largest cell in z (mid-height plateau) [mm]
    h_side=0.15,          # first cell at the side walls y = +-23 [mm]
    q_far=1.2,            # max ratio in the up/downstream H-blocks
    h_solid_wall=0.06,    # last solid-ring cell at the surface [mm]. Above
                          # the caps (layer B) these blocks are fluid and meet
                          # the ring's first cell (0.035-0.05 mm around the
                          # circumference) across r = R: ratio 1.2-1.7, inside
                          # the D-I transition limit 2.0 (target 1.5)
    n=dict(               # medium counts per class; (*) = the smallest count
                          # that keeps every ratio <= Q_MAX at medium
        TH=8,             # per 45-deg ring block (arc) / per 10 mm tile half-side
        R=17,             # (*) fluid ring, wall -> tile edge (5.14 mm diagonal)
        S=10,             # (*) solid ring, core -> wall (ratio <= 1.31: the
                          # radial direction is wall-parallel above the caps
                          # and inside the solid below)
        ST=15,            # (*) side strip, tile edge -> side wall (3 mm)
        U2=30, D2=30,     # (*) module <-> port column (50 mm each)
        U1=3, D1=3,       # port columns (10 mm each), uniform; the U2/D2
                          # plateau is matched to their cell (ratio <= 1.2)
        ZA=48,            # (*) z in [0, 65]
        ZB=23,            # (*) z in [65, 72]
    ),
    zlever=False,         # per-edge z-grading lever (housing only; see
                          # DESIGN_HOUSING). False = one z-mapping everywhere
)

CONFIGS = ("duct", "housing")

# ---- D-G housing (WP 2.5b) -------------------------------------------------
# Geometry beyond the duct box [mm]. The D-G text fixes: closed housing
# 46 x 72 mm, x in [-90, 90]; full-width 10 mm bottom inlet slot fed by a
# horizontal inlet channel; full-width 10 mm top outlet slot into a vertical
# riser; start length L = 10 D_h ~ 164 mm (doubling by gate). Choices made
# here (K0-design.md §10): the slots sit in the existing 10 mm port columns
# at the housing ends (x in [-90,-80] / [80,90]); the channel is 10 mm high
# (46 x 10 mm, the same D_h = 16.43 mm as the slot/riser, which is what makes
# the D-G "10 D_h ~ 164 mm" consistent) and runs below the floor with its
# top at z = 0; L is the STRAIGHT length upstream of the slot column (inlet
# patch at x = -90 - L) and the riser length above the lid (outlet patch at
# z = 72 + L); the slot column itself (the elbow under the inlet slot) is
# not counted in L.
GEOM_HOUSING = dict(
    h_ch=10.0,            # inlet channel height (z in [-h_ch, 0])
    port_Dh_mult=10.0,    # D-G: start L = 10 D_h
    y_far_strip=5.0,      # side strip width on the far-field vertex lines
                          # (port columns, channel, riser): 5 mm makes the
                          # far y-distribution smooth (max ratio 1.11 at
                          # medium) instead of the module's 3 mm strip, whose
                          # 0.26 mm tile-end cell meets 1.25 mm cells (4.8)
)

# Medium design of the housing: the duct DESIGN plus the housing classes.
# XL (channel length) and ZR (riser length) depend on the port length; their
# committed medium counts per port multiple are in n_port.
DESIGN_HOUSING = json.loads(json.dumps(DESIGN))
DESIGN_HOUSING["n"].update(
    U1=8, D1=8,           # port columns carry the slot / riser flow now and
                          # must follow the family (6/8/11/15, ratios
                          # 1.33/1.38/1.36; K0-design §2.3); 1.25 mm cells
    U2=33, D2=33,         # (*) module -> port column, two-sided down to the
                          # port cell, plateau <= h_far_max
    ZC=13,                # (*) inlet channel height, 10 mm
)
DESIGN_HOUSING.update(
    zlever=True,          # far z-mapping on |x| > 30 mm vertex lines
    h_far_max=3.75,       # U2/D2 plateau cap [mm] (3 x the port cell; the
                          # 2.5a duct plateau was 3.7 mm)
    h_run=4.0,            # streamwise plateau in the channel and riser [mm]
                          # (developing laminar duct flow; aspect ratio 3.2
                          # against the 1.25 mm cross-stream cells)
    h_ch_wall=0.5,        # first cell at the channel floor and ceiling [mm]
                          # (= h_floor, the housing floor cell it meets
                          # across the slot)
    h_ch_max=1.25,        # channel cross-stream (z) plateau [mm]
    n_port={              # (*) medium XL / ZR per port multiple: the
                          # smallest counts whose split keeps every level's
                          # ratio (with the neighbouring slot-column / lid
                          # cell) within FAM_TOL of the family
        "1": dict(XL=45, ZR=51),
        "2": dict(XL=86, ZR=92),
        "4": dict(XL=168, ZR=175),
    },
    # module-interior z-lever (K0-design §10.14): an OPTION, off by default
    # (K0_ZLEVER=module / --zlever module switches it on; the GDT1 Delta-t
    # ladder decides whether it is adopted). The cylinder lines (axis,
    # core, cylinder surface: every cap face column) keep the module
    # mapping with its 0.05 mm cap-plane cells (D-I cap film and Stokes
    # layer unchanged); the tile-edge lines (gap throats G, tile corners C)
    # and the side-wall lines of the strips take coarser cap-plane cells,
    # so the plumes that rise in the fluid ring cross taller cells. The
    # ring blocks blend the two lines linearly in the radial length
    # fraction (blockCreate.C), so the z-lines tilt across the ring by
    # D / L_ring (D the offset of the two lines' points, L_ring 1 mm at
    # the gap throat, 5.14 mm at the corner). Counts are unchanged.
    zlever_module=False,
    h_cap_gap=0.17,       # cap-plane cell of the gap-throat (G) lines [mm]
                          # at medium: the largest (0.01 mm steps) whose
                          # offset D keeps the gap-ray tilt atan(D / 1 mm)
                          # <= 30 deg (k0_targets.ZLEVER_TILT_MAX_DEG;
                          # 0.17: 29.6 deg, 0.18: 32.9)
    h_cap_corner=0.0,     # cap-plane cell of the corner (C) lines [mm];
                          # 0 = the far layer-B cell 7 mm / ZB (0.304 mm,
                          # uniform above the cap plane; tilt 19.6 deg over
                          # the 5.14 mm diagonal). It sets the plume band
                          # beside the cells and the added non-orthogonality
                          # in the G-C ring blocks near the tile sides
                          # (in-plane skew ~40 deg there). Coarse sweep at
                          # h_cap_gap 0.17 (K0-design §10.14): 0.20 / 0.25 /
                          # 0.304 / 0.40 mm -> plume band (d >= 1 mm) cell
                          # 1.74 / 1.96 / 2.18 / 2.59 x the cap cell,
                          # checkMesh max non-orthogonality 44.2 / 45.8 /
                          # 47.5 / 52.7 deg (lever off: 43.25)
)


def port_geometry(g, port_mult):
    """Port hydraulic diameters and lengths [mm] (D-G: L = 10 D_h)."""
    W = 2 * g["y_half"]
    dh_slot = 2 * W * g["x_port"] / (W + g["x_port"])      # riser 10 x 46
    dh_ch = 2 * W * g["h_ch"] / (W + g["h_ch"])            # channel 46 x h
    L1 = g["port_Dh_mult"] * max(dh_slot, dh_ch)
    return dict(Dh_slot=dh_slot, Dh_channel=dh_ch, L1=L1,
                L=float(port_mult) * L1, port_mult=float(port_mult))


def load_geometry():
    """GEOM with the cell D and H from properties.yaml (single source of
    truth; a missing or unreadable file is a hard failure, never a silent
    fallback to the nominal values)."""
    import yaml
    g = dict(GEOM, **GEOM_HOUSING)
    cell = yaml.safe_load(PROPS.read_text())["cell"]
    g["D"] = float(cell["diameter"]) * 1e3
    g["H"] = float(cell["height"]) * 1e3
    return g


def check_manifest(g, manifest):
    """Hard failure unless the make_geometry.py manifest (the parameter line
    makePhysics parses for the cell centres, d, h) describes the geometry
    this generator meshes. K0 meshes the true cell extent, so z0 must be 0
    (the template's snappy manifest, z0 = -1 mm / h = 66 mm, is rejected)."""
    line = Path(manifest).read_text().splitlines()[0]
    params = dict(kv.split("=", 1) for kv in line.lstrip("# ").split())
    want = dict(rows=g["rows"], cols=g["cols"], gap=g["gap"] * 1e-3,
                d=g["D"] * 1e-3, h=g["H"] * 1e-3, z0=0.0)
    bad = []
    for k, v in want.items():
        if k not in params:
            bad.append(f"{k} missing")
        elif not math.isclose(float(params[k]), float(v), rel_tol=1e-9,
                              abs_tol=1e-12):
            bad.append(f"{k}: manifest {params[k]} vs K0 {v:g}")
    if bad:
        raise SystemExit(f"{manifest}: geometry disagrees with "
                         f"k0_blockmesh.GEOM/properties.yaml: "
                         + "; ".join(bad))
    return params


# ============================================================================
# 1-D edge mappings (piecewise geometric), designed at medium, sampled at N
# ============================================================================

def _geo_len(h0, q, n):
    """Length of n cells starting at h0 with ratio q."""
    if n <= 0:
        return 0.0
    if abs(q - 1.0) < 1e-12:
        return h0 * n
    return h0 * (q ** n - 1.0) / (q - 1.0)


def _solve_q(L, n, h0):
    """Ratio q such that n cells from first cell h0 span L (q may be < 1)."""
    if n <= 0:
        raise ValueError("n must be >= 1")
    if abs(L - n * h0) < 1e-14 * L:
        return 1.0
    lo, hi = (1.0, 10.0) if L > n * h0 else (1e-6, 1.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _geo_len(h0, mid, n) > L:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


class Mapping:
    """Piecewise-geometric point distribution on [0, 1].

    sections: list of (f, nu, kappa): length fraction, cell-count fraction
    (at the design level) and kappa = n_design * ln(q) (signed; 0 = uniform).
    realize(N) samples the SAME mapping with N cells: n_i = round(nu_i N)
    (total fixed on the largest section, as blockMesh does), q_i =
    exp(kappa_i / n_i), blockMesh expansion E_i = q_i**(n_i - 1).
    """

    def __init__(self, sections, name=""):
        s = [(float(f), float(nu), float(k)) for f, nu, k in sections]
        ft, nt = sum(x[0] for x in s), sum(x[1] for x in s)
        self.sections = [(f / ft, nu / nt, k) for f, nu, k in s]
        self.name = name

    def inverted(self):
        return Mapping([(f, nu, -k) for f, nu, k in reversed(self.sections)],
                       self.name + "^-1")

    def counts(self, N):
        nus = [nu for _, nu, _ in self.sections]
        n = [max(1, int(nu * N + 0.5)) for nu in nus]
        n[max(range(len(n)), key=lambda i: nus[i])] += N - sum(n)
        if min(n) < 1:
            raise ValueError(f"{self.name}: section count < 1 at N={N}: {n}")
        return n

    def realize(self, N):
        """[(f, n_i, E_i)] in blockMesh multi-grading form."""
        out = []
        for (f, _, k), n in zip(self.sections, self.counts(N)):
            q = math.exp(k / n)
            out.append((f, n, q ** (n - 1)))
        return out

    def cells(self, N, L):
        """Cell sizes (length L) exactly as blockMesh's lineDivide makes them."""
        sizes = []
        for f, n, E in self.realize(N):
            q = E ** (1.0 / (n - 1)) if n > 1 else 1.0
            h0 = f * L * ((q - 1.0) / (q ** n - 1.0) if abs(q - 1) > 1e-12
                          else 1.0 / n)
            sizes += [h0 * q ** i for i in range(n)]
        return sizes

    def grading(self, N):
        """blockMesh grading token for this edge (in the edge direction)."""
        r = self.realize(N)
        if len(r) == 1:
            return f"{r[0][2]:.10g}"
        return "(" + " ".join(f"({f:.10g} {n / N:.10g} {E:.10g})"
                              for f, n, E in r) + ")"


UNIFORM = Mapping([(1, 1, 0)], "uniform")


class LeverMapping:
    """z-mapping of a module-interior vertex line under the module z-lever
    (K0_ZLEVER=module; K0-design §10.14).

    `base` is the module z-mapping of the same layer (length L, N_design
    cells at medium). The cap plane z = H is the START of the edge
    (layer B, cap -> lid) or, with flip=True, its END (layer A, floor ->
    cap). The line's distance from the cap plane is a fixed profile T of
    the base's continuous computational coordinate sigma (0 at the cap
    plane, N_design at the far end), designed at medium; at every level
    the line's points are T(sigma(z)) of the base's realised points z,
    with sigma(z) the inverse of the continuous (piecewise-exponential)
    base mapping. Profiles (h = the cap-plane cell at medium):

      track       T = h sigma                      sigma <= sigma*
                  T = T(sigma*) + s (x(sigma) - x(sigma*))   beyond:
                  uniform h cells until the base's cells catch up, then the
                  base's own cells scaled by s <= 1; sigma* and s make
                  T(N) = L and T C1 (h = s x'(sigma*)).
      geometric   T = h (q^sigma - 1) / (q - 1): one geometric run of
                  cells h, h q, h q^2 ... from the cap plane to the far end
                  (q from T(N) = L; q = 1 is uniform).

    So the offset D(z) = T(sigma(z)) - z between this line and the
    cylinder (module) lines is one function of physical position at every
    level: the tilt of the z-lines across the fluid ring is the same at
    every level, and the lever refines like the rest of the family (§3).
    Per-level section rounding cannot split the two lines apart (a
    track_scaled design with its own sections did: 2.8 mm offsets at
    coarse against 0.6 at medium). blockMesh gets the points exactly: one
    multi-grading section (f 1 1) per cell (lineDivide.C: each section
    rounds to one division, placed at the cumulative fraction)."""

    PROFILES = ("track", "geometric")

    def __init__(self, base, L, N_design, h, name, flip=False,
                 profile="track"):
        if profile not in self.PROFILES:
            raise ValueError(f"{name}: profile {profile!r}")
        self.base, self.L, self.N, self.h = base, float(L), int(N_design), h
        self.name, self.flip, self.profile = name, flip, profile
        self.sections = None           # not a piecewise-geometric mapping
        self._cache = {}
        secs = []                      # (f, n, q) of the design realisation
        for f, nn, E in base.realize(self.N):
            secs.append((f, nn, E ** (1.0 / (nn - 1)) if nn > 1 else 1.0))
        self._secs = secs
        if h < self.dxc(0.0) * (1 - 1e-9):
            raise ValueError(f"{name}: h {h} is below the base's cap cell")
        if profile == "track":
            self._design_track()
        else:
            self._design_geometric()

    def _design_track(self):
        h, L = self.h, self.L
        # sigma*: the first root of F = h - s(sig) x'(sig) (C1 junction)
        F = lambda sg: h - (L - h * sg) / (L - self.xc(sg)) * self.dxc(sg)
        lo, hi, step = 0.0, None, 0.05
        sg = step
        while sg < self.N - step:
            if F(sg) <= 0:
                lo, hi = sg - step, sg
                break
            sg += step
        if hi is None:
            raise ValueError(f"{self.name}: no C1 junction for h {h}")
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if F(mid) > 0:
                lo = mid
            else:
                hi = mid
        self.sigma_star = 0.5 * (lo + hi)
        self.z1 = self.xc(self.sigma_star)
        self.p1 = h * self.sigma_star
        self.s = (L - self.p1) / (L - self.z1)
        if not 0 < self.s <= 1 + 1e-12 or self.p1 >= L:
            raise ValueError(f"{self.name}: scale {self.s} (h {h} too large)")

    def _design_geometric(self):
        h, L, N = self.h, self.L, self.N
        if abs(N * h - L) < 1e-12 * L:
            self.q = 1.0
            return
        self.q = _solve_q(L, N, h)

    def T(self, sg):
        """Distance from the cap plane at base index sg (design level)."""
        if self.profile == "geometric":
            if abs(self.q - 1) < 1e-12:
                return self.h * sg
            return self.h * (self.q ** sg - 1) / (self.q - 1)
        if sg <= self.sigma_star:
            return self.h * sg
        return self.p1 + self.s * (self.xc(sg) - self.z1)

    # continuous base mapping, oriented from the cap plane
    def _x(self, sg):
        """Base position (edge direction) at computational index sg."""
        x, a = 0.0, 0
        last = len(self._secs) - 1
        for i, (f, nn, q) in enumerate(self._secs):
            if sg <= a + nn or i == last:
                t = min(max(sg - a, 0.0), nn)
                if abs(q - 1) < 1e-12:
                    return x + f * self.L * t / nn
                return x + f * self.L * (q ** t - 1) / (q ** nn - 1)
            x += f * self.L
            a += nn
        return self.L

    def _dx(self, sg):
        a = 0
        last = len(self._secs) - 1
        for i, (f, nn, q) in enumerate(self._secs):
            if sg <= a + nn or i == last:
                t = min(max(sg - a, 0.0), nn)
                if abs(q - 1) < 1e-12:
                    return f * self.L / nn
                return f * self.L * math.log(q) * q ** t / (q ** nn - 1)
            a += nn
        raise AssertionError

    def xc(self, sg):
        """Base distance from the cap plane at index sg counted from it."""
        return self.L - self._x(self.N - sg) if self.flip else self._x(sg)

    def dxc(self, sg):
        return self._dx(self.N - sg) if self.flip else self._dx(sg)

    def sigma(self, z):
        lo, hi = 0.0, float(self.N)
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if self.xc(mid) < z:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def phi(self, z):
        """Line position (from the cap plane) for base position z."""
        return self.T(self.sigma(z))

    def cells(self, N, L):
        if abs(L - self.L) > 1e-9 * self.L:
            raise ValueError(f"{self.name}: designed for L {self.L}, not {L}")
        if N not in self._cache:
            c = self.base.cells(N, self.L)
            if self.flip:
                c = c[::-1]
            pts = [0.0]
            for x in c:
                pts.append(pts[-1] + x)
            ph = [0.0] + [self.phi(p) for p in pts[1:-1]] + [self.L]
            out = [b - a for a, b in zip(ph[:-1], ph[1:])]
            if min(out) <= 0:
                raise ValueError(f"{self.name}: non-positive cell at N {N}")
            self._cache[N] = out[::-1] if self.flip else out
        return list(self._cache[N])

    def grading(self, N):
        """One (fraction 1 1) section per cell, in the edge direction."""
        return "(" + " ".join(f"({c / self.L:.12g} 1 1)"
                              for c in self.cells(N, self.L)) + ")"


def one_sided(L, N, h_start, name):
    """Single geometric section, first cell h_start (at the edge start)."""
    q = _solve_q(L, N, h_start)
    return Mapping([(1, 1, N * math.log(q))], name)


def fixed_ratio(N, E, name):
    """Single section with blockMesh expansion E (last/first) at N cells."""
    q = E ** (1.0 / (N - 1))
    return Mapping([(1, 1, N * math.log(q))], name)


def max_ratio(cells):
    return max([max(b / a, a / b) for a, b in zip(cells[:-1], cells[1:])]
               or [1.0])


def _kappa_c1(h0, hu, n):
    """kappa of an n-cell geometric section whose discrete first cell is h0
    and whose continuous end slope (per cell) equals hu, i.e. the mapping
    joins a plateau of hu-cells (or a mirrored section) with a continuous
    derivative. Solves n*u + ln(u/(e^u - 1)) = ln(hu/h0), kappa = n*u."""
    c = math.log(hu / h0)
    if c <= 1e-12:
        return 0.0
    u = max(c / n, 1e-9)
    for _ in range(60):
        eu = math.exp(u)
        F = n * u + math.log(u / (eu - 1.0)) - c
        dF = n + 1.0 / u - eu / (eu - 1.0)
        du = F / dF
        u = max(u - du, 1e-12)
        if abs(du) < 1e-15:
            break
    return n * u


def _sec_len(h0, kappa, n):
    """Length of the n-cell section (first cell h0, total log-growth kappa)."""
    if abs(kappa) < 1e-12:
        return h0 * n
    q = math.exp(kappa / n)
    return h0 * (math.exp(kappa) - 1.0) / (q - 1.0)


MIN_SEC = 4    # smallest section at the design level: a 1-3 cell section
               # rounds badly at the other levels (similarity breaks)


def _piecewise(L, N, h_a, h_b, name, q_max, h_max, grow_b, family=False,
               lead=None, fam_tol=None):
    """Search integer splits (n_a grow from h_a, n_u plateau, n_b grow from
    h_b toward the plateau if grow_b) with C1 junctions of the continuous
    mapping, so that sampling the same mapping at another level converges
    to a smooth distribution (junction ratio -> 1 under refinement).
    Returns the smoothest candidate (smallest max ratio) with plateau cell
    <= h_max; for the one-ended case the smallest plateau with ratio <=
    q_max (fastest growth allowed). family=True: among the candidates
    with ratio <= q_max at N, the one whose worst ratio over the whole
    family, normalised to the design level (ln q_k * N_k / N, so that an
    exactly similar family scores ln q_N at every level), is smallest:
    integer rounding of a short section can otherwise break a C1 junction
    at one level (WP 2.5b, U2/D2, channel, riser). lead = {level: cell}:
    the neighbouring cell before the edge start at each level, included in
    the family score (the channel's first cell against the slot-column
    cell, the riser's against the lid cell). fam_tol: a candidate is
    feasible only if that family score is <= ln(q_max) (1 + fam_tol), so
    the smallest feasible count leaves room for a level-robust split."""
    best = None
    nb_range = range(1, N) if grow_b else [0]
    for n_a in range(1, N + 1):
        for n_b in nb_range:
            n_u = N - n_a - n_b
            if n_u < 0 or n_a < MIN_SEC or 0 < n_u < MIN_SEC \
                    or 0 < n_b < MIN_SEC:
                continue

            def total(hu):
                ka = _kappa_c1(h_a, hu, n_a)
                kb = _kappa_c1(h_b, hu, n_b) if n_b else 0.0
                return (_sec_len(h_a, ka, n_a)
                        + (_sec_len(h_b, kb, n_b) if n_b else 0.0)
                        + n_u * hu), ka, kb

            h0 = max(h_a, h_b) if n_b else h_a
            if total(h0)[0] > L:
                continue
            lo, hi = h0, L
            for _ in range(80):
                hu = 0.5 * (lo + hi)
                if total(hu)[0] > L:
                    hi = hu
                else:
                    lo = hu
            hu = 0.5 * (lo + hi)
            if h_max is not None and hu > h_max * (1 + 1e-9):
                continue
            _, ka, kb = total(hu)
            la = _sec_len(h_a, ka, n_a)
            lb = _sec_len(h_b, kb, n_b) if n_b else 0.0
            secs = [(la / L, n_a / N, ka)]
            if n_u:
                secs.append((n_u * hu / L, n_u / N, 0.0))
            if n_b:
                secs.append((lb / L, n_b / N, -kb))
            m = Mapping(secs, name)
            q = max_ratio(m.cells(N, L))
            if family:
                vals = []
                for lv in LEVELS:
                    Nk = level_count(N, lv)
                    c = m.cells(Nk, L)
                    if lead is not None:
                        c = [lead[lv]] + c
                    vals.append(math.log(max_ratio(c)) * Nk / N)
                bad = q > q_max + 1e-9 or (
                    fam_tol is not None
                    and max(vals) > math.log(q_max) * (1 + fam_tol) + 1e-12)
                key = (bad, max(vals))
            elif grow_b:
                key = (q > q_max + 1e-9, q)
            else:
                key = (q > q_max + 1e-9, hu)
            if best is None or key < best[0]:
                best = (key, m)
    if best is None or best[0][0]:
        raise ValueError(f"{name}: N={N} too small for q<={q_max} "
                         f"(L={L}, ends {h_a}/{h_b}, h_max={h_max})")
    return best[1]


def grow_then_uniform(L, N, h_start, name, q_max=Q_MAX):
    """Geometric growth from h_start, then a uniform plateau (C1 junction);
    the smallest plateau cell whose max ratio stays <= q_max."""
    if N * h_start >= L:
        return one_sided(L, N, h_start, name)
    return _piecewise(L, N, h_start, None, name, q_max, None, grow_b=False)


def two_sided(L, N, h_a, h_b, name, q_max=Q_MAX, h_max=None, family=False):
    """Grow from first cell h_a, optional plateau (<= h_max), shrink to the
    last cell h_b; C1 junctions; the smoothest split with ratio <= q_max
    (family=True: smoothest over the level family, see _piecewise)."""
    return _piecewise(L, N, h_a, h_b, name, q_max, h_max, grow_b=True,
                      family=family)


def track(gmap, h_s, L, N, name, q_max=Q_MAX):
    """An edge that starts with a coarser cell h_s but keeps its points
    close to gmap's (index by index): uniform h_s cells until gmap's cells
    have grown to h_s (index m), then a fresh grow-then-uniform design from
    h_s over the remaining length and cells. The lead over gmap is at most
    sum_{i<m}(h_s - g_i) and all ratios stay <= Q_MAX."""
    cells = gmap.cells(N, L)
    if cells[0] >= h_s:
        return gmap
    m = next(i for i, c in enumerate(cells) if c >= h_s)
    rest = grow_then_uniform(L - m * h_s, N - m, h_s, name + "_rest",
                             q_max=q_max)
    fr = (L - m * h_s) / L
    nr = (N - m) / N
    sec = [(m * h_s / L, m / N, 0.0)]
    sec += [(f * fr, nu * nr, k) for f, nu, k in rest.sections]
    return Mapping(sec, name)


def grow_to_plateau(L, N, h_start, h_max, name, q_max=Q_MAX, family=False,
                    lead=None, fam_tol=None):
    """Geometric growth from h_start, then a uniform plateau of at most
    h_max (C1 junction); the smallest plateau with every ratio <= q_max
    (family=True: the split smoothest over the family, with the lead
    cells, see _piecewise)."""
    return _piecewise(L, N, h_start, None, name, q_max, h_max, grow_b=False,
                      family=family, lead=lead, fam_tol=fam_tol)


def track_scaled(gmap, h_s, L, N, name):
    """track() for a two-sided gmap (one that must END at a prescribed
    cell, e.g. the port column's): uniform h_s cells until gmap's cells have
    grown to h_s (index m), then gmap's own remaining cells, section by
    section, scaled by s = (L - m h_s) / (their length) <= 1. The edge ends
    where gmap ends, with gmap's shape and ratios (last cell s x gmap's), so
    the lines joining equal indices converge at the far end."""
    cells = gmap.cells(N, L)
    if cells[0] >= h_s:
        return gmap
    m = next(i for i, c in enumerate(cells) if c >= h_s)
    rest, a = [], 0
    for (_, _, k), n in zip(gmap.sections, gmap.counts(N)):
        b = a + n
        lo = max(a, m)
        if b > lo:
            nr = b - lo
            if nr < MIN_SEC:
                raise ValueError(f"{name}: tracking leaves a {nr}-cell "
                                 f"section (< {MIN_SEC})")
            rest.append((sum(cells[lo:b]), nr, nr * (k / n)))
        a = b
    s = (L - m * h_s) / sum(r[0] for r in rest)
    if not 0 < s <= 1 + 1e-12:
        raise ValueError(f"{name}: tracking scale {s}")
    secs = [(m * h_s / L, m / N, 0.0)]
    secs += [(ln * s / L, nr / N, kr) for ln, nr, kr in rest]
    return Mapping(secs, name)


FAM_TOL = 0.05   # channel / riser: family score <= 1.05 ln(Q_MAX)


def port_leads(g, n, d=None):
    """Cells before the start of the channel (slot-column cell x_port/U1)
    and of the riser (the D1 column's lid cell: the far layer-B cell with
    the z-lever, else the module mapping's) at every level: the lead cells
    of their family-aware designs."""
    Lzb = g["z_top"] - g["H"]
    zlever = True if d is None else bool(d.get("zlever"))
    if not zlever:
        zb = two_sided(Lzb, n["ZB"], d["h1_cap"], d["h_lid"], "zb")
    return dict(
        XL={lv: g["x_port"] / level_count(n["U1"], lv) for lv in LEVELS},
        ZR={lv: (Lzb / level_count(n["ZB"], lv) if zlever else
                 zb.cells(level_count(n["ZB"], lv), Lzb)[-1])
            for lv in LEVELS})


def level_count(n_medium, level):
    return max(1, int(round(n_medium * R_REFINE ** LEVELS[level])))


# ============================================================================
# 2-D multi-block topology
# ============================================================================

class Topology:
    def __init__(self, g, design, config="duct", port_mult=1):
        if config not in CONFIGS:
            raise ValueError(f"config {config!r} not in {CONFIGS}")
        self.g, self.d, self.config = g, design, config
        self.housing = config == "housing"
        self.zlever = bool(design.get("zlever", False))
        if self.zlever and not self.housing:
            raise ValueError("zlever is a housing-design option")
        self.zlever_module = bool(design.get("zlever_module", False))
        if self.zlever_module and not self.zlever:
            raise ValueError("zlever_module (the module-interior z-lever) "
                             "is a housing option on top of the far "
                             "z-lever (zlever)")
        self.vkind = {}          # tile-edge / strip-wall vertex -> G or C
        self.vcyl = set()        # cylinder vertices (axis, core, surface)
        # medium counts: the design classes (+ the port-length classes of
        # the housing, committed per port multiple in n_port)
        self.n_medium = dict(design["n"])
        self.port = None
        if self.housing:
            self.port = port_geometry(g, port_mult)
            key = f"{float(port_mult):g}"
            if key not in design["n_port"]:
                raise ValueError(f"port_mult {port_mult}: no committed "
                                 f"XL/ZR counts (have {sorted(design['n_port'])};"
                                 f" explore with --minimize)")
            self.n_medium.update(design["n_port"][key])
        self.v2d = []            # [(x, y)]
        self._vid = {}
        self.edges = {}          # (a, b) -> (class, Mapping) in a->b sense
        self.arcs = []           # (a, b, (cx, cy))
        self.blocks = []         # dict(v=(4 ids), ci, cj, zoneA, role, layers)
        # patch faces: name -> [("side", a, b, layer) | ("top", v, layer)]
        self.patches = {"inlet": [], "outlet": []}
        self.cyl = []            # per cylinder: dict of vertex ids
        self._build()

    # -- vertices / edges ----------------------------------------------------
    def vid(self, x, y):
        key = (round(x, 9), round(y, 9))
        if key not in self._vid:
            self._vid[key] = len(self.v2d)
            self.v2d.append((float(x), float(y)))
        return self._vid[key]

    def edge(self, a, b, cls, mapping):
        if (a, b) in self.edges:
            c0, m0 = self.edges[(a, b)]
            if c0 != cls or m0.sections != mapping.sections:
                raise ValueError(f"edge {a}->{b}: conflicting definitions "
                                 f"{c0}/{m0.name} vs {cls}/{mapping.name}")
            return
        if (b, a) in self.edges:
            c0, m0 = self.edges[(b, a)]
            inv = mapping.inverted()
            if c0 != cls or any(abs(p - q) > 1e-12 for s0, s1 in
                                zip(m0.sections, inv.sections)
                                for p, q in zip(s0, s1)):
                raise ValueError(f"edge {a}->{b}: conflicting definitions "
                                 f"{c0}/{m0.name} vs {cls}/{mapping.name}")
            return
        self.edges[(a, b)] = (cls, mapping)

    def edge_def(self, a, b):
        """(class, mapping in the a->b direction)."""
        if (a, b) in self.edges:
            return self.edges[(a, b)]
        cls, m = self.edges[(b, a)]
        return cls, m.inverted()

    def block(self, v, zoneA, role, layers=("A", "B")):
        v0, v1, v2, v3 = v
        ci = self.edge_def(v0, v1)[0]
        cj = self.edge_def(v1, v2)[0]
        for a, b, c in ((v3, v2, ci), (v0, v3, cj)):
            if self.edge_def(a, b)[0] != c:
                raise ValueError(f"block {v}: parallel edges differ in class")
        # right-handed in x-y (bottom face CCW seen from +z)
        (x0, y0), (x1, y1), (x3, y3) = (self.v2d[v0], self.v2d[v1],
                                        self.v2d[v3])
        if (x1 - x0) * (y3 - y0) - (y1 - y0) * (x3 - x0) <= 0:
            raise ValueError(f"block {v} ({role}) is left-handed")
        self.blocks.append(dict(v=tuple(v), ci=ci, cj=cj, zoneA=zoneA,
                                role=role, layers=tuple(layers)))

    # -- edge mapping designs (medium level) --------------------------------
    def _designs(self):
        g, d, n = self.g, self.d, self.n_medium
        R = g["D"] / 2
        P = (g["D"] + g["gap"]) / 2                  # tile half-width
        L_gap = P - R                                # 1 mm
        L_diag = math.sqrt(2) * P - R                # 5.142 mm
        m = {}
        m["ring_gap"] = one_sided(L_gap, n["R"], d["h1_wall"], "ring_gap")
        m["ring_diag"] = one_sided(L_diag, n["R"], d["h1_wall"], "ring_diag")
        # solid ring, core -> cylinder surface: the last cell (at the
        # surface) is h_solid_wall. Above the caps (layer B) these blocks are
        # fluid and meet the ring's first cell across r = R, so the ratio
        # h_solid_wall / h1_wall is the transition ratio there (D-I target
        # 1.5, limit 2.0); in layer A it resolves the solid side of the
        # interface.
        R_M = R - g["rho_core_mid"] * R
        R_K = R - g["rho_core_corner"] * R
        hs = d["h_solid_wall"]
        m["solid_M"] = one_sided(R_M, n["S"], hs, "solid_M").inverted()
        m["solid_K"] = one_sided(R_K, n["S"], hs, "solid_K").inverted()
        # last ring cell at the tile edge (G = mid-side, C = corner; at C the
        # size normal to the tile side is the diagonal cell * cos 45)
        hG = m["ring_gap"].cells(n["R"], L_gap)[-1]
        hC = m["ring_diag"].cells(n["R"], L_diag)[-1] / math.sqrt(2)
        h_side = P * 1.0 / n["TH"]                   # tile-side spacing
        self.h_edge = dict(G=hG, C=hC, side=h_side)
        L_st = g["y_half"] - g["rows"] * P
        m["strip_G"] = two_sided(L_st, n["ST"], hG * Q_MAX, d["h_side"],
                                 "strip_G")
        # at the tile corners the strip has more cells than it needs: one
        # geometric section from the wall cell, coarsening toward the tile
        m["strip_C"] = one_sided(L_st, n["ST"], d["h_side"],
                                 "strip_C").inverted()
        m["strip_C"].name = "strip_C"
        x_mod = g["cols"] * P
        L_up = (g["x_max"] - g["x_port"]) - x_mod
        # up/downstream x-edges, defined from the module (fine) end outward.
        # G (tile mid-sides, 1 mm from a cylinder) needs the fine start; the
        # C (tile corner) and W (side wall) edges start coarser but TRACK
        # the G distribution, so the y-lines joining equal indices stay
        # within ~6 mm of each other (measured non-orthogonality 31 deg at
        # medium, against 62 deg with independent gradings).
        h_port = g["x_port"] / n["U1"]
        h_c0 = min(hC * Q_MAX, h_side)
        if not self.housing:
            m["up_G"] = grow_then_uniform(L_up, n["U2"], hG * Q_MAX, "up_G",
                                          q_max=d["q_far"])
            # the far-field plateau runs into the uniform port column
            # (x_port / U1 cells): keep that junction within Q_MAX as well
            plat = m["up_G"].cells(n["U2"], L_up)[-1]
            if max(plat / h_port, h_port / plat) > Q_MAX + 1e-9:
                raise ValueError(f"U2/U1 junction ratio {plat / h_port:.3f} "
                                 f"(plateau {plat:.3f} mm, port "
                                 f"{h_port:.3f} mm)")
            m["up_C"] = track(m["up_G"], h_c0, L_up, n["U2"], "up_C",
                              q_max=d["q_far"])
        else:
            # housing: the port column is finer (U1 = 8, 1.25 mm) and now
            # carries the slot / riser flow, so U2/D2 grow from the module
            # to a plateau <= h_far_max and shrink back to the port cell
            m["up_G"] = two_sided(L_up, n["U2"], hG * Q_MAX, h_port, "up_G",
                                  q_max=d["q_far"], h_max=d["h_far_max"],
                                  family=True)
            plat = max(m["up_G"].cells(n["U2"], L_up))
            m["up_C"] = track_scaled(m["up_G"], h_c0, L_up, n["U2"], "up_C")
        self.h_edge.update(port=h_port, plateau=plat)
        m["up_W"] = m["up_C"]
        m["za"] = two_sided(g["H"], n["ZA"], d["h_floor"], d["h1_cap"], "za",
                            h_max=d["dz_max"])
        m["zb"] = two_sided(g["z_top"] - g["H"], n["ZB"], d["h1_cap"],
                            d["h_lid"], "zb")
        self.L = dict(gap=L_gap, diag=L_diag, strip=L_st, up=L_up,
                      za=g["H"], zb=g["z_top"] - g["H"])
        if self.housing:
            self._housing_designs(m, h_port)
        return m

    def _housing_designs(self, m, h_port):
        """Edge mappings of the D-G housing (medium design, sampled per
        level like every other mapping)."""
        g, d, n = self.g, self.d, self.n_medium
        P = (g["D"] + g["gap"]) / 2
        W = g["y_half"]
        # far-field y-distribution (vertex lines x <= -80, x >= 80: port
        # columns, channel, riser): side strip of width w (ST cells, wall
        # cell h_side) then the outer TH segment graded from the strip end
        # to the uniform tile-side spacing P / TH of the inner segments
        w = g["y_far_strip"]
        L_th = W - w - (g["rows"] - 1) * P            # 18 -> 10 mm: 8 mm
        m["strip_far"] = one_sided(w, n["ST"], d["h_side"],
                                   "strip_far").inverted()
        m["strip_far"].name = "strip_far"
        s_end = m["strip_far"].cells(n["ST"], w)[0]    # at the interior end
        h_in = P / n["TH"]
        best = None
        for i in range(1, 400):                        # first TH cell
            h0 = s_end * (0.5 + i / 200.0)
            try:
                mp = one_sided(L_th, n["TH"], h0, "th_far")
            except ValueError:
                continue
            c = mp.cells(n["TH"], L_th)
            q = max_ratio(list(reversed(
                m["strip_far"].cells(n["ST"], w))) + c + [h_in])
            if best is None or q < best[0] - 1e-12:
                best = (q, mp)
        m["th_far"] = best[1]
        self.far_y_ratio = best[0]
        # module-end strip lines x = +-30 mm (the U2/D2 blocks' module
        # face): the 3 mm strip shrinks toward the side wall at exactly
        # Q_MAX, which makes its tile-end cell as large as the ratio allows
        # (0.53 mm at medium, against 0.26 mm for strip_C) and cuts the jump
        # to the 1.25 mm tile-side cells from 4.8 to 2.4 (K0-design §7
        # item 4); the wall cell becomes 0.042 mm on these two lines only
        L_st = g["y_half"] - g["rows"] * P
        m["strip_Cend"] = fixed_ratio(n["ST"], Q_MAX ** (n["ST"] - 1),
                                      "strip_Cend").inverted()
        m["strip_Cend"].name = "strip_Cend"
        self.L.update(strip_Cend=L_st)
        # z: the per-edge lever. The far z-mapping has the same counts per
        # layer (conformal multi-block) but no cap-plane layer: layer B is
        # uniform (7 mm / ZB), layer A grows from the floor cell and shrinks
        # to that uniform cell at z = 65 mm (C1 at the layer junction).
        h_zb_far = (g["z_top"] - g["H"]) / n["ZB"]
        m["za_far"] = two_sided(g["H"], n["ZA"], d["h_floor"], h_zb_far,
                                "za_far", h_max=d["dz_max"], family=True)
        m["zb_far"] = Mapping([(1, 1, 0)], "zb_far")
        if self.zlever_module:
            # module-interior lever (§10.14): gap-throat lines h_cap_gap on
            # both sides of the cap plane; corner lines the uniform far
            # layer-B cell above and the lever mapping from it below (C1
            # across z = H at medium)
            Lzb = g["z_top"] - g["H"]
            m["za_gap"] = LeverMapping(m["za"], g["H"], n["ZA"],
                                       d["h_cap_gap"], "za_gap", flip=True)
            m["zb_gap"] = LeverMapping(m["zb"], Lzb, n["ZB"], d["h_cap_gap"],
                                       "zb_gap")
            h_c = d.get("h_cap_corner") or h_zb_far
            m["za_corner"] = LeverMapping(m["za"], g["H"], n["ZA"], h_c,
                                          "za_corner", flip=True)
            m["zb_corner"] = LeverMapping(m["zb"], Lzb, n["ZB"], h_c,
                                          "zb_corner", profile="geometric")
        # port column (and riser foot) cells at the lid / floor
        lid = (h_zb_far if self.zlever else
               m["zb"].cells(n["ZB"], g["z_top"] - g["H"])[-1])
        floor = d["h_floor"]
        L = self.port["L"]
        # channel height: floor/ceiling first cell h_ch_wall; the ceiling
        # (z = 0) cell meets the housing floor cell across the slot
        m["zc"] = two_sided(g["h_ch"], n["ZC"], d["h_ch_wall"], floor, "zc",
                            h_max=d["h_ch_max"], family=True)
        # channel length (x): starts at the slot-column cell (ratio 1 at
        # medium: the U1 and XL counts round differently per level, and a
        # 1.2 start leaves no margin, e.g. 1.23 at fine), grows upstream to
        # the streamwise plateau h_run; split smoothest over the family
        # including the junction with the slot-column cell at every level
        m["chan_x"] = grow_to_plateau(
            L, n["XL"], h_port, d["h_run"], "chan_x", family=True,
            lead=port_leads(g, n)["XL"], fam_tol=FAM_TOL)
        # riser length (z, from the lid up): starts at the lid cell of the
        # D1 column (ratio 1 at medium), grows to h_run; smoothest over the
        # family including the lid cell at every level
        m["zr"] = grow_to_plateau(L, n["ZR"], lid, d["h_run"], "zr",
                                  family=True, fam_tol=FAM_TOL,
                                  lead=port_leads(g, n, d)["ZR"])
        self.L.update(far_strip=w, far_th=L_th, ch=g["h_ch"], port=L)

    # -- build --------------------------------------------------------------
    def _build(self):
        g = self.g
        R = g["D"] / 2
        P = (g["D"] + g["gap"]) / 2
        pitch = 2 * P
        rows, cols = g["rows"], g["cols"]
        m = self.maps = self._designs()
        rM, rK = g["rho_core_mid"] * R, g["rho_core_corner"] * R
        x_mod, y_mod = cols * P, rows * P       # module half-extents (tiles)
        W = g["y_half"]

        centres = [(-(cols - 1) * pitch / 2 + c * pitch,
                    -(rows - 1) * pitch / 2 + r * pitch)
                   for r in range(rows) for c in range(cols)]
        self.centres = centres
        for k, (cx, cy) in enumerate(centres):
            ang = [math.radians(45 * s) for s in range(8)]
            O = self.vid(cx, cy)
            B = [self.vid(cx + (rM if s % 2 == 0 else rK) * math.cos(a),
                          cy + (rM if s % 2 == 0 else rK) * math.sin(a))
                 for s, a in enumerate(ang)]
            A = [self.vid(cx + R * math.cos(a), cy + R * math.sin(a))
                 for a in ang]
            T = []
            for s, a in enumerate(ang):
                c, sn = round(math.cos(a), 12), round(math.sin(a), 12)
                if s % 2 == 0:
                    T.append(self.vid(cx + P * round(c), cy + P * round(sn)))
                else:
                    T.append(self.vid(cx + P * math.copysign(1, c),
                                      cy + P * math.copysign(1, sn)))
            self.cyl.append(dict(O=O, B=B, A=A, T=T, centre=(cx, cy)))
            self.vcyl.update([O] + B + A)
            for s, v in enumerate(T):      # mid-side = gap throat G
                self.vkind[v] = "G" if s % 2 == 0 else "C"
            zone = f"cell{k}"
            # edges
            for s in range(8):
                t = (s + 1) % 8
                self.edge(B[s], B[t], "TH", UNIFORM)                 # core side
                self.edge(B[s], A[s], "S",                           # solid radial
                          m["solid_M"] if s % 2 == 0 else m["solid_K"])
                self.edge(A[s], A[t], "TH", UNIFORM)                 # arc
                self.arcs.append((A[s], A[t], (cx, cy)))
                self.edge(A[s], T[s], "R",
                          m["ring_gap"] if s % 2 == 0 else m["ring_diag"])
                self.edge(T[s], T[t], "TH", UNIFORM)                 # tile side
            for q in range(4):
                self.edge(O, B[2 * q], "TH", UNIFORM)
            # blocks: core quadrants, solid ring, fluid ring
            for q in range(4):
                self.block((O, B[2 * q], B[2 * q + 1], B[(2 * q + 2) % 8]),
                           zone, f"core{k}")
            for s in range(8):
                t = (s + 1) % 8
                self.block((B[s], A[s], A[t], B[t]), zone, f"solidring{k}")
            for s in range(8):
                t = (s + 1) % 8
                self.block((A[s], T[s], T[t], A[t]), "fluid", f"ring{k}")

        # side strips (tile edge y = +-y_mod -> side wall y = +-W)
        xs = [-x_mod + i * P for i in range(2 * cols + 1)]    # G/C alternate
        kind = {round(x, 9): ("C" if i % 2 == 0 else "G")
                for i, x in enumerate(xs)}
        for sgn in (-1, 1):
            yt, yw = sgn * y_mod, sgn * W
            for x in xs:
                a, b = self.vid(x, yt), self.vid(x, yw)
                self.vkind[b] = kind[round(x, 9)]    # side-wall line
                key = f"strip_{kind[round(x, 9)]}"
                if self.housing and abs(abs(x) - x_mod) < 1e-9:
                    key = "strip_Cend"      # module-end lines (U2/D2 face)
                self.edge(a, b, "ST", m[key])
            for x0, x1 in zip(xs[:-1], xs[1:]):
                self.edge(self.vid(x0, yw), self.vid(x1, yw), "TH", UNIFORM)
                if sgn < 0:
                    v = (self.vid(x0, yw), self.vid(x1, yw),
                         self.vid(x1, yt), self.vid(x0, yt))
                else:
                    v = (self.vid(x0, yt), self.vid(x1, yt),
                         self.vid(x1, yw), self.vid(x0, yw))
                self.block(v, "fluid", "strip")

        # up/downstream H-blocks: y breakpoints = walls, strips, tile half-sides
        ys = [-W] + [-y_mod + i * P for i in range(2 * rows + 1)] + [W]
        ytype = {}
        for i, y in enumerate(ys):
            if abs(abs(y) - W) < 1e-9:
                ytype[i] = "W"
            else:
                j = i - 1            # index within the tile breakpoints
                ytype[i] = "C" if j % 2 == 0 else "G"
        ycls = ["ST" if (abs(abs(y0) - W) < 1e-9 or abs(abs(y1) - W) < 1e-9)
                else "TH" for y0, y1 in zip(ys[:-1], ys[1:])]
        # far-field vertex lines (x <= -80, x >= 80): the duct keeps the
        # module breakpoints; the housing moves the strip boundary to
        # +-(W - y_far_strip) and grades the outer TH segment (the U2/D2
        # blocks bridge the two with straight, slightly slanted edges)
        if self.housing:
            wf = g["y_far_strip"]
            ys_far = [-W, -W + wf] + ys[2:-2] + [W - wf, W]
        else:
            ys_far = ys
        self.ys, self.ys_far = ys, ys_far
        xa_up, xb_up = g["x_min"], g["x_min"] + g["x_port"]
        xa_dn, xb_dn = g["x_max"] - g["x_port"], g["x_max"]
        x_in = xa_up - self.port["L"] if self.housing else None

        def y_at(x, i):          # row i breakpoint on the vertex line x
            return ys[i] if abs(abs(x) - x_mod) < 1e-9 else ys_far[i]

        for side in ("up", "dn"):
            if side == "up":
                cols_x = [(xa_up, xb_up, "U1"), (xb_up, -x_mod, "U2")]
                far_x = [xa_up, xb_up]
                if self.housing:          # inlet channel (layer C only)
                    cols_x.insert(0, (x_in, xa_up, "XL"))
                    far_x.insert(0, x_in)
            else:
                cols_x = [(x_mod, xa_dn, "D2"), (xa_dn, xb_dn, "D1")]
                far_x = [xa_dn, xb_dn]
            # y-edges on the far x-lines (port column / channel faces)
            for x in far_x:
                for j, ((y0, y1), c) in enumerate(zip(
                        zip(ys_far[:-1], ys_far[1:]), ycls)):
                    if c == "ST":   # strip segment: interior end -> wall
                        yt, yw = (y1, y0) if y0 < 0 else (y0, y1)
                        self.edge(self.vid(x, yt), self.vid(x, yw), "ST",
                                  m["strip_far" if self.housing
                                    else "strip_C"])
                    elif self.housing and j in (1, len(ycls) - 2):
                        # outer TH segment: fine end at the strip boundary
                        yf, yi = (y0, y1) if y0 < 0 else (y1, y0)
                        self.edge(self.vid(x, yf), self.vid(x, yi), "TH",
                                  m["th_far"])
                    else:
                        self.edge(self.vid(x, y0), self.vid(x, y1), "TH",
                                  UNIFORM)
            for xa, xb, cls in cols_x:
                for i in range(len(ys)):
                    a = self.vid(xa, y_at(xa, i))
                    b = self.vid(xb, y_at(xb, i))
                    if cls in ("U1", "D1"):
                        self.edge(a, b, cls, UNIFORM)
                    elif cls == "XL":      # fine end at the slot column
                        self.edge(b, a, cls, m["chan_x"])
                    elif cls == "U2":      # fine end at xb (= -x_mod)
                        self.edge(b, a, cls, m[f"up_{ytype[i]}"])
                    else:                  # D2: fine end at xa (= +x_mod)
                        self.edge(a, b, cls, m[f"up_{ytype[i]}"])
                if cls == "XL":
                    layers = ("C",)
                elif self.housing and cls == "U1":
                    layers = ("C", "A", "B")      # elbow under the slot
                elif self.housing and cls == "D1":
                    layers = ("A", "B", "R")      # riser above the slot
                else:
                    layers = ("A", "B")
                for i in range(len(ys) - 1):
                    v = (self.vid(xa, y_at(xa, i)), self.vid(xb, y_at(xb, i)),
                         self.vid(xb, y_at(xb, i + 1)),
                         self.vid(xa, y_at(xa, i + 1)))
                    self.block(v, "fluid", f"{side}_{cls}", layers)
                    if not self.housing:
                        if cls == "U1":
                            # x = x_min face, written (a_bot a_top b_top
                            # b_bot): outward (-x) normal as in the
                            # template's inlet
                            for lay in ("A", "B"):
                                self.patches["inlet"].append(
                                    ("side", v[0], v[3], lay))
                        if cls == "D1":
                            for lay in ("A", "B"):       # +x normal
                                self.patches["outlet"].append(
                                    ("side", v[2], v[1], lay))
                    elif cls == "XL":
                        # channel entry x = x_min - L (-x normal)
                        self.patches["inlet"].append(("side", v[0], v[3], "C"))
                    elif cls == "D1":
                        # riser top z = z_top + L (+z normal: the 2-D block
                        # is counter-clockwise seen from +z)
                        self.patches["outlet"].append(("top", v, "R"))

    # -- counts ---------------------------------------------------------------
    def counts(self, level):
        return {c: level_count(n, level) for c, n in self.n_medium.items()}

    # -- z-levels, layers and per-vertex z-mappings ---------------------------
    def zlevels(self):
        g = self.g
        zs = [0.0, g["H"], g["z_top"]]
        if self.housing:
            zs = [-g["h_ch"]] + zs + [g["z_top"] + self.port["L"]]
        return zs

    def layers(self):
        """[dict(name, cls, k0, k1)]: k0/k1 index zlevels()."""
        if not self.housing:
            return [dict(name="A", cls="ZA", k0=0, k1=1),
                    dict(name="B", cls="ZB", k0=1, k1=2)]
        return [dict(name="C", cls="ZC", k0=0, k1=1),
                dict(name="A", cls="ZA", k0=1, k1=2),
                dict(name="B", cls="ZB", k0=2, k1=3),
                dict(name="R", cls="ZR", k0=3, k1=4)]

    def is_far(self, v):
        """Vertex line outside the module + strips (|x| > x_mod): carries
        the far z-mapping when the z-lever is on."""
        g = self.g
        x_mod = g["cols"] * (g["D"] + g["gap"]) / 2
        return abs(self.v2d[v][0]) > x_mod + 1e-9

    def zmap(self, layer, v):
        """z-mapping of the vertical edge on 2-D vertex v in a layer."""
        if layer == "C":
            return self.maps["zc"]
        if layer == "R":
            return self.maps["zr"]
        key = {"A": "za", "B": "zb"}[layer]
        if self.zlever and self.is_far(v):
            return self.maps[key + "_far"]
        if self.zlever_module and v not in self.vcyl:
            # module-interior lever: tile-edge and strip-wall lines (a
            # KeyError here is a topology bug: every non-cylinder line of
            # the module region has a kind)
            return self.maps[key + {"G": "_gap", "C": "_corner"}[
                self.vkind[v]]]
        return self.maps[key]

    def zline_kind(self, v):
        """'far', 'gap', 'corner' or 'module': the z-mapping family of the
        vertical line on 2-D vertex v (layers A and B)."""
        if self.zlever and self.is_far(v):
            return "far"
        if self.zlever_module and v not in self.vcyl:
            return {"G": "gap", "C": "corner"}[self.vkind[v]]
        return "module"


# ============================================================================
# blockMeshDict writer + report
# ============================================================================

HEX_EDGES = [(0, 1), (3, 2), (7, 6), (4, 5),      # x (i)
             (0, 3), (1, 2), (5, 6), (4, 7),      # y (j)
             (0, 4), (1, 5), (2, 6), (3, 7)]      # z (k)


def render(topo, level):
    g = topo.g
    n = topo.counts(level)
    zs = topo.zlevels()
    layers = topo.layers()
    lay_by = {lay["name"]: lay for lay in layers}
    # global vertex ids for the (2-D vertex, z-level) pairs that blocks use,
    # z-level major (the duct uses every pair: its numbering is the 2.5a one)
    used = set()
    for b in topo.blocks:
        for ln in b["layers"]:
            lay = lay_by[ln]
            for v in b["v"]:
                used.update({(v, lay["k0"]), (v, lay["k1"])})
    gid = {}
    for k in range(len(zs)):
        for i in range(len(topo.v2d)):
            if (i, k) in used:
                gid[(i, k)] = len(gid)
    out = []
    w = out.append
    w("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
      "    class       dictionary;\n    object      blockMeshDict;\n}\n")
    head = f"// K0 structured O/H family (plan D-I, D-C), level: {level}\n"
    if topo.housing:
        head = (f"// K0 structured O/H family (plan D-I, D-C) on the D-G "
                f"housing, level: {level}\n"
                f"// ports: L = {topo.port['L']:.6g} mm "
                f"(= {topo.port['port_mult']:g} x 10 D_h, D_h = "
                f"{topo.port['Dh_slot']:.6g} mm); z-lever "
                f"{'on' if topo.zlever else 'off'}"
                + (f"; module-interior z-lever on (gap-throat cap-plane "
                   f"cell {topo.d['h_cap_gap']:g} mm at medium)"
                   if topo.zlever_module else "") + "\n")
    w(f"// Generated by scripts/mesh/k0_blockmesh.py - do not edit.\n"
      + head +
      f"// counts: " + " ".join(f"{c}={v}" for c, v in n.items()) + "\n\n")
    w("scale   0.001;\n\nvertices\n(\n")
    for (i, k), gi in sorted(gid.items(), key=lambda t: t[1]):
        x, y = topo.v2d[i]
        w(f"    ({x:.9g} {y:.9g} {zs[k]:.9g})  // {gi}\n")
    w(");\n\nblocks\n(\n")
    ncells = {}
    for lay in layers:
        for b in topo.blocks:
            if lay["name"] not in b["layers"]:
                continue
            v = b["v"]
            hexv = [gid[(x, lay["k0"])] for x in v] + \
                [gid[(x, lay["k1"])] for x in v]
            zone = b["zoneA"] if lay["name"] == "A" else "fluid"
            ni, nj, nk = n[b["ci"]], n[b["cj"]], n[lay["cls"]]
            ncells[zone] = ncells.get(zone, 0) + ni * nj * nk
            grads = []
            for e, (p0, p1) in enumerate(HEX_EDGES):
                if e < 8:        # in-plane edge: 2-D ids
                    a, bb = v[p0 % 4], v[p1 % 4]
                    cls, mp = topo.edge_def(a, bb)
                    grads.append(mp.grading(n[cls]))
                else:            # vertical edge on 2-D vertex v[e - 8]
                    grads.append(topo.zmap(lay["name"], v[e - 8]).grading(nk))
            w(f"    hex ({' '.join(map(str, hexv))}) {zone} ({ni} {nj} {nk})\n"
              f"        edgeGrading ({' '.join(grads)})"
              f"  // {b['role']} {lay['name']}\n")
    w(");\n\nedges\n(\n")
    for a, b, (cx, cy) in topo.arcs:
        for k, z in enumerate(zs):
            if (a, k) in gid and (b, k) in gid:
                w(f"    arc {gid[(a, k)]} {gid[(b, k)]} origin "
                  f"({cx:.9g} {cy:.9g} {z:.9g})\n")
    w(");\n\n")
    if topo.housing:
        w("defaultPatch\n{\n    name    walls;     // floor, lid, side and "
          "end walls, channel and riser walls, cell bottom caps "
          "(adiabatic)\n    type    wall;\n}\n\n")
    else:
        w("defaultPatch\n{\n    name    walls;     // floor, lid, side walls, "
          "cell bottom caps (adiabatic)\n    type    wall;\n}\n\n")
    w("boundary\n(\n")
    for name in ("inlet", "outlet"):
        w(f"    {name}\n    {{\n        type patch;\n        faces\n        (\n")
        for f in topo.patches[name]:
            if f[0] == "side":
                _, a, b, ln = f
                k0, k1 = lay_by[ln]["k0"], lay_by[ln]["k1"]
                ids = (gid[(a, k0)], gid[(a, k1)], gid[(b, k1)], gid[(b, k0)])
            else:
                _, v, ln = f
                k1 = lay_by[ln]["k1"]
                ids = tuple(gid[(x, k1)] for x in v)
            w("            (" + " ".join(map(str, ids)) + ")\n")
        w("        );\n    }\n")
    w(");\n\nmergePatchPairs\n(\n);\n")
    return "".join(out), ncells


def design_report(topo, level, nu=None):
    """Design (not measured) numbers for one level. Near-wall counts use the
    D-I thresholds of k0_targets: film at its 0.37 mm lower bound, Stokes
    layers sqrt(nu/(pi f)) at nu (default mu/rho from properties.yaml) and
    the D-I gate thickness (thinner of that and the quoted 0.18 mm at 4 Hz).
    Counts are at the design block edges; blockMesh's mid-sector cells are
    thinner, so the measured minima (k0_measure) can only be higher."""
    nu = T.coolant_nu() if nu is None else nu
    stokes = {f: T.stokes_thickness(f, nu) * 1e3 for f in T.REPORT_FREQS_HZ}
    sgate = {f: T.stokes_gate_thickness(f, nu) * 1e3
             for f in T.REPORT_FREQS_HZ}
    n = topo.counts(level)
    L = topo.L
    rep = dict(level=level, counts=n)
    rep["nearwall_targets"] = dict(
        nu_m2s=nu, film_mm=T.FILM_MM,
        stokes_mm={f"{f:g}Hz": v for f, v in stokes.items()},
        stokes_gate_mm={f"{f:g}Hz": v for f, v in sgate.items()})
    rep["ratio_to_coarser"] = {}
    if LEVELS[level] > -1:
        prev = [k for k, v in LEVELS.items() if v == LEVELS[level] - 1][0]
        npv = topo.counts(prev)
        rep["ratio_to_coarser"] = {c: n[c] / npv[c] for c in n}

    def edge_stats(mp, N, Lm):
        cells = mp.cells(N, Lm)
        qs = [max(b / a, a / b) for a, b in zip(cells[:-1], cells[1:])]

        def within(dist):
            cum, k = 0.0, 0
            for h in cells:
                cum += h
                if cum <= dist * (1 + 1e-9):
                    k += 1
            return k

        def end_of(k):          # wall distance where the k-th cell ends
            return sum(cells[:k])
        return dict(first=cells[0], last=cells[-1], max_ratio=max(qs or [1]),
                    cells_in_film=within(T.FILM_MM),
                    film_cellN_end=end_of(T.FILM_MIN_CELLS),
                    cells_in_stokes={f"{f:g}Hz": within(d)
                                     for f, d in stokes.items()},
                    cells_in_stokes_gate={f"{f:g}Hz": within(d)
                                          for f, d in sgate.items()},
                    stokes_cellN_end=end_of(T.STOKES_MIN_CELLS))

    rep["ring_gap"] = edge_stats(topo.maps["ring_gap"], n["R"], L["gap"])
    rep["ring_diag"] = edge_stats(topo.maps["ring_diag"], n["R"], L["diag"])
    rep["cells_across_gap"] = 2 * n["R"]
    rep["cap_zb"] = edge_stats(topo.maps["zb"], n["ZB"], L["zb"])
    za = topo.maps["za"].cells(n["ZA"], L["za"])
    rep["za"] = dict(first_floor=za[0], last_cap=za[-1], max_cell=max(za),
                     max_ratio=max(max(b / a, a / b)
                                   for a, b in zip(za[:-1], za[1:])))
    for k in ("strip_G", "strip_C", "up_G", "up_C", "up_W"):
        N = n["ST"] if k.startswith("strip") else n["U2"]
        Lk = L["strip"] if k.startswith("strip") else L["up"]
        c = topo.maps[k].cells(N, Lk)
        rep[k] = dict(first=c[0], last=c[-1], max_cell=max(c),
                      max_ratio=max(max(b / a, a / b)
                                    for a, b in zip(c[:-1], c[1:])))
    for k in ("solid_M", "solid_K"):
        R = topo.g["D"] / 2
        rho = topo.g["rho_core_mid" if k == "solid_M" else "rho_core_corner"]
        c = topo.maps[k].cells(n["S"], R - rho * R)
        rep[k] = dict(first_core=c[0], last_wall=c[-1],
                      max_ratio=max_ratio(c))
    # radial transition across r = R above the caps (layer B, both fluid)
    rep["capfoot_transition"] = rep["solid_M"]["last_wall"] / min(
        rep["ring_gap"]["first"], rep["ring_diag"]["first"])
    rep["arc_spacing"] = (math.pi / 4) * topo.g["D"] / 2 / n["TH"]
    rep["tile_side_spacing"] = (topo.g["D"] + topo.g["gap"]) / 2 / n["TH"]
    _, cells = render(topo, level)
    rep["cells"] = cells
    rep["cells_total"] = sum(cells.values())
    if topo.housing:
        rep["housing"] = housing_report(topo, level)
    return rep


def housing_report(topo, level):
    """Design numbers of the D-G housing additions at one level."""
    g, n, m, L = topo.g, topo.counts(level), topo.maps, topo.L

    def stats(mp, N, Lm):
        c = mp.cells(N, Lm)
        return dict(first=c[0], last=c[-1], max_cell=max(c),
                    max_ratio=max_ratio(c))
    rep = dict(port=dict(topo.port), zlever=topo.zlever,
               y_far_strip_mm=g["y_far_strip"])
    rep["chan_x"] = stats(m["chan_x"], n["XL"], L["port"])
    rep["zr"] = stats(m["zr"], n["ZR"], L["port"])
    rep["zc"] = stats(m["zc"], n["ZC"], g["h_ch"])
    for k in ("up_G", "up_C"):
        rep[k] = stats(m[k], n["U2"], L["up"])
    # U2/D2 tracking: largest distance between the G- and C-row points of
    # equal index (the y-lines joining them tilt by this over 10 mm)
    xg = [0.0]
    for c in m["up_G"].cells(n["U2"], L["up"]):
        xg.append(xg[-1] + c)
    xc = [0.0]
    for c in m["up_C"].cells(n["U2"], L["up"]):
        xc.append(xc[-1] + c)
    rep["up_track_max_mm"] = max(abs(a - b) for a, b in zip(xg, xc))
    # far-line y-distribution, wall -> centre, and its junctions
    s = list(reversed(m["strip_far"].cells(n["ST"], L["far_strip"])))
    t = m["th_far"].cells(n["TH"], L["far_th"])
    h_in = (g["D"] + g["gap"]) / 2 / n["TH"]
    rep["far_y"] = dict(wall=s[0], strip_end=s[-1], th_first=t[0],
                        th_last=t[-1], inner=h_in,
                        max_ratio=max_ratio(s + t + [h_in]))
    se = m["strip_Cend"].cells(n["ST"], L["strip_Cend"])
    rep["strip_Cend"] = dict(tile_end=se[0], wall=se[-1],
                             max_ratio=max_ratio(se),
                             jump_to_tile_side=h_in / se[0])
    # z: module vs far mapping at the cap plane and the layer ends
    za = m["za"].cells(n["ZA"], g["H"])
    zb = m["zb"].cells(n["ZB"], g["z_top"] - g["H"])
    zaf = m["za_far"].cells(n["ZA"], g["H"])
    zbf = m["zb_far"].cells(n["ZB"], g["z_top"] - g["H"])
    rep["capplane"] = dict(module=min(za[-1], zb[0]),
                           far=min(zaf[-1], zbf[0]),
                           za_far_max_ratio=max_ratio(zaf + zbf),
                           za_far_max_cell=max(zaf))
    rep["capplane"]["lever_ratio"] = rep["capplane"]["far"] / \
        rep["capplane"]["module"]
    # outflow column (D1 + riser): smallest z-cell; Co of the mean and of
    # the developed-laminar peak Q 500 outflow (k0_targets): the mean at
    # dt*, the peak at the level's Q 500 time step (the gated number), and
    # the largest dt that keeps the peak at Co 1
    col_z = (zaf + zbf) if topo.zlever else (za + zb)
    col_z = col_z + m["zr"].cells(n["ZR"], L["port"])
    w = T.CO_OUTFLOW_Q_M3S / (g["x_port"] * 2 * g["y_half"] * 1e-6)
    dz = min(col_z) * 1e-3
    dt_ref = T.co_outflow_dt(level)
    rep["outflow_column"] = dict(
        dz_min_mm=min(col_z), w_mean_mm_s=1e3 * w,
        co_mean_outflow=w * T.CO_OUTFLOW_DT_S / dz,
        peak_factor=T.CO_OUTFLOW_PEAK_FACTOR, dt_ref_s=dt_ref,
        co_peak_outflow=T.CO_OUTFLOW_PEAK_FACTOR * w * dt_ref / dz,
        dt_max_peak_s=T.CO_OUTFLOW_MAX * dz / (T.CO_OUTFLOW_PEAK_FACTOR * w))
    # cells: channel (+ slot column below the floor) and riser
    ycount = 2 * g["rows"] * n["TH"] + 2 * n["ST"]
    rep["cells_channel"] = ycount * n["ZC"] * (n["XL"] + n["U1"])
    rep["cells_riser"] = ycount * n["D1"] * n["ZR"]
    rep["zlever_module"] = zlever_module_report(topo, level)
    return rep


def _cum(c):
    out = [0.0]
    for x in c:
        out.append(out[-1] + x)
    return out


def zlever_module_report(topo, level):
    """Design numbers of the module-interior z-lever at one level (with
    the lever off: the same numbers for the uniform cap layer, the P00b
    baseline). Per line family (gap throat G, corner C): the cap-plane
    cells below / above z = H, the largest offset D of the line's points
    from the cylinder (module) line's, counted from the cap plane, and
    the tilt atan(D / L_ring) of the z-lines across the fluid ring (L_ring
    1 mm on the gap ray, 5.14 mm on the diagonal). Plume stations: the
    cap-plane cell at wall distance d on the two rays, from the linear
    radial blend of blockCreate.C, (1 - d/L_ring) h_cyl + (d/L_ring) h_line
    (the cell centred at d; k0_measure measures the mesh), with the P00b
    plume's Courant number (k0_targets.plume_courant)."""
    g, n, m = topo.g, topo.counts(level), topo.maps
    H, Lzb = g["H"], g["z_top"] - g["H"]
    R, P = g["D"] / 2, (g["D"] + g["gap"]) / 2
    ring = dict(gap=P - R, corner=math.sqrt(2) * P - R)
    za = m["za"].cells(n["ZA"], H)
    zb = m["zb"].cells(n["ZB"], Lzb)
    cap = min(za[-1], zb[0])
    rep = dict(on=topo.zlever_module, module_cap_mm=cap,
               h_cap_gap_design_mm=(topo.d.get("h_cap_gap")
                                    if topo.zlever_module else None),
               tilt_max_design_deg=T.ZLEVER_TILT_MAX_DEG, lines={})
    for kind in ("gap", "corner"):
        if topo.zlever_module:
            a = m[f"za_{kind}"].cells(n["ZA"], H)
            b = m[f"zb_{kind}"].cells(n["ZB"], Lzb)
        else:
            a, b = za, zb
        dA = max(abs(p - q) for p, q in zip(_cum(a[::-1]), _cum(za[::-1])))
        dB = max(abs(p - q) for p, q in zip(_cum(b), _cum(zb)))
        D = max(dA, dB)
        rep["lines"][kind] = dict(
            cap_below_mm=a[-1], cap_above_mm=b[0], offset_max_mm=D,
            offset_max_layer_A_mm=dA, offset_max_layer_B_mm=dB,
            ring_mm=ring[kind],
            tilt_deg=math.degrees(math.atan(D / ring[kind])),
            max_ratio=max_ratio(a + b), min_cell_mm=min(a + b),
            lid_mm=b[-1])
    st = []
    for kind, ds in (("gap", (0.25, 0.5, 1.0)),
                     ("corner", (0.5, 1.0, 2.0, 3.0, ring["corner"]))):
        ln = rep["lines"][kind]
        h_line = min(ln["cap_below_mm"], ln["cap_above_mm"])
        for dd in ds:
            lam = min(dd / ring[kind], 1.0)
            dz = (1 - lam) * cap + lam * h_line
            st.append(dict(ray=kind, d_mm=dd, **T.plume_courant(dz)))
    rep["plume_stations"] = st
    return rep


def table(topo, nu=None):
    reps = {lv: design_report(topo, lv, nu) for lv in LEVELS}
    rows = []
    rows.append("level   " + "".join(f"{c:>6}" for c in topo.n_medium)
                + "   fluid     solid(6)    total")
    for lv, r in reps.items():
        solid = sum(v for k, v in r["cells"].items() if k != "fluid")
        rows.append(f"{lv:<8}" + "".join(f"{r['counts'][c]:>6}"
                                          for c in topo.n_medium)
                    + f"  {r['cells']['fluid']:>8} {solid:>10} "
                      f"{r['cells_total']:>9}")
    rows.append("")
    t = reps["medium"]["nearwall_targets"]
    rows.append(f"near-wall counts at the block edges (diag / cap): film "
                f"{t['film_mm']:.2f} mm; Stokes 2 Hz {t['stokes_mm']['2Hz']:.4f}"
                f" mm; 4 Hz D-I gate {t['stokes_gate_mm']['4Hz']:.4f} mm "
                f"(computed {t['stokes_mm']['4Hz']:.4f}); nu {t['nu_m2s']:.4e}")
    rows.append("level    h1_gap  h1_diag  q_diag  h1_cap q_cap  film  S2Hz "
                "S4Hz   arc   q_za   q_strip q_up   r_N(fluid) r_N(cell0)")
    prev = None
    for lv, r in reps.items():
        rN_f = rN_c = float("nan")
        if prev:
            rN_f = (r["cells"]["fluid"] / prev["cells"]["fluid"]) ** (1 / 3)
            rN_c = (r["cells"]["cell0"] / prev["cells"]["cell0"]) ** (1 / 3)
        d, c = r["ring_diag"], r["cap_zb"]
        rows.append(
            f"{lv:<8} {r['ring_gap']['first']:.4f}  {d['first']:.4f}"
            f"  {d['max_ratio']:.3f}   {c['first']:.4f} {c['max_ratio']:.3f}"
            f"  {d['cells_in_film']}/{c['cells_in_film']}"
            f"   {d['cells_in_stokes_gate']['2Hz']}/{c['cells_in_stokes_gate']['2Hz']}"
            f"  {d['cells_in_stokes_gate']['4Hz']}/{c['cells_in_stokes_gate']['4Hz']}"
            f"   {r['arc_spacing']:.3f}  {r['za']['max_ratio']:.3f}  "
            f"{max(r['strip_G']['max_ratio'], r['strip_C']['max_ratio']):.3f}"
            f"   {max(r['up_G']['max_ratio'], r['up_C']['max_ratio'], r['up_W']['max_ratio']):.3f}"
            f"   {rN_f:.4f}     {rN_c:.4f}")
        prev = r
    return "\n".join(rows), reps


def min_counts(g, d, q_solid=1.31, config="duct"):
    """Smallest medium counts of the graded classes that satisfy the design
    targets (every ratio <= Q_MAX, q_far in the H-blocks, q_solid in the
    solid ring). TH, U1, D1 are free design choices and are kept. For the
    housing also U2/D2 (two-sided to the port cell, plateau <= h_far_max),
    ZC, and XL/ZR for every port multiple in d["n_port"] (returned under
    "n_port")."""
    R = g["D"] / 2
    P = (g["D"] + g["gap"]) / 2
    out = {}

    def first_ok(fn, lo=2, hi=400):
        for N in range(lo, hi):
            try:
                if fn(N):
                    return N
            except ValueError:
                pass
        raise ValueError("no count found")

    L_gap, L_diag = P - R, math.sqrt(2) * P - R
    out["R"] = first_ok(lambda N: max_ratio(one_sided(
        L_diag, N, d["h1_wall"], "d").cells(N, L_diag)) <= Q_MAX + 1e-9)
    Ls = [R - g["rho_core_mid"] * R, R - g["rho_core_corner"] * R]
    out["S"] = first_ok(lambda N: all(max_ratio(one_sided(
        L, N, d["h_solid_wall"], "s").cells(N, L)) <= q_solid for L in Ls))
    hG = one_sided(L_gap, out["R"], d["h1_wall"], "g").cells(
        out["R"], L_gap)[-1]
    hC = one_sided(L_diag, out["R"], d["h1_wall"], "d").cells(
        out["R"], L_diag)[-1] / math.sqrt(2)
    L_st = g["y_half"] - g["rows"] * P
    out["ST"] = first_ok(lambda N: bool(two_sided(
        L_st, N, hG * Q_MAX, d["h_side"], "st")))
    L_up = (g["x_max"] - g["x_port"]) - g["cols"] * P
    h_side = P / d["n"]["TH"]

    h_port = g["x_port"] / d["n"]["U1"]

    def up_ok(N):
        gm = grow_then_uniform(L_up, N, hG * Q_MAX, "u", q_max=d["q_far"])
        track(gm, min(hC * Q_MAX, h_side), L_up, N, "c", q_max=d["q_far"])
        plat = gm.cells(N, L_up)[-1]
        return max(plat / h_port, h_port / plat) <= Q_MAX + 1e-9

    def up_ok_housing(N):
        gm = two_sided(L_up, N, hG * Q_MAX, h_port, "u", q_max=d["q_far"],
                       h_max=d["h_far_max"], family=True)
        track_scaled(gm, min(hC * Q_MAX, h_side), L_up, N, "c")
        return True
    out["U2"] = out["D2"] = first_ok(up_ok_housing if config == "housing"
                                     else up_ok)
    out["ZA"] = first_ok(lambda N: bool(two_sided(
        g["H"], N, d["h_floor"], d["h1_cap"], "za", h_max=d["dz_max"])),
        lo=10)
    out["ZB"] = first_ok(lambda N: bool(two_sided(
        g["z_top"] - g["H"], N, d["h1_cap"], d["h_lid"], "zb")))
    if config != "housing":
        return out
    out["ZC"] = first_ok(lambda N: bool(two_sided(
        g["h_ch"], N, d["h_ch_wall"], d["h_floor"], "zc",
        h_max=d["h_ch_max"], family=True)))
    Lzb = g["z_top"] - g["H"]
    lid = (Lzb / d["n"]["ZB"] if d.get("zlever") else two_sided(
        Lzb, d["n"]["ZB"], d["h1_cap"], d["h_lid"], "zb").cells(
            d["n"]["ZB"], Lzb)[-1])
    out["n_port"] = {}
    leads = port_leads(g, d["n"], d)
    for key in d["n_port"]:
        L = port_geometry(g, float(key))["L"]
        out["n_port"][key] = dict(
            XL=first_ok(lambda N: bool(grow_to_plateau(
                L, N, h_port, d["h_run"], "x", family=True,
                lead=leads["XL"], fam_tol=FAM_TOL)), hi=2000),
            ZR=first_ok(lambda N: bool(grow_to_plateau(
                L, N, lid, d["h_run"], "r", family=True,
                lead=leads["ZR"], fam_tol=FAM_TOL)), hi=2000))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("level", nargs="?", choices=list(LEVELS))
    ap.add_argument("-o", "--output", type=Path,
                    help="blockMeshDict path (default: stdout summary only)")
    ap.add_argument("--report", type=Path, help="write design report JSON")
    ap.add_argument("--manifest", type=Path,
                    help="make_geometry.py MANIFEST.sha256 to cross-check "
                         "GEOM against (hard failure on disagreement)")
    ap.add_argument("--nu", type=float,
                    help="coolant kinematic viscosity [m2/s] for the Stokes "
                         "counts (default mu/rho from properties.yaml)")
    ap.add_argument("--table", action="store_true",
                    help="print the count/grading table for all levels")
    ap.add_argument("--minimize", action="store_true",
                    help="replace the (*) medium counts by the smallest "
                         "that meet the design targets (exploration)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL",
                    help="override a DESIGN entry for exploration (e.g. "
                         "n.TH=9, dz_max=3); never used for claim-bearing "
                         "meshes, which take the committed DESIGN")
    ap.add_argument("--config", choices=CONFIGS, default="duct",
                    help="duct (2.5a open duct, default) or housing (D-G "
                         "closed housing with inlet channel and riser)")
    ap.add_argument("--port-mult", type=float, default=1.0,
                    help="housing port length in units of the D-G start "
                         "length L = 10 D_h (1, 2 or 4: the port-gate "
                         "variants L/2L/4L)")
    ap.add_argument("--zlever", choices=T.ZLEVERS, default=None,
                    help="housing z-lever: 'far' (the committed design: the "
                         "0.05 mm cap-plane layer only on |x| <= 30 mm) or "
                         "'module' (the OPTION of K0-design §10.14: also "
                         "coarser cap-plane cells on the module's tile-edge "
                         "and strip-wall lines; off by default)")
    a = ap.parse_args(argv)
    if a.config != "housing" and a.port_mult != 1.0:
        ap.error("--port-mult needs --config housing")
    if a.config != "housing" and a.zlever:
        ap.error("--zlever needs --config housing")
    design = json.loads(json.dumps(DESIGN_HOUSING if a.config == "housing"
                                   else DESIGN))
    if a.zlever:
        design["zlever_module"] = a.zlever == "module"
    for kv in a.set:
        k, v = kv.split("=", 1)
        tgt = design
        *path, last = k.split(".")
        for p in path:
            tgt = tgt[p]
        if last not in tgt:
            ap.error(f"unknown design key {k}")
        tgt[last] = type(tgt[last])(float(v)) if isinstance(tgt[last], int) \
            else float(v)
    geom = load_geometry()
    if a.manifest:
        check_manifest(geom, a.manifest)
    if a.minimize:
        if a.config == "housing" and f"{a.port_mult:g}" not in \
                design["n_port"]:
            design["n_port"][f"{a.port_mult:g}"] = dict(XL=0, ZR=0)
        mc = min_counts(geom, design, config=a.config)
        print("minimal medium counts:", mc)
        if "n_port" in mc:
            design["n_port"].update(mc.pop("n_port"))
        design["n"].update(mc)
    topo = Topology(geom, design, a.config, a.port_mult)
    if a.table:
        txt, reps = table(topo, a.nu)
        print(txt)
        if a.report:
            a.report.write_text(json.dumps(reps, indent=1))
        return 0
    if not a.level:
        ap.error("LEVEL required (or --table)")
    text, cells = render(topo, a.level)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text)
    rep = design_report(topo, a.level, a.nu)
    if a.report:
        a.report.write_text(json.dumps(rep, indent=1))
    print(f"K0 {a.level}: " + " ".join(f"{c}={v}" for c, v in
                                        rep["counts"].items()))
    print(f"   cells: " + ", ".join(f"{k} {v}" for k, v in cells.items())
          + f"; total {sum(cells.values())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
