#!/usr/bin/env python3
"""Parametric geometry for the 6-cell 18650 immersion module.

Generates one closed ASCII STL per cell (cell0.stl ... cellN.stl) into
geometry/stl/, plus MANIFEST.sha256 (the STLs are gitignored; the manifest is
committed and is the reproducibility record).

The duct / validation box is NOT an STL: it is the blockMesh background box,
whose dimensions the case generator reads from data/reference/properties.yaml
(validation) or the case matrix (production duct).

Cell dimensions and layout come from properties.yaml + CLI overrides. The
default layout is a 2 x 3 grid at pitch = diameter + gap, cells axis-vertical
(z), matching the module of Liu et al. (ATE 233 (2023) 121184).
TODO(validation phase): confirm Liu's exact arrangement (2x3 vs 1x6 vs 3x2)
and in-box position against the full text / figures before the validation run.

Usage:
    python3 geometry/make_geometry.py                  # defaults from YAML
    python3 geometry/make_geometry.py --rows 1 --cols 6 --gap 0.002
    python3 geometry/make_geometry.py --segments 128   # finer tessellation
"""

import argparse
import hashlib
import math
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("pyyaml required: pip install -r requirements.txt")

REPO = Path(__file__).resolve().parent.parent
STL_DIR = REPO / "geometry" / "stl"
PROPS = REPO / "data" / "reference" / "properties.yaml"


def cylinder_facets(cx, cy, z0, z1, r, n):
    """Closed cylinder as triangle facets: side + two end-cap fans."""
    facets = []  # (normal, v1, v2, v3)
    ring = [(cx + r * math.cos(2 * math.pi * i / n),
             cy + r * math.sin(2 * math.pi * i / n)) for i in range(n)]
    for i in range(n):
        (xa, ya), (xb, yb) = ring[i], ring[(i + 1) % n]
        # outward normal of the side panel (mean of the two vertex normals)
        nx, ny = (xa + xb) / 2 - cx, (ya + yb) / 2 - cy
        m = math.hypot(nx, ny)
        nrm = (nx / m, ny / m, 0.0)
        # two triangles per side panel, outward-facing (CCW seen from outside)
        facets.append((nrm, (xa, ya, z0), (xb, yb, z0), (xb, yb, z1)))
        facets.append((nrm, (xa, ya, z0), (xb, yb, z1), (xa, ya, z1)))
        # bottom cap (normal -z): CCW seen from below
        facets.append(((0, 0, -1), (cx, cy, z0), (xb, yb, z0), (xa, ya, z0)))
        # top cap (normal +z)
        facets.append(((0, 0, 1), (cx, cy, z1), (xa, ya, z1), (xb, yb, z1)))
    return facets


def write_ascii_stl(path, name, facets):
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for nrm, v1, v2, v3 in facets:
            f.write(f"  facet normal {nrm[0]:.9e} {nrm[1]:.9e} {nrm[2]:.9e}\n")
            f.write("    outer loop\n")
            for v in (v1, v2, v3):
                f.write(f"      vertex {v[0]:.9e} {v[1]:.9e} {v[2]:.9e}\n")
            f.write("    endloop\n  endfacet\n")
        f.write(f"endsolid {name}\n")


def main():
    props = yaml.safe_load(PROPS.read_text())
    cell = props["cell"]

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rows", type=int, default=2)
    p.add_argument("--cols", type=int, default=3)
    p.add_argument("--gap", type=float, default=0.002,
                   help="inter-cell surface gap [m] (D-02: 2 mm)")
    p.add_argument("--diameter", type=float, default=cell["diameter"])
    p.add_argument("--height", type=float, default=cell["height"])
    p.add_argument("--z0", type=float, default=0.0,
                   help="cell base height [m] above origin")
    p.add_argument("--segments", type=int, default=96,
                   help="tessellation segments around the circumference")
    p.add_argument("--out", type=Path, default=STL_DIR,
                   help="output directory (default geometry/stl; cases pass "
                        "their own constant/triSurface so concurrent case "
                        "builds never race on shared files)")
    a = p.parse_args()

    r = a.diameter / 2
    pitch = a.diameter + a.gap
    # centre the grid on (0, 0)
    x0 = -(a.cols - 1) * pitch / 2
    y0 = -(a.rows - 1) * pitch / 2

    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    idx = 0
    for row in range(a.rows):
        for col in range(a.cols):
            cx, cy = x0 + col * pitch, y0 + row * pitch
            name = f"cell{idx}"
            path = out / f"{name}.stl"
            write_ascii_stl(path, name,
                            cylinder_facets(cx, cy, a.z0, a.z0 + a.height,
                                            r, a.segments))
            h = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest.append(f"{h}  {name}.stl")
            print(f"  {name}.stl  centre=({cx*1000:+.1f},{cy*1000:+.1f}) mm"
                  f"  sha256={h[:12]}...")
            idx += 1

    params = (f"# rows={a.rows} cols={a.cols} gap={a.gap} d={a.diameter} "
              f"h={a.height} z0={a.z0} segments={a.segments}")
    (out / "MANIFEST.sha256").write_text(
        params + "\n" + "\n".join(manifest) + "\n")
    print(f"{idx} cells, pitch {pitch*1000:.1f} mm; manifest written.")


if __name__ == "__main__":
    main()
