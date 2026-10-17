#!/usr/bin/env python3
"""Write a setExprFieldsDict that initialises p hydrostatically (NEW-14 prep).

Why (K0-design.md §5.1 b): the template starts p = 101325 Pa uniform.
chtMultiRegionFoam then forces p_rgh = p - rho*gh (createFluidFields.H,
"Force p_rgh to be consistent with p"), which leaves a hydrostatic mismatch
of up to rho*g*z (about 1 kPa in the duct, 3.3 kPa up the 164 mm riser)
against the outlet's fixedValue p_rgh. The first pressure solves remove it
as a spurious start-up flow (mm/s on the K0 mesh). With

    p = p_rgh,out + rho0 * ((g & C) - ghRef)

the start-up p_rgh equals the outlet value everywhere, i.e. the rest state
is in hydrostatic balance. rho0 is the fluid's own equation of state at the
initial temperature, so the balance is exact for the uniform initial T.

Inputs, all read from the case (hard failure if any is missing or not of
the supported form; never a silent default):
  constant/<region>/thermophysicalProperties  icoPolynomial rhoCoeffs<8>
  0/<region>/T                                internalField uniform T0
  constant/g                                  value (gx gy gz)
  constant/<region>/hRef                      optional (else 0), as the
                                              solver reads it
  0/<region>/p_rgh                            <patch> fixedValue uniform p0,
                                              or (D-G.1, OUTLET_P=p)
                                              prghPressure p uniform P: then
                                              p0 = P - rho0 (g & C_out -
                                              ghRef), C_out at the outlet
                                              height of system/caseMeta.yaml
                                              (ports.z_outlet_m), g along z;
                                              a value entry must equal p0
Output: system/<region>/setExprFieldsDict (setExprFields -region <region>
reads it there; syntax as the v2606 tutorial
compressible/rhoPimpleFoam/RAS/TJunctionAverage).

ghRef follows createFluidFields.H: ghRef = g & (cmptMag(g)/mag(g)) * hRef.

Usage: hydrostatic_init.py CASE [--region fluid] [--patch outlet]
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path


def _strip(txt):
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    return re.sub(r"//[^\n]*", "", txt)


def read_rho0(case, region, T0):
    txt = _strip((case / "constant" / region /
                  "thermophysicalProperties").read_text())
    if not re.search(r"equationOfState\s+icoPolynomial\s*;", txt):
        raise SystemExit("hydrostatic_init: equationOfState is not "
                         "icoPolynomial (only that form is supported)")
    m = re.search(r"rhoCoeffs<8>\s*\(([^)]*)\)", txt)
    if not m:
        raise SystemExit("hydrostatic_init: rhoCoeffs<8> not found")
    a = [float(v) for v in m.group(1).split()]
    return sum(c * T0 ** i for i, c in enumerate(a))


def read_uniform_internal(path):
    txt = _strip(Path(path).read_text())
    m = re.search(r"internalField\s+uniform\s+([^;\s]+)\s*;", txt)
    if not m:
        raise SystemExit(f"hydrostatic_init: {path}: internalField is not "
                         f"uniform (hydrostatic balance needs uniform T)")
    return float(m.group(1))


def read_g(case):
    txt = _strip((case / "constant" / "g").read_text())
    m = re.search(r"value\s*\(\s*([^)]*)\)", txt)
    return [float(v) for v in m.group(1).split()]


def read_href(case, region):
    p = case / "constant" / region / "hRef"
    if not p.exists():
        return 0.0
    txt = _strip(p.read_text())
    m = re.search(r"value\s+([^;\s]+)\s*;", txt)
    return float(m.group(1))


def read_patch_entry(path, patch):
    """(type, value, p) of a fixedValue or prghPressure p_rgh patch (p None
    for fixedValue, value None if absent); anything else is refused."""
    txt = _strip(Path(path).read_text())
    bf = txt[txt.index("boundaryField"):]
    m = re.search(rf"\b{re.escape(patch)}\s*\{{([^}}]*)\}}", bf)
    if not m:
        raise SystemExit(f"hydrostatic_init: patch {patch} not in {path}")
    body = m.group(1)
    t = re.search(r"type\s+(\w+)\s*;", body)
    kind = t.group(1) if t else None
    if kind not in ("fixedValue", "prghPressure"):
        raise SystemExit(f"hydrostatic_init: {path} {patch} is not "
                         f"fixedValue or prghPressure (no pressure reference)")
    v = re.search(r"\bvalue\s+uniform\s+([^;\s]+)\s*;", body)
    p = re.search(r"\bp\s+uniform\s+([^;\s]+)\s*;", body)
    if kind == "prghPressure" and not p:
        raise SystemExit(f"hydrostatic_init: {path} {patch}: prghPressure "
                         "without a uniform p")
    if kind == "fixedValue" and not v:
        raise SystemExit(f"hydrostatic_init: {path} {patch}: fixedValue "
                         "without a uniform value")
    return (kind, float(v.group(1)) if v else None,
            float(p.group(1)) if p else None)


def read_patch_value(path, patch):
    """The uniform value of a fixedValue p_rgh patch (refuses other types)."""
    kind, v, _ = read_patch_entry(path, patch)
    if kind != "fixedValue":
        raise SystemExit(f"hydrostatic_init: {path} {patch} is not "
                         f"fixedValue (no pressure reference)")
    return v


def read_patch_height(case):
    """z of the riser-top outlet face [m] from system/caseMeta.yaml
    (housing_physics ports.z_outlet_m); refused without it."""
    import yaml
    f = Path(case) / "system" / "caseMeta.yaml"
    if not f.exists():
        raise SystemExit("hydrostatic_init: a prghPressure outlet needs the "
                         "outlet height (system/caseMeta.yaml ports.z_outlet_m)")
    meta = yaml.safe_load(f.read_text()) or {}
    try:
        return float(meta["ports"]["z_outlet_m"])
    except (KeyError, TypeError) as e:
        raise SystemExit(f"hydrostatic_init: caseMeta has no ports.z_outlet_m ({e!r})")


def build(case, region="fluid", patch="outlet"):
    case = Path(case)
    T0 = read_uniform_internal(case / "0" / region / "T")
    rho0 = read_rho0(case, region, T0)
    g = read_g(case)
    mag = math.sqrt(sum(c * c for c in g))
    href = read_href(case, region)
    gh_ref = sum(c * abs(c) / mag for c in g) * href if mag > 1e-300 else 0.0
    kind, value, p_out = read_patch_entry(case / "0" / region / "p_rgh", patch)
    if kind == "fixedValue":
        p0 = value
    else:
        # D-G.1 (OUTLET_P=p): prghPressure fixes the static p on the face, so
        # the start-up p_rgh there is p - rho0 (g & C_out - ghRef) (v2606
        # prghPressureFvPatchScalarField::updateCoeffs); the face is the
        # horizontal riser top at z_out, g along z
        if abs(g[0]) > 0 or abs(g[1]) > 0:
            raise SystemExit("hydrostatic_init: a prghPressure outlet needs g "
                             "along z (the outlet height alone sets g & C)")
        p0 = p_out - rho0 * (g[2] * read_patch_height(case) - gh_ref)
        if value is not None and abs(value - p0) > 1e-9 * max(1.0, abs(p0)):
            raise SystemExit(f"hydrostatic_init: {patch} prghPressure value "
                             f"{value:.17g} is not its start-up p_rgh {p0:.17g} "
                             "(p - rho0 g.r_out): the render is inconsistent")
    gx, gy, gz = (f"{c:.17g}" for c in g)
    expr = (f"{p0:.17g} + {rho0:.17g}*({gx}*pos().x() + {gy}*pos().y()"
            f" + {gz}*pos().z() - ({gh_ref:.17g}))")
    body = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      setExprFieldsDict;
}}
// Generated by scripts/hydrostatic_init.py - do not edit.
// Hydrostatic initial p (NEW-14 prep, K0-design.md §5.1 b):
//   p = p_rgh,out + rho0*((g & C) - ghRef), so that the solver's start-up
//   p_rgh = p - rho*gh equals the {patch} {kind} start-up p_rgh {p0:.10g} Pa everywhere.
//   rho0 = icoPolynomial(T0 = {T0:.10g} K) = {rho0:.10g} kg/m3,
//   g = ({gx} {gy} {gz}), hRef = {href:.10g} m -> ghRef = {gh_ref:.10g}.

expressions
(
    p
    {{
        field       p;
        dimensions  [1 -1 -2 0 0 0 0];
        expression
        #{{
            {expr}
        #}};
    }}
);
"""
    out = case / "system" / region / "setExprFieldsDict"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(body)
    return dict(T0=T0, rho0=rho0, g=g, hRef=href, ghRef=gh_ref, p0=p0,
                outlet_type=kind, p_outlet=p_out,
                expression=expr, dict=str(out))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", type=Path)
    ap.add_argument("--region", default="fluid")
    ap.add_argument("--patch", default="outlet",
                    help="fixedValue (or prghPressure, D-G.1) p_rgh patch "
                         "that sets the level")
    a = ap.parse_args(argv)
    r = build(a.case, a.region, a.patch)
    print(f"hydrostatic p: rho0 {r['rho0']:.6g} kg/m3 at T0 {r['T0']:.6g} K, "
          f"p_rgh level {r['p0']:.10g} Pa ({a.patch}, {r['outlet_type']}), "
          f"ghRef {r['ghRef']:g}"
          f"\n  p = {r['expression']}\n  -> {r['dict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
