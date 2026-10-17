#!/usr/bin/env python3
"""INIT_FROM: start a housing case's heated stage from a shared
pre-circulated state (plan NEW-14, NEW-25; WP 2.5c).

A pre-circulated state is the reconstructed t_pre time directory of a
RUN_STAGE=precirc housing case (one per waveform, Qbar, grid, dt, scheme,
coupling, port length, properties). A member case (RUN_STAGE=heat, any heat
model and seed) starts from it by the RESTART path: its own Allrun.pre
builds the mesh, then this script

  1. checks the state key (housing_physics.state_key: the source's caseMeta
     must equal this case's on every key that shapes the isothermal flow,
     and the source must be a precirc case);
  2. checks mesh identity: the sha256 of every region's polyMesh files
     (points faces owner neighbour boundary) must equal the source's - the
     state is copied, never mapped;
  3. copies the source's latest time directory, which must be t_pre, into
     this case (all regions and uniform/, so the time index - and with it
     the phase lock of every step-cadenced monitor - continues);
  4. writes init_from.json: source path, caseMeta sha256, the copied time
     and a sha256 per copied file (the state hash of NEW-14's provenance).

Then Allrun.pre runs ./Allrun.heat (tilt seed, heat-on render).

Usage: init_from.py CASE SOURCE      (SOURCE relative paths: to the repo)
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import housing_physics as HP  # noqa: E402

MESH_FILES = ("points", "faces", "owner", "neighbour", "boundary")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def regions(case: Path):
    return [d.name for d in sorted((case / "constant").iterdir())
            if (d / "polyMesh").is_dir()]


def mesh_hashes(case: Path) -> dict:
    out = {}
    for r in regions(case):
        for f in MESH_FILES:
            p = case / "constant" / r / "polyMesh" / f
            if p.exists():
                out[f"{r}/{f}"] = sha256(p)
    return out


def resolve_source(src: str) -> Path:
    p = Path(src)
    return p if p.is_absolute() else HP.REPO / p


def init_from(case: Path, source: Path) -> dict:
    case, source = case.resolve(), source.resolve()
    if not (source / "system/caseMeta.yaml").exists():
        raise SystemExit(f"init_from: {source} has no system/caseMeta.yaml")
    mine, smeta = HP.load_meta(case), HP.load_meta(source)
    if smeta.get("run_stage") != "precirc":
        raise SystemExit(f"init_from: {source} is not a pre-circulation state "
                         f"(run_stage {smeta.get('run_stage')!r})")
    k_mine, k_src = HP.state_key(mine), HP.state_key(smeta)
    diff = {k: (k_src[k], k_mine[k]) for k in k_mine if k_src[k] != k_mine[k]}
    if diff:
        raise SystemExit("init_from: state key mismatch (source, this case): "
                         + json.dumps(diff, default=str))
    hm, hs = mesh_hashes(case), mesh_hashes(source)
    if not hm or hm != hs:
        bad = sorted(k for k in set(hm) | set(hs) if hm.get(k) != hs.get(k))
        raise SystemExit(f"init_from: mesh differs from the source's "
                         f"({len(bad)} file(s), e.g. {bad[:3]}); a state is "
                         "copied, never mapped")
    name, t = HP.latest_time(source)
    if name is None or abs(t - smeta["t_pre"]) > 0.5 * smeta["lattice"]["dt"]:
        raise SystemExit(f"init_from: the source's latest time {name} is not "
                         f"t_pre = {smeta['t_pre']:g} s (reconstructed state "
                         "missing?)")
    if (case / name).exists():
        raise SystemExit(f"init_from: {case / name} exists already")
    shutil.copytree(source / name, case / name)
    files = {str(p.relative_to(case / name)): sha256(p)
             for p in sorted((case / name).rglob("*")) if p.is_file()}
    state = hashlib.sha256(json.dumps(files, sort_keys=True).encode()
                           ).hexdigest()
    try:
        rel = str(source.relative_to(HP.REPO.resolve()))
    except ValueError:
        rel = None
    rec = dict(source=str(source), source_rel=rel,
               source_caseMeta_sha256=sha256(
        source / "system/caseMeta.yaml"), time=name, t_pre=smeta["t_pre"],
        state_sha256=state, files=files, mesh_sha256=hs,
        state_key=k_src)
    (case / "init_from.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print(__doc__)
        return 2
    rec = init_from(Path(argv[0]), resolve_source(argv[1]))
    print(f"init_from: copied {rec['source']}/{rec['time']} "
          f"({len(rec['files'])} files, state sha256 {rec['state_sha256'][:16]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
