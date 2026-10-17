#!/usr/bin/env python3
"""Compare a fresh reproduction with the published data of this package.

    python3 scripts/check_reproduction.py RESULTS_DIR [DATA_DIR]

RESULTS_DIR is the AISEE_OUT directory of a run of scripts/post/
aisee_results.py (scripts/reproduce.sh writes work/results); DATA_DIR is the
package's data/ (default). Checks, all hard (exit 1 on any failure):

1. RESULTS_DIR/aisee_results.csv against DATA_DIR/aisee_results.csv: the same
   cases in the same order, the same columns, every cell equal (numbers
   compared as floats, the rest as text); the largest absolute and relative
   difference of every numeric column is printed.
2. The scorer's gate record of every case (RESULTS_DIR/aisee_cache/<case>.json)
   against DATA_DIR/gates/<key>.json: gate names, verdicts and details, and the
   info block, after the same normalisation the published files carry
   (gate_record below: the solver-log stamp dropped, the case path reduced to
   the case directory name, non-finite numbers as null).
3. RESULTS_DIR/timeseries/<key>.csv.gz (scripts/make_timeseries.py --out)
   against DATA_DIR/timeseries/<key>.csv.gz: identical decompressed bytes.
4. The three figures of the paper exist and are non-empty PDF files:
   RESULTS_DIR/figures/fig_waveforms.pdf, fig_pareto.pdf, fig_histories.pdf.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIGS = ("fig_waveforms.pdf", "fig_pareto.pdf", "fig_histories.pdf")


def _finite(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_finite(v) for v in o]
    return o


def gate_record(key: str, case_dir: str, cache: dict) -> dict:
    """The published form of one scorer gate record (aisee_results.score)."""
    info = json.loads(json.dumps(cache["info"]))
    pre = info.get("precirc")
    if isinstance(pre, dict) and isinstance(pre.get("case"), str):
        pre["case"] = Path(pre["case"]).name
    return _finite({"key": key, "case_dir": case_dir,
                    "gates": [{"name": g["name"], "ok": bool(g["ok"]), "detail": g["detail"]}
                              for g in cache["gates"]],
                    "info": info})


def _canon(o) -> str:
    return json.dumps(o, sort_keys=True, allow_nan=False)


def check_csv(new: Path, ref: Path, fails: list):
    rn = list(csv.DictReader(new.open()))
    rr = list(csv.DictReader(ref.open()))
    cols_n = list(rn[0]) if rn else []
    cols_r = list(rr[0]) if rr else []
    if cols_n != cols_r:
        fails.append(f"CSV columns differ: {cols_n} vs {cols_r}")
        return
    kn, kr = [r["key"] for r in rn], [r["key"] for r in rr]
    if kn != kr:
        fails.append(f"CSV cases differ: {kn} vs {kr}")
        return
    ndiff, ncells = 0, 0
    print(f"{'column':<20} {'max |diff|':>12} {'max rel':>10}")
    for c in cols_r:
        dmax, rmax, numeric = 0.0, 0.0, True
        for a, b in zip(rn, rr):
            x, y = a[c], b[c]
            ncells += 1
            try:
                fx, fy = float(x), float(y)
            except ValueError:
                numeric = False
                if x != y:
                    ndiff += 1
                    fails.append(f"CSV {a['key']} {c}: {x!r} vs {y!r}")
                continue
            if math.isnan(fx) and math.isnan(fy):
                continue
            d = abs(fx - fy)
            dmax = max(dmax, d)
            rmax = max(rmax, d / abs(fy) if fy else (0.0 if d == 0 else math.inf))
            if fx != fy:
                ndiff += 1
                fails.append(f"CSV {a['key']} {c}: {x!r} vs {y!r}")
        if numeric:
            print(f"{c:<20} {dmax:>12.3g} {rmax:>10.3g}")
        else:
            print(f"{c:<20} {'(text)':>12} {'':>10}")
    print(f"CSV: {len(rr)} cases x {len(cols_r)} columns, {ncells} cells, {ndiff} different")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        sys.exit(__doc__)
    res = Path(argv[0])
    data = Path(argv[1]) if len(argv) > 1 else REPO / "data"
    fails = []
    check_csv(res / "aisee_results.csv", data / "aisee_results.csv", fails)

    sys.path.insert(0, str(REPO / "scripts" / "post"))
    import aisee_results as AR
    keys = {c["key"]: c for c in AR.cases()}
    ng = nt = 0
    for gf in sorted((data / "gates").glob("*.json")):
        pub = json.loads(gf.read_text())
        key, case_dir = pub["key"], pub["case_dir"]
        cf = res / "aisee_cache" / f"{case_dir}.json"
        if not cf.exists():
            fails.append(f"gates {key}: no fresh record {cf}")
            continue
        fresh = gate_record(key, case_dir, json.loads(cf.read_text()))
        if _canon(fresh) != _canon(pub):
            fails.append(f"gates {key}: fresh record differs from {gf.name}")
        else:
            ng += 1
        tf, tr = res / "timeseries" / f"{key}.csv.gz", data / "timeseries" / f"{key}.csv.gz"
        if not tf.exists():
            fails.append(f"timeseries {key}: {tf} missing")
        elif gzip.decompress(tf.read_bytes()) != gzip.decompress(tr.read_bytes()):
            fails.append(f"timeseries {key}: differs from {tr.name}")
        else:
            nt += 1
        if key not in keys:
            fails.append(f"{key}: not a case of aisee_results.py")
    print(f"gate records identical: {ng}; time series identical: {nt}")
    for f in FIGS:
        p = res / "figures" / f
        ok = p.is_file() and p.stat().st_size > 0 and p.read_bytes()[:5] == b"%PDF-"
        print(f"figure {f}: {'ok' if ok else 'MISSING'}" + (f" ({p.stat().st_size} bytes)" if ok else ""))
        if not ok:
            fails.append(f"figure {f} missing or not a PDF")
    if fails:
        print("\nREPRODUCTION CHECK FAILED:")
        for f in fails:
            print("  " + f)
        return 1
    print("\nREPRODUCTION CHECK PASSED: table, gate records, time series and figures match.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
