#!/usr/bin/env python3
"""Write the derived temperature series of every case to <out>/<key>.csv.gz.

Reads the extracted case directories (AISEE_CASE_ROOT, as for
scripts/post/aisee_results.py) through the scorer's own readers, so every
number is the one the scorer uses:

  t_s                    simulation time (s); the heat source is on for
                         t > t_pre = 30 s (system/caseMeta.yaml)
  Tmax_cell<i>_C         largest temperature in cell i (volFieldValue max,
                         monitor Tmax_cell<i>), deg C
  TvolAverage_cell<i>_C  volume-mean temperature of cell i (monitor
                         TvolAverage_cell<i>), deg C
  Tmax_module_C          max over the six cells of Tmax_cell<i>, deg C
                         (aisee_results.series)
  dTcell_K               max minus min of the six volume means, K
  Tmodule_mean_C         volume-weighted mean of the six volume means, deg C
                         (weights: the meshed cell volumes of the monitor
                         headers, as in aisee_results.series)

Usage: python3 scripts/make_timeseries.py [--out DIR]   (default data/timeseries)
The files are gzip with a zero timestamp, so the same input gives the same
bytes.
"""
from __future__ import annotations

import argparse
import gzip
import io
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts" / "post"))
import aisee_results as AR      # noqa: E402
import housing_io as HIO        # noqa: E402

K0 = 273.15
FMT = "%.10g"


def table(case: Path):
    t, smax, dtc, mod = AR.series(case)
    cols = {"t_s": t}
    for i in range(AR.N_CELLS):
        ts, v = HIO.column(case, f"cell{i}", f"Tmax_cell{i}")
        if not np.array_equal(ts, t):
            raise SystemExit(f"{case.name}: Tmax_cell{i} not on the series time base")
        cols[f"Tmax_cell{i}_C"] = v - K0
    for i in range(AR.N_CELLS):
        ts, v = HIO.column(case, f"cell{i}", f"TvolAverage_cell{i}")
        if not np.array_equal(ts, t):
            raise SystemExit(f"{case.name}: TvolAverage_cell{i} not on the series time base")
        cols[f"TvolAverage_cell{i}_C"] = v - K0
    cols["Tmax_module_C"] = smax - K0
    cols["dTcell_K"] = dtc
    cols["Tmodule_mean_C"] = mod - K0
    return cols


def write(cols: dict, path: Path):
    names = list(cols)
    data = np.column_stack([cols[n] for n in names])
    buf = io.StringIO()
    buf.write(",".join(names) + "\n")
    np.savetxt(buf, data, delimiter=",", fmt=FMT)
    raw = buf.getvalue().encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        with gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0, compresslevel=9) as g:
            g.write(raw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=REPO / "data" / "timeseries")
    a = ap.parse_args(argv)
    n = 0
    for c in AR.cases():
        if c["path"] is None:
            print(f"{c['key']}: case directory not found under {AR.ARCHIVE}", file=sys.stderr)
            continue
        cols = table(c["path"])
        write(cols, a.out / f"{c['key']}.csv.gz")
        n += 1
        print(f"{c['key']:<24} {cols['t_s'].size} samples, t {cols['t_s'][0]:g}..{cols['t_s'][-1]:g} s")
    print(f"wrote {n} series to {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
