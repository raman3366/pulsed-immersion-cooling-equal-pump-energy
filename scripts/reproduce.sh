#!/usr/bin/env bash
# Reproduce the paper's results table, the per-case numerical checks (the
# scorer's gates) and the paper's figures from the case archives in data/cases/.
#
#   bash scripts/reproduce.sh [WORK_DIR]        (default: work/ in this repository)
#
# 1. extracts every data/cases/*.tar.xz into WORK_DIR/cases (fresh copies);
# 2. scores every case with scripts/post/aisee_results.py (all gates
#    recomputed, --refresh) and writes WORK_DIR/results/aisee_results.csv, the
#    gate records WORK_DIR/results/aisee_cache/<case>.json and the figures
#    WORK_DIR/results/figures/fig_{geometry,waveforms,pareto,matched,histories}.pdf;
# 3. writes the derived time series to WORK_DIR/results/timeseries/;
# 4. prints the derived numbers of the paper (scripts/paper_numbers.py, also
#    saved to WORK_DIR/results/paper_numbers.txt);
# 5. compares the table, gate records and time series with data/ and checks that
#    the five figures exist
#    (scripts/check_reproduction.py) and exits non-zero on any difference.
#
# Six cases fail one or two gates and are flagged in the table (column
# gates_failed), as in the paper; they are named in --allow-fail below. Any
# other gate failure stops the scorer.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$REPO/work}"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
PY="${PYTHON:-python3}"

echo ">>> extracting data/cases/*.tar.xz into $WORK/cases"
mkdir -p "$WORK/cases"
n=0
for t in "$REPO"/data/cases/*.tar.xz; do
    name="$(basename "$t" .tar.xz)"
    rm -rf "${WORK:?}/cases/$name"
    tar -xJf "$t" -C "$WORK/cases"
    [ -f "$WORK/cases/$name/case.env" ] || { echo "extraction of $t failed"; exit 1; }
    n=$((n + 1))
done
echo "    $n cases"

export AISEE_CASE_ROOT="$WORK/cases"
export AISEE_OUT="$WORK/results"
mkdir -p "$AISEE_OUT"

echo ">>> scoring (results table, gates, figures)"
"$PY" "$REPO/scripts/post/aisee_results.py" --refresh --allow-fail \
    C1_Q90 steady_Q360 steady_Q360_dt6.25 steady_Q135 steady_Q360_cand2 steady_Q180

echo ">>> time series"
"$PY" "$REPO/scripts/make_timeseries.py" --out "$AISEE_OUT/timeseries"

echo ">>> derived numbers of the paper (energy ratios, checks, comparisons)"
"$PY" "$REPO/scripts/paper_numbers.py" "$AISEE_OUT" | tee "$AISEE_OUT/paper_numbers.txt"

echo ">>> comparison with data/"
"$PY" "$REPO/scripts/check_reproduction.py" "$AISEE_OUT" "$REPO/data"
