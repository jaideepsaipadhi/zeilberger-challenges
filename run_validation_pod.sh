#!/bin/bash
# Heavy-piece validation for the zeilberger-challenges repo — run on the pod.
# Usage: bash run_validation_pod.sh   (from the repo root; ~1-3h total, CPU only)
# Sends nothing anywhere: writes validation_pod_report.txt for you to paste back.
set -u
R=validation_pod_report.txt
echo "=== zeilberger-challenges pod validation $(date -u) ===" > $R

run() { # run <label> <timeout-sec> <workdir> <cmd...>
  local label=$1 t=$2 wd=$3; shift 3
  echo "--- $label ---" | tee -a $R
  ( cd "$wd" && timeout "$t" "$@" ) >> $R 2>&1
  echo "[exit $?]" >> $R
}

pip install --quiet mpmath sympy numpy 2>/dev/null || pip install --quiet --break-system-packages mpmath sympy numpy

# 1. KZ constants: full recompute + verification lemmas (the load-bearing numerics)
run kz_constants        7200 kz-constants python3 kz_constants.py
run c1_recompute        7200 kz-constants python3 c1_recompute.py
run c1_tail             3600 kz-constants python3 c1_tail.py
for f in kz-constants/verification/*.py; do
  run "verif/$(basename $f)" 3600 kz-constants/verification python3 "$(basename $f)"
done

# 2. Cone exponent enclosures (interval-arithmetic heavy)
for f in cone-exponent/*.py; do
  run "cone/$(basename $f)" 7200 cone-exponent python3 "$(basename $f)"
done

# 3. Restricted permutations: matsuo reproductions (long)
run matsuo_b            7200 restricted-permutations python3 matsuo_b.py
run matsuo_general      7200 restricted-permutations python3 matsuo_general.py
echo "--- diff regenerated matsuo results vs committed ---" | tee -a $R
for j in matsuo_b_results.json matsuo_general_results.json; do
  if diff <(python3 -m json.tool restricted-permutations/$j 2>/dev/null) \
          <(python3 -m json.tool restricted-permutations/results/$j) >/dev/null 2>&1; then
    echo "$j: MATCHES committed" >> $R
  else
    echo "$j: DIFFERS from committed (or not regenerated)" >> $R
  fi
done

# 4. Chomp: direct re-verification of the 10x42 triple (needs ~big RAM; ~hours)
#    Both engines independently. Comment out if pod < 64GB.
( cd chomp && gcc -O3 -fopenmp -o chomp3 chomp3.c && gcc -O3 -fopenmp -o chomp2 chomp2.c ) >> $R 2>&1
run chomp3_10x42        86400 chomp ./chomp3 10 42 10 42
run chomp2_10x42        86400 chomp ./chomp2 10 42

echo "=== done $(date -u) ===" >> $R
echo "Report written to $R — paste it back to Claude."
