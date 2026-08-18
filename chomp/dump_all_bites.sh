#!/bin/bash
set -e
export OMP_NUM_THREADS=13
mkdir -p data
: > data/bites_all.txt
for bar in "6 13" "6 93" "7 29" "7 30" "7 57" "8 10" "8 22" "8 23" "9 10" "9 26" \
           "10 14" "10 29" "10 33" "10 35" "11 18" "12 13" "13 16" "14 16" "10 42"; do
  a=$(echo $bar|cut -d' ' -f1); b=$(echo $bar|cut -d' ' -f2)
  echo "=== $a x $b ===" | tee -a data/bites_all.txt
  ./chomp3 $a $b $a $b | grep -E "BAR $a x $b |winning bite" | tee -a data/bites_all.txt
done
echo "DONE"
