#!/bin/bash
# Chomp Ch.2 hunt v4: triple search + TB track dumps for the rotation-number program.
set -e
export OMP_NUM_THREADS=13
LOG=chomp_results.log
gcc -O3 -fopenmp -o chomp3 chomp3.c
gcc -O2 -o chomp chomp.c
python3 cmp.py | tail -1 | tee -a $LOG
d1=$(./chomp 12 12 | grep P-pos); d3=$(./chomp3 12 12 | grep P-pos)
[ "$d1" = "$d3" ] && echo "v1/v3 agree: $d1" | tee -a $LOG || { echo MISMATCH-STOP; exit 1; }
run(){ echo "=== BOX $1 x $2 (RAM $3, ~$4 @13c) ===" | tee -a $LOG
       ./chomp3 $1 $2 | grep -E "BAR .* [2-9]+ winning|BOX|P-pos" | tee -a $LOG; }
runtb(){ echo "=== BOX $1 x $2 TB (RAM $3, ~$4) ===" | tee -a $LOG
       ./chomp3 $1 $2 TB > tb_$1x$2.out
       grep -E "BAR .* [2-9]+ winning|BOX|P-pos" tb_$1x$2.out | tee -a $LOG; }
# --- track-data boxes (TB dumps feed fit_tracks.py) ---
runtb 4 300  "44MB"  "2min"     # finish the 4-row tracks
runtb 5 150  "90MB"  "3min"
runtb 6 80   "60MB"  "2min"
runtb 7 110  "1.4GB" "4min"     # the a=7 rotation numbers
runtb 8 60   "1.1GB" "3min"
runtb 9 55   "3.4GB" "10min"
runtb 10 45  "3.7GB" "12min"
# --- pure hunt boxes ---
run 17 17  "0.3GB" "1min"
run 18 18  "1.2GB" "5min"
run 19 19  "4.5GB" "20min"
runtb 8 90   "23GB"  "45min"    # RAM-gated from here down: check free -g
runtb 9 70   "13GB"  "35min"
runtb 10 55  "22GB"  "50min"
runtb 7 150  "40GB"  "1h"
run 20 20  "17.3GB" "1.5h"
echo "=== SUMMARY ===" | tee -a $LOG
grep -h "winning moves" $LOG | sort -u | tee -a $LOG
echo "Now: python3 fit_tracks.py tb_*.out   -> rotation numbers + triple-candidate widths"
