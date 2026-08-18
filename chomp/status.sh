#!/bin/bash
echo "=== PROCESS ==="
ps -p 932 -o pid,etime,%cpu 2>/dev/null || echo "PID 932 done/dead -- check if queue finished below"
echo; echo "=== QUEUE POSITION (last 3 boxes) ==="
grep "=== BOX" console.log | tail -3
tail -2 console.log
echo; echo "=== ALL FINDINGS SO FAR ==="
grep -h "winning moves" chomp_results.log 2>/dev/null | sort -u
T=$(grep -c "3 winning moves\|[3-9] winning moves" chomp_results.log 2>/dev/null)
[ "$T" != "0" ] && echo "*** CHECK ABOVE: possible 3+ line! ***"
echo; echo "=== TRACK DUMPS AVAILABLE ==="
ls -la tb_*.out 2>/dev/null | awk '{print $NF, $5}'
echo; echo "=== DONE? ==="
grep -q "=== SUMMARY ===" console.log && echo "MAIN QUEUE COMPLETE -> run: python3 fit_tracks.py | tee fits.txt" || echo "still running"
