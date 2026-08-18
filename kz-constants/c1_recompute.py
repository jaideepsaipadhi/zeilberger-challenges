#!/usr/bin/env python3
"""
c1_recompute.py — settle the value of C_1 by computation.

THE QUESTION. Two values are in circulation:

    (A)  C_1 = 0.5212860373...        derived 10 Aug: exact DP for G(n) to
                                      n = 110, Richardson extrapolation.
                                      Recorded in the session transcript
                                      together with its method.
    (B)  C_1 = 0.52128605909203      present in later notes with NO recorded
                                      derivation, and used as the basis of a
                                      claim that the published Kauers-Zeilberger
                                      value is wrong in its eighth digit.

These DISAGREE at the eighth significant digit (0.521286_037_ vs
0.521286_059_). A longer extrapolation sharpens a value; it does not move it
by 2.2e-8. So at most one is right, and (A) is the one with a method attached.

This script recomputes from scratch and reports which the data supports.

THE OBJECT. G(n) = number of standard Young tableaux of shape [n,n,n] in which
every maximal run in every row has length >= 2. Equivalently (lattice-word
bijection) the number of words on {0,1,2} with n letters each such that every
prefix has #0 >= #1 >= #2 (the ballot condition) and every maximal run of
equal letters has length >= 2.

Conjecture 2a: G(n) ~ C_1 * 8^n / n^4.

THE DP. State: (c1, c2, d, r) at total length s, where c0 = s - c1 - c2 is the
count of letter 0, d in {0,1,2} is the last letter, and r in {1,2} records
whether the current run has length 1 or has already reached 2. Transitions
append a letter; a change of letter is legal only if r = 2. The ballot
condition is enforced as an invariant on (c0,c1,c2). All arithmetic is exact
Python integer arithmetic.

VALIDATION (runs first, and the script aborts if it fails). The published
initial terms of G are checked. If the DP does not reproduce them, nothing
below is meaningful.

EXTRAPOLATION. Set A_N = G(N) * N^4 / 8^N, so A_N -> C_1. Assuming the
standard expansion A_N = C_1 (1 + a_1/N + a_2/N^2 + ...), iterated Richardson
    A^{(k+1)}_N = ( (N+1) A^{(k)}_{N+1} - (N-k) A^{(k)}_N ) / (k+1)
accelerates convergence. We report the table of successive orders so the
stabilised digits are visible rather than asserted, and we print the distance
from each candidate value.

Run:  python3 c1_recompute.py --nmax 120
      python3 c1_recompute.py --nmax 200 --dps 40      (slower)
"""
import argparse, json, time
from mpmath import mp, mpf, nstr

ap = argparse.ArgumentParser()
ap.add_argument('--nmax', type=int, default=120)
ap.add_argument('--dps', type=int, default=40)
ap.add_argument('--orders', type=int, default=8)
args = ap.parse_args()
mp.dps = args.dps
T0 = time.time()

# published initial terms of G(n) (runs >= 2), from Kauers-Zeilberger
KNOWN = [1, 0, 1, 1, 5, 15, 69, 304, 1518, 7807, 42314, 236621,
         1364570, 8062975, 48680547, 299388670, 1871463427]

def G_series(nmax):
    """Exact G(n) for n = 0..nmax by layered DP.

    dp[(c1, c2, d, r)] = number of ballot words of length s with the given
    letter counts, last letter d, and run-state r (1 = current run has length
    1, 2 = current run already has length >= 2).  c0 = s - c1 - c2.
    """
    out = [0]*(nmax+1)
    out[0] = 1
    # start: first letter must be 0 (ballot), run length 1
    dp = {(0, 0, 0, 1): 1}
    s = 1
    total = 3*nmax
    # record G(n) when s == 3n and state is complete
    while s < total:
        ndp = {}
        for (c1, c2, d, r), v in dp.items():
            c0 = s - c1 - c2
            for e in (0, 1, 2):
                if e != d and r == 1:
                    continue            # cannot leave a run of length 1
            
                n0, n1, n2 = c0, c1, c2
                if e == 0:   n0 += 1
                elif e == 1: n1 += 1
                else:        n2 += 1
                if n0 > nmax or n1 > nmax or n2 > nmax:
                    continue
                if not (n0 >= n1 >= n2):
                    continue
                nr = 2 if e == d else 1
                key = (n1, n2, e, nr)
                ndp[key] = ndp.get(key, 0) + v
        dp = ndp
        s += 1
        if s % 3 == 0:
            n = s // 3
            if n <= nmax:
                tot = 0
                for (c1, c2, d, r), v in dp.items():
                    c0 = s - c1 - c2
                    if c0 == n and c1 == n and c2 == n and r == 2:
                        tot += v
                out[n] = tot
                if n <= 10 or n % 10 == 0:
                    print("   G(%d) = %s   [%.0fs]"
                          % (n, tot if n <= 8 else "%d digits" % len(str(tot)),
                             time.time()-T0), flush=True)
    return out

print("computing G(n) exactly to n = %d ..." % args.nmax, flush=True)
G = G_series(args.nmax)

print("\nvalidation against published terms:")
ok = True
for n, want in enumerate(KNOWN):
    if n > args.nmax: break
    got = G[n]
    flag = "ok" if got == want else "MISMATCH"
    if got != want: ok = False
    print("   n=%2d  computed %-12d published %-12d  %s" % (n, got, want, flag))
if not ok:
    raise SystemExit("\nDP does not reproduce the published terms — aborting; "
                     "nothing downstream would be meaningful.")
print("   all published terms reproduced.")

# ---------------------------------------------------------------- extrapolate
print("\nRichardson extrapolation of A_N = G(N)*N^4/8^N ...", flush=True)
NS = [n for n in range(args.nmax//2, args.nmax+1)]
A = {n: mpf(G[n]) * mpf(n)**4 / mpf(8)**n for n in NS}

cols = [[A[n] for n in NS]]
for k in range(args.orders):
    prev = cols[-1]
    nxt = []
    for i in range(len(prev)-1):
        N = NS[i]
        nxt.append(((N+1)*prev[i+1] - (N-k)*prev[i]) / (k+1))
    cols.append(nxt)

print("\n   order   last value")
for k, col in enumerate(cols):
    if col:
        print("   %5d   %s" % (k, nstr(col[-1], 20)))

best = cols[-1][-1] if cols[-1] else cols[0][-1]
# stability: difference between the last two entries of the deepest column
stab = None
if len(cols[-1]) >= 2:
    stab = abs(cols[-1][-1] - cols[-1][-2])

CAND_A = mpf('0.5212860373')
CAND_B = mpf('0.52128605909203')

print("\n--- comparison ---")
print("   extrapolated C_1        =", nstr(best, 20))
if stab is not None:
    print("   stability (last diff)   =", nstr(stab, 6))
print("   candidate A 0.5212860373              distance =",
      nstr(abs(best-CAND_A), 6))
print("   candidate B 0.52128605909203    distance =",
      nstr(abs(best-CAND_B), 6))

verdict = ("A" if abs(best-CAND_A) < abs(best-CAND_B) else "B")
print("\n   VERDICT: the data supports candidate %s." % verdict)
print("   (Compare each distance against the stability figure: a candidate")
print("    further away than the stability estimate is excluded.)")

json.dump({
  "nmax": args.nmax,
  "G_first": [str(G[n]) for n in range(min(11, args.nmax+1))],
  "validation_passed": ok,
  "extrapolated": nstr(best, 25),
  "stability": nstr(stab, 8) if stab is not None else None,
  "candidate_A": "0.5212860373", "dist_A": nstr(abs(best-CAND_A), 8),
  "candidate_B": "0.52128605909203", "dist_B": nstr(abs(best-CAND_B), 8),
  "verdict": verdict,
  "orders": [nstr(c[-1], 20) for c in cols if c],
}, open("c1_recompute_results.json", "w"), indent=1)
print("\nwrote c1_recompute_results.json  [%.0fs]" % (time.time()-T0))
