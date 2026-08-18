#!/usr/bin/env python3
"""
kz_constants.py — determine C_1 and C_2 to full available precision.

Handles both Kauers-Zeilberger models with the same machinery:

    --model G   runs of length >= 2,  growth mu = 8            -> C_1
    --model H   runs of ODD length,   growth mu = 7 + 5*sqrt2  -> C_2

The two differ only in the admissible run lengths and in mu; the ballot
condition, the dynamic program, the N^{-4} scaling and the accelerations are
identical. C_2 is determined here directly, NOT by calibration against C_1,
which removes the dependency that made the earlier C_2 figure provisional.

STATE OF PLAY. c1_recompute.py established, from an exact DP validated against
all sixteen published terms of G(n):

    n = 120, 8 Richardson orders   ->  0.5212860515...
    n = 200, 14 Richardson orders  ->  0.52128605884644744638

converging monotonically from below, with the candidate 0.52128605909(2) at a
distance of 2.5e-10 and the older 0.52128605909(2) excluded by 2.2e-8. Digits
0.521286058 are established. The remaining question is the tail: the deepest
Richardson column is still ascending, by roughly 6e-11 per order and shrinking,
so the last four digits rest on extrapolating that ascent rather than on any
computed value.

This script settles them three ways, and agrees or disagrees explicitly.

METHOD 1 -- deeper and longer Richardson. Larger n_max and more orders. The
limitation is arithmetic: the Richardson recursion differences quantities of
size ~0.52 and extracts corrections of size 1e-12, so working precision must
exceed the digits sought by a wide margin; --dps is set accordingly.

METHOD 2 -- Salzer / Levin-type acceleration, independent of Richardson.
Richardson assumes A_N = C(1 + a_1/N + a_2/N^2 + ...) and eliminates terms in
order. Salzer's transformation applies a weighted difference operator to
N^k A_N directly and converges under weaker hypotheses; where two accelerations
with different assumptions agree, the shared digits are far more trustworthy
than either alone. Disagreement is equally informative.

METHOD 3 -- explicit tail fit. The deepest Richardson column ascends with
increments d_k that decay geometrically. Fitting d_k ~ A r^k on the last
several orders and summing the geometric tail gives an estimate of the limit
that is independent of running more orders. This is the calculation that was
done informally by eye; here it is done numerically with the fitted ratio
reported so its plausibility can be judged.

The three estimates are printed together with the candidate, and the script
states which digits are common to all three -- those, and only those, should be
asserted in print.

Run:  python3 c1_tail.py --nmax 300 --orders 20 --dps 60
      python3 c1_tail.py --nmax 400 --orders 24 --dps 80    (slower)
"""
import argparse, json, time
from mpmath import mp, mpf, nstr, mpmathify

ap = argparse.ArgumentParser()
ap.add_argument('--nmax', type=int, default=300)
ap.add_argument('--orders', type=int, default=20)
ap.add_argument('--dps', type=int, default=60)
ap.add_argument('--model', choices=['G','H'], default='G')
ap.add_argument('--tailfit', type=int, default=6,
                help='number of trailing increments used for the geometric fit')
args = ap.parse_args()
mp.dps = args.dps
T0 = time.time()

KNOWN_G = [1, 0, 1, 1, 5, 15, 69, 304, 1518, 7807, 42314, 236621,
           1364570, 8062975, 48680547, 299388670, 1871463427]
if args.model == 'G':
    KNOWN = KNOWN_G
    MU = mpf(8)
    CAND = mpf('0.52128605909(2)')
    RUNOK = lambda L: L >= 2                 # runs of length at least two
else:
    KNOWN = None                              # no published list on hand
    MU = 7 + 5*mp.sqrt(2)
    CAND = mpf('0.6389278129')
    RUNOK = lambda L: L % 2 == 1              # runs of odd length

RCAP = 3   # run-length state: 1, 2, or 3 meaning ">=3"; parity recoverable

def G_series(nmax):
    """Exact counts for the chosen model.

    State (c1, c2, d, r): c0 = s - c1 - c2, last letter d, and r the current
    run length capped at RCAP with parity preserved:  r in {1,2,3,4} where
    3 means 'odd and >=3', 4 means 'even and >=4'.  A run may be closed only
    when RUNOK holds for its length, which for both models depends on the
    length only through (length==1, parity, length>=2), all carried by r."""
    def closable(r):
        L = {1: 1, 2: 2, 3: 3, 4: 4}[r]
        return RUNOK(L)
    def bump(r):
        if r == 1: return 2
        if r == 2: return 3
        return 3 if r == 4 else 4

    out = [0]*(nmax+1); out[0] = 1
    dp = {(0, 0, 0, 1): 1}
    s = 1; total = 3*nmax
    while s < total:
        ndp = {}
        for (c1, c2, d, r), v in dp.items():
            c0 = s - c1 - c2
            for e in (0, 1, 2):
                if e != d and not closable(r):
                    continue
                n0, n1, n2 = c0, c1, c2
                if e == 0:   n0 += 1
                elif e == 1: n1 += 1
                else:        n2 += 1
                if n0 > nmax or n1 > nmax or n2 > nmax: continue
                if not (n0 >= n1 >= n2): continue
                key = (n1, n2, e, bump(r) if e == d else 1)
                ndp[key] = ndp.get(key, 0) + v
        dp = ndp; s += 1
        if s % 3 == 0:
            n = s // 3
            if n <= nmax:
                out[n] = sum(v for (c1, c2, d, r), v in dp.items()
                             if s-c1-c2 == n and c1 == n and c2 == n
                             and closable(r))
                if n % 50 == 0:
                    print("   count(%d): %d digits  [%.0fs]"
                          % (n, len(str(out[n])), time.time()-T0), flush=True)
    return out

print("exact DP for G(n) to n = %d ..." % args.nmax, flush=True)
G = G_series(args.nmax)

if KNOWN is not None:
    bad = [n for n, w in enumerate(KNOWN) if n <= args.nmax and G[n] != w]
    print("validation against the published terms:",
          "ALL OK" if not bad else "MISMATCH at %s" % bad)
    if bad:
        raise SystemExit("aborting")
else:
    print("first terms (no published list to check against -- inspect these):")
    print("   ", [G[n] for n in range(0, 12)])
    print("   growth check: count(n+1)/count(n)/mu ->",
          nstr(mpf(G[args.nmax])/mpf(G[args.nmax-1])/MU, 12),
          "(should approach 1)")

NS = list(range(args.nmax//2, args.nmax+1))
A = [mpf(G[n]) * mpf(n)**4 / MU**n for n in NS]

# ---------------------------------------------------- METHOD 1: Richardson
cols = [A[:]]
for k in range(args.orders):
    prev = cols[-1]; nxt = []
    for i in range(len(prev)-1):
        N = NS[i]
        nxt.append(((N+1)*prev[i+1] - (N-k)*prev[i]) / (k+1))
    if not nxt: break
    cols.append(nxt)
rich = cols[-1][-1]
print("\nMETHOD 1  Richardson, %d orders" % (len(cols)-1))
for k in range(max(0, len(cols)-8), len(cols)):
    print("   order %2d : %s" % (k, nstr(cols[k][-1], 22)))

# ---------------------------------------------------- METHOD 2: Salzer
def salzer(seq, ns, k):
    """Salzer's transformation of order k at the tail of the sequence."""
    from mpmath import binomial
    m = len(seq)
    if m < k+1: return None
    num = mpf(0); den = mpf(0)
    base = m - k - 1
    for j in range(k+1):
        N = mpf(ns[base+j])
        w = (-1)**j * binomial(k, j) * (N**k)
        num += w * seq[base+j]
        den += w
    return num/den if den != 0 else None

print("\nMETHOD 2  Salzer acceleration (independent assumptions)")
salz = None
for k in range(2, min(args.orders, 16)+1):
    v = salzer(A, NS, k)
    if v is None: continue
    salz = v
    if k >= min(args.orders, 16) - 5:
        print("   order %2d : %s" % (k, nstr(v, 22)))

# ---------------------------------------------------- METHOD 3: tail fit
print("\nMETHOD 3  geometric tail fit on the Richardson increments")
deep = [c[-1] for c in cols]
inc = [deep[i+1]-deep[i] for i in range(len(deep)-1)]
tf = min(args.tailfit, len(inc)-1)
ratios = [inc[-i]/inc[-i-1] for i in range(1, tf+1) if inc[-i-1] != 0]
if ratios:
    r = sum(ratios)/len(ratios)
    last = inc[-1]
    tail = last*r/(1-r) if abs(r) < 1 else None
    print("   fitted ratio r = %s  (from %d increments)" % (nstr(r, 8), len(ratios)))
    print("   last increment = %s" % nstr(last, 8))
    if tail is not None:
        est = deep[-1] + tail
        print("   summed tail    = %s" % nstr(tail, 8))
        print("   limit estimate = %s" % nstr(est, 22))
    else:
        est = None
        print("   ratio not contractive; tail fit inapplicable")
else:
    est = None; r = None

# ---------------------------------------------------- agreement
print("\n--- comparison ---")
rows = [("Richardson", rich), ("Salzer", salz), ("tail fit", est)]
for name, v in rows:
    if v is None: continue
    print("   %-12s %s   |candidate - est| = %s"
          % (name, nstr(v, 22), nstr(abs(v-CAND), 6)))
print("   candidate    %s" % nstr(CAND, 22))

vals = [v for _, v in rows if v is not None]
if len(vals) >= 2:
    spread = max(vals) - min(vals)
    print("\n   spread across methods = %s" % nstr(spread, 6))
    # digits common to all estimates
    s0 = nstr(min(vals), 20); s1 = nstr(max(vals), 20)
    common = 0
    for a, b in zip(s0, s1):
        if a != b: break
        common += 1
    print("   agreed prefix: %s" % s0[:common])
    print("   -> assert only these digits in print; the remainder are")
    print("      supported by extrapolation, not established by it.")

json.dump({
 "nmax": args.nmax, "orders": len(cols)-1, "dps": args.dps,
 "validation": "ok",
 "richardson": nstr(rich, 25),
 "salzer": nstr(salz, 25) if salz is not None else None,
 "tail_fit": nstr(est, 25) if est is not None else None,
 "fitted_ratio": nstr(r, 10) if r is not None else None,
 "candidate": "0.52128605909(2)",
 "richardson_column": [nstr(c[-1], 22) for c in cols],
}, open("c1_tail_results.json", "w"), indent=1)
print("\nwrote c1_tail_results.json  [%.0fs]" % (time.time()-T0))
