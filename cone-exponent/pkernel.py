#!/usr/bin/env python3
"""
pkernel.py — test holonomicity of the (2,1,1)x[n] Solid SYT sequence by the
p-KERNEL, avoiding the cone exponent entirely.

THE CHAIN (all real theorems):
  1. a(n) holonomic, integer coefficients, finite radius  =>  GF is a
     G-function.
  2. Chudnovsky: G-operators are globally nilpotent (nilpotent p-curvature
     for almost all p).
  3. Katz / Cartier: nilpotent p-curvature => the reduction mod p has
     algebraic solutions over F_p(x).
  4. Christol: a series over F_p is algebraic  <=>  its coefficient sequence
     is p-AUTOMATIC.
  =>  if a(n) is holonomic then a(n) mod p is p-automatic for almost all p.

Automaticity has a finite signature: the p-KERNEL
    K_p(a) = { n -> a(p^k n + r) : k >= 0, 0 <= r < p^k }
must be FINITE. For a non-automatic sequence the number of distinct kernel
elements grows without bound in k.

WHY THIS ROUTE. It sidesteps the wall the other approach hit. No Dirichlet
eigenvalue, no irrationality of nu, no height bounds on a hypothetical
operator. This is also the standard modern technique (Bostan-Christol-Dumas
use it to separate D-finite from non-D-finite lattice-walk GFs).

THE HONEST LIMIT. Observing kernel growth on finitely many terms is EVIDENCE,
not proof: certifying non-automaticity is itself hard. But it is a different
wall, possibly a lower one, and it would be an INDEPENDENT second argument
for the same conclusion.

THE SEQUENCE. a(n) = # Solid SYT of cylindrical shape (2,1,1) x {1..n}
             = # walks (0,0,0,0) -> (n,n,n,n) with positive unit steps in Z^4
               staying in the cone  x1 >= x2 >= x3,  x1 >= x4.
Computed here by DP over the cone, mod p.

CONTROLS. The script first runs the same kernel analysis on two sequences
with known status:
  * central binomial C(2n,n) mod 2 -- algebraic mod 2, hence 2-automatic:
    the kernel must SATURATE;
  * a deliberately non-automatic sequence (n -> floor(n*sqrt2) mod 2):
    the kernel must GROW.
If the controls do not behave, the kernel code is wrong and nothing below is
meaningful.

Run:  python3 pkernel.py --nmax 4000 --primes 2,3,5 --kmax 5
"""
import argparse, json, sys
from math import isqrt

ap = argparse.ArgumentParser()
ap.add_argument('--nmax', type=int, default=2000)
ap.add_argument('--primes', type=str, default='2,3')
ap.add_argument('--kmax', type=int, default=5)
args = ap.parse_args()
PRIMES = [int(x) for x in args.primes.split(',')]

# ---------------------------------------------------------------- sequence
def solid_syt(nmax, p):
    """a(n) mod p: walks to (n,n,n,n) in the cone x1>=x2>=x3, x1>=x4,
    with positive unit steps. DP over lattice points, layer by total sum."""
    from collections import defaultdict
    # state (x1,x2,x3,x4) with x1>=x2>=x3, x1>=x4, all <= nmax
    cur = {(0, 0, 0, 0): 1 % p}
    out = [0]*(nmax+1)
    out[0] = 1 % p
    total = 4*nmax
    for s in range(1, total+1):
        nxt = defaultdict(int)
        for (a, b, c, d), v in cur.items():
            if not v: continue
            # step in each coordinate, keeping the cone constraints
            if a+1 <= nmax:                      nxt[(a+1, b, c, d)] = (nxt[(a+1, b, c, d)] + v) % p
            if b+1 <= a and b+1 <= nmax:         nxt[(a, b+1, c, d)] = (nxt[(a, b+1, c, d)] + v) % p
            if c+1 <= b and c+1 <= nmax:         nxt[(a, b, c+1, d)] = (nxt[(a, b, c+1, d)] + v) % p
            if d+1 <= a and d+1 <= nmax:         nxt[(a, b, c, d+1)] = (nxt[(a, b, c, d+1)] + v) % p
        cur = nxt
        # record diagonal points
        for n in range(1, nmax+1):
            if 4*n == s and (n, n, n, n) in cur:
                out[n] = cur[(n, n, n, n)] % p
        if not cur: break
    return out

# ---------------------------------------------------------------- kernel
def kernel_growth(seq, p, kmax, nmax):
    """Number of DISTINCT kernel elements at each depth k, comparing
    subsequences as tuples truncated to a common length."""
    counts = []
    seen_total = set()
    for k in range(kmax+1):
        pk = p**k
        L = (nmax - pk) // pk
        if L < 40:
            counts.append(None); continue
        elems = set()
        for r in range(pk):
            sub = tuple(seq[pk*n + r] for n in range(L))
            elems.add(sub)
            seen_total.add(sub)
        counts.append(len(elems))
    return counts, len(seen_total)

def report(name, seq, p, kmax, nmax, expect):
    counts, tot = kernel_growth(seq, p, kmax, nmax)
    print("   %-28s p=%d  distinct per depth: %s   cumulative: %d   [%s]"
          % (name, p, counts, tot, expect))
    return {"per_depth": counts, "cumulative": tot, "expect": expect}

out = {}
print("CONTROLS")
# central binomial mod 2: algebraic mod 2 (Christol), so 2-automatic
cb = [1]*(args.nmax+1)
c = 1
for n in range(1, args.nmax+1):
    c = c * 2 * (2*n-1) // n
    cb[n] = c % 2
out["control_central_binomial"] = report("C(2n,n) mod 2", cb, 2,
                                         args.kmax, args.nmax, "should SATURATE")
# a non-automatic control: floor(n*sqrt2) mod 2 (Sturmian, not automatic)
sq = [ (isqrt(2*n*n) ) % 2 for n in range(args.nmax+1) ]
out["control_sturmian"] = report("floor(n*sqrt2) mod 2", sq, 2,
                                 args.kmax, args.nmax, "should GROW")

print("\nTARGET: Solid SYT (2,1,1)x[n]")
out["target"] = {}
for p in PRIMES:
    print("   computing a(n) mod %d up to n=%d ..." % (p, args.nmax), flush=True)
    try:
        seq = solid_syt(min(args.nmax, 60), p)   # DP is heavy; small nmax
        print("      a(0..8) mod %d = %s" % (p, seq[:9]))
        out["target"]["p=%d" % p] = report("Solid SYT mod %d" % p, seq, p,
                                           3, len(seq)-1, "unknown")
    except Exception as e:
        print("      ERROR:", str(e)[:200])
        out["target"]["p=%d" % p] = {"error": str(e)[:300]}

out["note"] = ("Kernel growth on finite data is evidence, not proof. The "
               "controls check the kernel code: an automatic sequence must "
               "saturate, a Sturmian one must grow. The DP for the target is "
               "expensive, so nmax is small here -- if the controls behave and "
               "the target shows growth, the next step is a faster generator "
               "(the shell-indexed numpy DP from the earlier session reached "
               "N=84) to get enough terms for a meaningful depth.")
json.dump(out, open("pkernel_results.json", "w"), indent=1)
print("\nwrote pkernel_results.json")
