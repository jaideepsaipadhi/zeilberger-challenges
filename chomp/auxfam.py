#!/usr/bin/env python3
"""
auxfam.py — the two auxiliary families needed to close the p=3 induction.

WHY THESE TWO. From D = (b^3, 1^q) the moves land in only four places:
    (1,j)          -> ((j-1)^3, 1^q)        same family, smaller b
    (i,1) / tail   -> (b^3, 1^q'), q'<q     same family, or a rectangle (Gale)
    (2,j)          -> (b, (j-1)^2, 1^q)  =  A(b, j-1, q)
    (3,j)          -> (b, b, j-1, 1^q)   =  B(b, j-1, q)
So the p=3 thin tails are ALMOST closed: exactly two auxiliary shapes appear,

    A(b,c,q) = (b, c, c, 1^q),      B(b,c,q) = (b, b, c, 1^q),    1 <= c < b

and nothing beyond them. (The unbounded escape hierarchy afflicts GENERAL
two-value diagrams; for this slice one extra layer closes it.)

To prove the observed law  q(b,3) = 2b-3  by induction on b, we need every
move from (b^3, 1^{2b-3}) to land in N, i.e.:
  (1) (b'^3, 1^q) is N for b' < b whenever q != 2b'-3   [induction + the
      proved uniqueness of thin-tail tracks];
  (2) A(b,c,q) and B(b,c,q) are N for all reachable (c,q).
This script tests (2) and looks for the laws governing A and B.

For each family it reports, per (b,c), the set of q making it a P-position --
by the same chain argument used for thin tails, moving q -> q' < q stays in
the family, so at most ONE q per (b,c) can be P.

Run:  python3 auxfam.py --bmax 8 --qmax 30
"""
import sys, json, argparse
from functools import lru_cache

ap = argparse.ArgumentParser()
ap.add_argument('--bmax', type=int, default=7)
ap.add_argument('--qmax', type=int, default=24)
args = ap.parse_args()
sys.setrecursionlimit(200000)

def norm(rows):
    return tuple(sorted((x for x in rows if x > 0), reverse=True))

@lru_cache(maxsize=None)
def is_P(d):
    if not d:
        return False
    for i in range(len(d)):
        for j in range(1, d[i] + 1):
            m = norm(list(d[:i]) + [min(d[k], j-1) for k in range(i, len(d))])
            if is_P(m):
                return False
    return True

print("sanity: L-shape lemma (P iff a==b)")
bad = [(a,b) for b in range(1,7) for a in range(1,7)
       if is_P(norm([b] + [1]*(a-1))) != (a == b)]
print("   ", "OK" if not bad else "FAILS %s" % bad[:4])
print("sanity: thin-tail p=3 law q = 2b-3")
bad3 = []
for b in range(2, 7):
    qs = [q for q in range(0, args.qmax+1) if is_P(norm([b]*3 + [1]*q))]
    if qs != [2*b-3] and not (b == 2 and qs == [1]):
        bad3.append((b, qs))
print("   ", "OK" if not bad3 else "MISMATCH %s" % bad3)

out = {"A": {}, "B": {}}

for name, build in (("A", lambda b, c, q: [b] + [c]*2 + [1]*q),
                    ("B", lambda b, c, q: [b]*2 + [c] + [1]*q)):
    print("\n=== family %s ===" % name)
    print(" b  c    q with P     (at most one, by the chain argument)")
    for b in range(2, args.bmax + 1):
        for c in range(1, b):
            qs = [q for q in range(0, args.qmax + 1)
                  if is_P(norm(build(b, c, q)))]
            out[name]["%d,%d" % (b, c)] = qs
            print(" %2d %2d   %s" % (b, c, qs if qs else "none"), flush=True)

# ---- do A and B ever have a P at all? that is what matters for the induction
print("\n--- what the induction needs ---")
for name in ("A", "B"):
    withP = {k: v for k, v in out[name].items() if v}
    print(" family %s: %d of %d (b,c) pairs have a P-position"
          % (name, len(withP), len(out[name])))
    if withP:
        print("    those pairs:", list(withP.items())[:14])
    else:
        print("    NONE -> every A/B shape is N, and requirement (2) is")
        print("    automatic; the induction would close on (1) alone.")

json.dump(out, open("auxfam_results.json", "w"), indent=1)
print("\nwrote auxfam_results.json")
