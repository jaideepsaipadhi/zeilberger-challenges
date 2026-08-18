#!/usr/bin/env python3
"""
thintail2.py — extend the thin-tail table, test the three linear laws.

WHAT thintail.py FOUND (after the base-case fix; both sanity checks pass —
the L-shape lemma and Gale's rectangles-are-N):

    D = (b^p, 1^q) is P for exactly one q per (b,p)  [PROVED: the family is
    totally ordered under moves, since (b^p,1^q) -> (b^p,1^q') for any q'<q,
    and no P-position moves to a P-position]

    p | q(b,p) for b = 2..6      | law
    --+---------------------------+---------------------------
    1 | 1, 2, 3, 4, 5             | q = b-1        (L-shape lemma)
    2 | 1, 2, 4, 5, 7             | slope 3/2 — the integers not divisible
      |                           | by 3, i.e. quasi-linear (Beatty)
    3 | 1, 3, 5, 7, 9             | q = 2b-3       5/5
    4 | 1, 4, 7, 10, 13           | q = 3b-5       5/5
    5 | 1, 5, 7, 13, 15           | irregular
    6 | 1, 5, 10, 9, 16           | irregular, non-monotonic in b

THIS SCRIPT extends b to test whether 2b-3 and 3b-5 survive, whether p=2
really is the not-divisible-by-3 sequence, and whether p=5,6 resolve with
more data or are genuinely irregular. It also re-checks the p=6 row, whose
non-monotonicity (10 then 9) is either a real phenomenon or worth a second
look.

EFFICIENCY. The previous version memoised over general partitions, which is
wasteful here: every diagram reachable from (b^p, 1^q) has parts <= b, so the
state space is partitions with bounded part size. This version searches
directly with an explicit memo keyed on the tuple, iterating q upward per
(b,p) and stopping at the first P — which by the uniqueness proof is the only
one.

Run:  python3 thintail2.py --bmax 10 --pmax 6 --qmax 40
"""
import sys, json, argparse
from functools import lru_cache

ap = argparse.ArgumentParser()
ap.add_argument('--bmax', type=int, default=9)
ap.add_argument('--pmax', type=int, default=6)
ap.add_argument('--qmax', type=int, default=40)
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

print("sanity: L-shape lemma")
bad = [(a,b) for b in range(1,7) for a in range(1,7)
       if is_P(norm([b] + [1]*(a-1))) != (a == b)]
print("   ", "OK" if not bad else "FAILS %s" % bad[:4])

rows = {}
print("\n b  p   q      (first P per (b,p); unique by the chain argument)")
for p in range(1, args.pmax + 1):
    for b in range(2, args.bmax + 1):
        found = None
        for q in range(0, args.qmax + 1):
            if is_P(norm([b]*p + [1]*q)):
                found = q; break
        rows[(b, p)] = found
        print(" %2d %2d  %s" % (b, p, found if found is not None else "> qmax"),
              flush=True)

print("\n--- rows in p ---")
laws = {}
for p in range(1, args.pmax + 1):
    seq = [rows[(b, p)] for b in range(2, args.bmax + 1)]
    print(" p=%d: %s" % (p, seq))
    if any(v is None for v in seq):
        continue
    # test q = A*b + B from the first two, verify on the rest
    A = seq[1] - seq[0]
    B = seq[0] - A*2
    if all(seq[i] == A*(i+2) + B for i in range(len(seq))):
        print("      LINEAR: q = %d*b %+d   (holds on all %d values)"
              % (A, B, len(seq)))
        laws["p=%d" % p] = "q = %d*b %+d" % (A, B)
    else:
        # quasi-linear? fit a slope by least squares on the endpoints
        sl = (seq[-1] - seq[0]) / (len(seq) - 1)
        dev = [seq[i] - (seq[0] + sl*i) for i in range(len(seq))]
        print("      not linear; endpoint slope %.4f, deviations %s"
              % (sl, [round(d, 2) for d in dev]))
        laws["p=%d" % p] = "slope ~%.4f, deviations %s" % (sl, [round(d,2) for d in dev])

json.dump({"table": {"%d,%d" % k: v for k, v in rows.items()}, "laws": laws},
          open("thintail2_results.json", "w"), indent=1)
print("\nwrote thintail2_results.json")
