#!/usr/bin/env python3
"""
thintail.py — P/N status of the thin-tail family, since the theory stalled.

WHY. We proved:
  * Staircase Theorem: winning bites of any bar form a strict antichain
    (i < i' => j > j'), hence at most one per row.
  * Transpose closure of two-value diagrams: (b^p, c^q)^T = (a^c, p^{b-c}),
    parameters (p,c) -> (c,p).
  * L-shapes: (b, 1^{a-1}) is P iff a = b; hence the (2,2) bite wins iff the
    bar is square; hence the number of winning moves on a square bar is ODD
    (one transpose-fixed point, everything else pairing).

The natural next family was the thin tails D = (b^p, 1^q). But that family is
NOT transpose-closed: its transpose is (p+q, p^{b-1}), a long row on top of a
rectangle, which leaves the family unless p = 1. So the mirroring argument
that settled the L-shapes does not generalise, and five hand-computed cases
were not enough to see a criterion:

    (b,p,q) = (2,1,1) P   (3,1,2) P   (2,2,1) P   (2,1,2) N   (2,2,2) N

(2,2,1) being P is a genuine three-parameter fact not implied by anything
proved. So: compute, then look.

THIS SCRIPT solves Chomp directly on the two families by memoised retrograde
search over general Young diagrams (partitions), which is exact but
exponential -- fine for the small sizes needed to see a pattern:

  1. D = (b^p, 1^q)          the thin tails
  2. D = (b^p, c^q)          the general two-value diagrams, small ranges,
                             printed as a (p,c) grid per (a,b) so the
                             ANTICHAIN structure is visible directly

and reports, for the thin tails, a table indexed by (b,p,q) plus an automatic
search for simple criteria (parity patterns, linear relations among b,p,q,
and whether P-ness depends only on some small combination).

Run:  python3 thintail.py            (pure python; seconds to a minute)
      python3 thintail.py --bmax 6 --pmax 6 --qmax 8
Writes thintail_results.json.
"""
import sys, json, argparse
from functools import lru_cache

ap = argparse.ArgumentParser()
ap.add_argument('--bmax', type=int, default=5)
ap.add_argument('--pmax', type=int, default=5)
ap.add_argument('--qmax', type=int, default=7)
ap.add_argument('--gridmax', type=int, default=7)
args = ap.parse_args()
sys.setrecursionlimit(100000)

def norm(rows):
    """Young diagram as a weakly decreasing tuple, zeros dropped."""
    r = tuple(sorted((x for x in rows if x > 0), reverse=True))
    return r

def moves(d):
    """All positions reachable by one bite. Biting cell (i,j) (1-indexed)
    removes every cell (i',j') with i' >= i and j' >= j."""
    out = set()
    for i in range(len(d)):
        for j in range(1, d[i] + 1):
            new = list(d[:i]) + [min(d[k], j - 1) for k in range(i, len(d))]
            out.add(norm(new))
    return out

@lru_cache(maxsize=None)
def is_P(d):
    """True if d is a previous-player win (P-position). The position with only
    the poison square, (1,), is a LOSS for the mover, i.e. N: the mover must
    take it. The empty diagram means the poison was taken -> previous player
    lost, so () is N for the previous player... we use the convention that the
    player to move at () has already won, hence () is a P-position for the
    mover's opponent. Standard: () is terminal; the player forced to take the
    poison loses. Model: a position is P (previous wins) iff every move leads
    to an N position; the empty position is P (the mover has no move and the
    opponent just ate the poison)."""
    if not d:
        return True
    for m in moves(d):
        if is_P(m):
            return False
    return True

out = {}

# ---------------------------------------------------------------- sanity
print("sanity checks")
# L-shapes: (b, 1^{a-1}) is P iff a == b
bad = []
for b in range(1, 8):
    for a in range(1, 8):
        d = norm([b] + [1]*(a-1))
        got = is_P(d)
        want = (a == b)
        if got != want:
            bad.append((a, b, got, want))
print("   L-shape lemma (P iff a==b):", "OK" if not bad else "FAILS %s" % bad[:5])
out["Lshape_lemma_ok"] = not bad
# full rectangles are N (Gale), except the single poison square
badr = []
for a in range(1, 6):
    for b in range(1, 6):
        if a == 1 and b == 1: continue
        if is_P(norm([b]*a)): badr.append((a, b))
print("   rectangles are N (Gale):", "OK" if not badr else "FAILS %s" % badr[:5])
out["gale_ok"] = not badr

# ---------------------------------------------------------------- thin tails
print("\nTHIN TAILS  D = (b^p, 1^q)")
print(" b  p  q   P/N")
rows = []
for b in range(2, args.bmax + 1):
    for p in range(1, args.pmax + 1):
        for q in range(0, args.qmax + 1):
            d = norm([b]*p + [1]*q)
            v = is_P(d)
            rows.append({"b": b, "p": p, "q": q, "P": bool(v)})
            if v:
                print(" %2d %2d %2d    P" % (b, p, q))
out["thin_tails"] = rows
Ps = [(r["b"], r["p"], r["q"]) for r in rows if r["P"]]
print("\n   P-positions found:", len(Ps), "of", len(rows))
print("  ", Ps[:40])

# ---- automatic hunt for a simple criterion
print("\n   looking for a simple criterion ...")
crit = {}
# does P depend only on (b, p+q)?
bybp = {}
ok1 = True
for r in rows:
    k = (r["b"], r["p"] + r["q"])
    if k in bybp and bybp[k] != r["P"]: ok1 = False
    bybp[k] = r["P"]
crit["depends_only_on_b_and_p_plus_q"] = ok1
# does P depend only on (p, q)?
bypq = {}
ok2 = True
for r in rows:
    k = (r["p"], r["q"])
    if k in bypq and bypq[k] != r["P"]: ok2 = False
    bypq[k] = r["P"]
crit["depends_only_on_p_and_q"] = ok2
# linear relations among b,p,q holding on all P-positions
for (cb, cp, cq, c0) in [(1,-1,0,0),(1,0,-1,0),(0,1,-1,0),(1,-1,-1,0),
                         (1,-1,-1,1),(1,-1,0,-1),(1,0,-1,-1),(1,-2,-1,0)]:
    if Ps and all(cb*b + cp*p + cq*q + c0 == 0 for (b,p,q) in Ps):
        crit["relation_%d_%d_%d_%d" % (cb,cp,cq,c0)] = True
        print("      all P satisfy  %d*b + %d*p + %d*q + %d = 0"
              % (cb,cp,cq,c0))
out["thin_tail_criteria"] = crit
for k, v in crit.items():
    if v is True and not k.startswith("relation"):
        print("      %s: TRUE" % k)

# ---------------------------------------------------------------- (p,c) grids
print("\nTWO-VALUE GRIDS  D = (b^p, c^q),  p+q = a  — the antichain picture")
grids = {}
for a in range(2, args.gridmax + 1):
    for b in range(2, args.gridmax + 1):
        cells = []
        for p in range(0, a):
            for c in range(0, b):
                q = a - p
                d = norm([b]*p + [c]*q)
                if is_P(d):
                    cells.append((p, c))
        if cells:
            grids["%dx%d" % (a, b)] = cells
            # check the antichain property
            anti = all(not (p1 <= p2 and c1 <= c2) and not (p2 <= p1 and c2 <= c1)
                       for i, (p1, c1) in enumerate(cells)
                       for (p2, c2) in cells[i+1:])
            print("   a=%d b=%d: P-cells (p,c) = %s   antichain: %s"
                  % (a, b, cells, "yes" if anti else "NO"))
out["two_value_grids"] = grids

json.dump(out, open("thintail_results.json", "w"), indent=1)
print("\nwrote thintail_results.json — send this back")
