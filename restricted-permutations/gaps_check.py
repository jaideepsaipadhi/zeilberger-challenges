#!/usr/bin/env python3
"""
gaps_check.py — close the two gaps in holonomicity-note.tex.

GAP 1 (the one that threatens the theorem). Step (3) of the holonomicity
proof claims: the excused positions and values vary LINEARLY with n within
each residue class of n modulo lcm(r,s). Sections of holonomic sequences are
holonomic, so this is what lets us specialise the multivariate count to the
diagonal. It was asserted, not verified.

  The argument, now checked here: for the interleaving iota_{n,k}, block c is
  {c, c+k, c+2k, ...} cap [1,n], of length floor((n-c)/k) + 1, and the seams
  sit at the cumulative sums of those lengths. Each length is linear in n on
  each residue class of n mod k, hence so is each seam index. Positions give
  classes mod r, values mod s, so jointly mod lcm(r,s).

  This script computes the seam indices directly for a range of n, splits by
  residue class mod lcm(r,s), and CHECKS each is an exact affine function of n
  on that class (fits from two points, verifies on the rest).

GAP 2 (cosmetic but the table should be complete). The b-case was not run for
(3,4) and a couple of other pairs. This fills them in, using the symmetrised
excusal rule of Proposition 3(ii).

Run:  python3 gaps_check.py       (pure python, seconds)
Writes gaps_check_results.json — SEND THIS BACK.
"""
from itertools import permutations
from math import gcd
import json

def lcm(a, b): return a*b//gcd(a, b)

def interleave(n, k):
    order = []
    for start in range(1, k+1):
        x = start
        while x <= n:
            order.append(x); x += k
    return order

def seam_indices(n, k):
    """Indices j in 1..n-1 where iota_{n,k}(j+1) - iota_{n,k}(j) != k."""
    o = interleave(n, k)
    return [j for j in range(1, n) if o[j] - o[j-1] != k]

out = {}

# ---------------------------------------------------------------- GAP 1
print("GAP 1 — are the seam indices affine in n on residue classes?\n")
gap1 = []
for (r, s) in [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (3, 4), (4, 3), (4, 4)]:
    L = lcm(r, s)
    ok_all = True
    detail = []
    for k, name in ((r, "position"), (s, "value")):
        # seams for interleaving with step k; expect k-1 of them
        for cls in range(L):
            ns = [n for n in range(max(2*L, k+2), max(2*L, k+2) + 8*L) if n % L == cls]
            if len(ns) < 4: continue
            seqs = [seam_indices(n, k) for n in ns]
            counts = set(len(x) for x in seqs)
            if counts != {k-1}:
                ok_all = False
                detail.append({"which": name, "class": cls,
                               "error": "seam count not constant: %s" % sorted(counts)})
                continue
            # each seam slot should be affine in n on this class
            for slot in range(k-1):
                ys = [seqs[i][slot] for i in range(len(ns))]
                # fit from the first two points: y = alpha*(n) + beta
                dn = ns[1] - ns[0]
                if (ys[1] - ys[0]) % dn != 0:
                    alpha = (ys[1] - ys[0]) / dn
                else:
                    alpha = (ys[1] - ys[0]) // dn
                beta = ys[0] - alpha*ns[0]
                bad = [ns[i] for i in range(len(ns)) if alpha*ns[i] + beta != ys[i]]
                if bad:
                    ok_all = False
                    detail.append({"which": name, "class": cls, "slot": slot,
                                   "error": "not affine", "failures": bad[:4]})
                else:
                    detail.append({"which": name, "class": cls, "slot": slot,
                                   "affine": "j = %s*n + %s" % (alpha, beta)})
    gap1.append({"r": r, "s": s, "lcm": L, "all_affine": ok_all,
                 "sample": detail[:6]})
    print("  (r,s)=(%d,%d)  lcm=%d  -> %s"
          % (r, s, L, "ALL AFFINE" if ok_all else "FAILURE (see json)"), flush=True)
    if ok_all:
        for d in detail[:2]:
            if "affine" in d:
                print("        %s seam: %s" % (d["which"], d["affine"]), flush=True)
out["gap1"] = gap1
out["gap1_all_affine"] = all(g["all_affine"] for g in gap1)
print("\n  GAP 1 CLOSED:", out["gap1_all_affine"])

# ---------------------------------------------------------------- GAP 2
print("\nGAP 2 — fill the missing b-case table entries\n")

def b_rs(n, r, s):
    c = 0
    for p in permutations(range(1, n+1)):
        if all(abs(p[i+r]-p[i]) != s for i in range(n-r)): c += 1
    return c

def a_rs(n, r, s):
    c = 0
    for p in permutations(range(1, n+1)):
        if all(p[i+r]-p[i] != s for i in range(n-r)): c += 1
    return c

def transformed(n, r, s, absolute):
    po = interleave(n, r); vo = interleave(n, s)
    val_of_new = {i+1: vo[i] for i in range(n)}
    ps = set(k for k in range(1, n) if po[k] - po[k-1] != r)
    vup = set(u for u in range(1, n) if val_of_new[u+1] - val_of_new[u] != s)
    vdn = set(u for u in range(2, n+1) if val_of_new[u-1] - val_of_new[u] != -s)
    c = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for k in range(1, n):
            if k in ps: continue
            a_, b_ = p[k-1], p[k]
            if absolute:
                if abs(b_-a_) != 1: continue
                if b_ == a_+1 and a_ in vup: continue
                if b_ == a_-1 and a_ in vdn: continue
            else:
                if b_-a_ != 1: continue
                if a_ in vup: continue
            ok = False; break
        if ok: c += 1
    return c, len(ps), len(vup), len(vdn)

gap2 = []
print("(r,s)   n   a_rs  transf  |   b_rs  transf  | exc_a exc_b  match")
for (r, s) in [(3, 4), (4, 3), (4, 4), (2, 2), (2, 3)]:
    for n in range(5, 9):
        A = a_rs(n, r, s); Bv = b_rs(n, r, s)
        tA, ps, vu, vd = transformed(n, r, s, False)
        tB, _, _, _ = transformed(n, r, s, True)
        ea, eb = ps+vu, ps+vu+vd
        gap2.append({"r": r, "s": s, "n": n, "a_rs": A, "transf_a": tA,
                     "b_rs": Bv, "transf_b": tB, "excusal_a": ea,
                     "excusal_b": eb, "match_a": A == tA, "match_b": Bv == tB,
                     "predicted_a": (r-1)+(s-1), "predicted_b": (r-1)+2*(s-1)})
        print("(%d,%d)  %2d %6d %7d  | %6d %7d  | %4d %5d   %s"
              % (r, s, n, A, tA, Bv, tB, ea, eb,
                 "YES" if (A == tA and Bv == tB) else "NO"), flush=True)
out["gap2"] = gap2
out["gap2_all_match"] = all(x["match_a"] and x["match_b"] for x in gap2)
out["gap2_sizes_as_predicted"] = all(
    x["excusal_a"] == x["predicted_a"] and x["excusal_b"] == x["predicted_b"]
    for x in gap2)
print("\n  GAP 2 counts all match:", out["gap2_all_match"])
print("  excusal sizes match (r-1)+(s-1) and (r-1)+2(s-1):",
      out["gap2_sizes_as_predicted"])

print("\nBOTH GAPS CLOSED:", out["gap1_all_affine"] and out["gap2_all_match"]
      and out["gap2_sizes_as_predicted"])
json.dump(out, open("gaps_check_results.json", "w"), indent=1)
print("wrote gaps_check_results.json — send this back")
