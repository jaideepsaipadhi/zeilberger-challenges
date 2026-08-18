#!/usr/bin/env python3
"""
matsuo_b.py — the remaining gap: does the generalisation cover b_{r,s} too?

CONTEXT. matsuo_general.py established (brute force, n = 4..8, 32 rows) that
for a_{r,s} the relabelling
    positions: interleave the r arithmetic progressions of step r
    values:    interleave the s progressions of step s
turns "pi_{i+r} - pi_i != s" into the gap-1 condition with an excusal set of
size EXACTLY (r-1) + (s-1), CONSTANT IN n, and the transformed count
reproduces a_{r,s}(n) exactly for (2,2), (2,3), (3,2), (3,3), (2,4), (4,2),
(3,4). Bounded excusal sets => finitely many catalytic variables => holonomic
by the closure argument already proved for a_{2,2}.

THE GAP. Challenge 3 asks about a_{r,s} AND b_{r,s}, where
    b_{r,s}(n) = #{pi : |pi_{i+r} - pi_i| != s}
forbids BOTH +s and -s. The value-side relabelling must now handle two signs
at once, so the seam count may differ — plausibly ~2(s-1), still bounded, but
that is a guess and today has been a lesson in not trusting guesses.

WHAT THIS SCRIPT DOES (same brute force, no cleverness):
  1. computes b_{r,s}(n) directly from the definition;
  2. applies the SAME two relabellings;
  3. counts permutations avoiding |sigma_{k+1} - sigma_k| = 1 except at the
     seams, under three candidate excusal rules:
       (A) position seams + value seams, as in the a-case
       (B) (A) plus the mirror value seams (value u such that u-1 crosses a
           block boundary downward) — the natural two-sign version
       (C) (B) plus excusing when the PREVIOUS value sits at a seam
  4. reports which rule (if any) reproduces b_{r,s}(n), and the size of the
     excusal set it needs, and whether that size is constant in n.

If some rule reproduces b_{r,s} with an n-independent excusal set, Challenge
3 is fully covered. If none does, the b-half needs its own analysis and the
claim must be restricted to a_{r,s}.

Run:  python3 matsuo_b.py        (pure python, seconds)
Writes matsuo_b_results.json — SEND THIS BACK.
"""
from itertools import permutations
import json, time

T0 = time.time()
out = {}

def b_rs(n, r, s):
    c = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for i in range(n - r):
            if abs(p[i+r] - p[i]) == s:
                ok = False; break
        if ok: c += 1
    return c

def interleave(n, k):
    order = []
    for start in range(1, k+1):
        x = start
        while x <= n:
            order.append(x); x += k
    return order

def seams(n, r, s):
    pos_order = interleave(n, r)
    val_order = interleave(n, s)
    val_of_new = {i+1: val_order[i] for i in range(n)}
    pos_seams = set(k for k in range(1, n) if pos_order[k] - pos_order[k-1] != r)
    # forward value seams: new value u where u -> u+1 is NOT an old +s step
    val_seams_up = set(u for u in range(1, n)
                       if val_of_new[u+1] - val_of_new[u] != s)
    # backward: new value u where u -> u-1 is NOT an old -s step
    val_seams_dn = set(u for u in range(2, n+1)
                       if val_of_new[u-1] - val_of_new[u] != -s)
    return pos_seams, val_seams_up, val_seams_dn

def count_rule(n, pos_seams, vs_up, vs_dn, rule):
    c = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for k in range(1, n):
            if k in pos_seams:
                continue
            a_, b_ = p[k-1], p[k]
            if abs(b_ - a_) != 1:
                continue
            excused = False
            if rule in ('A', 'B', 'C') and a_ in vs_up and b_ == a_ + 1:
                excused = True
            if rule in ('B', 'C') and a_ in vs_dn and b_ == a_ - 1:
                excused = True
            if rule == 'C' and (a_ in vs_up or a_ in vs_dn):
                excused = True
            if not excused:
                ok = False; break
        if ok: c += 1
    return c

rows = []
print("(r,s)   n    b_rs(n)   pos  vup  vdn |   ruleA    ruleB    ruleC  | match")
for (r, s) in [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2)]:
    for n in range(4, 9):
        if n <= max(r, s): continue
        Bv = b_rs(n, r, s)
        ps, vu, vd = seams(n, r, s)
        cA = count_rule(n, ps, vu, vd, 'A')
        cB = count_rule(n, ps, vu, vd, 'B')
        cC = count_rule(n, ps, vu, vd, 'C')
        which = [nm for nm, val in (('A', cA), ('B', cB), ('C', cC)) if val == Bv]
        rows.append({"r": r, "s": s, "n": n, "b_rs": Bv,
                     "pos_seams": len(ps), "val_up": len(vu), "val_dn": len(vd),
                     "ruleA": cA, "ruleB": cB, "ruleC": cC,
                     "rules_matching": which})
        print("(%d,%d)  %2d  %8d   %3d  %3d  %3d | %8d %8d %8d  | %s"
              % (r, s, n, Bv, len(ps), len(vu), len(vd), cA, cB, cC,
                 ",".join(which) if which else "NONE"), flush=True)
out["rows"] = rows

print("\n--- ANALYSIS ---")
summary = []
by = {}
for row in rows:
    by.setdefault((row["r"], row["s"]), []).append(row)
for (r, s), rr in sorted(by.items()):
    common = set(['A', 'B', 'C'])
    for row in rr:
        common &= set(row["rules_matching"])
    sizes = sorted(set(x["pos_seams"] + x["val_up"] + x["val_dn"] for x in rr))
    summary.append({"r": r, "s": s, "rules_working_for_all_n": sorted(common),
                    "excusal_sizes_over_n": sizes,
                    "bounded_in_n": len(sizes) == 1})
    print("(r,s)=(%d,%d): rules working for all n = %s | excusal sizes = %s | "
          "bounded: %s" % (r, s, sorted(common) or "NONE", sizes,
                           len(sizes) == 1), flush=True)
out["summary"] = summary
out["some_rule_works_everywhere"] = all(x["rules_working_for_all_n"] for x in summary)
out["all_bounded"] = all(x["bounded_in_n"] for x in summary)

print("\nA SINGLE RULE REPRODUCES b_rs FOR EVERY (r,s), ALL n:",
      out["some_rule_works_everywhere"])
print("EXCUSAL SETS BOUNDED IN n:                            ", out["all_bounded"])
if out["some_rule_works_everywhere"] and out["all_bounded"]:
    print("\n=> b_{r,s} is covered too: Challenge 3 ($300) reduces IN FULL to")
    print("   the closure argument already proved for a_{2,2}.")
else:
    print("\n=> b_{r,s} is NOT covered by these rules. Either a different")
    print("   relabelling is needed for the two-sign case, or the b-half")
    print("   requires separate analysis — restrict any claim to a_{r,s}")
    print("   until this is resolved.")
json.dump(out, open("matsuo_b_results.json", "w"), indent=1)
print("\nwrote matsuo_b_results.json — send this back")
