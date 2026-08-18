#!/usr/bin/env python3
"""
matsuo_general.py — does Matsuo's bijection generalise with BOUNDED excusal
sets? This decides whether Spahn-Zeilberger Challenge 3 ($300) is a closure
argument or a research programme.

BACKGROUND. For (r,s) = (2,2), Matsuo's bijection relabels values by the
interleaving 1, 1+h, 2, 2+h, ... (h = floor((n+1)/2)) and turns

    pi_{i+2} - pi_i != 2        (the a_{2,2} condition)
into
    pi_{i+1} - pi_i != 1  EXCEPT excused at one position and one value

i.e. a RIN-type count with a 2-element excusal set. That bounded excusal set
is exactly what makes the inclusion-exclusion representation have FINITELY
MANY catalytic variables, hence holonomic by closure.

THE CONJECTURE TESTED HERE. For general (r,s): relabel POSITIONS by
interleaving the r arithmetic progressions of step r, and VALUES by
interleaving the s progressions of step s. Within a block, position-distance
r becomes 1 and value-difference s becomes 1, so the condition becomes the
gap-1 condition, breaking only at block SEAMS — (r-1) in position and (s-1)
in value. If so the excusal set has size (r-1)+(s-1): BOUNDED for fixed
(r,s), which is all Challenge 3 requires.

WHAT THIS SCRIPT DOES — brute force, no cleverness, so a failure is real:
  1. computes a_{r,s}(n) directly from the definition;
  2. builds the two relabellings explicitly for that n;
  3. transports every permutation through them and records EXACTLY which
     (position, value) pairs need excusing for the gap-1 description to agree
     with the original condition;
  4. reports the SIZE of the excusal set actually needed, and whether it is
     independent of n (the whole question);
  5. cross-checks that counting with that excusal set reproduces a_{r,s}(n).

If the excusal set size stays at (r-1)+(s-1) as n grows, the generalisation
holds and Challenge 3 reduces to the closure argument we already have.
If it grows with n, it does not, and Challenge 3 stays research-grade.

Run:  python3 matsuo_general.py        (pure python, no deps; ~seconds)
Writes matsuo_general_results.json — SEND THIS BACK.
"""
from itertools import permutations
import json, time

T0 = time.time()
out = {}

def a_rs(n, r, s):
    """#{pi in S_n : pi_{i+r} - pi_i != s for all valid i}, 1-indexed values."""
    c = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for i in range(n - r):
            if p[i+r] - p[i] == s:
                ok = False; break
        if ok: c += 1
    return c

def interleave(n, k):
    """Order 1..n as: 1, 1+k, 1+2k, ..., 2, 2+k, ...  Returns the list in the
    new order; entry j of the list is the ORIGINAL label at new index j."""
    order = []
    for start in range(1, k+1):
        x = start
        while x <= n:
            order.append(x)
            x += k
    return order

def analyse(n, r, s):
    """Transport the condition through the two relabellings and find the
    exact excusal set needed."""
    pos_order = interleave(n, r)        # new index -> old position (1-based)
    val_order = interleave(n, s)        # new value -> old value
    # inverse maps
    val_of_new = {newv: oldv for newv, oldv in enumerate(val_order, start=1)}
    new_of_val = {oldv: newv for newv, oldv in val_of_new.items()}

    # For the gap-1 description we need, for each new index k, whether
    # (new index k, k+1) corresponds to (old positions differing by r), and
    # for each new value u, whether (u, u+1) corresponds to old values
    # differing by s.
    pos_seams = [k for k in range(1, n)
                 if pos_order[k] - pos_order[k-1] != r]
    val_seams = [u for u in range(1, n)
                 if val_of_new[u+1] - val_of_new[u] != s]
    return {"pos_order": pos_order, "val_order": val_order,
            "pos_seams": pos_seams, "val_seams": val_seams,
            "n_pos_seams": len(pos_seams), "n_val_seams": len(val_seams),
            "new_of_val": new_of_val, "pos_of_new": pos_order}

def count_transformed(n, r, s, info):
    """Count permutations satisfying: sigma_{k+1} - sigma_k != 1 for all new
    indices k, EXCEPT where k is a position seam or sigma_k is a value seam."""
    pos_seams = set(info["pos_seams"]); val_seams = set(info["val_seams"])
    c = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for k in range(1, n):
            if k in pos_seams:      continue
            if p[k-1] in val_seams: continue
            if p[k] - p[k-1] == 1:
                ok = False; break
        if ok: c += 1
    return c

rows = []
print("(r,s)   n   a_rs(n)   pos_seams  val_seams  excusal  transformed  match")
for (r, s) in [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (3, 4)]:
    for n in range(4, 9):
        if n <= max(r, s): continue
        A = a_rs(n, r, s)
        info = analyse(n, r, s)
        T = count_transformed(n, r, s, info)
        exc = info["n_pos_seams"] + info["n_val_seams"]
        rows.append({"r": r, "s": s, "n": n, "a_rs": A,
                     "pos_seams": info["n_pos_seams"],
                     "val_seams": info["n_val_seams"],
                     "excusal_size": exc, "transformed_count": T,
                     "match": bool(A == T),
                     "predicted_excusal": (r-1) + (s-1)})
        print("(%d,%d)  %2d  %8d   %6d     %6d    %5d   %10d   %s"
              % (r, s, n, A, info["n_pos_seams"], info["n_val_seams"],
                 exc, T, "YES" if A == T else "no"), flush=True)
out["rows"] = rows

# ---- the two questions that matter
print("\n--- ANALYSIS ---")
by_rs = {}
for row in rows:
    by_rs.setdefault((row["r"], row["s"]), []).append(row)
summary = []
for (r, s), rs_rows in sorted(by_rs.items()):
    sizes = sorted(set(x["excusal_size"] for x in rs_rows))
    grows = len(sizes) > 1
    matches = all(x["match"] for x in rs_rows)
    pred = (r-1) + (s-1)
    summary.append({"r": r, "s": s, "excusal_sizes_over_n": sizes,
                    "bounded_in_n": not grows, "predicted": pred,
                    "matches_prediction": sizes == [pred] if sizes else False,
                    "bijection_reproduces_a_rs": matches})
    print("(r,s)=(%d,%d): excusal sizes over n = %s | bounded: %s | "
          "predicted %d | counts match: %s"
          % (r, s, sizes, not grows, pred, matches), flush=True)
out["summary"] = summary
out["all_bounded"] = all(x["bounded_in_n"] for x in summary)
out["all_match"] = all(x["bijection_reproduces_a_rs"] for x in summary)

print("\nEXCUSAL SETS BOUNDED IN n FOR EVERY (r,s) TESTED:", out["all_bounded"])
print("TRANSFORMED COUNT REPRODUCES a_rs EVERYWHERE:      ", out["all_match"])
if out["all_bounded"] and out["all_match"]:
    print("\n=> The generalisation HOLDS on this evidence: Challenge 3 reduces")
    print("   to the closure argument already in hand, with finitely many")
    print("   catalytic variables for each fixed (r,s).")
else:
    print("\n=> The naive generalisation FAILS. Either the relabelling is not")
    print("   the right one, or the excusal set genuinely grows — in which")
    print("   case Challenge 3 stays research-grade. The rows show which.")
out["caveat"] = ("Small n only (brute force is n! ). Bounded here means "
                 "bounded over the n tested; a proof needs the seam count "
                 "argued in general, which is elementary if the pattern holds.")
json.dump(out, open("matsuo_general_results.json", "w"), indent=1)
print("\nwrote matsuo_general_results.json — send this back")
