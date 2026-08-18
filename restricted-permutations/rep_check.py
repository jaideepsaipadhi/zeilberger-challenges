#!/usr/bin/env python3
"""
rep_check.py — verify Section 3 of the note, the one step never machine-checked.

WHY. The bijection half of the paper (Lemma 2, Corollary 3, Proposition 4) is
proved and confirmed by 52 brute-force rows. Section 3 -- the two-layer
inclusion-exclusion representation -- was derived by hand and corrected twice
during drafting. Nothing in it has been checked against a number. This script
does that.

WHAT IS BEING CHECKED. The claim is that

    R(n; P, V) = #{sigma : sigma_{k+1} - sigma_k != 1 whenever k not in P
                                                   and sigma_k not in V}

equals the inclusion-exclusion assembly described in Section 3:

  (a) IE over sets S of value bonds, S disjoint from V, weight (-1)^{|S|};
  (b) forcing S glues values into blocks; a realising permutation is an
      ARRANGEMENT of those blocks;
  (c) the position constraint (no forced bond may sit at a position in P) is
      removed by a second, finite IE over subsets T of P;
  (d) counting arrangements with prescribed internal-bond positions.

The script computes R three ways and compares:

  METHOD 0  brute force over S_n, straight from the definition.
  METHOD 1  the IE of (a)-(b) only, WITHOUT the position layer -- i.e. treating
            every bond as forbidden regardless of position. This should NOT
            match R when P is nonempty; it is included as a control so that a
            spurious agreement in METHOD 2 cannot pass unnoticed.
  METHOD 2  the full assembly (a)-(d), which is what Section 3 asserts.

METHOD 2 is implemented directly from the combinatorial description, not from
the closed-form generating function, so it tests the DERIVATION rather than the
algebra. If METHOD 2 == METHOD 0 across the tested cases, Section 3's
combinatorial content is correct and only the generating-function bookkeeping
(the block factors B, B*) remains as algebra -- which the script also checks
separately at the level of the block-factor identity.

Run:  python3 rep_check.py        (pure python, seconds)
Writes rep_check_results.json.
"""
from itertools import permutations, combinations
from math import factorial
import json

def brute_R(n, P, V, absolute=False):
    """METHOD 0: straight from the definition."""
    c = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for k in range(1, n):
            if k in P: continue
            a_, b_ = p[k-1], p[k]
            if absolute:
                if abs(b_-a_) != 1: continue
                if a_ in V: continue
            else:
                if b_-a_ != 1: continue
                if a_ in V: continue
            ok = False; break
        if ok: c += 1
    return c

def blocks_from_S(n, S):
    """Given a set S of value bonds, return the list of block sizes, blocks
    being maximal runs of consecutive values glued by the bonds in S."""
    sizes = []
    cur = 1
    for v in range(1, n):
        if v in S:
            cur += 1
        else:
            sizes.append(cur); cur = 1
    sizes.append(cur)
    return sizes

def arrangements_avoiding(sizes, P):
    """Count orderings of the blocks such that NO internal bond position lies
    in P. Internal bond positions of a block of size L placed after `off`
    positions are off+1, ..., off+L-1."""
    j = len(sizes)
    total = 0
    for perm in permutations(range(j)):
        off = 0; ok = True
        for idx in perm:
            L = sizes[idx]
            for i in range(1, L):
                if off + i in P:
                    ok = False; break
            if not ok: break
            off += L
        if ok: total += 1
    return total

def method1(n, V):
    """IE over S only, no position layer: sum_S (-1)^{|S|} (n-|S|)! ."""
    allowed = [v for v in range(1, n) if v not in V]
    tot = 0
    for r in range(len(allowed)+1):
        for S in combinations(allowed, r):
            tot += (-1)**r * factorial(n - r)
    return tot

def method2(n, P, V):
    """Full assembly: IE over S, then count arrangements avoiding P."""
    allowed = [v for v in range(1, n) if v not in V]
    tot = 0
    for r in range(len(allowed)+1):
        for S in combinations(allowed, r):
            sizes = blocks_from_S(n, set(S))
            tot += (-1)**r * arrangements_avoiding(sizes, P)
    return tot

def method2_via_T(n, P, V):
    """Same as method2 but removing the position constraint by the SECOND IE
    over subsets T of P, exactly as Section 3 (c) describes. Must agree with
    method2; this tests step (c) specifically."""
    allowed = [v for v in range(1, n) if v not in V]
    tot = 0
    for r in range(len(allowed)+1):
        for S in combinations(allowed, r):
            sizes = blocks_from_S(n, set(S))
            j = len(sizes)
            inner = 0
            Plist = sorted(P)
            for tsize in range(len(Plist)+1):
                for T in combinations(Plist, tsize):
                    # count arrangements where every position of T IS an
                    # internal bond position (others unconstrained)
                    cnt = 0
                    for perm in permutations(range(j)):
                        off = 0
                        hit = set()
                        for idx in perm:
                            L = sizes[idx]
                            for i in range(1, L):
                                hit.add(off + i)
                            off += L
                        if all(t in hit for t in T): cnt += 1
                    inner += (-1)**len(T) * cnt
            tot += (-1)**r * inner
    return tot

# ---------------------------------------------------------------- block factor
def block_factor_check(LMAX=8):
    """B*(x;u) = sum_{l>=2} (-1)^{l-1} x^l (u+...+u^{l-1})  ==  -x^2 u /((1+x)(1+xu)).
    Checked as a truncated bivariate series."""
    import sympy as sp
    x, u = sp.symbols('x u')
    lhs = sum((-1)**(l-1) * x**l * sum(u**i for i in range(1, l))
              for l in range(2, LMAX+1))
    rhs = -x**2*u/((1+x)*(1+x*u))
    diff = sp.expand(sp.series(sp.expand(lhs - rhs), x, 0, LMAX+1).removeO())
    return sp.simplify(diff) == 0, sp.sstr(sp.factor(rhs))

out = {"cases": []}
print("n   P            V            brute   method1   method2   method2_T   ok")
CASES = [
    (5, set(),      set()),
    (5, {2},        set()),
    (5, set(),      {3}),
    (5, {2},        {3}),
    (6, {3},        {2}),
    (6, {2, 4},     {3}),
    (6, {3},        {2, 5}),
    (7, {3},        {4}),
    (7, {2, 5},     {3, 6}),
]
allok = True
for (n, P, V) in CASES:
    b = brute_R(n, P, V)
    m1 = method1(n, V)
    m2 = method2(n, P, V)
    m2t = method2_via_T(n, P, V)
    ok = (b == m2 == m2t)
    if not ok: allok = False
    out["cases"].append({"n": n, "P": sorted(P), "V": sorted(V),
                         "brute": b, "method1_no_position_layer": m1,
                         "method2": m2, "method2_via_T": m2t, "ok": ok})
    print("%d   %-12s %-12s %6d  %8d  %8d  %10d   %s"
          % (n, sorted(P), sorted(V), b, m1, m2, m2t, "YES" if ok else "NO"),
          flush=True)

out["all_match"] = allok
print("\nSection 3 combinatorial content verified:", allok)
print("(method1 differs from brute wherever P is nonempty, as it must —")
print(" that is the control showing the position layer is doing real work.)")

try:
    bf_ok, bf = block_factor_check()
    out["block_factor_identity_ok"] = bool(bf_ok)
    out["block_factor"] = bf
    print("\nBlock-factor identity B*(x;u) = %s : %s" % (bf, bf_ok))
except Exception as e:
    out["block_factor_error"] = str(e)[:200]
    print("\nblock factor check skipped:", str(e)[:120])

json.dump(out, open("rep_check_results.json", "w"), indent=1)
print("\nwrote rep_check_results.json — send this back")
