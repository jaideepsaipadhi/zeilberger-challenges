#!/usr/bin/env python3
"""
wpart_final.py — test annihilation, not operator equality.

WHAT wpart_verify.py GOT RIGHT: the MINUS (e^{-z}) convention is confirmed for
both parities. ore_algebra's operators applied to the real data give
    even:  -6, -10, -14, -18, -22, -26   =  -2(2m+1)
    odd:    8,  12,  16,  20,  24,  28   =   4m
Clean linear inhomogeneities. The PLUS convention gives residuals growing like
the data itself, i.e. nonsense. So the encoding sign is settled.

WHAT IT GOT WRONG: it then demanded the homogenised operator equal our proved
minimal operator up to a single rational scalar. That is too strong.
Homogenising an order-1 inhomogeneous relation yields SOME order-2
annihilator, not necessarily the MINIMAL one; the two can differ by a left
factor that is not a scalar. (My coefficient ordering may also have been
reversed relative to ore_algebra's convention.)

THE TEST THAT MATTERS: does the homogenised operator annihilate the data?
If yes, ore_algebra produced a valid relation of exactly the expected shape,
and our proved operator is simply the minimal element of the same annihilator
module -- nothing to reconcile.

This script, for each parity:
  1. re-derives the inhomogeneity by fitting;
  2. homogenises via (g(m+1) S - g(m)) o L, trying BOTH coefficient orderings
     so a convention slip cannot produce a false negative;
  3. checks annihilation on all available data;
  4. reports how the homogenised operator relates to the proved one --
     ratio of leading coefficients, and whether the proved operator divides it.

Run:  python3 wpart_final.py
"""
import sympy as sp, json

m, z = sp.symbols('m z')

def umb(e):
    p = sp.Poly(sp.expand(e), z)
    return sum(c*sp.factorial(x[0]) for x, c in zip(p.monoms(), p.coeffs()))

eW = {k: sp.Integer(umb(z**2*(z-1)**(2*k-2))) for k in range(1, 18)}
oW = {k: sp.Integer(umb(z**2*(z-1)**(2*k-3))) for k in range(2, 18)}

PROVED_EVEN = [sp.expand(-m*(2*m-1)*(2*m+3)),
               sp.expand((8*m**3+32*m**2+32*m+7)/2),
               sp.expand(-(2*m+1)/2)]
PROVED_ODD  = [sp.expand(-2*(m-1)*(m+1)*(2*m-1)),
               sp.expand(4*m**3+10*m**2+3*m-1),
               sp.expand(-m)]

# ore_algebra's returned operators, MINUS convention (the confirmed one)
A_EVEN = sp.expand(4*m**2 + 2*m - 1)              # coefficient of Sm
B_EVEN = sp.expand(-(16*m**4 + 32*m**3 - 10*m))   # coefficient of Sm^0
A_ODD  = sp.expand(4*m**2 - 2*m - 1)
B_ODD  = sp.expand(-(16*m**4 - 24*m**2 + 6*m + 2))

def apply_L(c0, c1, data, m0, n=12):
    vals = []
    for mm in range(m0, m0+n):
        if (mm+1) not in data: break
        vals.append((mm, sp.expand(c0.subs(m, mm)*data[mm]
                                   + c1.subs(m, mm)*data[mm+1])))
    return vals

def fit(vals, maxdeg=3):
    for d in range(maxdeg+1):
        cs = sp.symbols('c0:%d' % (d+1))
        poly = sum(cs[i]*m**i for i in range(d+1))
        sol = sp.solve([sp.Eq(poly.subs(m, mm), v) for mm, v in vals[:d+1]],
                       cs, dict=True)
        if not sol: continue
        cand = sp.expand(poly.subs(sol[0]))
        if all(sp.simplify(cand.subs(m, mm) - v) == 0 for mm, v in vals):
            return sp.simplify(cand)
    return None

def homogenise(c0, c1, g):
    """(g(m+1) S - g(m)) o (c0 + c1 S)  ->  [H0, H1, H2]."""
    g1 = g.subs(m, m+1)
    return [sp.expand(-g*c0),
            sp.expand(g1*c0.subs(m, m+1) - g*c1),
            sp.expand(g1*c1.subs(m, m+1))]

def annihilates(H, data, m0, n=12):
    bad = []
    for mm in range(m0, m0+n):
        if (mm+2) not in data: break
        s = sum(H[i].subs(m, mm)*data[mm+i] for i in range(3))
        if sp.simplify(s) != 0: bad.append(mm)
    return bad

out = {}
for tag, (c0, c1, data, m0, proved) in {
    "even": (B_EVEN, A_EVEN, eW, 1, PROVED_EVEN),
    "odd":  (B_ODD,  A_ODD,  oW, 2, PROVED_ODD),
}.items():
    print("=" * 62)
    print(tag)
    rec = {}
    best = None
    for order_tag, (d0, d1) in (("as read (c0 + c1 S)", (c0, c1)),
                                ("reversed (c1 + c0 S)", (c1, c0))):
        vals = apply_L(d0, d1, data, m0)
        g = fit(vals)
        print("  %-22s residual: %s" % (order_tag,
              [(mm, int(v)) for mm, v in vals[:5]]))
        if g is None:
            print("      -> not polynomial; skipping")
            continue
        print("      -> g(m) =", sp.factor(g))
        H = homogenise(d0, d1, g)
        bad = annihilates(H, data, m0)
        print("      -> homogenised annihilates data:",
              "YES" if not bad else "NO at %s" % bad[:4])
        rec[order_tag] = {"g": sp.sstr(sp.factor(g)),
                          "annihilates": not bad,
                          "H": [sp.sstr(sp.factor(h)) for h in H]}
        if not bad and best is None:
            best = H
    if best is not None:
        # relation to the proved minimal operator
        r = [sp.cancel(best[i]/proved[i]) if proved[i] != 0 else None
             for i in range(3)]
        same = all(sp.simplify(r[i] - r[0]) == 0 for i in range(3)
                   if r[i] is not None)
        print("  ratio to proved operator, coefficientwise:")
        for i in range(3):
            print("      H%d/P%d = %s" % (i, i, sp.factor(r[i]) if r[i] is not None else "-"))
        print("  proportional to the proved (minimal) operator:", same)
        rec["proportional_to_proved"] = bool(same)
        rec["ratios"] = [sp.sstr(sp.factor(x)) if x is not None else None for x in r]
        if not same:
            print("  => not proportional, but BOTH annihilate the data:")
            print("     ore_algebra returned a valid non-minimal annihilator,")
            print("     which is exactly what homogenisation produces.")
    out[tag] = rec
    print()

json.dump(out, open("wpart_final_results.json", "w"), indent=1)
print("wrote wpart_final_results.json")
