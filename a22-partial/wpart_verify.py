#!/usr/bin/env python3
"""
wpart_verify.py — the gate my last script got wrong.

WHAT HAPPENED. wpart_gate.sage asked ore_algebra for a telescoper and then
demanded it equal our proved HOMOGENEOUS order-2 operator. It never could:
ct() returns an operator together with the BOUNDARY term, i.e. an
INHOMOGENEOUS relation, because the certificate does not vanish at z = 0.
(Our own solver imposes N(0) = 0 by hand to kill that; ore_algebra instead
hands the boundary back.) Same situation as the Hertzsprung control, where
ct gave an order-1 differential operator plus c(x) and we homogenised to
reach the order-7 recurrence.

Reading the even output by hand already shows it worked:
    L = (4m^2+2m-1) Sm - (16m^4+32m^3-10m)
    L applied to e_W:  m=1: -6,  m=2: -10,  m=3: -14   =  -2(2m+1)
so ore_algebra returned the valid inhomogeneous order-1 relation
    (4m^2+2m-1) e_W(m+1) - (16m^4+32m^3-10m) e_W(m) = -2(2m+1).

THIS SCRIPT does the check properly, for every operator the four runs
returned:
  1. apply it to the real data;
  2. fit the residual as a polynomial in m (that is the inhomogeneity);
  3. homogenise -- if L f = g with g a nonzero polynomial, then
     (g(m+1) S - g(m)) L  annihilates f -- and compare the result with the
     proved order-2 operator up to a rational factor;
  4. report which sign convention gives the cleaner inhomogeneity.

PROVED TARGETS (exact, certificates verified symbolically):
  even:  -m(2m-1)(2m+3),  (8m^3+32m^2+32m+7)/2,  -(2m+1)/2
  odd:   -2(m-1)(m+1)(2m-1),  4m^3+10m^2+3m-1,  -m

Run:  python3 wpart_verify.py       (sympy only, seconds)
"""
import sympy as sp, json

m, z = sp.symbols('m z')

def umb(expr):
    p = sp.Poly(sp.expand(expr), z)
    return sum(c*sp.factorial(e[0]) for e, c in zip(p.monoms(), p.coeffs()))

# ---------------------------------------------------------------- data
eW = {k: sp.Integer(umb(z**2*(z-1)**(2*k-2))) for k in range(1, 16)}
oW = {k: sp.Integer(umb(z**2*(z-1)**(2*k-3))) for k in range(2, 16)}
print("e_W(1..4) =", [int(eW[k]) for k in range(1, 5)])
print("o_W(2..5) =", [int(oW[k]) for k in range(2, 6)])

TARGET_EVEN = [sp.expand(-m*(2*m-1)*(2*m+3)),
               sp.expand((8*m**3+32*m**2+32*m+7)/2),
               sp.expand(-(2*m+1)/2)]
TARGET_ODD  = [sp.expand(-2*(m-1)*(m+1)*(2*m-1)),
               sp.expand(4*m**3+10*m**2+3*m-1),
               sp.expand(-m)]

# operators returned by ore_algebra, as [coeff of Sm^0, coeff of Sm^1]
RUNS = {
 "even MINUS": (eW, 1, [sp.expand(-(16*m**4+32*m**3-10*m)), sp.expand(4*m**2+2*m-1)], TARGET_EVEN),
 "even PLUS":  (eW, 1, [sp.expand(-(16*m**4-2*m)),          sp.expand(4*m**2-6*m+3)], TARGET_EVEN),
 "odd MINUS":  (oW, 2, [sp.expand(-(16*m**4-24*m**2+6*m+2)), sp.expand(4*m**2-2*m-1)], TARGET_ODD),
 "odd PLUS":   (oW, 2, [sp.expand(16*m**4-32*m**3+24*m**2-10*m+2), sp.expand(-4*m**2+10*m-7)], TARGET_ODD),
}

def apply_op(coeffs, data, m0, nmax=10):
    """coeffs[i] multiplies f(m+i). Returns [(m, value)]."""
    vals = []
    for mm in range(m0, m0+nmax):
        if any((mm+i) not in data for i in range(len(coeffs))):
            break
        s = sum(sp.Integer(coeffs[i].subs(m, mm)) * data[mm+i]
                for i in range(len(coeffs)))
        vals.append((mm, sp.nsimplify(s)))
    return vals

def fit_poly(vals, maxdeg=4):
    """Fit residual as a polynomial in m; return it or None."""
    for d in range(maxdeg+1):
        cs = sp.symbols('c0:%d' % (d+1))
        poly = sum(cs[i]*m**i for i in range(d+1))
        eqs = [sp.Eq(poly.subs(m, mm), v) for mm, v in vals[:d+1]]
        sol = sp.solve(eqs, cs, dict=True)
        if not sol:
            continue
        cand = sp.expand(poly.subs(sol[0]))
        if all(sp.simplify(cand.subs(m, mm) - v) == 0 for mm, v in vals):
            return sp.simplify(cand)
    return None

out = {}
for tag, (data, m0, coeffs, target) in RUNS.items():
    print("\n" + "="*60)
    print(tag)
    vals = apply_op(coeffs, data, m0)
    print("   L applied to data:", [(mm, int(v)) for mm, v in vals[:6]])
    g = fit_poly(vals)
    rec = {"operator": [sp.sstr(c) for c in coeffs],
           "applied": [(mm, str(v)) for mm, v in vals[:8]]}
    if g is None:
        print("   residual is NOT a low-degree polynomial in m -> suspect")
        rec["inhomogeneity"] = None
        out[tag] = rec
        continue
    print("   inhomogeneity g(m) =", sp.factor(g))
    rec["inhomogeneity"] = sp.sstr(sp.factor(g))

    # homogenise: (g(m+1) S - g(m)) o L
    gm, gm1 = g, g.subs(m, m+1)
    L = coeffs                      # L = c0 + c1 S
    # (g(m+1) S - g(m)) (c0 + c1 S) = -g(m)c0 + (g(m+1)c0(m+1) - g(m)c1) S
    #                                  + g(m+1) c1(m+1) S^2
    H = [sp.expand(-gm*L[0]),
         sp.expand(gm1*L[0].subs(m, m+1) - gm*L[1]),
         sp.expand(gm1*L[1].subs(m, m+1))]
    print("   homogenised order-2 operator obtained")
    rec["homogenised"] = [sp.sstr(sp.factor(h)) for h in H]

    # compare with the proved target up to a rational factor in m
    ratio = None; ok = True
    for a_, b_ in zip(H, target):
        if b_ == 0:
            if a_ != 0: ok = False; break
            continue
        r_ = sp.cancel(a_/b_)
        if ratio is None: ratio = r_
        elif sp.simplify(r_ - ratio) != 0: ok = False; break
    rec["matches_proved"] = bool(ok)
    rec["scalar"] = sp.sstr(sp.factor(ratio)) if ratio is not None else None
    print("   matches proved operator:", ok,
          ("(factor %s)" % sp.factor(ratio)) if ok and ratio is not None else "")
    out[tag] = rec

print("\n" + "="*60)
good = [k for k, v in out.items() if v.get("matches_proved")]
print("runs reproducing the proved operator:", good if good else "NONE")
out["verdict"] = good
json.dump(out, open("wpart_verify_results.json", "w"), indent=1)
print("wrote wpart_verify_results.json")
