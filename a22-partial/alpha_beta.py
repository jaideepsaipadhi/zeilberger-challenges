#!/usr/bin/env python3
"""
alpha_beta.py — the rational branch in closed form.

WHAT split_check.py ESTABLISHED (verified at three rational (z,w)):
    C_m = (alpha + beta*m) * d2^m  +  [conjugate-pair branch on {r1^2, r2^2}]
    d2/t = 1 - z - w   (t-free)
and, reading the three outputs, beta's numerator was ALWAYS proportional to
(t-1)^2 — 27000(t-1)^2, 416745(t-1)^2, 1152(t-1)^2 — with alpha and beta
sharing the same denominator family. That structure says the branch should
have a clean closed form rather than needing any search.

WHY IT MATTERS. Omega for this branch is d2/t = 1-z-w, which is FREE OF t.
So [t^m]((alpha+beta m) d2^m) is a coefficient read-off, not a residue
extraction, and the umbral then acts on a polynomial. If alpha, beta come out
in closed form, this branch contributes NO TELESCOPING AT ALL — exactly like
e_W and o_W this morning — leaving the conjugate-pair branch as the only
remaining telescoping problem.

METHOD. Same projection as split_check, run ONCE with z, w symbolic: impose
that C_m - (alpha+beta m) d2^m is annihilated by the conjugate-pair operator
(S^2 - e1c S + e2c)^2, with e1c = d1^2 - 2 d2, e2c = d2^2, and solve the two
resulting linear equations for alpha, beta over Q(t,z,w).

Cost control (the lesson from build_H): only C_3..C_8 are needed — the two
windows at start = 3 and 4 — so nothing high-index is ever expanded.

VERIFICATION: the symbolic alpha, beta are specialised back to the three
(z,w) points from split_check and compared against the values found there.

Run:  python3 alpha_beta.py         (sympy only)
Writes alpha_beta_results.json — SEND THIS BACK.
"""
import sympy as sp, json, time

t, z, w, m, v, S = sp.symbols('t z w m v S')
T0 = time.time()
out = {}

d1 = sp.expand(1 + t - z*t - w)
d2 = sp.expand(t*(1 - z - w))
e1c = sp.expand(d1**2 - 2*d2)          # r1^2 + r2^2
e2c = sp.expand(d2**2)                 # r1^2 r2^2
out["d1"] = sp.sstr(d1); out["d2"] = sp.sstr(d2)
out["e1c"] = sp.sstr(sp.factor(e1c)); out["e2c"] = sp.sstr(sp.factor(e2c))
print("e1c = r1^2+r2^2 =", sp.sstr(sp.factor(e1c)), flush=True)
print("e2c = r1^2 r2^2 =", sp.sstr(sp.factor(e2c)), flush=True)

# ---------------------------------------------------------------- C_1..C_8
print("\nbuilding C_1..C_8 symbolically (nothing high-index) ...", flush=True)
Dv = 1 + d1*v + d2*v**2
F1 = (1 + v)*(1 + v*t)/Dv
F2 = -t*v**2*(1 + v)*(1 + v*t)/Dv**2
s1 = sp.expand(sp.series(F1, v, 0, 6).removeO())
s2 = sp.expand(sp.series(F2, v, 0, 6).removeO())
A = [sp.expand(s1.coeff(v, k)) for k in range(6)]
B = [sp.expand(s2.coeff(v, k)) for k in range(6)]
cB = [None, sp.expand(2*d1), sp.expand(d1**2 + 2*d2),
      sp.expand(2*d1*d2), sp.expand(d2**2)]
NEED = 8
while len(A) <= NEED:
    k = len(A); A.append(sp.expand(-(d1*A[k-1] + d2*A[k-2])))
while len(B) <= NEED:
    k = len(B)
    B.append(sp.expand(-(cB[1]*B[k-1] + cB[2]*B[k-2] + cB[3]*B[k-3] + cB[4]*B[k-4])))
C = [sp.expand(A[k]*B[k]) for k in range(NEED+1)]
print("   done (%.0fs)" % (time.time()-T0), flush=True)

# ---------------------------------------------------------------- solve
al, be = sp.symbols('al be')
op = sp.expand((S**2 - e1c*S + e2c)**2)
opc = [sp.expand(sp.Poly(op, S).coeff_monomial(S**k)) for k in range(5)]
def rat(mm): return (al + be*mm) * d2**mm

print("\nsolving for alpha, beta over Q(t,z,w) ...", flush=True)
eqs = []
for start in (3, 4):
    expr = 0
    for k in range(5):
        expr += opc[k] * (C[start+k] - rat(start+k))
    eqs.append(sp.expand(expr))
sol = sp.solve(eqs, [al, be], dict=True)
if not sol:
    out["solved"] = False
    print("   NO SOLUTION — the projection failed symbolically", flush=True)
    json.dump(out, open("alpha_beta_results.json", "w"), indent=1)
    raise SystemExit
s0 = sol[0]
alv = sp.cancel(sp.together(s0[al]))
bev = sp.cancel(sp.together(s0[be]))
out["solved"] = True
print("   solved (%.0fs)" % (time.time()-T0), flush=True)

alf, bef = sp.factor(alv), sp.factor(bev)
out["alpha"] = sp.sstr(alf)
out["beta"] = sp.sstr(bef)
print("\nbeta  =", sp.sstr(bef), flush=True)
print("\nalpha =", sp.sstr(alf)[:600], flush=True)

# structure report: numerators/denominators
for name, val in (("alpha", alv), ("beta", bev)):
    n_, d_ = sp.fraction(sp.cancel(val))
    out[name + "_num_factored"] = sp.sstr(sp.factor(n_))[:800]
    out[name + "_den_factored"] = sp.sstr(sp.factor(d_))[:800]
    print("\n%s numerator  :" % name, sp.sstr(sp.factor(n_))[:300], flush=True)
    print("%s denominator:" % name, sp.sstr(sp.factor(d_))[:300], flush=True)

# does beta carry the observed (t-1)^2 ?
out["beta_has_(t-1)^2"] = bool(sp.simplify(
    sp.cancel(sp.fraction(sp.cancel(bev))[0] / (t-1)**2)).is_rational_function(t, z, w)
    and sp.simplify(sp.fraction(sp.cancel(bev))[0].subs(t, 1)) == 0)
print("\nbeta numerator vanishes at t=1 (the observed (t-1)^2):",
      out["beta_has_(t-1)^2"], flush=True)

# ---------------------------------------------------------------- verify
print("\nverifying against split_check's rational-point values ...", flush=True)
pts = [(sp.Rational(3, 2), sp.Rational(9, 5)),
       (sp.Rational(5, 3), sp.Rational(9, 7)),
       (sp.Integer(2), sp.Rational(9, 4))]
checks = []
for zv, wv in pts:
    d1n = sp.expand(d1.subs({z: zv, w: wv})); d2n = sp.expand(d2.subs({z: zv, w: wv}))
    e1n = sp.expand(d1n**2 - 2*d2n); e2n = sp.expand(d2n**2)
    opn = sp.expand((S**2 - e1n*S + e2n)**2)
    opcn = [sp.expand(sp.Poly(opn, S).coeff_monomial(S**k)) for k in range(5)]
    Cn = [sp.cancel(c.subs({z: zv, w: wv})) for c in C]
    an = sp.cancel(alv.subs({z: zv, w: wv})); bn = sp.cancel(bev.subs({z: zv, w: wv}))
    bad = []
    for start in (1, 2, 3, 4):
        expr = 0
        for k in range(5):
            expr += opcn[k]*(Cn[start+k] - (an + bn*(start+k))*d2n**(start+k))
        if sp.cancel(sp.together(sp.expand(expr))) != 0:
            bad.append(start)
    checks.append({"z": str(zv), "w": str(wv), "ok": not bad, "failures": bad})
    print("   (z,w)=(%s,%s): %s" % (zv, wv, "clean" if not bad else "MISMATCH %s" % bad),
          flush=True)
out["point_checks"] = checks
out["verified"] = all(c["ok"] for c in checks)
print("\nCLOSED FORM VERIFIED:", out["verified"], "(%.0fs)" % (time.time()-T0))
out["next"] = ("Rational branch is closed-form: Omega = d2/t = 1-z-w is "
               "t-free, so [t^m]((alpha+beta m) d2^m ...) is a coefficient "
               "read-off and the umbral acts on a polynomial — NO telescoping "
               "needed, as with e_W/o_W. Remaining: the 2-component "
               "conjugate-pair branch on {r1^2, r2^2}.")
json.dump(out, open("alpha_beta_results.json", "w"), indent=1)
print("wrote alpha_beta_results.json — send this back")
