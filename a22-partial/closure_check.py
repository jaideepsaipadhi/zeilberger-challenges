#!/usr/bin/env python3
"""
closure_check.py — gate the three-direction reduction before building it.

THE REDUCTION. e_X(m) = <k!j!> Res_t [ t^(-m-1) * 2 C_m ], where C_m = A_m B_m
is polynomial in (t,z,w) and satisfies the PROVED degree-6 recurrence with
coefficients constant in m. Eliminating x this way keeps everything
polynomial, so the umbral stays valid (this is what killed the branch-split
plan: alpha, beta were rational in z,w, outside the umbral's domain).

THE CLAIM IT RESTS ON. Telescoping the vector V_m = (C_m, ..., C_{m+5})
requires the module spanned by V_m over Q(t,z,w) to be CLOSED under d_t, d_z,
d_w — otherwise the certificates cannot be expressed in that vector and the
construction collapses.

WHY IT SHOULD HOLD: C_m is a combination of rho^m for rho in {r1^2, r1 r2,
r2^2} (each doubled), and d rho^m = m rho^(m-1) rho', which stays in the same
span with coefficients rational in t and polynomial in m.

WHY CHECK ANYWAY: "there is a reason to believe it" is exactly what was said
about the branch split, which then failed on a domain technicality. This is a
five-second test.

WHAT IS TESTED, at random rational (z,w) so everything is in Q(t):
  for each m, solve   d_v C_m  =  sum_{i=0..5} c_i(t) C_{m+i}
  over Q(t), for v = t, z, w, and report whether a solution exists.
  Also reports whether the c_i look polynomial in m (compared across m).

If any direction fails to close, the coupled reduction is not available and
we stop rather than building on it.

Run:  python3 closure_check.py       (sympy only, seconds)
Writes closure_check_results.json — SEND THIS BACK.
"""
import sympy as sp, json, random, time

t, z, w, v = sp.symbols('t z w v')
T0 = time.time()
out = {}

d1s = sp.expand(1 + t - z*t - w)
d2s = sp.expand(t*(1 - z - w))

def C_seq(nmax, d1, d2):
    Dv = 1 + d1*v + d2*v**2
    F1 = (1 + v)*(1 + v*t)/Dv
    F2 = -t*v**2*(1 + v)*(1 + v*t)/Dv**2
    s1 = sp.expand(sp.series(F1, v, 0, 6).removeO())
    s2 = sp.expand(sp.series(F2, v, 0, 6).removeO())
    A = [sp.cancel(s1.coeff(v, k)) for k in range(6)]
    B = [sp.cancel(s2.coeff(v, k)) for k in range(6)]
    cB = [None, sp.expand(2*d1), sp.expand(d1**2 + 2*d2),
          sp.expand(2*d1*d2), sp.expand(d2**2)]
    while len(A) <= nmax:
        k = len(A); A.append(sp.expand(-(d1*A[k-1] + d2*A[k-2])))
    while len(B) <= nmax:
        k = len(B)
        B.append(sp.expand(-(cB[1]*B[k-1] + cB[2]*B[k-2]
                             + cB[3]*B[k-3] + cB[4]*B[k-4])))
    return [sp.expand(A[k]*B[k]) for k in range(nmax+1)]

NM = 14
rnd = random.Random(23)
trials = []
for trial in range(2):
    zv = sp.Rational(rnd.randint(2, 9), rnd.randint(2, 7))
    wv = sp.Rational(rnd.randint(2, 9), rnd.randint(2, 7))
    print("\n(z,w) = (%s, %s)" % (zv, wv), flush=True)
    # C in the full (t,z,w) then substitute, so d_z, d_w are meaningful
    Cfull = C_seq(NM, d1s, d2s)
    rec = {"z": str(zv), "w": str(wv), "directions": {}}
    for dirname, sym in (("t", t), ("z", z), ("w", w)):
        results_dir = []
        for mm in (2, 3, 4, 5):
            target = sp.cancel(sp.diff(Cfull[mm], sym).subs({z: zv, w: wv}))
            basis = [sp.cancel(Cfull[mm+i].subs({z: zv, w: wv})) for i in range(6)]
            cs = sp.symbols('c0:6')
            expr = sp.expand(target - sum(cs[i]*basis[i] for i in range(6)))
            # solve over Q(t): collect powers of t
            polyexpr = sp.Poly(sp.expand(sp.numer(sp.together(expr))), t)
            eqs = polyexpr.coeffs()
            sol = sp.solve(eqs, cs, dict=True)
            ok = bool(sol)
            results_dir.append({"m": mm, "closes": ok,
                                "c": [sp.sstr(sp.cancel(sol[0].get(cs[i], 0)))[:80]
                                      for i in range(6)] if ok else None})
            print("   d_%s at m=%d: %s" % (dirname, mm,
                                           "closes" if ok else "DOES NOT CLOSE"),
                  flush=True)
        rec["directions"][dirname] = results_dir
    trials.append(rec)

out["trials"] = trials
allok = all(d["closes"] for tr in trials for dl in tr["directions"].values()
            for d in dl)
out["module_closed_under_all_three"] = allok
print("\nMODULE CLOSED UNDER d_t, d_z, d_w:", allok, "(%.0fs)" % (time.time()-T0))
if allok:
    print("=> the coupled three-direction reduction is available:")
    print("   V_m = (C_m..C_{m+5}), Omega = 1/t, certificates are 6-component,")
    print("   Theta = t^b only (C_m is polynomial, so no Q denominators),")
    print("   and only z,w carry boundary kills.")
else:
    print("=> NOT available. Do not build the coupled formulation.")
out["next"] = ("If closed: build the coupled telescoping with V_m, the "
               "companion matrix for the shift, and the three direction "
               "matrices measured here. If not: stop and reassess.")
json.dump(out, open("closure_check_results.json", "w"), indent=1)
print("wrote closure_check_results.json — send this back")
