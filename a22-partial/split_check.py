#!/usr/bin/env python3
"""
split_check.py — verify the claimed decomposition BEFORE building anything on it.

THE CLAIM (derived, not yet checked). C_m = A_m B_m is C-finite in m with
characteristic roots {r1^2, r1 r2, r2^2}, each of multiplicity 2, where r1, r2
are the roots of Y^2 + d1 Y + d2 = 0, d1 = 1+t-zt-w, d2 = t(1-z-w). Since
    r1 r2 = d2      (RATIONAL in t,z,w)
    r1^2 + r2^2 = d1^2 - 2 d2,   r1^2 * r2^2 = d2^2   (both RATIONAL)
the sequence should split canonically as

    C_m = [ (alpha + beta*m) * d2^m ]            <- RATIONAL branch
        + [ conjugate-pair branch on {r1^2, r2^2} ]   <- 2-component

with alpha, beta rational in (t,z,w). Two consequences were claimed:
  (i)  the rational branch needs only a scalar treatment, no field extension;
  (ii) Omega for that branch is d2/t = 1-z-w, which is FREE OF t, so its
       t-residue is a plain coefficient read-off, possibly closed-form.

WHY CHECK FIRST. This is four layers of restructuring deep, all derived in
conversation, none of it machine-checked. Everything downstream depends on it.

WHAT THIS SCRIPT DOES (all exact, at random rational (z,w) so it is fast):
  1. computes C_m directly for m = 1..14 from the verified A, B recurrences;
  2. extracts the rational branch by the standard projection: the part of C_m
     annihilated by (S - d2)^2, obtained as the residue-style component, and
     fits (alpha + beta m) d2^m to the first few terms;
  3. subtracts it and checks the REMAINDER is annihilated by the degree-4
     operator with characteristic polynomial (Y^2 - (d1^2-2d2) Y + d2^2)^2
     -- i.e. that the remainder really is the conjugate-pair branch;
  4. reports whether alpha, beta are rational functions (they must be) and
     prints them;
  5. checks claim (ii): that d2/t simplifies to 1-z-w exactly.

If step 3 fails, the split is wrong and nothing should be built on it.

Run:  python3 split_check.py        (sympy only, seconds)
Writes split_check_results.json — SEND THIS BACK.
"""
import sympy as sp, json, random, time

t, z, w, m, v, Y = sp.symbols('t z w m v Y')
T0 = time.time()
out = {}

d1s = sp.expand(1 + t - z*t - w)
d2s = sp.expand(t*(1 - z - w))

# ---------------------------------------------------------------- claim (ii)
omega_rat = sp.simplify(d2s/t)
out["d2_over_t"] = sp.sstr(omega_rat)
out["d2_over_t_is_1_minus_z_minus_w"] = bool(sp.simplify(omega_rat - (1 - z - w)) == 0)
print("claim (ii): d2/t =", sp.sstr(omega_rat),
      "-> t-free:", out["d2_over_t_is_1_minus_z_minus_w"], flush=True)

def C_seq(nmax, d1, d2):
    """C_m = A_m B_m from the verified recurrences (short series start)."""
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
        k = len(A); A.append(sp.cancel(-(d1*A[k-1] + d2*A[k-2])))
    while len(B) <= nmax:
        k = len(B)
        B.append(sp.cancel(-(cB[1]*B[k-1] + cB[2]*B[k-2]
                             + cB[3]*B[k-3] + cB[4]*B[k-4])))
    return [sp.cancel(A[k]*B[k]) for k in range(nmax+1)]

NM = 14
rnd = random.Random(11)
trials = []
for trial in range(3):
    zv = sp.Rational(rnd.randint(2, 9), rnd.randint(2, 7))
    wv = sp.Rational(rnd.randint(2, 9), rnd.randint(2, 7))
    d1 = sp.expand(d1s.subs({z: zv, w: wv}))
    d2 = sp.expand(d2s.subs({z: zv, w: wv}))
    C = C_seq(NM, d1, d2)
    print("\n(z,w) = (%s, %s)" % (zv, wv), flush=True)

    # ---- fit alpha, beta in the rational branch (alpha + beta m) d2^m
    al, be = sp.symbols('al be')
    # project: the rational branch is the unique (al+be*m) d2^m such that
    # C_m - that is annihilated by (S^2 - (d1^2-2d2) S + d2^2)^2.
    # Solve for al, be by imposing that annihilation on two windows.
    e1c = sp.expand(d1**2 - 2*d2)          # r1^2 + r2^2
    e2c = sp.expand(d2**2)                 # r1^2 * r2^2
    # conjugate-branch operator, squared: (S^2 - e1c S + e2c)^2
    S = sp.Symbol('S')
    op = sp.expand((S**2 - e1c*S + e2c)**2)
    opc = [sp.expand(sp.Poly(op, S).coeff_monomial(S**k)) for k in range(5)]

    def rat(mm):  return (al + be*mm) * d2**mm
    eqs = []
    for start in (3, 4):
        expr = 0
        for k in range(5):
            expr += opc[k] * (C[start+k] - rat(start+k))
        eqs.append(sp.expand(sp.cancel(expr)))
    sol = sp.solve(eqs, [al, be], dict=True)
    rec = {"z": str(zv), "w": str(wv)}
    if not sol:
        rec["solved"] = False
        print("   could not solve for alpha,beta", flush=True)
        trials.append(rec); continue
    s0 = sol[0]
    alv, bev = sp.cancel(s0[al]), sp.cancel(s0[be])
    rec["solved"] = True
    rec["alpha"] = sp.sstr(alv); rec["beta"] = sp.sstr(bev)
    print("   alpha =", sp.sstr(alv)[:120], flush=True)
    print("   beta  =", sp.sstr(bev)[:120], flush=True)

    # ---- check: remainder annihilated by the conjugate-pair operator
    bad = []
    for start in range(1, NM-4):
        expr = 0
        for k in range(5):
            expr += opc[k] * (C[start+k] - (alv + bev*(start+k))*d2**(start+k))
        if sp.cancel(sp.together(sp.expand(expr))) != 0:
            bad.append(start)
    rec["remainder_is_conjugate_branch"] = (not bad)
    rec["failures"] = bad[:6]
    print("   remainder annihilated by (S^2-e1c S+e2c)^2:",
          "YES" if not bad else "NO at %s" % bad[:5], flush=True)

    # ---- check the rational branch alone satisfies (S-d2)^2
    bad2 = []
    for start in range(1, NM-2):
        r_ = lambda mm: (alv + bev*mm)*d2**mm
        val = sp.cancel(r_(start+2) - 2*d2*r_(start+1) + d2**2*r_(start))
        if sp.simplify(val) != 0: bad2.append(start)
    rec["rational_branch_satisfies_(S-d2)^2"] = (not bad2)
    print("   rational branch killed by (S-d2)^2:",
          "YES" if not bad2 else "NO at %s" % bad2[:5], flush=True)
    trials.append(rec)

out["trials"] = trials
out["split_verified"] = all(tr.get("remainder_is_conjugate_branch") and
                            tr.get("rational_branch_satisfies_(S-d2)^2")
                            for tr in trials if tr.get("solved"))
print("\nSPLIT VERIFIED:", out["split_verified"], "(%.0fs)" % (time.time()-T0))
if not out["split_verified"]:
    print("The decomposition is NOT as claimed — do not build the coupled")
    print("formulation on it. Report this and we rethink.")
out["next"] = ("If verified: the conjugate branch is a 2-component coupled "
               "telescoping with Omega tied to {r1^2, r2^2} (rational "
               "symmetric functions e1c = d1^2-2d2, e2c = d2^2), and the "
               "rational branch has Omega = d2/t = 1-z-w, t-free, so its "
               "t-extraction is a coefficient read-off.")
json.dump(out, open("split_check_results.json", "w"), indent=1)
print("wrote split_check_results.json — send this back")
