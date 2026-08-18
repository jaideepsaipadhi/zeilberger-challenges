#!/usr/bin/env python3
"""
factor_recs.py — the factor recurrences, derived exactly and verified.

STRUCTURE BEING EXPLOITED. In the grand GF's X-part
    X = F2(s;t)F1(sig;t) + F1(s;t)F2(sig;t),
    F1(v;t) = (1+v)(1+v t)/D(v;t),
    F2(v;t) = -v^2 t (1+v)(1+v t)/D(v;t)^2,
    D(v;t) = 1 + (1+t-z t-w) v + t(1-z-w) v^2
the variable sig appears ONLY inside F1(sig;t) or F2(sig;t) — a
single-variable rational function of sig with (t,z,w) as parameters. D is
QUADRATIC in v, so the coefficient sequences

    A_m = [v^m] F1(v;t),        B_m = [v^m] F2(v;t)

satisfy linear recurrences in m with coefficients in Q(t,z,w) that follow
directly from the denominator — no telescoping search required. Writing
D = 1 + d1 v + d2 v^2  (d1 = 1+t-zt-w, d2 = t(1-z-w)) and N1 = (1+v)(1+vt)
= 1 + (1+t)v + t v^2:

    A_m + d1 A_{m-1} + d2 A_{m-2} = [v^m] N1      (zero for m >= 3)
so A is CONSTANT-COEFFICIENT of order 2 for m >= 3. For F2, denominator D^2
= 1 + 2d1 v + (d1^2+2d2) v^2 + 2 d1 d2 v^3 + d2^2 v^4, numerator
N2 = -t v^2 (1+v)(1+vt), giving a constant-coefficient order-4 recurrence
for m >= 5.

WHAT THIS SCRIPT DOES (all exact, no sampling):
  1. builds the recurrences symbolically from the denominators;
  2. verifies each against DIRECT SERIES EXPANSION of F1, F2 in v to high
     order — the coefficients are compared as rational functions of
     (t,z,w), i.e. exactly, not at sample points;
  3. reports the exact starting index from which each recurrence is
     homogeneous (the inhomogeneity dies once the numerator runs out);
  4. records the operators in a form the next step (symmetric product for
     the Hadamard-type diagonal) can consume.

WHY IT MATTERS: if these hold, the sigma- and s-extractions are FREE, and
the entire remaining difficulty is the symmetric-product combination plus
the [t^m] extraction and the two umbral directions — which is exactly where
the measured order-6 complexity should come from.

Run:  python3 factor_recs.py      (plain python3; sympy only)
Writes factor_recs_results.json — SEND THIS BACK. Should be quick.
"""
import sympy as sp, json

v, t, z, w, m = sp.symbols('v t z w m')

d1 = 1 + t - z*t - w
d2 = t*(1 - z - w)
D = 1 + d1*v + d2*v**2
N1 = (1 + v)*(1 + v*t)
N2 = -t*v**2*(1 + v)*(1 + v*t)
F1 = N1/D
F2 = N2/D**2

results = {"d1": sp.sstr(sp.expand(d1)), "d2": sp.sstr(sp.expand(d2))}

def series_coeffs(expr, upto):
    """Exact [v^k] for k = 0..upto, as rational functions of (t,z,w)."""
    ser = sp.series(expr, v, 0, upto + 1).removeO()
    ser = sp.expand(ser)
    return [sp.cancel(sp.together(ser.coeff(v, k))) for k in range(upto + 1)]

def check_recurrence(name, expr, den_poly, num_poly, order, upto=14):
    """den_poly = list of coefficients [c0, c1, ...] of the denominator in v.
    The claim: sum_j c_j * X_{k-j} = [v^k] num_poly  for all k (X_i = 0, i<0),
    hence homogeneous once k exceeds deg(num)."""
    out = {"target": name, "order": order}
    X = series_coeffs(expr, upto)
    nump = sp.Poly(sp.expand(num_poly), v)
    numc = {e[0]: c for e, c in zip(nump.monoms(), nump.coeffs())}
    degnum = nump.degree()
    out["numerator_degree"] = int(degnum)
    out["homogeneous_from_index"] = int(degnum) + 1
    bad = []
    for k in range(len(X)):
        lhs = 0
        for j, cj in enumerate(den_poly):
            if k - j >= 0:
                lhs += cj * X[k - j]
        rhs = numc.get(k, 0)
        if sp.simplify(sp.cancel(sp.together(lhs - rhs))) != 0:
            bad.append(k)
    out["identity_holds_for_all_k_tested"] = (not bad)
    out["failures"] = bad[:6]
    out["tested_up_to_index"] = upto
    out["denominator_coeffs"] = [sp.sstr(sp.expand(c)) for c in den_poly]
    # homogeneous operator, written as sum_j c_j X_{k-j} = 0 for k > degnum
    out["homogeneous_operator"] = (
        " + ".join("(%s)*X_{k-%d}" % (sp.sstr(sp.factor(c)), j)
                   for j, c in enumerate(den_poly) if c != 0) + " = 0")
    return out, X

# ---------------------------------------------------------------- F1
print("F1 = (1+v)(1+vt)/D,  D quadratic in v")
Dc = [sp.Integer(1), d1, d2]
r1, X1 = check_recurrence("A_m = [v^m]F1", F1, Dc, N1, 2)
print("   identity holds:", r1["identity_holds_for_all_k_tested"],
      "| homogeneous from m =", r1["homogeneous_from_index"])
print("   ", r1["homogeneous_operator"])
results["F1"] = r1

# ---------------------------------------------------------------- F2
print("\nF2 = -t v^2 (1+v)(1+vt)/D^2,  D^2 quartic in v")
D2poly = sp.Poly(sp.expand(D**2), v)
D2c = [sp.expand(D2poly.coeff_monomial(v**k)) for k in range(5)]
r2, X2 = check_recurrence("B_m = [v^m]F2", F2, D2c, N2, 4)
print("   identity holds:", r2["identity_holds_for_all_k_tested"],
      "| homogeneous from m =", r2["homogeneous_from_index"])
print("   ", r2["homogeneous_operator"])
results["F2"] = r2

# ------------------------------------------------- sanity: first few coeffs
results["A_first"] = [sp.sstr(sp.factor(c)) for c in X1[:4]]
results["B_first"] = [sp.sstr(sp.factor(c)) for c in X2[:5]]
print("\nA_0..A_3 =", results["A_first"])
print("B_0..B_4 =", results["B_first"])

# ------------------------------------------------- what the next step needs
results["next_step_note"] = (
    "Both factor sequences satisfy CONSTANT-COEFFICIENT (in m) recurrences "
    "over Q(t,z,w): A of order 2, B of order 4. The diagonal "
    "a22(2m) = [s^m sig^m t^m] is a Hadamard-type product in m of an "
    "A/B-sequence on the s-side with a B/A-sequence on the sig-side, so its "
    "operator is obtained by SYMMETRIC PRODUCT of the two factor operators "
    "(order <= 2*4 = 8 before the t-extraction). Then [t^m] and the two "
    "umbral directions remain.")
results["both_verified"] = bool(r1["identity_holds_for_all_k_tested"]
                                and r2["identity_holds_for_all_k_tested"])
print("\nBOTH FACTOR RECURRENCES VERIFIED EXACTLY:", results["both_verified"])

json.dump(results, open("factor_recs_results.json", "w"), indent=1)
print("wrote factor_recs_results.json — send this back")
