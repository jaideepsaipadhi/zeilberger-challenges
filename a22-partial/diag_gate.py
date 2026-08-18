#!/usr/bin/env python3
"""
diag_gate.py — END-TO-END GATE on the whole chain built today.

CLAIM UNDER TEST:
    e_X(m) = <k! j!> [t^m] ( 2 * A_m(t,z,w) * B_m(t,z,w) )
where A_m = [sigma^m]F1, B_m = [sigma^m]F2 (both C-finite, verified), and
<k!j!> is the double umbral z^k w^j -> k! j!.

KNOWN EXACT VALUES (independently established earlier in the campaign, used
as gates for the 400-term data run):
    e_X(1) = 0,   e_X(2) = -4,   e_X(3) = -48
    o_X(2) = -1,  o_X(3) = -11
so this checks the composition of EVERY piece: the factor recurrences, the
Hadamard product, the diagonal framing, and the umbral convention.

WHY THIS BEFORE THE TELESCOPING: the three-direction telescoping is the one
structurally untested step left. Running it against a chain that has a sign
or convention slip anywhere would produce an empty search that tells us
nothing — exactly the failure mode that cost this campaign three separate
multi-session detours. C_m is tiny for small m, so this costs seconds.

ODD PARITY: a22(2m-1) uses [s^m sigma^(m-1) t^m], giving
    o_X(m) = <k!j!> [t^m] ( B_m A_{m-1} + A_m B_{m-1} ).

Run:  python3 diag_gate.py       (plain python3; sympy only)
Writes diag_gate_results.json — SEND THIS BACK. Should take seconds.
"""
import sympy as sp, json

v, t, z, w = sp.symbols('v t z w')

d1 = sp.expand(1 + t - z*t - w)
d2 = sp.expand(t*(1 - z - w))
D = 1 + d1*v + d2*v**2
F1 = (1 + v)*(1 + v*t)/D
F2 = -t*v**2*(1 + v)*(1 + v*t)/D**2

MMAX = 4
print("expanding F1, F2 in sigma to order %d ..." % MMAX, flush=True)
s1 = sp.expand(sp.series(F1, v, 0, MMAX+1).removeO())
s2 = sp.expand(sp.series(F2, v, 0, MMAX+1).removeO())
A = [sp.expand(sp.cancel(s1.coeff(v, k))) for k in range(MMAX+1)]
B = [sp.expand(sp.cancel(s2.coeff(v, k))) for k in range(MMAX+1)]

results = {"A": [sp.sstr(sp.factor(a)) for a in A],
           "B": [sp.sstr(sp.factor(b)) for b in B]}
print("A_m:", results["A"])
print("B_m:", results["B"])

def umbral(expr):
    """<k! j!> : z^k w^j -> k! j!.  expr must be polynomial in z, w."""
    e = sp.expand(expr)
    p = sp.Poly(e, z, w)
    tot = 0
    for (kk, jj), c in zip(p.monoms(), p.coeffs()):
        tot += c * sp.factorial(kk) * sp.factorial(jj)
    return sp.expand(tot)

def coeff_t(expr, m):
    """[t^m] of a polynomial in t (C_m is polynomial in t here)."""
    e = sp.expand(expr)
    return sp.expand(e.coeff(t, m))

# ---------------------------------------------------------------- even
print("\nEVEN: e_X(m) = <k!j!>[t^m](2 A_m B_m)")
KNOWN_E = {1: 0, 2: -4, 3: -48}
even = {}
for m in range(1, MMAX+1):
    C = sp.expand(2 * A[m] * B[m])
    ct = coeff_t(C, m)
    val = umbral(ct)
    even[m] = sp.sstr(val)
    tag = ""
    if m in KNOWN_E:
        tag = "  expected %d  -> %s" % (KNOWN_E[m],
                                        "MATCH" if sp.simplify(val - KNOWN_E[m]) == 0
                                        else "MISMATCH")
    print("   e_X(%d) = %s%s" % (m, val, tag), flush=True)
results["e_X"] = even
results["e_X_expected"] = {str(k): v for k, v in KNOWN_E.items()}
results["e_X_all_match"] = all(
    sp.simplify(sp.sympify(even[m]) - KNOWN_E[m]) == 0 for m in KNOWN_E)

# ---------------------------------------------------------------- odd
print("\nODD: o_X(m) = <k!j!>[t^m](B_m A_{m-1} + A_m B_{m-1})")
KNOWN_O = {2: -1, 3: -11}
odd = {}
for m in range(2, MMAX+1):
    C = sp.expand(B[m]*A[m-1] + A[m]*B[m-1])
    ct = coeff_t(C, m)
    val = umbral(ct)
    odd[m] = sp.sstr(val)
    tag = ""
    if m in KNOWN_O:
        tag = "  expected %d  -> %s" % (KNOWN_O[m],
                                        "MATCH" if sp.simplify(val - KNOWN_O[m]) == 0
                                        else "MISMATCH")
    print("   o_X(%d) = %s%s" % (m, val, tag), flush=True)
results["o_X"] = odd
results["o_X_expected"] = {str(k): v for k, v in KNOWN_O.items()}
results["o_X_all_match"] = all(
    sp.simplify(sp.sympify(odd[m]) - KNOWN_O[m]) == 0 for m in KNOWN_O)

results["chain_verified"] = bool(results["e_X_all_match"] and results["o_X_all_match"])
print("\nCHAIN VERIFIED END TO END:", results["chain_verified"])
if not results["chain_verified"]:
    print("A mismatch localises the slip: if e_X(1)=0 holds but e_X(2) is off")
    print("by a sign or factor, suspect the 2*A_m*B_m normalisation or the")
    print("W/X sign convention; if everything is off, suspect the umbral or")
    print("the [t^m] index alignment.")

results["next"] = ("If clean: run the three-direction telescoping "
                   "sum_i p_i(m) Omega^i R = (d_u + m lam_u)G_u + (d_z-1)G_z "
                   "+ (d_w-1)G_w on the Furstenberg diagonal integrand, "
                   "gated against the 200-term e_X data mod 46337.")
json.dump(results, open("diag_gate_results.json", "w"), indent=1)
print("wrote diag_gate_results.json — send this back")
