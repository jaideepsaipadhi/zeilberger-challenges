#!/usr/bin/env python3
"""
build_H2.py — same construction, sized correctly this time.

WHAT WENT WRONG IN v1 (my error): P(x) has degree <= 6, so it needs only
C_1..C_6. v1 computed C up to index 12 symbolically in (t,z,w) and then
expanded the FULL product Q*H — the high-index C_m are enormous rational
functions and the expansion never finished. Fixed: compute C only to index 8
(6 for P, 7-8 for the tail check), and never expand more than needed.

ORDER OF WORK (cheap first, so a failure is caught in seconds not minutes):
  PHASE 1  everything at RANDOM RATIONAL (z,w): build Q, P, and check that
           P/Q reproduces C_m for m = 1..12. All arithmetic in Q(t): fast.
           If the construction is wrong, this fails immediately.
  PHASE 2  only if phase 1 is clean: build Q and P symbolically in (t,z,w)
           from C_1..C_6, and confirm ONE symbolic identity — the x^7
           coefficient of Q*H must vanish. That is the recurrence itself, so
           it is the meaningful symbolic check.

THE OFF-BY-ONE (do not reintroduce): the pure exponential form of A, B holds
for m >= 1, not m = 0 — that is why the degree-6 operator failed at k = 6 and
was clean from k = 7. So H = sum_{m>=1} C_m x^m, C_0 excluded (and C_0 = 0
anyway since B_0 = 0; the script asserts it).

Run:  python3 build_H2.py        (plain python3; sympy only)
Writes build_H2_results.json — SEND THIS BACK. Phase 1 should take seconds.
"""
import sympy as sp, json, time, random

x, t, z, w, Y, v, r1 = sp.symbols('x t z w Y v r1')
T0 = time.time()
results = {}

d1_s = sp.expand(1 + t - z*t - w)
d2_s = sp.expand(t*(1 - z - w))

def char_coeffs(d1, d2):
    """Degree-6 characteristic coefficients (descending) for C = A*B,
    reduced to the base field via e1 = -d1, e2 = d2."""
    e1, e2 = sp.expand(-d1), sp.expand(d2)
    r2s = e1 - r1
    minp = sp.Poly(r1**2 - e1*r1 + e2, r1)
    charY = sp.expand((Y - r1**2)**2 * (Y - r1*r2s)**2 * (Y - r2s**2)**2)
    out = []
    for c in sp.Poly(charY, Y).all_coeffs():
        rem = sp.rem(sp.Poly(sp.expand(c), r1), minp).as_expr()
        out.append(sp.cancel(sp.expand(rem)))
    assert not any(r1 in c.free_symbols for c in out), "reduction left r1"
    return out

def gen_C(nmax, d1, d2):
    """C_m = A_m*B_m for m = 0..nmax, in whatever field d1,d2 live in."""
    Dv = 1 + d1*v + d2*v**2
    F1 = (1 + v)*(1 + v*t)/Dv
    F2 = -t*v**2*(1 + v)*(1 + v*t)/Dv**2
    s1 = sp.expand(sp.series(F1, v, 0, 5).removeO())
    s2 = sp.expand(sp.series(F2, v, 0, 5).removeO())
    A = [sp.cancel(s1.coeff(v, k)) for k in range(5)]
    B = [sp.cancel(s2.coeff(v, k)) for k in range(5)]
    cB = [None, sp.expand(2*d1), sp.expand(d1**2 + 2*d2),
          sp.expand(2*d1*d2), sp.expand(d2**2)]
    while len(A) <= nmax:
        k = len(A); A.append(sp.cancel(-(d1*A[k-1] + d2*A[k-2])))
    while len(B) <= nmax:
        k = len(B)
        B.append(sp.cancel(-(cB[1]*B[k-1] + cB[2]*B[k-2]
                             + cB[3]*B[k-3] + cB[4]*B[k-4])))
    return [sp.cancel(A[k]*B[k]) for k in range(nmax+1)]

ORDER = 6

def build_PQ(d1, d2, nC=8):
    """Q(x) = sum_j red[j] x^j ; P(x) = truncation of Q*H to degree ORDER."""
    red = char_coeffs(d1, d2)
    Qx = sum(red[j]*x**j for j in range(ORDER+1))
    C = gen_C(nC, d1, d2)
    assert sp.cancel(C[0]) == 0, "C_0 must vanish (B_0 = 0)"
    H = sum(C[mm]*x**mm for mm in range(1, nC+1))
    prod = sp.Poly(sp.expand(sp.cancel(Qx*H)), x)
    P = sum(sp.cancel(prod.coeff_monomial(x**k))*x**k for k in range(ORDER+1))
    tail = [sp.cancel(prod.coeff_monomial(x**k)) for k in (ORDER+1, ORDER+2)]
    return Qx, P, tail, C, red

# ---------------------------------------------------------------- PHASE 1
print("PHASE 1 — construction and check at random rational (z,w)", flush=True)
rnd = random.Random(7)
checks = []
phase1_ok = True
for trial in range(3):
    zv = sp.Rational(rnd.randint(2, 9), rnd.randint(2, 7))
    wv = sp.Rational(rnd.randint(2, 9), rnd.randint(2, 7))
    d1n = sp.expand(d1_s.subs({z: zv, w: wv}))
    d2n = sp.expand(d2_s.subs({z: zv, w: wv}))
    Qn, Pn, tailn, Cn, _ = build_PQ(d1n, d2n, nC=8)
    ser = sp.expand(sp.series(sp.cancel(Pn/Qn), x, 0, 13).removeO())
    Cbig = gen_C(12, d1n, d2n)
    bad = [mm for mm in range(1, 13)
           if sp.cancel(sp.together(ser.coeff(x, mm) - Cbig[mm])) != 0]
    tail_ok = all(sp.cancel(c) == 0 for c in tailn)
    checks.append({"z": str(zv), "w": str(wv), "mismatches": bad,
                   "tail_vanishes": bool(tail_ok)})
    if bad or not tail_ok:
        phase1_ok = False
    print("   z=%s w=%s: %s | tail vanishes: %s (%.0fs)"
          % (zv, wv, "clean" if not bad else "MISMATCH at %s" % bad[:5],
             tail_ok, time.time()-T0), flush=True)
    results["point_checks"] = checks
    json.dump(results, open("build_H2_results.json", "w"), indent=1)

results["phase1_ok"] = phase1_ok
if not phase1_ok:
    print("\nPHASE 1 FAILED — construction is wrong; not attempting the")
    print("symbolic phase. Mismatch pattern says where: from m=1 means the")
    print("initial-value handling; only high m means the tail degree.")
    json.dump(results, open("build_H2_results.json", "w"), indent=1)
    raise SystemExit

# ---------------------------------------------------------------- PHASE 2
print("\nPHASE 2 — symbolic Q, P in (t,z,w) from C_1..C_6", flush=True)
Qs, Ps, tails, Cs, reds = build_PQ(d1_s, d2_s, nC=ORDER+2)
results["Q_x"] = sp.sstr(sp.factor(sp.expand(Qs)))
results["char_coeffs_desc"] = [sp.sstr(sp.factor(c)) for c in reds]
print("   Q, P built (%.0fs)" % (time.time()-T0), flush=True)
json.dump(results, open("build_H2_results.json", "w"), indent=1)

print("   checking the x^7 coefficient of Q*H vanishes symbolically ...",
      flush=True)
sym_ok = bool(sp.cancel(sp.together(tails[0])) == 0)
results["symbolic_x7_vanishes"] = sym_ok
print("   -> %s (%.0fs)" % (sym_ok, time.time()-T0), flush=True)

results["P_x"] = sp.sstr(sp.factor(sp.expand(Ps)))
results["H_verified"] = bool(phase1_ok and sym_ok)
print("\nH(x,t) = P/Q VERIFIED:", results["H_verified"])
results["next"] = ("Diagonal [x^m t^m] of the bivariate rational H = P/Q "
                   "(Furstenberg residue) -> algebraic function; then the "
                   "two umbral directions (d_z-1),(d_w-1) with z=0/w=0 "
                   "boundary kills; then parity recombination and "
                   "right-division against KK's order-8/degree-11 operator.")
json.dump(results, open("build_H2_results.json", "w"), indent=1)
print("wrote build_H2_results.json — send this back")
