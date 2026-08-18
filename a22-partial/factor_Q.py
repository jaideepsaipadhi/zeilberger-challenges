#!/usr/bin/env python3
"""
factor_Q.py — find the RIGHT certificate denominator shape (idea 3).

WHY. Every search so far used Theta = x^a t^b Q^c with a flat monomial box on
top, and every shape returned tele = 0 — including U~2500 with cross-terms.
That pattern is more consistent with a WRONG DENOMINATOR SHAPE than with a
box that is merely too small: if the true certificate has, say, low numerator
degree over Q1^2 * Q2 (Q1, Q2 irreducible factors of Q with different
exponents), no ansatz of the form N/(x^a t^b Q^c) can express it, and no
amount of extra U will ever find it.

Certificate denominators in creative telescoping are built from the SINGULAR
STRUCTURE of the integrand, i.e. from the irreducible factors of the
denominator, each with its own exponent — not from the denominator as a
single block.

THIS SCRIPT DOES NO SEARCHING. It reports, exactly:
  1. the irreducible factorisation of Q over Q(t,z,w)[x] and over Q[x,t,z,w];
  2. the same for the full denominator of R = H/(x t) (i.e. x*t*Q);
  3. which factors involve which variables — this decides which factors can
     appear in which direction's certificate (a factor free of z cannot help
     the z-direction's boundary behaviour, and a factor vanishing at z=0
     would BREAK the N_z(0)=0 boundary kill, which is exactly the error that
     produced spurious solutions in omega_step v1);
  4. whether each factor vanishes at z=0 or w=0 — anything that does must be
     excluded from Theta, or the umbral boundary argument is invalid;
  5. the discriminant/repeated-factor structure, since repeated factors are
     what force higher exponents in the certificate.

Output feeds directly into a redesigned Theta of the form
    Theta = x^a t^b * prod_i Q_i^{c_i}
with per-factor, per-direction exponents.

Run:  python3 factor_Q.py        (uses H_cache.pkl; sympy only)
Writes factor_Q_results.json — SEND THIS BACK.
"""
import sympy as sp, json, os, pickle, time

x, t, z, w = sp.symbols('x t z w')
T0 = time.time()
out = {}

if not os.path.exists("H_cache.pkl"):
    raise SystemExit("H_cache.pkl missing — run final_step3/4/5/6 once first")
Ps, Qs, ORDER = pickle.load(open("H_cache.pkl", "rb"))
Px, Qx = sp.sympify(Ps), sp.sympify(Qs)
print("H loaded (%.0fs); Q has x-degree %d" % (time.time()-T0, ORDER), flush=True)
out["order"] = ORDER

def describe(fac, tag):
    """Report a factorisation as a list of (factor, multiplicity, info)."""
    items = []
    for f_, e_ in fac:
        fs = sp.sstr(sp.expand(f_))
        info = {
            "factor": fs if len(fs) < 400 else fs[:400] + " ...(truncated)",
            "multiplicity": int(e_),
            "vars": sorted(str(s) for s in f_.free_symbols),
            "degree_x": int(sp.Poly(f_, x).degree()) if x in f_.free_symbols else 0,
            "vanishes_at_z0": bool(sp.simplify(f_.subs(z, 0)) == 0),
            "vanishes_at_w0": bool(sp.simplify(f_.subs(w, 0)) == 0),
        }
        items.append(info)
        print("   [%s] mult %d, vars %s, deg_x %d, z0=%s w0=%s"
              % (tag, info["multiplicity"], ",".join(info["vars"]),
                 info["degree_x"], info["vanishes_at_z0"], info["vanishes_at_w0"]),
              flush=True)
        print("        %s" % (fs if len(fs) < 200 else fs[:200] + " ..."), flush=True)
    return items

# ---------------------------------------------------------------- 1. Q
print("\nfactoring Q over Q[x,t,z,w] ...", flush=True)
fQ = sp.factor_list(sp.expand(Qx))
out["Q_content"] = sp.sstr(fQ[0])
print("   content:", sp.sstr(fQ[0]), flush=True)
out["Q_factors"] = describe(fQ[1], "Q")
print("   (%.0fs)" % (time.time()-T0), flush=True)

# ---------------------------------------------------------------- 2. Q as poly in x
print("\nfactoring Q as a polynomial in x over Q(t,z,w) ...", flush=True)
try:
    Qpx = sp.Poly(sp.expand(Qx), x)
    fQx = sp.factor_list(Qpx.as_expr(), x)
    out["Q_in_x_content"] = sp.sstr(fQx[0])
    out["Q_in_x_factors"] = describe(fQx[1], "Q|x")
except Exception as e:
    out["Q_in_x_error"] = str(e)[:300]
    print("   failed:", str(e)[:200], flush=True)

# ---------------------------------------------------------------- 3. P
print("\nfactoring P ...", flush=True)
try:
    fP = sp.factor_list(sp.expand(Px))
    out["P_content"] = sp.sstr(fP[0])
    out["P_factors"] = describe(fP[1], "P")
except Exception as e:
    out["P_error"] = str(e)[:300]
    print("   failed:", str(e)[:200], flush=True)

# ---------------------------------------------------------------- 4. shared
print("\ndo P and Q share factors (would cancel, changing the true denominator)?",
      flush=True)
try:
    g = sp.gcd(sp.expand(Px), sp.expand(Qx))
    out["gcd_P_Q"] = sp.sstr(sp.factor(g))
    print("   gcd(P,Q) =", sp.sstr(sp.factor(g))[:200], flush=True)
except Exception as e:
    out["gcd_error"] = str(e)[:200]

# ---------------------------------------------------------------- 5. verdict
zsafe = [f for f in out.get("Q_factors", [])
         if not f["vanishes_at_z0"] and not f["vanishes_at_w0"]]
zbad = [f for f in out.get("Q_factors", [])
        if f["vanishes_at_z0"] or f["vanishes_at_w0"]]
out["factors_usable_in_Theta"] = len(zsafe)
out["factors_forbidden_in_Theta"] = len(zbad)
print("\nVERDICT")
print("   factors safe for Theta (nonzero at z=0 and w=0): %d" % len(zsafe))
print("   factors that MUST be excluded (vanish at z=0 or w=0): %d" % len(zbad))
if zbad:
    print("   -> any Theta containing those breaks the N_z(0)=N_w(0)=0 boundary")
    print("      kill, which is exactly how omega_step v1 produced spurious")
    print("      solutions. The previous searches used Q^c wholesale, so if any")
    print("      forbidden factor is inside Q, EVERY search so far was malformed.")
out["next"] = ("Redesign Theta = x^a t^b * prod_i Q_i^{c_i} over the usable "
               "factors, with per-factor exponents and per-direction choices, "
               "then re-run the search over that (much smaller, correctly "
               "shaped) space.")
json.dump(out, open("factor_Q_results.json", "w"), indent=1)
print("\nwrote factor_Q_results.json — send this back")
