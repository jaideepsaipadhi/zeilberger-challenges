#!/usr/bin/env python3
"""
settle_q0.py — are ore_algebra's telescoper and ours both valid?

THE SITUATION. With the multivariate path now working, ore_algebra's ct() on
our Hertzsprung control returns

    telescoper  T_them = (x^5 - x^4 - x^3 + x^2) Dx + (x^4 - 6x^2 + 1)
    second out  (x^3 + 3x^2 + x - 1)   [printed as -x^3-3x^2-x+1]

Our proved answer (Session 12, certificate verified symbolically) is

    T_ours = q1 Dx + q0,  q1 = -(x^5 - x^4 - x^3 + x^2),
                          q0 = -(x^4 + 2x^3 + 2x - 1)
    with    q1 F' + q0 F = c(x) = -(x+1)(x^2+2x-1) = -x^3-3x^2-x+1.

The q1's agree up to sign. The q0's do NOT: theirs is x^4-6x^2+1, ours is
-(x^4+2x^3+2x-1). Yesterday I blamed my own mishandling of ct's return pair;
that explanation is now dead, because the pair is unpacked correctly here and
the difference persists.

WHAT THIS SCRIPT DOES. Applies BOTH operators to the actual A002464 series
and reports the right-hand side each produces:

    F(x) = sum b(n) x^n,   b = 1,1,0,0,2,14,90,646,5242,...

  * a valid telescoper for this integral gives a POLYNOMIAL RHS (the boundary
    term), so all coefficients past its degree must vanish;
  * if BOTH give polynomial RHS, both are valid elements of the telescoper
    ideal and there is nothing to reconcile -- ct simply returned a different
    generator;
  * if only ours does, ore_algebra's orientation or sign convention differs
    from the integral I encoded, and that must be understood BEFORE trusting
    its output on a_{2,2}.

It also checks whether the two operators are proportional over Q(x) (they are
not, unless the q0 discrepancy is illusory) and whether their difference
annihilates F, which is the precise statement that both lie in the same
telescoper ideal.

Run:  python3 settle_q0.py        (sympy only, seconds)
"""
import sympy as sp, json

x = sp.Symbol('x')
NT = 30

# ---------------------------------------------------------------- data
b = [1, 1, 0, 0, 2, 14, 90, 646, 5242, 47622, 479306, 5296790, 63779034]
while len(b) < NT + 6:
    n = len(b) - 4
    b.append((n+5)*b[n+3] - (n+2)*b[n+2] - (n-1)*b[n+1] + (n+1)*b[n])
F = [sp.Integer(v) for v in b[:NT]]
print("A002464:", [int(v) for v in F[:9]])

def deriv(seq):
    return [sp.Integer(i)*seq[i] for i in range(1, len(seq))] + [sp.Integer(0)]

Fp = deriv(F)

def apply_op(q1, q0, label):
    """Compute the series of q1*F' + q0*F and report the RHS."""
    c1 = sp.Poly(sp.expand(q1), x).all_coeffs()[::-1]   # ascending
    c0 = sp.Poly(sp.expand(q0), x).all_coeffs()[::-1]
    res = [sp.Integer(0)]*NT
    for i, cc in enumerate(c1):
        for j in range(NT - i):
            res[i+j] += cc * Fp[j]
    for i, cc in enumerate(c0):
        for j in range(NT - i):
            res[i+j] += cc * F[j]
    nz = [(i, v) for i, v in enumerate(res) if v != 0]
    tail = [(i, v) for i, v in nz if i >= 8]
    print("\n%s" % label)
    print("   q1 =", sp.sstr(sp.factor(q1)))
    print("   q0 =", sp.sstr(sp.factor(q0)))
    print("   RHS coefficients (nonzero):", [(i, int(v)) for i, v in nz[:10]])
    print("   coefficients from index 8 on:", len(tail),
          "nonzero" if tail else "none  -> POLYNOMIAL RHS")
    return {"label": label, "q1": sp.sstr(q1), "q0": sp.sstr(q0),
            "rhs_nonzero": [(i, str(v)) for i, v in nz[:12]],
            "tail_count": len(tail),
            "polynomial_rhs": not tail}

out = {}

# ours (proved)
q1_ours = -(x**5 - x**4 - x**3 + x**2)
q0_ours = -(x**4 + 2*x**3 + 2*x - 1)
out["ours"] = apply_op(q1_ours, q0_ours, "OURS (proved Session 12)")

# theirs, as returned
q1_them = x**5 - x**4 - x**3 + x**2
q0_them = x**4 - 6*x**2 + 1
out["theirs"] = apply_op(q1_them, q0_them, "ORE_ALGEBRA (as returned)")

# theirs with overall sign flipped (q1 matches ours after negation)
out["theirs_negated"] = apply_op(-q1_them, -q0_them,
                                 "ORE_ALGEBRA (overall sign flipped)")

# ---------------------------------------------------------------- relation
print("\n--- relation between the two operators ---")
diff_q1 = sp.expand(q1_ours + q1_them)     # ours = -theirs on q1
diff_q0 = sp.expand(q0_ours + q0_them)
print("   q1_ours + q1_theirs =", sp.sstr(sp.factor(diff_q1)))
print("   q0_ours + q0_theirs =", sp.sstr(sp.factor(diff_q0)))
out["q1_sum"] = sp.sstr(diff_q1)
out["q0_sum"] = sp.sstr(diff_q0)

# does the DIFFERENCE of the two operators annihilate F?
# (ours) - (-theirs) = (q0_ours + q0_them) as a multiplication operator
res = [sp.Integer(0)]*NT
dc = sp.Poly(sp.expand(diff_q0), x).all_coeffs()[::-1]
for i, cc in enumerate(dc):
    for j in range(NT - i):
        res[i+j] += cc * F[j]
nz = [(i, int(v)) for i, v in enumerate(res) if v != 0]
print("   (q0_ours+q0_theirs)*F  nonzero coefficients:", nz[:8],
      "..." if len(nz) > 8 else "")
out["difference_times_F"] = nz[:12]

print("\n--- verdict ---")
if out["ours"]["polynomial_rhs"] and (out["theirs"]["polynomial_rhs"]
                                      or out["theirs_negated"]["polynomial_rhs"]):
    print("   BOTH valid: ct returned a different generator of the same")
    print("   telescoper ideal. ore_algebra can be trusted on a_{2,2}.")
    out["verdict"] = "both valid"
elif out["ours"]["polynomial_rhs"]:
    print("   ONLY OURS is valid for this integral. ore_algebra's output")
    print("   corresponds to a different problem than the one encoded --")
    print("   the annihilator construction in doc_examples.sage must be")
    print("   re-derived before trusting ct on a_{2,2}.")
    out["verdict"] = "ours only — encoding differs"
else:
    print("   NEITHER gives a polynomial RHS — the series or the convention")
    print("   in this script is wrong. Stop and re-check before proceeding.")
    out["verdict"] = "neither — script convention wrong"

json.dump(out, open("settle_q0_results.json", "w"), indent=1)
print("\nwrote settle_q0_results.json")
