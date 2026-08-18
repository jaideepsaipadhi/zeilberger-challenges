# a22_ct.sage — encode the a_{2,2} diagonal for ore_algebra, gated first.
#
#   sage a22_ct.sage
#
# THE TARGET
#     e_X(m) = <k!j!> [x^m t^m] H(x,t;z,w),      H = P/Q  (built and verified)
# Realise both extractions as residues and the umbral functional as an
# integral against e^{-z-w}:
#     F = H(x,t;z,w) * x^{-m-1} * t^{-m-1} * e^{-z-w}
# and telescope x, t, z, w away, leaving an operator in Sm.
#
# ANNIHILATORS (recipe confirmed on the Hertzsprung and W-part controls; the
# e^{-z} sign is the MINUS convention, which cost two days to pin down):
#     Sm  acts as multiplication by 1/(x t)      ->  (x t) Sm - 1
#     d/dx log F = P_x/P - Q_x/Q - (m+1)/x
#         -> ann_x = x P Q Dx - [ x(P_x Q - P Q_x) - (m+1) P Q ]
#     d/dt log F   analogous
#     d/dz log F = P_z/P - Q_z/Q - 1
#         -> ann_z = P Q Dz - [ (P_z Q - P Q_z) - P Q ]
#     d/dw log F   analogous
#
# ct() takes ONE operator at a time and returns (telescopers, certificates)
# with the telescopers a plain LIST. eliminate() is limited to
# zero-dimensional ideals, so the plan is ITERATED ct: kill Dx, then Dt, then
# Dz, then Dw, rebuilding the ideal from the returned telescopers each time.
#
# EXPECT AN INHOMOGENEOUS RELATION at the end, as on both controls, and
# homogenise with   (g(m) S - g(m+1)) o L   [note the ordering; the reverse is
# wrong and leaves g(m+1)^2 - g(m)^2].
#
# STAGE 0 IS A GATE, not the target: the central binomial diagonal
#     [x^m t^m] 1/(1-x-t) = C(2m,m),   (m+1) a(m+1) = 2(2m+1) a(m)
# exercises exactly the Sm-as-1/(xt) encoding plus two residue directions, and
# the answer is elementary. If stage 0 fails, nothing below it is meaningful.
# STAGE 1 adds one umbral variable to the same gate. STAGE 2 is a_{2,2}.
#
# Stages run in order and each writes results as it finishes, so a run that is
# killed partway is still informative.

import json, time
from ore_algebra import *

T0 = time.time()
out = {}

def elapsed():
    return "%.0fs" % (time.time() - T0)

# ============================================================ STAGE 0
print("=" * 62)
print("STAGE 0 — gate: central binomial diagonal, known recurrence")
print("   [x^m t^m] 1/(1-x-t) = C(2m,m),  (m+1)a(m+1) = 2(2m+1)a(m)")
try:
    R0.<m, x, t> = ZZ[]
    A0.<Sm, Dx, Dt> = OreAlgebra(R0)
    D = 1 - x - t
    # F = 1/D * x^{-m-1} t^{-m-1}
    #   d/dx log F = 1/D - (m+1)/x   ->  x D Dx - [x - (m+1) D]  ... careful:
    #   d/dx log(1/D) = -D_x/D = 1/D  (since D_x = -1)
    ann_x = x*D*Dx - (x - (m+1)*D)
    ann_t = t*D*Dt - (t - (m+1)*D)
    ann_m = (x*t)*Sm - 1
    print("   ann_m =", ann_m)
    print("   ann_x =", ann_x)
    print("   ann_t =", ann_t)
    I = A0.ideal([ann_m, ann_x, ann_t])
    T1, C1 = I.ct(Dx)
    print("   after ct(Dx):", list(T1), "[%s]" % elapsed())
    A1.<Sm1, Dt1> = OreAlgebra(PolynomialRing(ZZ, 'm, t'))
    # rebuild in the smaller algebra
    gens = []
    for g in T1:
        gens.append(A1(str(g).replace('Sm', 'Sm1').replace('Dt', 'Dt1')))
    I1 = A1.ideal(gens)
    T2, C2 = I1.ct(Dt1)
    print("   after ct(Dt):", list(T2), "[%s]" % elapsed())
    out["stage0"] = {"after_Dx": [str(g) for g in T1],
                     "after_Dt": [str(g) for g in T2],
                     "expected": "(m+1) S - 2(2m+1)  up to scalar/inhomogeneity"}
except Exception as e:
    print("   ERROR:", str(e)[:400])
    out["stage0"] = {"error": str(e)[:600]}
json.dump(out, open("a22_ct_results.json", "w"), indent=1)

# ============================================================ STAGE 1
print()
print("=" * 62)
print("STAGE 1 — gate + one umbral variable (tests e^{-z} with a diagonal)")
try:
    R1.<m, x, t, z> = ZZ[]
    B.<Sm, Dx, Dt, Dz> = OreAlgebra(R1)
    D = 1 - x - t - z*x*t          # a z-dependent kernel so Dz is nontrivial
    Dx_ = -1 - z*t
    Dt_ = -1 - z*x
    Dz_ = -x*t
    ann_x = x*D*Dx - (-x*Dx_ - (m+1)*D)
    ann_t = t*D*Dt - (-t*Dt_ - (m+1)*D)
    ann_z = D*Dz - (-Dz_ - D)        # MINUS D  <- the e^{-z}
    ann_m = (x*t)*Sm - 1
    I = B.ideal([ann_m, ann_x, ann_t, ann_z])
    print("   ideal built [%s]" % elapsed())
    Tz, Cz = I.ct(Dx)
    print("   after ct(Dx): %d generator(s) [%s]" % (len(list(Tz)), elapsed()))
    out["stage1"] = {"after_Dx": [str(g) for g in Tz]}
except Exception as e:
    print("   ERROR:", str(e)[:400])
    out["stage1"] = {"error": str(e)[:600]}
json.dump(out, open("a22_ct_results.json", "w"), indent=1)

# ============================================================ STAGE 2
print()
print("=" * 62)
print("STAGE 2 — the real target: a_{2,2} via H = P/Q")
try:
    import pickle, os
    if not os.path.exists("H_cache.pkl"):
        raise RuntimeError("H_cache.pkl not found — rebuild it with "
                           "final_step3.py (it caches P, Q) before this stage")
    import sage.all as sg
    Ps, Qs, ORDER = pickle.load(open("H_cache.pkl", "rb"))
    print("   H loaded, Q has x-degree", ORDER, "[%s]" % elapsed())
    R2.<m, x, t, z, w> = ZZ[]
    Cc.<Sm, Dx, Dt, Dz, Dw> = OreAlgebra(R2)
    P = R2(str(sg.sympify(Ps)).replace('**', '^'))
    Q = R2(str(sg.sympify(Qs)).replace('**', '^'))
    print("   P, Q converted [%s]" % elapsed())
    ann_x = x*P*Q*Dx - ( x*(P.derivative(x)*Q - P*Q.derivative(x)) - (m+1)*P*Q )
    ann_t = t*P*Q*Dt - ( t*(P.derivative(t)*Q - P*Q.derivative(t)) - (m+1)*P*Q )
    ann_z = P*Q*Dz - ( (P.derivative(z)*Q - P*Q.derivative(z)) - P*Q )
    ann_w = P*Q*Dw - ( (P.derivative(w)*Q - P*Q.derivative(w)) - P*Q )
    ann_m = (x*t)*Sm - 1
    print("   annihilators built [%s]" % elapsed())
    I = Cc.ideal([ann_m, ann_x, ann_t, ann_z, ann_w])
    print("   ideal built [%s] — running ct(Dx), this is the expensive step"
          % elapsed())
    Tx, Cx = I.ct(Dx)
    print("   after ct(Dx): %d generator(s) [%s]" % (len(list(Tx)), elapsed()))
    out["stage2"] = {"after_Dx": [str(g)[:2000] for g in Tx],
                     "elapsed": elapsed()}
except Exception as e:
    print("   ERROR/STOP:", str(e)[:500])
    out["stage2"] = {"error": str(e)[:800], "elapsed": elapsed()}

json.dump(out, open("a22_ct_results.json", "w"), indent=1)
print()
print("wrote a22_ct_results.json — send this back (even if stage 2 was killed)")
