# a22_continue.sage — finish the telescoping: Dt, then Dz, then Dw.
#
#   sage a22_continue.sage
#
# WHERE WE ARE. a22_ct.sage established:
#   * STAGE 0 gate passed exactly — the diagonal encoding (Sm as multiplication
#     by 1/(xt)) plus iterated ct reproduces the central binomial recurrence
#     (m+1)a(m+1) = 2(2m+1)a(m) with no fuss;
#   * STAGE 2 ran ct(Dx) on the real a_{2,2} ideal: 5 generators, 1028s.
# That was the step the pole-order bound said a dense solver could not reach
# (~70,000 unknowns). Gröbner machinery did it in 17 minutes.
#
# THIS SCRIPT does the remaining three eliminations, checkpointing after each
# so a kill at any point leaves usable state:
#     ct(Dt)  -> operator in Sm, Dz, Dw
#     ct(Dz)  -> operator in Sm, Dw
#     ct(Dw)  -> operator in Sm alone   <- the target
#
# WHAT TO EXPECT AT THE END. As on both controls, ct returns an INHOMOGENEOUS
# relation; the boundary term is not killed automatically. Homogenise with
#     (g(m) S - g(m+1)) o L
# (note the ordering — the reverse leaves g(m+1)^2 - g(m)^2 and is wrong), then
# check against the 200-term e_X data and the independently measured shape,
# order 6 and degree 16.
#
# COST. Each step works in a smaller algebra but on larger coefficients, so
# ct(Dt) is the real indicator: if it lands in an hour or so the route is
# viable end to end; if it runs overnight, that is still information and the
# checkpoint file will say how far it got.
#
# The ideal is rebuilt from scratch here rather than resumed, because ct's
# output was not serialised last time. Rebuilding costs ~17 min (the ct(Dx)
# rerun) before the new work starts.

import json, time, sys, threading, pickle, os
sys.setrecursionlimit(200000)
threading.stack_size(512*1024*1024)
from ore_algebra import *

T0 = time.time()
def el(): return "%.0fs" % (time.time() - T0)
out = {}

def save():
    with open("a22_continue_results.json", "w") as f:
        json.dump(out, f, indent=1)

# ---------------------------------------------------------------- rebuild
if not os.path.exists("H_cache.pkl"):
    raise SystemExit("H_cache.pkl not found")
import sympy
Ps, Qs, ORDER = pickle.load(open("H_cache.pkl", "rb"))
R.<m, x, t, z, w> = ZZ[]
A.<Sm, Dx, Dt, Dz, Dw> = OreAlgebra(R)
P = R(str(sympy.sympify(Ps)).replace('**', '^'))
Q = R(str(sympy.sympify(Qs)).replace('**', '^'))
print("P, Q loaded [%s]" % el(), flush=True)

ann_m = (x*t)*Sm - 1
ann_x = x*P*Q*Dx - ( x*(P.derivative(x)*Q - P*Q.derivative(x)) - (m+1)*P*Q )
ann_t = t*P*Q*Dt - ( t*(P.derivative(t)*Q - P*Q.derivative(t)) - (m+1)*P*Q )
ann_z = P*Q*Dz - ( (P.derivative(z)*Q - P*Q.derivative(z)) - P*Q )
ann_w = P*Q*Dw - ( (P.derivative(w)*Q - P*Q.derivative(w)) - P*Q )

I = A.ideal([ann_m, ann_x, ann_t, ann_z, ann_w])
print("ideal built [%s] — ct(Dx) (~17 min, known to complete)" % el(), flush=True)
Tx, Cx = I.ct(Dx)
gens_x = list(Tx)
print("ct(Dx): %d generators [%s]" % (len(gens_x), el()), flush=True)
out["ct_Dx"] = {"n_generators": len(gens_x), "elapsed": el(),
                "generators": [str(g)[:1500] for g in gens_x]}
save()

# ---------------------------------------------------------------- ct(Dt)
print("\nct(Dt) — the viability indicator", flush=True)
try:
    R2.<m, t, z, w> = ZZ[]
    A2.<Sm, Dt, Dz, Dw> = OreAlgebra(R2)
    g2 = [A2(str(g)) for g in gens_x]
    I2 = A2.ideal(g2)
    print("   ideal rebuilt in Sm,Dt,Dz,Dw [%s]" % el(), flush=True)
    Tt, Ct = I2.ct(Dt)
    gens_t = list(Tt)
    print("   ct(Dt): %d generators [%s]" % (len(gens_t), el()), flush=True)
    out["ct_Dt"] = {"n_generators": len(gens_t), "elapsed": el(),
                    "generators": [str(g)[:1500] for g in gens_t]}
    save()
except Exception as e:
    print("   ERROR/STOP:", str(e)[:400], flush=True)
    out["ct_Dt"] = {"error": str(e)[:800], "elapsed": el()}
    save()
    raise SystemExit

# ---------------------------------------------------------------- ct(Dz)
print("\nct(Dz)", flush=True)
try:
    R3.<m, z, w> = ZZ[]
    A3.<Sm, Dz, Dw> = OreAlgebra(R3)
    g3 = [A3(str(g)) for g in gens_t]
    I3 = A3.ideal(g3)
    Tz, Cz = I3.ct(Dz)
    gens_z = list(Tz)
    print("   ct(Dz): %d generators [%s]" % (len(gens_z), el()), flush=True)
    out["ct_Dz"] = {"n_generators": len(gens_z), "elapsed": el(),
                    "generators": [str(g)[:1500] for g in gens_z]}
    save()
except Exception as e:
    print("   ERROR/STOP:", str(e)[:400], flush=True)
    out["ct_Dz"] = {"error": str(e)[:800], "elapsed": el()}
    save()
    raise SystemExit

# ---------------------------------------------------------------- ct(Dw)
print("\nct(Dw) — the last one; result should be an operator in Sm alone",
      flush=True)
try:
    R4.<m, w> = ZZ[]
    A4.<Sm, Dw> = OreAlgebra(R4)
    g4 = [A4(str(g)) for g in gens_z]
    I4 = A4.ideal(g4)
    Tw, Cw = I4.ct(Dw)
    gens_w = list(Tw)
    print("   ct(Dw): %d generators [%s]" % (len(gens_w), el()), flush=True)
    for g in gens_w:
        print("   ", str(g)[:2000], flush=True)
    out["ct_Dw"] = {"n_generators": len(gens_w), "elapsed": el(),
                    "generators": [str(g) for g in gens_w],
                    "certificates": [str(c)[:1500] for c in Cw]}
    save()
    print("\n>>> THIS IS THE TARGET OPERATOR (expect inhomogeneous;")
    print(">>> homogenise with (g(m) S - g(m+1)) o L, then check against")
    print(">>> the e_X data and the measured order 6 / degree 16.)")
except Exception as e:
    print("   ERROR/STOP:", str(e)[:400], flush=True)
    out["ct_Dw"] = {"error": str(e)[:800], "elapsed": el()}
    save()

save()
print("\nwrote a22_continue_results.json — send this back whatever happened")
