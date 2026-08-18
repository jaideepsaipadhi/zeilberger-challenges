# hertz_fix.sage — get the encoding right, then trust the tool.
#
#   sage hertz_fix.sage
#
# WHERE WE ARE. The patched ore_algebra passes its own three doctests and runs
# ct() on our Hertzsprung control. But settle_q0.py showed its answer is not a
# telescoper for OUR integral:
#
#     q1 agrees exactly (theirs = -ours)
#     q0 differs by exactly  -2(x+1)(x^2+2x-1) = 2c(x)
#
# i.e. off by twice the boundary term. That is not a random disagreement; it
# is the signature of a mis-specified e^{-z} factor, since the boundary
# contribution scales with precisely that term. The suspect is the line I
# wrote in doc_examples.sage:
#
#     ann_z = D*Dz + (x^2 - x) - D          # "- D" for the e^{-z}
#
# DERIVATION, done properly. Let K = (1+x)/D with D = 1 + x - z x + z x^2, and
# let P = K * e^{-z} be the actual integrand. Then
#
#     d/dz log P = d/dz log K - 1 = -(d/dz D)/D - 1,     d/dz D = x^2 - x.
#
# So P satisfies  D * P_z + (x^2 - x) * P + D * P = 0, i.e. the annihilator is
#
#     ann_z = D*Dz + (x^2 - x) + D          # PLUS D
#
# and my "- D" had the wrong sign. Flipping it shifts q0 by a multiple of the
# boundary term, which is exactly the observed discrepancy.
#
# For x there is no e^{-z} factor, so d/dx log K = 1/(1+x) - (d/dx D)/D with
# d/dx D = 1 - z + 2 z x, giving
#
#     ann_x = (1+x)*D*Dx - [ D - (1+x)(1 - z + 2 z x) ]
#
# THE TARGET (proved, Session 12, certificate verified symbolically):
#     q1 = -x^2 (x-1)^2 (x+1)
#     q0 = -(x^2+1)(x^2+2x-1)
#     boundary  c(x) = -x^3 - 3x^2 - x + 1
# A run is only trusted if ct returns this (up to an overall scalar).
#
# The script tries BOTH sign conventions for the e^{-z} term and reports which
# reproduces the proved answer, so the result does not depend on my having got
# the sign right in this comment.

import json
from ore_algebra import *

R.<x, z> = ZZ[]
A.<Dx, Dz> = OreAlgebra(R)

D  = 1 + x - z*x + z*x^2
Dz_D = x^2 - x            # d/dz D
Dx_D = 1 - z + 2*z*x      # d/dx D

# proved answer, for comparison
q1_target = -x^5 + x^4 + x^3 - x^2
q0_target = -(x^2+1)*(x^2+2*x-1)

out = {}
print("target  q1 =", factor(q1_target))
print("target  q0 =", factor(expand(q0_target)))
print()

for tag, sgn in (("PLUS  D  (derived)", +1), ("MINUS D  (what I wrote)", -1)):
    print("=" * 60)
    print(tag)
    ann_x = (1+x)*D*Dx - ( D - (1+x)*Dx_D )
    ann_z = D*Dz + Dz_D + sgn*D
    print("   ann_x =", ann_x)
    print("   ann_z =", ann_z)
    try:
        I = A.ideal([ann_x, ann_z])
        T, C = I.ct(Dz)
        gens = list(T.gens())
        print("   telescoper gens:", gens)
        print("   certificates   :", list(C))
        rec = {"ann_x": str(ann_x), "ann_z": str(ann_z),
               "telescoper": [str(g) for g in gens],
               "certificates": [str(c) for c in C]}
        # compare with the proved operator up to an overall rational scalar
        matched = False
        for g in gens:
            cf = g.coefficients(sparse=False)
            if len(cf) != 2:
                continue
            g0, g1 = R(cf[0]), R(cf[1])
            # g = g1*Dx + g0 ; want g1/q1_target == g0/q0_target as rationals
            num = expand(g1*q0_target - g0*q1_target)
            if num == 0:
                matched = True
                scal = None
                try:
                    scal = (g1/q1_target)
                except Exception:
                    pass
                print("   >>> MATCHES the proved operator (scalar %s)" % scal)
        rec["matches_proved"] = matched
        if not matched:
            print("   does NOT match the proved operator")
        out[tag] = rec
    except Exception as e:
        print("   ERROR:", str(e)[:300])
        out[tag] = {"error": str(e)[:400]}
    print()

with open("hertz_fix_results.json", "w") as f:
    json.dump(out, f, indent=1)
print("wrote hertz_fix_results.json — send this back")
