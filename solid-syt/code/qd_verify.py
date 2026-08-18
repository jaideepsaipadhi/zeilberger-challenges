import sympy as sp
from functools import lru_cache
import sys, json
sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

x,t=sp.symbols('x t')
N=17

Tser=sp.Integer(0)
for _ in range(N+3):
    Tser=sp.expand(sp.series(t*(2+Tser**3), t, 0, N+3).removeO())

# corrected R00
R00=sp.expand(sp.series((4*Tser-Tser**4)/(8*t), t,0,N).removeO())
R00_dp=sum(RK(m,0,0)*t**m for m in range(N))
print("R00 = (4T-T^4)/(8t):", sp.expand(R00-R00_dp)==0)

# the derived closed form for Qd_rev
sq=sp.sqrt(1-x*Tser**2).series(t,0,N).removeO()
sq=sp.expand(sq)
Qd_formula = sp.expand(sp.series( (Tser/(2*t) + (sq - 1 + x*Tser**2/2)/(t*x**2)) / sq, t, 0, N-1).removeO())
Qd_dp = sum(RK(m,i,i)*x**i*t**m for m in range(N) for i in range(m//2+1))
d=sp.expand(Qd_formula - Qd_dp)
ok=all(sp.simplify(d.coeff(t,m))==0 for m in range(N-1))
print("Qd_rev closed form matches DP diagonal to t^%d:"%(N-2), ok)
if not ok:
    for m in range(N-1):
        dm=sp.simplify(d.coeff(t,m))
        if dm!=0: print("  first mismatch t^%d:"%m, dm); break

# Lemma B chain: T(n) = [t^{3n}] Qd_rev(t;t)
Qdt=sp.expand(Qd_formula.subs(x,t))
Qdt=sp.expand(sp.series(Qdt,t,0,N-2).removeO())
P=json.load(open('pieces.json'))
Tvals=P['T']
okT=all(Qdt.coeff(t,3*n)==Tvals[n] for n in range((N-2)//3+1) if 3*n<N-2)
print("T(n) = [t^{3n}] Qd_rev(t;t):", okT, "  (n up to %d)"%((N-3)//3))
