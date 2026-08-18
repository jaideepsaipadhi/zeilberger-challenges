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
N=15

# T-series (her T = t(2+T^3)) -- same as BM's W
Tser=sp.Integer(0)
for _ in range(N+2):
    Tser=sp.expand(sp.series(t*(2+Tser**3), t, 0, N+2).removeO())

# R00 = Q6(0,0;t) = (4T - T^2)/(8t)
R00_formula=sp.expand(sp.series((4*Tser-Tser**2)/(8*t), t,0,N).removeO())
R00_dp=sum(RK(m,0,0)*t**m for m in range(N))
print("R00 = (4T-T^2)/(8t):", sp.expand(R00_formula-R00_dp)==0)

# Table version: S(x,t) = [(-2x/T)(1 - T^2/(2x)) + 1/(tx)]*sqrt(U)/2 + (1 - tx - t/x^2)*x/(2t)
U=1 - x*Tser*(1+Tser**3/4) + x**2*Tser**2/4
sqU=sp.sqrt(U).series(t,0,N).removeO()
S_tab=sp.expand(sp.series(((-2*x/Tser)*(1-Tser**2/(2*x)) + 1/(t*x))*sqU/2 + (1-t*x-t/x**2)*x/(2*t), t,0,N-1).removeO())
# candidate identification: S(x,t) = t*Q(x,0) - (t/2)*Q(0,0)
Qx0=sum(RK(m,i,0)*x**i*t**m for m in range(N) for i in range(m+1))
cand=sp.expand(sp.series(t*Qx0-(t*sp.Rational(1,2))*R00_dp, t,0,N-1).removeO())
d=sp.expand(S_tab-cand)
ok=all(sp.simplify(d.coeff(t,m))==0 for m in range(N-1))
print("Mishna Thm2.2/table: S(x,t) = t*Q(x,0) - (t/2)*Q(0,0):", ok)
if not ok:
    for m in range(N-1):
        dm=sp.simplify(d.coeff(t,m))
        if dm!=0: print("  first mismatch t^%d:"%m, dm); break
