import sympy as sp
from functools import lru_cache
import sys, json
sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

t,x,T,sig,F,s=sp.symbols('t x T sigma F s')

# ---------- Part 1: verify the mirror-extracted axis GF ----------
N=14
Tser=sp.Integer(0)
for _ in range(N+3):
    Tser=sp.expand(sp.series(t*(2+Tser**3), t, 0, N+3).removeO())
c_ser=sp.expand(sp.series(t*(4*Tser-Tser**4)/(8*t), t,0,N).removeO())  # c = t*R00
Dm = 1 - (1/x)*Tser*(1+Tser**3/4) + (1/x**2)*Tser**2/4   # Delta_minus(xbar) with xbar=1/x
sqDm = sp.sqrt(Dm).series(t,0,N).removeO()
R_formula = sp.expand(sp.series((c_ser + x + 1/(x*t) - 1/x**2 - (x + 2/(x*Tser) - Tser)*sqDm)/2, t,0,N-1).removeO())
# R(xbar) means: this should equal t*Q(1/x... i.e. compare with t*sum RK x^{-i}: R(xbar)=t*Q(xbar,0)
R_dp = sp.expand(t*sum(RK(m,i,0)*x**(-i)*t**m for m in range(N) for i in range(m+1)))
d=sp.expand(R_formula-R_dp)
okR=all(sp.simplify(d.coeff(t,m))==0 for m in range(N-1))
print("Axis GF (mirror extraction) verified:", okR)

# ---------- Part 2: minimal polynomial of A(t)=Qd_rev(t;t) ----------
# sigma^2 = 1 - t*T^2 (x=t), and 2t^3*F*sig = t^2*T + 2*sig - 2 + t*T^2  =>  sig*(2t^3F-2)=t^2T+tT^2-2
P1 = (1-t*T**2)*(2*t**3*F-2)**2 - (t**2*T+t*T**2-2)**2   # eliminate sigma
P2 = t*T**3 - T + 2*t                                     # T minimal
Pmin = sp.resultant(sp.Poly(P1,T), sp.Poly(P2,T))
Pmin = sp.factor(sp.expand(Pmin))
print("resultant factors:", )
facs = sp.factor_list(Pmin)[1]
# find the factor annihilating the actual series
A_ser = None
# build A(t) series from verified formula
sq=sp.sqrt(1-t*Tser**2).series(t,0,N).removeO()
A_ser=sp.expand(sp.series((Tser/(2*t) + (sq-1+t*Tser**2/2)/(t**3))/sq, t,0,N-2).removeO())
good=[]
for fac,mult in facs:
    if F not in fac.free_symbols: continue
    val=sp.expand(sp.series(fac.subs(F,A_ser), t, 0, 8).removeO())
    if val==0: good.append(fac)
print("annihilating factors:", len(good))
P=good[0]
print("minimal polynomial degree in F:", sp.degree(P,F), " degree in t:", sp.degree(P,t))
# is it a polynomial in t^3? 
in_t3 = all((m[0]%3==0) for m in sp.Poly(P,t,F).monoms())
print("polynomial in t^3:", in_t3)
json.dump(sp.srepr(P), open('minpoly_A.json','w'))
print("P =", sp.expand(P))
