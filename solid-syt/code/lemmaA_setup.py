import sympy as sp, json
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

x,t=sp.symbols('x t')
N=14
Tser=sp.Integer(0)
for _ in range(N+3):
    Tser=sp.expand(sp.series(t*(2+Tser**3), t, 0, N+3).removeO())
c_ser=sp.expand(sp.series((4*Tser-Tser**4)/8, t,0,N).removeO())  # c = t*R00 = t*(4T-T^4)/(8t)

# ---- (1) corrected axis GF ----
Dm = 1 - (1/x)*Tser*(1+Tser**3/4) + (1/x**2)*Tser**2/4
sqDm = sp.sqrt(Dm).series(t,0,N).removeO()
R_formula = sp.expand(sp.series((c_ser - x + 1/(x*t) - 1/x**2 + (x + Tser - 2/(x*Tser))*sqDm)/2, t,0,N-1).removeO())
R_dp = sp.expand(t*sum(RK(m,i,0)*x**(-i)*t**m for m in range(N) for i in range(m+1)))
d=sp.expand(R_formula-R_dp)
okR=all(sp.simplify(d.coeff(t,m))==0 for m in range(N-1))
print("CORRECTED axis GF R(xbar) verified:", okR)
if okR:
    json.dump("R(xbar) = (c - x + 1/(t x) - 1/x^2 + (x + T - 2/(xT))*sqrt(Dminus(xbar)))/2", open('axisGF.json','w'))

# ---- (2) the S pairing: S-GF(t) = CT_x[ C(t^2/x) * Q(x, t/x*C(t^2/x); t) ] ----
# Catalan GF C(z) = 1 + z C^2
z=sp.symbols('z')
Cat=sp.Integer(1)
for _ in range(N+2):
    Cat=sp.expand(sp.series(1+z*Cat**2, z, 0, (N+3)//2).removeO())
Cx = Cat.subs(z, t**2/x)   # C(t^2 xbar)
ystar = sp.expand(t/x * Cx)
# Q(x,y) from DP up to needed order
# CT_x [ Cx * Q(x, ystar) ] : expand Q(x,y)=sum RK x^i y^j t^m, substitute y->ystar
# Do it coefficient-wise in t to control blowup
SGF_target = {n: None for n in range(0, (N-1)//3+1)}
import json as js
P=js.load(open('pieces.json')); Sv=P['S']
# assemble series: sum over m,i,j RK(m,i,j) x^i ystar^j t^m ; then * Cx ; then CT_x ; then [t^{3n}]
expr=sp.Integer(0)
for m in range(N):
    for i in range(m+1):
        for j in range(m+1):
            v=RK(m,i,j)
            if v: expr += v * x**i * ystar**j * t**m
expr=sp.expand(expr*Cx)
expr=sp.expand(sp.series(expr, t, 0, N-1).removeO())
ok=True
for n in range(0,(N-2)//3+1):
    ct=expr.coeff(t,3*n).coeff(x,0) if expr.coeff(t,3*n)!=0 else 0
    # coeff(x,0) of a Laurent polynomial
    cf=expr.coeff(t,3*n)
    ct=sp.expand(cf).coeff(x,0)
    if ct!=Sv[n]: ok=False; print("pairing mismatch at n=%d: got %s want %d"%(n,ct,Sv[n])); break
print("S-GF pairing  S(n)=CT_x[C(t^2/x)*Q(x, (t/x)C(t^2/x))] at t^{3n}:", ok)
