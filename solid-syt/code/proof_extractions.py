"""Pin down the exact Wiener-Hopf splitting behind Prop Qd and Prop R, every sign machine-checked."""
import sympy as sp
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)
@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)
x,t=sp.symbols('x t')
N=13
T=sp.Integer(0)
for _ in range(N+3): T=sp.expand(sp.series(t*(2+T**3),t,0,N+3).removeO())
c=sp.expand(sp.series((4*T-T**4)/8,t,0,N).removeO())
Qd=sum(RK(m,i,i)*x**(-i)*t**m for m in range(N) for i in range(m//2+1))      # Qd(xbar;t)
R=sp.expand(t*sum(RK(m,i,0)*x**i*t**m for m in range(N) for i in range(m+1)))  # R(x)
Rb=sp.expand(t*sum(RK(m,i,0)*x**(-i)*t**m for m in range(N) for i in range(m+1)))  # R(xbar)
Delta=(1-t*x)**2-4*t**2/x
sqD=sp.sqrt(sp.expand(Delta*x)/x)
def ser(e,n=N-1): return sp.expand(sp.series(e,t,0,n).removeO())
sqDs=ser(sp.sqrt(1-2*t*x+t**2*x**2-4*t**2/x))
# (dagger): sqrt(D)*(x/t - xbar*Qd(xbar)) = c - xbar + x/t - x^2 - 2R(x)
lhs=ser(sqDs*(x/t - Qd/x)); rhs=ser(c - 1/x + x/t - x**2 - 2*R)
d=sp.expand(lhs-rhs)
print("(dagger) with R(x):", all(sp.simplify(d.coeff(t,m))==0 for m in range(-1,N-2)))
# now x -> xbar version: sqrt(D(xbar))*(xbar/t - x*Qd(x)) = c - x + xbar/t - xbar^2 - 2R(xbar)
Qdx=sum(RK(m,i,i)*x**(i)*t**m for m in range(N) for i in range(m//2+1))
sqDb=ser(sp.sqrt(1-2*t/x+t**2/x**2-4*t**2*x))
lhs2=ser(sqDb*(1/(x*t) - x*Qdx)); rhs2=ser(c - x + 1/(x*t) - 1/x**2 - 2*Rb)
d2=sp.expand(lhs2-rhs2)
print("(dagger-bar) with R(xbar):", all(sp.simplify(d2.coeff(t,m))==0 for m in range(-1,N-2)))
# factorization of D(xbar): candidates
D0=4*t**2/T**2
Dp=1-x*T**2
Dm=1-(1/x)*T*(1+T**3/4)+(1/x**2)*T**2/4
prod=ser(sp.expand(D0*Dp*Dm))
Db=sp.expand(1-2*t/x+t**2/x**2-4*t**2*x)
print("D(xbar) = D0*Dplus(x)*Dminus(xbar):", sp.expand(ser(prod)-ser(Db))==0)
# THE SPLITTING: divide (dagger-bar) by sqrt(D0*Dminus(xbar)):
sqDm=ser(sp.sqrt(Dm)); sqDp=ser(sp.sqrt(1-x*T**2)); sq0=ser(sp.sqrt(D0))
LHS3=ser(sqDp*(1/(x*t) - x*Qdx))            # sqrt(Dplus(x)) * (xbar/t - x Qd(x))
RHS3=ser(sp.expand((c - x + 1/(x*t) - 1/x**2 - 2*Rb))/(sq0*sqDm))
print("split identity LHS3=RHS3:", all(sp.simplify(sp.expand(LHS3-RHS3).coeff(t,m))==0 for m in range(-1,N-3)))
# classify supports: which x-powers can each unknown block occupy?
def _unused(e,n):
    e=sp.expand(e); pows=set()
    for m in range(-1,n):
        cm=sp.expand(e.coeff(t,m))
        p=sp.Poly(cm.subs(x,sp.Symbol('X')), sp.Symbol('X'), 1/sp.Symbol('X')) if cm!=0 else None
    return None
# unknown on LHS3: sqrt(Dplus)*x*Qd(x): powers >= +1? check min power
U1=ser(sqDp*x*Qdx)
mins=set(); 
for m in range(0,N-3):
    cm=sp.expand(U1.coeff(t,m))
    if cm==0: continue
    pl=sp.Poly(cm,x)
    mins.add(min(k[0] for k in pl.monoms()))
print("x-support min of sqrt(Dplus)*x*Qd(x):", min(mins), "(expect >=1: unknown-Qd block strictly positive)")
# unknown on RHS3: 2R(xbar)/(sqrt(D0*Dm)): powers <= 0?
U2=ser(sp.expand(2*Rb/(sq0*sqDm)))
maxs=set()
for m in range(-1,N-3):
    cm=sp.expand(U2.coeff(t,m))
    if cm==0: continue
    # laurent in x
    fr=sp.fraction(sp.together(cm)); num,den=fr
    pn=sp.Poly(num,x); pd=sp.Poly(den,x)
    maxs.add(pn.degree()-pd.degree())
print("x-support max of 2R(xbar)/sqrt(D0*Dm):", max(maxs), "(expect <=0: unknown-R block nonpositive)")
# known part: k(x) := (c - x + xbar/t - xbar^2)/(sqrt(D0*Dm)) + sqrt(Dplus)*xbar/t : compute its >0 and <=0 parts
KN=ser(sqDp/(x*t) - sp.expand((c - x + 1/(x*t) - 1/x**2)/(sq0*sqDm)))
def pos_part(e,n):
    out=0
    for m in range(-1,n):
        cm=sp.expand(e.coeff(t,m))
        if cm==0: continue
        pl=sp.Poly(sp.expand(cm*x**40), x)   # shift to make polynomial
        for (k,),co in zip(pl.monoms(),pl.coeffs()):
            if k-40>0: out+=co*x**(k-40)*t**m
    return sp.expand(out)
POSK=pos_part(KN,N-3)
POSL=pos_part(ser(sqDp*x*Qdx),N-3)
print("positive-part balance: [x Qd sqrt(Dplus)]^{>0} == [known]^{>0}:", sp.expand(POSL-POSK)==0)
# and derive Qd formula shape: xQd*sqrt(Dplus) is ITSELF all-positive => equals POSK => Qd = POSK/(x*sqrt(Dplus))
Qd_derived=ser(sp.expand(POSK/(x*sqDp)),N-4)
print("derived Qd == DP Qd(x):", sp.expand(Qd_derived-ser(Qdx,N-4))==0)

# ---- Prop R from the SAME splitting via the nonpositive part ----
def neg_part(e,n):
    out=0
    for m in range(-1,n):
        cm=sp.expand(e.coeff(t,m))
        if cm==0: continue
        pl=sp.Poly(sp.expand(cm*x**40), x)
        for (k,),co in zip(pl.monoms(),pl.coeffs()):
            if k-40<=0: out+=co*x**(k-40)*t**m
    return sp.expand(out)
A_known_neg = neg_part(ser(sp.expand((c - x + 1/(x*t) - 1/x**2)/(sq0*sqDm))), N-3)
B_known_neg = neg_part(ser(sqDp/(x*t)), N-3)
R_derived = ser(sp.expand((sq0*sqDm/2)*(A_known_neg - B_known_neg)), N-4)
print("Prop R from the <=0 part of the SAME splitting == DP R(xbar):",
      all(sp.simplify(sp.expand(R_derived-ser(Rb,N-4)).coeff(t,m))==0 for m in range(N-4)))
