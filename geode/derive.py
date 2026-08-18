import sympy as sp
from math import factorial, comb
from fractions import Fraction
from collections import defaultdict
x,z=sp.symbols('x z')

# ---------- reference: verified exact value ----------
KS=[2,3,4,5,6]
def C5(m):
    top=sum(k*mk for k,mk in zip(KS,m)); bot=1+sum((k-1)*mk for k,mk in zip(KS,m))
    v=factorial(top)//factorial(bot)
    for mk in m: v//=factorial(mk)
    return v
def G_exact(M):
    Mv=(M,)*5
    tot=0
    for a3 in range(M+1):
     for a4 in range(M+1):
      for a5 in range(M+1):
       for a6 in range(M+1):
        i=a3+a4+a5+a6
        coef=factorial(i)//(factorial(a3)*factorial(a4)*factorial(a5)*factorial(a6))
        tot+=(-1)**i*coef*C5([M+1+i,M-a3,M-a4,M-a5,M-a6])
    return tot

# ---------- STEP-BY-STEP CHECKS of the new derivation ----------
M=3
# (a) N' = sum of shifted m_k is CONSTANT = 5M+1, independent of alpha
Ns=set()
for a3 in range(M+1):
 for a4 in range(M+1):
  for a5 in range(M+1):
   for a6 in range(M+1):
    i=a3+a4+a5+a6
    Ns.add((M+1+i)+(M-a3)+(M-a4)+(M-a5)+(M-a6))
print("(a) N' constant?", Ns, " expected {%d}"%(5*M+1))
# (b) K' = 20M+2 - W
okK=True
for a3 in range(M+1):
 for a4 in range(M+1):
  for a5 in range(M+1):
   for a6 in range(M+1):
    i=a3+a4+a5+a6; W=a3+2*a4+3*a5+4*a6
    m=[M+1+i,M-a3,M-a4,M-a5,M-a6]
    if sum(k*mk for k,mk in zip(KS,m))!=20*M+2-W: okK=False
print("(b) K' = 20M+2-W identically:", okK)
# (c) a_W = (20M+2-W)!/(15M+2-W)! = (5M)! * C(20M+2-W, 5M)
okA=all(factorial(20*M+2-W)//factorial(15*M+2-W)==factorial(5*M)*comb(20*M+2-W,5*M) for W in range(0,10*M+1))
print("(c) a_W = (5M)!*C(20M+2-W,5M):", okA)
# (d) generating function identity: sum_{i,W} T[i,W] u^i y^W = prod_{k=1..4} (1+u y^k)^M
u,y=sp.symbols('u y')
T=defaultdict(int); T[(0,0)]=1
for idx in range(4):
    wgt=idx+1; T2=defaultdict(int)
    for (i,W),v in T.items():
        for a in range(M+1): T2[(i+a,W+wgt*a)]+=v*comb(M,a)
    T=T2
lhs=sum(v*u**i*y**W for (i,W),v in T.items())
rhs=sp.expand(sp.prod([(1+u*y**k)**M for k in range(1,5)]))
print("(d) T-generating function = prod (1+u y^k)^M:", sp.expand(lhs-rhs)==0)
# (e) Beta identity  i!/(M+1+i)! = (1/M!) int_0^1 x^i (1-x)^M dx
okB=all(sp.Rational(factorial(i),factorial(M+1+i))==sp.Rational(1,factorial(M))*sp.integrate(x**i*(1-x)**M,(x,0,1)) for i in range(0,6))
print("(e) Beta identity:", okB)

# ---------- THE CLOSED FORMULA ----------
def G_formula(M):
    integrand=(1-x)**M*sp.prod([((1+z)**k-x)**M for k in range(1,5)])
    I=sp.integrate(sp.expand(integrand),(x,0,1))
    F=sp.expand((1+z)**(10*M+2)*I)
    c=sp.Poly(F,z).coeff_monomial(z**(5*M))
    return sp.Rational(factorial(5*M),factorial(M)**5)*c

for Mt in [1,2,3]:
    gf=G_formula(Mt); ge=G_exact(Mt)
    print("M=%d: formula=%s  exact=%s  MATCH=%s"%(Mt,gf,ge,sp.simplify(gf-ge)==0))
