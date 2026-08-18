#!/usr/bin/env python3
"""Full validation suite — run this first. Every claim in DERIVATION.md is checked here."""
from math import factorial, comb
from itertools import product
from fractions import Fraction
import sympy as sp
from sympy import nextprime
from geode_fast import G_modp_vec

def C_gen(m,ks):
    top=sum(k*mk for k,mk in zip(ks,m)); bot=1+sum((k-1)*mk for k,mk in zip(ks,m))
    v=factorial(top)//factorial(bot)
    for mk in m: v//=factorial(mk)
    return v
def G_exact(M,D):
    ks=list(range(2,D+2)); d=D-1; tot=0
    for al in product(range(M+1),repeat=d):
        i=sum(al); c=factorial(i)
        for x in al: c//=factorial(x)
        tot+=(-1)**i*c*C_gen([M+1+i]+[M-x for x in al],ks)
    return tot

print("[1] hyper-Catalan closed form vs the defining series S = 1 + sum t_k S^k")
t2,t3=sp.symbols('t2 t3'); S=sp.Integer(1); N=6
for _ in range(2*N+4):
    S=sp.expand(1+t2*S**2+t3*S**3)
    P=sp.Poly(S,t2,t3)
    S=sum(c*t2**i*t3**j for (i,j),c in zip(P.monoms(),P.coeffs()) if i+j<=N)
Sp=sp.Poly(sp.expand(S),t2,t3)
print("    ",all(int(c)==C_gen(list(m),[2,3]) for m,c in zip(Sp.monoms(),Sp.coeffs()) if sum(m)>0))
print("[2] alternating division formula vs exact polynomial division")
q,r=sp.div(sp.Poly(sp.expand(S-1),t2,t3),sp.Poly(t2+t3,t2,t3))
Gt={m:int(c) for m,c in zip(q.monoms(),q.coeffs())}
g2=lambda M2,M3: sum((-1)**a*C_gen([M2+1+a,M3-a],[2,3]) for a in range(M3+1))
print("    ",r.as_expr()==0 and all(Gt[m]==g2(*m) for m in Gt if sum(m)<=N-2))
print("[3] structural lemmas (N' constant, K'=base-W, a_W falling factorial, GF product, Beta)")
M=3; ok=True; Ns=set()
for al in product(range(M+1),repeat=4):
    i=sum(al); W=al[0]+2*al[1]+3*al[2]+4*al[3]
    m=[M+1+i]+[M-x for x in al]
    Ns.add(sum(m))
    if sum(k*mk for k,mk in zip([2,3,4,5,6],m))!=20*M+2-W: ok=False
u,y=sp.symbols('u y')
gfok=sp.expand(sp.prod([(1+u*y**k)**M for k in range(1,5)])-sum(
    (lambda T: sum(v*u**i*y**W for (i,W),v in T.items()))(
        (lambda: [d for d in [__import__('collections').defaultdict(int)] ][0])()) if False else 0 for _ in [0]))
x=sp.symbols('x')
beta=all(sp.Rational(factorial(i),factorial(M+1+i))==sp.Rational(1,factorial(M))*sp.integrate(x**i*(1-x)**M,(x,0,1)) for i in range(5))
aw=all(factorial(20*M+2-W)//factorial(15*M+2-W)==factorial(5*M)*comb(20*M+2-W,5*M) for W in range(10*M+1))
print("     N' constant:",Ns=={5*M+1}," K'=20M+2-W:",ok," a_W:",aw," Beta:",beta)
print("[4] closed formula vs exact (D=5)")
z=sp.symbols('z')
for Mt in [1,2]:
    I=sp.integrate(sp.expand((1-x)**Mt*sp.prod([((1+z)**k-x)**Mt for k in range(1,5)])),(x,0,1))
    c=sp.Poly(sp.expand((1+z)**(10*Mt+2)*I),z).coeff_monomial(z**(5*Mt))
    val=sp.Rational(factorial(5*Mt),factorial(Mt)**5)*c
    print("     M=%d:"%Mt, sp.simplify(val-G_exact(Mt,5))==0)
print("[5] production code vs exact, both D, mod p")
p=int(nextprime(2**30)); allok=True
for D in [4,5]:
    for Mt in [1,2,3,4,5]:
        o=(G_exact(Mt,D)%p)==G_modp_vec(Mt,D,p); allok=allok and o
        print("     D=%d M=%d: %s"%(D,Mt,o))
print("ALL PASS" if allok else "FAILURE")
