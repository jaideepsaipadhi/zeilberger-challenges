import sympy as sp
from math import factorial, comb
from collections import defaultdict
from fractions import Fraction

# ---------- 2-variable validation: t2,t3 ----------
t2,t3=sp.symbols('t2 t3')
KS2=[2,3]
def C2(m):
    top=sum(k*mk for k,mk in zip(KS2,m)); bot=1+sum((k-1)*mk for k,mk in zip(KS2,m))
    v=factorial(top)//factorial(bot)
    for mk in m: v//=factorial(mk)
    return v
NN=7
S=sp.Integer(1)
for _ in range(2*NN+4):
    S=sp.expand(1+t2*S**2+t3*S**3)
    P=sp.Poly(S,t2,t3)
    S=sum(c*t2**a*t3**b for (a,b),c in zip(P.monoms(),P.coeffs()) if a+b<=NN)
Sp=sp.Poly(sp.expand(S),t2,t3)
okC=all(int(c)==C2(list(m)) for m,c in zip(Sp.monoms(),Sp.coeffs()) if sum(m)>0)
print("2-var closed form OK:",okC)
q,r=sp.div(sp.Poly(sp.expand(S-1),t2,t3),sp.Poly(t2+t3,t2,t3))
print("2-var exact division, remainder 0:", r.as_expr()==0)
Gtrue={m:int(c) for m,c in zip(q.monoms(),q.coeffs())}
# alternating formula in 2 vars: G[M] = sum_{a3} (-1)^a3 C[M2+1+a3, M3-a3]
def G2_alt(M2,M3):
    return sum((-1)**a*C2([M2+1+a,M3-a]) for a in range(M3+1))
bad=[(m,Gtrue[m],G2_alt(*m)) for m in Gtrue if sum(m)<=NN-2 and Gtrue[m]!=G2_alt(*m)]
print("2-var alternating formula matches TRUE division:", not bad, bad[:3])

# ---------- 5-variable: recursion (independent) vs alternating vs collapsed ----------
KS=[2,3,4,5,6]
def C5(m):
    top=sum(k*mk for k,mk in zip(KS,m)); bot=1+sum((k-1)*mk for k,mk in zip(KS,m))
    v=factorial(top)//factorial(bot)
    for mk in m: v//=factorial(mk)
    return v
memo={}
def G_rec(m):
    m=tuple(m)
    if min(m)<0: return 0
    if m in memo: return memo[m]
    mp=list(m); mp[0]+=1
    val=C5(mp)
    for k in range(1,5):
        mm=list(mp); mm[k]-=1
        if min(mm)>=0: val-=G_rec(mm)
    memo[m]=val
    return val
def G_alt(Mv):
    M2,M3,M4,M5,M6=Mv; tot=0
    for a3 in range(M3+1):
     for a4 in range(M4+1):
      for a5 in range(M5+1):
       for a6 in range(M6+1):
        i=a3+a4+a5+a6
        coef=factorial(i)//(factorial(a3)*factorial(a4)*factorial(a5)*factorial(a6))
        tot+=(-1)**i*coef*C5([M2+1+i,M3-a3,M4-a4,M5-a5,M6-a6])
    return tot
def G_collapsed(Mv):
    M2,M3,M4,M5,M6=Mv; Ms=[M3,M4,M5,M6]
    T=defaultdict(int); T[(0,0)]=1
    for idx,Mk in enumerate(Ms):
        wgt=idx+1; T2=defaultdict(int)
        for (i,W),v in T.items():
            for a in range(Mk+1): T2[(i+a,W+wgt*a)]+=v*comb(Mk,a)
        T=T2
    base_k=2*(M2+1)+3*M3+4*M4+5*M5+6*M6
    base_k1=(M2+1)+2*M3+3*M4+4*M5+5*M6
    pref=1
    for Mk in Ms: pref*=factorial(Mk)
    tot=Fraction(0)
    for (i,W),v in T.items():
        tot+=Fraction((-1)**i*factorial(i)*v*factorial(base_k-W), factorial(1+base_k1-W)*factorial(M2+1+i)*pref)
    assert tot.denominator==1
    return int(tot)
print("\n5-var cross-checks (recursion is independent of the closed-form route):")
for Mv in [(1,1,1,1,1),(2,2,2,2,2),(2,1,2,1,2),(3,2,1,2,3),(1,3,2,2,1)]:
    r_,a_,c_=G_rec(Mv),G_alt(Mv),G_collapsed(Mv)
    print("  M=%s  rec=%s  alt=%s  collapsed=%s  ALL MATCH=%s"%(Mv,r_,a_,c_,r_==a_==c_))
# cost comparison at the challenge size
print("\nCOST at M=(1000,)*5:")
print("  naive alpha-sum terms      : 1001^4 = %.3e"%(1001**4))
print("  collapsed (i,W) table size : %d x %d = %.3e"%(4001,10001,4001*10001))
print("  reduction factor           : %.0f x"%(1001**4/(4001*10001)))
