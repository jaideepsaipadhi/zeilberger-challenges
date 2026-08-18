from math import factorial, comb
from itertools import product
from collections import defaultdict
KS=[2,3,4,5,6]

def C_closed(m):   # m: list of m_k aligned with KS
    top=sum(k*mk for k,mk in zip(KS,m))
    bot=1+sum((k-1)*mk for k,mk in zip(KS,m))
    v=factorial(top)//factorial(bot)
    for mk in m: v//=factorial(mk)
    return v

# --- verify closed form against series via dict-based truncated expansion ---
N=2
def mul(A,B):
    R=defaultdict(int)
    for ma,ca in A.items():
        for mb,cb in B.items():
            m=tuple(x+y for x,y in zip(ma,mb))
            if max(m)<=N: R[m]+=ca*cb
    return R
S={(0,0,0,0,0):1}
for it in range(14):
    # S = 1 + sum_k t_k S^k
    new=defaultdict(int); new[(0,0,0,0,0)]+=1
    P={(0,0,0,0,0):1}
    for i,k in enumerate(KS):
        # S^k
        Sk={(0,0,0,0,0):1}
        for _ in range(k): Sk=mul(Sk,S)
        for m,c in Sk.items():
            mm=list(m); mm[i]+=1; mm=tuple(mm)
            if max(mm)<=N: new[mm]+=c
    S=dict(new)
ok=True
for m,c in S.items():
    if sum(m)==0: continue
    if c!=C_closed(list(m)): ok=False; print("MISMATCH",m,c,C_closed(list(m)));break
print("hyper-Catalan closed form verified vs series (exponents<=%d):"%N, ok)

# --- G coefficients from S via the defining relation C[m] = sum_k G[m-e_k] ---
# compute G by solving recursively over the box
G={}
def getG(m):
    m=tuple(m)
    if min(m)<0: return 0
    if m in G: return G[m]
    if sum(m)==0: 
        G[m]=0; return 0
    # C[m] = sum_k G[m-e_k]  => but we need an ordering; use lexicographic solve on t2 direction:
    raise RuntimeError
# simpler: build G directly by dividing series: G = (S-1)/(t2+...+t6) using ordered division
def divide(S):
    # returns G dict with (t2+..+t6)*G = S-1, on the truncated box
    Sm={m:c for m,c in S.items() if sum(m)>0}
    Gd=defaultdict(int)
    for m in sorted(Sm.keys(), key=lambda z:(sum(z),z)):
        # coefficient of m in (t2+...)*G = sum_k G[m-e_k]
        acc=0
        for i in range(5):
            mm=list(m); mm[i]-=1
            if min(mm)>=0: acc+=Gd[tuple(mm)]
        need=Sm[m]-acc
        # assign the deficit to the lexicographically-first available predecessor... 
        # correct approach: G[m - e_2] gets it if consistent
        mm=list(m); mm[0]-=1
        if min(mm)>=0: Gd[tuple(mm)]+=need
        else:
            if need!=0: print("division inconsistency at",m,need)
    return Gd
Gd=divide(S)
# verify (t2+...+t6)*G == S-1 on the box
bad=0
for m,c in S.items():
    if sum(m)==0: continue
    acc=0
    for i in range(5):
        mm=list(m); mm[i]-=1
        if min(mm)>=0: acc+=Gd[tuple(mm)]
    if acc!=c: bad+=1
print("division verified on box:", bad==0)

# --- alternating-sum formula ---
def G_alt(Mv):
    M2,M3,M4,M5,M6=Mv
    tot=0
    for a3 in range(M3+1):
     for a4 in range(M4+1):
      for a5 in range(M5+1):
       for a6 in range(M6+1):
        i=a3+a4+a5+a6
        m=[M2+1+i,M3-a3,M4-a4,M5-a5,M6-a6]
        coef=factorial(i)//(factorial(a3)*factorial(a4)*factorial(a5)*factorial(a6))
        tot+=(-1)**i*coef*C_closed(m)
    return tot

# --- (i,W)-collapsed formula ---
def G_collapsed(Mv):
    M2,M3,M4,M5,M6=Mv
    Ms=[M3,M4,M5,M6]
    T=defaultdict(int); T[(0,0)]=1
    for idx,Mk in enumerate(Ms):
        wgt=idx+1
        T2=defaultdict(int)
        for (i,W),v in T.items():
            for a in range(Mk+1):
                T2[(i+a,W+wgt*a)]+=v*comb(Mk,a)
        T=T2
    base_k =2*(M2+1)+3*M3+4*M4+5*M5+6*M6
    base_k1=(M2+1)+2*M3+3*M4+4*M5+5*M6
    pref=1
    for Mk in Ms: pref*=factorial(Mk)
    tot=0
    for (i,W),v in T.items():
        sum_k  = base_k  + 2*i - (W+2*i)
        sum_k1 = base_k1 +   i - (W+i)
        term=factorial(i)*v*factorial(sum_k)
        den=factorial(1+sum_k1)*factorial(M2+1+i)*pref
        assert term%den==0 or True
        tot+= (-1)**i * term / den if False else 0
    # do it with exact rationals
    from fractions import Fraction
    tot=Fraction(0)
    for (i,W),v in T.items():
        sum_k  = base_k  - W
        sum_k1 = base_k1 - W
        tot+= Fraction((-1)**i*factorial(i)*v*factorial(sum_k), factorial(1+sum_k1)*factorial(M2+1+i)*pref)
    assert tot.denominator==1, "not integer: %s"%tot
    return int(tot)

for Mv in [(1,1,1,1,1),(2,2,2,2,2),(2,1,2,1,2),(1,2,2,1,1)]:
    if max(Mv)<=N:
        gs=Gd[tuple(Mv)]
    else: gs=None
    ga=G_alt(Mv); gc=G_collapsed(Mv)
    print("M=%s  series=%s  alt=%s  collapsed=%s  alt==collapsed:%s  series==alt:%s"%(Mv,gs,ga,gc,ga==gc, (gs==ga) if gs is not None else "n/a"))
