import numpy as np, time
from math import lgamma, log10, comb, factorial
from fractions import Fraction
from collections import defaultdict

# ---- digit estimate for the 5D target ----
def lg10fact(n): return lgamma(n+1)/np.log(10)
M=1000
base_k=2*(M+1)+3*M+4*M+5*M+6*M
base_k1=(M+1)+2*M+3*M+4*M+5*M
approx=lg10fact(base_k)-lg10fact(1+base_k1)-lg10fact(M+1)-4*lg10fact(M)
print("G(1000^5): leading-term size ~ 10^%.0f  => about %.0f digits"%(approx,approx))
print("  (base_k=%d, base_k1=%d)"%(base_k,base_k1))

# ---- mod-p pipeline prototype (exactly the production algorithm), validated at small M ----
def G_modp(M, p):
    Ms=[M]*4
    Imax=4*M; Wmax=(1+2+3+4)*M
    T=np.zeros((Imax+1, Wmax+1), dtype=np.int64); T[0,0]=1
    binom=[comb(M,a)%p for a in range(M+1)]
    for idx in range(4):
        wgt=idx+1
        NT=np.zeros_like(T)
        for a in range(M+1):
            b=binom[a]
            if b==0: continue
            NT[a:, wgt*a:] = (NT[a:, wgt*a:] + b*T[:Imax+1-a, :Wmax+1-wgt*a]) % p
        T=NT
    # weights
    fact=[1]*(base_of(M)+2)
    return T

def base_of(M): return 2*(M+1)+3*M+4*M+5*M+6*M

def G_exact_modp(M,p):
    Imax=4*M; Wmax=10*M
    bk=2*(M+1)+3*M+4*M+5*M+6*M
    bk1=(M+1)+2*M+3*M+4*M+5*M
    N=max(bk, bk1+1, M+1+Imax)+2
    f=[1]*(N+1)
    for i in range(1,N+1): f[i]=f[i-1]*i%p
    def inv(x): return pow(int(x),p-2,p)
    T=np.zeros((Imax+1,Wmax+1),dtype=np.int64); T[0,0]=1
    binom=[comb(M,a)%p for a in range(M+1)]
    t0=time.time()
    for idx in range(4):
        wgt=idx+1
        NT=np.zeros_like(T)
        for a in range(M+1):
            b=binom[a]
            if b==0: continue
            NT[a:, wgt*a:] = (NT[a:, wgt*a:] + b*T[:Imax+1-a, :Wmax+1-wgt*a]) % p
        T=NT
    conv_t=time.time()-t0
    # separable weights: u_i = (-1)^i i! / ((M+1+i)! * (M!)^4) ; a_W = (bk-W)!/((1+bk1-W)!)
    pref=pow(f[M],4,p)
    u=np.array([ (p-1)**(i%2) if False else 0 for i in range(Imax+1)],dtype=np.int64)
    uu=[]
    for i in range(Imax+1):
        val=f[i]*inv(f[M+1+i])%p*inv(pref)%p
        if i%2: val=(p-val)%p
        uu.append(val)
    u=np.array(uu,dtype=np.int64)
    aw=[]
    for W in range(Wmax+1):
        if bk-W<0 or 1+bk1-W<0: aw.append(0); continue
        aw.append(f[bk-W]*inv(f[1+bk1-W])%p)
    a=np.array(aw,dtype=np.int64)
    tot=0
    for i in range(Imax+1):
        row=T[i]
        if not row.any(): continue
        s=int((row*a % p).sum() % p)
        tot=(tot+u[i]*s)%p
    return tot, conv_t

# validate mod p against exact small values
def G_exact(Mv):
    M2,M3,M4,M5,M6=Mv; Ms=[M3,M4,M5,M6]
    T=defaultdict(int); T[(0,0)]=1
    for idx,Mk in enumerate(Ms):
        wgt=idx+1; T2=defaultdict(int)
        for (i,W),v in T.items():
            for aa in range(Mk+1): T2[(i+aa,W+wgt*aa)]+=v*comb(Mk,aa)
        T=T2
    bk=2*(M2+1)+3*M3+4*M4+5*M5+6*M6; bk1=(M2+1)+2*M3+3*M4+4*M5+5*M6
    pref=1
    for Mk in Ms: pref*=factorial(Mk)
    tot=Fraction(0)
    for (i,W),v in T.items():
        tot+=Fraction((-1)**i*factorial(i)*v*factorial(bk-W), factorial(1+bk1-W)*factorial(M2+1+i)*pref)
    return int(tot)

p=(1<<30)-35   # prime? check
from sympy import isprime, nextprime
p=int(nextprime(2**30))
for Mtest in [2,3,5]:
    ex=G_exact((Mtest,)*5)%p
    md,_=G_exact_modp(Mtest,p)
    print("M=%d: exact mod p = %d ; pipeline mod p = %d ; MATCH=%s"%(Mtest,ex,md,ex==md))

# timing extrapolation
for Mt in [40,80]:
    t0=time.time(); _,ct=G_exact_modp(Mt,p); tt=time.time()-t0
    print("M=%d: conv %.2fs, total %.2fs  (table %dx%d)"%(Mt,ct,tt,4*Mt+1,10*Mt+1))
