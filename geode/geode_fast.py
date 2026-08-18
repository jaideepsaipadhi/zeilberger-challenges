"""Geode diagonal G(M,...,M) via the 1-D integral representation.
   D = number of variables (D=4: t2..t5 [Rubine]; D=5: t2..t6 [the $200 challenge]); d=D-1.
      G[M^D] = ((dM)!/(M!)^d) * [z^{dM}] (1+z)^E * int_0^1 h(x,z)^M dx
      h(x,z) = (1-x) * prod_{k=1}^{d}((1+z)^k - x)     [x-deg d+1, z-deg sd=d(d+1)/2]
      E = (sum_{k=2}^{D+1} k - sd)*M + 2
   Per prime: O(M^2) modular mults, O(M) memory."""
import numpy as np, time
from math import comb

def G_modp_vec(M, D, p):
    d=D-1; sd=d*(d+1)//2
    n=(d+1)*M+1                      # quadrature nodes: h^M has x-degree (d+1)M
    Zdeg=D*M                         # want [z^{DM}]  (a_W = ff(K', DM))
    E=(sum(range(2,D+2))-sd)*M+2
    inv=lambda a: pow(int(a)%p,p-2,p)
    xs=(np.arange(n,dtype=np.int64)+2)%p
    # ---- Pi(x)=prod(x-x_l) ----
    Pi=np.zeros(n+1,dtype=np.int64); Pi[0]=1
    for deg,xl in enumerate(xs):
        Pi[1:deg+2]=(Pi[1:deg+2]+Pi[0:deg+1])%p       # shift (multiply by x)
        Pi[0:deg+1]=(Pi[0:deg+1]-0)%p
        # subtract xl*old  -> need old before shift; redo properly:
    Pi=np.zeros(n+1,dtype=np.int64); Pi[0]=1
    for deg,xl in enumerate(xs):
        old=Pi.copy()
        Pi[:]=0
        Pi[1:deg+2]=old[0:deg+1]
        Pi[0:deg+1]=(Pi[0:deg+1]-old[0:deg+1]*int(xl))%p
        Pi%=p
    # ---- quadrature weights, all nodes at once ----
    invi=np.array([inv(i+1) for i in range(n)],dtype=np.int64)
    rem=np.full(n,int(Pi[n])%p,dtype=np.int64)
    integ=(rem*invi[n-1])%p
    for i in range(n-2,-1,-1):
        rem=(int(Pi[i+1])+rem*xs)%p
        integ=(integ+rem*invi[i])%p
    dv=np.ones(n,dtype=np.int64)
    for l in range(n):
        diff=(xs-int(xs[l]))%p; diff[l]=1
        dv=(dv*diff)%p
    Wq=(integ*np.array([inv(v) for v in dv],dtype=np.int64))%p
    # ---- g_j(z) coefficients for all nodes ----
    poly=np.zeros((1,n),dtype=np.int64); poly[0]=(1-xs)%p
    for k in range(1,d+1):
        f=[comb(k,i)%p for i in range(k+1)]
        new=np.zeros((poly.shape[0]+k,n),dtype=np.int64)
        for a in range(poly.shape[0]):
            new[a]=(new[a]+poly[a]*((1-xs)%p))%p        # b=0 term: f_0 = 1-x_j
            for b in range(1,k+1):
                new[a+b]=(new[a+b]+poly[a]*f[b])%p
        poly=new
    gc=np.zeros((sd+1,n),dtype=np.int64); gc[:poly.shape[0]]=poly
    g0=gc[0].copy()
    assert (g0!=0).all()
    ig0=np.array([inv(v) for v in g0],dtype=np.int64)
    # ---- P=g^M coefficients by the order-sd recurrence, contracted on the fly ----
    Cb=np.zeros(Zdeg+1,dtype=np.int64); Cb[0]=1
    for s in range(1,Zdeg+1): Cb[s]=Cb[s-1]*((E-s+1)%p)%p*inv(s)%p
    Pbuf=np.zeros((sd+1,n),dtype=np.int64)
    Pbuf[0]=np.array([pow(int(v),M,p) for v in g0],dtype=np.int64)
    acc=(Wq*Pbuf[0]%p*int(Cb[Zdeg]))%p
    invr=[inv(r+1) for r in range(Zdeg+1)]
    for r in range(0,Zdeg):
        s_acc=np.zeros(n,dtype=np.int64)
        for s in range(1,min(sd,r+1)+1):
            cm=(M*s-(r+1)+s)%p
            if cm==0: continue
            s_acc=(s_acc+gc[s]*cm%p*Pbuf[s-1])%p
        Pnew=s_acc*ig0%p*invr[r]%p
        Pbuf[1:]=Pbuf[:-1]; Pbuf[0]=Pnew
        acc=(acc+Wq*Pnew%p*int(Cb[Zdeg-(r+1)]))%p
    total=int(acc.sum()%p)
    fd=1
    for i in range(1,D*M+1): fd=fd*i%p
    fM=1
    for i in range(1,M+1): fM=fM*i%p
    return total*fd%p*inv(pow(fM,D,p))%p

if __name__=="__main__":
    from math import factorial
    from itertools import product
    from sympy import nextprime
    def C_gen(m,ks):
        top=sum(k*mk for k,mk in zip(ks,m)); bot=1+sum((k-1)*mk for k,mk in zip(ks,m))
        v=factorial(top)//factorial(bot)
        for mk in m: v//=factorial(mk)
        return v
    def G_exact(M,D):
        ks=list(range(2,D+2)); d=D-1; tot=0
        for al in product(range(M+1),repeat=d):
            i=sum(al); coef=factorial(i)
            for a in al: coef//=factorial(a)
            tot+=(-1)**i*coef*C_gen([M+1+i]+[M-a for a in al],ks)
        return tot
    p=int(nextprime(2**30))
    print("p =",p)
    allok=True
    for D in [4,5]:
        for M in [1,2,3,4,5]:
            ex=G_exact(M,D)%p; got=G_modp_vec(M,D,p)
            ok=ex==got; allok=allok and ok
            print("  D=%d M=%d: MATCH=%s"%(D,M,ok))
    print("ALL VALIDATIONS PASS:",allok)
    for M in [100,200,400]:
        t0=time.time(); G_modp_vec(M,5,p); print("  D=5 M=%d: %.2fs/prime"%(M,time.time()-t0))
