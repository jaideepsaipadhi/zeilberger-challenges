import numpy as np, time
from math import factorial, comb
from sympy import nextprime

# G[M] = ((5M)!/(M!)^5) * [z^{5M}] (1+z)^{10M+2} * int_0^1 h(x,z)^M dx,
#   h(x,z) = (1-x) * prod_{k=1..4} ((1+z)^k - x)   [degree 5 in x, degree 10 in z]
# ALGORITHM (per prime p):
#   - exact interpolatory quadrature in x with 5M+1 nodes  -> int = sum_j w_j h(x_j,z)^M
#   - for each node j: g_j(z) := h(x_j,z) is a degree-10 polynomial; its M-th power's
#     coefficients p_r satisfy an order-10 linear recurrence (from g P' = M g' P) -> O(10) per coeff
#   - [z^{5M}] of (1+z)^{10M+2} * sum_j w_j g_j^M  =  sum_j w_j sum_r C(10M+2, 5M-r) p_r^{(j)}
# Cost: O(M^2) modular mults per prime, O(M) memory.  (vs O(M^3) time, O(M^2) memory for the 2D table)

def G_modp(M,p):
    n=5*M+1                       # number of quadrature nodes
    xs=[(j+2)%p for j in range(n)]  # avoid x=1 (would make g_j(0)=0)
    # ---- quadrature weights: w_j = int_0^1 prod_{l!=j}(x-x_l)/(x_j-x_l) dx ----
    Pi=np.zeros(n+1,dtype=object); Pi[0]=1; deg=0
    for xl in xs:                                   # build prod (x - x_l)
        new=np.zeros(n+1,dtype=object)
        for i in range(deg+1):
            new[i+1]=(new[i+1]+Pi[i])%p
            new[i]=(new[i]-Pi[i]*xl)%p
        Pi=new; deg+=1
    inv=lambda a:pow(int(a)%p,p-2,p)
    invj=[inv(j+1) for j in range(n+1)]             # 1/(i+1) for integration
    W=[]
    for j,xj in enumerate(xs):
        # synthetic division Pi/(x-xj) -> coefficients q (degree n-1), and Pi'(xj)=prod_{l!=j}(xj-xl)
        q=np.zeros(n,dtype=object); rem=0
        for i in range(n-1,-1,-1):                  # descending synthetic division
            rem=(Pi[i+1]+rem*xj)%p if i<n-1 else Pi[n]%p
            q[i]=rem
        # derivative value
        d=1
        for l,xl in enumerate(xs):
            if l!=j: d=d*((xj-xl)%p)%p
        integ=0
        for i in range(n): integ=(integ+q[i]*invj[i])%p
        W.append(integ*inv(d)%p)
    # ---- per node: g_j(z) coefficients (degree 10), then M-th power by recurrence ----
    C_top=[0]*(5*M+1)                                # C(10M+2, 5M-r)
    # compute binomials C(10M+2, s) for s=0..5M
    Cb=[0]*(5*M+1); Cb[0]=1
    for s in range(1,5*M+1): Cb[s]=Cb[s-1]*((10*M+2-s+1)%p)%p*inv(s)%p
    total=0
    for j,xj in enumerate(xs):
        # g = (1-x) * prod_k ((1+z)^k - x)  as polynomial in z
        g=np.zeros(11,dtype=object); g[0]=(1-xj)%p
        cur=[ (1-xj)%p ]                             # start scalar
        poly=[(1-xj)%p]
        for k in range(1,5):
            # (1+z)^k - xj  coefficients
            f=[comb(k,i)%p for i in range(k+1)]; f[0]=(f[0]-xj)%p
            new=[0]*(len(poly)+k)
            for a,ca in enumerate(poly):
                if ca==0: continue
                for b,cb in enumerate(f): new[a+b]=(new[a+b]+ca*cb)%p
            poly=new
        g=poly+[0]*(11-len(poly))
        if g[0]==0: raise RuntimeError("g0=0 at node %d"%j)
        ig0=inv(g[0])
        # p_r coefficients of g^M via  g P' = M g' P
        P=[0]*(5*M+1); P[0]=pow(int(g[0]),M,p)
        for r in range(0,5*M):
            acc=0
            for s in range(1,11):
                if r+1-s<0: break
                if g[s]==0: continue
                acc=(acc+g[s]*((M*s-(r+1)+s)%p)%p*P[r+1-s])%p
            P[r+1]=acc%p*ig0%p*inv(r+1)%p
        sj=0
        for r in range(5*M+1): sj=(sj+Cb[5*M-r]*P[r])%p
        total=(total+W[j]*sj)%p
    # prefactor (5M)!/(M!)^5
    f5=1
    for i in range(1,5*M+1): f5=f5*i%p
    fM=1
    for i in range(1,M+1): fM=fM*i%p
    pref=f5*inv(pow(fM,5,p))%p
    return total*pref%p

# ---- VALIDATION against exact values ----
KS=[2,3,4,5,6]
def C5(m):
    top=sum(k*mk for k,mk in zip(KS,m)); bot=1+sum((k-1)*mk for k,mk in zip(KS,m))
    v=factorial(top)//factorial(bot)
    for mk in m: v//=factorial(mk)
    return v
def G_exact(M):
    tot=0
    for a3 in range(M+1):
     for a4 in range(M+1):
      for a5 in range(M+1):
       for a6 in range(M+1):
        i=a3+a4+a5+a6
        coef=factorial(i)//(factorial(a3)*factorial(a4)*factorial(a5)*factorial(a6))
        tot+=(-1)**i*coef*C5([M+1+i,M-a3,M-a4,M-a5,M-a6])
    return tot
p=int(nextprime(2**61))
print("prime p =",p)
for M in [1,2,3,4,5]:
    ex=G_exact(M)%p; got=G_modp(M,p)
    print("M=%d: exact mod p=%d  fast mod p=%d  MATCH=%s"%(M,ex,got,ex==got))
for M in [30,60]:
    t0=time.time(); G_modp(M,p); print("M=%d: %.2fs per prime"%(M,time.time()-t0))
