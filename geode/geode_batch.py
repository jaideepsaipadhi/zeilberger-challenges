"""Batched CPU backend: processes B primes simultaneously as (B,n) arrays.
   Speedups vs geode_fast.py:
     - primes batched into the vector dimension (amortizes numpy per-op overhead)
     - chunked accumulation: products < 2^62 accumulated 4 at a time -> ~2x fewer '%' ops
     - closed form for the Lagrange denominators at equispaced nodes:
         prod_{l!=j}(x_j - x_l) = (-1)^{n-1-j} * j! * (n-1-j)!   (x_j = j+2)
   Exactness guards: REQUIRE p < 2^30.5 = 1518500249: each product coef*P < p^2 < 2^61,
   and chunks of 4 such products stay < 2^63. Asserted at entry; primes from ~2^30 comply."""
import numpy as np

def inv_table(m, p):
    """inverses of 1..m mod p via inv[i] = -(p//i)*inv[p%i]"""
    inv=np.zeros(m+1,dtype=np.int64); inv[1]=1
    for i in range(2,m+1):
        inv[i]=(-(p//i)*inv[p%i])%p
    return inv[1:]

def G_modp_batch(M, D, primes):
    primes=np.asarray(primes,dtype=np.int64)
    assert (primes<1518500249).all() and (primes>2).all(), 'need p < 2^30.5 for int64 chunk safety'
    B=len(primes); P=primes.reshape(B,1)
    d=D-1; sd=d*(d+1)//2; n=D*M+1; Zdeg=D*M
    E=(sum(range(2,D+2))-sd)*M+2
    xs=np.arange(2,n+2,dtype=np.int64)          # nodes; < p for all our primes
    # ---- per-prime small tables ----
    invs=np.stack([inv_table(max(n,Zdeg+1,E%int(p)+1 if False else n), int(p)) for p in primes])  # inverses 1..n
    # need inverses of 1..max(n, Zdeg+1); n = Zdeg+1 so n suffices
    # factorials mod p up to n-1 for the denominator closed form
    fact=np.ones((B,n),dtype=np.int64)
    for i in range(1,n):
        fact[:,i]=fact[:,i-1]*i % P[:,0].reshape(B) if False else (fact[:,i-1]*i)%primes
    # ---- Pi(x) = prod_l (x - x_l) mod p, shape (B, n+1) ----
    Pi=np.zeros((B,n+1),dtype=np.int64); Pi[:,0]=1
    for deg in range(n):
        xl=int(xs[deg])
        old=Pi[:,:deg+1].copy()
        Pi[:,1:deg+2]=(Pi[:,1:deg+2]+old)%P if False else (Pi[:,1:deg+2]*0 + 0)
        # correct two-step update:
        Pi[:,1:deg+2]=old
        Pi[:,0:deg+1]=(-old*xl)%P
        if deg>0:
            pass
    # the loop above is wrong when deg>0 (loses previous low terms); rebuild properly:
    Pi=np.zeros((B,n+1),dtype=np.int64); Pi[:,0]=1
    for deg in range(n):
        xl=int(xs[deg])
        old=Pi[:,:deg+2].copy()
        upd=np.zeros_like(old)
        upd[:,1:]=old[:,:-1]
        upd=(upd-old*xl)%P
        Pi[:,:deg+2]=upd
    # ---- quadrature weights: integ_j = int_0^1 Pi(x)/(x-x_j) dx via synthetic division ----
    invi=invs                                  # (B, n): inverses of 1..n
    rem=np.repeat(Pi[:,n:n+1],n,axis=1)%P      # (B, n) start: leading coeff
    integ=(rem*invi[:,n-1:n])%P
    for i in range(n-2,-1,-1):
        rem=(Pi[:,i+1:i+2]+rem*xs)%P
        integ=(integ+rem*invi[:,i:i+1])%P
    # denominators: prod_{l!=j}(x_j-x_l) = (-1)^{n-1-j} j! (n-1-j)!
    j=np.arange(n)
    dv=(fact[:,j]*fact[:,n-1-j])%P
    sign=np.where((n-1-j)%2==1, -1, 1)
    dv=(dv*sign)%P
    # modular inverse of dv via Fermat (vectorized pow not available -> per element; use
    # dv^{-1} = dv^{p-2}: do with python pow per prime over n elems is 5e3*B... instead
    # invert via building inverse from invs of factorials: (j!)^{-1} = invfact
    invfact=np.ones((B,n),dtype=np.int64)
    for i in range(1,n):
        invfact[:,i]=(invfact[:,i-1]*invs[:,i-1])%primes
    idv=(invfact[:,j]*invfact[:,n-1-j])%P
    idv=(idv*sign)%P
    Wq=(integ*idv)%P
    # ---- g coefficients per node, per prime: (B, sd+1, n) ----
    from math import comb
    poly=np.zeros((B,1,n),dtype=np.int64); poly[:,0,:]=(1-xs)%P
    for k in range(1,d+1):
        f=[comb(k,i) for i in range(k+1)]
        new=np.zeros((B,poly.shape[1]+k,n),dtype=np.int64)
        for a in range(poly.shape[1]):
            new[:,a,:]=(new[:,a,:]+poly[:,a,:]*((1-xs)%P))%P
            for b in range(1,k+1):
                new[:,a+b,:]=(new[:,a+b,:]+poly[:,a,:]*f[b])%P
        poly=new
    gc=np.zeros((B,sd+1,n),dtype=np.int64); gc[:,:poly.shape[1],:]=poly
    g0=gc[:,0,:].copy()
    # ig0 via Fermat with per-prime pow (python loop over B x n... use pow on ints, B*n calls)
    ig0=np.zeros_like(g0)
    for bi in range(B):
        p=int(primes[bi])
        ig0[bi]=np.array([pow(int(v),p-2,p) for v in g0[bi]],dtype=np.int64)
    P0=np.zeros_like(g0)
    for bi in range(B):
        p=int(primes[bi])
        P0[bi]=np.array([pow(int(v),M,p) for v in g0[bi]],dtype=np.int64)
    # ---- Cb (B, Zdeg+1) ----
    Cb=np.ones((B,Zdeg+1),dtype=np.int64)
    for s in range(1,Zdeg+1):
        Cb[:,s]=(Cb[:,s-1]*((E-s+1)%primes))%primes
        Cb[:,s]=(Cb[:,s]*invs[:,s-1])%primes
    # ---- the march, chunked accumulation ----
    Pbuf=np.zeros((sd+1,B,n),dtype=np.int64); Pbuf[0]=P0
    acc=(Wq*P0)%P
    acc=(acc*Cb[:,Zdeg:Zdeg+1])%P
    invr=invs                                  # inverses of 1..n cover r+1 <= Zdeg = n-1... n-1+1=n OK
    for r in range(Zdeg):
        s_acc=np.zeros((B,n),dtype=np.int64)
        chunk=np.zeros((B,n),dtype=np.int64); c_in=0
        smax=min(sd,r+1)
        for s in range(1,smax+1):
            cm=(M*s-(r+1)+s)
            coef=(gc[:,s,:]*cm)%P              # < 2^31 * 4e4 < 2^47, one mod
            chunk=chunk+coef*Pbuf[s-1]         # each term < 2^62/4-ish; accumulate 4
            c_in+=1
            if c_in==4:
                s_acc=(s_acc+chunk)%P; chunk[:]=0; c_in=0
        if c_in: s_acc=(s_acc+chunk)%P
        Pn=(s_acc*ig0)%P
        Pn=(Pn*invr[:,r:r+1])%P
        Pbuf[1:]=Pbuf[:-1]; Pbuf[0]=Pn
        t=(Wq*Pn)%P
        acc=(acc+t*Cb[:,Zdeg-(r+1):Zdeg-r])%P
    tot=acc.sum(axis=1)%primes
    # prefactor
    out=[]
    for bi in range(B):
        p=int(primes[bi]); fd=1
        for i in range(1,D*M+1): fd=fd*i%p
        fM=1
        for i in range(1,M+1): fM=fM*i%p
        out.append(int(tot[bi])*fd%p*pow(pow(fM,D,p),p-2,p)%p)
    return out
