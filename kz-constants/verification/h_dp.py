"""H(n): tandem excursions of length 3n, quarter plane, all run lengths ODD.
   State: (x, y, dir, par) with par = current-run-length parity (1=odd => may end run).
   Scaled float DP (scale 1/(1+sqrt2) per step) + brute-force validation + Richardson for C2."""
import numpy as np
from itertools import product
V=[(1,0),(-1,1),(0,-1)]
def brute(n):
    # direct enumeration over direction sequences with odd runs, length 3n, cone, excursion
    from functools import lru_cache
    N=3*n
    total=0
    def rec(x,y,d,par,steps):
        nonlocal total
        if steps==N:
            if x==0 and y==0 and par==1: total+=1
            return
        # continue current run
        nx,ny=x+V[d][0],y+V[d][1]
        if nx>=0 and ny>=0: rec(nx,ny,d,1-par,steps+1)
        # or end run (allowed if par odd) and switch
        if par==1:
            for nd in range(3):
                if nd!=d:
                    nx,ny=x+V[nd][0],y+V[nd][1]
                    if nx>=0 and ny>=0: rec(nx,ny,nd,1,steps+1)
    for d0 in range(3):
        x,y=V[d0]
        if x>=0 and y>=0: rec(x,y,d0,1,1)
    return total
def dp(n,scale):
    N=3*n; L=N+2
    # arr[dir,par,x,y]
    arr=np.zeros((3,2,L,L))
    for d in range(3):
        x,y=V[d]
        if x>=0 and y>=0: arr[d,1,x,y]=scale
    for step in range(2,N+1):
        new=np.zeros_like(arr)
        for d in range(3):
            dx,dy=V[d]
            # continue: from (d,par) -> (d,1-par)
            for par in range(2):
                src=arr[d,par]
                sh=np.zeros_like(src)
                xs=slice(max(0,dx),L if dx>=0 else L+dx)
                # simpler: use np.roll with zeroing
                t=np.roll(np.roll(src,dx,axis=0),dy,axis=1)
                if dx==1: t[0,:]=0
                if dx==-1: t[-1,:]=0
                if dy==1: t[:,0]=0
                if dy==-1: t[:,-1]=0
                new[d,1-par]+=t
            # switch: from (d2,1) -> (d,1) for d2 != d
            for d2 in range(3):
                if d2==d: continue
                src=arr[d2,1]
                t=np.roll(np.roll(src,dx,axis=0),dy,axis=1)
                if dx==1: t[0,:]=0
                if dx==-1: t[-1,:]=0
                if dy==1: t[:,0]=0
                if dy==-1: t[:,-1]=0
                new[d,1]+=t
        # cone: negative coords impossible by construction (rolls zeroed); but v1=(-1,1): x can't go <0 handled
        arr=new*scale
    return arr[:,1,0,0].sum()
# validate small n vs brute (unscaled)
ok=True
for n in range(1,5):
    b=brute(n)
    d=dp(n,1.0)
    ok &= abs(d-b)<1e-6
    print("n=%d brute=%d dp=%.1f"%(n,b,d))
print("VALIDATION:","PASS" if ok else "FAIL")
mu=7+5*np.sqrt(2)
sc=1/(1+np.sqrt(2))   # per-step scale => per n: mu scaled out exactly
vals={}
import sys
for n in range(2,97,2):
    vals[n]=dp(n,sc)     # = H(n) * mu^{-n}
# c(n) = H(n) mu^{-n} n^4 -> C2; Richardson in 1/n
ns=sorted(vals)
c=[vals[n]*n**4 for n in ns]
# successive Richardson (orders 1..4) on the tail
seq=list(zip(ns,c))
for order in range(1,5):
    seq=[(n2,(c2*n2-c1*n1)/(n2-n1)) if order==1 else (n2, c2 + (c2-c1)*n1/(n2-n1)*1) for (n1,c1),(n2,c2) in zip(seq,seq[1:])]
print("last raw c(n):", ["%.8f"%x for _,x in list(zip(ns,c))[-3:]])
# clean Richardson: iterate x_k -> (n_{k+1} x_{k+1} - n_k x_k)/(n_{k+1}-n_k)
def rich(pairs):
    return [ (b_n, (b_n*b_v - a_n*a_v)/(b_n - a_n)) for (a_n,a_v),(b_n,b_v) in zip(pairs,pairs[1:]) ]
p=list(zip(ns,c))
for it in range(4):
    p=rich(p)
    print("Richardson pass %d tail: %.10f  %.10f"%(it+1,p[-2][1],p[-1][1]))
print("C2 estimate: %.10f"%p[-1][1])
