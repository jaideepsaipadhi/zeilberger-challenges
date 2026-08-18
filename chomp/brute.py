import sys
from functools import lru_cache
def solve(A,B):
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def wincnt(lam):
        c=0
        for i in range(len(lam)):
            if lam[i]==0: break
            lo=1 if i==0 else 0
            for t in range(lo,lam[i]):
                mu=list(lam)
                for k in range(i,len(lam)): mu[k]=min(mu[k],t)
                if wincnt(tuple(mu))==0: c+=1
        return c
    res={}
    for a2 in range(1,A+1):
        for b2 in range(1,B+1):
            res[(a2,b2)]=wincnt(tuple([b2]*a2+[0]*(A-a2)))
    # P count over the whole box
    allp=0
    def gen(k,mx,cur):
        nonlocal allp
        if k==A:
            if wincnt(tuple(cur))==0: allp+=1
            return
        for v in range(0,mx+1): gen(k+1,v,cur+[v])
    gen(0,B,[])
    return res,allp
A,B=int(sys.argv[1]),int(sys.argv[2])
sys.setrecursionlimit(100000)
res,allp=solve(A,B)
for (a2,b2),c in sorted(res.items()):
    print("BRUTE %d x %d : %d"%(a2,b2,c))
print("BRUTE P-positions:",allp)
