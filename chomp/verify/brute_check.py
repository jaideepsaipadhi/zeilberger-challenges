import sys
from functools import lru_cache
sys.setrecursionlimit(10000)
@lru_cache(None)
def P(l):
    if l[0]==1 and all(x==0 for x in l[1:]): return True
    for r in range(len(l)):
        for c in range(l[r]):
            if r==0 and c==0: continue
            ch=tuple(min(x,c) if k>=r else x for k,x in enumerate(l))
            if P(ch): return False
    return True
A,B=int(sys.argv[1]),int(sys.argv[2])
for a in range(2,A+1):
    for b in range(2,B+1):
        l=tuple([b]*a)
        n=sum(P(tuple(min(x,c) if k>=r else x for k,x in enumerate(l))) for r in range(a) for c in range(b) if (r,c)!=(0,0))
        if n!=1: print("%dx%d: %d winning moves"%(a,b,n))
print("brute done")
