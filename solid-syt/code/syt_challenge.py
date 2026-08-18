# Zeilberger's 1st Rigorous Challenge ($0.01): prove the empirical 2nd-order recurrence
# for g(n) = # Solid SYT of shape [[n,n],[n,1]]
# = linear extensions = lattice walks 0 -> (n,n,n,1) in coords (x11,x12,x21,x22),
#   staying in x11>=x12>=x22, x11>=x21>=x22, with x22 <= 1
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)

def g(n):
    # DP over (x11,x12,x21,x22), sum layers
    from collections import defaultdict
    f={(0,0,0,0):1}
    for s in range(1, 3*n+2):
        nf={}
        for x11 in range(min(n,s)+1):
            for x12 in range(min(x11,s)+1):
                for x21 in range(min(x11,s-x11-x12 if False else x11)+1):
                    x22 = s - x11 - x12 - x21
                    if x22<0 or x22>1: continue
                    if x12<x22 or x21<x22: continue
                    tot=0
                    for k,(a,b,c,d) in enumerate([(x11-1,x12,x21,x22),(x11,x12-1,x21,x22),
                                                   (x11,x12,x21-1,x22),(x11,x12,x21,x22-1)]):
                        if a<0 or b<0 or c<0 or d<0: continue
                        if a<b or a<c or b<d or c<d: continue
                        v=f.get((a,b,c,d))
                        if v: tot+=v
                    if tot: nf[(x11,x12,x21,x22)]=tot
        f=nf if s<3*n+1 else f|nf
        if s<3*n+1: pass
    return f.get((n,n,n,1),0)

# cleaner: full dict accumulate
def gseq(N):
    from collections import defaultdict
    all_f={(0,0,0,0):1}
    maxs=3*N+1
    # enumerate states in order of sum
    states=defaultdict(list)
    for x11 in range(N+1):
        for x12 in range(x11+1):
            for x21 in range(x11+1):
                for x22 in range(0, min(x12,x21,1)+1):
                    states[x11+x12+x21+x22].append((x11,x12,x21,x22))
    for s in range(1,maxs+1):
        for st in states[s]:
            x11,x12,x21,x22=st
            tot=0
            for (a,b,c,d) in [(x11-1,x12,x21,x22),(x11,x12-1,x21,x22),
                              (x11,x12,x21-1,x22),(x11,x12,x21,x22-1)]:
                if a<0 or b<0 or c<0 or d<0: continue
                if a<b or a<c or b<d or c<d: continue
                v=all_f.get((a,b,c,d))
                if v: tot+=v
            if tot: all_f[st]=tot
    return [all_f.get((n,n,n,1),0) for n in range(N+1)]

seq=gseq(40)
print("g(n) for n=1..12:", seq[1:13])
with open('sytchal.txt','w') as fo:
    for n,v in enumerate(seq): fo.write(f"{n} {v}\n")
