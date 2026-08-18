from math import comb
from collections import defaultdict
import json

N=46
# Q[m][(u,v)] = # quadrant walks steps {(1,1),(-1,0),(0,-1)} from (u,v) to (0,0), length m
# recurrence on first step
Mmax=3*N
Umax=2*N+2
prev={(0,0):1}
# we need Q for all m up to Mmax; store per m the dict
Qs=[prev]
for m in range(1,Mmax+1):
    cur={}
    for u in range(0,Umax+1):
        row_prev=prev
        for v in range(0,Umax+1):
            tot=row_prev.get((u+1,v+1),0)
            if u>=1: tot+=row_prev.get((u-1,v),0)
            if v>=1: tot+=row_prev.get((u,v-1),0)
            if tot: cur[(u,v)]=tot
    Qs.append(cur); prev=cur

def ballot(a,c):
    if c>a or c<0: return 0
    return comb(a+c,c)*(a-c+1)//(a+1)

K=[Qs[3*n].get((0,0),0) for n in range(N+1)]
T=[sum(Qs[3*n-a].get((a,a),0) for a in range(n+1)) for n in range(N+1)]
S=[sum(ballot(a,c)*Qs[3*n-a-c].get((a,a-c),0) for a in range(n+1) for c in range(a+1)) for n in range(N+1)]

seq={}
for line in open('sytchal55.txt'):
    n,v=line.split(); seq[int(n)]=int(v)
ok=all((3*n+1)*K[n]-2*S[n]+T[n]==seq[n] for n in range(1,N+1))
print("identity holds for all n <= %d:"%N, ok)
json.dump({'K':K,'T':T,'S':S}, open('pieces.json','w'))
print("T(1..6):", T[1:7])
print("S(1..6):", S[1:7])
