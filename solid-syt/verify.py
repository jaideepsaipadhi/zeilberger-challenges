#!/usr/bin/env python3
"""Self-contained verification: the order-2 degree-12 operator in data/recurrence.json
annihilates the 55 exact terms in data/sytchal55.txt. Tries all shift/order conventions
and reports which holds; PASS requires one convention to annihilate ALL testable n."""
import json, re
from fractions import Fraction

terms={}
for line in open('data/sytchal55.txt'):
    m=re.findall(r'-?\d+', line)
    if len(m)>=2: terms[int(m[0])]=int(m[1])
assert len(terms)>=50, "expected >=50 terms, got %d"%len(terms)
N=max(terms)
rec=json.load(open('data/recurrence.json'))
ps=[rec[k] for k in sorted(rec.keys())] if isinstance(rec,dict) else rec
def ev(coeffs,n):
    return sum(Fraction(c)*n**k for k,c in enumerate(coeffs))
r=len(ps)-1
def test(shift_dir, rev):
    P=ps[::-1] if rev else ps
    bad=0; tested=0
    for n in range(1+r, N+1):
        idx=[n-shift_dir*i for i in range(r+1)]
        if not all(1<=j<=N for j in idx): continue
        tot=sum(ev(P[i],n)*terms[idx[i]] for i in range(r+1))
        tested+=1
        if tot!=0: bad+=1
    return tested,bad
results={}
for sd in (1,-1):
    for rev in (False,True):
        t,b=test(sd,rev)
        results[(sd,rev)]=(t,b)
        print("convention shift=%+d rev=%s: %d tested, %d nonzero"%(sd,rev,t,b))
winners=[k for k,(t,b) in results.items() if t>=45 and b==0]
print("VERIFY:", "PASS (convention %s)"%str(winners[0]) if winners else "FAIL")
