#!/usr/bin/env python3
"""Extra validations: (A) our exact values H(1..6) vs Rubine's printed values (object identity,
arXiv 2512.21785); (B) production code vs exact for D=6..10 at small M; (C) fresh-prime check
of any reconstructed G_D*_M1000.txt files present."""
import sys, os, glob
sys.set_int_max_str_digits(0)
from math import factorial, comb
from fractions import Fraction
from collections import defaultdict
from itertools import product
from sympy import nextprime, prevprime
from geode_batch import G_modp_batch

def C_gen(m,ks):
    top=sum(k*x for k,x in zip(ks,m)); bot=1+sum((k-1)*x for k,x in zip(ks,m))
    v=factorial(top)//factorial(bot)
    for x in m: v//=factorial(x)
    return v
def G_exact(M,D):
    ks=list(range(2,D+2)); tot=0
    for al in product(range(M+1),repeat=D-1):
        i=sum(al); c=factorial(i)
        for a in al: c//=factorial(a)
        tot+=(-1)**i*c*C_gen([M+1+i]+[M-a for a in al],ks)
    return tot

ok=True
print("[A] H(1..6) vs Rubine's printed values (arXiv 2512.21785):")
rub=[12344,2408941884,894971463204720,446324644841317281200,
     263656050352833337510832640,173882340006327290808417397911384]
for n in range(1,7):
    m = G_exact(n,4)==rub[n-1]; ok&=m
    print("   H(%d): %s"%(n,"MATCH" if m else "FAIL"))
print("[B] production vs exact, D=6..10, M=1,2:")
p=int(nextprime(2**30))
for D in [6,7,8,9,10]:
    for M in [1,2]:
        m = G_exact(M,D)%p==G_modp_batch(M,D,[p])[0]; ok&=m
        print("   D=%d M=%d: %s"%(D,M,"MATCH" if m else "FAIL"))
print("[C] fresh-prime end-to-end checks of reconstructed values present in this directory:")
ps=[int(prevprime(2**30)), int(prevprime(2**30-10**6))]
for f in sorted(glob.glob("G_D*_M1000.txt")):
    D=int(f.split('_')[1][1:])
    G=int(open(f).read())
    fresh=G_modp_batch(1000,D,ps)
    m=all(G%q==int(r) for q,r in zip(ps,fresh)); ok&=m
    print("   %s (%d digits): %s"%(f,len(str(G)),"MATCH" if m else "FAIL"))
print("VALIDATE_EXTRA:", "ALL PASS" if ok else "FAILURE")
