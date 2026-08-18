import json, sympy as sp
P=json.load(open('pieces.json'))
R=json.load(open('piece_recs.json'))
n=sp.symbols('n')

for name in ['K','T','S']:
    d=R[name]; o,D,co=d['order'],d['degree'],d['coeffs']
    p0=sum(co[j]*n**j for j in range(D+1))
    p1=sum(co[D+1+j]*n**j for j in range(D+1))
    ratio=sp.factor(-p1/p0)   # seq(n)/seq(n-1)
    print(f"{name}(n)/{name}(n-1) =", ratio)

# closed-form hunt: verify against data via Gamma products
from math import comb
K,T,S=P['K'],P['T'],P['S']
import fractions
def check_formula(name, data, f, rng=range(1,30)):
    ok=all(f(m)==data[m] for m in rng)
    print(f"  {name} closed form verified on n in [1,30):", ok)

# K(n) = 4^n (3n)! / ((n+1)! (2n+1)!)  (Kreweras)
import math
check_formula("K", K, lambda m: 4**m*math.factorial(3*m)//(math.factorial(m+1)*math.factorial(2*m+1)))
