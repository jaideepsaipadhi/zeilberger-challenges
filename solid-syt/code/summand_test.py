# Summands:
#  Lemma B: t_B(n,a) = W(n, n, n-a)   = # quadrant Kreweras walks with n E, n N, (n-a) D steps
#  Lemma A: t_A(n,a,c) = Ballot(a,c) * W(n, n-c, n-a)
# W(alpha,beta,gamma) = k(alpha+beta+gamma; alpha-gamma, beta-gamma), steps E=(1,0),N=(0,1),D=(-1,-1), stay in quadrant.
# TEST: is t_B proper hypergeometric, i.e. are t_B(n,a+1)/t_B(n,a) and t_B(n+1,a)/t_B(n,a) rational in (n,a)?
from fractions import Fraction
from collections import defaultdict

# forward DP: count walks from (0,0) by (#E,#N,#D) -> endpoint determined; track position & step counts
# W(alpha,beta,gamma): DP over (e,nn,d) with position (e-d, nn-d) >= 0 throughout.
# order matters; count sequences: f[e][nn][d] = # valid orderings
NMAX=14
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)
@lru_cache(maxsize=None)
def W(e,nn,d):
    # walks using e E-steps, nn N-steps, d D-steps in some order, staying in quadrant
    if min(e,nn,d)<0: return 0
    if e-d<0 or nn-d<0: return 0   # endpoint must be in quadrant; also intermediate handled recursively
    if e==0 and nn==0 and d==0: return 1
    # remove last step
    tot=0
    # last step E: before it, position (e-1-d, nn-d) must be >=0 (it is if current valid)
    if e>0: tot+=W(e-1,nn,d)
    if nn>0: tot+=W(e,nn-1,d)
    if d>0:
        # last step D: position before = (e-d+1, nn-d+1); need current pos >=0 which holds; prior walk must be valid with counts (e,nn,d-1): its endpoint (e-(d-1), nn-(d-1)) >= (1,1) required? endpoint before D-step is (e-d+1, nn-d+1) >= (1,1) automatically since e>=d... e-d>=0 => e-d+1>=1 ok
        tot+=W(e,nn,d-1)
    return tot

# sanity: W(n,n,n) = Kreweras excursions A006335
print("W(n,n,n):", [W(n,n,n) for n in range(1,7)], " (expect 2,16,192,2816,46592,835584)")

# t_B table
tB={}
for n in range(1,NMAX+1):
    for a in range(0,n+1):
        tB[(n,a)]=W(n,n,n-a)

# ratio in a: fit rational function r(n,a) = t(n,a+1)/t(n,a) of small degree via exact linear algebra
import itertools
def fit_ratio(pairs, dmax=4):
    # pairs: list of ((n,a), value_num, value_den) with ratio = num/den as Fraction
    # ansatz: P(n,a) - r*Q(n,a) = 0, deg <= dmax in each
    monos=[(i,j) for i in range(dmax+1) for j in range(dmax+1) if i+j<=dmax]
    rows=[]
    for (n,a),r in pairs:
        row=[Fraction(n)**i*Fraction(a)**j for (i,j) in monos]+[-r*Fraction(n)**i*Fraction(a)**j for (i,j) in monos]
        rows.append(row)
    ncols=2*len(monos)
    if len(rows)<ncols+2: return None
    A=[row[:] for row in rows]
    piv=[]; rr=0
    for c in range(ncols):
        pr=next((k for k in range(rr,len(A)) if A[k][c]!=0), None)
        if pr is None: continue
        A[rr],A[pr]=A[pr],A[rr]
        pv=A[rr][c]; A[rr]=[x/pv for x in A[rr]]
        for k in range(len(A)):
            if k!=rr and A[k][c]!=0:
                f=A[k][c]; A[k]=[x-f*y for x,y in zip(A[k],A[rr])]
        piv.append(c); rr+=1
        if rr==len(A): break
    free=[c for c in range(ncols) if c not in piv]
    if not free: return None
    sol=[Fraction(0)]*ncols; sol[free[0]]=Fraction(1)
    for i,c in enumerate(piv):
        sol[c]=-sum(A[i][j]*sol[j] for j in free)
    return monos, sol

pairs=[]
for n in range(2,NMAX+1):
    for a in range(0,n):
        if tB[(n,a)]:
            pairs.append(((n,a), Fraction(tB[(n,a+1)],tB[(n,a)])))
res=fit_ratio(pairs, dmax=4)
print("\nt_B(n,a+1)/t_B(n,a) rational of joint degree <=4:", "YES" if res else "NO")
if res:
    monos,sol=res
    import sympy as sp
    nn_,aa_=sp.symbols('n a')
    m2=len(monos)
    P=sum(sp.Rational(sol[i])*nn_**monos[i][0]*aa_**monos[i][1] for i in range(m2))
    Q=sum(sp.Rational(sol[m2+i])*nn_**monos[i][0]*aa_**monos[i][1] for i in range(m2))
    print("  ratio =", sp.factor(sp.cancel(P/Q)))
