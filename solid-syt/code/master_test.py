from fractions import Fraction
from functools import lru_cache
import sys, sympy as sp
sys.setrecursionlimit(100000)
@lru_cache(maxsize=None)
def W(e,nn,d):
    if min(e,nn,d)<0 or e-d<0 or nn-d<0: return 0
    if e==0 and nn==0 and d==0: return 1
    tot=0
    if e>0: tot+=W(e-1,nn,d)
    if nn>0: tot+=W(e,nn-1,d)
    if d>0: tot+=W(e,nn,d-1)
    return tot

NMAX=13
# ratio in c of  t(n,a,c) = W(n, n-c, n-a),  0<=c<=a<=n
pairs=[]
for n in range(2,NMAX+1):
    for a in range(0,n+1):
        for c in range(0,a):
            v0=W(n,n-c,n-a); v1=W(n,n-c-1,n-a)
            if v0: pairs.append(((n,a,c), Fraction(v1,v0)))

def fit_ratio3(pairs, dmax=4):
    monos=[(i,j,k) for i in range(dmax+1) for j in range(dmax+1) for k in range(dmax+1) if i+j+k<=dmax]
    rows=[]
    for (n,a,c),r in pairs:
        row=[Fraction(n)**i*Fraction(a)**j*Fraction(c)**k for (i,j,k) in monos] \
           +[-r*Fraction(n)**i*Fraction(a)**j*Fraction(c)**k for (i,j,k) in monos]
        rows.append(row)
    ncols=2*len(monos)
    if len(rows)<ncols+3: return None
    import random
    random.seed(1); random.shuffle(rows)
    A=[row[:] for row in rows[:ncols+6]]
    piv=[]; rr=0
    for col in range(ncols):
        pr=next((k for k in range(rr,len(A)) if A[k][col]!=0), None)
        if pr is None: continue
        A[rr],A[pr]=A[pr],A[rr]
        pv=A[rr][col]; A[rr]=[x/pv for x in A[rr]]
        for k in range(len(A)):
            if k!=rr and A[k][col]!=0:
                f=A[k][col]; A[k]=[x-f*y for x,y in zip(A[k],A[rr])]
        piv.append(col); rr+=1
        if rr==len(A): break
    free=[c for c in range(ncols) if c not in piv]
    if not free: return None
    sol=[Fraction(0)]*ncols; sol[free[0]]=Fraction(1)
    for i,col in enumerate(piv):
        sol[col]=-sum(A[i][j]*sol[j] for j in free)
    # verify on ALL pairs
    n_,a_,c_=sp.symbols('n a c')
    m2=len(monos)
    P=sum(sp.Rational(sol[i])*n_**monos[i][0]*a_**monos[i][1]*c_**monos[i][2] for i in range(m2))
    Q=sum(sp.Rational(sol[m2+i])*n_**monos[i][0]*a_**monos[i][1]*c_**monos[i][2] for i in range(m2))
    ratio=sp.cancel(P/Q)
    for (n,a,c),r in pairs:
        val=ratio.subs({n_:n,a_:a,c_:c})
        if val!=r: return ('FAIL_VERIFY', ratio, (n,a,c), r, val)
    return ('OK', sp.factor(ratio))

res=fit_ratio3(pairs, dmax=4)
print("W(n,n-c-1,n-a)/W(n,n-c,n-a):", res[0])
if res[0]=='OK': print("  ratio =", res[1])
elif res[0]=='FAIL_VERIFY': print("  candidate", res[1], "fails at", res[2], ":", res[3], "vs", res[4])
