from fractions import Fraction
import json
P=json.load(open('pieces.json'))

def fit(order, D, seq, n0=1, spare=2):
    rows=[]; N=len(seq)-1
    for n in range(n0+order, N+1):
        row=[]
        for i in range(order+1):
            for j in range(D+1):
                row.append(Fraction(n)**j * seq[n-i])
        rows.append(row)
    ncols=(order+1)*(D+1)
    if len(rows)<ncols+spare: return None
    A=[r[:] for r in rows]
    piv=[]; r=0
    for c in range(ncols):
        pr=next((rr for rr in range(r,len(A)) if A[rr][c]!=0), None)
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        pv=A[r][c]; A[r]=[x/pv for x in A[r]]
        for rr in range(len(A)):
            if rr!=r and A[rr][c]!=0:
                f=A[rr][c]; A[rr]=[x-f*y for x,y in zip(A[rr],A[r])]
        piv.append(c); r+=1
        if r==len(A): break
    free=[c for c in range(ncols) if c not in piv]
    if not free: return None
    sol=[Fraction(0)]*ncols; sol[free[0]]=Fraction(1)
    for i,c in enumerate(piv):
        sol[c]=-sum(A[i][j]*sol[j] for j in free)
    from math import gcd, lcm
    L=1
    for x in sol: L=lcm(L,x.denominator)
    ints=[int(x*L) for x in sol]
    G=0
    for v in ints: G=gcd(G,v)
    ints=[v//G for v in ints] if G else ints
    for n in range(order+1, len(seq)):
        tot=0; idx=0
        for i in range(order+1):
            for j in range(D+1):
                tot+=ints[idx]*(n**j)*seq[n-i]; idx+=1
        if tot!=0: return None
    return ints

results={}
for name in ['K','T','S']:
    seq=P[name]
    done=False
    for order in range(1,6):
        for D in range(1,15):
            r=fit(order,D,seq)
            if r:
                print(f"{name}: minimal found -> order {order}, degree {D}")
                results[name]=dict(order=order,degree=D,coeffs=r)
                done=True; break
        if done: break
    if not done: print(f"{name}: nothing up to order 5 / degree 14")
json.dump(results, open('piece_recs.json','w'))
