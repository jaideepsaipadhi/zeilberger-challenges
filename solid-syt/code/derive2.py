import json, sympy as sp
from fractions import Fraction
P=json.load(open('pieces.json'))
K,T,S=P['K'],P['T'],P['S']
A=[(3*m+1)*K[m] for m in range(len(K))]

# rho(n) = S(n)/A(n): should be rational in n. Fit rho = p(n)/q(n) by exact linear algebra.
n=sp.symbols('n')
vals=[Fraction(S[m],A[m]) for m in range(len(A))]
# ansatz degrees up to 3/3
for dp in range(0,4):
    for dq in range(0,4):
        # p(n) - rho*q(n) = 0 -> linear in coeffs
        rows=[]; 
        for m in range(0,len(vals)):
            row=[Fraction(m)**j for j in range(dp+1)]+[-vals[m]*Fraction(m)**j for j in range(dq+1)]
            rows.append(row)
        import itertools
        ncols=dp+dq+2
        Amat=sp.Matrix([[sp.Rational(x) for x in r] for r in rows[:ncols+2]])
        ns=Amat.nullspace()
        if ns:
            v=ns[0]
            pnum=sum(v[j]*n**j for j in range(dp+1))
            pden=sum(v[dp+1+j]*n**j for j in range(dq+1))
            rho=sp.cancel(pnum/pden)
            # verify on all values
            if all(sp.Rational(S[m],A[m])==rho.subs(n,m) for m in range(len(A))):
                print("rho(n) = S(n)/((3n+1)K(n)) =", sp.factor(rho))
                sp.pickle=rho
                import pickle
                pickle.dump(sp.srepr(rho), open('rho.pkl','wb'))
                raise SystemExit
print("no rational rho up to degree 3/3")
