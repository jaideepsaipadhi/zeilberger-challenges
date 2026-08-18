import sympy as sp, json
from fractions import Fraction
t=sp.symbols('t')
n,N=sp.symbols('n N')
v=[sp.sympify(x) for x in json.load(open('ode_A.json'))]
q_const, qA, qA1, qA2, qA3 = v
qs=[qA,qA1,qA2,qA3]

# sanity: verify ODE on the series of A
import sympy as sp2
T=sp.Integer(0)
NT=26
for _ in range(NT+3):
    T=sp.expand(sp.series(t*(2+T**3), t, 0, NT+3).removeO())
sq=sp.sqrt(1-t*T**2).series(t,0,NT).removeO()
A=sp.expand(sp.series((T/(2*t) + (sq-1+t*T**2/2)/t**3)/sq, t,0,NT-2).removeO())
expr=q_const + qA*A + qA1*sp.diff(A,t) + qA2*sp.diff(A,t,2) + qA3*sp.diff(A,t,3)
expr=sp.expand(sp.series(expr,t,0,NT-6).removeO())
print("ODE verified on series to t^%d:"%(NT-7), expr==0)

# recurrence: [t^N] of q_const + sum_j q_j(t) A^{(j)} = 0
# [t^N](t^k A^{(j)}) = ff(N-k+j, j) * h(N-k+j), ff = falling factorial product (N-k+1)...(N-k+j)
terms={}   # shift m -> coefficient polynomial in N, where index = N+m
maxdeg=max(sp.degree(q,t) for q in qs)
for j,q in enumerate(qs):
    for k in range(sp.degree(q,t)+1):
        c=q.coeff(t,k)
        if c==0: continue
        m=j-k  # index N-k+j = N+m
        ff=sp.prod([ (N-k+1+i) for i in range(j) ]) if j>0 else sp.Integer(1)
        terms[m]=sp.expand(terms.get(m,0)+c*ff)
print("shifts present:", sorted(terms))
# support: h(N)=T(n) iff N=3n. Set N=3n; contributions only from shifts m ≡ 0 mod 3
rec={}
for m,cf in terms.items():
    if m%3==0:
        rec[m//3]=sp.expand(cf.subs(N,3*n))
    else:
        # must vanish identically when restricted to supported indices... they multiply h at non-multiples of 3 => auto-zero terms; fine
        pass
shifts=sorted(rec)  # in units of n; h(3n+3i)=T(n+i)
print("recurrence shifts (in n):", shifts, " => order", max(shifts)-min(shifts))
# recurrence: sum_i rec[i]*T(n+i) = 0 for 3n > deg(q_const) (=3) i.e. n>=2 conservatively
# verify against data
P=json.load(open('pieces.json')); Tv=P['T']
ok=True
for nn in range(2, 40):
    ssum=0
    for i in shifts:
        if 0<=nn+i<len(Tv):
            ssum+=int(rec[i].subs(n,nn))*Tv[nn+i]
        else: ssum=None; break
    if ssum not in (0,None): ok=False; print("fail at n=",nn); break
print("recurrence verified on data:", ok)

# THE IDENTITY: plug the closed-form ratio rho(m)=T(m)/T(m-1)
rho=lambda m: 6*(2*m-1)*(6*m-1)*(6*m+1)*(7*m+5)/((m+1)*(4*m+3)*(4*m+5)*(7*m-2))
lo=min(shifts)
ident=sp.Integer(0)
for i in shifts:
    prod=sp.Integer(1)
    for step in range(i-lo):
        prod=prod*rho(n+lo+1+step)
    ident+=rec[i]*prod
ident=sp.simplify(sp.together(ident))
print("closed form satisfies the recurrence identically:", ident==0)

# initial conditions + leading coeff check
lead=rec[max(shifts)]
roots=sp.solve(sp.Eq(lead,0),n)
print("leading coeff integer roots (n>=2)?:", [r for r in roots if r.is_integer and r>=2])
print("initial values match (n=0..5):", all(Fraction(Tv[m])==Fraction(int(Tv[m])) for m in range(6)))
