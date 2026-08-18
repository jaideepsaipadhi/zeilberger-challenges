import json, sympy as sp
n=sp.symbols('n')
rA=lambda m: 6*(3*m-1)*(3*m+1)/((m+1)*(2*m+1))       # A(m)/A(m-1), A=(3n+1)K
rT=lambda m: 6*(2*m-1)*(6*m-1)*(6*m+1)*(7*m+5)/((m+1)*(4*m+3)*(4*m+5)*(7*m-2))
c =lambda m: (8*m**2+m-24)/(4*(m+2)*(2*m+3))

# coordinates of g(n+i) in basis (A(n), T(n)):
w=[sp.Matrix([c(n), 1]),
   sp.Matrix([c(n+1)*rA(n+1), rT(n+1)]),
   sp.Matrix([c(n+2)*rA(n+1)*rA(n+2), rT(n+1)*rT(n+2)])]
det=lambda a,b: sp.together(a[0]*b[1]-a[1]*b[0])
q=[det(w[1],w[2]), -det(w[0],w[2]), det(w[0],w[1])]   # q0 g(n)+q1 g(n+1)+q2 g(n+2)=0
den=sp.lcm([sp.denom(x) for x in q])
qp=[sp.expand(sp.cancel(x*den)) for x in q]
g0=sp.gcd(sp.gcd(qp[0],qp[1]),qp[2])
qp=[sp.expand(sp.cancel(x/g0)) for x in qp]
print("derived degrees:", [sp.degree(x,n) for x in qp])

polys=json.load(open('recurrence.json'))
m=sp.symbols('m')
p=[sum(cf*m**j for j,cf in enumerate(pol)) for pol in polys]
guess=[sp.expand(p[2].subs(m,n+2)), sp.expand(p[1].subs(m,n+2)), sp.expand(p[0].subs(m,n+2))]
ok=all(sp.simplify(guess[i]*qp[j]-guess[j]*qp[i])==0 for i in range(3) for j in range(i+1,3))
print("derived operator == guessed degree-12 operator (up to scalar):", ok)

seq={}
for line in open('sytchal55.txt'):
    a,v=line.split(); seq[int(a)]=int(v)
fq=[sp.lambdify(n,x) for x in qp]
print("numeric check on n in [1,53]:",
      all(sum(int(round(fq[i](N)))*seq[N+i] for i in range(3))==0 for N in range(1,54)))
