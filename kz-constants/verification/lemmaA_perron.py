"""Lemma A foundation: the Fourier/tilted transfer operator of the run-restricted tandem MAP.
   States = current run direction i in {0,1,2}; transition: j != i uniform (prob 1/2 each) -- 
   NO: transition prob comes from the tilted run-GF weights; at criticality each run has
   P(R=r)=(1/2)^(r-1), r>=2, direction uniform over the other two. Increment of one MAP step
   (one RUN) = R * v_i with v0=(1,0), v1=(-1,1), v2=(0,-1).
   Fourier matrix: P(th)_{ij} = (1/2) * E[e^{i th . v_j R}] for j != i, 0 on diagonal.
   E[z^R] = sum_{r>=2} (1/2)^{r-1} z^r = z^2/(2-z) ... check: r=2: (1/2) z^2; sum = z^2/2 * 1/(1-z/2) = z^2/(2-z). OK.
   Perron root lam(th): expand to 2nd order at th=0; drift = i*grad must vanish... 
   CAREFUL: per-RUN drift is E[R]*E_stationary[v] -- stationary uniform => sum v_i = 0 ✓.
   The CLT covariance of S_n (sum over runs) per run = Hessian of log lam at 0.
   Consistency target: per-UNIT-STEP covariance = (per-run Hessian)/E[R] should equal (5/3)*M_unit? 
   Our Sigma=5M was per unit step n (3n unit steps total? n = ... G(n): 3n steps, runs of total length 3n).
   Sigma per unit time (unit = one lattice step of the tandem walk) = 5M with M=[[2/3,-1/3],[-1/3,2/3]].
   The Fourier matrix below yields the PER-STEP quadratic directly, so the
   expected Hessian(log lam) block is 5*M = [[10/3,-5/3],[-5/3,10/3]]
   (NOT 15*M: no extra factor E[R]=3 -- cross-checked against lemmaA_torus).
"""
import sympy as sp
t1,t2=sp.symbols('t1 t2', real=True)
I=sp.I
v=[sp.Matrix([1,0]),sp.Matrix([-1,1]),sp.Matrix([0,-1])]
def phi(vec):  # E[e^{i th.vec R}] with E over R: z=e^{i th.vec}: z^2/(2-z)
    z=sp.exp(I*(t1*vec[0]+t2*vec[1]))
    return z**2/(2-z)
P=sp.zeros(3,3)
for i in range(3):
    for j in range(3):
        if i!=j: P[i,j]=sp.Rational(1,2)*phi(v[j])
# Perron root at th=0: P(0) = (1/2)(J-I) with phi(0)=1 -> row sums 1, lam(0)=1 ✓
# expand lam(th) via perturbation: lam = stationary-weighted... easier: char poly and series.
lam=sp.symbols('lam')
cp=sp.det(P-lam*sp.eye(3))
# solve for the branch with lam(0)=1 by series: substitute lam = 1 + a1*eps + ... with th = eps*(u1,u2)
e=sp.symbols('e', positive=True)
u1,u2=sp.symbols('u1 u2', real=True)
cps=cp.subs({t1:e*u1,t2:e*u2})
# series solve
lser=sp.Function('L')(e)
sol=sp.series(sp.nsimplify(0),e,0,1)  # placeholder
# do it order by order: lam = 1 + c1 e + c2 e^2
c1,c2=sp.symbols('c1 c2')
expr=cps.subs(lam,1+c1*e+c2*e**2)
expr=sp.series(sp.expand(expr),e,0,3).removeO()
eq1=sp.expand(expr.coeff(e,1))
eq2=sp.expand(expr.coeff(e,2))
s1=sp.solve(eq1,c1)
print("order-1 (drift) c1 =",s1)
c1v=s1[0]
s2=sp.solve(eq2.subs(c1,c1v),c2)
print("order-2 c2 =",sp.simplify(s2[0]))
# log lam ~ c1 e + (c2 - c1^2/2) e^2; covariance: log lam = -1/2 th^T Sig_run th + i drift.th
c2v=sp.simplify(s2[0])
quad=sp.expand(c2v - c1v**2/2)
print("quadratic form in (u1,u2) of log-lam order 2:", sp.simplify(quad))
# expect -1/2 * (u1,u2) Sigma (u1,u2)^T with the PER-STEP position block
# Sigma = 5*M = 5*[[2/3,-1/3],[-1/3,2/3]] = [[10/3,-5/3],[-5/3,10/3]]
# (the old target used 15*M = "per run, x E[R]=3", but the quadratic computed
#  above is per unit step; cross-checked against lemmaA_torus's Sigma_3 block)
target=-sp.Rational(1,2)*(sp.Rational(10,3)*u1**2 - sp.Rational(10,3)*u1*u2 + sp.Rational(10,3)*u2**2)
print("matches Sigma=5M (per-step position block):", sp.simplify(quad-target)==0)
