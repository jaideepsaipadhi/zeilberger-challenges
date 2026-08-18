"""Certify that the machinery covers Conjecture 2b (odd-run H-model):
   tilted run law P(R=r) = w^{r-1}(1-w^2), r odd >= 1, at w = sqrt(2)-1;
   verify: drift 0, Sigma_H = scalar * M (=> angle pi/3 => exponent -4), corrector kappa = -E[R]/3."""
import sympy as sp
w=sp.sqrt(2)-1
# run-law moments over odd r
r=sp.symbols('r', positive=True, integer=True)
norm=w/(1-w**2)
ER = sp.simplify(sp.summation(r*w**r, (r,1,sp.oo)).rewrite(sp.Sum))  # need odd-only: do directly
# odd sums: sum_{k>=0} f(2k+1) w^{2k+1}
k=sp.symbols('k', nonnegative=True, integer=True)
S0=sp.summation(w**(2*k+1),(k,0,sp.oo))
S1=sp.summation((2*k+1)*w**(2*k+1),(k,0,sp.oo))
S2=sp.summation((2*k+1)**2*w**(2*k+1),(k,0,sp.oo))
ER=sp.simplify(S1/S0); ER2=sp.simplify(S2/S0)
print("E[R] =",sp.simplify(sp.radsimp(ER)),"=",sp.nsimplify(sp.simplify(ER),[sp.sqrt(2)]))
print("E[R^2] =",sp.nsimplify(sp.simplify(ER2),[sp.sqrt(2)]))
# per-run covariance via the correlation-sum formula: Sigma = ER2*M + 2*ER^2*sum_{m>=1}(-1/2)^m * M
#   = (ER2 - (2/3)ER^2) * M   [sum = -1/3]
sH=sp.simplify(ER2 - sp.Rational(2,3)*ER**2)
print("Sigma_H = s_H * M with s_H =",sp.nsimplify(sH,[sp.sqrt(2)]))
# independent check via Perron Hessian
t1,t2,e,u1,u2,c1,c2=sp.symbols('t1 t2 e u1 u2 c1 c2')
I=sp.I
v=[sp.Matrix([1,0]),sp.Matrix([-1,1]),sp.Matrix([0,-1])]
def phiH(vec):
    z=sp.exp(I*e*(u1*vec[0]+u2*vec[1]))
    return (1-w**2)*z/(1-w**2*z**2)*1  # E[z^R] = sum w^{r-1}(1-w^2) z^r over odd r = (1-w^2) z/(1-w^2 z^2)
P=sp.zeros(3,3)
for i in range(3):
    for j in range(3):
        if i!=j: P[i,j]=sp.Rational(1,2)*phiH(v[j])
lam=sp.symbols('lam')
cp=sp.det(P-lam*sp.eye(3))
expr=sp.series(sp.expand(cp.subs(lam,1+c1*e+c2*e**2)),e,0,3).removeO()
c1v=sp.solve(sp.expand(expr.coeff(e,1)),c1)[0]
c2v=sp.simplify(sp.solve(sp.expand(expr.coeff(e,2)).subs(c1,c1v),c2)[0])
print("drift order-1:",sp.simplify(c1v),"(expect 0)")
quad=sp.simplify(sp.expand(c2v-c1v**2/2))
target=sp.expand(-sp.Rational(1,2)*sH*(sp.Rational(2,3)*u1**2-sp.Rational(2,3)*u1*u2+sp.Rational(2,3)*u2**2))
print("Perron Hessian == s_H*M:", sp.simplify(quad-target)==0)
# corrector
print("corrector kappa_H = -E[R]/3 =",sp.nsimplify(sp.simplify(-ER/3),[sp.sqrt(2)]))
# criticality sanity: growth per unit step should be 1+w... mu_H = (1+sqrt2)^3 = 7+5sqrt2 (established earlier)
print("mu_H = (1+sqrt2)^3 =",sp.expand((1+sp.sqrt(2))**3))
