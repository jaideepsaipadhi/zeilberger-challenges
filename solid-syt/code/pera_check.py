"""Prove the per-a diagonal closed form: k(3n-a;a,a) from [x^a]Qd + Lagrange, as a Gamma identity."""
import sympy as sp
from functools import lru_cache
a,n=sp.symbols('a n', positive=True, integer=True)
# [x^a] Qd = (1/t)[ binom(2a,a) T^{2a+1}/(2*4^a)  -  binom(2a+2,a+1)*(a+1)/(2(a+2)*4^{a+1}) * T^{2a+4} ]
# derivation: Qd = T/(2t)*(1-u)^{-1/2} + (1/(t x^2)) [1 - (1-u/2)(1-u)^{-1/2}], u = x T^2
# [u^j](1-u)^{-1/2} = binom(2j,j)/4^j ; [u^j] of bracket = -binom(2j-2,j-1)/4^{j-1} * (j-1)/(2j)  (j>=1), 0 at j=0
j=sp.symbols('j', positive=True, integer=True)
lhs_j=-(sp.binomial(2*j,j)/4**j - sp.Rational(1,2)*sp.binomial(2*j-2,j-1)/4**(j-1))
target_j=-sp.binomial(2*j-2,j-1)/4**(j-1)*(j-1)/(2*j)
print("bracket coefficient identity:", sp.simplify(sp.combsimp(lhs_j-target_j))==0)
# Lagrange: L(N,k) = [t^N] T^k = (k/N) binom(N,(N-k)/3) 2^{N-(N-k)/3} when N=k mod 3
def L(N,k): return sp.together(k/sp.sympify(N))*sp.binomial(N,(N-k)/3)*2**(N-(N-k)/3)
# k(3n-a;a,a) = [t^{3n-a}] [x^a]Qd = [t^{3n-a+1}] of the T-part:
expr = sp.binomial(2*a,a)/(2*4**a)*L(3*n-a+1,2*a+1) - sp.binomial(2*a+2,a+1)*(a+1)/(2*(a+2)*4**(a+1))*L(3*n-a+1,2*a+4)
# target: K(n) * (a+1)! (2a+1)! n! (3n-a)! / (4^a (a!)^3 (n-a)! (3n)!) with K(n)=4^n(3n)!/((n+1)!(2n+1)!)
Kn=4**n*sp.factorial(3*n)/(sp.factorial(n+1)*sp.factorial(2*n+1))
target=Kn*sp.factorial(a+1)*sp.factorial(2*a+1)*sp.factorial(n)*sp.factorial(3*n-a)/(4**a*sp.factorial(a)**3*sp.factorial(n-a)*sp.factorial(3*n))
ratio=sp.simplify(sp.combsimp(sp.gammasimp(expr/target)))
print("per-a identity ratio (want 1):", ratio)
# numeric spot check independent of simplification
from functools import lru_cache
@lru_cache(maxsize=None)
def RK(m,i,jj):
    if i<0 or jj<0 or m<0: return 0
    if m==0: return 1 if (i==0 and jj==0) else 0
    return RK(m-1,i-1,jj)+RK(m-1,i,jj-1)+RK(m-1,i+1,jj+1)
ok=True
for nn in range(1,9):
    for aa in range(0,nn+1):
        v1=int(expr.subs({n:nn,a:aa}))
        v2=RK(3*nn-aa,aa,aa)
        if v1!=v2: ok=False; print("mismatch",nn,aa,v1,v2)
print("numeric check n<=8, all a:", ok)
