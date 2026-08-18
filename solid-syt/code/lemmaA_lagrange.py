import sympy as sp, json
from fractions import Fraction
from math import comb, factorial
n=sp.symbols('n', positive=True, integer=True)

# L(N,k) = [t^N] T^k = (k/N) * C(N, (N-k)/3) * 2^(N-(N-k)/3)   [Lagrange inversion, T = t(2+T^3)]
def Lsym(N,k):
    j=(N-k)/3
    return (sp.Integer(k)/N)*sp.binomial(N, j)*2**(N-j)

# claim: [t^{3n}] S-GF = (1/256)( -18*L(3n+4,1) + 108*L(3n+2,2) - L(3n+5,2) )
lhs = sp.Rational(1,256)*(-18*Lsym(3*n+4,1) + 108*Lsym(3*n+2,2) - Lsym(3*n+5,2))
# target: rho(n)*(3n+1)*K(n), K(n)=4^n (3n)!/((n+1)!(2n+1)!)
K = 4**n*sp.factorial(3*n)/(sp.factorial(n+1)*sp.factorial(2*n+1))
rhs = sp.Rational(3,8)*(9*n+16)/((n+2)*(2*n+3))*(3*n+1)*K

diff=sp.simplify(sp.together(sp.expand_func(lhs/rhs)))
print("ratio lhs/rhs simplified:", diff)
print("identity holds symbolically:", sp.simplify(diff-1)==0)

# numeric double-check
def Lnum(N,k):
    j=(N-k)//3
    assert (N-k)%3==0
    return Fraction(k,N)*comb(N,j)*2**(N-j)
P=json.load(open('pieces.json')); Sv=P['S']
ok=True
for m in range(0,40):
    v=Fraction(1,256)*(-18*Lnum(3*m+4,1)+108*Lnum(3*m+2,2)-Lnum(3*m+5,2))
    if v!=Sv[m]: ok=False; print("numeric fail at",m); break
print("numeric check n<40:", ok)
