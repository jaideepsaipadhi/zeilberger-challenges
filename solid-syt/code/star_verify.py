import sympy as sp
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

x,t=sp.symbols('x t')
N=13  # t-order

# series from DP
R  = sp.expand(t*sum(RK(m,i,0)*x**i*t**m for m in range(N) for i in range(m+1)))
c  = sp.expand(t*sum(RK(m,0,0)*t**m for m in range(N)))
Qdbar = sp.expand(sum(RK(m,i,i)*x**(-i)*t**m for m in range(N) for i in range(m//2+1)))

# delta and Y0
delta = (1-t*x)**2 - 4*t**2/x
sq = sp.sqrt(delta).series(t,0,N).removeO()
Y0 = sp.expand(sp.series(((1-t*x) - sq)/(2*t), t, 0, N-1).removeO())

LHS = sp.expand(sp.series(-Qdbar/x*sq, t, 0, N-1).removeO())
RHS = sp.expand(sp.series(c - 1/x - 2*R + 2*x*Y0, t, 0, N-1).removeO())
diff = sp.expand(LHS-RHS)
# truncate consistently: compare coefficients of t^m for m < N-1
ok=True
for m in range(N-1):
    d=sp.expand(diff.coeff(t,m))
    # ignore x-degrees beyond DP truncation range? DP series are complete per t-order, so d must vanish
    if sp.simplify(d)!=0:
        ok=False; print("mismatch at t^%d:"%m, sp.simplify(d)); break
print("Diagonal equation (STAR) verified to t^%d:"%(N-2), ok)
