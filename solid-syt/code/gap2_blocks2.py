import sympy as sp, json
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)

@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

x,t=sp.symbols('x t')
N=23

def tser(e):  # truncate to t-order < N, expanded
    return sp.expand(sp.series(sp.expand(e), t, 0, N).removeO())

z=sp.symbols('z')
Cat=sp.Integer(1)
for _ in range(N):
    Cat=sp.expand(sp.series(1+z*Cat**2, z, 0, N//2+2).removeO())
ystar=tser(t/x*Cat.subs(z,t**2/x**2))

# ground truth check first: t^2*Cat(t^2/x^2... note weight is Cat(t^2*xbar) NOT Cat(t^2/x^2)!
# CAREFUL: pairing weight was Cat(t^2/x) (argument t^2*xbar), while ystar = (t/x)*Cat(t^2/x^2)?? 
# Session 7: ystar = (t/x)Cat(t^2/x) with Cat arg t^2*xbar. y* solves t*y^2 - x*y + t=0 => y*=(x-sqrt(x^2-4t^2))/(2t)
# => y* = (t/x)*sum Cat_k (t/x)^{2k} = (t/x)Cat(t^2/x^2). And the WEIGHT was Cat(t^2/x) [arg t^2 xbar].
# Both appeared; recheck which is right: weight W = C(t^2 xbar), ystar_pair = t*xbar*C(t^2 xbar).
# ystar_pair solves x*y = t(1+y^2)? y=t*xbar*C: x*y = t*C; t(1+y^2)=t(1+t^2*xbar^2*C^2); C=1+ (t^2 xbar)C^2??
# Cat identity: C(z)=1+zC(z)^2 with z=t^2*xbar: C = 1 + t^2*xbar*C^2 => t*C = t + t^3*xbar*C^2 = t(1+ (t*xbar*C)^2 * x)...
# t(1+y^2) = t + t*(t^2 xbar^2 C^2) = t + xbar*(t^3 xbar C^2) = t + xbar*t*(C-1)*... since t^2 xbar C^2 = C-1:
# = t + xbar*t*(C-1) . And x*y = t*C. Equal iff tC = t + t*xbar*(C-1) iff C-1 = xbar(C-1) iff false.
# => ystar_pair does NOT solve xy=t(1+y^2); the quadratic root is y*=(t/x)Cat(t^2/x^2). My session-7 "KEY" note conflated them!
# Determine the correct pairing directly here.
Wt = tser(Cat.subs(z, t**2/x))            # C(t^2 xbar)
ys_pair = tser(t/x*Wt)                     # t*xbar*C(t^2*xbar)
ys_quad = ystar                            # (t/x)*C(t^2/x^2) : root of ty^2-xy+t=0

Q=sp.Integer(0)
ypows={}
def build_pows(base):
    d={0:sp.Integer(1)}
    for i in range(1,N):
        d[i]=tser(d[i-1]*base)
        if d[i]==0: break
    return d
Pp=build_pows(ys_pair)
for m in range(N):
    for i in range(m+1):
        for j in range(min(m+1, max(Pp)+1)):
            v=RK(m,i,j)
            if v and j in Pp: Q+= v*x**i*t**m*Pp[j]
Q=tser(Q)
G=tser(Wt*Q)
P=json.load(open('pieces.json')); Sv=P['S']
ok=all(sp.expand(G.coeff(t,3*n)).coeff(x,0)==Sv[n] for n in range(0,(N-1)//3))
print("ground truth (session-7 pairing, ys_pair): S(n)=CT[C(t^2 xbar) Q(x, t*xbar*C(t^2 xbar))]:", ok)
