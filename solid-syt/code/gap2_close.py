import sympy as sp, json
from functools import lru_cache
import sys
sys.setrecursionlimit(100000)
@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

# (i) the new combinatorial identity: S(n) = RK(3n+1; 1,0)
P=json.load(open('pieces.json')); Sv=P['S']
ok1=all(RK(3*n+1,1,0)==Sv[n] for n in range(0,13))
print("(i) S(n) = #RK-walks length 3n+1 to (1,0), n<=12:", ok1)

# (ii) [x^1] of the proven axis formula equals t^2 * S-GF-closed-form, as a rational identity in T
T=sp.symbols('T')
t=T/(2+T**3)          # rational uniformization t(T)
beta=1+T**3/4
# sqrt(Dminus(x)) = 1 + d1 x + d2 x^2 + ... : d1=-T*beta/2, d2=(T^2/8)(1-beta^2)
d1=-T*beta/2
d2=sp.Rational(1,8)*T**2*(1-beta**2)
# [x^1]R = (1/2)[ 1/t - 2/T + T*d1 + d2 ]
Rx1=sp.Rational(1,2)*(1/t - 2/T + T*d1 + d2)
# t^2 * Shat = t^2*[40t^2-18tT+(108t^3-1)T^2]/(256 t^5) = [40t^2-18tT+(108t^3-1)T^2]/(256 t^3)
lhs=sp.simplify(Rx1)
rhs=sp.simplify((40*t**2-18*t*T+(108*t**3-1)*T**2)/(256*t**3))
diff=sp.simplify(lhs-rhs)
print("(ii) [x^1]R == t^2 * S-GF closed form (rational identity in T):", diff==0)

# (iii) independent check: [x^1]R formula against DP series
x,tt=sp.symbols('x tt')
N=16
Ts=sp.Integer(0)
for _ in range(N+3):
    Ts=sp.expand(sp.series(tt*(2+Ts**3), tt, 0, N+3).removeO())
val=sp.expand(sp.series(Rx1.subs(T,Ts).subs(t,tt).doit() if False else (sp.Rational(1,2)*( (2+Ts**3)/Ts - 2/Ts + Ts*(-Ts*(1+Ts**3/4)/2) + sp.Rational(1,8)*Ts**2*(1-(1+Ts**3/4)**2))), tt,0,N).removeO())
dpv=sum(RK(m,1,0)*tt**(m+1) for m in range(N))
d3=sp.expand(val-sp.expand(sp.series(dpv/tt,tt,0,N).removeO()))
# val should equal sum RK(m,1,0) t^m ... [x^1]R = t*sum RK(m,1,0)t^m / t = coefficient series: t*Q1: [x^1]R(x)=t*R1
# recompute cleanly: [x^1]R should equal t*sum_m RK(m,1,0)t^m
d3=sp.expand(val - sp.expand(sp.series(tt*sum(RK(m,1,0)*tt**m for m in range(N)),tt,0,N).removeO()))
print("(iii) [x^1]R formula matches DP (t-series):", all(sp.simplify(d3.coeff(tt,m))==0 for m in range(N-1)))
