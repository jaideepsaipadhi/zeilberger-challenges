from functools import lru_cache
import sympy as sp, json
import sys
sys.setrecursionlimit(100000)

# reverse-Kreweras (our completions, steps E,N,D from origin) by endpoint & length
@lru_cache(maxsize=None)
def RK(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    # last step: E from (i-1,j); N from (i,j-1); D from (i+1,j+1)
    return RK(m-1,i-1,j)+RK(m-1,i,j-1)+RK(m-1,i+1,j+1)

# BM Kreweras (steps NE, W, S from origin)
@lru_cache(maxsize=None)
def KW(m,i,j):
    if i<0 or j<0 or m<0: return 0
    if m==0: return 1 if (i==0 and j==0) else 0
    # last step: NE from (i-1,j-1); W from (i+1,j); S from (i,j+1)
    return KW(m-1,i-1,j-1)+KW(m-1,i+1,j)+KW(m-1,i,j+1)

# orientation test on the diagonal
same=all(RK(m,a,a)==KW(m,a,a) for a in range(0,7) for m in range(0,16))
print("diagonal counts equal between Kreweras and reverse-Kreweras:", same)
if not same:
    for a in range(0,4):
        for m in range(0,10):
            if RK(m,a,a)!=KW(m,a,a): print("  first diff:",m,a,RK(m,a,a),KW(m,a,a)); break

# BM Theorem 2: t*Qd(x;t) = 2tZ - (1/x)*sqrt(1 - x t Z(1+Z) + x^2 t^2 Z^2) + 1/x,  Z=1+4t^3 Z^3
t,x=sp.symbols('t x')
NT=28
Z=sp.Integer(1)
for _ in range(NT):
    Z=sp.expand(sp.series(1+4*t**3*Z**3, t, 0, NT).removeO())
    Z=sp.Poly(Z,t).as_expr()
inner=1 - x*t*Z*(1+Z) + x**2*t**2*Z**2
sq=sp.series(sp.sqrt(inner), x, 0, 8).removeO()   # expand sqrt in x first
Qd=sp.expand((2*t*Z - (sq-1)/x)/t)
Qd=sp.series(Qd, t, 0, NT-2).removeO()
# check Qd coefficients against KW diagonal
ok=True
for a in range(0,6):
    for m in range(0,NT-4):
        cf=Qd.coeff(x,a).coeff(t,m)
        if cf!=KW(m,a,a): ok=False; print("Qd mismatch at a,m:",a,m,cf,KW(m,a,a)); break
print("BM Theorem 2 series matches Kreweras diagonal counts:", ok)

# Lemma B test: T(n) = [t^{3n}] Qd(t;t)
Qdt=sp.expand(Qd.subs(x,t))
Qdt=sp.series(Qdt,t,0,NT-4).removeO()
P=json.load(open('pieces.json'))
T=P['T']
okT=all(Qdt.coeff(t,3*n)==T[n] for n in range(0,(NT-5)//3))
print("T(n) = [t^{3n}] Qd(t;t) for available n:", okT,
      "   tested n range:", list(range(0,(NT-5)//3)))
