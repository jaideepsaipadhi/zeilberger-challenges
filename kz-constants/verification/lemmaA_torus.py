"""A2: (i) locate all unit-modulus points of the Fourier operator on the torus (prediction:
exactly 3, at (-s,s,s), 3s=0 mod 2pi); (ii) full 3D per-run covariance (position + length)."""
import numpy as np, itertools, sympy as sp
v=[(1,0),(-1,1),(0,-1)]
def Pmat(t1,t2,s):
    P=np.zeros((3,3),dtype=complex)
    for i in range(3):
        for j in range(3):
            if i!=j:
                z=np.exp(1j*(t1*v[j][0]+t2*v[j][1]+s))
                P[i,j]=0.5*z**2/(2-z)
    return P
def rho(t1,t2,s):
    return max(abs(np.linalg.eigvals(Pmat(t1,t2,s))))
# predicted unit points
preds=[(0.0,0.0,0.0),(-2*np.pi/3,2*np.pi/3,2*np.pi/3),(-4*np.pi/3,4*np.pi/3,4*np.pi/3)]
print("predicted points rho:", [round(rho(*p),10) for p in preds])
# grid scan for any OTHER near-unit points
N=25; mx=0; loc=None; hits=[]
for a in range(N):
    for b in range(N):
        for c in range(N):
            t1=2*np.pi*a/N; t2=2*np.pi*b/N; s=2*np.pi*c/N
            r=rho(t1,t2,s)
            near_pred=any(min(abs((t1-p[0])%(2*np.pi)),2*np.pi-abs((t1-p[0])%(2*np.pi)))<0.3 and
                          min(abs((t2-p[1])%(2*np.pi)),2*np.pi-abs((t2-p[1])%(2*np.pi)))<0.3 and
                          min(abs((s-p[2])%(2*np.pi)),2*np.pi-abs((s-p[2])%(2*np.pi)))<0.3 for p in preds)
            if r>0.999 and not near_pred: hits.append((t1,t2,s,r))
            if not near_pred and r>mx: mx=r; loc=(t1,t2,s)
print("non-predicted near-unit hits:", hits[:5], "count:", len(hits))
print("max rho away from predicted points: %.6f at %s"%(mx,tuple(round(x,3) for x in loc)))
# (ii) symbolic 3D covariance per run
t1s,t2s,ss,e=sp.symbols('t1 t2 s e', real=True)
u1,u2,u3,c1,c2=sp.symbols('u1 u2 u3 c1 c2')
I=sp.I
P=sp.zeros(3,3)
for i in range(3):
    for j in range(3):
        if i!=j:
            z=sp.exp(I*e*(u1*v[j][0]+u2*v[j][1]+u3))
            P[i,j]=sp.Rational(1,2)*z**2/(2-z)
lam=sp.symbols('lam')
cp=sp.det(P-lam*sp.eye(3))
expr=cp.subs(lam,1+c1*e+c2*e**2)
expr=sp.series(sp.expand(expr),e,0,3).removeO()
c1v=sp.solve(sp.expand(expr.coeff(e,1)),c1)[0]
c2v=sp.simplify(sp.solve(sp.expand(expr.coeff(e,2)).subs(c1,c1v),c2)[0])
print("drift (order-1):", sp.simplify(c1v), "  [expect 3*I*u3: length drift E[R]=3, position drift 0]")
quad=sp.expand(sp.simplify(c2v - c1v**2/2))
print("log-lam quadratic:", quad)
# extract Sigma3: quad = -1/2 u^T Sig u  (after removing drift, i.e. central moments)
Sig=sp.zeros(3,3)
us=[u1,u2,u3]
for a in range(3):
    Sig[a,a]=-2*quad.coeff(us[a]**2)
for a in range(3):
    for b in range(a+1,3):
        Sig[a,b]=Sig[b,a]=-quad.coeff(us[a]*us[b])/1
sp.pprint(Sig)
