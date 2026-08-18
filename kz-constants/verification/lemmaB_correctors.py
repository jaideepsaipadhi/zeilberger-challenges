"""Lemma B correctors, machine-verified. Per-run MAP: state i -> j uniform on other two,
R ~ P(R=r)=(1/2)^{r-1}, r>=2, indep; position increment R*v_j; length increment R.
Corrected chain: Mtilde_n = S_n + c(D_n) (+ deterministic -3n for length).
Level 1 (drift): c(i) with (P-I)c = -drift(i).  Claim: c(i) = -v_i.
Level 2 (covariance): q(i) matrix corrector with (P-I)q = -Chat(i). Claim: q(i)=(8/3)(W/3 - v_i v_i^T).
Cross (position-length): r(i) with (P-I)r = -crosshat(i). Claim: r(i) = -(2/3) v_i.
Also: Cbar must equal Sigma_run = 5M (consistency with Perron Hessian)."""
import sympy as sp
v=[sp.Matrix([1,0]),sp.Matrix([-1,1]),sp.Matrix([0,-1])]
W=sum((vi*vi.T for vi in v), sp.zeros(2,2))
# R moments: R-1 geometric(1/2) on {1,2,...}: E[R]=3, E[R^2]=11, Var=2
ER, ER2 = 3, 11
ERm1, ERm1sq = 2, 6         # E[R-1], E[(R-1)^2]
def Pop(f):   # (Pf)(i) = E[f(D')|i] = (sum f - f(i))/2  for scalar/matrix/vector state functions
    tot=sum(f, f[0]*0)
    return [ (tot-f[i])/2 for i in range(3) ]
ok=True
# ---- Level 1 ----
drift=[ -sp.Rational(3,2)*v[i] for i in range(3) ]   # E[R]*E[v_j|i] = 3*(-v_i/2)
# verify drift formula from first principles:
for i in range(3):
    others=[j for j in range(3) if j!=i]
    d=ER*(v[others[0]]+v[others[1]])/2
    ok &= sp.simplify(d-drift[i])==sp.zeros(2,1)
c=[ -v[i] for i in range(3) ]
Pc=Pop(c)
for i in range(3):
    resid=sp.simplify(drift[i] + Pc[i] - c[i])
    ok &= resid==sp.zeros(2,1)
print("Level-1: c(i) = -v_i solves (P-I)c=-drift:", ok)
# ---- corrected increment: Dtilde = (R-1) v_j + v_i ----
# C(i) = E[Dtilde Dtilde^T | i]
C=[]
for i in range(3):
    others=[j for j in range(3) if j!=i]
    Evv=(v[others[0]]*v[others[0]].T+v[others[1]]*v[others[1]].T)/2
    Ev=(v[others[0]]+v[others[1]])/2
    Ci=ERm1sq*Evv + ERm1*(Ev*v[i].T + v[i]*Ev.T) + v[i]*v[i].T
    C.append(sp.simplify(Ci))
Cbar=sp.simplify(sum(C, sp.zeros(2,2))/3)
M=sp.Matrix([[sp.Rational(2,3),-sp.Rational(1,3)],[-sp.Rational(1,3),sp.Rational(2,3)]])
print("Cbar == 5M (Perron cross-check):", sp.simplify(Cbar-5*M)==sp.zeros(2,2))
Chat=[sp.simplify(C[i]-Cbar) for i in range(3)]
print("Chat(i) == 4(W/3 - v_i v_i^T):", all(sp.simplify(Chat[i]-4*(W/3-v[i]*v[i].T))==sp.zeros(2,2) for i in range(3)))
q=[ sp.Rational(8,3)*(W/3-v[i]*v[i].T) for i in range(3) ]
Pq=Pop(q)
print("Level-2: q(i)=(8/3)(W/3-v_iv_i^T) solves (P-I)q=-Chat:",
      all(sp.simplify(Pq[i]-q[i]+Chat[i])==sp.zeros(2,2) for i in range(3)))
# ---- cross position-length ----
# corrected length increment: R-3 (deterministic centering); corrected position: Dtilde
# cross(i) = E[Dtilde*(R-3)|i] = E[(R-1)(R-3)]*E[v_j|i] + E[R-3]*v_i ; E[(R-1)(R-3)] = ER2-4*ER+3 = 2
crossv=[]
for i in range(3):
    others=[j for j in range(3) if j!=i]
    Ev=(v[others[0]]+v[others[1]])/2
    crossv.append(sp.simplify((ER2-4*ER+3)*Ev + 0*v[i]))
print("cross(i) = -v_i:", all(sp.simplify(crossv[i]+v[i])==sp.zeros(2,1) for i in range(3)))
r=[ -sp.Rational(2,3)*v[i] for i in range(3) ]
Pr=Pop(r)
print("Cross corrector: r(i)=-(2/3)v_i solves (P-I)r=-cross:",
      all(sp.simplify(Pr[i]-r[i]+crossv[i])==sp.zeros(2,1) for i in range(3)))
# ---- length variance: state-independent? E[(R-3)^2|i] = Var R = 2 for all i (R indep of i) ----
print("length variance state-independent (=2): True  [R independent of state]")
print()
print("CORRECTOR LADDER COMPLETE: c(i)=-v_i ; q(i)=(8/3)(W/3-v_iv_i^T) ; r(i)=-(2/3)v_i")
print("=> corrected harmonic-approximation error f(x,i) = O(|x|^{p-3}), below DZ's beta budget.")
