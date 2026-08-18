import sympy as sp
from math import factorial, comb
from itertools import product

# variables t2..t6 (5D challenge)
KS=[2,3,4,5,6]
ts=sp.symbols('t2 t3 t4 t5 t6')

# (1) verify hyper-Catalan closed form against the series S = 1 + sum_k t_k S^k
def C_closed(m):   # m = dict k->m_k
    top=sum(k*m[k] for k in m)
    bot=1+sum((k-1)*m[k] for k in m)
    v=factorial(top)//factorial(bot)
    for k in m: v//=factorial(m[k])
    return v

N=3   # truncation: total degree in each variable
S=sp.Integer(1)
for it in range(12):
    S=sp.expand(1+sum(ts[i]*S**KS[i] for i in range(5)))
    S=sp.Poly(S,*ts).as_expr()
    # truncate
    P=sp.Poly(S,*ts)
    S=sum(c*sp.prod([ts[i]**mon[i] for i in range(5)]) for mon,c in zip(P.monoms(),P.coeffs()) if max(mon)<=N)
Spoly=sp.Poly(sp.expand(S),*ts)
ok=True
for mon,c in zip(Spoly.monoms(),Spoly.coeffs()):
    if sum(mon)==0: continue
    if max(mon)>N: continue
    m={KS[i]:mon[i] for i in range(5)}
    if int(c)!=C_closed(m): ok=False; print("MISMATCH",mon,c,C_closed(m)); break
print("hyper-Catalan closed form verified vs series (all monomials, each exponent <=%d): %s"%(N,ok))

# (2) G = (S-1)/(t2+...+t6): get G coefficients by exact polynomial division on the truncation
num=sp.expand(S-1); den=sum(ts)
q,r=sp.div(sp.Poly(num,*ts),sp.Poly(den,*ts))
print("division remainder zero:", r.as_expr()==0)
Gp=q
def G_series(mvec):
    return int(Gp.coeff_monomial(sp.prod([ts[i]**mvec[i] for i in range(5)])))

# (3) MY COLLAPSED FORMULA
# G[M] = (1/ (M!)^4 ... ) sum over alpha in prod{0..M_k}, k=3..6:
#   (-1)^{|a|} * |a|!/(prod a_k!) * C[M + (1+|a|)e2 - sum a_k e_k]
def G_alt(Mv):
    M2,M3,M4,M5,M6=Mv
    tot=0
    for a3 in range(M3+1):
     for a4 in range(M4+1):
      for a5 in range(M5+1):
       for a6 in range(M6+1):
        i=a3+a4+a5+a6
        m={2:M2+1+i,3:M3-a3,4:M4-a4,5:M5-a5,6:M6-a6}
        coef=factorial(i)//(factorial(a3)*factorial(a4)*factorial(a5)*factorial(a6))
        tot+=(-1)**i*coef*C_closed(m)
    return tot

# (4) THE (i,W) COLLAPSE: summand depends on alpha only via i=|a| and W=a3+2a4+3a5+4a6
def G_collapsed(Mv):
    M2,M3,M4,M5,M6=Mv
    Ms=[M3,M4,M5,M6]
    # table T[i][W] = sum over alpha with |a|=i, W(a)=W of prod C(M_k, a_k)
    from collections import defaultdict
    T=defaultdict(int); T[(0,0)]=1
    for idx,Mk in enumerate(Ms):
        wgt=idx+1
        T2=defaultdict(int)
        for (i,W),v in T.items():
            for a in range(Mk+1):
                T2[(i+a,W+wgt*a)]+=v*comb(Mk,a)
        T=T2
    tot=0
    for (i,W),v in T.items():
        top=sum(k*M for k,M in zip(KS,Mv)) + 2*(1+i) - W - 2*i   # recompute cleanly below
        # aggregates: sum k m_k = [2(M2+1+i) + 3(M3-a3)+4(M4-a4)+5(M5-a5)+6(M6-a6)]
        base_k=2*(M2+1)+3*M3+4*M4+5*M5+6*M6
        sum_k = base_k + 2*i - (W + 2*i)   # since 3a3+4a4+5a5+6a6 = W + 2i  (weights 3,4,5,6 = (1,2,3,4)+2)
        base_k1=(M2+1)+2*M3+3*M4+4*M5+5*M6
        sum_k1 = base_k1 + i - (W + i)     # 2a3+3a4+4a5+5a6 = W + i
        numer=factorial(sum_k)
        denom=factorial(1+sum_k1)*factorial(M2+1+i)
        # prod over k of (M_k - a_k)! : we folded 1/(a_k!(M_k-a_k)!) = C(M_k,a_k)/M_k!
        pref=1
        for Mk in Ms: pref*=factorial(Mk)
        tot+=(-1)**i*factorial(i)*v*numer//1 * 1
        # careful: assemble exactly
    return None  # replaced below

for Mv in [(1,1,1,1,1),(2,2,2,2,2),(2,1,2,1,2)]:
    print("M=%s: series G=%s, alternating-sum G=%s, match=%s"%(Mv,G_series(Mv),G_alt(Mv),G_series(Mv)==G_alt(Mv)))
