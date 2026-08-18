import sympy as sp, json
from fractions import Fraction
t=sp.symbols('t')
P=json.load(open('pieces.json')); Sv=P['S']
NS=len(Sv)  # 47 values: S(0..46)
NT=3*(NS-1)+1  # t-order available

# S-GF series
SGF=sum(Sv[n]*t**(3*n) for n in range(NS))

# T and sigma series to matching order
T=sp.Integer(0)
for _ in range(60):
    T=sp.expand(sp.series(t*(2+T**3), t, 0, min(NT,60)).removeO())
NTT=min(NT,60)
sg=sp.expand(sp.sqrt(1-t*T**2).series(t,0,NTT).removeO())

# ansatz: q(t)*SGF = sum_{i<3,j<2} p_ij(t) T^i sg^j, deg(q)<=D, deg(p)<=D
def tryfit(D):
    unk=[]
    q=[sp.Symbol('q%d'%k) for k in range(D+1)]
    ps={}
    for i in range(3):
        for j in range(2):
            ps[(i,j)]=[sp.Symbol('p%d%d_%d'%(i,j,k)) for k in range(D+1)]
    qpoly=sum(q[k]*t**k for k in range(D+1))
    rhs=sum(sum(ps[(i,j)][k]*t**k for k in range(D+1))*T**i*sg**j for i in range(3) for j in range(2))
    expr=sp.expand(sp.series(qpoly*SGF-rhs, t, 0, NTT-1).removeO())
    eqs=[expr.coeff(t,m) for m in range(NTT-1)]
    allunk=q+[u for v in ps.values() for u in v]
    sol=sp.linsolve(eqs, allunk)
    if not sol: return None
    solv=list(sol)[0]
    if all(v==0 for v in solv): return None
    # check nontrivial q
    subs=dict(zip(allunk,solv))
    free=[s for s in solv if s.free_symbols]
    # pick a particular solution: set free params
    rep={}
    fs=set()
    for v in solv: fs|=v.free_symbols
    fs=list(fs)
    if not fs: 
        vals=solv
    else:
        rep={f:(1 if k==0 else 0) for k,f in enumerate(fs)}
        vals=[v.subs(rep) for v in solv]
    if all(v==0 for v in vals): return None
    return dict(zip(allunk,vals))

for D in range(2,13):
    r=tryfit(D)
    if r:
        print("TOWER FIT FOUND at degree D=%d"%D)
        json.dump({str(k):sp.srepr(sp.nsimplify(v)) for k,v in r.items()}, open('sgf_tower.json','w'))
        # display
        qpoly=sum(r[sp.Symbol('q%d'%k)]*t**k for k in range(D+1))
        print("q(t) =", sp.factor(qpoly))
        for i in range(3):
            for j in range(2):
                pp=sum(r[sp.Symbol('p%d%d_%d'%(i,j,k))]*t**k for k in range(D+1))
                if pp!=0: print("p_%d%d(t) =: "%(i,j), sp.factor(pp))
        break
else:
    print("no fit with {1,T,T^2}x{1,sigma} up to deg 12 — S-GF not in this tower (or needs bigger basis)")
