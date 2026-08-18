import sympy as sp, json
t,T,sg=sp.symbols('t T sigma')

# reduction: T^3 = (T-2t)/t ; sigma^2 = 1 - t*T^2
def reduce_expr(e):
    e=sp.expand(e)
    # reduce powers of T >= 3 and sigma >= 2 iteratively
    changed=True
    while changed:
        changed=False
        p=sp.Poly(e,T,sg)
        out=sp.Integer(0)
        for (i,j),cf in zip(p.monoms(),p.coeffs()):
            ii,jj=i,j
            fac=sp.sympify(cf)
            while ii>=3:
                fac=fac*(T-2*t)/t; ii-=3; changed=True
            while jj>=2:
                fac=fac*(1-t*T**2); jj-=2; changed=True
            out+=sp.expand(fac)*T**ii*sg**jj
        e=sp.expand(out)
    return e

Tp=sp.cancel((2+T**3)/(1-3*t*T**2))           # T'
def D(e):
    de=sp.diff(e,t)+sp.diff(e,T)*Tp+sp.diff(e,sg)*(-(T**2+2*t*T*Tp)/(2*sg))
    # clear sigma in denominator: multiply num/den... handle via together and replace 1/sg = sg/(1-tT^2)
    de=sp.together(de)
    num,den=sp.fraction(de)
    # den may contain sg and T powers; rationalize sigma: multiply num&den by sg if den has odd sg
    dpoly=sp.Poly(den, sg)
    if any(m[0]%2==1 for m in dpoly.monoms()):
        num=sp.expand(num*sg); den=sp.expand(den*sg)
    den=reduce_expr(den)   # now den should be sigma-free polynomial in T,t
    num=reduce_expr(num)
    # den still contains T; rationalize T: multiply by conjugates? Instead solve: we want num/den in basis; 
    # invert den in Q(t)[T]/(tT^3-T+2t): find inv via linear solve (3x3)
    dp=sp.Poly(den, T)
    assert sp.Poly(den,sg).degree()==0 or den.free_symbols.isdisjoint({sg}), "sigma left in denominator"
    # multiplication-by-den matrix in basis 1,T,T^2 over Q(t)
    M=sp.zeros(3,3)
    for j in range(3):
        col=reduce_expr(den*T**j)
        cp=sp.Poly(col,T)
        for i in range(3):
            M[i,j]=col.coeff(T,i)
    Minv=M.inv()
    # num = n0(sg-free part) + n1*sg ; invert den against each
    res=sp.Integer(0)
    for jsg in range(2):
        comp=sp.expand(num.coeff(sg,jsg))
        v=sp.Matrix([comp.coeff(T,i) for i in range(3)])
        w=Minv*v
        res+=sum(sp.cancel(w[i])*T**i for i in range(3))*sg**jsg
    return sp.expand(res)

# A = [T/(2t) + (sg - 1 + t T^2/2)/t^3] / sg  -> rationalize: * sg/sg
Araw=(T/(2*t) + (sg-1+t*T**2/2)/t**3)*sg  # = A*sg^2 = A*(1-tT^2)... careful
# A = N/sg with N = T/(2t)+(sg-1+tT^2/2)/t^3 ; A = N*sg/(1-tT^2)
N=T/(2*t)+(sg-1+t*T**2/2)/t**3
A=sp.cancel(sp.together(N*sg/(1-t*T**2)))
num,den=sp.fraction(A)
num=reduce_expr(sp.expand(num)); den=reduce_expr(sp.expand(den))
# put A in basis via same inversion trick
def tobasis(num,den):
    dpoly=sp.Poly(den,sg)
    if any(m[0]%2==1 for m in dpoly.monoms()):
        num=sp.expand(num*sg); den=sp.expand(den*sg)
    den=reduce_expr(den); num=reduce_expr(num)
    M=sp.zeros(3,3)
    for j in range(3):
        col=reduce_expr(den*T**j)
        for i in range(3): M[i,j]=col.coeff(T,i)
    Minv=M.inv()
    res=sp.Integer(0)
    for jsg in range(2):
        comp=sp.expand(num.coeff(sg,jsg))
        v=sp.Matrix([comp.coeff(T,i) for i in range(3)])
        w=Minv*v
        res+=sum(sp.cancel(w[i])*T**i for i in range(3))*sg**jsg
    return sp.expand(res)
A0=tobasis(num,den)
print("A in basis computed")
A1=D(A0); print("A' done")
A2=D(A1); print("A'' done")

def vec(e):
    return [sp.cancel(sp.expand(e).coeff(T,i).coeff(sg,j)) for i in range(3) for j in range(2)]
one=vec(sp.Integer(1)); vA=vec(A0); vA1=vec(A1); vA2=vec(A2)
Mm=sp.Matrix([[one[i],vA[i],vA1[i],vA2[i]] for i in range(6)])
ns=Mm.nullspace()
print("order-2 (with const) dependency:", len(ns)>0)
if ns:
    v=ns[0]
    den=sp.lcm([sp.denom(sp.together(z)) for z in v])
    v=[sp.expand(sp.cancel(sp.together(z)*den)) for z in v]
    g0=sp.gcd([z for z in v if z!=0]); v=[sp.expand(sp.cancel(z/g0)) for z in v]
    json.dump([sp.srepr(z) for z in v], open('ode_A.json','w'))
    for k,lbl in enumerate(["1","A","A'","A''"]):
        print(" coeff[%s] = %s"%(lbl, sp.factor(v[k])))
else:
    A3=D(A2); print("A''' done")
    vA3=vec(A3)
    Mm=sp.Matrix([[one[i],vA[i],vA1[i],vA2[i],vA3[i]] for i in range(6)])
    ns=Mm.nullspace()
    print("order-3 dependency:", len(ns)>0)
    if ns:
        v=ns[0]
        den=sp.lcm([sp.denom(sp.together(z)) for z in v])
        v=[sp.expand(sp.cancel(sp.together(z)*den)) for z in v]
        g0=sp.gcd([z for z in v if z!=0]); v=[sp.expand(sp.cancel(z/g0)) for z in v]
        json.dump([sp.srepr(z) for z in v], open('ode_A.json','w'))
        for k,lbl in enumerate(["1","A","A'","A''","A'''"]):
            print(" coeff[%s] deg %s"%(lbl, sp.degree(v[k],t) if v[k]!=0 else "-"))
