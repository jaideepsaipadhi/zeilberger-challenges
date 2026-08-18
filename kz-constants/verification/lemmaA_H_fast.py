"""Optimized completion of lemmaA_H.py's Perron cross-check (odd-run H model).

The original script computes det(P - lam*I) with entries containing
exp(I e (u.v)) and sqrt(2) radicals, then series-expands: the determinant of
symbolic exponentials is the bottleneck (hours). Here we expand each entry to
O(e^3) FIRST, so the determinant is a polynomial computation (seconds), and
the mathematics checked is identical:

  tilted odd-run law P(R=r) = w^{r-1}(1-w^2), r odd, w = sqrt(2)-1;
  E[z^R] = (1-w^2) z / (1 - w^2 z^2);
  transfer matrix P(e)_{ij} = (1/2) phi_H(v_j) for i != j;
  Perron root lam(e) = 1 + c1 e + c2 e^2 + O(e^3):
    * c1 = 0 (zero position drift),
    * quadratic form c2 - c1^2/2 must equal -(1/2) s_H (u^T M u) with
      s_H = E[R^2] - (2/3) E[R]^2 = 5/3 and M = [[2/3,-1/3],[-1/3,2/3]],
      i.e. -(1/2)*(5/3)*(2/3)*(u1^2 - u1 u2 + u2^2).
"""
import sympy as sp

w = sp.sqrt(2) - 1
e, u1, u2 = sp.symbols('e u1 u2', real=True)
I = sp.I
v = [sp.Matrix([1, 0]), sp.Matrix([-1, 1]), sp.Matrix([0, -1])]

def phiH_series(vec):
    t = u1 * vec[0] + u2 * vec[1]
    z = sp.exp(I * e * t)
    f = (1 - w**2) * z / (1 - w**2 * z**2)
    # truncate the ENTRY to O(e^3) before any matrix work
    return sp.expand(sp.series(f, e, 0, 3).removeO())

P = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        if i != j:
            P[i, j] = sp.Rational(1, 2) * phiH_series(v[j])

lam, c1, c2 = sp.symbols('lam c1 c2')
cp = sp.det(P - lam * sp.eye(3))
expr = sp.expand(cp.subs(lam, 1 + c1 * e + c2 * e**2))
expr = sp.expand(sp.series(expr, e, 0, 3).removeO())

c1v = sp.solve(sp.expand(expr.coeff(e, 1)), c1)[0]
c1v = sp.simplify(sp.radsimp(c1v))
print("drift order-1:", c1v, "(expect 0)")

c2v = sp.solve(sp.expand(expr.coeff(e, 2)).subs(c1, c1v), c2)[0]
quad = sp.simplify(sp.radsimp(sp.expand(c2v - c1v**2 / 2)))
print("Perron quadratic:", sp.nsimplify(quad, [sp.sqrt(2)]))

sH = sp.Rational(5, 3)   # = E[R^2] - (2/3) E[R]^2 = 3 - (2/3)*2, run-law moments
target = sp.expand(-sp.Rational(1, 2) * sH *
                   (sp.Rational(2, 3) * u1**2 - sp.Rational(2, 3) * u1 * u2 +
                    sp.Rational(2, 3) * u2**2))
print("Perron Hessian == s_H*M (s_H=5/3):", sp.simplify(quad - target) == 0)
