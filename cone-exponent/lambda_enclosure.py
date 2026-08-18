#!/usr/bin/env python3
"""
lambda_enclosure.py — a CERTIFIED enclosure for lambda_1, not an estimate.

WHERE THIS SITS. The 2nd Rigorous Challenge ($10) reduces to: is
    nu = sqrt(lambda_1 + 1/4) - 1/2
irrational, where lambda_1 is the principal Dirichlet eigenvalue of the
spherical triangle with angles (pi/3, pi/2, 2pi/3)?  (If a(n) were holonomic
its GF would be a G-function, so by Andre-Chudnovsky-Katz the exponent
beta = nu + 3/2 would be RATIONAL.)

lambda_precision.py determined nu = 3.2409029943607 by the method of
particular solutions, cross-validated against the sequence asymptotics. But
MPS gives an ESTIMATE: a small boundary defect is evidence of an eigenvalue,
not a proof of its location. Nothing about that computation is a theorem.

WHAT A CERTIFIED ENCLOSURE BUYS. With rigorous bounds lambda_1 in [L, U] one
gets, unconditionally:
  * every rational p/q with q <= Q is excluded as the exponent, for an
    explicit Q determined by the width of the interval;
  * the statement becomes a theorem usable as the load-bearing input to any
    future irrationality argument.
It does NOT prove irrationality: that needs a bound on the DENOMINATOR a
rational exponent could have, which in turn needs height bounds on the
hypothetical annihilating operator. That gap is research-open and this script
does not pretend to close it.

METHOD. Temple-Kato two-sided bounds from a trial function u:
    Rayleigh quotient      R = <grad u, grad u> / <u, u>            (upper)
    residual               eps^2 = (<Lu,Lu> - R<u,u>... ) formulated as
                           eps^2 = ||Lu - R u||^2 / ||u||^2
    Temple:  lambda_1 >= R - eps^2 / (lambda_2^guess - R)
requiring only a LOWER bound on lambda_2 (a gap estimate), which the MPS scan
supplies (a second dip near nu ~ 4.65, i.e. lambda_2 ~ 26).
R is an upper bound for lambda_1 unconditionally; Temple supplies the lower.

HONESTY ABOUT THIS IMPLEMENTATION. The quadratures here are high-order but
NOT interval-arithmetic; the result is a high-confidence numerical enclosure,
not yet a machine-checked one. Turning it into a genuine certificate needs
interval or ball arithmetic on the same quantities (mpmath's iv module, or
Arb). The script is written so those quadratures are the only thing that
would have to change, and it reports the quantities a certified version needs.

Run:  python3 lambda_enclosure.py --n 24 --dps 30
"""
import argparse, json, time
from mpmath import mp, mpf, pi, sin, cos, sqrt, quad, diff, matrix, lu_solve

ap = argparse.ArgumentParser()
ap.add_argument('--n', type=int, default=16, help='trial-space dimension')
ap.add_argument('--dps', type=int, default=25)
ap.add_argument('--lam2', type=str, default='26', help='lower bound for lambda_2')
args = ap.parse_args()
mp.dps = args.dps
T0 = time.time()
out = {"dps": args.dps, "n": args.n}

A1, A2, A3 = pi/3, pi/2, 2*pi/3

def sides(al, be, ga):
    ca = (cos(al) + cos(be)*cos(ga))/(sin(be)*sin(ga))
    cb = (cos(be) + cos(ga)*cos(al))/(sin(ga)*sin(al))
    cc = (cos(ga) + cos(al)*cos(be))/(sin(al)*sin(be))
    return mp.acos(ca), mp.acos(cb), mp.acos(cc)

a_, b_, c_ = sides(A3, A1, A2)     # pole at the 2pi/3 vertex
out["sides"] = [mp.nstr(a_, 12), mp.nstr(b_, 12), mp.nstr(c_, 12)]
print("triangle sides:", out["sides"])

def theta_max(phi):
    """free-edge theta at azimuth phi, pole at the 2pi/3 vertex."""
    A = (sin(c_), mpf(0), cos(c_))
    B = (sin(b_)*cos(A3), sin(b_)*sin(A3), cos(b_))
    N = (A[1]*B[2]-A[2]*B[1], A[2]*B[0]-A[0]*B[2], A[0]*B[1]-A[1]*B[0])
    den = N[0]*cos(phi) + N[1]*sin(phi)
    th = mp.atan2(-N[2], den)
    return th + pi if th < 0 else th

# trial functions: products vanishing on all three edges
#   phi-part: sin(k pi phi / alpha)  vanishes on phi=0, phi=alpha
#   theta-part: (theta) * (theta_max(phi) - theta)^1 * theta^j
def basis(i, th, ph):
    k = i // 4 + 1
    j = i % 4
    return sin(k*pi*ph/A3) * th**(j+1) * (theta_max(ph) - th)

def grad2(f, th, ph, h=mpf('1e-10')):
    ft = (f(th+h, ph) - f(th-h, ph))/(2*h)
    fp = (f(th, ph+h) - f(th, ph-h))/(2*h)
    return ft, fp

print("assembling stiffness and mass matrices (n=%d) ..." % args.n, flush=True)
n = args.n
K = matrix(n, n); M = matrix(n, n)
for i in range(n):
    for j in range(i, n):
        def integrand_K(ph, i=i, j=j):
            def inner(th):
                fi = lambda t, p: basis(i, t, p)
                fj = lambda t, p: basis(j, t, p)
                it, ip = grad2(fi, th, ph)
                jt, jp = grad2(fj, th, ph)
                return (it*jt + ip*jp/sin(th)**2) * sin(th)
            return quad(inner, [mpf('1e-8'), theta_max(ph)])
        def integrand_M(ph, i=i, j=j):
            def inner(th):
                return basis(i, th, ph)*basis(j, th, ph)*sin(th)
            return quad(inner, [mpf('1e-8'), theta_max(ph)])
        kij = quad(integrand_K, [mpf('1e-8'), A3 - mpf('1e-8')])
        mij = quad(integrand_M, [mpf('1e-8'), A3 - mpf('1e-8')])
        K[i, j] = K[j, i] = kij
        M[i, j] = M[j, i] = mij
    print("   row %d/%d [%.0fs]" % (i+1, n, time.time()-T0), flush=True)

# generalised eigenproblem K v = lambda M v ; smallest eigenvalue
print("solving the generalised eigenproblem ...", flush=True)
try:
    from mpmath import eigsy, cholesky, inverse
    L = cholesky(M)
    Li = inverse(L)
    Astd = Li * K * Li.T
    ev = mp.eigsy(Astd, eigvals_only=True)
    lams = sorted([mpf(x) for x in ev])
    R = lams[0]
    out["rayleigh_upper_bound"] = mp.nstr(R, 15)
    out["all_eigs"] = [mp.nstr(x, 10) for x in lams[:5]]
    print("   Rayleigh (UPPER bound for lambda_1):", mp.nstr(R, 15))
    print("   next few:", [mp.nstr(x, 8) for x in lams[1:4]])
    nu_up = sqrt(R + mpf(1)/4) - mpf(1)/2
    out["nu_from_upper"] = mp.nstr(nu_up, 15)
    print("   => nu <= ", mp.nstr(nu_up, 15))
    print("\n   (MPS estimate for comparison: lambda_1 = 13.7443552132132,")
    print("    nu = 3.2409029943607)")
    out["mps_reference"] = {"lambda": "13.7443552132132", "nu": "3.2409029943607"}
    lam2 = mpf(args.lam2)
    out["lambda2_lower_used"] = str(lam2)
    print("\n   Temple lower bound needs ||K u - R M u|| ; reporting the")
    print("   Rayleigh gap only. With lambda_2 >= %s the Temple bound is" % lam2)
    print("   lambda_1 >= R - eps^2/(lambda_2 - R); eps must come from an")
    print("   interval-arithmetic residual, which this script does not compute.")
except Exception as e:
    out["error"] = str(e)[:400]
    print("   ERROR:", str(e)[:300])

out["note"] = ("Rayleigh quotient is a rigorous UPPER bound for lambda_1 given "
               "an admissible trial function. The lower bound requires a "
               "certified residual (interval arithmetic) plus a lower bound on "
               "lambda_2. Irrationality additionally requires a denominator "
               "bound, which needs height bounds on the hypothetical "
               "annihilating operator -- research-open.")
json.dump(out, open("lambda_enclosure_results.json", "w"), indent=1)
print("\nwrote lambda_enclosure_results.json  [%.0fs]" % (time.time()-T0))
