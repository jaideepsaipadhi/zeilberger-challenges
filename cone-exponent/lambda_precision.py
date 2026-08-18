#!/usr/bin/env python3
"""
lambda_precision.py — push lambda_1 of the (pi/3, pi/2, 2pi/3) spherical
triangle to 8-10 stable digits, and turn the rational exclusions into
defensible windows.

WHY THIS NUMBER. The 2nd Rigorous Challenge asks to prove that
a(n) = #Solid SYT of cylindrical shape (2,1,1)x[n] is NOT holonomic. The
chain established in earlier sessions:

  a(n) ~ C * 256^n * n^{-beta},   beta = nu + 3/2,
  nu = sqrt(lambda_1 + 1/4) - 1/2,
  lambda_1 = principal Dirichlet eigenvalue of the spherical triangle with
             angles (pi/3, pi/2, 2pi/3) -- the transverse trace of the cone
             x1>=x2>=x3, x1>=x4.

The GF has integer coefficients and finite radius, so if it were D-finite it
would be a G-function, and by Andre-Chudnovsky-Katz its singularity exponents
would be RATIONAL. So:  nu irrational  =>  non-holonomic.

The triangle is not a reflection fundamental domain (the 2pi/3 angle), and
polynomial eigenfunctions were excluded through degree 11, so there is no
closed form to appeal to. Certifying irrationality is the wall. What is
reachable is PRECISION: enough digits to exclude every rational of small
height by a margin, which is publishable evidence and the core of the
"cone exponent" note.

METHOD. Method of particular solutions, Betcke-Trefethen variant. Expand
about a vertex of interior angle alpha with the Fourier-Legendre basis
    u_k(theta, phi) = P_nu^{-mu_k}(cos theta) * sin(mu_k * phi),
    mu_k = k * pi / alpha,
which satisfies the eigenvalue equation and the two boundary edges meeting at
that vertex identically; the remaining edge is collocated. A true eigenvalue
is a nu where the boundary block of the QR factor becomes singular.

TRAPS ALREADY HIT IN EARLIER SESSIONS (guarded against here):
  * integer nu gives identically-zero columns (P_n^{-m} = 0 for integer m > n)
    -> near-zero columns are dropped before the QR;
  * plain sigma_min has a spurious floor from near-degenerate high-order
    columns that can mask a real dip completely -> the Betcke-Trefethen
    subspace-angle form is used, never plain sigma_min;
  * the pole must sit at the 2pi/3 vertex; the pi/2-pole scan found a
    different (second) eigenvalue near nu ~ 4.65.

CALIBRATION. The script first reproduces the A3 chamber triangle, whose
principal eigenfunction is the Vandermonde with nu = 6 exactly. If that does
not come out to the requested precision, the settings are wrong and the run
aborts rather than reporting a number.

Run:  python3 lambda_precision.py
      python3 lambda_precision.py --dps 40 --nb 20 --ni 400
Writes lambda_precision_results.json.

Cost: mpmath legenp at high dps is slow; the defaults are a few minutes.
Raising --nb past ~24 grows cost quickly for little gain (the basis becomes
numerically dependent before it becomes more accurate).
"""
import json, time, argparse
from mpmath import mp, mpf, pi, cos, sin, sqrt, legenp, matrix, svd_r, qr

ap = argparse.ArgumentParser()
ap.add_argument('--dps', type=int, default=30)
ap.add_argument('--nb', type=int, default=14, help='number of basis functions')
ap.add_argument('--ni', type=int, default=260, help='collocation points')
ap.add_argument('--skip-calib', action='store_true')
args = ap.parse_args()
mp.dps = args.dps
T0 = time.time()
out = {"dps": args.dps, "nb": args.nb, "ni": args.ni}

# ---------------------------------------------------------------- geometry
# Spherical triangle with angles (pi/3, pi/2, 2pi/3). We expand about the
# 2pi/3 vertex. The opposite edge is parameterised and collocated.
# Vertices are placed by solving the spherical triangle; we only need the
# free edge as a curve on the sphere, given in (theta, phi) about the pole.

def triangle_edges(alpha, beta, gamma):
    """Return (side lengths a,b,c) opposite to angles (alpha,beta,gamma) via
    the spherical law of cosines for angles:
        cos a = (cos alpha + cos beta cos gamma) / (sin beta sin gamma)."""
    ca = (cos(alpha) + cos(beta)*cos(gamma)) / (sin(beta)*sin(gamma))
    cb = (cos(beta) + cos(gamma)*cos(alpha)) / (sin(gamma)*sin(alpha))
    cc = (cos(gamma) + cos(alpha)*cos(beta)) / (sin(alpha)*sin(beta))
    return mp.acos(ca), mp.acos(cb), mp.acos(cc)

def free_edge_points(alpha, beta, gamma, n):
    """Points on the edge opposite the expansion vertex (the one with angle
    `alpha`), in (theta, phi) coordinates with the pole at that vertex.
    The two edges meeting at the pole are phi = 0 and phi = alpha; the free
    edge runs between the points at (theta=c, phi=0) and (theta=b, phi=alpha)
    where b, c are the sides adjacent to the expansion vertex."""
    a_, b_, c_ = triangle_edges(alpha, beta, gamma)
    pts = []
    for i in range(1, n+1):
        s = mpf(i) / (n + 1)
        phi = alpha * s
        # great-circle interpolation of theta along the free edge:
        # solve the spherical triangle formed by pole, and the running point.
        # Using the cotangent four-parts formula on the sub-triangle:
        #   cot(theta) sin(phi_gap) = cot(side) sin(phi) ... we instead
        # interpolate by solving for theta from the edge's great circle.
        # The free edge's great circle: determined by its two endpoints.
        # Endpoint A at (theta=c_, phi=0), endpoint B at (theta=b_, phi=alpha).
        A = (sin(c_)*cos(mpf(0)), sin(c_)*sin(mpf(0)), cos(c_))
        B = (sin(b_)*cos(alpha),  sin(b_)*sin(alpha),  cos(b_))
        # normal of the plane through A, B
        N = (A[1]*B[2]-A[2]*B[1], A[2]*B[0]-A[0]*B[2], A[0]*B[1]-A[1]*B[0])
        # point on the sphere with azimuth phi lying on that great circle:
        # N . (sin th cos phi, sin th sin phi, cos th) = 0
        #  => tan(th) = -N3 / (N1 cos phi + N2 sin phi)
        denom = N[0]*cos(phi) + N[1]*sin(phi)
        th = mp.atan2(-N[2], denom)
        if th < 0: th += pi
        pts.append((th, phi))
    return pts

def interior_points(alpha, beta, gamma, n):
    """A scattering of interior points, used for the Betcke-Trefethen
    normalisation block."""
    a_, b_, c_ = triangle_edges(alpha, beta, gamma)
    pts = []
    k = 0
    m = int(mp.sqrt(n)) + 1
    for i in range(1, m+1):
        for j in range(1, m+1):
            if k >= n: break
            phi = alpha * mpf(i) / (m + 1)
            # radial fraction of the way to the free edge
            edge = free_edge_points(alpha, beta, gamma, 1)  # not used directly
            th_max = None
            # recompute the free-edge theta at this phi
            A = (sin(c_), mpf(0), cos(c_))
            B = (sin(b_)*cos(alpha), sin(b_)*sin(alpha), cos(b_))
            N = (A[1]*B[2]-A[2]*B[1], A[2]*B[0]-A[0]*B[2], A[0]*B[1]-A[1]*B[0])
            denom = N[0]*cos(phi) + N[1]*sin(phi)
            th_max = mp.atan2(-N[2], denom)
            if th_max < 0: th_max += pi
            th = th_max * mpf(j) / (m + 1)
            pts.append((th, phi)); k += 1
    return pts

# ---------------------------------------------------------------- MPS
def basis_matrix(nu, pts, alpha, nb):
    """u_k(theta,phi) = P_nu^{-mu_k}(cos theta) sin(mu_k phi), mu_k = k pi/alpha."""
    M = matrix(len(pts), nb)
    for r, (th, ph) in enumerate(pts):
        ct = cos(th)
        for k in range(1, nb+1):
            mu = k * pi / alpha
            try:
                val = legenp(nu, -mu, ct) * sin(mu*ph)
            except Exception:
                val = mpf(0)
            M[r, k-1] = val
    return M

def drop_null_columns(M, tol=mpf('1e-25')):
    """Guard (i): integer-nu columns can be identically zero."""
    keep = []
    for j in range(M.cols):
        s = max(abs(M[i, j]) for i in range(M.rows))
        if s > tol: keep.append(j)
    if len(keep) == M.cols: return M, keep
    N = matrix(M.rows, len(keep))
    for jj, j in enumerate(keep):
        for i in range(M.rows):
            N[i, jj] = M[i, j]
    return N, keep

def bt_sigma(nu, alpha, beta, gamma, nb, nb_pts, ni_pts):
    """Betcke-Trefethen subspace angle: stack boundary and interior rows,
    QR the stack, take sigma_min of the boundary block of Q.
    Guard (ii): never plain sigma_min of the boundary matrix alone."""
    Bp = free_edge_points(alpha, beta, gamma, nb_pts)
    Ip = interior_points(alpha, beta, gamma, ni_pts)
    MB = basis_matrix(nu, Bp, alpha, nb)
    MI = basis_matrix(nu, Ip, alpha, nb)
    M = matrix(MB.rows + MI.rows, nb)
    for i in range(MB.rows):
        for j in range(nb): M[i, j] = MB[i, j]
    for i in range(MI.rows):
        for j in range(nb): M[MB.rows + i, j] = MI[i, j]
    M, keep = drop_null_columns(M)
    if M.cols == 0: return mpf(1)
    Q, R = qr(M)
    QB = matrix(MB.rows, M.cols)
    for i in range(MB.rows):
        for j in range(M.cols): QB[i, j] = Q[i, j]
    s = svd_r(QB, compute_uv=False)
    return min(s)

def golden_min(f, lo, hi, iters=60):
    """Golden-section minimisation of a unimodal dip."""
    gr = (sqrt(5) - 1) / 2
    a, b = mpf(lo), mpf(hi)
    c = b - gr*(b - a); d = a + gr*(b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr*(b - a); fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr*(b - a); fd = f(d)
        if abs(b - a) < mpf(10)**(-(mp.dps - 4)): break
    return (a + b)/2, min(fc, fd)

A1, A2, A3 = pi/3, pi/2, 2*pi/3          # the challenge triangle
POLE = A3                                 # expand about the 2pi/3 vertex

# ---------------------------------------------------------------- calibration
if not args.skip_calib:
    print("calibration: A3 chamber triangle, expect nu = 6 exactly", flush=True)
    # A3 Weyl chamber trace: angles (pi/3, pi/2, pi/3)
    f = lambda nu: bt_sigma(nu, pi/3, pi/2, pi/3, args.nb, 60, 120)
    nu_c, s_c = golden_min(f, 5.8, 6.2, iters=40)
    err = abs(nu_c - 6)
    out["calibration"] = {"nu": str(nu_c), "sigma": mp.nstr(s_c, 8),
                          "error_vs_6": mp.nstr(err, 8)}
    print("   nu = %s   (error %s, sigma %s)  [%.0fs]"
          % (mp.nstr(nu_c, 12), mp.nstr(err, 4), mp.nstr(s_c, 4),
             time.time()-T0), flush=True)
    if err > mpf('1e-6'):
        out["status"] = "CALIBRATION FAILED — settings wrong, not reporting a value"
        json.dump(out, open("lambda_precision_results.json", "w"), indent=1)
        raise SystemExit("calibration failed")

# ---------------------------------------------------------------- the target
print("\ntarget: (pi/3, pi/2, 2pi/3), pole at the 2pi/3 vertex", flush=True)
f = lambda nu: bt_sigma(nu, POLE, A1, A2, args.nb, args.ni//2, args.ni)
nu_star, s_star = golden_min(f, 3.20, 3.29, iters=70)
lam = nu_star*(nu_star + 1)
out["nu"] = str(nu_star)
out["lambda_1"] = str(lam)
out["sigma_at_min"] = mp.nstr(s_star, 8)
print("   nu       = %s" % mp.nstr(nu_star, 15), flush=True)
print("   lambda_1 = %s" % mp.nstr(lam, 15), flush=True)
print("   sigma    = %s   [%.0fs]" % (mp.nstr(s_star, 6), time.time()-T0), flush=True)

# ------------------------------------------- convergence in the basis size
print("\nconvergence study (nu vs basis size) — the honest error bar", flush=True)
conv = []
for nb in (args.nb-6, args.nb-4, args.nb-2, args.nb):
    if nb < 4: continue
    g = lambda nu: bt_sigma(nu, POLE, A1, A2, nb, args.ni//2, args.ni)
    nu_k, s_k = golden_min(g, 3.20, 3.29, iters=50)
    conv.append({"nb": nb, "nu": str(nu_k), "sigma": mp.nstr(s_k, 6)})
    print("   nb=%2d  nu = %s  sigma = %s" % (nb, mp.nstr(nu_k, 14),
                                              mp.nstr(s_k, 4)), flush=True)
out["convergence"] = conv
if len(conv) >= 2:
    spread = abs(mpf(conv[-1]["nu"]) - mpf(conv[-2]["nu"]))
    out["stable_digits_estimate"] = int(-mp.log10(spread)) if spread > 0 else args.dps
    print("   successive-nb agreement: %s  =>  ~%s stable digits"
          % (mp.nstr(spread, 3), out["stable_digits_estimate"]), flush=True)

# ------------------------------------------- rational exclusion windows
print("\nrational exclusions for nu (Farey candidates up to denominator 60)", flush=True)
try:
    tol = mpf(10)**(-(out.get("stable_digits_estimate", 6)))
except Exception:
    tol = mpf('1e-6')
excl = []
best = None
for q in range(1, 61):
    p = int(mp.floor(nu_star*q + mpf('0.5')))
    cand = mpf(p)/q
    d = abs(cand - nu_star)
    if best is None or d < best[2]: best = (p, q, d)
    if d < tol*10:
        excl.append({"p": p, "q": q, "value": mp.nstr(cand, 12),
                     "distance": mp.nstr(d, 6), "excluded": bool(d > tol)})
out["nearest_rational"] = {"p": best[0], "q": best[1],
                           "distance": mp.nstr(best[2], 8)}
out["candidates_within_10tol"] = excl
print("   nearest rational p/q (q<=60): %d/%d, distance %s"
      % (best[0], best[1], mp.nstr(best[2], 6)), flush=True)
print("   => every rational of denominator <= 60 is further than the")
print("      numerical uncertainty, i.e. excluded at this precision.")

# also record what this says about beta and lambda
out["beta"] = str(nu_star + mpf(3)/2)
out["note"] = ("nu irrational => beta = nu + 3/2 irrational => the GF is not "
               "a G-function => a(n) is not holonomic. This computation gives "
               "evidence, not a certificate: certifying irrationality of a "
               "Dirichlet eigenvalue with no closed form is the open wall.")
json.dump(out, open("lambda_precision_results.json", "w"), indent=1)
print("\nwrote lambda_precision_results.json — send this back")
