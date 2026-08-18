#!/usr/bin/env python3
"""
final_step7.py — the ansatz reshaped from the factorisation of Q.

WHAT factor_Q.py FOUND (exact):
    Q = Q1^2 * Q2^2      with     gcd(P, Q) = 1
    Q1 = t*w*x + t*x*z - t*x + 1                       (degree 1 in x)
    Q2 = a degree-2-in-x irreducible                    (degree 2 in x)
  and NEITHER factor vanishes at z = 0 or w = 0, so the boundary-kill
  argument N_z(0) = N_w(0) = 0 was valid all along — that failure mode is
  ruled out.

WHY EVERY PREVIOUS SEARCH WAS SHAPED WRONG. All of them used
    Theta = x^a t^b Q^c   =   x^a t^b Q1^{2c} Q2^{2c}
which locks the two exponents together AND forces them EVEN. The natural
certificate denominators — Q1^3 Q2^2, Q1^2 Q2, Q1 Q2^2, Q1^4 Q2^2, ... —
were never in the search space at all. Not "too small": inexpressible. That
is a far better explanation of the persistent tele = 0 than box size.

THIS VERSION:
    Theta = x^a * t^b * Q1^{c1} * Q2^{c2}     c1, c2 INDEPENDENT
keeping the per-direction monomial caps from v6 (which demonstrably helped:
comparable reach fell from U=1547 to U=467) and the blocked elimination with
its rank self-check.

Also: P carries a factor x^2, so the numerator vanishes to second order at
x = 0; the ladder therefore includes small x-exponents a, which the old
Theta could not sensibly reach either.

Mathematics otherwise unchanged: four directions (residues x,t with the
-m/v twist; umbral z,w with the boundary kills), telescoper coefficients in
Q[m] only, self-generated e_X gate (asserts 0,-4,-48), nonzero-telescoper
diagnostic, per-shape checkpointing.

Run:  python3 final_step7.py
      python3 final_step7.py --umax 12000
Writes final_step7_results.json after every shape — SEND IT BACK.
"""
import numpy as np, sympy as sp, json, time, random, argparse, itertools, os, pickle

ap = argparse.ArgumentParser()
ap.add_argument('--umax', type=int, default=6000)
ap.add_argument('--block', type=int, default=256)
ap.add_argument('--nocheck', action='store_true')
args = ap.parse_args()

P = 46337
assert P < 2**16
x, t, z, w, m, v, Y, r1 = sp.symbols('x t z w m v Y r1')
T0 = time.time()
results = {"prime": P, "umax": args.umax}
inv = lambda a: pow(int(a) % P, P-2, P)

# ---------------------------------------------------------------- e_X gate
print("generating e_X data (gate) ...", flush=True)
d1 = sp.expand(1 + t - z*t - w); d2 = sp.expand(t*(1 - z - w))
Dv = 1 + d1*v + d2*v**2
F1 = (1 + v)*(1 + v*t)/Dv
F2 = -t*v**2*(1 + v)*(1 + v*t)/Dv**2
MD = 12
s1 = sp.expand(sp.series(F1, v, 0, MD+1).removeO())
s2 = sp.expand(sp.series(F2, v, 0, MD+1).removeO())
Aa = [sp.expand(s1.coeff(v, k)) for k in range(MD+1)]
Bb = [sp.expand(s2.coeff(v, k)) for k in range(MD+1)]
def umbral(e):
    p = sp.Poly(sp.expand(e), z, w)
    return sum(c*sp.factorial(kk)*sp.factorial(jj)
               for (kk, jj), c in zip(p.monoms(), p.coeffs()))
eX = {mm: sp.Integer(umbral(sp.expand(2*Aa[mm]*Bb[mm]).coeff(t, mm)))
      for mm in range(1, MD+1)}
print("   e_X(1..6) =", [int(eX[k]) for k in range(1, 7)], flush=True)
assert [int(eX[k]) for k in (1, 2, 3)] == [0, -4, -48], "e_X gate broken"
eXm = {k: int(vv) % P for k, vv in eX.items()}

# ---------------------------------------------------------------- H, factors
if not os.path.exists("H_cache.pkl"):
    raise SystemExit("H_cache.pkl missing — run an earlier final_step once")
Ps, Qs, ORDER = pickle.load(open("H_cache.pkl", "rb"))
Px, Qx = sp.sympify(Ps), sp.sympify(Qs)
fac = sp.factor_list(sp.expand(Qx))[1]
fac_sorted = sorted(fac, key=lambda fe: sp.Poly(fe[0], x).degree())
Q1 = sp.expand(fac_sorted[0][0])
Q2 = sp.expand(fac_sorted[1][0])
print("Q1 (deg_x %d):" % sp.Poly(Q1, x).degree(), sp.sstr(Q1), flush=True)
print("Q2 (deg_x %d): [%d terms]" % (sp.Poly(Q2, x).degree(),
                                     len(sp.Add.make_args(Q2))), flush=True)
assert sp.simplify(sp.expand(Q1**2 * Q2**2 - sp.expand(Qx))) == 0, \
    "Q != Q1^2 Q2^2 — factorisation assumption wrong, stop"
print("verified Q = Q1^2 * Q2^2 (%.0fs)" % (time.time()-T0), flush=True)
results["Q1"] = sp.sstr(Q1)

BASE = (x, t, z, w)
f_P  = sp.lambdify(BASE, Px, modules='math')
f_Q1 = sp.lambdify(BASE, Q1, modules='math')
f_Q2 = sp.lambdify(BASE, Q2, modules='math')
d_Q1 = {s: sp.lambdify(BASE, sp.diff(Q1, s), modules='math') for s in BASE}
d_Q2 = {s: sp.lambdify(BASE, sp.diff(Q2, s), modules='math') for s in BASE}
print("evaluators compiled (%.0fs)" % (time.time()-T0), flush=True)

def layout(r, dm, boxes, dmc):
    L = []
    for tag in ('x', 't', 'z', 'w'):
        cx, ct, cz, cw = boxes[tag]
        for e1_, e2_, e3_, e4_ in itertools.product(
                range(cx+1), range(ct+1), range(cz+1), range(cw+1)):
            if tag == 'z' and e3_ == 0: continue
            if tag == 'w' and e4_ == 0: continue
            for l in range(dmc+1):
                L.append(('c', tag, (e1_, e2_, e3_, e4_), l))
    ncert = len(L)
    for i in range(r+1):
        for l in range(dm+1):
            L.append(('p', i, l))
    return L, ncert

def row_at(pt, L, a, b, c1, c2, out):
    xv, tv, zv, wv, mv = pt
    q1 = f_Q1(xv, tv, zv, wv) % P
    q2 = f_Q2(xv, tv, zv, wv) % P
    if q1 == 0 or q2 == 0: raise ZeroDivisionError
    iq1, iq2 = inv(q1), inv(q2)
    Qv = pow(q1, 2, P)*pow(q2, 2, P) % P
    Hv = f_P(xv, tv, zv, wv) % P * inv(Qv) % P
    ix, it_, iz, iw = inv(xv), inv(tv), inv(zv), inv(wv)
    # dlog Theta / dv  for Theta = x^a t^b Q1^c1 Q2^c2
    dlt = {}
    for s, sym, base_exp, ivv in (('x', x, a, ix), ('t', t, b, it_),
                                  ('z', z, 0, iz), ('w', w, 0, iw)):
        val = (base_exp*ivv
               + c1*(d_Q1[sym](xv, tv, zv, wv) % P)*iq1
               + c2*(d_Q2[sym](xv, tv, zv, wv) % P)*iq2) % P
        dlt[s] = val
    Th = pow(xv, a, P)*pow(tv, b, P) % P * pow(q1, c1, P) % P * pow(q2, c2, P) % P
    if Th == 0: raise ZeroDivisionError
    iTh = inv(Th); ixt = ix*it_ % P
    lam = {'x': (-ix) % P, 't': (-it_) % P}
    dv = {'x': ix, 't': it_, 'z': iz, 'w': iw}
    mp = [pow(mv, l, P) for l in range(48)]
    for idx, lab in enumerate(L):
        if lab[0] == 'p':
            _, i, l = lab
            out[idx] = mp[l] * pow(ixt, i+1, P) % P * Hv % P
        else:
            _, tag, (e1_, e2_, e3_, e4_), l = lab
            mu = pow(xv, e1_, P)*pow(tv, e2_, P) % P \
                 * pow(zv, e3_, P) % P * pow(wv, e4_, P) % P
            ev = {'x': e1_, 't': e2_, 'z': e3_, 'w': e4_}[tag]
            inner = (ev * mu % P * dv[tag] - mu*dlt[tag]) % P
            if tag in ('x', 't'):
                inner = (inner + mv*lam[tag] % P * mu) % P
            else:
                inner = (inner - mu) % P
            out[idx] = (-(mp[l] * iTh % P * inner)) % P

def rref_slow(M):
    A = M.copy() % P; R_, Cn = A.shape
    piv = [-1]*Cn; rr = 0
    for cc in range(Cn):
        nz = np.nonzero(A[rr:, cc])[0]
        if len(nz) == 0: continue
        i0 = rr + nz[0]; A[[rr, i0]] = A[[i0, rr]]
        A[rr] = (A[rr]*inv(A[rr, cc])) % P
        colv = A[:, cc].copy(); colv[rr] = 0
        A = (A - np.outer(colv, A[rr])) % P
        piv[cc] = rr; rr += 1
        if rr == R_: break
    return A, piv, rr

def rref_blocked(M, B):
    A = M.copy() % P; R_, Cn = A.shape
    piv = [-1]*Cn; rr = 0; cc = 0
    while cc < Cn and rr < R_:
        hi = min(cc + B, Cn)
        for c2_ in range(cc, hi):
            if rr >= R_: break
            nz = np.nonzero(A[rr:, c2_])[0]
            if len(nz) == 0: continue
            i0 = rr + nz[0]
            if i0 != rr: A[[rr, i0]] = A[[i0, rr]]
            A[rr, cc:] = (A[rr, cc:]*inv(A[rr, c2_])) % P
            colv = A[:, c2_].copy(); colv[rr] = 0
            nzr = np.nonzero(colv)[0]
            if len(nzr):
                A[nzr, cc:] = (A[nzr, cc:] - np.outer(colv[nzr], A[rr, cc:])) % P
            piv[c2_] = rr; rr += 1
        cc = hi
    return A, piv, rr

FIRST = [True]

def run_shape(r, dm, a, b, c1, c2, hi, lo, dmc):
    boxes = {}
    for tag in ('x', 't', 'z', 'w'):
        cx = hi if tag == 'x' else lo
        ct = hi if tag == 't' else lo
        cz = hi if tag == 'z' else lo
        cw = hi if tag == 'w' else lo
        boxes[tag] = (cx, ct, cz, cw)
    L, ncert = layout(r, dm, boxes, dmc); U = len(L)
    if U > args.umax: return None, U
    t1 = time.time(); rows = U + 15
    rnd = random.Random(1)
    M = np.zeros((rows, U), dtype=np.int64); buf = [0]*U
    f = 0; guard = 0
    while f < rows and guard < 80*rows:
        guard += 1
        pt = tuple(rnd.randrange(2, P-2) for _ in range(5))
        try: row_at(pt, L, a, b, c1, c2, buf)
        except ZeroDivisionError: continue
        M[f] = buf; f += 1
    if f < rows: return {"U": U, "error": "insufficient points"}, U
    tasm = time.time() - t1
    t2 = time.time(); A, piv, rank = rref_blocked(M, args.block)
    telim = time.time() - t2
    if FIRST[0] and not args.nocheck and U <= 1500:
        _, _, rank2 = rref_slow(M)
        print("   SELF-CHECK: blocked %d vs reference %d -> %s"
              % (rank, rank2, "AGREE" if rank == rank2 else "DISAGREE"), flush=True)
        assert rank == rank2
        FIRST[0] = False
    free = [cc for cc in range(U) if piv[cc] < 0]
    sols = []
    for fc in free[:12]:
        vv = np.zeros(U, dtype=np.int64); vv[fc] = 1
        for cc in range(U):
            if piv[cc] >= 0: vv[cc] = (-A[piv[cc], fc]) % P
        sols.append(vv)
    hit = None; ntel = 0
    for vv in sols:
        ps = [sp.Integer(0)]*(r+1)
        for k, lab in enumerate(L):
            if lab[0] == 'p': ps[lab[1]] += sp.Integer(int(vv[k])) * m**lab[2]
        ps = [sp.expand(p_) for p_ in ps]
        if all(p_ == 0 for p_ in ps): continue
        ntel += 1
        if all(sum(int(ps[i].subs(m, mm)) % P * eXm[mm+i] for i in range(r+1)) % P == 0
               for mm in range(1, MD-r)):
            hit = ps; break
    rec = {"U": U, "n_cert": ncert, "r": r, "dm": dm, "a": a, "b": b,
           "c1": c1, "c2": c2, "hi": hi, "lo": lo, "dmc": dmc,
           "rank": int(rank), "nullity": U - int(rank),
           "nonzero_telescoper": ntel, "assembly_s": round(tasm, 1),
           "elim_s": round(telim, 1), "hit": bool(hit)}
    if hit: rec["telescoper"] = [sp.sstr(sp.factor(p_)) for p_ in hit]
    return rec, U

# ---- ladder over the NEW degrees of freedom (c1, c2) — odd and mixed
# exponents that x^a t^b Q^c could never express
SHAPES = [   # (r, dm, a, b, c1, c2, hi, lo, dmc)
    (6, 4, 1, 1, 1, 1, 4, 1, 2),
    (6, 4, 1, 1, 2, 1, 4, 1, 2),
    (6, 4, 1, 1, 1, 2, 4, 1, 2),
    (6, 4, 1, 1, 3, 2, 4, 1, 2),
    (6, 4, 1, 1, 2, 3, 4, 1, 2),
    (6, 5, 2, 1, 3, 2, 5, 1, 2),
    (6, 5, 2, 2, 3, 3, 5, 1, 2),
    (6, 5, 2, 2, 4, 2, 5, 1, 2),
    (8, 5, 2, 2, 3, 2, 6, 2, 2),
    (8, 5, 2, 2, 4, 3, 6, 2, 2),
    (8, 6, 2, 2, 4, 3, 6, 2, 3),
]
records, found = [], None
for sh in SHAPES:
    rec, U = run_shape(*sh)
    if rec is None:
        print("stopping: shape needs U=%d > umax=%d" % (U, args.umax)); break
    records.append(rec); results["records"] = records
    json.dump(results, open("final_step7_results.json", "w"), indent=1)
    if "error" in rec:
        print("   U=%-5d %s" % (U, rec["error"]), flush=True); continue
    print("   U=%-5d Th=x^%d t^%d Q1^%d Q2^%d r=%d dm=%d hi=%d lo=%d dmc=%d: "
          "null=%d tele=%d %s [asm %.0fs elim %.0fs]"
          % (rec["U"], rec["a"], rec["b"], rec["c1"], rec["c2"], rec["r"],
             rec["dm"], rec["hi"], rec["lo"], rec["dmc"], rec["nullity"],
             rec["nonzero_telescoper"], "<<< HIT" if rec["hit"] else "",
             rec["assembly_s"], rec["elim_s"]), flush=True)
    if rec["hit"]:
        found = rec; break

results["found"] = found
print("\nFINAL STEP:", "HIT" if found else "no hit in this ladder")
if not found:
    print("nonzero-telescoper vectors across all shapes:",
          sum(r_.get("nonzero_telescoper", 0) for r_ in records))
json.dump(results, open("final_step7_results.json", "w"), indent=1)
print("wrote final_step7_results.json — send this back")
