#!/usr/bin/env python3
"""
known_p.py — stop hunting for the telescoper; we already know it.

THE OVERSIGHT THIS FIXES. hunt_ex.py fitted the eX operator from 200 terms:
ORDER 6, DEGREE 16, coefficients mod 46337 (in ex_results.json). Every search
today hunted for p_i(m) from scratch with dm = 4..7. A degree-16 telescoper
CANNOT BE REPRESENTED at dm <= 7 — so tele = 0 was guaranteed regardless of
Theta, regardless of U. That likely explains every negative result today.

THE BETTER FORMULATION. With p known, the identity

    sum_i p_i(m) Omega^i R = (d_x + m lam_x)G_x + (d_t + m lam_t)G_t
                             + (d_z - 1)G_z + (d_w - 1)G_w

becomes a LINEAR SOLVE WITH KNOWN RIGHT-HAND SIDE: build b from the known p,
build the certificate columns M_c, and ask whether b lies in the column space
of M_c. Solvable  <=>  rank(M_c) == rank([M_c | b]).

That is smaller (no telescoper block), better conditioned, and gives a
DEFINITIVE yes/no per Theta shape instead of a search. If solvable, the
recovered certificate together with the known p IS the proof object.

INPUT: ex_results.json (from hunt_ex.py) — needs eX.shape = [6,16] and
eX.opvec_modp, the operator coefficients mod 46337 laid out as
opvec[k*(d+1)+j] = coefficient of m^j in p_k(m).

GATE: the loaded p is first checked against the e_X data generated here from
the verified chain — if it does not annihilate it, the script stops.

Run:  python3 known_p.py
      python3 known_p.py --exfile ex_results.json
Writes known_p_results.json — SEND THIS BACK.
"""
import numpy as np, sympy as sp, json, time, random, argparse, itertools, os, pickle

ap = argparse.ArgumentParser()
ap.add_argument('--exfile', default='ex_results.json')
ap.add_argument('--umax', type=int, default=8000)
ap.add_argument('--block', type=int, default=256)
args = ap.parse_args()

P = 46337
x, t, z, w, m, v = sp.symbols('x t z w m v')
T0 = time.time()
results = {"prime": P}
inv = lambda a: pow(int(a) % P, P-2, P)

# ---------------------------------------------------------------- e_X data
print("generating e_X data ...", flush=True)
d1 = sp.expand(1 + t - z*t - w); d2 = sp.expand(t*(1 - z - w))
Dv = 1 + d1*v + d2*v**2
F1 = (1 + v)*(1 + v*t)/Dv
F2 = -t*v**2*(1 + v)*(1 + v*t)/Dv**2
MD = 30
s1 = sp.expand(sp.series(F1, v, 0, 6).removeO())
s2 = sp.expand(sp.series(F2, v, 0, 6).removeO())
A = [sp.expand(s1.coeff(v, k)) for k in range(6)]
B = [sp.expand(s2.coeff(v, k)) for k in range(6)]
cB = [None, sp.expand(2*d1), sp.expand(d1**2 + 2*d2),
      sp.expand(2*d1*d2), sp.expand(d2**2)]
while len(A) <= MD:
    k = len(A); A.append(sp.expand(-(d1*A[k-1] + d2*A[k-2])))
while len(B) <= MD:
    k = len(B)
    B.append(sp.expand(-(cB[1]*B[k-1] + cB[2]*B[k-2] + cB[3]*B[k-3] + cB[4]*B[k-4])))
def umbral(e):
    p = sp.Poly(sp.expand(e), z, w)
    return sum(c*sp.factorial(kk)*sp.factorial(jj)
               for (kk, jj), c in zip(p.monoms(), p.coeffs()))
eX = {}
for mm in range(1, MD+1):
    eX[mm] = int(umbral(sp.expand(2*A[mm]*B[mm]).coeff(t, mm))) % P
print("   e_X(1..5) mod p =", [eX[k] for k in range(1, 6)], flush=True)

# ---------------------------------------------------------------- known p
if not os.path.exists(args.exfile):
    raise SystemExit("%s not found — it holds the fitted eX operator" % args.exfile)
ex = json.load(open(args.exfile))
shape = ex["eX"]["shape"]; opv = ex["eX"]["opvec_modp"]
r_kn, d_kn = int(shape[0]), int(shape[1])
print("known eX operator: order %d, degree %d (%d coefficients)"
      % (r_kn, d_kn, len(opv)), flush=True)
assert len(opv) == (r_kn+1)*(d_kn+1), "opvec length mismatch"
pk = [[int(opv[k*(d_kn+1)+j]) % P for j in range(d_kn+1)] for k in range(r_kn+1)]

def p_at(k, mm):
    s = 0; mp = 1
    for j in range(d_kn+1):
        s = (s + pk[k][j]*mp) % P
        mp = mp*mm % P
    return s

bad = [mm for mm in range(1, MD-r_kn)
       if sum(p_at(k, mm)*eX[mm+k] for k in range(r_kn+1)) % P != 0]
print("   annihilates generated e_X data:", "YES" if not bad else "NO at %s" % bad[:5],
      flush=True)
results["known_p_validates"] = (not bad)
assert not bad, "the fitted operator does not annihilate our e_X — stop"

# ---------------------------------------------------------------- H, factors
if not os.path.exists("H_cache.pkl"):
    raise SystemExit("H_cache.pkl missing")
Ps, Qs, ORDER = pickle.load(open("H_cache.pkl", "rb"))
Px, Qx = sp.sympify(Ps), sp.sympify(Qs)
fac = sorted(sp.factor_list(sp.expand(Qx))[1], key=lambda fe: sp.Poly(fe[0], x).degree())
Q1, Q2 = sp.expand(fac[0][0]), sp.expand(fac[1][0])
BASE = (x, t, z, w)
f_P = sp.lambdify(BASE, Px, modules='math')
f_Q1 = sp.lambdify(BASE, Q1, modules='math')
f_Q2 = sp.lambdify(BASE, Q2, modules='math')
dQ1 = {s: sp.lambdify(BASE, sp.diff(Q1, s), modules='math') for s in BASE}
dQ2 = {s: sp.lambdify(BASE, sp.diff(Q2, s), modules='math') for s in BASE}
print("evaluators ready (%.0fs)" % (time.time()-T0), flush=True)

def cert_layout(boxes, dmc):
    L = []
    for tag in ('x', 't', 'z', 'w'):
        cx, ct, cz, cw = boxes[tag]
        for e1_, e2_, e3_, e4_ in itertools.product(
                range(cx+1), range(ct+1), range(cz+1), range(cw+1)):
            if tag == 'z' and e3_ == 0: continue
            if tag == 'w' and e4_ == 0: continue
            for l in range(dmc+1):
                L.append((tag, (e1_, e2_, e3_, e4_), l))
    return L

def row_and_rhs(pt, L, a, b, c1, c2):
    xv, tv, zv, wv, mv = pt
    q1 = f_Q1(xv, tv, zv, wv) % P; q2 = f_Q2(xv, tv, zv, wv) % P
    if q1 == 0 or q2 == 0: raise ZeroDivisionError
    iq1, iq2 = inv(q1), inv(q2)
    Qv = pow(q1, 2, P)*pow(q2, 2, P) % P
    Hv = f_P(xv, tv, zv, wv) % P * inv(Qv) % P
    ix, it_, iz, iw = inv(xv), inv(tv), inv(zv), inv(wv)
    dlt = {}
    for s, sym, be, ivv in (('x', x, a, ix), ('t', t, b, it_),
                            ('z', z, 0, iz), ('w', w, 0, iw)):
        dlt[s] = (be*ivv + c1*(dQ1[sym](xv, tv, zv, wv) % P)*iq1
                  + c2*(dQ2[sym](xv, tv, zv, wv) % P)*iq2) % P
    Th = pow(xv, a, P)*pow(tv, b, P) % P * pow(q1, c1, P) % P * pow(q2, c2, P) % P
    if Th == 0: raise ZeroDivisionError
    iTh = inv(Th); ixt = ix*it_ % P
    lam = {'x': (-ix) % P, 't': (-it_) % P}
    dv = {'x': ix, 't': it_, 'z': iz, 'w': iw}
    mp = [pow(mv, l, P) for l in range(48)]
    # RHS: sum_i p_i(mv) * Omega^i * R,  R = H/(x t),  Omega = 1/(x t)
    rhs = 0
    for i in range(r_kn+1):
        rhs = (rhs + p_at(i, mv) * pow(ixt, i+1, P) % P * Hv) % P
    row = []
    for tag, (e1_, e2_, e3_, e4_), l in L:
        mu = pow(xv, e1_, P)*pow(tv, e2_, P) % P \
             * pow(zv, e3_, P) % P * pow(wv, e4_, P) % P
        ev = {'x': e1_, 't': e2_, 'z': e3_, 'w': e4_}[tag]
        inner = (ev*mu % P * dv[tag] - mu*dlt[tag]) % P
        if tag in ('x', 't'):
            inner = (inner + mv*lam[tag] % P * mu) % P
        else:
            inner = (inner - mu) % P
        row.append(mp[l] * iTh % P * inner % P)
    return row, rhs

def rank_of(M, B_):
    A_ = M.copy() % P; R_, Cn = A_.shape
    rr = 0
    for cc in range(Cn):
        nz = np.nonzero(A_[rr:, cc])[0]
        if len(nz) == 0: continue
        i0 = rr + nz[0]
        if i0 != rr: A_[[rr, i0]] = A_[[i0, rr]]
        A_[rr] = (A_[rr]*inv(A_[rr, cc])) % P
        colv = A_[:, cc].copy(); colv[rr] = 0
        nzr = np.nonzero(colv)[0]
        if len(nzr):
            A_[nzr] = (A_[nzr] - np.outer(colv[nzr], A_[rr])) % P
        rr += 1
        if rr == R_: break
    return rr

def solvable(a, b, c1, c2, hi, lo, dmc):
    boxes = {}
    for tag in ('x', 't', 'z', 'w'):
        boxes[tag] = (hi if tag == 'x' else lo, hi if tag == 't' else lo,
                      hi if tag == 'z' else lo, hi if tag == 'w' else lo)
    L = cert_layout(boxes, dmc); U = len(L)
    if U > args.umax: return None, U
    t1 = time.time(); rows = U + 20
    rnd = random.Random(3)
    M = np.zeros((rows, U), dtype=np.int64); bvec = np.zeros(rows, dtype=np.int64)
    f = 0; guard = 0
    while f < rows and guard < 80*rows:
        guard += 1
        pt = tuple(rnd.randrange(2, P-2) for _ in range(5))
        try: row, rhs = row_and_rhs(pt, L, a, b, c1, c2)
        except ZeroDivisionError: continue
        M[f] = row; bvec[f] = rhs; f += 1
    if f < rows: return {"U": U, "error": "insufficient points"}, U
    rk = rank_of(M, None)
    aug = np.concatenate([M, bvec.reshape(-1, 1)], axis=1)
    rka = rank_of(aug, None)
    rec = {"U": U, "a": a, "b": b, "c1": c1, "c2": c2, "hi": hi, "lo": lo,
           "dmc": dmc, "rank_Mc": int(rk), "rank_aug": int(rka),
           "SOLVABLE": bool(rk == rka), "secs": round(time.time()-t1, 1)}
    return rec, U

SHAPES = [   # (a, b, c1, c2, hi, lo, dmc)
    (1, 1, 1, 1, 4, 1, 2), (1, 1, 2, 1, 4, 1, 2), (1, 1, 1, 2, 4, 1, 2),
    (1, 1, 3, 2, 4, 1, 2), (2, 2, 3, 2, 5, 1, 2), (2, 2, 3, 3, 5, 1, 2),
    (2, 2, 4, 3, 6, 2, 2), (2, 2, 3, 2, 6, 2, 3), (3, 3, 4, 3, 6, 2, 3),
]
records, hit = [], None
for sh in SHAPES:
    rec, U = solvable(*sh)
    if rec is None:
        print("skipping U=%d > umax" % U); continue
    records.append(rec); results["records"] = records
    json.dump(results, open("known_p_results.json", "w"), indent=1)
    if "error" in rec:
        print("   U=%-5d %s" % (U, rec["error"]), flush=True); continue
    print("   U=%-5d Th=x^%d t^%d Q1^%d Q2^%d hi=%d lo=%d dmc=%d: "
          "rank %d vs aug %d -> %s [%.0fs]"
          % (rec["U"], rec["a"], rec["b"], rec["c1"], rec["c2"], rec["hi"],
             rec["lo"], rec["dmc"], rec["rank_Mc"], rec["rank_aug"],
             "SOLVABLE <<<" if rec["SOLVABLE"] else "not solvable", rec["secs"]),
          flush=True)
    if rec["SOLVABLE"]:
        hit = rec; break

results["solvable_shape"] = hit
print("\nRESULT:", "SOLVABLE — certificate exists in this family" if hit
      else "no Theta in this ladder makes it solvable")
json.dump(results, open("known_p_results.json", "w"), indent=1)
print("wrote known_p_results.json — send this back")
