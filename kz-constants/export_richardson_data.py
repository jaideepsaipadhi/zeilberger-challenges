"""Export the exact sequences and Richardson tables behind C1 and C2
(Kauers-Zeilberger Conjectures 2a/2b).

Outputs (in --out, default ./richardson_data):
  G_exact.csv          n, G(n)            exact integers, runs of length >= 2
  H_exact.csv          n, H(n)            exact integers, all runs odd
  G_richardson.csv     N, order_0..order_K   A_N = G(N) N^4 / 8^N and its Richardson columns
  H_richardson.csv     N, order_0..order_K   A_N = H(N) N^4 / (7+5*sqrt2)^N, same scheme
  summary.json         parameters, validation results, last entry of each column

Richardson step (identical to c1_recompute.py):
  A^{(k+1)}_N = ((N+1) A^{(k)}_{N+1} - (N-k) A^{(k)}_N) / (k+1)
which cancels one further term of the 1/N expansion at each order.

Usage:
  python3 export_richardson_data.py --gmax 500 --hmax 200 --orders 28 --dps 120
Quick test:
  python3 export_richardson_data.py --gmax 60 --hmax 40 --orders 8
"""
import argparse, csv, json, os, time
from mpmath import mp, mpf, sqrt, nstr

ap = argparse.ArgumentParser()
ap.add_argument('--gmax', type=int, default=500)
ap.add_argument('--hmax', type=int, default=200)
ap.add_argument('--orders', type=int, default=28)
ap.add_argument('--dps', type=int, default=120)
ap.add_argument('--out', default='richardson_data')
args = ap.parse_args()
mp.dps = args.dps
os.makedirs(args.out, exist_ok=True)
T0 = time.time()

G_PUBLISHED = [1, 0, 1, 1, 5, 15, 69, 304, 1518, 7807, 42314, 236621,
               1364570, 8062975, 48680547, 299388670, 1871463427]

# ------------------------------------------------------------------ G(n)
def G_series(nmax):
    """Ballot words in 3 letters, every run of every letter of length >= 2.
    State (c1, c2, last letter, r) with r=1 if the current run has length 1."""
    out = [0]*(nmax+1); out[0] = 1
    dp = {(0, 0, 0, 1): 1}
    s = 1
    while s < 3*nmax:
        ndp = {}
        for (c1, c2, d, r), v in dp.items():
            c0 = s - c1 - c2
            for e in (0, 1, 2):
                if e != d and r == 1:
                    continue
                n0, n1, n2 = c0 + (e == 0), c1 + (e == 1), c2 + (e == 2)
                if n0 > nmax or not (n0 >= n1 >= n2):
                    continue
                key = (n1, n2, e, 2 if e == d else 1)
                ndp[key] = ndp.get(key, 0) + v
        dp = ndp; s += 1
        if s % 3 == 0:
            n = s // 3
            out[n] = sum(v for (c1, c2, d, r), v in dp.items()
                         if c1 == n and c2 == n and r == 2)
            if n % 50 == 0:
                print("  G(%d): %d digits  [%.0fs]" % (n, len(str(out[n])), time.time()-T0), flush=True)
    return out

# ------------------------------------------------------------------ H(n)
V = [(1, 0), (-1, 1), (0, -1)]

def H_brute(n):
    N = 3*n; tot = 0
    def rec(x, y, d, par, k):
        nonlocal tot
        if k == N:
            tot += (x == 0 and y == 0 and par == 1); return
        nx, ny = x+V[d][0], y+V[d][1]
        if nx >= 0 and ny >= 0: rec(nx, ny, d, 1-par, k+1)
        if par == 1:
            for e in range(3):
                if e != d:
                    nx, ny = x+V[e][0], y+V[e][1]
                    if nx >= 0 and ny >= 0: rec(nx, ny, e, 1, k+1)
    for d in range(3):
        x, y = V[d]
        if x >= 0 and y >= 0: rec(x, y, d, 1, 1)
    return tot

def H_series(nmax):
    """Quarter-plane excursions with steps (1,0),(-1,1),(0,-1), every run odd.
    Exact-integer port of verification/h_dp.py.  State (x, y, dir, parity),
    parity 1 = current run length odd (run may end)."""
    N = 3*nmax
    out = [0]*(nmax+1); out[0] = 1
    dp = {}
    for d in range(3):
        x, y = V[d]
        if x >= 0 and y >= 0: dp[(x, y, d, 1)] = 1
    for k in range(1, N):
        rem = N - k - 1
        ndp = {}
        for (x, y, d, par), v in dp.items():
            moves = [(d, 1-par)]
            if par == 1:
                moves += [(e, 1) for e in range(3) if e != d]
            for e, np_ in moves:
                nx, ny = x+V[e][0], y+V[e][1]
                if nx < 0 or ny < 0 or nx+ny > rem:   # x+y drops by at most 1 per step
                    continue
                key = (nx, ny, e, np_)
                ndp[key] = ndp.get(key, 0) + v
        dp = ndp
        if (k+1) % 3 == 0:
            n = (k+1)//3
            out[n] = sum(v for (x, y, d, par), v in dp.items() if x == 0 and y == 0 and par == 1)
            if n % 50 == 0:
                print("  H(%d): %d digits  [%.0fs]" % (n, len(str(out[n])), time.time()-T0), flush=True)
    return out

# ------------------------------------------------------------------ Richardson
def richardson(seq, mu, orders):
    NS = list(range(len(seq)//2, len(seq)))
    cols = [[mpf(seq[n]) * mpf(n)**4 / mu**n for n in NS]]
    for k in range(orders):
        prev = cols[-1]
        cols.append([((NS[i]+1)*prev[i+1] - (NS[i]-k)*prev[i]) / (k+1)
                     for i in range(len(prev)-1)])
    return NS, cols

def write_table(path, NS, cols):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['N'] + ['order_%d' % k for k in range(len(cols))])
        for i, N in enumerate(NS):
            w.writerow([N] + [nstr(c[i], mp.dps - 10) if i < len(c) else '' for c in cols])

summary = {"gmax": args.gmax, "hmax": args.hmax, "orders": args.orders, "dps": args.dps}

print("G(n) to n = %d ..." % args.gmax, flush=True)
G = G_series(args.gmax)
assert G[:len(G_PUBLISHED)] == G_PUBLISHED[:args.gmax+1], "G does not match published terms"
summary["G_validation"] = "all %d published terms reproduced" % min(len(G_PUBLISHED), args.gmax+1)

print("H(n) to n = %d ..." % args.hmax, flush=True)
H = H_series(args.hmax)
for n in range(1, 5):
    assert H[n] == H_brute(n), "H(%d) disagrees with brute force" % n
summary["H_validation"] = "n=1..4 match independent brute-force enumeration"
summary["H_first"] = [str(h) for h in H[:8]]

with open(os.path.join(args.out, 'G_exact.csv'), 'w') as f:
    f.write('n,G(n)\n' + ''.join('%d,%d\n' % (n, g) for n, g in enumerate(G)))
with open(os.path.join(args.out, 'H_exact.csv'), 'w') as f:
    f.write('n,H(n)\n' + ''.join('%d,%d\n' % (n, h) for n, h in enumerate(H)))

for name, seq, mu in (('G', G, mpf(8)), ('H', H, 7 + 5*sqrt(2))):
    NS, cols = richardson(seq, mu, args.orders)
    write_table(os.path.join(args.out, '%s_richardson.csv' % name), NS, cols)
    summary['%s_last_by_order' % name] = [nstr(c[-1], 25) for c in cols if c]
    print("%s: deepest Richardson value %s" % (name, nstr(cols[-1][-1] if cols[-1] else cols[0][-1], 20)))

json.dump(summary, open(os.path.join(args.out, 'summary.json'), 'w'), indent=1)
print("wrote %s/  [%.0fs]" % (args.out, time.time()-T0))
