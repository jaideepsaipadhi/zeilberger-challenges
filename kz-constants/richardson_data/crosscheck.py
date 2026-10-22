"""Three independent accelerations of A_N for C1 and C2, read from the exact CSVs.

  Richardson (order K)   assumes A_N = C (1 + a1/N + a2/N^2 + ...)
  Salzer     (order k)   weighted differences of N^k A_N; different hypotheses
  Tail fit               fits the deepest Richardson increments d_k ~ a r^k and
                         sums the geometric remainder

Only digits shared by all three should be asserted.  Writes crosscheck.json.
Run from this folder:  python3 crosscheck.py
"""
import csv, json
from mpmath import mp, mpf, sqrt, binomial, nstr
mp.dps = 120

def load(name):
    with open(name) as f:
        r = csv.reader(f); next(r)
        return [int(v) for _, v in r]

def richardson(A, NS, K):
    cols = [A[:]]
    for k in range(K):
        p = cols[-1]
        cols.append([((NS[i]+1)*p[i+1] - (NS[i]-k)*p[i]) / (k+1) for i in range(len(p)-1)])
    return [c[-1] for c in cols]

def salzer(A, NS, k):
    base = len(A) - k - 1
    num = den = mpf(0)
    for j in range(k+1):
        w = (-1)**j * binomial(k, j) * mpf(NS[base+j])**k
        num += w*A[base+j]; den += w
    return num/den

def analyse(seq, mu, K=28, salz_orders=(12, 14, 16), fit=6):
    NS = list(range(len(seq)//2, len(seq)))
    A = [mpf(seq[n]) * mpf(n)**4 / mu**n for n in NS]
    deep = richardson(A, NS, K)
    inc = [deep[i+1]-deep[i] for i in range(len(deep)-1)]
    ratios = [inc[-i]/inc[-i-1] for i in range(1, fit+1)]
    r = sum(ratios)/len(ratios)
    tail = deep[-1] + inc[-1]*r/(1-r)
    sz = {k: salzer(A, NS, k) for k in salz_orders}
    vals = [deep[-1], sz[max(salz_orders)], tail]
    return {
        "richardson_order_%d" % K: nstr(deep[-1], 16),
        "salzer": {str(k): nstr(v, 16) for k, v in sz.items()},
        "tail_fit": nstr(tail, 16),
        "tail_fit_ratio": nstr(r, 6),
        "increments_last_5": [nstr(d, 3) for d in inc[-5:]],
        "spread_of_three": nstr(max(vals)-min(vals), 3),
    }

out = {"C1 (G, mu=8)": analyse(load("G_exact.csv"), mpf(8)),
       "C2 (H, mu=7+5sqrt2)": analyse(load("H_exact.csv"), 7 + 5*sqrt(2))}
json.dump(out, open("crosscheck.json", "w"), indent=1)
print(json.dumps(out, indent=1))
