Richardson data for Kauers-Zeilberger Conjectures 2a/2b
Jaideep Sai Padhi, arXiv:2608.11290 (Section on restricted runs)

G_exact.csv        G(n), n = 0..500, exact (every run length >= 2). Reproduces all 17 published terms.
H_exact.csv        H(n), n = 0..200, exact (every run length odd). Checked against brute force for n <= 4.
G_richardson.csv   A_N = G(N) N^4 / 8^N for N = 250..500, and Richardson columns of order 0..28.
H_richardson.csv   A_N = H(N) N^4 / (7+5*sqrt2)^N for N = 100..200, orders 0..28.
summary.json       parameters and the last entry of every column.
export_richardson_data.py   regenerates everything:  python3 export_richardson_data.py --gmax 500 --hmax 200 --orders 28 --dps 120

Richardson step: A^(k+1)_N = ((N+1) A^(k)_(N+1) - (N-k) A^(k)_N) / (k+1).
Deepest values: C1 ~ 0.52128605909, C2 ~ 0.63892781296.
Note: the higher-order columns for H use only 101 points, so treat the deepest
H orders with care; compare adjacent orders for stability.
