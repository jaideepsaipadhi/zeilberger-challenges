Richardson data for Kauers-Zeilberger Conjectures 2a/2b
Jaideep Sai Padhi, arXiv:2608.11290

RESULT (digits common to three independent accelerations)
  C1 = 0.52128605909   (spread of the three estimates 7.8e-12)
  C2 = 0.6389278129    (spread 2.1e-11)

  CSV values are printed to ~25 digits for reproducibility only. Do NOT read
  them as accurate beyond the figures above.

FILES
  G_exact.csv         G(n), n = 0..500, exact (every run length >= 2).
                      Reproduces all 17 published Kauers-Zeilberger terms.
  H_exact.csv         H(n), n = 0..200, exact (every run length odd).
                      Matches independent brute-force enumeration for n <= 4.
  G_richardson.csv    A_N = G(N) N^4 / 8^N, N = 250..500, Richardson orders 0..28.
  H_richardson.csv    A_N = H(N) N^4 / (7+5 sqrt2)^N, N = 100..200, orders 0..28.
  summary.json        parameters and the deepest entry of every order.
  crosscheck.py       Richardson vs Salzer vs geometric tail fit, from the CSVs.
  crosscheck.json     its output.
  ../export_richardson_data.py   regenerates the CSVs (about 3 minutes):
      python3 export_richardson_data.py --gmax 500 --hmax 200 --orders 28 --dps 120

RICHARDSON STEP
  A^(k+1)_N = ((N+1) A^(k)_(N+1) - (N-k) A^(k)_N) / (k+1)

HOW WELL IT CONVERGES (change in deepest value when adding one order)

  order        4         10        20        28
  C1        4.6e-9    1.7e-11   7.6e-13   1.9e-13
  C2        6.7e-8    4.7e-11   2.2e-12   5.7e-13

  Convergence is monotone but slow: increments shrink by only ~15% per order
  at depth, and the column is still rising at order 28. Richardson alone
  under-estimates the limit. This is the failure mode behind our earlier,
  wrong value C1 = 0.5212860373 (order 6, n = 110), which looked stable and
  was wrong in the eighth digit.

THREE INDEPENDENT ESTIMATES (crosscheck.json)

                  Richardson-28        Salzer-16            tail fit
  C1     0.5212860590904      0.5212860590836      0.5212860590914
  C2     0.6389278129613      0.6389278129437      0.6389278129646

  Salzer is also still rising with its order (C1: 12 -> 0.52128605907082,
  16 -> 0.52128605908360), so it approaches from below like Richardson; the
  tail fit (ratio ~0.85) sums the remaining geometric increments. Only the
  digits shared by all three are asserted above.
