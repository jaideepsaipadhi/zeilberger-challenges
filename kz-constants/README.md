# A local limit theorem for excursions of Markov-modulated random walks in cones

Proving **Kauers-Zeilberger Conjectures 2a and 2b** (arXiv:2006.10205):
G(n) ~ C1*8^n/n^4 and H(n) ~ C2*(7+5sqrt2)^n/n^4, with C2 = 0.6389278129(4)
computed here for the first time.

The theorem transports the Denisov-Wachtel-Zhang harmonic-function program to
Markov-additive walks via an explicit two-level Poisson corrector ladder
(c(i) = -v_i, q(i) = (8/3)(W/3 - v_i v_i^T)). Two structural identities fall out:
criticality of the tilt is exactly k-independence of path weights (2W(w_c)=1),
and each constant is a universal factor times the *square* of the discrete
cone-harmonic function's apex value.

`paper/` — the manuscript. `verification/` — machine checks of every
model-specific identity: Perron drift/covariance (lemmaA_perron.py), torus
aperiodicity scan (lemmaA_torus.py), the corrector ladder (lemmaB_correctors.py),
Monte Carlo certification of the exit exponent and cone profile (mc_exit.py),
the 2b constants (lemmaA_H.py), and the C2 extraction with calibration against
known C1 (h_dp.py).

Rigor note (also in the paper's introduction): three steps are transparent
citations to Denisov-Zhang rather than rewritten pages — B.1's potential
absorption, stretches of B.3's replay, D3's truncation import. Every
model-specific input to those arguments is verified here.
