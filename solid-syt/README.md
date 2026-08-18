# Proof of Zeilberger's recurrence for solid SYT of shape [[n,n],[n,1]]
(First Rigorous Challenge, solid SYT project)

## Contents
- paper/   : LaTeX source + compiled PDF (9pp, final) + auto-generated appendices
- code/    : all verification and derivation scripts (pure Python 3 + sympy)
- data/    : sequence data and exact symbolic objects

## Key files
- data/sytchal55.txt      : g(n) for n = 1..55 (two independent computations agree)
- data/recurrence.json    : the order-2, degree-12 operator (p0, p1, p2 coefficient lists)
- data/minpoly_A.json     : minimal polynomial of the diagonal series A(t) = Qd(t;t)
- data/ode_A.json         : exact order-3 inhomogeneous ODE for A (derived in Q(t,T,sigma))
- code/derive3.py         : Cramer derivation of the operator from the 2-dim module (matches recurrence.json exactly)
- code/lemmaB_*.py        : Lemma B chain (minpoly -> ODE -> recurrence -> ratio identity)
- code/lemmaA_*.py, code/gap2_close.py : Lemma A chain (pairing, kernel collapse, Lagrange, closing checks)
- code/star_verify.py, code/qd_verify.py : verification of the diagonal equation and the Qd closed form

## Quick verification (polish-pass additions)
- code/proof_extractions.py : machine-checks EVERY step of the joint Wiener-Hopf proof of
  Props Qd and R (the dagger identities, the factorization, the split form, both support
  classifications, and that the >0 / <=0 extractions reproduce Qd and R exactly)
- code/pera_check.py : symbolic proof of the per-a diagonal closed form (bracket coefficient
  identity, Gamma-ratio == 1, numeric check n<=8)

## Quick verification
    python3 verify.py
checks that the operator in data/recurrence.json annihilates all 55 exact terms
(convention: sum_{i=0}^{2} p_i(n) g(n-i) = 0). The deeper machine-verification trail
(functional equation, the (STAR) diagonal equation, the Qd closed form, both Lemma chains)
is reproducible from the scripts in code/ -- see the paper's appendix for the map.

## Status (Aug 2026)

Complete. The paper (paper/paper.pdf) contains full proofs of every step; all identities
are machine-verified by the scripts in code/ (see the two Quick verification sections above).
