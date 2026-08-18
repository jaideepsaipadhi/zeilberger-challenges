# Holonomicity of the restricted permutation counts $a_{r,s}$ and $b_{r,s}$

An answer to the third challenge of Spahn and Zeilberger,
[arXiv:2211.02550](https://arxiv.org/abs/2211.02550) (also ECA 3 (2023),
S2A10), which asks whether

$$a_{r,s}(n) = \#\{\pi \in S_n : \pi_{i+r} - \pi_i \neq s\},\qquad
  b_{r,s}(n) = \#\{\pi \in S_n : |\pi_{i+r} - \pi_i| \neq s\}$$

are holonomic for all $r,s > 1$, noting that no general theory was known.

**They are, for every fixed $(r,s)$.**

## The idea

Matsuo's bijection, used by Spahn and Zeilberger only for $(r,s) = (2,2)$,
generalises. Relabel *positions* by interleaving the $r$ arithmetic
progressions of step $r$, and *values* by interleaving the $s$ progressions of
step $s$. Within a block, position-distance $r$ becomes $1$ and value-difference
$s$ becomes $1$, so the condition turns into the adjacency condition
$\sigma_{k+1} - \sigma_k \neq 1$ — excused only where the interleaving breaks.

The excused set has size exactly

* $(r-1) + (s-1)$ for $a_{r,s}$,
* $(r-1) + 2(s-1)$ for $b_{r,s}$ (both signs must be excused),

**independently of $n$**. That boundedness is the whole point: it makes the
two-layer inclusion–exclusion representation carry finitely many catalytic
variables for each fixed $(r,s)$, and holonomicity then follows from standard
closure theorems.

## Contents

```
holonomicity-note.tex     the paper
holonomicity-note.pdf     compiled

matsuo_general.py         verifies the bijection for a_{r,s}
matsuo_b.py               verifies it for b_{r,s} (two-sign excusal rule)
gaps_check.py             seam indices are affine in n on residue classes;
                          completes the verification table

results/
  matsuo_general_results.json
  matsuo_b_results.json
  gaps_check_results.json
```

## Reproducing

```bash
python3 matsuo_general.py     # a-case,  n = 4..8, seven (r,s) pairs
python3 matsuo_b.py           # b-case,  three candidate excusal rules
python3 gaps_check.py         # affine seam check + completed table
```

Pure Python, no dependencies, seconds each. All three brute-force the counts
directly from the definitions and compare against the transformed counts — no
result below is asserted without a machine check.

### What the scripts establish

`matsuo_general.py` — for $(r,s) \in \{(2,2),(2,3),(3,2),(3,3),(2,4),(4,2),(3,4)\}$
and $4 \le n \le 8$, the transformed count equals $a_{r,s}(n)$ exactly in all 32
rows, with excusal size exactly $(r-1)+(s-1)$ and constant in $n$.

`matsuo_b.py` — the two-sign case needs the symmetrised rule (excuse an ascent
at a forward value seam *and* a descent at a backward one). That rule
reproduces $b_{r,s}(n)$ exactly in every row; the one-sided rule undercounts and
a more permissive rule overcounts. Excusal size exactly $(r-1)+2(s-1)$.

`gaps_check.py` — the seam indices are exact affine functions of $n$ on each
residue class modulo $\mathrm{lcm}(r,s)$, with slopes $c/r$ and $c/s$ as the
block-length formula $\lfloor (n-c)/k \rfloor + 1$ predicts. This is what
licenses the sectioning step of the proof. Also completes the table with
$(3,4)$, $(4,3)$, $(4,4)$ for both families.

## Scope of the verification

Brute force is $O(n!)$, so the checks run to $n = 8$. They confirm a pattern,
not a theorem — but the theorem is proved in the note (Lemma 2 and
Proposition 4), and the computations are there to catch an error in the
statement, not to substitute for the proof.

## Status

Not yet submitted; not yet externally reviewed.

## Author

Jaideep Sai Padhi — jpadhi@purdue.edu
