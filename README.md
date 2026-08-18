# Zeilberger challenge problems — solutions and partial results

Jaideep Sai Padhi · jpadhi@purdue.edu

Solutions to five challenge problems posed by Doron Zeilberger and his
collaborators, together with substantial partial progress on two more, and the
combined paper describing both the results and the method by which they were
obtained.

This repository consolidates work previously held in five separate
repositories, which it supersedes.

---

## Status

| Item | Prize | Status | Directory |
|---|---|---|---|
| Chomp — 2nd Computational Challenge | $100 | **Solved.** 10×42 has three winning moves | `chomp/` |
| Spahn–Zeilberger Challenge 3 | $300 | **Solved.** All a(r,s), b(r,s) holonomic | `restricted-permutations/` |
| Solid SYT — 1st Rigorous Challenge | — | **Solved.** Order-2 recurrence proved | `solid-syt/` |
| Geode Challenge, D = 5 | $200 | **Solved.** Also D = 6…10 | `geode/` |
| Kauers–Zeilberger 2a, 2b | 2 × $200 | **Solved.** LLT for excursions in cones | `kz-constants/` |
| Spahn–Zeilberger Challenge 1 | $100 | Partial — W-half proved; operator open | `a22-partial/` |
| Solid SYT — 2nd Rigorous Challenge | $10 | Partial — cone exponent bracketed | `cone-exponent/` |

---

## Headline results

**Chomp.** The bar 10×42 has exactly three winning opening bites, at (5,36),
(7,30) and (8,26), answering Ekhad and Zeilberger's Second Computational
Challenge. With it: eighteen previously unrecorded doubles, a Staircase
Theorem (winning bites form a strict antichain), a parity theorem for square
bars, and Beatty-type conjectures with rotation numbers in Q(√2).

**Holonomicity.** For every fixed (r,s), both a(r,s) and b(r,s) are holonomic.
Matsuo's bijection generalises: interleaving arithmetic progressions converts
the (r,s)-condition into the adjacency condition, excused on a set of size
exactly (r−1)+(s−1), resp. (r−1)+2(s−1) — bounded independently of n.

**Solid SYT.** g(n) = c(n)(3n+1)K(n) + T(n), giving the conjectured order-2,
degree-12 recurrence. By-products: an explicit algebraic generating function
for reverse-Kreweras diagonal walks, and a closed form for diagonal-endpoint
counts.

**Geode.** A one-dimensional integral representation reduces the computation
from O(M³) time and O(M²) memory to O(M²) and O(M), making G(1000⁵) — 8367
digits — a matter of minutes.

**Kauers–Zeilberger.** A local limit theorem for excursions of
Markov-modulated random walks in cones proves both 2a and 2b. Each constant is
a universal factor times the *square* of a discrete cone-harmonic function's
apex value. Refined numerically:

    C1 = 0.52128605909(2)      (published estimate: "close to 0.521286")
    C2 = 0.6389278129(4)       (published estimate: "close to 0.63892")

Both determined independently, by exact enumeration plus three mutually
independent accelerations agreeing to ~10⁻¹¹. See `kz-constants/`.

---

## Layout

```
paper/                    the combined paper (source + PDF)
chomp/                    note + thin-tail track computations
restricted-permutations/  note + four verification scripts
solid-syt/                paper + appendices
geode/                    note (integral representation)
kz-constants/             LLT paper + constant determination
a22-partial/              the a(2,2) reduction and its obstruction
cone-exponent/            eigenvalue determination and bounds
```

Engines, logs and bulk data from the original repositories (Chomp C solvers,
Geode CRT driver, solid-SYT derivation trail, KZ verification suite) are
retained in the git history of those repositories.

---

## Method

All of this work was produced in collaboration with an AI system (Claude). The
division of labour is set out in the paper: the system generated the ideas,
proofs and code; the author directed the investigation, ran the computations,
enforced verification, and caught the errors. The paper's Section 3 describes
the working method and Section 12 records the failures, which we regard as part
of the report.

The governing protocol, stated once because it did most of the work: **check
against a value independently known to be correct, as early and as cheaply as
possible.** Internal consistency is not validation.

## Licence

MIT.
