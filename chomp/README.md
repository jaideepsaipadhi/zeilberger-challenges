# The 10x42 Chomp bar has three winning moves

Answering the **Second Computational Challenge** of S. B. Ekhad and D. Zeilberger
(chompc, Aug 2018): find a and b such that the a x b Chomp bar has at least THREE
winning moves. Ken Thompson's 8x10 (two winning moves) had stood as the record since
the 1970s; the 2018 search to 14x14 found further doubles but no triple.

## The answer

**The 10 x 42 bar has exactly three winning bites: (5,36), (7,30), (8,26).**

They form a strict staircase (rows increasing, columns decreasing) -- which is forced:

**Staircase Theorem.** The winning bites of any bar form a strict antichain: if (i,j) and
(i',j') are both winning with i < i', then j > j'. *Proof:* the (i,j)-child is one move from
the (i',j')-child whenever i < i' and j <= j' (bite row i, column j); two P-positions one
move apart is impossible. Corollary: no winning bite lies in row 1 (row-1 bites leave full
rectangles, which are N-positions by Gale's strategy-stealing).

## Verification

1. `brute.py`/`cmp.py`: engines validated against brute force on 6 boxes (all bar counts + P-census).
2. History reproduced exactly: Thompson's 8x10 = 2 moves; Ekhad-Zeilberger's 6x13, 10x14 = 2.
3. The triple found by solving the 10x45 box (chomp3, area-layered bit-packed engine), then
   re-derived from the DIFFERENT 10x42 box (triple_verify.log), then counter-signed by the
   independent v2 engine (different iteration order, no early exit; triple_v2.log).

## Extended census (all first computations beyond the 2018 table)

Exactly eighteen bars with exactly two winning moves in the searched region (a<=b):
6x13, 6x93, 7x29, 7x30, 7x57, 8x10, 8x22, 8x23, 9x10, 9x26, 10x14, 10x29, 10x33, 10x35,
11x18, 12x13, 13x16, 14x16 -- plus transposes. Squares through 20 have a unique winning
move. Notable: doubles cluster (7x29/30, 8x22/23, 13x16/14x16 adjacent pairs), and the
multiplicity tracks of fixed-split families follow quasi-periodic laws with rotation
numbers in Q(sqrt(2)) for the 3-row lockstep families (see fit_tracks.py) -- a short note
is in preparation.

## Reproduce

    gcc -O3 -fopenmp -o chomp3 chomp3.c
    ./chomp3 10 42 10 42        # solves the box, prints the bar counts and the three bites

Engines: chomp.c (v1 reference), chomp2.c (v2 independent), chomp3.c (production).
Logs: chomp_results.log (the campaign), triple_verify.log, triple_v2.log (verification).
Track dumps: tb_*.out; analysis: fit_tracks.py.
