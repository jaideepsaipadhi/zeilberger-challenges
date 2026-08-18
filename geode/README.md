# Diagonal Geode numbers: G(1000^D) for D = 4,...,10
Target: G(1000,1000,1000,1000,1000). Rubine (arXiv 2512.21785) claimed the 4D case and noted the
5D case needs "new innovations". DERIVATION.md contains the new idea: a 1-D integral representation
reducing the computation from O(M^3) time / O(M^2) memory to O(M^2) time / O(M) memory.


## Results (all values in this repository, each verified against fresh primes)

| D | digits | file |
|---|--------|------|
| 4 | 6036 | G_D4_M1000.txt (reproduces Rubine, arXiv 2512.21785) |
| 5 | 8367 | G_D5_M1000.txt (the challenge) |
| 6 | 10866 | G_D6_M1000.txt |
| 7 | 13506 | G_D7_M1000.txt |
| 8 | 16267 | G_D8_M1000.txt |
| 9 | 19135 | G_D9_M1000.txt |
| 10 | 22097 | G_D10_M1000.txt |

The full derivation is in note.pdf and DERIVATION.md. Any value can be re-verified modulo
primes of your choosing in seconds: python3 -c "from geode_batch import G_modp_batch; print(G_modp_batch(1000, 5, [1073741789]))"
then compare with int(open('G_D5_M1000.txt').read()) % 1073741789.

## Run commands (adjust --workers to your real core count)
    python3 validate_all.py                                            # must print ALL PASS
    python3 run_geode_cpu.py --D 4 --M 1000 --workers 256 --batch 4    # acceptance, ~30-40s wall
    python3 run_geode_cpu.py --D 5 --M 1000 --workers 256 --batch 4    # the challenge, ~60s wall
Then verify the acceptance run:
    - G_D4_M1000.txt digit count should be ~6040-6110 (rigorous bound: 6264)
    - first/last 20 digits must equal Rubine's printout: 14060489925985310384 ... 41512407713629184000
      (rubine_G4_partial.txt has his printed value's start and end; the middle of that file is an
       incomplete transcription - do not diff the whole file. NOTE: his prose says "6303 digits",
       which is inconsistent with the growth rate of his own smaller values; trust the digits.)
    - independent residues: pick 2 primes NOT in state_D4_M1000.json, run
      python3 -c "from geode_batch import G_modp_batch; print(G_modp_batch(1000,4,[<p1>,<p2>]))"
      and check int(open('G_D4_M1000.txt').read()) % p matches each.
H(1)..H(6) already verified equal to Rubine's printed small values (object identity confirmed).

## Order of operations (single machine / reference)
1. `python3 validate_all.py`                             # must print ALL PASS
2. `python3 run_geode_cpu.py --D 4 --M 1000 --workers N` # ACCEPTANCE TEST vs Rubine's number
3. `python3 run_geode_cpu.py --D 5 --M 1000 --workers N` # the challenge
Resumable (state_D*_M*.json). Requirements: python3, numpy, sympy.
MEASURED at M=1000, D=5 (batched backend geode_batch.py, primes ~2^30): 2.64 s/prime.
927+5 primes => ~41 CPU-core-minutes total; divide by physical cores for wall clock
(prime batches are fully independent). CRT reconstruction adds seconds.
Old single-prime backend (geode_fast.py) is kept as the cross-check reference (~12 s/prime).

## GPU
See CUDA_DESIGN.md for the parallelisation analysis (short version: parallelise over
(prime x node), NOT over the recurrence, which is sequential and too small per step).
    python3 geode_cuda.py --selftest        # REQUIRED: GPU vs validated CPU, both D, several M
    python3 geode_cuda.py --D 4 --M 1000    # acceptance test
    python3 geode_cuda.py --D 5 --M 1000    # the challenge
geode_cuda.py is a cupy RawKernel (Montgomery arithmetic, 31-bit primes, ~4.6M threads).
Its per-thread control flow has been emulated in Python and matches the validated CPU code
exactly, but the CUDA itself has NOT been run on a GPU here - the selftest is the gate.

## Deployment notes (learned the hard way)
- The driver sets OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1 itself; at high --workers this is
  essential (an OpenBLAS thread explosion at 256 workers corrupted numpy imports on RunPod).
- Cloud "core counts" lie: check the cgroup quota (cat /sys/fs/cgroup/cpu.max; cores = quota/period)
  and set --workers to the REAL core count. Oversubscription is strictly slower.
- Python 3.11+ caps int->str at 4300 digits; the driver lifts it (sys.set_int_max_str_digits(0)).
- Runs are resumable; run long ladders under tmux/nohup. Progress ETAs stabilize after ~2 waves.

## Files
DERIVATION.md   the mathematics, with every verified step and the known pitfalls
geode_fast.py   the algorithm (validated, vectorized, D=4 and D=5)
validate_all.py the validation suite
run_geode.py    multi-prime CRT driver (resumable)
fast.py         independent scalar reference implementation (slower, for cross-checking)
derive.py verify.py explore2.py proto.py   the derivation/verification trail
