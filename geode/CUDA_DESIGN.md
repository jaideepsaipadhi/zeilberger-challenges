# GPU design notes

## The trap: do NOT parallelize the inner loop
The per-node work is the recurrence  p_{r+1} = (1/(g_0(r+1))) sum_{s=1}^{s_d} g_s (Ms-r-1+s) p_{r+1-s}.
It is SEQUENTIAL in r (5000 steps at M=1000) and each step touches only s_d=10 values.
A kernel launch per r-step gives ~50k element-ops per launch against ~5 us launch overhead:
the GPU would idle and the port would be SLOWER than numpy.

## The right mapping: one thread per (prime, node)
Both outer dimensions are embarrassingly parallel and there is no communication between them:
  - primes: ~930 independent runs (30-bit primes for the ~8400-digit answer)
  - quadrature nodes: DM+1 = 5001 per prime
=> 930 x 5001 = 4.65M independent threads, each marching its own 5000-step recurrence
   entirely in registers (rolling buffer of s_d=10 values + accumulator).
Total work 930 * 5001 * 5000 * ~10 = 2.3e11 modular multiplications, zero synchronisation.
Per-thread state: ~30 registers. No shared memory needed. This is close to an ideal CUDA workload.

## Arithmetic
Use 31-bit primes and Montgomery multiplication (u32 operands, u64 intermediate).
`%` on the device is an integer division and is slow; REDC is ~4 cheap ops.
  m = (u32)t * pinv ;  u = (t + (u64)m*p) >> 32 ;  if (u>=p) u-=p
with pinv = -p^{-1} mod 2^32. Keep every value in the Montgomery domain for the whole march;
convert once in, once out. 62-bit primes would halve the prime count but need 128-bit
intermediates (__umul64hi + Barrett) - slower per op; 31-bit is the sweet spot.

## Per-prime host precompute (small, O(n) each)
  invr[r] = (r+1)^{-1} mod p     via inv[i] = -(p/i)*inv[p mod i] mod p
  Cb[s]   = C(E, s) mod p        by the ratio recurrence
  w_j     = quadrature weights   O(n^2) per prime - THE REMAINING CPU BOTTLENECK
Weights: nodes are fixed (x_j = j+2), so w_j depends on p only through reduction. Options,
in increasing effort: (a) numpy per prime (~0.05 s => ~1 min total, fine); (b) one CUDA block
per prime building the node polynomial cooperatively; (c) compute the weights once as exact
rationals and reduce mod each prime.

## Expected runtime at M=1000, D=5
Kernel: 2.3e11 modmuls. At a realistic 2e10 Montgomery-modmul/s on one modern card that is
~10-20 s; even ten times pessimistic it is a few minutes. Host precompute ~1 min.
Memory: partials 930*5001*4 B = 19 MB; Cb+invr+w 3*930*5001*4 B = 56 MB. Nothing is large.

## Correctness protocol
`geode_cuda.py --selftest` runs the GPU path and the validated CPU path (geode_fast.py) on the
same small (M, prime) set and compares. Do not trust a full run that has not passed the selftest,
and still run the D=4 acceptance test against Rubine's published number.
