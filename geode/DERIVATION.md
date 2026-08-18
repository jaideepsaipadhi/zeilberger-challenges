# G(M,...,M) for the Geode: a 1-D integral representation

Notation: D = number of variables t_2..t_{D+1}; d = D-1; s_d = d(d+1)/2.
S = 1 + sum_{k>=2} t_k S^k ; S - 1 = (t_2+...+t_{D+1}) G. Target: G[M,...,M].

## 1. Hyper-Catalan closed form (classical; re-derived by Lagrange inversion, verified in code)
    C[m] = K! / ( (1+K-N)! prod_k m_k! ),  K = sum_k k m_k,  N = sum_k m_k.

## 2. Division by (t_2+...+t_{D+1}) (geometric expansion in t_2)
    G[M] = sum_{alpha} (-1)^{|a|} multinomial(|a|;a) C[M + (1+|a|)e_2 - sum_k a_k e_k],
    alpha ranging over {0..M}^d.  Naive size at M=1000, D=5: 1001^4 = 1.0e12 terms.

## 3. Three structural facts (each verified independently in validate_all.py)
 (a) N' = sum_k m'_k = D*M + 1 is CONSTANT in alpha.
     Hence, with W := sum_{k=1}^{d} k*a_{k+2}, K' = base_k - W and
        a_W := K'!/(1+K'-N')! = ff(K', DM) = (DM)! * C(K', DM),
     a POLYNOMIAL in W of degree DM.  [The |a|-dependence cancels identically.]
 (b) C(K',DM) = [z^{DM}](1+z)^{K'}, and the two-statistic table is a PRODUCT:
        sum_{i,W} T[i,W] u^i y^W = prod_{k=1}^{d} (1 + u y^k)^M     (binomial theorem)
     where T[i,W] = sum_{|a|=i, W(a)=W} prod_k C(M,a_k).
 (c) i!/(M+1+i)! = (1/M!) int_0^1 x^i (1-x)^M dx  (Beta), so the i-sum collapses to u = -x.

## 4. Result
    G[M^D] = ((DM)!/(M!)^D) * [z^{DM}] (1+z)^E * int_0^1 h(x,z)^M dx
    h(x,z) = (1-x) * prod_{k=1}^{d} ((1+z)^k - x)      [x-degree D, z-degree s_d]
    E = ( sum_{k=2}^{D+1} k - s_d ) * M + 2            [D=4: 8M+2 ; D=5: 10M+2]

## 5. Algorithm
 - exact interpolatory quadrature in x with DM+1 nodes (h^M has x-degree DM); nodes must avoid x=1
   (g_j(0) = (1-x_j)^D must be invertible).
 - per node, g_j(z) = h(x_j,z) has degree s_d; the coefficients of g_j^M satisfy the order-s_d
   recurrence from g P' = M g' P:  p_{r+1} = (1/(g_0 (r+1))) sum_{s=1}^{s_d} g_s (Ms - r - 1 + s) p_{r+1-s}.
 - contract on the fly against C(E, DM-r) and the quadrature weights.
 Cost per prime: O(M^2) modular multiplications, O(M) memory.
 (The earlier 2D-table method was O(M^3) time and O(M^2) memory; this is ~500x fewer ops at M=1000
  and ~10^4x less memory. The naive alpha-sum is 1e12 terms.)

## 6. Validation performed
 - closed form vs the defining series (2-var and 5-var)
 - alternating division formula vs exact polynomial division, and vs an independent recursion
   G[m] = C[m+e_2] - sum_{k>=3} G[m+e_2-e_k]
 - all structural lemmas (a),(b),(c) checked individually
 - the final formula vs exact values (D=5, M=1,2,3)
 - the production mod-p code vs exact values, D=4 AND D=5, M=1..5
 - measured timing is cleanly O(M^2)

## 7. Known pitfalls (cost real debugging time)
 - the z-extraction degree and the prefactor use D*M, NOT d*M
 - quadrature node x=1 is forbidden
 - keep p < 2^31 so int64 products stay exact
