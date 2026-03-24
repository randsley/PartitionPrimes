# PartitionPrimes

[![DOI](https://zenodo.org/badge/1175396956.svg)](https://doi.org/10.5281/zenodo.19118619)

Julia implementation of partition-theoretic prime detection via the quasi-shuffle algebra,
extending the framework of Craig, van Ittersum & Ono [1]. This repository implements six
prime-detecting expressions $E_1$–$E_6$ built from MacMahon partition functions, where for
$n \geq 2$: $E_i(n) = 0$ if and only if $n$ is prime — connecting additive number theory
(partitions) to multiplicative number theory (primes) via quasimodular forms.

The three companion papers [2, 3, 4] extend the original framework by deriving $E_5$ and $E_6$,
establishing the weight-12 cusp obstruction, and proving that $E_6$ is recoverable via
tau-cancellation between $M_6$ and $M_7$.

## Repository Layout

```
partition_primes.jl         Standalone module: σ, M1–M3, E1–E6, primality test

QuasiShuffleAlgebra/        Julia package (full implementation)
  src/
    QuasiShuffleAlgebra.jl  Module entry point
    util.jl                 Word, ZqElem types and helpers
    bernoulli.jl            Exact Bernoulli numbers (cached, Rational{BigInt})
    diamond.jl              Diamond product (eq. 4.3)
    quasishuffle.jl         Memoized quasi-shuffle product
    macmahon.jl             σ, M1–M7, M_direct, M_macmahonesque
    d_operator.jl           D operator via exact RREF over Q
    symmetrisation.jl       Symmetrisation (Theorem 4.4)
    prime_detection.jl      E1–E6, is_prime_partition, verify_range
    conjecture.jl           Computational test of the open conjecture
  test/                     906 tests, ~2 minutes
  examples/demo.jl          Annotated demo of every feature

notebooks/                  Jupyter notebooks (Julia kernel)
  PartitionPrimes_Tutorial  Introductory tutorial (uses partition_primes.jl)
  paper1_obstruction        Companion to Paper 1: Weight-12 Cusp Obstruction
  paper2_recovery           Companion to Paper 2: Cusp Cancellation & E₆
  paper3_computational      Companion to Paper 3: Computational Study
  E5_Exploration            E₅ deep dive (uses QuasiShuffleAlgebra)
  QuasiShuffleAlgebra_Guide Package API guide

scripts/                    Standalone derivation and verification scripts
  extract_e5.jl             E₅ extraction via prime evaluation matrix
  find_e5_direct.jl         Direct null-space search for E₅
  fit_macmahon.jl           Closed-form fitting for M₄, M₅, M₆
  verify_e5.jl              Independent verification of E₅

Paper/                      LaTeX sources for the three papers
Extend.md                   Implementation notes and mathematical background
```

## Quick Start

```julia
using Pkg
Pkg.develop(path="QuasiShuffleAlgebra")
using QuasiShuffleAlgebra

# Partition-theoretic primality test
is_prime_partition(97)    # true
is_prime_partition(100)   # false

filter(is_prime_partition, 2:30)
# [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

# Prime-detecting expression E1 (Theorem 1.1)
# E1(n) = (n²-3n+2)·M1(n) - 8·M2(n) = 0  iff  n is prime
E1(13)   # 0//1
E1(12)   # positive rational

# Quasi-shuffle product (algebra Z_q)
quasishuffle_words([1], [1, 1])
# 3·[1,1,1] + (1/6)·[3,1] + (1/6)·[1,3] - (1/3)·[1,1]

# Computational conjecture test
test_conjecture(3, 4; N=150, verbose=false)
# (holds=true, dim_basis=16, dim_prime_vanishing=6, dim_table1_span=16, ...)
```

## Key Mathematical Objects

**MacMahon functions** $M_a(n)$ — weighted sums over strict $a$-part partitions of $n$:

$$M_a(n) = \sum_{\substack{0 < s_1 < \cdots < s_a \\ m_i \geq 1,\; \sum m_i s_i = n}} m_1 m_2 \cdots m_a$$

**Prime-detecting expressions** (Theorem 1.1 of [1], extended in [2, 3, 4]):

| Expression | Basis | Degree | Source | Non-negative? |
|---|---|---|---|---|
| $E_1(n)$ | $\{M_1\}$ | $d=1$ | CIO [1] | Yes |
| $E_2(n)$ | $\{M_1, M_2\}$ | $d=2$ | CIO [1] | Yes |
| $E_3(n)$ | $\{M_1, \ldots, M_3\}$ | $d=2$ | CIO [1] | Yes |
| $E_4(n)$ | $\{M_1, \ldots, M_4\}$ | $d=3$ | CIO [1] | Yes |
| $E_5(n)$ | $\{M_1, \ldots, M_5\}$ | $d=2$ | This work [2] | No (shifted by $+856 E_4$) |
| $E_6(n)$ | $\{M_1, \ldots, M_7\}$ | $d=4$ | This work [3, 4] | No |

All six vanish if and only if $n$ is prime for $n \geq 2$. Note that $E_5$ skips $M_6$
due to the weight-12 cusp obstruction [3], and $E_6$ requires $M_7$ to cancel the
Ramanujan $\tau$-function contribution from $M_6$ [4].

**Bounded completeness:** No independent prime-vanishing expressions beyond $E_1$–$E_6$
were found in an exhaustive search through $a_{\max} = 12$ (weight 24) and polynomial
degree $d \leq 8$. See `Extend.md` and Paper 3 [2] for details.

## Implementation Notes

- All algebra is exact: `Rational{BigInt}` throughout, no floating point.
- **Weight convention**: in `M_macmahonesque(vec_a, n)`, `vec_a[1]` is the exponent of the
  *largest* part — matching every explicit formula in the paper, despite the definition
  notation going smallest-to-largest.
- The D operator solver uses RREF on `[A | b]` over $\mathbb{Q}$, returning a particular solution
  when the system is underdetermined (MacMahonesque basis functions can be linearly dependent).
- The conjecture test checks the $\mathbb{Q}[n]$-span by including $n^j \cdot E_i(n)$ generators
  for $j = 0, \ldots, d$, evaluated directly rather than via truncated coefficient vectors.

## References

[1] Craig, W., van Ittersum, J.-W., & Ono, K. (2024).
*Integer Partitions Detect the Primes.*
Proc. Natl. Acad. Sci. USA **121**, e2409816121.
[arXiv:2405.06451v2](https://arxiv.org/abs/2405.06451)

[2] Randsley, N. (2025).
*Extending Partition-Theoretic Prime Detection: A Computational Study of the Weight-12
Barrier, the Expression $E_5$, and a Basis-Dependent Recovery of $E_6$.*
(Paper 3 — main computational paper)

[3] Randsley, N. (2025).
*The Weight-12 Cusp Obstruction in Partition-Theoretic Prime Detection.*
(Paper 1 — obstruction theorem)

[4] Randsley, N. (2025).
*Cusp Cancellation and the First Recovery of Prime-Vanishing Relations Beyond Weight 12.*
(Paper 2 — recovery of $E_6$)

### Additional references

- Kang, S.-Y., Matsusaka, T., & Shin, S. (2025). Quasi-modularity in MacMahon partition variants and prime detection. *Ramanujan J.* **67**, 12.
- van Ittersum, J.-W., Mauth, N., Ono, K., & Singh, A. (2025). Quasimodular forms that detect primes are Eisenstein. Preprint.
- Kaneko, M. & Zagier, D. (1995). A generalized Jacobi theta function and quasimodular forms. *The Moduli Space of Curves*, Birkhäuser.
- Zagier, D. (2008). Elliptic modular forms and their applications. *The 1-2-3 of Modular Forms*, Springer.

## License

MIT License — Copyright 2026 Nigel Randsley
