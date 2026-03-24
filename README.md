# PartitionPrimes

[![DOI](https://zenodo.org/badge/1175396956.svg)](https://doi.org/10.5281/zenodo.19118619)

Julia implementation of the quasi-shuffle algebra for partition-theoretic prime detection, based on:

> **"Integer Partitions Detect the Primes"**
> William Craig, Jan-Willem van Ittersum & Ken Ono
> [arXiv:2405.06451v2](https://arxiv.org/abs/2405.06451) (July 2024)

The paper proves that for $n \geq 2$, certain non-negative expressions built from MacMahon
partition functions vanish **if and only if $n$ is prime** — connecting additive number
theory (partitions) to multiplicative number theory (primes) via quasimodular forms.

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

**Prime-detecting expressions** (Theorem 1.1, Table 1):

| Expression | Formula | Vanishes iff |
|---|---|---|
| $E_1(n)$ | $(n^2-3n+2)M_1 - 8M_2$ | $n$ is prime |
| $E_2(n)$ | $(3n^3-13n^2+18n-8)M_1 + \cdots - 960M_3$ | $n$ is prime |
| $E_3(n)$ | polynomial coefficients up to $M_4$ | $n$ is prime |
| $E_4(n)$ | polynomial coefficients up to $M_5$ | $n$ is prime |
| $E_5(n)$ | polynomial coefficients up to $M_5$ (derived computationally) | $n$ is prime |

**Open conjecture:** Any non-negative prime-vanishing expression in $\mathbb{Q}[n] \otimes \{M_a\}$
is a $\mathbb{Q}[n]$-linear combination of Table 1 entries. Tested computationally in
`src/conjecture.jl`; confirmed for polynomial degree $\leq 3$ and weight $\leq 5$ using $E_1$–$E_5$.
See `Extend.md` for the derivation of $E_5$.

## Implementation Notes

- All algebra is exact: `Rational{BigInt}` throughout, no floating point.
- **Weight convention**: in `M_macmahonesque(vec_a, n)`, `vec_a[1]` is the exponent of the
  *largest* part — matching every explicit formula in the paper, despite the definition
  notation going smallest-to-largest.
- The D operator solver uses RREF on `[A | b]` over $\mathbb{Q}$, returning a particular solution
  when the system is underdetermined (MacMahonesque basis functions can be linearly dependent).
- The conjecture test checks the $\mathbb{Q}[n]$-span by including $n^j \cdot E_i(n)$ generators
  for $j = 0, \ldots, d$, evaluated directly rather than via truncated coefficient vectors.

## Reference

Craig, W., van Ittersum, J.-W., & Ono, K. (2024).
*Integer Partitions Detect the Primes.*
[arXiv:2405.06451v2](https://arxiv.org/abs/2405.06451)

## License

MIT License — Copyright 2026 Nigel Randsley
