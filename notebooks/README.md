# Notebooks

Julia notebooks for exploring and verifying the mathematics in the PartitionPrimes project.

## Foundational

| Notebook | Description |
|----------|-------------|
| [PartitionPrimes_Tutorial.ipynb](PartitionPrimes_Tutorial.ipynb) | Introduction to partition-theoretic prime detection. Walks through divisor sums σ_k(n), MacMahon partition functions M_1–M_3, prime-detecting expressions E_1–E_3, and the partition primality test. Depends on `partition_primes.jl`. |
| [QuasiShuffleAlgebra_Guide.ipynb](QuasiShuffleAlgebra_Guide.ipynb) | API guide for the `QuasiShuffleAlgebra.jl` package. Covers the type system (Word, ZqElem), MacMahon function computation, prime detection, and algebraic operations. |

## Paper Companions

Each paper has a dedicated notebook that implements its key results with exact rational arithmetic, allowing readers to follow the mathematical arguments computationally.

| Notebook | Paper | Key content |
|----------|-------|-------------|
| [paper3_computational.ipynb](paper3_computational.ipynb) | Paper 1: *Computational Study (main paper)* | Systematic search for all prime-vanishing directions up to a_max=10, d=6. Reproduces E_1–E_4 from CIO (2024), derives E_5 (d=2 in {M_1,…,M_5}) and E_6 (d=4 in {M_1,…,M_7}), and confirms bounded completeness. Self-contained. |
| [E5_Exploration.ipynb](E5_Exploration.ipynb) | Paper 1 (supplement) | Detailed exploration of E_5: closed-form M_4/M_5 implementations, sign structure analysis, non-negativity remediation (E_5 + 856·E_4), and the M_6 pivot discovery. Depends on `QuasiShuffleAlgebra.jl`. |
| [paper1_obstruction.ipynb](paper1_obstruction.ipynb) | Paper 2: *The Weight-12 Cusp Obstruction* | Demonstrates that the Ramanujan τ-function and the Bernoulli prime 691 force a structural obstruction at weight 12. Verifies that no prime-vanishing expression in {M_1,…,M_6} can involve M_6. Self-contained (defines all functions inline). |
| [paper2_recovery.ipynb](paper2_recovery.ipynb) | Paper 3: *Cusp Cancellation and Recovery of E_6* | Shows that adjoining M_7 makes the weight-12 obstruction cancellable. Verifies that E_6 first appears at polynomial degree d=4 and the obstruction space has dimension d+2. Self-contained. |

## Dependencies

- **Tutorial notebook**: requires `partition_primes.jl` (in repo root)
- **Paper companion notebooks** (paper1, paper2, paper3): self-contained — all functions defined inline
- **E5_Exploration**: requires the `QuasiShuffleAlgebra` package (in repo)
- **QuasiShuffleAlgebra_Guide**: requires the `QuasiShuffleAlgebra` package (in repo)
- All notebooks require Julia 1.6+ and the `Plots` package

## References

- Craig, W., van Ittersum, J.-W., & Ono, K. (2024). *Integer Partitions Detect the Primes.* Proc. Natl. Acad. Sci. USA **121**.
- Randsley, N. (2025). *Extending Partition-Theoretic Prime Detection* (Paper 1 — computational study).
- Randsley, N. (2025). *The Weight-12 Cusp Obstruction* (Paper 2).
- Randsley, N. (2025). *Cusp Cancellation and Recovery of E6* (Paper 3).
