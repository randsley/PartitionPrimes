# Changelog

## v0.1.0 — March 2026

### Repository Organization
- Consolidated all notebooks into `notebooks/` with README
- Moved utility scripts into `scripts/`
- Added companion notebooks for all three papers
- Updated README with current repository layout

### E6 Discovery and Extended Search
- Derived E6: new prime-vanishing expression via M6+M7 tau-cancellation
- Added M8, M9, M10-M12 and cusp forms at weights 16-24
- Searched exhaustively through a_max=12, d=8: no E7-E9 found
- E1-E6 confirmed complete through weight 24
- Three-paper package: obstruction, recovery, and computational study

### E5 Discovery
- Derived E5 at polynomial degree d=2 in {M1,...,M5}
- Added closed-form M4, M5, M6 using divisor power sums and Ramanujan tau
- Established E5+856*E4 as the minimal non-negative shift
- M6 closed form requires tau(n) due to the weight-12 cusp form Delta

### Bug Fixes
- Implemented E5 and E6 in standalone `partition_primes.jl`
- Added `verify_nonnegativity` for E1-E4 non-negativity checking
- Fixed `compare_expressions` to include E4
- Fixed sigma(k, 0) to throw ArgumentError
- Fixed `@printf` and `is_prime_trial` issues in notebooks
- Fixed incorrect convolution identity in Tutorial notebook
- Fixed `macmahon_table` DP bugs in paper companion notebooks
- Fixed anonymous function return type syntax for Julia 1.12

### Foundation
- Implemented quasi-shuffle algebra Z_q with exact rational arithmetic
- MacMahon functions M1-M3 (closed form) and M_direct (combinatorial)
- Prime-detecting expressions E1-E4 from Craig-van Ittersum-Ono (2024)
- D operator, symmetrisation, diamond product
- 906 tests covering all algebraic operations
- Computational conjecture test framework
