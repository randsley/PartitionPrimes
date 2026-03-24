# Scripts

Standalone Julia scripts used during the derivation and verification of E5 and the closed-form MacMahon functions. These are research artifacts documenting the computational workflow, not library code.

## Files

| Script | Purpose |
|--------|---------|
| `fit_macmahon.jl` | Derive closed-form expressions for M4, M5, M6 as linear combinations of divisor power sums n^k * sigma_{2j-1}(n) and (for M6) the Ramanujan tau function. |
| `extract_e5.jl` | Extract E5 from the null space of the prime evaluation matrix using closed-form M4/M5/M6. Uses a larger prime set (N=600) to ensure M6 columns are independent. |
| `find_e5_direct.jl` | Direct null-space search for E5 at (d=6, a_max=5). All logic wrapped in functions to avoid Julia soft-scope issues in scripts. |
| `verify_e5.jl` | Independent verification of E5 using hard-coded coefficients. Activates the QuasiShuffleAlgebra package and checks E5 against its known formula. |

## Running

Each script is self-contained. Run from the repository root:

```bash
julia scripts/fit_macmahon.jl
julia scripts/extract_e5.jl
julia scripts/find_e5_direct.jl
julia scripts/verify_e5.jl      # requires QuasiShuffleAlgebra
```

Note: `verify_e5.jl` uses a hard-coded path to activate the QuasiShuffleAlgebra package. Adjust the `Pkg.activate(...)` path if your checkout location differs.
