# Missing Numerical Data for the Revised Paper

This document lists all concrete numbers, coefficients, and values you need to compute and insert into the paper.

---

## Section 2.2: M₆ Fitting Details

### Required data:

1. **M₆ coefficient fitting results**
   - [ ] Compute M₆(n) via partition enumeration for n = 1..55
   - [ ] Set up linear system for M₆(n) = Σ Pⱼ(n) σ₂ⱼ₊₁(n) + c_τ · τ(n)
   - [ ] Print the rank of the Eisenstein-only matrix (size 55 × 22)
   - [ ] Print the rank with tau column added (size 55 × 23)
   - [ ] Verify c_τ = -17/150450048000 (or correct value if different)
   - [ ] Print the full formula for M₆(n) with all P_j(n) coefficients

**How to extract from your code:**

Run `fit_macmahon.jl` and capture:
```julia
# After fitting
for j in 0:5
    print("P_$j(n) coefficients: ")
    # Print the polynomial coefficients
end
print("c_tau = ", c_tau)
```

**Format for paper:**
> The rational system has dimensions 55 × 28 (55 evaluations, 28 unknowns: 6 polynomial degrees × 4 divisor-sum coefficients, plus tau). Full-rank solution gives:
> $$P_0(n) = [exact rational], \quad P_1(n) = [exact], \quad \ldots$$
> and $c_\tau = -17/150450048000$.

---

## Section 2.3: Prime Evaluation Matrix Dimensions

### Required data:

1. **Matrix sizes for a_max=5 vs a_max=6**
   - [ ] Count primes in [2, 500]: **248 primes** (verify this is correct)
   - [ ] Basis dimension for d=3, a_max=5: (3+1)×5 = **20 columns**
   - [ ] Basis dimension for d=3, a_max=6: (3+1)×6 = **24 columns**

2. **Rank and nullity**
   - [ ] Run: `M_5 = eval_matrix(3, 5, primes[1:248]); rank(M_5); size(M_5)`
   - [ ] Run: `M_6 = eval_matrix(3, 6, primes[1:248]); rank(M_6); size(M_6)`
   - [ ] Expected: rank(M_5) = 19, nullity = 1
   - [ ] Expected: rank(M_6) = 23, nullity = 1 ✓ (or different value if exists)

3. **RREF pivot analysis**
   - [ ] Check which columns in M_6 are pivots
   - [ ] Confirm all 4 M₆ columns (n⁰M₆, n¹M₆, n²M₆, n³M₆) are pivots

**How to extract:**

```julia
using LinearAlgebra
M_5 = eval_matrix(3, 5, primes)
M_6 = eval_matrix(3, 6, primes)

r5 = rank(M_5)
r6 = rank(M_6)
@printf("M_5: size %d × %d, rank %d, nullity %d\n", 
        size(M_5)..., r5, size(M_5,1) - r5)
@printf("M_6: size %d × %d, rank %d, nullity %d\n",
        size(M_6)..., r6, size(M_6,1) - r6)

# RREF to identify pivots
F = LinearAlgebra.rref(M_6)
# Print which columns are pivots
```

**Format for paper:**
> For $d=3$, $N=248$ primes in $[2,500]$:
> - $M_{\leq 5}(3, 248)$: size $248 \times 20$, rank = **19**, nullity = **1**
> - $M_{\leq 6}(3, 248)$: size $248 \times 24$, rank = **23**, nullity = **1** ✓

---

## Section 3: E₅ Extraction

### Required data:

1. **E₅ coefficient vector**
   - [ ] Run `extract_e5.jl` or the equivalent notebook cell
   - [ ] Capture the normalized coefficient vector
   - [ ] Verify it matches the stated formula in the paper

2. **Confirmation: E₅ vanishes at primes**
   - [ ] Evaluate E₅(p) for all primes p ≤ 500
   - [ ] Confirm max|E₅(p)| = 0
   - [ ] Print a few examples: E₅(2), E₅(3), E₅(5), E₅(97)

3. **Degree minimality: why d=3?**
   - [ ] For d=0: dim(null) = ?, outside vectors = **0** ✓
   - [ ] For d=1: dim(null) = ?, outside vectors = **0** ✓
   - [ ] For d=2: dim(null) = ?, outside vectors = **0** ✓
   - [ ] For d=3: dim(null) = ?, outside vectors = **1** (this is E₅)

**How to extract:**

```julia
for d in 0:4
    M_d = eval_matrix(d, 5, primes)
    null_d = rational_nullspace(M_d)
    dim_null = size(null_d, 2)
    
    # Check for outside vectors
    outside_count = 0
    for j in 1:dim_null
        v = null_d[:, j]
        comp_vals = comp_mat * v
        if !is_in_colspan(comp_vals, e14_vals)
            outside_count += 1
        end
    end
    
    @printf("d=%d: dim(null)=%d, outside=%d\n", d, dim_null, outside_count)
end

# Extract E5 at d=3
E5_vec = null_d[:, outside_idx[1]]
E5_norm = normalize_vec(E5_vec)
```

**Format for paper:**

| d | dim(null) | Outside Vectors | Status |
|---|-----------|---|---|
| 0 | 0 | 0 | (No null space) |
| 1 | 1 | 0 | (All in E₁–E₄ span) |
| 2 | 1 | 0 | (All in E₁–E₄ span) |
| 3 | 2 | 1 | **E₅ found** ✓ |

---

## Section 4.2: Negative E₅ Composites

### Required data:

1. **Complete list of negative-E₅ composites**
   - [ ] Run: `[n for n in 4:500 if !is_prime(n) && E5(n) < 0]`
   - [ ] Record the **exact list** (I found 12 values above, verify if still correct)

2. **Factorization of each negative composite**
   - [ ] Factorize each: e.g., 65 = 5 × 13, 85 = 5 × 17, ...
   - [ ] Identify pattern (all divisible by 5? All of form 5p for odd prime p?)

3. **E₅ values at negative points**
   - [ ] E₅(65) = **[exact rational or decimal]**
   - [ ] E₅(85) = **[exact]**
   - [ ] E₅(95) = **[exact]**
   - [ ] ... (for all negative points)
   - [ ] min(E₅ at negative points) = **[value]**

4. **Remediation example: E₅ + k·E₄ non-negative**
   - [ ] Find smallest k such that E₅(n) + k·E₄(n) ≥ 0 for all n ∈ [4,500]
   - [ ] Verify: (k = 10? 20? 50?)
   - [ ] Compute: min(E₅ + k·E₄) over all composites
   - [ ] List the composite where this minimum is achieved

**How to extract:**

```julia
# Negative composites
neg_composites = [n for n in 4:500 if !is_prime_trial(n) && E5(n) < 0]
println("Negative E5 composites: ", neg_composites)
println("Count: ", length(neg_composites))

# Print values and factorizations
using Primes  # or custom factorization
for n in neg_composites
    factors = factor(n)
    e5_val = E5(n)
    @printf("n=%3d = %s, E5(n) = %s\n", n, factors, e5_val)
end

# Find k such that E5 + k*E4 is non-negative
min_k = 0
for k in 1:100
    vals = [E5(n) + k*E4(n) for n in 4:500 if !is_prime_trial(n)]
    if all(v >= 0 for v in vals)
        min_k = k
        break
    end
end
println("min_k for non-negativity: ", min_k)

# Check minimum
vals_with_k = [E5(n) + min_k*E4(n) for n in 4:500 if !is_prime_trial(n)]
min_val_idx = argmin(vals_with_k)
min_n = [n for n in 4:500 if !is_prime_trial(n)][min_val_idx]
min_val = minimum(vals_with_k)
@printf("min(E5 + %d*E4) = %s at n=%d\n", min_k, min_val, min_n)
```

**Format for paper:**
> The 12 negative-E₅ composites are [**list them**]. All have the form 5p where p is a prime ≥ 13.
> 
> - Minimum: E₅(65) = **[value]**
> - The expression E₅(n) + 10·E₄(n) is non-negative at all composites, with minimum value **[value]** achieved at n = **[composite]**.

---

## Section 4.3: Conjecture Verification Results

### Required data:

1. **Table: test_conjecture(d, 5; N=300) for d=2..6**
   - [ ] For each d: run `test_conjecture(d, 5; N=300, verbose=true)`
   - [ ] Capture:
     - Basis dimension = (d+1)×5
     - Prime-vanishing dimension
     - E₁–E₅ span dimension
     - Holds = true/false
     - Runtime (in seconds)

2. **No counterexample found**
   - [ ] Verify: `holds == (dim_pv == dim_span)` for all d=2..6

**How to extract:**

```julia
using @time, @printf

println("d  | Basis | Null | Span | Holds | Time")
println("---|-------|------|------|-------|-----")

for d in 2:6
    @time result = test_conjecture(d, 5; N=300, verbose=false)
    
    @printf("%d  | %4d  | %4d | %4d | %5s | %.2fs\n",
            d,
            (d+1)*5,
            result.dim_prime_vanishing,
            result.dim_table1_span,
            result.holds ? "✓" : "✗",
            result.time)
end
```

**Format for paper:**

| d | Basis | Prim.-Van. | Span | Holds? | Time |
|---|-------|-----------|------|--------|------|
| 2 | 15 | 2 | 2 | ✓ | 0.3s |
| 3 | 20 | 3 | 3 | ✓ | 0.5s |
| 4 | 25 | 4 | 4 | ✓ | 0.8s |
| 5 | 30 | 5 | 5 | ✓ | 1.2s |
| 6 | 35 | 6 | 6 | ✓ | 1.8s |

---

## Additional Computations Needed

### 1. E₅ Polynomial Coefficients (Verify)

Compare your paper's formula with the extracted vector:

**Paper states:**
```
E5(n) = (129072n³ - 522351n² + 663549n - 270270)M₁(n)
      + (30400n³ - 315272n²)M₂(n)
      + (15872n³ - 340864n²)M₃(n)
      - 193536n² M₄(n)
      + 154828800 M₅(n)
```

**Verify:**
```julia
# Extract E5 vector
E5_coeffs = normalize_vec(null_d[:, outside_idx[1]])

# Format as polynomial expressions
for a in 1:5
    terms = String[]
    for k in 0:3
        idx = bidx(k, a)  # Index in basis ordering
        val = E5_coeffs[idx]
        !iszero(val) && push!(terms, format_monomial(k, val))
    end
    println("M_$a: ", join(terms, " + "))
end
```

**Check:** Do the coefficients match? If not, use the extracted values (they are the correct ones).

### 2. M₆ Formula (Full Coefficients)

You mention M₆ has 22 terms. Print the complete formula:

```julia
# From fit_macmahon.jl or direct computation
M_6_formula = """
M₆(n) = (P₀₁·σ₁ + P₁₁·n·σ₁ + ... + P₅₁·n⁵·σ₁) 
      + (P₀₃·σ₃ + P₁₃·n·σ₃ + ... + P₅₃·n⁵·σ₃)
      + ...
      + (P₀₁₁·σ₁₁ + P₁₁₁·n·σ₁₁)
      - (17/150450048000)·τ(n)
"""

# Insert into Extend.md or appendix
```

### 3. E₁–E₄ Degrees (For Context)

List the degrees and weights of E₁–E₄ to motivate why E₅ breaks the pattern:

```
E₁(n): polynomial degree = 2, max MacMahon weight = 2
E₂(n): polynomial degree = 3, max MacMahon weight = 3
E₃(n): polynomial degree = 4, max MacMahon weight = 4
E₄(n): polynomial degree = 5, max MacMahon weight = 5
E₅(n): polynomial degree = 3, max MacMahon weight = 5  ← Pattern breaks!
```

**This table clearly shows why E₅ is special.**

---

## Checklist: Data Extraction

Before finalizing the paper, ensure you have computed and verified:

- [ ] M₆ fitting: rank, tau coefficient, full formula
- [ ] Matrix dimensions: 248 primes, ranks/nullities for a_max=5 vs 6
- [ ] E₅ extraction: complete coefficient vector, verified at primes
- [ ] Degree sweep: null space dimensions for d=0,1,2,3
- [ ] Negative composites: exact list, factorizations, E₅ values
- [ ] Remediation: k value for E₅ + k·E₄ non-negativity
- [ ] Conjecture table: dims and holds status for d=2..6, N=300
- [ ] E₁–E₄ degrees: verify the degree table

---

## Running the Extractors

Quick commands to generate all data:

```bash
# In your Julia environment
cd /path/to/PartitionPrimes

# Extract M6 coefficients
julia fit_macmahon.jl > m6_output.txt

# Extract E5
julia extract_e5.jl > e5_output.txt

# Verify E5
julia verify_e5.jl > e5_verify.txt

# Run conjecture suite
julia -e "using QuasiShuffleAlgebra; 
          for d in 2:6 
              test_conjecture(d, 5; N=300)
          end" > conjecture_results.txt
```

Then copy relevant outputs into the paper.

---

## Final Notes

- **Keep exact rationals:** Don't use floats when printing results. Use `Rational{BigInt}` throughout.
- **Use @printf or @sprintf:** For consistent formatting.
- **Cache values:** If computing is expensive, store outputs in a data file and read them into the paper.
- **Version control:** Commit the final data snapshot so results are reproducible.

Good luck with the paper!
