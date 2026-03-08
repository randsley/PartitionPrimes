# QuasiShuffleAlgebra.jl

A Julia package implementing the quasi-shuffle algebra for MacMahonesque partition
functions and partition-theoretic prime detection, based on the paper:

> **"Integer Partitions Detect the Primes"**
> William Craig, Jan-Willem van Ittersum, Ken Ono
> [arXiv:2405.06451v2](https://arxiv.org/abs/2405.06451) (July 2024)

The paper proves that primes are exactly the integers n >= 2 where certain linear
combinations of MacMahon partition functions vanish. This package provides:

- Exact computation of MacMahon and MacMahonesque partition functions
- The quasi-shuffle algebra Z_q (formal multiplication of generating series)
- The D operator (expressing n*M_a(n) as constant-coefficient combinations)
- Symmetrisation for recovering quasimodular forms
- Five prime-detecting expressions (E1--E5) from Theorem 1.1 and Table 1

All arithmetic is exact, using `Rational{BigInt}` throughout the algebra layers.

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/QuasiShuffleAlgebra")
```

**Dependencies:** Combinatorics.jl (v1.x), plus the standard library packages
LinearAlgebra, Printf, and Test. Requires Julia 1.6 or later.

## Quick Start

```julia
using QuasiShuffleAlgebra

# Partition-theoretic primality test
is_prime_partition(97)   # true
is_prime_partition(100)  # false

# Find primes in a range
filter(is_prime_partition, 2:50)
# [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Verify against trial division
verify_range(2, 200)
# Passed: perfect agreement on [2, 200], 46 primes detected.

# Compute quasi-shuffle product
quasishuffle_words([1], [1, 1])
# 3*[1,1,1] + (1/6)*[3,1] + (1/6)*[1,3] - (1/3)*[1,1]

# D operator
d_operator([1])
# Expresses n*M_(1)(n) as a linear combination of M_{word}(n)
```

## Mathematical Background

### Divisor Power Sums

The divisor power sum function is:

    sigma_k(n) = sum of d^k over all divisors d of n

Computed by `sigma(k, n)` (exported as the Unicode `σ`).

### MacMahon Partition Functions

The MacMahon partition function M_a(n) counts weighted partitions of n into
a parts with strictly increasing part sizes (equation 1.2 of the paper):

    M_a(n) = sum over 0 < s_1 < s_2 < ... < s_a, m_i > 0, sum(m_i * s_i) = n
             of m_1 * m_2 * ... * m_a

The package provides closed-form evaluations for a = 1, 2, 3 via Fourier
expansions of the generating functions U_a(q):

| Function | Formula |
|----------|---------|
| `M1(n)` | sigma_1(n) |
| `M2(n)` | ((-2n+1)*sigma_1(n) + sigma_3(n)) / 8 |
| `M3(n)` | ((40n^2-100n+37)*sigma_1(n) - 10(3n-5)*sigma_3(n) + 3*sigma_5(n)) / 1920 |

Closed-form evaluations are also provided for a = 4, 5, 6 via fitting to divisor power
sums (and the Ramanujan tau function for a = 6):

| Function | Formula structure |
|----------|-------------------|
| `M4(n)` | linear combination of σ₁, σ₃, σ₅, σ₇ with polynomial-in-n prefactors |
| `M5(n)` | linear combination of σ₁, σ₃, σ₅, σ₇, σ₉ with polynomial-in-n prefactors |
| `M6(n)` | linear combination of σ₁, σ₃, σ₅, σ₇, σ₉, σ₁₁, and τ(n) with polynomial-in-n prefactors |

The M₆ formula requires τ(n) because the weight-12 quasimodular form generating M₆ has
a cusp-form component Δ(q) = q∏(1−q^k)²⁴, whose n-th coefficient is the Ramanujan tau
function. Computed by `ramanujan_tau(n)` (cached, exact BigInt).

For arbitrary a, `M_direct(a, n)` computes M_a(n) by direct partition enumeration.

### MacMahonesque Partition Functions

The MacMahonesque generalisation M_{(v_1,...,v_a)}(n) replaces the uniform
multiplicity product with a monomial weight (equation 1.3):

    M_{(v_1,...,v_a)}(n) = sum of m_1^v_1 * m_2^v_2 * ... * m_a^v_a

over the same summation domain as M_a(n). Note that M_a(n) = M_{(1,...,1)}(n)
where the vector has a ones.

Computed by `M_macmahonesque(vec_a, n)` where `vec_a` is a vector of
non-negative integer exponents. Entries of 0 are valid: m^0 = 1.

**Weight convention:** Despite definition (1.3) assigning v_1 to m_1 (the
multiplicity of the *smallest* part), every formula in the paper — the D
operator example, Ψ₁, Ψ₂, and the Table 1 entries — uses the *opposite*
convention: v_1 is the exponent of the *largest* part. The implementation
follows the paper's formula convention (vec_a[1] = exponent of largest part),
which is the convention needed to make all identities check out numerically.
This was confirmed by cross-checking the paper's formula for nM_{(1,1)}(n)
at n = 4: the formula evaluates correctly only under this convention.

### Prime-Detecting Expressions

For n >= 2, each expression below is non-negative and equals zero **if and only
if n is prime** (Theorem 1.1):

**E1(n)** -- corresponds to quasimodular form 6*H_6:

    (n^2 - 3n + 2)*M_1(n) - 8*M_2(n) >= 0

**E2(n)** -- corresponds to 36*H_8:

    (3n^3 - 13n^2 + 18n - 8)*M_1(n) + (12n^2 - 120n + 212)*M_2(n) - 960*M_3(n) >= 0

**E3(n)** -- corresponds to 90*H_10 (uses M_4 via closed form).

**E4(n)** -- corresponds to 90*H_12 (uses M_4 and M_5 via closed forms).

**E5(n)** -- derived computationally as the unique (up to scaling) prime-vanishing
expression in the degree-3, weight-5 basis outside the Q[n]-span of E1–E4:

    (-270270 + 663549n - 522351n² + 129072n³)·M₁(n)
    + (-315272n² + 30400n³)·M₂(n)
    + (-340864n² + 15872n³)·M₃(n)
    + (-193536n²)·M₄(n)
    + 154828800·M₅(n)

This expression vanishes at all primes and is linearly independent of E1–E4 over Q[n].
Note: E5 can be negative at some composites (e.g. n=65, 85, 95), so it satisfies only
the prime-vanishing condition, not the full non-negativity claimed in Theorem 1.1.
See `Extend.md` for the derivation.

### The Algebra Z_q

The space Z_q is the Q-vector space generated by all MacMahonesque generating
series U_{vec_a}(q) = sum_{n>=1} M_{vec_a}(n) * q^n. It is closed under:

1. **Multiplication** (quasi-shuffle product) -- Theorem 4.2 (Bachmann-Kuhn)
2. **Differentiation** (D = q * d/dq operator)

This means any product or derivative of MacMahonesque series can be expressed as
a finite linear combination of MacMahonesque series with rational coefficients.

### The Quasi-Shuffle Product

The product structure on Z_q is defined recursively on words (equation after 4.3).
If u, v are words with first letters x, y and tails u', v':

    u * v = x(u' * v) + y(u * v') + (x diamond y)(u' * v')

where the **diamond product** combines two letters into a linear combination of
letters using Bernoulli numbers (equation 4.3).

### The D Operator

The D operator expresses n*M_{vec_a}(n) as a constant-coefficient linear
combination of MacMahonesque functions (Theorem 1.4):

    n * M_alpha(n) = sum of c_{beta} * M_beta(n)

where the sum runs over words beta with weight at most |alpha| + length(alpha) + 1.
This allows reducing polynomial-coefficient expressions (like the E_i) to
constant-coefficient form.

### Symmetrisation

For a word vec_a with all-odd entries, the symmetrised series:

    U^sym_{vec_a}(q) = sum over unique permutations sigma of U_{sigma(vec_a)}(q)

is a quasimodular form (Theorem 4.4). These generate all quasimodular forms
(Theorem 1.5).

## API Reference

### Types

```julia
const Word   = Vector{Int}                    # e.g. [1,3] represents exponents (v_1, v_2) = (1, 3)
const ZqElem = Dict{Word, Rational{BigInt}}   # formal Q-linear combination of words
```

A `ZqElem` represents an element of Z_q: a finite linear combination of
MacMahonesque generating series, with exact rational coefficients. For example,
`Dict([1,1] => 2//1, [3] => 1//6)` represents `2*U_{(1,1)} + (1/6)*U_{(3)}`.

### Arithmetic on Z_q Elements

| Function | Signature | Description |
|----------|-----------|-------------|
| `zq_add(a, b)` | `ZqElem, ZqElem -> ZqElem` | Add two elements |
| `zq_scale(c, z)` | `Number, ZqElem -> ZqElem` | Scale by a rational constant |
| `zq_multiply(a, b)` | `ZqElem, ZqElem -> ZqElem` | Multiply via quasi-shuffle product |
| `cleanup!(z)` | `ZqElem -> ZqElem` | Remove zero entries in place |
| `evaluate_zq(z, n)` | `ZqElem, Int -> Rational{BigInt}` | Evaluate at integer n: sum of c_w * M_w(n) |

### Bernoulli Numbers

```julia
bernoulli(n::Int) -> Rational{BigInt}
```

Compute the n-th Bernoulli number exactly. Uses the standard recurrence relation
with caching. Exploits B_n = 0 for odd n > 1.

**Examples:**

```julia
bernoulli(0)   # 1//1
bernoulli(1)   # -1//2
bernoulli(12)  # -691//2730
bernoulli(7)   # 0//1
```

### Diamond Product

```julia
diamond_coeffs(i::Int, j::Int) -> Dict{Int, Rational{BigInt}}
diamond(x::Int, y::Int) -> ZqElem
```

The diamond product combines two single letters into a linear combination of
letters. `diamond_coeffs(i, j)` returns the raw coefficient dictionary mapping
letter index m to its coefficient. `diamond(x, y)` returns the result as a
`ZqElem` with single-letter words.

The coefficients are (adapted from equation 4.3 with MacMahonesque indexing):

- Top coefficient: `c_{i+j+1} = i! * j! / (i+j+1)!`
- Lower coefficients: `c_m = [(-1)^i * C(i,m) + (-1)^j * C(j,m)] * B_{i+j+1-m} / (i+j+1-m)`

**Example:**

```julia
diamond(1, 1)
# Dict([3] => 1//6, [1] => -1//6)
# i.e. e_1 diamond e_1 = (1/6)*e_3 - (1/6)*e_1
```

### Quasi-Shuffle Product

```julia
quasishuffle_words(u::Word, v::Word) -> ZqElem
```

Compute the quasi-shuffle product of two words. The result is a `ZqElem` giving
the linear combination of words that results from multiplying the corresponding
generating series. Results are memoized for performance.

```julia
clear_qs_cache!()
```

Clear the memoization cache (useful for benchmarking or memory management).

**Example:**

```julia
quasishuffle_words([1], [1, 1])
# Dict([1,1,1] => 3, [3,1] => 1//6, [1,3] => 1//6, [1,1] => -1//3)
```

This means U_{(1)} * U_{(1,1)} = 3*U_{(1,1,1)} + (1/6)*U_{(3,1)} + (1/6)*U_{(1,3)} - (1/3)*U_{(1,1)},
which can be verified numerically: for each n, the convolution
sum_{i+j=n} M_{(1)}(i) * M_{(1,1)}(j) equals the right-hand side evaluated at n.

### MacMahon and MacMahonesque Functions

```julia
σ(k::Int, n::Int) -> Int
```

Sum of k-th powers of divisors of n. `σ(1, n)` is the ordinary sum of
divisors; `σ(0, n)` counts divisors.

```julia
M1(n::Int) -> Int           # = σ(1, n)
M2(n::Int) -> Rational      # closed-form via σ₁, σ₃
M3(n::Int) -> Rational      # closed-form via σ₁, σ₃, σ₅
M4(n::Int) -> Rational      # closed-form via σ₁, σ₃, σ₅, σ₇
M5(n::Int) -> Rational      # closed-form via σ₁, σ₃, σ₅, σ₇, σ₉
M6(n::Int) -> Rational      # closed-form via σ₁, σ₃, σ₅, σ₇, σ₉, σ₁₁, τ(n)
ramanujan_tau(n::Int) -> BigInt   # Ramanujan tau function (cached)
```

Closed-form MacMahon functions derived from the Fourier expansions of U_a(q).
M4–M6 use divisor power sums with polynomial-in-n prefactors; M6 additionally
requires the Ramanujan tau function because its generating series has a weight-12
cusp form component.

```julia
M_direct(a::Int, n::Int) -> Int
```

Compute M_a(n) by direct combinatorial enumeration of partitions into a
distinct part sizes with positive multiplicities. Exact but exponential in a;
practical for a <= 5 and n up to a few hundred.

```julia
M_macmahonesque(vec_a::Vector{Int}, n::Int, [T::Type = Int]) -> T
```

Compute M_{vec_a}(n) for arbitrary non-negative integer exponent vectors.
The type parameter T controls the accumulation type:

- `Int` (default): fast integer arithmetic, suitable for prime detection
- `Rational{BigInt}`: exact rational arithmetic, required for D operator
  and quasi-shuffle algebra computations

**Examples:**

```julia
M_macmahonesque([1], 6)          # 12 (= σ_1(6))
M_macmahonesque([2, 1], 6)       # 26
M_macmahonesque([1, 1], 6)       # 15 (= M_direct(2, 6))
M_macmahonesque([3, 0], 10)      # with a zero exponent: m^0 = 1
M_macmahonesque([1], 6, Rational{BigInt})  # 12//1 (exact rational)
```

### D Operator

```julia
d_operator(vec_a::Word) -> ZqElem
```

Express n * M_{vec_a}(n) as a constant-coefficient linear combination of
MacMahonesque functions. Uses coefficient matching: evaluates both sides at
sufficiently many values of n and solves the resulting linear system exactly
over Q using Gaussian elimination.

The weight bound for basis words is |vec_a| + length(vec_a) + 1, and the
maximum word length is length(vec_a) + 1 (from Theorem 1.4).

**Performance note:** The D operator solves a linear system whose size grows
with the weight and length of the input word. For `[1]` it takes under a
second; for `[1,1]` it takes approximately 2 minutes due to the ~83 basis
elements requiring exact rational Gaussian elimination.

**Example:**

```julia
result = d_operator([1])
# Expresses n*M_{(1)}(n) in MacMahonesque basis

# Verify numerically:
for n in 1:10
    lhs = Rational{BigInt}(n) * M_macmahonesque([1], n, Rational{BigInt})
    rhs = evaluate_zq(result, n)
    @assert lhs == rhs
end
```

### Symmetrisation

```julia
symmetrise(vec_a::Word) -> ZqElem
```

Compute the symmetrised series U^sym_{vec_a} by summing over all unique
permutations of vec_a. Each unique permutation appears with coefficient 1.

For vectors with all-odd entries, the result represents a quasimodular form
(Theorem 4.4).

**Example:**

```julia
symmetrise([1, 3])
# Dict([1, 3] => 1, [3, 1] => 1)

symmetrise([1, 1, 3])
# Dict([1, 1, 3] => 1, [1, 3, 1] => 1, [3, 1, 1] => 1)
```

### Conjecture Testing

```julia
test_conjecture(d::Int, a_max::Int; N::Int=300, verbose::Bool=true) -> NamedTuple
```

Tests whether every prime-vanishing expression in the basis
`{n^k · M_a(n) : 0 ≤ k ≤ d, 1 ≤ a ≤ a_max}` lies in the Q[n]-span of the
Table 1 entries E1–E5. Returns a NamedTuple with fields:
`holds` (Bool), `dim_prime_vanishing`, `dim_table1_span`, and
`counterexample` (a coefficient vector, or `nothing`).

```julia
scan_conjecture(d_max::Int, a_max::Int; N::Int=300)
```

Sweeps all `(d, a)` pairs with `1 ≤ d ≤ d_max`, `1 ≤ a ≤ a_max`, printing
a result line for each and stopping at the first counterexample.

**Example:**

```julia
test_conjecture(2, 2; N=100)   # quick check covering E1
test_conjecture(5, 5; N=300)   # covers all Table 1 entries
scan_conjecture(6, 6)          # systematic sweep
```

### Prime Detection

```julia
E1(n::Int) -> Rational     # (n^2-3n+2)*M_1(n) - 8*M_2(n)
E2(n::Int) -> Rational     # uses M_1, M_2, M_3
E3(n::Int) -> Rational     # uses M_1, M_2, M_3, M_4 (closed form)
E4(n::Int) -> Rational     # uses M_1, M_2, M_3, M_4, M_5 (closed forms)
E5(n::Int) -> Rational     # uses M_1, M_2, M_3, M_4, M_5 (closed forms; derived computationally)
```

Prime-detecting expressions from Theorem 1.1 and Table 1. For n >= 2:
- E1–E4 satisfy E_i(n) >= 0 for all n, with E_i(n) = 0 iff n is prime
- E5 vanishes precisely at primes but may be negative at some composites

All five use fast closed-form MacMahon functions (no direct enumeration).

```julia
is_prime_partition(n::Int) -> Bool
```

Returns true if n is prime, using the partition-theoretic criterion E1(n) = 0.
Valid for n >= 2; returns false for n < 2.

```julia
is_prime_trial(n::Int) -> Bool
```

Simple trial-division primality test, provided for cross-verification.

```julia
verify_range(lo::Int, hi::Int; verbose=true) -> Vector{Int}
```

Cross-check `is_prime_partition` against `is_prime_trial` over [lo, hi].
Returns a vector of mismatches (empty if all agree). With `verbose=true`,
prints a summary.

### Word Enumeration

```julia
all_words_up_to_weight(w; max_length=w+1, min_length=1) -> Vector{Word}
```

Generate all words (vectors of non-negative integers) with entries summing to
at most w, with length between min_length and max_length. Used internally by
the D operator to construct the MacMahonesque basis.

## File Structure

```
QuasiShuffleAlgebra/
  Project.toml              # Package metadata and dependencies
  src/
    QuasiShuffleAlgebra.jl  # Module definition, includes, exports
    util.jl                 # Word/ZqElem types, arithmetic helpers, word enumeration
    bernoulli.jl            # Cached exact Bernoulli numbers via recurrence
    diamond.jl              # Diamond product coefficients (eq. 4.3)
    quasishuffle.jl         # Memoized recursive quasi-shuffle product
    macmahon.jl             # σ, ramanujan_tau, M1-M6, M_direct, M_macmahonesque
    d_operator.jl           # D operator via exact rational RREF over Q
    symmetrisation.jl       # Symmetrised series (Theorem 4.4)
    prime_detection.jl      # E1-E5, is_prime_partition, verify_range
    conjecture.jl           # Computational test of the open conjecture
  test/
    runtests.jl             # Test runner
    test_bernoulli.jl       # B_0 through B_20 against known values
    test_quasishuffle.jl    # Paper's explicit example + convolution identities
    test_d_operator.jl      # Numerical verification + paper's explicit decomposition
    test_prime_detection.jl # E1 vs trial division, M_direct consistency
    test_conjecture.jl      # Infrastructure tests for conjecture machinery
  examples/
    demo.jl                 # Interactive demonstration script
```

## Testing

Run the full test suite (906 tests, approximately 2 minutes):

```bash
cd QuasiShuffleAlgebra
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Test Summary

| Test Suite | Tests | What it verifies |
|------------|-------|------------------|
| Bernoulli numbers | 21 | B_0 through B_20 against known values; B_odd = 0 for n > 1 |
| Quasi-shuffle product | 54 | Explicit paper result for [1]*[1,1]; convolution identity for n=1:30; [1]*[1] matches power series convolution for n=1:20 |
| D operator | 150 | d_operator([1]) and [1,1] verified numerically for n=1:50; paper's explicit decomposition of nM_{(1,1)}(n) verified for n=1:50 |
| Prime detection | 637 | E1 = 0 iff prime for n=2:200; E1 >= 0 for n=2:100; E5 vanishes at primes 2–23; E5 >= 0 at small composites; M1 = sigma_1; M_direct = M_macmahonesque for a=1,2,3 and n=1:30 |
| Conjecture infrastructure | 44 | Basis evaluation, RREF, null space, rank, test_conjecture at (d=2,a=2) and (d=3,a_max=5) |

### Running Individual Tests

```julia
using QuasiShuffleAlgebra, Test

# Quick Bernoulli check
@test bernoulli(12) == -691//2730

# Quasi-shuffle ground truth
result = quasishuffle_words([1], [1, 1])
@test result[[1,1,1]] == 3
@test result[[3,1]] == 1//6

# Prime detection
@test is_prime_partition(97)
@test !is_prime_partition(100)
```

## Running the Demo

```bash
cd QuasiShuffleAlgebra
julia --project=. examples/demo.jl
```

The demo script displays MacMahon function values, prime-detecting expressions,
quasi-shuffle products, convolution verification, D operator output, and
symmetrisation results.

## Design Decisions

### Exact Arithmetic

All algebra layers (Bernoulli numbers, diamond product, quasi-shuffle, D operator)
use `Rational{BigInt}` throughout. This avoids silent overflow that would occur
with `Rational{Int}` at weight >= 10, where Bernoulli number denominators become
large (e.g., B_12 = -691/2730).

The `M_macmahonesque` function accepts a type parameter so it can return either
`Int` (fast, for prime detection) or `Rational{BigInt}` (for algebra computations).

### Words with Zero Entries

Words like `[3, 0]` are distinct from `[3]`. The former is a two-part
MacMahonesque function where the *largest* part's multiplicity is cubed and the
smallest part's multiplicity appears to the power 0 (= 1). The D operator
produces such words naturally.

Under the weight convention described above, `M_macmahonesque([3, 0], n)` sums
m_largest^3 · m_smallest^0 = m_largest^3 over 2-strict-part partitions of n.

### Quasi-Shuffle Memoization

The recursive quasi-shuffle computation is memoized in a module-level cache.
Without memoization, the recursion is exponential in word length. The cache key
is `Tuple{Word, Word}`. Note that `quasishuffle(a, b) != quasishuffle(b, a)` in
general, so pairs are not sorted.

### D Operator: Coefficient Matching

The D operator uses "Option B" from the implementation plan: rather than
symbolically deriving the expansion via Ramanujan's differential equations, it
evaluates both sides numerically at enough points and solves the resulting
overdetermined linear system exactly. This is more robust and easier to verify,
at the cost of being computationally expensive for larger words.

The system is solved using custom RREF-based Gaussian elimination over
`Rational{BigInt}`. Julia's built-in `\` does not produce exact results for
rational matrices. The solver handles underdetermined systems (which arise because
some MacMahonesque basis functions are linearly dependent as integer sequences):
it runs RREF on the augmented matrix `[A | b]`, checks consistency, and returns
a particular solution with free variables set to zero. This means d_operator finds
*a* valid decomposition, which may differ in its word representation from the
specific decomposition listed in the paper — but both are numerically equivalent.

### Diamond Product Index Convention

The diamond product formula from equation 4.3 of the paper uses an index
convention where the top letter is `i + j + 1` (not `i + j`). This shift
arises from matching MacMahonesque exponents to the Bachmann-Kuhn algebra
generators. The implementation was validated against the paper's explicit
quasi-shuffle example: `[1] * [1,1] = 3*[1,1,1] + (1/6)*[3,1] + (1/6)*[1,3] - (1/3)*[1,1]`.

## Verification Identities

These identities from the paper serve as ground truth for testing:

**Quasi-shuffle (p.12):**

    U_{(1)} * U_{(1,1)} = 3*U_{(1,1,1)} + (1/6)*U_{(3,1)} + (1/6)*U_{(1,3)} - (1/3)*U_{(1,1)}

**Convolution (Theorem 1.4.1):**

    sum_{i+j=n} M_{(1)}(i) * M_{(1)}(j) = (1/6)*M_{(3)}(n) + 2*M_{(1,1)}(n) - (1/6)*M_{(1)}(n)

**Prime detection (Theorem 1.1.1):**

    (n^2 - 3n + 2)*M_1(n) - 8*M_2(n) = 0  iff  n is prime  (n >= 2)

All three are verified exactly in the test suite.

## Open Conjecture

The paper states but does not prove:

> Let P̃(x) = (p_1(x), …, p_a(x)) ∈ ℤ[x]^a be relatively prime integer polynomials.
> For n ≥ 2, suppose E(n) = p_1(n)M_1(n) + ⋯ + p_a(n)M_a(n) ≥ 0 and vanishes
> precisely on the primes. Then E(n) is a ℚ[n]-linear combination of the entries
> in Table 1.

This is a finite-dimensional linear algebra question for each fixed degree bound d
(polynomial coefficients) and weight bound a_max (MacMahon functions). The function
`test_conjecture(d, a_max)` tests it computationally by:

1. Building the basis `{n^k · M_a(n) : 0 ≤ k ≤ d, 1 ≤ a ≤ a_max}`
2. Finding the prime-vanishing subspace (null space of the evaluation matrix at
   primes, computed exactly over ℚ)
3. Checking whether every basis vector of that null space lies in the ℚ-span of
   Table 1 entries, using rational row reduction

**Status:** Confirmed holds at `(d=3, a_max=5)` using E1–E5 (verified in the test suite).
Larger ranges require E5 in the Table 1 span; without E5 there is always a counterexample
at `a_max ≥ 5`. A counterexample that persists even with E5 would be a significant
mathematical finding.

**Key computational discovery:** M₆ columns are always pivot columns in the prime
evaluation matrix RREF — meaning no prime-vanishing expression at any tested `(d, a_max)`
involves M₆. The E5 formula consequently involves only M₁–M₅, not following the pattern
of E1–E4 each introducing the next MacMahon function.

## References

- Craig, van Ittersum, Ono: [Integer Partitions Detect the Primes](https://arxiv.org/abs/2405.06451) (2024)
- Bachmann, Kuhn: *The algebra of generating functions for multiple divisor sums and applications to multiple zeta values*. Ramanujan J. 40 (2016), no. 3, 605--648
- Bachmann: [Lecture notes on multiple zeta values and modular forms](https://www.henrikbachmann.com/uploads/7/7/6/3/77634444/mzv_mf_2020_v_5_4.pdf)
- Andrews, Rose: *MacMahon's sum-of-divisors functions, Chebyshev polynomials, and quasi-modular forms*. J. Reine Angew. Math. 676 (2013), 97--103
- Zagier: *The 1-2-3 of Modular Forms*. In: Bruinier et al., Universitext, Springer 2008
