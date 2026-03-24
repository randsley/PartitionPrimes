"""
partition_primes.jl

Julia implementation of prime-detecting partition functions from:
  "Integer Partitions Detect the Primes"
  Craig, van Ittersum & Ono (arXiv:2405.06451v2, 2024)

Implements:
  - σ_k(n)           : divisor power sums
  - M_a(n)           : MacMahon partition functions (closed form & direct)
  - M_vec_a(n)       : MacMahonesque partition functions
  - E_i(n)           : prime-detecting expressions from Theorem 1.1 & Table 1
  - is_prime_partition: partition-theoretic primality test

All arithmetic is exact (uses Int / Rational as appropriate).
"""

module PartitionPrimes

export σ, M1, M2, M3, M_direct, M_macmahonesque
export E1, E2, E3, E4, E5, E6
export compare_expressions
export is_prime_partition, verify_range, verify_nonnegativity

# ──────────────────────────────────────────────────────────────────────────────
# 1.  Divisor power sums  σ_k(n) = Σ_{d | n} d^k
# ──────────────────────────────────────────────────────────────────────────────

"""
    σ(k, n) -> Int

Sum of k-th powers of divisors of n.
σ(1, n) is the ordinary sum-of-divisors; σ(0, n) counts divisors.
"""
function σ(k::Int, n::Int)::Int
    n < 1 && throw(ArgumentError("σ(k, n) requires n ≥ 1, got n = $n"))
    s = 0
    for d in 1:isqrt(n)
        if n % d == 0
            s += d^k
            if d != n ÷ d
                s += (n ÷ d)^k
            end
        end
    end
    return s
end

# ──────────────────────────────────────────────────────────────────────────────
# 2.  MacMahon functions via closed forms  (Section 3 of the paper)
#
#     From the Fourier expansions of U_a(q) given in the paper:
#
#     M1(n) = σ_1(n)
#
#     M2(n) = (1/8) [ (-2n+1)·σ_1(n) + σ_3(n) ]
#
#     M3(n) = (1/1920) [ (40n²-100n+37)·σ_1(n)
#                        - 10(3n-5)·σ_3(n)
#                        + 3·σ_5(n) ]
# ──────────────────────────────────────────────────────────────────────────────

"""
    M1(n) -> Int   [= σ_1(n)]
"""
M1(n::Int)::Int = σ(1, n)

"""
    M2(n) -> Rational{Int}
"""
function M2(n::Int)
    return ((-2n + 1) * σ(1, n) + σ(3, n)) // 8
end

"""
    M3(n) -> Rational{Int}
"""
function M3(n::Int)
    return ((40n^2 - 100n + 37) * σ(1, n)
            - 10(3n - 5) * σ(3, n)
            + 3 * σ(5, n)) // 1920
end

# ──────────────────────────────────────────────────────────────────────────────
# 3.  M_a(n) by direct combinatorial computation  (equation 1.2)
#
#     M_a(n) = Σ  m₁·m₂·…·mₐ
#     over all 0 < s₁ < s₂ < … < sₐ,  mᵢ > 0,  with Σ mᵢ sᵢ = n
#
#     We enumerate via recursive descent on the parts.
# ──────────────────────────────────────────────────────────────────────────────

"""
    M_direct(a, n) -> Int

Compute M_a(n) combinatorially.  Exact but exponential in a; practical for a ≤ 5,
n ≤ a few hundred.
"""
function M_direct(a::Int, n::Int)::Int
    a == 0 && return n == 0 ? 1 : 0
    n <= 0 && return 0
    acc = Ref(0)
    _macmahon_rec!(a, n, 0, 1, 1, acc)
    return acc[]
end

# Recursive helper: accumulates into `acc`
# depth   = number of parts chosen so far
# target  = remaining sum
# min_s   = minimum next part size (enforces strict increase)
# prod_so_far = product of multiplicities chosen so far
function _macmahon_rec!(a, target, depth, min_s, prod_so_far, acc)
    if depth == a - 1
        # must use exactly one more part s with multiplicity m = target/s
        for s in min_s:target
            if target % s == 0
                acc[] += prod_so_far * (target ÷ s)
            end
        end
        return
    end
    # choose next part s, multiplicity m ≥ 1
    remaining_parts = a - depth
    for s in min_s:(target ÷ remaining_parts)
        for m in 1:(target ÷ s - (remaining_parts - 1))
            contrib = target - m * s
            contrib > 0 || continue
            _macmahon_rec!(a, contrib, depth + 1, s + 1, prod_so_far * m, acc)
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# 4.  MacMahonesque functions  M_vec_a(n)  (equation 1.3)
#
#     M_vec_a(n) = Σ  m₁^v₁ · m₂^v₂ · … · mₐ^vₐ
#     same summation domain as M_a(n) but with monomial weights given by vec_a.
# ──────────────────────────────────────────────────────────────────────────────

"""
    M_macmahonesque(vec_a, n) -> Int

Compute M_{vec_a}(n) for vec_a = (v₁, v₂, …, vₐ) ∈ ℕᵃ.

Note: M_a(n) = M_{(1,1,…,1)}(n),  M_{(k)}(n) = Σ_{d|n} (n/d)^k  (single-part case).
"""
function M_macmahonesque(vec_a::Vector{Int}, n::Int)::Int
    a = length(vec_a)
    a == 0 && return n == 0 ? 1 : 0
    acc = Ref(0)
    _macmahonesque_rec!(vec_a, a, n, 0, 1, 1, acc)
    return acc[]
end

function _macmahonesque_rec!(vec_a, a, target, depth, min_s, wprod, acc)
    if depth == a - 1
        v = vec_a[depth + 1]
        for s in min_s:target
            if target % s == 0
                m = target ÷ s
                acc[] += wprod * m^v
            end
        end
        return
    end
    v = vec_a[depth + 1]
    remaining_parts = a - depth
    for s in min_s:(target ÷ remaining_parts)
        for m in 1:(target ÷ s - (remaining_parts - 1))
            contrib = target - m * s
            contrib > 0 || continue
            _macmahonesque_rec!(vec_a, a, contrib, depth + 1, s + 1, wprod * m^v, acc)
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# 5.  Prime-detecting expressions  (Theorem 1.1 and Table 1 of the paper)
#
#     Each E_i(n) satisfies:
#       E_i(n) ≥ 0  for all n ≥ 1
#       E_i(n) = 0  ⟺  n is prime   (for n ≥ 2)
# ──────────────────────────────────────────────────────────────────────────────

"""
    E1(n) -> Rational{Int}

Theorem 1.1 (1):  (n²-3n+2)·M₁(n) - 8·M₂(n)
Corresponds to quasimodular form 6H₆.
"""
function E1(n::Int)
    return (n^2 - 3n + 2) * M1(n) - 8 * M2(n)
end

"""
    E2(n) -> Rational{Int}

Theorem 1.1 (2):  (3n³-13n²+18n-8)·M₁(n) + (12n²-120n+212)·M₂(n) - 960·M₃(n)
Corresponds to quasimodular form 36H₈.
"""
function E2(n::Int)
    return (3n^3 - 13n^2 + 18n - 8) * M1(n) +
           (12n^2 - 120n + 212) * M2(n) -
           960 * M3(n)
end

"""
    E3(n) -> Rational{Int}

Table 1, row 3 (quasimodular form 90H₁₀):
(25n⁴-171n³+423n²-447n+170)·M₁(n)
+ (300n³-3554n²+12900n-14990)·M₂(n)
+ (2400n²-60480n+214080)·M₃(n)
- 725760·M₄(n)
"""
function E3(n::Int)
    m4 = M_direct(4, n)
    return (25n^4 - 171n^3 + 423n^2 - 447n + 170) * M1(n) +
           (300n^3 - 3554n^2 + 12900n - 14990) * M2(n) +
           (2400n^2 - 60480n + 214080) * M3(n) -
           725760 * m4
end

"""
    E4(n) -> Rational{Int}

Table 1, row 4 (quasimodular form 90H₁₂).
"""
function E4(n::Int)
    m4 = M_direct(4, n)
    m5 = M_direct(5, n)
    return (126n^5 - 1303n^4 + 5073n^3 - 9323n^2 + 8097n - 2670) * M1(n) +
           (3024n^4 - 48900n^3 + 288014n^2 - 737100n + 695490) * M2(n) +
           (60480n^3 - 1510080n^2 + 10644480n - 23496480) * M3(n) +
           (725760n^2 - 36288000n + 218453760) * m4 -
           580608000 * m5
end

"""
    E5(n) -> Rational

Appendix A.5: unique new prime-vanishing direction at degree 2 in {M₁,…,M₅}.
Normalized to primitive integer coefficients.
Uses M_direct for M₄ and M₅ (slow for large n).
"""
function E5(n::Int)
    m4 = M_direct(4, n)
    m5 = M_direct(5, n)
    bn = big(n)
    return (
        (-450450 + 675675*bn - 225225*bn^2) * M1(n)
        + (960960*bn - 120120*bn^2) * M2(n)
        + (2534912*bn - 166016*bn^2) * M3(n)
        + (7999488*bn - 322560*bn^2) * m4
        + 258048000 * m5
    )
end

"""
    E6(n) -> Rational

Appendix A.6: unique new prime-vanishing direction at degree 4 in {M₁,…,M₇}.
Normalized to primitive integer coefficients.
Uses M_direct for M₄–M₇ (very slow for large n; practical for n ≤ ~50).
"""
function E6(n::Int)
    m4 = M_direct(4, n)
    m5 = M_direct(5, n)
    m6 = M_direct(6, n)
    m7 = M_direct(7, n)
    bn = big(n)
    return (
        (105367732470 - 277943208789*bn + 289613157079*bn^2
         - 145583618619*bn^3 + 28545937859*bn^4) * M1(n)
        + (-51324522800*bn^3 + 4292382480*bn^4) * M2(n)
        + (-40313554176*bn^3 + 1776519808*bn^4) * M3(n)
        + (-32888346624*bn^3 + 770273280*bn^4) * m4
        + (-7741440000*bn^3 - 154828800*bn^4) * m5
        + (-483548921856000 + 37196070912000*bn) * m6
        + 892705701888000 * m7
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# 6.  Partition-theoretic primality test
# ──────────────────────────────────────────────────────────────────────────────

"""
    is_prime_partition(n) -> Bool

Returns true iff n is prime, using the partition-theoretic criterion of
Theorem 1.1(1): n is prime ⟺ (n²-3n+2)·M₁(n) - 8·M₂(n) = 0.

Valid for n ≥ 2 (returns false for n < 2, including n = 0 and n = 1).
"""
function is_prime_partition(n::Int)::Bool
    n < 2 && return false
    return iszero(E1(n))
end

# ──────────────────────────────────────────────────────────────────────────────
# 7.  Verification utilities
# ──────────────────────────────────────────────────────────────────────────────

"""
    is_prime_trial(n) -> Bool

Simple trial-division primality test for verification.
"""
function is_prime_trial(n::Int)::Bool
    n < 2 && return false
    n == 2 && return true
    n % 2 == 0 && return false
    for d in 3:2:isqrt(n)
        n % d == 0 && return false
    end
    return true
end

"""
    verify_range(lo, hi; verbose=true)

Cross-check is_prime_partition against trial division over [lo, hi].
Prints mismatches and a summary.
"""
function verify_range(lo::Int, hi::Int; verbose::Bool = true)
    mismatches = Int[]
    for n in lo:hi
        p_part = is_prime_partition(n)
        p_trial = is_prime_trial(n)
        if p_part != p_trial
            push!(mismatches, n)
        end
    end
    if verbose
        if isempty(mismatches)
            println("✓ Perfect agreement on [$lo, $hi]:  ",
                    count(is_prime_trial, lo:hi), " primes detected.")
        else
            println("✗ Mismatches at: ", mismatches)
        end
    end
    return mismatches
end

"""
    verify_nonnegativity(lo, hi; expressions=[E1, E2, E3, E4], verbose=true)

Verify that each E_i(n) ≥ 0 for all composite n in [lo, hi].
Returns a Dict mapping expression names to lists of violating n values.
"""
function verify_nonnegativity(lo::Int, hi::Int;
                               expressions::Vector = [E1, E2, E3, E4],
                               verbose::Bool = true)
    violations = Dict{String,Vector{Int}}()
    names = [string(e) for e in expressions]
    for (i, E) in enumerate(expressions)
        violations[names[i]] = Int[]
    end
    composites = filter(n -> !is_prime_trial(n), lo:hi)
    for n in composites
        for (i, E) in enumerate(expressions)
            val = E(n)
            if val < 0
                push!(violations[names[i]], n)
            end
        end
    end
    if verbose
        all_ok = all(isempty, values(violations))
        if all_ok
            println("✓ Non-negativity verified for $(join(names, ", ")) ",
                    "on $(length(composites)) composites in [$lo, $hi].")
        else
            for (name, vs) in violations
                isempty(vs) || println("✗ $name is negative at: ", vs)
            end
        end
    end
    return violations
end

"""
    compare_expressions(n)

Evaluate all implemented prime-detecting expressions at n and display results.
Useful for exploring the relationships between expressions.
"""
function compare_expressions(n::Int)
    println("n = $n  ($(is_prime_trial(n) ? "prime" : "composite"))")
    println("  E1(n) = $(E1(n))")
    println("  E2(n) = $(E2(n))")
    println("  E3(n) = $(E3(n))")
    println("  E4(n) = $(E4(n))")
    println()
    println("  Ratio E2/E1: ",
            iszero(E1(n)) ? "0/0 (both zero at primes)" : E2(n) // E1(n))
end

end  # module PartitionPrimes

# ──────────────────────────────────────────────────────────────────────────────
# Demo / quick-start
# ──────────────────────────────────────────────────────────────────────────────

import .PartitionPrimes: σ, M1, M2, M3, M_direct, M_macmahonesque,
    E1, E2, E3, E4, E5, E6,
    is_prime_partition, verify_range, verify_nonnegativity, compare_expressions

println("=" ^ 60)
println("Partition-theoretic prime detection")
println("Craig, van Ittersum & Ono (2024)")
println("=" ^ 60)

println("\n── MacMahon function values ──")
for n in [2, 3, 4, 5, 6, 7, 8, 10, 12, 13]
    println("  n=$(lpad(n,2))  M₁=$(lpad(M1(n),4))  M₂=$(M2(n))  M₃=$(M3(n))")
end

println("\n── Prime-detecting expression E1(n) ──")
println("  (should be 0 iff n is prime)")
for n in 2:20
    val = E1(n)
    marker = iszero(val) ? "  ← PRIME" : ""
    println("  n=$(lpad(n,2))  E1=$(val)$(marker)")
end

println("\n── Partition primality test vs trial division, n ∈ [2, 100] ──")
verify_range(2, 100)

println("\n── Detailed comparison at selected n ──")
for n in [11, 12, 13, 97, 100]
    compare_expressions(n)
end

println("\n── Primes in [2, 50] via partition test ──")
primes_found = filter(is_prime_partition, 2:50)
println("  ", primes_found)

println("\n── MacMahonesque M_{(2,1)}(n) vs M_direct(2,n) for n=6..10 ──")
for n in 6:10
    m21 = M_macmahonesque([2, 1], n)
    m2  = M_direct(2, n)
    println("  n=$n  M_{(2,1)}=$m21  M_{(1,1)}=$m2")
end
