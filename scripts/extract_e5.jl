"""
extract_e5.jl - Phase 3

Use larger N for prime evaluation to ensure M_6 columns are independent.
Only build the prime evaluation matrix (fast: closed forms for M4/M5/M6).
"""

import Pkg
Pkg.activate("/Users/nigelrandsley/GitHub/PartitionPrimes/QuasiShuffleAlgebra")
using QuasiShuffleAlgebra
using Printf

d     = 6
a_max = 6
N     = 500

ns        = collect(2:N)
primes_ns = filter(is_prime_trial, ns)
@printf("Primes in [2,%d]: %d; Basis dim: %d\n", N, length(primes_ns), (d+1)*a_max)

bidx(k, a) = (a-1)*(d+1) + k + 1

function show_vector(c, label)
    println("\n$label:")
    for a in 1:a_max
        parts = String[]
        for k in 0:d
            val = c[bidx(k, a)]
            iszero(val) && continue
            push!(parts, k == 0 ? "$val" : k == 1 ? "$(val)·n" : "$(val)·n^$k")
        end
        isempty(parts) || println("  M_$a: ", join(parts, " + "))
    end
    ma = [a for a in 1:a_max if any(!iszero(c[bidx(k,a)]) for k in 0:d)]
    println("  M_a present: $ma")
end

function normalize_vec(c)
    dens = [denominator(v) for v in c if !iszero(v)]
    isempty(dens) && return c
    lcm_d = reduce(lcm, dens)
    c_int = [numerator(v * lcm_d) for v in c]
    nonzero = c_int[c_int .!= 0]
    gcd_n = reduce(gcd, nonzero)
    return c_int .// gcd_n
end

println("Building prime evaluation matrix (N=$N)...")
@time prime_mat = eval_matrix(d, a_max, primes_ns)
@printf("Matrix size: %d × %d\n", size(prime_mat)...)

println("Computing null space...")
@time null_vecs = rational_nullspace(prime_mat)
dim_pv = size(null_vecs, 2)
@printf("Prime-vanishing subspace dim: %d\n", dim_pv)

# Check which vectors involve M_6
println("\nChecking M_6 presence in all null vectors...")
has_m6_any = [j for j in 1:dim_pv if any(!iszero(null_vecs[bidx(k,6), j]) for k in 0:d)]
@printf("%d null vectors involve M_6\n", length(has_m6_any))

# E1-E4 span check using composites 4..100
composites = filter(n -> !is_prime_trial(n), 4:100)
println("Building composite eval matrix ($(length(composites)) points)...")
@time comp_mat = eval_matrix(d, a_max, composites)

E_funcs = [E1, E2, E3, E4]
e14_vals = zeros(Rational{BigInt}, length(composites), 4*(d+1))
for (i, n) in enumerate(composites)
    for (fi, Ef) in enumerate(E_funcs)
        val = Rational{BigInt}(Ef(n))
        for j in 0:d
            e14_vals[i, j*4 + fi] = Rational{BigInt}(big(n)^j) * val
        end
    end
end

println("Finding null vectors outside E1-E4 span...")
outside_idx = Int[]
@time for j in 1:dim_pv
    coeff_vec = null_vecs[:, j]
    expr_vals = comp_mat * coeff_vec
    if !is_in_colspan(expr_vals, e14_vals)
        push!(outside_idx, j)
    end
end
@printf("%d null vectors outside E1-E4 Q[n]-span\n", length(outside_idx))

norm_vecs = [normalize_vec(null_vecs[:, j]) for j in outside_idx]

println("\nAll outside vectors:")
for (i, c) in enumerate(norm_vecs)
    show_vector(c, "#$i")
    has_m6 = any(!iszero(c[bidx(k,6)]) for k in 0:d)
    println("  Has M_6: $has_m6")
    cvals = comp_mat * c
    println("  Min/max at composites 4..100: $(minimum(cvals)) / $(maximum(cvals))")
end

# Find valid E5: non-negative at composites, not in E1-E4 span
# Try -c if c is all negative
println("\nSearching for non-negative candidate...")
for (i, c) in enumerate(norm_vecs)
    cvals = comp_mat * c
    if all(v >= 0 for v in cvals)
        println("✓ Vec #$i is non-negative at composites 4..100")
        show_vector(c, "E5 candidate (vec #$i)")
    elseif all(v <= 0 for v in cvals)
        println("✓ -Vec #$i is non-negative at composites 4..100")
        show_vector(-c, "E5 candidate (-vec #$i)")
    else
        println("  Vec #$i is mixed sign — need linear combination")
    end
end
