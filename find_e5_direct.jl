"""
find_e5_direct.jl

Extract E5 from (d=6, a_max=5) null space.
All wrapped in functions to avoid Julia soft-scope issues.
"""

import Pkg
Pkg.activate("/Users/nigelrandsley/GitHub/PartitionPrimes/QuasiShuffleAlgebra")
using QuasiShuffleAlgebra
using Printf

function normalize_vec(c)
    dens = [denominator(v) for v in c if !iszero(v)]
    isempty(dens) && return c
    lcm_d = reduce(lcm, dens)
    c_int = [numerator(v * lcm_d) for v in c]
    nonzero = c_int[c_int .!= 0]
    gcd_n = reduce(gcd, nonzero)
    return c_int .// gcd_n
end

function find_e5()
    d = 6; a_max = 5
    N = 300
    primes_ns = filter(is_prime_trial, 2:N)
    @printf("d=%d, a_max=%d, N=%d: %d primes, basis dim=%d\n",
            d, a_max, N, length(primes_ns), (d+1)*a_max)

    bidx(k, a) = (a-1)*(d+1) + k + 1

    println("Building prime evaluation matrix...")
    @time prime_mat = eval_matrix(d, a_max, primes_ns)
    println("Computing null space...")
    @time null_vecs = rational_nullspace(prime_mat)
    dim_pv = size(null_vecs, 2)
    @printf("Null space dim: %d\n", dim_pv)

    # Find outside-E1..E4 part
    composites = filter(n -> !is_prime_trial(n), 4:120)
    comp_mat = eval_matrix(d, a_max, composites)
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

    outside_idx = [j for j in 1:dim_pv
                   if !is_in_colspan(comp_mat * null_vecs[:, j], e14_vals)]
    @printf("Outside E1-E4 span: %d vectors\n", length(outside_idx))

    m5_const_col = bidx(0, 5)
    println("M_5 terms in outside vectors:")
    for j in outside_idx
        v = null_vecs[:, j]
        m5t = [(k, v[bidx(k, 5)]) for k in 0:d if !iszero(v[bidx(k, 5)])]
        println("  #$j: $m5t")
    end

    # Find base E5: outside vector with M_5 constant (k=0) nonzero
    e5_idx = findfirst(j -> !iszero(null_vecs[m5_const_col, j]), outside_idx)
    if e5_idx === nothing
        println("ERROR: No vector with M_5 constant!")
        return nothing
    end
    j_e5 = outside_idx[e5_idx]
    c = normalize_vec(null_vecs[:, j_e5])
    if c[m5_const_col] > 0
        c = -c
    end
    @printf("E5 base: null vector #%d, M_5(k=0) = %s\n", j_e5, c[m5_const_col])

    println("\nE5 expression:")
    for a in 1:a_max
        parts = String[]
        for k in 0:d
            val = c[bidx(k, a)]
            iszero(val) && continue
            push!(parts, k == 0 ? "$val" : k == 1 ? "$(val)·n" : "$(val)·n^$k")
        end
        isempty(parts) || println("  M_$a: ", join(parts, " + "))
    end

    # Verify prime vanishing at primes up to 500
    println("\nVerifying prime vanishing at primes in [2,500]...")
    big_primes = filter(is_prime_trial, 2:500)
    big_pm = eval_matrix(d, a_max, big_primes)
    pv = big_pm * c
    max_pv = maximum(abs.(pv))
    println("Max |E5(p)| = $max_pv", iszero(max_pv) ? " ✓" : " ✗")

    # Non-negativity at composites
    println("Checking composites 4..300...")
    all_comp = filter(n -> !is_prime_trial(n), 4:300)
    cm = eval_matrix(d, a_max, all_comp)
    cvals = cm * c
    neg_idx = findall(v -> v < 0, cvals)
    if isempty(neg_idx)
        println("✓ E5 ≥ 0 at all composites in [4,300]")
    else
        println("✗ NEGATIVE at $(length(neg_idx)) composites:")
        for i in neg_idx[1:min(5, end)]
            println("  n=$(all_comp[i]): $(cvals[i])")
        end
    end

    return c, d, a_max, bidx
end

result = find_e5()
if result !== nothing
    c, d, a_max, bidx = result
    println("\n" * "="^60)
    println("Julia E5 function:")
    println()
    println("function E5(n::Int)")
    for a in 1:a_max
        has_any = any(!iszero(c[bidx(k, a)]) for k in 0:d)
        has_any || continue
        ma = a==1 ? "M1(n)" : a==2 ? "M2(n)" : a==3 ? "M3(n)" : a==4 ? "M4(n)" : "M5(n)"
        println("    m$a = $ma")
    end
    println("    return (")
    for a in 1:a_max
        for k in 0:d
            v = c[bidx(k, a)]
            iszero(v) && continue
            n_str = k == 0 ? "" : k == 1 ? " * big(n)" : " * big(n)^$k"
            println("        + $v$n_str * m$a")
        end
    end
    println("    )")
    println("end")
end
