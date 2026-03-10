"""
extract_e8_e9.jl — search for prime-vanishing expressions E7, E8, E9
in the extended {M1,...,M9} basis.

Key structure:
  S_16 = span{Δ·E₄}   → cusp_form_16(n) in M8
  S_18 = span{Δ·E₆}   → cusp_form_18(n) in M9
  S_20 has dim 2       → M10 would introduce two new cusp form contributions

Each new cusp-form dimension opens a potential τ-cancellation channel.
Within {M1,...,M9} we ask: does {E1,...,E6} still span the full
prime-vanishing subspace, or do M8/M9 admit new directions?

Performance note: M_direct(8, n) for n up to 90 takes ~5-20 min.
                  M_direct(9, n) for n up to 110 takes ~30-90 min.
Run from PartitionPrimes root:
  julia --project=QuasiShuffleAlgebra Paper/extract_e8_e9.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "QuasiShuffleAlgebra"))
using QuasiShuffleAlgebra
using Printf

sep() = println("\n" * "="^60)

const A_LIST_8  = [1, 2, 3, 4, 5, 6, 7, 8]
const A_LIST_9  = [1, 2, 3, 4, 5, 6, 7, 8, 9]
const PRIMES_500 = filter(is_prime_trial, 2:500)
const COMPS_100  = filter(n -> !is_prime_trial(n), 4:100)
const COMPS_500  = filter(n -> !is_prime_trial(n), 4:500)

# ── Helpers ───────────────────────────────────────────────────────────────────

function normalize_vec(c::Vector{Rational{BigInt}})
    nonzero = filter(!iszero, c)
    isempty(nonzero) && return c
    lcm_d  = reduce(lcm, denominator.(nonzero))
    c_int  = [numerator(v * lcm_d) for v in c]
    nz_int = filter(!iszero, c_int)
    gcd_n  = reduce(gcd, nz_int)
    return c_int .// gcd_n
end

function eval_matrix_alist(d::Int, a_list::Vector{Int},
                           ns::AbstractVector{Int})::Matrix{Rational{BigInt}}
    basis = [(k, a) for a in a_list for k in 0:d]
    mat   = zeros(Rational{BigInt}, length(ns), length(basis))
    for (i, n) in enumerate(ns)
        cache = Dict{Int, Rational{BigInt}}()
        for (j, (k, a)) in enumerate(basis)
            ma = get!(cache, a) do; eval_Ma(a, n) end
            mat[i, j] = big(n)^k * ma
        end
    end
    return mat
end

function build_known_span(d::Int, extra_fns::Vector)
    rows  = length(COMPS_100)
    ncols = (d + 1) * (6 + length(extra_fns))
    S     = zeros(Rational{BigInt}, rows, ncols)
    col   = 1
    for j in 0:d
        for Ef in [E1, E2, E3, E4, E5, E6]
            for (r, n) in enumerate(COMPS_100)
                S[r, col] = Rational{BigInt}(big(n)^j) * Rational{BigInt}(Ef(n))
            end
            col += 1
        end
        for fn in extra_fns
            for (r, n) in enumerate(COMPS_100)
                S[r, col] = Rational{BigInt}(big(n)^j) * fn(n)
            end
            col += 1
        end
    end
    return S
end

# ── Step 1: Verify M8 and M9 ──────────────────────────────────────────────────
function step1_verify()
    sep()
    println("STEP 1: Fit and verify M8 and M9")
    println()

    print("  Fitting M8 (M_direct(8,n) for n=1..90) ... ")
    t0 = time()
    M8(1)
    @printf("done in %.1f s\n", time() - t0)
    errs = sum(M8(n) != Rational{BigInt}(M_direct(8, n)) for n in 1:20)
    @printf("  M8 vs M_direct n=1..20: %s\n", errs == 0 ? "✓" : "✗ $errs mismatches")

    print("  Fitting M9 (M_direct(9,n) for n=1..110) ... ")
    t0 = time()
    M9(1)
    @printf("done in %.1f s\n", time() - t0)
    errs = sum(M9(n) != Rational{BigInt}(M_direct(9, n)) for n in 1:20)
    @printf("  M9 vs M_direct n=1..20: %s\n", errs == 0 ? "✓" : "✗ $errs mismatches")
end

# ── Step 2: Nullity table ─────────────────────────────────────────────────────
function step2_nullity_table()
    sep()
    println("STEP 2: Nullity table — a_max = 7 vs 8 vs 9")
    println()
    println("   d  | null(1-7) | null(+M8) | null(+M9) | gain vs base")
    println("  ----|-----------|-----------|-----------|-------------")
    for d in 2:6
        mat7  = eval_matrix_alist(d, A_LIST_8[1:7], PRIMES_500)
        mat8  = eval_matrix_alist(d, A_LIST_8,      PRIMES_500)
        mat9  = eval_matrix_alist(d, A_LIST_9,      PRIMES_500)
        nv7   = size(rational_nullspace(mat7),  2)
        nv8   = size(rational_nullspace(mat8),  2)
        nv9   = size(rational_nullspace(mat9),  2)
        gain8 = nv8 - nv7
        gain9 = nv9 - nv7
        @printf("  d=%d |    %3d    |    %3d    |    %3d    | +M8:%d  +M9:%d  %s\n",
                d, nv7, nv8, nv9, gain8, gain9,
                (gain8 > 0 || gain9 > 0) ? "← NEW!" : "")
    end
end

# ── Step 3: Extract new expressions ───────────────────────────────────────────
function step3_extract()
    sep()
    println("STEP 3: Degree sweep — finding directions outside Q[n]·span(E1–E6)")
    println()
    println("   d  | basis_dim | null_dim | new dirs outside E1–E6 span")
    println("  ----|-----------|----------|-----------------------------")

    found_fns   = Function[]
    found_exprs = NamedTuple[]

    for d in 0:8
        basis = [(k, a) for a in A_LIST_9 for k in 0:d]
        pm    = eval_matrix_alist(d, A_LIST_9, PRIMES_500)
        nv    = rational_nullspace(pm)

        if size(nv, 2) == 0
            @printf("  d=%d |     %3d   |    0     |\n", d, length(basis))
            continue
        end

        cm   = eval_matrix_alist(d, A_LIST_9, COMPS_100)
        span = build_known_span(d, found_fns)

        outside_idx = Int[]
        for j in 1:size(nv, 2)
            is_in_colspan(cm * nv[:, j], span) || push!(outside_idx, j)
        end

        @printf("  d=%d |     %3d   |   %3d    | %d  %s\n",
                d, length(basis), size(nv, 2), length(outside_idx),
                length(outside_idx) > 0 ? "← NEW!" : "")

        for idx in outside_idx
            raw = nv[:, idx]
            vec = normalize_vec(raw)

            # Orient: prefer positive leading high-a constant
            for a_try in [9, 8, 7]
                oi = findfirst(((k, a),) -> a == a_try && k == 0, basis)
                if oi !== nothing && !iszero(vec[oi])
                    vec[oi] < 0 && (vec = -vec)
                    break
                end
            end
            if all(iszero, vec)
            elseif first(filter(!iszero, vec)) < 0
                vec = -vec
            end

            basis_c = copy(basis)
            vec_c   = copy(vec)
            eval_fn = function(n::Int)
                result = Rational{BigInt}(0)
                cache  = Dict{Int, Rational{BigInt}}()
                for (j2, (k, a)) in enumerate(basis_c)
                    ma = get!(cache, a) do; eval_Ma(a, n) end
                    result += vec_c[j2] * big(n)^k * ma
                end
                return result
            end

            push!(found_fns,   eval_fn)
            push!(found_exprs, (name  = "E$(6 + length(found_exprs) + 1)",
                                d     = d,
                                basis = basis_c,
                                vec   = vec_c))
        end
    end

    return found_exprs, found_fns
end

# ── Step 4: Print formulas ────────────────────────────────────────────────────
function step4_formulas(found_exprs)
    sep()
    if isempty(found_exprs)
        println("RESULT: No new prime-vanishing directions found for d ≤ 8.")
        println("        {E1,...,E6} spans the full prime-vanishing subspace")
        println("        for a_max ≤ 9, within these bounds.")
        return
    end

    println("STEP 4: Formulas for new prime-vanishing expressions")
    println()
    for expr in found_exprs
        println("  $(expr.name)(n)  [degree d=$(expr.d)]")
        println()
        by_a = Dict{Int, Vector{Tuple{Int, Rational{BigInt}}}}()
        for ((k, a), c) in zip(expr.basis, expr.vec)
            iszero(c) && continue
            push!(get!(by_a, a, Tuple{Int, Rational{BigInt}}[]), (k, c))
        end
        for a in sort(collect(keys(by_a)))
            terms = sort(by_a[a]; by = first)
            parts = String[]
            for (k, c) in terms
                cn = numerator(c)
                push!(parts, k == 0 ? "$cn" : k == 1 ? "$(cn)n" : "$(cn)n^$k")
            end
            poly = length(parts) == 1 ? parts[1] : "($(join(parts, " + ")))"
            println("    + $poly · M$a(n)")
        end
        println()
    end
end

# ── Step 5: Verify ────────────────────────────────────────────────────────────
function step5_verify(found_exprs, found_fns)
    isempty(found_exprs) && return
    sep()
    println("STEP 5: Verification")
    println()
    for (expr, fn) in zip(found_exprs, found_fns)
        prime_ok = all(iszero(fn(p)) for p in PRIMES_500)
        n_pos    = count(fn(n) > 0 for n in COMPS_500)
        n_neg    = count(fn(n) < 0 for n in COMPS_500)
        n_zero   = count(iszero(fn(n)) for n in COMPS_500)
        @printf("  %s vanishes at all %d primes in [2,500]: %s\n",
                expr.name, length(PRIMES_500), prime_ok ? "✓" : "✗")
        @printf("  %s composites in [4,500]:  +%d  -%d  zero:%d\n",
                expr.name, n_pos, n_neg, n_zero)
        println()
    end
end

# ── Step 6: Conjecture check ──────────────────────────────────────────────────
function step6_conjecture(found_exprs, found_fns)
    sep()
    names = isempty(found_exprs) ? "" : "," * join(getfield.(found_exprs, :name), ",")
    println("STEP 6: Does {E1,...,E6$names} span the full pv subspace?")
    println("        (a_list=1..9, N=300 primes)")
    println()
    println("   d  | pv_dim | span_dim | holds?")
    println("  ----|--------|----------|-------")
    N300 = filter(is_prime_trial, 2:2000)[1:300]
    for d in 2:6
        pm   = eval_matrix_alist(d, A_LIST_9, N300)
        nv   = rational_nullspace(pm)
        cm   = eval_matrix_alist(d, A_LIST_9, COMPS_100)
        span = build_known_span(d, found_fns)
        outside = 0
        for j in 1:size(nv, 2)
            is_in_colspan(cm * nv[:, j], span) || (outside += 1)
        end
        @printf("  d=%d |   %3d  |   %3d    | %s\n",
                d, size(nv, 2), rank_over_Q(span),
                outside == 0 ? "✓" : "✗ ($outside outside)")
    end
end

# ── Main ──────────────────────────────────────────────────────────────────────
step1_verify()
step2_nullity_table()
found_exprs, found_fns = step3_extract()
step4_formulas(found_exprs)
step5_verify(found_exprs, found_fns)
step6_conjecture(found_exprs, found_fns)
sep()
println("DONE.")
