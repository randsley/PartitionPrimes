"""
extract_e6_e7.jl — search for prime-vanishing expressions E6, E7
in the extended {M1, M2, M3, M4, M5, M6, M7} basis.

Key discovery: M7(n) contains both τ(n) and n·τ(n) terms because the depth-1
component of U_7(q) (as a quasimodular form) lives in M_12 = span{E_12, Δ},
and the Δ piece contributes via E_2·Δ. Individually, M6 and M7 columns are
pivots in the prime evaluation matrix. However, TOGETHER, they admit linear
combinations that cancel the τ terms (one equation cancels τ(p), another
cancels p·τ(p)), leaving a polynomial in p — potentially giving new null-space
directions that correspond to E6 and E7.

Run from PartitionPrimes root:
  julia --project=QuasiShuffleAlgebra Paper/extract_e6_e7.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "QuasiShuffleAlgebra"))
using QuasiShuffleAlgebra
using Printf

sep() = println("\n" * "="^60)

# Use all a up to 7 (include both M6 and M7)
const A_LIST      = [1, 2, 3, 4, 5, 6, 7]
const PRIMES_500  = filter(is_prime_trial, 2:500)
const COMPS_100   = filter(n -> !is_prime_trial(n), 4:100)
const COMPS_500   = filter(n -> !is_prime_trial(n), 4:500)

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
    ncols = (d + 1) * (5 + length(extra_fns))
    S     = zeros(Rational{BigInt}, rows, ncols)
    col   = 1
    for j in 0:d
        for Ef in [E1, E2, E3, E4, E5]
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

# ── Step 1: Verify M7 ─────────────────────────────────────────────────────────
function step1_verify_M7()
    sep()
    println("STEP 1: Verify M7 closed form (triggers lazy fit on first call)")
    println()
    print("  Fitting ... ")
    t0 = time()
    M7(1)
    @printf("done in %.2f s\n", time() - t0)

    local errs = 0
    for n in 1:30
        M7(n) != Rational{BigInt}(M_direct(7, n)) && (errs += 1)
    end
    @printf("  M7 vs M_direct n=1..30: %s\n", errs == 0 ? "✓" : "✗ $errs mismatches")
    @printf("  M7(28) = %s  (expected 1)\n", M7(28))

    println()
    println("  Note: M7 contains BOTH τ(n) and n·τ(n) terms.")
    println("  This is because the depth-1 piece of U_7(q) lives in M_12 = span{E_12, Δ},")
    println("  and E_2·Δ contributes τ-type terms to Fourier coefficients.")
end

# ── Step 2: Nullity comparison ────────────────────────────────────────────────
function step2_nullity_table()
    sep()
    println("STEP 2: Nullity across basis choices")
    println()
    println("   d  | null(M1-5) | null(+M6) | null(+M7) | null(+M6+M7) | new from M6+M7?")
    println("  ----|-----------|-----------|-----------|-------------|----------------")
    for d in 2:6
        mat5  = eval_matrix_alist(d, [1,2,3,4,5],   PRIMES_500)
        mat56 = eval_matrix_alist(d, [1,2,3,4,5,6], PRIMES_500)
        mat57 = eval_matrix_alist(d, [1,2,3,4,5,7], PRIMES_500)
        mat567= eval_matrix_alist(d, A_LIST,         PRIMES_500)
        nv5   = size(rational_nullspace(mat5),   2)
        nv56  = size(rational_nullspace(mat56),  2)
        nv57  = size(rational_nullspace(mat57),  2)
        nv567 = size(rational_nullspace(mat567), 2)
        gain  = nv567 - nv5
        @printf("  d=%d |    %3d    |    %3d    |    %3d    |     %3d     | %s\n",
                d, nv5, nv56, nv57, nv567,
                gain > 0 ? "+$gain ← NEW!" : "0")
    end
end

# ── Step 3: Extract new expressions ──────────────────────────────────────────
function step3_extract()
    sep()
    println("STEP 3: Degree sweep — finding directions outside Q[n]·span(E1–E5)")
    println()
    println("   d  | basis_dim | null_dim | new dirs outside E1–E5 span")
    println("  ----|-----------|----------|-----------------------------")

    found_fns   = Function[]
    found_exprs = NamedTuple[]

    for d in 0:8
        basis = [(k, a) for a in A_LIST for k in 0:d]
        pm    = eval_matrix_alist(d, A_LIST, PRIMES_500)
        nv    = rational_nullspace(pm)

        if size(nv, 2) == 0
            @printf("  d=%d |     %3d   |    0     |\n", d, length(basis))
            continue
        end

        cm   = eval_matrix_alist(d, A_LIST, COMPS_100)
        span = build_known_span(d, found_fns)

        local outside_idx = Int[]
        for j in 1:size(nv, 2)
            is_in_colspan(cm * nv[:, j], span) || push!(outside_idx, j)
        end

        @printf("  d=%d |     %3d   |   %3d    | %d  %s\n",
                d, length(basis), size(nv, 2), length(outside_idx),
                length(outside_idx) > 0 ? "← NEW!" : "")

        for idx in outside_idx
            raw = nv[:, idx]
            vec = normalize_vec(raw)

            # Orient: prefer positive M7 constant, then M6 constant, then first nonzero
            orient_idx = findfirst(((k,a),) -> a == 7 && k == 0, basis)
            orient_idx === nothing && (orient_idx = findfirst(((k,a),) -> a == 6 && k == 0, basis))
            if orient_idx !== nothing && !iszero(vec[orient_idx]) && vec[orient_idx] < 0
                vec = -vec
            elseif all(iszero, vec)
                # skip
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
            push!(found_exprs, (name  = "E$(5 + length(found_exprs) + 1)",
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
        println("        {E1,...,E5} spans the full prime-vanishing subspace")
        println("        for a_max ≤ 7, within these bounds.")
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
            push!(get!(by_a, a, Tuple{Int,Rational{BigInt}}[]), (k, c))
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
        local prime_ok = all(iszero(fn(p)) for p in PRIMES_500)
        local n_pos    = count(fn(n) > 0 for n in COMPS_500)
        local n_neg    = count(fn(n) < 0 for n in COMPS_500)
        local n_zero   = count(iszero(fn(n)) for n in COMPS_500)
        @printf("  %s vanishes at all %d primes in [2,500]: %s\n",
                expr.name, length(PRIMES_500), prime_ok ? "✓" : "✗")
        @printf("  %s composites in [4,500]:  +%d  -%d  zero:%d\n",
                expr.name, n_pos, n_neg, n_zero)
        local pos_cs = filter(n -> fn(n) > 0, COMPS_500)
        if !isempty(pos_cs)
            @printf("  All positive composites odd: %s\n", all(isodd, pos_cs) ? "Yes" : "No")
        end
        println()
    end
end

# ── Step 6: Extended conjecture check ─────────────────────────────────────────
function step6_conjecture(found_exprs, found_fns)
    isempty(found_exprs) && return
    sep()
    names = join(getfield.(found_exprs, :name), ",")
    println("STEP 6: Conjecture check — does {E1,...,E5,$names} span")
    println("        the full prime-vanishing subspace? (a_list=$(A_LIST), N=300)")
    println()
    println("   d  | pv_dim | span_dim | holds?")
    println("  ----|--------|----------|-------")
    N300 = filter(is_prime_trial, 2:2000)[1:300]
    for d in 2:6
        pm   = eval_matrix_alist(d, A_LIST, N300)
        nv   = rational_nullspace(pm)
        cm   = eval_matrix_alist(d, A_LIST, COMPS_100)
        span = build_known_span(d, found_fns)
        local outside = 0
        for j in 1:size(nv, 2)
            is_in_colspan(cm * nv[:, j], span) || (outside += 1)
        end
        @printf("  d=%d |   %3d  |   %3d    | %s\n",
                d, size(nv, 2), rank_over_Q(span),
                outside == 0 ? "✓" : "✗ ($outside outside)")
    end
end

# ── Main ──────────────────────────────────────────────────────────────────────
step1_verify_M7()
step2_nullity_table()
found_exprs, found_fns = step3_extract()
step4_formulas(found_exprs)
step5_verify(found_exprs, found_fns)
step6_conjecture(found_exprs, found_fns)
sep()
println("DONE.")
