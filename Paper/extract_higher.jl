"""
extract_higher.jl — search for prime-vanishing expressions beyond E6,
extending to {M1,...,M12} (weights up to 24).

Key motivation: S_24 has dimension 2 (spanned by Δ·E₁₂ and Δ²), making
M12 the first weight where two independent cusp forms coexist — the same
structural condition that generated E6 at weight 12/14.

Performance: M_direct(a, n) fitting cost scales sharply with a.
  M_direct(10,n) for n=1..120  ~30-90 s
  M_direct(11,n) for n=1..140  ~5-30 min
  M_direct(12,n) for n=1..160  ~30-120 min

Run from PartitionPrimes root:
  julia --project=QuasiShuffleAlgebra Paper/extract_higher.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "QuasiShuffleAlgebra"))
using QuasiShuffleAlgebra
using Printf

sep() = println("\n" * "="^60)

const PRIMES_500 = filter(is_prime_trial, 2:500)
const COMPS_100  = filter(n -> !is_prime_trial(n), 4:100)
const COMPS_500  = filter(n -> !is_prime_trial(n), 4:500)

# Known expressions (E1–E6); update if more are found
const KNOWN_E_FUNS = Function[E1, E2, E3, E4, E5, E6]
const N_KNOWN      = length(KNOWN_E_FUNS)

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

function normalize_vec(c::Vector{Rational{BigInt}})
    nonzero = filter(!iszero, c)
    isempty(nonzero) && return c
    lcm_d  = reduce(lcm, denominator.(nonzero))
    c_int  = [numerator(v * lcm_d) for v in c]
    nz_int = filter(!iszero, c_int)
    gcd_n  = reduce(gcd, nz_int)
    return c_int .// gcd_n
end

function build_known_span(d::Int, extra_fns::Vector)
    rows  = length(COMPS_100)
    ncols = (d + 1) * (N_KNOWN + length(extra_fns))
    S     = zeros(Rational{BigInt}, rows, ncols)
    col   = 1
    for j in 0:d
        for Ef in [KNOWN_E_FUNS; extra_fns]
            for (r, n) in enumerate(COMPS_100)
                S[r, col] = Rational{BigInt}(big(n)^j) * Rational{BigInt}(Ef(n))
            end
            col += 1
        end
    end
    return S
end

# ── Step 1: Fit and verify M10–M12 ────────────────────────────────────────────
function step1_verify()
    sep()
    println("STEP 1: Fit and verify M10, M11, M12")
    println()

    for (a, Ma, label) in [(10, M10, "n=1..120"), (11, M11, "n=1..140"), (12, M12, "n=1..160")]
        print("  Fitting M$a ($label) ... ")
        flush(stdout)
        t0 = time()
        Ma(1)  # triggers lazy fit
        @printf("done in %.1f s\n", time() - t0)
        errs = sum(Ma(n) != Rational{BigInt}(M_direct(a, n)) for n in 1:15)
        @printf("  M%d vs M_direct n=1..15: %s\n", a, errs == 0 ? "✓" : "✗ $errs mismatches")
    end
end

# ── Step 2: Nullity table ─────────────────────────────────────────────────────
function step2_nullity_table()
    sep()
    println("STEP 2: Nullity table — incremental gains from M10, M11, M12")
    println()
    println("   d  | null(1-9) | null(+M10) | null(+M11) | null(+M12) | gains")
    println("  ----|-----------|------------|------------|------------|------")
    for d in 2:6
        mat9   = eval_matrix_alist(d, collect(1:9),  PRIMES_500)
        mat10  = eval_matrix_alist(d, collect(1:10), PRIMES_500)
        mat11  = eval_matrix_alist(d, collect(1:11), PRIMES_500)
        mat12  = eval_matrix_alist(d, collect(1:12), PRIMES_500)
        nv9    = size(rational_nullspace(mat9),   2)
        nv10   = size(rational_nullspace(mat10),  2)
        nv11   = size(rational_nullspace(mat11),  2)
        nv12   = size(rational_nullspace(mat12),  2)
        gains  = filter(!iszero, [nv10-nv9, nv11-nv9, nv12-nv9])
        @printf("  d=%d |    %3d    |     %3d    |     %3d    |     %3d    | %s\n",
                d, nv9, nv10, nv11, nv12,
                isempty(gains) ? "" : join(["+$(g)" for g in [nv10-nv9, nv11-nv9, nv12-nv9]], " ") * "  ← NEW!" )
    end
end

# ── Step 3: Extract new expressions ───────────────────────────────────────────
function step3_extract()
    sep()
    println("STEP 3: Degree sweep d=0..8 with a_max=12")
    println("        Finding directions outside Q[n]·span(E1–E6)")
    println()
    println("   d  | basis_dim | null_dim | new dirs")
    println("  ----|-----------|----------|----------")

    found_fns   = Function[]
    found_exprs = NamedTuple[]
    a_list      = collect(1:12)

    for d in 0:8
        basis = [(k, a) for a in a_list for k in 0:d]
        pm    = eval_matrix_alist(d, a_list, PRIMES_500)
        nv    = rational_nullspace(pm)

        if size(nv, 2) == 0
            @printf("  d=%d |     %3d   |    0     |\n", d, length(basis))
            continue
        end

        cm   = eval_matrix_alist(d, a_list, COMPS_100)
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
            for a_try in reverse(a_list)
                oi = findfirst(((k, a),) -> a == a_try && k == 0, basis)
                if oi !== nothing && !iszero(vec[oi])
                    vec[oi] < 0 && (vec = -vec)
                    break
                end
            end
            if !all(iszero, vec) && first(filter(!iszero, vec)) < 0
                vec = -vec
            end

            basis_c = copy(basis); vec_c = copy(vec)
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
            push!(found_exprs, (name  = "E$(N_KNOWN + length(found_exprs) + 1)",
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
        println("RESULT: No new prime-vanishing directions found (d ≤ 8, a_max ≤ 12).")
        println("        {E1,...,E6} spans the full prime-vanishing subspace")
        println("        within these bounds.")
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

# ── Step 6: Span check ────────────────────────────────────────────────────────
function step6_span(found_exprs, found_fns)
    sep()
    names = isempty(found_exprs) ? "" : "," * join(getfield.(found_exprs, :name), ",")
    println("STEP 6: Does {E1,...,E6$names} span the full pv subspace?")
    println("        (a_max=12, N=300 primes)")
    println()
    println("   d  | pv_dim | span_dim | holds?")
    println("  ----|--------|----------|-------")
    N300 = filter(is_prime_trial, 2:2000)[1:300]
    for d in 2:6
        pm      = eval_matrix_alist(d, collect(1:12), N300)
        nv      = rational_nullspace(pm)
        cm      = eval_matrix_alist(d, collect(1:12), COMPS_100)
        span    = build_known_span(d, found_fns)
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
t_start = time()
step1_verify()
step2_nullity_table()
found_exprs, found_fns = step3_extract()
step4_formulas(found_exprs)
step5_verify(found_exprs, found_fns)
step6_span(found_exprs, found_fns)
sep()
@printf("DONE. Total elapsed: %.1f s\n", time() - t_start)
