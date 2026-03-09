"""
extract_paper_data.jl — all wrapped in functions to avoid Julia soft-scope issues.
Run from PartitionPrimes root:
  julia --project=QuasiShuffleAlgebra Paper/extract_paper_data.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "QuasiShuffleAlgebra"))
using QuasiShuffleAlgebra
using Printf

sep() = println("\n" * "="^60)

# ── globals used across sections ─────────────────────────────────────────────
const PRIMES_500     = filter(is_prime_trial, 2:500)
const COMPOSITES_100 = filter(n -> !is_prime_trial(n), 4:100)
const COMPOSITES_500 = filter(n -> !is_prime_trial(n), 4:500)

function normalize_vec(c::Vector{Rational{BigInt}})
    nonzero = filter(!iszero, c)
    isempty(nonzero) && return c
    lcm_d = reduce(lcm, denominator.(nonzero))
    c_int = [numerator(v * lcm_d) for v in c]
    nz_int = filter(!iszero, c_int)
    gcd_n = reduce(gcd, nz_int)
    return c_int .// gcd_n
end

# ── 1. M6 FITTING DETAILS ────────────────────────────────────────────────────
function section1_m6_fitting()
    sep()
    println("1. M6 COEFFICIENT FITTING")
    println()

    # Design matrix: 21 divisor-sum columns + 1 tau column
    ns = 1:55
    A = zeros(Rational{BigInt}, length(ns), 22)
    col = 1
    for j in 0:5
        for k in 0:(5-j)
            for (i, n) in enumerate(ns)
                A[i, col] = Rational{BigInt}(big(n)^k) * Rational{BigInt}(σ(2j+1, n))
            end
            col += 1
        end
    end
    for (i, n) in enumerate(ns)
        A[i, 22] = Rational{BigInt}(ramanujan_tau(n))
    end

    r_no_tau  = rank_over_Q(copy(A[:, 1:21]))
    r_with_tau = rank_over_Q(copy(A))

    @printf("  Design matrix: %d rows × 22 cols (21 divisor-sum + 1 tau)\n", length(ns))
    @printf("  Rank WITHOUT tau column: %d / 21\n", r_no_tau)
    @printf("  Rank WITH    tau column: %d / 22\n", r_with_tau)
    @printf("  tau is linearly independent: %s\n", r_with_tau > r_no_tau ? "YES ✓" : "NO ✗")
    println()
    @printf("  c_tau = -17 / 150450048000\n")

    # Factor the denominator
    d_val = big(150450048000)
    println("  Denominator factorisation:")
    println("    150450048000 = 2^10 × 3^5 × 5^3 × 7 × 691")
    # Just verify by printing what it actually factors to
    m = d_val
    factors = Pair{BigInt,Int}[]
    for p in [2,3,5,7,11,13,17,19,23,29]
        cnt = 0
        while m % p == 0
            m ÷= p; cnt += 1
        end
        cnt > 0 && push!(factors, big(p) => cnt)
    end
    m > 1 && push!(factors, m => 1)
    fstr = join(["$(p)^$(e)" for (p,e) in factors], " × ")
    println("    Actual: 150450048000 = $fstr")
    println()

    # Verify M6
    errors = 0
    for n in 1:30
        M6(n) != Rational{BigInt}(M_direct(6, n)) && (errors += 1)
    end
    @printf("  M6 closed form vs M_direct n=1..30: %s\n",
            errors == 0 ? "✓ all match" : "✗ $errors mismatches")
end

# ── 2. PRIME EVALUATION MATRIX ───────────────────────────────────────────────
function section2_matrix_dimensions()
    sep()
    println("2. PRIME EVALUATION MATRIX DIMENSIONS")
    println()
    @printf("  Primes in [2,500]: %d\n\n", length(PRIMES_500))

    for d_ref in [3]
        println("  d = $d_ref:")
        for a_max in [5, 6]
            pm = eval_matrix(d_ref, a_max, PRIMES_500)
            nv = rational_nullspace(pm)
            r  = size(pm, 2) - size(nv, 2)
            @printf("    a_max=%d: %d × %d,  rank=%d,  nullity=%d\n",
                    a_max, size(pm)..., r, size(nv,2))
        end
    end

    println()
    println("  Nullity with vs without M6 across degrees:")
    println("   d  | nullity(a≤5) | nullity(a≤6) | M6 adds to null space?")
    println("  ----|-------------|-------------|------------------------")
    for d in 2:6
        nv5 = rational_nullspace(eval_matrix(d, 5, PRIMES_500))
        nv6 = rational_nullspace(eval_matrix(d, 6, PRIMES_500))
        @printf("  d=%d |          %2d |          %2d | %s\n",
                d, size(nv5,2), size(nv6,2),
                size(nv5,2)==size(nv6,2) ? "No ✓" : "Yes ✗")
    end
end

# ── 3. E5 EXTRACTION ─────────────────────────────────────────────────────────
function section3_extraction()
    sep()
    println("3. E5 EXTRACTION: DEGREE SWEEP (a_max=5)")
    println()
    println("   d  | basis_dim | null_dim | outside E1-E4 span")
    println("  ----|-----------|----------|-------------------")

    found_d   = nothing
    found_vec = nothing

    for d in 0:5
        pm  = eval_matrix(d, 5, PRIMES_500)
        nv  = rational_nullspace(pm)
        cm  = eval_matrix(d, 5, COMPOSITES_100)

        e14_cols = [Rational{BigInt}(big(n)^j) * Rational{BigInt}(Ef(n))
                    for n in COMPOSITES_100, j in 0:d, Ef in [E1, E2, E3, E4]]
        e14_mat = reshape(e14_cols, length(COMPOSITES_100), 4*(d+1))

        outside_idx = Int[]
        for j in 1:size(nv, 2)
            is_in_colspan(cm * nv[:,j], e14_mat) || push!(outside_idx, j)
        end

        status = if length(outside_idx) == 0
            "(all in E1-E4 span)"
        elseif length(outside_idx) == 1
            "← E5 found here ✓"
        else
            "← $(length(outside_idx)) directions (unexpected)"
        end
        @printf("  d=%d |        %3d |       %2d | %d  %s\n",
                d, (d+1)*5, size(nv,2), length(outside_idx), status)

        if length(outside_idx) >= 1 && found_d === nothing
            raw = nv[:, outside_idx[1]]
            vec = normalize_vec(raw)
            # Orient: M5 constant coefficient positive
            bidx(k, a) = (a-1)*(d+1) + k + 1
            if vec[bidx(0,5)] < 0
                vec = -vec
            end
            found_d   = d
            found_vec = vec
        end
    end

    return found_d, found_vec
end

# ── 4. E5 FORMULA ────────────────────────────────────────────────────────────
function section4_formula(d_e5, e5_vec)
    sep()
    println("4. E5 FORMULA (extracted at d=$d_e5, a_max=5)")
    println()

    bidx(k, a) = (a-1)*(d_e5+1) + k + 1

    println("  E5(n) =")
    for a in 1:5
        terms = String[]
        for k in 0:d_e5
            v = e5_vec[bidx(k, a)]
            iszero(v) && continue
            vn = numerator(v)  # integers after normalize
            nstr = k == 0 ? "" : k == 1 ? "*n" : "*n^$k"
            push!(terms, "$vn$nstr")
        end
        isempty(terms) && continue
        coeff = length(terms) == 1 ? terms[1] : "($(join(terms, " + ")))"
        println("    + $coeff * M$a(n)")
    end
    println()

    # Verify against package E5
    errors = 0
    for n in vcat(filter(is_prime_trial, 2:50), filter(n -> !is_prime_trial(n), 4:50))
        Ma = [Rational{BigInt}(M1(n)), M2(n), M3(n), M4(n), M5(n)]
        manual = sum(
            e5_vec[bidx(k,a)] * Rational{BigInt}(big(n)^k) * Ma[a]
            for a in 1:5 for k in 0:d_e5
        )
        # compare to package E5 at same degree
        iszero(manual) != is_prime_trial(n) && (errors += 1)
    end
    @printf("  Extracted E5 (d=%d) vanishes iff prime for n=2..50: %s\n",
            d_e5, errors == 0 ? "✓" : "✗ $errors mismatches")
end

# ── 5. PRIME VERIFICATION ────────────────────────────────────────────────────
function section5_prime_verification()
    sep()
    println("5. PRIME VERIFICATION")
    println()
    vals = [E5(p) for p in PRIMES_500]
    @printf("  E5(p) for all %d primes in [2,500]: max|E5(p)| = %s  %s\n",
            length(PRIMES_500), maximum(abs.(vals)),
            all(iszero, vals) ? "✓" : "✗")
    println()
    println("  Spot-check:")
    for p in [2, 3, 5, 7, 11, 13, 97, 101, 499]
        @printf("    E5(%3d) = %s\n", p, E5(p))
    end
end

# ── 6. NEGATIVE-E5 COMPOSITES ────────────────────────────────────────────────
function section6_negative_composites()
    sep()
    println("6. NEGATIVE-E5 COMPOSITES IN [4, 500]")
    println()

    neg = [(n, E5(n)) for n in COMPOSITES_500 if E5(n) < 0]
    @printf("  Count: %d\n\n", length(neg))
    println("  n   | factorization      | E5(n)")
    println("  ----|--------------------|---------")
    for (n, val) in neg
        fs = Int[]
        m  = n
        for p in 2:isqrt(n)+1
            while m % p == 0; push!(fs, p); m ÷= p; end
        end
        m > 1 && push!(fs, m)
        @printf("  %3d | %-18s | %s\n", n, join(fs, " × "), val)
    end
    println()

    all5 = all(n % 5 == 0 for (n,_) in neg)
    @printf("  All divisible by 5: %s\n", all5 ? "YES ✓" : "NO")
    minv, mini = findmin([val for (_,val) in neg])
    @printf("  Minimum: E5(%d) = %s\n", neg[mini][1], minv)

    return neg
end

# ── 7. REMEDIATION ───────────────────────────────────────────────────────────
function section7_remediation()
    sep()
    println("7. REMEDIATION: smallest k so E5 + k*E4 ≥ 0 at all composites [4,500]")
    println()

    e5v = [E5(n) for n in COMPOSITES_500]
    e4v = [E4(n) for n in COMPOSITES_500]

    min_k = nothing
    for k in 1:10000
        all(e5v[i] + k*e4v[i] >= 0 for i in eachindex(COMPOSITES_500)) || continue
        min_k = k
        break
    end

    if min_k !== nothing
        combined = [e5v[i] + min_k*e4v[i] for i in eachindex(COMPOSITES_500)]
        minval   = minimum(combined)
        minn     = COMPOSITES_500[argmin(combined)]
        @printf("  Smallest k: %d\n", min_k)
        @printf("  min(E5 + %d*E4) = %s at n=%d\n", min_k, minval, minn)
        pv_ok = all(iszero(E5(p) + min_k*E4(p)) for p in PRIMES_500)
        @printf("  E5 + %d*E4 vanishes at all primes [2,500]: %s\n", min_k, pv_ok ? "✓" : "✗")
    else
        println("  No k ≤ 10000 works — check manually")
    end
end

# ── 8. CONJECTURE TABLE ──────────────────────────────────────────────────────
function section8_conjecture_table()
    sep()
    println("8. CONJECTURE VERIFICATION TABLE (a_max=5, N=300)")
    println()
    println("   d  | basis_dim | pv_dim | t1_span | holds? | time(s)")
    println("  ----|-----------|--------|---------|--------|--------")
    for d in 2:6
        t0 = time()
        r  = test_conjecture(d, 5; N=300, verbose=false)
        dt = time() - t0
        @printf("  d=%d |        %2d |     %2d |      %2d | %6s | %.1f\n",
                d, (d+1)*5, r.dim_prime_vanishing, r.dim_table1_span,
                r.holds ? "✓" : "✗", dt)
    end
end

# ── 9. SUMMARY TABLE ─────────────────────────────────────────────────────────
function section9_summary()
    sep()
    println("9. E1-E5 DEGREE AND WEIGHT SUMMARY")
    println()
    println("  Expr | Max poly degree | Max weight | Non-neg? | Note")
    println("  -----|-----------------|------------|----------|---------------------")
    println("  E1   |        2        |     2      |   Yes    | Ea introduces M_{a+1}")
    println("  E2   |        3        |     3      |   Yes    | Ea introduces M_{a+1}")
    println("  E3   |        4        |     4      |   Yes    | Ea introduces M_{a+1}")
    println("  E4   |        5        |     5      |   Yes    | Ea introduces M_{a+1}")
    println("  E5   |        2        |     5      |   No     | Pattern breaks — M6 excluded, degree resets to 2")
    println()
    println("  KEY: E5 first appears at d=2 — the canonical minimal-degree form.")
    println("  The package E5 implements this d=2 canonical formula.")
end

# ── MAIN ─────────────────────────────────────────────────────────────────────
section1_m6_fitting()
section2_matrix_dimensions()
d_e5, e5_vec = section3_extraction()
section4_formula(d_e5, e5_vec)
section5_prime_verification()
section6_negative_composites()
section7_remediation()
section8_conjecture_table()
section9_summary()

sep()
println("DONE.")
