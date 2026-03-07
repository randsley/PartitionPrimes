using QuasiShuffleAlgebra
using Printf

println("=" ^ 65)
println("Quasi-Shuffle Algebra for Partition-Theoretic Prime Detection")
println("Craig, van Ittersum & Ono — arXiv:2405.06451v2 (2024)")
println("=" ^ 65)

# ---------------------------------------------------------------------------
println("\n── MacMahon function values ──")
println("  M_a(n) = Σ m₁m₂⋯mₐ over strict a-part partitions")
println()
@printf "  %4s  %6s  %8s  %10s\n" "n" "M₁(n)" "M₂(n)" "M₃(n)"
println("  " * "-"^34)
for n in [2, 3, 4, 5, 6, 7, 8, 10, 12, 13]
    @printf "  %4d  %6d  %8s  %10s\n" n M1(n) M2(n) M3(n)
end

# ---------------------------------------------------------------------------
println("\n── Prime-detecting expression E₁(n) = (n²-3n+2)M₁(n) - 8M₂(n) ──")
println("  Zero iff n is prime (Theorem 1.1)")
println()
for n in 2:25
    val = E1(n)
    marker = iszero(val) ? "  ← prime" : ""
    @printf "  n=%2d  E₁ = %s%s\n" n val marker
end

# ---------------------------------------------------------------------------
println("\n── All four prime-detecting expressions at n=12 (composite) ──")
println("  Each should be strictly positive at composites")
println()
n = 12
@printf "  E₁(12) = %s\n" E1(n)
@printf "  E₂(12) = %s\n" E2(n)
@printf "  E₃(12) = %s\n" E3(n)
@printf "  E₄(12) = %s\n" E4(n)
println()
println("  At a prime (n=13):")
n = 13
@printf "  E₁(13) = %s\n" E1(n)
@printf "  E₂(13) = %s\n" E2(n)
@printf "  E₃(13) = %s\n" E3(n)
@printf "  E₄(13) = %s\n" E4(n)

# ---------------------------------------------------------------------------
println("\n── Partition primality test vs trial division, n ∈ [2, 200] ──")
verify_range(2, 200)

println("\n── Primes in [2, 60] via partition test ──")
println("  ", filter(is_prime_partition, 2:60))

# ---------------------------------------------------------------------------
println("\n── MacMahonesque functions M_{vec_a}(n) ──")
println("  Weight convention: vec_a[1] = exponent of LARGEST part")
println()
@printf "  %4s  %10s  %10s  %10s\n" "n" "M_{(3,0)}" "M_{(2,2)}" "M_{(1,3)}"
println("  " * "-"^38)
for n in [4, 5, 6, 8, 10]
    v30 = M_macmahonesque([3,0], n)
    v22 = M_macmahonesque([2,2], n)
    v13 = M_macmahonesque([1,3], n)
    @printf "  %4d  %10d  %10d  %10d\n" n v30 v22 v13
end

# ---------------------------------------------------------------------------
println("\n── Quasi-shuffle product: U₍₁₎ * U₍₁,₁₎ ──")
println("  Paper result (p.12): 3·U₍₁,₁,₁₎ + (1/6)·U₍₃,₁₎ + (1/6)·U₍₁,₃₎ - (1/3)·U₍₁,₁₎")
println()
qs = quasishuffle_words([1], [1, 1])
for (w, c) in sort(collect(qs), by = p -> (length(p.first), p.first))
    println("    $c · M_$w")
end

# ---------------------------------------------------------------------------
println("\n── Convolution identity verification (Theorem 1.4.1) ──")
println("  Σᵢ M₍₁₎(i)·M₍₁₎(n-i) = (1/6)M₍₃₎(n) + 2M₍₁,₁₎(n) - (1/6)M₍₁₎(n)")
println()
conv_formula = quasishuffle_words([1], [1])  # U_(1) * U_(1)
let all_ok = true
    for n in 1:15
        lhs = sum(M_macmahonesque([1], i, Rational{BigInt}) *
                  M_macmahonesque([1], n - i, Rational{BigInt}) for i in 1:n-1;
                  init = Rational{BigInt}(0))
        rhs = evaluate_zq(conv_formula, n)
        status = lhs == rhs ? "ok" : "FAIL"
        all_ok &= (lhs == rhs)
        @printf "  n=%2d  lhs=%8s  rhs=%8s  [%s]\n" n lhs rhs status
    end
    println(all_ok ? "  All match." : "  MISMATCHES FOUND.")
end

# ---------------------------------------------------------------------------
println("\n── D operator: n·M₍₁₎(n) as constant-coefficient combination ──")
d1 = d_operator([1])
println("  n·M₍₁₎(n) =")
for (w, c) in sort(collect(d1), by = p -> (length(p.first), p.first))
    println("    $c · M_$w")
end
println()
println("  Verifying for n = 1..10:")
let all_ok = true
    for n in 1:10
        lhs = Rational{BigInt}(n) * M_macmahonesque([1], n, Rational{BigInt})
        rhs = evaluate_zq(d1, n)
        status = lhs == rhs ? "ok" : "FAIL"
        all_ok &= (lhs == rhs)
        @printf "  n=%2d  n·M₁=%6s  formula=%6s  [%s]\n" n lhs rhs status
    end
    println(all_ok ? "  All match." : "  MISMATCHES FOUND.")
end

# ---------------------------------------------------------------------------
println("\n── D operator paper formula: n·M₍₁,₁₎(n) (p.4) ──")
println("  (1/22)(−21M₍₃,₁₎ + 72M₍₂,₂₎ − 9M₍₁,₃₎ + 24M₍₃,₀₎")
println("         − 24M₍₃,₀,₀₎ − 24M₍₂,₁,₀₎ + 24M₍₂,₀,₁₎ − 72M₍₁,₁,₁₎)")
println()
paper_d11 = ZqElem(
    [3, 1]    => -21//big(22),
    [2, 2]    =>  72//big(22),
    [1, 3]    =>  -9//big(22),
    [3, 0]    =>  24//big(22),
    [3, 0, 0] => -24//big(22),
    [2, 1, 0] => -24//big(22),
    [2, 0, 1] =>  24//big(22),
    [1, 1, 1] => -72//big(22),
)
let all_ok = true
    for n in 1:10
        lhs = Rational{BigInt}(n) * M_macmahonesque([1, 1], n, Rational{BigInt})
        rhs = evaluate_zq(paper_d11, n)
        status = lhs == rhs ? "ok" : "FAIL"
        all_ok &= (lhs == rhs)
        @printf "  n=%2d  n·M₂=%8s  paper formula=%8s  [%s]\n" n lhs rhs status
    end
    println(all_ok ? "  All match." : "  MISMATCHES FOUND.")
end

# ---------------------------------------------------------------------------
println("\n── Symmetrisation (Theorem 4.4) ──")
println("  U^sym_{vec_a} = Σ_σ U_{σ(vec_a)} is quasimodular for odd-entry vectors")
println()
for vec in [[1, 3], [1, 1, 3], [3, 5]]
    sym = symmetrise(vec)
    println("  symmetrise($vec):")
    for (w, c) in sort(collect(sym), by = p -> p.first)
        println("    $c · U_$w")
    end
end

# ---------------------------------------------------------------------------
println("\n── Open conjecture test ──")
println("  Is every prime-vanishing expression in Q[n]-span of Table 1?")
println()
println("  E1..E4 use M₁..M₅.  E5 is the unknown fifth Table 1 entry.")
println("  Conjecture holds at (d, a_max ≤ 4): E1..E4 span the prime-vanishing subspace.")
println("  At a_max = 5, M₅ introduces new prime-vanishing directions — E5 is needed.")
println()
for (d, a) in [(2, 2), (3, 4), (4, 4)]
    r = test_conjecture(d, a; N=150, verbose=false)
    holds_str = r.holds ? "holds ✓" : "COUNTEREXAMPLE FOUND"
    @printf "  (d=%d, a_max=%d): basis=%d, pv_dim=%d, t1_rank=%d  → %s\n" d a r.dim_basis r.dim_prime_vanishing r.dim_table1_span holds_str
end
println()
println("  With a_max=5 (M₅ included), E5 is required:")
let r = test_conjecture(3, 5; N=150, verbose=false)
    @printf "  (d=3, a_max=5): pv_dim=%d, t1_rank=%d  → E5 needed (unknown)\n" r.dim_prime_vanishing r.dim_table1_span
end
println()
println("  For a systematic sweep: scan_conjecture(4, 4)")
println()
println("=" ^ 65)
println("Done.")
