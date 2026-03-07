"""
test_conjecture.jl

Tests for the conjecture computational machinery.

Sanity checks (not the conjecture itself — that is the research question):
  1. Table 1 vectors are correctly encoded: E1(p) = 0 at primes, E1(n) > 0 at composites
  2. rational_nullspace gives a true null space (M * v = 0)
  3. rank_over_Q is consistent with nullspace dimension
  4. test_conjecture(2, 2) runs without error and reports a non-trivial prime-vanishing subspace
"""

using Test

@testset "Conjecture infrastructure" begin

    @testset "Table 1 coefficient vectors: E1 sanity" begin
        # E1 = (n²-3n+2)M₁ - 8M₂ vanishes at primes, positive at composites
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]
            @test iszero(E1(p))
        end
        for c in [4, 6, 8, 9, 10, 12, 14, 15, 16]
            @test E1(c) > 0
        end
    end

    @testset "eval_matrix: basis (k=0,a=1) column equals M1" begin
        ns = collect(2:20)
        d, a_max = 1, 2
        mat = eval_matrix(d, a_max, ns)
        # Column 1 = (k=0, a=1) = M₁(n) = σ₁(n)
        for (i, n) in enumerate(ns)
            @test mat[i, 1] == Rational{BigInt}(M1(n))
        end
    end

    @testset "rational_nullspace: M * v = 0 for all null vectors" begin
        M = Rational{BigInt}[1 2 3; 4 5 6; 7 8 9]
        N = rational_nullspace(M)
        @test size(N, 2) > 0
        tol = zeros(Rational{BigInt}, 3)
        for j in 1:size(N, 2)
            @test M * N[:, j] == tol
        end
    end

    @testset "rank_over_Q: known ranks" begin
        # Full rank
        M = Rational{BigInt}[1 0; 0 1; 0 0]
        @test rank_over_Q(M) == 2
        # Rank-1 matrix
        M2 = Rational{BigInt}[1 2; 2 4; 3 6]
        @test rank_over_Q(M2) == 1
    end

    @testset "test_conjecture(2, 2) runs and finds non-trivial subspace" begin
        result = test_conjecture(2, 2; N=100, verbose=false)
        # E1 lies in the (d=2, a=2) basis, so prime-vanishing dim ≥ 1
        @test result.dim_prime_vanishing >= 1
        # Table 1 rank should be ≥ 1 within these bounds
        @test result.dim_table1_span >= 1
    end

end
