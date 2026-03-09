@testset "Prime detection" begin
    # E1 matches trial division for n in 2:200
    @testset "E1 vs trial division" begin
        for n in 2:200
            @test is_prime_partition(n) == is_prime_trial(n)
        end
    end

    # E1 non-negative for all n >= 2
    @testset "E1 non-negativity" begin
        for n in 2:100
            @test E1(n) >= 0
        end
    end

    # E2 matches trial division (Theorem 1.1(2) / Table 1)
    @testset "E2 vs trial division" begin
        for n in 2:100
            @test iszero(E2(n)) == is_prime_trial(n)
        end
    end

    # Verify M1 = sigma_1
    @testset "M1 consistency" begin
        for n in 1:50
            @test M1(n) == σ(1, n)
            @test M1(n) == M_macmahonesque([1], n)
        end
    end

    # E5 vanishes precisely at primes (minimal d=2 canonical form)
    @testset "E5 vanishes at primes" begin
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]
            @test iszero(E5(p))
        end
    end
    @testset "E5 non-zero at composites 4..64" begin
        # E5 is not universally non-negative (may be negative at larger composites),
        # but all composites in [4,64] are positive.
        for n in filter(n -> !is_prime_trial(n), 4:64)
            @test !iszero(E5(n))
        end
    end
    @testset "E5 + 856*E4 non-negative at composites 4..100" begin
        # E5 alone can be negative (219 composites in [4,500]).
        # The combination E5 + 856*E4 is non-negative on all composites in [4,500];
        # 856 is the smallest such integer constant.
        for n in filter(n -> !is_prime_trial(n), 4:100)
            @test E5(n) + 856*E4(n) >= 0
        end
    end

    # M_direct consistency with M_macmahonesque
    @testset "M_direct vs M_macmahonesque" begin
        for n in 1:30
            @test M_direct(1, n) == M_macmahonesque([1], n)
            @test M_direct(2, n) == M_macmahonesque([1, 1], n)
            @test M_direct(3, n) == M_macmahonesque([1, 1, 1], n)
        end
    end

    # Zero-entry words: M_{(0)}(n) counts partitions into one part size s|n
    # with any multiplicity m, each contributing m^0 = 1. So M_{(0)}(n) = d(n) = σ(0,n).
    @testset "Zero-entry word M_{(0)}(n) == d(n)" begin
        for n in 1:30
            @test M_macmahonesque([0], n) == σ(0, n)
        end
    end
end
