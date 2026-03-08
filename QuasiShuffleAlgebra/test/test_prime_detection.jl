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

    # E5 vanishes at primes, non-negative at small composites
    @testset "E5 vanishes at primes" begin
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]
            @test iszero(E5(p))
        end
    end
    @testset "E5 non-negative at small composites" begin
        for n in [4, 6, 8, 9, 10, 12, 14, 16, 18, 20, 25]
            @test E5(n) >= 0
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
