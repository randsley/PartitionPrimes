@testset "Bernoulli numbers" begin
    # Known values
    @test bernoulli(0) == 1
    @test bernoulli(1) == -1//2
    @test bernoulli(2) == 1//6
    @test bernoulli(4) == -1//30
    @test bernoulli(6) == 1//42
    @test bernoulli(8) == -1//30
    @test bernoulli(10) == 5//66
    @test bernoulli(12) == -691//2730
    @test bernoulli(14) == 7//6
    @test bernoulli(16) == -3617//510
    @test bernoulli(18) == 43867//798
    @test bernoulli(20) == -174611//330

    # Odd Bernoulli numbers are zero for n > 1
    for n in [3, 5, 7, 9, 11, 13, 15, 17, 19]
        @test bernoulli(n) == 0
    end
end
