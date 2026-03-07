@testset "D operator" begin
    @testset "d_operator([1]) numerical" begin
        result = d_operator([1])
        for n in 1:50
            lhs = Rational{BigInt}(n) * M_macmahonesque([1], n, Rational{BigInt})
            rhs = evaluate_zq(result, n)
            @test lhs == rhs
        end
    end

    @testset "d_operator([1,1]) numerical" begin
        result = d_operator([1, 1])
        for n in 1:50
            lhs = Rational{BigInt}(n) * M_macmahonesque([1, 1], n, Rational{BigInt})
            rhs = evaluate_zq(result, n)
            @test lhs == rhs
        end
    end

    # Explicit coefficients from paper p.4:
    # n*M_{(1,1)}(n) = (1/22)( -21*M_{(3,1)} + 72*M_{(2,2)} - 9*M_{(1,3)}
    #                          + 24*M_{(3,0)} - 24*M_{(3,0,0)} - 24*M_{(2,1,0)}
    #                          + 24*M_{(2,0,1)} - 72*M_{(1,1,1)} )
    @testset "d_operator([1,1]) explicit coefficients (paper p.4)" begin
        result = d_operator([1, 1])
        @test get(result, [3, 1],    0//1) == -21//22
        @test get(result, [2, 2],    0//1) ==  72//22
        @test get(result, [1, 3],    0//1) ==  -9//22
        @test get(result, [3, 0],    0//1) ==  24//22
        @test get(result, [3, 0, 0], 0//1) == -24//22
        @test get(result, [2, 1, 0], 0//1) == -24//22
        @test get(result, [2, 0, 1], 0//1) ==  24//22
        @test get(result, [1, 1, 1], 0//1) == -72//22
        # Verify no unexpected non-zero words beyond those in the paper
        expected_words = Set([[3,1],[2,2],[1,3],[3,0],[3,0,0],[2,1,0],[2,0,1],[1,1,1]])
        for (w, c) in result
            if w ∉ expected_words
                @test iszero(c)
            end
        end
    end
end
