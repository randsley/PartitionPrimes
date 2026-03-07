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

    # Paper (p.4) explicit decomposition of n*M_{(1,1)}(n) (Theorem 1.4 example).
    # NOTE on convention: in M_macmahonesque(vec_a, n), vec_a[1] is the exponent
    # of the LARGEST part (deepest recursion), matching the paper's formula notation
    # where subscript v_1 corresponds to the multiplicity of the largest part size.
    @testset "d_operator([1,1]) paper decomposition is numerically equivalent" begin
        paper_result = ZqElem(
            [3, 1]    => -21//big(22),
            [2, 2]    =>  72//big(22),
            [1, 3]    =>  -9//big(22),
            [3, 0]    =>  24//big(22),
            [3, 0, 0] => -24//big(22),
            [2, 1, 0] => -24//big(22),
            [2, 0, 1] =>  24//big(22),
            [1, 1, 1] => -72//big(22),
        )
        for n in 1:50
            lhs = Rational{BigInt}(n) * M_macmahonesque([1, 1], n, Rational{BigInt})
            rhs = evaluate_zq(paper_result, n)
            @test lhs == rhs
        end
    end
end
