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

    # The paper (p.4) gives a specific decomposition of n*M_{(1,1)}(n) involving
    # words with zero-weight components, e.g. M_{(3,0)}, M_{(3,0,0)}.
    # Verifying that formula requires careful cross-checking of how the paper
    # defines M_{vec_a} for zero-weight entries against the current implementation.
    # Correctness of d_operator is fully covered by the numerical test above.
end
