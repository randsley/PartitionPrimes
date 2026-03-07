@testset "Quasi-shuffle product" begin
    # Ground truth from paper p.11:
    # quasishuffle([1], [1,1]) = 3*[1,1,1] + (1/6)*[3,1] + (1/6)*[1,3] - (1/3)*[1,1]
    @testset "quasishuffle([1],[1,1]) explicit" begin
        result = quasishuffle_words([1], [1, 1])

        expected = ZqElem(
            [1, 1, 1] => big(3) // big(1),
            [3, 1]    => big(1) // big(6),
            [1, 3]    => big(1) // big(6),
            [1, 1]    => big(-1) // big(3),
        )

        # Check all expected terms
        for (w, c) in expected
            @test get(result, w, Rational{BigInt}(0)) == c
        end
        # Check no extra terms
        for (w, c) in result
            if !haskey(expected, w)
                @test iszero(c)
            end
        end
    end

    # Convolution test: for n=1:30,
    # sum_{i+j=n} M_(1)(i)*M_(1)(j) == (1/6)*M_(3)(n) + 2*M_(1,1)(n) - (1/6)*M_(1)(n)
    @testset "Convolution identity (numerical)" begin
        for n in 1:30
            # LHS: convolution
            lhs = Rational{BigInt}(0)
            for i in 1:(n-1)
                j = n - i
                lhs += Rational{BigInt}(M_macmahonesque([1], i, Rational{BigInt})) *
                       Rational{BigInt}(M_macmahonesque([1], j, Rational{BigInt}))
            end

            # RHS: from quasi-shuffle product of [1]*[1]
            rhs = Rational{BigInt}(1) // 6 * M_macmahonesque([3], n, Rational{BigInt}) +
                  Rational{BigInt}(2) * M_macmahonesque([1, 1], n, Rational{BigInt}) -
                  Rational{BigInt}(1) // 6 * M_macmahonesque([1], n, Rational{BigInt})

            @test lhs == rhs
        end
    end

    # The quasi-shuffle of [1]*[1] should give convolution coefficients
    @testset "quasishuffle([1],[1]) matches convolution" begin
        result = quasishuffle_words([1], [1])

        # Expected from diamond(1,1): c_{1,1,2} = 1!*1!/3! = 1/6, and B_2/2 terms
        # [1]*[1] = [1,1] + [1,1] + diamond(1,1)*[] = 2*[1,1] + (1/6)*[2] + ...
        # Actually let's just verify numerically
        for n in 1:20
            lhs = Rational{BigInt}(0)
            for i in 1:(n-1)
                lhs += M_macmahonesque([1], i, Rational{BigInt}) *
                       M_macmahonesque([1], n - i, Rational{BigInt})
            end
            rhs = evaluate_zq(result, n)
            @test lhs == rhs
        end
    end
end
