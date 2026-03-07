"""
    rational_solve(A, b) -> Vector{Rational{BigInt}}

Solve A*x = b exactly over Rational{BigInt} using Gaussian elimination.
"""
function rational_solve(A::Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    m, n = size(A)
    aug = hcat(copy(A), reshape(copy(b), m, 1))

    pivot_row = 1
    pivot_cols = Int[]
    for col in 1:n
        found = 0
        for row in pivot_row:m
            if !iszero(aug[row, col])
                found = row
                break
            end
        end
        found == 0 && continue

        if found != pivot_row
            for j in 1:(n+1)
                aug[pivot_row, j], aug[found, j] = aug[found, j], aug[pivot_row, j]
            end
        end

        push!(pivot_cols, col)
        piv = aug[pivot_row, col]
        for row in (pivot_row + 1):m
            if !iszero(aug[row, col])
                factor = aug[row, col] // piv
                for j in col:(n + 1)
                    aug[row, j] -= factor * aug[pivot_row, j]
                end
            end
        end
        pivot_row += 1
    end

    rank = length(pivot_cols)
    for row in (rank + 1):m
        if !iszero(aug[row, n + 1])
            error("No exact solution (inconsistent)")
        end
    end

    rank == n || error("Underdetermined system: rank $rank < $n unknowns; basis may be incomplete")

    x = zeros(Rational{BigInt}, n)
    for k in rank:-1:1
        col = pivot_cols[k]
        s = aug[k, n + 1]
        for j in (k + 1):rank
            s -= aug[k, pivot_cols[j]] * x[pivot_cols[j]]
        end
        x[col] = s // aug[k, col]
    end

    return x
end

"""
    d_operator(vec_a::Word) -> ZqElem

Compute the D operator: express n*M_{vec_a}(n) as a constant-coefficient
linear combination of MacMahonesque functions via coefficient matching.
"""
function d_operator(vec_a::Word)::ZqElem
    weight = sum(vec_a)
    len = length(vec_a)
    max_weight = weight + len + 1

    basis = all_words_up_to_weight(max_weight; max_length = len + 1, min_length = 1)
    num_basis = length(basis)
    # Use num_basis + 5 equations (overdetermined) so the inconsistency check
    # in rational_solve can catch a wrong basis before back-substitution.
    num_eqs = num_basis + 5

    V = Matrix{Rational{BigInt}}(undef, num_eqs, num_basis)
    rhs = Vector{Rational{BigInt}}(undef, num_eqs)

    for n in 1:num_eqs
        rhs[n] = Rational{BigInt}(n) * M_macmahonesque(vec_a, n, Rational{BigInt})
        for j in 1:num_basis
            V[n, j] = M_macmahonesque(basis[j], n, Rational{BigInt})
        end
    end

    c = rational_solve(V, rhs)

    result = ZqElem()
    for j in 1:num_basis
        if !iszero(c[j])
            result[basis[j]] = c[j]
        end
    end

    return result
end
