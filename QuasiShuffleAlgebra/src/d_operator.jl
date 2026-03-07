"""
    rational_solve(A, b) -> Vector{Rational{BigInt}}

Solve A*x = b exactly over Q using RREF on the augmented matrix [A | b].
Handles overdetermined and underdetermined systems:
  - Overdetermined: checks consistency, then solves.
  - Underdetermined: returns a particular solution with free variables set to zero.
Throws an error if the system is inconsistent (no solution exists).
"""
function rational_solve(A::Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    m, n = size(A)
    aug = hcat(copy(A), reshape(copy(b), m, 1))

    # RREF on augmented matrix
    pivot_cols = Int[]
    pivot_row = 1
    for col in 1:(n + 1)
        piv = 0
        for row in pivot_row:m
            if !iszero(aug[row, col])
                piv = row
                break
            end
        end
        piv == 0 && continue

        # If pivot is in the RHS column, the system is inconsistent
        col == n + 1 && error("No exact solution (inconsistent system)")

        aug[[pivot_row, piv], :] = aug[[piv, pivot_row], :]
        aug[pivot_row, :] ./= aug[pivot_row, col]
        for row in 1:m
            row == pivot_row && continue
            iszero(aug[row, col]) && continue
            aug[row, :] .-= aug[row, col] .* aug[pivot_row, :]
        end
        push!(pivot_cols, col)
        pivot_row += 1
        pivot_row > m && break
    end

    # Check consistency: any row [0 ... 0 | nonzero] would have been caught above.
    # Read off particular solution (free variables = 0).
    x = zeros(Rational{BigInt}, n)
    for (r, col) in enumerate(pivot_cols)
        x[col] = aug[r, n + 1]
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
