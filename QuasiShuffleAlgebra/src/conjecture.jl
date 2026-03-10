"""
conjecture.jl

Computational test of the open conjecture (Extend.md, §Open Conjecture):

  Any prime-detecting expression E(n) = Σ pᵢ(n)·Mᵢ(n) ≥ 0
  (vanishing precisely at primes) is a Q[n]-linear combination
  of the five entries in Table 1.

Strategy (fixed degree bound d, weight bound a_max):
  1. Build basis B = { n^k · M_a(n) : 0 ≤ k ≤ d, 1 ≤ a ≤ a_max }
     Each Table 1 entry Eᵢ is a coefficient vector in this basis.
  2. Evaluate B at enough primes to determine the prime-vanishing subspace
     (the null space of the prime evaluation matrix, over Q).
  3. Check: does every basis vector of that null space lie in span{E1,...,E5}?
     If yes → conjecture holds at (d, a_max).
     If no  → coefficient vector of a counterexample is returned.

Note: this tests the linear part of the conjecture (vanishing at primes).
The non-negativity condition (E(n) ≥ 0 at composites) requires additional
positivity analysis beyond linear algebra.
"""

# ---------------------------------------------------------------------------
# Basis indexing
# ---------------------------------------------------------------------------

"""
    build_basis(d, a_max) -> Vector{Tuple{Int,Int}}

Ordered list of (k, a) pairs: the monomial basis { n^k · M_a(n) }.
Ordering: a varies in outer loop, k in inner — matches column layout of
eval_matrix so that Table 1 coefficient vectors are easy to assemble.
"""
function build_basis(d::Int, a_max::Int)::Vector{Tuple{Int,Int}}
    return [(k, a) for a in 1:a_max for k in 0:d]
end

"""
    basis_index(d, k, a) -> Int

Column index of (k, a) in the basis produced by build_basis(d, a_max).
"""
basis_index(d::Int, k::Int, a::Int) = (a - 1) * (d + 1) + k + 1

# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

"""
    eval_Ma(a, n) -> Rational{BigInt}

Evaluate M_a(n) using closed forms for a ≤ 6, M_direct otherwise.
M4–M6 use divisor-sum formulas (fast); M6 also uses ramanujan_tau.
"""
function eval_Ma(a::Int, n::Int)::Rational{BigInt}
    if a == 1
        return Rational{BigInt}(M1(n))
    elseif a == 2
        return M2(n)
    elseif a == 3
        return M3(n)
    elseif a == 4
        return M4(n)
    elseif a == 5
        return M5(n)
    elseif a == 6
        return M6(n)
    elseif a == 7
        return M7(n)
    elseif a == 8
        return M8(n)
    elseif a == 9
        return M9(n)
    elseif a == 10
        return M10(n)
    elseif a == 11
        return M11(n)
    elseif a == 12
        return M12(n)
    else
        return Rational{BigInt}(M_direct(a, n))
    end
end

"""
    eval_matrix(d, a_max, ns) -> Matrix{Rational{BigInt}}

Rows = values of n, columns = basis elements { n^k · M_a(n) }.
"""
function eval_matrix(d::Int, a_max::Int, ns::AbstractVector{Int})::Matrix{Rational{BigInt}}
    basis = build_basis(d, a_max)
    mat = Matrix{Rational{BigInt}}(undef, length(ns), length(basis))
    for (i, n) in enumerate(ns)
        Ma_cache = Dict{Int,Rational{BigInt}}()
        for (j, (k, a)) in enumerate(basis)
            ma = get!(Ma_cache, a) do; eval_Ma(a, n) end
            mat[i, j] = big(n)^k * ma
        end
    end
    return mat
end

# ---------------------------------------------------------------------------
# Table 1 coefficient vectors
# ---------------------------------------------------------------------------

"""
    table1_coeffs(d, a_max) -> Matrix{Rational{BigInt}}

Returns a matrix whose columns are the coefficient vectors of E1..E5
in the basis { n^k · M_a }.

E1 = (n²-3n+2)·M₁ - 8·M₂
E2 = (3n³-13n²+18n-8)·M₁ + (12n²-120n+212)·M₂ - 960·M₃
E3 = (25n⁴-171n³+423n²-447n+170)·M₁
     + (300n³-3554n²+12900n-14990)·M₂
     + (2400n²-60480n+214080)·M₃
     - 725760·M₄
E4 = (126n⁵-1303n⁴+5073n³-9323n²+8097n-2670)·M₁
     + (3024n⁴-48900n³+288014n²-737100n+695490)·M₂
     + (60480n³-1510080n²+10644480n-23496480)·M₃
     + (725760n²-36288000n+218453760)·M₄
     - 580608000·M₅
E5 = (-450450 + 675675n - 225225n²)·M₁ + (960960n - 120120n²)·M₂
     + (2534912n - 166016n²)·M₃ + (7999488n - 322560n²)·M₄ + 258048000·M₅
(Minimal polynomial degree d=2; vanishes iff prime; may be negative at composites.)
"""
function table1_coeffs(d::Int, a_max::Int)::Matrix{Rational{BigInt}}
    dim = (d + 1) * a_max
    T = Matrix{Rational{BigInt}}(undef, dim, 5)
    fill!(T, Rational{BigInt}(0))

    idx(k, a) = basis_index(d, k, a)

    # Helper: set coefficient if (k,a) is within bounds
    function set!(col, k, a, val)
        (k ≤ d && a ≤ a_max) || return
        T[idx(k, a), col] = Rational{BigInt}(val)
    end

    # E1: (n²-3n+2)·M₁ - 8·M₂
    set!(1, 2, 1,  1)
    set!(1, 1, 1, -3)
    set!(1, 0, 1,  2)
    set!(1, 0, 2, -8)

    # E2: (3n³-13n²+18n-8)·M₁ + (12n²-120n+212)·M₂ - 960·M₃
    set!(2, 3, 1,   3)
    set!(2, 2, 1, -13)
    set!(2, 1, 1,  18)
    set!(2, 0, 1,  -8)
    set!(2, 2, 2,   12)
    set!(2, 1, 2, -120)
    set!(2, 0, 2,  212)
    set!(2, 0, 3, -960)

    # E3
    set!(3, 4, 1,   25)
    set!(3, 3, 1, -171)
    set!(3, 2, 1,  423)
    set!(3, 1, 1, -447)
    set!(3, 0, 1,  170)
    set!(3, 3, 2,    300)
    set!(3, 2, 2,  -3554)
    set!(3, 1, 2,  12900)
    set!(3, 0, 2, -14990)
    set!(3, 2, 3,    2400)
    set!(3, 1, 3,  -60480)
    set!(3, 0, 3,  214080)
    set!(3, 0, 4, -725760)

    # E4
    set!(4, 5, 1,     126)
    set!(4, 4, 1,   -1303)
    set!(4, 3, 1,    5073)
    set!(4, 2, 1,   -9323)
    set!(4, 1, 1,    8097)
    set!(4, 0, 1,   -2670)
    set!(4, 4, 2,     3024)
    set!(4, 3, 2,   -48900)
    set!(4, 2, 2,   288014)
    set!(4, 1, 2,  -737100)
    set!(4, 0, 2,   695490)
    set!(4, 3, 3,     60480)
    set!(4, 2, 3,  -1510080)
    set!(4, 1, 3,  10644480)
    set!(4, 0, 3, -23496480)
    set!(4, 2, 4,     725760)
    set!(4, 1, 4,  -36288000)
    set!(4, 0, 4,  218453760)
    set!(4, 0, 5, -580608000)

    # E5: minimal-degree (d=2) canonical form
    set!(5, 0, 1,  -450450)
    set!(5, 1, 1,   675675)
    set!(5, 2, 1,  -225225)
    set!(5, 1, 2,   960960)
    set!(5, 2, 2,  -120120)
    set!(5, 1, 3,  2534912)
    set!(5, 2, 3,  -166016)
    set!(5, 1, 4,  7999488)
    set!(5, 2, 4,  -322560)
    set!(5, 0, 5, 258048000)

    return T
end

# ---------------------------------------------------------------------------
# Rational linear algebra
# ---------------------------------------------------------------------------

"""
    rational_rref!(M) -> Vector{Int}

In-place rational row-reduced echelon form over Q. Returns pivot column indices.
"""
function rational_rref!(M::Matrix{Rational{BigInt}})::Vector{Int}
    m, n = size(M)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        piv = findfirst(!iszero, view(M, row:m, col))
        piv === nothing && continue
        piv += row - 1
        M[[row, piv], :] = M[[piv, row], :]
        scale = M[row, col]
        M[row, :] ./= scale
        for r in 1:m
            r == row && continue
            iszero(M[r, col]) && continue
            M[r, :] .-= M[r, col] .* M[row, :]
        end
        push!(pivot_cols, col)
        row += 1
        row > m && break
    end
    return pivot_cols
end

"""
    rank_over_Q(M) -> Int
"""
function rank_over_Q(M::Matrix{Rational{BigInt}})::Int
    return length(rational_rref!(copy(M)))
end

"""
    rational_nullspace(M) -> Matrix{Rational{BigInt}}

Null space of M over Q. Columns are basis vectors.
"""
function rational_nullspace(M::Matrix{Rational{BigInt}})::Matrix{Rational{BigInt}}
    _, n = size(M)
    A = copy(M)
    pivot_cols = rational_rref!(A)
    free_cols = setdiff(1:n, pivot_cols)
    isempty(free_cols) && return Matrix{Rational{BigInt}}(undef, n, 0)
    null_vecs = zeros(Rational{BigInt}, n, length(free_cols))
    for (j, fc) in enumerate(free_cols)
        null_vecs[fc, j] = Rational{BigInt}(1)
        for (r, pc) in enumerate(pivot_cols)
            null_vecs[pc, j] = -A[r, fc]
        end
    end
    return null_vecs
end

"""
    is_in_colspan(v, S) -> Bool

Test whether vector v lies in the column span of S over Q.
"""
function is_in_colspan(v::Vector{Rational{BigInt}},
                       S::Matrix{Rational{BigInt}})::Bool
    aug = hcat(S, reshape(v, :, 1))
    return rank_over_Q(S) == rank_over_Q(aug)
end

# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

"""
    test_conjecture(d, a_max; N=300, verbose=true) -> NamedTuple

Test whether every prime-vanishing expression in basis B = {n^k·M_a : k≤d, a≤a_max}
lies in the Q[n]-span of Table 1 entries E1..E5 (i.e., Q-linear combinations of
n^j·Eᵢ for j=0..d, i=1..5).

Keyword arguments:
  N       — evaluate at n = 2..N (more = more reliable null space)
  verbose — print progress summary

Returns NamedTuple with fields:
  holds               :: Bool
  dim_basis           :: Int  (total basis dimension)
  dim_prime_vanishing :: Int  (dimension of prime-vanishing subspace)
  dim_table1_span     :: Int  (rank of Table 1 vectors, within bounds)
  counterexample      :: Union{Nothing, Vector{Rational{BigInt}}}
"""
function test_conjecture(d::Int, a_max::Int;
                         N::Int = 300,
                         verbose::Bool = true)
    ns = collect(2:N)
    primes_in_range = filter(is_prime_trial, ns)

    verbose && @printf("Basis: d=%d, a_max=%d → %d elements\n",
                       d, a_max, (d+1)*a_max)
    verbose && @printf("Evaluating at %d primes in [2,%d]…\n",
                       length(primes_in_range), N)

    # Step 1: prime evaluation matrix
    prime_mat = eval_matrix(d, a_max, primes_in_range)

    # Step 2: null space = prime-vanishing subspace
    null_vecs = rational_nullspace(prime_mat)
    dim_pv = size(null_vecs, 2)
    verbose && @printf("Prime-vanishing subspace dimension: %d\n", dim_pv)

    if dim_pv == 0
        verbose && println("Trivial: only the zero expression vanishes at all primes.")
        return (holds=true, dim_basis=(d+1)*a_max,
                dim_prime_vanishing=0, dim_table1_span=0,
                counterexample=nothing)
    end

    # Step 3: Table 1 Q[n]-span: include n^j · Eᵢ(n) for j=0..d, i=1..5
    # The conjecture asserts Q[n]-linear combinations of Table 1 entries,
    # so n·E1, n²·E1, n·E2, … are all valid generators.
    #
    # IMPORTANT: evaluate E1..E4 directly (not via truncated coefficient vectors).
    # table1_coeffs silently drops terms with k>d or a>a_max, so using
    # all_mat * t1c would give wrong values when E_i has degree > d.
    all_mat = eval_matrix(d, a_max, ns)      # (length(ns) × dim_basis)
    E_funcs = [E1, E2, E3, E4, E5]
    t1_base = zeros(Rational{BigInt}, length(ns), 5)
    for (i, n) in enumerate(ns)
        for (fi, Ef) in enumerate(E_funcs)
            t1_base[i, fi] = Rational{BigInt}(Ef(n))
        end
    end

    # Build Q[n]-span: columns are n^j · Eᵢ(n) for j=0..d, i=1..5
    ns_big = Rational{BigInt}[big(n) for n in ns]
    n_powers = [ns_big .^ j for j in 0:d]   # d+1 vectors, each length(ns)
    t1_qn_cols = Vector{Rational{BigInt}}[]
    for j in 0:d
        pow = n_powers[j+1]
        for i in 1:5
            push!(t1_qn_cols, pow .* t1_base[:, i])
        end
    end
    t1_vals = hcat(t1_qn_cols...)            # (length(ns) × 5*(d+1))

    dim_t1 = rank_over_Q(copy(t1_vals))
    verbose && @printf("Table 1 Q[n]-span rank (within bounds): %d\n", dim_t1)

    # Check each null space vector
    counterexample = nothing
    for j in 1:dim_pv
        coeff_vec = null_vecs[:, j]
        expr_vals = all_mat * coeff_vec      # values of expression over all n
        if !is_in_colspan(expr_vals, t1_vals)
            counterexample = coeff_vec
            break
        end
    end

    holds = isnothing(counterexample)
    verbose && println(holds ? "✓ Conjecture holds at (d=$d, a_max=$a_max)." :
                               "✗ Counterexample found at (d=$d, a_max=$a_max).")

    return (holds=holds, dim_basis=(d+1)*a_max,
            dim_prime_vanishing=dim_pv,
            dim_table1_span=dim_t1,
            counterexample=counterexample)
end

"""
    scan_conjecture(d_max, a_max; N=300) -> Nothing

Run test_conjecture for all (d, a) with 1 ≤ d ≤ d_max, 1 ≤ a ≤ a_max.
Stops and reports immediately if a counterexample is found.
"""
function scan_conjecture(d_max::Int, a_max::Int; N::Int = 300)
    for a in 1:a_max, d in 1:d_max
        result = test_conjecture(d, a; N=N, verbose=false)
        status = result.holds ? "✓" : "✗ COUNTEREXAMPLE"
        @printf("(d=%d, a=%d): pv_dim=%d, t1_rank=%d  %s\n",
                d, a, result.dim_prime_vanishing, result.dim_table1_span, status)
        if !result.holds
            println("\nConjecture FAILS at (d=$d, a=$a). Counterexample coefficient vector:")
            display(result.counterexample)
            return
        end
    end
    println("\nAll cases up to (d=$d_max, a=$a_max): conjecture holds.")
end

export build_basis, eval_matrix, table1_coeffs
export rational_rref!, rational_nullspace, rank_over_Q, is_in_colspan
export test_conjecture, scan_conjecture
