"""
fit_macmahon.jl

Derive closed-form expressions for M_4, M_5, M_6 as linear combinations of:
  n^k · σ_{2j-1}(n)  and (for M_6 only) τ(n) (Ramanujan tau function)

At weight 12 (M_6), the cusp form Δ(q) contributes τ(n) = n-th coefficient of Δ.
Polynomial degrees: σ_{2j-1} gets degree a-j in n (matches M_2, M_3 pattern).
"""

import Pkg
Pkg.activate("/Users/nigelrandsley/GitHub/PartitionPrimes/QuasiShuffleAlgebra")
using QuasiShuffleAlgebra

function sigma(k::Int, n::Int)::Rational{BigInt}
    s = big(0)
    for d in 1:isqrt(n)
        if n % d == 0
            s += big(d)^k
            if d != n ÷ d
                s += big(n ÷ d)^k
            end
        end
    end
    return Rational{BigInt}(s)
end

"""
Compute Ramanujan tau for n=1..N via  q * prod(1-q^k)^24.
Multiply by (1-q^k) 24 times for each k.
"""
function ramanujan_tau(N::Int)::Vector{BigInt}
    coeff = zeros(BigInt, N + 1)
    coeff[1] = 1  # coefficient of q^0
    for k in 1:N
        for _ in 1:24
            for n in N:-1:k
                coeff[n+1] -= coeff[n+1-k]
            end
        end
    end
    # tau(n) = coeff of q^n in q*prod = coeff[n] (= coeff of q^{n-1} in prod)
    return [coeff[n] for n in 1:N]
end

function rational_rref!(M::Matrix{Rational{BigInt}})
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

function fit_Ma(a::Int, n_range; tau_vec=nothing)
    n_unknowns = sum(a - j + 1 for j in 1:a)
    use_tau = (tau_vec !== nothing)
    if use_tau; n_unknowns += 1; end

    n_pts = length(n_range)
    @assert n_pts >= n_unknowns "Need ≥ $n_unknowns points, got $n_pts"

    A = Matrix{Rational{BigInt}}(undef, n_pts, n_unknowns)
    b = Vector{Rational{BigInt}}(undef, n_pts)

    for (i, n) in enumerate(n_range)
        b[i] = Rational{BigInt}(M_direct(a, n))
        col = 0
        for j in 1:a
            sig = sigma(2j - 1, n)
            for k in 0:(a - j)
                col += 1
                A[i, col] = Rational{BigInt}(big(n)^k) * sig
            end
        end
        if use_tau
            A[i, col+1] = Rational{BigInt}(tau_vec[n])
        end
    end

    aug = hcat(A, reshape(b, :, 1))
    pivot_cols = rational_rref!(aug)
    n_col = n_unknowns + 1

    (n_col in pivot_cols) && error("System is inconsistent!")
    (length(pivot_cols) < n_unknowns) && error("Underdetermined: rank=$(length(pivot_cols)) < $n_unknowns")

    x = zeros(Rational{BigInt}, n_unknowns)
    for (r, pc) in enumerate(pivot_cols)
        x[pc] = aug[r, n_col]
    end
    return x
end

function verify_fit(a::Int, x::Vector{Rational{BigInt}}, n_range; tau_vec=nothing)
    println("Verifying M_$a on n=$(first(n_range))..$(last(n_range))...")
    use_tau = (tau_vec !== nothing)
    for n in n_range
        expected = Rational{BigInt}(M_direct(a, n))
        computed = Rational{BigInt}(0)
        col = 0
        for j in 1:a
            sig = sigma(2j - 1, n)
            for k in 0:(a - j)
                col += 1
                computed += x[col] * Rational{BigInt}(big(n)^k) * sig
            end
        end
        if use_tau
            col += 1
            computed += x[col] * Rational{BigInt}(tau_vec[n])
        end
        if expected != computed
            println("  MISMATCH at n=$n: expected=$expected, got=$computed")
            return false
        end
    end
    println("  ✓ All match.")
    return true
end

function print_formula(a::Int, x::Vector{Rational{BigInt}}; has_tau=false)
    println("M_$a formula:")
    col = 0
    for j in 1:a
        terms = String[]
        for k in 0:(a - j)
            col += 1
            c = x[col]
            iszero(c) && continue
            s = k == 0 ? "$c" : k == 1 ? "$(c)*n" : "$(c)*n^$k"
            push!(terms, s)
        end
        isempty(terms) || println("  σ_{$(2j-1)}: ", join(terms, " + "))
    end
    if has_tau
        col += 1
        c = x[col]
        iszero(c) || println("  τ(n): $c")
    end
end

# ---- Precompute tau ----
println("Computing τ(n) for n=1..65...")
tau = ramanujan_tau(65)
@assert tau[1] == 1 && tau[2] == -24 && tau[3] == 252 && tau[5] == 4830
println("  Sanity check ✓  τ(1)=$(tau[1]), τ(2)=$(tau[2]), τ(5)=$(tau[5])")

# ---- M_4 ----
println("\n" * "="^55)
println("Fitting M_4 ($(4*5÷2) unknowns)...")
x4 = fit_Ma(4, 1:30)
print_formula(4, x4)
verify_fit(4, x4, 1:55)

# ---- M_5 ----
println("\n" * "="^55)
println("Fitting M_5 ($(5*6÷2) unknowns)...")
x5 = fit_Ma(5, 1:35)
print_formula(5, x5)
verify_fit(5, x5, 1:55)

# ---- M_6 (with tau) ----
println("\n" * "="^55)
println("Fitting M_6 ($(6*7÷2)+1=$(6*7÷2+1) unknowns, including τ)...")
x6 = fit_Ma(6, 1:50; tau_vec=tau)
print_formula(6, x6; has_tau=true)
verify_fit(6, x6, 1:60; tau_vec=tau)

println("\n✓ All fits complete!")
println()
println("="^60)
println("Julia function definitions for macmahon.jl:")
println()

function julia_fn(a, x; has_tau=false)
    lines = String[]
    push!(lines, "function M$(a)_closed(n::Int)")
    col = 0
    for j in 1:a
        has_term = any(!iszero(x[(a-j+1 > 0 ? (col+1) : 1):min(col+a-j+1, length(x))]) for _ in 1:1)
        # recount
        base = col
        has_nonzero = any(!iszero(x[base+k+1]) for k in 0:(a-j))
        if has_nonzero
            push!(lines, "    s$(2j-1) = sigma($(2j-1), n)")
        end
        col += a - j + 1
    end
    push!(lines, "    return (")
    col = 0
    for j in 1:a
        for k in 0:(a-j)
            col += 1
            c = x[col]
            iszero(c) && continue
            num, den = numerator(c), denominator(c)
            nterm = k == 0 ? "" : k == 1 ? " * big(n)" : " * big(n)^$k"
            push!(lines, "        + ($num // big($den))$nterm * s$(2j-1)")
        end
    end
    if has_tau
        col += 1
        c = x[col]
        if !iszero(c)
            num, den = numerator(c), denominator(c)
            push!(lines, "        + ($num // big($den)) * tau_n  # Ramanujan τ(n)")
        end
    end
    push!(lines, "    )")
    push!(lines, "end")
    return join(lines, "\n")
end

println(julia_fn(4, x4))
println()
println(julia_fn(5, x5))
println()
println(julia_fn(6, x6; has_tau=true))
