"""
compute_m6_coefficients.jl — Print the full closed-form polynomial coefficients for M6(n).

Run from PartitionPrimes root:
  julia --project=QuasiShuffleAlgebra Paper/compute_m6_coefficients.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "QuasiShuffleAlgebra"))
using QuasiShuffleAlgebra

function solve_m6()
    ns = 1:55

    # Build design matrix: columns ordered as (j,k) for j=0..5, k=0..(5-j), then tau
    A = zeros(Rational{BigInt}, length(ns), 22)
    col = 1
    for j in 0:5
        for k in 0:(5-j)
            for (i, n) in enumerate(ns)
                A[i, col] = Rational{BigInt}(big(n)^k) * Rational{BigInt}(σ(2j+1, n))
            end
            col += 1
        end
    end
    for (i, n) in enumerate(ns)
        A[i, 22] = Rational{BigInt}(ramanujan_tau(n))
    end

    b = Rational{BigInt}[Rational{BigInt}(M_direct(6, n)) for n in ns]

    # Gaussian elimination over Q
    Ab = hcat(A, reshape(b, :, 1))
    m = size(Ab, 1)
    pivot_row = 1
    for c in 1:22
        pr = findfirst(!iszero, Ab[pivot_row:end, c])
        pr === nothing && continue
        pr += pivot_row - 1
        Ab[pivot_row, :], Ab[pr, :] = Ab[pr, :], Ab[pivot_row, :]
        Ab[pivot_row, :] ./= Ab[pivot_row, c]
        for r in 1:m
            r == pivot_row && continue
            iszero(Ab[r, c]) && continue
            Ab[r, :] -= Ab[r, c] * Ab[pivot_row, :]
        end
        pivot_row += 1
        pivot_row > m && break
    end

    return Ab[1:22, end]
end

function print_m6_formula(sol)
    println("M₆(n) = Σⱼ Pⱼ(n)·σ_{2j+1}(n)  −  (17/150450048000)·τ(n)")
    println()

    col = 1
    for j in 0:5
        terms = String[]
        for k in 0:(5-j)
            v = sol[col]
            if !iszero(v)
                npart = k == 0 ? "" : k == 1 ? "*n" : "*n^$k"
                push!(terms, "($v)$npart")
            end
            col += 1
        end
        poly = isempty(terms) ? "0" : join(terms, " + ")
        println("  P_$j(n) = $poly")
        println("           [coefficient of σ_{$(2j+1)}(n)]")
        println()
    end
    println("  c_τ = $(sol[22])")
end

function verify_m6(sol)
    errors = 0
    for n in 1:30
        manual = begin
            acc = zero(Rational{BigInt})
            c = 1
            for j in 0:5
                for k in 0:(5-j)
                    acc += sol[c] * Rational{BigInt}(big(n)^k) * Rational{BigInt}(σ(2j+1, n))
                    c += 1
                end
            end
            acc += sol[22] * Rational{BigInt}(ramanujan_tau(n))
            acc
        end
        if manual != Rational{BigInt}(M_direct(6, n))
            errors += 1
            println("MISMATCH at n=$n")
        end
    end
    println("Verification M6(n=1..30): $(errors == 0 ? "✓ all match" : "✗ $errors mismatches")")
end

sol = solve_m6()
print_m6_formula(sol)
verify_m6(sol)
