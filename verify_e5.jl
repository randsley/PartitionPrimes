import Pkg
Pkg.activate("/Users/nigelrandsley/GitHub/PartitionPrimes/QuasiShuffleAlgebra")
using QuasiShuffleAlgebra

# E5 from (d=3, a_max=5) counterexample — closes conjecture as a generator
function E5_d3(n::Int)
    m1 = M1(n); m2 = M2(n); m3 = M3(n); m4 = M4(n); m5 = M5(n)
    bn = big(n)
    return (
        -270270 * m1
        + 663549 * bn * m1
        + (-522351) * bn^2 * m1
        + 129072 * bn^3 * m1
        + (-315272) * bn^2 * m2
        + 30400 * bn^3 * m2
        + (-340864) * bn^2 * m3
        + 15872 * bn^3 * m3
        + (-193536) * bn^2 * m4
        + 154828800 * m5
    )
end

function conj_test(d, a_max, E_funcs; N_prime=200)
    composites = [n for n in 4:100 if !is_prime_trial(n)]
    primes_ns = filter(is_prime_trial, 2:N_prime)
    pm = eval_matrix(d, a_max, primes_ns)
    null_vecs = rational_nullspace(pm)
    cm = eval_matrix(d, a_max, composites)
    ne = length(E_funcs)
    e_vals = zeros(Rational{BigInt}, length(composites), ne*(d+1))
    for (i, n) in enumerate(composites)
        for (fi, Ef) in enumerate(E_funcs)
            v = Rational{BigInt}(Ef(n))
            for j in 0:d
                e_vals[i, j*ne + fi] = Rational{BigInt}(big(n)^j) * v
            end
        end
    end
    n_out = count(j -> !is_in_colspan(cm * null_vecs[:,j], e_vals), 1:size(null_vecs,2))
    return (holds = n_out == 0, n_outside = n_out, dim_pv = size(null_vecs,2))
end

E14 = [E1, E2, E3, E4]
E15 = [E1, E2, E3, E4, E5_d3]

println("Conjecture tests at a_max=5:")
for d in 2:6
    r4 = conj_test(d, 5, E14)
    r5 = conj_test(d, 5, E15)
    println("  d=$d: E1-E4 outside=$(r4.n_outside)/$(r4.dim_pv), E1-E5 outside=$(r5.n_outside)/$(r5.dim_pv)  $(r5.holds ? "✓" : "✗")")
end

println("\nConjecture at a_max=4 (should hold with E1-E4 only):")
for d in 2:5
    r = conj_test(d, 4, E14)
    println("  d=$d: outside=$(r.n_outside)/$(r.dim_pv)  $(r.holds ? "✓" : "✗")")
end
