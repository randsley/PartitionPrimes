using QuasiShuffleAlgebra
using Printf

println("=" ^ 60)
println("Quasi-Shuffle Algebra for Partition-Theoretic Prime Detection")
println("Craig, van Ittersum & Ono (2024)")
println("=" ^ 60)

println("\n-- MacMahon function values --")
for n in [2, 3, 4, 5, 6, 7, 8, 10, 12, 13]
    @printf "  n=%2d  M1=%4d  M2=%s  M3=%s\n" n M1(n) M2(n) M3(n)
end

println("\n-- Prime-detecting expression E1(n) --")
println("  (should be 0 iff n is prime)")
for n in 2:20
    val = E1(n)
    marker = iszero(val) ? "  <- PRIME" : ""
    @printf "  n=%2d  E1=%s%s\n" n val marker
end

println("\n-- Partition primality test vs trial division, n in [2, 100] --")
verify_range(2, 100)

println("\n-- Primes in [2, 50] via partition test --")
primes_found = filter(is_prime_partition, 2:50)
println("  ", primes_found)

println("\n-- Quasi-shuffle product: [1] * [1,1] --")
qs = quasishuffle_words([1], [1, 1])
for (w, c) in sort(collect(qs), by = p -> (length(p.first), p.first))
    println("  $(c) * M_$(w)")
end

println("\n-- Convolution verification for n=1:10 --")
for n in 1:10
    lhs = sum(M_macmahonesque([1], i, Rational{BigInt}) *
              M_macmahonesque([1], n - i, Rational{BigInt}) for i in 1:n-1; init=Rational{BigInt}(0))
    rhs = evaluate_zq(quasishuffle_words([1], [1]), n)
    status = lhs == rhs ? "OK" : "FAIL"
    @printf "  n=%2d  convolution=%s  [%s]\n" n lhs status
end

println("\n-- D operator on [1,1] --")
d_result = d_operator([1, 1])
println("  n*M_(1,1)(n) = ")
for (w, c) in sort(collect(d_result), by = p -> (length(p.first), p.first))
    println("    $(c) * M_$(w)")
end

println("\n-- Symmetrisation of [1,3] --")
sym = symmetrise([1, 3])
for (w, c) in sort(collect(sym), by = p -> p.first)
    println("  $(c) * U_$(w)")
end
