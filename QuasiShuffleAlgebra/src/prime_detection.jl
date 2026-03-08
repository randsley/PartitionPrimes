# Prime-detecting expressions (Theorem 1.1 and Table 1)

function E1(n::Int)
    return (n^2 - 3n + 2) * M1(n) - 8 * M2(n)
end

function E2(n::Int)
    return (3n^3 - 13n^2 + 18n - 8) * M1(n) +
           (12n^2 - 120n + 212) * M2(n) -
           960 * M3(n)
end

function E3(n::Int)
    m4 = M4(n)
    return (25n^4 - 171n^3 + 423n^2 - 447n + 170) * M1(n) +
           (300n^3 - 3554n^2 + 12900n - 14990) * M2(n) +
           (2400n^2 - 60480n + 214080) * M3(n) -
           725760 * m4
end

function E4(n::Int)
    m4 = M4(n)
    m5 = M5(n)
    return (126n^5 - 1303n^4 + 5073n^3 - 9323n^2 + 8097n - 2670) * M1(n) +
           (3024n^4 - 48900n^3 + 288014n^2 - 737100n + 695490) * M2(n) +
           (60480n^3 - 1510080n^2 + 10644480n - 23496480) * M3(n) +
           (725760n^2 - 36288000n + 218453760) * m4 -
           580608000 * m5
end

function E5(n::Int)
    m3 = M3(n)
    m4 = M4(n)
    m5 = M5(n)
    bn = big(n)
    return (
        (-270270 + 663549*bn - 522351*bn^2 + 129072*bn^3) * M1(n)
        + (-315272*bn^2 + 30400*bn^3) * M2(n)
        + (-340864*bn^2 + 15872*bn^3) * m3
        + (-193536*bn^2) * m4
        + 154828800 * m5
    )
end

function is_prime_trial(n::Int)::Bool
    n < 2 && return false
    n == 2 && return true
    n % 2 == 0 && return false
    for d in 3:2:isqrt(n)
        n % d == 0 && return false
    end
    return true
end

function is_prime_partition(n::Int)::Bool
    n < 2 && return false
    return iszero(E1(n))
end

function verify_range(lo::Int, hi::Int; verbose::Bool = true)
    mismatches = Int[]
    for n in lo:hi
        if is_prime_partition(n) != is_prime_trial(n)
            push!(mismatches, n)
        end
    end
    if verbose
        if isempty(mismatches)
            println("Passed: perfect agreement on [$lo, $hi], ",
                    count(is_prime_trial, lo:hi), " primes detected.")
        else
            println("FAILED: mismatches at ", mismatches)
        end
    end
    return mismatches
end
