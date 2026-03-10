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
    # Minimal-degree (d=2) canonical form extracted from the null space of the
    # prime evaluation matrix at (a_max=5, N=95 primes in [2,500]).
    # Vanishes iff n is prime; first appears at polynomial degree d=2 (not d=3).
    m3 = M3(n)
    m4 = M4(n)
    m5 = M5(n)
    bn = big(n)
    return (
        (-450450 + 675675*bn - 225225*bn^2) * M1(n)
        + (960960*bn - 120120*bn^2) * M2(n)
        + (2534912*bn - 166016*bn^2) * m3
        + (7999488*bn - 322560*bn^2) * m4
        + 258048000 * m5
    )
end

function E6(n::Int)
    # Degree d=4 expression from null space of prime evaluation matrix at a_max=7.
    # Vanishes iff n is prime. Involves M6 and M7 via τ-cancellation mechanism.
    m4  = M4(n)
    m5  = M5(n)
    m6  = M6(n)
    m7  = M7(n)
    bn  = big(n)
    return (
        (105367732470 - 277943208789*bn + 289613157079*bn^2
         - 145583618619*bn^3 + 28545937859*bn^4) * M1(n)
        + (-51324522800*bn^3 + 4292382480*bn^4) * M2(n)
        + (-40313554176*bn^3 + 1776519808*bn^4) * M3(n)
        + (-32888346624*bn^3 + 770273280*bn^4) * m4
        + (-7741440000*bn^3 - 154828800*bn^4) * m5
        + (-483548921856000 + 37196070912000*bn) * m6
        + 892705701888000 * m7
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
