const _bernoulli_cache = Dict{Int, Rational{BigInt}}(
    0 => big(1) // big(1),
    1 => big(-1) // big(2),
)

"""
    bernoulli(n) -> Rational{BigInt}

Compute the n-th Bernoulli number using the standard recurrence.
B_n = 0 for odd n > 1. Results are cached.
"""
function bernoulli(n::Int)::Rational{BigInt}
    n < 0 && throw(ArgumentError("n must be non-negative"))
    haskey(_bernoulli_cache, n) && return _bernoulli_cache[n]

    # B_n = 0 for odd n > 1
    if n > 1 && isodd(n)
        _bernoulli_cache[n] = Rational{BigInt}(0)
        return Rational{BigInt}(0)
    end

    # Recurrence: B_n = -1/(n+1) * sum_{k=0}^{n-1} C(n+1,k) * B_k
    s = Rational{BigInt}(0)
    for k in 0:(n-1)
        s += binomial(big(n + 1), big(k)) * bernoulli(k)
    end
    val = -s // (n + 1)
    _bernoulli_cache[n] = val
    return val
end
