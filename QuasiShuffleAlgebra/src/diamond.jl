"""
    diamond_coeffs(i, j) -> Dict{Int, Rational{BigInt}}

Compute the diamond product coefficients for letters i and j.
Returns a dictionary mapping output letter m -> coefficient.

The diamond product of letters i and j is:
  i diamond j = sum_m c_m * e_m

Based on eq. 4.3, with the convention that letter k corresponds to U_{(k)}.
The top letter is i+j+1 (not i+j) due to the index convention matching
MacMahonesque exponents to the Bachmann-Kuhn algebra generators.
"""
function diamond_coeffs(i::Int, j::Int)::Dict{Int, Rational{BigInt}}
    coeffs = Dict{Int, Rational{BigInt}}()

    # Top coefficient at letter i+j+1: i!*j! / (i+j+1)!
    top = factorial(big(i)) * factorial(big(j)) // factorial(big(i + j + 1))
    if !iszero(top)
        coeffs[i + j + 1] = top
    end

    # Lower coefficients at letters k = 1..i+j
    # c_k = [(-1)^i * C(i,k) + (-1)^j * C(j,k)] * B_{i+j+1-k} / (i+j+1-k)
    # (k is the output letter index; this is eq. 4.3 of the paper with k = m+1)
    # k=0 is always zero: same-parity gives B_{i+j+1}=0 (odd index>1); different-parity
    # gives (-1)^i + (-1)^j = 0. So the loop starts at k=1.
    for m in 1:(i + j)
        bern_idx = i + j + 1 - m
        B = bernoulli(bern_idx)
        iszero(B) && continue
        term1 = iseven(i) ? binomial(big(i), big(m)) : -binomial(big(i), big(m))
        term2 = iseven(j) ? binomial(big(j), big(m)) : -binomial(big(j), big(m))
        val = (term1 + term2) * B // bern_idx
        if !iszero(val)
            coeffs[m] = get(coeffs, m, Rational{BigInt}(0)) + val
        end
    end

    # Clean zeros
    filter!(p -> !iszero(p.second), coeffs)
    return coeffs
end

"""
    diamond(x::Int, y::Int) -> ZqElem

Diamond product of single letters x and y.
Returns a ZqElem (linear combination of single-letter words).
"""
function diamond(x::Int, y::Int)::ZqElem
    coeffs = diamond_coeffs(x, y)
    result = ZqElem()
    for (m, c) in coeffs
        result[[m]] = c
    end
    return result
end
