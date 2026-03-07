# Type aliases
const Word = Vector{Int}
const ZqElem = Dict{Word, Rational{BigInt}}

"""Remove zero entries from a ZqElem in place."""
function cleanup!(z::ZqElem)
    filter!(p -> !iszero(p.second), z)
    return z
end

"""Add two ZqElem, returning a new ZqElem."""
function zq_add(a::ZqElem, b::ZqElem)::ZqElem
    result = copy(a)
    for (w, c) in b
        result[w] = get(result, w, Rational{BigInt}(0)) + c
    end
    cleanup!(result)
    return result
end

"""Scale a ZqElem by a rational coefficient."""
function zq_scale(c::Rational{BigInt}, z::ZqElem)::ZqElem
    iszero(c) && return ZqElem()
    result = ZqElem()
    for (w, v) in z
        result[w] = c * v
    end
    cleanup!(result)
    return result
end

zq_scale(c, z::ZqElem) = zq_scale(Rational{BigInt}(c), z)

"""Subtract: a - b."""
function zq_sub(a::ZqElem, b::ZqElem)::ZqElem
    return zq_add(a, zq_scale(Rational{BigInt}(-1), b))
end

"""
    all_words_up_to_weight(w; max_length=w+1, min_length=1)

Generate all words with non-negative integer entries summing to at most w,
with length between min_length and max_length.
"""
function all_words_up_to_weight(w::Int; max_length::Int = w + 1, min_length::Int = 1)
    words = Word[]
    for len in min_length:max_length
        _enum_words!(words, Int[], len, w)
    end
    return words
end

function _enum_words!(words, prefix, remaining_len, remaining_weight)
    if remaining_len == 0
        push!(words, copy(prefix))
        return
    end
    max_val = remaining_weight
    for v in 0:max_val
        push!(prefix, v)
        _enum_words!(words, prefix, remaining_len - 1, remaining_weight - v)
        pop!(prefix)
    end
end

"""
    evaluate_zq(z::ZqElem, n::Int) -> Rational{BigInt}

Evaluate a ZqElem at integer n by computing the linear combination of M_{word}(n).
"""
function evaluate_zq(z::ZqElem, n::Int)::Rational{BigInt}
    result = Rational{BigInt}(0)
    for (w, c) in z
        result += c * M_macmahonesque(w, n, Rational{BigInt})
    end
    return result
end
