const _qs_cache = Dict{Tuple{Word, Word}, ZqElem}()

"""
    quasishuffle_words(u::Word, v::Word) -> ZqElem

Compute the quasi-shuffle product of two words.
Uses memoization for performance.

Recursion (eq. after 4.3):
  1 * w = w * 1 = w
  xu * yv = x(u*yv) + y(xu*v) + (x diamond y)(u*v)
"""
function quasishuffle_words(u::Word, v::Word)::ZqElem
    key = (u, v)
    haskey(_qs_cache, key) && return _qs_cache[key]

    result = _quasishuffle_impl(u, v)
    _qs_cache[key] = result
    return result
end

function _quasishuffle_impl(u::Word, v::Word)::ZqElem
    # Base cases: empty word acts as identity
    if isempty(u)
        return isempty(v) ? ZqElem(Int[] => Rational{BigInt}(1)) : ZqElem(v => Rational{BigInt}(1))
    end
    if isempty(v)
        return ZqElem(u => Rational{BigInt}(1))
    end

    x = u[1]
    u_tail = u[2:end]
    y = v[1]
    v_tail = v[2:end]

    result = ZqElem()

    # Term 1: x(u_tail * v)
    inner1 = quasishuffle_words(u_tail, v)
    for (w, c) in inner1
        new_word = vcat([x], w)
        result[new_word] = get(result, new_word, Rational{BigInt}(0)) + c
    end

    # Term 2: y(u * v_tail)
    inner2 = quasishuffle_words(u, v_tail)
    for (w, c) in inner2
        new_word = vcat([y], w)
        result[new_word] = get(result, new_word, Rational{BigInt}(0)) + c
    end

    # Term 3: (x diamond y)(u_tail * v_tail)
    dia = diamond(x, y)
    inner3 = quasishuffle_words(u_tail, v_tail)
    for (dw, dc) in dia
        # dw is a single-letter word [m]
        for (w, c) in inner3
            new_word = vcat(dw, w)
            result[new_word] = get(result, new_word, Rational{BigInt}(0)) + dc * c
        end
    end

    cleanup!(result)
    return result
end

"""
    zq_multiply(a::ZqElem, b::ZqElem) -> ZqElem

Multiply two ZqElem using the quasi-shuffle product (bilinear extension).
"""
function zq_multiply(a::ZqElem, b::ZqElem)::ZqElem
    result = ZqElem()
    for (wa, ca) in a
        for (wb, cb) in b
            prod = quasishuffle_words(wa, wb)
            for (w, c) in prod
                result[w] = get(result, w, Rational{BigInt}(0)) + ca * cb * c
            end
        end
    end
    cleanup!(result)
    return result
end

"""Clear the quasi-shuffle cache."""
function clear_qs_cache!()
    empty!(_qs_cache)
end
