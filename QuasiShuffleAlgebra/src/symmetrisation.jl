"""
    symmetrise(vec_a::Word) -> ZqElem

Compute the symmetrised series U^sym_{vec_a} = sum over unique permutations.
For all-odd weight vectors, the result is quasimodular (Theorem 4.4).
"""
function symmetrise(vec_a::Word)::ZqElem
    result = ZqElem()
    seen = Set{Word}()
    for perm in permutations(vec_a)
        w = Word(perm)
        if w in seen
            continue
        end
        push!(seen, w)
        result[w] = Rational{BigInt}(1)
    end
    return result
end
