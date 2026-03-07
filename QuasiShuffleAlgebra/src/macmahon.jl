# Divisor power sums
"""
    σ(k, n) -> BigInt

Sum of k-th powers of divisors of n. Uses BigInt arithmetic to avoid overflow:
for k >= 5 and n >= ~8000 individual terms d^k exceed Int64 range.
"""
function σ(k::Int, n::Int)::BigInt
    s = big(0)
    for d in 1:isqrt(n)
        if n % d == 0
            s += big(d)^k
            if d != n ÷ d
                s += big(n ÷ d)^k
            end
        end
    end
    return s
end

# Closed-form MacMahon functions
M1(n::Int) = σ(1, n)

function M2(n::Int)
    return ((-2n + 1) * σ(1, n) + σ(3, n)) // 8
end

function M3(n::Int)
    return ((40n^2 - 100n + 37) * σ(1, n) - 10(3n - 5) * σ(3, n) + 3 * σ(5, n)) // 1920
end

# Direct combinatorial M_a(n)
"""
    M_direct(a, n) -> Int

Compute M_a(n) combinatorially (eq. 1.2).
"""
function M_direct(a::Int, n::Int)::Int
    a == 0 && return n == 0 ? 1 : 0
    n <= 0 && return 0
    acc = Ref(0)
    _macmahon_rec!(a, n, 0, 1, 1, acc)
    return acc[]
end

function _macmahon_rec!(a, target, depth, min_s, prod_so_far, acc)
    if depth == a - 1
        for s in min_s:target
            if target % s == 0
                acc[] += prod_so_far * (target ÷ s)
            end
        end
        return
    end
    remaining_parts = a - depth
    for s in min_s:(target ÷ remaining_parts)
        for m in 1:(target ÷ s - (remaining_parts - 1))
            contrib = target - m * s
            contrib > 0 || continue
            _macmahon_rec!(a, contrib, depth + 1, s + 1, prod_so_far * m, acc)
        end
    end
end

# MacMahonesque functions M_{vec_a}(n) with type parameter
"""
    M_macmahonesque(vec_a, n, [T=Int]) -> T

Compute M_{vec_a}(n) for vec_a = (v1,...,va).
T controls accumulation type: Int (fast) or Rational{BigInt} (exact algebra).
"""
function M_macmahonesque(vec_a::Vector{Int}, n::Int, ::Type{T} = Int) where T
    a = length(vec_a)
    a == 0 && return n == 0 ? T(1) : T(0)
    n <= 0 && return T(0)
    acc = Ref(T(0))
    _macmahonesque_rec!(vec_a, a, n, 0, 1, T(1), acc)
    return acc[]
end

function _macmahonesque_rec!(vec_a, a, target, depth, min_s, wprod::T, acc) where T
    if depth == a - 1
        v = vec_a[depth + 1]
        for s in min_s:target
            if target % s == 0
                m = target ÷ s
                acc[] += wprod * T(m)^v
            end
        end
        return
    end
    v = vec_a[depth + 1]
    remaining_parts = a - depth
    for s in min_s:(target ÷ remaining_parts)
        for m in 1:(target ÷ s - (remaining_parts - 1))
            contrib = target - m * s
            contrib > 0 || continue
            _macmahonesque_rec!(vec_a, a, contrib, depth + 1, s + 1, wprod * T(m)^v, acc)
        end
    end
end
