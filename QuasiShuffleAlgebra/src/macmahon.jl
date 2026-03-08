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

"""
    ramanujan_tau(n) -> BigInt

Ramanujan tau function: coefficient of q^n in Δ(q) = q · ∏_{k≥1}(1-q^k)^24.
Computed on-demand; results are cached.
"""
const _tau_cache = BigInt[]

function ramanujan_tau(n::Int)::BigInt
    n <= length(_tau_cache) && return _tau_cache[n]
    N = max(n, length(_tau_cache) + 64)
    coeff = zeros(BigInt, N + 1)
    coeff[1] = 1
    for k in 1:N
        for _ in 1:24
            for m in N:-1:k
                coeff[m+1] -= coeff[m+1-k]
            end
        end
    end
    resize!(_tau_cache, N)
    for i in 1:N
        _tau_cache[i] = coeff[i]  # τ(i) = coeff of q^{i-1} in ∏ = coeff[i]
    end
    return _tau_cache[n]
end

# Closed-form MacMahon functions
M1(n::Int) = σ(1, n)

function M2(n::Int)
    return ((-2n + 1) * σ(1, n) + σ(3, n)) // 8
end

function M3(n::Int)
    return ((40n^2 - 100n + 37) * σ(1, n) - 10(3n - 5) * σ(3, n) + 3 * σ(5, n)) // 1920
end

"""
    M4(n) -> Rational{BigInt}

Closed-form M_4(n) via divisor power sums (derived by fitting to M_direct).
Formula: M_4 = P_1(n)·σ_1 + P_2(n)·σ_3 + P_3(n)·σ_5 + c·σ_7
where polynomial degrees follow the pattern deg(σ_{2j-1}) = a - j for a=4.
"""
function M4(n::Int)::Rational{BigInt}
    s1 = σ(1, n); s3 = σ(3, n); s5 = σ(5, n); s7 = σ(7, n)
    bn = big(n)
    return (
        (3229 // big(967680)) * s1
        + (-47 // big(4608)) * bn * s1
        + (7 // big(1152)) * bn^2 * s1
        + (-1 // big(1152)) * bn^3 * s1
        + (47 // big(9216)) * s3
        + (-7 // big(1536)) * bn * s3
        + (1 // big(1280)) * bn^2 * s3
        + (7 // big(15360)) * s5
        + (-1 // big(7680)) * bn * s5
        + (1 // big(193536)) * s7
    )
end

"""
    M5(n) -> Rational{BigInt}

Closed-form M_5(n) via divisor power sums (derived by fitting to M_direct).
"""
function M5(n::Int)::Rational{BigInt}
    s1 = σ(1, n); s3 = σ(3, n); s5 = σ(5, n); s7 = σ(7, n); s9 = σ(9, n)
    bn = big(n)
    return (
        (10679 // big(17203200)) * s1
        + (-1571 // big(774144)) * bn * s1
        + (133 // big(92160)) * bn^2 * s1
        + (-1 // big(3072)) * bn^3 * s1
        + (1 // big(46080)) * bn^4 * s1
        + (1571 // big(1548288)) * s3
        + (-133 // big(122880)) * bn * s3
        + (3 // big(10240)) * bn^2 * s3
        + (-1 // big(46080)) * bn^3 * s3
        + (133 // big(1228800)) * s5
        + (-1 // big(20480)) * bn * s5
        + (1 // big(215040)) * bn^2 * s5
        + (1 // big(516096)) * s7
        + (-1 // big(3096576)) * bn * s7
        + (1 // big(154828800)) * s9
    )
end

"""
    M6(n) -> Rational{BigInt}

Closed-form M_6(n) via divisor power sums and Ramanujan τ(n) (derived by fitting).
At weight 12, the cusp form Δ(q) contributes a τ(n) term.
"""
function M6(n::Int)::Rational{BigInt}
    s1  = σ(1, n);  s3  = σ(3, n);  s5  = σ(5, n)
    s7  = σ(7, n);  s9  = σ(9, n);  s11 = σ(11, n)
    tau = ramanujan_tau(n)
    bn  = big(n)
    return (
        (550499 // big(4541644800)) * s1
        + (-153617 // big(371589120)) * bn * s1
        + (2159 // big(6635520)) * bn^2 * s1
        + (-67 // big(737280)) * bn^3 * s1
        + (11 // big(1105920)) * bn^4 * s1
        + (-1 // big(2764800)) * bn^5 * s1
        + (153617 // big(743178240)) * s3
        + (-2159 // big(8847360)) * bn * s3
        + (67 // big(819200)) * bn^2 * s3
        + (-11 // big(1105920)) * bn^3 * s3
        + (1 // big(2580480)) * bn^4 * s3
        + (2159 // big(88473600)) * s5
        + (-67 // big(4915200)) * bn * s5
        + (11 // big(5160960)) * bn^2 * s5
        + (-1 // big(10321920)) * bn^3 * s5
        + (67 // big(123863040)) * s7
        + (-11 // big(74317824)) * bn * s7
        + (1 // big(111476736)) * bn^2 * s7
        + (11 // big(3715891200)) * s9
        + (-1 // big(3096576000)) * bn * s9
        + (1 // big(268995133440)) * s11
        + (-17 // big(150450048000)) * tau
    )
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
        v = vec_a[a - depth]   # vec_a[1]: weight of largest part (depth a-1)
        for s in min_s:target
            if target % s == 0
                m = target ÷ s
                acc[] += wprod * T(m)^v
            end
        end
        return
    end
    v = vec_a[a - depth]       # vec_a[a]: weight of smallest part (depth 0)
    remaining_parts = a - depth
    for s in min_s:(target ÷ remaining_parts)
        for m in 1:(target ÷ s - (remaining_parts - 1))
            contrib = target - m * s
            contrib > 0 || continue
            _macmahonesque_rec!(vec_a, a, contrib, depth + 1, s + 1, wprod * T(m)^v, acc)
        end
    end
end
