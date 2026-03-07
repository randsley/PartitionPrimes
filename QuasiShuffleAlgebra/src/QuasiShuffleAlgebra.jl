"""
QuasiShuffleAlgebra.jl

Quasi-shuffle algebra for MacMahonesque partition functions and prime detection.

Based on: "Integer Partitions Detect the Primes"
  Craig, van Ittersum & Ono (arXiv:2405.06451v2, 2024)
"""
module QuasiShuffleAlgebra

using Combinatorics
using LinearAlgebra
using Printf

include("util.jl")
include("bernoulli.jl")
include("diamond.jl")
include("quasishuffle.jl")
include("macmahon.jl")
include("d_operator.jl")
include("symmetrisation.jl")
include("prime_detection.jl")

export Word, ZqElem
export zq_add, zq_scale, zq_multiply, cleanup!, evaluate_zq, all_words_up_to_weight
export bernoulli
export diamond_coeffs, diamond
export quasishuffle, quasishuffle_words
export σ, M1, M2, M3, M_direct, M_macmahonesque
export d_operator
export symmetrise
export E1, E2, E3, E4, is_prime_partition, verify_range, is_prime_trial

end
