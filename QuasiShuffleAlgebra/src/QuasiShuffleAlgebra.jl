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
include("conjecture.jl")

export Word, ZqElem
export zq_add, zq_scale, zq_multiply, cleanup!, evaluate_zq, all_words_up_to_weight
export bernoulli
export diamond_coeffs, diamond
export quasishuffle, quasishuffle_words
export σ, ramanujan_tau, M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, M_direct, M_macmahonesque
export cusp_form_16, cusp_form_18, cusp_form_20, cusp_form_22, cusp_form_24a, cusp_form_24b
export d_operator
export symmetrise
export E1, E2, E3, E4, E5, E6, is_prime_partition, verify_range, is_prime_trial
export build_basis, eval_matrix, eval_Ma, table1_coeffs
export rational_rref!, rational_nullspace, rank_over_Q, is_in_colspan
export test_conjecture, scan_conjecture

end
