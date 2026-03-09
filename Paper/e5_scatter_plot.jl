"""
e5_scatter_plot.jl — Generate Figure 1: E₅(n) at composite numbers n ∈ [4, 200].

Run from PartitionPrimes root:
  julia --project=QuasiShuffleAlgebra Paper/e5_scatter_plot.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "QuasiShuffleAlgebra"))
using QuasiShuffleAlgebra
using CairoMakie

# ── Data ─────────────────────────────────────────────────────────────────────

const RANGE = 4:200

composites = filter(n -> !is_prime_trial(n), RANGE)
primes_in_range = filter(is_prime_trial, RANGE)

# Compute E5 values as Float64 for plotting (exact BigInt → Float64)
e5_vals = [Float64(E5(n)) for n in composites]

# Sign-classify each composite
pos_idx = findall(v -> v > 0, e5_vals)
neg_idx = findall(v -> v < 0, e5_vals)

pos_n   = composites[pos_idx]
neg_n   = composites[neg_idx]
pos_v   = e5_vals[pos_idx]
neg_v   = e5_vals[neg_idx]

println("Composites in [4,200]:  $(length(composites))")
println("E5 > 0: $(length(pos_idx))   E5 < 0: $(length(neg_idx))")
println("Range of E5 values: [$(minimum(e5_vals)), $(maximum(e5_vals))]")

# ── Plot ──────────────────────────────────────────────────────────────────────

# Symlog-like display: sign(v) * log10(1 + |v|) keeps zero at zero,
# compresses the huge range while preserving sign information.
symlog(v) = sign(v) * log10(1 + abs(v))

pos_sy = symlog.(pos_v)
neg_sy = symlog.(neg_v)

fig = Figure(size = (900, 520), fontsize = 13)

ax = Axis(fig[1, 1],
    title  = "Figure 1.  E₅(n) at composite integers n ∈ [4, 200]",
    xlabel = "n  (composite)",
    ylabel = "sign(E₅) · log₁₀(1 + |E₅(n)|)   [symlog scale]",
    titlesize  = 14,
    xlabelsize = 12,
    ylabelsize = 12,
)

# Zero reference line
hlines!(ax, [0.0], color = :black, linewidth = 0.8, linestyle = :dash)

# Positive composites — teal upward stems
for (n, sy) in zip(pos_n, pos_sy)
    lines!(ax, [n, n], [0.0, sy], color = (:teal, 0.6), linewidth = 1.0)
end
scatter!(ax, pos_n, pos_sy,
    color      = :teal,
    marker     = :circle,
    markersize = 6,
    label      = "E₅(n) > 0  ($(length(pos_idx)) composites)")

# Negative composites — crimson downward stems
for (n, sy) in zip(neg_n, neg_sy)
    lines!(ax, [n, n], [0.0, sy], color = (:crimson, 0.5), linewidth = 1.0)
end
scatter!(ax, neg_n, neg_sy,
    color      = :crimson,
    marker     = :circle,
    markersize = 6,
    label      = "E₅(n) < 0  ($(length(neg_idx)) composites)")

# Annotate a few notable positive composites
notable = [25, 35, 49, 77]
for n in notable
    idx = findfirst(==(n), composites)
    idx === nothing && continue
    sy = symlog(e5_vals[idx])
    text!(ax, n + 1.5, sy,
        text = "n=$n",
        fontsize = 9,
        color = :teal,
        align = (:left, :center))
end

# Legend
axislegend(ax, position = :lt, framevisible = true, labelsize = 11)

# Annotation box
text!(ax, 155, minimum(neg_sy) * 0.6,
    text = "Canonical E₅ (degree d = 2)\n349 of 404 composites in [4,500] are negative.\nE₅ + 856·E₄ ≥ 0 on all composites in [4,500].",
    fontsize  = 9,
    color     = :gray30,
    align     = (:left, :center),
    justification = :left)

save(joinpath(@__DIR__, "e5_scatter.png"), fig, px_per_unit = 2)
println("\nSaved: Paper/e5_scatter.png")
