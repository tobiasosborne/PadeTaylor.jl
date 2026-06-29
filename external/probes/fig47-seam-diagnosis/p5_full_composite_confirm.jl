# P5 — FULL-COMPOSITE CONFIRMATION (bead vwgl): the acceptance test for #1.
#
# P4(a) showed a SINGLE window's local pole-count seed-variance collapses (|Δ|
# 9->2 near, 13->1 far). But a production windowed composite must (i) drop the
# GLOBAL |Δ| (monolithic = 172) toward single digits, and (ii) not reintroduce
# the seam at TILE BOUNDARIES. This builds a throwaway 5x5 overlapping composite
# and measures both.
#
# Scheme (emulates a production windowed_path_network_solve):
#   - 25 window centers on a 20-spaced 5x5 grid over [-40,40]^2.
#   - Each window solved WIDE: core [c±10] + margin 6 = [c±16], clamped to domain,
#     from the IC z=0, at a per-window seed threaded from the global seed
#     (window_seed = global_seed*1000 + window_index) so two global seeds genuinely
#     re-randomise every window (non-gameable).
#   - Keep each window's poles only in its Voronoi CORE (nearest-center == this
#     window) => hard partition, no double counting; margin gives interior nodes
#     to resolve core-edge poles. Union over windows = composite pole field.

using PadeTaylor, Printf
pI(z, u, up) = 6u^2 + z
const HALF = 50.0; const N = 101
xs = range(-HALF, HALF; length = N); ys = range(-HALF, HALF; length = N)
prob = PadeTaylorProblem(pI, (-0.1875, 0.3049), (0.0, HALF); order = 30)

const CENTERS = [(cx, cy) for cy in -40.0:20.0:40.0 for cx in -40.0:20.0:40.0]  # 25
const CORE_H = 10.0; const MARGIN = 6.0

nearest_center(p) = CENTERS[argmin([abs(complex(c[1], c[2]) - p) for c in CENTERS])]
restrict(ps, lo, hi) = filter(p -> lo <= real(p) <= hi && lo <= imag(p) <= hi, ps)

function poleset_agreement(A, B; tol = 0.5)
    isempty(A) && return (0.0, NaN, NaN)
    matched = 0; nn = Float64[]
    for p in A
        d = Inf; for q in B; dd = abs(p - q); dd < d && (d = dd); end
        push!(nn, d); d <= tol && (matched += 1)
    end
    sort!(nn); (matched / length(A), nn[(length(nn)+1)÷2], maximum(nn))
end

function composite_poles(seed)
    out = ComplexF64[]
    for (wi, (cx, cy)) in enumerate(CENTERS)
        xlo = max(-HALF, cx - CORE_H - MARGIN); xhi = min(HALF, cx + CORE_H + MARGIN)
        ylo = max(-HALF, cy - CORE_H - MARGIN); yhi = min(HALF, cy + CORE_H + MARGIN)
        zc = ComplexF64[complex(x, y) for y in ys for x in xs if xlo<=x<=xhi && ylo<=y<=yhi]
        pn = path_network_solve(prob, zc; h = 0.5, order = 30, rng_seed = seed*1000 + wi, extrapolate = true)
        for p in extract_poles(pn)
            nc = nearest_center(p)
            (nc[1] == cx && nc[2] == cy) && push!(out, p)
        end
    end
    out
end

# tile-boundary residual: poles within `band` of any inter-core line x or y in {-30,-10,10,30}
const BLINES = (-30.0, -10.0, 10.0, 30.0)
near_boundary(p; band = 1.5) = any(b -> abs(real(p)-b) < band, BLINES) || any(b -> abs(imag(p)-b) < band, BLINES)

println("="^72); println("MONOLITHIC baseline (|Δ|=172 expected)"); println("="^72); flush(stdout)
mono = Dict{Int,Vector{ComplexF64}}()
for seed in (0, 42)
    t0 = time()
    pn = path_network_solve(prob, ComplexF64[complex(x,y) for y in ys for x in xs];
                            h = 0.5, order = 30, rng_seed = seed, extrapolate = true)
    mono[seed] = restrict(extract_poles(pn), -49.0, 49.0)
    @printf("  mono seed=%d poles(in[-49,49])=%d  %.1fs\n", seed, length(mono[seed]), time()-t0); flush(stdout)
end

println("\n" * "="^72); println("WINDOWED COMPOSITE (5x5, overlap margin 6, per-window seed)"); println("="^72); flush(stdout)
comp = Dict{Int,Vector{ComplexF64}}()
for seed in (0, 42)
    t0 = time()
    comp[seed] = restrict(composite_poles(seed), -49.0, 49.0)
    @printf("  composite seed=%d poles(in[-49,49])=%d  %.1fs\n", seed, length(comp[seed]), time()-t0); flush(stdout)
end

mΔ = abs(length(mono[0]) - length(mono[42])); mmatch, mmed, mmax = poleset_agreement(mono[0], mono[42])
cΔ = abs(length(comp[0]) - length(comp[42])); cmatch, cmed, cmax = poleset_agreement(comp[0], comp[42])

# boundary vs interior breakdown of the composite seed-variance
cb0 = filter(near_boundary, comp[0]); cb42 = filter(near_boundary, comp[42])
ci0 = filter(!near_boundary, comp[0]); ci42 = filter(!near_boundary, comp[42])
bmatch, _, bmax = poleset_agreement(cb0, cb42); imatch, _, imax = poleset_agreement(ci0, ci42)

println("\n" * "="^72); println("RESULT — GLOBAL pole-set seed-invariance"); println("="^72)
@printf("MONOLITHIC : count s0=%d s42=%d  |Δ|=%d   match@0.5=%.1f%%  medianNN=%.3f  maxNN=%.2f\n",
        length(mono[0]), length(mono[42]), mΔ, 100mmatch, mmed, mmax)
@printf("COMPOSITE  : count s0=%d s42=%d  |Δ|=%d   match@0.5=%.1f%%  medianNN=%.3f  maxNN=%.2f\n",
        length(comp[0]), length(comp[42]), cΔ, 100cmatch, cmed, cmax)
println("\n-- TILE-BOUNDARY RESIDUAL (composite; near inter-core lines vs interior) --")
@printf("boundary band poles: s0=%d s42=%d  match@0.5=%.1f%%  maxNN=%.2f\n", length(cb0), length(cb42), 100bmatch, bmax)
@printf("interior     poles: s0=%d s42=%d  match@0.5=%.1f%%  maxNN=%.2f\n", length(ci0), length(ci42), 100imatch, imax)
println("""
READOUT:
  #1 CONFIRMED if composite |Δ| << 172 AND match >> monolithic AND the boundary
     band is NOT dramatically worse than the interior (overlap-core handled it).
  CAUTION if the boundary band match is much worse than interior => tile-boundary
     seam is the residual to engineer (more overlap / consensus assignment).
""")
