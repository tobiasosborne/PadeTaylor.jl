# B4 Stage-2 interpolation audition (bead padetaylor-0ln.37.8).
#
# Throwaway probe: measures, on the ACTUAL P_I^(2) tritronquee wedge
# underlay, the audition metrics for
#   (a) nearest-node Voronoi   (current `_stage2_fill`)
#   (b) distance-weighted blend of all overlapping in-disc Pade evals.
#
# Metrics reported:
#   - overlap multiplicity distribution over covered grid pixels;
#   - C0 jump at former-Voronoi boundaries: |u| difference between
#     adjacent grid pixels assigned to DIFFERENT nearest nodes, both
#     covered, under (a) vs (b);
#   - overlap-disagreement magnitude: at a multiply-covered pixel, the
#     spread max-min of the participating nodes' Pade evaluations.
#
# Run: julia --project=. external/probes/stage2-blend-audition/audition.jl

using PadeTaylor
using PadeTaylor.VectorPathNetwork: vector_path_network_solve
using PadeTaylor.VectorPathNetworkStage2: _validity_radius
using Statistics

# ---------------------------------------------------------------------------
# Build the P_I^(2) tritronquee wedge walk, seeded directly from the KKG
# asymptotic IC at z = -3 (the SURF_Z_SEED of the figure helper).  We use
# the package-default B2 adaptive :max_q_root walk, order 24, h ceiling 0.1
# — the same numerics surf_wedge_fill drives.
# ---------------------------------------------------------------------------
const T_PARAM = 0.0
f = painleve_hierarchy(:I, 2; t = T_PARAM)
z_seed = -3.0 + 0.0im
y_seed = ComplexF64.(pI2_tritronquee_ic(z_seed; t = T_PARAM, n_terms = 2))
prob = VectorPadeTaylorProblem(f, y_seed, (z_seed, 20.0 + 0.0im); order = 24)

# The B3 extended threading fan: 19 radial shells x 9 angles in the wedge.
targets = ComplexF64[r * cis(θ) for r in 2.0:1.0:20.0
                                for θ in range(-0.5, 0.5; length = 9)]

println("Running P_I^(2) wedge walk ...")
walk = vector_path_network_solve(prob, targets; order = 24, h = 0.1,
                                 tol = 1.0e-8)
nnode = length(walk.visited_z)
println("  visited nodes: ", nnode)

# A Cartesian grid over the wedge region the nodes actually populate.
# A grid fine enough to RESOLVE the R_gate discs (med radius ~0.04): we
# zoom on a representative pole-rich sub-window of the wedge and sample at
# ~0.012 spacing so several pixels sit across each disc — only then is an
# "adjacent pixel" C0 jump a real boundary measurement.
zr = [real(z) for z in walk.visited_z]
zi = [imag(z) for z in walk.visited_z]
x0 = quantile(zr, 0.35); x1 = quantile(zr, 0.65)
y0 = quantile(zi, 0.30); y1 = quantile(zi, 0.70)
xs = range(x0, x1; length = 380)
ys = range(y0, y1; length = 380)
grid = ComplexF64[x + y*im for y in ys, x in xs]   # row-major (ys, xs)
gridv = vec(grid)
println("  grid: ", length(ys), " x ", length(xs), " = ", length(gridv), " pixels")

# ---------------------------------------------------------------------------
# Recompute every node's B1 R_gate (the same gate _stage2_fill applies).
# ---------------------------------------------------------------------------
tol = 1.0e-8
Rgate = Float64[ _validity_radius(walk.visited_jets[k],
                                  walk.visited_denominator[k],
                                  real(walk.visited_h[k]), tol)
                 for k in 1:nnode ]
println("  R_gate: med=", round(median(Rgate); digits=4),
        " p10=", round(quantile(Rgate,0.1); digits=4),
        " p90=", round(quantile(Rgate,0.9); digits=4))

# Per-node shared-Q Pade evaluation at z (low-to-high Horner).
function evalpoly_lh(c, t)
    s = zero(ComplexF64)
    @inbounds for k in length(c):-1:1
        s = s*t + c[k]
    end
    s
end
function node_eval(k, z)
    h_v = real(walk.visited_h[k])
    t   = (z - walk.visited_z[k]) / h_v
    q   = evalpoly_lh(walk.visited_denominator[k], t)
    # component 1 only (u = V_0) — the figure underlay channel.
    evalpoly_lh(walk.visited_numerators[k][1], t) / q
end

# For each grid pixel: the list of nodes whose verified disc covers it.
covering = Vector{Vector{Int}}(undef, length(gridv))
for (i, z) in enumerate(gridv)
    lst = Int[]
    for k in 1:nnode
        if abs(z - walk.visited_z[k]) <= Rgate[k]
            push!(lst, k)
        end
    end
    covering[i] = lst
end
mult = length.(covering)
covered = mult .> 0
println()
println("=== OVERLAP MULTIPLICITY (covered pixels) ===")
ncov = count(covered)
println("  covered pixels: ", ncov, " / ", length(gridv),
        "  (", round(100*ncov/length(gridv); digits=1), " %)")
if ncov > 0
    cm = mult[covered]
    println("  multiplicity: mean=", round(mean(cm); digits=2),
            " med=", median(cm), " max=", maximum(cm))
    for m in 1:min(maximum(cm), 8)
        println("    mult=", m, ": ", count(==(m), cm),
                " (", round(100*count(==(m),cm)/ncov; digits=1), " %)")
    end
    println("    mult>=2: ", round(100*count(>=(2),cm)/ncov; digits=1), " %")
end

# ---------------------------------------------------------------------------
# Scheme (a) — Voronoi: nearest node (lexicographic tiebreak as Stage 1).
# ---------------------------------------------------------------------------
function nearest_node(z)
    idx, best = 1, abs(walk.visited_z[1] - z)
    for k in 2:nnode
        d = abs(walk.visited_z[k] - z)
        if d < best
            idx, best = k, d
        end
    end
    idx
end
nn = [nearest_node(z) for z in gridv]
u_voronoi = Vector{ComplexF64}(undef, length(gridv))
for i in eachindex(gridv)
    # Voronoi underlay: assigned node must itself cover the pixel.
    if abs(gridv[i] - walk.visited_z[nn[i]]) <= Rgate[nn[i]]
        u_voronoi[i] = node_eval(nn[i], gridv[i])
    else
        u_voronoi[i] = ComplexF64(NaN, NaN)
    end
end

# ---------------------------------------------------------------------------
# Scheme (b) — distance-weighted blend.  Smooth partition-of-unity bump
# w_k(z) = (1 - (d_k/R_k)^2)^2  for d_k < R_k, else 0 — a C1 bump that
# vanishes (with zero derivative) at the disc edge.  Blend = sum w_k u_k
# / sum w_k over the COVERING nodes only (each strictly in its own disc).
# ---------------------------------------------------------------------------
function bump(d, R)
    x = d / R
    x >= 1 ? 0.0 : (1 - x*x)^2
end
u_blend = Vector{ComplexF64}(undef, length(gridv))
disagree = Float64[]            # spread among participating nodes
for i in eachindex(gridv)
    lst = covering[i]
    if isempty(lst)
        u_blend[i] = ComplexF64(NaN, NaN)
        continue
    end
    z = gridv[i]
    wsum = 0.0
    acc  = zero(ComplexF64)
    vals = ComplexF64[]
    for k in lst
        d = abs(z - walk.visited_z[k])
        w = bump(d, Rgate[k])
        # a pixel exactly at a disc edge gets w=0 from every node; fall
        # back to inverse-distance so it is still blended honestly.
        ev = node_eval(k, z)
        push!(vals, ev)
        wsum += w
        acc  += w * ev
    end
    if wsum > 0
        u_blend[i] = acc / wsum
    else
        # all-edge degenerate: equal-weight mean of in-disc evals.
        u_blend[i] = sum(vals) / length(vals)
    end
    if length(lst) >= 2
        mags = abs.(vals)
        push!(disagree, maximum(mags) - minimum(mags))
    end
end

println()
println("=== OVERLAP DISAGREEMENT (mult>=2 pixels, ||u_k| spread) ===")
if !isempty(disagree)
    println("  n=", length(disagree),
            " med=", round(median(disagree); sigdigits=3),
            " mean=", round(mean(disagree); sigdigits=3),
            " p90=", round(quantile(disagree,0.9); sigdigits=3),
            " max=", round(maximum(disagree); sigdigits=3))
    # relative to local |u| scale
    relmed = median(disagree) / median(filter(isfinite, abs.(u_blend)))
    println("  med disagreement / med |u| = ", round(relmed; sigdigits=3))
end

# ---------------------------------------------------------------------------
# C0 jump at former-Voronoi boundaries.  For each pair of horizontally /
# vertically adjacent grid pixels assigned to DIFFERENT nearest nodes, with
# BOTH covered, record ||u_left| - |u_right|| under (a) and under (b).
# ---------------------------------------------------------------------------
ny, nx = size(grid)
lidx(r, c) = (c - 1) * ny + r        # vec() is column-major over (ny,nx)
jump_vor = Float64[]
jump_bld = Float64[]
for r in 1:ny, c in 1:nx
    i = lidx(r, c)
    for (dr, dc) in ((0,1), (1,0))
        r2, c2 = r + dr, c + dc
        (r2 <= ny && c2 <= nx) || continue
        j = lidx(r2, c2)
        nn[i] == nn[j] && continue                 # same Voronoi cell
        # both pixels covered under Voronoi:
        if isfinite(u_voronoi[i]) && isfinite(u_voronoi[j])
            push!(jump_vor, abs(abs(u_voronoi[i]) - abs(u_voronoi[j])))
        end
        if isfinite(u_blend[i]) && isfinite(u_blend[j])
            push!(jump_bld, abs(abs(u_blend[i]) - abs(u_blend[j])))
        end
    end
end

println()
println("=== C0 JUMP at former-Voronoi boundaries (||u| adjacent diff) ===")
println("  (a) Voronoi : n=", length(jump_vor),
        " med=", round(median(jump_vor); sigdigits=3),
        " p90=", round(quantile(jump_vor,0.9); sigdigits=3),
        " max=", round(maximum(jump_vor); sigdigits=3))
println("  (b) blend   : n=", length(jump_bld),
        " med=", round(median(jump_bld); sigdigits=3),
        " p90=", round(quantile(jump_bld,0.9); sigdigits=3),
        " max=", round(maximum(jump_bld); sigdigits=3))

# Baseline: jump between SAME-cell adjacent covered pixels (the intrinsic
# field gradient at this grid spacing — the noise floor a boundary jump
# must be compared against).
jump_same = Float64[]
for r in 1:ny, c in 1:nx
    i = lidx(r, c)
    for (dr, dc) in ((0,1), (1,0))
        r2, c2 = r + dr, c + dc
        (r2 <= ny && c2 <= nx) || continue
        j = lidx(r2, c2)
        nn[i] == nn[j] || continue
        if isfinite(u_voronoi[i]) && isfinite(u_voronoi[j])
            push!(jump_same, abs(abs(u_voronoi[i]) - abs(u_voronoi[j])))
        end
    end
end
if !isempty(jump_same)
    println("  baseline (same-cell adjacent, intrinsic gradient):")
    println("    med=", round(median(jump_same); sigdigits=3),
            " p90=", round(quantile(jump_same,0.9); sigdigits=3))
    println("  => Voronoi boundary jump / intrinsic = ",
            round(median(jump_vor)/median(jump_same); sigdigits=3), "x (med)")
    println("  => blend   boundary jump / intrinsic = ",
            round(median(jump_bld)/median(jump_same); sigdigits=3), "x (med)")
end
# ---------------------------------------------------------------------------
# Genuine Voronoi-seam count: a C0 jump can only be SMOOTHED by a blend if
# the boundary lies inside an OVERLAP (both nodes cover both pixels).  Count
# adjacent-different-node covered pairs where the overlap is non-empty.
# ---------------------------------------------------------------------------
function count_seams()
    seam_in_overlap = 0
    seam_total      = 0
    for r in 1:ny, c in 1:nx
        i = lidx(r, c)
        for (dr, dc) in ((0,1), (1,0))
            r2, c2 = r + dr, c + dc
            (r2 <= ny && c2 <= nx) || continue
            j = lidx(r2, c2)
            nn[i] == nn[j] && continue
            (isfinite(u_voronoi[i]) && isfinite(u_voronoi[j])) || continue
            seam_total += 1
            a, b = nn[i], nn[j]
            in_both(p) = abs(gridv[p]-walk.visited_z[a]) <= Rgate[a] &&
                         abs(gridv[p]-walk.visited_z[b]) <= Rgate[b]
            if in_both(i) && in_both(j)
                seam_in_overlap += 1
            end
        end
    end
    return seam_total, seam_in_overlap
end
seam_total, seam_in_overlap = count_seams()
println()
println("=== BLENDABLE SEAMS (boundary inside a true overlap) ===")
println("  covered Voronoi-boundary adjacent pairs : ", seam_total)
println("  ... of which lie inside a node-pair overlap (blendable): ",
        seam_in_overlap,
        seam_total > 0 ? "  ("*string(round(100*seam_in_overlap/seam_total;
                                            digits=1))*" %)" : "")
println()
println("audition complete.")
