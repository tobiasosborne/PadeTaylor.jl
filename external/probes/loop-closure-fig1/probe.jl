# external/probes/loop-closure-fig1/probe.jl
#
# Loop-closure probe — Figure 1 of FFW 2017 (P_III three-sheet pole spiral).
# Bead `padetaylor-e3h`.  Investigation only; no `src/` changes.
#
# ## Hypothesis
#
# The visible "seam" in `figures/output/ffw2017_fig_1.png` (white wedges
# between the spiral halves on sheets ±1 — and, to a lesser extent on
# sheet 0) might be caused by **path-network truncation drift**: the
# Stage-1 walker builds a *tree*, so any node `k` is reached by a unique
# IC→k path.  Two nodes A and B that are geometrically adjacent (a
# Delaunay edge in ζ-space) but tree-distant (LCA deep in the IC) carry
# *independent* accumulated truncation error.  When Stage 2 evaluates
# `w(M)` at the edge midpoint by snapping to the nearest visited node,
# the choice of "nearest" can flip across the Voronoi cell wall, and
# the discontinuity = |u_A(M) − u_B(M)| is what we measure here.
#
# Concretely, for every non-tree Delaunay edge (A, B) on the principal
# sheet:
#
#     M   = (z_A + z_B) / 2
#     t_A = (M − z_A) / h_A,   u_A = _evaluate_pade(visited_pade[A], t_A)
#     t_B = (M − z_B) / h_B,   u_B = _evaluate_pade(visited_pade[B], t_B)
#     ΔP  = |u_A − u_B|,   ΔP_rel = ΔP / (|u_A| + |u_B| + ε)
#
# If `ΔP_rel ≫ adaptive_tol (1e-10)` and the high-disagreement edges
# **cluster at the rendered seam locations** in the z-plane projection,
# the hypothesis is **Confirmed** — graph-consensus Stage-2 (averaging
# overlapping Padé patches with weights based on tree distance and Padé
# disc geometry) is the right v2 direction.
#
# If high disagreement exists but doesn't cluster at the seam, the
# verdict is **Partially confirmed** — consensus signal exists but the
# seam itself is from something else (bilinear interp at the sheet lift,
# Cartesian-grid undersampling near pole-dense neighbourhoods, etc.).
#
# If all `ΔP_rel ≲ adaptive_tol`, the verdict is **Null** — the walker
# tree is already self-consistent to the controller's tolerance, and the
# seam is a downstream rendering artefact (the Stage-1/2 boundary is
# not where the error lives).
#
# All three outcomes are useful.  This probe does not bias toward
# "Confirmed".  See `REPORT.md` in the same directory for the verdict.
#
# ## Algorithm
#
# 1. Reproduce `figures/ffw2017_fig_1.jl`'s Stage-1 solve verbatim.
#    Use the same R(ζ), ζ-window, adaptive controller, tolerance, and
#    `build_stage1_targets` rectangular-raster constructor as the
#    figure script.  Cite: `figures/ffw2017_fig_1.jl:130-241`.
# 2. Filter `visited_z` to the principal sheet (sheet [0]).  Because
#    the Fig 1 solve runs *without* `branch_points`,
#    `sol.visited_sheet[k]` is the empty `Int[]` for every k — so we
#    fall back to the ζ-strip predicate `-2π < imag(visited_z[k]) ≤ 2π`
#    (FFW md:103 strip definition).
# 3. Build the Delaunay triangulation on the sheet-0 ζ-coordinates.
#    `DelaunayTriangulation.each_edge` returns unique undirected edges
#    as `Tuple{Int, Int}` referring to the input-points array.
# 4. Tree edges: `(visited_parent[k], k)` for every k ≥ 2 (skipping
#    `visited_parent[k] == 0`).  We restrict to edges where BOTH
#    endpoints are sheet 0 (mapping the tree-edge indices through the
#    sheet-0 sub-indexing).
# 5. Non-tree edges = `each_edge(tri) \ tree_edges_set`.  For each
#    non-tree edge, evaluate the two-sided Padé disagreement at the
#    midpoint.  Record `(A, B, ΔP_abs, ΔP_rel, tree_dist, extrap_max)`.
# 6. Render the histogram and spatial map; write the verdict to
#    `REPORT.md`.
#
# Re-run:
#     julia --project=external/probes/loop-closure-fig1 \
#           external/probes/loop-closure-fig1/probe.jl
#
# Single Julia invocation (CLAUDE.md Rule 7).  Wall ~ 1–3 min depending
# on machine; the bulk is the Stage-1 walk.

using PadeTaylor
using PadeTaylor.CoordTransforms: pIII_transformed_rhs, pIII_z_to_ζ
using DelaunayTriangulation
using CairoMakie
using Statistics
using Printf

const PROBE_DIR = @__DIR__
const OUTPUT_DIR = joinpath(PROBE_DIR, "output")
mkpath(OUTPUT_DIR)

# ──────────────────────────────────────────────────────────────────────
# Step 1 — reproduce Fig 1's Stage-1 solve verbatim.
# Source: figures/ffw2017_fig_1.jl:130-241.
# ──────────────────────────────────────────────────────────────────────
const α, β, γ, δ = -0.5, -0.5, 1.0, -1.0
const Z0, U0, UP0 = 1.0 + 0.0im, 0.25 + 0.0im, 1.0 + 0.0im
const ORDER       = 30
const ADAPT_TOL   = 1.0e-10
const K_CONS      = 1.0e-3
const MAX_RESC    = 50

const Re_LO, Re_HI = -2.0, 5.0
const Im_LO, Im_HI = -6π + 0.5, 6π - 0.5

R(ζ) = max(0.1, (8.0 - real(ζ)) / 20.0)

function build_stage1_targets()
    tgs = ComplexF64[]
    rev = Re_LO
    while rev ≤ Re_HI + 1e-12
        rh = R(rev + 0.0im)
        imv = Im_LO
        while imv ≤ Im_HI + 1e-12
            push!(tgs, complex(rev, imv))
            imv += rh
        end
        rev += rh
    end
    return tgs
end

ζ0, w0, wp0 = pIII_z_to_ζ(Z0, U0, UP0)
println("=" ^ 78)
println("Loop-closure probe — FFW 2017 Fig 1 (bead padetaylor-e3h)")
println("=" ^ 78)
@printf("IC: (z, u, u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("    → (ζ, w, w') = (%s, %s, %s)\n", ζ0, w0, wp0)
@printf("R(ζ) floor at 0.1; ζ-window Re ∈ [%g, %g], Im ∈ [%g, %g]\n",
        Re_LO, Re_HI, Im_LO, Im_HI)

f = pIII_transformed_rhs(α, β, γ, δ)
prob = PadeTaylorProblem(f, (w0, wp0),
                          (ζ0, complex(Re_HI, Im_HI));
                          order = ORDER)

stage1_targets = build_stage1_targets()
@printf("Stage-1 targets: %d  (FFW md:95 reports 3041)\n",
        length(stage1_targets))

t_walk0 = time()
sol = path_network_solve(prob, stage1_targets;
                          h = R(ζ0),
                          node_separation = R,
                          step_size_policy = :adaptive_ffw,
                          adaptive_tol = ADAPT_TOL,
                          k_conservative = K_CONS,
                          max_rescales = MAX_RESC,
                          max_steps_per_target = 8000)
walk_wall = time() - t_walk0
N = length(sol.visited_z)
@printf("Stage-1 walk: %.2f s; visited tree nodes = %d (worklog 037: ~2664)\n",
        walk_wall, N)

# ──────────────────────────────────────────────────────────────────────
# Step 2 — sheet-0 filter via ζ-strip predicate.  Branched solve is OFF
# (no branch_points), so visited_sheet[k] == Int[] for every k.  We use
# the FFW md:103 strip definition: sheet 0 ⇔ Im ζ ∈ (-2π, 2π].
# ──────────────────────────────────────────────────────────────────────
function sheet0_mask(visited_z)
    return [(-2π < imag(z) ≤ 2π) for z in visited_z]
end

mask0 = sheet0_mask(sol.visited_z)
sheet0_idx = findall(mask0)            # indices into visited_*
n_sheet0 = length(sheet0_idx)
@printf("Sheet-0 nodes (Im ζ ∈ (-2π, 2π]): %d / %d\n", n_sheet0, N)
if n_sheet0 < 200
    @warn "Sheet-0 set < 200; widening predicate per honest-deviation latitude"
    mask0 = [(-2π - 0.5 < imag(z) ≤ 2π + 0.5) for z in sol.visited_z]
    sheet0_idx = findall(mask0)
    n_sheet0 = length(sheet0_idx)
    @printf("After widening: %d sheet-0 nodes\n", n_sheet0)
end

# Map global k ∈ 1..N ↦ sheet-0 local idx (0 if not in sheet 0).
global_to_local = zeros(Int, N)
for (loc, glob) in enumerate(sheet0_idx)
    global_to_local[glob] = loc
end

ζ_sheet0 = sol.visited_z[sheet0_idx]
pts2D    = [(real(z), imag(z)) for z in ζ_sheet0]

# ──────────────────────────────────────────────────────────────────────
# Step 3 — Delaunay triangulation on sheet-0 ζ-coordinates.
# ──────────────────────────────────────────────────────────────────────
@printf("Triangulating %d sheet-0 nodes ...\n", n_sheet0)
t_tri0 = time()
tri = triangulate(pts2D)
@printf("DelaunayTriangulation built in %.2f s\n", time() - t_tri0)

# Normalise undirected edges to (min, max) so set membership is symmetric.
norm_edge(a, b) = a < b ? (a, b) : (b, a)

delaunay_edges = Set{Tuple{Int,Int}}()
n_ghost_ref = Ref(0)
for e in each_edge(tri)
    a, b = e
    # DelaunayTriangulation.jl uses negative indices for "ghost"
    # vertices on the convex-hull boundary (the unbounded-face
    # representative).  Skip any edge incident to a ghost vertex.
    if a < 1 || b < 1
        n_ghost_ref[] += 1
        continue
    end
    push!(delaunay_edges, norm_edge(a, b))
end
n_ghost = n_ghost_ref[]
@printf("Unique Delaunay edges (sheet 0): %d  (ghost edges skipped: %d; expected ~%d ≈ 3N)\n",
        length(delaunay_edges), n_ghost, 3 * n_sheet0)

# ──────────────────────────────────────────────────────────────────────
# Step 4 — tree edges, restricted to sheet 0.  visited_parent[k] = 0 is
# the IC sentinel.  We keep only edges where both endpoints are sheet 0
# (after mapping to local sheet-0 indices).
# ──────────────────────────────────────────────────────────────────────
tree_edges_local = Set{Tuple{Int,Int}}()
for k in 2:N
    p = sol.visited_parent[k]
    p == 0 && continue
    lk = global_to_local[k]
    lp = global_to_local[p]
    (lk == 0 || lp == 0) && continue
    push!(tree_edges_local, norm_edge(lk, lp))
end
@printf("Tree edges (sheet 0 ∩ sheet 0): %d\n", length(tree_edges_local))

# ──────────────────────────────────────────────────────────────────────
# Step 5 — non-tree edges = Delaunay minus tree.
# ──────────────────────────────────────────────────────────────────────
nontree_edges = setdiff(delaunay_edges, tree_edges_local)
@printf("Non-tree Delaunay edges: %d\n", length(nontree_edges))

# ──────────────────────────────────────────────────────────────────────
# Helper — tree-distance via LCA on parent chains.  visited_parent is a
# Vector{Int} indexed by GLOBAL node index; root has parent = 0.  We
# build a depth array once and use the classical LCA-via-equalise-depth.
# ──────────────────────────────────────────────────────────────────────
function build_depths(parent::Vector{Int})
    n = length(parent)
    depth = fill(-1, n)
    # Roots: any node with parent = 0.
    for k in 1:n
        parent[k] == 0 && (depth[k] = 0)
    end
    # Iteratively fill (parent's depth + 1) when parent already filled.
    progress = true
    while progress
        progress = false
        for k in 1:n
            depth[k] >= 0 && continue
            p = parent[k]
            if p >= 1 && depth[p] >= 0
                depth[k] = depth[p] + 1
                progress = true
            end
        end
    end
    return depth
end

depth = build_depths(sol.visited_parent)

function tree_path_distance(parent, depth, a::Int, b::Int)
    da, db = depth[a], depth[b]
    steps = 0
    while da > db
        a = parent[a]; da -= 1; steps += 1
    end
    while db > da
        b = parent[b]; db -= 1; steps += 1
    end
    while a != b
        a = parent[a]; b = parent[b]; steps += 2
    end
    return steps
end

# ──────────────────────────────────────────────────────────────────────
# Step 6 — evaluate the loop-closure disagreement on each non-tree edge.
# NO disc-cap: we want raw |u_A − u_B|, even past |t| = 1 (where the
# Padé is technically extrapolating).
# ──────────────────────────────────────────────────────────────────────
const EPS_FLOOR = 1e-300

mutable struct EdgeReport
    A_local::Int          # local sheet-0 index
    B_local::Int
    A_global::Int         # global visited-tree index
    B_global::Int
    ΔP_abs::Float64
    ΔP_rel::Float64
    tree_dist::Int
    extrap_max::Float64
    M::ComplexF64         # midpoint in ζ
end

@printf("Evaluating %d non-tree edges ...\n", length(nontree_edges))
t_eval0 = time()

results = EdgeReport[]
n_blowups_ref = Ref(0)
for (la, lb) in nontree_edges
    ga = sheet0_idx[la]
    gb = sheet0_idx[lb]
    zA = sol.visited_z[ga]
    zB = sol.visited_z[gb]
    hA = sol.visited_h[ga]
    hB = sol.visited_h[gb]
    M  = (zA + zB) / 2
    tA = (M - zA) / hA
    tB = (M - zB) / hB
    uA = try
        PadeTaylor.PathNetwork._evaluate_pade(sol.visited_pade[ga], tA)
    catch err
        n_blowups_ref[] += 1
        complex(NaN, NaN)
    end
    uB = try
        PadeTaylor.PathNetwork._evaluate_pade(sol.visited_pade[gb], tB)
    catch err
        n_blowups_ref[] += 1
        complex(NaN, NaN)
    end
    (isfinite(real(uA)) && isfinite(real(uB))) || continue
    ΔP_abs = abs(uA - uB)
    denom  = abs(uA) + abs(uB) + EPS_FLOOR
    ΔP_rel = ΔP_abs / denom
    td     = tree_path_distance(sol.visited_parent, depth, ga, gb)
    em     = max(abs(tA), abs(tB))
    push!(results, EdgeReport(la, lb, ga, gb, ΔP_abs, ΔP_rel, td, em, M))
end
n_blowups = n_blowups_ref[]
eval_wall = time() - t_eval0
@printf("Eval done in %.2f s; %d good edges, %d blowups (Padé denom=0)\n",
        eval_wall, length(results), n_blowups)

# ──────────────────────────────────────────────────────────────────────
# Aggregate statistics.
# ──────────────────────────────────────────────────────────────────────
ΔP_rels   = [r.ΔP_rel for r in results]
ΔP_abss   = [r.ΔP_abs for r in results]
tree_dst  = [r.tree_dist for r in results]
extrap_mx = [r.extrap_max for r in results]

med  = median(ΔP_rels)
p90  = quantile(ΔP_rels, 0.90)
p99  = quantile(ΔP_rels, 0.99)
maxR = maximum(ΔP_rels)

# Pearson correlations.
function pearson(xs, ys)
    mx, my = mean(xs), mean(ys)
    num = sum((xs .- mx) .* (ys .- my))
    den = sqrt(sum((xs .- mx).^2) * sum((ys .- my).^2))
    return den == 0 ? NaN : num / den
end
corr_td  = pearson(Float64.(tree_dst), log10.(max.(ΔP_rels, 1e-300)))
corr_em  = pearson(extrap_mx,          log10.(max.(ΔP_rels, 1e-300)))

@printf("\n--- ΔP_rel statistics over %d edges ---\n", length(results))
@printf("median:  %.3e\n", med)
@printf("p90:     %.3e\n", p90)
@printf("p99:     %.3e\n", p99)
@printf("max:     %.3e\n", maxR)
@printf("Pearson r(tree_dist,   log10 ΔP_rel) = %+.3f\n", corr_td)
@printf("Pearson r(extrap_max,  log10 ΔP_rel) = %+.3f\n", corr_em)

# ──────────────────────────────────────────────────────────────────────
# Step 7 — histogram (log-x) of ΔP_rel.
# ──────────────────────────────────────────────────────────────────────
fig1 = Figure(size = (900, 600))
ax1  = Axis(fig1[1, 1];
            xlabel = "log10(ΔP_rel)",
            ylabel = "count",
            title  = "Loop-closure ΔP_rel — Fig 1 sheet 0, $(length(results)) non-tree edges",
            subtitle = @sprintf("median=%.2e  p90=%.2e  p99=%.2e  max=%.2e",
                                med, p90, p99, maxR))
logΔ = log10.(max.(ΔP_rels, 1e-20))
hist!(ax1, logΔ; bins = 80, color = (:steelblue, 0.8),
      strokecolor = :black, strokewidth = 0.4)
# Reference lines.
for (val, label, col) in (
        (log10(1e-10), "adaptive_tol = 1e-10", :red),
        (log10(1e-6),  "FFW Exp-2 ~ 1e-6",    :orange),
        (log10(4e-15), "sheet-0 conj-sym median 4e-15", :darkgreen),
        (log10(med),   "this run's median",   :purple))
    vlines!(ax1, [val]; color = col, linewidth = 1.5, linestyle = :dash,
            label = label)
end
axislegend(ax1; position = :rt, labelsize = 9)
save(joinpath(OUTPUT_DIR, "histogram.png"), fig1)
@printf("Wrote %s\n", joinpath(OUTPUT_DIR, "histogram.png"))

# ──────────────────────────────────────────────────────────────────────
# Step 8 — spatial map.  Left axis: ζ-plane.  Right axis: z-plane via
# z = exp(ζ/2) (PIII transform).  Sheet-0 scatter dots; non-tree edges
# colored by log10(ΔP_rel); top-10 worst edges highlighted red.
# ──────────────────────────────────────────────────────────────────────
# Order edges by descending ΔP_rel for stable top-K extraction.
order_by_rel = sortperm(ΔP_rels; rev = true)
top10_idx    = order_by_rel[1:min(10, length(order_by_rel))]

# Edge segments for linesegments!.
function edge_segments_ζ(results, idxs)
    segs = Point2f[]
    cols = Float64[]
    for i in idxs
        r = results[i]
        zA = sol.visited_z[r.A_global]
        zB = sol.visited_z[r.B_global]
        push!(segs, Point2f(real(zA), imag(zA)))
        push!(segs, Point2f(real(zB), imag(zB)))
        push!(cols, log10(max(r.ΔP_rel, 1e-20)))
    end
    return segs, cols
end

function edge_segments_z(results, idxs)
    segs = Point2f[]
    cols = Float64[]
    for i in idxs
        r = results[i]
        # z = exp(ζ/2), but on sheet 0 we keep the principal-branch
        # representative (Im ζ ∈ (-2π, 2π] ⇒ arg z ∈ (-π, π]).
        zA_ζ = sol.visited_z[r.A_global]
        zB_ζ = sol.visited_z[r.B_global]
        zA_z = exp(zA_ζ / 2)
        zB_z = exp(zB_ζ / 2)
        push!(segs, Point2f(real(zA_z), imag(zA_z)))
        push!(segs, Point2f(real(zB_z), imag(zB_z)))
        push!(cols, log10(max(r.ΔP_rel, 1e-20)))
    end
    return segs, cols
end

all_idx = 1:length(results)
ζ_segs, ζ_cols = edge_segments_ζ(results, all_idx)
z_segs, z_cols = edge_segments_z(results, all_idx)
ζ_top, _       = edge_segments_ζ(results, top10_idx)
z_top, _       = edge_segments_z(results, top10_idx)

fig2 = Figure(size = (1500, 800))
Label(fig2[0, 1:2],
      "Loop-closure disagreement (sheet 0). Color = log10(ΔP_rel). Red = top-10 worst.";
      fontsize = 14)

axζ = Axis(fig2[1, 1];
           xlabel = "Re ζ", ylabel = "Im ζ",
           title = "ζ-plane (sheet 0)",
           aspect = DataAspect(),
           limits = (Re_LO - 0.2, Re_HI + 0.2, -2π - 0.3, 2π + 0.3))

# Sheet-0 node scatter as a context layer.
scatter!(axζ,
         Float64.(real.(ζ_sheet0)),
         Float64.(imag.(ζ_sheet0));
         color = (:black, 0.25), markersize = 1.5)

# Colormap range — clamp to [-15, 0] for readability.
crange = (-15.0, max(0.0, maximum(ζ_cols)))
linesegments!(axζ, ζ_segs;
              color = repeat(ζ_cols, inner = 2),
              colormap = :plasma, colorrange = crange,
              linewidth = 0.6)
linesegments!(axζ, ζ_top; color = :red, linewidth = 2.0)

# z-plane bounds: |z|_max = exp(Re_HI / 2) ≈ 12.2; pad 5%.
const Z_HALF = 1.05 * exp(Re_HI / 2)
axz = Axis(fig2[1, 2];
           xlabel = "Re z", ylabel = "Im z",
           title = "z-plane via z = exp(ζ/2)  (sheet 0)",
           aspect = DataAspect(),
           limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))

z_sheet0_pts = [exp(z / 2) for z in ζ_sheet0]
scatter!(axz,
         Float64.(real.(z_sheet0_pts)),
         Float64.(imag.(z_sheet0_pts));
         color = (:black, 0.25), markersize = 1.5)

linesegments!(axz, z_segs;
              color = repeat(z_cols, inner = 2),
              colormap = :plasma, colorrange = crange,
              linewidth = 0.6)
linesegments!(axz, z_top; color = :red, linewidth = 2.0)

Colorbar(fig2[1, 3];
         colormap = :plasma, limits = crange,
         label = "log10(ΔP_rel)")

save(joinpath(OUTPUT_DIR, "spatial_map.png"), fig2)
@printf("Wrote %s\n", joinpath(OUTPUT_DIR, "spatial_map.png"))

# ──────────────────────────────────────────────────────────────────────
# Step 9 — verdict.  Cutoffs (per probe plan):
#   • Null:                 max(ΔP_rel) ≤ adaptive_tol (1e-10).
#   • Confirmed:            top edges (above 1e-6) spatially cluster
#                           near rendered seam.
#   • Partially confirmed:  high disagreement exists (above adaptive_tol)
#                           but doesn't cluster.
#   • Inconclusive:         ambiguous spatial pattern.
#
# The "spatially cluster at seam" judgement requires visual inspection
# of `spatial_map.png` against `figures/output/ffw2017_fig_1.png`; this
# script computes the *quantitative* part and surfaces enough numbers
# for the operator (and the report) to make the call.  In particular we
# print: how many edges exceed 1e-6 (FFW Exp-2 reference) and 1e-3
# (visibly bad), and how those edges distribute in ζ-space (binned
# by Im ζ strip half and by Re ζ near vs far from the high-pole-density
# right boundary).
# ──────────────────────────────────────────────────────────────────────
n_above_tol = count(>(ADAPT_TOL), ΔP_rels)
n_above_1e6 = count(>(1e-6), ΔP_rels)
n_above_1e3 = count(>(1e-3), ΔP_rels)
@printf("\n--- edge-count tallies ---\n")
@printf("ΔP_rel > 1e-10 (adaptive_tol):  %5d / %d (%.2f%%)\n",
        n_above_tol, length(results), 100*n_above_tol/length(results))
@printf("ΔP_rel > 1e-6 (FFW Exp-2):     %5d / %d (%.2f%%)\n",
        n_above_1e6, length(results), 100*n_above_1e6/length(results))
@printf("ΔP_rel > 1e-3 (visibly bad):   %5d / %d (%.2f%%)\n",
        n_above_1e3, length(results), 100*n_above_1e3/length(results))

# Spatial clustering: of top-K worst edges, where do they sit?
@printf("\n--- top-10 worst non-tree edges ---\n")
@printf("%4s %4s %12s %12s %10s %10s %20s\n",
        "A_g", "B_g", "ΔP_abs", "ΔP_rel", "tree_d", "extrap_max", "midpoint ζ")
for k in top10_idx
    r = results[k]
    @printf("%4d %4d %12.3e %12.3e %10d %10.3f %20s\n",
            r.A_global, r.B_global,
            r.ΔP_abs, r.ΔP_rel, r.tree_dist, r.extrap_max,
            string(round(r.M; digits = 3)))
end

# Re ζ and Im ζ distribution of "visibly bad" edges (> 1e-6).
bad = filter(r -> r.ΔP_rel > 1e-6, results)
if !isempty(bad)
    re_mids = [real(r.M) for r in bad]
    im_mids = [imag(r.M) for r in bad]
    @printf("\nBad-edge (ΔP_rel > 1e-6) midpoint distribution (N=%d):\n",
            length(bad))
    @printf("  Re ζ:  min=%.2f  median=%.2f  max=%.2f\n",
            minimum(re_mids), median(re_mids), maximum(re_mids))
    @printf("  Im ζ:  min=%.2f  median=%.2f  max=%.2f\n",
            minimum(im_mids), median(im_mids), maximum(im_mids))
    @printf("  Fraction with Re ζ > 3.0 (near pole-dense boundary): %.2f%%\n",
            100 * count(>(3.0), re_mids) / length(bad))
    @printf("  Fraction with |Im ζ| > π (near sheet boundary):      %.2f%%\n",
            100 * count(>(Float64(π)), abs.(im_mids)) / length(bad))
end

verdict = if maxR ≤ ADAPT_TOL
    "Null"
elseif n_above_1e6 == 0
    "Null-ish (no edges above 1e-6; worst $(round(maxR; sigdigits=2)))"
elseif n_above_1e3 > 0
    "Confirmed-pending-spatial-check (≥1 edge above 1e-3; see spatial map)"
else
    "Partially confirmed (edges above 1e-6 but none above 1e-3)"
end

println()
println("=" ^ 78)
@printf("VERDICT (preliminary, numbers-only): %s\n", verdict)
@printf("Total wall: %.2f s (walk %.2f s + tri+eval %.2f s)\n",
        time() - t_walk0, walk_wall, eval_wall)
println("=" ^ 78)

# Persist a one-line summary for the REPORT.md writer.
open(joinpath(OUTPUT_DIR, "summary.txt"), "w") do io
    @printf(io, "verdict_numbers_only: %s\n", verdict)
    @printf(io, "n_visited:            %d\n", N)
    @printf(io, "n_sheet0:             %d\n", n_sheet0)
    @printf(io, "n_delaunay_edges:     %d\n", length(delaunay_edges))
    @printf(io, "n_tree_edges:         %d\n", length(tree_edges_local))
    @printf(io, "n_nontree:            %d\n", length(nontree_edges))
    @printf(io, "n_good_evals:         %d\n", length(results))
    @printf(io, "n_pade_blowups:       %d\n", n_blowups)
    @printf(io, "median_dP_rel:        %.6e\n", med)
    @printf(io, "p90_dP_rel:           %.6e\n", p90)
    @printf(io, "p99_dP_rel:           %.6e\n", p99)
    @printf(io, "max_dP_rel:           %.6e\n", maxR)
    @printf(io, "n_above_1e-10:        %d\n", n_above_tol)
    @printf(io, "n_above_1e-6:         %d\n", n_above_1e6)
    @printf(io, "n_above_1e-3:         %d\n", n_above_1e3)
    @printf(io, "corr_treedist:        %.4f\n", corr_td)
    @printf(io, "corr_extrapmax:       %.4f\n", corr_em)
    @printf(io, "wall_walk_s:          %.2f\n", walk_wall)
    @printf(io, "wall_eval_s:          %.2f\n", eval_wall)
    if !isempty(bad)
        @printf(io, "bad_re_median:        %.3f\n", median([real(r.M) for r in bad]))
        @printf(io, "bad_im_median:        %.3f\n", median([imag(r.M) for r in bad]))
        @printf(io, "bad_frac_re_gt_3:     %.4f\n",
                count(>(3.0), [real(r.M) for r in bad]) / length(bad))
        @printf(io, "bad_frac_abs_im_gt_pi: %.4f\n",
                count(>(Float64(π)), [abs(imag(r.M)) for r in bad]) / length(bad))
    end
end
@printf("Wrote %s\n", joinpath(OUTPUT_DIR, "summary.txt"))
