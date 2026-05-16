"""
    PadeTaylorDiagnosticsExt

Package-extension providing the Delaunay-backed implementation of
`Diagnostics.quality_diagnose` per ADR-0016.  Activated when both
`PadeTaylor` and `DelaunayTriangulation` are loaded:

```julia
using PadeTaylor, DelaunayTriangulation

sol    = path_network_solve(prob, grid; diagnose = true)
report = sol.diagnostics            # ::DiagnosticReport
println(report)                     # categorised loop-closure summary
```

Equivalent post-hoc:

```julia
sol    = path_network_solve(prob, grid)        # no diagnostics attached
report = quality_diagnose(sol)                 # compute on demand
```

## Algorithm (verbatim from `external/probes/loop-closure-fig1/probe.jl`)

The probe at `external/probes/loop-closure-fig1/probe.jl:156-329` is
the algorithmic source; we port lines 156-329 here with two narrow
adjustments: (a) the sheet-0 predicate is generalised to honour
`visited_sheet[k] == [0]` when branched (the probe ran with no
branches, so `visited_sheet[k]` was always empty there); (b) the
top-N extraction is parameterised via `n_worst` rather than the
probe's hard-coded 10.

  1. **Sheet-0 filter.**  When `visited_sheet[k]` is empty for every
     visited node (the no-branches case, FFW Fig 1), filter to the
     ζ-strip `-2π < imag(z) ≤ 2π` (FFW md:103).  When non-empty,
     filter to `visited_sheet[k] == [sheet]`.  Sheet defaults to `0`
     (the principal sheet).
  2. **Delaunay triangulate** the sheet-0 nodes via
     `DelaunayTriangulation.triangulate` on 2-tuple coordinates.
     Filter ghost edges (negative vertex indices for the unbounded
     face's representatives — probe.jl:196-201).
  3. **Tree edges.**  Edges `(visited_parent[k], k)` for `k ≥ 2`,
     mapped through the sheet-0 sub-indexing.  An edge crossing
     off-sheet (one endpoint not in sheet 0) is dropped.
  4. **Non-tree edges = Delaunay edges − tree edges.**  These are the
     loop-closure population: every one corresponds to a cycle in the
     visited-node graph that the Stage-1 tree omits.
  5. **Per-edge ΔP_rel.**  For each non-tree edge `(A, B)`:
     `M = (z_A + z_B) / 2`, `t_X = (M - z_X) / visited_h[X]`,
     `u_X = _evaluate_pade(visited_pade[X], t_X)`, then
     `ΔP_rel = |u_A − u_B| / (|u_A| + |u_B| + ε)`.  Padé denominator
     blow-ups (probe.jl:310-321) are caught and the edge skipped.
  6. **Categorise** per `Diagnostics`'s module-docstring thresholds:
     `:well_closed`, `:noisy`, `:extrap_driven`, `:depth_driven`,
     `:branch_cut` (reserved, v2).
  7. **Aggregate** quantiles, top-N worst edges, and the centroid of
     "bad" midpoints (`ΔP_rel > tol_bad`).

## Tree-distance via LCA on parent chains

`visited_parent` is a tree on the global indices; `visited_parent[root]
= 0`.  We compute per-node depth once, then for each non-tree edge
walk both endpoints up to equal depth and continue until they meet
(probe.jl:236-274).  This is `O(depth)` per edge — fine for the
edge counts we see in practice (~2000 edges, depth a few hundred).

## Why this is an extension, not core

`DelaunayTriangulation.jl` carries `ExactPredicates`, `AdaptivePredicates`,
and `EnumX` as transitive deps.  Loading them eagerly would more than
double `using PadeTaylor`'s precompile time for users who only want a
Stage-1 walk.  ADR-0003 documents the precedent (Arblib, CommonSolve,
Makie); ADR-0016 ties this extension to it.

## References

  - ADR-0016 — `docs/adr/0016-diagnostics-extension.md` (this design).
  - Probe — `external/probes/loop-closure-fig1/probe.jl` (algorithmic source).
  - Probe verdict — `external/probes/loop-closure-fig1/REPORT.md:79-98`.
  - FFW 2017 §2.1.2 — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-103`.
"""
module PadeTaylorDiagnosticsExt

using PadeTaylor: PathNetworkSolution
using PadeTaylor.Diagnostics: DiagnosticReport, EdgeReport, quality_diagnose
using PadeTaylor.PathNetwork: _evaluate_pade
using DelaunayTriangulation: triangulate, each_edge
using Statistics: median, quantile, mean

import PadeTaylor.Diagnostics: quality_diagnose

const _EPS_FLOOR = 1e-300

# Sheet-0 mask: honour visited_sheet[k] == [sheet] when branched;
# fall back to the FFW ζ-strip predicate when not.  The probe at
# probe.jl:156-158 used only the strip branch (Fig 1 ran without
# branch points); this generalisation lets the same diagnostic
# work on branched walks without changing the caller-facing API.
function _sheet_mask(visited_z::AbstractVector,
                     visited_sheet::AbstractVector,
                     sheet::Int)
    if isempty(visited_sheet) || all(isempty, visited_sheet)
        return [(-2π < imag(z) ≤ 2π) for z in visited_z]
    end
    target = Int[sheet]
    return [s == target for s in visited_sheet]
end

# Depths via iterative parent-chain relaxation (probe.jl:236-257).
# `visited_parent[root] = 0`; we seed roots at depth 0 and fill the
# rest by `depth[k] = depth[parent[k]] + 1` until no more progress.
function _build_depths(parent::AbstractVector{Int})
    n = length(parent)
    depth = fill(-1, n)
    for k in 1:n
        parent[k] == 0 && (depth[k] = 0)
    end
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

# Tree-distance via LCA on parent chains (probe.jl:261-274).
function _tree_path_distance(parent::AbstractVector{Int},
                             depth::AbstractVector{Int},
                             a::Int, b::Int)
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

# Bucket an edge into one of the four v1 categories (`:branch_cut`
# is reserved for multi-sheet, never emitted here).
function _categorise(ΔP_rel::Float64, extrap_max::Float64,
                     tol_well::Float64, tol_bad::Float64)
    ΔP_rel ≤ tol_well   && return :well_closed
    ΔP_rel ≤ tol_bad    && return :noisy
    extrap_max > 1.0    && return :extrap_driven
    return :depth_driven
end

function quality_diagnose(sol::PathNetworkSolution;
                          sheet::Int    = 0,
                          tol_well::Real = 1e-10,
                          tol_bad::Real  = 1e-6,
                          n_worst::Integer = 10)
    sheet == 0 || throw(ArgumentError(
        "quality_diagnose: v1 supports sheet 0 only (got sheet=$sheet). " *
        "Multi-sheet diagnostics deferred to bead padetaylor-8py; track " *
        "there for the cut-aware Delaunay step that lifts this restriction."))

    N = length(sol.visited_z)
    mask  = _sheet_mask(sol.visited_z, sol.visited_sheet, sheet)
    s_idx = findall(mask)
    n_s   = length(s_idx)

    global_to_local = zeros(Int, N)
    for (loc, glob) in enumerate(s_idx)
        global_to_local[glob] = loc
    end

    # Delaunay on sheet-`sheet` ζ-coordinates.  Each point is a 2-tuple
    # (Re, Im); `each_edge` returns undirected (a, b) pairs.  Ghost
    # vertices on the convex hull come back as negative indices —
    # filter per probe.jl:196-201.
    pts2D = [(Float64(real(sol.visited_z[g])), Float64(imag(sol.visited_z[g])))
             for g in s_idx]
    tri = triangulate(pts2D)
    norm_edge(a, b) = a < b ? (a, b) : (b, a)
    delaunay_edges = Set{Tuple{Int,Int}}()
    for e in each_edge(tri)
        a, b = e
        (a < 1 || b < 1) && continue
        push!(delaunay_edges, norm_edge(a, b))
    end

    # Tree edges restricted to sheet `sheet` (both endpoints in-sheet).
    tree_edges = Set{Tuple{Int,Int}}()
    for k in 2:N
        p = sol.visited_parent[k]
        p == 0 && continue
        lk = global_to_local[k]
        lp = global_to_local[p]
        (lk == 0 || lp == 0) && continue
        push!(tree_edges, norm_edge(lk, lp))
    end

    nontree_edges = setdiff(delaunay_edges, tree_edges)

    depth = _build_depths(sol.visited_parent)
    tol_w = Float64(tol_well)
    tol_b = Float64(tol_bad)

    edges = EdgeReport[]
    sizehint!(edges, length(nontree_edges))
    for (la, lb) in nontree_edges
        ga, gb = s_idx[la], s_idx[lb]
        zA, zB = sol.visited_z[ga], sol.visited_z[gb]
        hA, hB = sol.visited_h[ga], sol.visited_h[gb]
        M  = (zA + zB) / 2
        tA = (M - zA) / hA
        tB = (M - zB) / hB
        uA = try
            _evaluate_pade(sol.visited_pade[ga], tA)
        catch
            continue
        end
        uB = try
            _evaluate_pade(sol.visited_pade[gb], tB)
        catch
            continue
        end
        (isfinite(real(uA)) && isfinite(imag(uA)) &&
         isfinite(real(uB)) && isfinite(imag(uB))) || continue
        ΔP_abs = Float64(abs(uA - uB))
        denom  = Float64(abs(uA)) + Float64(abs(uB)) + _EPS_FLOOR
        ΔP_rel = ΔP_abs / denom
        td     = _tree_path_distance(sol.visited_parent, depth, ga, gb)
        em     = Float64(max(abs(tA), abs(tB)))
        cat    = _categorise(ΔP_rel, em, tol_w, tol_b)
        push!(edges, EdgeReport(ga, gb, ΔP_abs, ΔP_rel, td, em,
                                ComplexF64(M), cat))
    end

    n_edges = length(edges)
    rels = Float64[e.ΔP_rel for e in edges]
    if n_edges == 0
        return DiagnosticReport(0, 0, 0, 0, 0, 0, NaN, NaN, NaN, NaN,
                                EdgeReport[], complex(NaN, NaN),
                                sheet, tol_w, tol_b)
    end

    n_well   = count(==(:well_closed),   (e.category for e in edges))
    n_noisy  = count(==(:noisy),         (e.category for e in edges))
    n_extr   = count(==(:extrap_driven), (e.category for e in edges))
    n_depth  = count(==(:depth_driven),  (e.category for e in edges))
    n_branch = 0   # reserved for v2; see ADR-0016 §"Consequences"

    med  = median(rels)
    p90  = quantile(rels, 0.90)
    p99  = quantile(rels, 0.99)
    maxR = maximum(rels)

    order_desc = sortperm(rels; rev = true)
    n_worst_eff = min(Int(n_worst), n_edges)
    worst = EdgeReport[edges[order_desc[i]] for i in 1:n_worst_eff]

    bad_mids = ComplexF64[e.midpoint for e in edges if e.ΔP_rel > tol_b]
    centroid = isempty(bad_mids) ? complex(NaN, NaN) :
                ComplexF64(mean(bad_mids))

    return DiagnosticReport(n_edges, n_well, n_noisy, n_extr, n_depth,
                            n_branch, med, p90, p99, maxR, worst,
                            centroid, sheet, tol_w, tol_b)
end

end # module PadeTaylorDiagnosticsExt
