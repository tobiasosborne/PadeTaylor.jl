"""
Delaunay geometry and tree-distance helpers for diagnostics (ADR-0016a).

Triangulation remains restricted to retained sheet-0 nodes, exactly as in
ADR-0016. Ghost edges are excluded; parent links outside the retained set
are counted separately from candidate non-tree edges. This distinction
keeps the midpoint-evaluation denominator truthful; allow_degenerate keeps
vector empty walks valid. Cross-sheet geometry remains padetaylor-8py work.
"""
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


function _diagnostic_geometry(sol, s_idx; allow_degenerate::Bool=false)
    N = length(sol.visited_z)
    global_to_local = zeros(Int, N)
    for (loc, glob) in enumerate(s_idx)
        global_to_local[glob] = loc
    end
    norm_edge(a, b) = a < b ? (a, b) : (b, a)
    delaunay_edges = Set{Tuple{Int,Int}}()
    if !allow_degenerate || length(s_idx) >= 3
        pts2D = [(Float64(real(sol.visited_z[g])), Float64(imag(sol.visited_z[g])))
                 for g in s_idx]
        tri = triangulate(pts2D)
        for (a, b) in each_edge(tri)
            (a < 1 || b < 1) && continue
            push!(delaunay_edges, norm_edge(a, b))
        end
    end
    tree_edges = Set{Tuple{Int,Int}}()
    n_off_sheet = 0
    for k in 2:N
        p = sol.visited_parent[k]
        p == 0 && continue
        lk, lp = global_to_local[k], global_to_local[p]
        if lk == 0 || lp == 0
            n_off_sheet += 1
            continue
        end
        push!(tree_edges, norm_edge(lk, lp))
    end
    return setdiff(delaunay_edges, tree_edges), n_off_sheet
end
