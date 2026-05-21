# ext/diagnostics_vector_adapter.jl
#
# The VC-7 **vector adapter** of `PadeTaylorDiagnosticsExt` — the method
# `quality_diagnose(::VectorPathNetworkSolution)`.  `include`d by
# `ext/PadeTaylorDiagnosticsExt.jl` at the end of the extension module,
# purely for the CLAUDE.md Rule 6 ≤200-effective-LOC cap: the scalar
# `quality_diagnose(::PathNetworkSolution)` plus this vector method in
# one file exceeded the cap, and the natural seam is scalar-vs-vector.
#
# This is NOT a separate module — it is `include`d *into*
# `PadeTaylorDiagnosticsExt`, so it shares that module's `using`s
# (`triangulate`, `each_edge`, `median`/`quantile`/`mean`, the
# `DiagnosticReport`/`EdgeReport` structs) and its private helpers
# (`_build_depths`, `_tree_path_distance`, `_categorise`, `_EPS_FLOOR`).
# A `module` here would break that sharing; the file is exposition that
# continues the extension module, the same `include`-a-fragment pattern
# `src/PadeTaylor.jl` uses for its own Rule-6 splits.
#
# ## What this adapter is — ADR-0025 Amendment 3 §"D-VC7 scope" / Amend. 7
#
# `quality_diagnose` (ADR-0016) — the loop-closure quality certificate —
# was typed for the scalar `PathNetworkSolution`.  The v0.2 vector
# path-network walk (the substrate of the headline `P_I⁽²⁾` tritronquée
# figure) produces a `VectorPathNetworkSolution`.  This file adds the
# additive vector method; the scalar method is byte-untouched (ADR-0001),
# so every v0.1 + v0.2 scalar test stays bit-identical.
#
# The vector method differs from the scalar one in exactly three ways
# (the A4 baseline probe `external/probes/v8b-baseline/REPORT.md` §2b
# catalogued the data-model gap and is the verified reference impl):
#
#   - **shared-`Q` evaluation.**  A scalar node stores one
#     `PadeApproximant`; a vector node stores `d` numerator polynomials
#     `visited_numerators[k]` over ONE denominator `visited_denominator[k]`
#     (ADR-0019).  The adapter evaluates component `i` as the Horner
#     ratio `Pᵢ(t)/Q(t)` — `_eval_shared_q` below — yielding a `d`-vector.
#   - **vector-norm `ΔP_rel`.**  `ΔP_rel` is the Euclidean 2-norm over
#     the `d` companion components:
#     `‖y_A − y_B‖ / (‖y_A‖ + ‖y_B‖ + ε)`.
#   - **no sheet mask.**  The `P_I⁽²⁾` companion system is meromorphic —
#     single-sheeted — so `VectorPathNetworkSolution` has no
#     `visited_sheet` field; every node participates, the report's
#     `sheet` is fixed at `0`, and there is no `sheet` kwarg.
#
# The Delaunay / tree-edge / non-tree / LCA tree-distance / categorise /
# aggregate machinery ports verbatim from the scalar method in
# `PadeTaylorDiagnosticsExt.jl`.
#
# See the `PadeTaylorDiagnosticsExt` module docstring (§"The vector
# adapter") and ADR-0025 Amendment 7 for the full design and the
# measured loop-closure numbers.

# Horner evaluation of the polynomial `Σ c[k+1] tᵏ` (coefficients
# low-to-high) — the convention the shared-`Q` numerators and
# denominator are stored in.  A private copy of
# `VectorPathNetworkStage2._eval_poly` / `VectorProblems._eval_poly`:
# copied, not imported, to keep the extension free of a reach into a
# sibling module's private symbol (the same module-boundary discipline
# the Stage-2 module follows for its own `_eval_poly` copy).
function _eval_poly(c::AbstractVector, t)
    s = zero(eltype(c)) * t
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

# Evaluate the `d`-vector shared-`Q` approximant of a node at the
# rescaled coordinate `t`: every component `Pᵢ(t)` over the *same*
# shared denominator `Q(t)` — the v0.2 keystone (ADR-0019) applied at
# evaluation time, exactly `VectorPathNetworkStage2._stage2_fill`'s
# inner loop.  Returns the `d`-vector state `y(t)`.
function _eval_shared_q(numerators::AbstractVector,
                        denominator::AbstractVector, t)
    q_t = _eval_poly(denominator, t)
    return [_eval_poly(num, t) / q_t for num in numerators]
end

"""
    quality_diagnose(sol::VectorPathNetworkSolution; tol_well=1e-10,
                     tol_bad=1e-6, n_worst=10) -> DiagnosticReport

Loop-closure quality certificate for a **vector** path-network walk —
the VC-7 criterion of ADR-0025 (Amendment 3 §"D-VC7 scope" / Amendment
7, bead `padetaylor-0ln.37.14`).  This is the additive vector analogue
of the scalar `quality_diagnose(::PathNetworkSolution)`; see the
`PadeTaylorDiagnosticsExt` module docstring's "vector adapter" section
and this file's header for the data-model differences and the
verbatim-ported machinery.

The procedure mirrors the scalar method:

  1. **No sheet mask.**  The `P_I⁽²⁾` companion system is meromorphic —
     single-sheeted — so `VectorPathNetworkSolution` carries no
     `visited_sheet` field; every visited node participates and the
     report's `sheet` is fixed at `0`.
  2. **Delaunay-triangulate** the visited-node cloud; drop ghost edges
     (negative vertex indices).
  3. **Tree edges** are `{(visited_parent[k], k) : k ≥ 2}`.
  4. **Non-tree edges = Delaunay − tree** — the loop-closure population.
  5. **Per-edge `ΔP_rel`.**  For non-tree edge `(A, B)`,
     `M = (z_A + z_B)/2`, `t_X = (M − z_X)/visited_h[X]`,
     `y_X = Pᵢ(t_X)/Q(t_X)` over node `X`'s shared-`Q` store, and
     `ΔP_rel = ‖y_A − y_B‖₂ / (‖y_A‖₂ + ‖y_B‖₂ + ε)`.
  6. **Categorise** with the `Diagnostics` module thresholds and
     **aggregate** the quantiles, top-`n_worst` edges, and bad-midpoint
     centroid.

Unlike the scalar method this adapter takes no `sheet` kwarg — there is
only one sheet.  A walk with fewer than 3 visited nodes, or with no
non-tree Delaunay edge, returns a well-formed empty `DiagnosticReport`
(`n_edges = 0`, `NaN` quantiles) — it does not throw, exactly as the
scalar method's empty-edge branch behaves.
"""
function quality_diagnose(sol::VectorPathNetworkSolution;
                          tol_well::Real = 1e-10,
                          tol_bad::Real  = 1e-6,
                          n_worst::Integer = 10)
    N = length(sol.visited_z)
    tol_w = Float64(tol_well)
    tol_b = Float64(tol_bad)

    # A walk with fewer than 3 visited nodes has no Delaunay
    # triangulation at all (a triangulation needs ≥ 3 non-collinear
    # points) — hence no non-tree edges and no loop-closure population.
    # This is a degenerate-but-valid input (a one- or two-target walk),
    # not an error: return a well-formed empty report, exactly as the
    # n_edges == 0 branch below does.  Without this guard `triangulate`
    # throws an `InsufficientPointsError` on a one-node walk.
    if N < 3
        return DiagnosticReport(0, 0, 0, 0, 0, 0, NaN, NaN, NaN, NaN,
                                EdgeReport[], complex(NaN, NaN),
                                0, tol_w, tol_b)
    end

    # No sheet mask: P_I⁽²⁾ in companion form is single-sheeted, so
    # every visited node participates and global index == local index.
    pts2D = [(Float64(real(z)), Float64(imag(z))) for z in sol.visited_z]
    tri = triangulate(pts2D)
    norm_edge(a, b) = a < b ? (a, b) : (b, a)
    delaunay_edges = Set{Tuple{Int,Int}}()
    for e in each_edge(tri)
        a, b = e
        (a < 1 || b < 1) && continue
        push!(delaunay_edges, norm_edge(a, b))
    end

    tree_edges = Set{Tuple{Int,Int}}()
    for k in 2:N
        p = sol.visited_parent[k]
        p == 0 && continue
        push!(tree_edges, norm_edge(k, p))
    end

    nontree_edges = setdiff(delaunay_edges, tree_edges)

    depth = _build_depths(sol.visited_parent)

    edges = EdgeReport[]
    sizehint!(edges, length(nontree_edges))
    for (a, b) in nontree_edges
        zA, zB = sol.visited_z[a], sol.visited_z[b]
        hA, hB = sol.visited_h[a], sol.visited_h[b]
        M  = (zA + zB) / 2
        tA = (M - zA) / hA
        tB = (M - zB) / hB
        yA = try
            _eval_shared_q(sol.visited_numerators[a],
                           sol.visited_denominator[a], tA)
        catch
            continue
        end
        yB = try
            _eval_shared_q(sol.visited_numerators[b],
                           sol.visited_denominator[b], tB)
        catch
            continue
        end
        (all(c -> isfinite(real(c)) && isfinite(imag(c)), yA) &&
         all(c -> isfinite(real(c)) && isfinite(imag(c)), yB)) || continue
        ΔP_abs = Float64(sqrt(sum(abs2, yA .- yB)))
        denom  = Float64(sqrt(sum(abs2, yA))) +
                 Float64(sqrt(sum(abs2, yB))) + _EPS_FLOOR
        ΔP_rel = ΔP_abs / denom
        td     = _tree_path_distance(sol.visited_parent, depth, a, b)
        em     = Float64(max(abs(tA), abs(tB)))
        cat    = _categorise(ΔP_rel, em, tol_w, tol_b)
        push!(edges, EdgeReport(a, b, ΔP_abs, ΔP_rel, td, em,
                                ComplexF64(M), cat))
    end

    n_edges = length(edges)
    rels = Float64[e.ΔP_rel for e in edges]
    if n_edges == 0
        return DiagnosticReport(0, 0, 0, 0, 0, 0, NaN, NaN, NaN, NaN,
                                EdgeReport[], complex(NaN, NaN),
                                0, tol_w, tol_b)
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
                            centroid, 0, tol_w, tol_b)
end
