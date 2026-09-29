"""
Midpoint evaluation and aggregation for the scalar and vector diagnostics.

ADR-0016a makes the analysed population explicit. Each candidate non-tree
edge is evaluated once or counted once as a pole exception / nonfinite
result. Unexpected exceptions propagate. No surviving edges after failed
evaluations is an error with a suggestion, not an empty clean certificate.
Quantiles retain their original meaning over successfully evaluated edges.
"""
_diagnostic_finite(x::Number) = isfinite(x)
_diagnostic_finite(x::AbstractVector) = all(isfinite, x)

function _diagnostic_edges(sol, s_idx, nontree_edges, evaluate, discrepancy,
                           tol_w, tol_b)
    depth = _build_depths(sol.visited_parent)
    edges = EdgeReport[]
    sizehint!(edges, length(nontree_edges))
    n_threw = n_nonfinite = 0
    for (la, lb) in nontree_edges
        ga, gb = s_idx[la], s_idx[lb]
        zA, zB = sol.visited_z[ga], sol.visited_z[gb]
        M = (zA + zB) / 2
        tA = (M - zA) / sol.visited_h[ga]
        tB = (M - zB) / sol.visited_h[gb]
        values = try
            (evaluate(sol, ga, tA), evaluate(sol, gb, tB))
        catch err
            err isa DomainError || rethrow()
            n_threw += 1
            continue
        end
        uA, uB = values
        if !(_diagnostic_finite(uA) && _diagnostic_finite(uB))
            n_nonfinite += 1
            continue
        end
        ΔP_abs, ΔP_rel = discrepancy(uA, uB)
        if !(isfinite(ΔP_abs) && isfinite(ΔP_rel))
            n_nonfinite += 1
            continue
        end
        td = _tree_path_distance(sol.visited_parent, depth, ga, gb)
        em = Float64(max(abs(tA), abs(tB)))
        cat = _categorise(ΔP_rel, em, tol_w, tol_b)
        push!(edges, EdgeReport(ga, gb, ΔP_abs, ΔP_rel, td, em,
                                ComplexF64(M), cat))
    end
    return edges, n_threw, n_nonfinite
end

function _diagnostic_report(edges, sampling, n_threw, n_nonfinite,
                            sheet, tol_w, tol_b, n_worst)
    n_edges = length(edges)
    n_dropped = n_threw + n_nonfinite
    if n_edges == 0 && n_dropped > 0
        throw(ArgumentError(
            "quality_diagnose: all $(sampling.n_edges_candidate) candidate " *
            "edges failed (threw=$n_threw, nonfinite=$n_nonfinite). " *
            "Suggestion: inspect the stored Padés and refine the node cloud " *
            "or shorten steps around poles before trusting a certificate."))
    end
    counts = (sampling.n_nodes, sampling.n_nodes_retained,
              sampling.n_nodes - sampling.n_nodes_retained,
              sampling.n_edges_candidate, n_dropped, n_threw, n_nonfinite,
              sampling.n_tree_edges_off_sheet)
    if n_edges == 0
        return DiagnosticReport(0, 0, 0, 0, 0, 0, NaN, NaN, NaN, NaN,
                                EdgeReport[], complex(NaN, NaN),
                                sheet, tol_w, tol_b, counts...)
    end
    rels = Float64[e.ΔP_rel for e in edges]
    n_well = count(e -> e.category == :well_closed, edges)
    n_noisy = count(e -> e.category == :noisy, edges)
    n_extr = count(e -> e.category == :extrap_driven, edges)
    n_depth = count(e -> e.category == :depth_driven, edges)
    order_desc = sortperm(rels; rev=true)
    worst = edges[order_desc[1:min(Int(n_worst), n_edges)]]
    bad_mids = ComplexF64[e.midpoint for e in edges if e.ΔP_rel > tol_b]
    centroid = isempty(bad_mids) ? complex(NaN, NaN) : ComplexF64(mean(bad_mids))
    return DiagnosticReport(n_edges, n_well, n_noisy, n_extr, n_depth, 0,
                            median(rels), quantile(rels, 0.90),
                            quantile(rels, 0.99), maximum(rels), worst,
                            centroid, sheet, tol_w, tol_b, counts...)
end
