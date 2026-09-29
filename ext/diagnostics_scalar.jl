"""
Scalar sheet selection and diagnostics entry point (ADR-0016a).

Empty sheet metadata means a branchless walk, not a coordinate frame.
Every such node participates. The ζ-strip from FFW 2017 md:103 cannot
be inferred from that metadata, so it is never applied implicitly.
Branched walks retain the existing [0] predicate and sheet-0 restriction.
Any node loss needs the caller's explicit retained-fraction budget.
"""
function _sheet_mask(visited_z::AbstractVector,
                     visited_sheet::AbstractVector, sheet::Int)
    if isempty(visited_sheet) || all(isempty, visited_sheet)
        return trues(length(visited_z))
    end
    length(visited_sheet) == length(visited_z) || throw(ArgumentError(
        "quality_diagnose: sheet metadata length does not match the nodes. " *
        "Suggestion: provide one visited_sheet entry per visited node."))
    target = Int[sheet]
    return [s == target for s in visited_sheet]
end

function _check_retained_nodes(N, n_s, min_retained_fraction)
    isfinite(min_retained_fraction) && 0 <= min_retained_fraction <= 1 ||
        throw(ArgumentError(
            "quality_diagnose: min_retained_fraction must be finite and in [0,1]. " *
            "Suggestion: use 1.0 for full coverage or explicitly budget a " *
            "known sheet-0 subset."))
    if N > 0 && n_s / N < min_retained_fraction
        throw(ArgumentError(
            "quality_diagnose: retained $n_s / $N nodes " *
            "($(100*n_s/N)%), below min_retained_fraction=$min_retained_fraction. " *
            "Suggestion: inspect visited_sheet; if this sheet-0 subset is " *
            "intentional, pass min_retained_fraction=$(n_s/N) explicitly."))
    end
end

_evaluate_scalar_node(sol, k, t) = _evaluate_pade(sol.visited_pade[k], t)
function _scalar_disagreement(uA, uB)
    ΔP_abs = Float64(abs(uA - uB))
    denom = Float64(abs(uA)) + Float64(abs(uB)) + _EPS_FLOOR
    return ΔP_abs, ΔP_abs / denom
end

function quality_diagnose(sol::PathNetworkSolution;
                          sheet::Int=0, tol_well::Real=1e-10,
                          tol_bad::Real=1e-6, n_worst::Integer=10,
                          min_retained_fraction::Real=1.0)
    sheet == 0 || throw(ArgumentError(
        "quality_diagnose: v1 supports sheet 0 only (got sheet=$sheet). " *
        "Multi-sheet diagnostics deferred to bead padetaylor-8py; track " *
        "there for the cut-aware Delaunay step that lifts this restriction."))
    N = length(sol.visited_z)
    mask = _sheet_mask(sol.visited_z, sol.visited_sheet, sheet)
    s_idx = findall(mask)
    n_s = length(s_idx)
    _check_retained_nodes(N, n_s, min_retained_fraction)
    nontree_edges, n_off_sheet = _diagnostic_geometry(sol, s_idx)
    tol_w, tol_b = Float64(tol_well), Float64(tol_bad)
    edges, n_threw, n_nonfinite = _diagnostic_edges(
        sol, s_idx, nontree_edges, _evaluate_scalar_node, _scalar_disagreement,
        tol_w, tol_b)
    sampling = (; n_nodes=N, n_nodes_retained=n_s,
                 n_edges_candidate=length(nontree_edges),
                 n_tree_edges_off_sheet=n_off_sheet)
    return _diagnostic_report(edges, sampling, n_threw, n_nonfinite,
                              sheet, tol_w, tol_b, n_worst)
end
