"""
Vector diagnostics adapter of PadeTaylorDiagnosticsExt (ADR-0025 Amend. 7).

Shared-Q patches evaluate each component Pᵢ(t)/Q(t), and midpoint
disagreement uses the vector 2-norm. Every node participates because
VectorPathNetworkSolution has no sheet metadata. ADR-0016a shares geometry,
edge accounting and aggregation with the scalar method: a failed edge is
counted once, and all-candidate failure throws with a suggestion. Valid
walks without candidate edges retain the existing empty-report convention.
"""
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


_evaluate_vector_node(sol, k, t) = _eval_shared_q(
    sol.visited_numerators[k], sol.visited_denominator[k], t)
function _vector_disagreement(yA, yB)
    ΔP_abs = Float64(sqrt(sum(abs2, yA .- yB)))
    denom = Float64(sqrt(sum(abs2, yA))) +
            Float64(sqrt(sum(abs2, yB))) + _EPS_FLOOR
    return ΔP_abs, ΔP_abs / denom
end

"""
    quality_diagnose(sol::VectorPathNetworkSolution; tol_well=1e-10,
                     tol_bad=1e-6, n_worst=10) -> DiagnosticReport

Diagnose every node using shared-Q midpoint evaluations and vector 2-norm
disagreement. Coverage and failed evaluations appear in the report. No
sheet keyword is accepted; fewer than three nodes have no candidate edges.
"""
function quality_diagnose(sol::VectorPathNetworkSolution;
                          tol_well::Real=1e-10, tol_bad::Real=1e-6,
                          n_worst::Integer=10)
    N = length(sol.visited_z)
    s_idx = collect(1:N)
    nontree_edges, n_off_sheet = _diagnostic_geometry(sol, s_idx; allow_degenerate=true)
    tol_w, tol_b = Float64(tol_well), Float64(tol_bad)
    edges, n_threw, n_nonfinite = _diagnostic_edges(
        sol, s_idx, nontree_edges, _evaluate_vector_node, _vector_disagreement,
        tol_w, tol_b)
    sampling = (; n_nodes=N, n_nodes_retained=N,
                 n_edges_candidate=length(nontree_edges),
                 n_tree_edges_off_sheet=n_off_sheet)
    return _diagnostic_report(edges, sampling, n_threw, n_nonfinite,
                              0, tol_w, tol_b, n_worst)
end
