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

## Algorithm (probe, corrected by ADR-0016a)

The probe at `external/probes/loop-closure-fig1/probe.jl:156-329` is
the original source. ADR-0016a corrects its frame-specific mask and
accounts for failed evaluations; `n_worst` parameterises top-N extraction.

  1. **Sheet-0 filter.** Branchless walks retain every node: empty sheet
     metadata does not identify a ζ-frame. Branched walks retain
     `visited_sheet[k] == [0]`. The default `min_retained_fraction=1.0`
     rejects any loss; a caller can explicitly budget a known subset.
  2. **Delaunay triangulate** the sheet-0 nodes via
     `DelaunayTriangulation.triangulate` on 2-tuple coordinates.
     Filter ghost edges (negative vertex indices for the unbounded
     face's representatives — probe.jl:196-201).
  3. **Tree edges.**  Edges `(visited_parent[k], k)` for `k ≥ 2`,
     mapped through the sheet-0 sub-indexing.  A parent link crossing
     off-sheet is counted as `n_tree_edges_off_sheet`, not a failed evaluation.
  4. **Non-tree edges = Delaunay edges − tree edges.**  These are the
     loop-closure population: every one corresponds to a cycle in the
     visited-node graph that the Stage-1 tree omits.
  5. **Per-edge ΔP_rel.**  For each non-tree edge `(A, B)`:
     `M = (z_A + z_B) / 2`, `t_X = (M - z_X) / visited_h[X]`,
     `u_X = _evaluate_pade(visited_pade[X], t_X)`, then
     `ΔP_rel = |u_A − u_B| / (|u_A| + |u_B| + ε)`.  Padé denominator
     pole exceptions and nonfinite values each count one dropped edge.
     Unexpected exceptions propagate; losing every candidate throws with
     a suggestion. The report displays both coverage and failure counts.
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

## The vector adapter — `quality_diagnose(::VectorPathNetworkSolution)`

ADR-0025 Amendment 3 (§"D-VC7 scope", bead `padetaylor-0ln.37.14`)
brings the same loop-closure certificate to the v0.2 *vector*
path-network walk — the substrate of the headline `P_I⁽²⁾` tritronquée
figure.  The scalar `quality_diagnose` above is typed for
`PathNetworkSolution`; the vector walk produces a
`VectorPathNetworkSolution`, whose data model differs in three ways
the A4 baseline probe (`external/probes/v8b-baseline/REPORT.md` §2b)
catalogued and whose §2b is the verified reference implementation this
adapter lifts:

  - **per-node approximant.**  A scalar node stores a single
    `PadeApproximant` and the `_evaluate_pade` evaluator returns a
    scalar.  A vector node stores a *shared-`Q`* approximant — `d`
    numerator polynomials `visited_numerators[k]` over **one**
    denominator `visited_denominator[k]` (ADR-0019).  The adapter
    evaluates component `i` as the Horner ratio `Pᵢ(t)/Q(t)`, the
    identical `_eval_poly` pattern `VectorPathNetworkStage2._stage2_fill`
    uses; the result is a `d`-vector `y`, not a scalar.
  - **`ΔP_rel` is a vector norm.**  The scalar `ΔP_rel` is generalised
    to the Euclidean 2-norm over the `d` companion components:
    `ΔP_rel = ‖y_A(M) − y_B(M)‖ / (‖y_A(M)‖ + ‖y_B(M)‖ + ε)`.
  - **no sheet mask.**  `PathNetworkSolution` carries a `visited_sheet`
    field and the scalar method filters to sheet 0.  The `P_I⁽²⁾`
    companion system is meromorphic — single-sheeted — so
    `VectorPathNetworkSolution` has **no** `visited_sheet` field at all
    (`VectorPathNetwork.jl` docstring).  The vector adapter therefore
    drops `_sheet_mask` entirely: every visited node participates,
    global index == local index, and the `sheet` field of the returned
    `DiagnosticReport` is fixed at `0` ("the single sheet").

Everything else — the Delaunay triangulation, the tree-edge
subtraction, the non-tree edge set, the midpoint `t = (M − z_X)/h_X`
rescaling, the `_build_depths` / `_tree_path_distance` LCA tree
distance, the `_categorise` thresholds, and the `DiagnosticReport`
aggregation — shares helpers with the scalar method. The
`EdgeReport` / `DiagnosticReport` structs are reused: they
store `ΔP_rel` / `ΔP_abs` as `Float64`, agnostic to a scalar-vs-vector
origin.  Both methods report sample counts and midpoint failures under
ADR-0016a; the vector method still retains every node.

The vector method itself lives in `ext/diagnostics_vector_adapter.jl`,
`include`d at the end of this module (a continuation of it, not a
submodule — it shares this module's `using`s and private helpers).
Geometry, midpoint accounting, scalar selection and the vector adapter
live in separate files, each below CLAUDE.md Rule 6's 200-line cap.

## Why this is an extension, not core

`DelaunayTriangulation.jl` carries `ExactPredicates`, `AdaptivePredicates`,
and `EnumX` as transitive deps.  Loading them eagerly would more than
double `using PadeTaylor`'s precompile time for users who only want a
Stage-1 walk.  ADR-0003 documents the precedent (Arblib, CommonSolve,
Makie); ADR-0016 ties this extension to it.

## References

  - ADR-0016 — `docs/adr/0016-diagnostics-extension.md` (this design).
  - ADR-0025 Amendment 3 §"D-VC7 scope" / Amendment 7 — the vector
    adapter spec and the measured before/after loop-closure numbers.
  - Probe — `external/probes/loop-closure-fig1/probe.jl` (scalar
    algorithmic source); `external/probes/v8b-baseline/probe.jl` §2b
    (the vector-adapter reference implementation).
  - Probe verdict — `external/probes/loop-closure-fig1/REPORT.md:79-98`;
    `external/probes/v8b-baseline/REPORT.md` §2b (the V8b baseline).
  - FFW 2017 §2.1.2 — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-103`.
"""
module PadeTaylorDiagnosticsExt

using PadeTaylor: PathNetworkSolution, VectorPathNetworkSolution
using PadeTaylor.Diagnostics: DiagnosticReport, EdgeReport, quality_diagnose
using PadeTaylor.PathNetwork: _evaluate_pade
using DelaunayTriangulation: triangulate, each_edge
using Statistics: median, quantile, mean

import PadeTaylor.Diagnostics: quality_diagnose

const _EPS_FLOOR = 1e-300

include("diagnostics_geometry.jl")
include("diagnostics_edges.jl")
include("diagnostics_scalar.jl")
include("diagnostics_vector_adapter.jl")

end # module PadeTaylorDiagnosticsExt
