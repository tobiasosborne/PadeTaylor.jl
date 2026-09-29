"""
    EdgeReport

This chapter defines the edge/report containers and their text display,
included in Diagnostics. ADR-0016a records the analysed population beside
its quality statistics; the full certificate follows this edge container.

One non-tree Delaunay edge of the visited-node cloud.  Fields:

  - `A`, `B`             — global indices into the parent
                            `PathNetworkSolution`'s `visited_*` arrays.
  - `ΔP_abs`             — `|P_A(M) - P_B(M)|` at the midpoint `M`.
  - `ΔP_rel`             — `ΔP_abs / (|P_A(M)| + |P_B(M)| + ε)`.
  - `tree_dist`          — number of tree edges on the path A↔B
                            (LCA-based; tree-distant edges are the
                            interesting loop-closure population).
  - `extrap_max`         — `max(|t_A|, |t_B|)` where
                            `t_X = (M - z_X) / visited_h[X]`.  Values
                            `> 1` indicate the edge midpoint sits
                            outside one (or both) endpoints' canonical
                            Padé disc.
  - `midpoint`           — `M = (z_A + z_B) / 2`.
  - `category`           — `:well_closed | :noisy | :extrap_driven |
                            :depth_driven | :branch_cut`; see the
                            module docstring for thresholds.
"""
struct EdgeReport
    A          :: Int
    B          :: Int
    ΔP_abs     :: Float64
    ΔP_rel     :: Float64
    tree_dist  :: Int
    extrap_max :: Float64
    midpoint   :: ComplexF64
    category   :: Symbol
end

"""
    DiagnosticReport

Loop-closure quality certificate for a `PathNetworkSolution`.

The five category tallies sum to `n_edges` on the analysed sheet. Quantiles
(`median_ΔP_rel`, `p90_ΔP_rel`, `p99_ΔP_rel`, `max_ΔP_rel`) are over
all evaluated non-tree edges.  `worst_edges` carries the top-N
offenders sorted by `ΔP_rel` descending; `n_worst` defaults to 10.

`bad_centroid` is the arithmetic mean of midpoints of edges with
`ΔP_rel > tol_bad`; `NaN+NaN·im` when none exist.  Useful for figure
scripts that want to circle the catastrophic region.

`sheet` records which sheet was analysed (`0` in v1; see the module
docstring's "Sheet 0 only" note).  `tol_well` and `tol_bad` echo the
thresholds the report was computed at, so a serialised report stays
self-describing.

Coverage is explicit: `n_nodes_retained + n_nodes_dropped == n_nodes`;
`n_edges + n_edges_dropped == n_edges_candidate`; and
`n_edges_dropped == n_edges_threw + n_edges_nonfinite`. Pole exceptions
count as `n_edges_threw`; nonfinite values or disagreement statistics
count as `n_edges_nonfinite`, once per edge. Unexpected errors propagate.
`n_tree_edges_off_sheet` counts parent links excluded by sheet selection;
these are not candidate non-tree edges and are outside that denominator.
No cross-sheet Delaunay graph is constructed in v1 (padetaylor-8py).
"""
struct DiagnosticReport
    n_edges         :: Int
    n_well_closed   :: Int
    n_noisy         :: Int
    n_extrap_driven :: Int
    n_depth_driven  :: Int
    n_branch_cut    :: Int
    median_ΔP_rel   :: Float64
    p90_ΔP_rel      :: Float64
    p99_ΔP_rel      :: Float64
    max_ΔP_rel      :: Float64
    worst_edges     :: Vector{EdgeReport}
    bad_centroid    :: ComplexF64
    sheet           :: Int
    tol_well        :: Float64
    tol_bad         :: Float64
    n_nodes         :: Int
    n_nodes_retained :: Int
    n_nodes_dropped  :: Int
    n_edges_candidate :: Int
    n_edges_dropped  :: Int
    n_edges_threw    :: Int
    n_edges_nonfinite :: Int
    n_tree_edges_off_sheet :: Int
end

# Compact text/plain summary.  Keep narrow (≤80 cols) so the default
# REPL print stays legible after `display(sol.diagnostics)`.
function Base.show(io::IO, ::MIME"text/plain", r::DiagnosticReport)
    println(io, "DiagnosticReport — sheet $(r.sheet), $(r.n_edges) non-tree Delaunay edges")
    println(io, "  retained nodes  : $(r.n_nodes_retained) / $(r.n_nodes) (dropped=$(r.n_nodes_dropped))")
    println(io, "  evaluated edges : $(r.n_edges) / $(r.n_edges_candidate) (dropped=$(r.n_edges_dropped))")
    println(io, "  failed edges    : threw=$(r.n_edges_threw), nonfinite=$(r.n_edges_nonfinite)")
    println(io, "  off-sheet tree  : $(r.n_tree_edges_off_sheet) excluded parent links")
    pct(n) = r.n_edges == 0 ? 0.0 : 100 * n / r.n_edges
    println(io, "  well_closed     : $(r.n_well_closed) ($(round(pct(r.n_well_closed); digits=1))%)  (ΔP_rel ≤ $(r.tol_well))")
    println(io, "  noisy           : $(r.n_noisy) ($(round(pct(r.n_noisy); digits=1))%)")
    println(io, "  extrap_driven   : $(r.n_extrap_driven) ($(round(pct(r.n_extrap_driven); digits=1))%)  (|t| > 1 at midpoint)")
    println(io, "  depth_driven    : $(r.n_depth_driven) ($(round(pct(r.n_depth_driven); digits=1))%)  (in-disc loop-closure failure)")
    println(io, "  branch_cut      : $(r.n_branch_cut)  (reserved — v1 sheet-0 only)")
    println(io, "  median ΔP_rel   : $(r.median_ΔP_rel)")
    println(io, "  p99 ΔP_rel      : $(r.p99_ΔP_rel)")
    print(io,   "  bad centroid    : $(r.bad_centroid)")
end
