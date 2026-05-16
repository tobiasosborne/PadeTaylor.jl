"""
    PadeTaylor.Diagnostics

Quality-certificate layer for `PathNetworkSolution` outputs.  Bead
`padetaylor-5t4`; promoted from the FFW 2017 Fig 1 loop-closure probe
shipped in `5a4d0a7` and documented at
`external/probes/loop-closure-fig1/REPORT.md:1-156`.

## What this module is for

`path_network_solve` walks a **tree** rooted at the IC (FW 2011 §3.1).
Two visited nodes that are geometrically adjacent in the ζ-plane — i.e.
neighbours under the Delaunay triangulation of the visited-node cloud —
but tree-distant (their LCA sits deep in the IC) carry *independent*
accumulated truncation error along their respective IC-to-node paths.
The local adaptive controller (`:adaptive_ffw`) sees only one step at a
time; it has no view of long-range loop closure.

The probe at `external/probes/loop-closure-fig1/probe.jl` measured the
per-edge midpoint disagreement
`ΔP_rel := |P_A(M) - P_B(M)| / (|P_A(M)| + |P_B(M)| + ε)` on every
non-tree Delaunay edge of FFW Fig 1's sheet-0 walk and reported a
**trimodal** distribution (REPORT.md:55-58): a machine-eps lobe, a
controller-tolerance lobe, and a ~6 % catastrophic tail clustered at
high-Re ζ where stored Padés must extrapolate past their canonical
disc to reach midpoints.  That probe was investigation-only.  This
module promotes those findings to a first-class quality certificate:
`quality_diagnose(sol)` returns a `DiagnosticReport` summarising the
loop-closure signal and flagging the worst offenders for the caller's
attention.

## Design (ADR-0016)

  - **Weak-dep extension.**  The Delaunay-triangulation step relies on
    `DelaunayTriangulation.jl`, which we wire in via a Julia 1.9+
    package extension (same precedent as `Arblib`, `CommonSolve`,
    `Makie` — see ADR-0003).  This module declares the data containers
    and an empty generic `quality_diagnose`; the heavy lifting lives
    in `ext/PadeTaylorDiagnosticsExt.jl` and is activated by `using
    DelaunayTriangulation` alongside `PadeTaylor`.  Calling
    `quality_diagnose` without that load surfaces a `MethodError`
    that the `path_network_solve(diagnose=true)` entry point catches
    and rethrows with an explicit `using DelaunayTriangulation`
    suggestion (CLAUDE.md Rule 1).
  - **Eager opt-in.**  `path_network_solve(...; diagnose=true)` calls
    `quality_diagnose(sol)` post-solve and attaches the
    `DiagnosticReport` to the returned solution's `diagnostics` field.
    The default `diagnose=false` keeps `diagnostics === nothing`,
    preserving every existing test invariant byte-for-byte.
  - **Sheet 0 only at v1.**  Mirrors the probe verbatim (REPORT.md:7-32):
    when `visited_sheet[k]` is empty (`branch_points = ()` solves),
    sheet 0 is defined by the FFW ζ-strip predicate `-2π < imag(z) ≤
    2π` (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
    FFW2017_painleve_riemann_surfaces_preprint.md:103`); when populated,
    sheet 0 is `visited_sheet[k] == [0]`.  Multi-sheet support (which
    requires explicit branch-cut bookkeeping for the Delaunay step)
    is deferred to bead `padetaylor-8py`; the `n_branch_cut` field
    on `DiagnosticReport` is reserved for it and is always `0` in v1.

## Edge categories

Each non-tree Delaunay edge is classified by its midpoint disagreement
and by whether *either* endpoint's stored Padé had to extrapolate past
its canonical disc (`|t| > 1`) to reach the midpoint:

  - `:well_closed`     — `ΔP_rel ≤ tol_well` (default `1e-10`).  Both
                         endpoints' Padés agree to controller tolerance
                         at the midpoint; loop closes.
  - `:noisy`           — `tol_well < ΔP_rel ≤ tol_bad` (default
                         `tol_bad = 1e-6`).  Loop closes only to coarser
                         tolerance; usually long tree paths or moderate
                         extrapolation.
  - `:extrap_driven`   — `ΔP_rel > tol_bad` AND `max(|t_A|, |t_B|) > 1`.
                         At least one endpoint extrapolated past its
                         disc; the disagreement combines honest
                         tree-divergence with Padé extrapolation
                         amplification.  Denser sampling (e.g. Poisson-
                         disk Stage-1 nodes, bead `padetaylor-zwh`) is
                         the right cure.
  - `:depth_driven`    — `ΔP_rel > tol_bad` AND both endpoints inside
                         their canonical discs.  The honest "graph
                         consensus" signal: two independently-walked
                         Padé patches disagree even though both are
                         in-disc.  A graph-consensus Stage-2 pass would
                         flag these as suspect.
  - `:branch_cut`      — reserved for v2 (multi-sheet); always `0`.

## References

  - FFW 2017 §2.1.2 — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-103`.
  - Probe — `external/probes/loop-closure-fig1/probe.jl` (the algorithmic source we promote here).
  - Probe verdict — `external/probes/loop-closure-fig1/REPORT.md:79-98`.
  - ADR-0016 — `docs/adr/0016-diagnostics-extension.md` (this design).
  - ADR-0003 — `docs/adr/0003-extensions-pattern.md` (the weak-dep precedent).
"""
module Diagnostics

export DiagnosticReport, EdgeReport, quality_diagnose

"""
    EdgeReport

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

Tallies (`n_*`) sum to `n_edges` on the analysed sheet.  Quantiles
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
end

"""
    quality_diagnose(sol::PathNetworkSolution; sheet=0, tol_well=1e-10,
                     tol_bad=1e-6, n_worst=10) -> DiagnosticReport

Compute a loop-closure quality certificate on `sol`.  Sheet 0 only at
v1: the function filters visited nodes to sheet 0 (via
`visited_sheet[k] == [0]` when the solve carried branch points, or the
FFW ζ-strip predicate `-2π < imag(z) ≤ 2π` when `visited_sheet[k]` is
empty), Delaunay-triangulates them, extracts non-tree edges, and
records `ΔP_rel` at each edge midpoint.  See the module docstring for
the category thresholds; bead `padetaylor-8py` tracks multi-sheet
support.

This generic is **method-less in core PadeTaylor**: the Delaunay-backed
implementation lives in `ext/PadeTaylorDiagnosticsExt.jl` and activates
when a package that loads `DelaunayTriangulation.jl` is present
alongside `PadeTaylor`.  Without that load, calling this function
surfaces a `MethodError`; the `path_network_solve(diagnose=true)`
entry point catches that and rethrows with an explicit `using
DelaunayTriangulation` suggestion.
"""
function quality_diagnose end

# Compact text/plain summary.  Keep narrow (≤80 cols) so the default
# REPL print stays legible after `display(sol.diagnostics)`.
function Base.show(io::IO, ::MIME"text/plain", r::DiagnosticReport)
    println(io, "DiagnosticReport — sheet $(r.sheet), $(r.n_edges) non-tree Delaunay edges")
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

end # module Diagnostics
