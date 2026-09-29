"""
    PadeTaylor.Diagnostics

Quality-certificate layer for `PathNetworkSolution` outputs.  Bead
`padetaylor-5t4`; promoted from the FFW 2017 Fig 1 loop-closure probe
shipped in `5a4d0a7` and documented at
`external/probes/loop-closure-fig1/REPORT.md:1-156`.

## What this module is for

`path_network_solve` walks a **tree** rooted at the IC (FW 2011 §3.1).
Two visited nodes that are geometrically adjacent in the solve's plane — i.e.
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
  - **Sheet 0 only at v1 (ADR-0016a).** Branchless walks retain every
    node. Empty `visited_sheet` metadata does not identify a ζ-frame:
    the strip in FFW 2017 md:103 belongs to the specific map z=exp(ζ/2).
    Populated metadata retains `visited_sheet[k] == [0]`. Any node loss
    throws unless explicitly budgeted by `min_retained_fraction`.
    `DiagnosticReport` records retained/dropped nodes, candidate/evaluated/
    dropped edges, pole exceptions, nonfinite results, and excluded
    off-sheet parent links. Multi-sheet Delaunay support remains deferred
    to `padetaylor-8py`; `n_branch_cut` is reserved and always 0 in v1.

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

include("diagnostics_reports.jl")

"""
    quality_diagnose(sol::PathNetworkSolution; sheet=0, tol_well=1e-10,
                     tol_bad=1e-6, n_worst=10,
                     min_retained_fraction=1.0) -> DiagnosticReport

Compute a loop-closure quality certificate on `sol`. Branchless walks
retain every visited node. Branched walks select `visited_sheet[k] == [0]`;
`sheet != 0` remains unsupported (padetaylor-8py). The default
`min_retained_fraction=1.0` throws on any node loss with a suggestion.
For an intentional sheet-0 subset, explicitly supply a lower fraction in
[0,1]; the report still records retained and dropped counts. No coordinate
strip is inferred. Delaunay-triangulate the retained cloud, extract
non-tree edges, and record midpoint disagreement. Pole exceptions and
nonfinite evaluations count as dropped edges; unexpected exceptions and
all-candidate failure throw. Quantiles describe evaluated edges only.
See the module docstring for category thresholds and ADR-0016a for coverage.

This generic is **method-less in core PadeTaylor**: the Delaunay-backed
implementation lives in `ext/PadeTaylorDiagnosticsExt.jl` and activates
when a package that loads `DelaunayTriangulation.jl` is present
alongside `PadeTaylor`.  Without that load, calling this function
surfaces a `MethodError`; the `path_network_solve(diagnose=true)`
entry point catches that and rethrows with an explicit `using
DelaunayTriangulation` suggestion.

## The vector method — `quality_diagnose(::VectorPathNetworkSolution)`

The extension also carries an **additive** method for the v0.2 vector
path-network walk (`VectorPathNetworkSolution`), the VC-7 criterion of
ADR-0025 (Amendment 3 §"D-VC7 scope"; Amendment 7).  It differs from
the scalar method in exactly three ways: it evaluates each node's
shared-`Q` approximant `Pᵢ(t)/Q(t)` (not a scalar `PadeApproximant`),
it generalises `ΔP_rel` to the vector 2-norm
`‖y_A − y_B‖ / (‖y_A‖ + ‖y_B‖ + ε)`, and it takes **no** `sheet` or
`min_retained_fraction` kwarg — the `P_I⁽²⁾` companion system is
single-sheeted, so there is no
`visited_sheet` field and the report's `sheet` is fixed at `0`.  The
Delaunay / tree-edge / LCA / categorise / aggregate machinery is
identical.  See the extension module's docstring for the full design.
"""
function quality_diagnose end

end # module Diagnostics
