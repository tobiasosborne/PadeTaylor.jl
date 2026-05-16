# ADR-0016 — Loop-closure diagnostics as a Delaunay-backed extension

**Status**: Accepted (2026-05-16) | **Bead**: `padetaylor-5t4` | **Worklog**: 048

## Decision

Promote the FFW 2017 Fig 1 loop-closure probe (shipped as
investigation-only in commit `5a4d0a7` at
`external/probes/loop-closure-fig1/probe.jl`) to a first-class
quality certificate `quality_diagnose(sol)::DiagnosticReport` on
`PathNetworkSolution`.  Wire it in via a Julia 1.9+ **weak-dep
package extension** (`DelaunayTriangulation` →
`ext/PadeTaylorDiagnosticsExt.jl`), per the precedent set in ADR-0003
for `Arblib`, `CommonSolve`, and `Makie`.  Add an **eager opt-in
kwarg** `diagnose::Bool = false` on `path_network_solve`: when `true`,
post-solve attach the `DiagnosticReport` to the returned solution's
new `diagnostics` field.  Sheet 0 only at v1; multi-sheet deferred to
bead `padetaylor-8py`.

## Context

The path-network walker (`src/PathNetwork.jl::path_network_solve`)
builds a **tree** rooted at the IC (FW 2011 §3.1).  Two visited nodes
that are geometrically adjacent (Delaunay neighbours in the ζ-plane)
but tree-distant — their LCA sits deep in the IC — carry
*independent* accumulated truncation error along their respective
IC-to-node paths.  The local adaptive controller (`:adaptive_ffw`)
sees only one step at a time; it has no view of long-range loop
closure.

The probe at `external/probes/loop-closure-fig1/probe.jl` (597 LOC)
measured the per-edge midpoint disagreement
`ΔP_rel := |P_A(M) - P_B(M)| / (|P_A(M)| + |P_B(M)| + ε)` on every
non-tree Delaunay edge of FFW Fig 1's sheet-0 walk.  The findings
(`external/probes/loop-closure-fig1/REPORT.md:34-77`):

  - 1689 non-tree Delaunay edges on the 850-node sheet-0 sub-walk.
  - `median(ΔP_rel) = 2.57e-12`; `p99 = 1.06e-01`; `max = 9.97e-01`.
  - **Trimodal** distribution: ~480-edge lobe at `1e-16…1e-13`
    (machine-eps closure), ~1000-edge lobe at `1e-13…1e-10`
    (controller-tolerance closure), ~200-edge tail at `1e-6…1e0`
    (catastrophic).
  - 6.3 % of loop closures are visibly catastrophic (`ΔP_rel > 1e-3`).
  - Catastrophic edges cluster in the high-Re ζ corner where `R(ζ)`
    already loses resolution; several top-10 worst edges have
    `extrap_max > 5` (midpoints outside both endpoints' canonical
    Padé discs).

The probe's verdict (REPORT.md:79-98): *partially confirmed* — the
consistency-check half of the original hypothesis is vindicated;
the correction-signal half is harder because at `|t| > 5` neither
Padé is trustworthy.

The investigation was useful enough on its own to want a permanent,
caller-accessible diagnostic — not just for the FFW Fig 1 walk, but
for every Painlevé solve.  Hence this ADR.

## Architecture

### Weak-dep extension

`DelaunayTriangulation.jl` is the algorithmic core of the diagnostic;
it carries `ExactPredicates`, `AdaptivePredicates`, `EnumX` as
transitive deps.  Loading it eagerly would more than double
`using PadeTaylor`'s precompile time for users who only want a
Stage-1 walk.  ADR-0003 sets the precedent (Arblib, CommonSolve,
Makie); we follow it verbatim:

```toml
[weakdeps]
DelaunayTriangulation = "927a84f5-c5f4-47a5-9785-b46e178433df"
[extensions]
PadeTaylorDiagnosticsExt = "DelaunayTriangulation"
[compat]
DelaunayTriangulation = "1"
```

Core `src/Diagnostics.jl` declares `DiagnosticReport`, `EdgeReport`,
and the empty generic `quality_diagnose`.  Method-on-`PathNetworkSolution`
lives in `ext/PadeTaylorDiagnosticsExt.jl`.

### Eager opt-in `diagnose::Bool = false` kwarg

`path_network_solve(...; diagnose = true)` calls
`quality_diagnose(sol)` post-solve and rebuilds the returned
solution with the report attached to a new
`diagnostics::Union{Nothing, DiagnosticReport}` field (field 11).
Default `false` keeps `sol.diagnostics === nothing` and preserves
every existing test invariant byte-for-byte.

If `diagnose=true` is passed but the extension is not loaded, the
`MethodError` from the empty generic is caught and rethrown as an
`ArgumentError` carrying the explicit `using DelaunayTriangulation`
suggestion (CLAUDE.md Rule 1).

### Sheet 0 only at v1

The probe ran on FFW Fig 1, which uses no branch points, so
`visited_sheet[k]` is always empty there.  The probe used the FFW
ζ-strip predicate `-2π < imag(z) ≤ 2π` to identify sheet 0
(`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:103`).  We generalise
slightly: when `visited_sheet[k]` is non-empty (branched solve), the
sheet-0 filter becomes `visited_sheet[k] == [0]`.

Multi-sheet diagnostics require a cut-aware Delaunay step (edges
that cross a branch cut must not be paired with same-sheet
midpoints).  That work is deferred to bead `padetaylor-8py`; the
`n_branch_cut` field on `DiagnosticReport` is reserved for it and
is always `0` in v1.  The `quality_diagnose` method throws
`ArgumentError` if called with `sheet != 0` until that bead lands.

### Edge categories

Each non-tree Delaunay edge is bucketed into one of four v1
categories (the fifth, `:branch_cut`, is reserved for v2):

| Category          | Predicate                                                  |
|-------------------|------------------------------------------------------------|
| `:well_closed`    | `ΔP_rel ≤ tol_well` (default `1e-10`)                      |
| `:noisy`          | `tol_well < ΔP_rel ≤ tol_bad` (default `tol_bad = 1e-6`)   |
| `:extrap_driven`  | `ΔP_rel > tol_bad` AND `max(|t_A|, |t_B|) > 1`             |
| `:depth_driven`   | `ΔP_rel > tol_bad` AND both endpoints in-disc              |
| `:branch_cut`     | reserved for v2; always 0 in v1                            |

The `:depth_driven` bucket is the honest "graph consensus" signal —
two independently-walked Padé patches disagree even though both are
in their canonical discs.  A future graph-consensus Stage-2 pass
(out of scope; tracked via bead `padetaylor-zwh` for the orthogonal
Poisson-disk node-placement direction) would target this bucket.

## Alternatives considered, rejected

  - **Hard dep on `DelaunayTriangulation`.**  Rejected: balloons
    `using PadeTaylor` precompile time for users who don't want
    diagnostics; violates the ADR-0003 lean-core principle.
  - **Lazy-only API (no `diagnose` kwarg; only `quality_diagnose(sol)`).**
    Rejected: defeats the "make the diagnostic visible at the call
    site" motivation.  The eager opt-in lets `figure scripts` and
    `bd-claim` runs report quality without extra boilerplate.
  - **k-NN edges instead of Delaunay.**  Rejected: the probe
    deliberately chose Delaunay because every non-tree edge
    corresponds to a *cycle* in the visited-node graph, and the
    `nontree = Delaunay \ tree` set is exactly the loop-closure
    population.  k-NN would over-count short edges and miss the
    convex-hull boundary edges that turned out to carry the
    catastrophic signal (REPORT.md:60-71).
  - **Symbol-typed `policy` instead of boolean `diagnose`.**
    Rejected for v1: only one policy exists today (compute the
    full report).  Future v2 policies (e.g. compute only the
    centroid for a fast spatial cue) can be added without a
    breaking change by introducing a `Symbol`-typed kwarg
    alongside the bool — the bool is opt-in, so deprecation is
    smooth.

## Consequences

### User-visible

  - New public API: `DiagnosticReport`, `EdgeReport`,
    `quality_diagnose`, all exported from `PadeTaylor`.
  - New `diagnostics` field on `PathNetworkSolution`; populated when
    `diagnose=true` is passed, `nothing` otherwise.
  - New kwarg `diagnose::Bool = false` on `path_network_solve`.
  - User friction: activating the diagnostic requires
    `using DelaunayTriangulation`.  The error path on the missing
    load is explicit per Rule 1.

### Test coverage

A separate sub-agent task ships the test suite (target `test/diagnose_test.jl`):
  - `diagnose=false` (default): `sol.diagnostics === nothing`.
  - `diagnose=true` without `using DelaunayTriangulation`: throws
    `ArgumentError` carrying the suggestion.
  - `diagnose=true` with extension loaded on a small fixture:
    `sol.diagnostics::DiagnosticReport`; categories sum to `n_edges`.
  - Post-hoc `quality_diagnose(sol)` matches the eager-opt-in result
    bit-for-bit.
  - Sheet-0 mask path: branched solve with `visited_sheet[k] == [0]`
    yields the same edge population as the equivalent ζ-strip path.
  - Mutation-proof: perturb `_categorise`'s `extrap_max` threshold
    or `_tree_path_distance`'s LCA loop and the test goes RED.

### Documentation

  - This ADR-0016.
  - `src/Diagnostics.jl` top-of-file literate docstring (Rule 10).
  - `ext/PadeTaylorDiagnosticsExt.jl` top-of-file literate docstring.
  - `src/PathNetwork.jl`: new `diagnose` kwarg documented inline;
    `PathNetworkSolution` docstring updated for the new field;
    backward-compat constructor commentary updated.
  - `src/PadeTaylor.jl`: architecture list updated (item 8); new
    re-exports and exports.
  - Worklog 048.

### Open follow-ups

  - **Bead `padetaylor-8py` — multi-sheet diagnostics.**  Lift the
    `sheet == 0` restriction.  Requires cut-aware Delaunay edge
    filtering: an edge whose two endpoints sit on different sheets,
    or whose midpoint sits on the opposite side of a branch cut from
    either endpoint, must be marked `:branch_cut` rather than
    evaluated.  Touches `_sheet_mask` and adds a cut-crossing
    predicate to the per-edge loop.

  - **Graph-consensus Stage-2 (longer-term).**  The `:depth_driven`
    bucket flags the population a consensus-Stage-2 pass would
    target.  Orthogonal to `padetaylor-zwh` (Poisson-disk node
    placement), which targets the `:extrap_driven` bucket by
    reducing extrapolation rather than averaging it out.

## References

  - **Probe** — `external/probes/loop-closure-fig1/probe.jl` (the
    algorithmic source we promote here).
  - **Probe verdict** — `external/probes/loop-closure-fig1/REPORT.md:79-98`.
  - **FFW 2017 §2.1.2** —
    `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-103`
    (Stage-1 adaptive controller spec and the ζ-strip sheet definition).
  - **ADR-0003** — `docs/adr/0003-extensions-pattern.md` (the
    weak-dep precedent we follow).
  - **ADR-0004** — `docs/adr/0004-path-network-architecture.md`
    (the path-network tree on which loop closure is measured).
  - **ADR-0013** — `docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md`
    (the `visited_sheet` field the sheet-0 mask consults when populated).
  - **CLAUDE.md** Rule 1 (fail loud with suggestion) — the
    extension-missing error path.
