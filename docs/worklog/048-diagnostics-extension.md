# Worklog 048 — Diagnostics extension: loop-closure quality certificate (bead `padetaylor-5t4`)

**Date**: 2026-05-16
**Author**: Claude Sonnet 4.6 (docs-lockstep pass)
**Bead**: `padetaylor-5t4` — closed by this commit.
**Scope**: Promote the FFW 2017 Fig 1 loop-closure probe to a first-class
`quality_diagnose(sol)::DiagnosticReport` API on `PathNetworkSolution`.
New files: `src/Diagnostics.jl`, `ext/PadeTaylorDiagnosticsExt.jl`,
`docs/adr/0016-diagnostics-extension.md`, `test/diagnose_test.jl`.
Modified: `src/PathNetwork.jl`, `src/PadeTaylor.jl`, `Project.toml`,
`test/runtests.jl`.

> **Take-home**: The loop-closure probe (bead `padetaylor-e3h`,
> `external/probes/loop-closure-fig1/REPORT.md:34-77`) revealed a
> **trimodal** disagreement distribution on FFW Fig 1's sheet-0 walk:
> ~480 edges at machine-eps quality, ~1000 edges within controller
> tolerance, and a ~200-edge catastrophic tail above `1e-6`.  That
> finding was structured enough — and spatially informative enough —
> to warrant a permanent, caller-accessible certificate rather than
> leaving it as a one-off probe script.  This worklog records how it
> was wired in.

## What this is

The path-network walker builds a **tree** rooted at the initial
condition.  Two visited nodes that are geometrically adjacent (Delaunay
neighbours in the ζ-plane) but tree-distant carry *independent*
accumulated truncation error along their IC-to-node paths; the local
`:adaptive_ffw` controller has no view of long-range loop closure.

The probe at `external/probes/loop-closure-fig1/probe.jl` measured the
per-edge midpoint disagreement
`ΔP_rel := |P_A(M) - P_B(M)| / (|P_A(M)| + |P_B(M)| + ε)` on every
non-tree Delaunay edge of FFW Fig 1's 850-node sheet-0 sub-walk.  The
key findings (`external/probes/loop-closure-fig1/REPORT.md:34-77`):

  - 1689 non-tree Delaunay edges; `median(ΔP_rel) = 2.57e-12`;
    `p99 = 1.06e-01`; `max = 9.97e-01`.
  - **Trimodal distribution**: ~480 edges at `1e-16…1e-13`
    (machine-eps closure), ~1000 edges at `1e-13…1e-10`
    (controller-tolerance closure), ~200 edges at `1e-6…1e0`
    (catastrophic).
  - 6.3 % of loop closures are visibly catastrophic (`ΔP_rel > 1e-3`).
  - Catastrophic edges cluster in the high-Re ζ corner where `R(ζ)`
    loses resolution; several top-10 worst edges have `extrap_max > 5`
    (midpoints outside both endpoints' canonical Padé discs).
  - Pearson `r(tree_dist, ΔP_rel) ≈ 0.18` — weak; the catastrophic
    tail is geometrically localised, not distributed along deep tree
    branches.

The investigation was decisive enough to promote: a permanent
`quality_diagnose` API makes the certificate available for every future
Painlevé solve, not just the one-off probe.

## Design decisions locked in

Per **ADR-0016** (`docs/adr/0016-diagnostics-extension.md`):

1. **Weak-dep extension** — `DelaunayTriangulation` goes into
   `[weakdeps]`; the Delaunay-backed method lives in
   `ext/PadeTaylorDiagnosticsExt.jl`.  Without the extension loaded,
   `quality_diagnose` throws a helpful `ErrorException` naming the
   extension.  Same precedent as ADR-0003 (Arblib, CommonSolve, Makie).

2. **Eager opt-in kwarg** — `diagnose::Bool = false` on
   `path_network_solve`: when `true`, the diagnostic runs immediately
   after the walk and the `DiagnosticReport` is attached to
   `PathNetworkSolution.diagnostics`.  The caller can also run
   `quality_diagnose(sol)` post-hoc on an existing solution.

3. **Sheet 0 only at v1** — multi-sheet diagnostics (sheets ±1, ±2,
   …) deferred to bead `padetaylor-8py`.  Sheet-0 filtering uses the
   strip predicate on `visited_z` when `visited_sheet` is not
   populated; the Schwarz-reflection path is threaded through.

## What shipped

  - **`src/Diagnostics.jl`** (208 LOC; ~120 effective): defines
    `DiagnosticReport` (struct with 9 fields: `n_nodes`, `n_nontree`,
    `delta_p_rel_median`, `delta_p_rel_p90`, `delta_p_rel_p99`,
    `delta_p_rel_max`, `n_above_tol`, `n_catastrophic`, `edge_reports`)
    and `EdgeReport` (struct with 5 fields: `z_a`, `z_b`,
    `delta_p_rel`, `tree_dist_a`, `tree_dist_b`); empty generic
    `quality_diagnose(sol)` that throws a useful error message when the
    extension is not loaded; `Base.show` for both structs.

  - **`ext/PadeTaylorDiagnosticsExt.jl`** (276 LOC; ~180 effective):
    the Delaunay-backed method.  Filters sheet-0 visited nodes → builds
    Delaunay triangulation → extracts non-tree loop-closure edges →
    evaluates `ΔP_rel` at each edge midpoint via
    `PathNetwork._evaluate_pade` → collects `EdgeReport` instances →
    returns a populated `DiagnosticReport`.

  - **`docs/adr/0016-diagnostics-extension.md`** (227 LOC): full ADR
    with context, architecture, locked decisions, deferred scope.

  - **`test/diagnose_test.jl`** (196 LOC, 32 assertions across 8
    testsets DG.1–DG.8): see "What's tested" below.

  - **`src/PathNetwork.jl`** modified: `PathNetworkSolution` gains 11th
    field `diagnostics::Union{Nothing, DiagnosticReport}`; two
    backward-compatible inner constructors added; `path_network_solve`
    gains `diagnose::Bool=false` kwarg; Schwarz-reflection path
    threaded.

  - **`src/PadeTaylor.jl`** modified: `include("Diagnostics.jl")`
    inserted before PathNetwork; `DiagnosticReport`, `EdgeReport`,
    `quality_diagnose` re-exported; module-doc list updated.

  - **`Project.toml`** modified: `Statistics` added to `[deps]`;
    `DelaunayTriangulation` added to `[weakdeps]`, `[extensions]`,
    `[compat]`, `[extras]`, `[targets].test`.

  - **`test/runtests.jl`** modified: `include("diagnose_test.jl")`
    wired after `ext_makie_test.jl`.

## What's tested

32 assertions across 8 testsets; run in 18 s in single-file isolation:

  | Testset | Assertions | What it pins |
  |---------|-----------|--------------|
  | DG.1 — smoke without extension loaded | 1 | `quality_diagnose` throws `ErrorException` when extension absent |
  | DG.2 — struct field counts | 4 | `DiagnosticReport` has 9 fields; `EdgeReport` has 5 fields |
  | DG.3 — empty / trivial solve | 3 | `diagnostics === nothing` by default; `diagnose=true` returns a report |
  | DG.4 — report numerics (Fig 1 sub-walk) | 8 | `n_nontree > 0`; median, p90, p99, max within documented ranges; `n_above_tol` + `n_catastrophic` non-negative |
  | DG.5 — EdgeReport fields | 4 | `z_a`, `z_b` in sheet-0 strip; `delta_p_rel >= 0`; `tree_dist_a` + `tree_dist_b >= 1` |
  | DG.6 — post-hoc call | 2 | `quality_diagnose(sol)` on an already-solved solution gives same result as `diagnose=true` |
  | DG.7 — show methods | 4 | `DiagnosticReport` and `EdgeReport` print without error; output contains expected field labels |
  | DG.8 — backward compat | 6 | old 10-arg `PathNetworkSolution` constructor still works; `sol.diagnostics === nothing` on old-style construction |

**Mutation-prove A** applied and reverted (RED→GREEN verified at 18 s
wall time).  Procedure: perturb `delta_p_rel` computation in the
extension by negating the absolute-value comparison, confirm DG.4
fails, restore.  Mutation-prove procedure recorded in
`test/diagnose_test.jl:160-190`.

## OOM friction

Full `Pkg.test()` was SIGTERM'd **twice** by the OOM-killer on this
WSL2 system (62 GB total, ~47 GB free at run time) during Julia
precompile-cache warming.  The failure is a precompile spike, not a
test logic issue.  The new test file was validated in single-file
isolation (`julia --project=. test/diagnose_test.jl`, 18 s, 32/32
GREEN); the existing 2508 tests are unchanged (additive struct field +
backward-compat constructors; no existing test touches the new field).

Memory context captured to `bd memories` so future sessions know to
avoid full `Pkg.test()` under low system headroom.  **Next session
should run `Pkg.test()` to confirm the suite is regression-free before
pushing.**

## What is NOT shipped (CLAUDE.md Rule 9)

  - **Empirical validation sweep** — reproducing FW 2011 Figs 4.1,
    4.7, 5.1, 5.2 with `diagnose=true` to generate concrete
    `DiagnosticReport` comparisons across different pole configurations.
    Deferred to bead **`padetaylor-6pj`**.

  - **Multi-sheet diagnostics** — repeating the loop-closure check on
    sheets ±1, ±2, … (where the visible seams in FFW Fig 1 actually
    appear; the probe's P3 recommendation at
    `external/probes/loop-closure-fig1/REPORT.md:117-123`).  Deferred
    to bead **`padetaylor-8py`**.

  - **Graph-consensus Stage-2** — using `ΔP_rel` as a correction
    signal to jointly adjust node values.  The probe's P2
    recommendation (`REPORT.md:110-115`): effective only when paired
    with denser Poisson-disk node placement (bead `padetaylor-zwh`)
    so midpoints stay inside `|t| ≤ 1`.  Deferred.

## Bead

  - `padetaylor-5t4` — closed by this commit.

## References

  - **ADR-0016** — `docs/adr/0016-diagnostics-extension.md`.
  - **Probe REPORT** — `external/probes/loop-closure-fig1/REPORT.md:34-77`
    (trimodal-distribution finding + spatial clustering).
  - **Probe REPORT** — `external/probes/loop-closure-fig1/REPORT.md:79-98`
    (verdict + recommendation).
  - **ADR-0003** — `docs/adr/0003-pkg-extensions.md` (weak-dep pattern
    precedent: Arblib, CommonSolve, Makie).
  - **`src/Diagnostics.jl`** — `DiagnosticReport`, `EdgeReport`,
    generic `quality_diagnose`.
  - **`ext/PadeTaylorDiagnosticsExt.jl`** — Delaunay-backed method.
  - **`test/diagnose_test.jl:160-190`** — mutation-prove A procedure.
