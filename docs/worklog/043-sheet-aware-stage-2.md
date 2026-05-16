# Worklog 043 — Sheet-aware Stage-2 (bead `padetaylor-hed` closed)

**Date**: 2026-05-16
**Author**: Claude Opus 4.7 (1M context)
**Bead**: `padetaylor-hed` — closed by this commit.
**Scope**: Step A5 of the 11-step FFW 2017 figure-reproduction plan.
**ADR**: none (additive extension to `path_network_solve`'s Stage-2
pathway, no new design pattern beyond what ADR-0013 already
introduced).

> **Take-home**: Two opt-in additions on top of the A4 substrate.
> First, `path_network_solve`'s Stage-2 nearest-visited lookup gains
> an opt-in `grid_sheet :: Vector{Vector{Int}}` kwarg that restricts
> the pool to visited nodes whose `visited_sheet[k]` equals the
> query's sheet.  Second, a new public accessor
> `eval_at_sheet(sol::PathNetworkSolution, z, sheet)` does the same
> per-point query post-hoc.  Default `grid_sheet = nothing`
> preserves byte-identical legacy behaviour; SA.1.1's 2 assertions
> verify the invariant.  Test count **2368 → 2398 GREEN** (+30
> assertions in `test/sheet_aware_stage2_test.jl`).  Three
> mutations (S1/S2/S3) bite distinct subsets (4/6/2 RED).  Unblocks
> B5 (FFW Figs 2/3/7).

## Ground truth read first

  - **ADR-0013** (`docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md`)
    §"Open follow-ups: A5 sheet-aware Stage-2" — the deferral this
    bead closes.  A4 ships `visited_sheet` population; A5 wires the
    Stage-2 lookup to honour it.

  - **Worklog 042 §"What is NOT shipped"** — `eval_at_sheet`
    accessor identified as needed for B5 figure scripts rendering
    per-sheet panels.

  - **`src/PathNetwork.jl` Stage-2 block** (~lines 587-625) — the
    existing Stage-2 loop that A5 extends.  `_nearest_visited`
    (lines 801-814) — the helper A5 specializes.

## Design synthesis

A5 has no new design choices — it's a natural sheet-restricted
extension of two existing primitives:

  1. **`_nearest_visited_on_sheet`** (~25 LOC) — same min-distance
     scan as `_nearest_visited`, with an early-continue on
     `visited_sheet[i] != sheet`.  Returns `0` if no matching-sheet
     visited node exists (sentinel for the fail-soft NaN path).
     Same lexicographic tiebreak (Re, then Im) for determinism.

  2. **`grid_sheet` kwarg** (~15 LOC validation + ~5 LOC dispatch
     in Stage-2 loop) — when supplied, each grid point's lookup
     calls `_nearest_visited_on_sheet`; when `nothing`, the
     unrestricted `_nearest_visited` is used.  Validation throws
     before the Stage-1 walk (which can take minutes).

  3. **`eval_at_sheet`** public accessor (~40 LOC + docstring) —
     same per-point Padé eval as Stage-2, with an explicit sheet
     argument.  Shape validation throws on length mismatch
     (`length(sheet) != length(visited_sheet[1])`).

The natural alternative was to make `grid_sheet` a struct field
on `PathNetworkSolution` so the sheet metadata travels with the
solution.  Rejected because:

  - It would be CALLER-side metadata, not WALKER-side fact — the
    walker doesn't pick grid_sheet; the caller specifies it per
    query.  Storing it in the solution conflates the two roles.

  - `eval_at_sheet` allows the caller to re-query a SOLVED tree
    with different sheets at zero extra walk cost, which a stored
    `grid_sheet` would forbid (you'd need to re-solve to change
    the sheet annotation).

  - Backward-compat: existing `PathNetworkSolution` construction
    sites (three external test fixtures + the Schwarz path) don't
    need any new positional arg.

## What shipped

  - **`src/PathNetwork.jl`** (+~85 LOC):
    - New kwarg `grid_sheet::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing`.
    - Validation block (length + per-element length checks; throws
      before Stage-1).
    - Stage-2 dispatch on `grid_sheet === nothing` to
      `_nearest_visited` (unchanged) vs `_nearest_visited_on_sheet`
      (new); fail-soft NaN output when `idx_v == 0`.
    - New helper `_nearest_visited_on_sheet` (~25 LOC).
    - New public accessor `eval_at_sheet(sol, z, sheet)` (~40 LOC
      + docstring); exported.
    - Docstring grows one kwarg paragraph for `grid_sheet`.

  - **`src/PadeTaylor.jl`** (+2 lines): `eval_at_sheet` re-export.

  - **`test/sheet_aware_stage2_test.jl`** (new, ~225 LOC): eight
    testsets SA.1.1–SA.1.8, 30 GREEN assertions; three mutations
    S1/S2/S3 documented in the footer with exact bite counts
    (4/6/2 RED).

  - **`test/runtests.jl`**: `include("sheet_aware_stage2_test.jl")`
    after the path_network_branch test.

  - **This worklog** (`docs/worklog/043-sheet-aware-stage-2.md`).

## Mutation-proof results

  | Mutation | Bite count (of 30 new assertions) | Where it hits |
  |----------|------------------------------------|----------------|
  | S1       | 4 RED                              | SA.1.3 (2), SA.1.8 (2) |
  | S2       | 6 RED                              | SA.1.3 (2), SA.1.5 (4) |
  | S3       | 2 RED                              | SA.1.5 only    |

S1 is the "helper itself is broken" case (drop sheet filter inside
`_nearest_visited_on_sheet`); S2 is the "Stage-2 dispatch ignores
the kwarg" case (always call the unrestricted helper); S3 is the
"eval_at_sheet is misrouted" case (drop sheet match in the
accessor).  Each bites a distinct surface area, validating the
three-pronged coverage.

## Frictions

### 1.  Wedge walker corners itself against cuts in refuse mode (again).

SA.1.6's first attempt constructed a 2-branch problem in refuse
mode (default), and the walker failed to reach the lower-half grid
point `-1 - 0.5im` from the IC at `z = 1` — same cornering as
PNB.1.2 / PNB.1.4 hit in worklog 042.  Resolution: switched the
test fixture to `cross_branch = true` so the walker can navigate
freely; the test's job is `eval_at_sheet` shape validation, not
walker behaviour under cuts.  Same lesson as worklog 042 §"Friction
1": test the contract, not the route.

### 2.  Solution is `exp(z - z₀)`, not `exp(z)`.

SA.1.4's first attempt asserted `sol_k.grid_u[i] ≈ exp(z)` — but
the IC is `(z₀, u₀, u'₀) = (1, 1, 1)`, so the analytic solution
is `u(z) = exp(z - z₀)` not `exp(z)`.  Caught by mutation-proof
sweep when the GREEN tests bit unexpectedly.  Fix: assert against
`exp(z - z₀)`.  Small lesson: even when the ODE is simple, the IC
shifts the solution and the test oracle must account for it.

### 3.  Tiebreak test fixture was non-degenerate.

SA.1.8's first attempt set up `visited_z = [0, 0.3, -0.3]` and
queried at `(0, -0.3i)` expecting tiebreak — but distance from
query to node 1 is `0.3`, to node 3 is `√(0.3² + 0.3²) = 0.424`.
Not equidistant.  Caught at first run.  Fix: 5-node fixture with
two same-sheet nodes equidistant from `(0, 0)`.  Documented inline
so future readers can see the equidistant geometry without solving
arithmetic in their head.

## Hard-won lesson

**Default kwargs as the backward-compat contract.**  The
`grid_sheet = nothing` default is the entire backward-compat story
— two assertions in SA.1.1 verify that a call with
`grid_sheet = nothing` is byte-identical to a call without the
kwarg.  This is a tighter pin than worklog 042's PNB.1.1
("branch_points = ()" byte-equivalence) because the alternative
dispatch path actually executes (the loop branches on
`grid_sheet === nothing`); equality means the OUTPUTS are
identical, not that the same code ran.  Mutation S2 (the natural
"forget the kwarg" mutation) does NOT bite SA.1.1 — because in
the `grid_sheet = nothing` case BOTH the legacy and "always call
unrestricted" paths produce identical output (they ARE the same
in that case).  S2 is caught by SA.1.3 and SA.1.5 instead.

The lesson: testing "default = legacy" is a weak signal that
nothing regressed; the load-bearing tests are the ones where the
non-default path produces a different output than the legacy
path.

## Empirical results

  | testset | assertions | wall (s) | mutations biting |
  |---------|-----------|----------|------------------|
  | SA.1.1  |  2        | 3.4      | (backward-compat invariant) |
  | SA.1.2  |  3        | 0.3      | (constructive sanity) |
  | SA.1.3  |  2        | 0.0      | S1, S2           |
  | SA.1.4  |  5        | 0.1      | (analytic value pin) |
  | SA.1.5  |  8        | 0.0      | S2, S3           |
  | SA.1.6  |  3        | 0.7      | (input validation) |
  | SA.1.7  |  2        | 0.2      | (input validation) |
  | SA.1.8  |  5        | 0.0      | S1               |
  | **Σ**   | **30**    | **~5**   | S1, S2, S3 all bite |

Full test suite: **2368 → 2398 GREEN** (+30 assertions, 5m 53s
wall on a single thread — +12 s vs A4-baseline).

## Beads

  - `padetaylor-hed` — closed by this commit.
  - To-be-filed next: B5 (FFW Figs 2/3/7 figure scripts) — three
    figures, three beads.  Now unblocked by A4 + A5.

## What is NOT shipped

  - **`PainleveProblem` forwarding** for `branch_points` /
    `cross_branch` / `grid_sheet` — same status as worklog 042
    (trivial; deferred until first B5 figure-script needs it).

  - **Sheet-aware `PathNetworkSolution.(::Complex)` callable
    sugar** — Julia's `(sol::Struct)(args...)` syntax for terse
    per-point eval.  `eval_at_sheet` covers the use case; sugar
    is v2 ergonomics.

  - **B5 figures** (FFW Figs 2, 3, 7) — three separate beads to
    file when starting.  Each needs its own A4 + A5 composition
    test against the figure's parameter set.

  - **Sheet-aware `PoleField.extract_poles`** — currently pulls
    poles from the full `visited_pade` without sheet awareness.
    For multi-sheet figures it makes sense to query
    `extract_poles(sol; sheet = [k])` returning only that sheet's
    poles.  Not yet bead-tracked.

## References

  - **ADR-0013** §"Open follow-ups: A5 sheet-aware Stage-2".
  - **Worklog 042** — A4 substrate that A5 wires into.
  - **`src/PathNetwork.jl`** — the file the new code lands in.
  - **`test/sheet_aware_stage2_test.jl`** SA.1.1–SA.1.8 —
    acceptance tests + mutation-proof footer.
