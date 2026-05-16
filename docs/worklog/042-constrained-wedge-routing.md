# Worklog 042 — Constrained-wedge routing + sheet bookkeeping (bead `padetaylor-bho` closed)

**Date**: 2026-05-16
**Author**: Claude Opus 4.7 (1M context)
**Bead**: `padetaylor-bho` — closed by this commit.
**ADR**: 0013 (Accepted in this commit).
**Scope**: Step A4 of the 11-step FFW 2017 figure-reproduction plan
— the **walker-side enforcement layer** for Riemann-sheet bookkeeping
during `PathNetwork.path_network_solve`.  Adds three opt-in kwargs
(`branch_points`, `branch_cut_angles`, `cross_branch`) and one new
internal module (`src/BranchTracker.jl`, ~150 LOC).  Together with
A3 (η-plane PVI RHS, worklog 041) this unblocks the substrate for
FFW Fig 2 column 2/3 and Fig 3.

> **Take-home**: A4 ships **both** the refuse mode (walker rejects
> any step crossing a branch cut) and the cross mode (walker is
> permitted to cross; per-branch sheet counter bumps accordingly)
> per the user's "A4a+A4b together" scope decision.  The geometric
> primitives live in a new `src/BranchTracker.jl` (~150 LOC) so
> PathNetwork.jl's growth is held to ~80 LOC of kwarg threading,
> per ADR-0013's split-vs-inline analysis.  Test count **2262 →
> 2368 GREEN** (+106 assertions across two new files,
> `test/branch_tracker_test.jl` and `test/path_network_branch_test.jl`);
> six mutations (3 in BranchTracker, 3 at the PathNetwork wiring
> level) each bite a distinct subset of the new tests.  Default
> kwargs (`branch_points = ()`) are byte-equivalent to the pre-A4
> walker — the load-bearing backward-compat invariant verified by
> PNB.1.1's 8 assertions.

## Ground truth read first

Per Law 1, foundational reads were FFW 2017 §2.2.2 + the relevant
SheetTracker primitives:

  - **FFW md:163-189** — the circumambulation algorithm.  md:178
    is the verbatim spec: "none of the paths overstep the branch
    cuts, indicated by solid lines, on a given sheet" (refuse mode);
    "the paths move through the branch cut by running only in a
    counterclockwise direction around ζ = 0" (cross mode).
    md:185-189 sheet parametrisation via `(θ_k, φ_k)` index ranges
    — the motivation for the per-branch TUPLE sheet representation
    rather than a fused scalar.

  - **`src/SheetTracker.jl`** (worklog 018) — the v1 primitives
    `winding_delta`, `accumulate_winding`, `sheet_index` that
    compute sheet indices AFTER a path-network walk has completed.
    These complement BranchTracker (walker-side enforcement) cleanly
    — same `(-π, π]`-normalised `winding_delta` reused inside
    `step_sheet_update`; no duplication.

  - **`src/PathNetwork.jl`** — the wedge-loop and visited-tree
    accumulation (now 818 LOC after the A4 additions; bead
    `padetaylor-qum` open for the orthogonal LOC-cap split).

  - **ADR-0013** (`docs/adr/0013-constrained-wedge-and-sheet-
    bookkeeping.md`) — written first in this session per Law 2,
    documents the half-line-cut decision, per-branch-tuple sheet
    choice, helper-module split rationale, and the v1 corner about
    Schwarz-reflection incompatibility.

## Design synthesis

### Layer 1: BranchTracker module (new, ~150 LOC)

Four public functions:

  - `segment_crosses_cut(z_cur, z_new, branch, cut_angle) -> Bool`
    — 2×2 parametric segment-ray intersection via Cramer's rule on
    complex imaginary parts.  Closed conventions: `0 < t < 1` and
    `s > 0` (open both ends on segment; closed on ray-direction
    `s = 0` excludes the branch point itself).  Parallel segments
    return `false`.

  - `any_cut_crossed(z_cur, z_new, branches, cut_angles) -> Bool`
    — short-circuit OR over `segment_crosses_cut` for multi-branch
    inputs.

  - `step_sheet_update(sheet_old, z_cur, z_new, branches, cut_angles)
    -> Vector{Int}` — for cross mode.  Per branch, if the segment
    crosses the cut, bump the per-branch counter by `sign(Δθ_k)`
    where `Δθ_k = winding_delta(z_cur, z_new, branch_k)` (reusing
    the SheetTracker primitive).  Returns a fresh `Vector{Int}` of
    the same length as `sheet_old`.

  - `resolve_cut_angles(branches, angles) -> NTuple{N, Float64}`
    — kwarg-resolution helper.  Scalar broadcasts to the right
    length; tuple/vector of matching length passes through;
    mismatched length throws (CLAUDE.md Rule 1).

### Layer 2: PathNetwork wiring (~80 LOC delta)

`path_network_solve` gains three new opt-in kwargs:

```julia
path_network_solve(prob, grid;
    branch_points     :: Tuple{Vararg{Complex}}        = (),
    branch_cut_angles :: Union{Real, Tuple{Vararg{Real}}} = π,
    cross_branch      :: Bool                            = false,
    initial_sheet     :: AbstractVector{<:Integer}      =
                            zeros(Int, length(branch_points)),
    # ... existing kwargs unchanged ...
)
```

Validation block (added inside the function body, before the main
loop) checks length consistency between `initial_sheet` and
`branch_points`, the orthogonality with `enforce_real_axis_symmetry`,
and the no-`cross_branch`-without-branches invariant.

The wedge inner loop gains two new branch-aware operations:

  1. **Filter forbidden candidates** (refuse mode only): a helper
     `_filter_forbidden_candidates(evals, z_cur, branches, cut_angles, CT)`
     replaces each forbidden candidate with the existing
     `(z_cur, CT(Inf, 0), CT(0, 0), nothing)` failure sentinel.  The
     downstream `_select_candidate` (both `:min_u` and
     `:steepest_descent`) naturally avoids `nothing`-Padé entries
     via the existing "all 5 wedge candidates failed" guard, which
     gains a contextual message about `cross_branch = true` when
     `branched` is true.

  2. **Update sheet** (after accepted step): for cross mode, call
     `step_sheet_update(visited_sheet[parent_idx], z_cur, z_new,
     branches, cut_angles)`; for refuse mode, copy the parent's
     sheet unchanged; for `branch_points = ()`, push `Int[]`.

`PathNetworkSolution` gains a 10th field `visited_sheet ::
Vector{Vector{Int}}`; a 9-arg backward-compat constructor preserves
existing call sites (three external test fixtures construct it
directly; their fixtures keep working unchanged).

### Layer 3: Module include order

`SheetTracker` was moved up in `src/PadeTaylor.jl`'s include order
(from after CoordTransforms to before PathNetwork) so BranchTracker
can `using ..SheetTracker: winding_delta` and PathNetwork can in
turn `using ..BranchTracker`.  No reverse-dependency conflicts —
SheetTracker is freestanding (no `using ..PathNetwork`).

## What shipped

  - **`docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md`**
    (~190 lines): the Decision, Why-these-choices, Consequences,
    Open follow-ups (5 v2 items: finite cut segments, A5/B5/PainleveProblem
    forwarding, PathNetwork LOC-cap split).

  - **`src/BranchTracker.jl`** (new, ~158 LOC): four public functions
    + module docstring chapter explaining the cut model, the
    intersection algorithm, the per-step sheet update law, and the
    caller contract.

  - **`src/PathNetwork.jl`** (+~80 LOC): new kwargs, validation,
    filter helper, wedge-loop wiring, struct field + 9-arg
    backward-compat constructor.  Docstring grows four kwarg
    paragraphs documenting the new API surface.

  - **`src/PadeTaylor.jl`** (+3 lines): include-order reshuffle
    (`SheetTracker` before `PathNetwork`); `BranchTracker` add;
    architecture-comment update.

  - **`test/branch_tracker_test.jl`** (new, ~180 LOC): seven testsets
    BT.1.1–BT.1.7, 30 assertions, three mutations (B1/B2/B3)
    documented in the footer with their exact bite counts.

  - **`test/path_network_branch_test.jl`** (new, ~225 LOC): seven
    testsets PNB.1.1–PNB.1.7, 75 assertions, three mutations
    (PB1/PB2/PB3) documented in the footer.

  - **`test/runtests.jl`**: include both new test files; add
    `BranchTracker` to the umbrella-loads assertion (+1).

  - **This worklog** (`docs/worklog/042-constrained-wedge-routing.md`).

## Mutation-proof results

All six mutations applied + restored per CLAUDE.md Rule 4.  Each
bites a distinct subset of the new tests; no mutation is "silent."

  | Mutation | Bite count (of 105 new assertions)         | Where it hits           |
  |----------|--------------------------------------------|-------------------------|
  | B1       | 5 RED                                      | BT.1.1 (3), BT.1.3 (1), BT.1.5 (1) |
  | B2       | 9 RED                                      | BT.1.5 (2), PNB.1.3 (7) |
  | B3       | 52 RED — broadest cascade                  | BT.1.5 (1), BT.1.6 (1), PNB.1.3 (50) |
  | PB1      | 2 RED                                      | PNB.1.5 only            |
  | PB2      | 8 RED                                      | PNB.1.3 only            |
  | PB3      | 1 RED                                      | PNB.1.6 wrong-length    |

B3's 52-RED cascade is the deepest signal — the
`segment_crosses_cut` gate inside `step_sheet_update` is load-bearing
on every single step's sheet update, so dropping it produces a
divergence cascade through the entire chain-walk verification in
PNB.1.3.  The other five mutations have tighter bites (2-9
assertions); each is documented with its specific failure mode in
the test file footers (`test/branch_tracker_test.jl` lines 161-185,
`test/path_network_branch_test.jl` lines 188-217).

## Frictions

### 1.  Wedge walker corners itself against cuts under `:min_u`.

The first cut of PNB.1.2 used `f(z,u,up) = up` (effectively `u =
exp(z)`) with grid `[1+1i, -1+0.5i, -1-0.5i]` and a branch at origin
with arg=π cut.  The `:min_u` selection drives the walker LEFT (smaller
`|exp(z)|`), and with goal at `-1+0.5i` the walker's wedge candidates
include `lower-left` directions.  At some intermediate step the
walker landed just below the cut (`z ≈ -0.99 - 0.015i`) and then
all 5 wedge candidates aimed up — all forbidden by refuse mode → fail.

Two responses considered:

  (a) Use `:steepest_descent` selection — but this also has
      degenerate corners and would just move the symptom.
  (b) Choose a test geometry the walker can navigate without
      cornering itself.

Per Rule 9 (senior-engineer-grade only), I went with (b): the test's
job is to prove the *filter* is load-bearing, not to demonstrate the
walker is omniscient.  Grid changed to all-upper-right
(`[1+1i, 0.5+0.5i, 0.5+1i]`) where the walker naturally stays above
the cut.  A separate direct unit test of the `_filter_forbidden_candidates`
helper proves the filter contract independent of walker behaviour.

This is the lesson worklog 040 named (FFW Fig 5 §"Friction 1"): "test
the contract, not the route the walker chose to honor it."

### 2.  Parallel-segment degeneracy in `segment_crosses_cut`.

My first PNB.1.5 attempt used IC = z=0, target = z=2, branch at z=1,
cut along arg = 0.  Walker steps along the real axis — but so does the
cut, so `det(d, u) = imag(d * conj(u)) = 0`, parallel, predicate
returns `false`, walker happily walks through the cut without
detecting.

This is *correct* behaviour per BT.1.2: parallel = degenerate = not
detected (FW md doesn't specify, and a tangent crossing is
geometrically ambiguous).  The test was just configured to exploit the
degeneracy unintentionally.  Resolution: redesign PNB.1.5 with a
cut perpendicular to the walker's direction
(IC=(-0.5,0), target=(-0.5,1), branch=(0,0.3), cut arg=π) where every
step direction unambiguously crosses.  All 5 candidates forbidden →
fail-loud message containing "cross_branch".

### 3.  Mutation PB1 is caught by PNB.1.5 only.

The natural test of "is `_filter_forbidden_candidates` actually
called?" is "in refuse mode on a geometry that requires the filter,
does the walker behave correctly?"  But the geometry where the
filter is OBSERVABLY load-bearing (i.e., the walker would otherwise
cross) is hard to construct without the corner-trap of Friction 1.

PNB.1.5 catches PB1 indirectly: with the filter disabled, the walker
finds a valid candidate (the cut-crossing one) and DOES NOT throw.
The expected `ErrorException` doesn't fire, `@test_throws` fails.

The direct unit test of the helper in PNB.1.2 does NOT bite PB1
(because the helper itself is unchanged — only its caller is
commented out).  This is the principled coverage split: helper tests
prove correctness; PNB.1.5 proves wiring.  Two distinct
responsibilities, two distinct tests.

### 4.  PathNetwork.jl now at 818 LOC (Rule 6 cap = 200).

This bead's additions move PathNetwork.jl from 744 → 818 LOC against
the Rule-6 cap of 200.  Bead `padetaylor-qum` is open for the
LOC-cap split (it was open BEFORE A4; A4 just makes it more
pressing).  ADR-0013 §"PathNetwork.jl LOC-cap split" calls this out
explicitly.  The split itself is orthogonal scope — A4 ships with
the existing module shape and the bead-tracked refactor follows.

## Hard-won lesson

**The "narrow but load-bearing" mutation hint.**  PB3 (drop
`initial_sheet` length check) bit only 1 RED out of 75 new
assertions.  A surface read of the bite count might suggest the
check is over-narrow ("only one test catches it — is it really
load-bearing?").  But the failure mode is *systemic*: without the
check, the wrong-length input cascades into a downstream `@assert`
inside `step_sheet_update` at solve time — a much harder-to-debug
error message ("AssertionError on line 95") replacing the
proximate-and-clear "initial_sheet length 1 must equal branch_points
length 2".  The single biting assertion in PNB.1.6 explicitly tests
that the *thrown exception type* is `ArgumentError`, not just that
*something* threw.  This is the difference between "the check exists"
and "the check fires at the right level with the right message" —
Rule 1's fail-fast-fail-loud is about the *quality* of the failure,
not just its existence.

**The single-assertion bite is the load-bearing one.**  Wide
mutation cascades (B3's 52 RED) reveal systemic dependencies; narrow
bites (PB3's 1 RED) reveal that a check has exactly one tight
responsibility.  Both are good.  A mutation that bites 0
assertions is the failure case — neither a narrow nor a wide bite.

## Empirical results

  | testset    | assertions | wall (s) | mutations biting |
  |------------|-----------|----------|------------------|
  | BT.1.1     |  7        | 0.0      | B1               |
  | BT.1.2     |  3        | 0.0      | (degenerate pin) |
  | BT.1.3     |  3        | 0.0      | B1               |
  | BT.1.4     |  4        | 0.0      | (multi-branch sanity) |
  | BT.1.5     |  4        | 0.1      | B1, B2, B3       |
  | BT.1.6     |  2        | 0.0      | B3               |
  | BT.1.7     |  7        | 0.1      | (kwarg validation) |
  | PNB.1.1    |  8        | 3.1      | (backward-compat invariant) |
  | PNB.1.2    |  8        | 0.4      | (filter helper unit + tree-edge check) |
  | PNB.1.3    | 51        | 0.1      | B2, B3, PB2      |
  | PNB.1.4    |  1        | 0.0      | (initial_sheet seed) |
  | PNB.1.5    |  3        | 0.2      | PB1              |
  | PNB.1.6    |  3        | 0.3      | PB3              |
  | PNB.1.7    |  1        | 0.3      | (Schwarz incompat) |
  | **Σ**      | **105**   | **~5**   | B1, B2, B3, PB1, PB2, PB3 all bite |

(Plus 1 umbrella-loads assertion for `BranchTracker`.)

Full test suite: **2262 → 2368 GREEN** (+106 assertions, 5m 41s
wall on a single thread — +9 s vs A3-baseline).

## Beads

  - `padetaylor-bho` — closed by this commit.
  - `padetaylor-qum` — open; A4's additions push PathNetwork.jl
    further over the LOC cap.  ADR-0013 §"PathNetwork.jl LOC-cap
    split" calls this out.
  - To-be-filed at next session start:
    - **A5** sheet-aware Stage-2 (no bead yet) — restrict barycentric
      lookup pool to matching `visited_sheet`; necessary before B5.
    - **B5** FFW Figs 2/3/7 (no beads yet) — gated on A5.

## What is NOT shipped

  - **`PainleveProblem` forwarding** for `branch_points` /
    `cross_branch` (trivial; deferred per ADR-0013 §"Open follow-ups"
    until B5 figure-scripts actually need it).

  - **Finite cut segments** (a `branch_cut_segments::Vector{Tuple{Complex,Complex}}`
    kwarg).  Per ADR-0013 §"Why these specific choices", the half-line
    representation covers every FFW 2017 figure that needs cuts; the
    finite-segment lift is v2 ergonomics, deferred.

  - **A5 sheet-aware Stage-2**.  Separate bead.  A4 populates
    `visited_sheet`; A5 makes the Stage-2 nearest-visited lookup
    *respect* the sheet (currently it picks by Euclidean distance only,
    ignoring sheet).  B5 figures need A5.

  - **B5 figure rendering** (FFW Figs 2, 3, 7).  Separate beads,
    gated on A5.

  - **Schwarz × branch composition**.  Currently fails-loud if
    `enforce_real_axis_symmetry = true` is combined with
    non-empty `branch_points` — the Schwarz mirror's
    upper-half-plane-only walk assumes simply-connected geometry,
    which conflicts with cuts.  A v2 lift would walk only the
    upper-half cut-respecting subset, but no caller has asked.

## References

  - **FFW 2017 §2.2.2** — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:163-189` (the circumambulation algorithm); md:178 (verbatim refuse + cross specification); md:185-189 (sheet parametrisation via `(θ_k, φ_k)`).
  - **ADR-0013** — `docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md`.
  - **Worklog 018** — `docs/worklog/018-sheet-tracker-pVI.md` §"What is NOT shipped" — the original deferral that this bead closes.
  - **Worklog 041 + A3** — η-plane PVI RHS; the immediately preceding A4 substrate.
  - **`src/BranchTracker.jl`** — the new helper module.
  - **`src/PathNetwork.jl`** — the file where the new kwargs land.
  - **`test/branch_tracker_test.jl` BT.1.1–BT.1.7** — primitive
    acceptance tests + mutation-proof footer.
  - **`test/path_network_branch_test.jl` PNB.1.1–PNB.1.7** —
    integration acceptance + mutation-proof footer.
  - **HANDOFF.md** §"Session 2026-05-15" §"What is NOT shipped" —
    A4 as the largest remaining infra piece for the FFW arc.
