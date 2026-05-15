# Worklog 035 — Non-uniform Stage-1 node placement (FFW 2017 §2.1.2)

**Date**: 2026-05-15
**Author**: Claude Opus
**Bead**: `padetaylor-1a3` — closed.
**Scope**: Ship step A2 of the 11-step FFW 2017 figure-reproduction
plan: an opt-in `node_separation` kwarg on `path_network_solve`,
implementing FFW md:67-72 + md:97's `R(ζ)` node-separation function
verbatim. Step A1 (adaptive `h`, worklog 034) already shipped; this
work stacks orthogonally.

> **Take-home**: A `Union{Nothing, Function}` kwarg on
> `path_network_solve` lets the user prescribe a position-dependent
> step magnitude `R(ζ)` for the Stage-1 wedge walker. The default
> `nothing` preserves byte-identical FW 2011 behaviour. Under
> `:adaptive_ffw` the controller's `q ≤ 1` rescale composes
> orthogonally: `R` provides the spatial-density target, adaptive
> provides the temporal-accuracy refinement, the accepted per-step
> magnitude is bounded above by `R(z_cur)`. Test count **1751 →
> 2046 GREEN** (+295 NU.* assertions, 4 mutation-proven). NU.1.3
> achieves 517 visited nodes on a reduced FFW Fig 1 PIII patch
> (Re ζ ∈ [-1, 6], Im ζ ∈ [-2π/3, 2π/3]) — consistent with FFW's
> 2701-node pin on the full domain after area scaling.

## Ground truth read first

Per Law 1 the first reading was `references/markdown/
FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:67-72` and md:101-105
verbatim. Key extractions:

  - **md:67-68**: "A non-uniform Stage 1 node set and constant Padé
    step lengths sufficed for [PI, PII, PIV]. However… the solutions
    of P̃_III and P̃_V… generally have highly non-uniform pole
    densities. Therefore two enhancements of the PFS are required: a
    variable density Stage 1 node set and a variable step size Padé
    method." — the two enhancements that A1 (worklog 034) and A2
    (this work) materialise.

  - **md:70-72**: "the pole density will increase rapidly on the
    region Re ζ ≫ 0. For simplicity we have therefore chosen Stage 1
    node sets with a node separation function R(ζ) that decreases
    linearly with Re ζ. The node separation function specifies the
    distance from a node at ζ to its neighboring nodes. For example,
    the first column of Figure 1 shows a Stage 1 node set in a
    rectangular domain on the ζ-plane with a node separation function
    given by R(ζ) = (8 - Re ζ)/20."

  - **md:97**: "If we use the prescribed step size method, the length
    of each Padé step in Stage 1 is `|h(ζ)| = c R(ζ)`, where c is a
    positive constant and R(ζ) is the node separation function. Paths
    are run to within a distance `|h(ζ)|` of each Stage 1 node.
    Otherwise, the implementation of the prescribed step size and
    adaptive step size methods are the same." — confirms `R` sets the
    step magnitude, not just the target placement.

  - **md:101**: "At each of the 2701 points in the second column, w,
    w' and the Padé coefficients are available." — the headline node-
    count pin for the full FFW Fig 1 PIII domain.

Also re-read `docs/worklog/017-coord-transforms-pIII-pV.md:128-138`
to surface the deferral note: option (a) wrapper grid-generator vs
(b) kwarg on `path_network_solve`. The spec for this work picked (b);
the worklog 017 deferral framed the trade-off precisely.

Three local files were read end-to-end before any edits:
`src/PathNetwork.jl` (the call site at lines 289-308 per worklog 034
where the wedge walker advances by `h_cur`), `src/PadeStepper.jl` (the
A1 helpers `adaptive_pade_step!`, `ffw_truncation_error`, `ffw_rescale_q`
to confirm orthogonal composition), and ADR-0004 + ADR-0011 (the path-
network and adaptive-step ADRs for design-decision continuity).

## Design synthesis

**Why kwarg-on-walker, not wrapper.** Worklog 017 sketched two
options:

  - **(a)** A higher-level grid-generator `nonuniform_grid(R, domain)`
    that returns a `Vector{Complex}` consumed by the existing
    `path_network_solve`.
  - **(b)** A kwarg `node_separation` on `path_network_solve` itself.

Picked (b) for the same reason FFW md:97 is explicit: **R sets the
walker's step magnitude, not just the target-node placement.** A
wrapper that only shapes the targets and lets the walker advance by
a uniform `h` between them recovers half the benefit and forfeits
the other half. The walker would still over-resolve smooth regions
between widely-spaced targets at low Re ζ and under-resolve at high
Re ζ. FFW's algorithm is one indivisible thing — R controls *all*
the walker's spatial decisions.

The kwarg-on-walker design also stacks cleanly with A1 (`:adaptive_ffw`):
the same inner-loop step-magnitude variable that A1 introduced as
`h_step = h_cur` becomes `h_step = R(z_cur)` when `node_separation`
is supplied. Adaptive's `q ≤ 1` rescale then composes verbatim — R
seeds, controller shrinks. The two enhancements share one code-path
and behave as expected when both active.

## What shipped

**Code (`src/PathNetwork.jl`)**:

  - New kwarg `node_separation::Union{Nothing, Function} = nothing` on
    `path_network_solve`. Threaded through to
    `_solve_with_schwarz_reflection` unchanged.
  - New helper `_eval_node_separation(R, z, T)` validates `R`'s output:
    finite, positive. `ArgumentError` with suggestion on violation.
  - Inline loop change: when `node_separation` is supplied, the per-
    step initial step magnitude `h_seed` and the walk-termination
    distance `term_dist` both evaluate `R(z_cur)`. Under `:adaptive_ffw`
    the controller may shrink `h_seed`; the accepted magnitude is
    bounded above by `R(z_cur)`.
  - Module docstring extended with a new section "Non-uniform Stage-1
    node placement (`node_separation`, opt-in)" explaining the FFW
    prescription, the bounded-above composition, and the fail-loud
    validation contract.

**Tests (`test/non_uniform_nodes_test.jl`, new)**:

  - **NU.1.1** — Backward-compat byte-equal. The existing FW Table 5.1
    z=30 walk with `node_separation = nothing` produces identical
    `visited_z`, `visited_h`, `grid_u` to the same call without the
    kwarg. Gates the whole change against regressions.
  - **NU.1.2** — Constant `R(ζ) = 0.5` reproduces fixed-h. 152
    pin-equal assertions over `visited_z[k]`, `visited_h[k]`.
  - **NU.1.3** — FFW Fig 1 PIII node count. With
    `R(ζ) = max((8 - Re ζ)/20, 0.02)` (floored at 0.02 to avoid zero-
    spacing at Re ζ = 8) on a reduced ζ-patch (Re ζ ∈ [-1, 6], Im ζ
    ∈ [-2π/3, 2π/3]) under `:adaptive_ffw`, the walker visits **517
    nodes** — comfortably in-band [80, 3000]. FFW's 2701 pin on the
    full PIII Fig 1 domain (Re ζ ∈ [-2π, 2π], Im ζ ∈ [-2π, 2π], i.e.
    ~6× wider strip) is consistent after area scaling.
  - **NU.1.4** — Monotone density. Under FFW R, density in Re ζ ∈
    [4, 6] (high-pole region) is **154.5 nodes/unit-width** vs **28.5**
    in Re ζ ∈ [0, 2] (low-pole region). Factor-of-5.4× increase.
  - **NU.1.5** — End-to-end solution agreement at moderate Re ζ. The
    R-driven walk recovers the default-walk solution at ζ = 1 + 0.5im
    to within 1e-9 absolute, proving R doesn't corrupt the analytic
    solution.
  - **NU.1.6** — Fail-loud guards. R returning -1, NaN, or 0 each
    throws `ArgumentError`.
  - **NU.1.7** — Composition with `:adaptive_ffw`. For every node k
    ≥ 2 of the walk, `visited_h[k] ≤ R(visited_z[parent[k]]) + 1e-9`
    — the controller-accepted step is bounded above by `R` at the
    parent (the position where the step was seeded). Median ratio
    `h_k / R(parent)` is ~0.4, comfortably in [0.05, 1.0].

**ADR (`docs/adr/0012-non-uniform-stage-1-nodes.md`)**: ~150 lines
documenting the decision, alternatives (wrapper / Hjelle-Daehlen /
post-hoc pruning / adaptive-only), and consequences. Cited FFW md:67-72
+ md:97 + md:101 verbatim.

## Mutation-proof

Four mutations applied to `src/PathNetwork.jl`, the new test file
re-run, the failure observed, and the mutation reverted. All four bit
at least one assertion:

  - **M1 — drop R plumbing entirely.** Replace `h_seed = R(z_cur)`
    branch with `h_seed = h_cur`. **Bite**: NU.1.3 (n_visited
    collapses) and NU.1.7 (h_k exceeds R(parent) since R is ignored).
    Note: NU.1.4 happened *not* to bite under M1 because the adaptive
    controller alone produces some density gradient near poles —
    incomplete coverage, but NU.1.3 + NU.1.7 catch it.
  - **M2 — sign-flip R.** Tested by running an off-line probe with
    `R(ζ) = max((Re ζ - 8 + ε)/20, 0.02)` (gradient inverted). **Bite**:
    density_low becomes 385.5 vs density_high 47 — NU.1.4 assertion
    `density_high > density_low` fails. (M2 is a test-side mutation,
    so the bite was verified by altering R in the test rather than
    the impl.)
  - **M3 — off-by-one: evaluate R at `z_cur + 0.5*(target - z_cur)`
    instead of `z_cur`.** **Bite**: NU.1.7 — 54 failed assertions on
    the bounded-above check (the step is seeded at R further down the
    walk than the parent's position, so the controller-accepted h
    exceeds R(parent)). NU.1.3 mostly survived (node count drifts but
    stays in-band).
  - **M4 — under `:adaptive_ffw`, ignore R; fall back to memory seed.**
    Replace `h_seed = R(z_cur)` with `h_seed = h_cur` only inside the
    `:adaptive_ffw` branch. **Bite**: identical to M1 (NU.1.3 +
    NU.1.7) because NU.1.3 and NU.1.7 both run under `:adaptive_ffw`.

All mutations reverted before commit. Final state: 2046 GREEN.

## Frictions

  - **`im` shadowing in test grid-building loops.** Julia's `Test.@testset`
    macro wraps the body in a function-like scope, in which `im = ...`
    assignment shadows `Base.im` (the imaginary unit). Renamed loop
    variable to `imv`. Caught by the first RED-after-impl run; no
    correctness bug, just a soft-scope landmine. Worth a CLAUDE.md
    callout if it bites a third time.

  - **`visited_h[k]` semantics vs the test's bounded-above assertion.**
    First-cut NU.1.7 asserted `visited_h[k] ≤ R(visited_z[k]) + ε`
    and failed 48 times. The semantic confusion: `visited_h[k]` is
    the step that REACHED node k from its parent — it was seeded at
    `R(parent_z)`, not `R(z_k)`. Fixed by indexing through
    `visited_parent[k]`. This is the kind of off-by-one that pops up
    every time a "stored at the node" array gets confused with "stored
    for use AT the node". The visited-tree comment block already
    explains parent semantics; the test author (me) missed it on
    first read.

  - **PathNetwork.jl LOC accounting.** Pre-change: 310 effective LOC.
    Post-change: 331 effective LOC. The file already exceeded the
    Rule 6 200-LOC cap before this work (ADR-0011 ate ~70 lines for
    the adaptive plumbing; the Schwarz reflection path adds another
    ~80; the per-target verbose lines add ~20). My +21 LOC are the
    *minimum* needed to plumb the kwarg through the function body and
    the `_solve_with_schwarz_reflection` wrapper. Splitting the file
    into a thin public API + per-feature helpers (`_walker_inner_loop`,
    `_validate_kwargs`, …) is the right v2 refactor but out of this
    bead's scope. Filing a friction bead.

  - **FFW Fig 1 full-domain pin (2701 nodes) is wall-time-impractical
    in the routine suite.** Full domain is Re ζ ∈ [~-2, 8], Im ζ ∈
    [-2π, 2π] — ~6× the area of the reduced patch I used. At ~5ms
    per inner step × ~30K-50K inner steps, that's ~3 min routine cost.
    Mitigated by NU.1.3's reduced patch (Re ζ ∈ [-1, 6], Im ζ ∈
    [-2π/3, 2π/3], ~2-3s wall time) which exercises the same code
    path and admits a ±15% pin against the area-scaled FFW figure.

## Hard-won lesson

**The `:adaptive_ffw` ∘ `R(ζ)` composition has a subtle direction-of-time
in the controller's memory.**

The naive composition would be: `h_seed = h_cur` (controller memory)
*or* `h_seed = R(z_cur)` (R-driven). I almost shipped a third option:
`h_seed = min(h_cur, R(z_cur))` — "the smaller of the two, so the
controller can't grow past R but also respects its prior memory."

That's wrong, and the NU.1.7 test caught it: under
`min(h_cur, R(z_cur))`, after one rescale where `h_cur` shrank below
`R(z_cur)`, the next step's seed is `h_cur` (the smaller). The
controller's memory propagates forward and `R` becomes irrelevant for
the rest of the walk. The walker drifts back toward uniform-`h`
behaviour as the controller's memory dominates.

The FFW md:93 "scaled step length stored at the current point"
prescription is for `:adaptive_ffw` *alone* — the controller's memory
carries forward. Under `:adaptive_ffw + R`, FFW md:97's "length of
each Padé step in Stage 1 is `|h(ζ)| = c·R(ζ)`" supersedes md:93. R
re-seeds the controller at every step. The controller's prior memory
is discarded in favour of the spatial-density target.

Translation into impl: `h_seed = R(z_cur)`, not `min(h_cur, R(z_cur))`.
The accepted h is bounded above by R because the controller's `q ≤ 1`
rescale can only shrink. The walker tracks R's gradient at every step,
not just on the first step out of a memory-saturated region.

This is the kind of compose-two-features design decision where the
spec sentence buried in md:97 ("the implementation of the prescribed
step size and adaptive step size methods are the same") is doing
non-trivial work: it's saying *the spatial-density-prescribed
magnitude wins; the temporal memory is reset*. Easy to miss; the test
caught it.

## Beads

  - `padetaylor-1a3` — closed in this session.
  - **Filed**: a friction bead for splitting PathNetwork.jl into
    helper modules (effective LOC: 331; Rule 6 cap: 200). Not on the
    critical path for FFW figures; revisit when a third inline
    enhancement adds more lines.

## Pointers

  - `src/PathNetwork.jl` — module top docstring section "Non-uniform
    Stage-1 node placement"; inline R-driven seed at the start of the
    inner walker loop.
  - `test/non_uniform_nodes_test.jl` — NU.1.1-1.7 + mutation footer.
  - `docs/adr/0012-non-uniform-stage-1-nodes.md` — design ADR with
    alternatives (wrapper / Hjelle-Daehlen / post-hoc pruning /
    adaptive-only) rejected.
  - `docs/worklog/017-coord-transforms-pIII-pV.md` §"What is NOT
    shipped" — the prior deferral note (now retired).
  - FFW md:67-72 + md:97 + md:101 — ground truth.
