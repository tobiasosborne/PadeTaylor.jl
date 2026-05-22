"""
    PadeTaylor.VectorPathNetwork

Minimal **vector** Stage-1 path-network walk: the substrate the v0.2
`A_4⁽¹⁾` and `P_I⁽²⁾` pole-field scatter figures (beads V8a / V8b)
stand on.  This module is the vector analogue of v0.1's
`PathNetwork.path_network_solve` — but deliberately *only* its Stage-1
tree construction, not the whole 1010-LOC v0.1 driver.

## Why this module is deliberately minimal

v0.1's `src/PathNetwork.jl` is large because it carries far more than
the Stage-1 walk: Stage-2 fine-grid barycentric extrapolation,
adaptive step-size control (`:adaptive_ffw`), non-uniform node
placement, Schwarz-reflection symmetry, branch-cut / Riemann-sheet
tracking, and a loop-closure diagnostics layer.  V7's mandate (v0.2
plan row V7) is the *minimal* tool the pole-field figures need: a
Stage-1-style wedge walk plus shared-`Q` pole extraction.  Everything
else in v0.1's `PathNetwork` is **explicitly out of scope** here and
registered as a follow-up bead, each with the forcing condition that
would pull it back in:

  - **Stage-2 fine-grid extrapolation** — bead `padetaylor-0ln.20`.
    Forced when a vector figure needs a *dense filled heatmap* of
    `‖y(z)‖` (not just a pole scatter); V8a/V8b only need the scatter.
  - **Edge detection / pole-field mask** — bead `padetaylor-0ln.21`.
    Forced when a vector figure needs the FW §3.2.2 pole-field
    boundary, not the pole locations themselves.
  - **Branch-cut / Riemann-sheet tracking** — bead `padetaylor-0ln.22`.
    Forced when a vector target has genuine branch points (the
    Noumi–Yamada `A_n⁽¹⁾` systems and `P_I⁽²⁾` in companion form are
    meromorphic — single-sheeted — so V7 needs no sheet machinery).
  - **Sheet-aware evaluation** — folded into `padetaylor-0ln.22`.

## Stage-1 — the wedge walk (mirrors v0.1 `PathNetwork`)

Read v0.1's `src/PathNetwork.jl` module docstring, "Stage 1 — path-tree
construction": initialise `visited = {z_0}` at the problem's IC; for
each `target` in `targets`, find the nearest already-visited node
`z_v`, and walk from `z_v` toward `target` via a 5-direction wedge
until within `h` of `target`.  Each landed point becomes a new visited
node.

## The wedge-selection criterion — principled `:max_q_root` (B2)

v0.1 picks the wedge direction that minimises `|u|` — a heuristic for
"step away from the next pole", since a scalar meromorphic solution's
`|u|` grows toward a pole.  V7 lifted that to **min-‖y‖**: a vector
meromorphic solution's components all blow up together at a shared
movable pole (ADR-0019), so `‖y‖` is a *proxy* for "distance to the
next pole".

But ‖y‖ is only a proxy, and ADR-0025's Phase-A4 baseline measured the
proxy failing — the shipped min-‖y‖ V8b wedge walk has ~70 %
catastrophic loop closures.  Bead `padetaylor-0ln.23` proposed the
principled refinement; B2 (`padetaylor-0ln.37.6`) discharges it: the
default wedge criterion is now **`:max_q_root`** — choose the candidate
whose landed node has the **largest distance to its nearest shared-`Q`
denominator root**, the most pole-free disc of validity.  This is exact,
not a proxy: the shared-`Q` roots ARE the system's local poles.  The
machinery lives in the sibling module `VectorWedgeStep` (split out for
the CLAUDE.md Rule 6 LOC cap; read its docstring for the full
rationale).  `:min_y` is kept as an opt-in — the verified V7 behaviour
— but `:max_q_root` is the default for vector walks.

## Adaptive step size — `adaptive = true` (B2)

The v1 walk used a single fixed `h`.  ADR-0025's Phase-A2 probe
(`external/probes/wedge-tractability/REPORT.md`) measured the fixed-`h`
walk **blocking past `|x| ≈ 8`**: a bridging walk between targets
crosses the pole field on a chord and a fixed-`h` chord step lands on a
pole, degenerating the shared-`Q` linear solve.  B2 makes `h` adaptive
(the default): each step is capped by the parent node's shared-`Q`
pole distance (reusing `StepControl.step_pade_root` per ADR-0021) so a
step never overshoots a pole, and `h` grows toward a maximum where the
pole field is sparse and shrinks where it is dense.  The controller is
`VectorWedgeStep._adaptive_h`; pass `adaptive = false` for the verified
V7 fixed-`h` behaviour.  The `h` kwarg is then the *maximum* step
(`h_max`); the walk floors `h` at `H_MIN_RATIO·h` (`h_min`) and throws
if the pole field forces a step below that (Rule 1 — an honest
wedged-walk signal).  C3 of ADR-0025 — the hand-tuned figure `h = 0.1`
— is retired by this adaptivity.

## The per-node *canonical* shared-`Q` store — why pole extraction is *cleaner*

Each visited node stores a shared-denominator approximant produced by
`VectorStepper.vector_pade_step_with_pade!`: `d` numerator polynomials
`numerators` over **one** denominator `denominator`.  This is the
keystone difference from v0.1's `PoleField`: there a scalar node
stores one `(m,n)` Padé and its single denominator; a *vector* node
under naïve per-component fitting would store `d` denominators — `d`
inconsistent pole sets agreeing only up to noise (ADR-0019).  The
shared `Q` collapses that to **one** polynomial per node whose roots
are *every* component's poles simultaneously.  `VectorPoleField`
exploits exactly this.

The stored approximant is the **canonical** one: centred *at the node
itself* and built with the *real* step magnitude `h`, in the rescaled
variable `t = (z - z_node)/h`.  This mirrors v0.1's
`PathNetwork._local_pade` design (`src/PathNetwork.jl`, "Per-node
canonical Padé"): the wedge-direction approximant a step produces is
centred at the *parent* and lives in a *complex*-`h` variable — its
roots would map to the z-plane via the wrong centre and the wrong
scale.  So after the wedge selects the next position, we take one
extra `vector_pade_step_with_pade!` from `(z_new, y_new)` with the
real step `h` to build the node's own canonical store.  `VectorPoleField`'s
`z = z_node + t·h` mapping is then exact.

## Target ordering — the FW/FFW double-run accuracy indicator (VC-10)

FW 2011 §3.1 (line 156) shuffles the Stage-1 targets into a *random*
order before walking; FFW 2017 (`references/markdown/FFW2017_painleve_
riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:
246-247`) makes that randomness a *diagnostic*: "the PFS method selects
the target nodes in Stage 1 in a random order ... if we compute the
same solution twice, different paths will be run, resulting in
solutions that should differ by approximately the numerical error."
The disagreement between two differently-ordered runs is a practical,
oracle-free accuracy estimate of the result — VC-10 of ADR-0025.

The walk's target processing order is therefore made **controllable**
by the optional `rng` kwarg.  With `rng === nothing` (the default) the
targets are processed in the *given* order — the V7/V8 behaviour,
reproduced bit-identically (ADR-0001 additive guarantee).  With `rng`
an `AbstractRNG`, the targets are `shuffle`d by it before the walk: a
*different* processing order builds a *different* path-network tree, so
the visited node sequence — and any pole field independently extracted
from it — differs by the walk's numerical error.  The walk for a
**fixed** `rng` is still fully deterministic (the same `rng` state ⇒
the same shuffle ⇒ the same tree, bit-identical); VC-10 is the
*complementary* fact that *different* orderings converge to the *same*
pole field within the reported accuracy.  Two seeded runs + a
nearest-neighbour pole-field match (the figure helper `surf_vc10_two_
run`) is the FW/FFW double-run indicator.

## The Stage-1 tree

`VectorPathNetworkSolution` carries the visited tree.
`visited_parent[k]` is the index of the node `visited_z[k]` was
reached from; the root (IC) node has parent `0`.  The edge set
`{(visited_parent[k], k) : k ≥ 2}` is the Stage-1 path tree — the
first node landed in a per-target walk chains off the nearest
already-visited node, later nodes of the same walk chain off their
immediate predecessor (the v0.1 `PathNetwork` convention verbatim).

## Stage-2 — the fine-grid fill (bead `padetaylor-0ln.20`)

V7 shipped only the Stage-1 tree, because the `A_4⁽¹⁾`/`P_I⁽²⁾`
*scatter* figures need only the pole locations — the roots of each
node's shared `Q`.  The headline `P_I⁽²⁾` tritronquée *surface* figure
needs more: a dense filled heatmap of `y(z)` over the pole-rich wedge.
That is Stage 2, the direct vector lift of v0.1's
`PathNetwork.path_network_solve` Stage-2 step (`src/PathNetwork.jl:29-32`,
"Stage 2 — fine-grid extrapolation"; the fill loop at `:669-693`).

When the caller passes a `fine_grid`, `vector_path_network_solve`, after
building the Stage-1 tree, evaluates the dense solution at every grid
point and stores the result in the `grid_z`/`grid_y` solution fields.
The fill mechanics — nearest-visited lookup, the `t = (z_f − z_v)/h_v`
rescaling with the real canonical step, the shared-`Q` `Pᵢ(t)/Q(t)`
evaluation, and the `extrapolate` disc-radius gate (ADR-0015) — live in
the sibling module `VectorPathNetworkStage2`, split out for the
CLAUDE.md Rule 6 LOC cap; read its docstring for the full procedure.
When no `fine_grid` is passed both Stage-2 fields are empty and the V7
scatter-figure call sites are byte-unaffected.

## Fail-fast contract (CLAUDE.md Rule 1)

Empty `targets`, non-positive `h`, a `wedge_angles` not of length 5,
and an `on_target_failure` outside `(:throw, :skip)` throw an
`ArgumentError` with a `Suggestion`.  An unreachable target, an
all-five-wedge-candidates-fail, and an adaptive-step collapse throw a
typed `VectorWalkError` carrying a `reason::Symbol` — and, under the
default `on_target_failure = :throw`, abort the whole solve.
Numerical breakdowns inside `vector_pade_step_with_pade!` (a step
landing on a pole — `Q(1) ≈ 0`) propagate unchanged; the wedge
selection is *meant* to steer around them, and a thrown `DomainError`
is the honest fail-loud signal that it could not.

## The resilient walk — `on_target_failure = :skip` (ADR-0026 D1)

The default `on_target_failure = :throw` aborts the whole run the
moment one target's walk fails.  For a dense space-filling target set
(the headline-figure wedge) that is too brittle: one unreachable
target costs every other target's coverage.  With
`on_target_failure = :skip` the per-target walk is wrapped so a
`VectorWalkError` is *caught*, recorded as a `VectorWalkFailure` in the
solution's tenth field `failed_targets`, and the walk moves on to the
next target.  This is fail-*loud-by-accounting*, never a silent
swallow (Rule 1): `failed_targets` is first-class data the caller must
surface.  An exception that is *not* a `VectorWalkError` — an
unrecognised bug — is re-thrown even under `:skip`; the walk never
skips an error it cannot classify.  Partial nodes a skipped walk laid
down before it failed are kept (the walk pushes to `visited_*` only on
an accepted step).  See `VectorWalkFailure` and the
`on_target_failure` kwarg doc.

## References

  - `src/PathNetwork.jl` — the v0.1 scalar path-network whose Stage-1
    wedge walk this module mirrors (`:21-27`, "Stage 1 — path-tree
    construction"; `:153`, the `DEFAULT_WEDGE`).
  - `src/VectorStepper.jl` — `vector_pade_step_with_pade!`, the
    per-step shared-`Q` primitive this walk loops.
  - `src/VectorPoleField.jl` — the shared-`Q` pole extractor that
    consumes `VectorPathNetworkSolution`.
  - `src/VectorProblems.jl` — the `_eval_poly` Horner evaluator and the
    `Pᵢ(t)/Q(t)` shared-`Q` dense-evaluation pattern Stage 2 mirrors.
  - `docs/v0p2_plan.md` — §"The 10 phases", V7 row.
  - `docs/adr/0019-shared-denominator-pade.md` — the shared `Q` whose
    roots are the pole field.
  - `docs/adr/0015-extrapolate-stage-2.md` — the `extrapolate` Stage-2
    disc-radius policy this module's Stage-2 fill honours.
"""
module VectorPathNetwork

using Random:               AbstractRNG, shuffle
using ..VectorProblems:     VectorPadeTaylorProblem
using ..VectorCoefficients: vector_taylor_coefficients
using ..VectorStepper:      VectorPadeStepperState, vector_pade_step_with_pade!
using ..VectorPathNetworkStage2: _stage2_fill
using ..VectorWedgeStep:    _select_wedge, _adaptive_h, WEDGE_STEP_POLICIES,
                            H_MIN_RATIO, VectorWalkError

export VectorPathNetworkSolution, vector_path_network_solve,
       VectorWalkFailure, VectorWalkError

# FW 2011 §3.1 ±22.5°/±45° wedge — the v0.1 `PathNetwork.DEFAULT_WEDGE`.
const DEFAULT_WEDGE = [-π / 4, -π / 8, 0.0, π / 8, π / 4]

# -----------------------------------------------------------------------------
# The skipped-target record (ADR-0026 D1)
# -----------------------------------------------------------------------------

"""
    VectorWalkFailure{T}

A first-class record of one target the resilient Stage-1 walk could
**not** reach.  ADR-0026 D1 makes `vector_path_network_solve` skip a
failed per-target walk and continue (`on_target_failure = :skip`)
rather than aborting the whole run.  A skip is *not* a silent swallow
— that would violate CLAUDE.md Rule 1.  It is fail-*loud-by-accounting*:
every skipped target is captured as one of these records and surfaced
to the caller in `VectorPathNetworkSolution.failed_targets`, so a
non-empty list is an unmissable signal the caller must report.

## Fields

  - `target::Complex{T}`       — the target point the walk was trying
                                 to reach;
  - `reason::Symbol`           — the `VectorWalkError.reason` of the
                                 failure (`:unreachable` |
                                 `:all_candidates_failed` |
                                 `:step_collapse`);
  - `z_stuck::Complex{T}`      — the position the walk was at when it
                                 failed (the `z_cur` in scope at the
                                 throw — read directly from the walk
                                 state, never parsed from the
                                 exception message);
  - `residual_dist::T`         — `|z_stuck − target|`, how far short of
                                 the target the walk stopped.

Partial nodes a failed walk laid down *before* it failed are kept in
the `visited_*` arrays (the walk pushes only on an accepted step, so
each is a valid solution node) — a `VectorWalkFailure` records "this
target was not reached," never "discard the filament."
"""
struct VectorWalkFailure{T <: AbstractFloat}
    target        :: Complex{T}
    reason        :: Symbol
    z_stuck       :: Complex{T}
    residual_dist :: T
end

# -----------------------------------------------------------------------------
# Solution container
# -----------------------------------------------------------------------------

"""
    VectorPathNetworkSolution{T}

The Stage-1 visited tree of a vector path-network walk.  For a tree of
`N` visited nodes, all arrays have length `N`:

  - `visited_z::Vector{Complex{T}}`            — node positions in the
                                                 complex plane;
  - `visited_y::Vector{Vector{Complex{T}}}`    — the `d`-vector state
                                                 `y` at each node;
  - `visited_h::Vector{Complex{T}}`            — the step length used
                                                 to reach each node
                                                 (the root carries the
                                                 nominal `h`);
  - `visited_numerators::Vector{Vector{Vector{Complex{T}}}}` — the `d`
                                                 shared-`Q` numerator
                                                 polynomials per node;
  - `visited_denominator::Vector{Vector{Complex{T}}}` — the *single*
                                                 shared denominator `Q`
                                                 per node (the pole
                                                 extractor roots this);
  - `visited_parent::Vector{Int}`              — the index of the node
                                                 each was reached from;
                                                 the root has parent 0;
  - `visited_jets::Vector{Vector{Vector{Complex{T}}}}` — the *raw*
                                                 (un-rescaled) order-`order`
                                                 vector Taylor jet at each
                                                 node: `d` coefficient
                                                 vectors `[c₀,…,c_order]`
                                                 in the `h'`-variable.
                                                 This is the **exact jet
                                                 the canonical shared-`Q`
                                                 Padé was built from**
                                                 (`_canonical_pade` already
                                                 computes it via
                                                 `vector_taylor_coefficients`);
                                                 caching it makes the B1
                                                 true-radius Stage-2 gate
                                                 cost nothing extra
                                                 (ADR-0025 Amendment 1 §A1,
                                                 "cache the order-24 jet").

The edge set `{(visited_parent[k], k) : k ≥ 2}` is the Stage-1 path
tree.  Element type is `Complex{T}` throughout — path-network steps
are complex-valued by construction even when the underlying problem is
real.

`visited_jets` is *additive* (ADR-0001): every v0.1 + v0.2 scalar test
is bit-identical, and the two pre-Stage-2 backward-compat constructors
(6-arg, 8-arg) default it to an empty vector.  When it is empty the
Stage-2 true-radius gate falls back to recomputing the jet from the
node `(z, y)` state — see `VectorPathNetworkStage2._stage2_fill`.

The two Stage-2 fields carry the fine-grid fill (bead `padetaylor-0ln.20`):

  - `grid_z::Vector{Complex{T}}`               — the fine-grid points
                                                 (the `fine_grid` kwarg
                                                 of `vector_path_network_solve`
                                                 verbatim);
  - `grid_y::Vector{Vector{Complex{T}}}`       — the dense interpolated
                                                 `d`-vector state `y(z)`
                                                 at each grid point; a
                                                 `d`-vector of `NaN + NaN·im`
                                                 marks a grid point not
                                                 covered by any visited
                                                 node's disc (fail-soft,
                                                 ADR-0015).

When `vector_path_network_solve` is called without a `fine_grid`, both
Stage-2 fields are empty — the V7 scatter-figure call sites are
byte-unaffected.

The tenth field carries the resilient-walk accounting (ADR-0026 D1):

  - `failed_targets::Vector{VectorWalkFailure{T}}` — one record per
                                                 target the walk could
                                                 not reach.  Empty
                                                 unless the solve was
                                                 called with
                                                 `on_target_failure =
                                                 :skip` *and* some
                                                 per-target walk
                                                 failed; a non-empty
                                                 list is the
                                                 fail-loud-by-accounting
                                                 signal the caller must
                                                 surface (see
                                                 `VectorWalkFailure`).

`failed_targets` is *additive* (ADR-0001) exactly as `visited_jets`
was: three backward-compat constructors (6-arg Stage-1-only, 8-arg
Stage-1+Stage-2, 9-arg pre-D1) default it — and the trailing Stage-2 /
jet fields — to empty vectors, so every hand-built test fixture and
every `on_target_failure = :throw` call site keeps its exact shape.
"""
struct VectorPathNetworkSolution{T <: AbstractFloat}
    visited_z           :: Vector{Complex{T}}
    visited_y           :: Vector{Vector{Complex{T}}}
    visited_h           :: Vector{Complex{T}}
    visited_numerators  :: Vector{Vector{Vector{Complex{T}}}}
    visited_denominator :: Vector{Vector{Complex{T}}}
    visited_parent      :: Vector{Int}
    grid_z              :: Vector{Complex{T}}
    grid_y              :: Vector{Vector{Complex{T}}}
    visited_jets        :: Vector{Vector{Vector{Complex{T}}}}
    failed_targets      :: Vector{VectorWalkFailure{T}}
end

# Backward-compat 6-arg constructor: a Stage-1-only solution (no
# fine-grid fill, no cached jets, no failed-target accounting).
# Defaults `grid_z`/`grid_y`/`visited_jets`/`failed_targets` to empty
# vectors so the V7 scatter-figure call sites and the hand-built test
# fixtures that predate Stage 2 (bead padetaylor-0ln.20) keep their
# exact shape.
VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar) where {T} =
    VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar,
                                 Complex{T}[], Vector{Complex{T}}[],
                                 Vector{Vector{Complex{T}}}[],
                                 VectorWalkFailure{T}[])

# Backward-compat 8-arg constructor: a Stage-1+Stage-2 solution built
# before the B1 jet cache (ADR-0025 Amendment 1) added `visited_jets`.
# Defaults `visited_jets` to empty — the Stage-2 true-radius gate then
# recomputes each node's jet from its `(z, y)` state (the documented
# fallback path in `VectorPathNetworkStage2._stage2_fill`) — and
# `failed_targets` to empty (the pre-ADR-0026 non-resilient walk).
VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar, gz, gy) where {T} =
    VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar, gz, gy,
                                 Vector{Vector{Complex{T}}}[],
                                 VectorWalkFailure{T}[])

# Backward-compat 9-arg constructor: a Stage-1+Stage-2+jet-cache
# solution built before ADR-0026 D1 added the resilient-walk
# `failed_targets` field.  Defaults it to empty — the non-resilient
# `on_target_failure = :throw` walk skips no target, so it has nothing
# to record.
VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar, gz, gy,
                             vjet) where {T} =
    VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar, gz, gy, vjet,
                                 VectorWalkFailure{T}[])

# -----------------------------------------------------------------------------
# Public driver
# -----------------------------------------------------------------------------

"""
    vector_path_network_solve(prob::VectorPadeTaylorProblem, targets;
                              order = prob.order, h = 0.5,
                              wedge_angles = DEFAULT_WEDGE,
                              step_policy = :max_q_root,
                              adaptive = true,
                              max_steps_per_target = 1000,
                              fine_grid = nothing,
                              extrapolate = false,
                              tol = 1e-8,
                              rng = nothing,
                              on_target_failure = :throw)
        -> VectorPathNetworkSolution

Build the minimal Stage-1 vector path-network covering `targets` — a
vector of complex points spanning the region of interest — from
`prob`'s initial condition; optionally also perform the Stage-2
fine-grid fill.

The walk seeds `visited = {z_0}` at `prob.zspan[1]` with the IC `y0`.
For each `target`, it locates the nearest already-visited node and
walks toward `target` via a 5-direction wedge, choosing the wedge
direction by the pole-avoidance criterion `step_policy` (the principled
`:max_q_root` — maximise the landed node's distance to its nearest
shared-`Q` root — by default; `:min_y` for the V7 proxy; see the module
docstring and `VectorWedgeStep`).  Each landed point becomes a visited
node storing the shared-`Q` approximant `(numerators, denominator)` from
`VectorStepper.vector_pade_step_with_pade!`.

Kwargs:
  - `order::Integer`   — Taylor truncation degree (defaults to
                         `prob.order`).
  - `h::Number`        — the wedge step magnitude.  With `adaptive =
                         true` (the default) it is the *maximum* step
                         `h_max`; the walk floors at `h_min =
                         H_MIN_RATIO·h` (`VectorWedgeStep.H_MIN_RATIO =
                         1e-3`) and adapts in between (see `adaptive`
                         below).  With `adaptive = false` it is the
                         single fixed step (the V7 behaviour).
  - `wedge_angles`     — five angle offsets relative to the goal
                         direction; must have exactly 5 entries.
  - `step_policy`      — `:max_q_root` (default — the principled
                         shared-`Q`-root-distance selector, ADR-0025 B2)
                         or `:min_y` (the V7 min-‖y‖ proxy).
  - `adaptive`         — `Bool` (default `true`).  When `true`, the step
                         magnitude is capped per-step by the parent
                         node's shared-`Q` pole distance and grows
                         toward `h_max = h` where the pole field is
                         sparse (`VectorWedgeStep._adaptive_h`); a step
                         never overshoots a pole.  When `false`, the
                         single fixed step `h` is used (V7 behaviour).
  - `max_steps_per_target` — cap on per-target walk length.
  - `fine_grid`        — `Union{Nothing, AbstractVector}` (default
                         `nothing`).  When a vector of complex points,
                         the solve performs the Stage-2 fine-grid fill
                         (see the module docstring, "Stage 2"): each
                         grid point is evaluated against the nearest
                         visited node's shared-`Q` Padé and the result
                         is stored in `grid_y`.  When `nothing`, the
                         Stage-2 fields stay empty (V7 scatter-figure
                         behaviour, byte-unaffected).
  - `extrapolate`      — `Bool` (default `false`, ADR-0015).  Gates the
                         Stage-2 validity check.  When `false`, a grid
                         point outside the nearest visited node's
                         **true validity disc** `R_gate` (the B1
                         true-radius gate, ADR-0025 Amendment 1)
                         receives a `d`-vector of `NaN + NaN·im`
                         (fail-soft per Rule 1).  When `true`, the check
                         is skipped and the nearest-node Padé is
                         evaluated regardless — gap-free figure panels
                         at the cost of evaluating outside the verified
                         disc.  Unused when `fine_grid === nothing`.
  - `tol`              — `Real` (default `1e-8`).  The truncation
                         tolerance the B1 true-radius gate (ADR-0025
                         Amendment 1) honours: the per-node honest
                         validity radius is the Jorba–Zou step at this
                         `tol`.  The calibrated safety factor `s(tol)`
                         is `0.34 / 0.36 / 0.30` for `tol = 1e-6 /
                         1e-8 / 1e-10`; values between the calibration
                         points clamp to the nearest.  Unused when
                         `fine_grid === nothing` or `extrapolate=true`.
  - `rng`              — `Union{Nothing, AbstractRNG}` (default
                         `nothing`).  Controls the Stage-1 target
                         *processing order* — the FW/FFW double-run
                         accuracy indicator (VC-10, ADR-0025; FFW 2017
                         `:246-247`).  With `nothing` the targets are
                         walked in the *given* order (the V7/V8
                         behaviour, reproduced bit-identically — ADR-0001
                         additive).  With an `AbstractRNG` the targets
                         are `shuffle`d by it first: a different order
                         builds a different path-network tree and a
                         different visited-node sequence, so a pole
                         field independently extracted from two
                         differently-`rng`-seeded runs differs by the
                         walk's numerical error — the practical accuracy
                         estimate FW/FFW exploit.  The walk for a
                         **fixed** `rng` state stays fully deterministic
                         (same shuffle ⇒ same tree, bit-identical).
  - `on_target_failure` — `Symbol` (default `:throw`).  The resilient-
                         walk policy (ADR-0026 D1).  With `:throw` (the
                         default) a per-target walk failure — an
                         unreachable target, an all-five-candidates-fail,
                         or an adaptive-step collapse — aborts the whole
                         solve by propagating a `VectorWalkError`; this
                         is byte-identical to the pre-ADR-0026 behaviour.
                         With `:skip` the failing per-target walk is
                         caught, recorded as a `VectorWalkFailure` in the
                         returned solution's `failed_targets`, and the
                         walk continues to the next target.  This is
                         *not* a silent swallow (Rule 1): every skip is
                         first-class accounting data the caller must
                         surface.  Partial nodes a skipped walk laid
                         down before it failed are kept — they are valid
                         solution nodes.  An exception that is *not* a
                         `VectorWalkError` (an unrecognised bug) is
                         re-thrown even under `:skip` — the walk never
                         skips an error it does not understand.

Throws `ArgumentError` (Rule 1) for empty `targets`, non-positive `h`,
a `wedge_angles` not of length 5, an unknown `step_policy`, an
`on_target_failure` not in `(:throw, :skip)`, or an empty `fine_grid`
(an empty fine grid is a caller mistake — pass `nothing` to skip
Stage 2).  With `on_target_failure = :throw`, a `VectorWalkError` if a
target is unreachable in `max_steps_per_target` steps, all five wedge
candidates fail, or (with `adaptive = true`) the pole field forces the
step below `h_min`; with `:skip` those become `failed_targets` records
instead.  A numerical breakdown inside the stepper (step landed on a
pole) propagates unchanged.
"""
function vector_path_network_solve(prob::VectorPadeTaylorProblem{F, CT},
                                   targets::AbstractVector;
                                   order::Integer = prob.order,
                                   h::Number = 0.5,
                                   wedge_angles::AbstractVector{<:Real} =
                                       DEFAULT_WEDGE,
                                   step_policy::Symbol = :max_q_root,
                                   adaptive::Bool = true,
                                   max_steps_per_target::Integer = 1000,
                                   fine_grid::Union{Nothing, AbstractVector} =
                                       nothing,
                                   extrapolate::Bool = false,
                                   tol::Real = 1e-8,
                                   rng::Union{Nothing, AbstractRNG} = nothing,
                                   on_target_failure::Symbol = :throw
                                   ) where {F, CT}
    isempty(targets) && throw(ArgumentError(
        "vector_path_network_solve: targets is empty — there is nothing " *
        "for the Stage-1 walk to reach.  Suggestion: pass a non-empty " *
        "vector of complex target points spanning the region of interest."))
    fine_grid !== nothing && isempty(fine_grid) && throw(ArgumentError(
        "vector_path_network_solve: fine_grid is empty — an empty Stage-2 " *
        "grid is a caller mistake.  Suggestion: pass `fine_grid = nothing` " *
        "to skip the Stage-2 fill, or a non-empty vector of grid points."))
    abs(h) > 0 || throw(ArgumentError(
        "vector_path_network_solve: h must be a non-zero step magnitude " *
        "(got $h).  Suggestion: pass a positive real h (FW 2011 default 0.5)."))
    length(wedge_angles) == 5 || throw(ArgumentError(
        "vector_path_network_solve: wedge_angles must have exactly 5 " *
        "entries (got $(length(wedge_angles))); FW 2011 §3.1 fixes the " *
        "5-direction wedge.  Suggestion: pass the default ±22.5°/±45°."))
    step_policy in WEDGE_STEP_POLICIES || throw(ArgumentError(
        "vector_path_network_solve: step_policy must be one of " *
        "$(WEDGE_STEP_POLICIES) (got :$step_policy).  :max_q_root is the " *
        "principled shared-Q-root-distance selector (ADR-0025 B2); :min_y " *
        "is the V7 min-‖y‖ proxy."))
    on_target_failure in (:throw, :skip) || throw(ArgumentError(
        "vector_path_network_solve: on_target_failure must be :throw or " *
        ":skip (got :$on_target_failure).  :throw (default) aborts the " *
        "solve on the first unreachable target; :skip records it in " *
        "failed_targets and continues (ADR-0026 D1).  Suggestion: pass " *
        "one of those two symbols."))

    # Working float type: real part of the problem's element type.
    T  = float(real(CT))
    C  = Complex{T}
    h_max = T(abs(h))
    # Adaptive walk: `h` is the ceiling, the floor is H_MIN_RATIO·h_max.
    # Fixed walk: h_min = h_max so `_adaptive_h` is never consulted.
    h_min = adaptive ? T(H_MIN_RATIO) * h_max : h_max

    # Seed the visited tree at the IC.  The root carries its own
    # canonical shared-Q approximant — centred at z0, real step h_max —
    # so the pole extractor treats every node uniformly.  The root's
    # `visited_h` records `h_max`: with adaptive `h` each node's step
    # may differ, so `visited_h` is the per-node canonical step the
    # pole extractor and the Stage-2 gate map roots through.
    z0 = C(prob.zspan[1])
    y0 = Vector{C}(prob.y0)
    num0, den0, jet0 = _canonical_pade(prob.f, z0, y0, Int(order), h_max, C)

    visited_z           = C[z0]
    visited_y           = Vector{C}[y0]
    visited_h           = C[C(h_max)]
    visited_numerators  = Vector{Vector{C}}[num0]
    visited_denominator = Vector{C}[den0]
    visited_parent      = Int[0]
    visited_jets        = Vector{Vector{C}}[jet0]

    # Stage-1 target processing order.  With `rng === nothing` the
    # targets are walked in the order given — the V7/V8 behaviour,
    # bit-identical (ADR-0001).  With an `rng` they are `shuffle`d by it
    # — a different processing order builds a different path-network
    # tree, the FW/FFW double-run accuracy indicator (VC-10, ADR-0025;
    # FW 2011 §3.1 line 156 shuffles; FFW 2017 `:246-247` makes the
    # cross-order disagreement the diagnostic).  Deterministic for a
    # fixed `rng` state.
    target_list = collect(C, targets)
    rng === nothing || (target_list = shuffle(rng, target_list))

    # The resilient-walk accounting (ADR-0026 D1): one VectorWalkFailure
    # per skipped target.  Stays empty under `on_target_failure = :throw`
    # (a walk failure aborts before anything is recorded) and under
    # `:skip` when every walk succeeds.
    failed_targets = VectorWalkFailure{T}[]

    for target in target_list
        # Skip a target we already sit on.
        any(z -> abs(z - target) ≤ 10 * eps(T), visited_z) && continue

        idx_v = _nearest_visited(visited_z, target)
        z_cur  = visited_z[idx_v]
        y_cur  = visited_y[idx_v]
        parent = idx_v
        n_steps = 0

        # The per-target wedge walk.  Under `on_target_failure = :skip`
        # it runs inside a try/catch: a `VectorWalkError` (an unreachable
        # target, an all-candidates-fail, or an adaptive-step collapse)
        # is recorded in `failed_targets` from the *live* `z_cur` — never
        # parsed from the exception message — and the walk continues to
        # the next target.  Any *other* exception is a bug we do not
        # understand, so it is re-thrown (Rule 1 — never skip an error
        # you cannot classify).  Under `:throw` (the default) no
        # try/catch wraps the walk: a `VectorWalkError` propagates and
        # aborts the solve, byte-identically to the pre-ADR-0026 walk.
        try
            # Loop guard re-reads the *current* step magnitude each
            # iteration — with adaptive `h` it shrinks/grows per step.
            while abs(z_cur - target) > real(visited_h[parent])
                n_steps += 1
                n_steps > max_steps_per_target && throw(VectorWalkError(
                    :unreachable,
                    "vector_path_network_solve: target $target unreachable " *
                    "in $max_steps_per_target steps; current z = $z_cur, " *
                    "|z - target| = $(abs(z_cur - target)).  Suggestion: " *
                    "raise max_steps_per_target, shorten h, or widen the " *
                    "wedge."))

                goal_dir = angle(target - z_cur)

                # Adaptive step (B2): the parent's shared-Q caps `h` to
                # the nearest-pole distance; `adaptive=false` pins it at
                # h_max.
                h_cur = adaptive ?
                    _adaptive_h(visited_denominator[parent],
                                real(visited_h[parent]), h_max, h_min) :
                    h_max

                z_new, y_new =
                    _select_wedge(prob.f, z_cur, y_cur, Int(order), h_cur,
                                  goal_dir, wedge_angles, step_policy)

                # The node's CANONICAL shared-Q store: centred at z_new,
                # real step h_cur, so its roots map to the z-plane
                # correctly.
                num, den, jet = _canonical_pade(prob.f, z_new, y_new,
                                                Int(order), h_cur, C)

                push!(visited_z, z_new)
                push!(visited_y, y_new)
                push!(visited_h, C(h_cur))
                push!(visited_numerators, num)
                push!(visited_denominator, den)
                push!(visited_parent, parent)
                push!(visited_jets, jet)

                parent = length(visited_z)   # next step chains off this node
                z_cur, y_cur = z_new, y_new
            end
        catch e
            # `:throw` mode (or a non-walk exception) → re-raise; `:skip`
            # mode catching a genuine `VectorWalkError` → record and
            # continue.  Partial nodes this walk already pushed stay in
            # `visited_*` — they are valid solution nodes.
            (on_target_failure === :skip && e isa VectorWalkError) || rethrow()
            push!(failed_targets,
                  VectorWalkFailure{T}(target, e.reason, z_cur,
                                       T(abs(z_cur - target))))
            continue
        end
    end

    # Stage 2 — the fine-grid fill (bead padetaylor-0ln.20).  When no
    # `fine_grid` is supplied the two Stage-2 grid fields stay empty (the
    # V7 scatter-figure behaviour).  `visited_jets` is populated either
    # way — it is the cached order-`order` jet of the Stage-1 walk, and
    # carrying it costs nothing (it was computed by `_canonical_pade`).
    if fine_grid === nothing
        return VectorPathNetworkSolution{T}(visited_z, visited_y, visited_h,
                                            visited_numerators,
                                            visited_denominator,
                                            visited_parent,
                                            Complex{T}[], Vector{Complex{T}}[],
                                            visited_jets, failed_targets)
    end

    grid_z = collect(C, fine_grid)
    grid_y = _stage2_fill(grid_z, visited_z, visited_h, visited_numerators,
                          visited_denominator, visited_jets, extrapolate,
                          real(tol), prob.f, visited_y, Int(order), C,
                          _nearest_visited)

    return VectorPathNetworkSolution{T}(visited_z, visited_y, visited_h,
                                        visited_numerators,
                                        visited_denominator, visited_parent,
                                        grid_z, grid_y, visited_jets,
                                        failed_targets)
end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Nearest visited node to `target` by Euclidean distance, with a
# lexicographic (Re, then Im) tiebreak for reproducibility — the v0.1
# `PathNetwork._nearest_visited` convention verbatim.
function _nearest_visited(visited_z::Vector{Complex{T}},
                          target::Complex{T}) where {T}
    idx, best = 1, abs(visited_z[1] - target)
    @inbounds for i in 2:length(visited_z)
        d = abs(visited_z[i] - target)
        if d < best || (d == best &&
                        (real(visited_z[i]) < real(visited_z[idx]) ||
                         (real(visited_z[i]) == real(visited_z[idx]) &&
                          imag(visited_z[i]) < imag(visited_z[idx]))))
            idx, best = i, d
        end
    end
    return idx
end

# The wedge-direction selector (`_select_wedge`) and the adaptive
# step-size controller (`_adaptive_h`) live in the sibling module
# `VectorWedgeStep` — split out for the CLAUDE.md Rule 6 200-LOC cap
# (B2 — bead `padetaylor-0ln.37.6`).  That module holds the principled
# `:max_q_root` shared-Q-root-distance criterion and the pole-capped
# adaptive `h`; read its literate docstring for the full design.

# Build the canonical shared-Q approximant for a node: `d` numerator
# polynomials over one denominator `Q`, centred at `(z, y)` and in the
# rescaled variable `t = (z' - z)/h` with REAL step `h`.  We run one
# `vector_pade_step_with_pade!` purely for its returned approximant and
# discard the landed state — exactly v0.1's `_local_pade` design.  A
# `DomainError` (the canonical step lands on a pole at t = 1) propagates
# unchanged: a node sitting on a pole is a fail-loud condition (Rule 1).
# It also returns the raw order-`order` vector Taylor jet — the jet the
# canonical Padé was built from — for the B1 jet cache (see the
# `visited_jets` field docstring on `VectorPathNetworkSolution`).
function _canonical_pade(f, z::Complex{T}, y::Vector{Complex{T}},
                         order::Int, h_mag::T,
                         ::Type{C}) where {T, C}
    state = VectorPadeStepperState{C}(z, y)
    _, numerators, denominator =
        vector_pade_step_with_pade!(state, f, order, C(h_mag))
    jets = vector_taylor_coefficients(f, C(z), Vector{C}(y), order)
    return numerators, denominator, jets
end

# The Stage-2 fine-grid fill (`_stage2_fill`), the wedge-step selector
# (`_select_wedge`) and the adaptive-h controller (`_adaptive_h`) live
# in the sibling modules `VectorPathNetworkStage2` / `VectorWedgeStep`,
# `include`d before this file (CLAUDE.md Rule 6 200-LOC split) and
# pulled in via the `using`s at the top; neither has a reverse
# dependency on this module.

end # module VectorPathNetwork
