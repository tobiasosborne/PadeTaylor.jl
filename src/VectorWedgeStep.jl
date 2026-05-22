"""
    PadeTaylor.VectorWedgeStep

The **principled wedge-step machinery** of the vector path-network walk:
the shared-`Q`-root direction selector and the adaptive step-size
controller.  This module is the discharge of bead `padetaylor-0ln.23`
and the substance of bead `padetaylor-0ln.37.6` (ADR-0025 Phase B2,
Lever 2).

## Why this module exists — the A4 baseline and the A2 failure

V7's `VectorPathNetwork` shipped a **min-‖y‖** wedge selector and a
**fixed** step `h`.  Both were honest v1 corners, and Phase-A measured
exactly how they fail:

  - **A2** (`external/probes/wedge-tractability/REPORT.md`): the
    fixed-`h = 0.1` production fan walk **blocks past `|x| ≈ 8`** — a
    bridging walk between two targets crosses the pole field on a chord
    and a fixed-`h` chord step lands on a pole, degenerating the
    shared-`Q` linear solve (`shared_denominator_pade: every singular
    value is below τ`).  The whole solve then aborts.  The probe also
    found a *single* ray reaches `|x| = 20` at any `h`, and a finer `h`
    pushes the completing frontier out — so the fixed step is the
    binding geometric limit, not the wedge itself.

  - **A4** (ADR-0025 Amendment 3): the shipped min-‖y‖ V8b walk has
    **~70 % catastrophic loop closures** — min-‖y‖ is only a *proxy*
    for "step away from the next pole" and the A4 baseline shows it
    steps the walk into trouble.

ADR-0025 Lever 2 retires both corners.  This module replaces the proxy
with the principled criterion and makes `h` adaptive.

## The shared-`Q`-root direction selector

The v0.1 scalar walk picks the wedge direction that minimises `|u|` — a
heuristic, since a scalar meromorphic solution's `|u|` grows toward a
pole.  V7 lifted that to **min-‖y‖** for the vector case.  But ‖y‖ is
only a proxy: a vector state can have a modest norm and still sit one
short step from a movable pole.

The principled criterion (`docs/adr/0025` Lever 2, bead `0ln.23`): of
the candidate wedge directions, choose the one whose **landed node has
the largest distance to its nearest shared-`Q` denominator root** — the
most pole-free disc of validity.  This is *exact*, not a proxy: the
shared-`Q` denominator's roots ARE the system's local poles (ADR-0019),
and the canonical shared-`Q` approximant at the landed node
(`VectorPathNetwork._canonical_pade`, centred at the node, real step
`h`, rescaled variable `t = (z' - z_node)/h`) puts a pole at z-plane
distance `h·|t*|` for each root `t*`.  So the pole-free disc radius of a
candidate's landed node is

    d_pole(candidate) = h · min |t* : Q_canonical(t*) = 0|

and `_select_wedge_max_q_root` returns the candidate maximising it.  A
candidate whose canonical-Padé build throws (it landed *on* a pole) is
sentinelled out — it has a zero pole-free disc by construction.

The min-‖y‖ selector is kept available (`step_policy = :min_y`) as the
verified V7 behaviour, but `:max_q_root` is the **default** for vector
walks: it is the principled method ADR-0025 commits to.

## Adaptive step size — the scale-derived, non-ratcheting law

Fixed `h` is the A2 failure.  The first adaptive controller (B2) was a
**geometric grow / geometric shrink** pair — and ADR-0026 Amendment 4's
S1 diagnostic proved it a *geometric sink*: its recurrence
`h_next = h_prev · min(GROW, POLE_SAFETY·min|t*|)` multiplied `h_prev`
in *both* branches, so once `h` shrank in a sustained pole field it
could never recover, `h = 0` was an attractor, and 86 % of walk nodes
collapsed below `h < 0.05`.  Worse, the FW-faithful nearest-node
chaining propagated a collapsed `h` into every later walk that chained
off the poisoned node.

The current controller (ADR-0026 Amendment 4, the S2 step law) replaces
the recurrence with a step that is an **absolute** function of the
local pole field, with *no multiplicative dependence on `h_prev`*:

    h = clamp( SAFETY · D_local ,  h_min ,  h_max )

  - `D_local = h_prev · min|t*|` is the **`h`-independent absolute
    z-plane distance to the nearest pole**.  The parent's canonical
    shared-`Q` lives in `t = (z' - z_parent)/h_prev`, so a denominator
    root `t*` is a pole at z-plane distance `h_prev·|t*|`; the
    *nearest* root gives `D_local`.  S1 EXP D verified `D_local` is
    invariant to `h_prev` to 5 significant figures — it is a property
    of the node and the field, not of the path that reached it.

  - `SAFETY < 1` sizes the step to ≲ the honest B1 validity disc, so
    adjacent path-network filaments tile the wedge gap-free (ADR-0026
    Amendment 3, bottleneck 1: a step larger than the honest disc opens
    gaps *along every filament*).

  - The pole measure `min|t*|` is **direction-agnostic** — the nearest
    pole in *any* direction, not the FW 2011 §3.1 `step_pade_root`
    forward projection ADR-0021 anticipated.  Deliberate B2 deviation:
    the `:max_q_root` selector already picks the step *direction* to
    dodge poles, so a directional cap would shrink `h` for a pole the
    walk steers around.

Because `h` jumps straight to `SAFETY·D_local` with no memory of how
small `h_prev` became, **a node leaving a dense pocket recovers to
`h_max` in a single step** — the geometric sink is gone by
construction.  The old `GROW` grow-rate limiter is deleted outright:
`SAFETY·D_local` is already a parent-safe absolute distance, so there
is no unsafe increase left for a rate-limiter to gate (Rule 9 — no
unused machinery).  `h` is never allowed below `h_min`: a step whose
`SAFETY·D_local` is genuinely sub-`h_min` means the *true* local pole
density exceeds what `h_min` can honestly resolve, and a fail-loud
throw is the honest signal (Rule 1).

## Fail-loud contract (Rule 1)

`_adaptive_h` throws `VectorWalkError(:step_collapse, …)` with a
`Suggestion` when the *genuine* local pole density forces
`SAFETY·D_local` below `h_min` — the walk cannot make honest progress
and must say so, not silently take a degenerate step.  Under the
non-ratcheting law this fires only on genuine density, never on a
self-inflicted geometric ratchet (the retired sink fired it
constantly).  `_select_wedge` throws
`VectorWalkError(:all_candidates_failed, …)` when *every* candidate
failed (each either threw in the stepper or threw building its
canonical `Q`) — the all-five-wedge-candidates-fail condition.

Both were bare `ErrorException`s before ADR-0026 D1; the typed
`VectorWalkError` (defined at the top of this module) lets the
resilient driver classify a caught failure by `reason::Symbol` — see
that struct's docstring.  An *uncaught* `VectorWalkError` (the
`on_target_failure = :throw` default) prints byte-identically to the
old `ErrorException` via the `showerror` method.

## References

  - `src/VectorPathNetwork.jl` — the Stage-1 walk driver that loops
    this machinery; the module docstring's "Stage 1 — the wedge walk".
  - `src/VectorStepper.jl` — `vector_pade_step_with_pade!`, the per-step
    shared-`Q` primitive each candidate direction runs.
  - `src/VectorPathNetworkStage2.jl` — the B1 Stage-2 pole-adjacency
    clamp, the same `0.5·h·min|t*|` measure the adaptive cap uses.
  - `docs/adr/0025-headline-figure-re-resolution.md` — Lever 2 (the
    dense-wedge walk); Amendment 2 (the walk threads the pole field for
    extraction); Amendment 3 (the A4 ~70 % catastrophic-loop baseline).
  - `docs/adr/0021-vector-step-control.md` — anticipated reusing the FW
    2011 §3.1 `step_pade_root` on the shared `Q`; B2's `_adaptive_h`
    deviates (direction-agnostic `min|t*|` cap — see "Adaptive step
    size" above for why the directional projection mis-couples with
    the `:max_q_root` selector).
  - `src/StepControl.jl` — `step_pade_root`, the FW 2011 §3.1
    directional pole-distance heuristic ADR-0021 anticipated.
  - `external/probes/wedge-tractability/REPORT.md` — A2: the fixed-`h`
    walk blocks past `|x| ≈ 8`; a finer/adaptive `h` reaches `|x| = 20`.
  - CLAUDE.md Rules 1, 4, 6, 9, 10.
"""
module VectorWedgeStep

using LinearAlgebra:        norm
using Polynomials:          Polynomial, roots
using ..VectorStepper:      VectorPadeStepperState, vector_pade_step_with_pade!

export _select_wedge, _adaptive_h, WEDGE_STEP_POLICIES, H_MIN_RATIO,
       SAFETY, VectorWalkError

# -----------------------------------------------------------------------------
# The typed walk-failure exception
# -----------------------------------------------------------------------------

"""
    VectorWalkError <: Exception

The typed exception a vector Stage-1 wedge walk throws when it cannot
honestly reach a target.  It is the *classified* successor to the bare
`ErrorException`s the walk used to throw — defined here, in the module
`include`d before `VectorPathNetwork`, so both the wedge machinery
(`_select_wedge`, `_adaptive_h`) and the driver
(`vector_path_network_solve`'s unreachable-target site) can construct it.

## Why a typed exception, not message-string parsing

ADR-0026 D1 makes the Stage-1 walk *resilient*: with
`on_target_failure = :skip` the driver catches a per-target walk
failure, records it as first-class data, and continues to the next
target.  To record *which kind* of failure occurred — and to be sure
the caught exception is genuinely a walk failure and not some
unrelated bug (Rule 1: never skip an error you do not understand) —
the driver must classify the exception.  Classifying by `occursin`-ing
substrings of the message string would be brittle and is exactly the
anti-pattern CLAUDE.md Rule 2 warns against.  A `reason::Symbol` field
is the genuinely correct design: the `catch` block reads `e.reason`
directly and `isa VectorWalkError` is the unambiguous "this is a walk
failure I understand" test.

## Fields

  - `reason::Symbol` — one of:
      * `:unreachable`           — the per-target walk exceeded
        `max_steps_per_target` without coming within `h` of the target;
      * `:all_candidates_failed` — every one of the five wedge
        candidates threw or had a degenerate canonical shared-`Q`
        store (`_select_wedge`);
      * `:step_collapse`         — the adaptive pole cap forced the
        step magnitude below `h_min` (`_adaptive_h`).
  - `msg::String` — the full human-readable message, retaining the
    `Suggestion:` remediation line (Rule 1 — a fail-loud throw still
    carries actionable advice even when it is later caught and skipped).

`Base.showerror` prints `msg` verbatim, so an *uncaught*
`VectorWalkError` (the `:throw` default) reads exactly as the old
`ErrorException` did — byte-identical operator-facing text.
"""
struct VectorWalkError <: Exception
    reason :: Symbol
    msg    :: String
end

Base.showerror(io::IO, e::VectorWalkError) = print(io, e.msg)

# The step-size safety factor for the scale-derived adaptive law.
#
# ADR-0026 Amendment 4 retired the old geometric-grow / geometric-shrink
# pair (`GROW`, `POLE_SAFETY`): that recurrence,
# `h_next = h_prev · min(GROW, POLE_SAFETY·min|t*|)`, made *both* terms
# proportional to `h_prev` — a genuine geometric sink with `h = 0` an
# attractor (S1 verdict).  The new law has no `h_prev` factor at all:
# `h` is set as an *absolute* fraction of the local pole-free disc
# `D_local`, so a node leaving a dense pocket recovers to `h_max` in a
# single step.  `GROW` is therefore deleted outright — there is nothing
# left for a grow rate-limiter to limit (S1's "may be retained" option
# is declined: `SAFETY·D_local` is already a parent-safe absolute
# distance, so an increase gate would be dead code, and Rule 9 forbids
# carrying unused machinery).
#
# `SAFETY` is the one free knob.  ADR-0026 Amendment 3 measured the
# honest B1 validity disc at ≈ 0.106 × the local pole *spacing*, and the
# walk's step must stay ≲ that honest disc so adjacent filaments tile
# without gaps.  `D_local` is the distance to the *nearest* pole, which
# for a mid-field node is roughly half the local pole spacing — so the
# honest disc is ≈ 0.106 × (2·D_local) ≈ 0.21 × D_local.  `SAFETY` ≈ 0.21
# sizes the step *at* the honest disc; a value below that gives a small
# tiling-overlap budget.
#
# S6b/S7 (ADR-0026 Amendment 5 / Amendment 6) — `SAFETY` lowered
# 0.25 → 0.10.  The S6c inner-wedge sweep
# (`external/probes/coverage-plateau-f4e/probe_s6_sweep.jl`) measured
# that a smaller `SAFETY` lifts honest coverage: at order 36, `SAFETY`
# 0.25→0.15→0.10 lifts inner-wedge coverage 0.232→0.289→0.314 (+9 pp) —
# sizing the step to the honest B1 disc rather than ≈ 2.4× it closes the
# along-filament gaps Amendment 3 flagged as bottleneck 1.
#
# Amendment 5's S6b *attempt* to lower `SAFETY` was reverted at the time
# because it regressed the `src/` pole-extraction suite (VPN.1.1 /
# VPN.3.1 / VPN.4.3 / VPO.4 went RED — empty pole fields).  That was NOT
# a `SAFETY` defect: it was the `radius_t` *scale-fixing heresy* in
# `extract_poles_shared_q`.  The old far-root filter kept denominator
# roots by `|t*| ≤ radius_t` — a window in the *rescaled* variable
# `t = Δz/h` — i.e. a window of a fixed number of *per-node steps*.  A
# smaller `SAFETY` shrinks `h`, so the *same physical pole* at a fixed
# z-distance maps to a larger `|t*| = D/h`, falls outside the window,
# and is discarded — fewer nodes root it, the `min_support ≥ 2` filter
# empties the field.  S7 (bead `padetaylor-gt1`, ADR-0026 Amendment 6)
# fixed that: `extract_poles_shared_q` now filters by the root's
# **z-plane distance** `h_node·|t*|` against a scale-STABLE z-radius
# `radius_t·h_max` (`h_max` = the walk's step ceiling, `h`-independent).
# The filter no longer shrinks with the per-node `h`, so lowering
# `SAFETY` keeps `min_support` voting healthy and the `src/` suite
# GREEN — the coupling Amendment 5 hit is gone.  `SAFETY = 0.10` is the
# measured coverage-best value, now unblocked.
#
# S7 also surfaced — and fixed — a *second* facet of the same
# `|t*|`-heresy: `extract_poles_shared_q` chose each cluster's
# representative as the smallest-`|t*|` candidate, which under the
# adaptive per-node `h` is a coarse-`h` *far* node, not the accurate
# z-plane-closest one.  The representative is now ordered by z-plane
# node-to-pole distance (see `VectorPoleField.extract_poles_shared_q`).
# With the finer `SAFETY = 0.10` step the `:max_q_root` selector steers
# the walk to wider-of-the-pole nodes, so the VPO.4 ℘-oracle's lone far
# target was additionally given a near-pole companion so the walk
# genuinely threads *past* the pole — the FW-faithful "walk past the
# pole to extract it" (the test's stated intent), not a tolerance
# relaxation.
const SAFETY = 0.10

# The adaptive step floor, as a fraction of `h_max`: `h_min = H_MIN_RATIO
# · h_max`.  A walk threading a dense pole field legitimately takes small
# steps — the floor is not a "typical step" bound but a genuine
# wedged-walk threshold: a step below `h_min` means the pole cap has
# collapsed `h` so far the walk can make no honest progress, and a
# fail-loud throw (Rule 1) is the honest signal.  `1/1000` lets the walk
# thread close to a pole (e.g. a `1e-4` step in a `|x| ≤ 20` window)
# while still firing on a genuine degeneration.
const H_MIN_RATIO = 1.0e-3

# The wedge-direction selection policies this module supports.  The V7
# `:min_y` proxy is retained; `:max_q_root` is the ADR-0025 principled
# default.  `VectorPathNetwork` validates the caller's `step_policy`
# against this set.
const WEDGE_STEP_POLICIES = (:max_q_root, :min_y)

# The pole-disc "clear" cap, as a multiple of the step magnitude: a
# candidate whose nearest shared-`Q` pole sits beyond `CLEAR_CAP · h` is
# pole-free *for this step* and its `:max_q_root` score is clamped to
# `CLEAR_CAP·h`.  Without the clamp the selector would chase the tiny,
# candidate-dependent *noise* roots a robust-Padé returns far out for a
# (near-)pole-free system, instead of letting the goal-alignment
# tie-break pick the least-detour direction.  `CLEAR_CAP = 10` — a pole
# ten step-lengths away constrains nothing the next single step does.
const CLEAR_CAP   = 10.0

# -----------------------------------------------------------------------------
# Candidate quality — the canonical pole-free disc radius
# -----------------------------------------------------------------------------

"""
    _candidate_pole_disc(f, z_new, y_new, order, h_mag, ::Type{C}) -> Real

The z-plane radius of the pole-free disc of a candidate's *landed* node
`(z_new, y_new)`: `h_mag · min|t*|` over the roots `t*` of the canonical
shared-`Q` denominator built at `(z_new, y_new)` with the real step
`h_mag`.

This is the principled `:max_q_root` selection score (see the module
docstring).  The canonical shared-`Q` lives in `t = (z' - z_new)/h_mag`,
so a denominator root `t*` is a pole at z-plane distance `h_mag·|t*|`;
the *nearest* root bounds the honest disc.

The score is **clamped at `CLEAR_CAP · h_mag`**: a candidate whose
nearest pole is beyond ten step-lengths is pole-free for the purpose of
the next step, and all such candidates tie — letting `_select_wedge`'s
goal-alignment tie-break choose the least-detour direction.  Without the
clamp the selector chases the tiny noise roots a robust-Padé returns far
out for a (near-)pole-free system.  A constant denominator
(`length(Q) ≤ 1`, no roots) is the limiting pole-free case, also scored
`CLEAR_CAP · h_mag`.

Returns `-Inf` (a candidate that loses every comparison) when the
canonical shared-`Q` build itself throws: that happens precisely when
the landed node sits *on* a pole (`Q(1) ≈ 0` inside
`vector_pade_step_with_pade!`), i.e. the worst possible candidate.
"""
function _candidate_pole_disc(f, z_new::Complex{T}, y_new::Vector{Complex{T}},
                              order::Int, h_mag::T, ::Type{C}) where {T, C}
    clear = T(CLEAR_CAP) * h_mag
    state = VectorPadeStepperState{C}(z_new, y_new)
    local denominator
    try
        _, _, denominator =
            vector_pade_step_with_pade!(state, f, order, C(h_mag))
    catch
        # The canonical step at the landed node threw — the node is on a
        # pole.  The worst candidate: lose every max comparison.
        return T(-Inf)
    end
    length(denominator) ≤ 1 && return clear         # constant Q — no pole
    rs = roots(Polynomial(collect(denominator)))
    isempty(rs) && return clear
    return min(clear, h_mag * T(minimum(abs, rs)))
end

# -----------------------------------------------------------------------------
# Wedge-direction selector
# -----------------------------------------------------------------------------

"""
    _select_wedge(f, z_cur, y_cur, order, h_mag, goal_dir, wedge_angles,
                  step_policy) -> (z_new, y_new)

Take one wedge step: evaluate the five candidate directions
`goal_dir + wedge_angles[k]` and return the landed state `(z_new,
y_new)` of the chosen one.

  - `step_policy = :max_q_root` (the ADR-0025 / bead-`0ln.23` default):
    choose the candidate whose landed node has the **largest pole-free
    disc** — `_candidate_pole_disc` — the principled "step toward the
    most pole-free region" criterion.
  - `step_policy = :min_y`: the V7 proxy — choose the candidate
    minimising the landed-state norm `‖y_new‖`.  Retained for the V7
    behaviour and as the mutation-proof control.

A candidate whose stepper call throws is sentinelled out (it landed on
a pole — the stepper's fail-loud `DomainError`).  Under `:max_q_root` a
candidate that *steps* cleanly but whose *canonical* `Q` then degenerates
(`_candidate_pole_disc` returns `-Inf`) is **also** excluded — the
driver rebuilds that same canonical store at the landed node, so a
candidate that cannot produce one is genuinely unusable and picking it
would crash the driver.

Throws `VectorWalkError(:all_candidates_failed, …)` (Rule 1) when
*every* candidate is unusable — no candidate stepped cleanly, or
(under `:max_q_root`) every clean candidate's canonical store
degenerates.  The all-five-wedge-candidates-fail condition.
"""
function _select_wedge(f, z_cur::Complex{T}, y_cur::Vector{Complex{T}},
                       order::Int, h_mag::T, goal_dir,
                       wedge_angles::AbstractVector{<:Real},
                       step_policy::Symbol) where {T}
    C   = Complex{T}
    n_w = length(wedge_angles)
    cand_z  = Vector{C}(undef, n_w)
    cand_y  = Vector{Vector{C}}(undef, n_w)
    cand_ok = falses(n_w)

    for k in 1:n_w
        θ = T(goal_dir) + T(wedge_angles[k])
        h_step = C(h_mag * cos(θ), h_mag * sin(θ))
        state  = VectorPadeStepperState{C}(z_cur, y_cur)
        try
            vector_pade_step_with_pade!(state, f, order, h_step)
            cand_z[k]  = state.z
            cand_y[k]  = state.y
            cand_ok[k] = true
        catch
            cand_ok[k] = false
        end
    end

    # Score each surviving candidate as a `(primary, tie_break)` tuple,
    # maximised lexicographically.  `primary` is the `:max_q_root`
    # pole-free disc radius (`_candidate_pole_disc`) or, for `:min_y`,
    # `-‖y‖` (minimise ‖y‖ as a maximised score).  `tie_break` is
    # `-|wedge_angle|`: when the primary scores tie — notably a pole-free
    # region where every candidate's disc is unbounded — the selector
    # prefers the most goal-aligned direction (smallest wedge offset),
    # the least-detour step.  Without this tie-break a uniformly
    # pole-free region would always pick wedge angle 1 and wander 45°
    # off the goal.
    #
    # Under `:max_q_root`, a `_candidate_pole_disc` of `-Inf` flags a
    # candidate whose canonical store cannot be built — the driver would
    # crash on it.  Such a candidate is *excluded* (not merely scored
    # low): `best_k` only ever advances onto a candidate with a finite
    # `primary`, so if every candidate degenerates `best_k` stays 0 and
    # the all-candidates-fail throw fires.
    best_k, best_score = 0, (T(-Inf), T(-Inf))
    for k in 1:n_w
        cand_ok[k] || continue
        primary = step_policy === :max_q_root ?
            _candidate_pole_disc(f, cand_z[k], cand_y[k], order, h_mag, C) :
            -norm(cand_y[k])
        isfinite(primary) || continue          # canonical store unbuildable
        score = (primary, -abs(T(wedge_angles[k])))
        if best_k == 0 || score > best_score
            best_k, best_score = k, score
        end
    end
    best_k == 0 && throw(VectorWalkError(:all_candidates_failed,
        "vector_path_network_solve: all 5 wedge candidates failed at " *
        "z = $z_cur; every candidate step threw or its canonical " *
        "shared-Q store degenerated (a pole of the local shared-Q " *
        "approximant).  Suggestion: shorten h or widen the wedge."))

    return cand_z[best_k], cand_y[best_k]
end

# -----------------------------------------------------------------------------
# Adaptive step size
# -----------------------------------------------------------------------------

"""
    _adaptive_h(denominator, h_prev, h_max, h_min) -> Real

The next step magnitude for the adaptive wedge walk: a **scale-derived,
non-ratcheting** step law (see the module docstring, "Adaptive step
size").

`denominator` is the *parent* node's canonical shared-`Q` polynomial
(low-to-high coefficients, in the rescaled variable `t = (z' -
z_parent)/h_parent`).  `h_prev` is the step that *reached* the parent;
`h_max` / `h_min` bound the adaptive range.

## Why the old geometric-grow / geometric-shrink law was retired

The previous `_adaptive_h` was a **geometric sink** (ADR-0026
Amendment 4, the S1 diagnostic verdict).  Its recurrence was

    h_next = h_prev · min(GROW, POLE_SAFETY·min|t*|)

— *both* branches multiply `h_prev`, so the step has no `h`-independent
floor.  Whenever `min|t*| < GROW/POLE_SAFETY = 3` was sustained (a
mid-density pole field — exactly the P_I⁽²⁾ tritronquée wedge), `h`
shrank geometrically and `h = 0` was an attractor it could never climb
back from.  Measured: 86 % of walk nodes collapsed to `h < 0.05`, and
the FW-faithful nearest-node chaining then *propagated* a collapsed `h`
into every later walk that chained off the poisoned node — a walk once
stalled to `h ≈ 3·10⁻⁶` and could not recover.  The S1 diagnostic
confirmed the cap `POLE_SAFETY·h_prev·min|t*|` *is* a genuine absolute
z-distance (the two prior audits were a false dichotomy: a true
absolute cap can still be a sink, because the cap is a fraction of the
*parent's* disc and the parent's disc shrank with `h_prev`).

## The scale-derived law

`h` is now an **absolute** function of the local pole field with **no
multiplicative dependence on `h_prev`**:

    h = clamp( SAFETY · D_local ,  h_min ,  h_max )

where `D_local = h_prev · min|t*|` is the `h`-independent absolute
z-plane distance to the nearest pole.  The canonical shared-`Q` lives
in `t = (z' - z_parent)/h_prev`, so a denominator root `t*` is a pole
at z-plane distance `h_prev·|t*|`; the *nearest* root gives `D_local`.
S1 EXP D verified `D_local` is invariant to `h_prev` to 5 significant
figures across a 16× span of `h` — it is a property of the node and the
field, not of the path that reached it.  Because `h` jumps straight to
`SAFETY·D_local` with no memory of how small `h_prev` became, **a node
leaving a dense pocket recovers to `h_max` in a single step** — the
geometric sink is gone by construction.

`SAFETY` (the one free knob) sizes the step to ≲ the honest B1 validity
disc so adjacent filaments tile gap-free — see the `const SAFETY`
comment for the disc-spacing reasoning.  The old `GROW` grow-rate
limiter is **deleted entirely**: `SAFETY·D_local` is already a
parent-safe absolute distance, so there is no unsafe increase for a
rate-limiter to gate.

The pole measure is **direction-agnostic** (`min|t*|`, the nearest pole
in *any* direction, not `step_pade_root`'s forward projection).
ADR-0021 anticipated a *directional* cap; B2 deviates deliberately:
`_select_wedge` picks the step *direction* to dodge a pole, so a
directional measure would shrink `h` for a pole the walk then steers
around.  Capping by the nearest pole in any direction keeps the step
short enough that whatever direction the selector chooses still lands
short.  This is the same `min|t*|` measure the B1 Stage-2
pole-adjacency clamp uses (`VectorPathNetworkStage2`).

When the shared `Q` has **no roots** — a constant `Q` (`length ≤ 1`),
or `roots` returns an empty set — there is no nearby pole, `D_local` is
unbounded, and `h` is set to `h_max`.

## Fail-loud (Rule 1)

Throws `VectorWalkError(:step_collapse, …)` when `SAFETY·D_local`
falls below `h_min` — the case where the *true* local pole density
genuinely exceeds what `h_min` can honestly resolve.  Under the old
geometric-sink law this fired constantly as a *self-inflicted* ratchet
artefact; under the non-ratcheting law it fires only on genuine
density, so it should be rare.  The ADR-0026 D1 `:skip` driver catches
it and records the point as a `failed_targets` record — a fail-loud
throw is still the honest signal, never a silent degenerate step.
"""
function _adaptive_h(denominator::AbstractVector{C}, h_prev::T,
                     h_max::T, h_min::T) where {T, C}
    # D_local — the h-independent absolute z-plane distance to the
    # nearest pole.  The canonical shared-Q lives in t = Δz/h_prev, so a
    # root t* is a pole at z-distance h_prev·|t*|.  No nearby pole (a
    # constant Q, or an empty root set) ⇒ the step is unconstrained from
    # below ⇒ h = h_max.
    h_target = h_max
    if length(denominator) > 1
        rs = roots(Polynomial(collect(denominator)))
        if !isempty(rs)
            D_local  = h_prev * T(minimum(abs, rs))
            # The scale-derived step: an *absolute* fraction of the
            # local pole-free disc, with no h_prev factor — so a node
            # leaving a dense pocket jumps straight back to h_max.
            h_target = T(SAFETY) * D_local
        end
    end

    h = clamp(h_target, h_min, h_max)

    # The :step_collapse throw now fires only for the genuine case:
    # SAFETY·D_local is itself below h_min, i.e. the true local pole
    # density exceeds what h_min can honestly resolve.  (Under clamp,
    # `h == h_min` whenever h_target ≤ h_min, so the genuine-collapse
    # test is on h_target, not the clamped `h`.)
    h_target < h_min && throw(VectorWalkError(:step_collapse,
        "vector_path_network_solve: adaptive step collapsed to " *
        "SAFETY·D_local = $h_target (< h_min = $h_min) — the shared-Q " *
        "pole field is genuinely too dense here for an honest step " *
        "(the nearest pole sits within $(h_target / SAFETY) of the node). " *
        " Suggestion: lower h_min if a finer walk is acceptable, or " *
        "accept this point as the walk's honest frontier (the :skip " *
        "driver records it as a failed target)."))
    return h
end

end # module VectorWedgeStep
