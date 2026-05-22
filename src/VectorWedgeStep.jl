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

## Adaptive step size

Fixed `h` is the A2 failure.  The adaptive controller threads the
pole field by two coupled rules:

  1. **Pole cap (never overshoot a pole).**  Before a step, the
     parent node already carries its canonical shared-`Q` (the walk
     recomputes it at every landed node).  Its roots `t*` are the
     local poles; the *nearest*, `min|t*|`, sits at z-plane distance
     `h_node·min|t*|` (the canonical `Q` lives in `t = (z' -
     z_node)/h_node`).  The step is capped at `POLE_SAFETY·h_node·
     min|t*|` — landing *short* of that pole.

     The cap is **direction-agnostic** — the nearest pole in *any*
     direction, not the FW 2011 §3.1 `step_pade_root` forward
     projection ADR-0021 anticipated.  This is a deliberate B2
     deviation: the `:max_q_root` selector picks the step *direction*
     to dodge poles, so a directional cap would shrink `h` for a pole
     the walk steers around — needlessly stalling the walk (it stalled
     the figure walk to `h ≈ 3·10⁻⁶` in testing).  Capping by the
     closest pole in any direction keeps the step short enough that
     whatever direction the selector chooses still lands short — the
     honest coupling of the cap with the selector.

  2. **Geometric grow / shrink.**  Where the pole field is sparse the
     cap is inert and `h` grows geometrically toward `h_max` (a clean
     step ⇒ `h ← min(h_max, GROW · h)`); where it is dense the cap
     binds and `h` is pulled down to it.  `h` is never allowed below
     `h_min` (a hard floor — a walk that needs a sub-`h_min` step to
     dodge a pole is genuinely wedged, and a fail-loud throw is the
     honest signal, Rule 1).

The net effect, measured against the A2 probe: a single filament
threads to `|x| = 20`, and an adaptive walk shrinks `h` exactly where
the fixed-`h = 0.1` walk landed on a pole — the wedge is threaded
without the `shared_denominator_pade` degeneration A2 documented.

## Fail-loud contract (Rule 1)

`_adaptive_h` throws `VectorWalkError(:step_collapse, …)` with a
`Suggestion` when the pole cap forces `h` below `h_min` — the walk
cannot make honest progress and must say so, not silently take a
degenerate step.  `_select_wedge` throws
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
       VectorWalkError

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

# The step-size growth / shrink constants and the pole-distance safety
# factor.  GROW > 1 lets `h` climb toward `h_max` in pole-sparse
# regions; POLE_SAFETY < 1 makes a capped step land *short* of the pole
# (FW 2011 §3.1 — never step onto the singularity the Padé flagged).
const GROW        = 1.5
const POLE_SAFETY = 0.5

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

The next step magnitude for the adaptive wedge walk (see the module
docstring, "Adaptive step size").

`denominator` is the *parent* node's canonical shared-`Q` polynomial
(low-to-high coefficients, in the rescaled variable `t = (z' -
z_parent)/h_parent`).  `h_prev` is the step that *reached* the parent;
`h_max` / `h_min` bound the adaptive range.

Two coupled rules:

  1. **Geometric grow.**  Start from `h_grow = min(h_max, GROW·h_prev)`
     — a clean previous step earns a larger next one, capped at `h_max`.
  2. **Pole cap.**  The shared `Q`'s roots `t*` are the system's local
     poles; the *nearest* — `min|t*|` — sits at z-plane distance
     `h_prev·min|t*|` (the canonical `Q` is scaled by `h_prev`).  The
     capped step is `POLE_SAFETY·h_prev·min|t*|`, landing *short* of
     that pole; `h = min(h_grow, pole_cap)`.

The cap is **direction-agnostic** (`min|t*|`, not `step_pade_root`'s
forward-projected distance).  ADR-0021 anticipated reusing the FW 2011
§3.1 *directional* `step_pade_root` on the shared `Q`; B2 deviates
deliberately and the reason is the coupling with the `:max_q_root`
selector.  `_select_wedge` picks the step *direction* to dodge a pole —
so capping `h` by the pole on the goal direction alone would shrink `h`
for a pole the walk then steers around, needlessly stalling the walk
(it stalled the figure walk to `h ≈ 3·10⁻⁶` in testing).  Capping by
the nearest pole in *any* direction is the honest coupling: it keeps
the step short enough that whatever direction `_select_wedge` chooses
still lands short of the closest pole, while letting `h` grow wherever
the closest pole is far.  This is the same `min|t*|` measure the B1
Stage-2 pole-adjacency clamp uses (`VectorPathNetworkStage2`).  A
constant `Q` (`length ≤ 1`) has no roots — the cap is then absent and
the grow rule governs.

Throws `VectorWalkError(:step_collapse, …)` (Rule 1) when the pole cap
forces `h` below `h_min`: a walk that needs a sub-`h_min` step to dodge
a pole is genuinely wedged, and a fail-loud throw is the honest signal
— never a silent degenerate step.
"""
function _adaptive_h(denominator::AbstractVector{C}, h_prev::T,
                     h_max::T, h_min::T) where {T, C}
    h_grow = min(h_max, T(GROW) * h_prev)

    # Pole cap: the nearest root of the canonical shared-Q (in the
    # t-variable) is a pole at z-distance h_prev·min|t*|.  Cap `h` short
    # of it by POLE_SAFETY so whatever wedge direction is chosen lands
    # short of the closest pole.
    h = h_grow
    if length(denominator) > 1
        rs = roots(Polynomial(collect(denominator)))
        if !isempty(rs)
            pole_cap = T(POLE_SAFETY) * h_prev * T(minimum(abs, rs))
            h = min(h, pole_cap)
        end
    end

    h < h_min && throw(VectorWalkError(:step_collapse,
        "vector_path_network_solve: adaptive step collapsed to h = $h " *
        "(< h_min = $h_min) — the shared-Q pole field is too dense here " *
        "for an honest step.  Suggestion: lower h_min if a finer walk is " *
        "acceptable, or accept this point as the walk's honest frontier."))
    return h
end

end # module VectorWedgeStep
