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

The one essential change from the scalar v0.1 walk is the
**wedge-selection criterion**.  v0.1 picks the wedge direction that
minimises `|u|` — a heuristic for "step away from the next pole",
since a scalar meromorphic solution's `|u|` grows toward a pole.  The
vector analogue minimises the **norm** `‖y‖` of the vector state: a
vector meromorphic solution's components all blow up together at a
shared movable pole (ADR-0019), so `‖y‖` is the natural scalar proxy
for "distance to the next pole".  We mirror v0.1's `:min_u` as
`:min_y` here.

A more principled refinement — picking the wedge candidate that
*maximises the distance from the nearest root of the local shared-`Q`
denominator* — is deferred as bead `padetaylor-0ln.23`; it would be
forced if the `min-‖y‖` heuristic is observed to step a vector walk
onto a pole in a production figure.  `min-‖y‖` is the V7 default
because it is the direct, verified vector lift of the v0.1 criterion
the figure pipeline already trusts.

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
an unreachable target, and all-five-wedge-candidates-fail all throw
with a `Suggestion`.  Numerical breakdowns inside
`vector_pade_step_with_pade!` (a step landing on a pole — `Q(1) ≈ 0`)
propagate unchanged; the wedge selection is *meant* to steer around
them, and a thrown `DomainError` is the honest fail-loud signal that
it could not.

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

using LinearAlgebra:    norm
using ..VectorProblems: VectorPadeTaylorProblem
using ..VectorStepper:  VectorPadeStepperState, vector_pade_step_with_pade!
using ..VectorPathNetworkStage2: _stage2_fill

export VectorPathNetworkSolution, vector_path_network_solve

# FW 2011 §3.1 ±22.5°/±45° wedge — the v0.1 `PathNetwork.DEFAULT_WEDGE`.
const DEFAULT_WEDGE = [-π / 4, -π / 8, 0.0, π / 8, π / 4]

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
                                                 the root has parent 0.

The edge set `{(visited_parent[k], k) : k ≥ 2}` is the Stage-1 path
tree.  Element type is `Complex{T}` throughout — path-network steps
are complex-valued by construction even when the underlying problem is
real.

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
byte-unaffected.  A 6-arg backward-compat constructor (Stage-1 fields
only) defaults `grid_z`/`grid_y` to empty vectors.
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
end

# Backward-compat 6-arg constructor: a Stage-1-only solution (no
# fine-grid fill).  Defaults `grid_z`/`grid_y` to empty vectors so the
# V7 scatter-figure call sites and the hand-built test fixtures that
# predate Stage 2 (bead padetaylor-0ln.20) keep their exact shape.
VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar) where {T} =
    VectorPathNetworkSolution{T}(vz, vy, vh, vnum, vden, vpar,
                                 Complex{T}[], Vector{Complex{T}}[])

# -----------------------------------------------------------------------------
# Public driver
# -----------------------------------------------------------------------------

"""
    vector_path_network_solve(prob::VectorPadeTaylorProblem, targets;
                              order = prob.order, h = 0.5,
                              wedge_angles = DEFAULT_WEDGE,
                              step_policy = :min_y,
                              max_steps_per_target = 1000,
                              fine_grid = nothing,
                              extrapolate = false)
        -> VectorPathNetworkSolution

Build the minimal Stage-1 vector path-network covering `targets` — a
vector of complex points spanning the region of interest — from
`prob`'s initial condition; optionally also perform the Stage-2
fine-grid fill.

The walk seeds `visited = {z_0}` at `prob.zspan[1]` with the IC `y0`.
For each `target`, it locates the nearest already-visited node and
walks toward `target` via a 5-direction wedge, choosing the wedge
direction by the pole-avoidance criterion `step_policy` (only `:min_y`
— minimise the new state's norm `‖y‖` — is supported at V7; see the
module docstring).  Each landed point becomes a visited node storing
the shared-`Q` approximant `(numerators, denominator)` from
`VectorStepper.vector_pade_step_with_pade!`.

Kwargs:
  - `order::Integer`   — Taylor truncation degree (defaults to
                         `prob.order`).
  - `h::Number`        — wedge step magnitude (FW 2011 default 0.5).
  - `wedge_angles`     — five angle offsets relative to the goal
                         direction; must have exactly 5 entries.
  - `step_policy`      — `:min_y` (the only V7 value).
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
                         Stage-2 disc-radius check.  When `false`, a
                         grid point more than `real(visited_h[idx])`
                         from its nearest visited node receives a
                         `d`-vector of `NaN + NaN·im` (fail-soft per
                         Rule 1).  When `true`, the check is skipped and
                         the nearest-node Padé is evaluated regardless
                         of `|t|` vs 1 — gap-free figure panels at the
                         cost of degraded accuracy past `|t| = 1`.
                         Unused when `fine_grid === nothing`.

Throws `ArgumentError` (Rule 1) for empty `targets`, non-positive `h`,
a `wedge_angles` not of length 5, an unknown `step_policy`, or an empty
`fine_grid` (an empty fine grid is a caller mistake — pass `nothing` to
skip Stage 2); `ErrorException` if a target is unreachable in
`max_steps_per_target` steps or all five wedge candidates fail.  A
numerical breakdown inside the stepper (step landed on a pole)
propagates unchanged.
"""
function vector_path_network_solve(prob::VectorPadeTaylorProblem{F, CT},
                                   targets::AbstractVector;
                                   order::Integer = prob.order,
                                   h::Number = 0.5,
                                   wedge_angles::AbstractVector{<:Real} =
                                       DEFAULT_WEDGE,
                                   step_policy::Symbol = :min_y,
                                   max_steps_per_target::Integer = 1000,
                                   fine_grid::Union{Nothing, AbstractVector} =
                                       nothing,
                                   extrapolate::Bool = false
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
    step_policy === :min_y || throw(ArgumentError(
        "vector_path_network_solve: step_policy must be :min_y (got " *
        ":$step_policy); V7 supports only the min-‖y‖ wedge criterion.  " *
        "The shared-Q root-distance refinement is deferred — bead " *
        "padetaylor-0ln.23."))

    # Working float type: real part of the problem's element type.
    T  = float(real(CT))
    C  = Complex{T}
    h_mag = T(abs(h))

    # Seed the visited tree at the IC.  The root carries its own
    # canonical shared-Q approximant — centred at z0, real step h_mag —
    # so the pole extractor treats every node uniformly.
    z0 = C(prob.zspan[1])
    y0 = Vector{C}(prob.y0)
    num0, den0 = _canonical_pade(prob.f, z0, y0, Int(order), h_mag, C)

    visited_z           = C[z0]
    visited_y           = Vector{C}[y0]
    visited_h           = C[C(h_mag)]
    visited_numerators  = Vector{Vector{C}}[num0]
    visited_denominator = Vector{C}[den0]
    visited_parent      = Int[0]

    for target in collect(C, targets)
        # Skip a target we already sit on.
        any(z -> abs(z - target) ≤ 10 * eps(T), visited_z) && continue

        idx_v = _nearest_visited(visited_z, target)
        z_cur  = visited_z[idx_v]
        y_cur  = visited_y[idx_v]
        parent = idx_v
        n_steps = 0

        while abs(z_cur - target) > h_mag
            n_steps += 1
            n_steps > max_steps_per_target && throw(ErrorException(
                "vector_path_network_solve: target $target unreachable in " *
                "$max_steps_per_target steps; current z = $z_cur, " *
                "|z - target| = $(abs(z_cur - target)).  Suggestion: raise " *
                "max_steps_per_target, shorten h, or widen the wedge."))

            goal_dir = angle(target - z_cur)
            z_new, y_new =
                _wedge_step(prob.f, z_cur, y_cur, Int(order), h_mag,
                            goal_dir, wedge_angles, step_policy)

            # Recompute the node's CANONICAL shared-Q store: centred at
            # z_new, real step h_mag.  The wedge-step approximant is
            # centred at the parent in a complex-h variable and would
            # mis-map its roots (see the module docstring).
            num, den = _canonical_pade(prob.f, z_new, y_new, Int(order),
                                       h_mag, C)

            push!(visited_z, z_new)
            push!(visited_y, y_new)
            push!(visited_h, C(h_mag))
            push!(visited_numerators, num)
            push!(visited_denominator, den)
            push!(visited_parent, parent)

            parent = length(visited_z)     # next step chains off this node
            z_cur, y_cur = z_new, y_new
        end
    end

    # Stage 2 — the fine-grid fill (bead padetaylor-0ln.20).  When no
    # `fine_grid` is supplied the two Stage-2 fields stay empty (the V7
    # scatter-figure behaviour).  Otherwise every grid point is
    # evaluated against its nearest visited node's shared-Q Padé.
    if fine_grid === nothing
        return VectorPathNetworkSolution{T}(visited_z, visited_y, visited_h,
                                            visited_numerators,
                                            visited_denominator,
                                            visited_parent)
    end

    grid_z = collect(C, fine_grid)
    grid_y = _stage2_fill(grid_z, visited_z, visited_h, visited_numerators,
                          visited_denominator, extrapolate, C,
                          _nearest_visited)

    return VectorPathNetworkSolution{T}(visited_z, visited_y, visited_h,
                                        visited_numerators,
                                        visited_denominator, visited_parent,
                                        grid_z, grid_y)
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

# Take one wedge step: evaluate the 5 candidate directions, pick the one
# that minimises ‖y_new‖ (the vector lift of v0.1's min-|u| criterion),
# and return (z_new, y_new).  A candidate whose stepper call throws is
# sentinelled (skipped) so it loses the min-‖y‖ selection — the same
# shape as v0.1's failed-candidate handling.
function _wedge_step(f, z_cur::Complex{T}, y_cur::Vector{Complex{T}},
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

    # min-‖y‖ wedge selection: the candidate with the smallest new-state
    # norm is the one stepping furthest from the (shared) movable pole.
    best_k, best_val = 0, T(Inf)
    for k in 1:n_w
        cand_ok[k] || continue
        v = norm(cand_y[k])
        if v < best_val
            best_k, best_val = k, v
        end
    end
    best_k == 0 && throw(ErrorException(
        "vector_path_network_solve: all 5 wedge candidates failed at " *
        "z = $z_cur; every candidate step threw (e.g. landed on a pole " *
        "of the local shared-Q approximant).  Suggestion: shorten h or " *
        "widen the wedge."))

    return cand_z[best_k], cand_y[best_k]
end

# Build the canonical shared-Q approximant for a node: `d` numerator
# polynomials over one denominator `Q`, centred at `(z, y)` and in the
# rescaled variable `t = (z' - z)/h` with REAL step `h`.  We run one
# `vector_pade_step_with_pade!` purely for its returned approximant and
# discard the landed state — exactly v0.1's `_local_pade` design.  A
# `DomainError` (the canonical step lands on a pole at t = 1) propagates
# unchanged: a node sitting on a pole is a fail-loud condition (Rule 1).
function _canonical_pade(f, z::Complex{T}, y::Vector{Complex{T}},
                         order::Int, h_mag::T,
                         ::Type{C}) where {T, C}
    state = VectorPadeStepperState{C}(z, y)
    _, numerators, denominator =
        vector_pade_step_with_pade!(state, f, order, C(h_mag))
    return numerators, denominator
end

# The Stage-2 fine-grid fill (`_stage2_fill`) and its Horner evaluator
# (`_eval_poly`) live in the sibling module `VectorPathNetworkStage2`,
# `include`d *before* this file (see `src/PadeTaylor.jl`) and pulled in
# via the `using` at the top.  Splitting them out keeps this file under
# the CLAUDE.md Rule 6 200 effective-LOC cap.  Stage 2 has no reverse
# dependency on this module — the nearest-visited scan it needs is
# passed to `_stage2_fill` as a function argument — so the include
# order is Stage 2 first, then this Stage-1 driver.

end # module VectorPathNetwork
