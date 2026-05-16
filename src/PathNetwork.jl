"""
    PadeTaylor.PathNetwork

Tier-2 implementation of the Fornberg-Weideman 2011 §3.1 path-network:
build a tree of Padé-Taylor segments in the complex plane that bridges
poles of the solution by deliberate off-axis detours, then evaluate
the local stored approximant at any fine-grid lattice node within
distance `h` of a tree node.

This module is a *sibling driver* to `Problems.solve_pade`, not an
extension of it.  Where `solve_pade` walks one segment after another
in fixed (real) `h` and returns a 1-D trajectory, `path_network_solve`
walks **five candidate complex directions per step**, picks the one
that minimises `|u|` (heuristic for "furthest from the next pole"),
and records the visited tree.  The data structure cannot be unified
with `PadeTaylorSolution` without losing the Stage-2 evaluation
capability — see ADR-0004 §"Decision".

## The algorithm (FW 2011 §3.1, `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-166`)

**Stage 1 — path-tree construction**.  Initialise `visited = {z_0}`.
For each target `z_t` in (shuffled) `grid`: find the nearest already-
visited node `z_v` (Euclidean, lexicographic tiebreak); from `z_v`
step toward `z_t` using the 5-direction wedge until within `h` of
`z_t`.  Each landed point becomes a new entry in `visited`, with the
local Padé approximant *centered at that point* (computed by a
canonical-direction step from there — see below).

**Stage 2 — fine-grid extrapolation** (FW 2011 line 166).  For each
`z_f` in `grid` find the nearest visited `z_v`; if `|z_f - z_v| ≤ h`,
evaluate `visited_pade[k]` at `t = (z_f - z_v) / h`.  Otherwise the
slot gets `NaN + NaN·im` — no silent extrapolation.

## Wedge convention

Default: FW 2011 ±22.5°/±45° (`FW2011...md:158-159`).  Override via
`wedge_angles` kwarg, e.g. `[-π/6, -π/12, 0, π/12, π/6]` for RF 2014's
±15°/±30° codification (`ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:148-153`).

## Per-node canonical Padé

At every visited node we store **one** Padé approximant centered at
that node, built with the canonical (real, positive) step length `h`.
Cost: one extra `Coefficients.taylor_coefficients_2nd` + one extra
`RobustPade.robust_pade` per visited node, on top of the 5 wedge-
direction Padés used for selection.  In exchange Stage 2 always
*interpolates* inside `|t| ≤ 1` (the disc of validity of the canonical
Padé), never extrapolates past `t = 1`.  For Tier-2 ship this
simplifies the dense-output code and avoids accuracy loss on off-
direction lattice nodes.

## Schwarz-reflection symmetry (opt-in, worklog 014 + bead `padetaylor-dtj`)

For ODEs that preserve complex conjugation — real-coefficient `f` AND
real IC on the real axis (`Im(z₀) = Im(u₀) = Im(u'₀) = 0`) — the
analytic solution satisfies `u(z̄) = ū(z)` globally (Schwarz reflection
of an analytic function across its real-axis trace).  The default path-
network is **not** numerically Schwarz-symmetric: the `shuffle(rng,
targets)` step (FW 2011 line 156) creates an asymmetric visited tree,
and the Stage-2 nearest-visited lookup at `z` and `conj(z)` can land on
non-conjugate visited nodes — cascading 4-5 orders of magnitude into
`|u|` at far-from-IC conjugate-pair grid cells.

The opt-in kwarg `enforce_real_axis_symmetry = true` cures this: walk
the unique upper-half + on-axis representatives only, then populate the
output `grid_u` / `grid_up` by indexing the upper-walk's outputs
(mirroring via `conj` for input cells with `Im(z) < 0`).  The result is
bit-exact symmetric across conjugate pairs and runs roughly 2× faster
than a full-plane walk.  The validity precondition (real coefficients
in `f`) cannot be checked at this layer and is the caller's promise;
the IC-on-real-axis precondition IS checked and throws an
`ArgumentError` if violated (CLAUDE.md Rule 1).  Default `false`
preserves the FW 2011 algorithm verbatim.

## Adaptive Padé step size (`:adaptive_ffw`, opt-in)

FFW 2017 §2.1.2 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:74-97`) replaces the FW
2011 fixed `h` with a controller that estimates per-step truncation
error `T(h) = |ε_{n+1} h^{n+1} / a(h)|` and rescales `h := q·h` with
`q = (k·Tol/T(h))^(1/(n+1))` until `T(h) ≤ Tol`.  Accepted `|q·h|`
seeds the next step's initial `h` — the controller has memory.

The `:adaptive_ffw` policy (opt-in via `step_size_policy` kwarg)
delegates to `PadeStepper.adaptive_pade_step!` for each wedge-
direction step in Stage 1.  The controller magnitude is the real
absolute step length; the wedge direction is preserved across
rescales (a real `q ∈ (0, 1]` multiplies a complex `h` without
rotation).

The Stage-2 evaluation (canonical-Padé interpolation) is **unchanged**:
each visited node stores its own real `visited_h[k]`, and Stage 2
indexes the local Padé at `t = (z_f - z_v) / visited_h[k]` — heterogeneous
node-by-node `h`s are naturally supported by the per-node-`h` storage
already in place (worklog 034).

Activate via `step_size_policy = :adaptive_ffw` and pass an
`adaptive_tol` (default `1e-12`) and optional `k_conservative`
(FFW default `1e-3`).

## Non-uniform Stage-1 node placement (`node_separation`, opt-in)

FFW 2017 §2.1.2 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:67-72` + md:97) observes
that the pole densities of `P̃_III` / `P̃_V` solutions grow
*exponentially* with `Re ζ` under the exponential transformations.  A
uniform Stage-1 step is then over-resolved at small `Re ζ` (wide pole
spacing, the walker wastes effort on smooth regions) and under-resolved
at large `Re ζ` (pole spacing collapses below `h` and the walker can no
longer reliably bridge them).

FFW prescribe a **node-separation function** `R(ζ)` that specifies the
distance from a node at `ζ` to its neighbours; for their Fig 1 PIII
solution `R(ζ) = (8 - Re ζ)/20` (md:72), a linearly decreasing magnitude.
We expose this via the opt-in kwarg `node_separation::Union{Nothing,
Function} = nothing` on `path_network_solve`.

When `node_separation === nothing` (the default), behaviour is byte-
identical to the FW 2011 fixed-`h` walker; backward compat for every
existing test is the load-bearing invariant.  When a function, the
per-step initial step magnitude is `R(z_cur)` rather than the FFW
"scaled step length stored at the current point" (md:93).  Under
`step_size_policy = :adaptive_ffw` the adaptive controller may then
shrink `R(z_cur)` via the `q ≤ 1` rescale; it cannot grow it.  Composition
is therefore **R provides the spatial-density target; adaptive
provides temporal-accuracy refinement**; the two stack orthogonally.

`R(ζ::Complex{T}) -> T` must return a positive finite real.  We check
this at every step and throw `ArgumentError` on non-positive or non-
finite output (CLAUDE.md Rule 1).

## Tier-3-to-5 deferrals (per ADR-0004)

  - **BVP composition (Tier 3)**: stub comment at end of this file.
  - **CoordTransforms (Tier 4)**: external wrapper, no change here.
  - **SheetTracker (Tier 5)**: external wrapper, no change here.
"""
module PathNetwork

using Base.Threads:  nthreads, @threads
using Printf:        @sprintf, @printf
using Random:        shuffle, MersenneTwister
using ..RobustPade:  PadeApproximant
using ..PadeStepper: PadeStepperState, pade_step_with_pade!,
                     adaptive_pade_step!, ffw_truncation_error, ffw_rescale_q,
                     _evaluate_pade, _evaluate_pade_deriv
using ..Problems:    PadeTaylorProblem
using ..BranchTracker: resolve_cut_angles, any_cut_crossed, step_sheet_update

export PathNetworkSolution, path_network_solve, eval_at, eval_at_sheet

const DEFAULT_WEDGE = [-π/4, -π/8, 0.0, π/8, π/4]   # FW 2011 ±22.5°/±45°

# -----------------------------------------------------------------------------
# Solution container
# -----------------------------------------------------------------------------

"""
    PathNetworkSolution{T}

Result of `path_network_solve`.  `visited_*` arrays are the Stage-1
tree; `grid_*` arrays are the Stage-2 fine-grid evaluation.  All
fields use `Complex{T}` for the spatial coordinate even when the
underlying problem is real, since path-network steps are complex-
valued by construction.

`visited_parent[k]` is the index (into the `visited_*` arrays) of the
node from which `visited_z[k]` was reached during Stage-1 tree
construction; the root (IC) node has parent `0`.  The set of edges
`{(visited_parent[k], k) : k ≥ 2}` is exactly the Stage-1 path tree —
this is what FW 2011 Fig 3.2 draws.  Because the FW 2011 algorithm
starts each per-target walk from the *nearest already-visited node*
(line 164), the parent of the first node landed in a new walk is that
nearest node, while subsequent nodes of the same walk chain off their
immediate predecessor.

A `NaN + NaN·im` entry in `grid_u` or `grid_up` signals that the
corresponding `grid` point was not within `h` of any visited node;
fail-loud per CLAUDE.md Rule 1.

`visited_sheet[k]` is a per-branch Riemann-sheet index tuple
(`Vector{Int}`) populated by `BranchTracker` when `branch_points` is
non-empty (ADR-0013).  For the default `branch_points = ()` call,
every entry is `Int[]` (zero-length).  When N branches are supplied,
`length(visited_sheet[k]) == N` for every visited node, with the
root carrying the `initial_sheet` kwarg's value.

A 9-arg convenience constructor (without `visited_sheet`) is
provided for backward-compat with pre-ADR-0013 callers (e.g.
existing test fixtures and external wrappers); it defaults
`visited_sheet` to the empty `Vector{Int}[]`.
"""
struct PathNetworkSolution{T <: AbstractFloat}
    visited_z      :: Vector{Complex{T}}
    visited_u      :: Vector{Complex{T}}
    visited_up     :: Vector{Complex{T}}
    visited_pade   :: Vector{PadeApproximant{Complex{T}}}
    visited_h      :: Vector{T}                # real (canonical) per node
    visited_parent :: Vector{Int}              # tree edge: parent index, 0 = root
    grid_z         :: Vector{Complex{T}}
    grid_u         :: Vector{Complex{T}}
    grid_up        :: Vector{Complex{T}}
    visited_sheet  :: Vector{Vector{Int}}      # per-branch sheet tuple per node (ADR-0013)
end

# Backward-compat 9-arg constructor: defaults visited_sheet to empty.
PathNetworkSolution{T}(vz, vu, vup, vp, vh, vpar, gz, gu, gup) where {T} =
    PathNetworkSolution{T}(vz, vu, vup, vp, vh, vpar, gz, gu, gup,
                           Vector{Int}[])

# -----------------------------------------------------------------------------
# Public driver
# -----------------------------------------------------------------------------

"""
    path_network_solve(prob, grid; h=0.5, order=prob.order,
                       wedge_angles=DEFAULT_WEDGE,
                       step_selection=:min_u,
                       step_size_policy=:fixed,
                       node_separation=nothing,
                       max_steps_per_target=1000,
                       rng_seed=0,
                       enforce_real_axis_symmetry=false,
                       verbose=false, progress_every=500) -> PathNetworkSolution

Build the FW 2011 §3.1 path-network covering `grid::AbstractVector{<:Complex}`
from `prob`'s IC, then return the Stage-1 tree + Stage-2 evaluations.

Required: `prob isa PadeTaylorProblem` with 2nd-order `y0 = (u_0, up_0)`.
The IC point `prob.zspan[1]` becomes `visited_z[1]`.

Kwargs:
  - `h::Real`            — canonical step length (FW default 0.5).
  - `order::Integer`     — Taylor truncation degree.  Defaults to
    `prob.order` (the order the `PadeTaylorProblem` was built with);
    pass explicitly only to override it.  FW 2011's standard is 30.
  - `wedge_angles`       — five angle offsets relative to goal direction.
  - `step_selection`     — `:min_u` (default; FW 2011) or `:steepest_descent`.
  - `step_size_policy`   — `:fixed` (default; FW 2011 constant `h`) or
    `:adaptive_ffw` (FFW 2017 §2.1.2 controller; uses `adaptive_tol`
    and `k_conservative` kwargs).
  - `adaptive_tol`       — controller acceptance threshold for
    `step_size_policy = :adaptive_ffw` (default `1e-12`).  Unused
    under `:fixed`.
  - `k_conservative`     — FFW 2017 conservative factor `k` in
    `q = (k·Tol/T(h))^(1/(n+1))` (default `1e-3` per md:91).  Unused
    under `:fixed`.
  - `max_rescales`       — per-step cap on FFW rescale iterations
    under `:adaptive_ffw` (default `50`).  Unused under `:fixed`.
  - `node_separation`    — `Union{Nothing, Function}`.  When `nothing`
    (default), the walker advances by a constant `h` (FW 2011) or by
    the controller-memory-seeded magnitude (`:adaptive_ffw`).  When a
    function `R::Complex{T} -> T`, the per-step initial step magnitude
    is `R(z_cur)` at every Stage-1 walker step.  Under `:adaptive_ffw`,
    the controller may shrink that initial seed via the `q ≤ 1` rescale
    but cannot grow it — so the per-step accepted `|h|` is bounded
    above by `R(z_cur)`.  `R` must return a positive finite real
    (`ArgumentError` thrown at step time otherwise — Rule 1).
    FFW Fig 1 prescription is `R(ζ) = (8 - Re ζ)/20`
    (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
    FFW2017_painleve_riemann_surfaces_preprint.md:72`).
  - `branch_points`      — `Tuple{Vararg{Complex}}` (ADR-0013).
    Empty `()` default = no branches, byte-equivalent to the pre-A4
    walker.  When non-empty, the wedge selector consults
    `BranchTracker.any_cut_crossed` on each candidate to enforce the
    Riemann-sheet structure.  Each branch point carries a half-line
    cut from itself along direction `branch_cut_angles[k]` (default
    `arg = π`, Julia `log` convention).
  - `branch_cut_angles`  — `Real` or `Tuple{Vararg{Real}}` of cut
    directions per branch (default `π`).  Scalar broadcasts to all
    branches; tuple must match `length(branch_points)`.
  - `cross_branch`       — `Bool` (default `false`).  In refuse mode
    (`false`) any wedge candidate whose step segment crosses a cut
    is forbidden; if all 5 are forbidden the walker throws (CLAUDE.md
    Rule 1).  In cross mode (`true`) candidates are allowed to cross
    and the per-branch sheet counter (`visited_sheet`) is updated
    via `sign(winding_delta(z_cur, z_new, b_k))` for each crossed
    branch.  Throws if combined with `branch_points = ()`.
  - `initial_sheet`      — `AbstractVector{<:Integer}` of length
    `length(branch_points)` (default zeros).  The root visited node
    inherits this sheet tuple; descendants accumulate from there.
    Defaults to `(0, 0, …)` — "this IC sits on the principal sheet
    of every branch".
  - `grid_sheet`         — `Union{Nothing, Vector{Vector{Int}}}`
    (default `nothing`).  When supplied (A5 / bead `padetaylor-hed`),
    each entry `grid_sheet[i]` specifies the sheet tuple that grid
    point `grid[i]` is queried against.  The Stage-2 nearest-visited
    lookup is restricted to visited nodes whose `visited_sheet[k]
    == grid_sheet[i]`; if no matching-sheet visited node lies within
    `visited_h[k]` of `grid[i]`, the Stage-2 output is `NaN + NaN·im`
    (same fail-soft pattern as the existing too-far-away case).  Must
    match `length(grid)`; each inner vector must have length
    `length(branch_points)`.  Throws on mismatch (Rule 1).  Default
    `nothing` preserves the existing sheet-agnostic Euclidean
    nearest-visited semantics.
  - `extrapolate`        — `Bool` (default `false`, ADR-0015).  When
    `false`, Stage-2 cells more than `visited_h[idx]` from the
    nearest visited node receive `NaN + NaN·im` — fail-soft per
    ADR-0004 + CLAUDE.md Rule 1.  When `true`, the disc-radius
    check is skipped and the nearest-node Padé is evaluated
    regardless of `|t|` vs 1 — matches FFW 2017 §2.1.1 (md:62)
    Stage-2 spec; eliminates the white-gap artefact in figure
    renders at the cost of degraded accuracy past `|t|=1`.
    Opt-in to preserve the default-on Rule 1 contract; figure
    scripts pass `extrapolate=true` to get filled panels.
  - `max_steps_per_target` — cap per-target Stage-1 walk length.
  - `rng_seed::Integer`  — for deterministic target shuffle (tests).
  - `enforce_real_axis_symmetry::Bool` — opt-in.  When `true`, walk only
    upper-half + on-axis targets (`Im(z) ≥ 0`) and populate the lower
    half via Schwarz reflection `u(z̄) = ū(z)`.  Correct iff the ODE
    coefficients are real *and* the IC sits on the real axis (`Im(z₀)
    = Im(u₀) = Im(u'₀) = 0`); we validate the latter and trust the
    caller on the former.  Eliminates the `shuffle(rng, targets)`-
    induced visited-tree asymmetry that costs 4-5 orders of magnitude
    on conjugate-pair cells far from the IC (worklog 014 §"Bug 1",
    bead `padetaylor-dtj`); roughly halves wall time on full-plane
    fills.  Default `false` preserves the FW 2011 algorithm verbatim.
  - `verbose::Bool` — when `true`, emit eager-flushed progress lines
    to `stdout`: one line per target start, one line every
    `progress_every` aggregate inner steps with elapsed wall, current
    `z_cur`, `|u_cur|`, distance to target, and steps/second.  Default
    `false` (silent).
  - `progress_every::Integer` — Stage-1-inner-step granularity for the
    `verbose` progress lines.  Defaults to 500.

## Threading

`_wedge_evaluations` runs the five wedge-direction candidate steps in
parallel via `Threads.@threads`.  Speedup is bounded by the wedge
count (5 candidates → up to ~5× per step) plus the canonical-Padé
build that follows selection (still single-threaded).  Realistic
end-to-end speedup with 5+ threads is ~3–4× on `BigFloat`-256 walks
where SVD dominates per-step cost.

`f` (the user RHS) must therefore be **thread-safe** — pure functions
are fine; closures over mutable shared state are not.  Set
`JULIA_NUM_THREADS=N` before launching Julia (or pass `-tN`); the
serial path is recovered automatically when `Threads.nthreads() == 1`.

Throws if Stage 1 cannot reach a target within `max_steps_per_target`
steps, if all 5 wedge candidates fail, or if `:adaptive_ffw` is
requested (Tier-4 deferral per ADR-0004).
"""
function path_network_solve(prob::PadeTaylorProblem,
                            grid::AbstractVector{<:Complex};
                            h::Real      = 0.5,
                            order::Integer = prob.order,
                            wedge_angles::AbstractVector{<:Real} = DEFAULT_WEDGE,
                            step_selection::Symbol  = :min_u,
                            step_size_policy::Symbol = :fixed,
                            adaptive_tol::Real = 1.0e-12,
                            k_conservative::Real = 1.0e-3,
                            max_rescales::Integer = 50,
                            node_separation::Union{Nothing, Function} = nothing,
                            branch_points::Tuple = (),
                            branch_cut_angles = π,
                            cross_branch::Bool = false,
                            initial_sheet::AbstractVector{<:Integer} =
                                zeros(Int, length(branch_points)),
                            grid_sheet::Union{Nothing, AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
                            extrapolate::Bool = false,
                            max_steps_per_target::Integer = 1000,
                            rng_seed::Integer = 0,
                            enforce_real_axis_symmetry::Bool = false,
                            verbose::Bool = false,
                            progress_every::Integer = 500)

    if enforce_real_axis_symmetry
        isempty(branch_points) || throw(ArgumentError(
            "path_network_solve: enforce_real_axis_symmetry=true combined " *
            "with non-empty branch_points is not supported (the Schwarz " *
            "mirror step assumes a simply connected upper half plane).  " *
            "Set one or the other.  See ADR-0013 §\"Open follow-ups\"."))
        return _solve_with_schwarz_reflection(
            prob, grid; h, order, wedge_angles, step_selection,
            step_size_policy, adaptive_tol, k_conservative, max_rescales,
            node_separation,
            max_steps_per_target, rng_seed, verbose, progress_every)
    end

    step_size_policy ∈ (:fixed, :adaptive_ffw) || throw(ArgumentError(
        "path_network_solve: step_size_policy must be :fixed or " *
        ":adaptive_ffw (got :$step_size_policy).  See ADR-0011 for " *
        ":adaptive_ffw semantics."))
    step_size_policy === :fixed || adaptive_tol > 0 || throw(ArgumentError(
        "path_network_solve: adaptive_tol must be positive (got $adaptive_tol)."))
    step_selection in (:min_u, :steepest_descent) || throw(ArgumentError(
        "path_network_solve: step_selection must be :min_u or " *
        ":steepest_descent (got $step_selection)."))
    length(wedge_angles) == 5 || throw(ArgumentError(
        "path_network_solve: wedge_angles must have exactly 5 entries " *
        "(got $(length(wedge_angles))).  FW 2011 §3.1 fixes this."))
    h > 0 || throw(ArgumentError(
        "path_network_solve: h must be positive real (got $h)."))

    # Branch / sheet input validation (ADR-0013).
    length(initial_sheet) == length(branch_points) || throw(ArgumentError(
        "path_network_solve: initial_sheet length $(length(initial_sheet)) " *
        "must equal branch_points length $(length(branch_points)).  " *
        "Default initial_sheet = zeros(Int, length(branch_points))."))
    cut_angles_resolved = resolve_cut_angles(branch_points, branch_cut_angles)
    sheet_root = collect(Int, initial_sheet)
    branched   = !isempty(branch_points)
    if !branched && cross_branch
        throw(ArgumentError(
            "path_network_solve: cross_branch=true requires non-empty " *
            "branch_points; got branch_points = ()."))
    end
    # A5 / bead padetaylor-hed: validate grid_sheet shape.  Done here so
    # malformed input throws BEFORE the Stage-1 walk (which can be
    # minutes-long).
    if grid_sheet !== nothing
        length(grid_sheet) == length(grid) || throw(ArgumentError(
            "path_network_solve: grid_sheet length $(length(grid_sheet)) " *
            "must equal grid length $(length(grid))."))
        for (i, s) in pairs(grid_sheet)
            length(s) == length(branch_points) || throw(ArgumentError(
                "path_network_solve: grid_sheet[$i] length $(length(s)) " *
                "must equal branch_points length $(length(branch_points))."))
        end
    end

    # Promote element type to Complex{T <: AbstractFloat}.
    T     = float(typeof(real(first(grid))))
    CT    = Complex{T}
    h_T   = T(h)

    # Initialise visited at IC.
    z_0   = CT(prob.zspan[1])
    u_0   = CT(prob.y0[1])
    up_0  = CT(prob.y0[2])
    pade_0 = _local_pade(prob.f, z_0, u_0, up_0, order, h_T)

    visited_z      = [z_0]
    visited_u      = [u_0]
    visited_up     = [up_0]
    visited_pade   = [pade_0]
    visited_h      = [h_T]
    visited_parent = [0]               # root has no parent
    visited_sheet  = branched ? Vector{Int}[copy(sheet_root)] :
                                 Vector{Int}[Int[]]

    # Stage 1: build path tree.  Random target shuffle (FW 2011 line 156).
    rng     = MersenneTwister(rng_seed)
    targets = shuffle(rng, collect(CT, grid))

    if verbose
        println("path_network_solve: Stage 1 begin — ",
                length(targets), " targets, h=", h_T, ", order=", order,
                ", threads=", nthreads())
        flush(stdout)
    end

    t_start_stage1 = time()
    total_n_steps  = 0
    # Per-target initial step magnitude.  Under :fixed this stays at
    # the constant `h_T`; under :adaptive_ffw FFW md:93 prescribes
    # "the initial step length is always the scaled step length stored
    # at the current point" — i.e. the previous accepted step's |q·h|,
    # threaded across the walk so the controller has memory.
    for (target_idx, target) in enumerate(targets)
        any(z -> isapprox(z, target; atol=10*eps(T)), visited_z) && continue

        idx_v = _nearest_visited(visited_z, target)
        z_cur, u_cur, up_cur = visited_z[idx_v], visited_u[idx_v], visited_up[idx_v]
        # Seed `h_cur` from the visited node's stored step length
        # (which is the FFW "scaled step length stored at the current
        # point" under adaptive).  Under :fixed all visited_h are h_T.
        h_cur = visited_h[idx_v]

        verbose && _verbose_target_start(target_idx, length(targets),
                                         target, z_cur, t_start_stage1)

        # The first node landed in this walk chains off the nearest
        # already-visited node (FW 2011 line 164); each later node of
        # the walk chains off its immediate predecessor.
        parent_idx = idx_v
        n_steps = 0
        # The walk-termination distance.  Under fixed/adaptive without
        # R it's `h_cur` (the prior accepted step length).  When R is
        # supplied we use R(z_cur) — the spatial density target — for
        # the same reason FFW prescribe R as the wedge-step magnitude:
        # the walker should approach the target with the local
        # density's resolution, not the controller's last memory.
        term_dist = node_separation === nothing ? h_cur :
                     _eval_node_separation(node_separation, z_cur, T)
        while abs(z_cur - target) > term_dist
            n_steps        += 1
            total_n_steps  += 1
            n_steps > max_steps_per_target && throw(ErrorException(
                "path_network_solve: target $target unreachable in " *
                "$max_steps_per_target steps from $(visited_z[idx_v]); " *
                "current z = $z_cur, |z - target| = $(abs(z_cur - target))."))

            goal_dir = angle(target - z_cur)

            # FFW md:67-72 / md:97 non-uniform Stage-1 nodes: when
            # `node_separation` is supplied, the per-step initial
            # magnitude is `R(z_cur)` regardless of `step_size_policy`.
            # Under :adaptive_ffw the controller may shrink R(z_cur)
            # but cannot grow it (q ≤ 1).  Under :fixed R(z_cur)
            # simply sets the magnitude — no rescale.
            h_seed = node_separation === nothing ? h_cur :
                      _eval_node_separation(node_separation, z_cur, T)

            # Adaptive rescale loop (FFW md:93): re-run the wedge,
            # re-select min-|u|, re-compute T(h); if T(h) > Tol shrink
            # h and try again.  Under :fixed this loop runs exactly
            # once with h_step = h_seed and skips the T(h) check.
            h_step    = h_seed
            n_resc    = 0
            z_new, u_new, up_new, pade_sel = z_cur, CT(NaN, NaN), CT(NaN, NaN), nothing
            while true
                evals = _wedge_evaluations(prob.f, z_cur, u_cur, up_cur,
                                           order, h_step, goal_dir, wedge_angles)
                # ADR-0013 refuse mode: nullify any candidate whose step
                # crosses a branch cut.  Stays a no-op in the default
                # (no-branches) case and in cross_branch=true.
                if branched && !cross_branch
                    evals = _filter_forbidden_candidates(
                        evals, z_cur, branch_points, cut_angles_resolved, CT)
                end
                idx_sel = _select_candidate(step_selection, evals, u_cur, up_cur,
                                            goal_dir, wedge_angles)
                z_new, u_new, up_new, pade_sel = evals[idx_sel]
                if !isfinite(abs(u_new)) || pade_sel === nothing
                    msg = "path_network_solve: all 5 wedge candidates failed at " *
                          "z=$z_cur (target=$target); shrink h or widen wedge."
                    if branched && !cross_branch
                        msg *= "  (Refuse mode active with " *
                               "$(length(branch_points)) branch_points: all " *
                               "candidates may be forbidden by the cut " *
                               "constraint — consider cross_branch=true or " *
                               "redirect branch_cut_angles.)"
                    end
                    throw(ErrorException(msg))
                end

                step_size_policy === :fixed && break

                # :adaptive_ffw — estimate T(h) at the selected wedge
                # direction.  The wedge angle θ_sel rotates the real
                # h_step magnitude; ffw_truncation_error consumes the
                # COMPLEX h that was actually taken.
                θ_sel = T(goal_dir) + T(wedge_angles[idx_sel])
                h_complex = Complex{T}(h_step * cos(θ_sel), h_step * sin(θ_sel))
                Th = ffw_truncation_error(prob.f, z_cur, u_cur, up_cur,
                                          Int(order), h_complex)
                Th ≤ adaptive_tol && break
                n_resc ≥ max_rescales && throw(ErrorException(
                    "path_network_solve: :adaptive_ffw max_rescales=" *
                    "$max_rescales exhausted at z=$z_cur, target=$target; " *
                    "T(h)=$Th > Tol=$adaptive_tol.  Suggestion: tighten " *
                    "max_rescales, relax adaptive_tol, or shorten initial h."))
                q = ffw_rescale_q(adaptive_tol, Th, Int(order); k = k_conservative)
                h_step *= q
                n_resc += 1
            end

            # Canonical (real-h) Padé centered at z_new, per ADR-0004's
            # design decision: Stage 2 always interpolates inside |t| ≤ 1
            # of a REAL-direction disc.  The wedge-direction `pade_sel`
            # above is consumed by step selection only; storing it as the
            # visited node's Padé would invalidate the `t = (z_f - z_v)/h`
            # Stage-2 evaluation (which assumes real h).  Cost: one extra
            # Coefficients + RobustPade per visited node.  Under
            # :adaptive_ffw the canonical Padé uses the accepted
            # h_step (FFW md:93: scaled length is stored at the point).
            pade_canonical = _local_pade(prob.f, z_new, u_new, up_new, order, h_step)
            # Update the per-walk seed for the next step.  Under :fixed
            # this stays at h_T (h_step == h_T == h_cur).  Under
            # :adaptive_ffw it carries the controller's accepted step
            # length forward.  When `node_separation` is supplied,
            # h_cur is irrelevant: the next iteration re-seeds h_seed
            # from R(z_new).
            h_cur = h_step

            push!(visited_z, z_new)
            push!(visited_u, u_new)
            push!(visited_up, up_new)
            push!(visited_pade, pade_canonical)
            push!(visited_h, h_step)
            push!(visited_parent, parent_idx)
            # Sheet bookkeeping (ADR-0013).  Refuse mode: copies parent's
            # sheet (no crossing happened by construction).  Cross mode:
            # checks each branch cut for crossing and bumps per
            # sign(winding_delta).  Default no-branches: pushes Int[].
            if branched
                push!(visited_sheet,
                      cross_branch ?
                          step_sheet_update(visited_sheet[parent_idx],
                                            z_cur, z_new,
                                            branch_points, cut_angles_resolved) :
                          copy(visited_sheet[parent_idx]))
            else
                push!(visited_sheet, Int[])
            end
            parent_idx = length(visited_z)   # next step chains off this node

            z_cur, u_cur, up_cur = z_new, u_new, up_new

            # Re-evaluate the termination distance at the new position.
            # Under R(ζ) this can shrink as Re ζ grows; without R it
            # stays at the controller-memory `h_cur` (which equals h_T
            # under :fixed).
            term_dist = node_separation === nothing ? h_cur :
                         _eval_node_separation(node_separation, z_cur, T)

            if verbose && total_n_steps % progress_every == 0
                _verbose_step(total_n_steps, target_idx, length(targets),
                              z_cur, u_cur, target, t_start_stage1)
            end
        end
    end

    if verbose
        elapsed = time() - t_start_stage1
        @printf(stdout, "path_network_solve: Stage 1 done — %d visited nodes, %d total inner steps, %.2f s (%.1f steps/s)\n",
                length(visited_z), total_n_steps, elapsed,
                total_n_steps / max(elapsed, eps()))
        flush(stdout)
    end

    # Stage 2: fine-grid extrapolation.  A5 / bead padetaylor-hed adds
    # sheet-aware nearest-visited when `grid_sheet` is supplied: each
    # grid point's pool is restricted to visited nodes whose
    # `visited_sheet[k] == grid_sheet[i]`.  When `grid_sheet === nothing`
    # the lookup is pure Euclidean nearest-visited (legacy behaviour).
    grid_z  = collect(CT, grid)
    grid_u  = Vector{CT}(undef, length(grid_z))
    grid_up = Vector{CT}(undef, length(grid_z))
    nan_CT  = CT(NaN, NaN)
    for (i, z_f) in enumerate(grid_z)
        idx_v = grid_sheet === nothing ?
                    _nearest_visited(visited_z, z_f) :
                    _nearest_visited_on_sheet(visited_z, visited_sheet,
                                              z_f, grid_sheet[i])
        if idx_v == 0
            # No matching-sheet visited node exists.
            grid_u[i]  = nan_CT
            grid_up[i] = nan_CT
            continue
        end
        z_v   = visited_z[idx_v]
        h_v   = visited_h[idx_v]
        # ADR-0015: extrapolate=true skips the disc-radius check
        # to match FFW md:62's Stage-2 spec (evaluate Padé at every
        # fine-grid cell regardless of |t| vs 1).
        if !extrapolate && abs(z_f - z_v) > h_v
            grid_u[i]  = nan_CT
            grid_up[i] = nan_CT
        else
            t = (z_f - z_v) / h_v
            grid_u[i]  = _evaluate_pade(visited_pade[idx_v], t)
            grid_up[i] = _evaluate_pade_deriv(visited_pade[idx_v], t) / h_v
        end
    end

    return PathNetworkSolution{T}(visited_z, visited_u, visited_up,
                                  visited_pade, visited_h, visited_parent,
                                  grid_z, grid_u, grid_up, visited_sheet)
end

# ADR-0013 refuse-mode filter: replace each forbidden candidate with the
# same shape as a Padé failure (Inf-|u| sentinel, nothing-Padé), so the
# downstream `_select_candidate` naturally skips it (min_u of Inf is
# large; steepest_descent's pade-nothing check follows the same path
# the Padé-failure already does in the existing "all 5 failed" guard).
function _filter_forbidden_candidates(evals, z_cur, branches, cut_angles, ::Type{CT}) where CT
    return map(evals) do e
        z_n, _, _, pade = e
        if pade !== nothing && any_cut_crossed(z_cur, z_n, branches, cut_angles)
            (z_cur, CT(Inf, 0), CT(0, 0), nothing)
        else
            e
        end
    end
end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Schwarz-reflection driver (worklog 014 §"Bug 1", bead `padetaylor-dtj`).
#
# For ODEs that preserve complex conjugation (real-coefficient `f`, real
# IC on the real axis), the analytic solution satisfies `u(z̄) = ū(z)`
# globally.  The default `shuffle(rng, targets)` step breaks this
# numerically — the visited-tree order is shuffle-dependent, so the
# Stage-2 nearest-visited lookup at `z` and `conj(z)` may land on
# non-conjugate visited nodes, blowing `|u|` apart by 4-5 orders of
# magnitude at far-from-IC conjugate pairs.
#
# Cure: walk the unique upper-half-or-on-axis representatives only, then
# build the requested grid by indexing the upper-walk's outputs (mirror
# via `conj` for input cells with `Im(z) < 0`).  Bit-exact symmetric;
# roughly half the wall-time of a full-plane walk.
#
# Correctness preconditions: the caller's `f` has real coefficients
# (cannot be checked at this layer; opt-in kwarg pushes the burden onto
# the caller) AND the IC `(z₀, u₀, u'₀)` lies on the real axis (we DO
# check this — fail loud per CLAUDE.md Rule 1).
function _solve_with_schwarz_reflection(prob::PadeTaylorProblem,
                                         grid::AbstractVector{<:Complex};
                                         h::Real,
                                         order::Integer,
                                         wedge_angles::AbstractVector{<:Real},
                                         step_selection::Symbol,
                                         step_size_policy::Symbol,
                                         adaptive_tol::Real,
                                         k_conservative::Real,
                                         max_rescales::Integer,
                                         node_separation::Union{Nothing, Function} = nothing,
                                         max_steps_per_target::Integer = 1000,
                                         rng_seed::Integer = 0,
                                         verbose::Bool = false,
                                         progress_every::Integer = 500)
    T  = float(typeof(real(first(grid))))
    CT = Complex{T}
    tol = 10 * eps(T)

    abs(imag(prob.zspan[1])) <= tol || throw(ArgumentError(
        "path_network_solve: enforce_real_axis_symmetry=true requires " *
        "imag(zspan[1]) ≈ 0 (got $(prob.zspan[1])).  Schwarz reflection " *
        "u(z̄) = ū(z) is exact only when the IC sits on the real axis. " *
        "Suggestion: shift the IC onto the real axis, or set the kwarg false."))

    u0, up0 = prob.y0[1], prob.y0[2]
    abs(imag(u0)) <= tol || throw(ArgumentError(
        "path_network_solve: enforce_real_axis_symmetry=true requires " *
        "imag(y0[1]) ≈ 0 (got $u0).  Suggestion: real ICs only, " *
        "or set the kwarg false."))
    abs(imag(up0)) <= tol || throw(ArgumentError(
        "path_network_solve: enforce_real_axis_symmetry=true requires " *
        "imag(y0[2]) ≈ 0 (got $up0).  Suggestion: real ICs only, " *
        "or set the kwarg false."))

    # Map each input cell to its upper-half-or-on-axis canonical
    # representative.  `complex(real(z), abs(imag(z)))` collapses the
    # ±0.0im signed-zero ambiguity to +0.0im for on-axis cells, so the
    # `idx_of` dict lookup below is deterministic.
    grid_input  = collect(CT, grid)
    upper_canon = CT[complex(real(z), abs(imag(z))) for z in grid_input]
    upper_unique = unique(upper_canon)

    sol_upper = path_network_solve(
        prob, upper_unique;
        h = h, order = order, wedge_angles = wedge_angles,
        step_selection = step_selection,
        step_size_policy = step_size_policy,
        adaptive_tol = adaptive_tol,
        k_conservative = k_conservative,
        max_rescales = max_rescales,
        node_separation = node_separation,
        max_steps_per_target = max_steps_per_target,
        rng_seed = rng_seed,
        enforce_real_axis_symmetry = false,
        verbose = verbose, progress_every = progress_every)

    idx_of = Dict{CT, Int}()
    for (i, z) in enumerate(sol_upper.grid_z)
        idx_of[z] = i
    end

    n = length(grid_input)
    grid_u  = Vector{CT}(undef, n)
    grid_up = Vector{CT}(undef, n)
    for (i, z) in enumerate(grid_input)
        j = idx_of[upper_canon[i]]
        if imag(z) < 0
            grid_u[i]  = conj(sol_upper.grid_u[j])
            grid_up[i] = conj(sol_upper.grid_up[j])
        else
            grid_u[i]  = sol_upper.grid_u[j]
            grid_up[i] = sol_upper.grid_up[j]
        end
    end

    # The visited tree (and hence `visited_parent`) is the upper-half
    # walk's verbatim; only the Stage-2 grid outputs are mirrored.
    return PathNetworkSolution{T}(
        sol_upper.visited_z, sol_upper.visited_u, sol_upper.visited_up,
        sol_upper.visited_pade, sol_upper.visited_h, sol_upper.visited_parent,
        grid_input, grid_u, grid_up)
end

# Compute the canonical local Padé centered at (z, u, up) with real step h.
# We clone the state so the IC isn't mutated.
function _local_pade(f, z::Complex{T}, u::Complex{T}, up::Complex{T},
                     order::Integer, h::T) where T
    s = PadeStepperState{Complex{T}}(z, u, up)
    _, pade = pade_step_with_pade!(s, f, Int(order), h)
    return pade
end

# Evaluate the user-supplied `node_separation` function `R(z)` and
# validate its output.  FFW md:72 prescribes a positive real step
# magnitude; non-finite or non-positive returns are fail-loud per
# CLAUDE.md Rule 1.  Result is coerced to `T` so the downstream
# step-magnitude arithmetic stays in the working float type.
function _eval_node_separation(R, z::Complex{T}, ::Type{T}) where T
    r = R(z)
    isfinite(r) || throw(ArgumentError(
        "path_network_solve: node_separation R($z) returned a non-finite " *
        "value ($r).  R must return a positive finite real (FFW 2017 §2.1.2 " *
        "md:72 prescribes R(ζ) ∈ ℝ_{>0})."))
    r > 0 || throw(ArgumentError(
        "path_network_solve: node_separation R($z) returned $r ≤ 0.  R " *
        "must return a positive finite real (FFW 2017 §2.1.2 md:72)."))
    return T(r)
end

# Nearest visited node to `target` by Euclidean distance.  Lexicographic
# tiebreak on (Re, Im) for reproducibility (spec gap 1, unified spec §12).
function _nearest_visited(visited_z::Vector{Complex{T}}, target::Complex{T}) where T
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

# A5 / ADR-0013 §"Open follow-ups" / bead padetaylor-hed.
# Same min-distance scan as `_nearest_visited`, restricted to visited
# nodes whose `visited_sheet[k]` matches `sheet`.  Returns `0` if no
# matching-sheet visited node exists — the caller treats that as
# fail-soft (NaN output in Stage 2), same as the existing
# "too-far-away" case (`abs(z_f - z_v) > visited_h[k]`).
#
# Lexicographic tie-break (real first, then imag) mirrors
# `_nearest_visited` so results are deterministic.
function _nearest_visited_on_sheet(visited_z::Vector{Complex{T}},
                                    visited_sheet::Vector{Vector{Int}},
                                    target::Complex{T},
                                    sheet::AbstractVector{<:Integer}) where T
    idx, best = 0, T(Inf)
    @inbounds for i in eachindex(visited_z)
        visited_sheet[i] == sheet || continue
        d = abs(visited_z[i] - target)
        if idx == 0 || d < best ||
           (d == best &&
            (real(visited_z[i]) < real(visited_z[idx]) ||
             (real(visited_z[i]) == real(visited_z[idx]) &&
              imag(visited_z[i]) < imag(visited_z[idx]))))
            idx, best = i, d
        end
    end
    return idx
end

# Try each of 5 wedge directions; return Vector of (z_new, u_new, up_new, pade).
# A failed candidate gets (z_cur, Inf+0im, 0+0im, nothing) so it loses min-|u|.
#
# Threading: each wedge index `k` does an independent PadeStepperState +
# pade_step_with_pade! call.  `evals[k]` is written by exactly one thread
# (its own `k`), and `f`, `z_cur`, `u_cur`, `up_cur`, `h`, `goal_dir`,
# `wedge_angles` are all read-only.  `pade_step_with_pade!` is pure
# w.r.t. `f` (the user RHS must itself be thread-safe; closures over
# mutable shared state break this contract — see `path_network_solve`
# docstring).  Result ordering is by `k`, not by thread completion
# order, so `argmin` downstream is deterministic.
function _wedge_evaluations(f, z_cur::Complex{T}, u_cur::Complex{T},
                            up_cur::Complex{T}, order::Integer, h::T,
                            goal_dir::T,
                            wedge_angles::AbstractVector{<:Real}) where T
    CT = Complex{T}
    n_w = length(wedge_angles)
    evals = Vector{Tuple{CT, CT, CT, Union{Nothing, PadeApproximant{CT}}}}(undef, n_w)
    @threads for k in 1:n_w
        θ = wedge_angles[k]
        h_step = CT(h * cos(T(goal_dir) + T(θ)),
                    h * sin(T(goal_dir) + T(θ)))
        s = PadeStepperState{CT}(z_cur, u_cur, up_cur)
        try
            _, pade = pade_step_with_pade!(s, f, Int(order), h_step)
            evals[k] = (s.z, s.u, s.up, pade)
        catch
            evals[k] = (z_cur, CT(Inf, 0), CT(0, 0), nothing)
        end
    end
    return evals
end

# Verbose-mode progress helpers — Float64-truncate BigFloat coordinates
# for human-readable single-line output.  Eager-flushed so progress
# remains visible under stdout buffering.
function _verbose_target_start(target_idx::Int, n_targets::Int,
                               target::Complex{T}, z_start::Complex{T},
                               t_start::Float64) where T
    elapsed = time() - t_start
    println(@sprintf("[%7.1fs] target %4d/%-4d  z=%s  start_z=%s",
                     elapsed, target_idx, n_targets,
                     _fmt(target), _fmt(z_start)))
    flush(stdout)
end

function _verbose_step(total_n_steps::Int, target_idx::Int, n_targets::Int,
                       z_cur::Complex{T}, u_cur::Complex{T},
                       target::Complex{T}, t_start::Float64) where T
    elapsed = time() - t_start
    sps     = total_n_steps / max(elapsed, eps())
    println(@sprintf("[%7.1fs]   step %7d  tgt %d/%d  z=%s  |u|=%.3e  |Δ|=%.3e  %.1f steps/s",
                     elapsed, total_n_steps, target_idx, n_targets,
                     _fmt(z_cur), Float64(abs(u_cur)),
                     Float64(abs(z_cur - target)), sps))
    flush(stdout)
end

# Truncate Complex{<:AbstractFloat} to a 5-digit fixed-format string —
# BigFloats would print ~75 digits per coordinate otherwise.
_fmt(z::Complex) = @sprintf("(%+.4e,%+.4e)", Float64(real(z)), Float64(imag(z)))

# Min-|u| (FW 2011 default) or steepest-descent (FW 2011 §5.4.1) selection.
function _select_candidate(selection::Symbol, evals,
                           u_cur::Complex{T}, up_cur::Complex{T},
                           goal_dir::T,
                           wedge_angles::AbstractVector{<:Real}) where T
    if selection === :min_u
        return argmin(abs(e[2]) for e in evals)
    else
        # :steepest_descent: θ_sd = arg(-u/u'); clip to wedge.
        θ_sd = abs(up_cur) > 0 ? angle(-u_cur / up_cur) : T(goal_dir)
        offsets = (T(goal_dir) + T(θ) for θ in wedge_angles)
        return argmin(abs(θ_sd - off) for off in offsets)
    end
end

# -----------------------------------------------------------------------------
# A5 public accessor: post-hoc per-point sheet-aware evaluation
# (bead padetaylor-hed, complement to the grid_sheet kwarg).
# -----------------------------------------------------------------------------

"""
    eval_at_sheet(sol::PathNetworkSolution, z::Complex, sheet::AbstractVector{<:Integer}) -> (u, up)

Evaluate the stored solution at a single point `z`, restricting the
nearest-visited lookup to visited nodes whose `visited_sheet[k]`
equals `sheet`.  Returns `(NaN+NaN·im, NaN+NaN·im)` if no
matching-sheet visited node lies within `visited_h[k]` of `z`
(fail-soft, same as Stage-2's too-far-away case).

Use this for post-hoc per-point queries: figure scripts rendering a
multi-sheet PVI heatmap call `eval_at_sheet(sol, z, [k_branch_1,
k_branch_2])` for each pixel.  When all queries are vectorisable
into a `grid` ahead of time, prefer the `grid_sheet` kwarg of
`path_network_solve` instead — it amortises the Stage-1 build cost
across all query points.

Throws `ArgumentError` if `length(sheet) != length(sol.visited_sheet[1])`
(only checked when there is at least one visited node with a non-empty
sheet tuple).
"""
function eval_at_sheet(sol::PathNetworkSolution{T},
                       z::Complex,
                       sheet::AbstractVector{<:Integer};
                       extrapolate::Bool = false) where T
    CT     = Complex{T}
    sheet_branched = !isempty(sol.visited_sheet) &&
                      !isempty(sol.visited_sheet[1])
    if sheet_branched
        length(sheet) == length(sol.visited_sheet[1]) || throw(ArgumentError(
            "eval_at_sheet: sheet length $(length(sheet)) must equal " *
            "visited_sheet inner length $(length(sol.visited_sheet[1]))."))
    elseif !isempty(sheet)
        throw(ArgumentError(
            "eval_at_sheet: solution has no branch points " *
            "(visited_sheet entries are empty), but caller passed " *
            "sheet of length $(length(sheet))."))
    end
    z_CT   = CT(z)
    sheet_typed = collect(Int, sheet)
    idx_v  = _nearest_visited_on_sheet(sol.visited_z, sol.visited_sheet,
                                        z_CT, sheet_typed)
    nan_CT = CT(NaN, NaN)
    idx_v == 0 && return (nan_CT, nan_CT)
    z_v    = sol.visited_z[idx_v]
    h_v    = sol.visited_h[idx_v]
    # ADR-0015: extrapolate=true skips disc-radius check (FFW md:62).
    # ADR-0015: extrapolate=true skips disc-radius check (FFW md:62).
    if !extrapolate && abs(z_CT - z_v) > h_v
        return (nan_CT, nan_CT)
    end
    t = (z_CT - z_v) / h_v
    u  = _evaluate_pade(sol.visited_pade[idx_v], t)
    up = _evaluate_pade_deriv(sol.visited_pade[idx_v], t) / h_v
    return (u, up)
end

# -----------------------------------------------------------------------------
# ADR-0015 sheet-blind public accessor: per-point Stage-2 eval without
# the A5 sheet restriction.  Replaces the local `stage2_eval_blind`
# helpers figure scripts have historically rolled themselves.
# -----------------------------------------------------------------------------

"""
    eval_at(sol::PathNetworkSolution, z::Complex; extrapolate=false) -> (u, up)

Evaluate the stored solution at a single point `z` using the nearest
visited node's local Padé approximant.  Sheet-blind: the lookup
ignores `visited_sheet`.  For sheet-restricted queries on multi-sheet
solutions, use `eval_at_sheet`.

When `extrapolate = false` (the default), returns `(NaN+NaN·im,
NaN+NaN·im)` if `z` lies more than `visited_h[idx_v]` from the
nearest visited node — fail-soft per ADR-0004 / CLAUDE.md Rule 1.

When `extrapolate = true`, the disc-radius check is skipped and the
nearest-node Padé is evaluated regardless; matches FFW 2017 §2.1.1
(md:62) Stage-2 behaviour.  Used by figure scripts that prefer
filled rendering panels over fail-soft NaN gaps.  See ADR-0015 for
the rationale.
"""
function eval_at(sol::PathNetworkSolution{T},
                 z::Complex;
                 extrapolate::Bool = false) where T
    CT     = Complex{T}
    z_CT   = CT(z)
    idx_v  = _nearest_visited(sol.visited_z, z_CT)
    z_v    = sol.visited_z[idx_v]
    h_v    = sol.visited_h[idx_v]
    nan_CT = CT(NaN, NaN)
    if !extrapolate && abs(z_CT - z_v) > h_v
        return (nan_CT, nan_CT)
    end
    t = (z_CT - z_v) / h_v
    u  = _evaluate_pade(sol.visited_pade[idx_v], t)
    up = _evaluate_pade_deriv(sol.visited_pade[idx_v], t) / h_v
    return (u, up)
end

# TIER-3 INTERFACE: BVP dispatcher (bead padetaylor-804) consumes
# PathNetworkSolution.{grid_z, grid_u, grid_up}.  No code here until
# Phase 11; see docs/adr/0004-path-network-architecture.md.

end # module PathNetwork
