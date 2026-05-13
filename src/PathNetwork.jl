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

## Tier-3-to-5 deferrals (per ADR-0004)

  - **BVP composition (Tier 3)**: stub comment at end of this file.
  - **CoordTransforms (Tier 4)**: external wrapper, no change here.
  - **SheetTracker (Tier 5)**: external wrapper, no change here.
  - **`:adaptive_ffw` step-size policy**: throws `ArgumentError` per
    fail-fast discipline; Tier-4 work item.
"""
module PathNetwork

using Random:        shuffle, MersenneTwister
using ..RobustPade:  PadeApproximant
using ..PadeStepper: PadeStepperState, pade_step_with_pade!,
                     _evaluate_pade, _evaluate_pade_deriv
using ..Problems:    PadeTaylorProblem

export PathNetworkSolution, path_network_solve

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

A `NaN + NaN·im` entry in `grid_u` or `grid_up` signals that the
corresponding `grid` point was not within `h` of any visited node;
fail-loud per CLAUDE.md Rule 1.
"""
struct PathNetworkSolution{T <: AbstractFloat}
    visited_z    :: Vector{Complex{T}}
    visited_u    :: Vector{Complex{T}}
    visited_up   :: Vector{Complex{T}}
    visited_pade :: Vector{PadeApproximant{Complex{T}}}
    visited_h    :: Vector{T}                # real (canonical) per node
    grid_z       :: Vector{Complex{T}}
    grid_u       :: Vector{Complex{T}}
    grid_up      :: Vector{Complex{T}}
end

# -----------------------------------------------------------------------------
# Public driver
# -----------------------------------------------------------------------------

"""
    path_network_solve(prob, grid; h=0.5, order=30,
                       wedge_angles=DEFAULT_WEDGE,
                       step_selection=:min_u,
                       step_size_policy=:fixed,
                       max_steps_per_target=1000,
                       rng_seed=0) -> PathNetworkSolution

Build the FW 2011 §3.1 path-network covering `grid::AbstractVector{<:Complex}`
from `prob`'s IC, then return the Stage-1 tree + Stage-2 evaluations.

Required: `prob isa PadeTaylorProblem` with 2nd-order `y0 = (u_0, up_0)`.
The IC point `prob.zspan[1]` becomes `visited_z[1]`.

Kwargs:
  - `h::Real`            — canonical step length (FW default 0.5).
  - `order::Integer`     — Taylor truncation degree (FW default 30).
  - `wedge_angles`       — five angle offsets relative to goal direction.
  - `step_selection`     — `:min_u` (default; FW 2011) or `:steepest_descent`.
  - `step_size_policy`   — `:fixed` (Tier-2 only); `:adaptive_ffw` throws.
  - `max_steps_per_target` — cap per-target Stage-1 walk length.
  - `rng_seed::Integer`  — for deterministic target shuffle (tests).

Throws if Stage 1 cannot reach a target within `max_steps_per_target`
steps, if all 5 wedge candidates fail, or if `:adaptive_ffw` is
requested (Tier-4 deferral per ADR-0004).
"""
function path_network_solve(prob::PadeTaylorProblem,
                            grid::AbstractVector{<:Complex};
                            h::Real      = 0.5,
                            order::Integer = 30,
                            wedge_angles::AbstractVector{<:Real} = DEFAULT_WEDGE,
                            step_selection::Symbol  = :min_u,
                            step_size_policy::Symbol = :fixed,
                            max_steps_per_target::Integer = 1000,
                            rng_seed::Integer = 0)

    step_size_policy === :fixed || throw(ArgumentError(
        "path_network_solve: :adaptive_ffw step-size policy is a Tier-4 " *
        "deferral (FFW 2017 §2.4, bead padetaylor-???); only :fixed is " *
        "supported in Tier-2.  See ADR-0004."))
    step_selection in (:min_u, :steepest_descent) || throw(ArgumentError(
        "path_network_solve: step_selection must be :min_u or " *
        ":steepest_descent (got $step_selection)."))
    length(wedge_angles) == 5 || throw(ArgumentError(
        "path_network_solve: wedge_angles must have exactly 5 entries " *
        "(got $(length(wedge_angles))).  FW 2011 §3.1 fixes this."))
    h > 0 || throw(ArgumentError(
        "path_network_solve: h must be positive real (got $h)."))

    # Promote element type to Complex{T <: AbstractFloat}.
    T     = float(typeof(real(first(grid))))
    CT    = Complex{T}
    h_T   = T(h)

    # Initialise visited at IC.
    z_0   = CT(prob.zspan[1])
    u_0   = CT(prob.y0[1])
    up_0  = CT(prob.y0[2])
    pade_0 = _local_pade(prob.f, z_0, u_0, up_0, order, h_T)

    visited_z    = [z_0]
    visited_u    = [u_0]
    visited_up   = [up_0]
    visited_pade = [pade_0]
    visited_h    = [h_T]

    # Stage 1: build path tree.  Random target shuffle (FW 2011 line 156).
    rng     = MersenneTwister(rng_seed)
    targets = shuffle(rng, collect(CT, grid))

    for target in targets
        any(z -> isapprox(z, target; atol=10*eps(T)), visited_z) && continue

        idx_v = _nearest_visited(visited_z, target)
        z_cur, u_cur, up_cur = visited_z[idx_v], visited_u[idx_v], visited_up[idx_v]

        n_steps = 0
        while abs(z_cur - target) > h_T
            n_steps += 1
            n_steps > max_steps_per_target && throw(ErrorException(
                "path_network_solve: target $target unreachable in " *
                "$max_steps_per_target steps from $(visited_z[idx_v]); " *
                "current z = $z_cur, |z - target| = $(abs(z_cur - target))."))

            goal_dir = angle(target - z_cur)
            evals = _wedge_evaluations(prob.f, z_cur, u_cur, up_cur,
                                       order, h_T, goal_dir, wedge_angles)

            idx_sel = _select_candidate(step_selection, evals, u_cur, up_cur,
                                        goal_dir, wedge_angles)

            z_new, u_new, up_new, pade_sel = evals[idx_sel]

            (!isfinite(abs(u_new)) || pade_sel === nothing) && throw(ErrorException(
                "path_network_solve: all 5 wedge candidates failed at " *
                "z=$z_cur (target=$target); shrink h or widen wedge."))

            push!(visited_z, z_new)
            push!(visited_u, u_new)
            push!(visited_up, up_new)
            push!(visited_pade, pade_sel)
            push!(visited_h, h_T)

            z_cur, u_cur, up_cur = z_new, u_new, up_new
        end
    end

    # Stage 2: fine-grid extrapolation.
    grid_z  = collect(CT, grid)
    grid_u  = Vector{CT}(undef, length(grid_z))
    grid_up = Vector{CT}(undef, length(grid_z))
    nan_CT  = CT(NaN, NaN)
    for (i, z_f) in enumerate(grid_z)
        idx_v = _nearest_visited(visited_z, z_f)
        z_v   = visited_z[idx_v]
        h_v   = visited_h[idx_v]
        if abs(z_f - z_v) > h_v
            grid_u[i]  = nan_CT
            grid_up[i] = nan_CT
        else
            t = (z_f - z_v) / h_v
            grid_u[i]  = _evaluate_pade(visited_pade[idx_v], t)
            grid_up[i] = _evaluate_pade_deriv(visited_pade[idx_v], t) / h_v
        end
    end

    return PathNetworkSolution{T}(visited_z, visited_u, visited_up,
                                  visited_pade, visited_h,
                                  grid_z, grid_u, grid_up)
end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Compute the canonical local Padé centered at (z, u, up) with real step h.
# We clone the state so the IC isn't mutated.
function _local_pade(f, z::Complex{T}, u::Complex{T}, up::Complex{T},
                     order::Integer, h::T) where T
    s = PadeStepperState{Complex{T}}(z, u, up)
    _, pade = pade_step_with_pade!(s, f, Int(order), h)
    return pade
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

# Try each of 5 wedge directions; return Vector of (z_new, u_new, up_new, pade).
# A failed candidate gets (z_cur, Inf+0im, 0+0im, nothing) so it loses min-|u|.
function _wedge_evaluations(f, z_cur::Complex{T}, u_cur::Complex{T},
                            up_cur::Complex{T}, order::Integer, h::T,
                            goal_dir::T,
                            wedge_angles::AbstractVector{<:Real}) where T
    CT = Complex{T}
    evals = Vector{Tuple{CT, CT, CT, Union{Nothing, PadeApproximant{CT}}}}(undef, 5)
    for (k, θ) in enumerate(wedge_angles)
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

# TIER-3 INTERFACE: BVP dispatcher (bead padetaylor-804) consumes
# PathNetworkSolution.{grid_z, grid_u, grid_up}.  No code here until
# Phase 11; see docs/adr/0004-path-network-architecture.md.

end # module PathNetwork
