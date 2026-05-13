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

## Tier-3-to-5 deferrals (per ADR-0004)

  - **BVP composition (Tier 3)**: stub comment at end of this file.
  - **CoordTransforms (Tier 4)**: external wrapper, no change here.
  - **SheetTracker (Tier 5)**: external wrapper, no change here.
  - **`:adaptive_ffw` step-size policy**: throws `ArgumentError` per
    fail-fast discipline; Tier-4 work item.
"""
module PathNetwork

using Base.Threads:  nthreads, @threads
using Printf:        @sprintf, @printf
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
                       rng_seed=0,
                       enforce_real_axis_symmetry=false,
                       verbose=false, progress_every=500) -> PathNetworkSolution

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
                            order::Integer = 30,
                            wedge_angles::AbstractVector{<:Real} = DEFAULT_WEDGE,
                            step_selection::Symbol  = :min_u,
                            step_size_policy::Symbol = :fixed,
                            max_steps_per_target::Integer = 1000,
                            rng_seed::Integer = 0,
                            enforce_real_axis_symmetry::Bool = false,
                            verbose::Bool = false,
                            progress_every::Integer = 500)

    if enforce_real_axis_symmetry
        return _solve_with_schwarz_reflection(
            prob, grid; h, order, wedge_angles, step_selection,
            step_size_policy, max_steps_per_target, rng_seed,
            verbose, progress_every)
    end

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

    if verbose
        println("path_network_solve: Stage 1 begin — ",
                length(targets), " targets, h=", h_T, ", order=", order,
                ", threads=", nthreads())
        flush(stdout)
    end

    t_start_stage1 = time()
    total_n_steps  = 0
    for (target_idx, target) in enumerate(targets)
        any(z -> isapprox(z, target; atol=10*eps(T)), visited_z) && continue

        idx_v = _nearest_visited(visited_z, target)
        z_cur, u_cur, up_cur = visited_z[idx_v], visited_u[idx_v], visited_up[idx_v]

        verbose && _verbose_target_start(target_idx, length(targets),
                                         target, z_cur, t_start_stage1)

        n_steps = 0
        while abs(z_cur - target) > h_T
            n_steps        += 1
            total_n_steps  += 1
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

            # Canonical (real-h) Padé centered at z_new, per ADR-0004's
            # design decision: Stage 2 always interpolates inside |t| ≤ 1
            # of a REAL-direction disc.  The wedge-direction `pade_sel`
            # above is consumed by step selection only; storing it as the
            # visited node's Padé would invalidate the `t = (z_f - z_v)/h`
            # Stage-2 evaluation (which assumes real h).  Cost: one extra
            # Coefficients + RobustPade per visited node.
            pade_canonical = _local_pade(prob.f, z_new, u_new, up_new, order, h_T)

            push!(visited_z, z_new)
            push!(visited_u, u_new)
            push!(visited_up, up_new)
            push!(visited_pade, pade_canonical)
            push!(visited_h, h_T)

            z_cur, u_cur, up_cur = z_new, u_new, up_new

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
                                         max_steps_per_target::Integer,
                                         rng_seed::Integer,
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

    return PathNetworkSolution{T}(
        sol_upper.visited_z, sol_upper.visited_u, sol_upper.visited_up,
        sol_upper.visited_pade, sol_upper.visited_h,
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

# TIER-3 INTERFACE: BVP dispatcher (bead padetaylor-804) consumes
# PathNetworkSolution.{grid_z, grid_u, grid_up}.  No code here until
# Phase 11; see docs/adr/0004-path-network-architecture.md.

end # module PathNetwork
