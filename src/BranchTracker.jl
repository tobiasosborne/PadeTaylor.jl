"""
    PadeTaylor.BranchTracker

Walker-side enforcement layer for Riemann-sheet bookkeeping during
`PathNetwork.path_network_solve`.  Complementary to the caller-side
`SheetTracker` primitives shipped in worklog 018 — where SheetTracker
analyses an already-walked path, BranchTracker decides per-step
whether a candidate wedge direction is allowed (cuts respected) and
updates the per-branch sheet counters when a step deliberately
crosses a cut.

## Why a separate module (ADR-0013 §"Why a new helper module")

PathNetwork.jl is already past the Rule-6 LOC cap.  This module
contains the geometric primitives (segment-ray intersection,
cut-crossing predicate, per-branch sheet update) in one cohesive
place, independently testable, while PathNetwork's wedge loop grows
by only ~50 LOC of kwarg threading and per-candidate calls.

## Branch cut model (ADR-0013 §"Decision")

Each branch point `b ∈ ℂ` carries a **half-line cut** running from
`b` along direction `arg = α`:

    C(b, α) = { b + s · exp(i · α) : s ≥ 0 }

The default `α = π` matches Julia's `log` convention — cuts run
from each branch point along the negative-real direction.  Per-branch
`α_k` lets callers steer cuts away from regions they intend the
walker to visit.

## Cut-crossing predicate

Given a step segment `S = [z_cur, z_new]` and a cut `C(b, α)`, the
predicate `segment_crosses_cut(z_cur, z_new, b, α)` returns `true`
iff the open segment crosses the open ray.  The intersection is the
solution `(t, s) ∈ ℝ²` of

    z_cur + t · (z_new - z_cur) = b + s · exp(i · α)

with `0 < t < 1` and `s > 0`.  Closed endpoint conventions:

  - `t = 0` and `t = 1` are NOT counted (lands exactly on the start
    or end of the step — a measure-zero numerical accident).
  - `s = 0` is NOT counted (the branch point itself is treated as
    crossed only when the ray side `s > 0` is crossed; `s = 0`
    represents the walker passing *through* the branch point,
    which is its own kind of degeneracy and is handled by the
    Padé stepper's fail-loud at singular `Q(z)`).

Parallel segments (`det = 0`) return `false`.

## Sheet update law (cross-mode)

In `cross_branch = true` mode, for each branch `b_k` whose cut is
crossed by the step `[z_cur, z_new]`, the sheet counter for `b_k`
is updated by `sign(Δθ_k)` where `Δθ_k = winding_delta(z_cur, z_new,
b_k)` (the SheetTracker primitive shipped in worklog 018,
normalised to `(-π, π]`).  Counter-clockwise crossing (positive
`Δθ`) bumps `+1`; clockwise bumps `-1`.  No crossing → no bump
(per-step `Δθ` may be nonzero — small rotation NOT enclosing the
branch — and that's the unambiguous signal to leave the counter
alone).

The per-step caller contract is the same as SheetTracker's: the
step magnitude must be small enough that `|Δθ_k| < π` at every
branch (else the principal-branch normalisation in `winding_delta`
hides one full revolution).  At the FW 2011 default `h = 0.5`
with branches at distance `≥ 1`, this holds with ~30° margin.

## Caller responsibilities

  - `branch_points` is a tuple of `Complex` (any element type
    compatible with the problem's `T`).  Empty tuple is the no-op.
  - `branch_cut_angles` is a scalar `Real` (broadcast to all
    branches) or a tuple of the same length as `branch_points`.
  - Sheet counters are integers (`Int`) per branch, threaded through
    `PathNetworkSolution.visited_sheet :: Vector{Vector{Int}}` —
    one inner vector per visited node, each inner vector of length
    `length(branch_points)`.

  Fail-loud predicate: in `refuse` mode, when all wedge candidates
  are forbidden, the caller (PathNetwork's wedge loop) throws an
  `ErrorException` whose message names `cross_branch = true` as the
  most common fix.  This module does not throw — it returns booleans
  and updated counters.
"""
module BranchTracker

using ..SheetTracker: winding_delta

export segment_crosses_cut,
       any_cut_crossed,
       step_sheet_update,
       resolve_cut_angles

# -----------------------------------------------------------------------------
# Cut-crossing predicate
# -----------------------------------------------------------------------------

"""
    segment_crosses_cut(z_cur, z_new, branch, cut_angle) -> Bool

Test whether the segment `[z_cur, z_new]` crosses the half-line cut
`{branch + s · exp(i·cut_angle) : s ≥ 0}`.  See module docstring for
endpoint conventions.

Algorithm: parametrise segment as `z_cur + t·(z_new - z_cur)` and
ray as `branch + s·exp(i·α)`.  Equate and solve the 2×2 linear
system via Cramer's rule on complex-valued imaginary parts.
"""
function segment_crosses_cut(z_cur, z_new, branch, cut_angle)
    d   = z_new - z_cur
    δ   = branch - z_cur
    u   = cis(cut_angle)                      # = exp(i·α)
    det = imag(d * conj(u))
    iszero(det) && return false               # parallel
    t = imag(δ * conj(u)) / det
    s = imag(δ * conj(d)) / det
    return 0 < t < 1 && s > 0
end

"""
    any_cut_crossed(z_cur, z_new, branches, cut_angles) -> Bool

True iff the step crosses ANY of the listed branch cuts.  `branches`
and `cut_angles` must have equal length (use `resolve_cut_angles` to
broadcast a scalar `cut_angles` to a tuple of the right length first).
"""
function any_cut_crossed(z_cur, z_new, branches, cut_angles)
    @assert length(branches) == length(cut_angles)
    @inbounds for k in eachindex(branches)
        segment_crosses_cut(z_cur, z_new, branches[k], cut_angles[k]) &&
            return true
    end
    return false
end

# -----------------------------------------------------------------------------
# Sheet-counter update
# -----------------------------------------------------------------------------

"""
    step_sheet_update(sheet_old, z_cur, z_new, branches, cut_angles) -> Vector{Int}

Per-step sheet-counter update for `cross_branch = true` mode.  Returns
a NEW `Vector{Int}` of the same length as `sheet_old` (and
`branches`); for every branch whose cut is crossed by `[z_cur,
z_new]`, the corresponding counter is updated by
`sign(winding_delta(z_cur, z_new, branch))`.  Branches whose cut is
NOT crossed are passed through unchanged.

Allocates a fresh vector each call (cheap: typically 1-3 Int entries);
PathNetwork's wedge loop calls this at most once per accepted step.
"""
function step_sheet_update(sheet_old::AbstractVector{<:Integer},
                           z_cur, z_new, branches, cut_angles)
    @assert length(sheet_old) == length(branches) == length(cut_angles)
    out = collect(Int, sheet_old)             # fresh allocation, Int-typed
    @inbounds for k in eachindex(branches)
        if segment_crosses_cut(z_cur, z_new, branches[k], cut_angles[k])
            Δθ = winding_delta(z_cur, z_new, branches[k])
            out[k] += Δθ > 0 ? 1 : (Δθ < 0 ? -1 : 0)
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Kwarg-resolution helper
# -----------------------------------------------------------------------------

"""
    resolve_cut_angles(branches, angles) -> NTuple{N, Float64}

Normalise the public `branch_cut_angles` kwarg to a tuple of `Float64`
matching `length(branches)`.  A scalar broadcasts.  A tuple/vector of
matching length passes through (with element type coercion to
`Float64`).  Mismatched length throws `ArgumentError` (CLAUDE.md
Rule 1).

Used at the top of `path_network_solve` to validate user input once
and pass the resolved tuple to `any_cut_crossed` /
`step_sheet_update` without per-call branching.
"""
function resolve_cut_angles(branches, angles::Real)
    return ntuple(_ -> Float64(angles), length(branches))
end

function resolve_cut_angles(branches, angles)
    length(angles) == length(branches) || throw(ArgumentError(
        "BranchTracker.resolve_cut_angles: branch_cut_angles length " *
        "$(length(angles)) does not match branch_points length " *
        "$(length(branches)).  Pass a scalar to broadcast, or a tuple " *
        "of matching length for per-branch control."))
    return ntuple(k -> Float64(angles[k]), length(branches))
end

end # module BranchTracker
