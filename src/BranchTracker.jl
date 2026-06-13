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

The per-step caller contract is that no single step may GRAZE a
tracked branch.  `winding_delta` itself is always correct — a
straight chord can never subtend `|Δθ| ≥ π` about an exterior
branch, so the principal-value result IS the chord's true winding
(bug `padetaylor-61um` corrected the earlier "loses a full
revolution" framing, which was geometrically false for the
straight steps the walker realises).  The genuine hazard is a step
whose chord passes very CLOSE to the branch relative to its length:
then the discrete two-node chord cannot resolve which side of the
branch the true continuation passed, and the `sign(Δθ_k)` bump is
ill-conditioned.  `step_sheet_update` ENFORCES this (Rule 1): a
crossed-cut step with grazing ratio `min_dist / |chord| < graze_tol`
(default `WINDING_GRAZE_TOL = 0.1`, calibrated on the in-suite
corpus; equivalently `|Δθ_k| → π`) throws a `DomainError` with a
"shorten h / densify nodes near the branch" suggestion rather than
silently mis-accumulating.  At the FW 2011 default `h = 0.5` with
branches at distance `≥ 1` the grazing ratio stays well above
`graze_tol`, so the guard never trips a well-conditioned walk.  See
ADR-0031.

## Caller responsibilities

  - `branch_points` is a tuple of `Complex` (any element type
    compatible with the problem's `T`).  Empty tuple is the no-op.
  - `branch_cut_angles` is a scalar `Real` (broadcast to all
    branches) or a tuple of the same length as `branch_points`.
  - Sheet counters are integers (`Int`) per branch, threaded through
    `PathNetworkSolution.visited_sheet :: Vector{Vector{Int}}` —
    one inner vector per visited node, each inner vector of length
    `length(branch_points)`.

  Fail-loud predicates: in `refuse` mode, when all wedge candidates
  are forbidden, the caller (PathNetwork's wedge loop) throws an
  `ErrorException` whose message names `cross_branch = true` as the
  most common fix.  `step_sheet_update` itself throws a `DomainError`
  on a winding-ambiguous (branch-grazing) crossed-cut step (the
  `graze_tol` guard above; bug `padetaylor-61um`).  `segment_crosses_cut`
  / `any_cut_crossed` / `resolve_cut_angles` do not throw on geometry —
  they return booleans / counters / the resolved tuple.
"""
module BranchTracker

using ..SheetTracker: winding_delta

export segment_crosses_cut,
       any_cut_crossed,
       step_sheet_update,
       resolve_cut_angles

# Default grazing-ratio threshold for the winding-ambiguity guard in
# `step_sheet_update` (bug padetaylor-61um).  A crossed-cut single step whose
# chord's closest approach to the branch is below this fraction of the chord
# length is refused as winding-ambiguous (the discrete two-node chord cannot
# resolve which side of the branch the true continuation passed).  CALIBRATED on
# the in-suite `cross_branch=true` corpus (2026-06-13, measured at every
# crossed-cut step across corpus_winding/branch_tracker/path_network_branch/
# sheet_aware_stage2/corpus_pathnet_winding/ffw_fig_2/3/7): the LONE ambiguous
# step (CPN.7.3 FINE) has ratio 0.0605; the NEAREST safe step (CPN.7.2 COARSE,
# which correctly returns [1]) has ratio 0.2045; every other in-suite walk step
# has ratio ≥ 0.71.  The default 0.1 catches FINE with 1.65× margin and clears
# COARSE by 2.05×.  See ADR-0031 for the full distribution and the
# |winding_delta| → π angular equivalence.  Override per-call via the
# `graze_tol` kwarg; threading it through `path_network_solve` is a deferred lift
# (bead padetaylor-61um notes; no in-suite walk needs it).
const WINDING_GRAZE_TOL = 0.1

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
    step_sheet_update(sheet_old, z_cur, z_new, branches, cut_angles;
                      graze_tol = WINDING_GRAZE_TOL) -> Vector{Int}

Per-step sheet-counter update for `cross_branch = true` mode.  Returns
a NEW `Vector{Int}` of the same length as `sheet_old` (and
`branches`); for every branch whose cut is crossed by `[z_cur,
z_new]`, the corresponding counter is updated by
`sign(winding_delta(z_cur, z_new, branch))`.  Branches whose cut is
NOT crossed are passed through unchanged.

Throws `DomainError` (Rule 1, bug `padetaylor-61um`) when a crossed-cut
step GRAZES the branch — chord closest-approach `min_dist < graze_tol ·
|chord|` — because then the discrete two-node chord cannot resolve which
side of the branch the true continuation passed and the `sign(Δθ)` bump
is ill-conditioned (equivalently `|winding_delta| → π`).  `graze_tol`
defaults to `WINDING_GRAZE_TOL = 0.1` (calibrated on the in-suite
`cross_branch` corpus; see the module docstring's "Sheet update law" and
ADR-0031); pass a smaller value to admit nearer-grazing steps or a larger
one to be stricter.  A well-conditioned step (`h` small, branches not
skimmed) never trips it.

Allocates a fresh vector each call (cheap: typically 1-3 Int entries);
PathNetwork's wedge loop calls this at most once per accepted step.
"""
function step_sheet_update(sheet_old::AbstractVector{<:Integer},
                           z_cur, z_new, branches, cut_angles;
                           graze_tol::Real = WINDING_GRAZE_TOL)
    @assert length(sheet_old) == length(branches) == length(cut_angles)
    out = collect(Int, sheet_old)             # fresh allocation, Int-typed
    @inbounds for k in eachindex(branches)
        if segment_crosses_cut(z_cur, z_new, branches[k], cut_angles[k])
            # Rule-1 fail-loud winding-ambiguity guard (bug padetaylor-61um).
            # A single STRAIGHT walker step whose chord GRAZES the branch is
            # winding-AMBIGUOUS: the discrete two-node chord cannot tell whether
            # the true ODE continuation passed the branch on the short or the
            # long side, so the sign(winding_delta) sheet bump is unreliable.
            # The signal is the chord's closest approach to the branch RELATIVE
            # to the chord length (the grazing ratio); equivalently
            # |winding_delta| → π.  Below `graze_tol` the bump is REFUSED loudly
            # rather than silently mis-accumulated into a corrupt sheet index.
            #
            # NOTE the defect was NEVER a "lost revolution": winding_delta is
            # mathematically CORRECT for the straight chord (it returns the
            # principal-value winding, which a straight chord can never push past
            # ±π).  The real bug was this UNGUARDED ill-conditioned step.  See
            # ADR-0031 and the SheetTracker.winding_delta docstring.
            d    = z_new - z_cur
            ts   = real(conj(d) * (branches[k] - z_cur)) / abs2(d)
            ts   = clamp(ts, zero(ts), oneunit(ts))
            mind = abs(z_cur + ts * d - branches[k])
            mind < graze_tol * abs(d) && throw(DomainError((z_cur, z_new, branches[k]),
                "step_sheet_update: the straight step from $z_cur to $z_new grazes " *
                "branch $k (at $(branches[k])): closest approach $mind over a chord of " *
                "length $(abs(d)) (ratio $(mind / abs(d)) < graze_tol=$graze_tol).  The " *
                "Riemann-sheet winding is AMBIGUOUS — the chord cannot resolve which " *
                "side of the branch the true analytic continuation passed (equivalently " *
                "|winding_delta| → π).  Suggestion: shorten h, or densify nodes near " *
                "the branch (e.g. node_separation R(z) that shrinks toward branch_points), " *
                "so each step clears the branch, then retry."))
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
