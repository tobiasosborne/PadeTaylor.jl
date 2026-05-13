"""
    PadeTaylor.SheetTracker

Tier-5 architectural wrapper (Phase 14, bead `padetaylor-grc`) for
PVI — the sixth Painlevé equation, whose solutions live on multi-sheet
Riemann surfaces because of fixed singularities at BOTH `z = 0` and
`z = 1`.  By Picard's theorem (FFW2017...md:137), no single coordinate
transform can map both branch points out of the finite plane; the
single-transform trick that worked for PIII/PV (Phase 13's
`CoordTransforms`) is fundamentally inadequate here.

This module ships the **ζ-plane PVI transform** (FFW 2017 eq. 3 at
md:144), which DOES map `z = 0` to `ζ = -∞` — but leaves the second
branch point at `z = 1` as an infinite lattice `ζ = 2π·i·k`, `k ∈ ℤ`.
The trade-off: an equation that integrates everywhere except on the
imaginary-axis lattice, where deliberate **circumambulation** (FFW
2017 §2.2.2, md:163-178) is required to walk *around* each branch
point without "overstepping" the branch cuts.

To support circumambulation, this module also ships **winding-number
primitives** that consume a path (e.g., `PathNetworkSolution.visited_z`)
plus a branch point and accumulate the signed angle traversed, then
expose the resulting Riemann-sheet index `s ∈ ℤ`.  Higher-level routing
(constraining the wedge selection so paths do NOT overstep cuts) is
deferred — the primitives in this module are the building blocks, but
the full PathNetwork-level routing is a follow-up.

## What's in v1

  - `pVI_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''`
    Factory returning the ζ-plane P̃_VI RHS closure (FFW 2017 eq. 3,
    md:144).  Plug into `PadeTaylorProblem`.  Note: the z↔ζ coordinate
    conversion is mathematically identical to PV's (FFW2017...md:146
    "obtained by setting u(z) = w(ζ), z = e^ζ in PVI"), so callers
    reuse `pV_z_to_ζ` / `pV_ζ_to_z` from `CoordTransforms`.

  - `winding_delta(z_old, z_new, branch) -> Float64`
    Signed angle change `Δθ = arg(z_new - branch) - arg(z_old - branch)`
    normalised to `(-π, π]`.  Caller must ensure path steps are small
    enough that no single step's `|Δθ| ≥ π` (else the normalisation
    masks the discontinuity; document at usage sites).

  - `accumulate_winding(path, branch) -> Vector{Float64}`
    Cumulative winding along a path with respect to a branch point.
    `out[i]` is the total signed angle traversed from `path[1]` to
    `path[i]`.  `out[1] == 0` by convention (the IC contributes no
    winding).

  - `sheet_index(total_winding) -> Int`
    Convert accumulated winding to integer sheet index `s = round(total
    / 2π)`.  A counterclockwise loop around `branch` adds `+2π` and
    increments `s` by `+1`; clockwise loop subtracts and decrements.

## What's NOT in v1 (deferral notes)

  - **η-plane PVI transform** (FFW 2017 eq. 5, md:154).  The
    second exponential `ζ = e^η` makes the branch-point-free region
    `Re η < log(2π)` more compact at the cost of nested `e^(e^η)`
    arithmetic.  Useful for FFW Fig 2 first column ONLY; the
    ζ-plane suffices for Figs 2 (column 2), 3, 7.  File a bead if a
    downstream caller actually needs the η-plane.

  - **Constrained-wedge PathNetwork variant** that REFUSES to overstep
    branch cuts and INCREMENTS the sheet counter on deliberate
    circumambulation.  This is the "production-grade circumambulation"
    referenced in the bead's `cross_branch=true` kwarg sketch.  Bigger
    architectural change in PathNetwork.jl; the v1 primitives here
    let a caller compute the sheet index AFTER a regular PathNetwork
    walk, but they don't ENFORCE branch-cut avoidance during the walk.
    File a bead if needed.

  - **Multi-branch winding** (winding around BOTH `ζ = 0` AND
    `ζ = 2π·i`).  The primitives here work for ONE branch at a time;
    callers wanting multi-branch sheet indexing compose two
    `accumulate_winding` calls.  Adequate for v1.

## Sheet semantics (FFW2017...md:180-189)

For PVI in the ζ-plane with branch points at `ζ = 2π·i·k`, FFW's sheet
parametrisation (eq. 6, md:187-189) maps sheet `k` to angular ranges
`(2k-1)π < θ_k ≤ (2k+1)π` for the `ζ = 0`-branch, with analogous
ranges for the `ζ = 1`-equivalent branch.  In `[-π, π]` convention,
sheet 0 is the principal sheet.  The integer returned by
`sheet_index` matches this convention.

## References

  - `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md`
    :29-37 (PVI equation), :135-160 (the ζ- and η-plane transforms),
    :163-189 (circumambulation + sheet index), :191-195 (Fig 2 caption
    with IC + parameters).
  - `docs/figure_catalogue.md §5` row T5.
  - `src/CoordTransforms.jl` — `pV_z_to_ζ` reused for PVI's coordinate
    conversion (mathematically identical map, distinct equation).
"""
module SheetTracker

export pVI_transformed_rhs
export winding_delta, accumulate_winding, sheet_index

# -----------------------------------------------------------------------------
# PVI ζ-plane transformed RHS (FFW 2017 eq. 3, md:144)
# -----------------------------------------------------------------------------

"""
    pVI_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''

Factory returning the ζ-plane PVI RHS closure.  Plug directly into
`PadeTaylorProblem(rhs, (w₀, w'₀), (ζ_start, ζ_end); order = ...)`.

The transformed equation (FFW2017...md:144):

    d²w/dζ² = (1/2)(1/w + 1/(w-1) + 1/(w - e^ζ)) (dw/dζ)²
              - (e^ζ/(e^ζ - 1) + e^ζ/(w - e^ζ)) (dw/dζ)
              + w(w-1)(w-e^ζ) / (e^ζ - 1)² ·
                (α + β e^ζ/w² + γ(e^ζ - 1)/(w-1)² + δ e^ζ(e^ζ - 1)/(w - e^ζ)²)

Requires `w ∉ {0, 1, e^ζ}` and `e^ζ ≠ 1` (i.e., `ζ ≠ 2π·i·k`).  These
are the fixed singular surfaces of the transformed equation, exactly
the branch-point lattice that motivates this module's existence.
"""
function pVI_transformed_rhs(α, β, γ, δ)
    return (ζ, w, wp) -> begin
        eζ = exp(ζ)
        eζ_m1 = eζ - 1
        w_m1  = w - 1
        w_meζ = w - eζ

        first  = (1 / w + 1 / w_m1 + 1 / w_meζ) * wp^2 / 2

        second = -(eζ / eζ_m1 + eζ / w_meζ) * wp

        param  = α +
                 β * eζ / w^2 +
                 γ * eζ_m1 / w_m1^2 +
                 δ * eζ * eζ_m1 / w_meζ^2
        third  = w * w_m1 * w_meζ / eζ_m1^2 * param

        return first + second + third
    end
end

# -----------------------------------------------------------------------------
# Winding-number primitives
# -----------------------------------------------------------------------------

"""
    winding_delta(z_old, z_new, branch) -> Float64

Signed angle change from `z_old` to `z_new` as seen from `branch`,
normalised to `(-π, π]`.  Caller must ensure no single path step has
`|Δθ| ≥ π` (else the normalisation hides the wrap and the cumulative
sum loses one full revolution at that step).

For path-network walks with step `h = 0.5` and branch points at
distance `> 1`, single-step `|Δθ|` is bounded by `arcsin(0.5/1) ≈ 0.52
≈ 30°` — well within the safe range.
"""
function winding_delta(z_old, z_new, branch)
    Δθ = angle(z_new - branch) - angle(z_old - branch)
    # Normalise to (-π, π].
    if Δθ ≤ -π
        Δθ += 2 * π
    elseif Δθ > π
        Δθ -= 2 * π
    end
    return Δθ
end

"""
    accumulate_winding(path::AbstractVector{<:Complex}, branch) -> Vector{Float64}

Cumulative winding angle (in radians) along `path` w.r.t. `branch`.
`out[1] == 0.0` (the starting point contributes no winding).  `out[i]`
for `i > 1` is the total signed `Δθ` summed from `path[1]` to `path[i]`.

For a closed counterclockwise loop enclosing `branch`, `out[end] ≈
2π`; clockwise gives `-2π`; non-enclosing gives `≈ 0`.
"""
function accumulate_winding(path::AbstractVector{<:Complex}, branch)
    n = length(path)
    out = zeros(Float64, n)
    @inbounds for i in 2:n
        out[i] = out[i-1] + winding_delta(path[i-1], path[i], branch)
    end
    return out
end

"""
    sheet_index(total_winding::Real) -> Int

Convert accumulated winding to integer Riemann-sheet index per FFW
2017's convention (md:187-189).  `s = round(total / 2π)`: +1 per
counterclockwise revolution, -1 per clockwise.
"""
sheet_index(total_winding::Real) = round(Int, total_winding / (2 * π))

end # module SheetTracker
