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

  - `pVI_eta_transformed_rhs(α, β, γ, δ) -> (η, v, vp) -> v''`
    Factory returning the η-plane PVI RHS closure (FFW 2017 eq. 5,
    md:154).  The second exponential `ζ = e^η` stacked atop `z = e^ζ`
    maps PVI's `ζ = 0` branch point (and thus `z = 0`) out of the
    finite plane, leaving the `z = 1` lattice at `η = log|2π·k| +
    i·arg(2π·i·k)` for `|k| ≥ 1` (md:148, eq. 4).  Region `Re η <
    log(2π) ≈ 1.838` is branch-point-free (md:157).

  - `pVI_z_to_η(z, u, up) -> (η, v, vp)`,
    `pVI_η_to_z(η, v, vp) -> (z, u, up)`
    Composition of the two exponential transforms (`ζ = log z`, `η =
    log ζ`), with their inverses.  `v(η) = u(z)` (the dependent
    variable is unchanged); the chain rule gives `vp = z · ζ · up`
    forward and `up = vp / (z · ζ)` back.

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

## η-plane double-exp transform (FFW 2017 §2.2.1, md:146-162)

The η-plane equation is the price-paid form of the ζ-plane equation:
nested `exp(exp(η))` arithmetic in every term, in exchange for a
*finite* representation of the `z = 0` half-line on the *finite*
η-plane.  The ζ-plane handles this only by walking off to `Re ζ = -∞`.

The branch-point-free region `Re η < log(2π) ≈ 1.83788` is rectangular
in `(Re η, Im η)`; FFW Fig 2 column 1 shows a PVI solution computed
there.  Outside this region (`Re η > log(2π)`) the lattice points
`η = log|2π·k| + i·arg(2π·i·k)` are branch points of the η-plane
equation — they correspond to `ζ = 2π·i·k` (which corresponded to
`z = 1`) — and circumambulation in the η-plane is required, no
different in spirit from §2.2.2 ζ-plane circumambulation, just on a
different parametrisation of the surface.

When to prefer the η-plane over the ζ-plane (md:157-162): the η-plane
gives a *compact* rectangular pole-free region containing infinitely
many ζ-sheets folded together; the ζ-plane gives a strip-per-sheet
unfolding of the same surface.  For producing a single figure of one
solution on one rectangular region (Fig 2 column 1), the η-plane is
the right tool; for tracking analytic continuation across sheets (Figs
2 column 2, 3), the ζ-plane is.

## What's NOT in v1 (deferral notes)

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
export pVI_eta_transformed_rhs, pVI_z_to_η, pVI_η_to_z
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
# PVI η-plane transformed RHS (FFW 2017 eq. 5, md:154)
# -----------------------------------------------------------------------------

"""
    pVI_eta_transformed_rhs(α, β, γ, δ) -> (η, v, vp) -> v''

Factory returning the η-plane PVI RHS closure (FFW 2017 eq. 5,
md:154).  Plug directly into
`PadeTaylorProblem(rhs, (v₀, v'₀), (η_start, η_end); order = ...)`.

The η-plane equation, written character-by-character from FFW 2017
md:154 (`E = exp(exp(η))` abbreviates the nested exponential):

    d²v/dη² = (1/2)(1/v + 1/(v-1) + 1/(v - E)) (dv/dη)²
              - (e^η·E/(E - 1) + e^η·E/(v - E) - 1) (dv/dη)
              + v(v-1)(v-E) · e^(2η) / (E - 1)² ·
                (α + β E/v² + γ(E-1)/(v-1)² + δ E(E-1)/(v - E)²)

The `-1` inside the `(dv/dη)` bracket is the only place where the
η-plane equation is NOT a structural copy of the ζ-plane equation
(eq. 3, md:144) with `e^ζ → E`: it is the chain-rule artefact from
stacking `ζ = e^η` on top of `z = e^ζ`.  Forgetting it is the easy
hand-derivation mistake.

Requires `v ∉ {0, 1, E}` and `E ≠ 1` (i.e., `e^η ∉ {2π·i·k : k ∈ ℤ}`,
which on the principal branch means `η ∉ {log|2π·k| + i·arg(2π·i·k) :
|k| ≥ 1}`).  These are the fixed singular surfaces of the η-plane
equation; the region `Re η < log(2π)` is free of them.
"""
function pVI_eta_transformed_rhs(α, β, γ, δ)
    return (η, v, vp) -> begin
        eη   = exp(η)
        E    = exp(eη)            # e^(e^η)
        E_m1 = E - 1
        v_m1 = v - 1
        v_mE = v - E

        first  = (1 / v + 1 / v_m1 + 1 / v_mE) * vp^2 / 2

        second = -(eη * E / E_m1 + eη * E / v_mE - 1) * vp

        param  = α +
                 β * E / v^2 +
                 γ * E_m1 / v_m1^2 +
                 δ * E * E_m1 / v_mE^2
        third  = v * v_m1 * v_mE * eη^2 / E_m1^2 * param

        return first + second + third
    end
end

# -----------------------------------------------------------------------------
# PVI z ↔ η coordinate transforms (composition of two exponentials)
# -----------------------------------------------------------------------------

"""
    pVI_z_to_η(z, u, up) -> (η, v, vp)

Convert a PVI state `(z, u(z), u'(z))` to the η-plane state `(η,
v(η), dv/dη)` via the composition `ζ = log z`, `η = log ζ` (principal
branches).  Since `v(η) = w(ζ) = u(z)` (the dependent variable is
unchanged through both transforms; FFW2017...md:146, md:151), the
chain rule gives `dv/dη = (dw/dζ)·(dζ/dη) = (z·u')·ζ = z·log(z)·u'`.

Singular at `z = 0` (where `ζ = -∞`) and at `z = 1` (where `ζ = 0`,
so `η = -∞`); both correspond to PVI's fixed singularities.
"""
function pVI_z_to_η(z, u, up)
    ζ  = log(z)
    η  = log(ζ)
    v  = u
    vp = z * ζ * up
    return (η, v, vp)
end

"""
    pVI_η_to_z(η, v, vp) -> (z, u, up)

Inverse of `pVI_z_to_η`.  `ζ = exp(η)`, `z = exp(ζ) = exp(exp(η))`,
`u = v`, `u' = vp / (z·ζ)`.  Round-trip with `pVI_z_to_η` is exact
to floating-point error on the principal branch.
"""
function pVI_η_to_z(η, v, vp)
    ζ  = exp(η)
    z  = exp(ζ)
    u  = v
    up = vp / (z * ζ)
    return (z, u, up)
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
