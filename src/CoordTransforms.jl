"""
    PadeTaylor.CoordTransforms

Tier-4 architectural wrapper (Phase 13, bead `padetaylor-bvh`) —
exponential coordinate transforms that map the **fixed branch point at
`z = 0`** in PIII and PV out of the finite complex plane.  In the
transformed `ζ`-plane the resulting equations P̃_III and P̃_V have
meromorphic solutions, so the existing pole-friendly Padé-Taylor IVP
machinery (PathNetwork, BVP, etc.) — built for single-valued
meromorphic ODEs in the FW 2011 tradition — applies unchanged.

The "Tier-4" name comes from `docs/figure_catalogue.md §5`: the FFW
2017 figures 1, 4, 5, 6 cannot be reproduced without these transforms.

## The transformations (FFW 2017 §2.1, lines 39-48)

  - **PIII**: `z = exp(ζ/2)`, `u(z) = exp(-ζ/2) · w(ζ)`.  Strip width
    `4π` in the ζ-plane corresponds to ONE sheet in the z-plane.  The
    transformed equation is `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:43`:

        w'' = (1/w)(w')² + (1/4)(α w² + γ w³ + β eᶻ + δ e^(2ζ) / w)

  - **PV**: `z = exp(ζ)`, `u(z) = w(ζ)` (no `u`-side rescaling).
    Strip width `2π` per z-sheet.  Transformed equation
    (FFW2017...md:47):

        w'' = (1/(2w) + 1/(w-1))(w')²
              + (w - 1)² (α w + β/w)
              + γ eᶻ w
              + δ e^(2ζ) w (w + 1) / (w - 1)

Both transformed equations are meromorphic in `w(ζ)` (FFW2017...md:43-49,
citing the meromorphy results of refs. [14,15]).  No branch point at
`ζ = -∞` because the exponential transform smooths it out.

## Public API

  - `pIII_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''`
    Factory returning the P̃_III RHS closure suitable for
    `PadeTaylorProblem`.  Captures the four PIII parameters.

  - `pV_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''`
    Same for P̃_V.

  - `pIII_z_to_ζ(z, u, up) -> (ζ, w, wp)`
    Coordinate + state conversion at a single point for PIII.  Branch
    convention: `ζ = 2·log(z)` uses Julia's principal `log`.  Callers
    interested in non-principal sheets can shift `ζ` by `4π·im·s` for
    `s ∈ ℤ`.

  - `pIII_ζ_to_z(ζ, w, wp) -> (z, u, up)`
    Inverse of `pIII_z_to_ζ`.

  - `pV_z_to_ζ(z, u, up) -> (ζ, w, wp)`,  `pV_ζ_to_z(ζ, w, wp) -> (z, u, up)`
    Same for PV; principal branch `ζ = log(z)`; non-principal sheets at
    `ζ + 2π·im·s`.

## What is NOT in v1

Per the bead `padetaylor-bvh` (Tier-4 v1 scope):

  - **Non-uniform Stage-1 nodes** (FFW2017...md:67-72): the FFW paper
    notes that PIII / PV pole densities are highly non-uniform across
    the ζ-plane, growing rapidly for `Re ζ ≫ 0`.  FFW uses a node-
    separation function `R(ζ)` that decreases linearly with `Re ζ`.
    The transforms in this module compose with the existing
    `path_network_solve` (uniform target grid) — accuracy in the
    high-`Re ζ` tail will degrade per FFW.  A separate bead is the
    place to add adaptive node placement; the transforms themselves
    are useful as-is for moderate-`|ζ|` solves.

  - **Adaptive Padé step size** (FFW2017...md:74-97): same story.  The
    transforms work with FW's fixed-`h` stepper today.  Adaptive `h`
    is independent — file a bead if a downstream user needs it.

  - **Sheet tracking** (PVI / Phase 14 / bead `padetaylor-grc`): PIII
    and PV branch points are only at `z = 0` and are removed by the
    transform.  No sheet-index bookkeeping needed at this layer.

## Derivation cross-checks

The IC transform is derived once in the worklog (`docs/worklog/017-…`).
Quick recap:

  - `w(ζ) = exp(ζ/2) · u(z(ζ))` for PIII, so `dw/dζ = (1/2) exp(ζ/2) u
    + exp(ζ/2) (du/dz)(dz/dζ) = (1/2) z u + (z²/2) u'`.  At `z = 1`:
    `(w, w') = (u, (u + u')/2)`.  Test CT.1.1 pins this at a specific
    sample.

  - Inverse: `u = w / z`, `u' = (2 w' - z u) / z² = (2 w' - w) / z²`.

  - For PV: `w = u`, `w' = z u'`; inverse `u = w`, `u' = w' / z`.

## References

  - `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
    FFW2017_painleve_riemann_surfaces_preprint.md` §2.1 lines 39-48
    (the transforms + transformed equations), :67-71 (non-uniform
    nodes deferral note).
  - `docs/figure_catalogue.md §5` (Tier-4 acceptance plan).
  - `src/PathNetwork.jl`, `src/Problems.jl` — the downstream consumers.
"""
module CoordTransforms

export pIII_transformed_rhs, pV_transformed_rhs
export pIII_z_to_ζ, pIII_ζ_to_z
export pV_z_to_ζ,   pV_ζ_to_z

# -----------------------------------------------------------------------------
# PIII — z = exp(ζ/2), u(z) = exp(-ζ/2) · w(ζ)
# -----------------------------------------------------------------------------

"""
    pIII_z_to_ζ(z, u, up) -> (ζ, w, wp)

Convert a PIII state `(z, u(z), u'(z))` to the P̃_III ζ-plane state
`(ζ, w(ζ), dw/dζ)` via `ζ = 2 log z`, `w = z u`, `w' = (z u + z² u') / 2`.
Principal-branch `log`; callers wanting sheet `s` add `4π·im·s` to ζ.
"""
function pIII_z_to_ζ(z, u, up)
    ζ  = 2 * log(z)
    w  = z * u
    wp = (z * u + z^2 * up) / 2
    return (ζ, w, wp)
end

"""
    pIII_ζ_to_z(ζ, w, wp) -> (z, u, up)

Inverse of `pIII_z_to_ζ`.  `z = exp(ζ/2)`, `u = w / z`, `u' = (2 w' - w) / z²`.
"""
function pIII_ζ_to_z(ζ, w, wp)
    z  = exp(ζ / 2)
    u  = w / z
    up = (2 * wp - w) / z^2
    return (z, u, up)
end

"""
    pIII_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''

Factory returning the P̃_III RHS closure.  Plug directly into a
`PadeTaylorProblem(rhs, (w₀, w'₀), (ζ_start, ζ_end); order = ...)`.

Requires `w ≠ 0` (else `1/w` blows up).  Throws on `w == 0`-evaluation
is delegated to Julia's IEEE infinity arithmetic — fail-loud per
CLAUDE.md Rule 1 (downstream Padé-Taylor stepper will throw on
non-finite Taylor coefficients).
"""
function pIII_transformed_rhs(α, β, γ, δ)
    return (ζ, w, wp) -> begin
        eζ  = exp(ζ)
        e2ζ = eζ * eζ
        return wp^2 / w + (α * w^2 + γ * w^3 + β * eζ + δ * e2ζ / w) / 4
    end
end

# -----------------------------------------------------------------------------
# PV — z = exp(ζ), u(z) = w(ζ)
# -----------------------------------------------------------------------------

"""
    pV_z_to_ζ(z, u, up) -> (ζ, w, wp)

Convert a PV state `(z, u(z), u'(z))` to the P̃_V ζ-plane state
`(ζ, w(ζ), dw/dζ)` via `ζ = log z`, `w = u`, `w' = z u'`.  Principal-
branch `log`; sheets at `+2π·im·s`.
"""
function pV_z_to_ζ(z, u, up)
    ζ  = log(z)
    w  = u
    wp = z * up
    return (ζ, w, wp)
end

"""
    pV_ζ_to_z(ζ, w, wp) -> (z, u, up)

Inverse of `pV_z_to_ζ`.  `z = exp(ζ)`, `u = w`, `u' = w' / z`.
"""
function pV_ζ_to_z(ζ, w, wp)
    z  = exp(ζ)
    u  = w
    up = wp / z
    return (z, u, up)
end

"""
    pV_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''

Factory returning the P̃_V RHS closure.  Requires `w ≠ 0` (the `β/w`
and `1/(2w)` terms) and `w ≠ 1` (the `1/(w-1)` and `(w+1)/(w-1)`
terms).
"""
function pV_transformed_rhs(α, β, γ, δ)
    return (ζ, w, wp) -> begin
        eζ  = exp(ζ)
        e2ζ = eζ * eζ
        return (1 / (2*w) + 1 / (w - 1)) * wp^2 +
               (w - 1)^2 * (α * w + β / w) +
               γ * eζ * w +
               δ * e2ζ * w * (w + 1) / (w - 1)
    end
end

end # module CoordTransforms
