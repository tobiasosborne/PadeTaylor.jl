"""
    PadeTaylor.EdgeDetector

Tier-2.5 pole-field edge-detection module (Phase 12.5, bead `padetaylor-c2p`)
— the 5-point Laplacian classifier from Fornberg–Weideman 2011 §3.2.2
(`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:202-208`).

This module is a **diagnostic primitive**, not a step in any solver
pipeline: callers supply a 2D array of `u(z)` values on an equispaced
complex lattice (typically the dense output of `PathNetwork.path_network_solve`
over a Cartesian grid), and the detector reports either the raw discrete
Laplacian residual `Δu` (eq. 3.3) or a Boolean bitmap classifying each
interior cell as "pole-field" or "smooth" by thresholding `log₁₀|Δu|`
against the level curve drawn in FW Fig. 3.3.

## The algorithm (FW 2011 §3.2.2 verbatim)

Writing `z = x + iy`, every meromorphic `u(z)` obeys Laplace's equation
`∂²u/∂x² + ∂²u/∂y² = 0` away from its poles.  The discrete 5-point
stencil

```
Δu ≈ ([ 1  ]
      [1 -4 1] · u) / h²
     ([ 1  ])
```

is an `O(h²)` approximation to `∇²u`.  By the harmonicity above, the
analytic value is zero away from poles; the discrete residual is small
(`~|u'''| · h²`) in smooth regions and large near singularities (where
`u` is *not* harmonic — analyticity breaks at a pole).  FW Fig. 3.3
displays `log₁₀|Δu|` as a colour field over the lattice; the contour at
level `0.001` (FW's empirical choice) separates pole-field cells from
smooth cells.

The level is in units of `log₁₀|Δu|`, i.e., the comparison is
`log₁₀|Δu| > level`, equivalently `|Δu| > 10^level`.  Default `level =
0.001` matches FW's published contour and selects cells where `|Δu| >
10^0.001 ≈ 1.0023` — a hair above the natural separation between
"essentially zero" (analytic smooth, `|Δu| ~ 10⁻¹⁰` for typical grids)
and "non-zero" (pole proximity, `|Δu| ≫ 1`).

## What the bead description got slightly wrong

`padetaylor-c2p`'s description reads "|∇²u|/h² > 0.001 threshold (FW's
empirical level)".  This conflates two things: FW's eq. (3.3) already
defines `Δu ≈ (stencil · u) / h²` (the `/h²` is inside `Δu`), and the
0.001 level is on `log₁₀|Δu|`, not on bare `|Δu|`.  See FW2011...md:208
verbatim: "Fig. 3.3 shows `log₁₀|Δu|` ... it is easy to select a
contour level (here 0.001)".  We follow the paper.

## v1 admission about manual classification

FW 2011 line 401 verbatim: "we implemented this 'manually' after having
inspected, in preliminary test calculations, where the pole field edges
appeared."  i.e., the published runs combined visual inspection with
the automatic detector.  Our v1 ships the **diagnostic primitive**;
downstream callers (Phase 12 v2, bead `padetaylor-k31`) compose it with
the dispatcher.  The user can also call it interactively to verify the
threshold is appropriate for their problem before relying on it.

## Boundary handling

Interior cells (rows 2..end-1, columns 2..end-1) receive the full
5-point stencil.  Boundary cells (first/last row, first/last column)
have no neighbours on one side; the residual there is undefined.  We
mark boundary entries of `laplacian_residual`'s output as `NaN + NaN*im`
and boundary entries of `pole_field_mask`'s output as `false`.  Callers
who need full-lattice classification should pad their input grid by one
cell on each side before calling.

## Public API

  - `laplacian_residual(u_grid::AbstractMatrix{Complex{T}}, h::Real) ->
    Matrix{Complex{T}}` — eq. (3.3) verbatim.  Returns a matrix of the
    same shape as `u_grid`; boundary entries are `NaN + NaN·im`.

  - `pole_field_mask(u_grid::AbstractMatrix{Complex{T}}, h::Real;
    level = 0.001) -> BitMatrix` — convenience wrapper that returns
    `log₁₀|Δu| > level` per cell; boundary entries `false`.

  - `pole_field_mask(residual::AbstractMatrix{Complex{T}}; level = 0.001)`
    — same, for pre-computed `Δu` from `laplacian_residual`.

## Fail-fast contract (CLAUDE.md Rule 1)

Throws `ArgumentError` on:
  - `size(u_grid)` with either dimension `< 3` (no interior cells).
  - `h ≤ 0`.

NaN/Inf entries in `u_grid` propagate to `Δu` per IEEE 754; no special
handling.  Callers concerned about NaNs upstream of the detector should
check before invoking.

## References

  - Fornberg & Weideman, *A numerical methodology for the Painlevé
    equations*, J. Comput. Phys. 230 (2011) 5957–5973, §3.2.2 +
    Fig. 3.3 + line 401.  `references/markdown/FW2011_*.md:202-208,
    395-410`.
"""
module EdgeDetector

export laplacian_residual, pole_field_mask

"""
    laplacian_residual(u_grid::AbstractMatrix{Complex{T}}, h::Real) where {T<:AbstractFloat}

Apply the 5-point cross stencil `[1, 1, -4, 1, 1] / h²` to `u_grid` per
FW 2011 eq. (3.3).  Returns a matrix the same shape as `u_grid` with
interior cells populated and boundary cells set to `NaN + NaN·im`.

The stencil indexing convention: `u_grid[i, j]` is at lattice
coordinate `z = x_j + i·y_i` if the matrix rows index `y` and columns
index `x` (i.e., the conventional `imshow` orientation).  The stencil
is rotationally symmetric, so the convention is for the caller's
convenience and does not affect the residual values.
"""
function laplacian_residual(
    u_grid::AbstractMatrix{Complex{T}},
    h::Real,
) where {T<:AbstractFloat}
    nrow, ncol = size(u_grid)
    nrow ≥ 3 && ncol ≥ 3 || throw(ArgumentError(
        "u_grid must be at least 3x3 (got $(size(u_grid))); " *
        "detail: stencil needs one cell of padding on each side"))
    h > 0 || throw(ArgumentError(
        "h must be positive (got $h); " *
        "detail: stencil denominator is h²"))

    h_T   = T(h)
    inv_h2 = inv(h_T^2)
    Δu    = fill(Complex{T}(NaN, NaN), nrow, ncol)
    @inbounds for j in 2:(ncol-1), i in 2:(nrow-1)
        Δu[i, j] = inv_h2 * (
            u_grid[i-1, j] + u_grid[i+1, j] +
            u_grid[i,   j-1] + u_grid[i,   j+1] -
            4 * u_grid[i, j]
        )
    end
    return Δu
end

"""
    pole_field_mask(u_grid::AbstractMatrix{Complex{T}}, h::Real; level = 0.001)

Bitmap classifier: returns a `BitMatrix` the same shape as `u_grid`
with `mask[i, j] = log₁₀|Δu[i, j]| > level` on interior cells and
`false` on boundary cells.

Default `level = 0.001` matches FW 2011 Fig. 3.3 (the level curve
labelled in the figure).  Pass `level = -3` (say) to detect a much
weaker non-harmonicity if the problem demands.

See module docstring for the FW reference and the boundary-handling
convention.
"""
function pole_field_mask(
    u_grid::AbstractMatrix{Complex{T}},
    h::Real;
    level::Real = 0.001,
) where {T<:AbstractFloat}
    Δu = laplacian_residual(u_grid, h)
    return pole_field_mask(Δu; level = level)
end

function pole_field_mask(
    Δu::AbstractMatrix{Complex{T}};
    level::Real = 0.001,
) where {T<:AbstractFloat}
    nrow, ncol = size(Δu)
    mask = falses(nrow, ncol)
    level_T = T(level)
    @inbounds for j in 1:ncol, i in 1:nrow
        z = Δu[i, j]
        # Boundary cells were set to NaN by `laplacian_residual`; we
        # leave the mask `false` there.  Other NaNs (e.g. from upstream
        # NaN in u_grid) also propagate to `false` via the `isnan` check
        # — defensive but cheap.
        isnan(real(z)) && continue
        if log10(abs(z)) > level_T
            mask[i, j] = true
        end
    end
    return mask
end

end # module EdgeDetector
