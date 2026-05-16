# figures/_ffw2017_fig_7_helpers.jl
#
# Pure-Julia helper kernel for `figures/ffw2017_fig_7.jl`.  Defines the
# bilinear-interpolation primitive `bilinear_w` and the load-bearing
# z-plane rendering invariant `zplane_log_modulus`, both kept in this
# separate file so `test/ffw_fig_7_test.jl` can include + exercise the
# EXACT same helper code path that the figure renders with (FF7.1.8
# mutation-prove M5: swap `log10(aw)` for `aw` here and the test goes
# RED).
#
# This file is `include()`d by both `figures/ffw2017_fig_7.jl` (which
# also `using`s CairoMakie + plot routines) and by
# `test/ffw_fig_7_test.jl` (which does NOT pull CairoMakie).  Keep
# this file Makie-free so the test environment can load it cleanly.
#
# Helpers take the ζ-plane render lattice as explicit arguments
# (`xs`, `ys`, `re_lo`, `re_hi`) rather than reading them from
# enclosing scope.  This keeps the kernel callable from any module
# scope (the test's @testset scope; the figure script's top-level
# scope).  Both call sites wrap the helpers in zero-arg shims that
# bind the in-scope grid constants — see the figure script and the
# test file for those shims.

"""
    bilinear_w(W::Matrix{ComplexF64}, ζ::Complex,
               xs::AbstractVector{Float64}, ys::AbstractVector{Float64})
        -> ComplexF64

Bilinear interpolation of `w(ζ)` from the ζ-plane render array `W`
over the lattice `xs × ys`.  Returns `NaN + NaN·im` if the query
point is outside the lattice OR if any of the four surrounding corner
cells is non-finite (NaN propagation: caller sees holes in the
visited tree).

`W` must be shaped `(length(xs), length(ys))`; `xs` and `ys` must be
evenly spaced.
"""
function bilinear_w(W::Matrix{ComplexF64}, ζ::Complex,
                    xs::AbstractVector{<:Real},
                    ys::AbstractVector{<:Real})
    NX = length(xs); NY = length(ys)
    rx, iy_ = real(ζ), imag(ζ)
    (rx < xs[1] || rx > xs[end] ||
     iy_ < ys[1] || iy_ > ys[end]) &&
        return complex(NaN, NaN)
    dx = xs[2] - xs[1]
    dy = ys[2] - ys[1]
    i = clamp(floor(Int, (rx - xs[1]) / dx) + 1, 1, NX - 1)
    j = clamp(floor(Int, (iy_ - ys[1]) / dy) + 1, 1, NY - 1)
    tx = (rx - xs[i]) / dx
    ty = (iy_ - ys[j]) / dy
    w00 = W[i, j];     w10 = W[i+1, j]
    w01 = W[i, j+1];   w11 = W[i+1, j+1]
    (isfinite(real(w00)) && isfinite(real(w10)) &&
     isfinite(real(w01)) && isfinite(real(w11))) ||
        return complex(NaN, NaN)
    return (1 - tx) * (1 - ty) * w00 +
           tx       * (1 - ty) * w10 +
           (1 - tx) * ty       * w01 +
           tx       * ty       * w11
end

"""
    zplane_log_modulus(W::Matrix{ComplexF64}, z::Complex,
                       sheet_lift::Integer,
                       xs::AbstractVector{Float64},
                       ys::AbstractVector{Float64},
                       re_lo::Real, re_hi::Real) -> Float64

z-plane rendering invariant: returns
`log10(abs(bilinear_w(W, ζ_lift, xs, ys)))` where
`ζ_lift = log|z| + i (arg z + 2π · sheet_lift)`.  Returns `NaN` when
`|z|` is outside the annulus `[exp(re_lo), exp(re_hi)]`, when the
lifted Im ζ is outside the rendered ζ-window, or when the bilinear
interp returns non-finite / zero.

This is the LOAD-BEARING differentiator vs Fig 3 (which renders
`arg u(z)`) and Fig 6 (which renders `|u(z)|`).  FFW md:281 prescribes
`log10|u|` for Fig 7 col 3 because PVI's z-plane residues scale with
`|z₀|` (Table 3) — pole spikes and the bulk solution coexist on one
heatmap only when log-scaled.  Mutation M5 (replace `log10(aw)` with
`aw`) flips test FF7.1.8 directly because the test queries this
helper.
"""
function zplane_log_modulus(W::Matrix{ComplexF64}, z::Complex,
                            sheet_lift::Integer,
                            xs::AbstractVector{<:Real},
                            ys::AbstractVector{<:Real},
                            re_lo::Real, re_hi::Real)
    absz = abs(z)
    absz < exp(re_lo) - 0.005 && return NaN
    absz > exp(re_hi) + 0.005 && return NaN
    ζ_lift = complex(log(absz), angle(z) + 2π * sheet_lift)
    (imag(ζ_lift) < ys[1] || imag(ζ_lift) > ys[end]) && return NaN
    w = bilinear_w(W, ζ_lift, xs, ys)
    isfinite(real(w)) || return NaN
    aw = abs(w)
    aw > 0 || return NaN
    return log10(aw)
end
