# figures/ffw2017_fig_5.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, **Figure 5** — the
# tronquée P_III solution on its 1.5-sheet pole-free sector, together
# with the relative-condition-number heatmap that motivates the
# IVP+BVP hybrid algorithm.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:220-265
#         (Figure 5 caption + §3.2 "A tronquée P_III solution").
#
# ## The problem (FFW md:222)
#
# > "The first and second columns of Figure 5 show a tronquée P_III
# >  solution that is pole-free on the region `-3π/4 < arg z < 9π/4`
# >  (column 1), which in the ζ-plane corresponds to the region
# >  `-3π/2 < Im ζ < 9π/2` (column 2).  In the z-plane the asymptotic
# >  behavior of the solution on the pole-free region is `u ~ z^{1/3},
# >  z → ∞`."
#
# Parameters (FFW md:240): `(α, β, γ, δ) = (1, -1/20, 0, -1)`.
#
# The pole-free sector in the ζ-plane is a 6π-wide horizontal strip —
# `1.5 × 4π` PIII sheets (FFW md:222) — far wider than the principal
# strip `(-2π, 2π]`.  FFW md:226: "if our enhanced PFS method is used
# without a BVP solver to compute the tronquée P_III solution, then the
# error is on the order of 10⁻¹"; the BVP step on the smooth sector
# brings the error down (FFW caption md:240) to per-strip `4e-6 / 3e-7
# / 2e-8 / 3e-9 / 4e-8` (from bottom to top across five horizontal
# strips of the sector).  Our figure renders the joint Float64 stack;
# the achieved errors and the gap to FFW's targets are documented in
# `test/ffw_fig_5_test.jl` and worklog 040.
#
# ## The two boundary asymptotic ICs (FFW md:240-244 verbatim)
#
# FFW publish two precomputed asymptotic ICs at points just inside the
# sector's curved boundaries in the z-plane:
#
#     z₁ = 30 e^{( 9π/4 - π/12) i}    (upper-boundary IC)
#     z₂ = 30 e^{(-3π/4 + π/12) i}    (lower-boundary IC)
#
#     u(z₁)  = -2.000735432319 + 2.376177147900 i
#     u(z₂)  =  2.384379236170 - 1.993845650158 i
#     u'(z₁) = -5.939523100e-3 + 3.402038641e-2 i
#     u'(z₂) =  6.050817704e-3 + 3.398020750e-2 i
#
# Both points have `|z| = 30`; in the ζ-frame (`ζ = 2 log z`):
#
#     ζ₁ = 2 log 30 + i · 13π/3 ≈ 6.802 + 13.614 i
#     ζ₂ = 2 log 30 - i ·  4π/3 ≈ 6.802 -  4.189 i
#
# Both fall just INSIDE the FFW sector `Im ζ ∈ (-3π/2, 9π/2)` — they
# are FFW's "two points close to the boundary" (md:232).
#
# ## Sheet-aware asymptotic-IC wrapper
#
# The shipped helper `PadeTaylor.pIII_asymptotic_ic(z; ...)` (FFW
# md:230 series) uses Julia's principal complex `z^{1/3}`.  For the
# upper IC at `arg z = 13π/6`, the principal cube root is
# `30^{1/3} · exp(i π/18)`, but the correct (sector-consistent) branch
# is `30^{1/3} · exp(i · 13π/18)` — a different value entirely.  We
# therefore build a thin sheet-aware wrapper at the figure layer
# (`tronquee_ic_sheet`) that computes `s = z^{1/3}` using the **caller's
# unfolded argument** instead of Julia's principal angle.  The
# coefficients `a_1 = -β/3` and `a_2 = ((4/9) + δ a_1²)/(2δ)` are
# lifted verbatim from `IVPBVPHybrid.pIII_asymptotic_ic`'s docstring
# (ADR-0014 / worklog 039) — no algorithmic change, just a different
# choice of cube-root branch.
#
# This wrapper is needed because the FFW Fig 5 sector spans 6π in
# `Im ζ` — wider than ONE PIII strip — so the principal-branch
# assumption made inside the shipped helper is insufficient.  The
# wrapper is figure-local; the shipped helper is byte-unchanged.
#
# ## The hybrid driver call (using IVPBVPHybrid from ADR-0014)
#
# We call `solve_pole_free_hybrid` on the FFW md:240-244 sector with
# `n_terms = 15`, `N = 50` collocation per slice, `adaptive_tol = 1e-10`
# PFS controller, `glue_tol = 1e-8`, n_slices = 12.  Per worklog 039
# §"What is NOT shipped", the in-tree `pIII_asymptotic_ic` hard-codes
# `a_3, ..., a_N = 0`; the `n_terms = 15` call therefore matches the
# `n_terms = 2` call numerically.  The honest gap to FFW's per-strip
# `4e-6 / 3e-7 / 2e-8 / 3e-9 / 4e-8` (md:240) is recorded in worklog
# 040's §"Empirical sweep".  The figure still renders the correct
# **qualitative** structure — pole-free sector + the FFW md:262
# exponential κ-growth — because the BVP solver achieves spectral-
# floor accuracy GIVEN whatever BCs the PFS walks supply, and the BCs
# themselves are accurate to ~10⁻³ in `u` at `|z| = 30` (FFW
# n_terms = 2 reproduces md:243 to four decimal places, verified in
# worklog 040 §"Ground-truth pin").
#
# ## Two panels
#
# **Panel A — |w(ζ)|**.  Modulus heatmap over the FFW Fig 5 sector
# `Re ζ ∈ [-1.5, 6.802]`, `Im ζ ∈ [-4π/3, 13π/3]` (the two-IC strip
# we cover at v1; the full FFW strip `-3π/2 < Im ζ < 9π/2` would
# require additional asymptotic-IC corner placements).  Should show a
# smooth dim plain in the sector and pole structure outside the FFW
# strip boundaries (rendered for context if the PFS walks reach those
# cells).  The two FFW md:243 asymptotic-IC points are annotated as
# crosses on the upper and lower curved boundaries.
#
# **Panel B — log₁₀(κ_r)**.  Closed-form relative condition number
# from FFW md:262:
#
#     κ_r(ζ) ~ (27/16) · e^{2 Re ζ / 3}    on the pole-free sector
#
# Rendered as a horizontal-stripe heatmap (the formula depends only on
# `Re ζ` — the sector is well-conditioned uniformly in `Im ζ`, but
# CONDITIONING DEGRADES exponentially as `Re ζ` increases).  At the
# FFW IC anchor `Re ζ = 2 log 30 ≈ 6.802`, `κ_r ≈ 157` (FFW md:264
# verbatim).  This panel is the "condition-number heatmap" that FFW
# Fig 5 column 3 reports.  An EMPIRICAL κ_r evaluation at sample
# points (using FFW eq. (7) md:254 with the hybrid-driver-evaluated
# w, w', w'') is logged at script bottom for cross-validation against
# the closed form.
#
# ## Acceptance
#
# Per `docs/figure_catalogue.md §5` (Fig 5 row, T4): the figure script
# runs with `n_terms = 15` (per the bead spec — though the helper's
# v1 caveat means this is equivalent to `n_terms = 2`), `N = 50`,
# `adaptive_tol = 1e-10`, `glue_tol = 1e-8`.  `test/ffw_fig_5_test.jl`
# `FF5.1.*` pins: IC round-trip exact (FF5.1.1), pole-free sector
# bounded |u| inside (FF5.1.2), per-strip self-cross-check (FF5.1.3-5),
# condition-number formula κ_r ≈ 157 at z = 30 (FF5.1.6), κ-monotonicity
# (FF5.1.7), and n_terms-convergence behavior (FF5.1.8 — documenting
# the helper's v1 cap of a_2).

using PadeTaylor
using PadeTaylor.CoordTransforms: pIII_transformed_rhs, pIII_z_to_ζ, pIII_ζ_to_z
using CairoMakie
using Printf

# =====================================================================
# Parameters (FFW md:240 verbatim)
# =====================================================================
const α, β, γ, δ = 1.0, -1/20, 0.0, -1.0

# Two FFW md:243 asymptotic-IC points in z-frame.  Verbatim copies.
const ARG_Z1 = 13π/6     # = 9π/4 - π/12, upper-boundary IC
const ARG_Z2 = -2π/3     # = -3π/4 + π/12, lower-boundary IC
const Z1_MOD = 30.0
const Z2_MOD = 30.0
const Z1     = Z1_MOD * cis(ARG_Z1)
const Z2     = Z2_MOD * cis(ARG_Z2)

# FFW md:243 published values (reference oracle for tests).
const U_Z1_FFW  = -2.000735432319 + 2.376177147900im
const U_Z2_FFW  =  2.384379236170 - 1.993845650158im
const UP_Z1_FFW = -5.939523100e-3 + 3.402038641e-2im
const UP_Z2_FFW =  6.050817704e-3 + 3.398020750e-2im

# In the ζ-frame: ζ = 2 log z; the upper IC at unfolded arg z = 13π/6
# has Im ζ = 13π/3; the lower at -4π/3.  Re ζ = 2 log 30 for both.
#
# **v1 sector choice (honest restriction)**.  The FULL FFW Fig 5 sector
# spans 6π in Im ζ — 1.5 PIII sheets — making the asymptotic-series
# cube root sheet-dependent across the BVP collocation nodes.  Our
# shipped `pIII_asymptotic_ic` uses Julia's principal `z^{1/3}`, which
# is valid only on ONE PIII sheet (Im ζ ∈ (-2π, 2π]).  Rather than
# build a per-node sheet-aware adapter (which would require modifying
# the IVPBVPHybrid driver's API — explicitly out of scope per CLAUDE.md
# Rule 6 / the bead's "compose, don't refactor" mandate), we render a
# PRINCIPAL-SHEET slice of the FFW sector.  Panel A's BVP solve lives
# inside this principal slice; Panel B's condition-number heatmap
# (closed form, Im-independent) extends across the full FFW range.
const RE_ANCHOR = 2 * log(Z1_MOD)            # ≈ 6.802
# Principal-sheet sub-strip: Im ζ ∈ [-3π/2 + 0.05, 2π - 0.05].
# Covers the lower-FFW-IC at arg z = -2π/3 (Im ζ = -4π/3 ≈ -4.19) —
# inside [-3π/2 + ε, 2π - ε].  The upper-FFW-IC at arg z = 13π/6
# (Im ζ = 13π/3 ≈ 13.61) is OUTSIDE this principal-sheet strip and is
# shown as an annotation only.
const IM_HI     = 2π - 0.05
const IM_LO     = -3π/2 + 0.05
const RE_EXTENT = 3.0                        # BVP slab Re ζ ∈ [3.8, 6.8]
                                              # (|z| ∈ [6.7, 30]) — inside the
                                              # asymptotic regime where the
                                              # n_terms=2 IC is reliable

# Full FFW Fig 5 sector range (for Panel B's closed-form κ heatmap).
const IM_HI_FFW = 9π/2
const IM_LO_FFW = -3π/2

# Hybrid-driver tuning (per task brief).  `n_terms = 15` is the
# bead-spec instruction; per worklog 039 §"What is NOT shipped" the
# v1 helper hard-codes `a_3 ... a_N = 0`, so this is numerically
# equivalent to `n_terms = 2` — documented in the empirical sweep.
const N_TERMS    = 15
const N_BVP_COLL = 50
const N_SLICES   = 12
const ADAPT_TOL  = 1.0e-10
const GLUE_TOL   = 1.0e-8
const K_CONS     = 1.0e-3

const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_5.png")

# =====================================================================
# Sheet-aware tronquée P_III asymptotic IC.
#
# The shipped `PadeTaylor.pIII_asymptotic_ic` uses Julia's principal
# `z^(1/3)`; for points with unfolded `arg z ∉ (-π, π]` this gives
# the wrong sheet.  The FFW Fig 5 upper IC at `arg z = 13π/6` needs
# the +2π sheet of the cube root.  We rebuild the helper's series
# (md:230) here using the CALLER-SUPPLIED unfolded argument; the
# coefficients `a_1 = -β/3`, `a_2 = ((4/9) + δ a_1²) / (2δ)` are
# verbatim from the shipped helper's docstring (ADR-0014).
#
# Returns `(u, u')` with `u(z) ≈ s · (1 + a_1 s^{-2} + a_2 s^{-4})`
# and the corresponding derivative; valid for `|z| ≫ 1` on any single
# sheet of `z^{1/3}` consistent with `arg_z`.
# =====================================================================
function tronquee_ic_sheet(z_modulus::Real, arg_z::Real;
                            n_terms::Integer = N_TERMS,
                            β::Real = β, δ::Real = δ)
    # Coefficients.  The shipped helper only knows a_1 and a_2 closed
    # form; a_3+ are zero (worklog 039 §"What is NOT shipped").
    a1 = -float(β) / 3
    a2 = n_terms ≥ 2 ? ((4.0/9.0) + float(δ) * a1^2) / (2 * float(δ)) : 0.0
    # Cube root with the caller's unfolded argument.
    s_mag = z_modulus^(1/3)
    s_arg = arg_z / 3
    s     = complex(s_mag * cos(s_arg), s_mag * sin(s_arg))
    # For low |z| the divergent asymptotic series gives garbage; use
    # the leading-order `u ≈ z^{1/3}` as the seed instead.  Threshold
    # `|z| < 2` is conservative — at |z| = 2 the a_1 z^{-1/3} correction
    # is ~10⁻² and a_2 z^{-1} is ~10⁻¹; below that the series is
    # unreliable as a Newton initial guess.
    if z_modulus < 2.0
        u  = s
        up = inv(3 * s^2)
        return (u, up)
    end
    s_inv2 = inv(s^2)
    u_sum  = one(s) + a1 * s_inv2 + a2 * s_inv2^2
    up_sum = one(s) - a1 * s_inv2 - 3 * a2 * s_inv2^2
    u  = s * u_sum
    up = up_sum / (3 * s^2)
    return (u, up)
end

# Principal-sheet adapter for the hybrid driver.  The sub-strip we
# render lives entirely in `Im ζ ∈ (-2π, 2π]` (mod the 0.05 nudges
# that keep us off the strip boundaries), so Julia's principal
# `z^{1/3}` is valid throughout — no sheet bookkeeping required.
# The closure simply forwards z's principal argument to
# `tronquee_ic_sheet`.
function make_principal_sheet_ic(sector::NamedTuple;
                                  n_terms = N_TERMS, β = β, δ = δ)
    return (z::Number) -> begin
        return tronquee_ic_sheet(abs(z), angle(z);
                                  n_terms = n_terms, β = β, δ = δ)
    end
end

# =====================================================================
# Build the PainleveProblem and call the hybrid driver.
#
# Per the bead spec: use `IVPBVPHybrid.solve_pole_free_hybrid` from
# ADR-0014.  The IC at z₂ (lower boundary) is the seed for the
# `PainleveProblem`; the driver's internal `_pfs_ray_walk` builds
# separate walks anchored at the upper and lower asymptotic-IC corners.
# =====================================================================
sector = (im_lo = IM_LO, im_hi = IM_HI,
           re_anchor = RE_ANCHOR, re_extent = RE_EXTENT)

asymptotic_ic_fn = make_principal_sheet_ic(sector;
                                            n_terms = N_TERMS, β = β, δ = δ)

# Build the PainleveProblem seeded at Z2 (the driver uses pp's frame
# maps for coordinate round-trips; the underlying problem.f is the
# only RHS used during the PFS walks).  Z2 has arg z = -2π/3 ∈ (-π, π],
# in the principal sheet — `asymptotic_ic_fn(Z2)` returns the
# sheet-correct value.
u0_seed, up0_seed = asymptotic_ic_fn(Z2)
pp = PainleveProblem(:III; α = α, β = β, γ = γ, δ = δ,
                       u0 = u0_seed, up0 = up0_seed,
                       zspan = (Z2, Z2 + 1.0), order = 30)

@printf("FFW 2017 Fig 5: (α, β, γ, δ) = (%.3f, %.4f, %.3f, %.3f)\n",
        α, β, γ, δ)
@printf("                |z| = %.1f for both ICs\n", Z1_MOD)
@printf("                Z1 (upper) arg z = %.4f rad (= 13π/6)\n", ARG_Z1)
@printf("                Z2 (lower) arg z = %.4f rad (= -2π/3)\n", ARG_Z2)
@printf("                ζ sector: Im ζ ∈ [%.4f, %.4f] (width %.4f ≈ %.2fπ)\n",
        IM_LO, IM_HI, IM_HI - IM_LO, (IM_HI - IM_LO) / π)
@printf("                Re ζ ∈ [%.4f, %.4f]\n",
        RE_ANCHOR - RE_EXTENT, RE_ANCHOR)
@printf("                tuning: n_terms = %d, N = %d, n_slices = %d, adaptive_tol = %.0e\n",
        N_TERMS, N_BVP_COLL, N_SLICES, ADAPT_TOL)

# Cross-check the sheet-aware IC formula against FFW md:243.  The
# wrapper takes (|z|, arg z) and uses the caller's unfolded argument
# for `z^{1/3}` — so we CAN reproduce FFW's published IC values at
# both z₁ (upper IC, arg = 13π/6) and z₂ (lower IC, arg = -2π/3) by
# evaluating with the FFW unfolded arguments.  The gap |Δ| ~ 10⁻³
# is the n_terms=2 truncation error vs FFW's optimally-truncated
# (n_terms ≥ 15) series.  Documented in worklog 040 §"Ground truth
# pin".
u1_computed, up1_computed = tronquee_ic_sheet(Z1_MOD, ARG_Z1;
                                                n_terms = N_TERMS,
                                                β = β, δ = δ)
u2_computed, up2_computed = tronquee_ic_sheet(Z2_MOD, ARG_Z2;
                                                n_terms = N_TERMS,
                                                β = β, δ = δ)
@printf("                tronquee(|z|=30, arg=13π/6) = %.6f %+.6f i  vs FFW %.6f %+.6f i  (|Δ| = %.2e)\n",
        real(u1_computed), imag(u1_computed),
        real(U_Z1_FFW), imag(U_Z1_FFW),
        abs(u1_computed - U_Z1_FFW))
@printf("                tronquee(|z|=30, arg=-2π/3) = %.6f %+.6f i  vs FFW %.6f %+.6f i  (|Δ| = %.2e)\n",
        real(u2_computed), imag(u2_computed),
        real(U_Z2_FFW), imag(U_Z2_FFW),
        abs(u2_computed - U_Z2_FFW))

# Hybrid solve.
t0 = time()
sol = solve_pole_free_hybrid(pp, sector, asymptotic_ic_fn;
        pfs_kwargs = (; h = 0.4,
                      step_size_policy = :adaptive_ffw,
                      adaptive_tol = ADAPT_TOL,
                      k_conservative = K_CONS,
                      max_rescales = 50,
                      max_steps_per_target = 800),
        bvp_kwargs = (; N = N_BVP_COLL, tol = 1e-10, maxiter = 30),
        n_slices   = N_SLICES,
        glue_tol   = GLUE_TOL)
@printf("Hybrid solve in %.2f s; %d BVP slices.\n",
        time() - t0, length(sol.bvp_slices))
for (i, slice) in enumerate(sol.bvp_slices)
    @printf("  slice %2d  Re ζ = %.4f  iters = %d  res = %.2e\n",
            i, real(slice.z_a), slice.iterations, slice.residual_inf)
end

# =====================================================================
# Panel A: |w(ζ)| heatmap over the BVP sector.
#
# Build a uniform Re × Im grid covering [re_lo, re_hi] × [im_lo, im_hi].
# For each cell, evaluate `sol(ζ)` — the IVPBVPSolution's callable
# (BVP linear-interp between bracketing slices in Re ζ).  Outside the
# sector we render NaN (transparent).
# =====================================================================
const NX = 220
const NY = 320
const RE_LO  = RE_ANCHOR - RE_EXTENT   # = -1.5
const RE_HI  = RE_ANCHOR                # ≈ 6.802

xs_pA = collect(range(RE_LO + 0.005, RE_HI - 0.005; length = NX))
ys_pA = collect(range(IM_LO + 0.01, IM_HI - 0.01; length = NY))

absW_pA = fill(NaN, NX, NY)
for jx in 1:NY, ix in 1:NX
    ζ = complex(xs_pA[ix], ys_pA[jx])
    try
        w, _ = sol(ζ)
        absW_pA[ix, jx] = isfinite(w) ? abs(w) : NaN
    catch e
        # DomainError outside sector → leave NaN.
        absW_pA[ix, jx] = NaN
    end
end
n_finite_pA = count(isfinite, absW_pA)
@printf("Panel A: %d / %d cells finite (%.1f%%)\n",
        n_finite_pA, length(absW_pA), 100 * n_finite_pA / length(absW_pA))

# =====================================================================
# Panel B: log₁₀(κ_r) heatmap via closed form (FFW md:262).
#
# κ_r(ζ) = (27/16) · exp(2 Re ζ / 3)   on the pole-free sector.
# We render on the FULL FFW Fig 5 ζ-window: Re ζ matching Panel A's
# BVP rectangle but Im ζ spanning the full FFW sector
# `(-3π/2, 9π/2)`.  Since the closed-form κ_r is Im ζ-independent,
# this is one-dimensional and renders as horizontal stripes.
# =====================================================================
const NY_FFW = 480
xs_pB = xs_pA
ys_pB = collect(range(IM_LO_FFW + 0.01, IM_HI_FFW - 0.01; length = NY_FFW))

κ_pB = fill(NaN, NX, NY_FFW)
for jx in 1:NY_FFW, ix in 1:NX
    re = xs_pB[ix]
    κ_pB[ix, jx] = (27.0/16.0) * exp(2 * re / 3)
end
log_κ_pB = log10.(κ_pB)
@printf("Panel B: log10(κ_r) range [%.3f, %.3f]  (FFW md:264: max ≈ log10(157) = 2.20)\n",
        minimum(log_κ_pB), maximum(log_κ_pB))

# =====================================================================
# Empirical κ_r at sample points — cross-validation against FFW
# eq. (7) md:254.  We pick a handful of interior ζ points, evaluate
# w, w', and approximate w'' via finite-difference along Im ζ;
# then κ_r = |1/w''| · |-(w')²/w + (1/4)(2αw² + 3γw³ - δe^{2ζ}/w)|.
# =====================================================================
function empirical_κ(sol, ζ; h = 0.002)
    w0, wp0 = sol(ζ)
    wp_p, _ = sol(ζ + complex(0, h))
    wp_m, _ = sol(ζ - complex(0, h))
    # Approximate w''(ζ) via FD on w' along Im axis.
    wpp = (wp_p - wp_m) / (2 * im * h)
    eζ  = exp(ζ); e2ζ = eζ * eζ
    numer = -wp0^2 / w0 + (2 * α * w0^2 + 3 * γ * w0^3 - δ * e2ζ / w0) / 4
    return abs(numer / wpp)
end

sample_ζs = ComplexF64[
    complex(0.0, 0.0),       # mid-sector, low Re
    complex(2.0, 0.0),
    complex(4.0, 0.0),
    complex(6.0, 0.0),       # high Re, near FFW IC
    complex(2.0, 2.0),
    complex(4.0, -1.0),
]
@printf("\nEmpirical-vs-formula κ_r at %d sample points:\n", length(sample_ζs))
@printf("  %-30s  %-12s  %-12s  ratio\n",
        "ζ",
        "κ_r (formula)",
        "κ_r (empirical)")
for ζs in sample_ζs
    κ_formula = (27.0/16.0) * exp(2 * real(ζs) / 3)
    try
        κ_emp = empirical_κ(sol, ζs)
        @printf("  (%6.3f, %6.3f)              %.4e    %.4e    %.3f\n",
                real(ζs), imag(ζs), κ_formula, κ_emp, κ_emp / κ_formula)
    catch e
        @printf("  (%6.3f, %6.3f)              %.4e    --(domain)         --\n",
                real(ζs), imag(ζs), κ_formula)
    end
end

# =====================================================================
# Render.  Two panels, ~1500 × 1100 px target.
# =====================================================================
# |w| cap for Panel A.  On the pole-free sector, |w| ~ |z|^{4/3};
# at |z| = exp(Re ζ / 2), the maximum is |w|_max ≈ exp(2 Re ζ / 3) ≈ 93
# at Re ζ = 2 log 30.  Cap at 120 to keep the heatmap dynamic-range
# legible across the rectangle.
const U_CAP = 120.0
const ABS_LO_PA = 0.0

fig = Figure(size = (1500, 1100))
Label(fig[0, 1:2],
      "FFW 2017 Fig 5 — tronquée P_III, (α, β, γ, δ) = (1, −1/20, 0, −1) — pole-free sector −3π/2 < Im ζ < 9π/2 (Panel A: principal-sheet sub-strip; Panel B: full sector via closed form)";
      fontsize = 14, padding = (0, 0, 10, 10))

# ---- Panel A: |w(ζ)| heatmap ----------------------------------------
axA_subtitle = @sprintf("Principal-sheet BVP rectangle: Re ζ ∈ [%.2f, %.2f], Im ζ ∈ [%.2f, %.2f]; width %.2fπ",
                         RE_ANCHOR - RE_EXTENT, RE_ANCHOR, IM_LO, IM_HI, (IM_HI - IM_LO) / π)
axA = Axis(fig[1, 1];
           xlabel = "Re ζ", ylabel = "Im ζ",
           title  = "A.  |w(ζ)| (tronquée P_III, BVP sector)",
           subtitle = axA_subtitle,
           titlesize = 13, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (RE_LO - 0.1, RE_HI + 0.1, IM_LO - 0.3, IM_HI + 0.3))

clipped_absW = clamp.(absW_pA, ABS_LO_PA, U_CAP)
heatmap!(axA, xs_pA, ys_pA, clipped_absW;
         colormap = :viridis, colorrange = (ABS_LO_PA, U_CAP),
         nan_color = :transparent)

# Dashed strip boundaries at the FFW md:240 horizontal strip dividers
# falling inside this sub-strip.  FFW reports 5 strips bottom→top with
# per-strip errors `4e-6, 3e-7, 2e-8, 3e-9, 4e-8` (md:240); strip
# width is 6π/5 = 1.2π.  Only ~3 dividers fall in our principal sheet.
for k in 1:4
    h = -3π/2 + k * (6π / 5)
    if IM_LO ≤ h ≤ IM_HI
        lines!(axA, [RE_LO, RE_HI], [h, h];
               color = (:white, 0.5), linewidth = 0.6, linestyle = :dash)
    end
end

# Driver-auto-placed boundary anchors (the actual asymptotic-IC
# sources for our principal-sheet BVP rectangle).
ε_drv = 0.001 * (IM_HI - IM_LO)
scatter!(axA, [Float64(RE_ANCHOR), Float64(RE_ANCHOR)],
              [Float64(IM_HI - ε_drv), Float64(IM_LO + ε_drv)];
        color = :red, marker = :xcross, markersize = 14, strokewidth = 1.5,
        strokecolor = :white)
text!(axA, Float64(RE_ANCHOR - 0.7), Float64(IM_HI - 0.4);
      text = "ζ_top (driver anchor)", color = :white, fontsize = 9)
text!(axA, Float64(RE_ANCHOR - 0.7), Float64(IM_LO + 0.25);
      text = "ζ_bot (driver anchor)", color = :white, fontsize = 9)

Colorbar(fig[1, 1, Right()],
         label = "|w(ζ)|  (capped at $U_CAP)",
         colormap = :viridis, colorrange = (ABS_LO_PA, U_CAP),
         labelsize = 10)

# ---- Panel B: log10(κ_r) heatmap (closed form, FFW md:262) ----------
axB_subtitle = @sprintf("FFW md:262 formula; max ≈ %.2f at Re ζ = %.2f (κ_r ≈ %.1f at z = 30, FFW md:264)",
                         maximum(log_κ_pB), RE_ANCHOR,
                         (27/16) * exp(2 * RE_ANCHOR / 3))
axB = Axis(fig[1, 2];
           xlabel = "Re ζ", ylabel = "Im ζ",
           title  = "B.  log₁₀ κ_r(ζ) = log₁₀((27/16) e^{2 Re ζ / 3})  — full FFW sector",
           subtitle = axB_subtitle,
           titlesize = 13, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (RE_LO - 0.1, RE_HI + 0.1, IM_LO_FFW - 0.3, IM_HI_FFW + 0.3))
heatmap!(axB, xs_pB, ys_pB, log_κ_pB;
         colormap = :inferno,
         colorrange = (minimum(log_κ_pB), maximum(log_κ_pB)))

# All FFW 5-strip boundaries on B.  Per md:240: bottom→top per-strip
# errors `4e-6, 3e-7, 2e-8, 3e-9, 4e-8`.  Strip width = 6π/5 = 1.2π.
for k in 1:4
    h = -3π/2 + k * (6π / 5)
    lines!(axB, [RE_LO, RE_HI], [h, h];
           color = (:white, 0.5), linewidth = 0.6, linestyle = :dash)
end

# Mark FFW md:243's two ASYMPTOTIC IC corners on Panel B (these are
# the exact FFW Fig 5 source points, not auto-placed; one is on a
# non-principal PIII sheet, so absent from Panel A's BVP solve).
ζ1_ffw = complex(RE_ANCHOR, 2 * ARG_Z1)   # ≈ 6.802 + 13.614 i
ζ2_ffw = complex(RE_ANCHOR, 2 * ARG_Z2)   # ≈ 6.802 -  4.189 i
scatter!(axB, [Float64(real(ζ1_ffw)), Float64(real(ζ2_ffw))],
              [Float64(imag(ζ1_ffw)), Float64(imag(ζ2_ffw))];
        color = :white, marker = :xcross, markersize = 14, strokewidth = 1.5,
        strokecolor = :black)
text!(axB, Float64(real(ζ1_ffw) - 1.2), Float64(imag(ζ1_ffw) + 0.25);
      text = "ζ₁ FFW (z₁ = 30 e^{(9π/4-π/12)i})",
      color = :white, fontsize = 9)
text!(axB, Float64(real(ζ2_ffw) - 1.2), Float64(imag(ζ2_ffw) - 0.7);
      text = "ζ₂ FFW (z₂ = 30 e^{(-3π/4+π/12)i})",
      color = :white, fontsize = 9)

Colorbar(fig[1, 2, Right()],
         label = "log₁₀ κ_r",
         colormap = :inferno,
         colorrange = (minimum(log_κ_pB), maximum(log_κ_pB)),
         labelsize = 10)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("\n  wrote %s\n", OUTPNG)
