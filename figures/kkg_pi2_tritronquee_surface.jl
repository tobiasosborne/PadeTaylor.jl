# figures/kkg_pi2_tritronquee_surface.jl
#
# THE HEADLINE FIGURE of the v0.2 vector solver (bead `padetaylor-0ln.28`,
# v0.2 plan row "V8b+", F4).  It renders the full complex-plane surface of
# the P_I^(2) tritronquée `V_0(x, 0)` over the square `[-20,20]²`,
# reproducing Kapaev–Klein–Grava (KKG) 2015 Figs 7.4 (`Re V_0`) and 7.5
# (`Im V_0`) — as **both 2D heatmaps and 3D surfaces**, in the visual
# idiom of the Fornberg–Weideman pole-field figures (`figutil.jl`).
#
# Where the V8b figure (`kkg_pi2_tritronquee_pole_field.jl`) shipped only
# a real-axis BVP slice + a pole *scatter*, F4 is the full Re/Im surface
# KKG actually plot: the deferred-corner promise of the V8b README is
# discharged here.
#
# ## The equation, the solution
#
# P_I^(2) is the fourth-order member of the Painlevé-I hierarchy
# (KKG 2015 eq. (1.1); pillar C §1,
# `docs/v0p2_pillarC_painleve_hierarchy_findings.md`), here at `t = 0`:
#
#     u'''' + 10 u'^2 + 20 u u'' + 40 (u^3 - 6 t u + 6 x) = 0.
#
# Its **tritronquée** `V_0` (KKG 2015 §7) is the distinguished separatrix
# solution: pole-free on a wide angular sector (~270°, centred on the
# negative real axis) and pole-rich only in the wedge straddling
# `arg x ≈ 0` (KKG eq. (1.4); pillar C §4).  On the negative real axis it
# is real, positive and monotone — the **ridge along the negative real
# axis is positive** (`u → +(6|x|)^{1/3}`, KKG eq. (1.3)); a negative
# ridge would signal the sign-bug regression guarded by `test/`.
#
# ## How the surface is computed — the F3 triple-method kernel
#
# Every load-bearing computation lives in the Makie-free kernel
# `figures/_kkg_pi2_surface_helpers.jl` (F3, bead `padetaylor-0ln.33`;
# the wedge half re-resolved by B3, bead `padetaylor-0ln.37.7`), shared
# verbatim with its acceptance test `figures/test_kkg_pi2_surface.jl`.
# This script ONLY renders — it `include`s the kernel and calls
# `kkg_pi2_surface()`.
#
# The kernel stitches two regions onto one `[-20,20]²` Cartesian grid:
#
#   * **Region 1 — the pole-free sector** (`arg x ∈ (36°,324°)`): `V_0`
#     is analytic, so `Re u`/`Im u` are harmonic and are computed THREE
#     algorithmically-disjoint ways — (1) a fan of radial `vector_bvp_solve`
#     rays, (2) the in-house 2D-Chebyshev Laplace solve on the conformal
#     `w = log x` rectangle, (3) the Gridap FEM Laplace solve on the same
#     rectangle — and **majority-voted per grid point** (the 3-sample
#     median: an outlier voter is discarded, the two that concur win;
#     ADR-0024).  The `spread` map records the max pairwise disagreement.
#
#   * **Region 2 — the pole-rich wedge** (`|arg x| < 36°`): the
#     validated pole *field* + an honest partial `|u|` surface underlay
#     (ADR-0025 Amendment 2).  The A2 tractability probe proved a filled
#     honest wedge surface numerically unreachable — honest coverage
#     saturates at ~8-18 % — so the wedge panel is rescoped: the B2
#     `:max_q_root` adaptive walk threads the B3 extended fan through
#     the whole wedge to `|x| = 20`, `extract_poles_shared_q`
#     (`min_support ≥ 2`, VC-6) extracts the validated pole *locations*
#     (the FW 2011 Fig 4.7 idiom, the primary deliverable), and the F2
#     Stage-2 fill — run `extrapolate = false`, B1 true-radius gate —
#     supplies an honest `|u|` underlay only where the verified node
#     discs cover.  The rest of the wedge is `NaN` (an honest gap).
#
#   The `±36°` Stokes-line strips are `NaN`-masked (the asymptotic seed
#   and both region solvers degrade there).  Cells outside the `|x| ≤ 20`
#   disc are `NaN`.  `NaN` cells render in a neutral grey.
#
# ## Rendering choices (this script's only decisions)
#
#   * **Layout.**  A 2×2 panel grid: Re-heatmap | Re-3D / Im-heatmap |
#     Im-3D.  Re row = KKG Fig 7.4; Im row = KKG Fig 7.5.
#
#   * **Display clamp.**  In the wedge `|u|` spikes to `O(10³)`; an
#     un-clamped colour/`z` range would wash the smooth ~270° sector to a
#     flat single colour.  Re/Im are therefore clamped to `±SURF_CLAMP`
#     for display, chosen as a few × the smooth-sector amplitude (the
#     negative-axis ridge `u(-20) ≈ +4.9`; the sector excursions stay
#     within ~±10).  `SURF_CLAMP = 15.0` keeps the smooth sector fully
#     resolved while the clamped wedge still shows the characteristic KKG
#     "jagged" pole structure (clamped peaks read as a sawtooth plateau).
#     This is a *display* clamp only — the kernel's matrices are untouched.
#
#   * **Colormap.**  A diverging `:RdBu` for the signed `Re`/`Im` fields
#     (zero at white) — the FW-family convention for signed quantities;
#     the unsigned `|u|` pole-field figures use viridis, but Re/Im are
#     signed and a diverging map reads the sign at a glance.
#
#   * **3D view.**  `Axis3` azimuth/elevation borrowed from
#     `figutil.jl`'s `pole_field_figure` (the MATLAB-`surf`-like FW view).
#
#   * **Pole overlay.**  The wedge panel's primary content is the
#     validated pole *field* (ADR-0025 Amendment 2): the extracted
#     `poles` are scattered as black dots over the honest partial `|u|`
#     underlay — the FW 2011 Fig 4.7 idiom.  Where the B1 gate found no
#     honest datum the underlay is `NaN`/grey: an explicit gap, not an
#     extrapolation.
#
# ## v1 corners (CLAUDE.md Rule 9) — inherited from the F3 kernel
#
#   1. **Inner-arc asymptotic seed.**  Methods 2/3's inner arc
#      (`|x| ≈ 2`) is the KKG `n_terms = 2` asymptotic series, accurate
#      only to `O(10⁻²)` at `|x| ≈ 2`.  The harmonic solve damps it and
#      the exact ray-fan voter anchors the vote.  *Forcing condition for
#      v2:* a figure pin inside `|x| ≲ 3` to better than `1e-2` needs a
#      small-`|x|` Taylor/BVP patch.  (The sector's corner — Phase C.)
#
#   2. **RETIRED — wedge Stage-2 extrapolation (C2).**  The shipped
#      figure evaluated the wedge Padé far outside its disc
#      (`extrapolate = true`).  ADR-0025 Amendment 1's B1 true-radius
#      gate retires it: the Stage-2 fill now runs `extrapolate = false`,
#      the honest gaps render `NaN`/grey, and the wedge panel is the
#      validated pole field + partial underlay (Amendment 2), not a
#      dishonestly-filled surface.
#
#   3. **RETIRED — hand-tuned wedge step `h = 0.1` (C3).**  ADR-0025
#      Amendment 4's B2 adaptive `:max_q_root` walk (the package
#      default) retires the hand-tune: the step dodges the nearest
#      shared-`Q` pole and `h` adapts to the local density, threading
#      the B3 extended fan to `|x| ≈ 18-20`.
#
#   4. **Stokes-strip masking.**  `±3°` of arc on each `±36°` Stokes line
#      is `NaN`-masked — both region solvers degrade at the boundary.
#      *Forcing condition for v2:* a uniform connection formula across the
#      Stokes line would close the gap (Phase E).
#
# The deferred fully-filled honest wedge *surface* (a 2D-fill
# architecture) is recorded as a deferred bead under ADR-0025
# Amendment 2 §Deferred.
#
# References:
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015: eq. (1.1) the ODE, eq. (1.3) the negative-axis
#     asymptotics, eq. (1.4) the pole-free sector, §7 Figs 7.4/7.5 the
#     Re/Im surfaces this script reproduces.
#   * `docs/adr/0025-headline-figure-re-resolution.md` — the
#     re-resolution ADR; Amendment 1 (B1 gate), Amendment 2 (wedge
#     rescope to pole field + honest partial underlay), Amendment 4
#     (B2 adaptive walk).
#   * `docs/adr/0024-laplace-harmonic-extension.md` — the triple-method
#     majority-vote decision; the conformal `w = log x` rectangle map.
#   * `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — §4 the sector
#     geometry, the sign-corrected tritronquée IC.
#   * `figures/_kkg_pi2_surface_helpers.jl` — the F3 compute kernel.
#   * `figures/figutil.jl` — the FW-family `Axis3` rendering convention.

using CairoMakie
using Printf

# The F3 kernel `include`s the Gridap-backed Laplace voter, whose
# `_kkg_pi2_gridap_helper.jl` uses an unqualified `Point`.  Pulled into
# the same `Main` as `using CairoMakie`, `Point` becomes ambiguous
# (CairoMakie/GeometryBasics vs Gridap.Fields both export it).  Isolating
# the kernel in its own module gives it a clean namespace — `using Gridap`
# there resolves `Point` unambiguously — and keeps CairoMakie's `Point`
# out of the kernel's scope.  The kernel file itself is untouched.
module KKGSurfaceKernel
include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))
end
using .KKGSurfaceKernel: kkg_pi2_surface, SURF_GRID_N, SURF_XY_LIM

const OUTPNG = joinpath(@__DIR__, "output",
                        "kkg_pi2_tritronquee_surface.png")

# Display clamp — a few × the smooth-sector amplitude (see header).  The
# kernel's matrices are NOT modified; this is applied only to display
# copies.  ±15 fully resolves the ~270° smooth sector while the clamped
# wedge plateau still reads as the KKG jagged pole structure.
const SURF_CLAMP = 15.0

# ----------------------------------------------------------------------
# Compute — the exact F3 pipeline the acceptance test reproduces.
# ----------------------------------------------------------------------
@printf("KKG 2015 Figs 7.4/7.5 — P_I^(2) tritronquée V₀(x,0) surface\n")
@printf("  triple-method majority-vote sector + path-network wedge (ADR-0024)\n")
@printf("  grid: %d×%d over [-%.0f,%.0f]²\n",
        SURF_GRID_N, SURF_GRID_N, SURF_XY_LIM, SURF_XY_LIM)

t0  = time()
res = kkg_pi2_surface()
@printf("  kernel: %.1f s — %s\n", time() - t0, res.message)

xs, ys = res.xs, res.ys

# Quantitative sanity readout — the negative-real-axis ridge must be
# positive (the sign-bug regression guard; KKG eq. (1.3)).
let iy0 = argmin(abs.(ys)),               # the y ≈ 0 row
    ix  = argmin(abs.(xs .- (-15.0)))     # x ≈ -15 column
    ridge = res.Re_u[ix, iy0]
    @printf("  ridge check: Re V₀(-15, 0) = %+.4f  (KKG (6·15)^(1/3) = +%.4f)\n",
            ridge, cbrt(90.0))
end

nfin = count(isfinite, res.Re_u)
@printf("  filled grid cells: %d / %d ; wedge poles: %d\n",
        nfin, length(res.Re_u), length(res.poles))

# ----------------------------------------------------------------------
# Display matrices — a clamped copy for the colour/z range.  NaN cells
# (outside the disc, masked strips, no datum) stay NaN: Makie renders
# them in `nan_color`.
# ----------------------------------------------------------------------
display_clamp(M) = [isfinite(v) ? clamp(v, -SURF_CLAMP, SURF_CLAMP) : NaN
                    for v in M]
Re_disp = display_clamp(res.Re_u)
Im_disp = display_clamp(res.Im_u)

const CRANGE   = (-SURF_CLAMP, SURF_CLAMP)
const CMAP     = :RdBu                       # diverging — signed Re/Im
const NANGREY  = RGBAf(0.82, 0.82, 0.82, 1)  # neutral for NaN cells

# Pole overlay coordinates (wedge pole field, KKG 7.4/7.5 jagged region).
pole_re = Float64.(real.(res.poles))
pole_im = Float64.(imag.(res.poles))

# ----------------------------------------------------------------------
# Render — a 2×2 panel grid.
#   row 1 (KKG Fig 7.4): Re V₀ — heatmap | 3D surface
#   row 2 (KKG Fig 7.5): Im V₀ — heatmap | 3D surface
# ----------------------------------------------------------------------
fig = Figure(size = (1320, 1180))

Label(fig[0, 1:2],
      "P_I⁽²⁾ tritronquée V₀(x, 0) over the complex x-plane [-20,20]²  " *
      "—  KKG 2015 Figs 7.4 (Re) / 7.5 (Im)";
      fontsize = 16, padding = (0, 0, 6, 10))

"""
    heatmap_panel(gp, M, title) -> Axis

A 2D heatmap of the clamped field `M` in an `Axis` at grid slot `gp`,
square data aspect, NaN cells in neutral grey, the wedge pole field
overlaid as small black dots.
"""
function heatmap_panel(gp, M, title)
    ax = Axis(gp; title = title, xlabel = "Re x", ylabel = "Im x",
              titlesize = 13, aspect = DataAspect(),
              limits = (-SURF_XY_LIM, SURF_XY_LIM,
                        -SURF_XY_LIM, SURF_XY_LIM))
    hm = heatmap!(ax, xs, ys, M; colormap = CMAP, colorrange = CRANGE,
                  nan_color = NANGREY)
    # the wedge pole field — KKG's jagged region, made explicit
    scatter!(ax, pole_re, pole_im; color = :black, markersize = 4.0)
    return ax, hm
end

"""
    surface_panel(gp, M, title) -> Axis3

A 3D `Axis3` surface of the clamped field `M`, the MATLAB-`surf`-like
FW view (azimuth/elevation from `figutil.jl`), NaN cells transparent.
"""
function surface_panel(gp, M, title)
    ax = Axis3(gp; title = title, xlabel = "Re x", ylabel = "Im x",
               zlabel = "value", titlesize = 13,
               azimuth = 1.30π, elevation = 0.16π,
               aspect = (1.0, 1.0, 0.55))
    surface!(ax, xs, ys, M; colormap = CMAP, colorrange = CRANGE,
             nan_color = RGBAf(0, 0, 0, 0))
    return ax
end

# ---- row 1 — Re V₀ (KKG Fig 7.4) ------------------------------------
axR_h, hmR = heatmap_panel(fig[1, 1], Re_disp,
    "A. Re V₀ — heatmap  [KKG Fig 7.4]")
axR_s      = surface_panel(fig[1, 2], Re_disp,
    "B. Re V₀ — surface  [KKG Fig 7.4]")

# ---- row 2 — Im V₀ (KKG Fig 7.5) ------------------------------------
axI_h, hmI = heatmap_panel(fig[2, 1], Im_disp,
    "C. Im V₀ — heatmap  [KKG Fig 7.5]")
axI_s      = surface_panel(fig[2, 2], Im_disp,
    "D. Im V₀ — surface  [KKG Fig 7.5]")

# ---- shared colourbar ------------------------------------------------
Colorbar(fig[1:2, 3], hmR;
         label = @sprintf("V₀ component  (display clamp ±%.0f)", SURF_CLAMP))

# A standing note: the clamp + the wedge rescope, on the figure itself.
const FIGNOTE = string(
    "Smooth ~270° sector: triple-method majority vote (ray-fan BVP / ",
    "2D-Chebyshev + Gridap FEM Laplace, ADR-0024).  Wedge |arg x|<36° ",
    "(ADR-0025 Amend. 2): validated pole field — ",
    @sprintf("%d extracted poles (black dots) — over an honest partial ",
             length(res.poles)),
    "|u| underlay, B1-gated (no Padé out of disc); display clamped ",
    @sprintf("to ±%.0f.  Grey: NaN — outside |x|≤20 disc, ±3° Stokes ",
             SURF_CLAMP),
    "strips, and the wedge gaps the honest B1 gate leaves unfilled.")
Label(fig[3, 1:3], FIGNOTE; fontsize = 9, padding = (0, 0, 8, 0))

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)

# Light render-sanity check (F4 is the render; the F3 kernel is already
# pinned by `figures/test_kkg_pi2_surface.jl`).  Assert the PNG was
# actually produced and is non-trivially sized — a save that silently
# emitted nothing must fail loud (CLAUDE.md Rule 1).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("  render-sanity: PNG present, %d bytes\n", bytes)
end
