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
#     validated pole *field* + a **dual-fill provenance** `|u|` surface
#     (ADR-0025 Amendment 14 / ADR-0026 Amendment 10).  The headline
#     figure's wedge is rendered FW-2011-style filled, with the
#     B1-honest certified core visually distinguished from the FW-style
#     extrapolated overlay:
#
#       - the `:max_q_root` resilient adaptive walk threads a dense
#         disc-spaced Cartesian lattice through the whole wedge to
#         `|x| = 20`, `extract_poles_shared_q` (`min_support ≥ 2`,
#         VC-6) extracts the validated pole *locations* (the FW 2011
#         Fig 4.7 idiom);
#       - the Stage-2 fine-grid fill runs **twice on the same walk** —
#         once with `extrapolate = false` (B1 true-radius gate, ~29 %
#         of the wedge is certified at order 48), once with
#         `extrapolate = true` (no validity gate, every cell finite,
#         identical to FW 2011's published rendering);
#       - this figure script composites the two: **certified cells at
#         full opacity, FW-style extrapolated cells at reduced alpha**.
#         The reader sees the full surface FW renders, and reads
#         provenance off the alpha channel.
#
#     ADR-0025's strict "no Padé out of disc" rule remains the package
#     default and applies to every non-headline figure; this dual-fill
#     is a per-figure override, justified by FW precedent (FW's own
#     pole-field figures evaluate each Padé over its full step with no
#     validity gate — FW-md:395-397) and honest *only when* the figure
#     marks provenance, which this one does.
#
#   A thin `±1°` strip straddling each `±36°` Stokes line is `NaN`-masked
#   — the band the sector fan does not reach and the wedge underlay does
#   not honestly cover (ADR-0025 Amendment 12's Phase-E1 audition
#   narrowed it `±3° → ±1°`).  Cells outside the `|x| ≤ 20` disc are
#   `NaN`.  `NaN` cells render in a neutral grey.
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
#   * **Pole overlay.**  The validated pole *field* (ADR-0025
#     Amendment 14 / ADR-0026 Amendment 10) is scattered as black dots
#     over the dual-fill surface — the FW 2011 Fig 4.7 idiom.  The
#     wedge underlay is filled to the wedge boundary (FW-2011-style);
#     full-opacity cells are B1-certified and reduced-alpha cells are
#     FW-style extrapolated past the verified disc.  The pole field
#     itself is unchanged from the previous certified-only rendering —
#     it is validated via VC-4 (dominant-balance per-pole certificate)
#     + VC-5 (conjugate-pairing) and is independent of the underlay
#     fill mode.
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
#   2. **RETIRED → DUAL-FILL PROVENANCE — wedge Stage-2 extrapolation
#      (C2).**  The shipped V8b figure evaluated the wedge Padé far
#      outside its disc (`extrapolate = true`) with no provenance marker
#      — dishonestly filled.  ADR-0025 Amendment 1's B1 true-radius gate
#      retired the blanket dishonest fill; ADR-0025 Amendment 14 (this
#      figure's `tf9` closeout) then re-introduces `extrapolate = true`
#      as a **per-figure, provenance-marked overlay**: the Stage-2 fill
#      runs twice on the same walk (certified `extrapolate = false`
#      core + FW-style filled overlay), and the figure script paints
#      certified cells at full opacity and extrapolated cells at
#      reduced alpha (α ≈ 0.50).  Justified by FW precedent — FW 2011's
#      own published pole-field figures FILL precisely because their
#      Stage-2 evaluates each Padé over its full step with no validity
#      gate (FW-md:395-397).  The strict B1 gate remains the package
#      default for every non-headline figure.
#
#   3. **RETIRED — hand-tuned wedge step `h = 0.1` (C3).**  ADR-0025
#      Amendment 4's B2 adaptive `:max_q_root` walk (the package
#      default) retires the hand-tune: the step dodges the nearest
#      shared-`Q` pole and `h` adapts to the local density, threading
#      the B3 extended fan to `|x| ≈ 18-20`.
#
#   4. **NARROWED (C4) — Stokes-strip masking.**  ADR-0025 Amendment 12's
#      Phase-E1 audition (bead `padetaylor-0ln.37.18`) narrowed the
#      masked strip `±3° → ±1°`: the audition measured the sector ODE
#      solve healthy right up to the 36° Stokes line and the
#      triple-method vote honest down to a `1°` fan margin.  The grey
#      strip shrinks 3× (`6° → 2°` total).  *Forcing condition for the
#      residual `±1°` mask:* a uniform connection formula across the
#      Stokes line would close the last `2°` (an Amendment-12 deferred
#      bead).
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
using .KKGSurfaceKernel: kkg_pi2_surface, SURF_GRID_N, SURF_XY_LIM,
                         SURF_STITCH_MASK_DEG, SURF_PN_ORDER

const OUTPNG = joinpath(@__DIR__, "output",
                        "kkg_pi2_tritronquee_surface.png")

# Display clamp — a few × the smooth-sector amplitude (see header).  The
# kernel's matrices are NOT modified; this is applied only to display
# copies.  ±15 fully resolves the ~270° smooth sector while the clamped
# wedge plateau still reads as the KKG jagged pole structure.
const SURF_CLAMP = 15.0

# Reduced alpha for the FW-style extrapolated wedge overlay (ADR-0025
# Amendment 14 / ADR-0026 Amendment 10).  Extrapolated cells — wedge
# cells the B1 honest gate leaves uncovered, filled here FW-2011-style
# (every Padé evaluated over its full step with no validity gate, exactly
# FW's published rendering) — are painted at this alpha to mark their
# provenance.  Certified cells (full opacity, alpha = 1.0) read as
# vivid; extrapolated cells read as washed-out but still legible.  `0.50`
# is the chosen value: at 0.45 the extrapolated overlay is too pale to
# read pole structure on the colour-clamped wedge, at 0.55 the visual
# distinction from full-opacity cells is too subtle on the 2D heatmap
# panels.  Picked alpha (rather than hatch or desaturation) because it
# composes natively in Makie's `heatmap!`/`surface!` and the alpha drop
# reads cleanly on both the 2D heatmap and the 3D surface panels in one
# layout pass.
const EXTRAP_ALPHA = 0.50

# ----------------------------------------------------------------------
# Compute — the exact F3 pipeline the acceptance test reproduces.
# ----------------------------------------------------------------------
@printf("KKG 2015 Figs 7.4/7.5 — P_I^(2) tritronquée V₀(x,0) surface\n"); flush(stdout)
@printf("  triple-method majority-vote sector + path-network wedge (ADR-0024)\n"); flush(stdout)
@printf("  grid: %d×%d over [-%.0f,%.0f]²\n",
        SURF_GRID_N, SURF_GRID_N, SURF_XY_LIM, SURF_XY_LIM); flush(stdout)

t0  = time()
@printf("  kernel: starting kkg_pi2_surface (order=%d, expect tens of minutes at order 48)...\n",
        SURF_PN_ORDER); flush(stdout)
res = kkg_pi2_surface()
@printf("  kernel: %.1f s — %s\n", time() - t0, res.message); flush(stdout)

xs, ys = res.xs, res.ys

# Quantitative sanity readout — the negative-real-axis ridge must be
# positive (the sign-bug regression guard; KKG eq. (1.3)).
let iy0 = argmin(abs.(ys)),               # the y ≈ 0 row
    ix  = argmin(abs.(xs .- (-15.0)))     # x ≈ -15 column
    ridge = res.Re_u[ix, iy0]
    @printf("  ridge check: Re V₀(-15, 0) = %+.4f  (KKG (6·15)^(1/3) = +%.4f)\n",
            ridge, cbrt(90.0)); flush(stdout)
end

nfin     = count(isfinite, res.Re_u)
nfin_ext = count(isfinite, res.Re_u_extrap)
@printf("  certified filled grid cells: %d / %d (%.1f%%)\n",
        nfin, length(res.Re_u), 100 * nfin / length(res.Re_u)); flush(stdout)
@printf("  FW-style filled grid cells:  %d / %d (%.1f%%)\n",
        nfin_ext, length(res.Re_u_extrap),
        100 * nfin_ext / length(res.Re_u_extrap)); flush(stdout)
@printf("  wedge poles: %d\n", length(res.poles)); flush(stdout)

# ----------------------------------------------------------------------
# Display matrices — clamped copies for the colour/z range.  NaN cells
# (outside the disc, masked strips, no datum) stay NaN: Makie renders
# them in `nan_color`.
#
# The wedge is rendered as a TWO-LAYER composite (ADR-0025 Amendment 14
# / ADR-0026 Amendment 10):
#
#   * `Re_certified` / `Im_certified` — the B1-honest cells (the
#     sector vote + the certified wedge core); painted at full opacity.
#     Wedge cells the B1 gate did not certify are `NaN` here.
#
#   * `Re_extrap_only` / `Im_extrap_only` — the FW-2011-style extrapolated
#     wedge cells *only* (`!wedge_covered`); painted at `EXTRAP_ALPHA`
#     reduced opacity.  Certified cells are `NaN` here so the alpha
#     layer composites cleanly underneath the full-opacity certified
#     layer rather than dimming it.
#
# At render time the two layers are drawn into the same axis (extrap
# layer first, certified layer on top).  The reader sees a fully filled
# wedge — FW-faithful — and reads provenance off the alpha channel.
# ----------------------------------------------------------------------
display_clamp(M) = [isfinite(v) ? clamp(v, -SURF_CLAMP, SURF_CLAMP) : NaN
                    for v in M]

# Build the certified layer — the kernel's `Re_u`/`Im_u` are already the
# B1-certified surface (sector vote + certified wedge core; uncovered
# wedge cells are NaN).  Just clamp for display.
Re_certified = display_clamp(res.Re_u)
Im_certified = display_clamp(res.Im_u)

# Build the extrap-only layer — start from the kernel's `Re_u_extrap`
# (the FW-style filled overlay) and NULL OUT cells that the certified
# layer also fills.  Result: the extrap layer carries values ONLY where
# the wedge gate uncovered them, which is exactly the set of cells the
# alpha drop should mark.  This avoids the alpha layer compositing on
# top of the certified layer (which would dim certified cells too).
function extrap_only(M_extrap, M_cert)
    n_x, n_y = size(M_extrap)
    out = fill(NaN, n_x, n_y)
    @inbounds for j in 1:n_y, i in 1:n_x
        # If the certified layer covers this cell, the extrap layer must
        # stay NaN here — the certified layer paints in full opacity.
        if !isnan(M_cert[i, j])
            continue
        end
        v = M_extrap[i, j]
        if isfinite(v)
            out[i, j] = clamp(v, -SURF_CLAMP, SURF_CLAMP)
        end
    end
    return out
end
Re_extrap_only = extrap_only(res.Re_u_extrap, res.Re_u)
Im_extrap_only = extrap_only(res.Im_u_extrap, res.Im_u)

const CRANGE   = (-SURF_CLAMP, SURF_CLAMP)
const CMAP     = :RdBu                       # diverging — signed Re/Im
const NANGREY  = RGBAf(0.82, 0.82, 0.82, 1)  # neutral for NaN cells
# An alpha-aware NaN colour for the extrap layer — the cells the extrap
# layer leaves NaN are either (a) cells the certified layer paints (so
# the certified layer covers them) or (b) cells outside the figure /
# masked strips (the certified layer leaves NaN too, and the underlying
# axis background will read).  Fully transparent so neither case dims.
const NANCLEAR = RGBAf(0, 0, 0, 0)

# Pole overlay coordinates (wedge pole field, KKG 7.4/7.5 jagged region).
pole_re = Float64.(real.(res.poles))
pole_im = Float64.(imag.(res.poles))

# Provenance accounting — certified vs extrapolated cell counts, surfaced
# on the figure caption and the stdout report.  `Re_extrap_only` is a
# `Matrix{Float64}` with NaN-elsewhere semantics; `!isnan` counts the
# finite reduced-alpha cells.
n_certified    = count(isfinite, res.Re_u)
n_extrapolated = count(!isnan, Re_extrap_only)

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
    heatmap_panel(gp, M_cert, M_extrap, title) -> Axis

A 2D heatmap of the clamped field at grid slot `gp`, rendered as a
two-layer dual-provenance composite (ADR-0025 Amendment 14 / ADR-0026
Amendment 10):

  * the **extrapolated layer** (`M_extrap`, FW-2011-style filled wedge
    cells the certified layer leaves NaN) is drawn first at
    `EXTRAP_ALPHA = $(EXTRAP_ALPHA)` reduced opacity;

  * the **certified layer** (`M_cert`, sector vote + B1-honest wedge
    cells) is drawn on top at full opacity.

The wedge pole field is then overlaid as small black dots.  Square data
aspect; cells outside `|x| ≤ 20` and inside the masked Stokes strips
render in `NANGREY` (the certified layer's neutral grey shows through
where both layers are NaN).
"""
function heatmap_panel(gp, M_cert, M_extrap, title)
    ax = Axis(gp; title = title, xlabel = "Re x", ylabel = "Im x",
              titlesize = 13, aspect = DataAspect(),
              limits = (-SURF_XY_LIM, SURF_XY_LIM,
                        -SURF_XY_LIM, SURF_XY_LIM))
    # extrapolated wedge cells — drawn at reduced alpha, with a fully-
    # transparent nan_color so the certified layer below paints through.
    heatmap!(ax, xs, ys, M_extrap; colormap = CMAP, colorrange = CRANGE,
             nan_color = NANCLEAR, alpha = EXTRAP_ALPHA)
    # certified cells (and the masked / outside-disc NaN background) —
    # full opacity, neutral grey for NaN.
    hm = heatmap!(ax, xs, ys, M_cert; colormap = CMAP, colorrange = CRANGE,
                  nan_color = NANGREY)
    # the wedge pole field — KKG's jagged region, made explicit
    scatter!(ax, pole_re, pole_im; color = :black, markersize = 4.0)
    return ax, hm
end

"""
    surface_panel(gp, M_cert, M_extrap, title) -> Axis3

A 3D `Axis3` surface of the clamped field, rendered as a two-layer
dual-provenance composite (ADR-0025 Amendment 14 / ADR-0026 Amendment 10):
the extrapolated wedge layer at `EXTRAP_ALPHA` reduced opacity, the
certified layer (sector + B1-honest wedge core) on top at full opacity.
NaN cells transparent.  The MATLAB-`surf`-like FW view (azimuth /
elevation from `figutil.jl`).
"""
function surface_panel(gp, M_cert, M_extrap, title)
    ax = Axis3(gp; title = title, xlabel = "Re x", ylabel = "Im x",
               zlabel = "value", titlesize = 13,
               azimuth = 1.30π, elevation = 0.16π,
               aspect = (1.0, 1.0, 0.55))
    # extrapolated wedge cells first, at reduced alpha
    surface!(ax, xs, ys, M_extrap; colormap = CMAP, colorrange = CRANGE,
             nan_color = NANCLEAR, alpha = EXTRAP_ALPHA)
    # certified cells (and everything else) on top, full opacity
    surface!(ax, xs, ys, M_cert; colormap = CMAP, colorrange = CRANGE,
             nan_color = NANCLEAR)
    return ax
end

# ---- row 1 — Re V₀ (KKG Fig 7.4) ------------------------------------
axR_h, hmR = heatmap_panel(fig[1, 1], Re_certified, Re_extrap_only,
    "A. Re V₀ — heatmap  [KKG Fig 7.4]")
axR_s      = surface_panel(fig[1, 2], Re_certified, Re_extrap_only,
    "B. Re V₀ — surface  [KKG Fig 7.4]")

# ---- row 2 — Im V₀ (KKG Fig 7.5) ------------------------------------
axI_h, hmI = heatmap_panel(fig[2, 1], Im_certified, Im_extrap_only,
    "C. Im V₀ — heatmap  [KKG Fig 7.5]")
axI_s      = surface_panel(fig[2, 2], Im_certified, Im_extrap_only,
    "D. Im V₀ — surface  [KKG Fig 7.5]")

# ---- shared colourbar ------------------------------------------------
Colorbar(fig[1:2, 3], hmR;
         label = @sprintf("V₀ component  (display clamp ±%.0f)", SURF_CLAMP))

# A standing note: the dual-provenance wedge rendering, on the figure
# itself.  Honesty marker — the reader must be told which cells are
# certified and which are extrapolated (ADR-0025 Amendment 14 /
# ADR-0026 Amendment 10).
const FIGNOTE = string(
    "Smooth ~270° sector: triple-method majority vote (ray-fan BVP / ",
    "2D-Chebyshev + Gridap FEM Laplace, ADR-0024).  Wedge |arg x|<36° ",
    "(ADR-0025 Amend. 14 / ADR-0026 Amend. 10): dual-fill provenance — ",
    @sprintf("FULL OPACITY cells are B1-honest (no Padé evaluated past " *
             "its verified disc, %.1f%% of grid); REDUCED ALPHA " *
             "(α = %.2f) cells are evaluated past their verified disc, " *
             "FW-2011-style (%.1f%% additional).  ",
             100.0 * n_certified / length(res.Re_u),
             EXTRAP_ALPHA,
             100.0 * n_extrapolated / length(res.Re_u)),
    @sprintf("%d validated wedge poles (black dots); display clamped " *
             "to ±%.0f.  Grey: NaN — outside |x|≤20 disc, ±%.0f° " *
             "Stokes strips.",
             length(res.poles), SURF_CLAMP, SURF_STITCH_MASK_DEG))
Label(fig[3, 1:3], FIGNOTE; fontsize = 9, padding = (0, 0, 8, 0))

@printf("  building Makie figure + saving PNG...\n"); flush(stdout)
mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG); flush(stdout)

# Light render-sanity check (F4 is the render; the F3 kernel is already
# pinned by `figures/test_kkg_pi2_surface.jl`).  Assert the PNG was
# actually produced and is non-trivially sized — a save that silently
# emitted nothing must fail loud (CLAUDE.md Rule 1).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("  render-sanity: PNG present, %d bytes\n", bytes); flush(stdout)
end
