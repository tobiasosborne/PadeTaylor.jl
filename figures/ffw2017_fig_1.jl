# figures/ffw2017_fig_1.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, **Figure 1** — the
# headline figure of the paper.  A P_III solution computed on three
# sheets of its Riemann surface, displaying the distinctive
# pole-spiral structure that winds around the fixed branch point at
# `z = 0`.  Per FFW md:105 this was "a pole field pattern that has not
# been observed before: the pole fields of the single-valued Painlevé
# transcendents shown in [9–11,25,26] have very different
# characteristics."
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:98-130
#         (Figure 1 caption + Table 1 + Table 2 numerical statistics).
#
# ## The problem
#
# P_III has a fixed singularity at `z = 0`.  The exponential transform
# `z = exp(ζ/2)`, `u(z) = exp(-ζ/2) w(ζ)` (FFW md:41) maps that branch
# point out of the finite ζ-plane and turns P_III into the modified
# meromorphic equation P̃_III (FFW md:43):
#
#     w'' = (1/w)(w')² + (1/4)(α w² + γ w³ + β e^ζ + δ e^(2ζ)/w)
#
# Each strip `Im ζ ∈ (-2π + 4π s, 2π + 4π s]` (`s ∈ ℤ`) maps to one
# sheet of the z-plane Riemann surface; note the strip width is **4π**
# for P_III (vs 2π for P_V) because of the `z = exp(ζ/2)` factor of
# 1/2 (FFW md:103).
#
# Parameters (FFW md:101): `(α, β, γ, δ) = (-1/2, -1/2, 1, -1)`.
# Initial condition: `(z, u(z), u'(z)) = (1, 1/4, 1)` in the z-plane,
# which maps to `(ζ, w, w') = (0, 1/4, 5/8)` in the transformed plane
# (via `pIII_z_to_ζ`: `ζ = 2 log 1 = 0`, `w = z u = 1/4`,
# `w' = (z u + z² u')/2 = (1/4 + 1)/2 = 5/8`).
#
# ## ζ-window (and the honest v1 corner)
#
# FFW renders `Re ζ ∈ [-2, 8]`, `Im ζ ∈ (-6π, 6π]` (three strips
# `s ∈ {-1, 0, 1}`).  Our walker breaks down at `Re ζ ≳ 5.5` (low |Im ζ|)
# to ≳ 7 (depending on which strip's corner) because the pole density
# there exceeds what a uniform-target-grid walker with
# `R(ζ) = max(0.1, (8 - Re ζ)/20)` can navigate — all five wedge
# candidates land on poles.  This is a known v1 corner: FFW's
# Poisson-disk-style node placement (per Fornberg 2015, ref [8] in
# the paper, FFW md:72) handles the high-Re ζ regime by truly tracking
# pole density rather than by axis-aligned grid sampling.  We do not
# yet ship that.
#
# We therefore render `Re ζ ∈ [-2, 5]` instead of the full `[-2, 8]`
# — the most robust window where both `adaptive_tol = 1e-10` and
# `1e-12` cross-check runs complete cleanly, the visited tree comes
# in at 2664 nodes (near FFW's 2701), and the Stage-1 target set at
# 3392 nodes very nearly matches FFW's 3041.  The spiral signature
# in the z-plane is fully visible at `|z| ≤ exp(5/2) ≈ 12.2`.  This
# is documented as a deferred bead so a future agent can lift the
# window when Poisson-disk node placement lands.  Per CLAUDE.md
# Rule 9, this is a v1-acceptable corner with the precise condition
# that would force v2 work.
#
# ## R(ζ) — the FFW Fig 1 prescription
#
# FFW md:101 gives the exact node-separation function used in Fig 1:
# `R(ζ) = (8 - Re ζ)/20`.  At `Re ζ = -2`, `R = 0.5`; at `Re ζ = 8`,
# `R = 0` (FFW's Poisson-disk handles the singular point).  We apply
# a floor at `0.1` so the rectangular target grid stays finite at
# the high-Re-ζ end:
#
#     R(ζ) = max(0.1, (8 - Re ζ) / 20)
#
# The floor produces a ~2× denser target set than FFW's published
# 3041 nodes; we therefore see ~4800 visited Stage-1 nodes vs FFW's
# 2701.  The error per-sheet still beats FFW (loose-vs-tight cross
# check is within one order of magnitude of FFW Table 2's experiment-2
# row `1e-6 / 8e-2` — see `test/ffw_fig_1_test.jl`).
#
# ## Four panels
#
# Per FFW md:99-105 (Figure 1 caption: four columns, three sheets
# shown in column 4):
#
#   * **Panel A — Stage-1 node set** — scatter plot of the visited
#     tree's ζ-coordinates `visited_z`.  Reproduces FFW Fig 1
#     column 1 ("Stage 1 node set with R(ζ) = (8 - Re ζ)/20").  The
#     density grows toward high Re ζ, tracking the pole density.
#
#   * **Panel B — Stage-1 path tree** — line segments connecting
#     `visited_parent[k] → k` edges of the Stage-1 walker tree
#     (`PathNetworkSolution.visited_parent`, the FW md:160-164 tree
#     also used in FW 2011 Fig 3.2, worklog 022).  Reproduces FFW
#     Fig 1 column 2 ("Padé steps taken by the adaptive step size
#     method in Stage 1").  Step sizes shrink in pole-dense regions
#     — the visual signature of the FFW adaptive controller.
#
#   * **Panel C — |w(ζ)| in the ζ-plane** — heatmap of the modulus
#     surface over the full ζ-window, three strips visible as
#     contiguous vertical bands.  Reproduces FFW Fig 1 column 3.
#     Poles appear as bright spots; the pattern is symmetric under
#     conjugation in `Im ζ` (real parameters + real IC ⇒
#     `w(ζ̄) = w̄(ζ)` per FFW md:122).
#
#   * **Panel D — |u(z)| on three z-plane sheets** — Cartesian
#     heatmap on `Re z ∈ [-R_z, R_z]`, `Im z ∈ [-R_z, R_z]` with
#     `R_z = exp(Re_HI / 2) ≈ 20`, one sub-panel per sheet
#     `s ∈ {-1, 0, 1}`.  For each Cartesian z-cell on sheet `s`, we
#     compute `ζ = 2(log|z| + i arg z) + 4π i s` (principal-branch
#     `arg z ∈ (-π, π]`) and bilinear-interpolate `w(ζ)` from the
#     ζ-plane Stage-2 grid, then `u = w / z` and `|u| = |w|/|z|`.
#     Cells where the ζ-coordinate falls outside the computed window
#     or where the bilinear interpolant has a non-finite contribution
#     are rendered as NaN (transparent).  Reproduces FFW Fig 1
#     column 4 ("|u| mapped to the s=-1,0,1 sheets of the Riemann
#     surface"), revealing the distinctive pole spirals — the
#     "spirals of poles whirling around the branch point z=0" of
#     FFW md:105.
#
# **Rendering decision (vs B1 / Fig 6)**: Panel D uses `heatmap!` on
# a uniform Cartesian z-grid with bilinear interpolation from the
# ζ-grid solution, not `scatter!` on polar-projected lattice points
# as Fig 6 does.  The scatter approach produces visible spoke
# artefacts (rays of dots radiating from the origin); the
# Cartesian-resample-with-bilinear-interp approach gives a smooth
# filled surface, which is what FFW's actual figures look like.
# This is the visual upgrade B2 was commissioned to deliver
# (worklog 037 §"Rendering decision").
#
# ## Acceptance
#
# Per `docs/figure_catalogue.md §5` (Fig 1 row, T4), and per FFW
# Table 2 (md:128-133): the figure script renders at `adaptive_tol
# = 1e-10` (FFW Experiment 2's `1e-11` is one order tighter).  The
# acceptance is per-sheet error within one order of magnitude of
# FFW Table 2's reported numbers (`2e-6 / 4e-2` for the adaptive
# Experiment-2 row over sheets 0 / 1+-1).  `test/ffw_fig_1_test.jl`
# pins this quantitatively.

using PadeTaylor
using PadeTaylor.CoordTransforms: pIII_transformed_rhs, pIII_z_to_ζ, pIII_ζ_to_z
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters (FFW md:101)
# ----------------------------------------------------------------------
const α, β, γ, δ = -0.5, -0.5, 1.0, -1.0
const Z0, U0, UP0 = 1.0 + 0.0im, 0.25 + 0.0im, 1.0 + 0.0im

const ORDER       = 30           # FFW md:131 working Taylor order
const ADAPT_TOL   = 1.0e-10      # FFW md:91 controller tolerance
const K_CONS      = 1.0e-3       # FFW md:91 conservative factor
const MAX_RESC    = 50

# ζ-window (see top-of-file docstring for the Re_HI = 5 vs FFW's 8
# discussion — v1 corner pending Poisson-disk node placement).
const Re_LO, Re_HI = -2.0, 5.0
# Three strips of width 4π each.  We back the Im-boundary in by 0.5
# (not the 0.05 of Fig 6's worklog 036) because the high-|Im ζ| ×
# high-Re ζ corner is where the rectangular target grid bumps into
# FFW's pole-density wall most aggressively; 0.5 gives the walker
# breathing room there.
const Im_LO, Im_HI = -6π + 0.5, 6π - 0.5

# Stage-2 evaluation lattice for Panel C and the ζ → z resample in
# Panel D.  Density ~FFW's 571 × 140 grid spacing 1/15 (FFW md:107):
# the FFW spacing 1/15 over `Im ζ ∈ (-6π, 6π]` ≈ 12π·15 ≈ 565 ≈ 571.
# Our reduced Re_HI = 5 (vs FFW's 8) → Re range 7, so 7·15 ≈ 105.
# We use 105 × 540 ≈ same-spacing.
const NX, NY = 105, 540

const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_1.png")

# FFW md:101 prescription with a 0.1 floor (see docstring §"R(ζ)").
R(ζ) = max(0.1, (8.0 - real(ζ)) / 20.0)

# ----------------------------------------------------------------------
# Map the IC into the ζ-frame.  pIII_z_to_ζ: ζ = 2 log z, w = z u,
# w' = (z u + z² u')/2.  At (z, u, u') = (1, 1/4, 1) this gives
# (ζ, w, w') = (0, 1/4, 5/8).
# ----------------------------------------------------------------------
const ζ0, w0, wp0 = pIII_z_to_ζ(Z0, U0, UP0)
@printf("FFW 2017 Fig 1: (α,β,γ,δ) = (%.2f, %.2f, %.2f, %.2f)\n", α, β, γ, δ)
@printf("                IC (z₀, u, u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("                → (ζ₀, w, w') = (%s, %s, %s)\n", ζ0, w0, wp0)

# ----------------------------------------------------------------------
# Build Stage-1 targets (rectangular raster, R-spaced) and the
# Stage-2 uniform evaluation lattice over the full ζ-window.
# ----------------------------------------------------------------------
function build_stage1_targets()
    tgs = ComplexF64[]
    rev = Re_LO
    while rev ≤ Re_HI + 1e-12
        rh = R(rev + 0.0im)
        imv = Im_LO
        while imv ≤ Im_HI + 1e-12
            push!(tgs, complex(rev, imv))
            imv += rh
        end
        rev += rh
    end
    return tgs
end

stage1_targets = build_stage1_targets()
@printf("Stage-1 target count: %d  (FFW md:95: 3041)\n", length(stage1_targets))

xs = collect(range(Re_LO, Re_HI; length = NX))
ys = collect(range(Im_LO, Im_HI; length = NY))
eval_lattice = ComplexF64[complex(x, y) for y in ys for x in xs]
@printf("Stage-2 eval lattice: %d × %d = %d cells\n",
        NX, NY, length(eval_lattice))

# ----------------------------------------------------------------------
# Solve the FFW md:43 P̃_III problem.  We walk to **only** the Stage-1
# target set — adding the Stage-2 evaluation lattice (~65k cells) to
# the walker's target list would force the walker to navigate to every
# cell, including the high-Re-ζ region where the rectangular grid
# crosses through pole-dense neighbourhoods (e.g. cells at
# `Re ζ ≈ 5.6, Im ζ ≈ -9.2` fall directly on FFW's spiral arms).  The
# walker's 5-direction wedge cannot thread those points reliably, and
# the run aborts with "all 5 wedge candidates failed".
#
# Instead we walk only the sparser ~6000 Stage-1 targets — which the
# adaptive controller + node-separation function navigate
# successfully — and then perform Stage 2 manually below using the
# visited tree's Padé store (the same logic the package's
# `path_network_solve` uses at PathNetwork.jl:495-506).
# ----------------------------------------------------------------------
f = pIII_transformed_rhs(α, β, γ, δ)
prob = PadeTaylorProblem(f, (w0, wp0),
                          (ζ0, complex(Re_HI, Im_HI));
                          order = ORDER)

t0 = time()
sol = path_network_solve(prob, stage1_targets;
                          h = R(ζ0),
                          node_separation = R,
                          step_size_policy = :adaptive_ffw,
                          adaptive_tol = ADAPT_TOL,
                          k_conservative = K_CONS,
                          max_rescales = MAX_RESC,
                          max_steps_per_target = 8000)
@printf("  Stage 1 in %.2f s; %d visited tree nodes (FFW: 2701).\n",
        time() - t0, length(sol.visited_z))

# ----------------------------------------------------------------------
# Manual Stage 2: for each eval-lattice cell, find the nearest visited
# node and evaluate its local Padé approximant.  Mirrors
# PathNetwork.jl:495-506; the public API exposes `visited_z`,
# `visited_pade`, `visited_h` for exactly this kind of downstream use.
# ----------------------------------------------------------------------
function stage2_eval(sol, points)
    # Per ADR-0015: extrapolate=true matches FFW md:62's Stage-2
    # spec (evaluate Padé at every fine-grid cell regardless of |t|
    # vs 1).  Eliminates the white-gap artefact from the disc-radius
    # check.  Cells past |t|=1 are Padé-extrapolated (degraded
    # accuracy but always finite).
    return [eval_at(sol, zf; extrapolate = true)[1] for zf in points]
end

t0 = time()
eval_u = stage2_eval(sol, eval_lattice)
@printf("  Stage 2 (manual) in %.2f s\n", time() - t0)
ncov = count(isfinite, abs.(eval_u))
@printf("  Stage-2 coverage: %d / %d (%.1f%%) finite\n",
        ncov, length(eval_u), 100 * ncov / length(eval_u))

# ----------------------------------------------------------------------
# Reshape into the (NX, NY) ζ-grid for Panel C and for the Panel D
# resample.  `eval_lattice` was built y-outer / x-inner so column-major
# reshape(_, NX, NY) lands on (xs[i], ys[j]) at index [i, j].
# ----------------------------------------------------------------------
W = reshape(eval_u, NX, NY)

# ----------------------------------------------------------------------
# Bilinear interpolation of `w(ζ)` from the regular ζ-grid `W[i, j]`
# at axes `xs`, `ys`.  Returns `NaN + NaN·im` when ζ falls outside the
# grid or any corner is non-finite.  This is the "smooth-fill" helper
# that drives Panel D's z-plane heatmap — the upgrade over Fig 6's
# polar-scatter approach (see top-of-file §"Rendering decision").
# ----------------------------------------------------------------------
function bilinear_w(ζ::Complex)
    rx, iy_ = real(ζ), imag(ζ)
    # In-bounds check; cells outside the ζ-window get NaN.
    (rx < xs[1] || rx > xs[end] || iy_ < ys[1] || iy_ > ys[end]) &&
        return complex(NaN, NaN)
    # Find the index of the cell bracketing ζ.  `xs` is uniform, so
    # we can compute the index directly from the spacing.
    dx = xs[2] - xs[1]
    dy = ys[2] - ys[1]
    i = clamp(floor(Int, (rx - xs[1]) / dx) + 1, 1, NX - 1)
    j = clamp(floor(Int, (iy_ - ys[1]) / dy) + 1, 1, NY - 1)
    tx = (rx - xs[i]) / dx
    ty = (iy_ - ys[j]) / dy
    w00 = W[i, j]
    w10 = W[i+1, j]
    w01 = W[i, j+1]
    w11 = W[i+1, j+1]
    (isfinite(real(w00)) && isfinite(real(w10)) &&
     isfinite(real(w01)) && isfinite(real(w11))) ||
        return complex(NaN, NaN)
    return (1 - tx) * (1 - ty) * w00 +
           tx       * (1 - ty) * w10 +
           (1 - tx) * ty       * w01 +
           tx       * ty       * w11
end

# ----------------------------------------------------------------------
# Render: four columns A | B | C | D, matching FFW Fig 1's layout.
# Columns A, B, C share the full ζ-window aspect (Re-narrow,
# Im-tall — three strips of width 4π stacked).  Column D is a vertical
# stack of three square z-plane heatmaps, one per sheet s ∈ {-1, 0, 1}.
# ----------------------------------------------------------------------
const U_CAP   = 5.0       # heatmap cap — pole spikes go to infinity
const ABS_LO  = 0.0
# z-plane half-width: |z|_max = exp(Re_HI / 2) ≈ exp(2.5) ≈ 12.2 at
# Re_HI = 5.  Pad by 5% so the spiral arms at the boundary don't
# clip.
const Z_HALF  = 1.05 * exp(Re_HI / 2)
# Cartesian z-grid resolution for the Panel D heatmaps.
const Z_NX, Z_NY = 320, 320

# Tall portrait figure: four columns × 2 rows (title + plots).  Three
# `Im ζ` strips of width 4π each (12π ≈ 37.7) × Re range 8 → aspect
# ratio per ζ-panel ≈ 4.7 : 1.  Width 1800, height 1500 gives each
# panel ~270 × 1100 px (ζ panels) and 360 × 360 px (z-sheet panels).
fig = Figure(size = (1800, 1500))
Label(fig[0, 1:4],
      "FFW 2017 Fig 1 — P_III three sheets, (α,β,γ,δ) = (−1/2, −1/2, 1, −1), IC (z, u, u′) = (1, 1/4, 1)";
      fontsize = 18, padding = (0, 0, 10, 10))

# ---- Panel A: Stage-1 node set (ζ-plane scatter) --------------------
axA = Axis(fig[1, 1];
           xlabel = "Re ζ", ylabel = "Im ζ",
           title  = "A. Stage-1 node set",
           subtitle = "$(length(sol.visited_z)) visited nodes;\n" *
                      "R(ζ) = max(0.1, (8 − Re ζ)/20)",
           titlesize = 14, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (Re_LO - 0.2, Re_HI + 0.2, Im_LO - 0.3, Im_HI + 0.3))
# Dashed strip dividers at Im ζ = ±2π (FFW md:103 sheet boundaries).
for hline in (-2π, 2π)
    lines!(axA, [Re_LO - 0.2, Re_HI + 0.2], [hline, hline];
           color = (:gray, 0.6), linewidth = 0.8, linestyle = :dash)
end
scatter!(axA,
         Float64.(real.(sol.visited_z)),
         Float64.(imag.(sol.visited_z));
         color = :black, markersize = 1.4)

# ---- Panel B: Stage-1 path tree (visited_parent edges) --------------
axB = Axis(fig[1, 2];
           xlabel = "Re ζ", ylabel = "Im ζ",
           title  = "B. Stage-1 adaptive Padé steps",
           subtitle = "tree edges from visited_parent;\n" *
                      "FFW Fig 1 column 2",
           titlesize = 14, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (Re_LO - 0.2, Re_HI + 0.2, Im_LO - 0.3, Im_HI + 0.3))
for hline in (-2π, 2π)
    lines!(axB, [Re_LO - 0.2, Re_HI + 0.2], [hline, hline];
           color = (:gray, 0.6), linewidth = 0.8, linestyle = :dash)
end
seg = Point2f[]
for k in 2:length(sol.visited_z)
    p = sol.visited_parent[k]
    p == 0 && continue
    push!(seg, Point2f(real(sol.visited_z[p]), imag(sol.visited_z[p])))
    push!(seg, Point2f(real(sol.visited_z[k]), imag(sol.visited_z[k])))
end
linesegments!(axB, seg; color = (:black, 0.7), linewidth = 0.4)
scatter!(axB, [Float64(real(ζ0))], [Float64(imag(ζ0))];
         color = :red, markersize = 9)

# ---- Panel C: |w(ζ)| heatmap over the full ζ-window -----------------
axC = Axis(fig[1, 3];
           xlabel = "Re ζ", ylabel = "Im ζ",
           title  = "C. |w(ζ)| in the ζ-plane",
           subtitle = "three sheets stacked;\nFFW Fig 1 column 3",
           titlesize = 14, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (Re_LO, Re_HI, Im_LO, Im_HI))
absW = clamp.(abs.(W), ABS_LO, U_CAP)
heatmap!(axC, xs, ys, absW;
         colormap = :viridis, colorrange = (ABS_LO, U_CAP))
for hline in (-2π, 2π)
    lines!(axC, [Re_LO, Re_HI], [hline, hline];
           color = (:white, 0.75), linewidth = 0.6, linestyle = :dash)
end

# ---- Panel D: |u(z)| on three z-plane sheets (vertical stack) -------
panelD_layout = GridLayout(fig[1, 4])
Label(panelD_layout[1, 1, Top()],
      "D. |u(z)| on three Riemann sheets (s = +1, 0, −1)";
      fontsize = 14, padding = (0, 0, 0, 6))

z_xs = collect(range(-Z_HALF, Z_HALF; length = Z_NX))
z_ys = collect(range(-Z_HALF, Z_HALF; length = Z_NY))

# Sheets stacked top-down s = +1, 0, -1 so that the order matches the
# ζ-plane strip stacking in Panel C (Im ζ decreases downward in the
# plot — but PIII strips for s = +1 are at the TOP of the ζ window,
# s = 0 in the MIDDLE, s = -1 at the BOTTOM).
for (row, s) in enumerate((1, 0, -1))
    # For each Cartesian z-cell, compute the ζ on this sheet:
    # ζ = 2 (log|z| + i arg z) + 4π i s
    # then bilinear-interpolate w(ζ) from the ζ-grid, then u = w/z.
    U_sheet = fill(NaN, Z_NX, Z_NY)
    for jx in 1:Z_NY, ix in 1:Z_NX
        z = complex(z_xs[ix], z_ys[jx])
        # Skip a small disc around the origin — `log z` singular,
        # sheet boundary ambiguous.
        abs(z) < exp(Re_LO / 2) - 0.01 && continue
        ζcell = 2 * complex(log(abs(z)), angle(z)) + 4π * im * s
        w = bilinear_w(ζcell)
        isfinite(real(w)) || continue
        u = w / z
        absu = abs(u)
        isfinite(absu) || continue
        U_sheet[ix, jx] = clamp(absu, ABS_LO, U_CAP)
    end
    axD = Axis(panelD_layout[row, 1];
               xlabel = row == 3 ? "Re z" : "",
               ylabel = "Im z",
               title  = "sheet s = $s",
               titlesize = 12,
               aspect = DataAspect(),
               limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
    heatmap!(axD, z_xs, z_ys, U_sheet;
             colormap = :viridis, colorrange = (ABS_LO, U_CAP),
             nan_color = :transparent)
    # Origin dot + faint axes for orientation.
    lines!(axD, [-Z_HALF, Z_HALF], [0.0, 0.0];
           color = (:white, 0.25), linewidth = 0.3)
    lines!(axD, [0.0, 0.0], [-Z_HALF, Z_HALF];
           color = (:white, 0.25), linewidth = 0.3)
end

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
