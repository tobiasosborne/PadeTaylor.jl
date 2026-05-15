# figures/ffw2017_fig_4.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, **Figure 4** — the
# tronquée P_V solution on three sheets of its Riemann surface.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:200-240
#         (Figure 4 caption + §3.1 "A tronquée P_V solution" + Table 4).
#
# ## The problem
#
# Per FFW md:209: "for a given set of parameters `α, β, γ, δ` with
# `δ ≠ 0`, there is a unique tronquée P_V solution with `u(z) ~ -1`,
# `z → ∞` for `-π < arg z < π`.  The solution in Figure 4 is such a
# solution and its ICs were obtained by substituting `u = Σ a_n z^{-n}`
# into P_V and evaluating the truncated expansion far out on the
# positive real axis."
#
# Per FFW md:236 (Figure 4 caption, verbatim):
# > A tronquée P_V solution in the ζ and z planes with parameters
# > `(α, β, γ, δ) = (1, 0, 1/4, -1/2)` and ICs `u(30) = -1.05294551349665`
# > and `u'(30) = 2.47019460566845e-3` in the z-plane.  The third column
# > is a phase portrait of the solution shown in the second column.
# > The symmetry-based error estimates for the solution on sheets 0-2
# > are **3e-10, 7e-7 and 1e-6**, respectively.
#
# The IC values are FFW's published precomputed asymptotic-series
# values; we use them verbatim (the asymptotic-series reconstruction is
# a separate concern — see the deferred bead at end of worklog 038).
#
# ## PV transform
#
# Per FFW md:45-47, the PV exponential transform `z = exp(ζ)`,
# `u(z) = w(ζ)` (no rescale, unlike PIII) sends the fixed singularity
# `z = 0` out of the finite plane and yields the modified meromorphic
# equation `P̃_V`:
#
#     w'' = (1/(2w) + 1/(w-1)) (w')² + (w-1)² (α w + β/w)
#         + γ e^ζ w + δ e^(2ζ) w (w+1) / (w-1)
#
# Each strip `Im ζ ∈ (-π + 2π s, π + 2π s]` maps to one sheet of the
# multi-valued z-plane Riemann surface (PV strip width 2π per FFW md:103).
#
# ## ζ-window
#
# The IC at `z = 30` maps to `ζ₀ = log 30 ≈ 3.4012`, far inside the
# high-Re-ζ pole-rich region (vs Fig 6's IC at `ζ₀ = 0`).  We render:
#
#     Re ζ ∈ [-1, 4]
#     Im ζ ∈ (-π + 0.01, 5π - 0.01]
#     sheet 0: Im ζ ∈ (-π,  π],   z-arg ∈ (-π,  π]
#     sheet 1: Im ζ ∈ ( π, 3π],   z-arg ∈ ( π, 3π]
#     sheet 2: Im ζ ∈ (3π, 5π],   z-arg ∈ (3π, 5π]
#
# corresponding to `|z| ∈ [e^{-1}, e^4] ≈ [0.37, 54.6]` in the z-plane,
# repeated three times around the branch point.  The IC at `Re ζ ≈ 3.4`
# is interior to the high end of this window — the walker propagates
# **inward** (toward decreasing `Re ζ`) into the pole-free sector.
#
# ## Tronquée signature
#
# Per FFW md:209: `u(z) ~ -1` as `z → ∞` on `-π < arg z < π` (sheet 0).
# The pole-free sector of the principal sheet is therefore the
# **smooth blue plain** on the right side of the ζ-plane heatmap
# (high `Re ζ`, sheet 0, `Im ζ ∈ (-π, π]`) where `|w(ζ)| ≈ 1`.  Sheets
# 1 and 2 (`Im ζ ∈ (π, 3π]` and `(3π, 5π]`) are not in the asymptotic
# sector and exhibit the dense pole/zero pattern characteristic of PV.
#
# Per FFW md:207 (Table 4 + caption): "all the zeros of the solution
# have double multiplicity.  They can be identified on the phase
# portraits in the third column as points around which each color is
# assumed twice in the order indicated by the color wheel above
# Figure 3 (red→yellow→green etc. for a counterclockwise traversal
# around the point)."  This double-multiplicity zero signature is the
# qualitative feature distinguishing PV tronquée from the generic PV
# solution of Figure 6.
#
# ## Method
#
#   * `PadeTaylorProblem` built directly on the FFW md:47 `P̃_V` RHS
#     via `CoordTransforms.pV_transformed_rhs(α, β, γ, δ)`.  The IC
#     `(z₀, u, u') = (30, -1.05294551349665, 2.47019460566845e-3)`
#     maps under `pV_z_to_ζ` to `(ζ₀, w, w')`:
#
#         ζ₀ = log 30
#         w  = u = -1.05294551349665
#         w' = z u' = 30 · 2.47019460566845e-3 = 0.0741058381700535
#
#   * Stage 1 walks via `path_network_solve` with
#     `step_size_policy = :adaptive_ffw` (ADR-0011, worklog 034) at
#     `adaptive_tol = 1e-10`.
#
#   * Stage-1 node placement uses the FFW md:72-style non-uniform
#     `node_separation` kwarg (ADR-0012, worklog 035):
#
#         R(ζ) = max(0.05, (4 - Re ζ) / 10)
#
#     — Fig 6's prescription, which gives `R(3.4) ≈ 0.06` at the IC
#     (the high-density side) and `R(-1) = 0.5` (the low-density side).
#     The walker propagates outward in **both** directions from the IC,
#     so the node-separation function being densest at the IC end is
#     fine — the FW md:160-164 tree fans out from the root.
#
#   * Stage-1 / Stage-2 split (worklog 037 pattern): walker walks ONLY
#     the rectangular Stage-1 target grid (~870 targets), and we
#     manually evaluate Stage-2 from the visited tree at a finer
#     uniform lattice (~120 × 220 = 26,400 cells).  This keeps the
#     walker away from cells that fall directly on pole spirals.
#
# ## Rendering — three columns, three sheets, nine panels
#
# Column layout matches FFW Fig 4's published layout:
#
#   * **Column 1 — `|w(ζ)|` ζ-plane heatmap** per sheet.  Modulus
#     surface over each `Im ζ` strip.  Sheet 0 shows the smooth
#     pole-free plain at high `Re ζ`; sheets 1 & 2 show dense pole
#     spikes.
#
#   * **Column 2 — `|u(z)|` z-plane modulus** per sheet.  Cartesian
#     `Re z, Im z ∈ [-|z|_max, |z|_max]` grid (~300×300), with
#     **bilinear interpolation from the ζ-grid** (the worklog 037
#     B2 rendering upgrade — no polar-scatter spoke artefacts).  For
#     each Cartesian z-cell on sheet `s`, compute
#     `ζ = log|z| + i arg z + 2π i s` (principal-branch `arg z`) and
#     interpolate `w(ζ)`.  PV has no rescale, so `u(z) = w(ζ)` directly.
#
#   * **Column 3 — `arg u(z)` z-plane phase portrait** per sheet.  Same
#     Cartesian z-grid; color-wheel-encoded `arg u ∈ [-π, π]` via HSV
#     hue.  Double-multiplicity zeros of the PV tronquée (FFW md:207)
#     appear as points where each color is traversed **twice** for a
#     single counterclockwise loop around the zero (vs once for a
#     simple zero / pole).
#
# ## Acceptance
#
# Per `docs/figure_catalogue.md §5` (Fig 4 row, T4) and FFW md:236:
# per-sheet error within one order of magnitude of FFW Table-2-style
# `3e-10 / 7e-7 / 1e-6` for sheets 0/1/2.  `test/ffw_fig_4_test.jl`
# pins this quantitatively (FF4.1.1 IC round-trip; FF4.1.2-4
# per-sheet loose-vs-tight self-cross-check; FF4.1.5 pole-free sector
# check; FF4.1.6 pole-density gradient).

using PadeTaylor
using PadeTaylor.CoordTransforms: pV_transformed_rhs, pV_z_to_ζ, pV_ζ_to_z
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters + IC (FFW md:236, verbatim)
# ----------------------------------------------------------------------
const α, β, γ, δ = 1.0, 0.0, 0.25, -0.5
const Z0  = 30.0 + 0.0im
const U0  = -1.05294551349665 + 0.0im
const UP0 =  2.47019460566845e-3 + 0.0im

const ORDER       = 30           # FW/FFW working Taylor order
const ADAPT_TOL   = 1.0e-10      # FFW md:91 controller tolerance
const K_CONS      = 1.0e-3       # FFW md:91 conservative factor
const MAX_RESC    = 50

# ζ-window: three contiguous strips of vertical width 2π stacked.
const Re_LO, Re_HI = -1.0, 4.0
const Im_LO, Im_HI = -π + 0.01, 5π - 0.01       # epsilon-shrunk to keep
                                                 # boundary lookups in-tree

# Stage-1 + Stage-2 lattice densities (Stage-2 finer than Stage-1 grid).
const NX, NY = 120, 220                          # ζ-grid Stage-2 cells
const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_4.png")

# FFW md:72-style node-separation function, scaled to our `Re ζ` range.
# At Re ζ = -1 → R = 0.5 (FW 2011 default); at Re ζ = 4 → R = 0.05
# (floor); typical step `R` correlates with local pole density.
R(ζ) = max(0.05, (4.0 - real(ζ)) / 10.0)

# ----------------------------------------------------------------------
# Map the IC into the ζ-frame.  pV_z_to_ζ: ζ = log z, w = u, w' = z u'.
# At (z, u, u') = (30, -1.0529..., 2.470e-3) this gives
# (ζ, w, w') = (log 30, -1.0529..., 30 · 2.470e-3 = 0.0741...).
# ----------------------------------------------------------------------
const ζ0, w0, wp0 = pV_z_to_ζ(Z0, U0, UP0)
@printf("FFW 2017 Fig 4: (α,β,γ,δ) = (%.2f, %.2f, %.2f, %.2f)\n", α, β, γ, δ)
@printf("                IC (z₀, u, u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("                → (ζ₀, w, w') = (%s, %s, %s)\n", ζ0, w0, wp0)

# ----------------------------------------------------------------------
# Build Stage-1 targets (rectangular raster, R-spaced) and the Stage-2
# uniform evaluation lattice over the full ζ-window.
# ----------------------------------------------------------------------
function build_stage1_targets()
    tgs = ComplexF64[]
    re = Re_LO
    while re ≤ Re_HI + 1e-12
        rh = R(re + 0.0im)
        imv = Im_LO
        while imv ≤ Im_HI + 1e-12
            push!(tgs, complex(re, imv))
            imv += rh
        end
        re += rh
    end
    return tgs
end

stage1_targets = build_stage1_targets()
@printf("Stage-1 target count: %d (R-spaced)\n", length(stage1_targets))

xs = collect(range(Re_LO, Re_HI; length = NX))
ys = collect(range(Im_LO, Im_HI; length = NY))
eval_lattice = ComplexF64[complex(x, y) for y in ys for x in xs]
@printf("Stage-2 eval lattice: %d × %d = %d cells\n",
        NX, NY, length(eval_lattice))

# ----------------------------------------------------------------------
# Solve the FFW md:47 `P̃_V` problem.  Walk only the Stage-1 target
# set (worklog 037 pattern): adding the full ~26k-cell eval lattice to
# the walker's target list forces it to navigate to every cell, which
# at sheets 1-2 means cells on pole spirals where the 5-direction
# wedge can't reliably thread.  Stage-2 then reads the visited tree
# manually via `stage2_eval` below — same logic as the walker's
# internal Stage-2 at PathNetwork.jl:495-506.
# ----------------------------------------------------------------------
f = pV_transformed_rhs(α, β, γ, δ)
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
@printf("  Stage 1 in %.2f s; %d visited tree nodes.\n",
        time() - t0, length(sol.visited_z))

# ----------------------------------------------------------------------
# Manual Stage 2 (worklog 037 helper).  Nearest-visited lookup + local
# Padé evaluation; out-of-disc cells return NaN.
# ----------------------------------------------------------------------
function stage2_eval(sol, points)
    nv = length(sol.visited_z)
    out = Vector{ComplexF64}(undef, length(points))
    for (i, zf) in enumerate(points)
        best_k = 1
        best_d = abs(zf - sol.visited_z[1])
        for k in 2:nv
            d = abs(zf - sol.visited_z[k])
            if d < best_d
                best_d = d
                best_k = k
            end
        end
        hv = sol.visited_h[best_k]
        if best_d > hv
            out[i] = complex(NaN, NaN)
        else
            t = (zf - sol.visited_z[best_k]) / hv
            out[i] = PadeTaylor.PathNetwork._evaluate_pade(
                sol.visited_pade[best_k], t)
        end
    end
    return out
end

t0 = time()
eval_u = stage2_eval(sol, eval_lattice)
@printf("  Stage 2 (manual) in %.2f s\n", time() - t0)
ncov = count(isfinite, abs.(eval_u))
@printf("  Stage-2 coverage: %d / %d (%.1f%%) finite\n",
        ncov, length(eval_u), 100 * ncov / length(eval_u))

# ----------------------------------------------------------------------
# Reshape into the (NX, NY) grid.  eval_lattice was built y-outer /
# x-inner, matching column-major reshape(_, NX, NY).
# ----------------------------------------------------------------------
W = reshape(eval_u, NX, NY)

# ----------------------------------------------------------------------
# Bilinear interpolation of `w(ζ)` from the regular ζ-grid `W[i, j]`.
# Returns NaN+NaN·im when ζ falls outside the grid or any corner is
# non-finite.  The "smooth-fill" helper from worklog 037 — drives the
# z-plane panel rendering, avoiding polar-scatter spoke artefacts.
# ----------------------------------------------------------------------
function bilinear_w(ζ::Complex)
    rx, iy_ = real(ζ), imag(ζ)
    (rx < xs[1] || rx > xs[end] || iy_ < ys[1] || iy_ > ys[end]) &&
        return complex(NaN, NaN)
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
# Sheet definitions.  PV strip width is 2π (vs PIII's 4π).
#   sheet 0: Im ζ ∈ (-π,  π],   z-arg ∈ (-π,  π]   (principal)
#   sheet 1: Im ζ ∈ ( π, 3π],   z-arg ∈ ( π, 3π]
#   sheet 2: Im ζ ∈ (3π, 5π],   z-arg ∈ (3π, 5π]
# ----------------------------------------------------------------------
sheet_im_ranges = [(-π + 2π*s, π + 2π*s) for s in 0:2]

# ----------------------------------------------------------------------
# Render: three columns × three sheets = nine panels.
# ----------------------------------------------------------------------
const U_CAP   = 5.0       # |w|/|u| cap — pole spikes go to infinity
const ABS_LO  = 0.0
# z-plane half-width: |z|_max = exp(Re_HI) ≈ exp(4) ≈ 54.6.
const Z_RMAX  = exp(Re_HI)
const Z_HALF  = 1.05 * Z_RMAX
# Cartesian z-grid resolution per sheet panel.
const Z_NX, Z_NY = 300, 300

fig = Figure(size = (1500, 1100))
Label(fig[0, 1:3],
      "FFW 2017 Fig 4 — tronquée P_V, (α,β,γ,δ) = (1, 0, 1/4, −1/2), IC (z, u, u′) = (30, −1.0529, 2.47e−3)";
      fontsize = 16, padding = (0, 0, 6, 6))

# Cartesian z-grid axes shared across sheet panels.
z_xs = collect(range(-Z_HALF, Z_HALF; length = Z_NX))
z_ys = collect(range(-Z_HALF, Z_HALF; length = Z_NY))

for (i, (im_lo, im_hi)) in enumerate(sheet_im_ranges)
    sheet_label = "sheet $(i-1)"
    sheet_idx   = i - 1

    # Indices of `ys` falling in this strip (for the ζ-plane heatmap).
    ym = findall(y -> im_lo ≤ y ≤ im_hi, ys)
    sub_ys = ys[ym]
    sub_W  = W[:, ym]
    sub_absW = abs.(sub_W)
    capped_absW = clamp.(sub_absW, ABS_LO, U_CAP)

    # ---------------- Column 1: |w(ζ)| in the ζ-plane ----------------
    axA = Axis(fig[i, 1];
               xlabel = i == 3 ? "Re ζ" : "",
               ylabel = "Im ζ",
               title  = "|w(ζ)|, $sheet_label",
               titlesize = 13,
               aspect = DataAspect(),
               limits = (Re_LO, Re_HI, im_lo, im_hi))
    heatmap!(axA, xs, sub_ys, capped_absW;
             colormap = :viridis, colorrange = (ABS_LO, U_CAP),
             nan_color = :transparent)

    # ---------------- Column 2: |u(z)| in the z-plane ----------------
    # Cartesian-resample bilinear-interp (worklog 037 B2 upgrade).
    # PV transform: u(z) = w(ζ) with ζ = log|z| + i arg z + 2π i s.
    U_sheet  = fill(NaN, Z_NX, Z_NY)
    Phi_sheet = fill(NaN, Z_NX, Z_NY)
    for jx in 1:Z_NY, ix in 1:Z_NX
        z = complex(z_xs[ix], z_ys[jx])
        # Skip a small disc around the origin: log(0) singular.
        abs(z) < exp(Re_LO) - 0.01 && continue
        # Skip cells outside |z| ≤ exp(Re_HI): outside the walked window.
        abs(z) > exp(Re_HI) + 0.01 && continue
        ζcell = complex(log(abs(z)), angle(z)) + 2π * im * sheet_idx
        w = bilinear_w(ζcell)
        isfinite(real(w)) || continue
        # PV: u = w, no rescale.
        u = w
        absu = abs(u)
        isfinite(absu) || continue
        U_sheet[ix, jx]   = clamp(absu, ABS_LO, U_CAP)
        Phi_sheet[ix, jx] = Float64(angle(u))
    end

    axB = Axis(fig[i, 2];
               xlabel = i == 3 ? "Re z" : "",
               ylabel = "Im z",
               title  = "|u(z)|, $sheet_label",
               titlesize = 13,
               aspect = DataAspect(),
               limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
    heatmap!(axB, z_xs, z_ys, U_sheet;
             colormap = :viridis, colorrange = (ABS_LO, U_CAP),
             nan_color = :transparent)
    # Faint origin axes for orientation.
    lines!(axB, [-Z_HALF, Z_HALF], [0.0, 0.0];
           color = (:white, 0.25), linewidth = 0.3)
    lines!(axB, [0.0, 0.0], [-Z_HALF, Z_HALF];
           color = (:white, 0.25), linewidth = 0.3)

    # ---------------- Column 3: arg u(z) phase portrait ---------------
    axC = Axis(fig[i, 3];
               xlabel = i == 3 ? "Re z" : "",
               ylabel = "Im z",
               title  = "arg u(z), $sheet_label",
               titlesize = 13,
               aspect = DataAspect(),
               limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
    heatmap!(axC, z_xs, z_ys, Phi_sheet;
             colormap = :hsv, colorrange = (-π, π),
             nan_color = :transparent)
    lines!(axC, [-Z_HALF, Z_HALF], [0.0, 0.0];
           color = (:black, 0.25), linewidth = 0.3)
    lines!(axC, [0.0, 0.0], [-Z_HALF, Z_HALF];
           color = (:black, 0.25), linewidth = 0.3)
end

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
