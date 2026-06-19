# figures/heun_closed_form_reduction.jl
#
# ============================================================================
# HEADLINE FIGURE 3 (P1 epic, bead padetaylor-e4z): a RIGOROUS validation
# portrait — PadeTaylor's `heun_g` reproducing a KNOWN CLOSED FORM as it
# navigates the Heun singular points z = 0, 1, a = 2.
#
# THE IDENTITY (the oracle)
# -------------------------
# On the exceptional parameter slice (a, q) = (2, 2) the general Heun
# function degenerates to a Gauss `₂F₁`, which in turn is elementary:
#
#     HeunG(a=2, q=2, α=1, β=2, γ=3/2, δ=1; z)
#         = ₂F₁(1/2, 1; 3/2; 2z − z²)
#         = atanh(√w) / √w,          w = 2z − z².
#
#   * The HeunG → ₂F₁ reduction is DLMF §31.7 (the exceptional-(a,q) class;
#     `docs/heun_research.md:385-406`).  The specific (a, q, α, β, γ, δ)
#     slice and its `2z − z²` argument substitution are pinned + verified
#     to |diff| < 1e-44 by the INDEPENDENT mpmath recurrence in
#     `test/_oracle_corpus_heun.jl:18` (CHN.1).
#   * The ₂F₁(1/2, 1; 3/2; w) = atanh(√w)/√w elementary closed form is
#     DLMF 15.4.x, transcribed and asserted (not merely pinned) in
#     `test/corpus_heun_test.jl:59-61` via `_f2f1_half1_3half`.
#
# This is the BEST oracle in the Heun corpus: fully closed-form, in
# standard `ComplexF64` arithmetic, with NO external dependency.  It lets
# us draw a true validation portrait — not "does heun_g render", but "does
# heun_g equal a function we can compute another way, everywhere."
#
# WHY THE FIGURE IS RIGOROUS (Rule 5 spirit)
# ------------------------------------------
# We do not merely render.  The script:
#   (1) re-asserts the CHN.1 pinned values at z = 0.3, 0.5 (rtol 1e-12) and
#       z = 0.7 (rtol 1e-9) against the `_oracle_corpus_heun.jl` constants,
#       failing LOUD (Rule 1) if heun_g drifts from the suite's pins;
#   (2) reports the MAX |heun_g − closed_form| over the principal disc
#       |z| ≤ 0.95 and asserts it is < 1e-10 (the walk floor on the disc);
#   (3) confirms the PNG is materially produced (> 50 KB).
# All three numbers are PRINTED.
#
# RENDERING
# ---------
# Reuses the Wegert phase-portrait (domain-colouring) idiom of
# `figures/heun_complex_portrait.jl`: hue = arg f, brightness modulated by a
# gentle log-magnitude ripple + zero-shading.  Singular points marked with a
# white X; Wolfram-convention cuts shown downward (matching the walker's
# `branch_cut_angles` / `enforce_real_axis_symmetry` principal-sheet choice).
#
# THE WALK
# --------
# `heun_g` per single point is the public API, but evaluating it at ~150²
# grid points one-at-a-time would re-walk from z₀ = 0.05 for every pixel.
# Instead — exactly as `heun_complex_portrait.jl` does — we run ONE shared
# `path_network_solve` over the whole grid, replicating verbatim what
# `heun_g` does internally for real parameters (Frobenius IC at z₀ = 0.05,
# `enforce_real_axis_symmetry = true`, same order / h / n_taylor defaults).
# The CHN.1 acceptance checks above call the PUBLIC `heun_g` directly, so the
# shared-walk grid is cross-checked against the canonical single-point path.
# ============================================================================

using PadeTaylor
using PadeTaylor.Heun: _heun_g_frobenius_at_0, _heun_g_rhs
using PadeTaylor: PadeTaylorProblem, path_network_solve
using CairoMakie
using Printf

# Pull the suite's pinned oracle constants (HEUNG_2F1_Z0p3/Z0p5/Z0p7) from
# the SAME file the test suite uses — never hardcode from memory (Law 1).
include(joinpath(@__DIR__, "..", "test", "_oracle_corpus_heun.jl"))

const OUTPNG = joinpath(@__DIR__, "output", "heun_closed_form_reduction.png")
mkpath(dirname(OUTPNG))

# ----------------------------------------------------------------------
# Parameters — the CHN.1 exceptional slice (match the call site exactly).
# test/corpus_heun_test.jl:77 calls heun_g(2,2,1,2,3//2,1; ...): γ is the
# RATIONAL 3//2 there.  Float64(3//2) == 1.5 exactly (dyadic), so the walk
# arithmetic is identical; we keep 3//2 in the public-API acceptance calls
# below so they mirror the suite byte-for-byte.
# ----------------------------------------------------------------------
const G_a, G_q, G_α, G_β, G_γ, G_δ = 2.0, 2.0, 1.0, 2.0, 1.5, 1.0

# Singular points of this HeunG: z = 0, z = 1, z = a = 2.
const SING = (0.0 + 0im, 1.0 + 0im, 2.0 + 0im)

# The CLOSED FORM (panel 2 + error oracle): standard ComplexF64 arithmetic.
#   ₂F₁(1/2, 1; 3/2; w) = atanh(√w) / √w,   w = 2z − z².
# At z = 0, w = 0 → 0/0; the removable-singularity limit is 1 (HeunG(0)=1).
function closed_form(z)
    w = 2z - z^2
    abs(w) < 1e-14 && return one(ComplexF64)        # atanh(√w)/√w → 1 as w→0
    sw = sqrt(complex(w))
    return atanh(sw) / sw
end

# ----------------------------------------------------------------------
# (1) ACCEPTANCE — re-assert CHN.1 against the suite's pins (Rule 5).
#     Calls the PUBLIC heun_g, single-point, exactly as the suite does.
# ----------------------------------------------------------------------
println("="^72)
println("CHN.1 reproduction — public heun_g(2,2,1,2,3//2,1; z) vs pinned oracle")
println("="^72)

const CHN1_POINTS = ((0.3, HEUNG_2F1_Z0p3, 1e-12),
                     (0.5, HEUNG_2F1_Z0p5, 1e-12),
                     (0.7, HEUNG_2F1_Z0p7, 1e-9))

chn1_relerrs = Float64[]
for (zr, oracle, rtol) in CHN1_POINTS
    val = real(heun_g(G_a, G_q, G_α, G_β, 3 // 2, G_δ; z = zr + 0im))
    rel = abs(val - oracle) / abs(oracle)
    push!(chn1_relerrs, rel)
    @printf("  z = %.1f : heun_g = %.15f  oracle = %.15f  rel-err = %.3e  (rtol %.0e)\n",
            zr, val, oracle, rel, rtol)
    if rel > rtol
        error("CHN.1 ACCEPTANCE FAILED at z=$zr: rel-err $rel exceeds rtol $rtol. " *
              "heun_g has drifted from the suite's pinned 2F1 reduction oracle — " *
              "this is a numerical breakdown, not a tolerance to relax (Rule 1, Rule 2).")
    end
end
println("  CHN.1 reproduction: PASS (all 3 points within pinned rtol).")

# ----------------------------------------------------------------------
# Grid + the single shared HeunG walk (replicates heun_g internals).
# ----------------------------------------------------------------------
const RE_MIN, RE_MAX = -1.5, 2.5
const IM_MIN, IM_MAX = -1.5, 1.5
const NX, NY = 180, 144
const ORDER, HSTEP, NTAYL, EPS_START = 30, 0.05, 60, 0.05

xs = range(RE_MIN, RE_MAX; length = NX)
ys = range(IM_MIN, IM_MAX; length = NY)
grid = ComplexF64[complex(x, y) for y in ys for x in xs]
@printf("\nGrid: %d × %d = %d targets over Re∈[%.1f,%.1f] Im∈[%.1f,%.1f]\n",
        NX, NY, length(grid), RE_MIN, RE_MAX, IM_MIN, IM_MAX)

@printf("Walking HeunG (single shared path_network_solve)...")
t0 = time()
z0 = ComplexF64(EPS_START)
u0, up0, _ = _heun_g_frobenius_at_0(complex(G_a), complex(G_q), complex(G_α),
                                    complex(G_β), complex(G_γ), complex(G_δ),
                                    z0, NTAYL)
f_g  = _heun_g_rhs(complex(G_a), complex(G_q), complex(G_α),
                   complex(G_β), complex(G_γ), complex(G_δ))
prob = PadeTaylorProblem(f_g, (u0, up0), (z0, ComplexF64(RE_MAX + 0.0im)); order = ORDER)
sol  = path_network_solve(prob, grid;
                          h = HSTEP, order = ORDER,
                          enforce_real_axis_symmetry = true,
                          extrapolate = true,
                          max_steps_per_target = 2000)
@printf(" %.1fs (%d visited)\n", time() - t0, length(sol.visited_z))

G_vals = reshape(sol.grid_u, NX, NY)

# Closed form on the same grid.
CF_vals = Matrix{ComplexF64}(undef, NX, NY)
for j in 1:NY, i in 1:NX
    CF_vals[i, j] = closed_form(grid[(j - 1) * NX + i])
end

# ----------------------------------------------------------------------
# (2) ACCEPTANCE — MAX |heun_g − closed_form| over the principal disc.
#     |z| ≤ 0.95: inside the Frobenius / first-continuation region, away
#     from the branch points at z=1, z=2; the two must agree to the walk
#     floor.  Assert < 1e-10 and PRINT the actual value (Rule 5).
# ----------------------------------------------------------------------
disc_max, disc_argz, n_disc = let dmax = 0.0, dz = 0.0 + 0im, n = 0
    for j in 1:NY, i in 1:NX
        z = grid[(j - 1) * NX + i]
        abs(z) ≤ 0.95 || continue
        g, c = G_vals[i, j], CF_vals[i, j]
        (isfinite(g) && isfinite(c)) || continue
        n += 1
        d = abs(g - c)
        if d > dmax
            dmax = d; dz = z
        end
    end
    (dmax, dz, n)
end
@printf("\nMAX |heun_g − closed_form| over |z| ≤ 0.95 (%d grid points): %.3e", n_disc, disc_max)
@printf("  (worst at z = %.3f%+.3fi)\n", real(disc_argz), imag(disc_argz))
if disc_max ≥ 1e-10
    error("DISC ACCEPTANCE FAILED: max |heun_g − closed_form| = $disc_max on " *
          "|z| ≤ 0.95 is ≥ 1e-10.  The shared walk disagrees with the elementary " *
          "closed form on the principal disc — a numerical breakdown (Rule 1).")
end
println("  Disc agreement: PASS (machine-precision identity inside |z| ≤ 0.95).")

# ----------------------------------------------------------------------
# Domain colouring (Wegert phase portrait) — reused from heun_complex_portrait.jl.
# ----------------------------------------------------------------------
function hsv2rgb(h::Real, s::Real, v::Real)
    h6 = mod(h, 1.0) * 6
    i = floor(h6); f = h6 - i
    p = v * (1 - s); q = v * (1 - s * f); t = v * (1 - s * (1 - f))
    r, g, b = if     i < 1; (v, t, p)
              elseif i < 2; (q, v, p)
              elseif i < 3; (p, v, t)
              elseif i < 4; (p, q, v)
              elseif i < 5; (t, p, v)
              else          (v, p, q)
              end
    return RGBf(clamp(r, 0, 1), clamp(g, 0, 1), clamp(b, 0, 1))
end

function color_grid(vals)
    img = Matrix{RGBf}(undef, NX, NY)
    for j in 1:NY, i in 1:NX
        v = vals[i, j]
        if !isfinite(v)
            img[i, j] = RGBf(0.10, 0.10, 0.10); continue
        end
        h = mod(angle(v) / (2π) + 1.0, 1.0)
        r = abs(v)
        L = log10(r + 1e-12)
        ripple = 1.0 - 0.08 * cos(2π * L * 2)
        zero_shade = r < 0.1 ? (r / 0.1) * 0.7 + 0.3 : 1.0
        img[i, j] = hsv2rgb(h, 1.0, clamp(1.0 * ripple * zero_shade, 0.0, 1.0))
    end
    return img
end

# log10 pointwise error grid (panel 3).  Mask non-finite / branch regions
# (where either side blew up) to NaN so viridis leaves them blank.
function log_err_grid(g_vals, c_vals; floor_lo = -16.0)
    out = Matrix{Float64}(undef, NX, NY)
    for j in 1:NY, i in 1:NX
        g, c = g_vals[i, j], c_vals[i, j]
        if !(isfinite(g) && isfinite(c))
            out[i, j] = NaN; continue
        end
        d = abs(g - c)
        out[i, j] = d == 0 ? floor_lo : clamp(log10(d), floor_lo, 2.0)
    end
    return out
end

img_g   = color_grid(G_vals)
img_c   = color_grid(CF_vals)
logerr  = log_err_grid(G_vals, CF_vals)

function mark_sing!(ax, z; label = nothing)
    scatter!(ax, [real(z)], [imag(z)]; marker = :xcross, markersize = 20,
             color = :white, strokecolor = :black, strokewidth = 2)
    label === nothing || text!(ax, real(z), imag(z); text = label,
        offset = (9, 9), color = :white, fontsize = 12,
        strokecolor = :black, strokewidth = 1)
end
mark_sing_dark!(ax, z) = scatter!(ax, [real(z)], [imag(z)];
    marker = :xcross, markersize = 18, color = :black,
    strokecolor = :white, strokewidth = 1.5)

function mark_cut!(ax, from_z)
    lines!(ax, [real(from_z), RE_MAX + 0.05], [imag(from_z), imag(from_z)];
           color = (:white, 0.65), linestyle = :dash, linewidth = 1.5)
end

# Unit-disc outline (|z| = 1) — the region where the identity is exact.
unit_circle() = ([cos(t) for t in range(0, 2π; length = 200)],
                 [sin(t) for t in range(0, 2π; length = 200)])

# ----------------------------------------------------------------------
# Render: 1×3 panels.
# ----------------------------------------------------------------------
fig = Figure(size = (1750, 660), fontsize = 14)
Label(fig[0, 1:4],
      "HeunG closed-form reduction — PadeTaylor.jl `heun_g` vs " *
      "₂F₁(½,1;³⁄₂; 2z−z²) = atanh(√w)/√w   (w = 2z−z²;  a=2, q=2, α=1, β=2, γ=³⁄₂, δ=1)";
      fontsize = 21, padding = (0, 0, 8, 0))

ucx, ucy = unit_circle()

# Panel 1 — PadeTaylor heun_g phase portrait.
ax1 = Axis(fig[1, 1];
           title = "PadeTaylor  heun_g(2,2,1,2,³⁄₂,1; z)   (Wegert phase portrait)",
           xlabel = "Re z", ylabel = "Im z", titlesize = 13,
           aspect = DataAspect(), limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
image!(ax1, (RE_MIN, RE_MAX), (IM_MIN, IM_MAX), img_g)
lines!(ax1, ucx, ucy; color = (:white, 0.55), linewidth = 1.0, linestyle = :dot)
mark_sing!(ax1, SING[1]; label = "0")
mark_sing!(ax1, SING[2]; label = "1")
mark_sing!(ax1, SING[3]; label = "a=2")
mark_cut!(ax1, SING[2]); mark_cut!(ax1, SING[3])

# Panel 2 — closed-form phase portrait (independent ComplexF64 arithmetic).
ax2 = Axis(fig[1, 2];
           title = "Closed form  atanh(√w)/√w,  w = 2z−z²   (independent ℂ arithmetic)",
           xlabel = "Re z", ylabel = "Im z", titlesize = 13,
           aspect = DataAspect(), limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
image!(ax2, (RE_MIN, RE_MAX), (IM_MIN, IM_MAX), img_c)
lines!(ax2, ucx, ucy; color = (:white, 0.55), linewidth = 1.0, linestyle = :dot)
mark_sing!(ax2, SING[1]; label = "0")
mark_sing!(ax2, SING[2]; label = "1")
mark_sing!(ax2, SING[3]; label = "a=2")
mark_cut!(ax2, SING[2]); mark_cut!(ax2, SING[3])

# Panel 3 — log10 pointwise error heatmap.
ax3 = Axis(fig[1, 3];
           title = "log₁₀ |heun_g − closed form|   (machine-ε on the unit disc)",
           xlabel = "Re z", ylabel = "Im z", titlesize = 13,
           aspect = DataAspect(), limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
hm = heatmap!(ax3, xs, ys, logerr; colormap = :viridis, colorrange = (-16.0, 0.0))
lines!(ax3, ucx, ucy; color = (:white, 0.7), linewidth = 1.2, linestyle = :dot)
mark_sing_dark!(ax3, SING[1]); mark_sing_dark!(ax3, SING[2]); mark_sing_dark!(ax3, SING[3])
mark_cut!(ax3, SING[2]); mark_cut!(ax3, SING[3])
Colorbar(fig[1, 4], hm; label = "log₁₀ |Δ|", width = 14)

# Annotation: machine-precision agreement on the unit disc + honest growth.
# (@sprintf needs a single literal format string; build the caption in parts.)
caption = string(
    @sprintf("Inside the unit disc (dotted) the two agree to machine precision: max |heun_g − closed form| over |z| ≤ 0.95 = %.2e (< 1e-10).  ", disc_max),
    "Error grows HONESTLY toward the branch points z=1, z=2 (× marks; dashed = downward Wolfram-convention cuts).  ",
    "CHN.1 pins reproduced at z=0.3,0.5 (rtol 1e-12) and z=0.7 (rtol 1e-9).  ",
    @sprintf("Grid %d×%d; single shared path_network_solve walk.", NX, NY))
Label(fig[2, 1:4], caption;
      fontsize = 11, color = :gray30, padding = (0, 0, 0, 4))

colsize!(fig.layout, 4, Relative(0.03))

save(OUTPNG, fig)
sz = filesize(OUTPNG)
@printf("\nSaved: %s  (%.1f KB)\n", OUTPNG, sz / 1024)
if sz < 50 * 1024
    error("PNG ACCEPTANCE FAILED: $OUTPNG is only $(sz) bytes (< 50 KB) — " *
          "the figure did not render to a full image (Rule 1).")
end

println("="^72)
println("ALL ACCEPTANCE CHECKS PASSED.")
@printf("  CHN.1 rel-errs:  z=0.3 → %.3e | z=0.5 → %.3e | z=0.7 → %.3e\n",
        chn1_relerrs[1], chn1_relerrs[2], chn1_relerrs[3])
@printf("  Disc max-error (|z|≤0.95): %.3e\n", disc_max)
@printf("  PNG: %.1f KB\n", sz / 1024)
println("="^72)
