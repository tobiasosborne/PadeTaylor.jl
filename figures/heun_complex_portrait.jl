# figures/heun_complex_portrait.jl
#
# Domain-colour portraits of HeunG and HeunC on the complex plane.
# The headline visual for the Heun arc (ADR-0018, worklog 051).
#
# **Domain colouring** is the classic visual for meromorphic functions:
#   - hue  = arg f(z)         (full circle of colours represents -π .. π)
#   - value/brightness = magnitude indicator
#   - we use brightness = 1 - 0.6 · exp(-|f|/scale) + light/dark contour
#     stripes at log-spaced magnitude levels (the "checkerboard" you see
#     in Trefethen's complex-variables textbook)
#
# This makes both the BRANCH STRUCTURE (sharp colour discontinuities
# across cuts) and the POLE / ZERO STRUCTURE (bullseye colour swirls
# centred on each pole or zero) immediately legible.
#
# Two panels (1×2), large.  Singular points marked; Wolfram-convention
# cuts shown.  Computed via single shared path-network walks per function.

using PadeTaylor
using PadeTaylor.Heun: _heun_g_frobenius_at_0, _heun_g_rhs,
                      _heun_c_frobenius_at_0, _heun_c_rhs
using PadeTaylor: PadeTaylorProblem, path_network_solve
using CairoMakie
using Printf

# Manual HSV → RGB conversion to avoid namespace fights with ColorTypes.
function hsv2rgb(h::Real, s::Real, v::Real)
    h6 = mod(h, 1.0) * 6
    i = floor(h6); f = h6 - i
    p = v * (1 - s); q = v * (1 - s*f); t = v * (1 - s*(1-f))
    r, g, b = if     i < 1; (v, t, p)
              elseif i < 2; (q, v, p)
              elseif i < 3; (p, v, t)
              elseif i < 4; (p, q, v)
              elseif i < 5; (t, p, v)
              else          (v, p, q)
              end
    return RGBf(clamp(r, 0, 1), clamp(g, 0, 1), clamp(b, 0, 1))
end

const OUTPNG = joinpath(@__DIR__, "output", "heun_complex_portrait.png")
mkpath(dirname(OUTPNG))

# ----------------------------------------------------------------------
# Parameters — chosen for dramatic phase structure
# ----------------------------------------------------------------------

# HeunG: parameters chosen for visible phase structure across the plane.
# (a=2, q=2.0, α=1.5, β=1.5, γ=0.5, δ=0.5)
# Local exponents:
#   z=0: {0, 1-γ} = {0, 0.5}        → square-root branch at the origin
#   z=1: {0, 1-δ} = {0, 0.5}        → square-root branch at z=1
#   z=a=2: {0, 1-ε} = {0, 1-3} = {0, -2}  → ε=α+β+1-γ-δ=3, simple-then-double pole
# These half-integer exponents produce dramatic phase swirls.
const G_a, G_q, G_α, G_β, G_γ, G_δ = 2.0, 2.0, 1.5, 1.5, 0.5, 0.5

# HeunC: same demo parameters as test suite, singularities at 0, 1.
const C_q, C_α, C_γ, C_δ, C_ε = 1.0, 1.0, 2.0, 3.0, -1.0

# Domain: tighter around singularities for visual drama.
const RE_MIN, RE_MAX = -1.5, 3.5
const IM_MIN, IM_MAX = -2.0, 2.0

# Higher grid resolution for smoothness.
const NX, NY = 180, 144

const ORDER  = 30
const HSTEP  = 0.05
const NTAYL  = 60
const EPS_START = 0.05

# ----------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------
xs = range(RE_MIN, RE_MAX; length = NX)
ys = range(IM_MIN, IM_MAX; length = NY)
grid = ComplexF64[complex(x, y) for y in ys for x in xs]
@printf("Grid: %d × %d = %d targets\n", NX, NY, length(grid))

# ----------------------------------------------------------------------
# Walk: HeunG
# ----------------------------------------------------------------------
@printf("Walking HeunG...")
t0 = time()
z0 = ComplexF64(EPS_START)
u0_g, up0_g, _ = _heun_g_frobenius_at_0(
    complex(G_a), complex(G_q), complex(G_α), complex(G_β),
    complex(G_γ), complex(G_δ), z0, NTAYL)
f_g = _heun_g_rhs(complex(G_a), complex(G_q), complex(G_α), complex(G_β),
                  complex(G_γ), complex(G_δ))
prob_g = PadeTaylorProblem(f_g, (u0_g, up0_g),
                           (z0, ComplexF64(RE_MAX + 0.0im));
                           order = ORDER)
sol_g = path_network_solve(prob_g, grid;
                            h = HSTEP, order = ORDER,
                            enforce_real_axis_symmetry = true,
                            extrapolate = true,
                            max_steps_per_target = 2000)
@printf(" %.1fs (%d visited)\n", time()-t0, length(sol_g.visited_z))

G_vals = reshape(sol_g.grid_u, NX, NY)

# Walk: HeunC
@printf("Walking HeunC...")
t0 = time()
u0_c, up0_c, _ = _heun_c_frobenius_at_0(
    complex(C_q), complex(C_α), complex(C_γ), complex(C_δ),
    complex(C_ε), z0, NTAYL)
f_c = _heun_c_rhs(complex(C_q), complex(C_α), complex(C_γ),
                  complex(C_δ), complex(C_ε))
prob_c = PadeTaylorProblem(f_c, (u0_c, up0_c),
                           (z0, ComplexF64(RE_MAX + 0.0im));
                           order = ORDER)
sol_c = path_network_solve(prob_c, grid;
                            h = HSTEP, order = ORDER,
                            enforce_real_axis_symmetry = true,
                            extrapolate = true,
                            max_steps_per_target = 2000)
@printf(" %.1fs (%d visited)\n", time()-t0, length(sol_c.visited_z))

C_vals = reshape(sol_c.grid_u, NX, NY)

# ----------------------------------------------------------------------
# Domain colouring: HSV with hue=arg, magnitude-modulated brightness.
# Adds gentle log-spaced "shaded contour" stripes for magnitude legibility.
# ----------------------------------------------------------------------
# Phase portrait with subtle log-magnitude contour shading.
#   - hue       = arg f   (full colour wheel; main signal)
#   - sat       = 1       (vivid)
#   - val       = 1 with subtle ~10% darkening on log-spaced |f| isolines
#                 and ~30% darkening as |f| → 0 (zeros visible as dark spots)
# This is the Wegert "phase portrait" style — every pixel is bright,
# all phase variation is legible, magnitude is read off the contour pattern.
function color_grid(vals)
    img = Matrix{RGBf}(undef, NX, NY)
    for j in 1:NY, i in 1:NX
        v = vals[i, j]
        if !isfinite(v)
            img[i, j] = RGBf(0.10, 0.10, 0.10)
            continue
        end
        h = mod(angle(v) / (2π) + 1.0, 1.0)
        r = abs(v)
        L = log10(r + 1e-12)
        # Log-magnitude contour stripes (every 0.5 decade): ±8% brightness ripple
        ripple = 1.0 - 0.08 * cos(2π * L * 2)
        # Zero-shade: darken near |f| = 0 so zeros stand out
        zero_shade = r < 0.1 ? (r / 0.1) * 0.7 + 0.3 : 1.0
        # High-magnitude saturation cap (don't oversaturate at poles)
        val_mag = 1.0 * ripple * zero_shade
        sat = 1.0
        img[i, j] = hsv2rgb(h, sat, clamp(val_mag, 0.0, 1.0))
    end
    return img
end

# Log-magnitude image for the companion panel (viridis).
function log_mag_grid(vals, lo = -2.0, hi = 3.0)
    out = Matrix{Float64}(undef, NX, NY)
    for j in 1:NY, i in 1:NX
        v = vals[i, j]
        out[i, j] = isfinite(v) ? clamp(log10(abs(v) + 1e-12), lo, hi) : NaN
    end
    return out
end

img_g     = color_grid(G_vals)
img_c     = color_grid(C_vals)
logmag_g  = log_mag_grid(G_vals)
logmag_c  = log_mag_grid(C_vals)

# ----------------------------------------------------------------------
# Render: 2×2 panels (phase portrait + log magnitude per function) + wheel
# ----------------------------------------------------------------------
fig = Figure(size = (1600, 1100), fontsize = 14)
Label(fig[0, 1:4], "Heun functions on the complex plane — PadeTaylor.jl";
      fontsize = 24, padding = (0, 0, 10, 0))

# Build a colour-wheel legend image (a unit disc in arg/abs space).
function wheel_image(R_max = 4.0, N = 128)
    xs = range(-R_max, R_max; length = N)
    ys = range(-R_max, R_max; length = N)
    img = Matrix{RGBf}(undef, N, N)
    for j in 1:N, i in 1:N
        z = complex(xs[i], ys[j])
        r = abs(z)
        if r > R_max
            img[i, j] = RGBf(1, 1, 1)
        else
            h = mod(angle(z) / (2π) + 1.0, 1.0)
            L = log10(r + 1e-9)
            v_mag = 0.55 + 0.45 * tanh(0.7 * L)
            contour = 0.85 + 0.15 * (sin(π * L) * 0.5 + 0.5)
            val = clamp(v_mag * contour, 0.0, 1.0)
            sat = clamp(1.0 - 0.4 * exp(-r/0.05) - 0.2 * (1 - exp(-r/200)), 0.5, 1.0)
            img[i, j] = hsv2rgb(h, sat, val)
        end
    end
    return xs, ys, img
end

# Helper: draw a singularity marker (white X with black outline)
function mark_sing!(ax, z; label = nothing)
    scatter!(ax, [real(z)], [imag(z)];
             marker = :xcross, markersize = 22,
             color = :white, strokecolor = :black, strokewidth = 2)
    if label !== nothing
        text!(ax, real(z), imag(z); text = label,
              offset = (10, 10), color = :white, fontsize = 13,
              strokecolor = :black, strokewidth = 1)
    end
end

# Helper: draw a cut from z to (RE_MAX, Im(z))
function mark_cut!(ax, from_z)
    lines!(ax, [real(from_z), RE_MAX + 0.05], [imag(from_z), imag(from_z)];
           color = (:white, 0.7), linestyle = :dash, linewidth = 2)
end

# --- Row 1: HeunG phase + magnitude
ax_gp = Axis(fig[1, 1];
             title = "arg HeunG(a=2, q=2, α=1.5, β=1.5, γ=0.5, δ=0.5; z)  (Wegert phase portrait)",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
image!(ax_gp, (RE_MIN, RE_MAX), (IM_MIN, IM_MAX), img_g)
mark_sing!(ax_gp, 0; label = "0")
mark_sing!(ax_gp, 1; label = "1")
mark_sing!(ax_gp, 2; label = "a=2")
mark_cut!(ax_gp, 1); mark_cut!(ax_gp, 2)

ax_gm = Axis(fig[1, 2];
             title = "log₁₀ |HeunG|",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
hm_g = heatmap!(ax_gm, range(RE_MIN, RE_MAX; length=NX),
                       range(IM_MIN, IM_MAX; length=NY),
                       logmag_g;
                colormap = :viridis, colorrange = (-2.0, 3.0))
mark_sing!(ax_gm, 0); mark_sing!(ax_gm, 1); mark_sing!(ax_gm, 2)
mark_cut!(ax_gm, 1); mark_cut!(ax_gm, 2)
Colorbar(fig[1, 3], hm_g; label = "log₁₀ |HeunG|", width = 14)

# --- Row 2: HeunC phase + magnitude
ax_cp = Axis(fig[2, 1];
             title = "arg HeunC(q=1, α=1, γ=2, δ=3, ε=−1; z)  (Wegert phase portrait)",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
image!(ax_cp, (RE_MIN, RE_MAX), (IM_MIN, IM_MAX), img_c)
mark_sing!(ax_cp, 0; label = "0")
mark_sing!(ax_cp, 1; label = "1")
mark_cut!(ax_cp, 1)

ax_cm = Axis(fig[2, 2];
             title = "log₁₀ |HeunC|",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
hm_c = heatmap!(ax_cm, range(RE_MIN, RE_MAX; length=NX),
                       range(IM_MIN, IM_MAX; length=NY),
                       logmag_c;
                colormap = :viridis, colorrange = (-2.0, 3.0))
mark_sing!(ax_cm, 0); mark_sing!(ax_cm, 1)
mark_cut!(ax_cm, 1)
Colorbar(fig[2, 3], hm_c; label = "log₁₀ |HeunC|", width = 14)

# Row labels (function identification)
Label(fig[1, 0], "HeunG\nhalf-int exp";
      tellheight = false, rotation = π/2, fontsize = 14)
Label(fig[2, 0], "HeunC\nstandard";
      tellheight = false, rotation = π/2, fontsize = 14)

# Colour-wheel legend (right column, spanning both rows)
ax_w = Axis(fig[1:2, 4];
            title = "colour wheel\n(hue = arg f)",
            xlabel = "Re f", ylabel = "Im f",
            aspect = DataAspect(),
            limits = (-4, 4, -4, 4))
wxs, wys, wimg = wheel_image()
image!(ax_w, (-4.0, 4.0), (-4.0, 4.0), wimg)
text!(ax_w, 3.7, -0.4; text = "0", color = :black, fontsize = 11)
text!(ax_w, 0, 3.7;    text = "+π/2", color = :black, fontsize = 11, align = (:center, :center))
text!(ax_w, -3.7, -0.4; text = "±π", color = :black, fontsize = 11)
text!(ax_w, 0, -3.7;   text = "−π/2", color = :black, fontsize = 11, align = (:center, :center))

colsize!(fig.layout, 0, Relative(0.05))
colsize!(fig.layout, 1, Relative(0.38))
colsize!(fig.layout, 2, Relative(0.38))
colsize!(fig.layout, 3, Relative(0.04))
colsize!(fig.layout, 4, Relative(0.15))

Label(fig[3, 1:4],
      "× = regular singular point   ⋯ = Wolfram-convention branch cut along [singular, +∞).  Grid: $(NX)×$(NY) = $(NX*NY) targets.  Single shared path_network_solve walk per function.";
      fontsize = 11, color = :gray35,
      padding = (0, 0, 0, 4))

save(OUTPNG, fig)
@printf("Saved: %s\n", OUTPNG)
