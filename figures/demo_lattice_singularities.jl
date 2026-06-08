# figures/demo_lattice_singularities.jl
#
# DEMO: plotting the *value* of a PadeTaylor solution across the complex
# plane, for two solutions whose singularities form a regular pattern —
#
#   (a) the logistic equation  u' = u(1-u)   (recast u'' = u(1-u)(1-2u)),
#       u(z) = 1/(1 + e^{-z}), whose poles form an infinite VERTICAL ROW
#       on the imaginary axis at z = iπ(2k+1), spacing 2πi; and
#
#   (b) the equianharmonic Weierstrass ℘,  u'' = 6 u^2,  u = ℘(z-1; 0, 2),
#       whose double poles form a full 2D hexagonal LATTICE.
#
# This is the path-network tier (FW 2011 §3): from one initial condition
# the solver grows a tree of safe paths across the plane, steering around
# the pole fields, and stores a local Padé approximant at every node so
# any grid point can be evaluated cheaply.  We then read |u| off the grid
# and overlay the poles that `extract_poles` reads straight back out of
# the stored Padé denominators — the package turning a computed solution
# into a literal MAP of its singularities.
#
# Run:  julia --project=figures figures/demo_lattice_singularities.jl
# Out:  figures/output/demo_lattice_singularities.png

using PadeTaylor
using CairoMakie
using Printf

const H_STEP = 0.5      # FW 2011 working step length
const ORDER  = 30       # FW 2011 Taylor order
const OUTPNG = joinpath(@__DIR__, "output", "demo_lattice_singularities.png")

# Rectangular complex lattice: reshape(grid, Nx, Ny)[i,j] == (xs[i], ys[j]).
function rect_grid(xlo, xhi, ylo, yhi, Nx, Ny)
    xs = range(xlo, xhi; length = Nx)
    ys = range(ylo, yhi; length = Ny)
    grid = ComplexF64[complex(x, y) for y in ys for x in xs]   # x inner, y outer
    return xs, ys, grid
end

# NOTE on `pole_kwargs`: the stock `extract_poles` defaults (cluster_atol=0.1,
# min_support=3, radius_t=5) OVER-SPLIT each physical pole into several duplicate
# clusters on a dense 2D path-network pole field — the per-node Padé estimates of
# one pole spread by ~0.1-0.2 (> cluster_atol) and fragment, while sparsely-seen
# poles miss the 3-node support threshold (bug padetaylor-fzse). We pass per-panel
# tuned params: a wider cluster_atol (< the inter-pole spacing) collapses the
# fragments to one rep per pole; min_support=2 / a wider radius_t recover the
# sparsely-sampled outer poles. These give ~one cross per true singularity.

# Solve `prob` over the rectangular window and return (xs, ys, |u|-matrix, poles).
function field_over_window(prob, xlo, xhi, ylo, yhi, Nx, Ny; label = "", pole_kwargs = NamedTuple())
    xs, ys, grid = rect_grid(xlo, xhi, ylo, yhi, Nx, Ny)
    @printf("  [%s] solving on %d×%d lattice over [%.1f,%.1f]×[%.1f,%.1f] ...\n",
            label, Nx, Ny, xlo, xhi, ylo, yhi); flush(stdout)
    t0 = time()
    sol = path_network_solve(prob, grid;
                             h = H_STEP,
                             max_steps_per_target = 4000,
                             enforce_real_axis_symmetry = true)
    absu = abs.(reshape(sol.grid_u, Nx, Ny))
    poles = extract_poles(sol; pole_kwargs...)
    ncov = count(isfinite, absu)
    @printf("  [%s] done in %.1f s; %d tree nodes; %d/%d cells finite; %d poles mapped.\n",
            label, time() - t0, length(sol.visited_z), ncov, Nx * Ny, length(poles));
    flush(stdout)
    return xs, ys, absu, poles
end

# One viridis |u| heatmap panel.  Overlays the ANALYTIC true poles (open black
# circles) and what `extract_poles` actually found (red ×), so the figure itself
# diagnoses the pole-finder: a red × off a circle is a spurious rep, a circle with
# no × is a missed pole (bug padetaylor-fzse).
function panel!(gp, xs, ys, absu, found, truep; title, cap)
    ax = Axis(gp; title = title, titlesize = 12,
              xlabel = "Re z", ylabel = "Im z", aspect = DataAspect(),
              limits = (first(xs), last(xs), first(ys), last(ys)))
    hm = heatmap!(ax, xs, ys, clamp.(absu, 0.0, cap);
                  colormap = :viridis, colorrange = (0.0, cap),
                  nan_color = RGBAf(1, 1, 1, 1))
    win(p) = first(xs) ≤ real(p) ≤ last(xs) && first(ys) ≤ imag(p) ≤ last(ys)
    tw = filter(win, truep); fw = filter(win, found)
    scatter!(ax, Float64.(real.(tw)), Float64.(imag.(tw));   # analytic truth
             color = :white, marker = :circle, markersize = 13, strokewidth = 1.3,
             strokecolor = :black)
    scatter!(ax, Float64.(real.(fw)), Float64.(imag.(fw));   # extract_poles
             color = (:red, 0.95), marker = :xcross, markersize = 8)
    return ax, hm, length(fw), length(tw)
end

# ----------------------------------------------------------------------
@printf("DEMO: lattice/row of singularities — value of u(z) across ℂ\n"); flush(stdout)

# (a) Logistic: vertical row of simple poles at z = iπ(2k+1).
flog(z, u, up) = u * (1 - u) * (1 - 2u)
prob_log = PadeTaylorProblem(flog, (0.5, 0.25), (0.0, 1.0); order = ORDER)
xs_l, ys_l, absu_l, poles_l =
    field_over_window(prob_log, -3.0, 3.0, -11.0, 11.0, 61, 161; label = "logistic",
                      pole_kwargs = (cluster_atol = 0.3, min_support = 2, radius_t = 8.0))

# (b) Equianharmonic Weierstrass ℘: 2D hexagonal lattice of double poles.
# IC for u = ℘(z-1; g2=0, g3=2): a double pole at z=1 (the v1 CPB.1 oracle).
fwp(z, u, up) = 6 * u^2
prob_wp = PadeTaylorProblem(fwp, (1.071822516416917, 1.710337353176786),
                            (0.0, 1.5); order = ORDER)
xs_w, ys_w, absu_w, poles_w =
    field_over_window(prob_wp, -4.0, 6.0, -5.0, 5.0, 121, 121; label = "weierstrass",
                      pole_kwargs = (cluster_atol = 0.4,))

# ----------------------------------------------------------------------
# Analytic true singularities (independent of the solver) for the overlay.
true_log = ComplexF64[im * π * (2k + 1) for k in -3:3]               # logistic row
const Ω  = 1.3630340904278908                                       # equianharmonic half-period (FW md:297)
true_wp  = ComplexF64[1 + 2Ω * (m + n/2) + im * (Ω * sqrt(3) * n)   # ℘(z-1;0,2) lattice
                      for m in -4:4 for n in -4:4]

@printf("phase: rendering 2-panel figure\n"); flush(stdout)
fig = Figure(size = (1180, 640))
Label(fig[0, 1:2], "PadeTaylor.jl — |u(z)| across ℂ.  black ○ = analytic poles,  red × = extract_poles";
      fontsize = 15, font = :bold, padding = (0, 0, 6, 0))

_, hm_l, nfl, ntl = panel!(fig[1, 1], xs_l, ys_l, absu_l, poles_l, true_log;
    title = "(a) logistic  u' = u(1-u):  vertical ROW of simple poles  (z = iπ(2k+1))",
    cap = 3.0)
Colorbar(fig[1, 1, Right()], hm_l; label = "|u|  (clamp ≤ 3)", width = 12)

_, hm_w, nfw, ntw = panel!(fig[1, 2], xs_w, ys_w, absu_w, poles_w, true_wp;
    title = "(b) Weierstrass ℘  u'' = 6u²:  2D hexagonal LATTICE of double poles",
    cap = 40.0)
Colorbar(fig[1, 2, Right()], hm_w; label = "|u|  (clamp ≤ 40)", width = 12)

@printf("  in-window: logistic found %d / true %d ; weierstrass found %d / true %d\n",
        nfl, ntl, nfw, ntw); flush(stdout)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)

# Render-sanity (Rule 1 — fail loud).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("phase: wrote %s (%d bytes)\n", OUTPNG, bytes); flush(stdout)
end
