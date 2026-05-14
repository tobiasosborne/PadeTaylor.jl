# figures/figutil.jl
#
# Shared helpers for the FW-family pole-field surface figures
# (FW 2011 Fig 3.1, 4.3, 4.4 and their relatives).  This is NOT a
# package module — it is `include`d directly by the figure scripts so
# each script stays a single self-contained piece of exposition while
# the lattice-construction and Makie-rendering boilerplate lives in
# one place.
#
# Two helpers:
#
#   * `pole_field_lattice(half, N)` builds the `N×N` square lattice
#     over `[-half, half]²` together with the flat `Vector{ComplexF64}`
#     of target points in the order that makes `reshape(_, N, N)[i, j]`
#     land on `(xs[i], ys[j])` — the orientation `surface!` expects.
#
#   * `pole_field_figure(xs, ys, Z; …)` builds the `Axis3` magnitude
#     surface (viridis fill + faint wireframe, MATLAB-`surf`-like view,
#     z compressed so the pole spikes do not flatten the field) and
#     returns `(fig, ax)` so the caller can overlay figure-specific
#     annotations before `save`-ing.

using CairoMakie

"""
    pole_field_lattice(half::Real, N::Integer) -> (xs, ys, grid)

Uniform `N×N` lattice over `[-half, half]²`.  `xs`, `ys` are the axis
ranges; `grid::Vector{ComplexF64}` is flattened so that
`reshape(grid, N, N)[i, j] == complex(xs[i], ys[j])`.
"""
function pole_field_lattice(half::Real, N::Integer)
    xs = range(-float(half), float(half); length = N)
    ys = range(-float(half), float(half); length = N)
    # y outer / x inner: matches a column-major reshape(_, N, N).
    grid = ComplexF64[complex(x, y) for y in ys for x in xs]
    return xs, ys, grid
end

"""
    pole_field_figure(xs, ys, Z; title, ucap = 50.0,
                      zlabel = "|u(z)|", size = (950, 720)) -> (fig, ax)

Render the magnitude surface `Z` (already clamped to `[0, ucap]`) as an
`Axis3` viridis surface with a faint wireframe overlay, in the
MATLAB-`surf`-like view FW uses for its pole-field figures.  Returns
`(fig, ax)`; the caller is responsible for any further overlays and
for `save`-ing.
"""
function pole_field_figure(xs, ys, Z;
                           title::AbstractString,
                           ucap::Real = 50.0,
                           zlabel::AbstractString = "|u(z)|",
                           aspect = (1.0, 1.0, 0.45),
                           size = (950, 720))
    fig = Figure(; size = size)
    ax  = Axis3(fig[1, 1];
                xlabel = "x", ylabel = "y", zlabel = zlabel, title = title,
                # Look down at ~30°, z compressed: the unbounded pole
                # spikes stay readable without flattening the xy field.
                # `aspect` defaults to a square xy footprint; pass a
                # custom tuple for a non-square domain (e.g. Fig 4.1).
                azimuth = 1.30π, elevation = 0.16π,
                aspect = aspect,
                zticks = 0:10:Int(ucap))
    surface!(ax, xs, ys, Z; colormap = :viridis, colorrange = (0.0, Float64(ucap)))
    wireframe!(ax, xs, ys, Z; color = (:black, 0.15), linewidth = 0.3)
    return fig, ax
end

"""
    pole_scatter_axis(gp, poles; title, xlim, ylim,
                      markersize = 1.8, dotcolor = :black) -> Axis

Scatter `poles::AbstractVector{<:Complex}` as a tight dot field in an
`Axis` placed at grid position `gp` (e.g. `fig[1, 1]`), square data
aspect, FW pole-location-plot framing (FW 2011 Fig. 4.7 / 4.8): black
dots over the window, a faint cross at the origin.  Poles outside the
`[xlim] × [ylim]` window are clipped.  Returns the `Axis` so the caller
can add panel-specific overlays.
"""
function pole_scatter_axis(gp, poles;
                           title::AbstractString = "",
                           xlim::Tuple = (-50.0, 50.0),
                           ylim::Tuple = (-50.0, 50.0),
                           markersize::Real = 1.8,
                           dotcolor = :black)
    ax = Axis(gp; title = title, xlabel = "x", ylabel = "y",
              titlesize = 13, aspect = DataAspect(),
              limits = (xlim[1], xlim[2], ylim[1], ylim[2]))
    inwin = filter(p -> xlim[1] ≤ real(p) ≤ xlim[2] &&
                        ylim[1] ≤ imag(p) ≤ ylim[2], poles)
    # faint origin cross, drawn under the dots
    lines!(ax, [xlim[1], xlim[2]], [0.0, 0.0]; color = (:gray, 0.35), linewidth = 0.4)
    lines!(ax, [0.0, 0.0], [ylim[1], ylim[2]]; color = (:gray, 0.35), linewidth = 0.4)
    scatter!(ax, Float64.(real.(inwin)), Float64.(imag.(inwin));
             color = dotcolor, markersize = markersize)
    return ax
end
