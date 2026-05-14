"""
    PadeTaylorMakieExt

Package-extension providing a Makie plot recipe for `PainleveSolution`,
per ADR-0003 + ADR-0007.  Loaded automatically when both `PadeTaylor`
and `Makie` (e.g. via `CairoMakie` / `GLMakie`) are present:

```julia
using PadeTaylor, CairoMakie

pp  = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0, zspan = (0.0, 3.0))
sol = path_network_solve(pp, grid)
fig = painleveplot(sol)          # complex-z-plane view + poles
save("pii.png", fig)
```

## Design (ADR-0003 §"Translation only — no algorithmic logic")

The extension is *presentation-only*.  It reads a `PainleveSolution`
through the public surface the wrapper already exposes —
`grid_values(sol)` for the `z`-frame point set and `poles(sol)` for the
`z`-frame pole locations (ADR-0007) — and renders them with Makie
primitives.  No solving, no pole-finding, no frame arithmetic lives
here: that is all upstream in `PadeTaylor` proper.  Consequently the
extension is automatically correct for `:transformed` problems too —
`grid_values` / `poles` have already mapped everything back to the
natural `z`-frame.

## The view

`painleveplot` draws one complex-`z`-plane `Axis` (`Re z` vs `Im z`) —
the signature Painlevé visualisation, the same framing FW 2011 uses for
its pole-field figures.  The raw-type dispatch:

  - a `PadeTaylorSolution` raw (a `solve_pade` trajectory) is drawn as
    a path with its initial condition marked;
  - a `PathNetworkSolution` raw (a `path_network_solve` grid) is drawn
    as a scatter of evaluation points coloured by `|u|`.

Extracted poles are overlaid as red diamonds (suppressible via
`show_poles = false`); `pole_kwargs` is forwarded to `poles`.

Reference: docs/adr/0007-painleve-solution-wrapper.md,
docs/adr/0003-extensions-pattern.md.
"""
module PadeTaylorMakieExt

using PadeTaylor: PainleveSolution, PadeTaylorSolution, PathNetworkSolution,
                  grid_values, poles, equation, parameters, solutionname
import PadeTaylor: painleveplot
using Makie: Makie, Figure, Axis, lines!, scatter!, axislegend, Point2f

const _EQUATION_NAMES = (I = "Painlevé I",   II = "Painlevé II",
                         III = "Painlevé III", IV = "Painlevé IV",
                         V = "Painlevé V",   VI = "Painlevé VI")

# Provenance title: equation name + parameters + named-solution tag.
function _title(sol::PainleveSolution)
    base = get(_EQUATION_NAMES, equation(sol), string(equation(sol)))
    p = parameters(sol)
    isempty(p) ||
        (base *= "  (" * join(("$k = $v" for (k, v) in pairs(p)), ", ") * ")")
    solutionname(sol) === nothing ||
        (base *= " — " * string(solutionname(sol)))
    return base
end

# Complex vector → Makie 2-D points in the (Re, Im) plane.
_points(zs) = Point2f.(real.(zs), imag.(zs))

# Raw-type dispatch: a callable trajectory is a path; a path-network is
# a coloured scatter of independent evaluations.
function _draw_raw!(ax, ::PadeTaylorSolution, gz, gu)
    lines!(ax, _points(gz); color = :steelblue, linewidth = 2,
           label = "trajectory")
    scatter!(ax, _points(gz[1:1]); color = :seagreen, markersize = 11,
             label = "initial condition")
    return nothing
end

function _draw_raw!(ax, ::PathNetworkSolution, gz, gu)
    scatter!(ax, _points(gz); color = abs.(gu), colormap = :viridis,
             markersize = 7, label = "grid |u|")
    return nothing
end

"""
    painleveplot(sol::PainleveSolution; size = (820, 620),
                 show_poles = true, pole_kwargs = (;)) -> Makie.Figure

See the `PadeTaylor.painleveplot` docstring.  `size` is the figure
size; `show_poles` toggles the pole overlay; `pole_kwargs` is forwarded
to `poles(sol; pole_kwargs...)`.
"""
function painleveplot(sol::PainleveSolution;
                      size = (820, 620),
                      show_poles::Bool = true,
                      pole_kwargs = (;))
    fig = Figure(; size = size)
    ax  = Axis(fig[1, 1]; title = _title(sol),
               xlabel = "Re z", ylabel = "Im z")

    gz, gu, _ = grid_values(sol)
    _draw_raw!(ax, sol.raw, gz, gu)

    if show_poles
        pls = poles(sol; pole_kwargs...)
        isempty(pls) || scatter!(ax, _points(pls);
                                 color = :red, marker = :diamond,
                                 markersize = 12, label = "poles")
    end

    axislegend(ax; position = :rb, framevisible = true)
    return fig
end

end # module PadeTaylorMakieExt
