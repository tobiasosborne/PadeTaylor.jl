# src/PainleveSolution.jl — included into `module Painleve` (ADR-0007).
#
# ## Why this file exists
#
# `PainleveProblem` (src/Painleve.jl, ADR-0006) unified problem
# *construction* for the six Painlevé equations.  This file is the
# output-side mirror: `PainleveSolution`, a thin self-describing wrapper
# around whichever of the six raw solve-output types the solver layer
# produced.
#
# The solver substrate is equation-agnostic by design (ADR-0001): a
# `PadeTaylorSolution` returned from `solve_pade` carries no record that
# it *is* a Painlevé-I solution, with which parameters, on which branch;
# a `PathNetworkSolution` is a bare grid scatter; none of the six raw
# types pretty-prints; pole extraction exists for one of them.
# `PainleveSolution` adds exactly the missing layer — provenance and a
# uniform `z`-frame access surface — without touching the substrate.
#
# ## What it carries
#
# The struct stores the Painlevé identity (`equation`, `params`, the
# named-solution tag `name`) and the coordinate frame (`frame` plus the
# `to_frame` / `from_frame` triple maps), copied verbatim from the
# `PainleveProblem` the solution came from, alongside the raw output
# `raw`.  `name` is `nothing` for a problem built from explicit ICs; it
# is reserved for the deferred named-transcendent constructors
# (`tritronquee`, `hastings_mcleod`, …, ADR-0006 "Deferred"), which will
# set it — storing the field now keeps that thread purely additive.
#
# ## The `z`-frame contract
#
# Everything `PainleveSolution` exposes — the callable `sol(z)`,
# `poles(sol)`, `grid_values(sol)` — speaks the *natural* `z`-frame.
# For a `:direct` problem (PI / PII / PIV) that is the solve frame and
# the maps are the identity.  For a `:transformed` problem (PIII / PV /
# PVI) the underlying `PadeTaylorProblem` was built in a `ζ`-frame; the
# wrapper does the `z ↔ ζ` mapping at every boundary so the caller never
# sees `ζ`.  The one honest leak: `sol.raw` for a `:transformed`
# solution holds `ζ`-frame internals (its Padé store, step lengths, and
# `grid_*` arrays are all `ζ`-frame).  `grid_values(sol)` is the
# `z`-frame view of those point values; the Padé store has no `z`-frame
# image and stays `ζ`-frame (ADR-0007, superseding ADR-0006 refinement
# #2).
#
# ## Why a wrapper, not a common representation
#
# `PainleveSolution` unifies *provenance and access*, not *internal
# representation*.  A `PadeTaylorSolution` is a dense callable
# trajectory; a `PathNetworkSolution` is a grid of independent
# evaluations.  Forcing a common internal form on the two would lose
# information either way, so the wrapper keeps `raw` as-is and dispatches
# behaviour on its type — `sol(z)` works for the dense trajectory and
# fails loud for the grid (CLAUDE.md Rule 1: a grid has no dense
# interpolant, and returning a nearest-grid value would be a
# truthful-looking lie).
#
# Reference: docs/adr/0007-painleve-solution-wrapper.md.

# -----------------------------------------------------------------------------
# The wrapper struct
# -----------------------------------------------------------------------------

"""
    PainleveSolution{S, TF, FF}

Self-describing wrapper around a raw solve-output `raw::S` for one of
the six Painlevé equations.  Fields: `equation::Symbol` (`:I`…`:VI`),
`params::NamedTuple`, `name::Union{Symbol,Nothing}` (a named-solution
tag such as `:tritronquee`, or `nothing` for a generic IC),
`frame::Symbol` (`:direct` | `:transformed`), the `to_frame` /
`from_frame` coordinate maps `(z,u,up) ↔ (ζ,w,wp)`, and `raw` (the
underlying `PadeTaylorSolution` / `PathNetworkSolution` / …).

Produced by `solve_pade(::PainleveProblem)` and
`path_network_solve(::PainleveProblem, grid)`; see those methods.  The
supported access surface — `sol(z)`, `poles(sol)`, `grid_values(sol)` —
is entirely in the natural `z`-frame.  Reaching into `sol.raw` for a
`:transformed` solution exposes `ζ`-frame internals.
"""
struct PainleveSolution{S, TF, FF}
    equation   :: Symbol
    params     :: NamedTuple
    name       :: Union{Symbol, Nothing}
    frame      :: Symbol
    to_frame   :: TF
    from_frame :: FF
    raw        :: S
end

# Build from a `PainleveProblem` + the raw solve output it produced.
# `name` stays `nothing` until the named-transcendent constructors land.
_painleve_solution(pp::PainleveProblem, raw, name::Union{Symbol,Nothing} = nothing) =
    PainleveSolution(pp.equation, pp.params, name, pp.frame,
                     pp.to_frame, pp.from_frame, raw)

# Coordinate-only projection of a frame map.  The `(z,u,up) ↔ (ζ,w,wp)`
# triple maps' first component depends *only* on the coordinate — true
# for `_ident3` and for every exp-transform in `CoordTransforms` (e.g.
# `ζ = 2 log z` ignores `u`, `u'`).  So feeding dummy solution-value
# arguments recovers just the coordinate image.  This is the same idiom
# `_build_transformed` already uses for `ζ_end`.
_coord(map, c) = map(c, zero(c), zero(c))[1]

# -----------------------------------------------------------------------------
# Callable — dense evaluation in the natural z-frame
# -----------------------------------------------------------------------------

# A `PathNetworkSolution` is a grid of independent evaluations, not a
# dense interpolant: calling it fails loud rather than faking a value.
_assert_callable(::PadeTaylorSolution) = nothing
_assert_callable(raw) = throw(ArgumentError(
    "PainleveSolution: this solution wraps a $(nameof(typeof(raw))), which " *
    "is a grid of independent evaluations, not a dense interpolant — there " *
    "is no `sol(z)` callable.  Suggestion: read the grid with " *
    "`grid_values(sol)`, or re-solve with `solve_pade` for a callable " *
    "trajectory."))

"""
    (sol::PainleveSolution)(z) -> (u, u')

Evaluate the solution at `z` in the natural `z`-frame.  For a `:direct`
problem this forwards transparently to the wrapped trajectory.  For a
`:transformed` problem it maps `z → ζ`, evaluates the `ζ`-frame
trajectory, and maps the result `(w, w') → (u, u')` back.

Throws if the wrapped `raw` is not a dense callable trajectory (e.g. a
`PathNetworkSolution` — use `grid_values` for grid output).
"""
function (sol::PainleveSolution)(z)
    _assert_callable(sol.raw)
    sol.frame === :direct && return sol.raw(z)
    ζ        = _coord(sol.to_frame, z)
    w, wp    = sol.raw(ζ)
    _, u, up = sol.from_frame(ζ, w, wp)
    return (u, up)
end

# -----------------------------------------------------------------------------
# Poles — z-frame pole locations
# -----------------------------------------------------------------------------

"""
    poles(sol::PainleveSolution; kwargs...) -> Vector{<:Complex}

Pole locations of the solution in the natural `z`-frame.  Forwards to
`PoleField.extract_poles` (which reads the per-node Padé store of a
`PathNetworkSolution` or a `PadeTaylorSolution`); `kwargs` are passed
through (`radius_t`, `min_residue`, `cluster_atol`, `min_support`).
For a `:transformed` problem the extracted `ζ`-frame poles are mapped
back to the `z`-frame.

Throws — rather than returning an empty vector — if the wrapped `raw`
has no pole-extraction route, so that an empty result always means "no
poles found", never "extraction not wired".
"""
function poles(sol::PainleveSolution; kwargs...)
    raw = sol.raw
    (raw isa PathNetworkSolution || raw isa PadeTaylorSolution) || throw(ArgumentError(
        "PainleveSolution: pole extraction is not wired for a " *
        "$(nameof(typeof(raw))) raw solution.  Suggestion: extract poles " *
        "from `sol.raw` directly, or solve with `solve_pade` / " *
        "`path_network_solve`, whose Padé store `extract_poles` reads."))
    ζpoles = extract_poles(raw; kwargs...)
    sol.frame === :direct && return ζpoles
    return [_coord(sol.from_frame, ζ) for ζ in ζpoles]
end

# -----------------------------------------------------------------------------
# Grid values — z-frame point values
# -----------------------------------------------------------------------------

# Per-raw-type point-value arrays, in the raw solve frame.
_raw_grid(r::PathNetworkSolution) = (r.grid_z, r.grid_u, r.grid_up)
_raw_grid(r::PadeTaylorSolution)  =
    (r.z, getindex.(r.y, 1), getindex.(r.y, 2))
_raw_grid(r) = throw(ArgumentError(
    "PainleveSolution: grid_values is not defined for a " *
    "$(nameof(typeof(r))) raw solution.  Suggestion: read `sol.raw` " *
    "directly, or use `sol(z)` if the solution is a callable trajectory."))

"""
    grid_values(sol::PainleveSolution) -> (z, u, u')

Point values of the solution in the natural `z`-frame: parallel vectors
of evaluation points and the state there.  For a `PathNetworkSolution`
raw this is the Stage-2 grid; for a `PadeTaylorSolution` raw it is the
segment breakpoints.  For a `:transformed` problem each `(ζ, w, w')`
triple is mapped back through `from_frame` — this is the `z`-frame view
of a transformed-frame solution, whose `raw.grid_*` are `ζ`-frame.
"""
function grid_values(sol::PainleveSolution)
    z, u, up = _raw_grid(sol.raw)
    sol.frame === :direct && return (z, u, up)
    mapped = sol.from_frame.(z, u, up)
    return (getindex.(mapped, 1), getindex.(mapped, 2), getindex.(mapped, 3))
end

# -----------------------------------------------------------------------------
# Accessors + provenance display
# -----------------------------------------------------------------------------

"""
    equation(sol::PainleveSolution) -> Symbol

The Painlevé equation symbol (`:I` … `:VI`) this solution belongs to.
"""
equation(sol::PainleveSolution) = sol.equation

"""
    parameters(sol::PainleveSolution) -> NamedTuple

The equation parameters (`(;)`, `(; α)`, or `(; α, β, γ, δ)`).
"""
parameters(sol::PainleveSolution) = sol.params

"""
    solutionname(sol::PainleveSolution) -> Union{Symbol, Nothing}

The named-solution tag (`:tritronquee`, …) or `nothing` for a solution
built from explicit initial conditions.
"""
solutionname(sol::PainleveSolution) = sol.name

const _EQUATION_NAMES = (I = "Painlevé I",   II = "Painlevé II",
                         III = "Painlevé III", IV = "Painlevé IV",
                         V = "Painlevé V",   VI = "Painlevé VI")

# Cheap, total structural summary of the wrapped raw output.  Never
# throws and never computes poles — extraction can be expensive and can
# itself fail loud, so it has no place in `show`.
_raw_summary(r::PadeTaylorSolution) =
    "$(length(r.h))-segment Padé trajectory"
_raw_summary(r::PathNetworkSolution) =
    "path network — $(length(r.visited_z)) visited nodes, " *
    "$(length(r.grid_z)) grid points"
_raw_summary(r) = string(nameof(typeof(r)))

# z-frame domain endpoints, or `nothing` when the raw type has no
# well-defined span (kept total so `show` never throws).
_raw_span(r::PadeTaylorSolution)  = (first(r.z), last(r.z))
_raw_span(r::PathNetworkSolution) =
    isempty(r.visited_z) ? nothing : (first(r.visited_z), nothing)
_raw_span(r) = nothing

function Base.show(io::IO, ::MIME"text/plain", sol::PainleveSolution)
    eqname = get(_EQUATION_NAMES, sol.equation, string(sol.equation))
    println(io, "PainleveSolution — ", eqname)
    if isempty(sol.params)
        println(io, "  parameters      : (none)")
    else
        println(io, "  parameters      : ",
                join(("$k = $v" for (k, v) in pairs(sol.params)), ", "))
    end
    println(io, "  solution branch : ",
            sol.name === nothing ? "(generic IC)" : repr(sol.name))
    println(io, "  frame           : ", sol.frame)
    span = _raw_span(sol.raw)
    if span !== nothing
        a, b = span
        az = sol.frame === :direct ? a : _coord(sol.from_frame, a)
        if b === nothing
            println(io, "  domain (z)      : IC at ", az)
        else
            bz = sol.frame === :direct ? b : _coord(sol.from_frame, b)
            println(io, "  domain (z)      : ", az, " → ", bz)
        end
    end
    print(io, "  representation  : ", _raw_summary(sol.raw))
end

function Base.show(io::IO, sol::PainleveSolution)
    eqname = get(_EQUATION_NAMES, sol.equation, string(sol.equation))
    tag    = sol.name === nothing ? "" : " " * string(sol.name)
    print(io, "PainleveSolution(", eqname, tag, ", ", sol.frame, ")")
end
