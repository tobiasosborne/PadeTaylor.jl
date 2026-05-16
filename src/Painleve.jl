"""
    PadeTaylor.Painleve

Per-equation problem-builder layer for the six Painlevé equations
(ADR-0006).  This module does **not** add a `painleve()` value-function
— the Painlevé transcendents have infinitely many solutions, so "which
solution" must stay an explicit caller choice (CLAUDE.md Rule 1).  What
it adds is a single discoverable constructor, `PainleveProblem`, that
assembles the correct `PadeTaylorProblem` for a chosen equation:

  - selects the right-hand side (a closure of the canonical form +
    parameters);
  - for the branch-point equations PIII / PV / PVI, applies the
    exponential coordinate transform — *wrapping* the existing
    `CoordTransforms` / `SheetTracker` factories, not duplicating them
    — so the underlying `PadeTaylorProblem` is meromorphic and the
    pole-friendly solver substrate applies unchanged;
  - bookkeeps the equation symbol + parameters so the problem is
    self-describing;
  - stores the forward/inverse coordinate maps so transformed-frame
    results can be carried back to the natural `z`-frame.

It does **not** solve anything (that stays `solve_pade` /
`path_network_solve` / `bvp_solve` / the dispatchers), it does **not**
default initial conditions (ICs are always required), and it does
**not** add automatic Riemann-sheet routing for PIII / PV / PVI — the
standing Tier-5 gap.

## The six canonical forms

  - **PI**   `u'' = 6u² + z`
    (`FW2011_painleve_methodology_JCP230.md`)
  - **PII**  `u'' = 2u³ + zu + α`
    (`FW2014_second_PII_exploration_FoCM14.md:47`, eq. 2)
  - **PIV**  `u'' = (u')²/(2u) + (3/2)u³ + 4zu² + 2(z² − α)u + β/u`
    (`ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:45`)
  - **PIII / PV / PVI** — the transformed `ζ`-plane RHS factories live
    in `src/CoordTransforms.jl` (PIII, PV; FFW 2017 §2.1) and
    `src/SheetTracker.jl` (PVI; FFW 2017 §2.2), each with the source
    equation cited verbatim in its docstring.  PVI reuses PV's
    `z = exp(ζ)` transform.

## Coordinate frames

`frame === :direct` (PI, PII, PIV): the problem is built directly in
`z`; `to_frame` / `from_frame` are the identity.  `frame ===
:transformed` (PIII, PV, PVI): the caller supplies ICs and `zspan` in
the natural `z`-frame; the constructor maps the IC into the `ζ`-frame
and builds the `PadeTaylorProblem` there.  `to_frame` /  `from_frame`
are the `(z,u,up) ↔ (ζ,w,wp)` maps.

## Forwarding methods + `PainleveSolution`

`solve_pade(pp; …)` and `path_network_solve(pp, grid; …)` accept a
`PainleveProblem` directly and return a `PainleveSolution` — the
self-describing solve-output wrapper defined in
`src/PainleveSolution.jl` (`include`d into this module; ADR-0007).
For `:direct` problems they forward transparently to `pp.problem`
(clean now that the solvers honour `prob.order`, bead
`padetaylor-9xf`).  For `:transformed` problems the two methods
differ: `path_network_solve` handles the complex `ζ`-domain, so it
maps the caller's `z`-frame grid into the `ζ`-frame, solves, and the
returned `PainleveSolution` presents `z`-frame poles / grid values;
`solve_pade` does real-axis stepping and therefore serves a
`:transformed` problem only when its `ζ`-domain is real-typed,
throwing a `path_network_solve`-pointing message otherwise.  ADR-0007
records the reasoning and supersedes ADR-0006 refinement #2 (the
earlier "`:transformed` forwarding throws" decision).

## References

  - `docs/adr/0006-painleve-problem-layer.md` — the governing ADR.
  - `src/CoordTransforms.jl`, `src/SheetTracker.jl` — the wrapped
    PIII / PV / PVI RHS factories + transforms.
  - `src/Problems.jl` — `PadeTaylorProblem`, the struct this builds.
"""
module Painleve

using ..Problems:        PadeTaylorProblem, PadeTaylorSolution
using ..PathNetwork:     PathNetworkSolution
using ..PoleField:       extract_poles
using ..CoordTransforms: pIII_transformed_rhs, pV_transformed_rhs,
                         pIII_z_to_ζ, pIII_ζ_to_z, pV_z_to_ζ, pV_ζ_to_z
using ..SheetTracker:    pVI_transformed_rhs, pVI_eta_transformed_rhs,
                         pVI_z_to_η, pVI_η_to_z

import ..Problems:    solve_pade
import ..PathNetwork: path_network_solve

export PainleveProblem, PainleveSolution
export poles, grid_values, equation, parameters, solutionname
export tritronquee, hastings_mcleod
export pii_rational, pii_airy, piv_entire

# -----------------------------------------------------------------------------
# Canonical RHS factories for the no-transform equations (PI, PII, PIV).
# PIII / PV / PVI factories are wrapped from CoordTransforms / SheetTracker.
# -----------------------------------------------------------------------------

# PI — u'' = 6u² + z.
_pI_rhs() = (z, u, up) -> 6 * u^2 + z

# PII — u'' = 2u³ + zu + α  (FW 2014 eq. 2).
_pII_rhs(α) = (z, u, up) -> 2 * u^3 + z * u + α

# PIV — u'' = (u')²/(2u) + (3/2)u³ + 4zu² + 2(z²−α)u + β/u  (RF 2014).
# The β/u and (u')²/(2u) terms are singular at u = 0; PIV solutions
# carry movable zeros, so the closure fails loud there rather than
# letting an `Inf` leak into the Taylor jet (CLAUDE.md Rule 1).
function _pIV_rhs(α, β)
    return (z, u, up) -> begin
        iszero(u) && throw(DomainError(u,
            "PainleveProblem(:IV): RHS is singular at u = 0 (the β/u and " *
            "(u')²/(2u) terms).  PIV solutions have movable zeros and the " *
            "trajectory has reached one.  Suggestion: choose ICs or a " *
            "domain whose solution does not pass exactly through u = 0."))
        return up^2 / (2 * u) + (3//2) * u^3 + 4 * z * u^2 +
               2 * (z^2 - α) * u + β / u
    end
end

# -----------------------------------------------------------------------------
# Wrapper struct
# -----------------------------------------------------------------------------

"""
    PainleveProblem

Self-describing wrapper around a `PadeTaylorProblem` for one of the six
Painlevé equations.  Fields: `equation::Symbol` (`:I`…`:VI`),
`params::NamedTuple`, `problem` (the underlying `PadeTaylorProblem`,
built in the solve-frame), `frame::Symbol` (`:direct` | `:transformed`),
the `to_frame` / `from_frame` coordinate maps `(z,u,up) ↔ (ζ,w,wp)`,
and `name::Union{Symbol,Nothing}` — a named-solution tag
(`:tritronquee`, `:hastings_mcleod`, …) set by the named-transcendent
constructors (ADR-0008), or `nothing` for a problem built from explicit
initial conditions.  `name` rides through the forwarding methods into
`PainleveSolution.name`.

Construct via `PainleveProblem(equation; kwargs...)` — see that method
— or via a named constructor such as `tritronquee` / `hastings_mcleod`.
"""
struct PainleveProblem{P, TF, FF}
    equation   :: Symbol
    params     :: NamedTuple
    problem    :: P
    frame      :: Symbol
    to_frame   :: TF
    from_frame :: FF
    name       :: Union{Symbol, Nothing}
end

# Identity (z,u,up) ↔ (z,u,up) map for the no-transform equations.
_ident3(a, b, c) = (a, b, c)

# -----------------------------------------------------------------------------
# Keyword validation — fail loud on missing / unexpected parameters.
# -----------------------------------------------------------------------------

function _validate_kw(kw::NamedTuple, eq::Symbol,
                      required::Tuple, optional::Tuple)
    allowed = (required..., optional...)
    for k in keys(kw)
        k in allowed || throw(ArgumentError(
            "PainleveProblem(:$eq): unexpected keyword `$k`.  Accepted: " *
            "$(join(string.(allowed), ", ")).  Did you pass a parameter " *
            "this equation does not take?"))
    end
    for r in required
        haskey(kw, r) || throw(ArgumentError(
            "PainleveProblem(:$eq): missing required keyword `$r`.  " *
            "Required: $(join(string.(required), ", ")).  See the " *
            "canonical form in the `Painleve` module docstring / ADR-0006."))
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Constructors
# -----------------------------------------------------------------------------

"""
    PainleveProblem(equation::Symbol; kwargs...) -> PainleveProblem

Build the `PadeTaylorProblem` for one of the six Painlevé equations.
The IC point is `zspan[1]`; ICs are always required (never defaulted).

    PainleveProblem(:I;   u0, up0, zspan, order = 30)
    PainleveProblem(:II;  α,          u0, up0, zspan, order = 30)
    PainleveProblem(:IV;  α, β,       u0, up0, zspan, order = 30)
    PainleveProblem(:III; α, β, γ, δ, u0, up0, zspan, order = 30)
    PainleveProblem(:V;   α, β, γ, δ, u0, up0, zspan, order = 30)
    PainleveProblem(:VI;  α, β, γ, δ, u0, up0, zspan, order = 30)

For `:III` / `:V` / `:VI`, `u0`, `up0`, `zspan` are given in the
natural `z`-frame; the constructor maps the IC into the `ζ`-frame.
`zspan[1]` must avoid the fixed branch point(s) (`z = 0` for III / V;
`z ∈ {0, 1}` for VI).  Throws `ArgumentError` on an unknown equation,
a missing required parameter, or a parameter the equation does not
take.
"""
function PainleveProblem(equation::Symbol; kwargs...)
    kw = NamedTuple(kwargs)
    equation === :I   && return _build_I(kw)
    equation === :II  && return _build_II(kw)
    equation === :III && return _build_III(kw)
    equation === :IV  && return _build_IV(kw)
    equation === :V   && return _build_V(kw)
    equation === :VI  && return _build_VI(kw)
    throw(ArgumentError(
        "PainleveProblem: unknown equation $(repr(equation)); expected one " *
        "of :I, :II, :III, :IV, :V, :VI."))
end

function _build_I(kw)
    _validate_kw(kw, :I, (:u0, :up0, :zspan), (:order,))
    prob = PadeTaylorProblem(_pI_rhs(), (kw.u0, kw.up0), kw.zspan;
                             order = get(kw, :order, 30))
    return PainleveProblem(:I, (;), prob, :direct, _ident3, _ident3, nothing)
end

function _build_II(kw)
    _validate_kw(kw, :II, (:α, :u0, :up0, :zspan), (:order,))
    prob = PadeTaylorProblem(_pII_rhs(kw.α), (kw.u0, kw.up0), kw.zspan;
                             order = get(kw, :order, 30))
    return PainleveProblem(:II, (; α = kw.α), prob, :direct,
                           _ident3, _ident3, nothing)
end

function _build_IV(kw)
    _validate_kw(kw, :IV, (:α, :β, :u0, :up0, :zspan), (:order,))
    prob = PadeTaylorProblem(_pIV_rhs(kw.α, kw.β), (kw.u0, kw.up0), kw.zspan;
                             order = get(kw, :order, 30))
    return PainleveProblem(:IV, (; α = kw.α, β = kw.β),
                           prob, :direct, _ident3, _ident3, nothing)
end

# --- transformed-frame equations (PIII, PV, PVI) -----------------------------
# Shared assembly: validate, branch-point guard, transform the IC + zspan
# endpoints into the ζ-frame, build the PadeTaylorProblem there.
function _build_transformed(kw, eq::Symbol, rhs, z_to_ζ, ζ_to_z,
                            branch_points;
                            frame_tag::Symbol = :transformed,
                            extra_optionals::Tuple = ())
    _validate_kw(kw, eq, (:α, :β, :γ, :δ, :u0, :up0, :zspan),
                 (:order, extra_optionals...))
    z0 = kw.zspan[1]
    any(b -> z0 == b, branch_points) && throw(ArgumentError(
        "PainleveProblem(:$eq): zspan[1] = $z0 is a fixed branch point " *
        "$(branch_points); the coordinate transform is singular there.  " *
        "Suggestion: start the problem off the branch point(s)."))
    ζ0, w0, wp0 = z_to_ζ(z0, kw.u0, kw.up0)
    ζ_end       = z_to_ζ(kw.zspan[2], kw.u0, kw.up0)[1]
    prob = PadeTaylorProblem(rhs, (w0, wp0), (ζ0, ζ_end);
                             order = get(kw, :order, 30))
    return PainleveProblem(eq, (; α = kw.α, β = kw.β, γ = kw.γ, δ = kw.δ),
                           prob, frame_tag, z_to_ζ, ζ_to_z, nothing)
end

_build_III(kw) = _build_transformed(
    kw, :III, pIII_transformed_rhs(kw.α, kw.β, kw.γ, kw.δ),
    pIII_z_to_ζ, pIII_ζ_to_z, (0,))

_build_V(kw) = _build_transformed(
    kw, :V, pV_transformed_rhs(kw.α, kw.β, kw.γ, kw.δ),
    pV_z_to_ζ, pV_ζ_to_z, (0,))

# PVI accepts a `frame` kwarg selecting the coordinate of the
# underlying `PadeTaylorProblem`.  `:transformed` (default) is the
# ζ-plane equation (FFW eq. 3, md:144); `:transformed_eta` is the
# η-plane equation (FFW eq. 5, md:154), reached by stacking a second
# exponential atop the ζ-plane to push the `z = 0` branch out of the
# finite plane.  The fixed branch points (`z ∈ {0, 1}`) are the same
# IC-guard in either frame; the difference is in the solve-frame RHS
# and the `z ↔ {ζ, η}` maps.
function _build_VI(kw)
    frame = get(kw, :frame, :transformed)
    if frame === :transformed
        return _build_transformed(
            kw, :VI, pVI_transformed_rhs(kw.α, kw.β, kw.γ, kw.δ),
            pV_z_to_ζ, pV_ζ_to_z, (0, 1);
            frame_tag = :transformed, extra_optionals = (:frame,))
    elseif frame === :transformed_eta
        return _build_transformed(
            kw, :VI, pVI_eta_transformed_rhs(kw.α, kw.β, kw.γ, kw.δ),
            pVI_z_to_η, pVI_η_to_z, (0, 1);
            frame_tag = :transformed_eta, extra_optionals = (:frame,))
    else
        throw(ArgumentError(
            "PainleveProblem(:VI): unknown frame $(repr(frame)); expected " *
            ":transformed (ζ-plane, FFW eq. 3) or :transformed_eta (η-plane, " *
            "FFW eq. 5).  See `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/" *
            "FFW2017_painleve_riemann_surfaces_preprint.md:144` (ζ) and `:154` (η)."))
    end
end

# -----------------------------------------------------------------------------
# The PainleveSolution wrapper (struct + z-frame access surface).
# Included here so the forwarding methods below can construct it; it is
# part of `module Painleve`, one conceptual namespace with the builder.
# -----------------------------------------------------------------------------

include("PainleveSolution.jl")

# -----------------------------------------------------------------------------
# Named-transcendent constructors (`tritronquee`, `hastings_mcleod`) —
# specific literature-pinned solutions, the third file of this module.
# -----------------------------------------------------------------------------

include("PainleveNamed.jl")

# -----------------------------------------------------------------------------
# Parametrised-family closed-form constructors (`pii_rational`,
# `pii_airy`, `piv_entire`) — the second kind of named Painlevé solutions
# in the ADR-0008 taxonomy.  ADR-0010 records the API decision.
# -----------------------------------------------------------------------------

include("PainleveClosedForm.jl")

# -----------------------------------------------------------------------------
# Forwarding methods — solver entry points that accept a `PainleveProblem`
# and return a self-describing `PainleveSolution` (ADR-0007).
# -----------------------------------------------------------------------------
# Both unwrap to the underlying `PadeTaylorProblem`, call the equation-
# agnostic solver, and re-wrap with the Painlevé identity + frame maps.
# `:direct` problems forward transparently.  `:transformed` problems
# differ by solver — see the module docstring's "Forwarding methods"
# section for the reasoning.

"""
    solve_pade(pp::PainleveProblem; kwargs...) -> PainleveSolution

Solve the Painlevé problem with the fixed-step Padé-Taylor integrator
and wrap the result in a `PainleveSolution`.  `:direct` problems (PI /
PII / PIV) forward transparently.  A `:transformed` problem (PIII / PV /
PVI) is served only when its `ζ`-domain is real-typed; a complex
`ζ`-domain throws — `solve_pade` does real-axis stepping — with a
suggestion to use `path_network_solve`.
"""
function solve_pade(pp::PainleveProblem; kwargs...)
    if pp.frame !== :direct && eltype(pp.problem.zspan) <: Complex
        throw(ArgumentError(
            "solve_pade(::PainleveProblem): the :$(pp.equation) problem is " *
            "solved in a transformed $(pp.frame) frame with a *complex* " *
            "domain, but solve_pade does fixed-step real-axis stepping " *
            "(`state.z < z_end` is undefined for complex z).  Suggestion: " *
            "use `path_network_solve(pp, grid)` for transformed-frame " *
            "problems — it handles the complex solve-frame domain and " *
            "returns a z-frame PainleveSolution.  (solve_pade does serve " *
            "a transformed problem whose solve-frame domain is genuinely " *
            "real, i.e. built from real-valued ICs and span.)"))
    end
    return _painleve_solution(pp, solve_pade(pp.problem; kwargs...))
end

"""
    path_network_solve(pp::PainleveProblem, grid; kwargs...) -> PainleveSolution

Solve the Painlevé problem with the FW 2011 §3.1 path-network over
`grid` and wrap the result in a `PainleveSolution`.  `grid` is given in
the natural `z`-frame: for a `:transformed` problem (PIII / PV / PVI) it
is mapped into the `ζ`-frame before solving, and the returned
solution's `poles` / `grid_values` map back to the `z`-frame.
"""
function path_network_solve(pp::PainleveProblem, grid; kwargs...)
    if pp.frame === :direct
        raw = path_network_solve(pp.problem, grid; kwargs...)
    else
        ζgrid = [_coord(pp.to_frame, z) for z in grid]
        raw   = path_network_solve(pp.problem, ζgrid; kwargs...)
    end
    return _painleve_solution(pp, raw)
end

end # module Painleve
