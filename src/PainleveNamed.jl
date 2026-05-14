# src/PainleveNamed.jl — included into `module Painleve` (ADR-0008).
#
# ## Why this file exists
#
# `PainleveProblem` (src/Painleve.jl) builds the problem for a chosen
# equation *given* its initial conditions.  But the interesting Painlevé
# transcendents — the tritronquée, Hastings–McLeod, … — are *specific*
# solutions, singled out from the equation's infinitely many by initial
# or boundary conditions that the literature has computed once, to high
# precision.  A user should not have to know those numbers.  This file
# is the discoverable home for them: a named constructor per such
# solution, returning the `PainleveProblem` with the cited IC and a
# `name` tag baked in.
#
# ## What v1 ships, and why only this
#
# A survey of every in-tree reference (worklog 031) sorted the named
# Painlevé solutions into three kinds.  Only one kind is cleanly
# shippable as a constructor:
#
#   - **point-IC named solutions** — the paper gives `(u, u')` at a
#     point to high precision.  Exactly two qualify with 16-digit data:
#     the **PI tritronquée** (FW 2011 §4.1) and the **PII
#     Hastings–McLeod** solution at `α = 0` (FW 2014).  Those are the
#     two this file ships.
#   - **asymptotic-BC named solutions** (general-`α` Hastings–McLeod,
#     the secondary HM, Ablowitz–Segur, …) — defined by `z → ±∞`
#     behaviour; turning that into a point IC is a connection problem,
#     not a constructor.  Out of scope (ADR-0008).
#   - **closed-form parametrised families** (PII rational / Airy, PIV
#     entire) — shippable, but a different API shape with exact
#     oracles; deferred to bead `padetaylor-icf`.
#
# ## How the tag travels
#
# Each constructor builds an ordinary `PainleveProblem` via the keyword
# constructor, then stamps it with `name` (`:tritronquee` /
# `:hastings_mcleod`) through `_with_name`.  `name` rides the
# `PainleveProblem` through the forwarding methods into
# `PainleveSolution.name` (ADR-0007 + ADR-0008), so a solved named
# transcendent reports its own provenance:
# `solutionname(solve_pade(tritronquee(:I)))  === :tritronquee`.
#
# Reference: docs/adr/0008-named-transcendent-constructors.md.

# A copy of `pp` with its named-solution tag set.  Used only by the
# named constructors below — every other code path leaves `name`
# `nothing`.
_with_name(pp::PainleveProblem, name::Symbol) =
    PainleveProblem(pp.equation, pp.params, pp.problem, pp.frame,
                    pp.to_frame, pp.from_frame, name)

# Shared zspan guard: the cited `(u, u')` are the solution's value at a
# *specific* point; placing them at a different `zspan[1]` would
# silently pose a different IVP (CLAUDE.md Rule 1).
function _require_ic_point(fn::AbstractString, zspan, ic_point)
    zspan[1] == ic_point || throw(ArgumentError(
        "$fn: the cited initial condition is defined at z = $ic_point; " *
        "zspan[1] must equal $ic_point (got $(zspan[1])).  Suggestion: " *
        "pass `zspan = ($ic_point, z_end)`."))
    return nothing
end

"""
    tritronquee(equation::Symbol = :I; zspan = (0.0, 10.0), order = 30)
        -> PainleveProblem

The Painlevé-I **tritronquée** solution — the solution of `u'' = 6u² +
z` with one pole-free sector, the one the FW 2011 figure programme is
built around.

Returns a `PainleveProblem` carrying the FW 2011 §4.1 initial
condition

    u(0)  ≈ -0.1875543083404949
    u'(0) ≈  0.3049055602612289

(`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:224-229` — computed by Fornberg &
Weideman with a 32-digit-precision Maple BVP solver over `[-20i, 20i]`,
accurate to better than `10⁻²⁰`) and `name = :tritronquee`, which
propagates into `PainleveSolution.name` when solved.

`zspan` is the integration window (`zspan[1]` must be `0.0` — the point
the IC is defined at); `order` is the Taylor truncation degree.

Throws an `ArgumentError` for any `equation` other than `:I`: only PI
has an in-tree tritronquée initial condition.  (The PII / PIV
"tronquée" solutions in the references are defined by asymptotic
boundary conditions and are out of scope — ADR-0008.)
"""
function tritronquee(equation::Symbol = :I;
                     zspan = (0.0, 10.0), order::Integer = 30)
    equation === :I || throw(ArgumentError(
        "tritronquee: only :I has an in-tree tritronquée initial condition " *
        "(FW 2011 §4.1).  The PII / PIV \"tronquée\" solutions in the " *
        "references are asymptotic-BC-defined and out of scope (ADR-0008).  " *
        "Got equation = $(repr(equation))."))
    _require_ic_point("tritronquee", zspan, 0.0)
    # FW 2011 §4.1 eq. (4.1) — verbatim.
    u0  = -0.1875543083404949
    up0 =  0.3049055602612289
    pp  = PainleveProblem(:I; u0 = u0, up0 = up0, zspan = zspan, order = order)
    return _with_name(pp, :tritronquee)
end

"""
    hastings_mcleod(; branch = :positive, zspan = (0.0, 10.0), order = 30)
        -> PainleveProblem

The **Hastings–McLeod** solution of Painlevé II at `α = 0` (`u'' = 2u³
+ zu`) — pole-free and non-oscillatory on the whole real axis; the PII
solution behind the Tracy–Widom distribution.

Returns a `PainleveProblem` carrying the FW 2014 initial condition

    (u(0), u'(0)) ≈ (±0.3670615515480784, ∓0.2953721054475501)

(`references/markdown/FW2014_second_PII_exploration_FoCM14/
FW2014_second_PII_exploration_FoCM14.md:252-258`) and `name =
:hastings_mcleod`.  The `α = 0` Hastings–McLeod solution has two
sign-symmetric copies that "only differ in sign" (the PII `α = 0`
symmetry `u → -u`); `branch` selects between them:

  - `:positive` → `u(0) > 0` (the `(+, -)` sign pattern above);
  - `:negative` → `u(0) < 0` (the `(-, +)` pattern).

`zspan` is the integration window (`zspan[1]` must be `0.0`); `order`
is the Taylor truncation degree.  Throws an `ArgumentError` for a
`branch` other than `:positive` / `:negative`.
"""
function hastings_mcleod(; branch::Symbol = :positive,
                         zspan = (0.0, 10.0), order::Integer = 30)
    branch in (:positive, :negative) || throw(ArgumentError(
        "hastings_mcleod: branch must be :positive or :negative (got " *
        "$(repr(branch))).  The PII α = 0 Hastings–McLeod solution has two " *
        "sign-symmetric copies; :positive selects u(0) > 0, :negative " *
        "selects u(0) < 0."))
    _require_ic_point("hastings_mcleod", zspan, 0.0)
    # FW 2014 — verbatim; the (±, ∓) sign pattern.
    s   = branch === :positive ? 1.0 : -1.0
    u0  =  s * 0.3670615515480784
    up0 = -s * 0.2953721054475501
    pp  = PainleveProblem(:II; α = 0.0, u0 = u0, up0 = up0,
                          zspan = zspan, order = order)
    return _with_name(pp, :hastings_mcleod)
end
