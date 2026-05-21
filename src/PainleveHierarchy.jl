"""
    PadeTaylor.PainleveHierarchy

Problem-builder layer for the **Painlevé-I hierarchy** — the v0.2
north-star vector target alongside the Noumi–Yamada systems.  Like the
scalar `Painleve` builder and the vector `NoumiYamada` builder, this
module adds no value-function (a Painlevé-hierarchy member has infinitely
many solutions, so "which solution" stays an explicit caller choice
through the initial condition — CLAUDE.md Rule 1); what it adds is a
single discoverable RHS factory `painleve_hierarchy(:I, m)` and a
self-describing builder `PainleveHierarchyProblem(m; …)` that assembles a
`VectorProblems.VectorPadeTaylorProblem` solvable by `vector_solve_pade`.

## The Painlevé-I hierarchy

The PI hierarchy is the sequence of ODEs arising from stationary
reductions of the KdV hierarchy; the `m`-th member `P_I^(m)` is a scalar
ODE of order `2m` (pillar C §3,
`docs/v0p2_pillarC_painleve_hierarchy_findings.md:137-168`).

  - **m = 1** — PI itself, the 2nd-order ODE `u_xx = 6u^2 + x`
    (KKG 2015,
    `references/tex/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41/tritronquee_coeff.tex:124`).
  - **m = 2** — `P_I^(2)`, the **fourth-order ODE** and the primary V6
    deliverable (KKG 2015 eq. (p12), `tritronquee_coeff.tex:125-130`,
    transcribed verbatim):

        u_xxxx + 10 u_x^2 + 20 u u_xx + 40(u^3 - 6t u + 6x) = 0,

    depending on the parameter `t ∈ ℂ`.

## Companion-form vectorisation — the `x` and `t` roles

A high-order scalar ODE is integrated by the v0.2 vector pipeline only
after reduction to a first-order vector system.  The companion form
takes `y = (u, u', u'', …, u^{(2m-1)})` as the state vector; the first
`2m-1` rows are the trivial chain `y_k' = y_{k+1}` and the last row is
the ODE solved for the top derivative.  The independent variable is `x`
(the variable the integrator steps in); the parameter `t` is *captured
into the RHS closure* — it is fixed for one problem, not a state
component.  For `P_I^(2)`, isolating `u_xxxx` gives the 4-vector RHS
(pillar C §2, `findings.md:86-91`)

    y' = ( y_2,
           y_3,
           y_4,
           -10 y_2^2 - 20 y_1 y_3 - 40(y_1^3 - 6t y_1 + 6x) ),

`y = (y_1, y_2, y_3, y_4) ∈ ℂ^4`.  For `m = 1` the companion form is the
2-vector `y' = (y_2, 6 y_1^2 + x)` — the verbatim first-order lift of
`u_xx = 6u^2 + x`.  This `m = 1` member uses the *same* normalisation as
v0.1 `PainleveProblem(:I)` (`Painleve.jl:101`, `_pI_rhs() = 6u^2 + z`):
no rescaling is needed between the hierarchy convention and v0.1's PI, so
a hierarchy-`m=1` trajectory and a v0.1 PI trajectory from the same IC
coincide (test PH.1.2 — the reduction/consistency anchor).

## Scope — m ∈ {1, 2} explicit; m ≥ 3 deferred

`m = 1` and `m = 2` are implemented **explicitly** — the companion RHS is
a verbatim transcription of the source equation, the honest
representation.  The general-`m` member would require the
Lenard–Gelfand–Dikii recursion (pillar C §3, `findings.md:114-135`:
`(d/dx) L_{N+1} = (∂_x^3 + 4u ∂_x + 2u_x) L_N`, `L_0 = 1/2`) to generate
the order-`2m` ODE symbolically — that is **not** a v0.2 acceptance item.
`painleve_hierarchy(:I, m)` therefore throws for `m ≥ 3` with a clear
message; the deferred corner is registered as follow-up bead
`padetaylor-qi0` (Rule 9) with the forcing condition recorded: *needed if
a v0.2+ target requires a PI-hierarchy member of order ≥ 6 (m ≥ 3)*.

## The tritronquée asymptotic initial condition

`pI2_tritronquee_ic(x0)` returns the asymptotic seed for the `P_I^(2)`
tritronquée at a large real point `x0` (pillar C §4,
`findings.md:172-281`; KKG 2015 §7).  At `t = 0` the tritronquée obeys
`u ~ -∛6 · x^{1/3}` as `|x| → ∞` (KKG `tritronquee_coeff.tex` `V0_sector`
multline).  Seeding at large *negative* real `x0`, write `r = |x0|` and
take the real branch `Y = -∛6·r^{1/3}`; for `x < 0` the chain rule is
`d/dx = -d/dr`.

`n_terms = 1` is the leading order `u(x) = -∛6·r^{1/3}`; differentiating
gives the four components

    y_1 = -∛6 · r^{1/3}
    y_2 = +(1/3)  ∛6 · r^{-2/3}
    y_3 = +(2/9)  ∛6 · r^{-5/3}
    y_4 = +(10/27) ∛6 · r^{-8/3}.

All four after `y_1` are positive — `u` decreases as `|x|` grows, so
`u' > 0` for `x < 0`, and successive derivatives stay positive.  Pillar C
§4 (`findings.md:274`) prints `y_3` with a leading *minus*; that is a
transcription sign-slip — differentiation gives `+(2/9)`, and the test
PH.1.4 oracle is the differentiation-derived value.

`n_terms = 2` (bead V8b) adds the **first non-zero correction** of the
KKG `u = Y + Σ c_n Y^{-n}` series (KKG 2015 §7, eq. (7.2),
`findings.md:196-202`).  At `t = 0` the coefficients `c_1 = 2t`, `c_2 =
c_3 = c_4 = 0`, `c_5 = -(8/3)t^3` all vanish and `c_6 = 1`, so the series
through the first correction is `u = Y + Y^{-6}`.  With `Y = -∛6·r^{1/3}`
the correction is `Y^{-6} = (∛6)^{-6}·r^{-2} = 6^{-2}·r^{-2} = r^{-2}/36`
(the `b_1 = 1/36` of `findings.md:257`).  Differentiating
`u = -∛6·r^{1/3} + 6^{-2}·r^{-2}` three times w.r.t. `x` (`d/dx = -d/dr`,
`x < 0`) the `n_terms = 2` seed adds, on top of the `n_terms = 1` values,

    Δy_1 = + 6^{-2} · r^{-2}     Δy_2 = + 2·6^{-2} · r^{-3}
    Δy_3 = + 6·6^{-2} · r^{-4}   Δy_4 = + 24·6^{-2} · r^{-5}.

`n_terms = 2` is supported only for `t = 0` — the general `c_n(t)` series
(`c_1 = 2t`, …) is out of v0.2 scope; `t ≠ 0` with `n_terms ≥ 2` throws
(Rule 9 v1 corner — forced if a `t ≠ 0` tritronquée figure is ever
required).  `n_terms ≥ 3` throws (the `c_7…` corrections are not
implemented).  `n_terms ≤ 0` throws.

## Fail-fast contract (CLAUDE.md Rule 1)

`painleve_hierarchy` throws for an equation symbol other than `:I`, for
`m < 1`, and for `m ≥ 3`.  `PainleveHierarchyProblem` additionally
validates `length(y0) == 2m` and a non-degenerate `xspan` (the
`VectorPadeTaylorProblem` it wraps re-checks `order` and the span).
`pI2_tritronquee_ic` throws for `n_terms ∉ {1, 2}` and for `n_terms = 2`
with `t ≠ 0`.  No silent NaN / zero.

## References

  - `references/tex/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41/tritronquee_coeff.tex`
    lines 124-130 — PI (`y_xx = 6y^2 + x`) and the `P_I^(2)` eq. (p12),
    the verbatim source of both companion RHSs.
  - `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — §1 (the explicit
    `P_I^(2)` equation), §2 (the 4-vector companion form), §3 (the PI
    hierarchy / Lenard recursion), §4 (the tritronquée asymptotic IC).
  - `docs/adr/0022-vector-problem-types.md` — the governing ADR (extended
    with the `PainleveHierarchyProblem` subsection).
  - `src/VectorProblems.jl` — `VectorPadeTaylorProblem` /
    `vector_solve_pade`, the driver this builder feeds.
  - `src/NoumiYamada.jl` — the sibling vector problem-builder whose
    wrapper-struct + forwarding shape this module mirrors.
  - `src/Painleve.jl` — v0.1 `PainleveProblem(:I)`, the `m = 1`
    reduction anchor.
"""
module PainleveHierarchy

using ..VectorProblems: VectorPadeTaylorProblem

import ..VectorProblems: vector_solve_pade
import ..VectorBVP: vector_bvp_solve

export painleve_hierarchy, painleve_hierarchy_jacobian,
       PainleveHierarchyProblem, pI2_tritronquee_ic

# -----------------------------------------------------------------------------
# Companion-form RHS closures (m ∈ {1, 2}).
# -----------------------------------------------------------------------------

# m = 1 — PI, `u_xx = 6u^2 + x`, in 2-vector companion form.
# Same normalisation as v0.1 `Painleve._pI_rhs` (6u^2 + z); no rescale.
_rhs_PI1() = (x, y) -> begin
    length(y) == 2 || throw(ArgumentError(
        "painleve_hierarchy(:I, 1): state vector has length $(length(y)) " *
        "but the m=1 companion system has 2 components (u, u'). " *
        "Suggestion: pass a 2-vector state."))
    y1, y2 = y
    return [y2, 6 * y1^2 + x]
end

# m = 2 — P_I^(2), the 4th-order ODE in 4-vector companion form.
# KKG 2015 eq. (p12): u_xxxx + 10 u_x^2 + 20 u u_xx + 40(u^3 - 6t u + 6x) = 0.
# Isolating u_xxxx and stacking the trivial chain (pillar C §2).
_rhs_PI2(t) = (x, y) -> begin
    length(y) == 4 || throw(ArgumentError(
        "painleve_hierarchy(:I, 2): state vector has length $(length(y)) " *
        "but the P_I^(2) companion system has 4 components " *
        "(u, u', u'', u'''). Suggestion: pass a 4-vector state."))
    y1, y2, y3, y4 = y
    return [y2,
            y3,
            y4,
            -10 * y2^2 - 20 * y1 * y3 - 40 * (y1^3 - 6 * t * y1 + 6 * x)]
end

"""
    painleve_hierarchy(equation::Symbol, m::Integer; t = 0) -> f

Build the first-order companion-form RHS closure `f(x, y) -> Vector` for
the `m`-th member of the Painlevé-I hierarchy.  `equation` must be `:I`;
`m` selects the member (`m = 1` ⇒ PI, a 2-vector system; `m = 2` ⇒
`P_I^(2)`, a 4-vector system).  The parameter `t` is captured into the
closure — it enters `P_I^(2)` parametrically (KKG 2015 eq. (p12)) and is
unused by `m = 1` (PI carries no `t`).

The returned closure takes the independent variable `x` and the state
vector `y` and returns `y'`; `y` may carry plain numbers or `Taylor1`
coefficients (the v0.2 vector RHS contract, ADR-0020).

Throws `ArgumentError` (Rule 1) for an `equation` other than `:I`, for
`m < 1`, and for `m ≥ 3` — the general-`m` Lenard recursion is deferred
(see the module docstring; bead recorded for the follow-up).
"""
function painleve_hierarchy(equation::Symbol, m::Integer; t = 0)
    equation === :I || throw(ArgumentError(
        "painleve_hierarchy: unknown equation $(repr(equation)); only the " *
        "Painlevé-I hierarchy (:I) is implemented for v0.2. " *
        "Suggestion: pass :I."))
    m ≥ 1 || throw(ArgumentError(
        "painleve_hierarchy(:I, m): m must be ≥ 1 (got $m); the m-th PI- " *
        "hierarchy member P_I^(m) is defined for m = 1, 2, …. " *
        "Suggestion: pass m ≥ 1."))
    m == 1 && return _rhs_PI1()
    m == 2 && return _rhs_PI2(t)
    throw(ArgumentError(
        "painleve_hierarchy(:I, $m): only m ∈ {1, 2} are implemented " *
        "explicitly for v0.2.  m ≥ 3 requires the general Lenard–Gelfand– " *
        "Dikii recursion to generate the order-2m ODE symbolically — that " *
        "is a deferred follow-up bead, not a v0.2 acceptance item. " *
        "Suggestion: use m = 1 (PI) or m = 2 (P_I^(2))."))
end

# -----------------------------------------------------------------------------
# Analytic companion-form Jacobian ∂f/∂y.
# -----------------------------------------------------------------------------

"""
    painleve_hierarchy_jacobian(equation::Symbol, m::Integer; t = 0) -> Jf

Build the closure `Jf(x, y) -> d×d matrix` returning the **exact analytic
Jacobian** `∂f/∂y` of the `P_I^(m)` companion-form RHS — the hand-derived
counterpart to the `Taylor1` autodiff Jacobian `VectorBVP.vector_bvp_solve`
computes internally.  Supplying it via `jacobian = …` lets the BVP Newton
solve skip the per-node autodiff evaluations.

The companion form makes the first `2m-1` rows the trivial chain
`y_k' = y_{k+1}`; their Jacobian rows are the constant identity shift
(a single `1` on the super-diagonal).  Only the last row carries the ODE.

  - **m = 1** (`d = 2`).  `f = (y₂, 6y₁² + x)`, so

        Jf = [ 0     1 ;
              12y₁   0 ].

  - **m = 2** (`d = 4`).  Row 4 of the RHS is
    `-10y₂² - 20y₁y₃ - 40(y₁³ - 6t·y₁ + 6x)`; differentiating,
    `∂/∂y₁ = -20y₃ - 120y₁² + 240t`, `∂/∂y₂ = -20y₂`, `∂/∂y₃ = -20y₁`,
    `∂/∂y₄ = 0`, so

        Jf = [ 0                       1       0      0 ;
               0                       0       1      0 ;
               0                       0       0      1 ;
              -20y₃-120y₁²+240t      -20y₂   -20y₁    0 ].

The returned matrix follows the element type of `y` (a `Taylor1` or
`BigFloat` state yields a matrix of that element type) — it is built with
`zero(eltype(y))` / `one(eltype(y))` so it composes with the generic-`T`
Newton solve.

Throws `ArgumentError` (Rule 1) for an `equation` other than `:I` and for
`m ≥ 3` / `m < 1`, mirroring `painleve_hierarchy` (the general-`m` member
needs the deferred Lenard recursion — see the module docstring).
"""
function painleve_hierarchy_jacobian(equation::Symbol, m::Integer; t = 0)
    equation === :I || throw(ArgumentError(
        "painleve_hierarchy_jacobian: unknown equation $(repr(equation)); " *
        "only the Painlevé-I hierarchy (:I) is implemented for v0.2. " *
        "Suggestion: pass :I."))
    m ≥ 1 || throw(ArgumentError(
        "painleve_hierarchy_jacobian(:I, m): m must be ≥ 1 (got $m). " *
        "Suggestion: pass m ≥ 1."))
    if m == 1
        return (x, y) -> begin
            length(y) == 2 || throw(ArgumentError(
                "painleve_hierarchy_jacobian(:I, 1): state vector has " *
                "length $(length(y)) but the m=1 companion system has 2 " *
                "components. Suggestion: pass a 2-vector state."))
            y1 = y[1]
            z = zero(y1); o = one(y1)
            return [z  o;
                    12 * y1  z]
        end
    end
    if m == 2
        return (x, y) -> begin
            length(y) == 4 || throw(ArgumentError(
                "painleve_hierarchy_jacobian(:I, 2): state vector has " *
                "length $(length(y)) but the P_I^(2) companion system has " *
                "4 components. Suggestion: pass a 4-vector state."))
            y1, y2, y3, _ = y
            z = zero(y1); o = one(y1)
            return [z  o  z  z;
                    z  z  o  z;
                    z  z  z  o;
                    (-20*y3 - 120*y1^2 + 240*t)  (-20*y2)  (-20*y1)  z]
        end
    end
    throw(ArgumentError(
        "painleve_hierarchy_jacobian(:I, $m): only m ∈ {1, 2} are " *
        "implemented explicitly for v0.2.  m ≥ 3 requires the general " *
        "Lenard–Gelfand–Dikii recursion — a deferred follow-up bead, not " *
        "a v0.2 acceptance item.  Suggestion: use m = 1 or m = 2."))
end

# -----------------------------------------------------------------------------
# Problem builder
# -----------------------------------------------------------------------------

"""
    PainleveHierarchyProblem

Self-describing wrapper around a `VectorPadeTaylorProblem` for a
Painlevé-I hierarchy member.  Fields:

  - `m::Int`                            — the hierarchy index (`P_I^(m)`,
                                           order `2m`);
  - `t`                                 — the parameter captured into the
                                           `P_I^(2)` RHS (`nothing`-free:
                                           stored as given, unused for
                                           `m = 1`);
  - `problem::VectorPadeTaylorProblem`   — the underlying first-order
                                           vector IVP, with the
                                           companion-form RHS, built in
                                           the natural `x`-frame.

Like `PainleveProblem` and `NoumiYamadaProblem`, the wrapper carries the
family metadata `(m, t)` so the problem is self-describing, and exposes
the underlying `VectorPadeTaylorProblem` as `problem`.

Construct via `PainleveHierarchyProblem(m; y0, xspan, t, order)`.
"""
struct PainleveHierarchyProblem{P, Tt}
    m::Int
    t::Tt
    problem::P
end

"""
    PainleveHierarchyProblem(m; y0, xspan, t = 0, order = 30)
        -> PainleveHierarchyProblem

Build the `VectorPadeTaylorProblem` for the `m`-th Painlevé-I hierarchy
member in companion form.  `y0` is the `2m`-vector initial condition
`(u, u', …, u^{(2m-1)})` at `x_0 = xspan[1]`; `xspan` is the integration
window; `t` is the `P_I^(2)` parameter (unused for `m = 1`); `order` is
the Taylor truncation degree (FW 2011 §5.1 default 30).

Throws `ArgumentError` (Rule 1) for `m ∉ {1, 2}`, a `y0` whose length is
not `2m`, or a degenerate `xspan`.  Solve the result with
`vector_solve_pade(prob; h, …)` (a forwarding method is provided) or
`vector_solve_pade(prob.problem; h, …)`.
"""
function PainleveHierarchyProblem(m::Integer; y0::AbstractVector,
                                  xspan::Tuple, t = 0,
                                  order::Integer = 30)
    # painleve_hierarchy validates `m` (≥ 1, ∈ {1,2}); call it first so the
    # m ≥ 3 / m < 1 messages are the single source of truth.
    rhs = painleve_hierarchy(:I, m; t = t)
    d   = 2 * m
    length(y0) == d || throw(ArgumentError(
        "PainleveHierarchyProblem: y0 has length $(length(y0)) but the " *
        "P_I^($m) companion system has order 2m = $d (state " *
        "(u, u', …, u^{(2m-1)})). " *
        "Suggestion: pass an initial-condition vector of length 2m = $d."))
    xspan[1] == xspan[2] && throw(ArgumentError(
        "PainleveHierarchyProblem: xspan endpoints coincide " *
        "($(xspan[1]) == $(xspan[2])). " *
        "Suggestion: provide a non-degenerate integration interval."))
    prob = VectorPadeTaylorProblem(rhs, y0, xspan; order = order)
    return PainleveHierarchyProblem{typeof(prob), typeof(t)}(Int(m), t, prob)
end

"""
    vector_solve_pade(prob::PainleveHierarchyProblem; kwargs...)
        -> VectorPadeTaylorSolution

Solve the Painlevé-hierarchy problem with the v0.2 vector Padé–Taylor
integrator.  A thin forwarding method — the pattern `NoumiYamada` uses —
that unwraps to the underlying `VectorPadeTaylorProblem` and calls the
equation-agnostic `VectorProblems.vector_solve_pade`.  All keyword
arguments (`h`, `max_steps`, `step_policy`) pass through unchanged; the
hierarchy member is integrated in its natural `x`-frame, so the returned
`VectorPadeTaylorSolution` needs no coordinate remapping.
"""
vector_solve_pade(prob::PainleveHierarchyProblem; kwargs...) =
    vector_solve_pade(prob.problem; kwargs...)

"""
    vector_bvp_solve(php::PainleveHierarchyProblem, B_a, B_b, g;
                     z_a = nothing, z_b = nothing, jacobian = nothing,
                     kwargs...)
        -> VectorBVPSolution

Boundary-value-solve a Painlevé-I hierarchy member with the
Chebyshev-collocation Newton solver `VectorBVP.vector_bvp_solve`.  A thin
forwarding convenience — the global-collocation analogue of the
`vector_solve_pade(php; …)` IVP forwarding above.

The companion-form RHS and the integration segment are taken from the
wrapped `PainleveHierarchyProblem`: the RHS is `php.problem.f` (the
`painleve_hierarchy(:I, php.m; t = php.t)` closure baked in at
construction), and the segment defaults to the problem's `zspan` —
`z_a = php.problem.zspan[1]`, `z_b = php.problem.zspan[2]`.  Either
endpoint may be overridden via the `z_a` / `z_b` keywords (the v0.2
`P_I^(2)` tritronquée figure pins the BVP on a sub-segment of the IVP
window).

The boundary condition is the general linear two-point form
`B_a·y(z_a) + B_b·y(z_b) = g` (`B_a, B_b` are `2m×2m`, `g` length `2m`).

By default `jacobian` is the **exact analytic companion-form Jacobian**
`painleve_hierarchy_jacobian(:I, php.m; t = php.t)` — the hand-derived
override path, which spares the Newton solve its per-node `Taylor1`
autodiff.  Pass `jacobian = nothing` to force the autodiff path, or any
other `(z, y) -> matrix` closure.  All remaining keywords (`N`, `tol`,
`maxiter`, `initial_guess`) pass straight through to the core solver.
"""
function vector_bvp_solve(php::PainleveHierarchyProblem, B_a, B_b, g;
                          z_a = nothing, z_b = nothing,
                          jacobian = nothing, kwargs...)
    za = z_a === nothing ? php.problem.zspan[1] : z_a
    zb = z_b === nothing ? php.problem.zspan[2] : z_b
    Jf = jacobian === nothing ?
         painleve_hierarchy_jacobian(:I, php.m; t = php.t) : jacobian
    return vector_bvp_solve(php.problem.f, za, zb, B_a, B_b, g;
                            jacobian = Jf, kwargs...)
end

# -----------------------------------------------------------------------------
# Tritronquée asymptotic initial condition for P_I^(2).
# -----------------------------------------------------------------------------

"""
    pI2_tritronquee_ic(x0; t = 0, n_terms = 1) -> Vector

Asymptotic-series initial condition `(y_1, y_2, y_3, y_4) =
(u, u', u'', u''')` for the `P_I^(2)` tritronquée at a large real point
`x0` (pillar C §4; KKG 2015 §7).  At `t = 0` the tritronquée obeys
`u ~ -∛6 · x^{1/3}` as `|x| → ∞`; with `r = |x0|` and `d/dx = -d/dr` for
`x < 0`:

  - `n_terms = 1` (leading order) — `u = -∛6·r^{1/3}`, giving
    `(y_1,…,y_4) = (-∛6·r^{1/3}, (1/3)∛6·r^{-2/3}, (2/9)∛6·r^{-5/3},
    (10/27)∛6·r^{-8/3})`.
  - `n_terms = 2` (bead V8b) — adds the first non-zero KKG correction
    `Y^{-6} = 6^{-2}·r^{-2}` (`c_6 = 1`, `t = 0`; KKG eq. (7.2)).  The
    seed becomes the `n_terms = 1` values plus
    `(6^{-2}r^{-2}, 2·6^{-2}r^{-3}, 6·6^{-2}r^{-4}, 24·6^{-2}r^{-5})`.

The element type follows `x0`'s (a `BigFloat` `x0` yields a `BigFloat`
seed).  `n_terms = 2` is supported only for `t = 0`; the general
`c_n(t)` series is out of v0.2 scope, so `t ≠ 0` with `n_terms ≥ 2`
throws (Rule 9 v1 corner — forced if a `t ≠ 0` tritronquée figure is
required).  `n_terms ≥ 3` and `n_terms ≤ 0` throw (Rule 1).

Throws `ArgumentError` for `n_terms ∉ {1, 2}`, and for `n_terms = 2`
with `t ≠ 0`.
"""
function pI2_tritronquee_ic(x0; t = 0, n_terms::Integer = 1)
    n_terms in (1, 2) || throw(ArgumentError(
        "pI2_tritronquee_ic: only n_terms ∈ {1, 2} are implemented for " *
        "v0.2 — n_terms = 1 is the leading order, n_terms = 2 adds the " *
        "first non-zero KKG correction (c_6 = 1) at t = 0.  The further " *
        "c_7… corrections (KKG 2015 eq. (7.2)) are not implemented. " *
        "Suggestion: pass n_terms = 1 or 2."))
    n_terms == 2 && !iszero(t) && throw(ArgumentError(
        "pI2_tritronquee_ic: n_terms = 2 is implemented only for t = 0 " *
        "(got t = $t); at t = 0 the KKG series simplifies to u = Y + " *
        "Y^{-6}, but for t ≠ 0 the general c_n(t) series (c_1 = 2t, …) is " *
        "out of v0.2 scope.  Suggestion: pass t = 0, or n_terms = 1."))
    c  = cbrt(6 * one(float(real(typeof(x0)))))
    r  = abs(x0)
    y  = [-c * r^(1//3),
           (1//3)  * c * r^(-2//3),
           (2//9)  * c * r^(-5//3),
           (10//27) * c * r^(-8//3)]
    n_terms == 1 && return y
    # n_terms = 2: add the c_6 = 1 correction Y^{-6} = 6^{-2}·r^{-2} and
    # its three x-derivatives (d/dx = -d/dr for x < 0).
    s = one(c) / 36          # 6^{-2}
    y[1] += s * r^(-2)
    y[2] += 2  * s * r^(-3)
    y[3] += 6  * s * r^(-4)
    y[4] += 24 * s * r^(-5)
    return y
end

end # module PainleveHierarchy
