"""
    PadeTaylor.NoumiYamada

Problem-builder layer for the even-parity **Noumi–Yamada** `A_{2n}^{(1)}`
higher-order Painlevé systems — the v0.2 north-star target.  This module
is the vector analogue of the scalar `Painleve` builder: it does not add
a value-function (a Noumi–Yamada system has infinitely many solutions,
so "which solution" stays an explicit caller choice — CLAUDE.md Rule 1),
but a single discoverable constructor `NoumiYamadaProblem(n; …)` that
assembles a `VectorProblems.VectorPadeTaylorProblem` solvable by
`vector_solve_pade`.

## The `A_{2n}^{(1)}` system

For `l = 2n` (`n = 1, 2, …`), Noumi & Yamada 1998 define a first-order
system for the `(2n+1)` unknowns `f_0, …, f_{2n}` with `(2n+1)`
parameters `α_0, …, α_{2n}` (all indices mod `2n+1`):

    A_{2n}^{(1)}:   f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1}
                                 −  Σ_{1≤r≤n} f_{j+2r} ) + α_j

— transcribed verbatim from
`references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`,
lines 85–88, equation `\\eqref{A2n}`; cross-checked against the explicit
small-`l` Appendix formulas in `docs/v0p2_pillarB_noumi_yamada_findings.md`
§1.4.  The bracket is the *odd-offset* partial sum
`f_{j+1} + f_{j+3} + … + f_{j+2n-1}` minus the *even-offset* partial sum
`f_{j+2} + f_{j+4} + … + f_{j+2n}`.  Concretely:

  - `A_2^(1)` (`n = 1`):  `f_0' = f_0(f_1 − f_2) + α_0`
    (⇒ scalar PIV; pillar B §2);
  - `A_4^(1)` (`n = 2`):  `f_0' = f_0(f_1 − f_2 + f_3 − f_4) + α_0`
    (the v0.2 `n ≥ 4` figure / acceptance target; pillar B §3.2).

All other components are obtained from `f_0'` by the cyclic index
rotation `j → j + 1` (NY1998 main.tex:107) — `NoumiYamada` builds the
RHS for *every* component with one index-rotated comprehension, so the
cyclic structure is enforced by construction rather than transcribed.

## The `Σf_j = t` constraint as a checked invariant — not eliminated

For the `k = 1` normalisation (`Σα_j = 1`, pillar B §1.2) the dependent
variables satisfy `f_0 + f_1 + … + f_{2n} = t` (pillar B §1.2,
Matsuda 2012 lines 228–233).  This is an **invariant of the flow**, not
an externally imposed algebraic constraint: summing the system gives

    Σ_j f_j' = Σ_j f_j ( bracket_j )  +  Σ_j α_j .

Each `f_j·bracket_j` term cancels in the total sum (every product
`f_a f_b` appears once with `+` and once with `−` across the cyclic
family), so `Σ_j f_j' = Σ_j α_j = 1`, hence `d/dt (Σf_j − t) = 0`.

We therefore integrate the **full `(2n+1)`-component system directly**
— no constraint elimination, no reduced canonical-pair coordinates.
The constraint `Σf_j = t` is used two ways: (1) `NoumiYamadaProblem`
*validates* that the initial condition satisfies `Σf0 = tspan[1]` so the
invariant reads as `Σf_j = t` (rather than `Σf_j = t + const`), throwing
if it does not; (2) it is the natural numerical self-test of a solved
trajectory (test NY.1.4).  Integrating the full system keeps the RHS a
verbatim transcription of `\\eqref{A2n}` — the honest representation —
and defers the canonical-pair reduction to the V6 hierarchy layer.
ADR-0022 records this decision.

## Scope — even parity only

This module implements **only** the even-parity `A_{2n}^{(1)}` system.
It covers every v0.2 acceptance item: `A_2^(1)` (the V5b PIV-reduction
target) and `A_4^(1)` (the V8a pole-field figure and the PRD `n ≥ 4`
acceptance).  The odd-parity `A_{2n+1}^{(1)}` system (NY1998
`\\eqref{A2n+1}`, main.tex:92–102) is the *quadratic-form* variant
(⇒ PV / `A_5`); it is not needed for any v0.2 acceptance item and is
deliberately deferred — see bead `padetaylor-0ln.18` (the follow-up
registered for it), which would be forced if a v0.2+ target requires a
PV-family `A_{2n+1}` system.

## References

  - `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`
    lines 85–88 — equation `\\eqref{A2n}`, the verbatim source.
  - `docs/v0p2_pillarB_noumi_yamada_findings.md` §1 (system, both
    parities), §1.2 (the `Σf_j = t` constraint), §1.4 (small-`l`
    explicit formulas), §5.4 / §7.2 (the `A_4` Type C oracle).
  - `docs/adr/0022-vector-problem-types.md` — the governing ADR.
  - `src/VectorProblems.jl` — `VectorPadeTaylorProblem` /
    `vector_solve_pade`, the driver this builder feeds.
  - `src/Painleve.jl` — the scalar per-equation builder pattern mirrored
    here for vector systems.
"""
module NoumiYamada

using ..VectorProblems: VectorPadeTaylorProblem

import ..VectorProblems: vector_solve_pade

export noumi_yamada_rhs, NoumiYamadaProblem

# -----------------------------------------------------------------------------
# RHS factory
# -----------------------------------------------------------------------------

"""
    noumi_yamada_rhs(n::Integer; α::AbstractVector) -> f

Build the right-hand-side closure for the even-parity Noumi–Yamada
system `A_{2n}^{(1)}`.  Returns `f(t, fvec) -> Vector` of length
`2n+1`, where `fvec = (f_0, …, f_{2n})` is the dependent-variable
vector and the returned vector is `(f_0', …, f_{2n}')` given by

    f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1}  −  Σ_{1≤r≤n} f_{j+2r} ) + α_j ,

with all indices wrapped mod `2n+1` (NY1998 `\\eqref{A2n}`,
main.tex:85–88).  `α` is the `(2n+1)`-vector of parameters; it must sum
to `1` (the `k = 1` normalisation — pillar B §1.2).

The closure's contract matches the v0.2 vector RHS contract (ADR-0020):
`fvec` may carry plain numbers or `Taylor1` coefficients, and the RHS is
ordinary array + scalar arithmetic, so it composes with the vector
Taylor-jet generator unchanged.

Throws `ArgumentError` (Rule 1) for `n < 1`, an `α` whose length is not
`2n+1`, or an `α` that does not sum to `1`.
"""
function noumi_yamada_rhs(n::Integer; α::AbstractVector)
    n ≥ 1 || throw(ArgumentError(
        "noumi_yamada_rhs: n must be ≥ 1 (got $n); A_{2n}^{(1)} is defined " *
        "for n = 1, 2, … (n = 1 ⇒ A_2, n = 2 ⇒ A_4). " *
        "Suggestion: pass n ≥ 1."))
    d = 2n + 1
    length(α) == d || throw(ArgumentError(
        "noumi_yamada_rhs: α has length $(length(α)) but A_{2n}^{(1)} with " *
        "n = $n has $d parameters (α_0, …, α_$(d - 1)). " *
        "Suggestion: pass an α of length 2n+1 = $d."))
    _check_alpha_sum(α)
    # Capture α and n in the closure.  Index helper: a 0-based component
    # index j maps to the 1-based Julia slot mod(j, d) + 1.  The bracket
    # for component j is Σ_{r=1}^{n} f_{j+2r-1} − Σ_{r=1}^{n} f_{j+2r}.
    return function (t, fvec)
        length(fvec) == d || throw(ArgumentError(
            "noumi_yamada_rhs(A_$(2n)): state vector has length " *
            "$(length(fvec)) but the system has $d components. " *
            "Suggestion: pass a state vector of length 2n+1 = $d."))
        slot(j) = mod(j, d) + 1
        return [begin
            bracket = sum(fvec[slot(j + 2r - 1)] - fvec[slot(j + 2r)]
                          for r in 1:n)
            fvec[slot(j)] * bracket + α[slot(j)]
        end for j in 0:(d - 1)]
    end
end

"""
    _check_alpha_sum(α)

Throw `ArgumentError` unless `Σα_j ≈ 1` — the `k = 1` normalisation
(pillar B §1.2; NY1998 sets `α_0 + … + α_l = k`, and the project uses
`k = 1`).  The tolerance scales with the working precision so a
`BigFloat` parameter set is held to a tighter standard than a `Float64`
one.
"""
function _check_alpha_sum(α)
    s = sum(α)
    tol = 16 * eps(float(real(eltype(α))))
    abs(s - 1) ≤ tol || throw(ArgumentError(
        "Noumi–Yamada: the parameters must satisfy Σα_j = 1 (the k = 1 " *
        "normalisation; pillar B §1.2), but Σα_j = $s. " *
        "Suggestion: rescale α so its components sum to exactly 1."))
    return nothing
end

# -----------------------------------------------------------------------------
# Problem builder
# -----------------------------------------------------------------------------

"""
    NoumiYamadaProblem

Self-describing wrapper around a `VectorPadeTaylorProblem` for an
even-parity Noumi–Yamada `A_{2n}^{(1)}` system.  Fields:

  - `n::Int`                     — the system index (`A_{2n}^{(1)}`);
  - `α::Vector`                  — the `(2n+1)` parameters;
  - `problem::VectorPadeTaylorProblem` — the underlying vector IVP, with
    the `A_{2n}^{(1)}` RHS, built in the natural `t`-frame.

Like the scalar `PainleveProblem`, this wrapper carries the system
metadata `(n, α)` so the problem is self-describing, and exposes the
underlying `VectorPadeTaylorProblem` as `problem` — `vector_solve_pade`
is called on `problem` directly (a `NoumiYamada` system is integrated in
its natural `t`-frame, so no coordinate transform / forwarding shim is
needed, unlike the branch-point Painlevé equations).

Construct via `NoumiYamadaProblem(n; α, f0, tspan, order)`.
"""
struct NoumiYamadaProblem{P, A}
    n::Int
    α::A
    problem::P
end

"""
    NoumiYamadaProblem(n; α, f0, tspan, order = 30) -> NoumiYamadaProblem

Build the `VectorPadeTaylorProblem` for the even-parity Noumi–Yamada
system `A_{2n}^{(1)}`.  `α` is the `(2n+1)`-vector of parameters
(must sum to `1`, the `k = 1` normalisation); `f0` is the
`(2n+1)`-vector initial condition `(f_0(t_0), …, f_{2n}(t_0))` at
`t_0 = tspan[1]`; `tspan` is the integration window; `order` is the
Taylor truncation degree (FW 2011 §5.1 default 30).

The IC must satisfy the `Σf_j = t` constraint at the start point —
`Σf0 = tspan[1]` — so the conserved quantity `Σf_j − t` reads as
`Σf_j = t` along the whole trajectory (the constraint is an invariant
of the flow; see the module docstring).  The constructor throws if it
does not, rather than silently integrating a `Σf_j = t + const`
trajectory.

Throws `ArgumentError` (Rule 1) for `n < 1`, an `α` or `f0` whose length
is not `2n+1`, an `α` not summing to `1`, or an `f0` violating
`Σf0 = tspan[1]`.

Solve the result with `vector_solve_pade(prob.problem; h, …)`.
"""
function NoumiYamadaProblem(n::Integer; α::AbstractVector,
                            f0::AbstractVector, tspan::Tuple,
                            order::Integer = 30)
    n ≥ 1 || throw(ArgumentError(
        "NoumiYamadaProblem: n must be ≥ 1 (got $n); A_{2n}^{(1)} is " *
        "defined for n = 1, 2, …. Suggestion: pass n ≥ 1."))
    d = 2n + 1
    length(α) == d || throw(ArgumentError(
        "NoumiYamadaProblem: α has length $(length(α)) but A_{2n}^{(1)} " *
        "with n = $n has $d parameters. " *
        "Suggestion: pass an α of length 2n+1 = $d."))
    length(f0) == d || throw(ArgumentError(
        "NoumiYamadaProblem: f0 has length $(length(f0)) but A_{2n}^{(1)} " *
        "with n = $n has $d components. " *
        "Suggestion: pass an initial-condition vector of length 2n+1 = $d."))
    _check_alpha_sum(α)
    # Constraint check: Σf0 = t0 so the invariant reads Σf_j = t.
    t0 = tspan[1]
    s0 = sum(f0)
    tol = sqrt(eps(float(real(promote_type(eltype(f0), typeof(t0))))))
    abs(s0 - t0) ≤ tol * max(1, abs(t0)) || throw(ArgumentError(
        "NoumiYamadaProblem: the initial condition must satisfy the " *
        "constraint Σf0 = tspan[1] (= $t0) so the flow invariant reads " *
        "Σf_j = t, but Σf0 = $s0. " *
        "Suggestion: choose an f0 whose components sum to tspan[1] " *
        "(e.g. the A_$(2n) Type C point f_j = t0/$d)."))
    rhs  = noumi_yamada_rhs(n; α = α)
    prob = VectorPadeTaylorProblem(rhs, f0, tspan; order = order)
    return NoumiYamadaProblem{typeof(prob), typeof(α)}(Int(n), α, prob)
end

# -----------------------------------------------------------------------------
# Forwarding method — solve a NoumiYamadaProblem directly.
# -----------------------------------------------------------------------------

"""
    vector_solve_pade(prob::NoumiYamadaProblem; kwargs...)
        -> VectorPadeTaylorSolution

Solve the Noumi–Yamada problem with the v0.2 vector Padé–Taylor
integrator.  This is a thin forwarding method — exactly the pattern the
scalar `Painleve` builder uses for `solve_pade` — that unwraps to the
underlying `VectorPadeTaylorProblem` and calls the equation-agnostic
`VectorProblems.vector_solve_pade`.  All keyword arguments (`h`,
`max_steps`, `step_policy`) pass through unchanged.  A `NoumiYamada`
system is integrated in its natural `t`-frame, so the returned
`VectorPadeTaylorSolution` needs no coordinate-frame remapping; it is
returned as-is.
"""
vector_solve_pade(prob::NoumiYamadaProblem; kwargs...) =
    vector_solve_pade(prob.problem; kwargs...)

end # module NoumiYamada
