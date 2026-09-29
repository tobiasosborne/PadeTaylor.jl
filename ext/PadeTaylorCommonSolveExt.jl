"""
    PadeTaylorCommonSolveExt

Package-extension adapter wiring `PadeTaylorProblem` + `PadeTaylorAlg`
into the `CommonSolve.jl` `init` / `step!` / `solve!` / `solve` interface
per ADR-0003 (Stage 1 design lock).  Loaded automatically when both
`PadeTaylor` and `CommonSolve` are present:

```julia
using PadeTaylor, CommonSolve

prob = PadeTaylorProblem(f, (u0, up0), (0.0, 1.5); order = 30)
alg  = PadeTaylorAlg(; h = 0.5)

# Convenience (CommonSolve's default `solve = solve! ∘ init`):
sol = solve(prob, alg)

# Streaming:
integ = init(prob, alg)
while !integ.done
    step!(integ)
end
sol = solve!(integ)
```

## Design (ADR-0003 §"What goes in each extension")

  - `PadeTaylorAlg` is declared in the main `PadeTaylor` module so it
    is reachable without a qualified name; this extension adds the
    `CommonSolve` methods on `(PadeTaylorProblem, PadeTaylorAlg)`.
  - `init(prob, alg; check_in_class = true)` — constructs a
    `PadeTaylorIntegrator` wrapping the inner `PadeStepperState` +
    accumulators.  Pre-pushes the IC onto the trajectory.  The
    `check_in_class` keyword is the SAME out-of-class guard `solve_pade`
    applies by default (bug `padetaylor-v1ub`, ADR-0033); CommonSolve's
    default `solve(prob, alg; kw...) = solve!(init(prob, alg; kw...))`
    forwards it, so `solve(prob, alg; check_in_class = false)` disables it.
  - `step!(integ)` — takes one step of signed length
    `dir · min(alg.h, |z_end − state.z|)` with `dir = sign(z_end − z_start)`,
    exactly as `solve_pade` does (bug `padetaylor-xhjw`), so a DESCENDING
    span is integrated leftward and the final step lands exactly on
    `z_end`.  With the guard on it calls `pade_step_with_defect!` +
    `check_in_class!`; with it off, the unchecked `pade_step_with_pade!`.
    Sets `integ.done = true` once `dir · (z_end − state.z) ≤ 0`.

Before bug `padetaylor-ncfa` this layer compared `state.z ≥ z_end`, so a
descending span such as `(0.0, -1.5)` produced an integrator born done and
`solve` silently returned the one-node initial condition while `solve_pade`
returned four nodes; and it always used the unchecked stepper.  Both
divergences are gone and are pinned by exact-equality tests in
`test/ext_commonsolve_test.jl`.  A COMPLEX span is not supported by either
path in v1: `solve_pade` and this layer both throw the same `MethodError`
from the real-ordered loop condition (`isless(::Complex, ::Complex)`);
that shared limitation is also pinned there.
  - `solve!(integ)` — drives `step!` in a loop until `integ.done`;
    returns `PadeTaylorSolution{T, Y, P}` assembled from the
    accumulators.
  - `solve(prob, alg)` — `CommonSolve` provides a default
    `solve = solve! ∘ init`; we don't override.

This is a **translation layer**.  The trajectory bytes from
`solve(prob, alg)` are bit-identical to
`solve_pade(prob; h = alg.h, max_steps = alg.max_steps,
check_in_class)` — the loop condition, the signed clamped step, the
choice of checked/unchecked stepper and the checker's state evolve
identically, so `z`, `y`, `h` and every dense-evaluation value agree
with `==` (tested on ascending and descending spans, guard on and off).

## Fail-fast contract

  - `init` enforces `alg.h > 0` and the same 2nd-order requirement
    as `solve_pade`.
  - `step!` on a `done` integrator is a no-op (returns the integrator
    unchanged); not an error.  Matches `OrdinaryDiffEq.jl` semantics.
  - Hitting `max_steps` mid-`step!` throws `ErrorException`.

## References

  - `docs/adr/0003-extensions-pattern.md` — design rationale.
  - `Project.toml` `[weakdeps]` + `[extensions]` — declares this ext.
  - `src/Problems.jl::solve_pade` — the reference driver this layer mirrors.
"""
module PadeTaylorCommonSolveExt

using PadeTaylor:             PadeTaylorProblem, PadeTaylorSolution, PadeTaylorAlg
using PadeTaylor.RobustPade:  PadeApproximant
using PadeTaylor.PadeStepper: PadeStepperState, pade_step_with_pade!
using PadeTaylor.OutOfClass:  OutOfClassChecker, pade_step_with_defect!, check_in_class!
import CommonSolve

# =============================================================================
# Streaming integrator
# =============================================================================

"""
    PadeTaylorIntegrator{F, T, Y, P, H}

Streaming integrator state.  Holds the problem, algorithm, inner
`PadeStepperState`, and the four parallel vectors that accumulate the
trajectory, plus the integration direction `dir = sign(z_end - z_start)`
and the optional out-of-class checker (`nothing` when the guard is off).
`done` flips `true` once `dir · (zspan[2] - state.z) ≤ 0`.

Not part of the public API; obtain via `init(prob, alg)` and consume
via `step!`/`solve!`.
"""
mutable struct PadeTaylorIntegrator{F, T, Y, P, H <: Real}
    prob     :: PadeTaylorProblem{F, T, Y}
    alg      :: PadeTaylorAlg{H}
    state    :: PadeStepperState{T}
    z_vec    :: Vector{T}
    y_vec    :: Vector{Y}
    h_vec    :: Vector{T}
    pade_vec :: Vector{P}
    dir      :: T
    checker  :: Union{Nothing, OutOfClassChecker}
    steps    :: Int
    done     :: Bool
end

# =============================================================================
# CommonSolve interface
# =============================================================================

function CommonSolve.init(prob::PadeTaylorProblem{F, T, Y},
                          alg::PadeTaylorAlg{H};
                          check_in_class::Bool = true) where {F, T, Y, H <: Real}
    alg.h > 0 || throw(ArgumentError(
        "PadeTaylorAlg: h must be positive (got $(alg.h)).  " *
        "Suggestion: pass a strictly-positive step length."))
    Y <: Tuple || error(
        "init(prob, ::PadeTaylorAlg): 1st-order (scalar y0) branch is not " *
        "implemented in v1.  Suggestion: rewrite as a 2nd-order system " *
        "with `y0 = (u0, up0)`, or file a bead requesting 1st-order support.")

    z_start, z_end = prob.zspan
    state   = PadeStepperState{T}(z_start, prob.y0[1], prob.y0[2])

    P_T      = PadeApproximant{T}
    z_vec    = T[z_start]
    y_vec    = Y[prob.y0]
    h_vec    = T[]
    pade_vec = P_T[]

    # Direction and guard exactly as `solve_pade` (src/Problems.jl): `dir` is
    # ±1 (the constructor rejects z_start == z_end, so never zero); the
    # checker is `nothing` when the guard is disabled.  `done` is the
    # negation of solve_pade's loop condition, so for any legal problem the
    # integrator is born NOT done (bug padetaylor-ncfa).
    dir     = sign(z_end - z_start)
    checker = check_in_class ? OutOfClassChecker() : nothing
    done    = !(dir * (z_end - state.z) > zero(T))

    return PadeTaylorIntegrator{F, T, Y, P_T, H}(
        prob, alg, state, z_vec, y_vec, h_vec, pade_vec, dir, checker, 0, done)
end

function CommonSolve.step!(integ::PadeTaylorIntegrator{F, T, Y, P, H}) where {F, T, Y, P, H}
    integ.done && return integ                      # no-op on done integrator
    integ.steps += 1
    integ.steps > integ.alg.max_steps && error(
        "step!: exceeded max_steps = $(integ.alg.max_steps) " *
        "(current z = $(integ.state.z), target z_end = $(integ.prob.zspan[2])).  " *
        "Suggestion: raise max_steps, shorten zspan, or increase h.")

    z_end  = integ.prob.zspan[2]
    h_T    = T(integ.alg.h)
    dir    = integ.dir
    h_step = dir * min(h_T, abs(z_end - integ.state.z))   # signed; |h_step| ≤ h

    if integ.checker === nothing
        _, P_u = pade_step_with_pade!(integ.state, integ.prob.f, integ.prob.order, h_step)
    else
        _, P_u, δ = pade_step_with_defect!(integ.state, integ.prob.f,
                                           integ.prob.order, h_step)
        check_in_class!(integ.checker, δ, integ.state.z)
    end

    push!(integ.z_vec, integ.state.z)
    push!(integ.y_vec, (integ.state.u, integ.state.up))
    push!(integ.h_vec, h_step)
    push!(integ.pade_vec, P_u)

    if !(dir * (z_end - integ.state.z) > zero(T))
        integ.done = true
    end
    return integ
end

function CommonSolve.solve!(integ::PadeTaylorIntegrator{F, T, Y, P, H}) where {F, T, Y, P, H}
    while !integ.done
        CommonSolve.step!(integ)
    end
    return PadeTaylorSolution{T, Y, P}(integ.z_vec, integ.y_vec,
                                       integ.h_vec, integ.pade_vec)
end

end # module PadeTaylorCommonSolveExt
