"""
    PadeTaylorCommonSolveExt

Package-extension adapter wiring `PadeTaylorProblem` + `PadeTaylorAlg`
into the `CommonSolve.jl` `init` / `step!` / `solve!` / `solve` interface
per ADR-0003 (Stage 1 design lock).  Loaded automatically when both
`PadeTaylor` and `CommonSolve` are present:

```julia
using PadeTaylor, CommonSolve

prob = PadeTaylorProblem(f, (u0, up0), (0.0, 1.5); order = 30)
alg  = PadeTaylorAlg(; h_max = 0.5)

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
  - `init(prob, alg)` — constructs a `PadeTaylorIntegrator` wrapping
    the inner `PadeStepperState` + accumulators.  Pre-pushes the IC
    onto the trajectory.
  - `step!(integ)` — calls `pade_step_with_pade!` once, sized as
    `min(alg.h_max, z_end - state.z)` so the final step lands exactly
    on `z_end`.  Mutates the integrator; returns it for chaining.
    Sets `integ.done = true` when `state.z ≥ z_end`.
  - `solve!(integ)` — drives `step!` in a loop until `integ.done`;
    returns `PadeTaylorSolution{T, Y, P}` assembled from the
    accumulators.
  - `solve(prob, alg)` — `CommonSolve` provides a default
    `solve = solve! ∘ init`; we don't override.

This is a **translation layer**.  The trajectory bytes from
`solve(prob, alg)` are bit-identical to
`solve_pade(prob; h_max = alg.h_max, max_steps = alg.max_steps)`
modulo trivial evaluation-order differences (none here —
`pade_step_with_pade!` is deterministic and the integrator's `step!`
calls it with identical inputs).

## Fail-fast contract

  - `init` enforces `alg.h_max > 0` and the same 2nd-order requirement
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
import CommonSolve

# =============================================================================
# Streaming integrator
# =============================================================================

"""
    PadeTaylorIntegrator{F, T, Y, P, H}

Streaming integrator state.  Holds the problem, algorithm, inner
`PadeStepperState`, and the four parallel vectors that accumulate the
trajectory.  `done` flips `true` when `state.z ≥ zspan[2]`.

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
    steps    :: Int
    done     :: Bool
end

# =============================================================================
# CommonSolve interface
# =============================================================================

function CommonSolve.init(prob::PadeTaylorProblem{F, T, Y},
                          alg::PadeTaylorAlg{H}) where {F, T, Y, H <: Real}
    alg.h_max > 0 || throw(ArgumentError(
        "PadeTaylorAlg: h_max must be positive (got $(alg.h_max)).  " *
        "Suggestion: pass a strictly-positive step length."))
    Y <: Tuple || error(
        "init(prob, ::PadeTaylorAlg): 1st-order (scalar y0) branch is not " *
        "implemented in v1.  Suggestion: rewrite as a 2nd-order system " *
        "with `y0 = (u0, up0)`, or file a bead requesting 1st-order support.")

    z_start = prob.zspan[1]
    state   = PadeStepperState{T}(z_start, prob.y0[1], prob.y0[2])

    P_T      = PadeApproximant{T}
    z_vec    = T[z_start]
    y_vec    = Y[prob.y0]
    h_vec    = T[]
    pade_vec = P_T[]

    # Degenerate-zspan guard: if z_start ≥ z_end the integrator is born done.
    done = state.z ≥ prob.zspan[2]

    return PadeTaylorIntegrator{F, T, Y, P_T, H}(
        prob, alg, state, z_vec, y_vec, h_vec, pade_vec, 0, done)
end

function CommonSolve.step!(integ::PadeTaylorIntegrator{F, T, Y, P, H}) where {F, T, Y, P, H}
    integ.done && return integ                      # no-op on done integrator
    integ.steps += 1
    integ.steps > integ.alg.max_steps && error(
        "step!: exceeded max_steps = $(integ.alg.max_steps) " *
        "(current z = $(integ.state.z), target z_end = $(integ.prob.zspan[2])).  " *
        "Suggestion: raise max_steps, shorten zspan, or increase h_max.")

    z_end   = integ.prob.zspan[2]
    h_max_T = T(integ.alg.h_max)
    h_step  = min(h_max_T, z_end - integ.state.z)

    _, P_u = pade_step_with_pade!(integ.state, integ.prob.f, integ.prob.order, h_step)

    push!(integ.z_vec, integ.state.z)
    push!(integ.y_vec, (integ.state.u, integ.state.up))
    push!(integ.h_vec, h_step)
    push!(integ.pade_vec, P_u)

    if integ.state.z ≥ z_end
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
