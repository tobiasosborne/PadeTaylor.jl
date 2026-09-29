# 089 — CommonSolve extension: descending span, out-of-class guard (bead ncfa)

Date: 2026-09-29. Branch `pace/w9-commonsolve`.

## Reproduction (measured before any change)

Script: `PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, -1.5); order = 30)`,
`fW(z,u,up) = 6u^2`, `h = 0.5`.

    DESC CS nodes=1 z=[0.0]  solve_pade nodes=4 z=[0.0, -0.5, -1.0, -1.5]
    COMPLEX CS threw MethodError: no method matching isless(::ComplexF64, ::ComplexF64)
    COMPLEX solve_pade threw MethodError: no method matching isless(::ComplexF64, ::ComplexF64)

So the bead's first claim is confirmed: `solve(prob, alg)` silently returned
the one-node initial condition on a legal descending span, because
`init` set `done = state.z ≥ z_end` (old ext line 122) and `step!` clamped
with `min(h_T, z_end - z)` (old line 138). The bead's second claim is also
confirmed by reading: the ext always called the unchecked
`pade_step_with_pade!`, while `solve_pade` defaults `check_in_class = true`
(`src/Problems.jl:212`, `:255`).

Complex spans: NEITHER path supports them. `solve_pade`'s loop condition
`dir * (z_end - state.z) > zero(T)` (`src/Problems.jl:258`) is real-ordered and
throws `MethodError(isless, ComplexF64, ComplexF64)`. The fix for complex spans
is in `solve_pade` itself (out of this bead's scope), so this change gives the
extension the SAME behaviour and pins it (CS.5.3); see "What remains".

## Change

`ext/PadeTaylorCommonSolveExt.jl`:
- `init(prob, alg; check_in_class::Bool = true)` (`:127`) computes
  `dir = sign(z_end - z_start)` and `checker = check_in_class ?
  OutOfClassChecker() : nothing` (`:150-151`), exactly as `solve_pade`
  (`src/Problems.jl:242`, `:255`). `done = !(dir*(z_end - z) > 0)` (`:152`) is
  the negation of `solve_pade`'s loop condition, so no legal problem is born done.
- `PadeTaylorIntegrator` gains `dir::T` and `checker::Union{Nothing,OutOfClassChecker}`.
- `step!` uses `h_step = dir * min(h_T, abs(z_end - z))` (`:169`) and
  dispatches to `pade_step_with_defect!` + `check_in_class!` when the checker
  is present (`:171-177`), else the unchecked stepper — the same branch as
  `src/Problems.jl:267-272`.
- CommonSolve's default `solve(prob, alg; kw...) = solve!(init(prob, alg; kw...))`
  forwards `check_in_class`, so no `PadeTaylorAlg` field was added (keeps the
  public struct `src/PadeTaylor.jl:325` unchanged).
- Module docstring rewritten: direction logic, guard, the ncfa history, and the
  bit-identity claim now stated precisely and backed by tests.

`test/ext_commonsolve_test.jl`:
- CS.2.2 rewritten: the "unreachable through public API" comment was false;
  it now asserts a descending integrator is born live and its first step is
  `z = [0.0, -0.5]`, `h = [-0.5]`.
- CS.5.1: descending span, h = 0.5 and 0.4 (clamped final step), exact `==`
  on `z`, `y`, `h`, and dense values at three points; `z[end] == -1.5`;
  literal `z == [0.0, -0.5, -1.0, -1.5]`.
- CS.5.2: essential-singularity ODE `u'' = u(1+2z)/z^4` on `(-1.0, -0.02)`,
  `h = 0.1` (setup from `test/corpus_out_of_class_test.jl:104-113`): both
  `solve_pade` and `CommonSolve.solve` throw `OutOfClassError` by default; with
  `check_in_class = false` both return and agree with `==` on `z`, `y`, `h`.
- CS.5.3: complex span, both throw `MethodError` (shared limitation pinned).
- Ascending bit-identity was already pinned by CS.1.1 (unchanged, still green).

## Tests

`test/ext_commonsolve_test.jl` (run with CommonSolve supplied by a stacked
scratch environment, since it is a test-only extra): 64 pass / 0 fail / 0 broken.

## Mutation proof

- M1 `done = state.z ≥ z_end` in `init`: 44 pass, 14 fail, 6 error (CS.2.2, CS.5.1).
- M2 drop `dir *` from `h_step`: 45 pass, 2 fail, 1 error (CS.2.2, CS.5.1).
- M3 `checker = nothing`: 62 pass, 2 fail (CS.5.2).
All restored; `cmp` against the pre-mutation copy confirmed byte-identity.

## What remains

- Complex spans are unsupported in `solve_pade` itself (MethodError from a
  real-ordered comparison). A fix belongs in `src/Problems.jl` (e.g. a
  path-parameter formulation along the segment) and then here; CS.5.3 will go
  RED in the extension-vs-driver sense only if one path changes, which is the
  intent. Suggested new bead.
