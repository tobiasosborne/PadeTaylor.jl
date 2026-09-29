# Worklog 083 — Zero-component seeds end to end (padetaylor-0o9)

2026-09-29; worktree `w3-zerocomp`, branch `pace/w3-zerocomp`, base `b9f4132`. Read
`CLAUDE.md` and `AGENTS.md` at session start. No tracker writes, commits,
pushes, branch changes, dependency changes, full suite, or network operations.

## Result and measured scope

Canonical Type A/A₂, Type A/A₄ and Type B/A₄ work through Taylor jets → direct
SharedPade → dispatched VectorStepper → Float64 vector_solve_pade → dense
output. Orders 6 and 30, starts 0 and 1, h=1/8, spans of length 1, fixed and
Jorba–Zou policies: every checked value agrees with its closed form.
Default-order-30 cyclic rotations also work end to end. **Do not close 0o9
as universally fixed:** zero-first blocks reveal a surviving unpivoted-QR
column-selection bug, with silent errors in direct cell A and some order-6
vector solves. No numerical implementation was changed in this measurement task.

Ground truth: Matsuda's explicit seeds and parameters are in
`references/tex/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53/main.tex:324`
and `:328`; Type A for A₂ follows by substitution into the system at
`references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex:85`.
The test uses literal independent slopes rather than implementation-derived values (`test/zero_component_end_to_end_test.jl:51`).
GGT explains rank reduction and leading-factor cancellation at
`references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:134`
and `:136`; the Chebfun implementation was read in the read-only main checkout,
`external/chebfun/padeapprox.m:92`, `:111`, `:123`. Governing ADRs were read first:
`docs/adr/0027-sharedpade-graceful-reduction.md:35` and
`docs/adr/0028-sharedpade-dual-construction-pareto-dispatch.md:474`.

## Numbers

The new regression has 2,507 passing assertions, 0 failures/errors/broken.
It covers 64 trajectories, 512 segments, 576 nodes (including ICs), and 512
segment-midpoint dense samples. Canonical primitive/stepper checks also use
ComplexF64 with real-valued data; the complex driver is excluded for the
separate bug below. All identically-zero components remain exactly zero.

| Seed | Trajectories | Segments | Worst one-step error | Worst node error | Worst dense error |
|---|---:|---:|---:|---:|---:|
| A₂ Type A | 16 | 128 | 0 | 2.220446049250313e-16 | 2.220446049250313e-16 |
| A₄ Type A | 24 | 192 | 0 | 2.220446049250313e-16 | 2.220446049250313e-16 |
| A₄ Type B | 24 | 192 | 0 | 1.1102230246251565e-16 | 2.220446049250313e-16 |

The broader diagnostic scan measured 104 Float64 trajectories (832 segments),
including all rotations at both orders. 32 low-order rotated trajectories fail
closed-form accuracy; the 52 order-30 trajectories agree within 2.22e-16.
For the geometric jet `live[k+1]=0.5^k`, m=3:

| Input | Returned numerator(s) | Returned Q | Error in live value at t=1 |
|---|---|---|---:|
| `[live, zeros(7)]` | `[[1],[0]]` | `[1,-0.4999999999999999]` | roundoff |
| `[zeros(7), live]` | `[[0],[1]]` | `[1]` | 1 (returned 1, correct 2) |
| live first among 5 | live `[1]`, other `[0]` | `[1,-0.4999999999999999]` | roundoff |
| live last among 5 | live `[1]`, other `[0]` | `[1]` | 1 |

At order 6, rotated Type A/A₄ `(0,t,0,0,0)` starting at t=0 returns the all-zero
endpoint at t=1 (correct second entry 1); maximum node error 1, dense error
0.9375. Rotated Type B starting at zero similarly has node error 1/3, dense
error 0.3125. Both policies give these same failures. For Type A/A₂ starting
at t=1, order-6 rotations give endpoint errors 0.0056657364288108525 and
0.09726490403811239: square-cell selection can improve cell A without fixing it.

## Surviving mechanism and exact reproducer

The +2 matching window remains (`src/SharedPade.jl:160`). The all-jet norm is nonzero (`:235`) and its all-zero guard (`:246`) does not fire. Rank
reduction correctly reaches m_cur=1 (`:275`, `:283`). The failure occurs at
unpivoted `qr(adjoint(A_full*D))` and `F.Q[:,m_cur+1]` (`:313`, `:314`).
For `[zeros(7), live]`, the reduced matrix is `[0 0; 0.25 0.5]`. Its known
null vector `[1,-0.5]` has residual exactly 0; QR's selected column is `[0,1]`
with residual 0.5. Reweighting by D does not fix its direction. The first QR
column is zero, so the first m_cur columns do not span the live constraints;
rank(A_full)=m_cur alone does not justify selecting that Q column.

For a linear live jet, the same error selects b proportional to z. Numerator
recovery (`src/SharedPade.jl:320`) truncates `(z0+h*z)*z` at degree 1. Cancellation
(`:340`, `:347`, `:348`) then leaves the constant z0 and loses the slope.
The stepper calls the dual dispatcher (`src/VectorStepper.jl:251`): at default
order 30 the square cell recovers the measured seeds. At d>m it falls back
(`src/SharedPadeCellB.jl:116`, `src/SharedPadeDispatch.jl:99`). At order 6 for A₂,
both candidates can be wrong and the defect selector still returns one (`:112`).
This is a QR/null-space recovery bug, not the deleted Q(0) throw or zero-norm guard.

The following was run as temporary `test/_zero_component_counterexample_test.jl`
with every mutation restored: **0 pass / 6 fail / 0 error / 0 broken**. The file
is removed from the shipped test set; the exact known-correct RED reproducer
is preserved here for the repair bead (save it under test/ and run the command below).

```julia
using Test, PadeTaylor
using PadeTaylor.SharedPade: shared_denominator_pade
using PadeTaylor.NoumiYamada: NoumiYamadaProblem
using PadeTaylor.VectorStepper: VectorPadeStepperState, vector_pade_step_with_pade!
using PadeTaylor.VectorProblems: vector_solve_pade
@testset "0o9 surviving QR counterexamples (known-correct oracles)" begin
    live = [0.5^k for k in 0:6]
    nums, den = shared_denominator_pade([zeros(7), live], 3)
    @test sum(nums[2]) / sum(den) ≈ 2.0 atol=1e-14 rtol=0
    prob = NoumiYamadaProblem(2; α=[0.,1.,0.,0.,0.], f0=zeros(5),
                             tspan=(0.,1.), order=6)
    st = VectorPadeStepperState{Float64}(0.,zeros(5))
    vector_pade_step_with_pade!(st,prob.problem.f,6,.125)
    @test st.y ≈ [0.,.125,0.,0.,0.] atol=1e-14 rtol=0
    for policy in (:fixed,:jorba_zou)
        sol = vector_solve_pade(prob;h=.125,step_policy=policy)
        @test sol.y[end] ≈ [0.,1.,0.,0.,0.] atol=1e-14 rtol=0
        @test sol(.0625) ≈ [0.,.0625,0.,0.,0.] atol=1e-14 rtol=0
    end
end
```

Representative RED output: `Evaluated: 1.0 ≈ 2.0`; then vectors of five zeroes
versus `[0,0.125,0,0,0]`, `[0,1,0,0,0]`, and `[0,0.0625,0,0,0]`.

## Line-number audit and documentation changes

Checked the actual base files and git log, including `127edae` and `ccd6bad`.
The seven caveat/guard sites are reconciled, with measured qualifications:

- `src/NoumiYamadaSymmetry.jl:89`: end-to-end question now measured; `:74` no
  longer claims only Type C can be a vector oracle. It describes the surviving QR bug.
- PIV base `test/noumi_yamada_piv_test.jl:87` and `:262`: replaced blanket solver
  degeneracy claims with independent scalar-anchor rationale and default-order
  vector coverage; the same stale requirement at base `:139` was also corrected.
- Figure base `test/noumi_yamada_a4_figure_test.jl:72` and its header `:19`:
  all-nonzero IC is a fixture choice. Same duplicate corrected in
  `figures/_noumi_yamada_a4_helpers.jl:97`.
- Symmetry base `test/noumi_yamada_symmetry_test.jl:89`: Type C is one oracle;
  Types A/B default-order tests are now referenced, including the top chapter.
- SharedPade SP.1.5 base `:303`: prose was already corrected by `ccd6bad`, but
  `test/shared_pade_test.jl:318` still had the obsolete guard/non-test name.
  Renamed it and corrected the inline comment; its original assertion is preserved.
- SharedPade SP.5: base `:646`/`:653` had shifted to `:654`; the 'four throws'
  phrase was already explicitly historical. Preserve that accurate history;
  correct the still-live deleted-guard description at `test/shared_pade_test.jl:683`.
- ADR-0028 base `:427` and `:482`: replace the unmeasured zero-jet item with the
  measured support plus surviving ordering/QR issue. `README.md:175` states the
  measured public scope. No new algorithmic decision; the existing ADR is annotated.

Cancellation line numbers shifted: findfirst is now `src/SharedPade.jl:340`, slices `:347`/`:348`, all-zero numerator `:374`;
`:357` is now denominator-guard commentary. The norm and guard remain `:235`/`:246`.

## Mutation proof

Every mutation ran only the new standalone regression. Each was restored before
any subsequent mutation; the final restored run is recorded below. Output columns
are pass / fail / error / broken, exactly as Julia reported.

| ID | Exact implementation perturbation | RED output | Load-bearing checks |
|---|---|---|---|
| M1 | SharedPade `iszero(cnorm)` → `any(jet -> all(iszero,jet), jets)` | 3 / 0 / 2 / 0 | Mixed-zero acceptance; canonical direct construction |
| M2 | Both SharedPade zero-numerator branches `T[zero(T)]` → `T[one(T)]` | 2091 / 19556 / 0 / 0 | Zero numerators, zero state legs, step/node/dense values |
| M3 | VectorStepper `_rescale_by_powers(jet,h_T)` → `_rescale_by_powers(jet,one(T))` | 1439 / 1068 / 0 / 0 | Known live slopes; step/node/dense values |
| M4 | VectorCoefficients bootstrap assignment → `y[i][k]=zero(T)` | 11 / 1 / 1 / 0 | Complete analytic Taylor jets, including live derivative |
| M5 | VectorProblems dense `t=(z_T-sol.z[k])/h_k` → `zero(T)` | 1995 / 512 / 0 / 0 | Every dense midpoint; node checks remain green |
| M6 | Rational-oracle returns `(α,f)` → `(2 .* α, t -> 2 .* f(t))` | 2501 / 6 / 0 / 0 | Independent exact parameter/formula anchors |
| M7 | SharedPade `Vector{T}(b ./ b1)` → `Vector{T}(b)` | 2261 / 246 / 0 / 0 | Known shared denominator and normalisation |
| M8 | VectorStepper `state.z = state.z + h_T` → `state.z = state.z + 2h_T` | 947 / 536 / 0 / 0 | Point advancement, segment counts, node/dense coordinates |

M2 changes the dynamics and thus the adaptive segment count, explaining its
larger assertion total; M8 halves segment counts. M1/M4 surface fail-loud errors; M3/M5 produce incorrect
closed-form values. No pre-existing assertion or tolerance was weakened/deleted.
The exploratory full-ordering test initially also required Q=[1] for square-cell
rotations; that representation is nonunique, so shipped rotation checks assert
closed-form values. Its genuine failing cases remain in the reproducer above.

## Validation and remaining work

All Julia calls used, from this worktree root:
`JULIA_DEPOT_PATH="$PWD/.depot:$HOME/.julia" julia --project=. test/<file>_test.jl`.
Julia 1.12.3; one top-level Julia process at a time. juliaup's launcher lock was
read-only in the sandbox; the same prescribed command ran through approved
sandbox escalation. Private depot stayed in this worktree. No weak deps required.
Final restored validation: all six files pass (3,033 assertions; 0 fail/error/broken).

| Standalone test file | Pass | Fail/error/broken |
|---|---:|---|
| `zero_component_end_to_end_test.jl` | 2507 | 0/0/0 |
| `shared_pade_test.jl` | 192 | 0/0/0 |
| `vector_stepper_test.jl` | 55 | 0/0/0 |
| `noumi_yamada_symmetry_test.jl` | 175 | 0/0/0 |
| `noumi_yamada_piv_test.jl` | 37 | 0/0/0 |
| `noumi_yamada_a4_figure_test.jl` | 67 | 0/0/0 |

Restoration verified byte-for-byte for all five perturbed source files; only the
intended NoumiYamadaSymmetry docstring differs from HEAD. `git diff --check` passes.

The original exploratory run stopped on `MethodError: isless(::Float64,::ComplexF64)`
in `src/VectorProblems.jl:307`: h_T is complex at `:263`, but gap_mag is real at
`:300`, so `min(h_T,gap_mag)` fails. Reproduce with a ComplexF64 canonical seed
problem and `vector_solve_pade(prob;h=0.125)`. This is outside 0o9; report/re-file
when complex vector IVPs are required. Complex primitive/stepper results do not
establish complex driver support.

Keep/re-scope 0o9 to “unpivoted QR returns a non-null denominator for zero-first
component blocks”. A repair is forced by arbitrary component ordering, lower
orders, or direct SharedPade correctness; use the six RED assertions above and
independent `[1,-0.5]` null vector as its acceptance oracle. The square-cell
rank-deficient recovery and accepting two bad defect scores also need analysis.
No bd/bead edits, per the task override. BigFloat, non-real complex steps, arbitrary
orders, and pole-field walks of zero-component seeds were not measured. Register only `test/zero_component_end_to_end_test.jl`.
