# 088 — Shared-Padé order-independent null vector

Date: 2026-09-29. Bead: `padetaylor-brbv`. Branch: `pace/w8-qrnull`.
Started at 21:35:46 CEST; final validation finished before 21:55; deadline 21:58.
Base commit: `420fb2f`. Changes left uncommitted; acceptance is partial.

## Ground truth before code

Read `CLAUDE.md` and `AGENTS.md` first; verified the clean worktree and git log.
The precise reproducer is `docs/worklog/083-zero-component-end-to-end.md:63-112`.
GGT Algorithm 2 step 6 uses the right singular vector:
`references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:228-236`.
The last paragraph explicitly identifies reweighted QR as an extra refinement.
Read the actual Chebfun file in the read-only main checkout:
`/home/tobiasosborne/Projects/PadeTaylor.jl/external/chebfun/padeapprox.m:106-117`.
It takes the SVD null vector at :109, then does scalar QR reweighting at :111-117.
Read ADR-0019:15-19,79-88 (unchanged stacked-QR assumption), ADR-0027:35-60
(degree reduction and cancellation), and ADR-0028:514-531 (cell B and dispatch).

## Root cause confirmed before changing the implementation

For `live = [0.5^k for k in 0:6]`, m=3, both `[live,zeros(7)]` and
`[zeros(7),live]` reduce correctly from m=3 to m=1. At m=3 the largest
singular values are respectively 0.6601446929461752 and 0.6601446929461753;
the zero-first stack's second singular value is only 1.5805203028961613e-17.
At m=1 both have singular values `[0.5590169943749475,0]` and rank 1.

| Order | Reduced A | Old unit QR b | ‖A b‖₂ | Old Q | Live value at t=1 |
|---|---|---|---:|---|---:|
| live first | `[.25 .5; 0 0]` | `[-.894427190999916,.44721359549995787]` | 5.551115123125783e-17 | `[1,-.4999999999999999]` | 1.9999999999999996 |
| zero first | `[0 0; .25 .5]` | `[0,1]` | 0.5 | `[1]` | 1 |

The SVD right vector has residual exactly 0 in both orders. The literal
known-correct vector `[1,-.5]` also has residual 0. The error is the final
unpivoted QR column, not rank reduction, trimming or the zero-norm guard.
For tall A, its first m rows need not span the live constraints; the chosen
Q column need not annihilate later rows. Full-column-rank stacks have no
null direction at all and need the smallest-singular-value vector.

Before the source edit, the new regression file gave **471 pass / 185 fail /
8 error / 0 broken** (664 assertions). Its exact worklog-083 subset gave
**0 pass / 6 fail / 0 error / 0 broken**. Representative RED output:

```text
Evaluated: 1.0 ≈ 2.0 (atol=1.0e-14, rtol=0)
Evaluated: [0,0,0,0,0] ≈ [0,0.125,0,0,0]
Evaluated: [0,0,0,0,0] ≈ [0,1,0,0,0]
Evaluated: [0,0,0,0,0] ≈ [0,0.0625,0,0,0]
```

## Smallest correct change

`src/SharedPade.jl:310-322` now takes `conj.(Vt[end,:])` from the existing
accepted SVD for d>1. Conjugation is essential because Vt is V', as already
documented in `src/LinAlg.jl:36-45`. For d=1, retain the original Chebfun
reweighting so the scalar oracle remains unchanged. No extra decomposition,
component filtering or ordering heuristic is needed. Degree reduction,
numerator recovery, cancellation and trimming are unchanged.
The direct reproducer now returns exactly Q=`[1,-.5]` and live value 2 for
both orders. Both fixed and Jorba–Zou order-6 rotated Type A solves recover
the exact step, endpoint and dense-output values at the original 1e-14 tolerance.

The module chapter and public function docstring explain this choice at
`src/SharedPade.jl:80-94,203-209`. The public API is updated at `README.md:178`;
ADR-0036 records the algorithmic decision. The existing zero-component test's
chapter now links to the repaired counterexamples (`test/zero_component_end_to_end_test.jl:18-23`).

## Permutation property and independent oracle

`test/shared_pade_order_invariance_test.jl:32-46` ports the six assertions
verbatim in substance, keeping every original tolerance.
`:49-100` checks all cyclic rotations of five families in Float64 and
ComplexF64: a geometric jet plus zero; five mixed/duplicate components;
two shared poles (genuinely complex in ComplexF64); polynomial jets; and a
well-conditioned full-column-rank least-squares stack. There are 16 rotations
per type, 32 total, evaluated at t=-.5, .25 and 1. Denominators and mapped
per-component values must agree with the unrotated result at `64eps()`.
Literal P/Q polynomials and formal division (`:20-30`) supply independent
known-correct values for the rational/polynomial families.

| Type | Maximum Q coefficient difference | Maximum permuted value difference | Maximum literal-oracle value error |
|---|---:|---:|---:|
| Float64 | 1.3322676295501878e-15 | 6.217248937900877e-15 | 4.884981308350689e-15 |
| ComplexF64 | 1.3322676295501878e-15 | 6.217248937900877e-15 | 1.7763568394002505e-15 |

The restored new file passes **664 / 0 fail / 0 error / 0 broken**.

## Mutation proof (restored)

M1: replace `if d == 1` at `src/SharedPade.jl:313` with
`if true # MUTATION M1: restore QR recovery for stacked blocks.`
This exactly restores the faulty denominator algorithm for every component
count (the absolute weights are unchanged by conjugation). Run the new file:
**471 pass / 185 fail / 8 error / 0 broken**. All six reproducer assertions
fail; the rotation property gives 471 pass / 179 fail / 8 error. The representative
`1.0 ≈ 2.0` RED output above is reproduced after an initial GREEN run.

M2: replace `b = Vector{T}(conj.(Vt[end, :]))` at :312 with
`b = Vector{T}(Vt[end, :]) # MUTATION M2: omit the complex conjugation.`
Run the new file: **645 pass / 19 fail / 0 error / 0 broken**. All six real
reproducers pass; the complex rational denominator and value oracle fail:

```text
Evaluated: [1, -0.49999999999999956+0.24999999999999986im,
            6.397344083151129e-17-0.12499999999999999im]
         ≈ [1, -0.5-0.25im, 0+0.125im]
```

After each mutation, restored the saved GREEN source. `cmp` confirmed exact
restoration. After both restores, reran the final new file: **664 pass**, and
reran the requested existing files. No mutation remains in source.

## Same idiom in Cell B: inspected and measured

Confirmed `src/SharedPadeCellB.jl:137-138` uses unpivoted QR too. Its shape
differs: A is m_eff×(m_eff+1), and its adjoint has only m_eff columns
(`:129-130`). The selected final Q column follows **all** constraints, so it
is null even with dependent/zero rows. Cell A's tall stack has additional
constraints after the selected column.
For the same geometric jets and m=3, cell B has m_eff=4 and rank 1. Both
orders produce correct live values: 2.000000000000002 and 2.0. Measured QR
residuals are 6.206335383118183e-17 and 0. The denominators differ (degrees
4 and 2 after cancellation), consistent with a four-dimensional null space.
This measurement does **not** establish a cell-B correctness bug. Cell B was
left unchanged, as requested; its degenerate reduction policy is separate.

## Tests and session constraints

All Julia runs were sequential standalone files from this worktree root:

```sh
JULIA_DEPOT_PATH="$PWD/.depot:$HOME/.julia" julia --project=. test/<file>_test.jl
```

Julia 1.12.3. The first sandboxed invocation failed before starting Julia
because juliaup could not create its read-only launcher lock. The same
prescribed command then ran with approved sandbox escalation, as in worklog
083:168-172. All precompile output used this worktree's private depot.
No weak dependencies, package operations, network requests, full suite,
quality gates, beads edits, commits or pushes. Existing tolerances/assertions
were not changed; only the zero-component file's documentation was updated.
Temporary measurement files are removed from the worktree after recording.

| File under test/ | Pass | Fail | Error | Broken |
|---|---:|---:|---:|---:|
| shared_pade_order_invariance_test.jl | 664 | 0 | 0 | 0 |
| shared_pade_test.jl | 192 | 0 | 0 | 0 |
| zero_component_end_to_end_test.jl | 2507 | 0 | 0 | 0 |
| vector_stepper_test.jl | 55 | 0 | 0 | 0 |
| shared_pade_dispatch_test.jl | 21 | 5 | 0 | 0 |
| noumi_yamada_test.jl | 102 | 0 | 0 | 0 |
| noumi_yamada_symmetry_test.jl | 175 | 0 | 0 | 0 |
| noumi_yamada_piv_test.jl | 37 | 0 | 0 | 0 |

## Dispatch acceptance remains incomplete (no assertion/tolerance changed)

The unchanged HEAD source passes `shared_pade_dispatch_test.jl` **26/26**.
With the correct SVD primitive, its expected `:square` choices become
`:diagonal` at `test/shared_pade_dispatch_test.jl:87,89,92,95`. Its harmonic
epsilon-perturbation pick changes from diagonal to square at :136. These
are five regressions against the existing diagnostic-selection assertions.
Do not close the bead as fully accepted or claim all requested files green.

Measured literal-solution errors at t=.17,.41,.73,1 (maximum over components):

| Fixture | Old cell-A error | SVD cell-A error | Cell-B error |
|---|---:|---:|---:|
| harmonic | 5.219115029042598e-9 | 4.3298697960381105e-15 | 5.10702591327572e-15 |
| exp / t·exp | 2.3402857707299773e-9 | 2.220446049250313e-16 | 4.440892098500626e-16 |
| Calogero–Moser | 9.656595661988732e-11 | 2.220446049250313e-16 | 4.440892098500626e-16 |
| tan companion | 113.68985860298797 | 8.781079330333341e-8 | 8.023860971206886e-7 |

The formulas are `(cos(ht),-sin(ht))`, `(exp(ht),ht·exp(ht))`,
`(q,-q,ht/(2q),-ht/(2q)), q=√(1+(ht)²/2)`, and `(tan(ht),sec²(ht))`,
verified by substitution into the fixture RHS at `test/shared_pade_dispatch_test.jl:86-95`.
Thus the four changed picks improve actual accuracy. The harmonic unperturbed
defects are A=5.1115995565334034e-14 and B=5.372190906474322e-14; after the
test's 50eps coefficient perturbation A=3.36764206077491e-14 and B wins.
The selector compares scores with a purely relative 100eps band
(`src/SharedPadeDispatch.jl:109-112`); an absolute rounding floor is separate work.
Additional weighted-SVD, QR-on-SVD-basis and pivoted-QR probes did not preserve
all existing picks. They were measured, then discarded; no variant was shipped.
Changing selection policy or its expected choices requires integration review.

## Remaining scope

Register `test/shared_pade_order_invariance_test.jl` in `test/runtests.jl` in
the orchestrator's integration pass; this lane was instructed not to edit it.
The invariance argument requires an isolated smallest singular direction.
A full-column-rank stack with a repeated smallest singular value has no
unique denominator; canonical handling is separate work if such an input
becomes an API requirement. Ill-conditioned denominators can amplify rounding
in coefficients and values; the property fixtures have separated directions.
The unrelated ComplexF64 solver-driver comparison remains outside this bead
(`docs/worklog/083-zero-component-end-to-end.md:187-191`).
