# ADR-0036 — Order-independent shared-Padé null vector

Status: Proposed (2026-09-29); primitive validated, dispatch acceptance pending.
Bead: `padetaylor-brbv`.

## Ground truth and failure

GGT Algorithm 2 step 6 obtains the denominator from the null right singular
vector (`references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:228-236`).
The following paragraph explicitly describes QR reweighting as an addition
to that algorithm. Chebfun first takes `V(:,n+1)` at
`external/chebfun/padeapprox.m:106-109`, then refines it at `:111-117`.

ADR-0019:79-88 assumed that this scalar refinement ports unchanged to a stack.
It does not: the stack can be tall, and its first `m_cur` rows need not span
the constraints. With `[zeros(7), [0.5^k for k in 0:6]]`, degree reduction
correctly reaches 1 but the reduced matrix is `[0 0; 0.25 0.5]`. The old
QR column is `[0,1]`, with residual 0.5, while `[1,-0.5]` has residual 0.
The resulting live value at t=1 is 1 instead of the exact 2. Measurement and
the six failing end-to-end assertions are recorded in worklogs 083 and 088.

## Decision

For `d>1`, recover `b = conj.(Vt[end,:])` directly from the SVD already used
for rank classification. This follows GGT Algorithm 2 and preserves the
smallest-singular-value least-squares solution for full-column-rank stacks
accepted by the current ADR-0027 implementation. A component permutation
is a row permutation, leaving `A'*A` invariant; an isolated smallest right
singular vector therefore changes only by phase and roundoff. Normalising
`Q(0)=1` removes that phase. Conjugation is required because `Vt` is `V'`
(`src/LinAlg.jl:36-45`).

Retain Chebfun's reweighting for `d=1` to preserve the existing scalar oracle.
Keep degree reduction, numerator recovery, common-factor cancellation and
trimming unchanged. This supersedes ADR-0019's unchanged stacked-QR port;
it does not redesign ADR-0027 or ADR-0028's selector.

## Validation and limits

`test/shared_pade_order_invariance_test.jl` includes the exact six-assertion
reproducer and every cyclic rotation of real/complex one-pole, two-pole,
polynomial and full-column-rank jets. Literal rational formulas validate the
values independently of the decomposition. Existing shared-Padé, vector
stepper and Noumi–Yamada tests pass with their original tolerances. The
dispatch test has five failures: four expected cell choices change because
cell A improves, and one near-roundoff choice flips under epsilon perturbation.
Its unchanged-HEAD baseline passes 26/26. Acceptance remains incomplete;
the selector and existing assertions are unchanged. Worklog 088 records
independent oracle measurements, counts and both restored mutations.

The invariance argument requires a unique smallest singular direction.
Repeated smallest singular values in a full-column-rank stack do not define
a unique denominator; choosing among them is outside this repair.
`src/SharedPadeCellB.jl:137-138` contains the same unpivoted QR idiom. Its
matrix is `m_eff × (m_eff+1)` (`:129-130`), so the chosen final QR column
follows all constraints and remains null, including with zero rows. This
differs from cell A's tall stack, where only the first `m_cur` constraints
precede the chosen column. The geometric reproducer at m=3 gives cell-B
residuals 6.21e-17 and 0 for the two orders and live value 2 in both cases.
The denominators differ because the rank-1 matrix has a four-dimensional
null space. Cell B's separate reduction policy is unchanged.
