# Worklog 021 — Classical Padé via Toeplitz `\` shipped as the F64 default (bead `padetaylor-txg`)

**Date**: 2026-05-13 (follow-on to worklog 020)
**Author**: Claude Opus
**Scope**: Implementation of bead `padetaylor-txg`'s acceptance (a)–(f).
Adds `classical_pade_diagonal` to `src/RobustPade.jl` and an
element-type-driven `method::Symbol` dispatch to `robust_pade`.
Tightens PN.2.2 + PN.2.3 long-range rtols by 1000× and 100,000×
respectively.  All four mutation-proofs bite.

> **Take-home**:  Worklog 020 diagnosed the F64 long-range accuracy
> gap as algorithmic (GGT 2013 SVD's robustness machinery is wasted
> on smooth ℘-trajectory blocks); this session ships the fix.
> `robust_pade(c, m, m)` now dispatches by `T`: `:classical` (FW 2011
> §5.1.4 Toeplitz `\`) for F32 / F64 / Complex variants, `:svd` (GGT
> 2013 Algorithm 2) for BigFloat / Arb / generic AbstractFloat.
> Off-diagonal `(m ≠ n)` and singular Toeplitz both auto-route to the
> SVD path.  Test suite **1314 → 1377 GREEN** (+63 CP.1.*).  Headline
> empirical at the full PathNetwork (Stage 1 + Stage 2): z=10⁴ F64
> rel-err `1.4·10⁻¹⁰` — beats FW Table 5.1's published `2.34·10⁻¹⁰`
> by 1.7× (and is ~50,000× tighter than the prior SVD-default
> `6.05·10⁻⁶`).

## What changed

### `src/RobustPade.jl`

  - **New** `classical_pade_diagonal(c, m)` (~90 LOC).  Pure FW 2011
    §5.1.4: pad `c` to length `2m+1`; build the `m × m` Toeplitz
    `T_mat[i,j] = cv[m+i-j+1]` (eq. 5.4); RHS `-cv[m+2 : 2m+1]`; solve
    via `lu(T_mat; check=false)` + `\`; throw `SingularException(0)`
    on `!issuccess(F)`.  Numerator via the compact convolution form
    of eq. (5.5): `a_k = Σ_{j=0..min(k,m)} b_full[j+1] · cv[k-j+1]`.
    No trim, no rank-counting, no QR-reweighting; well-conditioned
    cases come out clean and the singular case throws.
  - **New** `_default_pade_method(::Type{T})` dispatcher: `:classical`
    for F32/F64/Complex variants, `:svd` everywhere else.
  - **`robust_pade`** gains `method::Symbol = _default_pade_method(T)`
    kwarg.  Validates against `(:classical, :svd)`.  When
    `method == :classical && m == n && n > 0`, tries
    `classical_pade_diagonal`; on `SingularException` falls through to
    the existing SVD path.  Off-diagonal `(m ≠ n)` or `n == 0` skip
    the classical branch entirely.  Other exception types rethrow.
  - **Exports**: `classical_pade_diagonal` added.
  - **Imports**: `lu`, `issuccess`, `SingularException` added from
    `LinearAlgebra`.
  - **Module docstring**: rewritten top-of-file to introduce the
    two-method architecture; cross-links FW 2011 lines 346 + 350 +
    eqs. 5.4/5.5, GGT 2013 §2 + Algorithm 2, ADR-0001 / 0002 / 0005,
    worklog 020.

### `test/classical_pade_test.jl` (new, ~260 LOC)

9 testsets / 63 assertions:

| testset | covers                                                                                          |
|---------|-------------------------------------------------------------------------------------------------|
| CP.1.1  | `exp(z)` Padé(2,2) closed form via classical + dispatch + `r_cls ≡ r_svd` sweep                |
| CP.1.2  | `exp(z)` Padé(15,15) classical retains full diagonal, function-value to ~eps at \|z\| ≤ 0.5    |
| CP.1.3  | `log(1.2 - z)` Padé(10,10) function-value match                                                 |
| CP.1.4  | `(1+z²)` Padé(1,1) singular Toeplitz → SVD fallback collapse to `(0,0)`                         |
| CP.1.5  | Off-diagonal `(3, 2)` Padé with `method=:classical` auto-routes to `:svd`                       |
| CP.1.6  | Element-type dispatch: F64 default → classical retains (20,20); F64 `:svd` reduces to (7,7);    |
|         | Complex{F64} default → classical; BigFloat default → SVD; BigFloat `:classical` override works  |
| CP.1.7  | `1/(1-z/2)` (rank-1 Toeplitz at high m) → SingularException → SVD fallback recovers `(0, 1)`    |
| CP.1.8  | Fail-fast on negative `m` and unknown `method`                                                  |
| CP.1.9  | `Complex{Float64}` classical Padé sanity                                                        |

Plus a documented mutation-proof block at the bottom (CP.1.M) listing
the four mutations and their verified bite counts.

### `test/runtests.jl`

Single line: `include("classical_pade_test.jl")` added between
`robustpade_test.jl` and `coefficients_test.jl`.

### `test/robustpade_test.jl`

Tests 2.1.2 (`exp(z)` 20→7 reduction), 2.1.3 (`log(1.2-z)` 20→10
reduction), 2.1.4 (`tan(z⁴)` Froissart removal), 2.1.6 (noisy
`1/(1-z)` tol-thresholded recovery) gained explicit `method = :svd`
kwargs.  These tests assert GGT 2013 Algorithm 2 specific behaviour
(diagonal hopping, QR reweighting, tol-thresholded rank counting)
that classical doesn't have.  Comments added explaining why.

### `test/pathnetwork_test.jl`

  - PN.2.2 z=30 F64: rtol `1e-9 → 1e-12`, imag bound `1e-9 → 1e-12`.
    Measured under classical: well below `1e-12` on the canonical run.
  - PN.2.3 z=10⁴ F64: rtol `5e-5 → 5e-10`, imag bound `5e-5 → 1e-8`.
    Measured: rel-err `1.4·10⁻¹⁰`, `|imag(u)| 2.6·10⁻⁹`.  Margins
    ~3.6× and ~3.8× respectively.

Both rtol tightenings are direct consequences of the dispatch change;
mutation P3 (revert F64 default to `:svd`) bites both as expected.

### `docs/adr/0005-classical-pade-default-at-float64.md` (new)

The dispatch decision (element-type-driven), the empirical evidence,
when SVD remains load-bearing (Froissart doublets / singular Toeplitz
/ arb-prec block boundaries), and the explicit deferral of FW 2011
line 346's row-removal min-norm fallback in favor of routing singular
F64 cases to `:svd`.

### `HANDOFF.md`

Refreshed test count (1314 → 1377), session summary, hard-won
lessons (39–41), open-bead state.

## Mutation-proof procedure

All four mutations bit then restored cleanly (1377/1377 GREEN
post-restore).

**Mutation P1 — flip sign of RHS in `classical_pade_diagonal`**.
Changed `rhs[i] = -cv[m + i + 1]` to `rhs[i] = cv[m + i + 1]`.
The classical Padé denominator polynomial gets reversed, so every
downstream `u(z+h)` evaluation is wrong.

  - Direct bite: CP.1.1 (coefficient checks + dispatch equality +
    5 `r_cls ≈ r_svd` checks), CP.1.3 (2 log function-value asserts),
    CP.1.9 (Complex{F64} exp).
  - Cascade bite: Phase 5 PadeStepper (8 fails — every step in the
    inner integrator is corrupted), Phase 6 Problems (9 fails across
    6.1.1–6.1.5), PathNetwork (PN.1.2 fail + PN.2.2 / PN.2.3 errors —
    long walks hit Inf/NaN), Phase 9 PI tritronquée (2 fails),
    LatticeDispatcher LD.1.2 (51 fails — cosh test fails at every
    BVP-fill cell).
  - **Total RED ~60+**.  Classical's correctness is load-bearing
    across the entire IVP stack.

**Mutation P2 — remove the singular-fallback `try`/`catch`**.
Changed the `try classical_pade_diagonal(c, m) catch e ... end` block
inside `robust_pade(... ; method = :classical)` to a direct return.

  - Bite: CP.1.4 (1f + 1e — `robust_pade(c=[1,0,1], 1, 1;
    method=:classical)` now throws instead of returning `(0,0)`
    collapse), CP.1.7 (1f + 1e — same on the rank-1 geometric series),
    RobustPade 2.1.5 (1 error — `(1+z²)` defect-1 test uses F64
    default which is now classical-only without fallback), Dispatcher
    DP.2.1 + DP.3.1 (1f + 2e — segment junctions hit singular cases
    in some configurations), LatticeDispatcher LD.2.1 (1f — fail-fast
    guard returns `SingularException` instead of expected
    `ArgumentError`).
  - **Total RED 5**.  Singular fallback is load-bearing.

**Mutation P3 — flip F64 / Complex{F64} default to `:svd`**.  Changed
the four `_default_pade_method` overloads for F32/F64/Complex variants
to return `:svd` instead of `:classical`.

  - Bite: CP.1.6 (4 of 10 — F64 default and Complex{F64} default
    `P.μ == 20` assertions fail because SVD reduces to (7, 7)),
    PN.2.2 (2 — F64 rel-err `6.6·10⁻¹²` now exceeds new `1e-12` rtol,
    imag bound also fails), PN.2.3 (2 — F64 rel-err `6.05·10⁻⁶` blows
    through `5·10⁻¹⁰`).
  - **Total RED 8**.  Proves BOTH the dispatch-default change AND the
    PN.2.2 / PN.2.3 rtol tightening are load-bearing on classical
    being the default for F64.

**Mutation P4 — drop the `m == n` guard so off-diagonal also routes
through `classical_pade_diagonal`**.

  - Bite: CP.1.5 (4 of 6 — `P_cls.μ == P_svd.μ` fails because
    `classical_pade_diagonal` always returns a `(m, m)` approximant,
    so the (3, 2) request comes back as (3, 3); function values at
    the four `z` samples disagree).
  - **Total RED 4**.  Off-diagonal routing is load-bearing.

Restoration: each mutation restored before the next was applied; full
suite GREEN at 1377/1377 after the last restore.

## Frictions surfaced

**F1.  CP.1.1's initial tolerance was confused about which error
mode it was checking.**  My first draft asserted `_eval_pade(P_cls,
z) ≈ exp(z)` to `1e-14` at `z ∈ {-0.3, -0.1, 0.1, 0.3, 0.5}`.  But
Padé(2,2) of `exp` has truncation error `~z⁵/120` ≈ `3·10⁻³` at
`|z| = 0.5` — the algorithm's intrinsic accuracy ceiling at order 2,
NOT floating-point roundoff.  10 out of 21 assertions failed for
this entirely-predictable reason.  Fix: removed the `exp(z)` checks
from CP.1.1; the `r_cls ≈ r_svd` check (which holds to ~eps
regardless of order) is the load-bearing assertion.  Lesson 41
captures the principle.

**F2.  CP.1.7 picked a degenerate input by accident.**  Setting up a
test for "well-conditioned high-`m` classical Padé", I used `c_k =
(1/2)^k` — the Taylor series of `1/(1 - z/2)`.  The resulting `m × m`
Toeplitz has `T[i, j] = α^(m+i-j)` with `α = 1/2`: every row `i` is
`α^(i-1) ·` row 1, so `rank(T) = 1` → exactly singular.
`classical_pade_diagonal` correctly threw `SingularException`, so my
test errored.  The "well-conditioned" property of a Taylor series
(convergence radius) is NOT the same as well-conditioning of the
Padé Toeplitz it builds.  Fix: pivoted CP.1.7 from a positive
no-reduction test to a *positive singular-fallback test* — `robust_pade
(c, 15, 15; method = :classical)` catches the throw and recovers
the true `(0, 1)` Padé.  This complements CP.1.4's `(1+z²)` defect-1
case (singular at low `m` from a polynomial input) with a different
mechanism (singular at high `m` from a smooth rational input).

**F3.  Worklog 020's wedge-walker rel-err was looser at full
PathNetwork.**  Worklog 020 measured `6.15·10⁻¹¹` at z=10⁴ F64 via
a custom 5-direction wedge walker (Stage 1 only — IC-to-target
single chain).  Production code's `path_network_solve` adds Stage 2
barycentric interpolation, which is one extra Padé eval per
target.  Measured under the production path: `1.4·10⁻¹⁰` — ~2.3×
looser, but still beating FW's published `2.34·10⁻¹⁰` by 1.7×.
Bead acceptance said "~1e-10"; setting rtol = `5e-10` with ~3.6×
margin matches the spirit while staying ahead of cross-platform
variance.  Lesson 40 captures.

## What is NOT shipped (explicit deferrals)

  - **FW 2011 line 346's row-removal min-norm fallback**.  We route
    singular F64 Toeplitz to `:svd` instead.  GGT 2013 was designed
    as the principled solution to the same condition; we already pay
    for the SVD path's existence (it's the BigFloat default).  Adding
    FW's specific min-norm path is a second code path with no clear
    advantage on the F64 cases we hit empirically.  If a downstream
    caller surfaces a case where the FW row-removal beats `:svd` on
    F64, file a follow-up bead.
  - **Off-diagonal classical**.  `classical_pade_diagonal` handles
    only `(m, m)` per FW 2011 §5.2.  General `(m, n)` admits the
    same Toeplitz construction with `T[i, j] = c_{m+i-j}` for `i,j
    ∈ 1:n`, but our use sites (PadeStepper, PathNetwork) all request
    diagonal Padé, so the off-diagonal case has no consumer.  Route
    to `:svd` covers the few existing tests that exercise off-diagonal.
  - **Adaptive `method` selection**.  We could detect rank deficiency
    *before* the LU and skip straight to SVD without the try/catch
    round-trip.  Empirically the throw-catch is negligible cost
    (~0.1% of `Pkg.test()` wall time), so not worth the complexity.
  - **Tightening PN.2.2 below `1e-12`**.  Cross-platform variance
    bound; could probe with a wider matrix of OS/LAPACK builds before
    committing to a tighter floor.  Not blocking.
  - **BF-256 sweep at z = 10⁴**.  Still deferred (worklog 019); the
    PN.2.3 BF-256 path is not in the routine suite due to multi-hour
    wall time at BigFloat-256 path-network walks.

## Pointers

  - `src/RobustPade.jl:155–245` — `classical_pade_diagonal` impl.
  - `src/RobustPade.jl:285–305` — `robust_pade` dispatch with
    `try`/`catch` fallback.
  - `src/RobustPade.jl:120–135` — `_default_pade_method`.
  - `test/classical_pade_test.jl` — CP.1.* testsets + mutation log.
  - `test/pathnetwork_test.jl:155–230` — PN.2.2 / PN.2.3 tightened.
  - `docs/adr/0005-classical-pade-default-at-float64.md` — design
    decision.
  - `docs/worklog/020-classical-pade-toeplitz-backslash.md` — the
    investigation worklog that triggered this implementation.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:330–350`
    — FW 2011 §5.1.4 ground truth (eqs. 5.4, 5.5, lines 346 + 350).

## Bead state

Closed in this session:

  - `padetaylor-txg` — Ship classical-Padé (Toeplitz \) as F64 default
    in RobustPade.  Acceptance (a)–(f) all satisfied.  Supersedes
    `padetaylor-u7o` (the four speculative attack vectors that worklog
    020 obsoleted).

## Hard-won lessons (for HANDOFF.md §"Hard-won")

**39.  Mutation cascades expose architectural load-bearing-ness.**
Mutation P1 (sign-flip on classical's RHS) was expected to bite ~10
CP.1.* assertions; it bit 60+ across Phase 5/6/9, PathNetwork, and
LatticeDispatcher.  The mutation count is a direct measure of how
deep an algorithm sits in the dependency graph; classical Padé is
consumed by every step in every IVP path-network walk.  When
predicting mutation scope, run the mutation rather than reasoning
top-down about the call graph.

**40.  Probe topology must match production code path.**  Worklog 020
probed via a custom Stage-1-only wedge walker (~`6.15·10⁻¹¹` z=10⁴
F64 rel-err); the production `path_network_solve` adds Stage 2
barycentric interpolation, which is one extra Padé eval per target
(~`1.4·10⁻¹⁰` measured).  ~2.3× looser, still well inside the
intended order of magnitude.  Lesson: when setting test rtols based
on a probe number, re-measure under the production code path before
committing.  Probes are correct for the architecture they model, not
for the architecture they intend to model.

**41.  Taylor-series well-conditioning ≠ Padé-Toeplitz
well-conditioning.**  CP.1.7's first attempt used `c_k = (1/2)^k`
(the Taylor series of `1/(1-z/2)` — a smooth analytic function with
infinite convergence radius) intending to demonstrate classical's
no-reduction property.  But the `m × m` Toeplitz built from this
geometric series has rank 1 (every row is a scaled copy of row 1),
so classical correctly threw `SingularException`.  The function's
analytic properties don't constrain the rank of the Toeplitz solver
matrix.  When choosing a non-degenerate test input, estimate the
Toeplitz's rank explicitly (or pick a function with provably
non-degenerate Padé table — e.g., a sum of multiple distant poles).
