# ADR-0005 — Classical Padé via Toeplitz `\` as the `Float64` default

**Status**: Accepted (2026-05-13)
**Supersedes**: portions of ADR-0002 that implied SVD is universally
the right Padé-conversion backend.
**Context**: Bead `padetaylor-txg`; worklog 020's empirical diagnosis
and the FW 2011 §5.1.4 close-read it triggered.

## Decision

`PadeTaylor.RobustPade.robust_pade` dispatches to one of two Padé
algorithms via the `method::Symbol` kwarg.  The default is
element-type-driven:

| element type                          | default method | rationale                                                          |
|---------------------------------------|----------------|---------------------------------------------------------------------|
| `Float32`, `Float64`                  | `:classical`   | FW 2011 §5.1.4 Toeplitz `\` (LU/QR); fast and accurate at F64.      |
| `Complex{Float32}`, `Complex{Float64}` | `:classical`   | Same; the path-network's wedge walker is the dominant consumer.    |
| `BigFloat`, generic `AbstractFloat`   | `:svd`         | GGT 2013 Algorithm 2; Jacobi SVD's relative-accuracy is required.   |
| `Arblib.Arb`                          | `:svd`         | Routes through `BigFloat` per ADR-0002.                             |

A new entry point `classical_pade_diagonal(c, m)` implements FW 2011
eqs. (5.4) + (5.5) for the diagonal `(m, m)` case.  Off-diagonal
requests `(m ≠ n)` with `method = :classical` transparently route to
the `:svd` path.  Exactly singular Toeplitz `lu` factorisations
(`!issuccess(F)`) also fall through to `:svd` — GGT 2013's principled
treatment of the same rank-deficiency FW 2011 line 346 handled by
row-removal + min-norm.

Callers can override via `method = :svd` (force the robust path) or
`method = :classical` (force the classical path).

## Why this matters — the experimental case

Worklog 020 probed both algorithms on the FW 2011 test ODE `u'' = 6
u²` with the equianharmonic-℘ trajectory, at `Float64`, `h = 0.5`,
`order = 30`, via a 5-direction wedge walker mirroring PathNetwork
Stage 1:

| target | method     | wall (s) | rel-err   | vs FW Table 5.1     |
|--------|------------|----------|-----------|---------------------|
| z=30   | `:svd`     | 5.84     | 6.6e-12   | 87× worse           |
| z=30   | `:classical` | 0.01   | 1.54e-13  | 2× worse            |
| z=10⁴  | `:svd`     | 17.76    | 6.05e-6   | 25,800× worse       |
| z=10⁴  | `:classical` | 3.45   | **6.15e-11**  | **3.8× BETTER**     |

Per-step trace (step crossing the `z = 1` lattice pole): SVD `6.2·10⁻¹⁰`
absolute error, classical `7.1·10⁻¹³` — **870× per-step accuracy
improvement**.  At full PathNetwork (Stage 1 + Stage 2), the
end-to-end `z = 10⁴` `Float64` rel-err under classical is `1.4·10⁻¹⁰`
(measured 2026-05-13) — still beating FW's published `2.34·10⁻¹⁰` by
1.7×, and ~50,000× better than the SVD path.

## Why classical wins the smooth-conditioned case

GGT 2013 Algorithm 2 (the `:svd` path):

1. SVD of the `(n+1) × n` Toeplitz `C̃`.  At `Float64`, LAPACK
   Demmel–Kahan via stdlib `svd`.
2. Rank counting against `tol · ‖c‖₂` — diagonal-hopping until the
   numerical rank matches the requested `n`.
3. Null right-singular-vector → denominator.
4. Chebfun-style QR-reweighting refinement to preserve exact
   zeros at block boundaries (`padeapprox.m:111–117`).
5. Numerator from upper-block + trim leading/trailing near-zeros.

Costs ~`n² × n_sweeps` Givens rotations during the Jacobi/QR step
plus tolerance comparisons across the diagonal-hopping loop; each
rotation contributes ~1 ulp of roundoff.  At `n = 16` (`order = 30`,
`m = n = 15`), total per-Padé roundoff is `~n² · eps ≈ 5·10⁻¹³`
absolute error in the coefficients.  In the path-network walk this
compounds non-linearly near `℘` poles (worklog 019 §"Per-step error
probe": ~65× amplification at the `z = 1` pole crossing).

FW 2011 §5.1.4 classical Padé (the `:classical` path):

1. Build the `m × m` Toeplitz `T` (FW eq. 5.4).
2. LU with partial pivoting: `b_tail = T \ rhs`.  Solve cost
   `n³ / 3 ≈ 1100` flops at `n = 15`, single sweep.
3. Numerator via convolution `a_k = Σ c_{k-j} · b_j` (FW eq. 5.5).
4. Done — no rank check, no QR-reweighting, no trim.

Roundoff is `~n · eps ≈ 3·10⁻¹⁵` per coefficient.  Both faster and
more accurate **for the case where SVD's robustness machinery is
unused** — smooth, well-conditioned ℘-trajectory blocks.

## When SVD remains load-bearing

Three regimes where `:classical` is incorrect or worse and the
default for arb-prec keeps `:svd`:

1. **Froissart doublets** — GGT 2013 §3 explicitly designed SVD to
   suppress spurious pole–zero pairs that `\` can produce.  Rare in
   path-network walks; possible for noisy or pathological input
   sequences.
2. **Singular Toeplitz** — `\` returns `Inf`/`NaN` (or `lu` reports
   `!issuccess`); the `:svd` path gracefully reduces to the correct
   block.  We handle this case via the auto-fallback inside
   `robust_pade(... ; method = :classical)`.
3. **Arbitrary precision** — at `BigFloat`-256 the dynamic range
   opens up; near a Padé table block boundary the gap between the
   smallest "signal" singular value and the largest "noise" singular
   value may span only a few orders of magnitude.
   `GenericLinearAlgebra`'s one-sided Jacobi (Demmel–Veselić) is
   load-bearing for the `σᵢ < tol · ‖c‖₂` threshold in this regime
   (ADR-0002).  LU's absolute-accuracy guarantees would mis-classify
   small genuine singular values as "noise".

For `Float64` the dynamic range is narrow, the rank-counting is
unambiguous, and the SVD robustness is wasted work.

## What about FW's own singular fallback?

FW 2011 line 346:

> In rare cases, it can happen that (5.4) becomes outright singular.
> [...] we remove the last row of (5.4). Matlab's `\` solver will then
> give the solution to the underdetermined system that has the
> smallest norm, thereby restoring numerical stability.

We *could* port this as a second-tier fallback (after LU rejects but
before SVD).  We chose not to: GGT 2013 was specifically designed as
the principled solution to the same problem, and we already pay for
the SVD path's existence (it's the default at `BigFloat` anyway).
Routing singular F64 cases to `:svd` rather than to FW's row-removal
min-norm is one fewer code path to maintain.  If empirical evidence
later shows row-removal min-norm gives meaningfully better results on
the singular F64 case, we'll revisit; the v1 acceptance does not need
it.

## Test coverage (CP.1.*)

`test/classical_pade_test.jl` — 9 testsets / 63 assertions:

| testset | covers                                                                                            |
|---------|---------------------------------------------------------------------------------------------------|
| CP.1.1  | `exp(z)` Padé(2,2) closed-form via both methods; coefficient + `r_cls ≡ r_svd` match.            |
| CP.1.2  | `exp(z)` Padé(15,15) classical retains full degree, function-value to ~eps.                       |
| CP.1.3  | `log(1.2 - z)` Padé(10,10) classical function-value (no SVD-style block reduction).               |
| CP.1.4  | Singular Toeplitz at `(1+z²)` Padé(1,1): SingularException + SVD fallback recovers `(0,0)`.       |
| CP.1.5  | Off-diagonal `(3, 2)` Padé with `method=:classical` auto-routes to `:svd`.                        |
| CP.1.6  | Element-type dispatch defaults: F64 → `:classical`, BigFloat → `:svd`, kwarg overrides work.      |
| CP.1.7  | Geometric series `1/(1-z/2)` → rank-1 Toeplitz at high `m` → SingularException → SVD fallback.    |
| CP.1.8  | Fail-fast on negative `m` and unknown `method`.                                                   |
| CP.1.9  | Complex{Float64} classical Padé sanity.                                                           |

Mutation-proof: four mutations (sign-flip on RHS, drop singular
fallback, flip F64 default, drop off-diagonal route) all bite as
predicted.  Full procedure + counts in `docs/worklog/021-classical-pade-default.md`
§"Mutation-proof".

## Tightened test acceptance

The `:classical` default also tightens two long-range path-network
test bounds (acceptance d in `padetaylor-txg`):

| test  | before    | after    | ratio    |
|-------|-----------|----------|----------|
| PN.2.2 z=30 F64  | rtol `1e-9`   | rtol `1e-12`  | 1,000×    |
| PN.2.3 z=10⁴ F64 | rtol `5e-5`   | rtol `5e-10`  | 100,000×  |

Both are direct consequences of the dispatch change.  Mutation P3
(reverting the F64 default to `:svd`) bites both as expected.

## References

  - **FW 2011 §5.1.4** — eqs. (5.4), (5.5), line 346 (singular
    fallback), line 350 ("method of choice ... Toeplitz approach
    and ... backslash operator").  `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:330–350`.
  - **GGT 2013 §2 + Algorithm 2** — `references/markdown/GGT2013_robust_pade_via_SVD_SIAMRev55/GGT2013_robust_pade_via_SVD_SIAMRev55.md`.
  - **ADR-0001** — four-layer architecture.
  - **ADR-0002** — bigfloat-SVD via GenericLinearAlgebra; explains
    why `:svd` is load-bearing at `BigFloat`.
  - **Worklog 020** — empirical diagnosis and FW 2011 close-read.
  - **Worklog 021** — implementation log + mutation-proof.
  - `src/RobustPade.jl` — `classical_pade_diagonal` + dispatch.
  - `test/classical_pade_test.jl` — CP.1.* tests.
