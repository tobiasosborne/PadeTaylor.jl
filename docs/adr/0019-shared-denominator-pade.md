# ADR-0019 — Shared-denominator (type-II Hermite–Padé) robust Padé

**Status**: Accepted (2026-05-20) | **Bead**: `padetaylor-0ln` epic
(`.2` V1a impl + tests, `.3/.4/.5` oracles V1b/c/d, `.6` V1e triple-oracle
suite + this ADR) | **Findings**: `docs/v0p2_pillarA_hermite_pade_findings.md`.

## Decision

`SharedPade.shared_denominator_pade(jets, m; tol)` is the PadeTaylor.jl
v0.2 keystone: it computes a **shared-denominator** (simultaneous /
type-II Hermite–Padé) approximant for a `d`-component vector of formal
power series — `d` numerator polynomials `P₁,…,P_d` over **one** shared
denominator `Q`, all normalised to `Q(0) = b[1] = 1`.

`Q` is recovered as the smallest-singular-value right null vector of a
**stacked block-Toeplitz matrix** `A_full` (`dm × (m+1)`), built by
vertically stacking one `m × (m+1)` Toeplitz block per component, then
refined by the column-reweighted QR step ported unchanged from Chebfun's
`padeapprox.m`. Each numerator is recovered from an upper Toeplitz block.
Four numerical-breakdown modes **throw** rather than return a NaN/zero
lie (Rule 1). The algorithm is **additive**: one new module
(`src/SharedPade.jl`, ≤ 200 LOC), no edit to the scalar `RobustPade`
path. Validation is **three independent oracles** wired into
`test/shared_pade_test.jl` (V1e).

## Context

A vector meromorphic ODE solution `f = (f₁,…,f_d)` has all `d` components
blowing up at the **same** movable poles — the singularities are a
property of the flow, not of the individual component
(`docs/v0p2_pillarA_hermite_pade_findings.md:46–68`, §2). Fitting each
component independently (e.g. with scalar `robust_pade`, or AAA per
component) yields `d` *different* pole sets that agree only up to noise.
The correct object is therefore one shared `Q` — a type-II Hermite–Padé
approximant in the sense of Mano–Tsuda 2017 eq. (1.2): `P^(0)` plays the
role of `Q`, and `P^(j)/P^(0) ≈ f_j/f₀`. The shared `Q` is what makes the
recovered pole set a single, consistent estimate of the system's
singularities — the prerequisite for lifting the PadeTaylor pole-tracking
pipeline from scalar to vector ODEs.

The scalar v0.1 path (`RobustPade.robust_pade`, `method=:svd`) already
solves the scalar matching condition `P(z) = f(z)Q(z) + O(z^{m+n+1})` as
a Toeplitz null-space problem (GGT 2013 Algorithm 2, ADR-0001 /
ADR-0005). v0.2 needs the vector generalisation without disturbing that
path.

## Architecture

### The stacked block-Toeplitz SVD

Mano–Tsuda 2017 §2.2 eq. (2.6)
(`docs/v0p2_pillarA_hermite_pade_findings.md:69–115`) generalises the
scalar GGT system directly. For `d` components, stack one `m × (m+1)`
Toeplitz block per component vertically:

```
        ┌ block₁ ┐   block_i[r,c] = c⁽ⁱ⁾_{m+r-c}
A_full =│   ⋮    │   from jets[i]
        └ block_d ┘
```

`A_full` is `dm × (m+1)`. Its smallest-σ right null vector is the shared
`Q`: every block contributes `m` matching equations against the *same*
unknown coefficient vector `b₀,…,bₘ` — exactly the constraint "all
components share `Q`". **One SVD of one rectangular matrix**; no coupled
or iterated solves (findings §3). The numerator step `aᵢ = upper · b`
uses the `(m+1)×(m+1)` upper Toeplitz block of each jet.

### `d = 1` reduces to GGT — the primary correctness oracle

For `d = 1`, `A_full` is a single `m × (m+1)` block, *identical* to the
scalar GGT `C̃` at `n = m`. The SVD, the QR-reweighting, the numerator
recovery and the `b[1]=1` normalisation are the verbatim scalar
operations. Hence `shared_denominator_pade([jet], m)` reproduces
`robust_pade(jet, m, m; method=:svd)` bit-for-bit modulo SVD sign
conventions. This reduction is the cheapest, tightest internal check and
is asserted directly in test SP.1.1.

### The QR-reweighting (ported unchanged from `padeapprox.m`)

After the SVD picks `b₀ = Vt[end,:]`, Chebfun's `padeapprox.m` lines
278–280 — *beyond* GGT 2013 Algorithm 2 itself — refine it:
`D = diag(|b₀| + √ε)`, `QR((A_full·D)ᵀ)`, `b = D·Q[:,end]`. The
column-reweighting better preserves the genuine exact zeros of `b` for
blocks at accuracy near the tolerance (findings §4). It ports without
change to the stacked `A_full`. The V1e mutation-proof (mutation M4)
confirms it is genuinely load-bearing — dropping it leaves the 1e-12
Float64 results intact but breaks the bit-exact `Rational{BigInt}` and
BigFloat-1e-30 oracle agreements.

### Degree reduction and the four defensive throws

A `while` loop (the vector analogue of `RobustPade`'s reduction) shrinks
`m` and rebuilds `A_full` whenever the smallest-σ null space is not
isolated (`ρ < m`). Four breakdown modes throw with a `suggestion` /
`detail` message (findings §7 "Failure modes", Rule 1):

1. **jet too short** — `length(jets[i]) < m+1`: cannot fill the block.
2. **all-zero jets** — every σ below `τ = tol·‖c‖` (`ρ = 0`): the input
   is indistinguishable from zero at tolerance.
3. **Q(0) ≈ 0** — recovered `b[1]` below tol: a pole sits at the
   expansion centre, so the `Q(0)=1` normalisation would divide by ~0.
4. **non-isolated null space** — more than one σ at/below `τ`: the
   shared denominator is not unique (Mano–Tsuda p. 12, "unique iff
   rank = m").

### Arb-precision routing

`Matrix{BigFloat}` `A_full` routes through `LinAlg.pade_svd`, which
dispatches to `GenericLinearAlgebra.svd` (one-sided Jacobi) — the same
arb-precision route ADR-0002 established for the scalar path. Test SP.3.6
exercises this at 256-bit precision.

## The three-oracle validation strategy

`shared_denominator_pade` recovers `Q` via an SVD + a QR-as-null-space
refinement. A single oracle that also used SVD/QR could hide a sign-flip,
a Toeplitz off-by-one or a wrong-column bug behind a matching
implementation. v0.2 therefore validates against **three structurally
independent** reference implementations (findings §6), all wired into
`test/shared_pade_test.jl` by bead V1e (108 assertions, GREEN):

- **Oracle 1 — Calgo 766** (`external/probes/calgo766-oracle/`). ACM TOMS
  Algorithm 766 (Cabay–Jones–Labahn 1997), FORTRAN, a Beckermann–Labahn
  striped-Sylvester recurrence — a *different algorithm* from the
  block-Toeplitz SVD. Double-precision build; outputs pinned as Julia
  literals in `test/_oracle_shared_pade.jl` (no live `ccall`, the v0.1
  `_oracles.jl` convention). Calgo's per-series degree vector forces a
  spurious common factor, so the comparison is by **function values**
  (the factor cancels in the reduced ratio) and **pole subset**.

- **Oracle 2 — block-Toeplitz determinant**
  (`external/probes/shared-pade-determinant-oracle/`). Mano–Tsuda 2017
  Proposition 2.3 bordered block-Toeplitz determinant — `Q`'s coefficients
  as a list of signed maximal minors, with **no SVD and no
  QR-as-null-space anywhere**. On `Rational{BigInt}` input the minors are
  exact, so this is a *tight* oracle: test SP.3.3 asserts **bit-exact**
  equality, SP.3.6 asserts BigFloat agreement to 1e-30.

- **Oracle 3 — AAA** (`external/probes/aaa-oracle/`).
  Nakatsukasa–Sète–Trefethen 2018, barycentric fit + Loewner-matrix SVD,
  fitting from **sampled function values** not the Taylor jet. AAA fits
  each component independently with no shared-denominator constraint; the
  cross-validation is the *coincidence* of the `d` per-component pole sets
  with the roots of the shared `Q` (diagnostic-grade, 1e-9 — findings
  §6 "AAA as a benchmark").

The V1e mutation-proof (4 mutations to `src/SharedPade.jl`, recorded in
the `SPM` block of `test/shared_pade_test.jl`) confirms the suite catches
regressions: wrong-column and transposed-Toeplitz mutations bite all three
oracles (72 / 71 of 108 RED); a numerator off-by-one bites only the
numerator assertions (50 RED, AAA stays GREEN — the suite localises the
fault); dropping the QR-reweighting bites only the exact / extended-
precision determinant-oracle asserts (38 RED) — proving each oracle and
each precision tier earns its place.

## Alternatives considered, rejected

- **`d` independent scalar `robust_pade` calls.** Rejected: produces `d`
  *different* denominators, hence `d` inconsistent pole sets — defeats the
  whole purpose (findings §2). It is kept only as the AAA-style diagnostic
  cross-check, not the primitive.

- **Coupled / iterated solve for the shared `Q`.** Rejected: the stacked
  block-Toeplitz SVD solves the simultaneous system in **one** SVD
  (findings §3). An iterative scheme would add convergence-failure modes
  and conditioning questions for no accuracy gain.

- **Minimal Mano–Tsuda system (`n` rows per block, `m = n(L-1)`).**
  Rejected for the impl: the GGT square-block convention (`n = m` rows per
  block) makes the `d=1` case bit-identical to the scalar path, which is
  the primary internal oracle. The minimal `m × (m+1)` system *is* used —
  inside the determinant oracle (Oracle 2), which reduces the taller
  `A_full` to a square subsystem by fraction-free row selection.

- **A single SVD-based oracle.** Rejected: not independent — a wrong-column
  or Toeplitz-layout bug could hide behind a matching SVD oracle (the
  motivation for Oracle 2's SVD-free determinant route).

- **Live FORTRAN `ccall` in the test suite.** Rejected: the suite must not
  depend on gfortran being installed (Rule 11, "quality gates run
  locally"). Calgo outputs are pinned as literals, mirroring the v0.1
  `_oracles.jl` Octave-pinning convention.

## Consequences

- A new public primitive `shared_denominator_pade` and an exported module
  `SharedPade`; the v0.1 scalar `RobustPade` path is untouched (additive
  architecture — Rule 9 "senior-grade", no v2 corners).
- `test/shared_pade_test.jl` grows to 108 assertions and `include`s two
  live probe oracles (`oracle.jl`, `aaa.jl`) via `@__DIR__`-relative
  paths; it remains runnable standalone and under `runtests.jl` with no
  external (gfortran) dependency.
- `test/_oracle_shared_pade.jl` joins `test/_oracles.jl` as a pinned-oracle
  literals file.
- The shared-`Q` primitive is the foundation for the v0.2 vector
  pole-tracking pipeline; downstream vector-ODE stepping work builds on
  this module.

## References

- `docs/v0p2_pillarA_hermite_pade_findings.md` — the spec: §2 block
  construction (`:46–135`), §3 scalar→shared-Q generalisation (`:151–192`),
  §4 Froissart suppression (`:193–247`), §6 the Calgo 766 + AAA oracles
  (`:323–378`), §7 algorithm steps 1–8 + failure modes (`:379–475`).
- `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`
  §2.2 eq. (2.6) and Proposition 2.3 — the block-Toeplitz null-space
  system and its determinantal representation.
- `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` — the scalar
  Algorithm 2 the `d=1` case reduces to.
- `external/chebfun/padeapprox.m` lines 278–280 — the QR-reweighting port.
- `src/SharedPade.jl` — the implementation (literate, ≤ 200 LOC).
- `test/shared_pade_test.jl` — the V1e triple-oracle suite + mutation
  record; `test/_oracle_shared_pade.jl` — pinned Calgo literals.
- `external/probes/calgo766-oracle/`, `…/shared-pade-determinant-oracle/`,
  `…/aaa-oracle/` — the three oracle probes (each with a README).
- ADR-0001 (four-layer architecture), ADR-0002 (BigFloat SVD via
  GenericLinearAlgebra), ADR-0005 (classical-Padé default at Float64).
- CLAUDE.md Rules 1, 4, 5, 6, 9, 10, 11.
