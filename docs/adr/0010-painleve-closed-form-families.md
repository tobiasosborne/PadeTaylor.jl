# ADR-0010 — Parametrised closed-form Painlevé families

**Status**: Accepted (2026-05-15) | **Bead**: `padetaylor-icf`
**Context**: Continues the named-transcendent arc (ADR-0008 shipped the
point-IC kind: `tritronquee(:I)`, `hastings_mcleod()`).  ADR-0008
explicitly deferred the parametrised-family kind to this bead.  ADR-0007
(the `PainleveSolution` wrapper) and ADR-0006 (the `PainleveProblem`
builder) supply the surrounding namespace; this ADR records the API
shape for the new families.

## Context

The ADR-0008 reconnaissance (worklog 031) sorted named Painlevé
solutions into three kinds:

  1. **Point-IC named solutions** — a 16-digit literature `(u, u')` at
     a specific point.  Shipped: `tritronquee`, `hastings_mcleod`.
  2. **Asymptotic-BC named solutions** — defined by `z → ±∞` behaviour;
     turning that into a point IC is a connection problem, not a
     constructor.  Out of scope.
  3. **Closed-form parametrised families** — exact formulas:
     PII rational `u_n` (FW2014 eq. 6), PII Airy `u_{n+1/2}` (FW2014
     eqs. 8-10), PIV entire (RF2014 md:91).  Different in shape from
     kind 1: a *family* indexed by `n` (or `kind`), with an *exact
     closed form* as the IC source and as the test oracle.  Deferred
     from v1 because the parametrised-family API + closed-form-oracle
     test pattern deserves its own treatment.

This ADR cashes the deferral in.

## Decision

Three new constructors in `src/PainleveClosedForm.jl` (a third file of
`module Painleve`, sibling to `PainleveNamed.jl`):

  - `pii_rational(n::Integer; zspan, order) -> PainleveProblem`
    — for `n ∈ {1, 2, 3}` (FW2014 eq. 6).
  - `pii_airy(n::Integer; θ = 0.0, zspan, order) -> PainleveProblem`
    — for `n ∈ {0, 1}` (FW2014 eqs. 8-9; `u_{5/2}` deferred).
  - `piv_entire(kind::Symbol; zspan, order) -> PainleveProblem`
    — for `kind ∈ {:minus_2z, :minus_two_thirds_z}` (RF2014 md:91).

Each returns a `PainleveProblem` consistent with the existing pattern
(ADR-0008 "rejected alternative: bypass the solver" still applies —
the constructor's job is to *pose* the problem; the closed form is
used to derive the IC and to validate, not to substitute for the
solve).

The `name` tag is a single Symbol per family — `:pii_rational`,
`:pii_airy`, `:piv_entire` — with the discriminating parameter stored
in `pp.params`:

  - `pii_rational(n)`:  `params = (; α = n)` — `α = n` is the PII
    equation parameter for `u_n`, so the discriminator IS the
    equation's free parameter; no extra field needed.
  - `pii_airy(n; θ)`:   `params = (; α = n + 1//2, θ)` — `θ` is a
    *solution-family* parameter (not an equation parameter), stored
    in `params` for full provenance.
  - `piv_entire(kind)`: `params = (; α, β)` — the `(α, β)` pair IS
    the discriminator (different `kind`s have different equation
    parameters); `kind` is recoverable from `(α, β)`.

The closed-form oracle helpers are internal-prefixed:

  - `_u_pii_rational(n, z) -> (u, u')`
  - `_u_pii_airy(n, θ, z) -> (u, u')`
  - `_u_piv_entire(kind, z) -> (u, u')`

They are not exported but are accessible via
`PadeTaylor.Painleve._u_pii_rational(...)` for tests and any caller
that wants the analytic ground truth.

### Pole guards (CLAUDE.md Rule 1)

`pii_rational(1)` with `zspan[1] = 0` would give `u_1(0) = -1/0`,
silently emitting `Inf` into the IVP state — unacceptable.  The
constructor fails-loud at any `zspan[1]` coinciding with a known
real-axis pole of the family member.  `piv_entire(kind)` with
`zspan[1] = 0` would give `u(0) = 0`, where the PIV RHS's `β/u` term
is singular; same fail-loud guard.

### Test idiom

The closed form IS the oracle.  Two layers of assertion per family:

  1. **Exact-equality IC check**: `pp.problem.y0 == (u_n(zspan[1]),
     u_n'(zspan[1]))` — no tolerance.  The constructor is *correct
     iff* it stamps the formula's exact value.
  2. **Approximate solver cross-check**: solve from the IC; evaluate
     `sol(z_test)`; assert `≈` the closed form `u_n(z_test)` at
     `rtol ≈ 1e-10`.  Validates that the solver reproduces the
     trajectory the formula predicts.

For linear closed forms (PIV entire) the second check is
`atol = 1e-10` — `u'' = 0` makes the Padé-Taylor IVP exact to
roundoff accumulation (~2e-12 over 20 steps).

The new `SpecialFunctions` dependency (for Airy `Ai/Bi/Ai'/Bi'`) is
justified by `pii_airy`'s closed form requiring evaluation at the IC
point.  Other paths in the package do not touch Airy functions.

## Consequences

**Positive**:

  - Users reaching for "the PII rational solution for `α = 3`" or
    "the PIV entire solution `u = -2z`" no longer have to derive
    the IC from the literature themselves.
  - Self-validating tests: 56 new assertions verify the
    constructor IC matches the formula exactly AND the solver
    reproduces the trajectory.  Mutation-proven (4 mutations bite at
    least one assertion apiece, with cross-test cascades).
  - The named-transcendent surface now covers both ADR-0008 kinds —
    a meaningful step toward the package's "useful Painlevé
    infrastructure" goal.

**Negative**:

  - New runtime dep: `SpecialFunctions.jl` (~30 transitive deps via
    `OpenSpecFun_jll`).  Mitigated by being a load-once-at-package-
    init cost; the existing path-network/BVP machinery does not call
    Airy.
  - `pii_airy` v1 caps at `n ∈ {0, 1}` (i.e., `u_{1/2}` and
    `u_{3/2}`).  `u_{5/2}` is in-tree (FW2014 eq. 10) but its
    derivative formula in terms of `(Φ, z, Φ')` requires careful
    hand-derivation.  Deferred to follow-up; recorded in this ADR's
    "Out of scope" section.

## Alternatives considered

**A. Return a richer `ExactPainleveSolution` type with the closed
form embedded.**  Rejected — the closed form is an *oracle*, not a
production access path.  Users who want `u_n(z)` exactly can call
the `_u_pii_*` helpers directly; users who want the *numerical*
trajectory go through `solve_pade(constructor(...))`.  The two paths
serve different needs; bundling them confuses the API.

**B. Expose `pii_rational(n)` as a value-function `u = pii_rational(n,
z)` returning the analytic value.**  Rejected per ADR-0006 §"Why no
painleve() value-function": the package is about *solving* the
equations, not implementing closed forms.  Closed forms exist only
for special cases; the value-function API would have to refuse most
problems, which is worse than not offering it.

**C. Ship all of FW2014 eq. 4's Bäcklund-recurrence-derived higher
rational solutions (`u_n` for `n ≥ 4`).**  Rejected for v1 — Bäcklund
recurrence is a separate API surface (it acts ON existing solutions),
deserves its own ADR.  Filed as a follow-up bead idea (not yet a
formal bead).

**D. Compute `pii_airy(n)` derivatives via finite-difference instead
of the closed-form `Φ' = -z/2 - Φ²` identity.**  Rejected — defeats
the "closed-form oracle" purpose; finite-difference introduces
roundoff that we are explicitly trying to *test against*.

## Out of scope (deferred)

  - **PII Airy `u_{5/2}`** — FW2014 eq. 10.  Hand-derivation of the
    `dN/dz`, `dD/dz` rational-in-Φ expressions is mechanical but
    error-prone; deferred to a follow-up bead.
  - **PII rational `u_n` for `n ≥ 4`** — via Bäcklund recurrence
    (FW2014 eq. 4).  Different API shape (recurrence on existing
    solutions); separate work item.
  - **PIV Generalized Hermite / Okamoto rationals** — RF2014 cites
    these (Noumi & Yamada ref [44]) but does not give the formulas
    in-line.  Out of scope per ADR-0008.
  - **PII rational complex-axis behaviour** — `u_3` has a pole near
    `z ≈ 1.508` on the real axis; the closed-form helpers handle
    this correctly (the formula doesn't care), but the IVP solver
    integrating across the pole needs the Padé bridge.  Covered by
    the existing `solve_pade` machinery; tests evaluate just below
    the pole to keep the IVP straightforward.
