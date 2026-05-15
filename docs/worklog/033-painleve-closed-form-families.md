# Worklog 033 — Closed-form Painlevé families

**Date**: 2026-05-15 (continues worklog 031)
**Author**: Claude Opus
**Beads**: `padetaylor-icf` — closed.
**Scope**: Ship the parametrised-family-with-exact-oracle kind of named
Painlevé solutions, deferred from worklog 031 / ADR-0008.

> **Take-home**: `pii_rational(n)`, `pii_airy(n; θ)`, `piv_entire(kind)`
> ship the seven new closed-form Painlevé constructors in v1 — three PII
> rational `u_n` for `n ∈ {1,2,3}`, two PII Airy `u_{n+½}` for `n ∈ {0,1}`,
> two PIV entire `u = -2z` and `u = -(2/3)z`.  The IC for each is derived
> from the closed-form analytic formula at `zspan[1]`; the same formula is
> the *oracle* for the self-validating tests.  ADR-0010 accepted; test
> suite **1652 → 1708 GREEN** (+56 CF.* assertions); new dep
> `SpecialFunctions.jl` for the Airy `Ai/Bi`.

## Ground truth read first — three Sonnet subagents in parallel

  1. **FW 2014 §2.2 + §2.3** (`FW2014_*.md:114-167`): the PII rational
     formulas (eq. 6, md:119-120) and PII Airy formulas (eqs. 7-10,
     md:130-162).  Critical confirmation: FW 2014 §1 eq. 2 (md:47)
     uses the form `u'' = 2u³ + zu + α` — *matches* our existing
     `_pII_rhs(α)` factory verbatim, so no sign-flip drama.  The
     rational `u_n` solves PII at `α = n`; the Airy `u_{n+½}` at
     `α = n + 1/2`.  FW 2014 gives the canonical IC for `u_3`
     (`(u_3(0), u_3'(0)) = (0, 0)`, md:218) — that drives the
     `pii_rational(3)` default `zspan = (0.0, 5.0)`.

  2. **RF 2014 §1-2** (`ReegerFornberg2014_*.md:60-130`): the two PIV
     entire solutions.  Critical gap: RF 2014 md:91 *names* `u = -2z`
     and `u = -(2/3)z` but DOES NOT state the `(α, β)` parameter pairs
     in-line.  The Noumi-Yamada reference [44] has them; we don't have
     that paper in-tree.  **Derived algebraically here**: substituting
     `u = -2z` into the PIV equation and matching coefficients gives
     `α = 0, β = -2`; substituting `u = -(2/3)z` gives `α = 0, β = -2/9`.
     Derivation recorded in `_u_piv_entire`'s docstring.

  3. **Audit of `PainleveNamed.jl` + `Painleve.jl` + ADR-0008**: the
     `_with_name(pp, sym)` helper rebuilds the `PainleveProblem` with a
     new name tag; works for `pii_rational` and `piv_entire` (params
     already correct).  For `pii_airy` the `θ` solution-family
     parameter must extend `params` beyond the default `(; α)`, so the
     constructor uses the inner `PainleveProblem(...)` form directly.
     The existing closed-form-oracle test idiom in
     `painleve_named_test.jl` (exact `==` on the IC tuple plus `≈` at
     downstream `z`) carries over unchanged.

The three agents ran in ~80 s wall-clock total in parallel.

## What shipped

### `src/PainleveClosedForm.jl` (new) — 210 lines

Three constructors + three closed-form oracle helpers.  Slotted as a
sibling of `PainleveNamed.jl` and `included` from `src/Painleve.jl`
after it.

The PII Airy implementation uses the `Φ' = -z/2 - Φ²` identity
(FW 2014 md:154) to fold Airy-function derivatives through.  This is
the load-bearing identity: it lets `u_{n+½}'(z)` be expressed as a
rational function of `Φ` and `z` (no additional Airy evals needed
beyond `Φ` itself).  For `u_{1/2}`: `u' = -Φ' = z/2 + Φ²`.  For
`u_{3/2}`: a quotient rule with `(6Φ² + z)Φ' + Φ` and `4ΦΦ' + 1` as
the numerator/denominator derivatives.

The PIV entire `(α, β)` derivation is recorded inline in
`_u_piv_entire`'s docstring — substituting `u = -2z` into the PIV RHS
and demanding the residual vanish identically gives a two-equation
system in `(α, β)` (the `1/z` and `z` terms must each vanish), which
yields `(0, -2)` uniquely.  Same for `u = -(2/3)z` → `(0, -2/9)`.

### `src/Painleve.jl` — three lines changed

  - `include("PainleveClosedForm.jl")` after `PainleveNamed.jl`.
  - `export pii_rational, pii_airy, piv_entire` alongside the
    existing `tritronquee, hastings_mcleod`.

### `Project.toml` — new dep

`SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"` in `[deps]`.
Justification: PII Airy solutions require evaluation of `Ai/Bi/Ai'/Bi'`
at the IC point.  No other code path in the package uses
SpecialFunctions; the dep is load-once, ~30 transitive packages via
`OpenSpecFun_jll`.

### `test/painleve_closed_form_test.jl` — 220 lines, 56 assertions

Four blocks (CF.1.*, CF.2.*, CF.3.*, CF.4.*).  The idiom:

  - **CF.X.1, CF.X.2, …**: assert `pp.problem.y0 == (u_exact, up_exact)`
    by exact equality (the constructor must stamp the formula's value).
  - **CF.X.k (the solver-cross-check)**: `solve_pade(pp)` then evaluate
    `sol(z_test)`; assert `≈` the closed form at `rtol ≈ 1e-10`.
  - **CF.4.***: fail-loud guards (out-of-range `n`, unknown `kind`,
    pole at `zspan[1]`).

For PIV entire the cross-check uses `atol = 1e-10` rather than rtol —
the closed form is linear, so the Padé-Taylor IVP accumulates only
~2e-12 of roundoff over 20 steps (linear-truncation theoretically zero;
the residual is float arithmetic noise).

### `docs/adr/0010-painleve-closed-form-families.md` (Accepted)

The decision record: context, the three families, the `name` /
`params` convention, the four rejected alternatives (richer wrapper
type; value-function; Bäcklund recurrence in v1; finite-difference
derivatives).  Out-of-scope section flags `u_{5/2}` and the
Generalized-Hermite/Okamoto PIV rationals for follow-up.

## Mutation-proof

Four mutations on `src/PainleveClosedForm.jl`, each verified RED
against the new CF.* tests then reverted.

  - **MC1** — sign flip on `_u_pii_rational(1)`: `-1/z` → `+1/z`.
    Bites 3 assertions: CF.1.1 (IC tuple mismatch) and CF.1.4 (solver
    cross-check at z=1.5 fails for both `u` and `u'`).
  - **MC2** — swap `(α, β)` for `piv_entire(:minus_2z)`: `(0, -2)` →
    `(0, -2/9)`.  Bites 3 assertions: CF.3.1 (`params` wrong) and
    CF.3.3 (solver cross-check explodes — wrong equation entirely).
  - **MC3** — drop the `zspan[1] == 0` pole guard in `pii_rational(1)`.
    Bites 1 assertion: CF.4.4 (the `@test_throws` no longer throws,
    and the call returns a PainleveProblem with `Inf` IC).
  - **MC4** — sign flip in `_airy_phi_prime`'s chain rule:
    `-_AIRY_C·(...)` → `+_AIRY_C·(...)`.  Bites 5 assertions across
    CF.2.* (every IC and cross-check assertion that consumes `Φ`
    is affected).

Each mutation applied, fast-loop tested via
`include("test/painleve_closed_form_test.jl")` (~3 s baseline,
~3 s under mutation), confirmed RED, reverted.  Procedure
recorded in the test file footer.

## Frictions

  - **`Irrational` vs `Float64` in `pp.params`.**  The first GREEN
    attempt failed CF.2.4: the test passed `θ = π` (Julia's
    `Irrational{:π}` constant) but expected `pp.params.θ ==
    Float64(π)`.  The constructor stored `θ` as-is, so the test's
    NamedTuple equality compared `Irrational{:π} == Float64`, which
    is false (different types, even if `==` is defined for the
    values).  Fix: `θ_stored = float(θ)` at the constructor entry
    to normalise any `Real` input to a concrete float type.
    Resurfaces a known Julia gotcha — `θ::Real` accepts
    `Irrational`, but `Irrational` doesn't compose with most
    Float64 operations cleanly.  Recorded as a small CLAUDE-style
    lesson below.

  - **RF 2014 doesn't state `(α, β)` for the entire solutions.**
    The agent flagged this immediately; the derivation by direct
    substitution into the PIV equation is straightforward (the
    residual must vanish *identically* in `z`, which gives two
    equations for `(α, β)`) but it's not trivial to spot without
    actually doing the algebra.  The derivation now lives in
    `_u_piv_entire`'s docstring with the algebraic steps, so a
    future reader doesn't have to redo it.  This is the right shape
    for "ground truth not literally in-tree but derivable from
    in-tree material" — record the derivation, not just the answer.

  - **PIV linear solution roundoff floor.**  The first GREEN attempt
    failed CF.3.3 at `atol = 1e-12` (actual error 2.5e-12 over 20
    steps).  Relaxed to `atol = 1e-10`.  Linear solutions of PIV
    have `u'' = 0` analytically, so the Padé-Taylor IVP truncates
    exactly at order 1 in pure algebra — but float arithmetic of
    `(u')²/(2u) + (3/2)u³ + …` evaluating to "exactly zero"
    introduces ~eps roundoff per step, accumulating linearly.
    Lesson: even for "exact" closed forms, the solver's roundoff
    floor over `N` steps is `~N · eps(T)`; pin tolerances
    accordingly.

## Hard-won lesson

**Closed-form derivations belong in the docstring, not in a separate
note.**  RF 2014's `(α, β)` derivation is six lines of algebra.
Without it, someone reading `_u_piv_entire`'s `(α, β) = (0, -2)` and
`(0, -2/9)` has to either trust the implementation or re-derive.  With
the derivation inline, the math is self-contained and the reviewer
can verify by reading.  The same principle applies to the Bäcklund
recurrence we deferred — when we ship `u_n` for `n ≥ 4`, the
recurrence relation should be derived inline, not cited externally.

This generalises Rule 1 (fail loud) to *documentation*: do not silently
require the reader to consult an external paper for the proof that an
in-tree formula is correct; provide the proof.  Especially when the
external paper itself "cites [44]" rather than stating the value
directly.

## Total scope

  - 1 new source file (`PainleveClosedForm.jl`, 210 lines).
  - 3 lines added to `src/Painleve.jl` (include + 1 export line).
  - 1 new dep in `Project.toml` (`SpecialFunctions`).
  - 1 new test file (`painleve_closed_form_test.jl`, 220 lines, 56
    assertions).
  - 1 new ADR (`docs/adr/0010-painleve-closed-form-families.md`).
  - 1 worklog (this file).
  - Test suite: 1652 → 1708 GREEN.

The `pii_airy(2)` (i.e. `u_{5/2}`) follow-up could be filed as a P3
bead alongside the v2 work on Bäcklund-recurrence-derived higher-`n`
rationals — both are extensions of this ADR rather than new ADRs.
Not filed in this session; the v1 scope is complete and shippable
as-is.
