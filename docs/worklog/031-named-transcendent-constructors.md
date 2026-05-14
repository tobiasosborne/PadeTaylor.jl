# Worklog 031 — Named Painlevé transcendent constructors

**Date**: 2026-05-14 (continues worklog 030)
**Author**: Claude Opus
**Beads**: `padetaylor-c86` — closed.  `padetaylor-icf` filed (P3,
deferred closed-form families).
**Scope**: Thread 2 of the Painlevé-infrastructure roadmap — the
named-transcendent constructors ADR-0006 deferred.

> **Take-home**: `tritronquee(:I)` and `hastings_mcleod()` build the
> `PainleveProblem` for those two literature-pinned solutions, IC and
> `name` tag baked in.  `solve_pade(tritronquee(:I))` now returns a
> `PainleveSolution` whose `solutionname` is `:tritronquee`.  Only two
> named transcendents ship in v1, and the reconnaissance below records
> exactly why.  Test suite **1590 → 1630 GREEN**; ADR-0008 accepted.

## Ground truth read first — the reconnaissance

A Sonnet subagent surveyed every in-tree reference for named special
solutions.  The finding that shaped the whole bead: named Painlevé
solutions come in three structurally different kinds, and only one is
cleanly shippable as a constructor.

  1. **Point-IC named solutions** — the paper gives `(u, u')` at a
     point to high precision.  **Exactly two qualify with 16-digit
     data:**
       - **PI tritronquée**: `u(0) ≈ -0.1875543083404949`,
         `u'(0) ≈ 0.3049055602612289` —
         `references/markdown/FW2011_painleve_methodology_JCP230/
         FW2011_painleve_methodology_JCP230.md:224-229`, verified
         verbatim at coding time (FW's 32-digit Maple BVP over
         `[-20i, 20i]`, "accurate to better than `10⁻²⁰`").
       - **PII Hastings–McLeod** (`α = 0`):
         `(u(0), u'(0)) ≈ (±0.3670615515480784, ∓0.2953721054475501)`
         — `references/markdown/FW2014_second_PII_exploration_FoCM14/
         FW2014_second_PII_exploration_FoCM14.md:252-258`, verified
         verbatim.  Two sign-symmetric copies that "only differ in
         sign".
  2. **Asymptotic-BC named solutions** — general-`α` Hastings–McLeod,
     the secondary HM, Ablowitz–Segur, the Boutroux `k = 0` family.
     Defined by `z → ±∞` behaviour; turning that into a point IC is a
     *connection problem* (a BVP solve), not a constructor.  Out of
     scope.
  3. **Closed-form parametrised families** — PII rational `uₙ`
     (`u₁ = -1/z`, `u₂`, `u₃` explicit), PII Airy solutions for
     `α = n + ½`, PIV entire solutions `-2z` / `-⅔z`.  Genuinely
     shippable, but a different API shape (parametrised families) with
     *exact closed-form oracles* — deferred to bead `padetaylor-icf`.

Per CLAUDE.md Law 1 the two v1 citations were re-read directly from
the markdown — the subagent's quotes matched exactly, but a literature
subagent's output is verified, not trusted (CLAUDE.md Rule 3).

## What shipped

### `src/PainleveNamed.jl` (new) — bead `padetaylor-c86`

Two named constructors, the third file of `module Painleve`:

  - `tritronquee(equation = :I; zspan, order)` — the PI tritronquée,
    FW 2011 §4.1 IC.
  - `hastings_mcleod(; branch, zspan, order)` — the PII `α = 0`
    Hastings–McLeod solution; `branch ∈ (:positive, :negative)`
    selects the sign-symmetric copy.

Each returns an ordinary `PainleveProblem` built via the keyword
constructor, then stamped with `name` (`:tritronquee` /
`:hastings_mcleod`) through the `_with_name` helper.  The caller still
chooses `zspan` and `order`; the IC is not the caller's to choose — it
*is* the named solution.

### `src/Painleve.jl` — the `name` field

`PainleveProblem` gained a final field
`name::Union{Symbol,Nothing}` — the named-solution tag, or `nothing`
for a problem built from explicit ICs.  The plain
`PainleveProblem(equation; …)` constructor and all four `_build_*`
helpers set `name = nothing`; the named constructors set it; and
`_painleve_solution` (in `PainleveSolution.jl`) now defaults its
`name` argument to `pp.name`, so the tag rides through the forwarding
methods into `PainleveSolution.name`.  ADR-0007 had already given
`_painleve_solution` a `name` parameter for exactly this — this bead
wired its source.

### `docs/adr/0008-named-transcendent-constructors.md` (new, Accepted)

Records the v1 scope, the three-kinds taxonomy, the `name`-field
decision, the fail-loud guards, and four rejected alternatives (a
value-function; bypassing the solver; defaulting the IC point away
from the user; storing `name` only on `PainleveSolution`).

### Tests — 1590 → 1630 GREEN

`test/painleve_named_test.jl` (new) — `NT.1`–`NT.3`:

  - **NT.1.*** — the built `PainleveProblem` carries the cited IC,
    asserted *verbatim* against the reference values (the IC is the
    oracle — CLAUDE.md Rule 5), plus the right equation / params /
    frame / `name`.
  - **NT.2.1** — `name` propagates through `solve_pade` /
    `path_network_solve` into `solutionname(sol)`.
  - **NT.2.2** — the Hastings–McLeod `:positive` and `:negative`
    branch solutions are exact negatives of each other.  PII at
    `α = 0` is `u'' = 2u³ + zu`, invariant under `u → -u`; this is a
    genuine cross-check that the branch sign is wired through *both*
    the IC and the integration, not a tautology.
  - **NT.2.3** — the HM `:positive` branch decays on `x > 0`
    (`u ~ Ai(x) → 0`), consistent with its defining pole-free
    non-oscillatory character.
  - **NT.3.1** — the three fail-loud guards.

Mutation-proof: three mutations (perturb the tritronquée IC's last
digit; revert the `name` plumbing; drop the HM branch sign flip), each
confirmed RED and reverted.  Procedures in the test file footer.

## Frictions

  - **The bead is small; the *scoping* was the work.**  Writing two
    constructors is ~40 LOC.  Deciding that it is *only* two —
    distinguishing a 16-digit point IC from an asymptotic BC that
    merely looks like one in the paper's prose — needed the full
    reference survey.  The honest deliverable here is as much the
    three-kinds taxonomy (ADR-0008, worklog) as the code.

  - **Adding a struct field touches every constructor.**  `name` on
    `PainleveProblem` meant updating the four `_build_*` positional
    constructions.  Small, but a reminder that a 7th field is not
    free — justified because `name` is genuinely part of a named
    problem's identity, and defaults cleanly to `nothing` otherwise.

  - **The IC values were already in-tree as folklore.**  The
    tritronquée IC appears as a bare literal in
    `phase9_tritronquee_test.jl`, `fw_fig_41_test.jl`, and
    `lattice_dispatcher_test.jl`.  It is now a *cited, tested*
    constant in `src/PainleveNamed.jl` — but the existing literals
    were deliberately left alone for this bead (retrofitting them to
    call `tritronquee(:I)` is a separate, optional cleanup, not part
    of shipping the constructor).

## What remains genuinely out of scope

  - **Closed-form named families** (`padetaylor-icf`, P3) — PII
    rational / Airy, PIV entire.  Parametrised families with exact
    oracles; their own ADR.
  - **Asymptotic-BC named solutions** — general-`α` Hastings–McLeod,
    secondary HM, Ablowitz–Segur, Boutroux `k = 0`.  A constructor
    cannot produce these from a point IC; the right shape is a
    `hastings_mcleod(α)` that runs a connection-problem BVP
    internally — a meaningfully larger piece of work, not filed yet.
  - **PIV Generalized Hermite / Okamoto rationals** — the formulas
    are not in-tree (the reference cites Clarkson 2006, absent).

## Hard-won lesson

**A literature subagent's job is to *sort*, not just to *find*.**  The
useful output of the reconnaissance was not "here are the named
solutions" — it was "here are the three *kinds*, and here is which
kind each one is".  The taxonomy is what made the v1 scope obvious and
defensible: ship the point-IC kind, defer the closed-form kind, reject
the asymptotic-BC kind.  When delegating literature work, ask the
subagent for the *classification*, not just the catalogue — the
classification is the decision.
