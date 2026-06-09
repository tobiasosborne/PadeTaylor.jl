# ADR-0030 — Formal-methods / verified-numerics frontier: DEFER (with a forcing condition)

**Status**: **accepted (decision: DEFER)** — 2026-06-09. Resolves the Tier-3
decision bead `padetaylor-krgy.13` (epic `padetaylor-krgy`, the test-hardening
sweep). Revisit only under the forcing condition in §Decision.
**Date**: 2026-06-09
**Beads**: `padetaylor-krgy.13`; methodology Tier-3 of
`docs/test_corpus/03_hardening_methodology.md` (the "frontier, flagged
aspirational" tier).

## Context

The test-hardening sweep's methodology synthesis (`03_hardening_methodology.md`,
Tier-3) surveyed machine-checked formal verification of numerical code:
proof assistants (Lean 4, Coq) for numerical kernels; SMT-backed verification;
proof-carrying / certified code; floating-point rounding proofs (Gappa / Flocq,
Melquiond). The research was explicit that **no Tier-3 claim reached a 3-0
verification vote** — these are high-assurance but high-cost, and the literature
support for them as a *practical* gate for a small team is weaker than for the
Tier-1/Tier-2 techniques. This ADR records the conscious decision rather than
letting the option lapse silently (Rule 9: a deferred corner is documented with
the exact condition that would force the work).

The relevant context that has CHANGED since the bead was filed: the sweep has now
shipped a deep, multi-technique assurance stack —

- **certified ball oracles** (ADR-0029, `krgy.5`): Arb midpoint-radius arithmetic
  gives *rigorous, machine-checked enclosures* of the closed-form oracle values
  (radius ~1e-74). This already delivers the single highest-value piece of
  "verified numerics": a proof-carrying golden the package is pinned against;
- **property-based** (`krgy.1`), **metamorphic** (`krgy.3`), **differential**
  (`krgy.14`), **convergence-order / MMS** (`krgy.6`), **automated mutation**
  (`krgy.4`), **accuracy-regression** (`krgy.8`), **static analysis** (`krgy.2`),
  and an **enforced provenance contract** (`krgy.9`).

Together these constrain the numerical core from many independent directions at a
maintenance cost a two-person, no-CI project can actually carry.

## Decision

**DEFER full formal verification of PadeTaylor's numerical kernels.** Do not
attempt Lean 4 / Coq proofs of the SVD, the GGT-2013 robust-Padé construction,
the Taylor recurrences, or the stepper — and do not adopt SMT-backed or
proof-carrying-code infrastructure — at this time.

**Rationale (the senior-engineer call):**

1. **Marginal assurance is small; marginal cost is large.** The contested
   correctness questions for this package have been edge-guard / fail-loud bugs
   (`xhjw`, `jznu`, the corpus's `q0yq`/`53tu`/`61um`/`v1ub`) and degenerate-regime
   divergences — NOT wrong core arithmetic. Across the whole sweep, no computed
   value was wrong; the numerical core is clean under property, metamorphic,
   differential, certified-ball, and convergence-order testing. A Lean 4 proof of
   the SVD-Padé kernel would defend a perimeter that is not where the bugs live,
   while imposing an ongoing proof-maintenance tax on every refactor.
2. **The highest-value "verified numerics" is already in hand** — certified ball
   oracles (ADR-0029) give rigorous error bounds where they matter most (the
   goldens), without proving the imperative kernel.
3. **Cost/maintainability for a small no-CI team.** Mechanized proofs of
   floating-point and linear-algebra kernels are a specialist, full-time effort
   (cf. the Flocq/Gappa and CompCert-FP bodies of work). They do not fit a
   discipline-heavy but small team, and — per Rule 11 — there is no CI to keep a
   proof corpus from rotting.

## Forcing condition (when to revisit — Rule 9)

Re-open this decision if ANY of the following becomes true:

- a **published headline claim** about a PadeTaylor result requires a
  machine-checked proof (e.g. a claimed pole-free sector or an exact identity
  that a referee/collaborator wants certified beyond testing); OR
- a specific numerical **kernel's correctness is contested** and is genuinely
  *unprovable by the existing testing stack* (no oracle, no metamorphic relation,
  no manufactured solution can pin it) — i.e. a true oracle-problem residue that
  only a proof can close; OR
- the project gains **sustained formal-methods capacity** (a contributor who owns
  Lean 4 / Flocq) for whom the maintenance cost is no longer prohibitive.

## The sanctioned cheapest entry, IF ever pursued

Should the forcing condition trigger, the lowest-cost foothold — NOT adopted now,
recorded so a future pilot starts in the right place — is a **single Lean 4 lemma
for one load-bearing algebraic identity**, not a kernel proof. Candidates: the
diagonal-`(m,m)` Padé order law `2(order÷2)+1` (the `krgy.6` finding); a
Bäcklund / affine-Weyl group relation already exercised algebraically as a
metamorphic relation (`NYS.1.6`); or a Taylor-recurrence closed form. The repo
already carries Lean 4 tooling (the `lean4` skills), so a narrow pilot is
mechanically cheap to start — the cost is the *ongoing* proof maintenance, which
is exactly what this ADR defers.

## Consequences

- `padetaylor-krgy.13` closed as **deferred-with-condition** (not abandoned).
- No Lean 4 / Coq / SMT dependency or proof corpus is added to the project.
- The forcing condition above is the trip-wire; until it fires, assurance rests
  on the shipped Tier-1/Tier-2 stack + certified ball oracles.
- TaylorModels.jl as a "turnkey verified integrator" remains REFUTED (0-3,
  `03_hardening_methodology.md`) and is explicitly not the back door to this tier.

## Citations

- `docs/test_corpus/03_hardening_methodology.md` — Tier-3 framing; the refuted
  Taylor-model claim; the certified-oracle (Tier-2 §2.2) that subsumes the
  highest-value verified-numerics need.
- `docs/adr/0029-certified-ball-oracles.md` — the rigorous-enclosure capability
  already shipped.
- Melquiond, *Gappa* / Boldo–Melquiond, *Flocq* — the reference point for
  machine-checked floating-point proofs (the cost basis for this deferral).
