# Worklog 072 — corpus-v2 pole-zoo: orchestrated serial build + the extract_poles finding

**Session 2026-06-08.** Designed and delivered a second ground-truth corpus epic
(`padetaylor-25og`) extending `padetaylor-p1v0` into the **pole-morphology** blind
spots the v1 corpus left, then — demonstrating the new path-network coverage on a
live plot — surfaced a genuine `extract_poles` bug (`padetaylor-fzse`).

## The arc

1. **Suggestion → catalogue (read-only fan-out).** From a request for "many more
   exactly-solvable ODEs hitting blind spots", produced the blind-spot matrix
   (`{pole order} × {arrangement} × {residue known}`) and fanned out **5 read-only
   derivation subagents** (sympy/mpmath, no Julia — Rule 7) to turn the menu into
   catalogue-grade records. Skepticism (Rule 3) paid: the fan-out caught **5 of my
   own brief errors** before any code (parabolic-cylinder has no real poles;
   `extract_poles` has no multiplicity field; algebraic-PVI `1/√z` branches at the
   PVI fixed singularity; the ℘ oracle needs dps≥80 for complex roots; …) plus a
   candidate library bug. Plan: `docs/test_corpus/02_corpus_extension_plan.md`.

2. **Filed the epic** — `padetaylor-25og` + 18 children (14 test files + 4
   follow-on), parented, priority-ordered.

3. **Orchestrated serial build.** One **Sonnet** recon pass produced a reusable
   BUILD CONTRACT (naming, header skeleton, `runtests.jl` idiom, mutation-proof
   convention, per-API signatures + minimal examples). Then **15 beads built
   serially, one Opus coding agent each** (Rule 7 — never two Julia processes), each
   doing oracle-capture → test → standalone-GREEN → mutation-proof → register, with
   the orchestrator reviewing every return, committing/pushing per bead, and
   raising/annotating beads on findings.

## What shipped

14 new corpus test files + the CVB.4 addition, ~370 assertions, all GREEN +
mutation-proven, `src/` untouched. **End-of-phase full `Pkg.test()`: 8859 pass /
10 broken / 0 fail** (18m32s). The +372 pass / +5 broken reconcile exactly to the
new work (the 5 new `@test_broken` markers on top of v1's 5).

- **1 new library bug** — `padetaylor-v1ub`: `solve_pade` silently bridges an
  essential singularity (`e^{1/z}`). CONFIRMED against the *real* solver, not just
  the mpmath emulation: relERR 2e-16 far from 0 → 1e+17 *and wrong sign* at dist
  0.02, finite throughout, no throw. Annotated on the bead with the full curve.
- **5 `@test_broken` auto-flip markers** for the known guard bugs — CBr.3 + CPN.7
  (61um, the latter the FIRST real-walk exercise), CFail.1 (v1ub), CBvx.4 scalar +
  vector (53tu). Each flips to "unexpectedly passes" the day its guard lands.
- **4 clean-verdict gap closures (no bug):** shared-Q keystone (#5, Jacobi-triple
  external oracle), D1-coupled BVP Jacobian (#6, orthopoly + Kummer), PVI ζ-RHS
  wiring (#9, external algebraic solution), coupling vector-BC τ-assembly (#10).
- **5 brief errata** (`ERRATA.md`, corpus-extension section) — recipe slips caught
  by independent verification, **package correct every time**. No wrong computed
  value anywhere; the numerical core stays clean.

Per maintainer decision the 3 remaining children were spun out standalone:
`v1ub` (fix — design-led, must not misfire on genuine pole-bridging), `90oh`
(extract_poles multiplicity API + ADR), `vimm` (blocked on reference acquisition).

## The friction — the `extract_poles` finding (`padetaylor-fzse`)

To show the new path-network coverage I plotted `|u(z)|` across ℂ for the logistic
(vertical pole *row*) and the equianharmonic ℘ (2D pole *lattice*) —
`figures/demo_lattice_singularities.jl`. The **field is computed flawlessly** (all
24 641 cells finite; bright spots exactly on the analytic poles). But overlaying
`extract_poles` exposed a real bug: some poles got 3–4 crosses, others none.

Probed against the exact lattice (FW half-period Ω=1.363, spacing ≈2.73):

> **default** (cluster_atol=0.1, min_support=3, radius_t=5) → 20 reps, **12 with
> nearest-neighbour < 0.6** (over-split duplicates, min-nn 0.109 *just* above
> cluster_atol). **tuned** cluster_atol=0.4 → 17 reps, 0 duplicates — but only
> **13 on true poles, 4 false positives** (0.93–1.24 from any true pole) and **4
> missed**. ~76 % precision / recall; tuning only trades duplicate-clumps for
> spurious+missed.

**Root cause.** (1) The greedy single-pass first-fit clustering
(`src/PoleField.jl:213`) cannot self-merge — a pole whose cross-node root estimates
spread wider than `cluster_atol` fragments into several reps. (2) `radius_t=5`
admits *far* nodes whose (15,15) Padé places the pole poorly (~2.5 away in z),
inflating that spread AND planting **ghost denominator-roots** ~1 away from real
poles; several far nodes agree on a ghost, so it clears `min_support=3` and the
`1e-8` Froissart gate → a coherent false positive.

**Lesson / why the corpus missed it.** `polefield_test` (and the new CPN cases)
pin a *single trajectory* (`min_support=1`, few close nodes) or a hand-placed grid;
the **dense generic 2D path-network with stock defaults** is a distinct, untested
regime. The reported pole *locations* are accurate — it is the **precision/recall
of the map** that fails. Fix is algorithmic (agglomerative/self-merging clustering
+ an agreement/validity gate + tighter far-root admission), not a knob; the bead
carries the fix menu and asks for a dense-elliptic-lattice corpus test asserting
~one rep per true pole at default params. Both `extract_poles` callers share
`_extract_poles_core`, so the fix lands for the scalar and vector paths together,
and `90oh`'s multiplicity work slots in alongside.

## Discipline notes

- One Julia process at a time across ~20 subagents (Rule 7), enforced by serial
  dispatch; full `Pkg.test()` only at end-of-phase. Mutation-proved every
  load-bearing assertion (Rule 4); `src/` byte-clean after each (M-impl restored).
- Per-bead commit/push cadence; the build-status table + ERRATA kept in lockstep
  (Law 2). `bd remember corpus-v2-expected-broken-count` records that a nonzero
  broken count is EXPECTED (the auto-flip markers), so a future agent doesn't read
  it as a regression — investigate only Fails.
