# Worklog 076 — v0.2 vector-epic closeout + the v1ub out-of-class guard

**Date**: 2026-06-19 (session after the P1 Heun epic, worklog/ HANDOFF "LATEST").
**Beads**: closed `padetaylor-d4a`, `padetaylor-zcz`, `padetaylor-0ln.40`
(vector-epic critical path) and `padetaylor-v1ub` (the only confirmed
silent-wrong-answer bug). ADR-0033 added. Suite **9333 pass / 2 broken / 0 fail**.

## Orchestration shape (Rule-7-safe)

Read-only Sonnet investigators in parallel (no `julia`) → orchestrator
synthesis + the load-bearing design decision → ONE serial Opus coding agent
(owns the single `julia` process) → orchestrator independent verification
(full `Pkg.test()`) + doc lockstep + commit/push. The orchestrator made the
algorithmic call (which discriminator) rather than delegating it.

## Part 1 — the v0.2 vector epic was secretly DONE (reconcile-before-build)

`bd ready`/in-progress flagged the P1 vector epic `0ln` (85 %), with `0ln.40`
("port the full FW path-network driver to the vector solver") IN_PROGRESS and
its child `d4a` (re-render the P_I^(2) headline figure + measure coverage)
also in-flight. Four untracked diagnostic probe dirs under `external/probes/`
mapped 1:1 onto `d4a`'s closed S-beads — the tell that recent work happened
here and never landed cleanly.

Three read-only scouts reconciled the beads against ground truth (worklog 060,
ADR-0026 Amendment 10, `figures/_kkg_pi2_surface_helpers.jl:1061`, the
2026-06-19 full-`Pkg.test`-green HEAD `9fc5451`). Finding: **d4a and zcz were
substantively complete but never closed** — the figure is wired
(`on_target_failure=:skip` + dense disc-spaced targets), rendered, certified
coverage measured (~29 % at order 48, saturated: order 24→36→48→72 =
13→23→29→29 %), VC suite GREEN (`942127e`), full `Pkg.test` GREEN. The cap is
the no-extrapolation honesty contract, not a walk defect (ADR-0026 Amd 9).

Closeout (commit `cb15ba0`): finalized the ADR-0026 status header
("Accepted — implemented & closed 2026-05-26"); gitignored the four diagnostic
scratch probes per ADR-0026:199's own "reproducible, gitignored" intent
(non-destructive — kept locally for the deferred `xw3`); closed d4a/zcz/0ln.40.

**Operational lesson (bd export).** `bd export` writes to **stdout** by default;
piping it never updates the tracked file, and the auto-export-on-close hook left
`.beads/issues.jsonl` PARTIAL (only `d4a` reached it; `zcz`/`0ln.40` lagged)
while `bd show` correctly showed all three CLOSED. `bd dolt commit` is not a real
subcommand (no-op). Reliable resync = **`bd export -o .beads/issues.jsonl`**.
The `bd-export-lag-working-set` memory was corrected to this.

## Part 2 — v1ub: fail loud on out-of-class input (ADR-0033)

`solve_pade` is meromorphic-only but had no guard: on `u'' = u(1+2z)/z⁴`
(exact `u = e^{1/z}`, essential singularity at z=0) it marched toward 0 and
returned finite, plausible, **wrong** values (relerr 4e-10 @|z|=0.1 →
**1.2e17 AND wrong sign @0.02**) with no throw — the only confirmed
silent-wrong-answer bug, a direct Rule-1 violation. Confirmed in
`test/corpus_out_of_class_test.jl` (CFail.1d was a `@test_broken` auto-flip
marker awaiting exactly this fix).

### The design call (orchestrator-locked) — Padé *convergence*, not *radius*

The tempting guard — watch the Jorba–Zou radius / coefficient growth and trip
when it collapses — **is wrong and would break pole bridging** (the package's
headline feature). The radius of convergence equals distance-to-nearest-
singularity for **both** a pole and an essential point, so it collapses
identically; it cannot tell "a pole I can bridge" from "an essential point I
cannot." The recon's first proposal (track h_JZ collapse) was overridden on
exactly this ground.

The sound discriminator is **Padé sequence convergence** (GGT 2013 §8;
de Montessus de Ballore / Nuttall–Pommerenke): meromorphic ⇒ successive-order
Padé approximants agree (*including while bridging a pole*); essential ⇒ no
finite-degree rational converges, so they disagree and the disagreement grows.

**Implementation** (`src/OutOfClass.jl`, literate, ~100 LOC):
- δ = two-order convergence defect — from the SAME scaled jet the step builds,
  form `[m,n]` and `[m-1,n-1]` Padé (truncate by two coeffs, no extra ODE),
  evaluate both at the step endpoint t=1, take `|P − P_lo|/(|P|+tiny)`.
- **History-gated throw**: `OutOfClassError` iff δ > τ AND δ grew monotonically
  over the last K transitions. The gate is the safety mechanism — it fires on
  the e^{1/z} approach (sustained growth) but **cannot** trip on the single-step
  across-0 bridge (CFail.1e — no ≥K history) or a one-off pole-bridge δ spike.
- Calibration DERIVED from a separation probe (`external/probes/out-of-class-
  guard/`), never fitted: max in-class δ across the whole pole-bridge corpus =
  **3.9e-7** (hardest case: a single big step bracketing a double pole), e^{1/z}
  fires at δ ≈ 0.99 → **τ = 1e-3** (~3.8 orders of margin), **K = 2**.
- Default-on with `check_in_class=false` escape hatch (bitwise-identical
  numerics when off). `PadeStepper.pade_step_with_pade!` left byte-for-byte
  unchanged (a SEPARATE checked stepper) → zero `path_network_solve` hot-path
  regression and the `@inferred` type-stability gate untouched.

CFail.1c reframed silent-lie → `@test_throws OutOfClassError` (+ an escape-hatch
test pinning the legacy behaviour); CFail.1d flipped `@test_broken` → `@test`.
Mutation-proven with two independent guard mutations (τ→1e10; flip the monotone
comparison), each reddening exactly CFail.1c+1d → both conjuncts load-bearing.

**Result**: full `Pkg.test` GREEN at **9333 pass / 2 broken / 0 fail** (28m07s);
no false positive across the pole-bridge corpus. The broken count 3→2 was
lockstepped into README, CLAUDE.md, `scripts/quality_gate.sh`, and the
`corpus-v2-expected-broken-count` bd memory.

### Follow-ons filed
- out-of-class guard for `path_network_solve` / the vector wedge walk (deferred
  until a wedge out-of-class probe is motivated; the δ formula is type-generic).
- CFail.2 Chazy natural-boundary fixture (same δ-runaway discriminator should
  catch it; needs a far-side oracle).
- `OutOfClass._evaluate_pade_deriv_via` DRY cleanup (import the PadeStepper twin
  rather than re-derive — P4, not worth re-running a green 28-min suite for).

## Lessons

1. **Reconcile before building.** The "next bead" in a long-running project can
   be secretly-done bookkeeping. Cheap read-only scouts against worklogs + git +
   the actual artifacts paid for themselves: the real consequential work was a
   soundness bug two priority tiers down, not the in-progress P1 epic.
2. **The discriminator subtlety is the whole fix.** A guard built on the obvious
   signal (radius/coefficient growth) would have passed the marker test and
   silently broken pole bridging. The radius is blind to singularity *type*;
   only Padé *convergence* sees it. Getting the math right up front (and
   overriding the recon) was the load-bearing decision.
3. **Soundness > features.** A confirmed silent-wrong-answer bug in a
   correctness-claiming library outranks an 85 %-done feature epic, even at P3.
