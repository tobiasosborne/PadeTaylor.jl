# ADR-0033 — `solve_pade`: fail loud on out-of-class (non-meromorphic) input

**Status**: **accepted** — 2026-06-19. Resolves bug `padetaylor-v1ub` (the only
confirmed silent-wrong-answer bug; epic `padetaylor-25og`).
**Date**: 2026-06-19
**Beads**: `padetaylor-v1ub`. Closes the auto-flip marker CFail.1d in
`test/corpus_out_of_class_test.jl`.

## Context

`solve_pade` is **meromorphic-only**: its entire value proposition is that a
diagonal Padé can see *past a pole* where plain Taylor truncation diverges. But
it had **no guard** against an input whose solution has a NON-pole singularity,
and on such an input it returned finite, plausible, *wrong* values with **no
throw and no NaN** — a direct Rule-1 (fail-loud) violation in a library whose
whole point is honest results.

The motivating out-of-class input (Family G,
`docs/test_corpus/02_corpus_extension_plan.md:292-306`) is the essential
singularity

    u'' = u·(1 + 2z) / z⁴        exact   u(z) = e^{1/z},

whose `z = 0` is an essential singularity (Laurent series `1 + 1/z + 1/(2z²) + …`
has infinitely many negative-power terms — NOT a pole). Marching the **real**
`solve_pade` from `z₀ = -1` toward `0` at fixed `h = 0.1` (oracle `e^{1/z}`,
mpmath dps=50):

    |z|=0.5  relerr 2.1e-16     |z|=0.1  relerr 4.0e-10
    |z|=0.05 relerr 4.9e-2      |z|=0.02 relerr 1.2e+17  AND WRONG SIGN

The value at `z = -0.02` is finite, has the wrong sign (`u ≈ -2.4e-5`, yet
`e^{1/z} > 0` for all real z), and relerr ≈ 1.2e17 — a confident, finite,
plausible-looking lie. Confirmed empirically in `test/corpus_out_of_class_test.jl`
(bead `padetaylor-v1ub`; the file shipped CFail.1d as a `@test_broken` auto-flip
marker awaiting exactly this fix).

## Why the discriminator is *Padé convergence*, not *radius* (the locked design)

The obvious-looking guard — watch the Taylor jet's coefficient growth, or the
Jorba–Zou radius-of-convergence `h_JZ`, and trip when it collapses — **DOES NOT
WORK**, and shipping it would have broken the package's headline feature (pole
bridging). The radius of convergence equals the distance to the nearest
singularity for **both** a pole and an essential point, so it collapses
identically in both cases. A radius/coefficient-growth guard cannot distinguish
"approaching a pole I can bridge" from "approaching an essential singularity I
cannot."

The correct discriminator is **Padé sequence convergence**, grounded in GGT 2013
§8 (`references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
GGT2013_robust_pade_via_SVD_SIREV55.md`, the "Modified Padé Approximant with
Pointwise Convergence" section, ~md:386):

  - For a **meromorphic** function the Padé approximants converge (in capacity /
    in measure — the **Nuttall–Pommerenke** theorem [1,21,26], and for the
    `m → ∞`, `n` fixed case **de Montessus de Ballore** [1,20]), so
    successive-order approximants AGREE. This holds *including while bridging a
    pole* — a diagonal Padé arranges its denominator to vanish at the pole, and
    the `[m,n]` and `[m-1,n-1]` approximants of the same jet land on the same
    value at the step endpoint. This convergence is precisely *why* the package
    can see past a pole.

  - For an **essential** singularity (GGT 2013 **Fig. 9** is exactly
    `exp((z+1.5)⁻²)`, an essential point at `z = -1.5`) NO finite-degree rational
    converges — the problem is outside the convergence theory entirely.
    Successive-order approximants DISAGREE, and the disagreement grows as the jet
    senses the essential singularity more strongly.

## Decision

**1. Discriminator = two-order Padé convergence defect.** Per step, from the SAME
scaled Taylor jet the step already builds (`c̃_k = h^k c_k`), form the working-
order Padé `P = [m,n]` (m = n = order ÷ 2) AND the reduced-order
`P_lo = [m-1,n-1]` (the same jet truncated by two coefficients — **no extra
ODE/jet evaluation**), evaluate both at the step endpoint `t = 1`, and take

    δ = |P(1) − P_lo(1)| / (|P(1)| + tiny).

**2. Throw policy = history-gated sustained monotone growth.** The driver keeps a
short ring buffer of δ across consecutive steps and throws `OutOfClassError`
**iff** `δ > τ` AND δ has grown monotonically across the last `K` step-to-step
transitions (K ≥ 2). The history gate is the safety mechanism (see below).

**3. Default-on, with an escape hatch.** `solve_pade(...; check_in_class = true)`
by default (Rule 1). `check_in_class = false` restores the legacy unguarded
behaviour for a user knowingly probing out-of-class input. Disabling is
byte-for-byte identical numerics (verified: guard-on and guard-off `u(1.05)` on
the ℘ pole bridge are bitwise equal).

**4. Typed error.** `OutOfClassError <: Exception` (idiom precedent
`VectorWalkError`, `src/VectorWedgeStep.jl:211`) carrying `z`, `δ`, the δ history,
and a rich message with a `Suggestion:` line. Catchable and self-documenting.

## Why the history gate (and not a single-step threshold)

A single-step `δ > τ` test would be unsafe. The history gate is what lets the
guard fire on the e^{1/z} **approach** (sustained monotone growth over the last
steps) while it **cannot** fire on:

  - the across-0 single large step (CFail.1e): one step leaves no ≥K-step history
    → the monotone-growth gate is never satisfied → the headline "bridge across
    z = 0" stays green;
  - an isolated single-segment pole bridge (a one-off δ spike): no sustained
    history, and the spike sits far below τ anyway.

## Measured τ/K calibration (DERIVED from data, never fitted to pass)

A blocking separation probe (`external/probes/out-of-class-guard/probe.jl`,
`probe2.jl`) measured δ at every step for the e^{1/z} approach and for ≥5
legitimate pole-bridge cases drawn from the regression corpus. The two
populations separate by ~6 orders of magnitude:

| case | type | max δ (all steps) | firing δ |
|---|---|---|---|
| **e^{1/z} approach** (z₀=-1→-0.02, h=0.1) | OUT-OF-CLASS | — | δ climbs `2.7e-14 → 3.2e-10 → 0.99` |
| ℘ double pole z=1, march 0→0.99 (h≤0.05) | in-class | 3.6e-16 | — |
| ℘ double pole z=1, single-seg bridge (h=1.5) | in-class | **3.9e-7** | — |
| tan simple pole π/2, single-seg bridge (h=3) | in-class | 1.5e-11 | — |
| tan march 0→1.56 (h=0.05) | in-class | 2.3e-16 | — |
| coth real pole z=0, bridge -1→1 (h=2) | in-class | 1.2e-9 | — |
| 1/(1-z) pole z=1 (exact rational), bridge (h=1.5) | in-class | 3.0e-15 | — |
| tanh complex pole iπ/2, march 0→1.5i (h=0.1i) | in-class | 2.5e-16 | — |
| sn complex pole, march 0→1.0i (h=0.1i) | in-class | 1.8e-16 | — |

**Max in-class δ over every step = 3.9e-7** (the hardest case — a single big step
bracketing a DOUBLE pole). The e^{1/z} guard must fire at `δ ≈ 0.99`.

  - **`τ = 1e-3`** sits in the gap with **~3.8 orders** of margin above the
    in-class ceiling (3.9e-7) and **~3 orders** below the firing value (0.99).
  - **`K = 2`**. The e^{1/z} approach gives ≥2 consecutive monotone-growth
    transitions ending at the firing step (`2.7e-14 < 3.2e-10 < 0.99`); the
    across-0 single step and the isolated single-segment bridges give none.

The populations separate cleanly (no overlap, the pole bridge does NOT produce
sustained-growing δ), so the guard is safe — it does not break pole bridging.

## Consequences

- **The bug is closed.** `solve_pade` to `z = -0.02` on the e^{1/z} recast now
  throws `OutOfClassError`; the CFail.1d auto-flip marker (`refuses_out_of_class`)
  returns `true` and its `@test_broken` is flipped to a passing `@test`.
- **Far-from-singularity accuracy is unchanged.** solve_pade to z=-0.5/-0.2/-0.1/
  -0.05 still returns the finite accurate/degraded values (δ at/below the floor;
  the guard is silent). CFail.1a/1b keep their finite pins.
- **No false positives** on the pole-bridging / near-pole regression corpus:
  GREEN across `corpus_scalar_pole_bridge` (27), `problems` (25), `padestepper`
  (16), `corpus_higher_order_pole` (43), `corpus_elliptic_lattice` (30),
  `corpus_periodic_pole` (12), `corpus_riccati_rational` (19),
  `corpus_riccati_special` (19), `pathnetwork` (120), `adaptive_step` (43),
  `fw_fig_41` (5), `polefield` (43).
- **Zero hot-path regression.** `PadeStepper.pade_step_with_pade!` is left
  **byte-for-byte unchanged** (`git diff src/PadeStepper.jl` empty), so the
  `path_network_solve` tight loop and the `@inferred`/allocation type-stability
  gate are untouched. The δ logic lives in a SEPARATE checked stepper in the new
  literate `src/OutOfClass.jl`; the default `solve_pade(check_in_class=false)`
  path remains on the unchecked stepper at the original cost.

## Mutation-proof (Rule 4)

CFail.1c (`@test_throws OutOfClassError`) and CFail.1d (the flipped `@test`) are
the load-bearing assertions. TWO independent guard mutations, each disabling
firing, were run RED and restored byte-for-byte:

  - **M-guard-τ**: bump `OUT_OF_CLASS_TAU` 1e-3 → 1e10 (no δ ever exceeds τ).
    Result: exactly CFail.1c + CFail.1d reddened (12 pass / 2 fail); 1a/1b/1e
    stayed GREEN. Restored; 14/14 GREEN.
  - **M-guard-monotone**: flip the monotone test `history[i] > history[i-1]` → `<`
    (a growing δ never satisfies the gate). Result: identical reddening — only
    CFail.1c + CFail.1d (12 pass / 2 fail). Restored; 14/14 GREEN.

Each mutation kills one conjunct of `δ > τ AND monotone-growth`, confirming both
halves are necessary and neither is dead. (`git status --porcelain
src/OutOfClass.jl` empty after each restore.)

## Open / deferred

- **CFail.2 Chazy (natural boundary)** remains DEFERRED (no far-side oracle past a
  movable natural boundary; plan Family G). The same δ-runaway discriminator
  should catch it — a natural boundary is also non-meromorphic — but it is not
  pinned here for lack of an oracle. Follow-on bead if/when an oracle is built.
- **Complex-direction `path_network_solve`** is NOT guarded by this change — the
  guard lives in the real-axis `solve_pade` driver. The δ formula
  (`two_order_defect`) is element-type-generic and would port to the wedge walk
  if an out-of-class probe there is ever needed; deferred until motivated.

## References

- GGT 2013 §8 (Nuttall–Pommerenke / de Montessus de Ballore; Fig. 9 essential-
  singularity example) — `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
  GGT2013_robust_pade_via_SVD_SIREV55.md` (~md:386).
- `src/OutOfClass.jl` — the literate guard module (δ formula, checked stepper,
  history-gated checker, `OutOfClassError`).
- `src/Problems.jl:230-252` — the `solve_pade` wire-in (`check_in_class` kwarg).
- `test/corpus_out_of_class_test.jl` — CFail corpus (CFail.1c throw, CFail.1d
  flip, the mutation-proof footer).
- `external/probes/out-of-class-guard/` — the re-runnable separation probe.
- `docs/test_corpus/02_corpus_extension_plan.md:292-306` — Family G recipe.
- `src/VectorWedgeStep.jl:211` — `VectorWalkError`, the typed-error idiom.
- bug `padetaylor-v1ub`.
