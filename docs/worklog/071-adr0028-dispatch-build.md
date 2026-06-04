# Worklog 071 — ADR-0028 dual-construction dispatch: audition → build → Froissart fix

**Session 2026-06-04.** Took ADR-0028 (shared-Padé dual-construction dispatch)
from *proposed* to *built & shipped*, via an empirical audition that corrected the
design, a calibrated build, and one genuine build-time discovery. Beads
`padetaylor-flnr` (build), `padetaylor-0ql3` (audition), `padetaylor-7nm5`
(discriminator research), the Froissart-consumer fix, `padetaylor-unk` (recovered).

## The arc

1. **Audition first (probe `external/probes/adr0028-dual-construction-audition/`).**
   Before building the gated design, prototyped cell A (shipped GGT diagonal `(m,m)`,
   asserted bit-identical to `shared_denominator_pade`) and three cell-B variants,
   and scored selection rules against ground truth. **Three corrections to the ADR**
   (Amendment 1):
   - cell B must be the **wide-square** Mano–Tsuda at degree `m_eff = ⌈m/d⌉·d`, NOT
     the ADR's `⌈m/d⌉`-keep-degree (over-determined when `d∤m`, ~1e9× worse on
     exp/CM/NY — the "source-silent" `d∤m` gap the Mano–Tsuda scout flagged, decisive);
   - the dispatch is justified but the win-boundary is `{genuine pole × precision}`,
     not entire-vs-meromorphic (cell A wins BF-256 Calogero–Moser ×900 — the FW
     Table 5.1 regime);
   - the `(R,g,K)` Pareto + conditioning gap `g` are unreliable (6/8) — `g` doesn't
     separate the regimes and isn't comparable across the tall-A/wide-B shapes.

2. **Discriminator research + the uniform selector (Amendment 2).** Surveyed
   defect/residual control (Enright, Higham, Corless–Kaya 2025), singularity-from-jet
   (Domb–Sykes, Mercer–Roberts), AAA. The **relative ODE defect**
   `‖ỹ'(t) − h·f(z₀+h·t, y(t))‖/‖h·f‖` is a *single uniform* selector — **18/18**
   across pole-free / pole-crossing / BigFloat / BF-crossing (probe `pole_crossing.jl`,
   with a genuine vector pole-crossing oracle: tan-companion `u''=2u+2u³` → `(tan z,
   sec²z)`). It needs no held-out point, no pole detection, no fallback. Cell A never
   wins a pole-*crossing* step; its only win is a *regular* high-precision near-pole step.

3. **The build (Amendment 3).** Three modules — `SharedPadeCellB.build_square_cell`,
   `SharedPadeDefect.relative_defect`, `SharedPadeDispatch.shared_pade_select`
   (d=1 → cell A; build A always; try B guarded → fall back to A on any degeneracy;
   A-default type-scaled tie-break) — wired at the single site `VectorStepper.jl:242`.
   `shared_denominator_pade` + the d=1 bit-identity (SP.1.1) untouched. 26 isolated
   unit tests; 4 mutants proven RED (cell-B window, always-A, always-B = the BF-CM
   A-win, the guard).

## The friction — the build-time discovery (Froissart consumers)

Wiring the stepper (Phase 2a) made every fixed-step test pass+improve but **collapsed
the ℘ FW Table 5.1 walk** (VPO.2/VPO.4 error). Root cause, and the lesson:

> **The audition tested isolated single steps; it could not see the value-vs-pole-
> fidelity conflict the production walk exposes.** Cell B optimises the step *value*
> (the defect), but on a regular stretch its wide-square `Q` can plant a **Froissart
> doublet** — a near-cancelling pole/zero pair, a spurious pole essentially *at* the
> node. Harmless for the value (the doublet cancels — the defect happily picks B),
> but the two consumers that read `Q`'s **roots** — the adaptive step controller
> (`_adaptive_h`) and the pole extractor (`extract_poles_shared_q`) — are corrupted.
> The **vector** root consumers lacked the residue Froissart filter the **scalar**
> `PoleField` already has.

A selector-level pole-fidelity veto was tried and **rejected** — it over-fires on
entire jets where the near-node root is harmless, throwing away the value gain; the
selector cannot know whether `Q`'s roots will be consumed. The fix belongs at the
consumer.

**Calibrated, not guessed** (the maintainer chose "investigate first"; probe
`external/probes/adr0028-froissart-consumer/`): across 342 ℘-walk nodes, genuine
℘-pole residues `≥ 2.1` (median 1.1e4) vs cell-B doublet residues `≤ 5.5e-8` (median
3e-13) — an **~8-order clean gap**. Threshold `POLE_MIN_RESIDUE = 1e-4` (the scalar
default `1e-8` is too low — a doublet reaches 5.5e-8). Lifted the scalar residue
criterion to the shared-Q form (`shared_q_residue = max_c |P_c(t)/Q'(t)|`) in
`_candidate_pole_disc`/`_adaptive_h` (nearest *genuine* pole) and
`extract_poles_shared_q`. **The filter touches pole-roots only — never the step
value** — so it is cell-agnostic and a strict correctness improvement. Mutation-proven:
filter off + dispatch wired ⇒ the ℘ walk collapses again.

## Re-baseline

The dispatch *evens* the shared-Q error across components. VS.1.2 (harmonic, order
30, h 0.7), measured: cell B → F64 cos 2.2e-16 / −sin 1.1e-16, BF-256 cos 1.76e-30 /
−sin 5.8e-31 — vs cell A's −sin ~5e-9 (F64) / ~1.6e-17 (BF), a ~7-/14-order win on
the worst component. Single-step entire/high-`d` tolerances tightened to measured ×
margin; accumulated-walk and figure bounds left as upper bounds.

## Lessons

- **Audition single steps, but a dispatch that changes the stored `Q` must be
  validated against everything that *consumes* `Q` — value AND roots.** The poles
  are this package's deliverable; a value-optimal cell can be pole-pessimal.
- **The selector cannot carry a concern it lacks the context for.** Pole fidelity is
  a consumer concern; fix it at the consumer, not by vetoing the value choice.
- **Calibrate thresholds on the real workload.** The scalar `min_residue = 1e-8`
  would have leaked the 5.5e-8 doublet; the ℘-walk data set it at 1e-4.

## Status

Full `Pkg.test()`: **7796 / 7796 GREEN** (19m42s, idle box — note a concurrent
heavy Julia session can OOM the suite, bd memory `full-pkg-test-got-sigterm`; the
first two attempts here were disrupted by a concurrent `MetricsQNM` session).
ADR-0028 → *accepted* (Amendment 3). `padetaylor-unk` recovered (entire-system
accuracy: −sin ~5e-9 → ~1e-16). Probes retained as the evidentiary record.
