# Worklog 086 — root cause of the default-order vector walk's far-wedge error (bead `padetaylor-rxbo`)

**Date**: 2026-09-29
**Scope**: Diagnosis only.  Find, by measurement, why the default-order
(`rng = nothing`, `h = 0.1`) vector path-network walk on the P_I^(2) fixture
is wrong by up to 2.18 relative for |z| ≥ 5.82 (worklog 082 §2-3,
`test/vector_field_seam_test.jl` VFSEAM.3).  **No `src/` change, no test
change**: VFSEAM.3 stays `@test_broken`.  Probes live in
`external/probes/vector-default-walk-rxbo/probe{1..7}.jl` (probe1 serialises
the three walks to `walks.jls`, ~6.8 MB, deleted before hand-off and not
committed; re-run probe1 (~1 min) before probes 2-7, which read it).

> Take-home: the error is NOT a per-step truncation defect.  Every edge of
> the worst chain has local error ~1e-12 (restarted from the reference state).
> The default target ORDER (`figures/_kkg_pi2_helpers.jl:163-164`: r = 2, 4, 6,
> 8 and within each ring θ = −0.5 first) makes the walk reach the far wedge
> first at 6·e^{−0.5i}.  To get there it squeezes past a pole at
> z ≈ 5.0-5.15 − 2.46i: ‖y‖ up to 7.06e5, adaptive h down to 0.0111.
> `_nearest_visited` (`src/VectorPathNetwork.jl:772`) then grows EVERY later
> far-wedge target from that subtree.  So all far nodes descend from the
> passage, and the upper far wedge is reached by a ~125-step creep along the
> |z| ≈ 5.8 arc (arg −0.33 → +0.43) in which the IVP amplifies relative error
> by 1e5-6e6.  Error on the chain: 3.7e-10 before the passage, 5.1e-7 after
> it (×1400), 1.9 at the end of the arc.  The seeded walks reach the upper far
> wedge directly: along their far chains ‖y‖ ≤ 5.41e4 and h ≥ 0.0188, and the
> end-to-end amplification is 1.5e5 instead of ≥ 1.4e13 (saturated).

## Ground truth opened

- `docs/worklog/082-vector-two-seed-field-gate.md` (the measurement), `test/vector_field_seam_test.jl:125-139` (VFSEAM.3).
- `docs/worklog/077-fig47-path-dependence-seam-rootcause.md` (the scalar analogue).
- `src/VectorPathNetwork.jl:723` (shuffle only when `rng !== nothing`), `:772` (each target starts from the nearest visited node), `:813-821` (adaptive h, `_select_wedge`).
- `src/VectorWedgeStep.jl:368-386` (`_candidate_pole_disc`), `:421-481` (`_select_wedge`), `:569-596` (`_adaptive_h`: `h = SAFETY·D_local`, `SAFETY = 0.10` at `:281`).
- `figures/_kkg_pi2_helpers.jl:157,163-164,293-296,339` (h = 0.1, target radii/angles, fan order, the figure's walk).

## Reference used

At each default-walk node the reference is the seed-0 walk (h = 0.1) evaluated by Stage 2
(`fine_grid = wd.visited_z`).  Where seed 0 does not cover the node, the default-order
h = 0.05 walk is used instead.  Both sided with each other to ≤ 7.7e-5 in worklog 082 §3.
Metric as in VFSEAM: `max_c |Δy_c| / max(1, |y_c|)`.  Default walk: 795 nodes, 473 covered,
126 above 1e-4 (the count differs from VFSEAM.3's 65 because VFSEAM.3 evaluates at seed-0
nodes and this counts default nodes).

## 1. No clean→bad edge: the error is accumulated, not injected by one step (probe1, probe3)

- probe1: of the bad nodes, **0** have a parent that agrees to ≤ 1e-6.  The earliest bad
  ancestor is node 415 (|z| = 5.810, rel 1.58e-4); its parent 414 is at 5.53e-5.
- probe3, local error: every edge on the chain of the worst node (493, z = 5.5407+2.5441i,
  rel 1.94, chain depth 326) was re-stepped from the REFERENCE parent state with the same
  Δz and order 24.  Local error is 1e-16…7e-12 on every edge.  The two exceptions (9.96e-6
  at 466→467 and 4.30e-5 at 478→479) fall where the reference's own Stage-2 value is noisy.
  Accumulated error along the chain:

| chain pos | node | |z| | accumulated rel |
|---|---|---|---|
| 150 | 290 | 5.248 | 3.68e-10 |
| 200 | 367 | 5.726 | 5.14e-7 |
| 240 | 407 | 5.788 | 2.64e-5 |
| 258 | 425 | 5.957 | 1.07e-3 |
| 298 | 465 | 6.032 | 1.69e-1 |
| 326 | 493 | 6.097 | 1.94 |

## 2. Where the chain goes, and the near-pole passage (probe2, probe5, probe7)

The chain of node 493 starts at z = −3 and runs through the LOWER half-plane
(0.69 − 0.99i, 2.76 − 1.82i, 4.41 − 1.81i).  It squeezes past a pole at
z ≈ 4.97-5.15 − 2.44…−2.49i, then creeps along |z| ≈ 5.75-6.05 from
5.41 − 1.93i up through the real axis to 5.54 + 2.54i.  Node values on it:

| node | z | ‖y‖ | h (adaptive) | nearest shared-Q pole |
|---|---|---|---|---|
| 301 | 4.914 − 2.413i | 3.26e4 | 0.0217 | 0.205 |
| 311 | 5.052 − 2.461i | 3.79e5 | 0.0137 | 0.132 |
| 316 | 5.105 − 2.490i | 6.99e5 | 0.0117 | — |
| 348 | 5.147 − 2.468i | 1.79e5 | 0.0133 | 0.142 |
| 368 | 5.410 − 1.934i | 2.75e2 | 0.0502 | 0.471 |

probe7 took every far default node (|z| ≥ 5.5) that the reference covers, bad AND good
(126 each), and recorded the maximum ‖y‖ and minimum h over its ancestor chain.  **All 252
chains pass through the same passage: chain max‖y‖ = 7.06e5 and min h = 0.0111 for every
one of them.**  Nodes 305, 308, 312, 316 and 349 each lie on all 126 bad chains.  For
comparison, the seed-0 far nodes (482): chain max‖y‖ has median 1.03e4 and max 5.41e4;
chain min h has min 0.0188.

Why the default order does this.  `kkg_target_wedge()` lists targets in the order
r ∈ (2, 4, 6, 8), and within each ring θ ∈ (−0.5, −0.25, 0, 0.25, 0.5)
(`figures/_kkg_pi2_helpers.jl:163-164,293-296`).  So the first far target is 6·e^{−0.5i}
= 5.27 − 2.88i, beyond the lower passage.  Every later far target is started from
`_nearest_visited` (`src/VectorPathNetwork.jl:772`), which is on that subtree.  The whole
far wedge is therefore one subtree grown tangentially from the lower-wedge branch.  A
shuffled order usually walks some upper far target first, straight from the inner wedge.

## 3. Amplification: default chain vs seed-0 chain (probe4, probe5, probe6)

Method: add a relative kick ε to the state at chain position i, re-step the SAME z-sequence
to the end, and take the ratio (end change)/ε.

- probe4, restarting from the reference state: from pos 150 (|z| = 5.248, before the
  passage) the end error is 1.04 and the amplification (ε = 1e-10) is 7.7e10.  From pos 200
  (|z| = 5.726, after the passage) the end error is 8.6e-5 and the amplification 7.5e5.
  That 8.6e-5 is about the reference's own ~1e-10 error at node 367 times 7.5e5.  Restarting
  from pos 300 gives 3.8e-10.  Separately, replaying the default chain from its own node-1
  state (Δz recomputed as a difference, i.e. rounding-level changes) lands 2.18 away from
  the stored end value: the chain is chaotic at roundoff level.
- probe5 (ε = 1e-13).  Default chain to node 493: amplification 1.4e13 from pos 1.  It is
  saturated at O(1)/1e-13 for every pos ≤ 171; then 3.4e11 at pos 181, 4.5e8 at pos 191,
  5.8e6 at pos 201, 2.7e5 at pos 241 and 1.2e4 at pos 271.  **Seed-0 chain to its node
  nearest z_493 (node 263, 0.04 away; depth 178, all in the UPPER half-plane): 1.5e5 from
  pos 1, and ≤ 660 from every position with |z| ≥ 1.1.**
- probe6 (ε = 1e-12, each edge split into s sub-steps):

| segment (pos) | |z| | s = 1 | s = 2 | s = 4 |
|---|---|---|---|---|
| 161→181 (passage) | 5.47→5.71 | 431 | 29 | 2.3 |
| 171→191 (passage exit) | 5.62→5.73 | 1.48e4 | 2.42e3 | 2.95e3 |
| 181→201 | 5.71→5.75 | 1.80e5 | 1.72e5 | 1.80e5 |
| 241→326 (arc creep) | 5.77→6.10 | 2.10e5 | 2.11e5 | 1.79e6 |

Inside the passage the amplification depends on the step: ×431 at the walk's own steps
against ×2.3 with 4 sub-steps.  That is a property of the discrete Padé map there, not of
the ODE.  After the passage, from ‖y‖ ≈ 1.8e5 down to 275 (181→201), and along the arc
(241→326), the amplification does not change with the step: that part is the ODE, measured
in the max(1,|y|) metric.

## 4. Verdict

This is an **under-guarded configuration of the walk, with a numerically sensitive
near-pole passage at its root**.  It is not a single wrong step, and not a
local-truncation bug.  Two factors multiply:

1. **A passage too close to a pole (walk-controllable).**  The default order routes the far
   wedge through z ≈ 5.0-5.15 − 2.46i at ‖y‖ up to 7.06e5 and h down to 0.0111.  Across it
   the chain error grows 3.7e-10 → 5.1e-7 (×1400), partly through step-dependent
   amplification of the discrete map (×431 vs ×2.3 at 4 sub-steps).
2. **A long tangential creep through a sensitive region (order-controllable).**  Because
   every far target starts from that subtree (`:772`), the upper far wedge is reached along
   the |z| ≈ 5.8 arc, which amplifies relative error by 2e5-6e6.  The seed walks reach the
   same points with 1.5e5 end-to-end.

**Invariant the walk fails to enforce.**  The walk has no bound on how close a Stage-1 chain
may pass to a pole: nothing caps ‖y‖ along a chain or floors h relative to h_max.
`_adaptive_h` (`src/VectorWedgeStep.jl:596`) only shrinks the step, down to 0.111·h_max
here; it never re-routes or flags.  `:max_q_root` (`:465-470`) scores only the landed
node's pole-free disc, one step at a time.  The walk also never checks that the error
propagated along a chain stays bounded.

**Candidate guard (measured separation, not yet implemented).**  A per-chain near-pole
signature separates the walks cleanly:

| walk | chain max‖y‖ over far nodes | chain min h over far nodes |
|---|---|---|
| default | 7.06e5 on every far chain | 0.0111 on every far chain |
| seed 0 | ≤ 5.41e4 | ≥ 0.0188 |

This makes it a walk-level flag, not a per-node one: in the default walk the bad and good
far nodes share the same passage.  A thresholded guard, e.g. "warn/throw when a chain passes
at ‖y‖ > 1e5 or h < 0.15·h_max", is uncalibrated beyond this one fixture.  The robust
remedy is the FW/FFW two-run check that VFSEAM already does: cross-evaluate two orderings'
fields and flag disagreement.

**Not measured (time box).**
- Which wedge direction `_select_wedge` picked at nodes 301-349, i.e. whether an off-goal
  wedge would have kept the chain farther from the pole.
- Whether the step dependence in the passage comes from shared-Q cell-A/B dispatch
  (`shared_pade_select`) or from Padé conditioning.
- Whether any `src/` change (e.g. an h floor with re-routing) cures the far wedge without
  breaking the other suites.

## 5. What was NOT changed, and why

- **VFSEAM.3's `@test_broken` is kept.**  No `src/` fix was made in the time box, and a fix
  in the walk's routing needs its own ADR: it changes the default-order behaviour pinned by
  ADR-0001 bit-identity.
- **The shipped figure is NOT switched to another ordering/step here.**  Worklog 082 §3 and
  this worklog support two candidate fixes at the figure level: (a) h = 0.05 in the default
  order, which agreed with the seeds to ≤ 7.65e-5 on the far wedge; (b) an upper-first or
  seeded order.  Either one changes the pinned counts in `test/kkg_pi2_figure_test.jl` and
  the published Panel B.  That is a separate, reviewed change.
- The stale "389 nodes, 21 poles" text (`figures/_kkg_pi2_helpers.jl:148`,
  `test/kkg_pi2_figure_test.jl:47`) is untouched (out of scope for a diagnosis).
