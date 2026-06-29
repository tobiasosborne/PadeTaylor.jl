# Worklog 078 — FW 2011 Fig 4.7 "seam": cure scoping, mechanism selection, and confirmation

**Date:** 2026-06-29. **Author:** Claude Opus 4.8 (1M).
**Beads:** `padetaylor-vwgl` (P0, the bug — cure direction now fixed to
bounded-window compositing), `padetaylor-sny7` (P1, the gate — now a
pole-set seed-invariance test). Cure/feature beads: `padetaylor-xwzf`
(new — `windowed_path_network_solve`), `padetaylor-ingn` (new — verify
the suppressed poles are spurious, not real), `padetaylor-mro9` (new —
per-window edge-gating with a window-local seed). Field-only siblings
re-scoped: `padetaylor-fe9`, `padetaylor-8dg`.
**Probes:** `external/probes/fig47-seam-diagnosis/p3_edgegated_seedinvariance.jl`,
`p4_mechanism_selection.jl`, `p4b_overlap_graph.jl`,
`p5_full_composite_confirm.jl`.
**Predecessor:** worklog 077 (the root cause). **Status:** **SCOPING
ONLY. No cure code applied** — the project remains frozen behind `vwgl`
per the maintainer directive. This worklog records the auditioned
solutions, the measured mechanism selection, and the recommended plan.

> **Take-home.** Five cure families were auditioned and adversarially
> verified; none is a clean monotherapy. Three measurements then settled
> the design. (1) The seam is a pole-field-**interior** grain boundary —
> edge-gating, already shipped and already in the figure, does **not**
> cure it. (2) A bounded-window composite (FW's own md:147 remedy) shrinks
> the pole-field seed-variance by construction; reconciliation's
> cross-branch constraints are too sparse (~1.9 % of disc overlaps) to be
> the primary cure. (3) A full 5×5 composite over `[-50,50]²` drops the
> two-seed pole-count mismatch from **151 to 12** (pole-set agreement
> **77 % → 99.4 %**, median pole displacement **→ 0.000**), with **no
> tile-boundary seam**, and runs **~2× faster** than the monolithic solve.
> Recommended cure: `windowed_path_network_solve` + a pole-set
> seed-invariance gate.

---

## 1. The auditioned solutions (5 families, all adversarially verified)

Each family was developed in depth against the real code and then
adversarially critiqued. Every one returned a **`partial`** verdict — the
useful result is *why*, and which regime each addresses.

| Family | Mechanism | Verdict — the decisive weakness |
|---|---|---|
| **Bounded-window composite** (FW md:147) | independent ≤20×20 tiles from the IC, composited | partial — *containment*; far tiles ungated → smooth-sector bloom unless windowed-edge-gated |
| **Edge-gated + 2D harmonic fill** (`fe9`) | confine IVP to pole field, fill smooth region by a global elliptic solve | partial — `extract_poles` reads the **unchanged** solve → never touches the pole scatter; no-op on dense fields |
| **Edge-gated + per-row BVP fill** (`8dg`) | FW md:190 per-row Chebyshev–Newton with an asymptotic outer BC | partial — same category error: cures the |u| **field** (Fig 4.1), not the pole **scatter** (Fig 4.7) |
| **Cross-branch Gauss–Newton consensus** | re-derive node states from disc-overlap constraints | partial — variational mode `δu''=12u·δu` is the near-null-space of overlap consistency; damping needed to converge = seed-independent **smoothing** (forbidden non-fix) |
| **Wildcard: walk→topology-only + global least-squares** | discard path-accumulated state, re-derive from IC+overlap+asymptotic | partial — the asymptotic far-field BC is load-bearing, likely sector-wrong (single-branch −√(−z/6) across PI's 5 sectors), and inapplicable to the edge-gated figures |

**Cross-cutting findings that reshaped the plan**

- The seam has **two regimes**: (A) smooth-sector field path-dependence
  (the dominant *magnitude*; the IVP-unstable sectors) and (B) the
  pole-lattice grain boundary (the *visible* Fig 4.7 artifact, since the
  figure plots pole locations, not the field).
- "Confine the IVP to the pole field" (FW md:401) is **already shipped**
  (`edge_gated_pole_field_solve`) and **already in `fig_4_7.jl:106`** —
  so the stabilization families' step 1 is done, and the seam survives it.
- The smooth-fill families (`fe9`/`8dg`) produce **no poles**, so they
  cannot move the pole scatter — they are **field-only** (valuable for
  Fig 4.1, not for `vwgl`).
- Any seed-invariance gate built only from self-consistency invariants is
  **gameable** (freeze-shuffle + `enforce_real_axis_symmetry` passes all
  three internal checks; over-damped reconciliation smooths to agreement).
  The gate must embed a **physical-truth anchor**.

## 2. P3 — the decisive measurement: regime B (pole-field interior)

`p3_edgegated_seedinvariance.jl` re-solves the seed-0 edge-gated
**confinement** (the figure's exact target set) at two seeds. Confinement
is what the figure already does; the question was whether it removes the
seed-dependence.

| metric (panel e, seeds 0/42) | A: plain full grid | B: edge-gated, confined |
|---|---|---|
| pole count | 8159 / 7987 | 8137 / 7966 |
| **pole-count seed-variance \|Δ\|** | **172** | **171** |
| pole-set match @0.5 | 76.6 % | 81.5 % |
| max pole NN displacement | 1.429 | 1.429 |

**Verdict: regime B.** Confining the walk to the pole field leaves the
pole-scatter seed-variance essentially unchanged (172 → 171). The grain
boundary is *inside* the pole field. (Nuance: worklog 077's "35 % of
cells, max rel 1.0" is near-pole-inclusive — the poles' huge-|u| halos;
off-pole the smooth background is mostly seed-stable. The bug is the pole
lattice jittering, confirming regime B a second way.) This **eliminates
`fe9`/`8dg` as cures for the bug** and forces the cure to act inside the
pole field.

## 3. P4 / P4b — mechanism selection (windowing vs reconciliation)

**(a) Windowing works — even for far/long-trunk tiles** (`p4_mechanism_selection.jl`):

| core region | monolithic \|Δ\| / match / maxNN | windowed \|Δ\| / match / maxNN |
|---|---|---|
| NEAR [-10,10]² (short trunk) | 9 / 97.9 % / 1.52 | **2 / 100 % / 0.38** |
| FAR [18,38]×[-10,10] (long trunk) | 13 / 100 % / 0.47 | **1 / 99 % / 0.68** |

The far window collapses to |Δ|=1 — refuting the "long trunk is as bad as
monolithic" objection. All cells in one window share the trunk
(common-mode), so their *differential* positions are seed-stable.

**(b) Reconciliation's cross-branch constraints are sparse** (`p4b_overlap_graph.jl`):

```
17776 nodes · 48529 disc-overlap edges · avg degree 5.5 · 99.9% of nodes ≥2 overlaps
CROSS-BRANCH overlaps (nearest common ancestor >16 hops): global 1.7%, seam annulus 1.9%
```

~98 % of disc overlaps are same-branch (already mutually consistent — the
1.4e-4 stitch result). Only ~1.9 % carry cross-branch information. A
global least-squares is dominated by the same-branch bulk (the
unstable-mode near-null-space), so the seam is weakly forced → robust
convergence needs regularization = the forbidden smoothing. **Windowing
eliminates the cross-branch differential by construction; reconciliation
tries to reconcile it after the fact against a 1.9 %-sparse signal.**

## 4. P5 — full-composite confirmation (the acceptance test for #1)

`p5_full_composite_confirm.jl` builds a throwaway 5×5 overlapping
composite over `[-50,50]²` (each window solved wide = core ± margin 6,
from the IC, at `window_seed = global_seed*1000 + window_index`; poles
kept only in each window's Voronoi core).

| metric (poles in [-49,49]², seeds 0/42) | Monolithic | **Windowed composite** |
|---|---|---|
| pole count | 7443 / 7292 | 7223 / 7211 |
| **pole-count seed-variance \|Δ\|** | **151** | **12** |
| pole-set match @0.5 | 77.1 % | **99.4 %** |
| median NN displacement | 0.076 | **0.000** |
| wall time (both seeds) | ~214 s | **~102 s** |

**Tile-boundary residual:** boundary-band match 98.9 % vs interior 99.2 %
— comparable, **no boundary seam**; the overlap-core scheme handled it.

#1 is confirmed: global |Δ| **151 → 12** (~12×), match **77 → 99.4 %**,
median pole displacement **→ 0.000**, boundaries clean, **~2× faster**.

## 5. Recommended Phase-2 plan (gated behind the freeze)

**Cure — `windowed_path_network_solve`** (`padetaylor-xwzf`), bounded-window
composite per FW md:147:
- New `src/WindowedComposite.jl`: tile into ≤20×20 windows, each solved
  wide (core + overlap margin) from the IC, with
  `window_seed = combine(rng_seed, window_idx)` (genuine per-seed
  re-randomisation → non-gameable). Composite by Voronoi-core assignment;
  `windowed_extract_poles` unions core poles + reuses `PoleField`
  diameter-capped merge for boundary dedup. Rewire `fw2011_fig_4_7.jl` and
  `fw2011_fig_4_8.jl` (the composites FW md:147 names).
- Per-window edge-gating (`padetaylor-mro9`) is a refinement (trims
  far-window smooth-sector bloom), deferred — plain windows already hit
  |Δ|=12. Needs a window-local seed because `edge_gated_pole_field_solve`
  seeds from cells near the IC (`EdgeGatedSolve.jl:305-318`).

**Gate — `test/field_seam_test.jl`** (`padetaylor-sny7`), truth-anchored,
land RED first: assert composite pole-set seed-invariance (`|Δcount|` +
Hausdorff) over two global seeds; RED on monolithic (|Δ|=151 / 77 %),
GREEN on composite (|Δ|=12 / 99.4 %); tolerance pinned at the measured
composite floor (**not** relaxed — Rule 2); cross-check the Weierstrass-℘
recall/precision oracle so it cannot be gamed by smoothing.

**Honest residuals carried forward (Rule 9):**
- **Residual |Δ|=12** (~0.6 % unmatched; a few poles at maxNN ~0.9–1.5):
  this is *containment*, not annihilation — matching FW's own posture. The
  `sny7` tolerance sits at this measured floor.
- **Composite finds ~220 fewer poles** than monolithic (7223 vs 7443):
  almost certainly the spurious *seed-dependent* poles correctly
  suppressed (the count gap *is* the bug), but `padetaylor-ingn` must
  confirm via the ℘ truth-anchor that these are spurious, not real poles
  dropped — a named verification, not an assumption.

## Reproduction

```
julia --project=. external/probes/fig47-seam-diagnosis/p3_edgegated_seedinvariance.jl   # regime B
julia --project=. external/probes/fig47-seam-diagnosis/p4_mechanism_selection.jl        # windowing (a)
julia --project=. external/probes/fig47-seam-diagnosis/p4b_overlap_graph.jl             # reconciliation (b)
julia --project=. external/probes/fig47-seam-diagnosis/p5_full_composite_confirm.jl     # full-composite confirm
```
