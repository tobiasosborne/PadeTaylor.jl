# ADR-0032 — `extract_poles`: disjoint-support self-merge for dense pole fields

**Status**: **accepted** — 2026-06-14. Resolves bug `padetaylor-fzse` (epic
`padetaylor-l7yt`, the bug-fix campaign).
**Date**: 2026-06-14
**Beads**: `padetaylor-fzse`; co-design note for `padetaylor-90oh` (multiplicity API).

## Context

`padetaylor-fzse` reported that `extract_poles` on a dense 2D path-network pole
field (the equianharmonic ℘ lattice of `figures/demo_lattice_singularities.jl`)
has "~76% precision / 76% recall": over-split duplicates, **4 interstitial false
positives**, and 4 missed poles. The bead and the read-only audition both
designed a fix around *killing the false positives* (an agreement-spread gate, a
residue bump, tighter `radius_t`).

## What the calibration probe established (the bead premise is wrong)

A blocking calibration probe (`external/probes/fzse-calibration/probe.jl`) on the
actual 121×121 field, classified against the **closed-form lattice including
out-of-window poles**, overturned the diagnosis (the same pattern as
`padetaylor-61um`):

- **At the DEFAULT params (`cluster_atol = 0.1, min_support = 3, radius_t = 5`)
  precision is 1.00**: all 32 reps lie within **0.208** of a real lattice point.
  There are **no interstitial false positives and no aliases.**
- The bead's "4 FPs at 0.93–1.24 from any true pole" were (a) genuine
  **out-of-window** lattice poles that boundary nodes legitimately see (mislabeled
  by measuring distance to in-window poles only), and (b) an **artifact of the
  bead's own tuning** to `cluster_atol = 0.4`, which lets node-local roots gather
  `min_support` and survive as aliases (the probe measured 4 alias reps at
  `0.4`, and 21 at `0.4, min_support = 2`). **Widening `cluster_atol` is the
  wrong cure** — it trades over-splitting for aliases.
- The genuine defect is **over-splitting only**: the greedy first-fit cluster
  cannot merge a *chain* of representatives wider than `cluster_atol`, so a pole
  whose cross-node estimates spread beyond `cluster_atol` fragments into 2–4
  duplicate reps. 4 in-window poles are duplicated at the default.
- The recall gap (5/17 missed) is on the **bottom-edge rows** and is mostly a
  **walk-coverage** limit (3 of 5 are seen by `< 2` nodes — no clustering change
  recovers them; lowering `min_support` only adds aliases). It is **orthogonal**
  to pole extraction.

## Decision

Add a **post-pass single-linkage self-merge** to `PoleField._extract_poles_core`
(`merge_atol` kwarg, default `h_max = maximum(|h|)`), run *after* the
`min_support` filter. Two reps are merged iff

    |rep_a − rep_b| ≤ merge_atol   AND   support(a) ∩ support(b) = ∅.

The merge keeps the **best-placed** rep (smallest z_dist ⇒ smallest index). No
agreement-spread gate, no residue bump, no `radius_t` change, no `min_support`
change — the calibration showed those are unnecessary (no FPs) or harmful (they
add aliases / mask coverage).

## Why the disjoint-support condition (distance alone is insufficient)

Distance cannot distinguish **over-split froth** (one pole seen by many nodes)
from a genuine **near-coalescent pole pair** (two distinct poles at small
separation). The **support sets** can: a node sees one pole at one `t*`, so it
contributes to exactly one froth rep — the froth reps of a single pole have
**disjoint** support. A coalescent pair gives each node two roots, so both reps
**share** that node. Merging only disjoint-support reps therefore collapses the
cross-node froth while **never** merging a coalescent pair.

This was found empirically: a distance-only merge at `merge_atol = h_max`
regressed `CPN.4` (a δ-sweep near-coalescent pair) and `CRic.3` (the Hermite
H₄′/H₄ zero-pair) and changed the `ffw_fig_4`/`ffw_fig_6` pole-counts; the
disjoint-support condition fixed all four while keeping the over-split fix.

A welcome consequence: because the disjoint-support condition (not the radius)
protects distinct poles, `merge_atol` need not be tuned tightly — it defaults to
`h_max`, a natural scale, with no upper-bound fragility.

## Calibration & results

`legacy (merge_atol = 0)` → 32 reps, 4 over-split dup-poles. Disjoint-support
self-merge → **24 reps, 2 dup-poles, precision 1.0**. The 2 residual dups are the
℘ **double** pole resolved by a single node into a doublet (two roots
~0.15–0.2 apart, shared support) — see below. Recall is unchanged (24 reps cover
the same lattice poles the walk visited).

## Consequences & scope

- **Default-on.** `extract_poles` now self-merges by default; `merge_atol = 0`
  recovers the exact legacy clustering. Well-separated fields (all existing
  PF.* / figure tests) are unaffected — distinct poles share support and/or sit
  beyond `merge_atol`. Verified GREEN across 12 in-suite scalar consumers.
- The discovery demo (`figures/demo_lattice_singularities.jl`) drops its
  `cluster_atol = 0.4` workaround for the ℘ panel and uses stock defaults.
- **Mutation-proven** (`polefield_test.jl` M7/M8): merge-off reddens PF.5.1;
  dropping `isdisjoint` reddens the coalescent-pair fixtures.

## Open / deferred

- **Double-pole multiplicity → `padetaylor-90oh`.** A 2nd-order pole that a
  *single* node splits into two roots (shared support) is deliberately **not**
  merged here — by distance it is indistinguishable from a coalescent pair. The
  multiplicity API should report it as one location of order 2 rather than two
  locations; the support sets the self-merge already tracks are the data it
  needs. This is the residual 2 dup-poles in the calibration.
- **Recall / walk coverage.** The bottom-edge misses are a path-network density
  matter (the walk visits those rows sparsely), not an extraction bug. A
  density-aware walk or per-figure node seeding is the lever; lowering the
  `min_support` default is not (it adds aliases). Out of scope for `fzse`.
- **Vector twin parity.** `VectorPoleField` uses a scale-derived *greedy*
  tolerance (`CLUSTER_FRAC·min(h)`) but no post-pass single-linkage self-merge;
  porting this disjoint-support post-pass there is a follow-on if a dense vector
  pole field over-splits.

## References

- `src/PoleField.jl` — `_extract_poles_core` (the self-merge + `_single_linkage_merge`).
- `test/polefield_test.jl` — PF.5.1 + the M7/M8 mutation footer.
- `external/probes/fzse-calibration/probe.jl` — the calibration (re-runnable).
- `figures/demo_lattice_singularities.jl` — the discovery demo (workaround removed).
