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

## Dense-field correction (bug `padetaylor-o0w4`, 2026-06-19)

**Status of this section**: accepted — extends the decision above. Resolves
`padetaylor-o0w4` (P1, the suite was RED on `test/ffw_fig_1_test.jl` FF1.1.6).

### The bug the original single-linkage introduced

The self-merge above was specified and shipped as **single-linkage**: a
union-find over the "close (`≤ merge_atol`) AND disjoint-support" edge graph.
Single-linkage **chains transitively** — `A~B` and `B~C` merge `{A,B,C}` even
when `|A − C| ≫ merge_atol`. On the equianharmonic-℘ calibration that is benign
(its froth sits within `merge_atol`; its genuinely-distinct lattice poles are
spaced `2ω ≈ 2.73 ≫ merge_atol = h_max`, so no chain bridges two real poles).

It is **not** benign on a *dense* pole field whose genuinely-distinct poles sit
**closer than `merge_atol`**. The FFW 2017 Fig 1 PIII three-sheet spiral
(`test/ffw_fig_1_test.jl`) has exactly that in its high-`Re ζ` region: FFW md:72
says "pole density will increase rapidly on the region `Re ζ ≫ 0`", so adjacent
distinct poles there are `< h_max` apart. Single-linkage chained a whole row of
distinct poles into a handful of representatives, **erasing the pole-density
gradient** the figure is built to show.

Measured (`extract_poles` at `cluster_atol = 0.15, min_support = 2`, sheet 0):

| variant | FF1.1.5 total | FF1.1.6 high `Re∈[3,5]` | FF1.1.6 low `Re∈[-1,1]` | FF1.1.6 |
|---|---|---|---|---|
| pre-fzse (no self-merge) | 2259 | 374 | 18 | PASS |
| fzse single-linkage (shipped) | 178 | **6** | **10** | **FAIL** (gradient inverted) |
| o0w4 diameter-capped (this fix) | 870 | 107 | 12 | PASS |

Single-linkage over-collapsed `2259 → 178` overall *and* inverted the gradient
(`6 < 10`). The `2259 → 870` de-frothing is still real and desirable; only the
dense-region over-merge was the defect.

### Decision: diameter-cap the single-linkage transitive closure

Keep the disjoint-support **single**-linkage edge graph — it is load-bearing and
must stay transitive (see below) — but **reject any union that would push the
merged component's diameter past `merge_atol`**. Edges are processed in
increasing distance; a union is accepted only if every member-to-member distance
in the would-be-merged component stays `≤ merge_atol`. This bounds each merged
group's diameter at `merge_atol`:

- **Over-split froth of one pole** has diameter `≪ merge_atol` (its cross-node
  spread is precisely what the merge targets — measured ≤ ~0.35 on the ℘
  lattice at `merge_atol = h_max = 0.5`), so the cap never blocks it: the
  calibration is **unchanged** (`n = 24, ndup = 2`, PF.5.1 still GREEN).
- **A chain of distinct dense poles** spans `> merge_atol` end-to-end, so the
  diameter test rejects the bridging union and the poles stay distinct: FFW
  Fig 1's gradient is preserved.

### Why NOT pure complete-linkage (the obvious alternative, and why it fails)

The first attempt was *complete*-linkage (merge a set only when **all** pairwise
distances `≤ merge_atol` AND all pairs disjoint-support). It capped the diameter
correctly and fixed FF1.1.6, but it **regressed PF.5.1** (`ndup(merged) = 3`,
the calibrated value is 2). Root cause, traced on the ℘ lattice
(`external/probes` ad-hoc, near `z = -0.363 + 2.36i`): a froth cluster's
**dominant** representative is seen by *many* nodes, so it **shares support**
with most of its own fragments. A fragment therefore reaches the dominant rep
**only through a disjoint intermediary fragment** (single-linkage's transitive
route, e.g. `15 ~disjoint~ 23 ~disjoint~ 24 ~disjoint~ 9`). Requiring *all-pairs*
disjoint orphans such a fragment, leaving a residual lattice duplicate. The
support-dimension linkage must stay **single** (transitive); only the
**distance** dimension takes the diameter cap. The shipped fix is therefore
*diameter-capped single-linkage*, not complete-linkage.

### Mutation-proof (Rule 4)

`test/ffw_fig_1_test.jl` FF1.1.6 is the load-bearing assertion. Removing the
diameter cap (forcing `within = true` ⇒ unbounded single-linkage) reproduces the
bug **exactly**: FF1.1.5 = 178, FF1.1.6 high/low = 6/10, both FF1.1.6 asserts
RED. Restoring the cap returns FF1.1.5 = 870, high/low = 107/12, GREEN. (The
prior M7/M8 mutations — self-merge OFF, disjoint-support dropped — remain valid
against PF.5.1 and the coalescent-pair fixtures; the diameter cap is a third,
orthogonal, load-bearing knob.)

### Scope of the change

`src/PoleField.jl` only: `_single_linkage_merge` → `_diam_capped_merge` (the
diameter-cap on the union step), the call site rename, and the paired docstring
updates. No test tolerance was relaxed (Rule 2). All FF1.* / PF.* / corpus
coalescent-pair fixtures (CPN.4, CRic.3, CRic.4) and the ffw_fig_4 / ffw_fig_6
pole-gradient tests stay GREEN. The disjoint-support guard and `merge_atol = 0`
legacy escape hatch are unchanged.

## References

- `src/PoleField.jl` — `_extract_poles_core` (the self-merge + `_diam_capped_merge`).
- `test/polefield_test.jl` — PF.5.1 + the M7/M8 mutation footer.
- `test/ffw_fig_1_test.jl` — FF1.1.5 / FF1.1.6 (the `padetaylor-o0w4` regression).
- `external/probes/fzse-calibration/probe.jl` — the calibration (re-runnable).
- `figures/demo_lattice_singularities.jl` — the discovery demo (workaround removed).
