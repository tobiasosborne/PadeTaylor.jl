# ADR-0034 — Bounded-window composite solve for the Fig 4.7 path-dependence seam

**Status**: **Accepted** — 2026-06-30.  Resolves bug `padetaylor-vwgl` (P0,
the Fig 4.7 seam); implementation bead `padetaylor-xwzf`; regression gate
`padetaylor-sny7`; verification bead `padetaylor-ingn`.
**Date**: 2026-06-30
**Worklogs**: 077 (root cause), 078 (scoping, mechanism selection, confirmation).

## Context

`path_network_solve` builds one monolithic walk tree from the IC across the
entire grid.  The Stage-1 target order is randomised (`rng_seed` shuffles
targets), so the branch structure — which cell is whose parent, which trunk a
far cell hangs off — is seed-dependent.  Two seeds produce two different
trees; where their subtrees meet, independently-accumulated IVP errors
disagree: a **grain-boundary seam** in the reconstructed pole lattice.
Worklog 077 root-caused it and showed the discontinuity covers ~35 % of the
`[-50,50]²` plane at up to O(1) relative error.

Bead `padetaylor-vwgl` (P0) — the seam visible in FW 2011 Fig 4.7 — is the
manifestation in the figures.  FW acknowledge the artefact in
`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:304-310`:

> "because of the random nature of our path selection algorithm … the
> difference between the upper and lower half-plane solutions gives a
> reasonable estimate for the numerical error."
> — FW2011 md:310

Their stated remedy — run twice and compare as an error indicator — confirms
they expected seed-to-seed variation; they do not claim the monolithic walk
is free of it.

**P3 measurement (probe `p3_edgegated_seedinvariance.jl`)** settled the
critical structural question — which regime is the seam — by re-solving the
edge-gated confined pole field at two seeds:

| metric (seeds 0 / 42) | plain full grid | edge-gated, confined |
|---|---|---|
| pole count | 8159 / 7987 | 8137 / 7966 |
| **pole-count seed-variance \|Δ\|** | **172** | **171** |
| pole-set match @0.5 | 76.6 % | 81.5 % |

**Verdict: Regime B** — a pole-field-*interior* grain boundary.  Confining
the walk to the pole field (`edge_gated_pole_field_solve`, already shipped
since ADR-0017, already active in `fw2011_fig_4_7.jl:106`) leaves |Δ| at
172 → 171: essentially unchanged.  The seam is *inside* the pole field; the
cure must act there.

## Decision

Adopt FW 2011's own construction (md:147) as the cure: tile the domain into
overlapping **bounded windows**, solve each **independently** from the shared
IC, and **composite by Voronoi-core ownership**.  FW describe their Fig 4.7
production method as:

> "each pole field shown in Fig. 4.7a–e is a 5×5 composite of 25 completely
> independent runs (each starting from the same ICs at z = 0, and covering
> areas of the extent 20×20 in the complex plane).  Similarly, Fig. 4.8 is
> a composite of 18 independent runs."
> — FW2011 md:147

This is productionised as the new module `src/WindowedComposite.jl` with
three public names:

- `windowed_path_network_solve(prob, xs, ys; window_extent, overlap, ...)`
  — the composite solver.
- `windowed_extract_poles(wsol; merge_atol, ...)` — composite pole field
  from the output.
- `WindowedCompositeSolution{T}` — the result struct carrying the composited
  field, per-window solves, Voronoi assignments, and per-window seeds.

The module is **additive**: it wraps `PathNetwork.path_network_solve`;
the default `path_network_solve` path is unchanged and all existing tests
stay green byte-for-byte.

Defaults reproduce FW's recipe: `window_extent = 20.0`, `overlap = 6.0`.

The P0 bug (`padetaylor-vwgl`) is gated by the seed-invariance regression
test `test/field_seam_test.jl` (bead `padetaylor-sny7`): assert composite
pole-set seed-invariance (|Δcount| + Hausdorff) over two global seeds; RED
on the monolithic walk (|Δ|=151 / 77 %), GREEN on the composite (|Δ|=12 /
99.4 %); tolerance pinned at the measured composite floor, not relaxed.

## Why the composite works — and not the alternatives

**The mechanism.** Every cell in one window hangs off the *same* short trunk
back to the shared IC, so the cells' *differential* positions are seed-stable
even though the trunk itself varies with the seed.  The cross-branch
differential that produced the seam never forms — there is no second subtree
to disagree with.  Worklog 078 §3a confirmed this even for far/long-trunk
windows:

| core region | monolithic \|Δ\| / match | windowed \|Δ\| / match |
|---|---|---|
| NEAR `[-10,10]²` (short trunk) | 9 / 97.9 % | **2 / 100 %** |
| FAR `[18,38]×[-10,10]` (long trunk) | 13 / 100 % | **1 / 99 %** |

The FAR window collapses to |Δ|=1, refuting the objection that long trunks
are as bad as the monolithic walk.  All cells in one window share the trunk
(common mode); only the *differential* matters, and within a window it is
zero.

**Why reconciliation is not the primary fix.**  Probe `p4b_overlap_graph.jl`
measured the disc-overlap graph (worklog 078 §3b):

```
17776 nodes · 48529 disc-overlap edges · avg degree 5.5
CROSS-BRANCH overlaps (nearest common ancestor >16 hops):
    global 1.7%    seam annulus 1.9%
```

~98 % of disc overlaps are same-branch — already mutually consistent.  Only
~1.9 % carry cross-branch information.  A global least-squares built on this
graph is dominated by the same-branch bulk (the unstable near-null-space); the
seam is weakly constrained and robust convergence requires regularisation —
which is smoothing, the forbidden non-fix (see Alternatives).

**Why smooth-region fill (`fe9`/`8dg`) is not the fix.**  Both families cure
the *field* (Fig 4.1), not the *pole scatter* (Fig 4.7): they produce no
poles and cannot move the pole lattice.  P3's regime-B verdict — the seam is
inside the pole field — eliminates them as cures for the seam.

## Per-window seeding is load-bearing (non-gameable gate)

Each window is solved at `_window_seed(rng_seed, wi)`:

```julia
_window_seed(rng_seed::Integer, wi::Integer) = Int(rng_seed) * 1_000_003 + wi
```

This is **arithmetic, not `hash`**: `hash` is not stable across Julia
versions and would silently break the package's bit-reproducibility contract.
The `1_000_003` multiplier (a prime above any realistic window count) keeps
`(rng_seed, wi)` pairs collision-free for sane inputs.

The non-gameability property is essential.  The `sny7` gate re-runs the
composite at two global seeds and asserts the pole field barely moves.  If
windows ignored the global seed (e.g. always seed 0), the gate would pass
trivially — the composite would be seed-independent because nothing
re-randomised.  Threading the global seed into every window means two global
seeds re-randomise every window's internal walk, so a gate that passes
certifies genuine path-independence, not a frozen pipeline.

## Confirmed numbers (P5, full-composite acceptance probe)

Probe `p5_full_composite_confirm.jl` — a 5×5 composite over `[-50,50]²`,
poles in `[-49,49]²`, two global seeds (worklog 078 §4):

| metric (seeds 0 / 42) | Monolithic | **Windowed composite** |
|---|---|---|
| pole count | 7443 / 7292 | 7223 / 7211 |
| **pole-count seed-variance \|Δ\|** | **151** | **12** |
| pole-set match @0.5 | 77.1 % | **99.4 %** |
| median NN displacement | 0.076 | **0.000** |
| wall time (both seeds) | ~214 s | **~102 s** |

Tile-boundary band match 98.9 % vs interior 99.2 % — no boundary seam; the
overlap-core partition handled it.

## Consequences

- **The seam is contained.** |Δ|≈12 is *containment*, not annihilation — a
  handful of poles near the noise floor still wobble.  This matches FW's own
  posture (they tolerate it and plot only poles); the `sny7` tolerance is
  pinned to the measured composite floor, not relaxed.

- **Composite finds ~220 fewer poles** than monolithic (7223 vs 7443).  These
  are almost certainly the spurious seed-dependent poles correctly suppressed —
  the count gap *is* the bug — but this is a **named verification**
  (`padetaylor-ingn`) via the Weierstrass-℘ truth anchor, not an assumption
  baked in.  Until `ingn` closes, the gap is an honest residual.

- **`window_extent` must be meaningfully smaller than the domain (several
  windows per axis).**  Gate calibration surfaced a hard usage caveat: at
  `window_extent ≈ domain` the composite is no better than monolithic
  (measured: HALF=25 gave composite 91 % ≈ monolithic 93 %; HALF=30 gave
  composite 97 % vs monolithic 78 %).  Figures use `[-50,50]` with the
  default 20.0 (5 per axis on each domain half), well inside the good regime.
  The `windowed_path_network_solve` docstring warns on this.

- **Default path unchanged.**  `path_network_solve` is not touched; the
  `@inferred`/allocation type-stability gates remain green; zero regression on
  the full test suite.

- **Rule 1 (fail loud) at two new callsites.** `windowed_path_network_solve`
  throws `ArgumentError` on malformed inputs (window captures no grid cell;
  non-increasing axes; `window_extent ≤ 0`; `overlap < 0`; fewer than 3
  nodes per axis) and throws `ErrorException` with window context if a per-
  window solve fails.  An unassigned cell after the composite loop throws
  (a coverage bug, not a silent NaN).

## Alternatives considered and rejected

**(a) Edge-gating alone.**  Already shipped (ADR-0017).  P3 measured
|Δ| 172 → 171: no cure.  The seam is Regime B (inside the pole field),
so surface-level confinement cannot remove it.

**(b) Smooth-region BVP/harmonic fill (`fe9`/`8dg`).**  Harmonically fill
the smooth regions from the pole-field boundary, or apply FW md:190 per-row
Chebyshev–Newton.  Both families cure the |u| *field* (Fig 4.1 use case);
neither produces poles, so neither can move the pole scatter.  Regime B
measurement eliminates them as cures for the seam (`padetaylor-fe9` /
`padetaylor-8dg` re-scoped to field-only utility).

**(c) Cross-branch overlap reconciliation / K-nearest consensus.**  Derive
node states from disc-overlap constraints (Gauss–Newton on the overlap graph)
or a K-nearest-path consensus at Stage-2.  Rejected: the cross-branch
constraint graph is only ~1.9 % of the overlap edges (p4b measurement).  A
global least-squares is dominated by the same-branch bulk (the unstable near-
null-space of the overlap consistency equations); robust convergence requires
regularisation, which is smoothing — a forbidden non-fix that masks the error
instead of removing it.

**(d) Freeze the shuffle / `enforce_real_axis_symmetry`.**  Make the walk
deterministic by sorting targets rather than shuffling, or enforce conjugate
symmetry post-hoc.  Rejected as forbidden non-fixes: they mask seed-
dependence rather than remove it.  A seed-invariance gate built on frozen
state passes trivially and certifies nothing about path-independence.

## References

- FW2011 md:147 — the 5×5 / 18-run composite construction this module
  productionises (the primary source the decision rests on).
- FW2011 md:208 — FW's Laplacian-residual edge detector; confirms the
  `"irregularities caused by different paths"` in smooth-region cells.
- FW2011 md:304-310 — FW's own admission of path-dependence and their error
  indicator (run twice, compare); confirms the monolithic walk was always
  known to be seed-variant.
  `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md`
- `docs/worklog/077-fig47-path-dependence-seam-rootcause.md` — root cause
  (the walk enforces no path-independence; ~35 % of `[-50,50]²` at O(1)
  rel error).
- `docs/worklog/078-fig47-seam-cure-scoping.md` — 5-family audition, P3/P4/
  P4b/P5 measurements, mechanism selection, confirmed numbers.
- `src/WindowedComposite.jl` — the implemented module (literate chapter;
  the ADR is the terser decision record).
- `test/field_seam_test.jl` — seed-invariance regression gate (bead
  `padetaylor-sny7`).
- `external/probes/fig47-seam-diagnosis/` — the four re-runnable probes
  (p3, p4, p4b, p5) whose measurements this ADR distils.
- ADR-0017 — edge-gating (`EdgeGatedSolve`) — the Regime-B precondition.
- Bead `padetaylor-vwgl` (P0, the bug); `padetaylor-xwzf` (impl);
  `padetaylor-sny7` (gate); `padetaylor-ingn` (verification: spurious vs
  real poles).
