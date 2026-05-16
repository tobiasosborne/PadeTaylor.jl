# Loop-closure probe — FFW 2017 Fig 1 (sheet [0])

**Bead**: `padetaylor-e3h`.  Probe directory: `external/probes/loop-closure-fig1/`.
**Date**: 2026-05-16.  **Author**: Claude Opus 4.7 (orchestrated session).

## Hypothesis

The path-network walker (`src/PathNetwork.jl::path_network_solve`) builds
a **tree** rooted at the IC.  Each visited node's stored canonical Padé
inherits accumulated truncation error along its unique IC-to-node path;
the local step controller (`:adaptive_ffw`) sees only one step at a
time.  The user's hypothesis (2026-05-16): *every loop in the visited-
node graph should agree, modulo branch cuts.*  If true, loop-closure
disagreements at Delaunay-edge midpoints supply both consistency checks
and a correction signal.

This probe tests that hypothesis empirically on FFW 2017 Fig 1's
principal sheet (sheet [0]) ζ-plane solve.

## Setup

Re-runs the Fig 1 walker verbatim (`figures/ffw2017_fig_1.jl` parameters,
PIII `(α, β, γ, δ) = (-1/2, -1/2, 1, -1)`, IC `(z, u, u') = (1, 1/4, 1)`,
ζ-window `Re ζ ∈ [-2, 5]`, `Im ζ ∈ (-6π+0.5, 6π-0.5]`, `R(ζ) = max(0.1,
(8-Re ζ)/20)`, `:adaptive_ffw` at `tol = 1e-10`, `order = 30`).  Filters
to sheet [0] via ζ-strip predicate `-2π < imag(visited_z[k]) ≤ 2π`
(Fig 1 doesn't pass `branch_points`, so `visited_sheet[k]` is empty —
verified by `smoke.jl`).  Delaunay-triangulates the sheet-0 z-positions;
extracts non-tree (loop-closure) edges; computes the midpoint
disagreement `ΔP_rel := |P_A(M) - P_B(M)| / (|P_A(M)| + |P_B(M)| + ε)`
for each non-tree edge using `PadeTaylor.PathNetwork._evaluate_pade` on
both endpoints' stored Padés.

## Numbers

| metric | value |
|---|---:|
| visited tree nodes | 2664 (worklog 037 ≈ 2664 — exact match) |
| sheet-0 nodes | 850 |
| Delaunay edges (sheet 0, ghost edges dropped) | 2526 |
| tree edges (sheet-0 ∩ sheet-0) | 840 |
| non-tree edges (loop-closure population) | **1689** |
| `median(ΔP_rel)` | `2.57e-12` |
| `p90(ΔP_rel)` | `2.22e-05` |
| `p99(ΔP_rel)` | `1.06e-01` |
| `max(ΔP_rel)` | `9.97e-01` |
| edges with `ΔP_rel > 1e-10` (above adaptive_tol) | 438 (25.9 %) |
| edges with `ΔP_rel > 1e-6` (above FFW Exp-2) | 202 (12.0 %) |
| edges with `ΔP_rel > 1e-3` (visibly catastrophic) | **107 (6.3 %)** |
| Pearson r(`tree_dist`, log₁₀ ΔP_rel) | `+0.176` |
| Pearson r(`extrap_max`, log₁₀ ΔP_rel) | `+0.137` |
| Padé denominator blowups | 0 |
| Wall (walker re-run / edge-eval / total) | 16.5 s / 0.16 s / ~65 s |

The histogram is **trimodal**: a ~480-edge lobe at `10⁻¹⁶…10⁻¹³`
(machine-eps closure, conjugate-symmetry-quality), the dominant
~1000-edge lobe at `10⁻¹³…10⁻¹⁰` (closure within controller tolerance),
and a ~200-edge tail at `10⁻⁶…10⁰` (catastrophic).

## Spatial finding

Of the 202 "bad" edges (`ΔP_rel > 1e-6`):
- **65 %** sit at `Re ζ > 3.0` (the high-Re ζ corner where `R(ζ)` floors
  at `0.1` — the sparsest part of the sampling).  Compared to ~30 % of
  the area, this is a 2× enrichment.
- **78 %** sit at `|Im ζ| > π` (vicinity of the sheet boundary
  `Im ζ = ±2π`).  Compared to ~25 % of the sheet-0 strip area, this is
  a 3× enrichment.
- Several top-10 worst edges sit at `Re ζ ≈ 5` with `extrap_max > 5` —
  i.e., midpoints are in the Padé **extrapolation** regime `|t| ≫ 1`,
  not the canonical disc.

The low Pearson correlations (`r ≈ 0.18` for tree_dist; `r ≈ 0.14` for
extrap_max) are a consequence of the bimodal structure: a global linear
fit averages the well-closed interior bulk against the catastrophic
boundary tail.  The **clustering** is the stronger signal than the
correlations.

## Verdict — Partially confirmed

Loop-closure disagreement is **real, large, and spatially structured**,
but **not in the way the original hypothesis predicted**.

  - The hypothesis that "tree-depth causes the seam via accumulated
    error along long IC-to-node paths" is **NOT supported on sheet 0**:
    Pearson `r(tree_dist, ΔP_rel) ≈ 0.18` is weak, and the catastrophic
    tail is geometrically localised (high-Re ζ + sheet boundary), not
    distributed along the deepest tree branches.
  - The catastrophic 6.3 % of loop closures *do* exist, *are* well above
    visible-artefact threshold, and *do* cluster in a structured way —
    but the cluster is in the **undersampled high-Re ζ corner** where
    `R(ζ)` already loses resolution and stored Padés must extrapolate
    past their disc to reach midpoints.
  - The "consistency check" half of the user's intuition is **vindicated**:
    a graph-consensus pass would flag these 107 edges as suspect.  The
    "correction signal" half is **harder**: at `|t| > 5` neither Padé
    is trustworthy, so naive averaging doesn't fix the issue — denser
    sampling does.

**Sheet 0 itself looks clean visually** in `figures/output/ffw2017_fig_1.png`;
the visible seams the parent conversation flagged (the white wedges on
sheets ±1) are **out of scope for this probe**.  Likely those seams have
a different cause (sheet-lift artefact, bilinear-interp at strip
boundaries, or — possibly — a tree-depth effect that genuinely materialises
only on the deeper sheets).  A follow-up probe at sheets ±1 is required
before drawing architectural conclusions about Fig 1's visible seam.

## Recommendations for follow-up beads

1. **P2 — "Graph-consensus Stage-2, jointly with denser sampling at high-Re ζ."**
   The signal is real but spatially localised to the region where `R(ζ)`
   already loses resolution.  Consensus alone is insufficient at
   `|t| > 5`; denser node placement (Poisson-disk per FFW Fornberg 2015,
   already a deferred bead `padetaylor-zwh`) must land jointly so that
   midpoints stay inside `|t| ≤ 1`.

2. **P3 — "Loop-closure probe on sheets ±1."**  Re-run with
   `branch_points = (0.0+0.0im,)` and `cross_branch = true` so
   `visited_sheet[k]` can be read directly.  Likely (untested here):
   the wedges on sheets ±1 may be a path-network-depth artefact — long
   walks deep into the spiral arms accumulating error — distinct from
   the boundary effect on sheet 0.  If confirmed, that *would* validate
   the original hypothesis on the figure's visible artefact and motivate
   architectural change.

## Caveats

  - Sheet [0] only.  Sheets ±1 deferred (monodromy bookkeeping).
  - No disc cap applied to `_evaluate_pade`.  Several top-10 edges have
    `|t| > 5`; the disagreement at those edges combines honest tree-
    divergence with Padé extrapolation amplification.  A consensus-Stage-2
    design wouldn't trust either endpoint at `|t| ≫ 1`, so this is the
    right population to flag, not a measurement artefact.
  - Convex-hull edges: some long red chords in `spatial_map.png` span
    across the high-Re hull boundary; the midpoint of such edges falls
    outside the support of the walker tree.
  - `_evaluate_pade` is an internal API.  Prior art for accessing it
    from outside the package: `test/ffw_fig_1_test.jl:160-162` already
    does so; verified stable as of HEAD `1c5e353`.
  - `DelaunayTriangulation.each_edge` returns ghost edges (negative
    vertex indices for convex-hull boundary representatives); 21 of
    these were filtered before counting.  First run crashed on this
    until handled.

## Reproduce

```
cd external/probes/loop-closure-fig1
julia --project=. probe.jl
# Outputs: output/{histogram.png, spatial_map.png, summary.txt}
# Wall ~65 s (Stage-1 walk 16 s + Delaunay precompile 13 s + plot ~35 s).
```

`smoke.jl` runs the three "architectural unknown" verifications from the
plan (visited_sheet empty? `_evaluate_pade` callable? DelaunayTriangulation
available?).  All three pass on HEAD `1c5e353`.
