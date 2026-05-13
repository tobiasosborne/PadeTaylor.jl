# Worklog 012 — Phase 9: PI tritronquée pole-field qualitative reproduction

**Date**: 2026-05-13
**Author**: Claude Opus
**Scope**: Phase 9 / bead `padetaylor-kvi` — Tier-C qualitative
reproduction of the PI tritronquée pole-field structure per DESIGN.md
§4 Phase 9 and FW 2011 §4.1.  "No new code" phase per DESIGN —
purely a composition test of PathNetwork (Phase 10) + EdgeDetector
(Phase 12.5).  1160/1160 GREEN (1148 prior + 12 new).

> Take-home: a clean 12-assertion qualitative test that pins all four
> defining invariants of Boutroux's tritronquée — 4-of-5 pole-free
> sectors, conjugate symmetry, leading-pole magnitude, and
> discriminative-vs-null-IC sensitivity — for under 3 seconds of
> suite time.  No new code; the test IS the deliverable.

## Why this works at all — what made Phase 9 "no new code"

DESIGN.md §4 was written before PathNetwork (Phase 10) and BVP
(Phase 11) shipped; at the time, Phase 9 was sketched as "Run PI
tritronquée at FW eq. (4.1) ICs, generate pole-field plot, compare
visually."  With Phase 10 + Phase 12.5 (EdgeDetector) now in place,
the same deliverable becomes a **composition test** with no novel
implementation: feed PathNetwork the tritronquée problem on a 2D
Cartesian grid, reshape its flat output to a matrix, run EdgeDetector
on the result, assert qualitative invariants.

The visual-comparison-to-FW step is replaced by the four numeric
invariants in §"Test plan" below — invariants that capture the
*algebraic* content of Boutroux's tritronquée classification (5-fold
sector structure, one pole-bearing sector centred at angle 0) and
the *physical* properties of real ICs (conjugate symmetry, leading
real-axis pole magnitude).

## A correction to my own initial assumption

My first orientation probe placed the pole-bearing sector around angle
180° (Boutroux's conventional pole-bearing direction = "the sector
opposite the asymptotic ridge").  This was wrong for **this**
formulation of PI.  Empirically:

  - For `u'' = 6u² + z` with the FW eq. (4.1) ICs, the leading pole
    on the real axis is at `z ≈ 2.07` (positive real), NOT at
    `z ≈ -2.something` (negative real).
  - The pole-bearing sector is centred on **angle 0** (positive real
    axis), width ≈ 60°-72° (one of the 5 sectors of width 2π/5).
  - The pole-free 4-sector contiguous arc spans ≈ 288° going through
    72°, 144°, 180°, 216°, 288°.

This is consistent with FW Table 5.1's pole locations (all on the
positive real axis, increasing in magnitude with z).  The Boutroux
classification labels which sector is pole-bearing depending on
ODE normalisation; the FW formulation places it on positive real.

The lesson: **when a published convention disagrees with empirical
output, probe before adjusting the test.**  Saved 30 minutes of
"my test is wrong" hand-wringing by dumping the actual mask as ASCII
with explicit (x, y) labels.

## Test plan

`test/phase9_tritronquee_test.jl` ships 5 testsets / 12 assertions:

  - **PT.1.1**: composition completeness + determinism.  No NaN sentinels
    in `grid_u` (full PathNetwork coverage); mask non-empty (`≥ 20`
    flagged cells out of 529 interior); mask `< (N-2)²/2` (not majority);
    second call with `rng_seed = 0` (default) returns the identical mask.
  - **PT.1.2**: 4-of-5 pole-free sectors property.  Histogram mask cells
    over 12 × 30° angular bins; at least 8 bins empty (≥ 4·72° = 288°
    pole-free, minus a bin or two of binning slack);
    `count[bin_around_0°] / total ≥ 0.75` (majority concentration in the
    pole-bearing sector centred at angle 0); `count[bin_at_180°] == 0`
    (the negative-real-axis bin must be empty).
  - **PT.1.3**: conjugate symmetry under `y → −y`.  Reflect the mask
    across the real axis; tolerate ≤ 4 mismatched cells (path-network's
    wedge walker can produce slightly different paths from `z₀` to
    conjugate-pair grid cells, leading to a few cells near the flagged-
    region boundary differing in their `log₁₀|Δu|` value).
  - **PT.1.4**: leading positive-real-axis pole magnitude.  The grid
    cell at `(x, y) ≈ (2.33, 0)` (closest to the pole at `z ≈ 2.07`,
    just past it) has `|u| > 100` (empirically ≈ 387); the corresponding
    mask cell is `true`.
  - **PT.2.1**: discriminative — null IC `(u(0), u'(0)) = (0, 0)`
    produces a *different* sector concentration than the tritronquée IC.
    For tritronquée: `centre_two_share ≥ 0.75`.  For null: `< 0.70`
    (empirically ≈ 0.46).  This pins the test as sensitive to the
    actual problem state, not just any IC.

Test infrastructure: 25×25 grid over `[-4, 4]²` at `h_grid = 1/3`,
PathNetwork at `h = 0.5`, `order = 30`.  Path-network produces 231
visited nodes; full grid (625 points) covered without NaN.  Total
wall ≈ 2.1 s (including the discriminative null-IC re-solve at PT.2.1).

## Mutation-proof (verified 2026-05-13)

**Mutation C** — flip sign of `z` in the ODE: `6u² - z` instead of
`6u² + z`.  This is a different equation (a real-`z`-shifted Painlevé
variant) with a rotated pole pattern.  Verified bite: **7 fails** out
of 12 assertions —

  - PT.1.1: the "non-empty pole field" assertion still passes (poles
    exist somewhere); the "second run equals first" still passes (still
    deterministic).  But the "no NaN" assertion FAILS — the perturbed
    ODE drives some grid cells to unreachable in the path-network's
    `max_steps_per_target` budget, generating NaN sentinels.
  - PT.1.2: all 4 sector assertions fail.  The centre-two bins (around
    angle 0) no longer dominate; `count[bin_at_180°] == 0` fails because
    the rotated pattern places poles near the negative real axis now.
  - PT.1.4: `|u(2.33, 0)|` drops below 100 (the pole has moved); the
    `mask[j_x, i_real] == true` assertion also fails.

The remaining 5 assertions (3 in PT.1.1, the symmetry test PT.1.3, the
null-IC discriminative test PT.2.1) survive because they capture
structural properties (full coverage, determinism, conjugate symmetry,
discriminative-vs-null-IC) that are not specific to the tritronquée
sector pattern.  This is a clean bite — 7 of the 12 assertions are
specifically tritronquée-shaped, and Mutation C breaks them all.

**Mutation D** (documented but does NOT bite, as expected) — perturb
the tritronquée IC by 1e-2 (1000× the FW eq. 4.1 uncertainty):
`u_tri = -0.1875543… - 0.01`.  This is a near-tritronquée (FW Fig 3.1
territory).  Within `[-4, 4]²`, the 4-of-5 pole-free property still
holds approximately — pole fields exist in the other 4 sectors per
FW 2011 line 139 but are pushed beyond the window radius.  The test
correctly passes here, confirming it asserts *qualitative tritronquée
structure*, not bitwise-exact ICs.  This non-bite IS the test's
robustness contract: the test catches "wrong equation" but tolerates
"slightly wrong tritronquée ICs" — exactly the discrimination
boundary you want for a Tier-C qualitative test.

Mutation C restored before commit per CLAUDE.md Rule 4.

## What changed

`test/phase9_tritronquee_test.jl` (165 LOC) — the new test file.
`test/runtests.jl` — `include("phase9_tritronquee_test.jl")` added.
`docs/figure_catalogue.md` — Fig 3.1 row updated to "PARTIAL — Phase 9
ships qualitative reproduction".  FW's full 161² Fig 3.1 acceptance
deferred to `padetaylor-k31` (Phase 12 v2 2D dispatcher).

No source-file changes (Phase 9 is "no new code" per DESIGN.md).

## Frictions surfaced

  - **F1. Grid-orientation confusion.**  My first probe used
    `[x + im*y for y in ys for x in xs]` (Julia comprehension nested-
    loop semantics; outer `y`, inner `x`).  Reshaping the flat result
    column-major gave `u_grid[i, j] = u(xs[i] + im·ys[j])` — row=x,
    col=y.  But my ASCII visualisation header read "top-row = high y"
    — a copy-paste from a previous attempt.  Wasted ~5 minutes
    interpreting the wrong orientation.  Fixed by switching to a
    Matrix comprehension `[xs[i] + im·ys[j] for i in 1:N, j in 1:N]`
    + `vec(...)`, which makes the row/col semantics explicit in the
    construction.  **Lesson**: when reshaping flat output from a
    Stage-2 driver to a 2D matrix, prefer an explicit Matrix
    comprehension + `vec` over nested `for ... for` flat comprehension.
  - **F2. The PadeTaylorProblem zspan non-degeneracy guard.**
    `path_network_solve` only needs `zspan[1]` for the IC; `zspan[2]`
    is unused (the grid drives Stage 1).  But `PadeTaylorProblem`'s
    constructor rejects `zspan[1] == zspan[2]`.  Workaround: pass a
    far-corner placeholder for `zspan[2]`.  This is mildly awkward
    API-wise but not blocking; bead-worthy candidate for a future
    "drop the zspan constraint for path-network-only callers"
    refinement.
  - **F3. The Boutroux convention disagreed with empirical output**.
    See §"A correction to my own initial assumption" above.  Caught by
    probing before writing assertions; saved by the simple expedient
    of dumping the mask with explicit (x, y) labels before deciding
    which sector "should" be pole-bearing.

## Pointers

  - [`test/phase9_tritronquee_test.jl`](../../test/phase9_tritronquee_test.jl)
    — 5 testsets + mutation-proof commentary at bottom.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    lines 52-54 (5-fold sectorial structure), 139-145 (Fig 3.1 +
    caption), 214-229 (§4.1 tritronquée + eq. 4.1 ICs).
  - [`docs/figure_catalogue.md`](../figure_catalogue.md) §1 row Fig 3.1
    — now marked "PARTIAL" with Phase 9 acceptance noted.
  - [`docs/worklog/011-edge-detector.md`](011-edge-detector.md) — the
    prerequisite EdgeDetector phase.

## Bead state

Closed in this session:
  - `padetaylor-kvi` — Phase 9 PI tritronquée qualitative GREEN at this
    commit.

Still open from prior sessions (no change):
  - `padetaylor-rgp` (next-up, scheduled to close as living-document
    tracking),
  - `padetaylor-k31` (next compositional phase: Phase 12 v2 2D dispatcher),
  - `padetaylor-bvh`, `padetaylor-grc`, `padetaylor-61j`, `padetaylor-8pi`.

## Hard-won lesson (for HANDOFF.md §"Hard-won")

**15. Probe orientation before writing test assertions on 2D-grid
output.**  My first attempt at Phase 9's sector test asserted the
pole-bearing sector was at angle 180° (the "Boutroux convention" I
remembered from PI theory).  Empirically it's at angle 0° for the
FW formulation `u'' = 6u² + z`.  The disagreement was caught in
seconds by dumping the mask as ASCII with explicit (x, y) axis labels,
but I'd have spent 30+ minutes patching test assertions to a wrong
mental model.  When 2D output disagrees with theory, dump first,
theorise second.
