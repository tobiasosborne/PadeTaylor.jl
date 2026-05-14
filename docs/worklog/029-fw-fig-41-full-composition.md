# Worklog 029 — FW 2011 Fig 4.1 full reproduction (steps ii + iii)

**Date**: 2026-05-14 (continues worklog 028)
**Author**: Claude Opus
**Bead**: `padetaylor-gky` — closed.
**Scope**: Compose the FW Fig 4.1 three-step BVP+IVP+BVP recipe into a
single figure. Step (i) — the imaginary-axis BVP — was pinned in
worklog 016 (`test/fw_fig_41_test.jl`, bead `padetaylor-0c3`); this
worklog adds steps (ii) the pole-field run-outs and (iii) the smooth-
band BVP fill.

> **Take-home**: `figures/fw2011_fig_4_1.jl` ships FW's "near-
> tritronquée case calculated by the steps (i)…(ii)…(iii)" (md:220).
> Step (ii) is `edge_gated_pole_field_solve` (worklog 028) run from
> each of the two BVP-derived ICs; step (iii) is `bvp_solve` per grid
> row, bridging each smooth run between its flanking anchors. All 119
> smooth-run BVPs converged with zero fallbacks. Test suite **1521 →
> 1526 GREEN** (+5 from `fw_fig_41_test.jl`'s new `FF.2.*` testsets).

## Ground truth read first

- `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:214-229` — §4.1 + the Fig 4.1
  caption: "(i) use a BVP solver over [-20i,20i] to obtain values for
  u and u' at z = 0 and z = 20i (ii) use this data to separately run
  out the two pole fields (iii) fill in the area between pole field
  edges by further BVP solutions." `:401` — smooth regions are
  unstable for any IVP solver (the reason step (ii) must be the
  edge-gated solve, not a plain `path_network_solve`). `:190` — the
  per-grid-line BVP fill pattern.
- The figure image `_page_7_Figure_2.jpeg` — a `|u(z)|` surface over
  `x ∈ [-10,10]`, `y ∈ [-10,30]`: a vast smooth plain, the upper pole
  field (back-left, run out from z=20i), the front-right pole field
  (run out from z=0), and the solid imaginary-axis line marking the
  step-(i) BVP segment.
- `test/fw_fig_41_test.jl` + `docs/worklog/016-fw-fig-41-pin.md` — the
  step-(i) recipe (`-√(-z/6)` negative-branch BCs, N=240) and the
  honest note that the quantitative tritronquée pole-*location* oracle
  is not in-tree.

## What shipped

### `figures/fw2011_fig_4_1.jl` (new)

The three steps, each mapped to a shipped component:

  - **(i)** `bvp_solve` on `[-20i, 20i]` — the worklog-016 recipe.
    Recovers `u(0) = -0.187554308340…` (FW eq. 4.1) and `u(20i)`,
    `u'(20i)`, plus the smooth solution along the imaginary axis.
  - **(ii)** two `edge_gated_pole_field_solve` runs (worklog 028) —
    one from the IC `(0, u(0), u'(0))`, one from `(20i, u(20i),
    u'(20i))`. The edge-gating is *load-bearing here*: FW md:401 says
    a plain IVP integrated across the near-tritronquée's smooth plain
    is corrupted; the gated solve keeps each run-out inside its pole
    field. Result: a 364-cell pole field from z=0, a 713-cell one from
    z=20i, on the 41×81 lattice.
  - **(iii)** per-grid-row `bvp_solve` of the smooth region. Each row's
    maximal smooth run is bridged between its flanking anchors — a
    pole-field-edge cell's IVP value, the step-(i) BVP spine value
    where the row crosses the imaginary axis, or the leading-term
    asymptote at a grid edge. **119 runs, all converged, zero
    leading-term fallbacks.**

The render is a `|u(z)|` surface (the shared `pole_field_figure`
helper, given a non-square `aspect` for the 1:2 domain) with the
imaginary-axis BVP segment drawn on top. Visual match to FW's Fig 4.1:
smooth plain + the two pole fields in the right places + the BVP
spine.

### `EdgeGatedSolution.u_grid` (new field)

`edge_gated_pole_field_solve` now also returns the final solve's `u`
scattered onto the `(nx, ny)` lattice (`NaN` outside the field+1-ring)
— the convenient 2D view the figure's step-(iii) stitching needs.
Additive; the worklog-028 tests are unaffected.

### `test/fw_fig_41_test.jl` — `FF.2.*` (new, 5 assertions)

There is no in-tree quantitative oracle for PI tritronquée pole
*locations* (FW Table 5.1 is the equianharmonic ℘, not PI — worklog
016). But the composition admits oracle-free cross-validations, and
those are what `FF.2.*` pin:

  - **FF.2.1** — the step-(ii) IVP run-out and the step-(i) BVP spine
    agree on the imaginary axis near `z = 0` to ≤ 1e-6. Two genuinely
    independent methods — spectral BVP vs. Taylor–Padé IVP — on the
    same near-tritronquée solution.
  - **FF.2.2** — the step-(ii) run-out from a real IC is conjugate-
    symmetric near `z = 0` (real ODE coefficients + real IC).
  - **FF.2.3** — the step-(iii) fill primitive: `bvp_solve` on an
    imaginary-axis sub-segment, given exact spine values as BCs,
    reproduces the spine in the interior to ≤ 1e-8.

Mutation-proven — M1 (shift the spine evaluation point off the
comparison), M2 (compare a cell to itself instead of its conjugate
reflection), M3 (feed the step-(iii) fill crude leading-term anchors
instead of exact spine values); each confirmed RED, reverted.
Procedure in the test file footer.

## Frictions

  - **The step-(iii) anchor logic, not the BVP, is the substance.**
    `bvp_solve` already exists and is tested; what was new was
    deciding *what BC each smooth run gets* — flanking IVP cell, spine
    crossing, or grid-edge leading term. This is figure-specific
    composition (two ICs, a vertical BVP spine, this geometry), so it
    lives in the figure script, not a `src/` module — consistent with
    worklog 016's "encapsulating it in a wrapper is premature
    abstraction; no downstream caller exists".
  - **Julia soft-scope.** The step-(iii) per-row loop counters tripped
    the top-level-`for` soft-scope ambiguity; wrapped the loop in a
    `let` block so the counters are loop-local.
  - **`bvp_solve`'s callable analytically continues off-segment.**
    Caught while writing the FF.2.1 mutation: `spine(im*y + 1.0)`
    does *not* throw — `real(t*)` stays in `[-1,1]` so the barycentric
    interpolant evaluates at a complex parameter (the analytic
    continuation off the axis). The mutation still bites (the
    continued value differs by ~0.1–0.5), but it is a reminder that
    the BVP callable's `DomainError` guard is on `real(t*)` only.

## What remains genuinely out of scope

The quantitative tritronquée pole-*location* pin. As worklog 016
recorded, FW Table 5.1 is the equianharmonic ℘, not PI; an honest PI
tritronquée pole oracle would need FW's Maple computation or an
independent recomputation, neither in-tree. Phase 9
(`test/phase9_tritronquee_test.jl`) covers tritronquée pole-field
*structure* qualitatively, and `FF.2.*` now cross-validate the Fig 4.1
*composition* — together that is the honest acceptance for this bead
(structural + cross-consistency), with the numeric pole-location pin
left explicitly unclaimed.

## Hard-won lesson

**Oracle-free cross-validation is available even when an oracle is
not.** The bead flagged "the hard part is the QUANTITATIVE oracle" and
worried steps (ii)+(iii) could only get "structural-only acceptance".
But a multi-method composition validates itself: the BVP spine and the
IVP run-out are independent computations of the *same* solution, so
their agreement where they overlap (FF.2.1) is a real quantitative
pin — no external oracle needed. FW themselves suggest exactly this
kind of self-consistency check (md:304-310). When you can't get an
oracle, look for two independent paths to the same number.
