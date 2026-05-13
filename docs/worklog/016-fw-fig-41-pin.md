# Worklog 016 — FW 2011 Fig 4.1 step-(i) BVP quantitative pin (bead `padetaylor-0c3`)

**Date**: 2026-05-13 (late evening, follow-up to worklog 015)
**Author**: Claude Opus
**Scope**: Reproduce the tritronquée FW eq. 4.1 reference values via
the FW Fig 4.1 step-(i) recipe — Chebyshev-Newton BVP for
`u'' = 6u² + z` on the vertical imaginary-axis segment `[-20i, +20i]`
with leading-term `u(z) = -√(-z/6)` Dirichlet BCs.  Test-only delivery;
no new `src/` code.  GREEN at 1250 → 1262.

> Take-home: the FW Fig 4.1 step (i) is a single `bvp_solve` call with
> the right BCs, initial guess, and N.  At N=240 the result pins
> `u(0)` to ≤ 3.5e-13 vs FW eq. 4.1 and `u'(0)` to ≤ 5.3e-11.  FW
> claims "better than 1e-20" at unspecified N; we hit the bead's 1e-10
> spec with N=240 and ≈100 ms wall time.

## What changed

  - `test/fw_fig_41_test.jl` (new, ~110 LOC) — 6 testsets / 12 assertions:
    - **FF.1.1** — Newton converged (≤10 iters, finite residual).
    - **FF.1.2** — `u(0)`, `u'(0)` match FW eq. 4.1 to ≤1e-10 (real
      parts); imag parts ≤1e-10 (Schwarz reflection on real axis).
    - **FF.1.3** — `u(+20i) = imposed BC` to ≤1e-13 (Dirichlet
      structural check on the barycentric callable).
    - **FF.1.4** — `u'(+20i)` sanity: finite, magnitude bounded
      (0.01 < |up| < 0.1), within 1e-3 of leading-term derivative
      `1/(12·√(-z/6))`.  Loose tolerance documents the o(1) correction
      explicitly — `u'(20i)` is NOT pinnable at 1e-10 against the
      leading-term derivative; an independent oracle would be needed
      for a tighter pin.
    - **FF.1.5** — Schwarz symmetry `u(-iy) = conj(u(iy))` on the BVP
      segment at 5 sample points to ≤1e-12.
    - **FF.1.6** — mutation-self-proof: the `+√(-z/6)` branch (wrong
      sign) drives Newton into a divergent basin (`@test_throws
      ErrorException`).  Embodies the structural BC-sign sensitivity.

  - `test/runtests.jl` — include the new test file.

  - `docs/figure_catalogue.md` row 4.1 — marked PARTIAL with the step-(i)
    pin shipped + the deferral note for steps (ii) and (iii).

  - No `src/` changes.  The recipe is the test itself; encapsulating
    it in a `fw_fig_41_axis_bvp(z_max; ...)` wrapper is premature
    abstraction (CLAUDE.md "Don't add ... abstractions beyond what the
    task requires") — no downstream caller exists, and the bead spec
    is a pin, not a feature.  The wrapper can be added later if needed.

## Spectral convergence probe

The N→err mapping at Float64 on this BVP (probe results, 2026-05-13):

| N | iter | residual_inf | err_u(0)       | err_u'(0)    |
|---|------|--------------|----------------|--------------|
|  20 |  8 | 3.1e-12     | 5.2e-2          | 1.1e-1       |
|  30 |  7 | 1.1e-11     | 2.5e-2          | 6.6e-2       |
|  40 |  7 | 2.3e-11     | 1.1e-2          | 3.6e-2       |
|  60 |  6 | 1.3e-10     | 1.3e-3          | 8.1e-3       |
|  80 |  6 | 4.3e-10     | 1.2e-4          | 1.3e-3       |
| 100 |  6 | 3.3e-10     | 9.3e-6          | 1.8e-4       |
| 120 |  6 | 2.0e-9      | 7.4e-7          | 2.4e-5       |
| 150 |  5 | 3.1e-9      | 1.8e-8          | 9.9e-7       |
| 180 |  5 | 1.2e-8      | 4.7e-10         | 3.9e-8       |
| 200 |  5 | 2.5e-8      | 4.2e-11         | 4.4e-9       |
| 220 |  5 | 2.0e-8      | 3.8e-12         | **4.9e-10**  |
| 240 |  5 | n/a         | **3.5e-13**     | **5.3e-11**  |
| 250 |  5 | 5.8e-8      | 1.0e-13         | 1.8e-11      |

N=240 is the smallest value at which BOTH `u(0)` and `u'(0)` meet the
≤1e-10 bead spec.  Wall time per BVP ≈ 100 ms (5-6 Newton iters,
each a 239×239 complex LU solve via LAPACK).

**Note the `residual_inf` floor that GROWS with N.**  This is the
spectral-truncation-error floor `cond(D₂) · eps(T) ≈ N² · eps(T)`
documented in worklog 006 §"Step-norm Newton".  At N=200, `N²·eps ≈
8.8e-12`, but `residual_inf` reports ~2.5e-8 — the floor reflects
the OVERALL residual at the converged Newton point (including the
non-trivial nonlinear part of the BVP), not just the linear-stage
truncation.  Newton CONVERGES (step-norm criterion ≤ `eps^(3/4)`);
the `residual_inf` is informational.  The test does not gate on it.

## Algorithmic findings

### BC sign drives basin selection

FW eq. 1.2 (FW2011...md:54-55) gives the leading asymptote as `u(z) =
±√(-z/6) + o(1)`.  The `±` corresponds to the FIVE asymptotic sectors
of PI (FW Fig 1.1, sector angles `2π/5`).  The tritronquée is unique:
pole-free in 4 of the 5 sectors, with sector boundaries at angles
`±π/5, ±3π/5, π` (the pole-bearing sector is centred on the positive
real axis).

For the `[-20i, +20i]` BVP, BOTH endpoints lie in pole-free sectors.
The TRITRONQUÉE branch at these points is `u(z) ~ -√(-z/6)` — verified
empirically: this branch's Newton converges in ≤8 iters and gives
`u(0) ≈ -0.188` (matching FW eq. 4.1's negative value); the `+√` branch
drives Newton into a divergent basin (`‖Δu‖ → 2.0+`, `‖R‖ → 4e4`)
within 10 iters.  FF.1.6 pins this.

The sign-choice is the load-bearing structural detail of the recipe.
Other PI solutions (tronquée cases, fully meromorphic solutions) have
different branch choices on the same axes — the BVP solver picks ONE
of them based on the BC sign + initial-guess basin.

### `u'(20i)` is not pinnable at 1e-10 against the leading term

The BVP is Dirichlet in `u`, not `u'`.  At `z = 20i`, `u(20i)` is
imposed (the BC); `u'(20i)` is the converged BVP's natural output and
differs from the leading-term derivative `1/(12·√(-z/6))` by `O(1)`
correction.  Empirical:

  - leading-term `u'(20i) = 1/(12·(1.291 - 1.291i)) ≈ 0.0323 + 0.0323i`
  - BVP-computed `u'(20i) ≈ 0.0324 + 0.0325i`
  - difference `≈ 2.5e-4` (the o(1) correction)

A 1e-10 pin against the leading term would FAIL.  An independent
oracle (mpmath, Mathematica `NDSolve` with high precision) could
deliver a tighter pin, but is out of scope for v1 — the bead's
canonical acceptance is `u(0), u'(0)` to ≤1e-10, and that is met.
FF.1.4 documents the gap explicitly via the looser ≤1e-3 assertion
and an order-of-magnitude bracket.

## What is NOT shipped (deferral notes)

**Step (ii) — outward pole-field IVPs from (0, u(0), u'(0)) and (20i,
u(20i), u'(20i)).**  These are `path_network_solve` calls with the
BVP-derived ICs.  Mechanically straightforward; the deferral is about
QUANTITATIVE pin criteria: the bead's "pole locations agree with FW
Table 5.1" mis-references (FW Table 5.1 tabulates the equianharmonic
℘-function, not PI tritronquée pole locations).  An honest tritronquée
pole-location oracle would come from FW's Maple computation or an
independent recomputation; neither is in-tree.  Phase 9
(`test/phase9_tritronquee_test.jl`) already covers tritronquée
pole-field STRUCTURE qualitatively (4-of-5 pole-free sectors).

**Step (iii) — interior BVP fill of the smooth band between pole field
edges.**  This is `lattice_dispatch_solve`'s line-190 algorithm (Phase
12 v2 v1).  Wiring it up for a SPECIFIC tritronquée 2D grid + showing
the band fill is a visualization deliverable, not a quantitative pin.

Neither (ii) nor (iii) gets a follow-up bead in this session.  Both
could be valuable for figure reproduction (Phase 12 v3, hypothetical)
but are not on the critical path for any current bead.

## Frictions surfaced

1. **The leading-term BC is NOT a tight pin at the endpoints.**  My
   first-instinct test plan was to pin `u(±20i)` AND `u'(±20i)` against
   the leading term `±√(-z/6)`.  This works for `u` (Dirichlet imposes
   it) but fails for `u'` at the 1e-10 level (o(1) correction).
   Documented the gap in FF.1.4 with a looser tolerance.

2. **`residual_inf` grows with N; Newton step-norm is the real signal.**
   First-cut test gated `sol.residual_inf ≤ 1e-11`; FAILED at N≥80.
   Removed the gate — the BVP uses step-norm Newton (worklog 006
   lesson 8), and the residual at convergence reflects truncation
   error not Newton stalling.

3. **N=240 is large but cheap.**  ~100 ms per BVP solve, dominated by
   one LAPACK LU.  No need to reach for BigFloat or N=400+ for the v1
   spec.  If a future bead asks for `≤1e-20` (FW's claim), bumping to
   BF-256 + N=300+ would deliver — empirical probe shows BF-256 and
   F64 give the SAME spectral convergence at small N (precision is not
   the bottleneck at this segment length / order).

## Beads

  - `padetaylor-0c3` — closed in this session.

## Pointers

  - `test/fw_fig_41_test.jl` — the canonical step-(i) recipe + 12 pins.
  - `src/BVP.jl` — the underlying Chebyshev-Newton solver; no changes.
  - `docs/figure_catalogue.md` §1 row 4.1 — PARTIAL acceptance marked.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    :48 (PI), :54-55 (leading-term asymptote), :214-222 (Fig 4.1
    composition), :226-227 (FW eq. 4.1: u(0), u'(0) reference).
  - `docs/worklog/006-phases-10-11-path-network-bvp.md` — BVP-solver
    Phase 11 implementation + step-norm Newton convergence rationale.
  - `docs/worklog/013-phase-12-v2-lattice.md` — original v1-scope
    decision that filed this bead.

## Hard-won lessons (for HANDOFF.md §"Hard-won")

23. **The leading-term Dirichlet BC pins `u` exactly but not `u'`**.
    For BVP step-(i) recipes that derive `u'` at the segment endpoints
    from the BVP (rather than imposing a Neumann BC), the converged
    `u'(z_endpoint)` differs from the analytical leading-term
    derivative by the o(1) correction.  At `z = ±20i` for tritronquée,
    this gap is ≈2.5e-4 — large enough to break a naive 1e-10 pin.
    Test the gap with a LOOSER tolerance documenting it explicitly,
    OR generate an independent oracle for a tighter pin.

24. **N-tuning probe before pinning rtol on a fresh BVP**.  Different
    segments (length, smoothness, orientation) require different N.
    `[-18, -14]` smooth band (Fig 4.6): N=20 suffices.  `[-20i, +20i]`
    long segment (Fig 4.1): need N=240 for ≤1e-10.  Probe the
    convergence empirically before committing the test's N parameter.
    The mapping is roughly `err ~ R^(-N)` for spectral convergence;
    `R` depends on the analyticity strip of the solution near the
    segment, which is problem-specific.

25. **`residual_inf` is NOT a convergence metric for spectral BVPs**.
    The BVP's step-norm Newton (worklog 006 lesson 8) converges in `‖Δu‖_∞
    ≤ eps^(3/4) ≈ 1.6e-12` while `‖R‖_∞` may sit at `N²·eps(T)` or
    higher.  Gating tests on `residual_inf ≤ TINY` is structurally
    wrong — Newton CAN converge with `‖R‖_∞ ~ 1e-8` at N=200 and
    still deliver `u(0)` to 1e-13.  Gate on what you actually want
    (per-eval error against an oracle) not on the residual.
