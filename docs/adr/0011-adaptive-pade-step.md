# ADR-0011 — Adaptive Padé step size (FFW 2017 §2.1.2)

**Status**: Accepted (2026-05-15) | **Bead**: `padetaylor-8ui` | **Worklog**: 034

## Context

The FW 2011 Padé-Taylor method ships in this repo with a *constant*
step size `h = 0.5` per `path_network_solve`'s default (ADR-0004).
That choice is FW 2011's own (`references/markdown/
FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:271-281`)
and it works for the single-valued Painlevé transcendents
(`P_I`, `P_{II}`, `P_{IV}`) on which FW 2011 / FW 2014 / FW 2015
/ RF 2014 ship their figures.

FFW 2017 turns its attention to the multivalued transcendents (`P_{III}`,
`P_V`, `P_{VI}`) on their Riemann surfaces, and observes
(`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:67-72`) that the pole
densities of these solutions grow exponentially with `Re ζ` under the
PIII/PV/PVI exponential transformations.  A uniform step `h = 0.5` is
then catastrophically wasteful at low `Re ζ` (orders-of-magnitude
over-resolution between widely-spaced poles) and dangerously
under-resolved at high `Re ζ` (where poles crowd to spacing ≪ h).

FFW 2017 §2.1.2 (md:74-97) introduces an adaptive step-size controller
specialised to the Padé-Taylor stepper.  Worklog 017 §"What is NOT
shipped" deferred this; bead `padetaylor-8ui` is the deferral
materialised, and this ADR records the design.  Step A1 of the 11-step
FFW-figure-reproduction plan; downstream work (non-uniform Stage-1
nodes, all seven FFW figures) depends on this.

## Decision

Ship an opt-in `step_size_policy = :adaptive_ffw` kwarg on
`PathNetwork.path_network_solve` that activates the FFW 2017 controller.
The default remains `:fixed` (backward compatible, ADR-0004 verbatim).
Three new helpers live in `src/PadeStepper.jl` as the implementation
core:

  - `ffw_truncation_error(f, z, u, up, n, h) → Real` — compute the
    FFW estimator `T(h) = |ε_{n+1} h^{n+1} / a(h)|` from one extra
    Taylor coefficient.  Algebra is in the rescaled variable
    `t = h'/h` (consistent with the rest of the stepper's `h^k`
    rescaling discipline, see `PadeStepper.jl` module docstring):
    `c̃_k = h^k c_k`, `ε̃_{n+1} = c̃_{n+1} + Σ_{r=1..ν} b_r c̃_{n+1-r}`,
    `a_rescaled(1) = a(h)`.  Return type is real-valued (the FFW
    controller magnitude is real even when `h` is complex).

  - `ffw_rescale_q(Tol, T_h, n; k = 1e-3) → Real` — the FFW eq. 2
    rescale factor `q = (k·Tol/T_h)^(1/(n+1))`.  Pure; defensive
    against `T_h = 0` (returns `Inf`).  FFW's recommended `k = 1e-3`
    is the default conservative factor (md:91).

  - `adaptive_pade_step!(state, f, n, h_init; adaptive_tol, ...)`
    `→ (state, P_u, meta)` — the full controller: iterate
    `T(h) > Tol ⇒ h := q·h` (FFW md:88-91) until acceptance or
    `max_rescales` exhausted, then run the standard
    `pade_step_with_pade!`.  Returns metadata `(h_used, h_step, T_h,
    n_rescales)` so callers can thread the controller's memory across
    steps.

The path-network driver implements the adaptive policy *inline* (not
via `adaptive_pade_step!`) because FFW md:93 specifies that under
rescale the wedge re-runs at the new `h_step` and the min-|u|
direction is re-selected.  `adaptive_pade_step!` is the simpler
contract appropriate for single-step callers (tests, direct API
consumers) and a future BVP/non-wedge integrator.

## Alternatives considered

**A.  Fixed `h = 0.5` (FW 2011 default; current shipped state).**
Rejected: degrades end-to-end on FFW 2017's PIII/PV figures and on
any path-network walk crossing a region of rapidly varying pole
density.  Worklog 017 already deferred adaptive `h` only because the
in-tree user list at the time (Phase 13 coord-transform tests) didn't
require it; bead `padetaylor-8ui` materialises the deferral now that
the FFW figure reproduction pipeline is the active work.

**B.  PI-controller on observed-residual.**  A classical
`OrdinaryDiffEq.jl`-style PI controller on the post-step residual of
the ODE (`u'' - f(z, u, u')` evaluated at `z+h`) would be more
general (no Taylor-coefficient hookup needed).  Rejected: it doesn't
use the structure of the Padé approximant — every step would recompute
a residual independent of the data the stepper *already has*.  FFW's
`T(h)` exploits the truncation series `q·w - p = Σ ε_k h^k` to a
quality-per-cost factor we should not abandon.  Worth revisiting if
this controller misfires on PVI (Tier-5 sheet-tracking).

**C.  Jorba-Zou order-h product (Jorba & Zou 2005 §3.2 eq. 3-8).**
The package's `StepControl` module already implements this; it was
*not* adopted as the path-network default because worklog 004
attempt-C showed it conflicts with FW's Padé-bridge paradigm (gives
overly conservative steps that miss the off-axis detour opportunity).
Rejected here for the same reason: the Jorba-Zou bound is a *Taylor*
truncation bound, not a *Padé* one, so it ignores the Padé
denominator's analytic-continuation effect.

**D.  Constant `q < 1` halving (FW-style step shrink).**  Replace
the FFW q-formula with `q = 0.5` (or any constant) and iterate.
Rejected: this is mutation M4 in the test file; loses Tol-monotonicity
(every Tol produces the same convergent `h`).  FFW's eq. 2 is the
minimum-viable controller; cheaper alternatives shed the headline
property.

**E.  Adaptive `h` exposed at `solve_pade`/the IVP layer too.**
Not in this ADR's scope; `path_network_solve` is the only consumer
the bead's downstream (FFW figures) requires.  `solve_pade`'s
existing real-axis stepping uses `StepControl.step_jorba_zou` for the
1-D regime; threading FFW adaptation through there is a separate
work item gated on a future user.

## Consequences

**Code**:
  - `src/PadeStepper.jl` gains the three helpers + an "Adaptive Padé
    step size" docstring section.  Code LOC remains under the 200-LOC
    cap (Rule 6).
  - `src/PathNetwork.jl`'s `path_network_solve` gains `adaptive_tol`,
    `k_conservative`, `max_rescales` kwargs and an inline adaptive
    branch.  The fixed-`h` default path is unchanged; the new path is
    isolated under `step_size_policy === :adaptive_ffw`.
  - `_solve_with_schwarz_reflection` threads the new kwargs through
    unchanged.

**Tests**:
  - `test/adaptive_step_test.jl` adds 7 testsets / 43 assertions
    (AS.1.1 - AS.1.7).  All mutation-proven (M1-M4); footer logs
    bite counts.
  - `test/pathnetwork_test.jl` PN.4.1 updated: the previous
    `@test_throws ArgumentError` for `:adaptive_ffw` is replaced
    with the same throw against `:bogus_policy` (the unknown-symbol
    fail-fast branch).
  - Aggregate count: 1708 → 1751 GREEN (+43).

**Figures unblocked**:
  - All seven FFW 2017 figures (Fig 1-4 PIII, Fig 5-7 PV) become
    achievable in principle.  In practice they additionally require:
    non-uniform Stage-1 nodes (FFW md:67-72; separate bead — Step A2
    of the 11-step plan), and `SheetTracker` mapping (Tier-5, already
    shipped via `src/SheetTracker.jl`).

**Performance**:
  - Per-step cost under `:adaptive_ffw` is *higher* than under `:fixed`
    by an additional `ffw_truncation_error` call (one extra
    `taylor_coefficients_2nd` at order `n+1`, one extra `robust_pade`)
    per rescale iteration.  Expected typical `n_rescales = 0-2` on
    smooth segments; up to ~5 near pole walls.  Net wall time is
    competitive or better than fixed because adaptive *takes fewer
    coarse steps* in smooth regions.

**Compatibility**:
  - Default `step_size_policy = :fixed` preserves byte-identical
    output on every existing test and figure.  Opt-in via kwarg.

## References

  - FFW 2017 §2.1.2 — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-97`.
  - FW 2011 §5.1 line 277 (fixed `h = 0.5` baseline) —
    `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:271-281`.
  - ADR-0004 — path-network architecture (consumed by this ADR).
  - `docs/worklog/017-coord-transforms-pIII-pV.md` §"What is NOT
    shipped" — the prior deferral note.
  - `docs/worklog/034-adaptive-pade-step.md` — implementation diary.
  - `src/PadeStepper.jl` — module docstring section "Adaptive Padé
    step size — FFW 2017 §2.1.2" + the three exported helpers.
  - `src/PathNetwork.jl` — docstring section "Adaptive Padé step size
    (`:adaptive_ffw`, opt-in)" + the inline rescale loop in
    `path_network_solve`.
  - `test/adaptive_step_test.jl` AS.1.1-AS.1.7 — acceptance tests.
