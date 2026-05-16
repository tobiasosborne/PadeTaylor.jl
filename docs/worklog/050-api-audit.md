# Worklog 050 — API audit pass: pre-v1.0 holistic review (bead `padetaylor-qur`)

**Date**: 2026-05-16
**Author**: Claude Sonnet 4.6 (docs-lockstep pass)
**Bead**: `padetaylor-qur` — closed by commit 9296731.
**Scope**: Holistic pre-v1.0 review of the ~39 exported symbols across
23 modules in `src/`.  Produces a normalisation recommendation list +
spawns 12 follow-up beads.  No source code modified.
New: `docs/api-review-2026-05-16.md` (743 LOC).

> **Take-home**: The public API surface is in better shape than it looks
> from the kwarg-name scattershot — every solver returns a `*Solution`
> struct, `!`-suffix discipline is clean, `order` threading is uniform,
> Symbol-valued kwargs are uniformly validated, and 74 % of throws carry
> `Suggestion:` text.  The three v1.0-blocking findings are precise and
> mechanical: two kwarg renames (`h` canonicalisation, `max_iter`
> normalisation) and one export-list omission (`pii_*` / `piv_entire`).
> Neither the design nor the implementation needs rethinking — just three
> targeted patches.

## What this is

A pre-v1.0 holistic read of every public-facing symbol in the package:
kwarg naming, return-type shape, error-message style, `!`-suffix
discipline, PascalCase/snake_case consistency, export-list completeness.
The output is an advisory document at `docs/api-review-2026-05-16.md`;
bead `padetaylor-qur` tracks the audit itself.  Actual normalisation
work lands in the 12 spawned follow-up beads.

The audit read 23 source files in full (all modules under `src/`) and
re-counted the exported symbol list from `src/PadeTaylor.jl:159-180`
(39 symbols, including `painleveplot` stub and `PadeTaylorAlg`).  It
did NOT run the test suite and did NOT modify any source file.

## Top findings (v1.0-blocking)

Three findings are **v1.0-blocking**; they are the canonical pre-v1.0
work items for the next ready agent.

### (a).1 — Step-size kwarg: three names for one concept

From `docs/api-review-2026-05-16.md:75-115`:

> `solve_pade(; h_max)` vs `path_network_solve(; h)` vs
> `lattice_dispatch_solve(; h_path)` — three names for the step length
> the Padé-Taylor stepper takes.  The naming asymmetry is load-bearing
> in `lattice_dispatch_solve`'s docstring
> (`src/LatticeDispatcher.jl:240-242`):
>
> > "`h_path` — IVP step size for the path-network walker. Forwarded as
> > `h_path` to `path_network_solve` (manual-fallback path) **and as
> > `h` to `edge_gated_pole_field_solve` (the kwarg name on edge-gated
> > is `h`, not `h_path`).**"
>
> Canonicalise to `h` everywhere.  `solve_pade`'s `h_max` is
> half-justified but misleading — `solve_pade` does no `min(h_max,
> h_adapt)` selection (`src/Problems.jl:202`); it simply uses `h_max`
> as the step length.

**Bead**: `padetaylor-xds`.  Safe migration path: `Base.@deprecate_binding`-style kwarg shims for one release.

### (a).2 — Max-iteration kwarg: six different names

From `docs/api-review-2026-05-16.md:119-145`:

> `bvp_solve(; maxiter)` vs `dispatch_solve(; bvp_maxiter)` vs
> `edge_gated_pole_field_solve(; max_iter)` vs `solve_pade(; max_steps)`
> vs `path_network_solve(; max_steps_per_target)` vs
> `path_network_solve(; max_rescales)`.  Not all the same concept, but
> `maxiter` vs `max_iter` is purely stylistic — the single un-snake-cased
> kwarg in the codebase.
>
> Standardise on `max_<noun>` with snake-cased noun.  Rename
> `bvp_solve`'s `maxiter` → `max_iter`; cascade through
> `dispatch_solve`'s `bvp_maxiter` → `bvp_max_iter`.

**Bead**: `padetaylor-0xn`.

### (e).3 — Export-list omission: `pii_rational`, `pii_airy`, `piv_entire`

From `docs/api-review-2026-05-16.md:519-538`:

> The `Painleve` module declares `pii_rational`, `pii_airy`, `piv_entire`
> in its own export list (`src/Painleve.jl:93`) but they are NOT
> re-exported from the top-level `PadeTaylor` module
> (`src/PadeTaylor.jl:175-178`).  Calling `using PadeTaylor; pii_rational(1)`
> would fail with `UndefVarError` despite the constructor existing in
> the package.  The closest analogues (`tritronquee`, `hastings_mcleod`)
> ARE re-exported on `src/PadeTaylor.jl:177`.  Single-line fix.

**Bead**: `padetaylor-gvz`.  One-line addition to `src/PadeTaylor.jl`.

## Positive findings

The audit validated these as already consistent (from
`docs/api-review-2026-05-16.md:685-719` "Things that are already consistent"):

1. **All solvers return `*Solution` structs** — no raw tuple/array
   returns from any public driver.  `eval_at` / `eval_at_sheet` return
   `Tuple{Complex, Complex}` (correct for per-point queries).
2. **`!`-suffix discipline is clean** — the three `!` functions
   (`pade_step!`, `pade_step_with_pade!`, `adaptive_pade_step!`) all
   mutate their first arg (`PadeStepperState`) and are not exported.
   No non-`!` function silently mutates exported state.
3. **PascalCase / snake_case perfect** for all 39 exported symbols —
   no camelCase outliers; `BVPSolution` PascalCase acronym is
   Julia-idiomatic (matches `LAPACK`, `FFT`, `SVD` precedents).
4. **`order` kwarg threading uniform** — every solver defaults to
   `prob.order`; `PadeTaylorProblem`'s default is `30` (FW 2011
   canonical); closed-form family constructors and named transcendents
   all carry `order = 30` through consistently.

## What is NOT in this commit

No source files were modified.  The audit document is advisory: it
catalogues inconsistencies and positive findings, assigns bead IDs,
but makes zero changes to `src/`, `test/`, or `ext/`.

No source changes → no tests were run for this commit.  The existing
2601 expected assertions (2508 base + 32 DG.* + 61 LD.X.*; all verified
in prior sessions) are unchanged.

The 12 spawned beads carry the actual normalisation work:

| Bead | Priority | Scope |
|------|----------|-------|
| `padetaylor-xds` | v1.0-blocking (P2) | Canonical `h` rename |
| `padetaylor-0xn` | v1.0-blocking (P2) | `max_iter` normalisation |
| `padetaylor-gvz` | v1.0-blocking (P2) | Export `pii_*` / `piv_entire` |
| `padetaylor-e0k` | P3 cosmetic | Bare `error()` replacement (4 sites) |
| `padetaylor-cr6` | P3 cosmetic | `Suggestion:` lines for terse throws |
| `padetaylor-1nz` | P3 cosmetic | `BVPProblem` constructor |
| `padetaylor-bbl` | P3 cosmetic | `pole_field_mask` `level` → `edge_level` rename |
| `padetaylor-s1q` | P3 cosmetic | `quality_diagnose` sheet API |
| `padetaylor-brf` | P4 defer | Symbol kwargs glossary |
| `padetaylor-289` | P4 defer | Painlevé constructor taxonomy in README |
| `padetaylor-kdx` | P4 defer | `grid_u` shape collision documentation |
| `padetaylor-egc` | P4 defer | `verbose` routing unification |

## Bead

  - `padetaylor-qur` — closed by commit 9296731.

## References

  - **Audit document** — `docs/api-review-2026-05-16.md` (743 LOC):
    §1 executive summary; §3 kwarg naming (findings (a).1–(a).6);
    §4 return-type consistency; §5 error-message style; §6 bang-suffix;
    §7 capitalisation + export-list omission; §8 additional checks;
    §9 proposed follow-up beads; §10 positive findings; §11 limitations.
  - **`src/PadeTaylor.jl:159-180`** — the export block audited for
    re-export completeness.
  - **`src/Painleve.jl:93`** — `pii_rational`, `pii_airy`, `piv_entire`
    present in `Painleve` module exports but absent from top-level re-export.
  - **`src/LatticeDispatcher.jl:240-242`** — the docstring that apologises
    for the `h_path` / `h` split, cited as the canonical "API smell" marker.
  - **CLAUDE.md Rules 1, 2, 10**.
