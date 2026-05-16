# Worklog 041 — η-plane PVI RHS + IC helpers (bead `padetaylor-riu` closed)

**Date**: 2026-05-16
**Author**: Claude Opus 4.7 (1M context)
**Bead**: `padetaylor-riu` — closed by this commit.
**Scope**: Step A3 of the 11-step FFW 2017 figure-reproduction plan
— the **η-plane PVI RHS** (FFW 2017 eq. 5, md:154) and the
`z ↔ η` IC helpers, plus a `:transformed_eta` frame for
`PainleveProblem(:VI; ...)`.  The substrate piece that unblocks
FFW Fig 2 column 1 (η-plane render).

> **Take-home**: A3 was carried by yesterday's long FFW arc session
> as far as the mutation-proof step before being stashed.  This
> session re-pops the stash, reverts the deliberate **M2 sign-flip
> marker** at `src/SheetTracker.jl:248`, re-verifies all three
> mutations (M1, M2, M3) bite their target assertions, and ships.
> Test count **2220 → 2262 GREEN** (+42 assertions, all in the
> previously-written `test/eta_pvi_test.jl`).  No source changes
> beyond restoring the marker; no new design decisions.  ADR-less
> additive extension to `SheetTracker.jl`.

## What this session did

The substantive work — `pVI_eta_transformed_rhs`,
`pVI_z_to_η`/`pVI_η_to_z`, the `:transformed_eta` frame branch in
`Painleve._build_VI`, the test file with 6 testsets and a
mutation-proof footer — was completed yesterday in the FFW 2017
figure-arc session (worklogs 034–040) and stashed mid-mutation-test
per HANDOFF.md (`git stash list` showed
`stash@{0}: On main: A3 η-plane PVI WIP (stopped mid-mutation-test;
has M2 marker)`).

This session's three steps:

  1. `git stash pop` — surfaced the WIP onto `main`.

  2. Reverted M2 at `src/SheetTracker.jl:248`:
     `vp = -z * ζ * up  # M2 mutation: sign flip` → `vp = z * ζ * up`.
     (Per the function's docstring chain rule `dv/dη = (dw/dζ)(dζ/dη)
     = (z·u')·ζ = z·log(z)·u'`, the unflipped form is the correct
     one.  Inverse `pVI_η_to_z` divides by `z·ζ`; with the unflipped
     forward, round-trip is exact.)

  3. **Re-verified all three mutations** per CLAUDE.md Rule 4 (the
     yesterday-session footer's bite-count claims were within ε but
     I re-measured the exact bite counts directly):

       - **M1** (drop `- 1` in η-RHS second term): **4 RED of 42**
         — ET.1.1 (1) + ET.1.1bis (1) on the hand-pinned closed-form
         RHS values; ET.1.3 (2) on the end-to-end direct-vs-η-
         transformed agreement.
       - **M2** (sign flip on `pVI_z_to_η`'s `vp`): **6 RED of 42**
         — ET.1.2 (3) on direct chain-rule check + both round-trip
         assertions; ET.1.3 (2) on end-to-end; ET.1.5 (1) on the
         PainleveProblem round-trip through `to_frame`/`from_frame`.
       - **M3** (dispatch `:transformed_eta` to ζ-plane factories):
         **2 RED of 42** — ET.1.5's `zspan[1] ≈ η_expected` and
         `y0[2] ≈ vp_expected` assertions, since `to_frame` now
         returns `log z, z·up` not `log log z, z·log(z)·up`.

     Updated `test/eta_pvi_test.jl`'s mutation-proof footer with the
     exact bite counts in place of yesterday's approximate description.

## What shipped (carried from the previous session's stash)

  - **`src/SheetTracker.jl`** (+~140 LOC): `pVI_eta_transformed_rhs`
    factory (FFW eq. 5 written character-by-character from md:154,
    with the chain-rule `- 1` in the second-term bracket — the only
    place where the η-plane equation is NOT a structural copy of
    the ζ-plane equation with `e^ζ → E = exp(exp(η))`), plus the
    coordinate helpers `pVI_z_to_η` / `pVI_η_to_z` (composition of
    `ζ = log z`, `η = log ζ`).  Top-of-file docstring chapter
    section added documenting the η-plane double-exp transform and
    its branch-point-free region `Re η < log(2π) ≈ 1.83788`
    (rectangular in the η-plane vs ζ-plane's left-half-infinity).

  - **`src/Painleve.jl`** (+~30 LOC): `_build_VI` now dispatches on
    a new `frame::Symbol` kwarg with values `:transformed` (default,
    ζ-plane, unchanged) or `:transformed_eta` (η-plane, new).
    `_build_transformed` extended with two optional kwargs
    `frame_tag` and `extra_optionals` so it can stamp the right
    frame string and validate the new optional `:frame` kwarg.  The
    `solve_pade(::PainleveProblem)` complex-zspan guard relaxed
    from `pp.frame === :transformed` to `pp.frame !== :direct` so
    `:transformed_eta` is caught with the same fail-loud message.

  - **`src/PainleveSolution.jl`** (1-line docstring update): the
    `frame` field's enumeration grows the third value
    `:transformed_eta`.

  - **`src/PadeTaylor.jl`**: three new exports
    (`pVI_eta_transformed_rhs`, `pVI_z_to_η`, `pVI_η_to_z`).

  - **`test/eta_pvi_test.jl`** (new, ~325 LOC): six testsets
    ET.1.1 / ET.1.1bis (hand-pinned RHS at α=β=γ=δ=0 and ={1,1,1,1}
    using the closed-form algebra inline in comments — `(2/3)(1 -
    log 2) + (125/24)(log 2)²` at the canonical sample point); ET.1.2
    (z↔η round-trip at FFW Fig 2's IC `z=10, u=0.4295…,
    u'=-1.617e-3`); ET.1.3 (end-to-end PVI direct-vs-η-transformed
    agreement at one Padé step from `z₀=5` to `z=5.05`); ET.1.4
    (branch-point-free region check on samples inside / outside / on
    the k=±1 lattice branch points); ET.1.5
    (`:transformed_eta` frame wiring through `PainleveProblem(:VI)`
    including unknown-frame fail-loud).  Wired into
    `test/runtests.jl`.

## Frictions

### 1.  The yesterday-session's footer counts were slightly low.

Yesterday's mutation-proof footer (written before stashing) claimed
"three independent assertions" for M1 and "ET.1.5 line 245 RED" for
M3.  The actual bites are M1 = 4 RED (not 3), M2 = 6 RED across
THREE testsets not two, M3 = 2 RED not 1.  Plausible cause: the
footer was drafted while writing the test file, before the mutation
loops were actually run — a forecast not a measurement.  Today's
measurement and footer-rewrite makes the doc claim match the
empirical reality.  **Lesson:** mutation footers should be filled
in *after* the mutation loop has actually been run; placeholder
text from the planning phase decays into approximation.

### 2.  `solve_pade`'s complex-zspan guard string was load-bearing.

`solve_pade(::PainleveProblem)` originally checked `pp.frame ===
:transformed` to gate its complex-zspan guard.  Without updating
the check to also cover `:transformed_eta`, an η-frame problem with
complex zspan (the typical case) would silently skip the guard and
hit a less-helpful downstream "state.z < z_end is undefined for
complex z" deep in the stepper.  Yesterday's session caught this
and relaxed the check to `pp.frame !== :direct`.  Confirmed today
via grep that this is the only `:transformed`-specific gate in the
solve-side code; the `painleveplot` Makie recipe and the
`figure_catalogue` documentation entries are frame-agnostic.

### 3.  No ADR for A3.

A3 is a pure additive extension to `SheetTracker.jl` + an opt-in
dispatch in `_build_VI`.  No new design decisions: the η-plane
equation is FFW 2017 eq. 5 verbatim, the `z → ζ → η` composition is
mathematically forced, and the `:transformed_eta` frame tag is a
parallel to the existing `:transformed` rather than a new pattern.
The HANDOFF and figure-catalogue describe A3 inline; no ADR file is
warranted.  (Contrast A1, A2, A6 which each got an ADR because each
made a genuine choice — controller policy, kwarg vs config struct,
hybrid-driver-vs-deepen-component, respectively.)

## Hard-won lesson

**The footer-as-spec antipattern.** I noticed yesterday's footer
described mutation bites in slightly idealised form
("three independent assertions catch the algebra error" — when
actually four did, of two qualitatively different shapes).  This
is the same antipattern as ADR-as-aspiration documented in
worklog 009: writing the doc *as if* the result, then never going
back to make it match what actually happened.  Concretely, the M1
bite is 4 RED not 3 because the closed-form pin at α=β=γ=δ=1 ALSO
bites (not just at α=β=γ=δ=0) — both share the `+ 1` constant
offset in the algebraic reduction.  The M2 bite is 6 not 3 because
the round-trip through `PainleveProblem.to_frame`/`from_frame` (a
different code path) ALSO catches it, in ET.1.5 — independently
load-bearing in the API surface, not just the helper.  The M3 bite
is 2 not 1 because the η-frame `y0[2]` (vp) is a *different* value
under ζ-plane factories, in addition to `zspan[1]` (η) — the
chain-rule factor `z·log(z)` vs `z` produces a different vp
component too.

**The right discipline: mutation loops first, footer-text last.**
Run the mutation, capture the exact pass/fail count, then write
the footer describing what bit.  This catches "wider-than-expected
bite" (load-bearing reach you didn't intend) and "narrower-than-
expected" (a brittle test the mutation slips past).  Both are
worth knowing.

## Empirical results

  | testset    | assertions | wall (s) | mutations biting |
  |------------|-----------|----------|------------------|
  | ET.1.1     |  2        | 0.1      | M1               |
  | ET.1.1bis  |  2        | 0.0      | M1               |
  | ET.1.2     | 10        | 5.0      | M2               |
  | ET.1.3     |  3        | 4.0      | M1, M2           |
  | ET.1.4     | 15        | 0.0      | (region pin; no mutation target) |
  | ET.1.5     | 10        | 0.7      | M2, M3           |
  | **Σ**      | **42**    | **~10**  | M1, M2, M3 all bite |

Full test suite: **2220 → 2262 GREEN** (+42 assertions, 5m 32s
wall on a single thread — +14 s vs yesterday's baseline).

## Beads

  - `padetaylor-riu` — closed by this commit.

## What is NOT shipped

  - **FFW Fig 2 column 1 rendering**.  That's a separate downstream
    step (B5).  This bead's job is the η-plane substrate; the
    figure is the next deliverable on top of it.  Per HANDOFF
    sequencing, B5 is gated on A4 (constrained-wedge routing) +
    A5 (sheet-aware Stage-2) too — A3 alone unblocks the
    *finite-rectangle* render but not the multi-sheet circumambulation
    in column 2/3.

  - **PVI sheet-counter primitive in η-plane**.  Yesterday-shipped
    `winding_delta` / `accumulate_winding` / `sheet_index` work
    against any branch point — they accept the branch-point coordinate
    as input, so they can be re-used against the η-plane lattice
    points `η = log|2π·k| + i·arg(2π·i·k)` without modification.  No
    code change needed; the friction is per-call, not API.

  - **`path_network_solve` against η-frame problems**.  The η-plane
    equation has `O(exp(exp(η)))` magnitudes at moderate `Re η`,
    which is mostly fine for the Padé arithmetic at `Re η < log(2π)`
    but stresses Float64 once outside.  No tests pin this; the
    figure-render bead (B5 / -bho) will exercise it and surface any
    issues.

## References

  - **FFW 2017 §2.2.1** —
    `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:146-162`
    (η-plane PVI: md:151 the `z = exp(exp(η))` map, md:154 eq. 5
    the equation, md:157 the branch-point-free region, md:148 eq. 4
    the lattice).  md:195 FFW Fig 2 caption (the IC point reused as
    ET.1.2's pin).
  - **Worklog 018** (`docs/worklog/018-sheet-tracker-pVI.md`)
    §"What is NOT shipped" — the original deferral that this bead
    closes.
  - **`src/SheetTracker.jl`** — the η-plane addition.
  - **`src/Painleve.jl` `_build_VI`** — the `:transformed_eta` frame
    branch.
  - **`test/eta_pvi_test.jl` ET.1.1–ET.1.5** — acceptance tests +
    mutation-proof footer.
  - **HANDOFF.md** §"Session 2026-05-15" §"What is NOT shipped from
    the 11-step plan" — the A3 stash description and resume
    instructions this session followed.
  - **Worklogs 034–040** — the FFW 2017 figure arc that A3 belongs
    to.
