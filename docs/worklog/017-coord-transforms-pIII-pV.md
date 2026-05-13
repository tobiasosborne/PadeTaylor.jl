# Worklog 017 — CoordTransforms: PIII / PV exponential maps (bead `padetaylor-bvh`)

**Date**: 2026-05-13 (late evening, follow-up to worklog 016)
**Author**: Claude Opus
**Scope**: Phase 13 / Tier-4 v1 — ship the exponential coordinate-map
helpers for the third and fifth Painlevé equations per FFW 2017 §2.1
lines 39-48.  New module `src/CoordTransforms.jl` (~100 LOC) + 6
testsets / 23 assertions in `test/coord_transforms_test.jl`.  Three
mutations bit as predicted.  GREEN at 1262 → 1285.

> Take-home: the PIII and PV transforms (`z = e^(ζ/2)`, `u = e^(-ζ/2) w`
> for PIII; `z = e^ζ`, `u = w` for PV) push the fixed branch point at
> `z = 0` out to `ζ = -∞`, leaving meromorphic equations in `ζ`.  This
> module ships ONLY the transform mechanics — IC conversion + RHS
> factory closures.  Non-uniform Stage-1 nodes + adaptive Padé `h` are
> still TODO for the full FFW figures.

## What changed

  - `src/CoordTransforms.jl` (new, ~100 LOC) — exports
    `pIII_transformed_rhs`, `pV_transformed_rhs`,
    `pIII_z_to_ζ`, `pIII_ζ_to_z`, `pV_z_to_ζ`, `pV_ζ_to_z`.

    The two RHS factories take the four PIII/PV parameters `(α, β, γ,
    δ)` and return a closure `(ζ, w, wp) -> w''` ready to drop into
    `PadeTaylorProblem`.  The four IC helpers map a state at a single
    point between the `(z, u, u')` and `(ζ, w, w')` representations
    (forward + inverse).

  - `src/PadeTaylor.jl` — includes the new module, re-exports its 6
    names; module-header architecture list gains row 12.

  - `test/coord_transforms_test.jl` (new, ~200 LOC) — 6 testsets:
    - **CT.1.1** — PIII IC round-trip at the FFW Fig 1 IC `(z=1, u=1/4,
      u'=1)`.  Verifies `(ζ, w, w') = (0, 1/4, 5/8)` (hand-derived;
      see module docstring) AND the inverse round-trip.
    - **CT.1.2** — PIII RHS evaluates to FFW eq. at `(0, 1/4, 5/8)`
      with FFW Fig 1 parameters.  Hand-computed value `w'' = 111/256`
      pinned to 1e-14.
    - **CT.1.3** — PIII end-to-end agreement: one PathNetwork Padé
      step from `z=1` to `z=1.05`, computed both DIRECTLY and via
      `transform → solve in ζ → invert`.  `|u_direct - u_recovered| ≤
      1e-10`.  Validates the transform is consistent with the direct
      equation (the load-bearing test).
    - **CT.1.4** — PV IC round-trip at `(z=1, u=1/2, u'=3/10)`.
    - **CT.1.5** — PV RHS evaluates to FFW eq. at the same point.
      Hand-computed `w'' = -0.1675` pinned to 1e-14.
    - **CT.1.6** — PV end-to-end agreement at `z=1.05`, same shape as
      CT.1.3.

    PIII test parameters mirror FFW 2017 Fig 1 caption verbatim:
    `α = β = -1/2, γ = 1, δ = -1`, `(u(1), u'(1)) = (1/4, 1)`.  PV
    parameters chosen to keep `w ∉ {0, 1}` along the step neighborhood.

  - `docs/figure_catalogue.md §6` row T4 — marked PARTIAL with the
    ExpCoords transform shipped + the deferral notes for non-uniform
    nodes and adaptive Padé `h` (each is independent of the other and
    of this work).

## Mutation-proof procedure

Three load-bearing mutations applied + restored.

**Mutation L** — in `pIII_transformed_rhs`, swap the `/4` divisor to
`/3` on the parametrised term group.  **Bit predicted**: CT.1.2 (the
RHS pin) and CT.1.3 (end-to-end agreement).  **Verified bite**: CT.1.2
line 64 RED (hand-computed `111/256 ≈ 0.4336` vs perturbed value);
CT.1.3 lines 110-111 RED (transformed solve diverges from direct ≫
1e-10).  CT.1.1 stays GREEN (IC round-trip independent of the RHS).

**Mutation M** — in `pIII_z_to_ζ`, flip the sign of `z² u'` in `wp`.
**Bit predicted**: CT.1.1 (IC pin) and CT.1.3 (transformed solve
starts from wrong `w'`).  **Verified bite**: CT.1.1 lines 38 + 44 RED
(`wp ≈ 5/8` and `up_recovered ≈ up₀` both fail — `wp` reads `-3/8`,
inverse round-trip yields `up = -1` instead of `+1`); CT.1.3 lines
110-111 RED.

**Mutation N** — in `pV_transformed_rhs`, swap `(w+1)/(w-1)` to
`(w-1)/(w+1)` (the `δ e^(2ζ)`-term denominator/numerator swap).  **Bit
predicted**: CT.1.5 (PV RHS pin) and CT.1.6 (PV end-to-end agreement).
**Verified bite**: CT.1.5 line 158 RED (hand-computed `-0.1675` vs
perturbed value); CT.1.6 lines 194-195 RED.

All three mutations restored before commit; full suite GREEN at 1285.

## Design decisions

**Module is helpers-only, not a driver.**  The transform layer
encapsulates the COORDINATE LOGIC (RHS closures + IC conversion)
without owning the integration loop.  Downstream callers compose:
`prob = PadeTaylorProblem(pIII_transformed_rhs(α, β, γ, δ), (w₀, w'₀),
(ζ_start, ζ_end))` then call `path_network_solve(prob, ...)` directly.
This matches the existing API surface (cf. `LatticeDispatcher` and
`Dispatcher` are composition layers; `RobustPade` and `BVP` are
helpers).  No new `Solution` struct, no driver function with a long
kwarg list.

**RHS factory closures over four parameters.**  Avoids re-evaluating
the parameter algebra on every step.  The closure captures `α, β, γ,
δ` once at construction; `exp(ζ)` and `exp(2ζ)` are recomputed per
step (the algebraic cost is small vs the Taylor-coefficient generation
that dominates).

**`PadeTaylorProblem` accepts the transformed RHS verbatim.**
`PadeTaylorProblem`'s 2nd-order branch already expects `f(z, u, u')`
where `z` is complex.  Renaming `z → ζ`, `u → w`, `u' → w'` is purely
notational — no scaffolding needed in `Problems.jl`.

**End-to-end agreement is the load-bearing test, not the symbolic
RHS pin.**  The hand-computed RHS values (CT.1.2, CT.1.5) lock in the
formula transcription, but they're not the spec.  The spec is "the
transform should not change the analytic solution," which CT.1.3
and CT.1.6 verify by comparing the direct PIII / PV solve against the
transformed solve at a single step.  If either RHS formula is wrong
in a way that the symbolic pin misses (e.g., a sign flip the
mutation-test framework couldn't catch with the chosen sample point),
end-to-end agreement still catches it.

**Single PathNetwork-step at `z = 1.05`.**  Chose `h = 0.05` in `z`
(so `h ≈ 0.0976` in ζ) — small enough that both direct PIII (with
its `1/u` term) and transformed P̃_III (with its `1/w` term) stay far
from the singular surfaces.  PathNetwork's intrinsic truncation error
at order 30, h=0.5 is around `1e-12` for smooth segments; the 1e-10
tolerance is conservative.  Longer-range agreement could be tested
with a multi-step grid, but adds wall-time + Stage-2 nearest-visited
bookkeeping that complicates the test.

## What is NOT shipped (deferral notes)

**Non-uniform Stage-1 node placement (FFW2017...md:67-72).**  FFW
notes that PIII/PV pole densities grow exponentially with `Re ζ`,
making uniform Stage-1 grids wasteful at small `Re ζ` and
under-resolved at large `Re ζ`.  FFW prescribes a node-separation
function `R(ζ) = (8 - Re ζ) / 20` for their Fig 1.  The current
`path_network_solve` takes a flat grid; adding `R(ζ)` would require
either (a) a higher-level grid-generator that wraps `path_network_solve`,
or (b) a kwarg on `path_network_solve` itself.  Independent of the
transform — file a bead if a downstream user needs it.

**Adaptive Padé step size (FFW2017...md:74-97).**  Similar story;
independent of the transform; not on critical path for shipping the
transform helpers.  The bead description noted this could be filed
separately "if substantial" — I'm not filing it (no downstream caller
asking for it yet).

**Multi-sheet output recovery.**  The transforms encode ONE
canonical branch (principal `log`).  Users wanting sheet `s` would
shift `ζ` by `4π·im·s` (PIII) or `2π·im·s` (PV) before solving.  The
module docstring mentions this; no API support is needed beyond what
exists.

## Frictions surfaced

  - **PV needs `w ≠ 1`**, not just `w ≠ 0`.  P̃_V has both `1/(w-1)`
    and `(w+1)/(w-1)` terms.  Initial PV test IC `(u=1, u'=...)` would
    have hit the singular surface immediately.  Picked `u=1/2` (so
    `w=1/2 ≠ 1`) instead.  Documented in the module docstring.

  - **Hand-derivation of `wp_0` for the PIII IC** caught a subtle
    factor on first cut.  `dw/dζ = (1/2) e^(ζ/2) u + e^(ζ/2) u' · (dz/dζ)`
    expands to `(z u + z² u') / 2`, NOT `(z u + z u') / 2`.  The extra
    `z` comes from `dz/dζ = (1/2) e^(ζ/2) = z/2`.  CT.1.1 catches sign
    AND factor errors via the hand-pinned `5/8`.  Mutation M (sign
    flip) bit; a factor-flip (using `z` instead of `z²`) would also
    bite.

## Beads

  - `padetaylor-bvh` — closed in this session.
  - **Not filed**: non-uniform Stage-1 nodes, adaptive Padé `h`.
    Both are sketched in FFW 2017 §2.1.2; useful for high-`Re ζ`
    figures (FFW Figs 1, 4 ramping past `Re ζ ≈ 8` are the canonical
    cases).  No downstream caller in-tree is blocked.  File when
    actually needed.

## Pointers

  - `src/CoordTransforms.jl` — module top docstring documents the
    transforms + derivation + deferrals.
  - `test/coord_transforms_test.jl` — 6 testsets + mutation-proof
    procedure at end-of-file.
  - `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md`
    :29-37 (PIII/PV/PVI equations), :39-48 (the transforms +
    transformed equations), :67-72 (non-uniform nodes note),
    :100-104 (Fig 1 caption with parameters + IC), :43 (P̃_III eq.),
    :47 (P̃_V eq.).
  - `docs/figure_catalogue.md §5` row T4 — PARTIAL acceptance marked.

## Hard-won lessons (for HANDOFF.md §"Hard-won")

26. **Hand-derive the IC-conversion factors on paper before coding**.
    `dw/dζ` for the PIII transform involves a chain-rule cascade:
    `w(ζ) = e^(ζ/2) u(z(ζ))` ⇒ `w' = (1/2) e^(ζ/2) u + e^(ζ/2) (du/dz)
    (dz/dζ) = (z u + z² u') / 2`.  The `z²` (not `z`) on the `u'` term
    is the load-bearing factor — easy to slip when transcribing.  Pin
    it in a test with a hand-computed value, then mutation-test by
    flipping the sign or the factor.

27. **End-to-end agreement is a stronger spec than a symbolic-RHS
    pin**.  Symbolic pins can miss errors that compensate algebraically
    at one sample point.  Direct-vs-transformed agreement at a single
    Padé step catches anything that DOES change the analytic solution.
    Both shipped — symbolic pins catch transcription errors with crisp
    error messages; end-to-end catches deeper algebra mistakes.

28. **Helpers-only module beats premature driver abstraction.**
    Phase 13's scope could have been a full `pIII_pV_solve(α, β, γ, δ,
    z_grid; ...)` driver wrapping PathNetwork.  Choosing
    helpers-only (RHS factory + IC conversion) lets callers compose
    naturally with the existing `PadeTaylorProblem` + `path_network_solve`
    API.  Net `+100 LOC` source vs ~`+250 LOC` for a driver wrapper;
    same downstream usability.
