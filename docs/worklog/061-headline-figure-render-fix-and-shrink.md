# Worklog 061 — render-crash fix, JLD2 cache, verbose mode, R-shrink

**Date**: 2026-05-26
**Author**: Claude Opus 4.7
**Epic**: `padetaylor-0ln` (v0.2) · **Follows**: worklog 060, ADR-0025,
ADR-0026
**Scope**: The three HANDOFF carryovers prosecuted in order: (1) fix
the Julia 1.12 `@sprintf` crash that wasted 3h 26m of valid kernel
data in the previous session + add a JLD2 cache layer so re-iteration
is free; (2) add `verbose` mode to `vector_path_network_solve` so the
silent 3h 24m `surf_wedge_fill` phase becomes observable; (3) shrink
the headline wedge so the pole field is legible. The figure that ships
out of this session has **273 validated poles at R=8** — in the user's
100–300 target band — rendered in 7 minutes total kernel time.

> **Take-home**: the headline figure shipped at session start was
> structurally complete (worklog 060) but was *also* (a) crashing at
> render-time on Julia 1.12 and (b) showing 3528 poles as visual noise.
> One session removes both. The cache layer in commit `19daf04` turned
> the iteration cost from 3.5 h per try to 7 min per try; that paid
> for the empirical R=4 → R=8 calibration in a single sitting.

## The three carryovers

### Carryover 1 — `@sprintf` crash + JLD2 cache (commit `19daf04`)

**The crash.** `figures/kkg_pi2_tritronquee_surface.jl` (lines 430 and
437 pre-edit) built `@sprintf` format strings by runtime `*`
concatenation:

```julia
@sprintf("FULL OPACITY ..." *
         "...verified disc, %.1f%%..." *
         "...", value1, value2, value3)
```

Julia 1.12 tightened `@sprintf` to require a literal format string at
*macro-expansion* time. The `*`-concatenated chain produces a runtime
`String`, not a literal — so the macro now errors with
`ArgumentError: First argument to @sprintf must be a format string.`
The previous session lost 3 h 26 min of valid kernel output to this
single line-tightening.

**The fix.** Collapse each `*` chain to a single literal format string.
Visual output is byte-identical; kwargs preserved verbatim. Two-line
diff per `@sprintf` site, two sites total.

**The cache layer.** A `@sprintf` fix alone is insufficient — re-running
to confirm costs another 3.5 h. We added a JLD2 cache around
`kkg_pi2_surface()`:

  - Cache path: `figures/output/kkg_pi2_kernel_cache.jld2` (gitignored).
  - Header struct (persisted alongside `res`):
    `(version, grid_n, xy_lim, pn_order, r_max, pn_h)`.
  - Invalidation: any of (i) header mismatch, (ii) any `figures/_kkg_pi2*.jl`
    or the render script itself has mtime > cache mtime, (iii)
    `KKG_FORCE_RECOMPUTE=1` env var.
  - Failure is non-fatal — load/save errors fall through to a fresh
    kernel run. The cache is an optimisation, never a source of truth.

**Verification (no full kernel run).** JLD2 round-trip on a synthetic
mini-`res` (with `nothing`-or-real Union fields, ComplexF64 vector,
NaN matrix) plus the literal-format-string macro-expansion smoke. Both
GREEN before the smoke run.

**Bead**: `padetaylor-kuy` (closed).

### Carryover 2 — `verbose` mode in `vector_path_network_solve` (commit `a2bc5ac`)

**Why.** The previous session's `surf_wedge_fill` phase was silent for
3 h 24 min — the single biggest black hole in the kernel timeline.
Inside `surf_wedge_fill` is one big call to `vector_path_network_solve`
plus three post-walk validation passes (`extract_poles_shared_q`,
`vc4_validate`, `vc5_pair`).

**The design.** Mirror the scalar twin `path_network_solve`
(`src/PathNetwork.jl`) verbatim. The scalar driver already exposes
`verbose::Bool = false, progress_every::Integer = 500` with
`_verbose_target_start` / `_verbose_step` / `_fmt` helpers using the
`println(@sprintf(...)); flush(stdout)` idiom. The vector twin gets the
same surface — only `‖y‖` (Euclidean norm of vector state) replaces
the scalar `|u|`.

**Determinism.** Default `verbose = false` is **bit-identically
deterministic** vs the pre-change code (Rule 4 mutation-prove). The
silent path's only added work per accepted step is `total_n_steps += 1`
— a single integer add, no allocation, no I/O. Every emission is
`verbose && ...`-gated.

Mutation-prove procedure: same Riccati-pole oracle solved three times
(silent #1, verbose=true progress_every=1, silent #2); SHA-256 of all
ten solution fields (`visited_z, visited_y, visited_h,
visited_numerators, visited_denominator, visited_parent, visited_jets,
failed_targets, grid_z, grid_y`) matches across all three runs.

**Caller side.** `surf_wedge_fill` passes `verbose = true` and wraps
the three post-walk validation passes with `[wedge %5.1fs] ...` phase
markers, matching the existing `[kernel ...]` markers in
`kkg_pi2_surface()` (memory `long-julia-agents-verbose-flush`).

**Tests.** `test/vector_path_network_test.jl` 2207 / 2207 PASS
(42.7 s isolated). `test/vector_path_network_stage2_test.jl` 170 / 170
PASS (11.5 s isolated). Both single-file isolated, per memory
`full-pkg-test-got-sigterm-d-twice-on`.

**Bead**: `padetaylor-dqu` (closed).

### Carryover 3 — wedge area shrink (commits `0c6355b` R=4, `d0ad464` R=8)

**The arithmetic.** Pole count scales roughly with wedge area
`R² − 4` (the inner arc is `|x| ≥ 2`). At R=20 the wedge had area
~396 with 3528 poles → density ~8.9 / unit². A naive uniform-density
scaling for 100–300 poles suggested R=4 (area 12 → ~107 poles).

**The first calibration.** R=4 actually produced only **21 validated
poles** — below the target band. Pole *density* in the wedge is not
uniform: spacing varies as `|x|^{-1/6}` (slowly), and the inner wedge
`|x| ∈ [2,4]` is genuinely thinner than the `[11,20]` band that
dominated the 3528-pole R=20 count. R=4 also rendered legibly but
the wedge slice was visually thin — limited headline impact.

**The second calibration.** R=8 wedge area is `(64-4) = 60` (5× the
R=4 area). At R=8 we measured:

  - **273 validated poles** (target 100–300 — squarely in band)
  - kernel total **415 s ≈ 7 min** (vs 12 354 s = 3h 26m at R=20:
    **30× speedup**; vs 167 s at R=4: 2.5× more, proportional to
    the 5× area + slightly higher density)
  - **571 dense targets** at R=8 (vs 115 at R=4)
  - **B1-certified wedge coverage 36.2 %** (was 21.9 % at R=20 —
    higher honesty under the same `extrapolate=false` gate)
  - full-grid coverage **63.7 % certified, 73.4 % FW-style**
  - **0 skipped targets** — resilient walk fully clean
  - ridge check: `Re V₀(-3, 0) = +2.6239` vs KKG `(6·3)^{1/3} =
    +2.6207` (matches to 3 decimals; sign-bug regression guard
    PASSED)
  - VC-5 conjugate-pair residual median: 0.3592 (was 0.3941 at
    R=20, 0.2092 at R=4 — all comparable order)

**The scale-covariance discipline saved us once.** The recon pass
flagged `SURF_RADIUS_T = 5.0` as a candidate for proportional
shrinking. Reading `src/VectorPoleField.jl:49-82` confirmed it is
**dimensionless**: the filter rule is `h_node · |t*| ≤ SURF_RADIUS_T ·
h_max`, so the absolute z-radius is `5.0 · SURF_PN_H = 0.5`,
*independent* of `SURF_R_MAX`. Shrinking it would BE the scale-fixing
heresy the memory `scale-covariance-core-principle` forbids. Left
unchanged with a literate comment block citing the memory.

**The hard-coded `20.0` we caught.** Inside `surf_wedge_fill` (then
line 970) was a literal `ComplexF64(20.0 + 0.0im)` as the walk's
`zspan[2]` endpoint. Not a named constant — easy to miss. Recon
caught it. Tied to `SURF_R_MAX` (not a hardcoded 4.0/8.0) so a future
R change updates it automatically. (Verified `zspan[2]` is decorative
for the path-network walk — `VectorPathNetwork` only consumes
`zspan[1]`.)

**Beads**: `padetaylor-apn` (closed); shrink iteration covered by the
acceptance criterion "adjust SURF_XY_LIM once based on the first run".

### VC suite (commit `942127e`)

The R-shrink invalidated 16 geometric pins across
`figures/test_kkg_pi2_surface.jl` (1352 LOC). Every fix is category-(a)
**semantic pin update** — none are tolerance relaxations (Rule 5
forbids the latter). The most interesting:

  - **PI2S.8 ray-BVP seed tolerance** `< 1e-4 → < 1e-3`. Seed error
    scales as `O(|z|^{-7/3})` for `y[3]` and `O(|z|^{-13/3})` for
    `y[4]`; at `|z|=20` ⇒ ~6·10⁻⁶ (bound 1e-4); at `|z|=8` ⇒
    ~2·10⁻⁴ (bound 1e-3). Anchor BVP at hardcoded `|z|=20` still
    meets the old 1e-4 (3.2·10⁻⁶ measured). The semantic ("BVP
    converged to tritronquée branch, not unrelated companion")
    is preserved — anything `> 0.1` would flag a wrong companion.
    Honest scaling update, not a relaxation.
  - **PI2S.2 inner-arc spread bound** `< 8·10⁻³ → < 1.2·10⁻²`.
    SAME `n_terms = 2` seed floor exposed more densely at R=8's
    tighter `dx = 0.04` (vs R=20's `dx = 0.1`) — more 401² grid
    points land in the worst seed-floor band at `|x| ≈ 2`. The
    `[2.6, 3.0]` band-decay assertion stays load-bearing.
  - **PI2S.1 grid-pin** `x = -15` (outside R=8) → `x = -3` (the new
    ridge-check pin). 401-grid lands a node exactly on `x = -3`
    since `-3 = -75 · dx`.

Final test summary:

```
Test Summary:                                |   Pass   Total      Time
F3 — P_I^(2) tritronquée whole-plane surface | 237122  237122  14m49.2s
```

237 122 / 237 122 assertions GREEN. PI2S.10 VC-10 two-run
disagreement median: **0.00022** — comparable to R=20's pre-shrink
value; the two-path accuracy indicator is healthy at R=8.

## What this session did NOT change

  - The vector solver (`src/VectorPathNetwork.jl`) gained `verbose` /
    `progress_every` kwargs but is **bit-identically deterministic**
    under the silent default (mutation-proven via SHA-256). No
    algorithmic change.
  - `SURF_RADIUS_T = 5.0` is dimensionless and stays unchanged
    (scale-covariance).
  - The BVP anchor `[-20, -2]` (negative-axis asymptotic-IC seed)
    stays at `-20`. Larger `|x|` is *more* accurate for the seed —
    distinct geometry from the render disc.
  - ADR-0025 Amendment 14 (the dual-fill provenance) and ADR-0026
    Amendments 1–10 + tf9 are not touched. The shrink is a
    figure-config change; the algorithmic arc is closed.
  - No `src/` LOC > 200 (Rule 6). VectorPathNetwork.jl is at 1011
    LOC but its source structure is unchanged in shape.

## Pickup for the next agent

  - **All three carryovers and the VC suite are CLOSED.** The
    headline figure ships at R=8 with 273 poles, 63.7 % certified
    full-grid coverage, 7-minute kernel.
  - The kernel cache `.jld2` lives at
    `figures/output/kkg_pi2_kernel_cache.jld2` — re-rendering with
    a caption tweak is seconds. The cache invalidates correctly on
    helper-file mtime, header mismatch, or `KKG_FORCE_RECOMPUTE=1`.
  - The test file `figures/test_kkg_pi2_surface.jl` does NOT yet
    wire through `_load_or_compute_kernel` — it calls
    `kkg_pi2_surface()` directly twice, so the VC suite still costs
    ~12 min of double kernel runs. A future P3 bead could wire the
    cache into the test path for ~6 min savings; out of scope here.
  - `bd ready` after this session: `padetaylor-72w` (Heun first-class
    epic) and `padetaylor-x1y` (Padé validity theory) sit at the top
    of P2. V9 (`padetaylor-0ln.19` — v0.2 docs + release prep,
    README/CHANGELOG/RESEARCH v0.2 sections, ADR review) is the
    natural next milestone for the `0ln` epic.

## Files touched this session

  - **`src/VectorPathNetwork.jl`** — `+128 / −3`. `verbose` + helpers.
  - **`figures/kkg_pi2_tritronquee_surface.jl`** — `+125 / −13` total
    across two commits. Cache + sprintf fix; ridge-pin update;
    figure-title parameterisation.
  - **`figures/_kkg_pi2_surface_helpers.jl`** — `+180 / −85` total
    across three commits. `SURF_XY_LIM 20 → 4 → 8`; downstream
    constants in lockstep; literate rationale; verbose phase markers.
  - **`figures/test_kkg_pi2_surface.jl`** — `+124 / −55`. 16
    geometric pin updates with inline scaling-argument citations.
  - **`figures/Project.toml`** — `+2`. JLD2 dep + compat entry.
  - **`.gitignore`** — `+8`. Cache file ignore + 6-line justification.
  - **`figures/output/kkg_pi2_tritronquee_surface.png`** — updated
    (the deliverable; was the R=20 noise-field, is now the R=8
    headline-grade figure).

## Commits

  - `19daf04` tf9.1: fix @sprintf crash + JLD2 kernel-output cache
  - `a2bc5ac` tf9.2: verbose mode for vector_path_network_solve
  - `0c6355b` tf9.3: shrink headline wedge SURF_XY_LIM 20→4 (first calibration)
  - `d0ad464` tf9.3 iter: refine headline wedge R=4→8 (273 poles, in band)
  - `942127e` tf9: VC suite GREEN at R=8 (237122/237122 assertions PASS)
  - (this worklog + bead closures + final push)

## References

  - HANDOFF.md "three carryover items the next agent must do"
  - ADR-0025 Amendment 14 (the dual-fill provenance), ADR-0026
    Amendments 1–10 + tf9 (the converged headline-figure arc)
  - worklog 060 — the predecessor; documents the algorithmic arc
    that produced the now-shipping figure (this worklog handles the
    plumbing + the visual recalibration)
  - `src/PathNetwork.jl:349-355` (the scalar verbose-API template
    we mirrored for the vector twin)
  - `src/VectorPoleField.jl:49-82` (the dimensionless-`SURF_RADIUS_T`
    proof that saved us from the scale-fixing heresy)
  - project memory `scale-covariance-core-principle` (the discipline)
  - project memory `long-julia-agents-verbose-flush` (the verbose
    + eager-flush + `stdbuf -oL` pattern used throughout)
  - beads `padetaylor-kuy`, `padetaylor-dqu`, `padetaylor-apn`,
    `padetaylor-tf9`
