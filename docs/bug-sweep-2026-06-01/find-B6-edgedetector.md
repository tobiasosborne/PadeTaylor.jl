# Bug sweep B6 — EdgeDetector + EdgeGatedSolve

## Area

Tier-2.5 pole-field edge detection and region-growing gate:

- `src/EdgeDetector.jl` — the 5-point discrete-Laplacian classifier
  (`laplacian_residual`, `pole_field_mask`, `_auto_level`) per
  FW 2011 §3.2.2 eq. (3.3).
- `src/EdgeGatedSolve.jl` — region growing (`_dilate`, `_erode`,
  `_open`, `_flood_fill`, `_solve_targets`,
  `edge_gated_pole_field_solve`) that confines `path_network_solve`
  to the pole field per FW 2011 md:401.

Special focus per the assignment: (a) the 5-point Laplacian stencil
coefficients + sign + `/h²`; (b) the h-aware threshold level
(ADR-0009); (c) the region-growing flood-fill neighbour connectivity
and boundary handling.

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:204-208`
  — §3.2.2, the 5-point stencil eq. (3.3) `Δu ≈ [stencil]·u/h²`, and the
  statement that the contour level `0.001` is on `log₁₀|Δu|`.
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:198`
  — Fig. 3.3 caption ("level curve 0.001 (solid line)").
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:304-310`
  — the random-path-selection accuracy-indicator passage (source of the
  documented run-to-run variation).
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:401`
  — "smooth regions are unstable regions for any IVP solver … force the
  path selection algorithm to complete pole fields before stepping into
  smooth regions" (the EdgeGatedSolve rationale).
- `docs/adr/0009-edge-detector-h-aware-level.md:42,55-60` — the
  h-aware level formula `LEVEL0 + 2·log₁₀(min(h, H0)/H0)`, anchor
  `(0.25, 0.001)`, clamp at `H0`.
- `docs/adr/0017-lattice-dispatcher-strict-mode.md:44-57` — the
  edge-gated IVP source / re-scatter contract.
- `docs/worklog/032-edge-detector-h-scaling.md:11-22,84-114,195-218`
  — derivation + four mutation-proofs of the h-scaling fix.
- `src/PathNetwork.jl:375-398,479-481` — `path_network_solve`
  signature, `rng_seed` default `0`, `shuffle(rng, …)` target order.
- `src/LatticeDispatcher.jl:323-349,447-462` — re-scatter index
  alignment and `_dilate_one`.
- `src/Problems.jl:108-132` — `PadeTaylorProblem.zspan` is the IC span;
  `zspan[1]` is the IC used as the seed centre.
- `test/edge_detector_test.jl:48-122,213-239` — stencil exact-zero,
  pole-residual, and `_auto_level` invariant assertions.

## Findings

No transcription bug (sign / off-by-one / normalisation / branch-cut /
aliasing) was found in either file. The stencil, the h-aware level,
and the morphological / flood-fill primitives all match ground truth.
The items below are LOW-confidence observations recorded for
completeness; none is a demonstrable mismatch against a cited reference
line.

### [LOW] Docstring coordinate convention is transposed vs the caller (cosmetic, not load-bearing)

- **Location**: `src/EdgeDetector.jl:166-170` vs
  `src/EdgeGatedSolve.jl:225,311` and `src/LatticeDispatcher.jl:304`.
- **Ground truth**: FW2011_*.md:204-206 — the eq. (3.3) cross stencil
  `[1; 1 -4 1; 1]/h²` is symmetric in the two lattice axes.
- **Code behaviour**: `laplacian_residual`'s docstring says
  "`u_grid[i, j]` is at `z = x_j + i·y_i` … rows index `y` and columns
  index `x`". Every actual caller (`EdgeGatedSolve._solve_targets`
  line 225, the seed loop line 311, `LatticeDispatcher` line 304) builds
  `grid[i,j] = xs[i] + im·ys[j]` — rows index `x`, columns index `y` —
  the transpose of the documented convention. Because the stencil
  (`src/EdgeDetector.jl:188-192`) is `u[i-1,j]+u[i+1,j]+u[i,j-1]+u[i,j+1]
  -4u[i,j]`, fully symmetric in `i↔j`, the residual value at every cell
  is identical under transposition. So the discrepancy changes nothing
  numerically.
- **Mechanism**: none — the docstring itself (line 170) states "the
  stencil is rotationally symmetric, so the convention is for the
  caller's convenience and does not affect the residual values." No
  intermittent (or any) discontinuity arises.
- **Intermittent?**: No.
- **Confidence**: 0.95 that this is harmless (i.e. NOT a bug). It is
  reported only because a future anisotropic stencil edit would silently
  inherit a wrong axis assignment.

### [LOW] Region-growth coverage is RNG/shuffle-dependent across passes (documented FW behaviour, not a transcription bug)

- **Location**: `src/EdgeGatedSolve.jl:326,348` (the two
  `_solve_targets` → `path_network_solve` calls) and
  `src/PathNetwork.jl:395,479-481`.
- **Ground truth**: FW2011_*.md:304-310 — "because of the random nature
  of our path selection algorithm … the difference between the two
  solutions can serve as a practical accuracy indicator." The
  run-to-run variation is an explicit, intended property of the FW
  algorithm.
- **Code behaviour**: `_solve_targets` calls `path_network_solve` with
  no `rng_seed`, so it uses the default seed `0`
  (`src/PathNetwork.jl:395`). The walk is therefore deterministic for a
  fixed target *list*, but the target list changes every region-growth
  pass (`targets = _dilate(field, grow_rings)` grows each iteration) and
  the final solve uses a different list again
  (`final_targets = _dilate(field, 1)`). `shuffle(rng, grid)`
  (`src/PathNetwork.jl:481`) reorders that list, so the visited-tree
  topology — and thus which Stage-2 cells fall within `visited_h` of a
  visited node vs return `NaN` (`extrapolate=false` default) — differs
  between passes. NaN cells are masked `false`
  (`src/EdgeDetector.jl:254`), so the field can grow slightly further on
  one pass than another.
- **Mechanism**: This affects how *far* the field grows, i.e. coverage,
  not the *values* of admitted cells. The mask is fail-safe in the only
  direction that matters: a NaN residual never produces a spurious
  `true` (the `isnan(real(z)) && continue` guard at
  `src/EdgeDetector.jl:254` precedes the threshold test), so a
  smooth-sector cell cannot leak into the field. A genuine
  discontinuity in the returned solution would require the gate to admit
  a smooth cell on some runs; the fail-safe direction prevents that.
- **Intermittent?**: Yes (shuffle/order-dependent), but the visible
  effect is field-extent jitter, not a value discontinuity — and it is
  the FW-documented design property, not a port error.
- **Confidence**: 0.8 that this is intended behaviour rather than a
  bug. If the maintainer's reported discontinuity is in EdgeGatedSolve
  output, the more likely culprit is *upstream* in `path_network_solve`
  / the Padé / step machinery (out of this area), not the gate logic.

### [LOW] `mask .|= field` permanently anchors the seed disc as pole-field regardless of its true Laplacian

- **Location**: `src/EdgeGatedSolve.jl:308-314,334,336`.
- **Ground truth**: FW2011_*.md:401 — the gate's job is to "complete
  pole fields before stepping into smooth regions"; the module docstring
  step 1 (`src/EdgeGatedSolve.jl:43-45`) asserts "the pole field
  originates at the IC."
- **Code behaviour**: the seed `field` is a Euclidean disc of radius
  `seed_r` around `z_ic` (line 311, inclusive `≤`). Every pass forces
  `mask .|= field` (lines 334, 336) before the flood-fill, so seed cells
  are treated as pole-field even if their discrete Laplacian is below
  threshold. This is by design (comment line 328-329: "the field is
  pole-field by definition").
- **Mechanism**: If the IC were placed such that the seed disc straddled
  a smooth region (e.g. a non-Painlevé problem, or an IC not at the
  pole-field origin), the forced-true seed could anchor the flood-fill
  into smooth terrain, producing inconsistent admitted regions. For the
  shipped Painlevé use case (IC at the pole-field origin per FW §3.1)
  this premise holds and the behaviour is correct.
- **Intermittent?**: No — deterministic given the IC and grid.
- **Confidence**: 0.6 that this is a deliberate, documented design
  choice rather than a latent bug. Recorded as a caller-precondition
  caveat, not a transcription error.

## Areas verified correct

- **5-point Laplacian stencil — coefficients, sign, and `/h²`**
  (`src/EdgeDetector.jl:188-192`). The code computes
  `inv_h2 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])`.
  This matches FW 2011 eq. (3.3) verbatim
  (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:206`:
  cross stencil `[1; 1 -4 1; 1]·u/h²`): four `+1` neighbours, a `-4`
  centre, divided by `h²`. Signs and the `1/h²` normalisation are
  correct. Cross-checked by `test/edge_detector_test.jl:48-89`
  (harmonic `z²` → residual ~roundoff) and `:91-122` (`1/(z-z₀)` →
  residual ≫ 100 near the pole). No sign flip, no missing/extra
  neighbour, no off-by-one in the `2:(nrow-1) × 2:(ncol-1)` interior
  loop (line 187).

- **h-aware threshold level `_auto_level`** (`src/EdgeDetector.jl:156-157`).
  Computes `T(_LEVEL0) + 2*log10(min(T(h), T(_H0))/T(_H0))` with
  `_H0 = 0.25`, `_LEVEL0 = 0.001` (lines 145-146). This matches
  ADR-0009:55-60 (`level(h) = LEVEL0 + 2·log₁₀(min(h, H0)/H0)`) and the
  worklog 032 mutation-proofs exactly. Sign is `+2·log₁₀` (M2 in
  worklog 032:199 proves a sign flip goes RED); the `min(h, H0)` clamp
  is present (M3 worklog 032:204 proves dropping it cascades to RED on
  EG.1.x). Confirmed by `test/edge_detector_test.jl:225-239`.

- **Threshold comparison is on `log₁₀|Δu|`, not bare `|Δu|`**
  (`src/EdgeDetector.jl:255`: `log10(abs(z)) > level_T`). Matches
  FW2011_*.md:208 ("Fig. 3.3 shows `log₁₀|Δu|` … select a contour level
  (here 0.001)") and the module docstring's explicit correction of the
  bead description (lines 75-82). The `Δu` already contains the `/h²`,
  so the level is applied to the normalised residual as FW intends.

- **NaN / boundary handling is fail-safe** (`src/EdgeDetector.jl:186`
  fills the result with `NaN+NaN·im`; only interior cells overwritten;
  `_mask_from_residual` line 254 `isnan(real(z)) && continue` before the
  threshold). A NaN residual — whether a boundary cell or an
  upstream-NaN-poisoned interior cell — yields `false`, never a spurious
  `true`. This is the safe direction: a NaN can only under-classify
  (cell stays smooth-side), never inject a phantom pole-field cell.
  Confirmed by `test/edge_detector_test.jl:139-143` (boundary cells
  false).

- **Fail-fast guards** (`src/EdgeDetector.jl:177-182`): rejects grids
  `< 3×3` and `h ≤ 0` with `ArgumentError` + detail text, per CLAUDE.md
  Rule 1. The 1-arg `pole_field_mask(Δu; level=:auto)` correctly throws
  because `:auto` is unresolvable without `h` in scope (lines 227-231),
  matching ADR-0009:85-87.

- **`issorted(xs; lt = ≤)` strict-increasing check**
  (`src/EdgeGatedSolve.jl:288`, `src/LatticeDispatcher.jl:276-279`).
  `issorted(v; lt)` tests `!lt(v[i+1], v[i])` on consecutive pairs; with
  `lt = ≤` this is `!(v[i+1] ≤ v[i])` ⟺ `v[i+1] > v[i]` ⟺ strictly
  increasing. The error message ("strictly increasing") matches the
  actual predicate — no off-by-one or wrong-direction comparison.

- **Uniform-spacing isotropy guard** (`src/EdgeGatedSolve.jl:298-303`):
  requires `Δx ≈ Δy` within `rtol=1e-10` before using a single `h_grid`,
  correct because the eq. (3.3) stencil is isotropic. `h_grid = Δx`
  (line 303) is the lattice spacing, which is the `h` the stencil needs
  — distinct from the IVP step `h` (line 326 passes the IVP `h` to
  `path_network_solve`, line 333 passes `h_grid` to `pole_field_mask`).
  The two are correctly kept separate.

- **`_dilate` / `_erode` — no in-place aliasing** (`src/EdgeGatedSolve.jl:140-181`).
  Each of the `r` single-ring passes does `prev = copy(out)` (lines 144,
  167) and reads `prev` while writing `out`, so a cell flipped earlier
  in a pass cannot influence a later cell in the same pass. `r`
  sequential passes give a true Chebyshev-distance-`r` (8-connected)
  morphological op. `_erode`'s grid-edge handling (line 173: a window
  hanging off the lattice erodes the cell away) is internally
  consistent with its docstring. No RNG, no order dependence in these
  primitives.

- **`_flood_fill` — correct 8-connected seeded fill** (`src/EdgeGatedSolve.jl:191-213`).
  Seeds `reached` from `seed ∩ mask`, marks `reached[ii,jj]=true` at
  push time (line 208, before pushing), so no cell is enqueued twice and
  the result is order-independent. Bounds check `1 ≤ ii ≤ nx && 1 ≤ jj ≤
  ny` (line 206) is correct. A stranded false-positive component not
  connected to the seed is correctly excluded. `_open = _dilate(_erode(…))`
  (line 186) is standard morphological opening.

- **Growth-loop ring arithmetic** (`src/EdgeGatedSolve.jl:290-293,325,339`).
  The `grow_rings ≥ 2` guard is correct: a pass solves `targets =
  dilate(field, grow_rings)`, but only cells whose full stencil
  neighbourhood lies inside `targets` get a non-NaN residual, i.e. the
  classifiable set is `targets` eroded by 1 = `dilate(field,
  grow_rings-1)`. With `grow_rings = 1` no new cell beyond `field` is
  ever classifiable; the guard's detail text (line 291-293) states this
  correctly. `new_cells = reachable .& targets .& .!field` (line 339)
  correctly clips any opening-spill back into `targets` and removes
  already-admitted cells.

- **Seed centre = IC** (`src/EdgeGatedSolve.jl:306` `z_ic =
  Complex{T}(prob.zspan[1])`). `prob.zspan[1]` is `z_start`, the IVP
  initial point (`src/Problems.jl:120-132`), which is where the Painlevé
  pole field originates (FW §3.1). Correct.

- **LatticeDispatcher re-scatter index alignment**
  (`src/LatticeDispatcher.jl:338-345,447-462`). The re-scatter iterates
  `_dilate_one(gated.field_mask)` in `for j in 1:ny, i in 1:nx`
  column-major order, matching `EdgeGatedSolve._solve_targets`'s loop
  order (`src/EdgeGatedSolve.jl:222`). `_dilate_one(mask)` (reads
  original `mask`, single pass) is equivalent to the final solve's
  `_dilate(field, 1)` (single ring pass), so the `k`-th
  `pn_solution.grid_up[k]` lands on the same `(i,j)` it was solved at.
  No index drift. (This is outside the two assigned files but was on
  the cross-module path; verified correct.)
