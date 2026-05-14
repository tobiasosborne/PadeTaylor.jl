# Worklog 027 — `PoleField.extract_poles`: pole locations from a solved path-network

**Date**: 2026-05-14 (continues worklog 026)
**Author**: Claude Opus
**Bead**: `padetaylor-xvf` — closed.
**Scope**: The capability gap flagged at the end of HANDOFF.md's
"Session 2026-05-14" entry: FW 2011 Fig 4.7/4.8 are *pole-location
plots*, and nothing in the codebase reads pole positions back out of a
solved path-network. `src/PoleField.jl` (new) closes that gap.

> **Take-home**: `extract_poles(sol::PathNetworkSolution)` reads the
> per-node Padé store, maps each stored denominator's roots back to the
> `z`-plane (`z = z_v + h·t*`), and clusters what survives two filters
> into one estimate per physical pole. Verified against the
> equianharmonic Weierstrass-℘ lattice (FW §5.1.1) — in-region poles
> recovered to ~1e-8, well inside the figure-catalogue's 1e-6 Float64
> spec. Test suite **1487 → 1509 GREEN** (+22 from `polefield_test.jl`).

## Ground truth read first

- `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:147` — Fig 4.7a–f is a 5×5
  composite of 25 independent path-network runs; Fig 4.8 a composite
  of 18. Each is a pole-location scatter plot.
- `:281-318` — the Weierstrass-℘ test problem `u'' = 6u²`, whose
  solution `u(z) = ℘(z + c₁; 0, c₂)` has analytically known
  second-order poles on a rhombic lattice; real-axis spacing
  `x = 1 + 2ωk`, `ω = Γ(1/3)³/(2^{13/6}π) ≈ 1.363`. This is the
  oracle — no pole-finder oracle is needed, the lattice is closed-form.
- `docs/adr/0004-path-network-architecture.md` + `src/PathNetwork.jl`
  module docstring — the per-node Padé store (`visited_pade`,
  `visited_h`, `visited_z`) `extract_poles` reads, and the
  `t = (z - z_v)/h` rescaled-variable convention.
- `src/StepControl.jl:240` `step_pade_root` — the existing
  `roots(Polynomial(P.b))` idiom; `extract_poles` follows it.

## `src/PoleField.jl` (new, ~180 LOC incl. docstrings)

A rational function's poles are its denominator's zeros. At visited
node `k` the stored Padé `P(t) = N(t)/D(t)` is in `t = (z - z_v)/h`, so
its poles in the `z`-plane are `z = z_v + h·t*` for each root `t*` of
`D`. `extract_poles` takes the union over the visited tree and filters:

  1. **`radius_t`** (default 5.0) — drop roots with `|t*|` beyond a few
     canonical steps; a local Padé does not place distant singularities.
  2. **`min_residue`** (default 1e-8) — drop Froissart doublets: a
     denominator root that is also nearly a numerator root has residue
     `|N(t*)/D'(t*)|` at the noise floor (the Float64 classical-Toeplitz
     path has no SVD rank guard — ADR-0005).
  3. **cross-node support** (`min_support`, default 3) — the decisive
     filter. Cluster surviving candidates greedily in increasing `|t*|`
     (best-placed node wins the representative slot); a cluster is a
     physical pole only if ≥ `min_support` *distinct* nodes
     independently land a root in it. Node-local artefacts never
     accrue cross-node support — this is the same "composite of
     independent runs" robustness FW's Fig 4.7 itself relies on.

## TDD record (spec-from-scratch, RED → GREEN → mutation-proven)

The two non-obvious findings, both surfaced by the tests going RED:

**Finding 1 — second-order poles split, and the dedup must cope.**
℘ has *double* poles; a finite-order Padé represents a double pole as
two near-coincident simple roots straddling the true location. The
split shrinks with `|t*|`: from the IC node at `t = 2` it is sub-1e-6
(z = 1 recovered to 4.6e-7), but from a far node at `t ≈ 4.7` the
off-axis poles split by ~0.06. The first cut (greedy keep-first dedup
at `atol = 1e-3`) kept one split half → 8e-4 errors. Fix: cluster at
`atol = 0.1` (wider than any split a *well-placed* node produces,
narrower than the 2.36 inter-pole spacing) and let the smallest-`|t*|`
candidate be the representative. The well-placed nodes' estimates are
already coalesced to ~1e-8, so the representative is accurate and the
straggling far-node halves merely add support.

**Finding 2 — "spurious" poles were genuine edge poles.**
With the support filter in, PF.1.2's "no spurious" check still failed:
8 reported poles sat 4e-6 … 0.05 off the lattice. Diagnosis: every one
is a *genuine* ℘ lattice pole lying just outside the covered grid box
(the `n = ±2` rows at `Im ≈ ±4.72`, the `m = ±2` columns) — correctly
identified as real (support 3–15) but placed poorly because only
distant nodes can see them. This is not an algorithm bug; it is the
honest behaviour, and matches FW (Fig 4.7 is "displayed over the
region [-50,50]²" — they clip to a window). Fix is test-side: the
accuracy assertions clip to the covered grid box. `extract_poles`
itself reports everything ≥ `min_support` nodes agree on; the caller
(a Fig 4.7 reproduction) clips to its display window.

**Tests** (`test/polefield_test.jl`, 22 assertions):
  - **PF.1.1** — single-node (IC-only) network, `min_support = 1`:
    surfaces the nearest pole `z = 1` to ≤ 1e-6 (the (15,15) double-
    pole floor at `t = 2`).
  - **PF.1.2** — 2D grid over a 7-pole region: (a) no spurious poles
    in-box, (b) all 7 expected lattice poles recovered to ≤ 1e-6
    (actually ~1e-8), (c) conjugate symmetry.
  - **PF.2.1** — degenerate network (constant denominators) → `[]`,
    no throw on the empty `roots` call.
  - **PF.2.2** — the cross-node support filter is load-bearing:
    `min_support = 1` yields strictly more, off-lattice "poles".

**Mutation-proven** — M1 (drop the `h` rescale in `z_v + h·t`), M2
(disable `radius_t`), M3 (disable the support filter). Each was applied
to `src/PoleField.jl`, the suite confirmed RED (3+2 / 43 / 49
failures respectively), the mutation reverted. Procedure in the test
file footer.

## Frictions

  - **`min_support = 3` default vs. single-node networks.** A network
    with one visited node (the IC-only grid in PF.1.1) cannot satisfy
    a cross-node support of 3 — `extract_poles` would return `[]`. The
    `min_support` kwarg is the escape hatch (`= 1` reads a single
    node's Padé directly); documented in the docstring. A genuine v1
    corner: poles near the edge of a *real* Fig 4.7 run will have low
    support and may need a lower `min_support` — left as a caller
    knob, not auto-tuned.

## Deferred / next

  - The Fig 4.7/4.8 *figure reproduction* itself (a `figures/`
    script consuming `extract_poles`) is the natural follow-up — it
    will need the per-IC PI runs and the display-window clip described
    above. Not filed as a bead yet; the HANDOFF "suggested next work"
    already names it.

## Hard-won lesson

**A RED test that looks like an algorithm bug can be a wrong
acceptance criterion.** Finding 2's "spurious poles" were the
algorithm working correctly — reporting real poles the network under-
resolved. CLAUDE.md Rule 2 ("all bugs are deep") cuts both ways: the
investigation showed the *test* assumed something false (every
reported pole is well-placed), not the code. The fix was to make the
acceptance match the physics (accuracy is a function of node
coverage), not to relax a tolerance. Cross-checking each "spurious"
location against the closed-form lattice — `℘_pole(m, n)` for larger
`m, n` — is what made the distinction visible.
