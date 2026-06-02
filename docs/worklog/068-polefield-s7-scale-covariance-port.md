# Worklog 068 — scalar `PoleField` ported to the §S7 scale-covariance fix (C3)

**Date**: 2026-06-02
**Author**: Claude Opus 4.8 (1M) — orchestrated (Opus coding subagent)
**Bead**: `padetaylor-bez` (sweep C3, P1 bug).
**Scope**: port the landed vector §S7 far-root fix to the scalar twin
`src/PoleField.jl`; the vector side was fixed in commit `19b48ba`
(ADR-0026 Amendment 6/7) but `PoleField.jl` was never updated.

> **Result**: scalar `_extract_poles_core` now filters and orders by
> z-plane distance against the `h`-independent step ceiling `h_max`,
> at parity with `VectorPoleField.extract_poles_shared_q`. New
> closed-form-oracle test **PF.4.1** RED→GREEN; both halves of the
> atomic fix mutation-proven (M5 filter, M6 sort key). `polefield_test.jl`
> **35/35 GREEN** (was 29).

## Root cause (sweep C3, worklog 065)

`src/PoleField.jl:144` filtered Padé denominator roots by the per-node,
scale-DEPENDENT magnitude `|t*| ≤ radius_t` (rescaled-`t` space). Since
`t = (z − z_ctr)/h`, a pole at a fixed z-distance `D` from a node sits at
`|t*| = D/h`; under a varying / adaptive `h`, shrinking `h` pushes a
genuine pole past the fixed `|t*|` window and it is silently dropped →
intermittently empty / sparse pole fields. The cluster **sort key**
(`:147`, `RT(abs(t))`) carried the same heresy: ordering by `|t*|` crowns
a coarse far node (small `|t*|`) over a fine near one as the cluster
representative — the worse-placed sighting wins. ADR-0026 Amendment 7
treats filter + sort key as **one atomic fix**.

## Fix (`src/PoleField.jl`, port of `VectorPoleField.jl:259-339`)

1. Recover `h_max = isempty(hs) ? zero(RT) : maximum(abs, hs)` — the
   walk's `h`-INDEPENDENT step ceiling; `radius_z = radius_t · h_max`.
2. Far-root filter on z-plane distance: `z_dist = abs(h*t); z_dist ≤ radius_z`
   (was `abs(t) ≤ radius`). For a near-uniform walk `h ≈ h_max`, so the
   test reduces EXACTLY to the legacy `|t*| ≤ radius_t` — backward
   compatible (PF.1.*/2.*/3.* unchanged).
3. Cluster sort key = `RT(z_dist)` (was `RT(abs(t))`) → the representative
   is the z-plane-closest node. `sort!(...; by = c -> c[2])` unchanged.

Net code-LOC ≈ +5; the rest is literate-programming prose. No public-API
or default-`radius_t` change.

## Verification

- **PF.4.1** (new): stitched coarse-`h`(0.5)/fine-`h`(0.05) `PadeTaylorSolution`
  on the equianharmonic-℘ problem (`u''=6u²`, closed-form lattice
  `℘_pole(0,0)=1`). At `min_support=4`: (a) the genuine pole is recovered;
  (b) the representative is the z-closest node (err 2.94e-9, not the
  min-`|t*|` far node's 4.57e-7); (c) no spurious in-`in_box` poles.
- **Mutation-proven**: M5 (revert filter → `abs(t) ≤ radius`) ⇒ PF.4.1(a)
  RED (cluster support 3 < 4, pole dropped); M6 (revert sort key →
  `abs(t)`) ⇒ PF.4.1(b) RED (representative drifts to 4.57e-7). Restored.
- `julia --project=. test/polefield_test.jl` → **Pass 35 / Total 35** (13.5 s).

## Note (not a bug)

On a 1D real-axis trajectory the widened z-radius lets far nodes also
"see" genuine off-region poles and place them poorly (`-0.22 ± 2.06i`-class
candidates) — the same geometry the file header and PF.1.2/2.2 already
handle via `in_box` clipping. PF.4.1(c) is therefore scoped to the covered
box, consistent with the rest of the suite. Inherent to 1D geometry, not
introduced by the fix; multi-node-fan path networks (PF.1.2) unaffected.

## References

- `src/PoleField.jl`; `src/VectorPoleField.jl:223-242,259-339`.
- `docs/adr/0026-vector-resilient-walk-dense-targets.md` (new Amendment 11).
- `docs/worklog/065-bug-sweep-verify-and-fix-plan.md` (C3 confirmation).
