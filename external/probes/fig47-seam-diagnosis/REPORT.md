# Probe — FW 2011 Fig 4.7 spurious "seam" (grain boundary): diagnosis & root cause

**Date:** 2026-06-29. **Beads:** `padetaylor-vwgl` (P0 root cause), `padetaylor-sny7`
(P1 regression gate). **Full write-up:** `docs/worklog/077-fig47-path-dependence-seam-rootcause.md`.

## What was found

Regenerating FW Fig 4.7 surfaced a spurious **discontinuity / grain boundary**
in the pole field (worst in panels **(e)** near-tritronquée and **(a)**
`u(0)=u'(0)=0`, present in all six): two lattice-orientation domains meet at a
sharp seam in a region that should be a single coherent quasi-lattice, plus
scattered spurious poles (incl. (f)'s forbidden sectors). It has been present
since project start (worklog 014, "y=0 row discontinuity") and never pinned.

## Root cause — path-DEPENDENCE of the monolithic walk (NOT the Stage-2 stitch)

`path_network_solve` enforces **no path-independence**. The walk
(`src/PathNetwork.jl:481` shuffled tree; `:605` per-node forward propagation)
reaches spatially-adjacent regions via *different* tree branches that accumulate
*different* IVP error in the smooth/unstable sectors (FW md:401). Where they
meet, the field jumps and the per-node pole lattices lose registration.

## Decisive evidence (run these scripts, `julia --project=. <file>`)

| script | what it measures | result |
|---|---|---|
| `p2_jumpmap_and_seeddiff.jl` | nearest-vs-2nd-nearest field jump **and** two-seed field diff, panel (e) `[-50,50]²` | jump max **1.37e-4** (stitch is FINE); two-seed diff **35% of cells >1%, max 1.0** (grossly path-dependent) |
| `seed_poles.jl` | extract poles at rng_seed 0 vs 42 | pole **count** seed-dependent (**8159 vs 7987**); the seam **moves** with seed |
| `pole_groundtruth_weierstrass.jl` | Weierstrass-℘ (known lattice) recall/precision/symmetry | all **1.000** — the extraction machinery is clean; the bug is the WALK, and ℘ is a misleading proxy (doubly-periodic, closes its loops) |
| `pi_fivefold_symmetry.jl` | PI(a) conjugate + 5-fold symmetry | 0.997 / 0.998–1.000 — the gross sector structure is **physical**; only the fine grain boundary is spurious |

Amplifier = window size: two-seed diff is **8.5e-3 at `[-25,25]²`** but **1.0 at
`[-50,50]²`** (`p2_jumpmap_and_seeddiff.jl` with `HALF=25` vs `HALF=50`). FW avoid grossness via a 5×5 composite of 25 independent 20×20
runs (FW md:147) and by only caring about pole fields (FW md:208 admits the
"ridges in flat areas from different paths").

## The diagnostic that catches this (the reusable lesson)

The **path-independence (two-seed) test**: re-solve with a different `rng_seed`
and require the field/poles to agree. Global stats, pole-recall, symmetry, and
proxy (℘) problems all MISS it; only interrogating the FIELD and testing the
path-independence INVARIANT exposes it.

## Status

Root-caused only — **no fix applied** (per directive). Cure = enforce
path-independence (smooth-region BVP fill / bounded-window composite / K-nearest
consensus, bead `padetaylor-6b5`) gated by the `sny7` seed-invariance test.
