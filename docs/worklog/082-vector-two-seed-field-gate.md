# Worklog 082 — vector two-seed FIELD gate (bead `padetaylor-lg4y`)

**Date**: 2026-09-29
**Scope**: Measure two-seed field agreement of `vector_path_network_solve` on
the shipped P_I^(2) wedge march, then port the scalar seam gate
(`test/field_seam_test.jl`) to the vector pipeline.  Test-only delivery: new
file `test/vector_field_seam_test.jl` (161 lines incl. header, GREEN at
12 pass + 1 broken, ~44 s wall).  No `src/` change.

> Take-home: two SEEDED orderings (0, 42, 7) agree field-wide to ≤ 5.8e-6
> relative, including 96-100 far-wedge nodes — no seam between random
> orderings.  But the DEFAULT ordering (`rng = nothing`) at `h = 0.1` — the
> exact walk the shipped P_I^(2) figure runs (`figures/_kkg_pi2_helpers.jl:339`)
> — disagrees with all of them by up to 2.18 relative on 65 far-wedge nodes
> (|z| ∈ [5.82, 8.0]).  Four independent walks (default@h=0.05, seed 7@h=0.05,
> seed 7@h=0.1, seed 42) side with seed 0 there.  This is a real vector-walk
> field defect; it is pinned as `@test_broken`, not papered over.

## Ground truth opened

- `test/field_seam_test.jl:1-130` — the scalar FSEAM gate (non-triviality +
  bidirectional pole-set match ≥ 0.90).
- `src/VectorPathNetwork.jl:451-470,642-662` (signature; `rng` kwarg),
  `:723` (the shuffle), `:819-838` (per-step `_select_wedge` + node push),
  `:401-412` (`VectorPathNetworkSolution` fields).
- `src/VectorPathNetworkStage2.jl:41-127` — the B1 true-radius gate: a
  `fine_grid` point outside the nearest node's honest disc is `NaN`.
- `test/vector_path_network_test.jl:853-885` — VPN.4.3 (Riccati, one pole).
- `figures/_kkg_pi2_helpers.jl:130-160,293-356` — the P_I^(2) fixture.
- `figures/test_kkg_pi2_surface.jl:938-1040` — PI2S.10, the existing VC-10
  two-run POLE check (median 0.35 disagreement, pole spacing 0.69).

## 1. Pole-set measurement (the scalar metric) — inconclusive on this fixture

Fixture: `kkg_pole_field`'s march (seed z = -3, order 24, h = 0.1, 20-target
fan); `extract_poles_shared_q(radius_t = 5, cluster_atol = 0.2, min_support = 2)`.

| run | visited nodes | poles | walk time |
|---|---|---|---|
| default (`rng = nothing`) | 795 | 85 | 12.0 s (incl. compile) |
| seed 0 | 848 | 88 | 6.3 s |
| seed 42 | 824 | 94 | 5.7 s |
| seed 7 | 812 | 84 | 5.9 s |

Pole-set match (fraction of A within tol of some pole of B, both ways):

| pair | @0.5 | @0.05 | @1e-4 | Hausdorff |
|---|---|---|---|---|
| 0 vs 42 | 0.807 / 0.809 | 0.057 / 0.053 | 0 / 0 | 1.18 |
| default vs 0 | 0.788 / 0.841 | 0.141 / 0.136 | 0 / 0 | 1.09 |
| default vs 42 | 0.847 / 0.819 | 0.141 / 0.128 | 0 / 0 | 0.929 |
| 0 vs 7 | 0.864 / 0.881 | 0.170 / 0.179 | 0 / 0 | 1.08 |
| 42 vs 7 | 0.915 / 0.905 | 0.170 / 0.190 | 0 / 0 | 1.08 |

Seeds 0 and 42 score 81 % — the scalar FSEAM "RED" regime (monolithic 77.6 %)
— yet (§2) their FIELDS agree to 4e-7.  So on this fixture the shared-Q
pole extraction jitter (0.05-0.5) swamps any walk signal: the pole-set metric
cannot discriminate a seam from extraction noise.  This is why the vector
gate compares the field instead.  (Side note, out of scope: the fixture's
docstrings, `figures/_kkg_pi2_helpers.jl:148` and
`test/kkg_pi2_figure_test.jl:47`, say "389 nodes, 21 poles"; measured today
795 nodes / 85 poles for the same default walk — stale text.)

## 2. Field measurement — the gate's metric

Run A with `fine_grid = B.visited_z`, compare A's Stage-2 value against B's
stored `visited_y`, relative to max(1, |y_B|), all four components; only
nodes inside A's honest discs (finite) count.

| A at B's nodes | covered | covered far (|z|≥5.5) | max rel (u, u', u'', u''') |
|---|---|---|---|
| seed 0 @ seed 42 | 158/824 | 96 | 3.79e-7, 2.30e-6, 1.43e-6, 5.84e-6 |
| seed 42 @ seed 0 | 164/848 | 100 | 3.73e-7, 3.36e-6, 7.48e-7, 3.42e-6 |
| default @ seed 0 | 215/848 | — | max 2.18; 65 nodes > 1e-4 |
| seed 0 @ default | 223/795 | — | u max 0.25; 60-63 nodes > 1e-4 |
| default @ seed 42 | 266/824 | — | u max 0.22; 84-89 nodes > 1e-4 |

Localisation of the default-vs-seed-0 disagreement: all 65 bad nodes have
|z| ∈ [5.82, 7.99], arg ∈ [-0.06, 0.52] (clusters near 5.3-5.8 + 1.5-2.9i,
6.0-7.0 + 3.5-3.8i, 7.6-7.7 + 2.0i, 5.8-6.0 - 0.3..0i, 7.97 - 0.03i).  The
inner wedge (|z| < 5.5, 77 covered nodes) agrees to 9.48e-8.

## 3. Which walk is wrong — independent cross-check

Evaluated at seed 0's far-wedge nodes (|z| ≥ 5.5), comparing an independent
walk ("oracle") against seed 0 and against the default h = 0.1 walk:

| oracle walk | nodes covered by both | max rel vs seed 0 | max rel vs default | n(vs default > 1e-4) |
|---|---|---|---|---|
| default order, h = 0.05 | 79 | 7.65e-5 | 2.94e-1 | 27 |
| seed 7, h = 0.05 | 78 | 2.14e-5 | 2.40e-1 | 24 |
| seed 7, h = 0.1 | 75 | 4.12e-7 | 2.40e-1 | 23 |

Every independent walk — including the SAME default order at a finer step —
sides with seed 0.  The outlier is the default-order h = 0.1 walk.  P_I^(2)
has the Painlevé property (single-valued meromorphic solution), so this is
walk error, not a sheet effect.  Root cause NOT investigated (time box); the
likely class is the one worklog 077 found for the scalar walk: a far-wedge
subtree that passes close to a pole and accumulates error.  Consequence: the
far-wedge part of the shipped Panel B (poles at |z| ≳ 5.8) is computed from a
field that is wrong at O(0.1-1) relative.

## 4. The gate (`test/vector_field_seam_test.jl`)

- VFSEAM.1 non-triviality: seeds 0/42 visit different nodes; > 100 nodes each.
- VFSEAM.2 both directions: covered ≥ 100, far-covered ≥ 50, max rel over all
  four components ≤ 1e-4 (measured ≤ 5.84e-6; seam regime 0.2-2.2).
- VFSEAM.3 known seam: `@test_broken` default-vs-seed-0 max rel ≤ 1e-4;
  pinned TRUE facts: covered ≥ 100, min |z_bad| ≥ 5.5, inner-wedge max rel
  ≤ 1e-4, bad count ≤ 100 (measured 65).

## 5. Mutation proof (executed; `src/` restored, `git diff src/` empty)

- **M-frozen-seed** — `src/VectorPathNetwork.jl:723` shuffle line replaced by
  a comment.  RED: VFSEAM.1 `w0.visited_z != w42.visited_z` fails (identical
  trees); VFSEAM.2 passes trivially at 4.9e-16 on 795/795 nodes (the gaming
  VFSEAM.1 exposes); VFSEAM.3's `@test_broken` unexpectedly passes (Error).
  Summary: 11 pass, 1 fail, 1 error.
- **M-drift** — after the `_select_wedge` call (`:819-821`) insert
  `y_new = y_new .* (1 + 1e-5)` (path-length-dependent error).  RED:
  VFSEAM.2 fails both directions (max rel 1.75e+2 and 8.11e+1 vs 1e-4);
  VFSEAM.3 localisation fails (min |z_bad| = 1.82; inner max rel 20.0).

## Follow-up needed (for the orchestrator to file)

- A bead for the default-order far-wedge field defect (§3): root-cause the
  default h = 0.1 walk's far-wedge subtree; decide whether the P_I^(2) figure
  should switch ordering / step or grow a two-run guard.  When fixed, delete
  the VFSEAM.3 `@test_broken` marker (it will report as an Error).
- The suite's expected `@test_broken` count rises from 2 to 3 once this file
  is registered in `test/runtests.jl`.
