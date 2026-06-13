# Worklog 074 — Bug-fix campaign (epic `padetaylor-l7yt`)

**Date:** 2026-06-13. **Outcome:** 4 of the 6 confirmed library bugs surfaced by
the corpus-v2 (`p1v0`/`25og`) + test-hardening (`krgy`) sweeps are FIXED, committed,
and pushed. The 2 hardest (`61um`, `fzse`) remain, with `61um`'s correct fix
direction now established (the audition's recommendation was proven wrong — see
below). Expected-broken count **10 → 6** (reconciled in `quality_gate.sh` +
`CLAUDE.md` + bd memory `corpus-v2-expected-broken-count`).

## Orchestration shape (Rule 7-safe)

Because Rule 7 forbids parallel Julia, the campaign ran as: **one read-only
fix-audition workflow** (6 Opus agents, one per bug — confirm bug, audition ≥2
fixes, design oracle + test battery; NO `julia`), then **strictly serial Julia
execution** in the main loop. A **minimal campaign env** at `/tmp/ptcampaign`
(core PadeTaylor + Supposition + GenericLinearAlgebra/Polynomials/SpecialFunctions/
TaylorSeries, NO Makie/DelaunayTriangulation) runs single test files in ~10–50 s
and dodges the SIGTERM-prone diagnostics-ext precompile (bd memory
`full-pkg-test-got-sigterm…`). To recreate it:
`julia -e 'using Pkg; Pkg.activate("/tmp/ptcampaign"); Pkg.develop(path="."); Pkg.add(["Supposition","GenericLinearAlgebra","Polynomials","SpecialFunctions","TaylorSeries"])'`.

## The four fixes

| Bug | Commit | Fix | Verification |
|-----|--------|-----|--------------|
| `53tu` | `8a3c3c5` | BVP+VectorBVP callables guard the documented pre-image disc `abs(t_star) ≤ 1+100eps` (was `real(t*)`-only → silent off-oblique-segment extrapolation) | 3 markers→`@test_throws`, 2 edge testsets, 3 Supposition props, **triple mutation-proof** (M1/M2 loci-independent + M3 fuzz), 7-file regression |
| `xhjw` | `b0227e3` | direction-aware `solve_pade` driver+callable (`dir=sign(z_end-z_start)`) — descending spans were a silent degenerate one-node trajectory (Rule 1) | wolframscript ℘ oracle, 25 assertions, **4-component mutation-proof** (loop/clamp/scan/guard), 8-consumer regression (624) |
| `jznu` | `3f96794` | FIX-3: scope the SharedPade↔scalar:svd "bit-for-bit" claim to the full-rank regime; pin both degenerate-regime counterexamples as SP.1.8 masters | docstring + tests only (no src logic), MT-A (fallback-specific) + MT-B (main path) mutation-proven |
| `q0yq` | `7712341` | inter-slice **coherence guard** in the IVPBVP callable — `jump > 3·max(min(\|w_a\|,\|w_b\|),1)` (a divergent/corrupt slice). Docstring promised a glue throw the body never performed | CBV.9 marker→`@test_throws` + GM2, MT1/MT2 mutation-proven, ffw_fig_5 (23) + ivp_bvp_hybrid (65) regression |

## Two audition errors caught (Rule 3 skepticism paid off)

1. **`xhjw` snap** — the audition (and my first draft) added a `state.z = z_end`
   snap "to kill float residue." A sweep of **570 spans** found ZERO where it
   changed behaviour (the clamp's `h_step == gap` already lands exactly), and
   MT-5 reddened nothing. It was **removed** as verified-redundant untestable
   code (Rule 9). The scan fix was ALSO invisible on smooth ℘ (any segment's
   (15,15) Padé extrapolates to ~2e-13 even at t≈9), so the load-bearing test is
   a descending **pole-bridge** at z=−1.8 where only the correct segment values
   it right.
2. **`q0yq` calibration** — the audition's `10·max(scale)` would **never fire**:
   CBV.9's +1000 injection inflates the LARGER bracketing value, so
   `max≈1078 → 10·max≈10780 > jump≈1006`. Measured the real ratios (honest
   `jump/|w|≈0.093` vs catastrophe `≈14`, a 150× gap) and used a **`min`-based**
   bound, `3×`, floored at 1.

## `61um` — the audition's fix is geometrically WRONG (proven)

The audition recommended chord-midpoint **bisection** to "recover the exact
continuous winding." Probe (`winding_delta` on CWD.5's geometry, branch 0,
z_old=0.1, z_new=0.1·e^{i·1.1π}):

```
winding_delta (2-arg principal value)      = -0.9π
chord-bisection  n = 2, 4, 16, 64          = -0.9π   (NEVER converges)
trajectory subdivision (sample the arc)    = +1.1π   ✓ recovers truth
```

A straight chord subtends ≤π about an exterior point, so bisecting the **chord**
returns the minor arc at every depth. A true subtended angle of +1.1π REQUIRES a
**curved** path — so the bug is *coarse sampling of the curved ODE trajectory*,
and the fix must sample the **actual trajectory** (reachable in the walker via
the Padé / Stage-2 grid), not the straight chord. The CWD.5 `@test_broken`
(asserting the 2-arg `winding_delta` returns +1.1π) is **information-theoretically
impossible** to satisfy and must be REFRAMED, not flipped.

**Correct fix direction (for the next session):** architectural — accumulate
sheet-winding along the densely-sampled trajectory in the `PathNetwork` walker /
`step_sheet_update`, OR fail loud when a step is too coarse to resolve. Inheritors:
`accumulate_winding`, `sheet_index`, `step_sheet_update`, and the walker
(`src/PathNetwork.jl:620-633`). Live real-walk fixtures: CWD.5
(`corpus_winding_test.jl:183`), CBr.3 (`corpus_elementary_branch_test.jl:158`),
CPN.7.2/7.3 (`corpus_pathnet_winding_test.jl`). Start by tracing where the walker
calls `step_sheet_update` and whether the Stage-2 trajectory is reachable there.

## Follow-on beads filed

- `padetaylor-rsln` (P3) — BVP `extrapolate=true` documented but unimplemented.
- `padetaylor-tqvz` (P3) — BVP disc-vs-strict-segment contract decision.
- `padetaylor-x0p0` (P2) — VectorProblems.jl:262 has xhjw's identical
  forward-only bug (distinct adaptive-stepper fix).
- `padetaylor-sn9a` (P3) — v2 IVPBVP PFS-vs-BVP derivative glue match (FFW
  md:247, needs BF-256).
- Plus `61um` and `fzse` remain OPEN on the campaign epic `padetaylor-l7yt`.

## Discipline notes

- Every fix carried an INDEPENDENT oracle (wolframscript ℘ regenerated for xhjw
  and the q0yq pole-past oracle; empirical capture for jznu's masters — the bead's
  "equal as power series" phrasing was itself wrong) — Rule 5, verified not copied.
- Every load-bearing assertion was mutation-proven (perturb → RED → restore
  byte-clean, `git diff src/` confirmed empty).
- `IVPBVPHybrid.jl` (797 LOC) is a PRE-EXISTING Rule 6 violation; the q0yq guard
  added to it (necessarily) does not address that — a module split is its own task.
