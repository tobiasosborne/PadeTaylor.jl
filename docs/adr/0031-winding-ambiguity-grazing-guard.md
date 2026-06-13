# ADR-0031 — Riemann-sheet winding: `winding_delta` is correct; fail loud on branch-grazing steps

**Status**: **accepted** — 2026-06-13. Resolves bug `padetaylor-61um` (epic
`padetaylor-l7yt`, the bug-fix campaign).
**Date**: 2026-06-13
**Beads**: `padetaylor-61um`. Supersedes the framing in the corpus headers of
`test/corpus_winding_test.jl` (CWD.5), `test/corpus_elementary_branch_test.jl`
(CBr.3), `test/corpus_pathnet_winding_test.jl` (CPN.7.3).

## Context

`padetaylor-61um` was filed as a Rule-1 (fail-loud) gap:
`SheetTracker.winding_delta(z_old, z_new, branch)` computes the per-step signed
angle as the difference of principal-value `arg()`s about a branch, normalised
to `(-π, π]`. Its documented precondition — that a single step subtend `|Δθ| < π`
— was unenforced. The corpus fixtures CWD.5 / CBr.3 / CPN.7.3 asserted (via
`@test_broken`) that for a step whose *true* subtended angle is `+1.1π` the
function "silently loses a full `+2π` revolution" (returning `−0.9π`), corrupting
every downstream Riemann-sheet index in `BranchTracker.step_sheet_update` →
`PathNetwork`'s walker (`src/PathNetwork.jl:627`). The bead proposed detecting a
single step whose true subtended angle is `≥ π` (via "the signed area / atan2 of
the relative chord") and throwing.

## What the audition + independent oracles established

A read-only audition (6 investigators + 4 adversarial verifiers, each using
`wolframscript` for independent ground truth) and an independent main-thread
`wolframscript` check **converged** on a correction to the bead's premise. The
load-bearing facts (all re-verified, re-runnable):

1. **The walker realises a STRAIGHT z-step.** `z_new = z_cur + h_step` is a
   single complex displacement (`src/PathNetwork.jl:944-949`); the sheet
   bookkeeping (`segment_crosses_cut`, `winding_delta`) operates only on the
   straight chord `[z_cur, z_new]`. There is no curved arc in the realised walk.

2. **A straight chord can NEVER subtend `|Δθ| ≥ π` about an exterior branch**
   (proven over 2000 random configs: `max |straight-chord winding| < π`). So
   `winding_delta` returns the chord's *exact* true winding and **loses
   nothing**. For CPN.7.3's realised FINE chord the straight-chord winding is
   `−0.94π`, equal to `winding_delta` to `< 1e-15`; the `sqrt(1−z²)` solution
   does **not** flip sheets along that chord. `visited_sheet = [−1]` was the
   geometrically-correct label for the realised straight step.

3. **The `+1.1π` the `@test_broken` markers expected is the MAJOR-arc angle of a
   CURVED path** the walker never realises, and is information-theoretically
   unrecoverable from two endpoints (which way a path winds between its
   endpoints is not a function of the endpoints). Chord-bisection never
   converges to it (worklog 074:47-74).

4. **The bead's proposed detector cannot fire.** "Subtended angle `≥ π`" computes
   the chord's own winding, which is `< π` in magnitude for every straight step
   (`0.94π` for FINE) — so it would *never* trigger on the case it was designed
   for.

The genuine defect is therefore **not** a lost revolution but **silent
ill-conditioning**: when a single straight step's chord *grazes* a tracked branch
(passes very close relative to its length), the discrete two-node chord cannot
resolve which side of the branch the true ODE continuation passed, and the
`sign(winding_delta)` sheet bump is unreliable — yet it was returned with no
diagnostic (a Rule-1 violation).

## Decision

1. **`winding_delta` stays a pure, unguarded principal-value primitive.** It is
   mathematically correct for the straight chord it sees; it is also consumed by
   `accumulate_winding` / `sheet_index`. Only its docstring is corrected (the
   "loses a revolution" precondition prose is replaced by the geometric truth).

2. **Add a fail-loud winding-ambiguity guard in
   `BranchTracker.step_sheet_update`** — the walker's sole sheet-bump site. When
   a step crosses a branch cut, compute the chord's closest approach to that
   branch relative to the chord length (the *grazing ratio*); if
   `min_dist / |chord| < graze_tol` throw a `DomainError` naming the branch, the
   distances, and a "shorten `h` / densify nodes near the branch" suggestion
   (Rule 1). Equivalently this is the `|Δθ| → π` zone.

3. **Reframe (do not flip) the three `@test_broken` markers.** Each premised a
   geometrically-impossible `+1.1π`. They become: a live pin that `winding_delta`
   returns the *correct* chord winding, plus a `@test_throws DomainError` on the
   grazing `step_sheet_update` (unit) and, for CPN.7.3, on the real walk.

## Why the grazing ratio (not the bead's subtended angle, not `|Δθ|`)

The detection signal must (a) fire on the genuinely ambiguous step and (b) never
trip a well-conditioned walk. Measured grazing ratio at **every** crossed-cut
step across the in-suite `cross_branch=true` corpus (`corpus_winding`,
`branch_tracker`, `path_network_branch`, `sheet_aware_stage2`,
`corpus_pathnet_winding`, `ffw_fig_2/3/7`; 38 steps total):

| step | grazing ratio | \|Δθ\|/π | verdict |
|---|---|---|---|
| CPN.7.3 FINE | **0.0605** | 0.917 | catch (ambiguous) |
| CPN.7.2 COARSE | **0.2045** | 0.556 | keep (safe, returns `[1]`) |
| ffw_fig_7 | 0.71, 0.81 | ≤0.37 | safe |
| sheet_aware / path_network rings | ≥0.98 | ≤0.13 | safe |
| ffw_fig_2/3 | ≥2.22 | ≤0.13 | safe |

The grazing ratio separates the lone ambiguous step from all safe steps by
**3.4×** (only one safe step, CPN.7.2, lies below 0.71). It is the root-cause
signal: a Padé built over a chord that grazes a branch cannot resolve that
branch, so the winding is ill-posed regardless of endpoint radii. The `|Δθ|`
proxy is radius-blind and separates these two fixtures by only 1.65×; the bead's
"subtended `≥ π`" never fires.

**The read-only audition's proposed `tau = 0.2` was wrong** — it estimated
CPN.7.2 COARSE's ratio as `0.79` from the *ideal* ring geometry, but the
*realised* walker step has ratio `0.2045`, leaving only a 2% margin. The
threshold was set on the **measured realised distribution**, not the ideal.

## Threshold

`WINDING_GRAZE_TOL = 0.1` (a module const in `BranchTracker`; per-call override
via the `graze_tol` kwarg on `step_sheet_update`). It catches CPN.7.3 FINE
(`0.0605`) with `1.65×` margin, clears the nearest safe step CPN.7.2 COARSE
(`0.2045`) by `2.05×`, and every other in-suite walk step by `≥ 7×`.
Mutation-proven both directions: `tol → 0` reddens the FINE/CBr.3 throws;
`tol → 0.3` (above CPN.7.2) reddens the safe COARSE leg.

## Consequences

- **`@test_broken` count drops 6 → 3** (CWD.5, CBr.3, CPN.7.3 reframed). The
  `corpus-v2-expected-broken-count` bd memory + `scripts/quality_gate.sh` are
  updated in lockstep.
- No shipped `cross_branch=true` walk is affected (measured: all clear `tau` by
  `≥ 7×`, including the adaptive multi-branch ffw_fig_7 ζ-walk).
- The guard is a *practical* safety net, not a proof: no finite threshold can
  *guarantee* correctness from two endpoints (a small-`|Δθ|` step could in
  principle hide added revolutions on a pathologically curved true path). For the
  small straight steps the walker realises, grazing is the operative hazard.

## Open / deferred

- **Cut-crossing vs analytic continuation.** The sheet law
  (`sign(winding_delta)` on `segment_crosses_cut`) is a *cut-crossing
  topological counter*, which can diverge from the multivalued solution's *actual*
  sheet when the chosen `cut_angle` does not match the function's true branch cut
  (e.g. `cut_angle = π` vs `sqrt(1−z²)`'s rightward cut from `z = +1`). Settling
  whether the package's bookkeeping should track one or the other is left to a
  follow-up (no current consumer relies on the distinction). Filed as a note on
  `padetaylor-61um`'s follow-on.
- **Threading `graze_tol` through `path_network_solve`.** Currently a power user
  must call `step_sheet_update` directly to override the default; no in-suite
  walk needs it. Deferred (PathNetwork is already over the Rule-6 LOC cap).
- **Walker-level auto-subdivision near branches** (insert an intermediate
  visited node instead of throwing) is a larger architectural change; the throw
  with a "shorten `h`" suggestion is the v1 behaviour.

## References

- `src/SheetTracker.jl` `winding_delta` (corrected docstring);
  `src/BranchTracker.jl` `step_sheet_update` + `WINDING_GRAZE_TOL` (the guard).
- `src/PathNetwork.jl:944-949` (straight z-step), `:627` (the bump call site).
- `test/corpus_winding_test.jl` CWD.5, `test/corpus_elementary_branch_test.jl`
  CBr.3, `test/corpus_pathnet_winding_test.jl` CPN.7.2/7.3 (the reframed
  fixtures + mutation-proof footers).
- `docs/worklog/074-bugfix-campaign.md:47-74` (the prior session's analysis that
  first established "reframe, not flip").
