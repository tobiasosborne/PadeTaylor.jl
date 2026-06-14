# Worklog 075 — Bug-fix campaign CLOSE (epic `padetaylor-l7yt`): final 3 bugs

**Date:** 2026-06-13 → 2026-06-14. **Outcome:** the bug-fix campaign is COMPLETE.
The 3 bugs left after worklog 074 — `x0p0` (vector descending span), `61um`
(winding revolution-loss), `fzse` (extract_poles precision/recall) — are FIXED,
committed, pushed, and closed. **All 6 epic children done** (53tu/xhjw/jznu/q0yq
in worklog 074 + 61um/fzse here); `x0p0` was an xhjw follow-on sibling.
Expected-broken **6 → 3** (lockstep: `quality_gate.sh` + `CLAUDE.md` + bd memory).

## Orchestration (Rule 7-safe)

One **read-only audition workflow** (11 Opus agents: 6 investigators → 4
adversarial verifiers → 1 synthesizer, each using `wolframscript` for independent
oracles; NO `julia`), then **strictly serial Julia** in the main loop, one bug at
a time, in the minimal campaign env `/tmp/ptcampaign`. Per-bug: claim → fix →
test → mutation-prove → regression → commit/push/close. The audition's structured
output is in the run transcript; its conclusions were checked against my own
independent `wolframscript` computations (Rule 3) before any code.

## The headline lesson: TWO of the three bugs were NOT what the bead said

`61um` and `fzse` both turned out to be **mis-framed** — the bead's premise was
geometrically/empirically false, and discovering that *was* the fix. This is the
campaign's recurring shape (Rule 3 skepticism + independent oracles pay off):

- **`61um`** — the bead said `winding_delta` "silently loses a full revolution"
  on a step subtending `|Δθ| ≥ π`. Independent `wolframscript` (mine + all 3
  investigators + both verifiers, converging): the walker realises a **straight**
  z-step, and a straight chord can **never** subtend `|Δθ| ≥ π` about an exterior
  branch — so `winding_delta` returns the chord's *exact* winding and loses
  nothing. The `+1.1π` the markers expected is a *curved-path* major arc,
  unrecoverable from two endpoints. The bead's proposed "subtended-angle ≥ π"
  detector provably never fires.
- **`fzse`** — the bead said ~76% precision with 4 interstitial false positives.
  The calibration probe (vs the closed-form lattice *including out-of-window*
  poles): precision is **1.00** at the default; the "FPs" were genuine
  out-of-window poles + an artifact of the bead's own `cluster_atol = 0.4`
  tuning. The real defect was over-splitting only.

## The three fixes

| Bug | Commit | Fix | Verification |
|-----|--------|-----|--------------|
| `x0p0` | `1caeaf3` | dir-aware `vector_solve_pade` driver+callable (mirrors scalar `xhjw`; signed-h clamp threads the vector Jorba–Zou ceiling) — descending spans were a silent degenerate one-node trajectory (Rule 1) | wolframscript equianharmonic-℘ companion oracle, +35 VP.2.* assertions, **4-locus mutation-proof** (loop/clamp/scan/guard), regression VP 72/72 + VPO 27/27 |
| `61um` | `9f5021a` (ADR-0031) | `winding_delta` left correct (docstring fixed); **fail-loud grazing guard** in `BranchTracker.step_sheet_update` (`min_dist/|chord| < graze_tol`, default 0.1); CWD.5/CBr.3/CPN.7.3 markers **reframed** (not flipped) | threshold **measured** at every in-suite crossed-cut step (FINE 0.0605 vs safe-min 0.2045 vs all-others ≥0.71); mutation-proven both directions; all `cross_branch` walks incl. ffw_fig_2/3/7 GREEN |
| `fzse` | `4d97319` (ADR-0032) | **disjoint-support single-linkage self-merge** post-pass in `extract_poles` (`merge_atol`, default `h_max`) — collapses cross-node over-split froth without merging near-coalescent pairs | calibration probe; +PF.5.1; mutation-proven M7 (merge-off) / M8 (drop `isdisjoint`); 12 in-suite scalar consumers GREEN |

## Two calibration corrections the read-only audition got wrong (Rule 3)

1. **`61um` threshold.** The audition proposed `tau = 0.2` from the *ideal* ring
   geometry. Instrumenting `step_sheet_update` and measuring the **realised**
   walker steps: CPN.7.2 COARSE (a *safe* step) has grazing ratio **0.2045**, so
   `tau = 0.2` gives only a 2% margin. Re-set `tau = 0.1` (1.65×/2.05× margins;
   all other walks ≥ 7×).
2. **`fzse` gate.** The audition designed an agreement-spread + residue gate to
   kill the "false positives." The calibration showed there are **no** FPs at the
   default → no gate needed. And a distance-only self-merge (a candidate design)
   regressed 4 real tests by merging near-coalescent pairs; the **disjoint-support**
   condition (froth has disjoint node support; a coalescent pair shares nodes) is
   what makes the merge correct.

## Discipline notes

- Every fix carried an **independent oracle** (wolframscript ℘ companion for
  x0p0; wolframscript chord-winding + sqrt-continuation for 61um; closed-form
  lattice for fzse) — verified, not copied (Rule 5).
- Every load-bearing assertion was **mutation-proven** (perturb → RED → restore
  byte-clean; `git diff src/` empty each time).
- `winding_delta` and `extract_poles` are now both better-documented primitives;
  the disproven "loses a revolution" / "interstitial FP" framings were corrected
  in their docstrings + the corpus headers (Law 2).

## Follow-ons filed

- 61um: `cut-crossing-vs-continuation` semantics decision; walker-level
  auto-subdivision near branches (Fix-B deferred).
- fzse: double-pole multiplicity → co-designed into `padetaylor-90oh` (the
  support sets the self-merge tracks are its data); recall-is-walk-coverage
  follow-on (the bottom-edge misses are a path-network density matter).

## Status

Epic `padetaylor-l7yt` COMPLETE (6/6 + the x0p0 sibling). Suite GREEN across all
affected files in the campaign env; a full `Pkg.test()` on an idle box should now
read **N pass / 3 broken / 0 fail** (the 3 = v1ub + cm-n2 + pi2-tritronquee).
