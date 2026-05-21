# Worklog 056 — vector BVP solver + V8b P_I^(2) tritronquée figure

**Date**: 2026-05-21
**Author**: Claude Opus 4.7 (orchestrator) + serial Opus coding subagents
**Epic**: `padetaylor-0ln` (v0.2 — vector / multi-component Painlevé-type systems)
**Scope**: Orchestrated, serial (one Julia process at a time, Rule 7) build of a
vector boundary-value-problem solver (`VectorBVP`) and the v0.2 epic's final
figure — V8b, the P_I^(2) tritronquée vs Kapaev–Klein–Grava 2015. Closes the
last open blocker (`padetaylor-r2m`) and the V8b figure (`padetaylor-0ln.18`).

> **Take-home**: the v0.2 vector stack now has a BVP solver, and the P_I^(2)
> tritronquée figure ships. The figure was blocked for most of the session by
> what three agents in succession misdiagnosed as a "separatrix isolation"
> problem requiring research-grade Newton globalisation + continuation. It was
> not. It was a **one-line sign error** in `pI2_tritronquee_ic`: the IC
> returned the negative (non-solution) branch of the tritronquée on the
> negative real axis. The leading balance `u³ = -6x` gives `u > 0` for `x < 0`
> (`u³ = 6|x| > 0`); the IC returned `u < 0`, which has P_I^(2) ODE residual
> ~-9600 at `x = -20` — not a solution at all. Every prior agent then pinned
> the BVP to non-solution boundary data and got non-tritronquée results.

## Orchestration method

Parent Opus orchestrator; each unit of work delegated to one Opus coding
subagent (or, for read-only ground-truth gathering, parallel Sonnet recon
subagents — read-only fan-out is allowed in parallel; only Julia processes are
serialised, Rule 7). After each coding phase the orchestrator **re-verified the
new test file standalone** (Rule 3 skepticism — never trust the subagent's
GREEN claim), confirmed the commit, raised friction beads, then dispatched the
next phase. Beads created before code (Rule 8); the orchestrator closed each
bead only after independent re-verification.

## What shipped — `src/VectorBVP.jl`, ADR-0023, the V8b figure, 1 bug fix

| Phase | Bead | Deliverable | Commit | Tests |
|---|---|---|---|---|
| VB1 | 0ln.24 | `src/VectorBVP.jl` — first-order vector-system Chebyshev-Newton BVP solver (`y'=f(z,y)`, general linear BC `B_a y_a + B_b y_b = g`, autodiff Jacobian via `Taylor1` + optional override, barycentric vector callable) + **ADR-0023** | `0fd7104` | 307 |
| VB2 | 0ln.25 | Rule-6 LOC split (207→200 effective) + nonlinear cross-validation: companion-form scalar PI ≡ `vector_bvp_solve` vs the oracle-validated scalar `bvp_solve` | `3e1297c` | +118 (425) |
| VB3 | 0ln.26 | Wire `VectorBVP` into the `PadeTaylor` umbrella + `painleve_hierarchy_jacobian` analytic Jacobian helper + `PainleveHierarchyProblem` BVP convenience | `4f1e467` | +187 (wire-in) |
| V8b-A | 0ln.18 | `pI2_tritronquee_ic` lifted to `n_terms=2` (KKG `c_6=1` correction at `t=0`) | `51bc21e` | +43 (PH.1.7) |
| Sign fix | 0ln.31, 0ln.30 | `pI2_tritronquee_ic` wrong-branch sign bug fixed; `PH.1.3`+`PH.1.4` tautological ODE-residual checks replaced with genuine ones | `bb1e208` | 159 PH GREEN |
| V8b | 0ln.18 | `figures/kkg_pi2_tritronquee_pole_field.jl` — KKG 2015 Figs 7.4/7.5: vector-BVP tritronquée on the negative real axis + 21-pole field extracted by the vector path network | `d38628e` | 27 (PI2F.*) |

**Full-suite verification**: `julia --project=. -e 'using Pkg; Pkg.test()'` →
**4472 / 4472 GREEN** (12m47s). The umbrella wire-in, the `pI2_tritronquee_ic`
sign fix, and the three new test files (`vector_bvp_test.jl`,
`vector_bvp_wirein_test.jl`, `kkg_pi2_figure_test.jl`) introduce no regression.

`VectorBVP` design (ADR-0023): a first-order vector-system collocation solver —
P_I^(2) plugs in as the `d=4` companion system the v0.2 stack already produces.
The autodiff Jacobian (one `Taylor1`-jet evaluation of the RHS per component,
no new dependency) removes a whole class of hand-derivative drift (Rule 2);
the analytic `painleve_hierarchy_jacobian` doubles as its cross-check oracle.

## The V8b misdiagnosis chain (Rule 2 / Rule 3 — the session's main lesson)

V8b's figure needs the P_I^(2) tritronquée, which is a separatrix — forward IVP
diverges (worklog 055 lesson 47), so a BVP is mandatory. The figure was blocked
through **three** agent passes, each diagnosis wrong, until the root cause
surfaced:

1. **First V8b build agent** — tried the 2+2-endpoint BVP, got a "wrong branch"
   (u positive where it pinned u negative); concluded the solver couldn't
   isolate the separatrix. Honestly refused to fabricate a figure (correct
   instinct) but the *diagnosis* was wrong.
2. **Separatrix probe** — confirmed `VectorBVP` is sound at d=4, but concluded
   the tritronquée is a genuine separatrix the 2+2 BVP cannot select, and that
   Armijo + continuation was the fix. Also wrong.
3. **Recipe-finding probe** — checked the dominant balance and plugged both
   sign branches into the actual ODE. Found: `u = -∛6|x|^{1/3}` has row-4
   residual ~-9600 at `x=-20`; `u = +∛6|x|^{1/3}` has ~0. The IC was returning
   a non-solution. With the correct sign, the **plain** `vector_bvp_solve`
   converges in 2–3 Newton iterations.

The orchestrator independently verified the arithmetic before acting:
`(-∛6|x|^{1/3})³ = -6|x|`, and for `x<0`, `6x = -6|x|`, so `u³+6x = -12|x| ≠ 0`;
`(+∛6|x|^{1/3})³ = +6|x| = -6x`, so `u³+6x = 0`. Decisive.

## Hard-won lessons

50. **The P_I^(2) tritronquée is positive on the negative real axis.** The
    leading balance `40(u³+6x)=0` gives `u³=-6x`; for `x<0`, `u³=6|x|>0`, so the
    real solution has `u>0`. KKG's `u ~ -∛6·x^{1/3}` uses the *real* cube root
    `x^{1/3}=-|x|^{1/3}` for `x<0`, so `u = -∛6·x^{1/3} = +∛6·|x|^{1/3} > 0`.
    `pI2_tritronquee_ic` had conflated `x^{1/3}` with `|x|^{1/3}`. The earlier
    `be9dc3c` "y_3 sign correction" was internally consistent with the *wrong*
    branch — it fixed a symptom inside a wrong frame. Always verify an IC by
    plugging it into the equation, not by checking it differentiates
    consistently with itself.
51. **A companion-form ODE residual is a tautology if you use the RHS.** For a
    companion system the highest-derivative row `f[d]` is *defined* as the
    negation of the rest of the equation, so `‖u^(d) - f[d]‖` is identically
    zero for any state. `PH.1.3` and `PH.1.4` both did this and asserted
    nothing (Rule 5 violation). A genuine residual check needs `u^(d)` derived
    independently — from a closed form or a finite difference of `u^(d-1)`.
52. **"Investigate before build" caught a 1-line fix masquerading as research.**
    The user-directed investigation pass turned a planned multi-phase
    "separatrix solver" build into a one-line sign correction. Cheap probes
    before expensive builds.
53. **`vector_bvp_solve` is sound at d=4 / 2+2 split** — proven by the probe's
    linear 4th-order sanity test (solErr ~1e-15, BC residual exactly 0.0).
    VB1's tests had only exercised d=2; the probe closed that gap.
54. **The tritronquée pole-field march needs a fine wedge step.** The vector
    path network marches from a BVP-anchored seed into the `arg x≈0` pole
    wedge; the FW-default `h=0.3` overshoots a pole onto a blown-up state.
    `h=0.1` completes (389 nodes, 21 poles). Hand-tuned — the principled fix is
    the shared-Q root-distance wedge criterion (`padetaylor-0ln.23`).

## Beads — this session

**Created**: `0ln.24`–`0ln.26` (VB1–VB3), `0ln.27` (ChebUtil de-duplication —
`_chebyshev_D1` copied into `VectorBVP.jl`), `0ln.28` (V8b+ — the full KKG
Re/Im surface, deferred behind Stage-2 fill `0ln.20`), `0ln.29` (separatrix
isolation — created then **closed as a misdiagnosis**), `0ln.30` (PH tautology),
`0ln.31` (the sign-fix bug).
**Closed**: `0ln.24`, `0ln.25`, `0ln.26`, `0ln.18` (V8b), `0ln.29`
(misdiagnosis), `0ln.30`, `0ln.31`, and **`padetaylor-r2m`** (the v0.2 BVP-gap
blocker — resolved by VB1–VB3).
**Open / deferred**: `0ln.19` (V9 — v0.2 docs/release prep, the last planned
epic phase), `0ln.20`–`0ln.23` (V7 deferred sophistications), `0ln.27`
(ChebUtil dedup), `0ln.28` (V8b+ full surface).

## v1 corners (Rule 9 — documented, not hidden)

- `VectorBVP`: the `x>0` branch of `pI2_tritronquee_ic` is unimplemented
  (fail-fast guard); 3+1 / 4+0 BC splits are ill-conditioned (the separatrix
  growing mode) — the negative-axis 2+2 split is the supported path.
- V8b figure: the wedge-march step `h` is hand-tuned (`0ln.23`); the full KKG
  Fig 7.4/7.5 Re/Im *surface* over `[-20,20]²` is deferred — it needs vector
  path-network Stage-2 fine-grid fill (`0ln.28` → `0ln.20`). V8b ships the BVP
  curve + the pole-location scatter.

## Pickup point for the next agent

**V9 (`padetaylor-0ln.19`) — v0.2 docs + release prep — is the last planned
epic phase.** README/CHANGELOG/RESEARCH v0.2 sections, ADR review (ADR-0023 is
new), HANDOFF refresh, the v0.2 epic close-out. After V9 the epic
`padetaylor-0ln` can close. The orchestration pattern that worked: serial Opus
coding subagents, parallel Sonnet recon, one red-green-mutation bead each,
orchestrator re-verifies standalone + closes the bead + raises friction beads.

Probe artefact at `external/probes/pi2-tritronquee-bvp/` (committed — the
`external/probes/` tree is tracked): the recipe-finding scripts `step1`–`step12`
and `RECIPE.md`. Note `RECIPE.md`'s `tri_seed` `Δy₄` term has a `-24s`/`+24s`
sign slip — harmless to its 2+2 BVP, which only pins `u,u'`.
