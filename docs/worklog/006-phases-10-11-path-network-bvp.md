# Worklog 006 — Phases 10 + 11: PathNetwork + BVP; the agentic workflow that shipped two phases in one session

**Date**: 2026-05-13
**Author**: Claude Opus (orchestrating); multiple Sonnet subagents (literature, oracle, architecture, marker, octave)
**Scope**: Tier-2 path-network module (Phase 10) + Tier-3 Chebyshev-Newton BVP module (Phase 11), both GREEN + mutation-proven, plus full algorithmic-spec + figure-catalogue + ADR + recipe + Octave oracle scaffolding.  Three commits land in this session: `910aab9` (BVP ground truth + oracle), `0ada60f` (Phase 10 PathNetwork), `cc7d8ca` (Phase 11 BVP).

> Take-home: 284/284 tests passing (218 prior + 66 new across PathNetwork
> + BVP).  Senior-grade Julia: generic in `T <: AbstractFloat` with
> Float64 + BigFloat-256 + Complex{T} all exercised in tests; analytic
> Jacobians; fail-fast guards; step-norm Newton convergence; literate
> top-of-file docstrings.  Both modules pass mutation-proof.  Tier-3
> figures from `docs/figure_catalogue.md §6` are now algorithmically
> reachable pending the Dispatcher (Phase 12, bead `padetaylor-8lk`).

## Orchestration pattern — Opus codes, Sonnet does the boring work

The session's working agreement was set early by the user
("DO dispatch sonnet subagents to do summarisation, research work for you,
and to do grunt work like marker scripts.  You code, they do your boring work").
Saved as feedback memory `feedback-delegate-grunt-work`.  The pattern that emerged:

  - **Sonnet subagents (in parallel where independent)**: literature
    survey (5 papers FW2011 §3.1 + §3.2 + FW2014 + FW2015 + RF2014 +
    FFW2017), Julia BVP ecosystem survey, FFI BVP options survey,
    BVP-ref acquisition + pdftotext, marker_single PDF→markdown
    conversion, Octave oracle script (`external/probes/bvp-oracle/capture.m`,
    440 LOC), and architecture-design draft (ADR-0004 + DESIGN fragment).
    Total: 8 Sonnet subagents dispatched in this session.
  - **Opus (parent)**: figure catalogue synthesis from the 5 paper
    surveys, unified path-network spec, ADR review, code writing for
    `src/PathNetwork.jl` and `src/BVP.jl`, test code, mutation-proof
    procedures, commit composition, this worklog.

Two Sonnet subagents hit idle timeout (BVP-ref acquisition at ~12 min;
the marker subagent kept going past).  In both cases the in-flight
work landed on disk before the timeout (PDFs downloaded, oracle
generated, etc.), and the only loss was the final report message.
**Mitigation pattern that worked**: for long-running grunt tasks (marker
conversion of a 783KB book, Octave oracle with N=24 PI Newton), spawn
the heavy invocation with `run_in_background: true` inside the
subagent — the subagent itself stays responsive to its idle-timeout
while the work proceeds.  Documented in the Octave-oracle subagent's
prompt for future use.

Per CLAUDE.md Rule 7 (no parallel Julia agents): one `julia` process
at a time across all subagents.  Honoured by sequencing the Octave-
oracle subagent's final `julia --project=. -e 'include("test/_oracle_bvp.jl")'`
syntax check before Opus ran `Pkg.test()` on the PathNetwork module.
Read-only literature subagents ran in parallel safely; Octave
(`/usr/bin/octave`) did not contend with Julia precompile.

## What changed — Phase 10 (PathNetwork)

`src/PathNetwork.jl` (306 lines incl. ~70 docstring, ~150 effective LOC).
Public API: `path_network_solve(prob, grid; kwargs...) -> PathNetworkSolution`
with kwargs `h, order, wedge_angles, step_selection, step_size_policy,
max_steps_per_target, rng_seed`.

Algorithm verbatim from FW 2011 §3.1 (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-166`):

  - **Stage 1**: for each target in shuffled `grid`, find nearest
    visited node, step toward target with 5-direction wedge until
    within `h`.  Default wedge `[-π/4, -π/8, 0, π/8, π/4]` (FW 2011
    ±22.5°/±45°); parametrisable to RF 2014's ±15°/±30°.
  - **Stage 2**: for each fine-grid node find nearest visited;
    evaluate visited's stored canonical Padé at `t = (z_f - z_v) / h`;
    NaN+NaN·im sentinel for uncovered (Rule 1).

**Design decision: canonical Padé per visited node.**  Each visited
node gets one Padé centered there with real `h`, in addition to the
wedge-direction Padé used for selection.  Cost: one extra Padé per
node (`Coefficients` + `RobustPade`).  Benefit: Stage 2 always
interpolates inside `|t| ≤ 1`, never extrapolates past `t = 1`.  The
alternative — store the wedge-selected Padé and extrapolate — would
make Stage 2's off-direction accuracy depend on which wedge direction
won at each node, an undesirable coupling.

**Signature widening**: `PadeStepper.pade_step_with_pade!(state, f,
order::Int, h::Real)` → `h::Number`.  Backwards compatible — existing
tests pass `Real h` with `Real`-typed state; new PathNetwork tests
pass `Complex h` (the wedge-direction displacement) with
`Complex`-typed state.  The internal `T(h)` coercion enforces type
matching at runtime.

`Random` moved from `[extras]` to `[deps]` in `Project.toml` (used by
PathNetwork for deterministic target shuffle via `MersenneTwister(rng_seed)`).

### Test plan + mutation-proof (PathNetwork)

4 testsets shipped (PN.1.1, PN.1.2, PN.2.1, PN.4.1).  PN.2.2 (FW Table
5.1 z=30 quantitative ≤1e-13) + PN.3.1 (:steepest_descent ≡ :min_u)
deferred to a follow-up commit per ADR-0004's v1/v2 split.  Filed as
bead `padetaylor-yt1` for the next session.

  - **Mutation A**: in `_select_candidate`, replace `argmin` →
    `argmax` for `:min_u` branch.  Steers TOWARD poles.  Bite: PN.1.2
    lines 79+83 (finite-|u| check + stay-far-from-pole check).
    Verified: 2 fails on PN.1.2.
  - **Mutation D**: invert Stage-2 coverage check `> h_v` → `< h_v`.
    Covered points NaN-out; uncovered get extrapolated.  Bite: PN.1.1
    lines 45-46+58-60, PN.1.2 lines 94-96, PN.2.1 line 109.
    Verified: 9 fails across PN.1.1 (5) + PN.1.2 (3) + PN.2.1 (1) —
    the NaN-sentinel logic is load-bearing across the entire Stage 2.

## What changed — Phase 11 (BVP)

`src/BVP.jl` (418 lines incl. ~280 docstring, 139 effective LOC).
Public API: `bvp_solve(f, ∂f_∂u, z_a, z_b, u_a, u_b; N, tol, maxiter,
initial_guess) -> BVPSolution`; callable `(sol)(z) -> (u, u')`.

Algorithm verbatim from `references/bvp_recipe.md` §§1-5:

  1. Chebyshev extrema nodes `t_j = cos(jπ/N)`, j=0..N (DMSUITE
     convention: t_0 = +1, t_N = -1).
  2. Affine map `[-1, 1] → [z_a, z_b]` via `z(t) = (z_a+z_b)/2 +
     (z_b-z_a)/2 · t`, with `t = -1 ↔ z_a` and `t = +1 ↔ z_b`.
  3. `D₁` from standard Trefethen SMIM / Weideman-Reddy formula;
     `D₂ = D₁ · D₁` (option (a), adequate for N ≤ 50).
  4. Newton on interior nodes 2..N (BC rows enforced by direct
     substitution).  Analytic Jacobian `J = D₂_int,int − ((z_b - z_a)²/4) ·
     diag(∂F/∂u(z_int, u_int))`.  Residual `R = (D₂ u)_int + bc_col −
     ((z_b - z_a)²/4) · F(z_int, u_int)`.
  5. Barycentric Berrut-Trefethen 2004 evaluation on Chebyshev-2
     weights `w_j = (-1)^j` with endpoint halving.

The `(z_b - z_a)²/4` affine-scale factor is **the** subtle place to
get the algorithm wrong.  FW 2011 line 190 says "the Jacobian is
trivial to calculate" without writing it.  Recipe §4 derives it once.
Mutation A (below) confirms its load-bearing status.

### Step-norm Newton convergence (the bug I hit on first cut)

My first cut used residual-norm convergence with default `tol =
100·eps(T) ≈ 2.2e-14`.  All BVP tests errored: Newton hit the
discrete-residual floor (~1e-12 for PI N=20, cond(D₂) ≈ N² · eps), and
the tolerance was structurally unachievable.

**Fix**: switch to step-norm convergence `‖Δu‖_∞ ≤ tol` with default
`eps(T)^(3/4)`.  This is the production-Newton standard (NLsolve.jl-
style).  Iteration stops when Newton stops making progress, not when
an irreducible spectral floor is hit.  The discrete residual is
recorded in the `BVPSolution.residual_inf` field for diagnostics but
is no longer the convergence test.

This matches the recipe §7 open-spec-gap entry that FW 2011 leaves
unspecified.  Lesson recorded in the `bvp_solve` docstring.

### Spec-drift catch by the Octave oracle subagent

The Group 4 PI BVP test was initially seeded with BCs from the
`problems-oracle/_oracle_problems.jl` (`u_at_0_5 = 4.0044…`).  But
that oracle is for `u'' = 6u²` (equianharmonic Weierstrass ℘), NOT
Painlevé I (`u'' = 6u² + z`).  The Octave subagent caught the spec
drift on its own, recomputed correct PI BCs via `mpmath.odefun` at
40 dps, and cross-validated `dup_at_zb = 16.2045…` against
`mpmath.up(0.5) = 16.2045…` to ~1e-10 absolute.  Documented inline
in `test/_oracle_bvp.jl` Group 4 header + in `references/bvp_recipe.md`.

This is the worklog-002 / worklog-005 spec-drift pattern recurring
**at the oracle layer**, not the impl layer.  Lesson: when one
oracle exists in the project, agents may grab it without
double-checking the underlying ODE.  Mitigation: oracle files should
state the ODE prominently at the top.

### Test plan + mutation-proof (BVP)

7 testsets (BV.1.1, BV.1.2, BV.1.3, BV.1.4, BV.2.1, BV.3.1, BV.4.1,
BV.5.1) covering D₂ primitives, linear BVP, two nonlinear PI BVP
problems, barycentric eval, callable + endpoint derivatives,
fail-fast guards, BigFloat-256 sanity.

  - **Mutation A**: drop the (1/4) affine-scale factor.  Bite: 19 fails
    across BV.1.2 (8/8), BV.1.3 (1/4), BV.1.4 (2/3), BV.2.1 (8/9),
    BV.3.1 (2/4), BV.5.1 (1/2) — universally load-bearing.
  - **Mutation B**: swap u_a/u_b in BC enforcement.  Bite: 10 fails on
    the asymmetric tests (BV.1.3, BV.1.4, BV.2.1, BV.3.1); symmetric
    linear u(±1)=1 case correctly not bitten.
  - **Mutation C**: drop the `(-1)^j` alternating sign in barycentric
    weights.  Bite: 8 fails on BV.2.1 (Newton unaffected; callable's
    interpolation breaks).

## Frictions surfaced

  - **F1. Tool path discovery**: `marker_single` lives in per-project
    `~/Projects/<project>/.venv/bin/`, not on system PATH.  User
    correction halfway through the BVP-acquisition subagent: "loads of
    tools incl marker are in venv and .venv subfolders of ~/Projects".
    Saved as user memory `user-python-tools-venv-pattern`.
  - **F2. TIB VPN access for paywalled refs**: user activated the
    institution VPN mid-session to fetch SIAM Review (Berrut-Trefethen
    2004) and ACM TOMS (Weideman-Reddy 2000), both paywalled.  Saved
    as reference memory `reference-tib-vpn-paywalls`.
  - **F3. Subagent idle timeouts on long downloads + marker**: as
    noted in §"Orchestration pattern" above.  `run_in_background: true`
    inside the subagent for the heavy invocations avoided re-running.
  - **F4. CLAUDE.md Rule 8 vs TaskCreate**: user override
    ("you may use TaskCreate, claude.md is too strong on this").
    Saved as feedback memory `feedback-taskcreate-allowed`.  Rule 8
    forbids `TodoWrite/TaskCreate` for cross-session work tracking;
    `bd` remains the only persistent tracker, but in-session
    TaskCreate is now permitted.

## Pointers

  - [`docs/figure_catalogue.md`](../figure_catalogue.md) — 79 figures
    tiered T0-T5; this session unblocks Tier 2 (path-network) and is
    one Dispatcher away from unblocking Tier 3.
  - [`docs/unified_path_network_spec.md`](../unified_path_network_spec.md)
    — 14-section algorithm spec with 8 open-spec-gap decisions.
  - [`docs/adr/0004-path-network-architecture.md`](../adr/0004-path-network-architecture.md)
    — ADR for Phase 10 path-network.
  - [`docs/design_section_6_path_network.md`](../design_section_6_path_network.md)
    — Phases 10-14 table for `DESIGN.md §6`.
  - [`references/bvp_recipe.md`](../../references/bvp_recipe.md) —
    canonical Chebyshev-Newton recipe in 8 sections.
  - [`src/PathNetwork.jl`](../../src/PathNetwork.jl),
    [`src/BVP.jl`](../../src/BVP.jl) — the modules.
  - [`test/pathnetwork_test.jl`](../../test/pathnetwork_test.jl),
    [`test/bvp_test.jl`](../../test/bvp_test.jl) — tests +
    mutation-proof commentary.
  - [`external/probes/bvp-oracle/`](../../external/probes/bvp-oracle/)
    — Octave oracle (DMSUITE-based) + `oracles.txt`.
  - [`test/_oracle_bvp.jl`](../../test/_oracle_bvp.jl) — 47 pinned
    constants from the Octave run.
  - [`references/markdown/{BerrutTrefethen2004,TrefethenSMIM_2000_book,WeidemanReddy2000_DMSUITE_ACMTOMS26}/...md`](../../references/markdown/)
    — marker-converted primary sources.

## Bead state

Closed in this session:
  - `padetaylor-804` — BVP module (Phase 11 GREEN at commit cc7d8ca).

Open + ready (filed in this session):
  - `padetaylor-8lk` — **P0** Phase 12: Dispatcher.  The next P0 work
    item; composes PathNetwork + BVP via the EdgeDetector.
  - `padetaylor-bvh` — **P2** Phase 13: CoordTransforms (PIII/PV exp
    coordinates).  Tier-4 deferral.
  - `padetaylor-grc` — **P2** Phase 14: SheetTracker (PVI multi-sheet).
    Tier-5 deferral.
  - `padetaylor-yt1` — **P1** PathNetwork PN.2.2 + PN.3.1 follow-up
    (FW Table 5.1 quantitative + :steepest_descent test).

Still open from prior sessions:
  - `padetaylor-1jf` — PathNetwork v1 shipped at commit 0ada60f but
    bead remains open until PN.2.2 + PN.3.1 (`padetaylor-yt1`) land.
  - `padetaylor-8cr` — v2 P0 umbrella for FW Table 5.1 long-range.
  - `padetaylor-rgp` — figure-acceptance catalogue (shipped as
    `docs/figure_catalogue.md` but bead stays open as the tracking
    layer until Tier-2-through-5 figures actually reproduce in tests).
  - `padetaylor-c2p` — Edge detector.  Tier-3 prerequisite for the
    Dispatcher.
  - `padetaylor-kvi`, `padetaylor-jhq`, `padetaylor-2vz`,
    `padetaylor-61j`, `padetaylor-8pi` — older items per prior worklogs.
