# Worklog 007 — Phase 12: Dispatcher v1 (1D ordered-chain IVP↔BVP composition)

**Date**: 2026-05-13
**Author**: Claude Opus (orchestrating); one Sonnet subagent (FW 2011 algorithmic dive)
**Scope**: Tier-3 Dispatcher module shipped on v1 acceptance — 1D
ordered-chain composition of `PathNetwork` (IVP) and `BVP` (Chebyshev-
Newton) segments per FW 2011 §4.4 + line 192's derivative-match
diagnostic.  329/329 tests GREEN (284 prior + 45 new across 9
dispatcher testsets).  v2 (2D lattice + automatic edge detection) filed
as bead `padetaylor-k31`.

> Take-home: the dispatcher is intentionally **algorithmically thin** —
> no new numerics, just orchestration discipline.  The contribution is
> the public-API shape + the junction-diagnostic contract.  The
> Sonnet algorithmic dive surfaced the critical scope insight:
> FW 2011 line 190 says "161 separate BVP solutions; one for each grid
> line" — they NEVER solve a 2D elliptic BVP.  The dispatcher composes
> 1D segments only.  This let me scope v1 honestly without committing
> to machinery that's not even in the source paper.

## Orchestration pattern — Sonnet algorithmic dive, Opus codes the rest

The 2026-05-13 sessions established Opus-as-coder, Sonnet-as-grunt
(memory `feedback-delegate-grunt-work`).  Phase 12 used a single
parallel Sonnet subagent:

  - **Sonnet (read-only Explore)**: deep read of FW 2011 §§3.1, 3.2,
    4.4, 5.1.2, 5.4.1, 5.5.  Reported (a) algorithmic specification of
    the 5-point Laplacian classifier, (b) the BVP solver's 1D-per-row
    structure (the killer scope insight), (c) FW2011...md:192's
    1e-7 / 1e-8 derivative-match tolerance, (d) line 401's "we did it
    manually" admission, (e) parameter choices (h=0.5, order=30).
    400-word return; Opus consumed in one pass.

  - **Opus (parent, sole julia process)**: scope decision, ADR-cross-ref,
    test plan, RED skeleton, GREEN impl, three mutations, this worklog,
    commit composition.

Per CLAUDE.md Rule 7 (no parallel Julia agents): the Sonnet subagent
was Explore-only (no julia invocation).  Opus held the single Julia
process throughout.

## v1 scope decision — and why FW 2011 §4.4 is "the algorithm"

The bead's original framing was "covering an arbitrary 2D complex-plane
domain" with per-grid-point region tagging.  The Sonnet dive surfaced
two facts that re-scoped v1:

  1. **FW 2011 never solved a 2D BVP**.  Line 190: "161 separate BVP
     solutions; one for each grid line."  The 2D lattice is filled by
     **independent 1D BVPs along rows or columns**, each composed with
     IVP path-network solutions at the row endpoints.
  2. **FW 2011 classified manually in their published runs**.  Line 401
     verbatim: "we implemented this 'manually' after having inspected,
     in preliminary test calculations, where the pole field edges
     appeared."  The 5-point Laplacian classifier of §3.2.2 is offered
     as a diagnostic that "can readily be fully automated," but FW
     did not automate it.

Together: the algorithmic essence is the **1D IVP↔BVP junction with
derivative-match diagnostic** (line 192 — "If the agreement is not to
within some tolerance, typically set to 10⁻⁷ or 10⁻⁸, an increase in N
is indicated").  The 2D lattice and the automatic edge detector are
*productivity layers*, not new algorithm.

**v1 ships the algorithmic essence**; v2 (bead `padetaylor-k31`) wraps
it for 2D lattice production with the 5-point Laplacian classifier
(formerly bead `padetaylor-c2p`, now consumed by v2).  This split is
faithful to FW and avoids over-committing to undocumented inference
about what the 2D dispatcher should be.

## What changed — `src/Dispatcher.jl` (322 lines, ≈150 effective LOC)

Public API:

  - `IVPSegment(z_end::Complex{T}, dense_grid::Vector{Complex{T}})` —
    one IVP segment of the chain.  Starts at the previous segment's
    terminus (or `prob_ivp.y0` for the first segment), walks the
    FW 2011 §3.1 5-direction wedge path-network to within `h` of
    `z_end`, evaluates Stage-2 at `dense_grid`.
  - `BVPSegment(z_end, u_b, dense_grid; initial_guess=nothing)` — BVP
    segment.  Left BC `u_a` is set by the dispatcher from the previous
    terminus; user supplies `u_b`.  Initial guess defaults to the
    linear ramp.
  - `DispatcherSolution{T}` — composed result.  Carries per-segment
    `:ivp` / `:bvp` kind, the underlying `PathNetworkSolution`s and
    `BVPSolution`s, junction coordinates + match tuples, concatenated
    dense output with `:ivp` / `:bvp` per-cell region tagging.
  - `dispatch_solve(prob_ivp, f_bvp, ∂f_∂u_bvp, segments; kwargs...)`
    — the public driver.  `derivative_match_tol = 1e-7` per FW line
    192; `derivative_match_strict = false` by default (diagnostic
    only, FW-faithful), `true` for hard-throw regression tests.

The signature **takes both the 2nd-order IVP `f` (via `prob_ivp.f`) and
the 1st-order BVP `f_bvp` + `∂f_∂u_bvp`** as separate arguments.  This
is honest about the algorithmic constraint that BVP applies to
`u'' = F(z, u)` (no `u'` dependence), while IVP path-network handles
general `u'' = f(z, u, u')`.  For autonomous-in-`u'` ODEs (the FW
target class), the user adapts trivially: `f_bvp(z, u) = f_ivp(z, u, _)`
ignoring the third arg.

Internal walk: `cur_(z, u, up)` accumulator threads state across
segments.  At an IVP segment k ≥ 1: build a fresh `PadeTaylorProblem`
with IC at `(cur_z, cur_u, cur_up)`, walk to `seg.z_end` plus the
user's `dense_grid`, take the terminus as the new accumulator.  At a
BVP segment: solve with `u_a = cur_u, u_b = seg.u_b`, eval barycentric
at `dense_grid` plus `seg.z_end`, propagate the BVP-recovered `(u, u')`
at the right endpoint as the new accumulator.

## The junction-match contract (FW 2011 line 192 verbatim)

At every interior junction between segments k and k+1:
  - **For X→BVP transitions** (where BVP starts at z = cur_z): record
    `(Δu, Δu')` where `Δu = |cur_u - u_BVP(cur_z)|` (zero by BC
    construction; mutation-proof target) and `Δu' = |cur_up -
    u'_BVP(cur_z)|` (the **diagnostic** — barycentric derivative of
    BVP at left endpoint vs IVP's analytic derivative there).
  - **For X→IVP transitions** (where IVP starts at z = cur_z with IC
    = `(cur_u, cur_up)`): record `(0, 0)` since the IC propagation is
    by construction exact.

The 1e-7 tolerance is a **diagnostic threshold**, not a convergence
test — FW says "an increase in N is indicated" if exceeded, leaving
the call to the user.  Default behaviour: record + don't throw.
`derivative_match_strict = true` opt-in for hard throw on violation;
used in regression tests where stability of the diagnostic IS the
contract.

## Test plan + mutation-proof

`test/dispatcher_test.jl` ships 9 testsets (5 outer, with DP.4.1
sub-fail-fast tests counting as 5 sub-tests).  45 new assertions:

  - **DP.1.1**: single `IVPSegment` ≡ direct `path_network_solve`.
    Equianharmonic ℘ ODE (Phase-10 test bed).  Passes through to the
    underlying path-network without alteration.
  - **DP.1.2**: single `BVPSegment` ≡ direct `bvp_solve`.  Linear ODE
    `u'' = u` on [-1, 1] (Phase-11 BV.1.2 reference setup).
  - **DP.2.1**: two-segment IVP→BVP with junction at z=0.3.  Linear
    ODE, IC `u(0)=1, u'(0)=0 → u(t)=cosh(t)`.  Δu ≤ 1e-12 (BC), Δu' ≤
    1e-7 (FW threshold).  Closed-form cross-check at sample points.
  - **DP.3.1**: three-segment IVP→BVP→IVP (the FW 2011 §4.4 pattern).
    Two interior junctions; cosh(0.75) within 1e-5 at the third
    segment's grid point.
  - **DP.4.1**: 5 fail-fast guards — empty segments, negative tol,
    negative h, N_bvp < 4, strict-mode tol violation.

Mutation procedure verified 2026-05-13:

  - **Mutation A** (drop u' in BVP→IVP propagation, set `zero(CT)`):
    bites DP.3.1 line 163 only — the downstream IVP starts with u'=0
    instead of sinh(0.6)≈0.6367, u(0.75)≈1.199 vs cosh(0.75)≈1.295.
    Isolates the BVP→IVP transition.
  - **Mutation B** (swap u_a/u_b in the BVP solve): bites 5 tests —
    DP.2.1 Δu and Δu' (the BC mismatch propagates), DP.3.1 first
    junction's Δu/Δu', and DP.3.1's cosh(0.75) downstream.  This is
    the Phase-11 BVP Mutation-B analogue at the dispatch layer; the
    `cur_u` ↔ `seg.u_b` distinction is load-bearing across every
    IVP→BVP-coupled test.
  - **Mutation C** (drop dense_grid push for segments k ≥ 2): bites
    DP.3.1 — `idx_075 = findfirst(...)` evaluates to `nothing`,
    downstream indexing operations error.  Confirms the
    grid-concatenation invariant.

## Frictions surfaced

  - **F1. PadeTaylorProblem zspan re-construction**: each IVP segment
    constructs a *fresh* `PadeTaylorProblem` with `zspan = (cur_z,
    seg.z_end)` and `y0 = (cur_u, cur_up)`, overriding any
    `prob_ivp.zspan[2]` or `prob_ivp.y0` the user passed in for
    segments beyond the first.  This is correct (each segment carries
    its own state) but surprising — documented in the dispatch_solve
    docstring.  Caught at design time, not impl.
  - **F2. The "no parallel Julia" rule + Sonnet subagent contention**:
    only one Sonnet subagent (Explore type, read-only) was dispatched
    for Phase 12.  Compare worklog 006 which dispatched 8 Sonnet
    subagents in parallel; Phase 12's algorithmic scope is narrower and
    the FW 2011 dive answered all open questions in one report.
  - **F3. The DP.4.1 strict-mode test**: chose N_bvp=4 +
    derivative_match_tol=1e-12 to deliberately exceed the spectral
    floor.  Empirically verified the strict-mode `ErrorException` does
    fire; could in principle be brittle if BVP at N=4 on a short
    interval happens to achieve <1e-12 by accident.  Recorded for
    future tightening if needed.

## Pointers

  - [`src/Dispatcher.jl`](../../src/Dispatcher.jl) — the module.
  - [`test/dispatcher_test.jl`](../../test/dispatcher_test.jl) — tests
    + mutation-proof commentary at the bottom.
  - [`docs/adr/0004-path-network-architecture.md`](../adr/0004-path-network-architecture.md)
    — ADR governing Tier-3 architecture (deferral plan now consumed by
    Phase 12 v1).
  - [`docs/unified_path_network_spec.md`](../unified_path_network_spec.md)
    §8 — dispatcher algorithm synthesis (now realised in v1's 1D form).
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    lines 155-208 (§3.1 + §3.2), 249-261 (§4.4 composition), 192
    (derivative-match tol), 401 (manual classification admission).

## Bead state

Closed in this session:
  - `padetaylor-8lk` — Phase 12 Dispatcher v1 (GREEN at this commit).

Filed in this session:
  - `padetaylor-k31` — **P1** Phase 12 v2: 2D lattice dispatcher with
    automatic edge detection.  Consumes the 5-point Laplacian
    classifier formerly tracked as `padetaylor-c2p`.

Still open from prior sessions (no change):
  - `padetaylor-1jf`, `padetaylor-8cr`, `padetaylor-yt1`, `padetaylor-rgp`,
    `padetaylor-c2p` (now redundant with `padetaylor-k31`; consider
    superseding), `padetaylor-kvi`, `padetaylor-jhq`, `padetaylor-2vz`,
    `padetaylor-bvh`, `padetaylor-grc`, `padetaylor-61j`, `padetaylor-8pi`.

## Hard-won lesson (for the lesson list in HANDOFF.md §"Hard-won")

**Read the source paper for SCOPE, not just for algorithm**.  FW 2011
line 190 ("161 separate BVP solutions; one for each grid line") and
line 401 ("we implemented this 'manually'") together re-scoped Phase
12 v1 from "2D dispatcher with edge detector" to "1D ordered-chain
dispatcher".  Reading the algorithm without reading what the authors
*actually did* leads to over-engineering the v1 scope.  Sonnet's
algorithmic dive caught both lines in one report; my pre-dive plan had
the wrong v1 scope.  Future agents: when scoping a phase, ask the
algorithmic-dive subagent specifically "what did the authors do in
practice, vs what does the algorithm permit?"
