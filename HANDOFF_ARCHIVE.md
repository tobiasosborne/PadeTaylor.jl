# HANDOFF_ARCHIVE.md — pre-2026-05-13-part-3 session entries

> This file accumulates historical session entries that have rotated off
> the front of `HANDOFF.md` to keep it readable.  Sessions are in
> reverse chronological order (most-recent-archived at the top), matching
> the convention `HANDOFF.md` uses for its retained recent sessions.
>
> When you need historical context — what an old bead did, what an older
> hard-won lesson said, the state of the test suite three months ago —
> read this file in full.  Otherwise, prefer `HANDOFF.md` for the
> current state.

### Session 2026-05-13 (post-v0.1.0, part 3) — classical-Padé Toeplitz `\` finding (worklog 020)

User-driven investigation: "did FW use arb prec?"  Answer: no.  FW
2011 line 350 explicitly describes their Padé method —
**Toeplitz approach + MATLAB `\` (LU/QR), no SVD anywhere**.  GGT
2013 (which we port via `src/RobustPade.jl`) was published two
years after FW 2011 and replaces `\` with SVD for robustness
against Froissart doublets.  For well-conditioned cases (℘ on the
equianharmonic trajectory), SVD's robustness is unused and costs
both per-step accuracy and wall time.

Empirical wedge-walker probe (5-direction `:min_u`, h=0.5, order=30,
F64):

| target | method     | wall (s) | rel-err   | vs FW Table 5.1   |
|--------|------------|----------|-----------|-------------------|
| z=30   | :svd       | 5.84     | 6.6e-12   | 87× worse         |
| z=30   | :classical | 0.01     | 1.54e-13  | 2× worse          |
| z=10⁴  | :svd       | 17.76    | 6.05e-6   | 25,800× worse     |
| z=10⁴  | :classical | 3.45     | 6.15e-11  | **3.8× BETTER**   |

Per-step probe (10-step trace along 22.5° wedge, F64 vs BF-256
truth): step 3 (closest to the z=1 lattice pole) — SVD 6.2e-10,
classical 7.1e-13 — **870× per-step accuracy improvement**.

**Beads filed (3)**:

  - `padetaylor-txg` **P1** — *Ship classical-Padé as F64 default*.
    Add `classical_pade_diagonal(c, m)` to `src/RobustPade.jl` per
    `FW2011_*.md:346-350`; dispatch by element type (classical for
    `Float64`/`Float32`/their complex, SVD-with-Jacobi for
    `BigFloat`/`Arb`); fallback to SVD on singular Toeplitz.
    Tightens PN.2.2 rtol `1e-9 → ~1e-12` and PN.2.3 rtol
    `5e-5 → ~1e-10`.  Acceptance + ADR-0005 sketched in the bead
    description.
  - `padetaylor-7zw` **P2** — *BF-256 tritronquée pin: rerun FW
    Fig 4.1 step-(i) BVP at BigFloat-256*.  BVP is already generic
    in `T <: AbstractFloat` (BV.5.1 confirms); test/fw_fig_41_test.jl
    just needs a parallel BF-256 testset.
  - `padetaylor-u7o` **superseded by padetaylor-txg** — the four
    speculative attack vectors are obsolete; classical-Padé is the
    answer.

No source changes committed yet — the implementation work is the
next session, blocking on the new P1 bead.

### Session 2026-05-13 (post-v0.1.0, part 2) — z=10⁴ investigation + PathNetwork enhancements (worklog 019, commit 29de074)

Pivot from naive "rerun PN.2.2 at z=10⁴ BF-256" (~4.5 h wall) to a
Float64 investigation question driven by what FW actually claimed.

  - **PN.2.3 testset** (`test/pathnetwork_test.jl`): routine Float64
    z=10⁴ regression test, `rtol = 5e-5` vs `u_at_10000_FW_ref =
    21.02530339471055`, plus `|imag(u)| < 5e-5` and `length(visited_z)
    > 20_000` (step-count regression detector).  Mutation E
    (canonical-Padé invariant) re-verified — bites PN.2.2 + PN.2.3 +
    cascade through Phase-9 / LatticeDispatcher.
  - **PathNetwork.path_network_solve enhancements**: new opt-in
    kwargs `verbose::Bool = false` + `progress_every::Integer = 500`
    for eager-flushed Stage-1 progress lines; `_wedge_evaluations`
    now runs the 5 wedge candidates under `Threads.@threads` (per-
    thread `PadeStepperState`, deterministic argmin by index).
    User RHS `f` must be thread-safe (documented in docstring).
    Measured 1.32× speedup at 8 threads on full `Pkg.test()`.
  - **F64 sweep over `(order, h, z)`** documented in worklog 019.
    Finding: `order` saturates by 30 in F64 long-range (truncation
    below roundoff); raising `order` does NOT close the FW gap.
  - **`external/probes/pathnetwork-long-range/capture.jl`**: offline
    BF-256 sweep probe over `z ∈ {30, 100, 500, 1000, 10000}`,
    ~3.3 h wall at 8 threads, re-runnable.  Not committed with full
    `oracles.txt` (multi-hour run; placeholder in commit).
  - `Project.toml`: add `Printf` (stdlib) for verbose-mode
    `@sprintf` / `@printf` helpers.

Test suite: 1311 → **1314 GREEN** (+3 PN.2.3 assertions).  Wall
~1m45s at 8 threads.

### Session 2026-05-13 (post-v0.1.0, part 1) — Documenter site + v0.1.0 release (commits 30b3298, 38a49ae, tag v0.1.0)

Two P1 beads closed:

  - **`padetaylor-36w`** — Documenter.jl docs site at `docs/{Project.toml,
    make.jl, src/*.md}`.  Self-bootstrapping build (`julia
    --project=docs docs/make.jl`); local-only per CLAUDE.md Rule 11
    (no `deploydocs`, no CI).  Sections: Home (overview + Phase-6
    pole-bridge headline + 14-tier status table), Architecture
    (synthesis of all four ADRs), API (per-module `@autodocs`
    blocks for all 14 modules + extensions), Figures (curated
    tier-by-tier from `docs/figure_catalogue.md`).
  - **`padetaylor-8ll`** — v0.1.0 release: Project.toml bumped
    `0.1.0-dev → 0.1.0`; new `CHANGELOG.md` with Keep-a-Changelog
    v0.1.0 entry covering all 14 phases + headline empirical
    results + known limitations; README.md status table refreshed
    (Phases Z-6 → all 14 tiers; tests 218 → 1311); annotated
    `git tag v0.1.0`.  **Not pushed** per CLAUDE.md Rule 6.

Local branch: 7 commits ahead of `origin/main`, tag `v0.1.0`,
nothing pushed.

### Beads filed this session (parts 1-3)

  - `padetaylor-36w`, `padetaylor-8ll`, `padetaylor-g9x` — all closed.
  - `padetaylor-txg` **P1** (classical-Padé F64 default) — open.
  - `padetaylor-7zw` **P2** (BF-256 tritronquée pin) — open.
  - `padetaylor-u7o` — superseded by `padetaylor-txg`.

### Open beads end-of-session (post worklog 020)

One P1 + several P2/P3:

  - **`padetaylor-txg`** P1 — classical-Padé default (next session).
  - `padetaylor-7zw` P2 — BF-256 tritronquée pin.
  - `padetaylor-g9x` (closed), `padetaylor-bhw`, `padetaylor-bho`,
    `padetaylor-8ui`, `padetaylor-1a3`, `padetaylor-gky`,
    `padetaylor-8pi`, `padetaylor-61j` — see `bd ready -n 30`.

### Hard-won lessons added this session (worklogs 019-020)

**32. Read the paper's wall times to constrain its arithmetic
precision.**  FW Table 5.1's 26.5 s for z=10⁴ Padé rules out
`vpa` — MATLAB's arbitrary-precision arithmetic is ~100× slower
than `double`.  Their numbers MUST be `Float64`.  When a paper
claims a specific rel-err number, the only-arithmetic-that-fits-
the-wall-time is often diagnostic of the underlying impl choices.

**33. Per-step roundoff amplifies non-linearly near pole crossings.**
10-step path-network probe at h=0.5 along 22.5° wedge sees two
algebraically-equivalent `h^k` rescaling schemes diverge by
2.5e-10 after just 10 steps — most concentrated in the single step
crossing the `z = 1` lattice pole (~65× amplification).  Intrinsic
to the flow near singularities, not a bug.

**34. `order` saturates at roundoff long before truncation in
Float64 long-range.**  Order 30, 40, 50 give identical rel-err at
z=10⁴ in Float64.  Raising `order` cannot close the FW gap; the
gap is roundoff-limited.

**35. When porting from a paper, port the actual algorithm, not
the "modern equivalent".**  We adopted GGT 2013 Algorithm 2 (Robust
Padé via SVD) — published two years AFTER FW 2011.  FW used
classical Padé via Toeplitz `\`.  Algorithmic progress moved the
recipe in a direction (SVD robustness) that costs accuracy on the
well-conditioned case.  When reproducing a paper's published
numbers, port the paper's algorithm.

**36. "Modern equivalent" sometimes means "added work the paper
didn't need".**  GGT 2013's robustness is real (Froissart
doublets, near-singular blocks).  On Painlevé pole-fields the
robustness machinery is unused; the bidiagonalization +
iterative diagonalization + tolerance-based degree-reduction cost
580–1000× per-step accuracy AND 5-580× wall time relative to
classical `\` on the smooth case.  Robustness has a price; spend
it only where needed.

**37. The right precision dispatch is element-type-driven.**  At
`Float64` the dynamic range is narrow and Demmel-Kahan SVD (or
LU) suffices.  At `BigFloat` it opens up and relative-accuracy
Jacobi SVD is the load-bearing tool.  Dispatch by `T`; default to
the right method per tier; expose alternatives as opt-in.

**38. Classical Padé via Toeplitz `\` is both faster AND more
accurate than GGT 2013 SVD for well-conditioned F64 Padé tables.**
Probe (z=10⁴, h=0.5, order=30, F64): SVD 17.8 s / 6.05e-6 rel-err;
classical 3.5 s / 6.15e-11 rel-err.  Classical beats FW's published
2.34e-10 by 3.8×.  See worklog 020 for full empirical table +
dispatch design (`padetaylor-txg`).

## Previous handoff content (pre-v0.1.0)

### Session 2026-05-13 late evening (part 4) — `padetaylor-grc` closed

One GREEN commit lands this session: Phase 14 / Tier-5 SheetTracker
for the sixth Painlevé equation per FFW 2017 §2.2.

  - `padetaylor-grc` **closed** — new module `src/SheetTracker.jl`
    (~150 LOC) exports:
    - `pVI_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''`
      Factory for the ζ-plane PVI RHS (FFW 2017 eq. 3, md:144).
      Removes the `z = 0` branch point; leaves `z = 1` as a
      `ζ = 2π·i·k` lattice on the imag axis.
    - `winding_delta`, `accumulate_winding`, `sheet_index` —
      path-side primitives for computing Riemann-sheet indices
      after a regular PathNetwork walk.  Normalisation to `(-π, π]`;
      sheet index via `round(total / 2π)`.
    - **PVI coordinate conversion = PV's** (FFW2017...md:146);
      callers reuse `pV_z_to_ζ` / `pV_ζ_to_z`.  No duplicate helper.
  - **26 new assertions** across 7 testsets (ST.1.1-ST.1.7):
    hand-pinned RHS at degenerate + non-degenerate parameters,
    end-to-end PVI direct-vs-transformed agreement at one Padé step
    (≤1e-10), winding primitives on closed CCW/CW loops + non-
    enclosing paths, sheet-index conversion, branch-lattice
    magnitude check at `ζ = 2π·im`.
  - Mutation-proven: Mutations O (RHS `/2 → /3`) + P (drop winding
    normalisation) + Q (`round` → `floor`) all bit as predicted.
  - `docs/figure_catalogue.md §6` row T5 marked PARTIAL: ζ-plane +
    winding primitives shipped; η-plane PVI eq. (Fig 2 column 1
    only) and constrained-wedge PathNetwork routing remain.

Test suite: 1285 → **1311 GREEN** (+25 ST.1.* + 1 umbrella `isdefined`).
Wall ~2m00s.

### Beads filed this session (part 4)

None.

### Open beads end-of-session (worklog 018)

Two P2 beads remain — both admin/non-code:
  - `padetaylor-61j` (Willers 1974 paper acquisition),
  - `padetaylor-8pi` (GLA piracy friction record).

No P0 / P1 / P2-code beads open.  All Tier-1 through Tier-5
architectural deliverables are at least PARTIAL acceptance.

### Hard-won lessons added this session (worklog 018)

29. **"Branch point in Float64" is a magnitude check, not isfinite**.
    At `ζ = 2π·im` exactly, `exp(ζ)` returns `1 + O(eps)·im` (the
    rounded `2π` constant is non-exact), so `(e^ζ - 1) = O(eps)` and
    the RHS blows up to `~1e30` — huge but still finite.  Test
    magnitude (`> 1e10`), not `isfinite`.  Downstream Padé-Taylor
    finite-Taylor-coefficient check is the actual fail-loud point.

30. **`atan2` angle differences need `(-π, π]` normalisation for
    winding**.  Raw `angle(z_new) - angle(z_old)` has `±2π` jumps
    at the branch cut.  Normalising per step to `(-π, π]` gives
    the "shortest signed path" interpretation, valid iff each path
    step is < π in angular extent.  Callers must walk with `h <<
    distance-to-branch` (worked out in the module docstring).

31. **Topological direction ≠ algorithm's same-side-wrap convention**.
    A step from `(-1, +0.01)` to `(-1, -0.01)` is geometrically
    downward (clockwise around origin) but normalised winding-delta
    reads `+0.02` (counterclockwise small).  The algorithm interprets
    every sub-π step as "going the short way around"; only actual
    circumambulation (loops > π built up over multiple steps)
    registers sheet changes.  Don't pin direction on individual
    cross-cut steps in tests.

## Last commit before this handoff (previous session)

CoordTransforms PIII/PV GREEN (bead `padetaylor-bvh` closed) — worklog 017.

### Session 2026-05-13 late evening (part 3) — `padetaylor-bvh` closed

One GREEN commit lands this session: Phase 13 / Tier-4 exponential
coordinate transforms for PIII and PV per FFW 2017 §2.1.

  - `padetaylor-bvh` **closed** — new module `src/CoordTransforms.jl`
    (~100 LOC) exports:
    - `pIII_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''` and
      same for PV — RHS factory closures ready for `PadeTaylorProblem`.
    - `pIII_z_to_ζ` / `pIII_ζ_to_z` and `pV_z_to_ζ` / `pV_ζ_to_z` —
      forward + inverse IC conversions at a single point.
    Module-header chapter cites FFW2017...md:39-48 verbatim (the two
    transformed equations P̃_III and P̃_V) and documents the
    deferrals (non-uniform Stage-1 nodes, adaptive Padé `h`, multi-
    sheet recovery).
  - **23 new assertions** across 6 testsets (CT.1.1-CT.1.6):
    IC round-trips, hand-pinned symbolic RHS values, end-to-end
    agreement direct-vs-transformed at one Padé step (≤1e-10).
  - Mutation-proven: Mutations L (PIII RHS `/4 → /3`) + M (PIII IC
    sign flip) + N (PV RHS `(w+1)/(w-1)` swap) all bit as predicted.
  - `docs/figure_catalogue.md §6` row T4 marked PARTIAL — ExpCoords
    shipped; non-uniform nodes + adaptive `h` are independent
    follow-ups, not filed as beads.

Test suite: 1262 → **1285 GREEN** (+22 CT.1.* + 1 umbrella `isdefined`).
Wall ~1m56s.

### Beads filed this session (part 3)

None.

### Open beads end-of-session (worklog 017)

Three P2 beads remain — same set as after worklog 016, minus `padetaylor-bvh`:
  - `padetaylor-grc` (Phase 14 SheetTracker — PVI, Tier-5),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction).

No P0 or P1 beads open.

### Hard-won lessons added this session (worklog 017)

26. **Hand-derive IC-conversion factors on paper before coding**.
    PIII's `w' = (z u + z² u') / 2` has a load-bearing `z²` (not `z`)
    on the `u'` term, from the `dz/dζ = z/2` factor in the chain rule.
    Pin the result in a test with a hand-computed value (CT.1.1's
    `5/8`); mutation-test the sign + factor.

27. **End-to-end agreement is a stronger spec than a symbolic RHS pin**.
    Symbolic pins (CT.1.2, CT.1.5) catch transcription errors with
    crisp messages; direct-vs-transformed agreement (CT.1.3, CT.1.6)
    catches algebra errors that compensate at one sample point.  Ship
    both.

28. **Helpers-only module beats premature driver abstraction**.  Phase
    13's scope could have been a full `pIII_pV_solve(...)` driver.
    Choosing helpers-only (RHS factory + IC conversion) lets callers
    compose with `PadeTaylorProblem` + `path_network_solve` naturally,
    ~`100 LOC` source vs ~`250 LOC` for a driver wrapper, same
    downstream usability.

## Last commit before this handoff (previous session)

FW 2011 Fig 4.1 step-(i) BVP quantitative pin GREEN (bead `padetaylor-0c3` closed) — worklog 016.

### Session 2026-05-13 late evening (part 2) — `padetaylor-0c3` closed

One GREEN commit lands this session: reproduce the tritronquée FW
eq. 4.1 reference values via the FW Fig 4.1 step-(i) recipe.

  - `padetaylor-0c3` **closed** — `test/fw_fig_41_test.jl` (new ~110
    LOC) demonstrates the canonical FW Fig 4.1 step (i): Chebyshev-
    Newton BVP for `u'' = 6u² + z` on `[-20i, +20i]` with leading-term
    `u(z) = -√(-z/6)` Dirichlet BCs.  At N=240, `u(0)` pins to ≤3.5e-13
    and `u'(0)` to ≤5.3e-11 vs FW eq. 4.1 — both well under the bead's
    1e-10 spec.  6 testsets / 12 assertions; structural mutation-self-
    proof at FF.1.6 (the `+√` branch diverges).
  - **No `src/` changes**.  The recipe is the test itself; encapsulating
    it in a `fw_fig_41_axis_bvp(z_max; ...)` wrapper would be premature
    abstraction (no downstream caller, the bead is a pin not a feature).
  - `docs/figure_catalogue.md` row 4.1 marked PARTIAL with the step-(i)
    pin shipped + deferral notes for steps (ii) and (iii).

Test suite: 1250 → **1262 GREEN** (+12 FF.1.*).  Wall ~1m49s.

### Beads filed this session (part 2)

None.

### Open beads end-of-session (worklog 016)

Four P2 beads remain — same set as after worklog 015, minus
`padetaylor-0c3`:
  - `padetaylor-bvh` (Phase 13 CoordTransforms),
  - `padetaylor-grc` (Phase 14 SheetTracker),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction).

No P0 or P1 beads open.

### Hard-won lessons added this session (worklog 016)

23. **The leading-term Dirichlet BC pins `u` exactly but not `u'`**.
    For BVP step-(i) recipes that derive `u'(endpoint)` from the BVP
    (no Neumann BC), the converged `u'` differs from the analytical
    leading-term derivative by an o(1) correction.  At `z = ±20i` for
    tritronquée: gap ≈ 2.5e-4.  Use a looser tolerance documenting
    the gap, or fetch an independent oracle for a tighter pin.

24. **N-tuning probe before pinning rtol on a fresh BVP**.  Different
    segments need different N.  `[-18, -14]` (Fig 4.6): N=20 suffices.
    `[-20i, +20i]` (Fig 4.1): N=240 for ≤1e-10.  Probe spectral
    convergence empirically (each +20 in N typically buys ~10× error
    reduction) BEFORE committing test's N parameter.

25. **`residual_inf` is not a convergence metric for spectral BVPs**.
    The BVP's step-norm Newton (worklog 006 lesson 8) converges in
    `‖Δu‖_∞ ≤ eps^(3/4)` while `‖R‖_∞` floats at `N²·eps(T)` or higher.
    Gating tests on `residual_inf ≤ TINY` is structurally wrong —
    Newton CAN converge with `‖R‖_∞ ~ 1e-8` at N=200 and still deliver
    `u(0)` to 1e-13.  Gate on per-eval error against the oracle, not
    on the residual.

## Last commit before this handoff (previous session)

PathNetwork `enforce_real_axis_symmetry` kwarg GREEN (bead `padetaylor-dtj` closed) — worklog 015.

### Session 2026-05-13 late evening — `padetaylor-dtj` closed

One GREEN commit lands this session: the library-level cure for the
y-asymmetry bug surfaced in worklog 014.

  - `padetaylor-dtj` **closed** — `PathNetwork.path_network_solve` gains
    opt-in `enforce_real_axis_symmetry::Bool` kwarg.  When `true`, walks
    only upper-half + on-axis canonical representatives (collapses ±0.0im
    signed zeros via `complex(real(z), abs(imag(z)))`) and mirrors lower-
    half outputs via `conj`.  Bit-exact `u(z̄) = ū(z)` invariant; roughly
    2× faster than full-plane walk.  Validates IC-on-real-axis precondition
    (throws `ArgumentError` if `Im(z₀) > tol`, `Im(u₀) > tol`, or
    `Im(u'₀) > tol`); trusts caller on real-coefficient `f`.  Default
    `false` preserves FW 2011 algorithm verbatim.
  - **4 new testsets / 13 assertions** PN.6.1-PN.6.4, mutation-proven
    (Mutations G + I bite as predicted; see worklog 015).
  - `examples/tritronquee_3d.jl` simplified — drops the user-space upper-
    half walk + manual mirror in favor of the kwarg.  Net -15 LOC user
    code, same PNG output.
  - `src/PathNetwork.jl` gains a new module-header chapter "Schwarz-
    reflection symmetry (opt-in)" documenting the algorithmic correctness
    preconditions and the design choices (signed-zero collapse, recursive
    call into the unkwarged path, on-axis output handling).

Test suite: 1237 → **1250 GREEN** (+12 PN.6.* + 1 extra in PN.6.2 for
`length(visited_z) > 1` walked-tree check).  Wall ~1m47s.

### Beads filed this session

None.

### Open beads end-of-session (worklog 015)

Five P2 beads remain — same set as before, minus `padetaylor-dtj`:
  - `padetaylor-bvh` (Phase 13 CoordTransforms),
  - `padetaylor-grc` (Phase 14 SheetTracker),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction),
  - `padetaylor-0c3` (Fig 4.1 quantitative pin).

No P0 or P1 beads open.

### Hard-won lessons added this session (worklog 015)

20. **Julia `isequal(+0.0, -0.0) == false`** — and the same for
    `+0.0im` vs `-0.0im`.  Dict-keying on `Complex{Float64}` with
    on-axis cells is a footgun.  `complex(real(z), abs(imag(z)))`
    canonicalizes to `+0.0im` and eliminates the ambiguity.  Caught
    in PN.6.1's first GREEN attempt; documented in worklog 015 + the
    `_solve_with_schwarz_reflection` helper's docstring.

21. **Visited-tree-property tests need grids that force walking**.
    PN.6.2's first cut used a grid entirely within `h` of the IC ⇒
    Stage-1 took 0 steps ⇒ `visited_z == [IC]` regardless of mutation
    ⇒ Mutation I didn't bite the test.  Rewrote with off-axis targets
    > `h` from IC to force multi-step walks.  Lesson: when asserting
    properties of the visited tree, design the grid to grow it.

22. **Promote validated user-space workarounds to library API**.
    Worklog 014's `examples/tritronquee_3d.jl` carried the Schwarz-
    reflection workaround as 15 LOC of grid-partitioning user code.
    This session lifts it to an opt-in kwarg — the workaround was
    already exercised on the canonical PI tritronquée case, so the
    library version inherits the empirical validation.  Test plan
    just adds the bit-exact-symmetry invariants the user-space code
    couldn't easily assert.

### Session 2026-05-13 evening — all 4 P1 beads closed

Four GREEN commits land in this session, knocking out every P1 bead:

  - `9d1cc4c` **EdgeDetector** (Phase 12.5, `padetaylor-c2p` closed):
    FW 2011 §3.2.2 5-point Laplacian pole-field classifier.  60-LOC
    module, 743 assertions, 2 mutations.  Threshold on `log₁₀|Δu|`
    per FW Fig 3.3 (not bare `|Δu|`; bead description had a misread,
    documented in worklog 011).
  - `1c7056e` **Phase 9 PI tritronquée qualitative** (`padetaylor-kvi`
    closed): 5 testsets / 12 assertions verifying 4-of-5 pole-free
    sectors, conjugate symmetry, leading-pole magnitude at 25×25
    [-4,4]².  Mutation-proven (Mutation C bites 7/12 by flipping the
    ODE z-sign).  No source changes per DESIGN.md §4 "no new code".
  - `cc8c804` **Figure-catalogue refresh** (`padetaylor-rgp` closed as
    living document): `docs/figure_catalogue.md` §6 + §8 brought
    current; PathNetwork/EdgeDetector/BVP/Dispatcher v1 marked
    shipped; Fig 3.1 row marked PARTIAL with Phase 9 acceptance.
  - `4db29b1` **Phase 12 v2 LatticeDispatcher** (`padetaylor-k31`
    closed): 2D-grid composition machinery (~115 LOC) — PathNetwork +
    EdgeDetector partition + per-row BVP fill per FW2011...md:190.
    4 testsets / 76 assertions, mutation-proven (Mutations E + F).
    Tested against cosh closed form to ≤ 1e-10.

Test suite: 404 → **1237 GREEN** (+833 assertions across the four
shipped phases).  Wall ~1m45s.

### Beads filed this session

  - `padetaylor-0c3` **P2** Phase 12 v2.1: FW Fig 4.1 quantitative pin
    (`u(20i)` ≤ 1e-10) for `lattice_dispatch_solve`.  Different
    compositional pattern (vertical BVP + two outward pole fields)
    than v1's line-190 horizontal-row algorithm.  See worklog 013
    §"v1 scope decision".
  - `padetaylor-dtj` **P2** PathNetwork: add `enforce_real_axis_symmetry`
    kwarg.  `path_network_solve`'s `shuffle(rng, targets)` breaks
    conjugate symmetry for real-coeff/real-IC problems.  Workaround
    lives in `examples/tritronquee_3d.jl`; library-level kwarg is the
    proper fix.  See worklog 014.

### Open beads end-of-session

All remaining beads are P2:
  - `padetaylor-bvh` (Phase 13 CoordTransforms),
  - `padetaylor-grc` (Phase 14 SheetTracker),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction),
  - `padetaylor-0c3` (Fig 4.1 quantitative pin, new this session).

No P0 or P1 beads remain open.

### Hard-won lessons added this session (worklogs 011-013)

14. **Trust the source paper over the bead description**.  `padetaylor-c2p`'s
    bead said `|∇²u|/h² > 0.001`; FW2011...md:208 verbatim is
    `log₁₀|Δu| > 0.001` (the `/h²` is inside `Δu` per eq. 3.3).  See
    worklog 011 §"Threshold interpretation".

15. **Probe orientation before writing test assertions on 2D-grid output**.
    For PI `u'' = 6u² + z`, the tritronquée's pole-bearing sector is
    centred on the *positive* real axis (angle 0), not negative.  The
    ASCII-dump-first heuristic caught this in seconds vs the
    theory-first approach that would have wasted 30+ min.  See worklog
    012 §"A correction to my own initial assumption".

16. **v2-bead acceptance criteria often blend "machinery" and "specific
    case"**.  `padetaylor-k31` fused the 2D-lattice composition layer
    with the FW Fig 4.1 acceptance pin.  These are different
    deliverables (different compositional patterns).  v1 ships the
    machinery; v2.1 ships the specific case.  See worklog 013 §"v1
    scope decision".

17. **`shuffle(rng, targets)` in PathNetwork breaks conjugate symmetry
    even when the underlying ODE preserves it**.  For real-coeff
    real-IC problems, the shuffle's random target order creates an
    asymmetric visited tree whose Stage-2 cascade can blow `|u|` apart
    by 4-5 orders of magnitude at conjugate-pair cells.  Fix: walk
    only upper-half + mirror via `u(z̄) = ū(z)`.  See worklog 014 §
    "Bug 1".

18. **Use even N for 2D Cartesian grids centred on a symmetry axis**.
    With odd N, the centre cell sits on the axis and gets walked
    along a special direction, producing path-dependent discontinuity
    from off-axis cells.  Even N avoids the special-case row.  See
    worklog 014 §"Bug 2".

19. **FW Fig 4.1 geometry ≠ line-190 horizontal-row algorithm**.
    FW Fig 4.1's vertical BVP with asymptotic BCs serves *open*
    smooth regions extending to infinity; line-190's horizontal-row
    BVP fill serves *interior* smooth runs bounded on both sides by
    pole fields.  Pure tritronquée has no bridgeable interior runs
    — the horizontal-row cure does not apply.  See worklog 014
    §"BVP-cure exploration".

### Examples directory shipped

`examples/` ships two working scripts + its own Project.toml/Manifest.toml:
  - `tritronquee_3d.jl` — 3D surface + 2D heatmap + EdgeDetector mask
    for PI tritronquée at 500×500 over `[-20, 20]²` in ~17 s.  Uses
    upper-half walk + conjugate mirror (Bug 1 workaround) + even N
    (Bug 2 workaround).  Reproduces the canonical FW Fig 3.1
    tritronquée pole-field qualitatively.
  - `tritronquee_bvp_compose.jl` — investigative reference for the
    "BVP-fill cures tritronquée artifact" claim (does not — 0 BVPs
    trigger on the pure tritronquée geometry).

PNGs are gitignored; regenerate via `julia --project=examples
examples/<name>.jl`.

Goodluck. Read CLAUDE.md again before you start.
