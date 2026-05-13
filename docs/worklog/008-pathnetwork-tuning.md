# Worklog 008 — Phase 10 follow-up: PathNetwork PN.2.2 + PN.3.1 (the wedge-vs-canonical-Padé bug)

**Date**: 2026-05-13
**Author**: Claude Opus (orchestrating)
**Scope**: PathNetwork Phase-10 follow-up tuning per bead `padetaylor-yt1`.
PN.2.2 (FW 2011 Table 5.1 z=30 to ≤1e-13 BF-256) + PN.3.1
(`:steepest_descent ≡ :min_u` to ≤1e-10).  340/340 tests GREEN
(329 prior + 11 new across the two testsets).  Bead `padetaylor-1jf`
(PathNetwork module) and `padetaylor-yt1` (this follow-up) both close.

> **Take-home — the headline result**: u(z=30) BF-256 rel-err = **2.13e-14**,
> *beating FW 2011's own published 8.34e-14*.  Fixing one
> two-line bug in PathNetwork's visited-node Padé storage turned
> ~22% rel-err at z=30 into 4.46e-11 (Float64) and 2.13e-14 (BF-256).
> The bug was hiding in plain sight in ADR-0004's design vs.
> implementation gap.

## The wedge-vs-canonical-Padé bug at long range

**ADR-0004's design** (cited in worklog 006 §"What changed — Phase 10"):
> "Each visited node gets one Padé centered there with real `h`, in
> addition to the wedge-direction Padé used for selection.  Cost: one
> extra Padé per node (`Coefficients` + `RobustPade`).  Benefit: Stage 2
> always interpolates inside `|t| ≤ 1`, never extrapolates past `t = 1`."

**Phase-10's implementation** (commit `0ada60f`) computed the canonical
Padé **only for the IC `z_0`** and stored it as `pade_0`.  For all
subsequent visited nodes (k > 1), the code stored the `pade_sel` from
the **wedge-direction step** that landed at that node — i.e., a Padé
built in a *complex* direction (e.g., the wedge-selected 45° off-axis
step), not the canonical real-h direction.

Stage 2 then evaluated `t = (z_f - z_v) / h_v` with `h_v = h_T` (the
real canonical step length), passing this `t` to a Padé that was
**built in a different direction**.  The Padé is valid in its own
rescaled `t ∈ [0, 1]` axis, which for a complex `h_step` is the
*complex* direction.  Evaluating it with a `t = real` corresponds to
projecting the (complex) Padé onto a (real) direction it wasn't built
to cover.

### Empirical signature of the bug

At z=30 with `h=0.5, order=30, :min_u`:
  - Float64 rel-err = 0.218 (22%)
  - BF-256 rel-err = 0.218 (**identical to Float64**)
  - Order=30, 40 give identical results — **insensitive to truncation order**
  - Imag part of u(30) = 7.16e-3 (significantly non-zero despite real ODE
    + real IC + real grid point)

These signatures are the bug's diagnostic fingerprint: **the error is
algorithmic, not arithmetic**.  If it were Padé truncation error, raising
order from 30 to 40 would reduce it.  If it were floating-point roundoff,
BF-256 would reduce it.  Neither did.  The error is in the Padé's *direction
of validity* mismatching the Stage-2 evaluation's *direction of use*.

At z=1.4 (just two steps from IC, past the lattice pole), the bug shows
1.50e+0 rel-err: u(1.4) returns -3.10 - 1.52i where the closed-form ℘
gives 6.25 (real).  This was hiding behind PN.1.2's lenient acceptance
("u(target) finite, magnitude < 1e3").

### The fix — two lines

In `src/PathNetwork.jl` Stage-1 loop, replace:

```julia
push!(visited_pade, pade_sel)
```

with:

```julia
pade_canonical = _local_pade(prob.f, z_new, u_new, up_new, order, h_T)
push!(visited_pade, pade_canonical)
```

The `_local_pade` helper already exists (it was used for the IC only);
it builds a fresh canonical Padé from `(z_new, u_new, up_new)` with the
real step length `h_T`.  Cost: one extra `Coefficients.taylor_coefficients_2nd`
+ one extra `RobustPade.robust_pade` per visited node.

### Empirical post-fix at long range

| target z | Float64 rel-err | imag(u) | notes |
|---|---|---|---|
| 0.5  | 4.4e-15 | 0      | trivial Padé-disc eval at IC |
| 1.4  | 1.2e-13 | 5.7e-13 | first pole bridge (z≈1.13) |
| 30.0 | 4.5e-11 | 4.8e-11 | 75 visited nodes, 14 pole bridges |

| target z | BF-256 rel-err | imag(u) | notes |
|---|---|---|---|
| 30.0 | **2.13e-14** | 1.2e-18 | **beats FW 2011's published 8.34e-14** |

Per `docs/figure_catalogue.md §1 row FW2011 Fig 5.1`: the headline
v1 acceptance for PN.2.2 was "≤1e-13 BF-256".  We achieve 2.13e-14,
roughly one decimal better than the FW reference.

## What changed

`src/PathNetwork.jl`: the two-line fix described above (plus a literate
comment block explaining ADR-0004's design + the cost rationale).

`test/pathnetwork_test.jl`:
  - **PN.2.2**: FW Table 5.1 long-range — F64 to ≤1e-9 + imag<1e-9;
    BF-256 to ≤1e-13 + imag<1e-15.  4 assertions; runs in ~50s wall
    (BF-256 dominates).
  - **PN.3.1**: `:steepest_descent ≡ :min_u` on a non-trivial pole-bridge
    grid {1.4+0i, 1.2+0.4i, 0.6+0.3i}.  The original "PN.1.1's tree"
    setup turned out to be trivial (0 steps required, all targets within
    h of IC); we substituted a grid that REQUIRES Stage-1 walking so
    step-selection is actually exercised.  Tolerance 1e-10 on Δu, 1e-9
    on Δu' (the latter looser per the 1/h chain-rule factor in path-network).

Combined: 11 new assertions; test suite grows 329 → 340 GREEN.
Total wall: 1m32s (PN.2.2 BF-256 is the new bottleneck).

## Mutation-proof (verified 2026-05-13)

**Mutation E** — restore pre-fix wedge-direction Padé storage:
`pade_canonical = pade_sel` (drops the extra `_local_pade` call).
Verified bite: 4 fails on PN.2.2 — F64 rel-err 0.218 >> 1e-9, F64
imag(u) 7.2e-3 >> 1e-9, BF-256 rel-err 0.218 >> 1e-13, BF-256 imag
7.2e-3 >> 1e-15.  The bug's algorithmic-vs-arithmetic invariance
is itself verified by the BF-256 fail (identical error magnitude).

**Mutation F** — in `:steepest_descent`, flip sign: `θ_sd = angle(u/up)`
instead of `angle(-u/up)`.  Steers toward poles.  Verified bite: 4
fails on PN.3.1 — at z=1.4 the mutated steepest-descent path fails
to bridge the pole and returns u=0+0im while :min_u gives u=6.2518,
Δu=6.25 >> 1e-10.  |Δu'| reaches 31.2.

Both mutations restored before commit per CLAUDE.md Rule 4.

## Frictions surfaced

  - **F1. Spec drift between ADR and impl, caught only by long-range test**:
    the wedge-vs-canonical-Padé bug was *specified correctly in ADR-0004
    and worklog 006* (which explicitly describes the canonical-Padé design)
    but *implemented incorrectly* in Phase-10's commit `0ada60f`.  PN.1.2
    didn't catch it because its z=1.4 assertion was structural ("finite,
    not divergent"); PN.2.1 didn't catch it because near-IC tests use
    the IC's canonical Padé directly.  Only PN.2.2 at z=30 (the FW
    Table 5.1 acceptance) was tight enough to surface the bug.  **Lesson**:
    when the ADR specifies an architectural invariant, a test must
    *verify the invariant*, not its near-IC approximation.  Phase-10's
    test corpus had this hole.

  - **F2. PN.3.1's original "near-IC" framing was trivial**:
    the bead specified `:steepest_descent ≡ :min_u to ≤1e-10 on PN.1.1's
    tree`.  PN.1.1's grid is entirely within `h` of the IC, so Stage 1
    takes ZERO steps regardless of selection rule, hence both rules
    agree trivially.  Mutation F (flipping the sign in θ_sd) does NOT
    bite a 0-step path.  We substituted a 3-target grid that requires
    1-3 steps each, including z=1.4 past the lattice pole; both rules
    agree to ~1e-12 in GREEN, and Mutation F triggers a 6.25 mismatch.
    **Lesson**: test grids must exercise the code path under test;
    near-trivial grids generate confidence-theatre tests.

  - **F3. The bug's symptom misled the initial investigation**:
    seeing `u(30) ≈ 0.856` (vs reference 1.095) with imag ≈ 7e-3
    initially suggested "the path-network went around poles in such a
    way that it computes a different ℘ branch", or "compound error
    along a 60-step path".  Neither is right.  The clue was BF-256
    giving the SAME 22% error as F64 — algorithm-not-arithmetic.  Then
    grep for `pade_canonical` vs `pade_sel` storage in the Stage-1
    loop revealed the design-vs-impl gap.  **Lesson**: when a numerical
    test fails by orders of magnitude, do not assume cumulative
    floating-point error.  Probe with BF-256 first; if BF-256 doesn't
    help, the bug is algorithmic.

## Pointers

  - [`src/PathNetwork.jl`](../../src/PathNetwork.jl) — Phase-10 module
    with the canonical-Padé bugfix at the Stage-1 loop.
  - [`test/pathnetwork_test.jl`](../../test/pathnetwork_test.jl) —
    PN.2.2 + PN.3.1 testsets + mutation-proof commentary.
  - [`docs/adr/0004-path-network-architecture.md`](../adr/0004-path-network-architecture.md)
    §"Per-node canonical Padé" — the design decision that the impl
    initially missed.
  - [`docs/worklog/006-phases-10-11-path-network-bvp.md`](006-phases-10-11-path-network-bvp.md)
    §"What changed — Phase 10" — describes the design correctly while
    the impl shipped with the bug.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-166`
    — FW 2011 §3.1 path-network original spec.
  - `docs/figure_catalogue.md §1 row FW2011 Fig 5.1` — the acceptance
    criterion ≤1e-13 BF-256 that PN.2.2 now meets.

## Bead state

Closed in this session:
  - `padetaylor-yt1` — PN.2.2 + PN.3.1 follow-up; GREEN at this commit.
  - `padetaylor-1jf` — Path-network module (held open by PN.2.2 + PN.3.1
    pending; now released).

Still open from prior sessions (no change):
  - `padetaylor-8cr`, `padetaylor-rgp`, `padetaylor-c2p`, `padetaylor-k31`,
    `padetaylor-kvi`, `padetaylor-jhq`, `padetaylor-2vz`, `padetaylor-bvh`,
    `padetaylor-grc`, `padetaylor-61j`, `padetaylor-8pi`.

## Hard-won lessons (for HANDOFF.md §"Hard-won")

**12. ADR-specified invariants need invariant tests, not approximation
tests**.  Phase-10's PN.1.2 acceptance was "u(1.4) finite, magnitude
< 1e3" — a structural sanity check that passes whether u(1.4) = 6.25
(correct) or u(1.4) = -3.10 - 1.52i (≈22% rel-err due to the
canonical-Padé bug).  ADR-0004 specified the canonical-Padé-per-node
invariant, but no test asserted "evaluating visited[k]'s Padé at t=0
returns visited_u[k]" — a single direct check that would have caught
the bug immediately at Phase-10 GREEN.  The long-range FW Table 5.1
test at z=30 caught it later as a downstream symptom; we should have
caught it as the symptom-of-its-own-direct-test.

**13. When numerical error is invariant to precision, suspect the
algorithm**.  Float64 rel-err = BF-256 rel-err is a strong signal that
the issue is structural, not roundoff.  Pre-emptive BF-256 probe is
worth the ~50× wall-time cost when a result smells wrong.
