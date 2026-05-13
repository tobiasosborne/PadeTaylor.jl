# Worklog 010 — Phase 8: PadeTaylorArblibExt package extension

**Date**: 2026-05-13
**Author**: Claude Opus
**Scope**: Phase 8 SVD-only Arb support (`PadeTaylorArblibExt`) per
ADR-0003 + ADR-0002.  Routes `pade_svd(::Matrix{Arb})` and
`pade_svd(::Matrix{Acb})` through BigFloat / Complex{BigFloat} +
GenericLinearAlgebra Jacobi.  404/404 tests GREEN (374 prior + 30 new
across 4 ArblibExt testsets).  Bead `padetaylor-jhq` closes.

> Take-home: a deliberately scope-bounded extension — only the SVD
> dispatch goes through, not the full IVP pipeline.  `Polynomials.roots(::Polynomial{Arb})`
> is the blocking friction (`padetaylor-8pi`) for end-to-end Arb
> integration, but the SVD layer alone is independently valuable for
> users doing arb-prec Padé approximation directly.

## Why SVD-only

`Arblib.jl` provides `Arb` (real ball) and `Acb` (complex ball) ball-
arithmetic types but does NOT ship an SVD primitive
(`RESEARCH.md §5.1` — confirmed by source inspection).  The natural
route is **BigFloat dispatch with radius discard** per ADR-0002:

  1. `BigFloat.(A)` extracts mid-points.
  2. `GenericLinearAlgebra.svd` runs Jacobi at current global
     `setprecision(BigFloat)`.
  3. Results lifted back to `Arb` / `Acb` (radius zero).

This is enough to make `RobustPade.robust_pade(::Vector{Arb})` work
when invoked directly — useful for users computing Padé approximants
in arb precision without doing IVP integration.

What's NOT included (deliberate v1 deferrals):

  - `Polynomials.roots(::Polynomial{Arb})` — blocked by
    `padetaylor-8pi`.  Affects `StepControl.step_pade_root` only;
    `step_jorba_zou` (the default) doesn't need it.  So `solve_pade`
    works at BigFloat but NOT at Arb if `step_pade_root` is used.
  - `default_tol(::Type{Arb})` — defer to user-supplied `tol` kwarg.
  - Per-Arb-value precision routing — the SVD runs at the global
    `setprecision(BigFloat)`; document this caveat.

## What changed

`Project.toml`: Arblib was being auto-promoted from `[weakdeps]` to
`[deps]` by a `Pkg.add` during probing.  Reverted; Arblib stays in
`[weakdeps]` + `[extensions]` per ADR-0003 (opt-in).  Also added
Arblib to `[extras]` + `[targets].test` so the test environment
resolves it.

`ext/PadeTaylorArblibExt.jl` (~110 LOC including literate docstring):
  - `pade_svd(A::AbstractMatrix{Arb}; full)` — BigFloat-routed.
  - `pade_svd(A::AbstractMatrix{Acb}; full)` — Complex{BigFloat}-routed.
  - Both return `Arb`/`Acb` outputs with radius zero (singular values
    are real, so `S` is always `Vector{Arb}`).

`test/ext_arblib_test.jl` (4 testsets, 30 new assertions):
  - AB.1.1: `pade_svd(::Matrix{Arb}) ≡ pade_svd(BigFloat.(A))`
    bit-identically (radii are zero by construction).
  - AB.1.2: same for `::Matrix{Acb}` vs `Complex{BigFloat}` path.
  - AB.2.1: SVD round-trip A ≈ U·diag(S)·Vt for symmetric Arb input;
    max reconstruction error < 1e-60 at BF-256.
  - AB.3.1: `full=true` returns n×n Vt for the GGT 2013 null-vector
    path (load-bearing for `RobustPade.robust_pade`).

## Mutation-proof (verified 2026-05-13)

**Mutation I** — downcast `A_bf` to `Matrix{Float64}` before the
Jacobi call.  This is THE correct precision-loss mutation per worklog
001 §"GenericLinearAlgebra svd! piracy": the naïve `s/Generic/Linear`
swap doesn't bite because GLA pirates `LinearAlgebra.svd!`; the bite
comes from forcing precision loss via Float64 downcast.

Verified bite: 7 fails across AB.1.1 (S/U/Vt mid-point mismatch at
Float64 precision ~1e-15 vs BigFloat ~1e-65), AB.2.1 (reconstruction
error 1.03e-15 >> 1e-60 BF-256 floor), AB.3.1 (Vt sign/precision flip
at 3 entries).  Confirms BF-256 dispatch is load-bearing.

**Mutation J** — replace `F.U` with `F.U'` (transpose) when lifting.
Breaks SVD identity `A = U·S·Vt`.

Verified bite: 2 fails — AB.1.1 (U mid-points don't match BigFloat
path), AB.2.1 (reconstruction error 3.71 vs 1e-60 target).

Both mutations restored before commit per CLAUDE.md Rule 4.

## Frictions surfaced

  - **F1. Pkg.add silently promotes weakdeps to deps**: my probe ran
    `Pkg.add("Arblib")` in the project env, which silently moved Arblib
    from `[weakdeps]` to `[deps]`.  This would have broken ADR-0003's
    opt-in design.  Reverted by re-editing `Project.toml`.  **Lesson**:
    when probing optional deps, use a temporary env
    (`Pkg.activate(; temp=true)`) or rely on `Pkg.test()` for resolution.
  - **F2. `@test all(BigFloat(s) for s in S_arb) .== S_bf` operator
    precedence**: the generator-`all` form returns a Bool, then `.==`
    tries `Bool .== Vector{BigFloat}`, raising `TypeError`.  Fixed by
    `BigFloat.(S_arb) == S_bf`.  Caught on first GREEN attempt.
  - **F3. Per-precision Arb vs global BigFloat precision**: Arb values
    carry their own precision; `BigFloat.(arb)` uses the global
    `setprecision(BigFloat)`, NOT the Arb's per-value precision.  For
    consistent results, set `setprecision(BigFloat, p)` to match.
    Documented in the ext docstring as a v1 caveat.

## Bead state

Closed in this session:
  - `padetaylor-jhq` — Phase 8 ArblibExt (SVD-only) GREEN at this commit.

Still open from prior sessions (no change):
  - `padetaylor-rgp`, `padetaylor-c2p`, `padetaylor-k31`,
    `padetaylor-kvi`, `padetaylor-bvh`, `padetaylor-grc`,
    `padetaylor-61j`, `padetaylor-8pi`.

`padetaylor-8pi` (Polynomials.roots-over-Arb friction) remains open
as the **blocker for end-to-end Arb IVP integration**.  Phase 8 v1
delivers SVD-only Arb support; full Arb IVP is a future v2 item once
that bead is resolved (e.g., by routing `step_pade_root` through
BigFloat at Arb input, or replacing with `step_jorba_zou`).

## Pointers

  - [`Project.toml`](../../Project.toml) — `[weakdeps]` + `[extensions]`
    + `[extras]` `Arblib` entries.
  - [`ext/PadeTaylorArblibExt.jl`](../../ext/PadeTaylorArblibExt.jl)
    — the extension.
  - [`test/ext_arblib_test.jl`](../../test/ext_arblib_test.jl)
    — 4 testsets + mutation-proof commentary.
  - [`docs/adr/0003-extensions-pattern.md`](../adr/0003-extensions-pattern.md)
    §"PadeTaylorArblibExt.jl" — design.
  - [`docs/adr/0002-bigfloat-svd-via-genericlinalg.md`](../adr/0002-bigfloat-svd-via-genericlinalg.md)
    — radius-discard caveat.
  - `RESEARCH.md §5.1` — Arblib SVD landscape (none).
