# ADR-0002 — Bigfloat SVD via GenericLinearAlgebra

**Status**: Accepted (2026-05-09)
**Context**: Stage 1 design lock; the GGT 2013 Padé routine requires
SVD of an `n × (n+1)` Toeplitz matrix at every Padé conversion. For
the arb-prec tier (`T ∈ {BigFloat, Arblib.Arb}`), Julia's stdlib
`LinearAlgebra.svd` throws `MethodError`. We need a path.

## Decision

`PadeTaylor.LinAlg.pade_svd(A::AbstractMatrix{T})` dispatches:

| `T` | Backend | Algorithm | Rationale |
|---|---|---|---|
| `Float64` | `LinearAlgebra.svd` | LAPACK Demmel-Kahan (`DGESVD`) | Matches Chebfun's `padeapprox.m` exactly; GGT's `tol = 1e-14` rank-counting tolerates DK at Float64 precision. |
| `Float32` | `LinearAlgebra.svd` | LAPACK Demmel-Kahan (`SGESVD`) | Same reasoning at `Float32`. |
| `BigFloat` | `GenericLinearAlgebra.svd` | Householder bidiagonalisation + bidiagonal QR (including Demmel-Kahan zero-shift iterations) | The default `tol = 2^(-p+10)` keeps the rank threshold above the backend's absolute error scale for the small GGT matrices used here. |
| `Arblib.Arb` | (extension) Convert to `BigFloat`, dispatch as above | As above | `Arblib.jl` ships no SVD whatsoever (verified by source inspection — see `RESEARCH.md §5.1`). The `Arb → BigFloat` conversion is precision-lossy on the radius but acceptable for the SVD step alone (see "Caveats" below). |
| `Complex{T}` for `T <: AbstractFloat` | `GenericLinearAlgebra.svd` (or stdlib if `T = Float64`) | as above | Same dispatch logic. |

## Why the generic SVD is sufficient for arb-prec

For Float64 and the GGT-typical regime (matrices of size `n ≤ 60`,
condition numbers up to `10¹²`), Demmel-Kahan implicit-shift QR is
the right choice — used by Chebfun's `padeapprox.m` line 93 (`svd`).

At arbitrary precision, GGT 2013 Algorithm 2 classifies a singular
value as zero if `σ_i < tol · ||c||_2`. The implementation chooses
`tol = 2^(-p+10)`, ten bits above the working-precision floor. For an
`n × (n+1)` GGT Toeplitz block, each coefficient occurs in at most `n`
entries, so `σ_max ≤ ||C̃||F ≤ √n · ||c||_2`. The backend's usual
`O(2⁻ᵖ · σ_max)` absolute singular-value error is therefore below the
threshold by a factor of order `2¹⁰/√n`; at the typical `n ≤ 60`, this
leaves ample margin for rank counting. CRP.4 pins the outcome on the
BigFloat corpus (`test/corpus_robust_pade_test.jl:190-201`).

Source inspection shows that `GenericLinearAlgebra.svd` first performs
Householder bidiagonalisation and then solves the bidiagonal problem by
QR, selecting the Demmel-Kahan zero-shift iteration when needed. It is
not a one-sided Jacobi implementation (`~/.julia/packages/
GenericLinearAlgebra/X90Kh/src/svd.jl:328-345,648-654,81-87,235-244`).

## Why `GenericLinearAlgebra.jl` is the chosen library

- **It is the only production-grade Julia library that provides
  generic SVD over `T <: AbstractFloat` including `BigFloat`** (see
  `RESEARCH.md §5.1`).
- **Active maintenance** — CI-tested against Julia 1.10+; used in
  arb-prec research workflows.
- **No FFI** — pure Julia; works on any platform Julia supports
  (load-bearing for the PRD's "phone deployment" forcing function
  if we ever target it).
- **Open issues at `κ > 10¹⁸`** — documented; not relevant for GGT
  matrices at our typical sizes. If such an instance appears, reassess
  the backend against a pinned high-precision oracle rather than
  assuming a different algorithm is already in use.

## Caveats

### Arb → BigFloat conversion is lossy on radius

Converting `Arb(mid ± rad)` to `BigFloat` discards `rad`. The SVD
returns `BigFloat` matrices whose entries are *no wider than* the
input mid-points; the Arb-radius information about uncertainty in
the input matrix is **lost** through the SVD step.

**Mitigation**: this is acceptable because the Padé routine's
downstream consumer is `RobustPade.robust_pade`, which uses the SVD
results to produce a `PadeApproximant{T}` whose interpretation
already rests on the GGT 2013 normalisation `||b||₂ = 1`. The
reported precision of the Padé approximant is set by the
*coefficient* arithmetic (where Arb radii are correctly tracked),
not by the SVD step itself.

If we ever want a fully Arb-rigorous SVD, we'd need to port a verified-
arithmetic SVD algorithm — well outside v1 scope. Documented as
deferred work in `RESEARCH.md §8`.

### `GenericLinearAlgebra.svd` over `Arb` element type — open question

`RESEARCH.md §8 Q4` flags this. The current path is to convert
`Matrix{Arb}` → `Matrix{BigFloat}` before calling
`GenericLinearAlgebra.svd`. **If empirical testing shows direct `Arb`
dispatch works** (the scalar interface required is `abs`, `sqrt`,
`/`, all of which `Arb` provides), we can skip the conversion in
the extension package. Verify in Phase 8.

## Citations

- `RESEARCH.md §5.1` — empirical confirmation of the Julia SVD
  landscape (no SVD in `Arblib.jl`; `GenericLinearAlgebra.jl` as the
  only path; `LinearAlgebra.svd` throws `MethodError` on `BigFloat`).
- `RESEARCH.md §5.3` — the absolute-accuracy margin supplied by the
  ten-bit gap in `default_tol`, applied to GGT's rank-counting threshold.
- `external/chebfun/padeapprox.m:93, 106` — Chebfun's MATLAB `svd`
  call (= LAPACK Demmel-Kahan); we deliberately diverge from this
  for `T ≠ Float64`.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
  GGT2013_robust_pade_via_SVD_SIREV55.md:213-217` — the
  `tol · ||c||₂` threshold used for rank counting.

## Consequences

- `LinAlg.jl` is ~60 LOC: dispatch table + `Arb → BigFloat` shim in
  the extension.
- `GenericLinearAlgebra.jl` becomes a hard `Project.toml` dep
  (alongside `LinearAlgebra` from stdlib).
- The `Arb`-element-type path lives in `PadeTaylorArblibExt` (per
  ADR-0003), pulled in only when the user has `Arblib.jl` loaded.
- Any future "rigorous Arb SVD" work is deferred and tracked in
  `RESEARCH.md §8`.
