# ADR-0002 — Bigfloat SVD via GenericLinearAlgebra (one-sided Jacobi)

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
| `BigFloat` | `GenericLinearAlgebra.svd` | One-sided Jacobi (Demmel-Veselić) | Relative-accuracy guarantee `c · 2⁻ᵖ · σᵢ` per SV — load-bearing for arb-prec rank counting near a Padé table block boundary. |
| `Arblib.Arb` | (extension) Convert to `BigFloat`, dispatch as above | One-sided Jacobi | `Arblib.jl` ships no SVD whatsoever (verified by source inspection — see `RESEARCH.md §5.1`). The `Arb → BigFloat` conversion is precision-lossy on the radius but acceptable for the SVD step alone (see "Caveats" below). |
| `Complex{T}` for `T <: AbstractFloat` | `GenericLinearAlgebra.svd` (or stdlib if `T = Float64`) | as above | Same dispatch logic. |

## Why one-sided Jacobi for arb-prec

For Float64 and the GGT-typical regime (matrices of size `n ≤ 60`,
condition numbers up to `10¹²`), Demmel-Kahan implicit-shift QR is
the right choice — used by Chebfun's `padeapprox.m` line 93 (`svd`).

**At arbitrary precision the regime changes.** GGT 2013 Algorithm 2
classifies a singular value as "zero" if `σ_i < tol · ||c||_2`. Near a
block boundary in the Padé table, the gap between the smallest
"signal" SV and the largest "noise" SV may span only a few orders of
magnitude — far less than the dynamic range of the working precision.

- **Demmel-Kahan** guarantees absolute error `c · 2⁻ᵖ · σ_max` on each
  singular value. Small-but-genuine SVs may be perturbed below the
  threshold and incorrectly classified as zero, causing spurious
  block reduction.
- **One-sided Jacobi (Demmel-Veselić 1992)** guarantees *relative*
  error `c · 2⁻ᵖ · σ_i` per singular value. Small SVs stay reliably
  above the threshold, regardless of `κ(A)`.

For GGT matrices at `n ≤ 60` the runtime overhead of Jacobi is
negligible (workbench's own SVD dispatcher routes `n ≤ 500` to Jacobi
precisely on this trade-off; see `RESEARCH.md §5.3`).

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
  matrices at our typical sizes. If we hit an instance, fall back
  to a hand-rolled Demmel-Veselić port; do not reach for Demmel-Kahan
  at arb-prec.

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
- `RESEARCH.md §5.3` — algorithm-choice argument: relative-accuracy
  Jacobi vs absolute-accuracy Demmel-Kahan, applied to GGT's
  rank-counting threshold.
- `external/chebfun/padeapprox.m:93, 106` — Chebfun's MATLAB `svd`
  call (= LAPACK Demmel-Kahan); we deliberately diverge from this
  for `T ≠ Float64`.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
  GGT2013_robust_pade_via_SVD_SIREV55.md:213–217` — the `tol · ||c||₂`
  thresholding that makes relative-accuracy SVD load-bearing at
  arb-prec.

## Consequences

- `LinAlg.jl` is ~60 LOC: dispatch table + `Arb → BigFloat` shim in
  the extension.
- `GenericLinearAlgebra.jl` becomes a hard `Project.toml` dep
  (alongside `LinearAlgebra` from stdlib).
- The `Arb`-element-type path lives in `PadeTaylorArblibExt` (per
  ADR-0003), pulled in only when the user has `Arblib.jl` loaded.
- Any future "rigorous Arb SVD" work is deferred and tracked in
  `RESEARCH.md §8`.
