# Block-Toeplitz-determinant shared-Q oracle (`shared-pade-determinant-oracle`)

An **independent** reference implementation of the type-II Hermite–Padé
(simultaneous-Padé / shared-denominator) approximant, for mutation-proving
`PadeTaylor.SharedPade.shared_denominator_pade` (`src/SharedPade.jl`). This
is the bead **V1c** oracle; bead V1e wires it into the project test suite.

## Why a second route

`shared_denominator_pade` recovers the shared denominator `Q` as the **SVD
null vector** of a stacked block-Toeplitz matrix `A_full`, then refines it
with a column-reweighted **QR** step ported from `padeapprox.m`. An oracle
that *also* took an SVD — or used QR as a null-space extractor — would not
be independent: a sign-flip, an off-by-one in the Toeplitz layout, or a
wrong-column bug in the SVD+QR path could hide behind a matching SVD+QR
oracle.

This oracle recovers `Q` by a route with **no singular-value decomposition
and no QR-as-null-space step anywhere** — the determinantal / exact-linear-
algebra route of Mano–Tsuda 2017.

## Files

| File         | Role |
|--------------|------|
| `oracle.jl`  | Self-contained module `SharedPadeDeterminantOracle` exporting `shared_q_via_determinant(jets, m) -> (numerators, denominator)`. No project dependency. **Committed.** |
| `verify.jl`  | Verification script: runs the oracle on the SP.1.1 / SP.1.2 inputs, cross-checks against the exactly-known `Q` (Rational input) and against `SharedPade.shared_denominator_pade` (Float64 input). **Committed.** |
| `README.md`  | This file. |

## Reproduce

```sh
julia --project=../../.. external/probes/shared-pade-determinant-oracle/verify.jl
```

(One Julia process — CLAUDE.md Rule 7. No build step; `oracle.jl` is pure
Julia and depends only on `LinearAlgebra.det`.)

## The algorithm — block-Toeplitz determinant (no SVD)

Ground truth: Mano–Tsuda 2017 §2.2 "Simultaneous Padé polynomials",
Proposition 2.3, read directly from the arXiv LaTeX source
(`references/tex/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285/hp_arXiv_final.tex:1484-1528`);
see also `docs/v0p2_pillarA_hermite_pade_findings.md:116-134` §2.

Proposition 2.3 gives the type-II denominator polynomial `P^(0)_0` (the
object playing the role of `Q`) as a **bordered block-Toeplitz determinant**:

```
P^(0)_0(w) = (1/NP^(0)) · det [ 1, w, w², …, wᵐ      ]
                              [ A¹_m(n, m+1)          ]
                              [        ⋮              ]
                              [ A^{L-1}_m(n, m+1)     ]
```

Expanding that `(m+1)×(m+1)` determinant along its **top row** (Laplace /
cofactor expansion): the coefficient of `wᵏ` in `Q` is the signed cofactor

```
b_k  =  (-1)^k · det( A_full with column (k+1) deleted )  =  (-1)^k · M_k
```

where `M_k` is the `k`-th **maximal minor** of the `m × (m+1)` stacked
block-Toeplitz matrix — an `m × m` determinant. The vector
`b = (b_0,…,b_m)` is, up to a global scale, the shared-denominator
coefficient vector. It is a right null vector of `A_full` *by construction*:
each `(A_full · b)_r` is the Laplace expansion of a determinant with two
equal rows, hence exactly zero. That is the same defining property
`SharedPade` extracts via SVD — here it is a list of `m+1` determinants.

Cramer's rule is the identical fact in different dress: fixing `b₀ = 1` and
solving `A_red · (b₁,…,b_m)ᵀ = -A_full[:,1]` gives
`b_j = det(A_red, col j ↦ -A_full[:,1]) / det(A_red)`, and these ratios are
the same signed minors. The oracle uses the minor-list form directly — it
treats every coefficient symmetrically and never needs `b₀ ≠ 0`.

**No SVD. No QR null-space step.** Only determinants and (for the row
selection below) fraction-free Gaussian elimination.

### The `d ≥ 2` determinant subtlety (V1e must know this)

Maximal minors are defined only for an `m × (m+1)` matrix. Mano–Tsuda's
*minimal* simultaneous-Padé system (eq. 2.6, `hp_arXiv_final.tex:1405-1424`)
**is** exactly that shape: with block height `n` and `m = n(L-1)` the stack
`A¹_m(n,m+1);…;A^{L-1}_m(n,m+1)` has `(L-1)·n = m` rows.

`SharedPade`, however, uses the **GGT square-block convention** `n = m` rows
per block (`src/SharedPade.jl:31-44`), so its `A_full` is the *taller*
`dm × (m+1)` matrix. For `d = 1` that is already `m × (m+1)` and the minors
apply directly. For `d ≥ 2` the stack has `dm > m` rows — the extra rows are
**redundant matching equations** (an exact `P/Q` satisfies more than the
minimal `m`).

The oracle therefore first reduces `A_full` to an `m × (m+1)` **square
subsystem** by selecting `m` linearly independent rows
(`_select_independent_rows`), using fraction-free Gaussian elimination with
column-pivot search — *not* an SVD or QR rank-revealing factorisation. For
`Rational{BigInt}` input the selection is exact (a row is accepted iff its
residual leading entry is a nonzero rational). Because an exact `P/Q`
satisfies *every* matching equation, any `m` independent rows yield the same
null vector — the selection is bookkeeping, not approximation. This is the
direct implementation of the bead-text instruction "select `m` independent
rows of the stacked matrix" and of Mano–Tsuda's square-selection in §2.2.

If fewer than `m` independent rows exist, `rank(A_full) < m` and no unique
shared `Q` exists — the oracle throws (Rule 1; Mano–Tsuda p. 12, "unique iff
rank = m").

## Exact / arbitrary-precision arithmetic — the tight-oracle property

A determinant of an `m × m` matrix is a polynomial in its entries, so for
`Rational{BigInt}` input the minors `M_k` are computed **exactly** by
Julia's generic rational `LinearAlgebra.det` — no roundoff at all. The
oracle output is then the *exact* rational `Q` and `Pᵢ`. For `BigFloat`
input the determinants are evaluated at extended precision. This is what
makes V1c a *tight* oracle: it is not itself a source of numerical error,
so any discrepancy with `SharedPade` is a real `SharedPade` bug, not oracle
noise.

`shared_q_via_determinant` matches `SharedPade`'s output convention exactly:
coefficients low-to-high, denominator normalised to `Q(0) = 1`, trailing
near-zeros trimmed — so V1e can compare the two tuples directly.

## Verification results

`verify.jl` run 2026-05-20, Julia 1.10, on the SP.1.1 (`d=1` scalar
rational) and SP.1.2 (`d=2` shared-`Q`) inputs of `test/shared_pade_test.jl`.
Each case is fed both as `Rational{BigInt}` jets (exact route) and as
`Float64` jets (vs `SharedPade`).

### Exact route — `Rational{BigInt}` input vs the exactly-known `Q`

| Case   | Denominator `Q`              | Numerators `Pᵢ`            |
|--------|------------------------------|----------------------------|
| SP.1.1 | `den == known Q` **exactly** | `P₁ == known P₁` **exactly** |
| SP.1.2 | `den == known Q` **exactly** | `P₁, P₂ == known` **exactly** |

Bit-for-bit equality: the determinant oracle on exact rational input
reproduces the construction `Q` (and numerators) with zero error.

### Float64 route — oracle vs `SharedPade` vs ground truth

| Case   | `|Δden|` vs truth / SharedPade | pole roots vs truth / SP | values vs truth / SP |
|--------|-------------------------------|--------------------------|----------------------|
| SP.1.1 | 1.1e-16 / 5.6e-17             | 2.2e-16 / 1.1e-16        | 2.5e-16 / 4.6e-16    |
| SP.1.2 | 1.1e-16 / 6.7e-16             | 4.6e-16 / 2.0e-15        | 1.4e-17 / 4.4e-16    |

All observed agreements are at **machine epsilon** (~1e-15 or below) — far
inside the `1e-12` floor V1e will assert (the same floor as the Calgo 766
oracle in `external/probes/calgo766-oracle/`). The failure-mode check
(jet too short) throws with an informative message, as required by Rule 1.

**Conclusion:** the determinant oracle is sound. It agrees with
`SharedPade.shared_denominator_pade` to machine epsilon on both the scalar
and the `d=2` shared-`Q` cases, and reproduces the exactly-known `Q`
bit-for-bit on rational input — with **no SVD anywhere** in its code path.

## References

- `references/tex/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285/hp_arXiv_final.tex:1365-1528`
  — §2.2 Simultaneous Padé polynomials: eq. (2.6) (the homogeneous
  block-Toeplitz system) and Proposition 2.3 (the bordered block-Toeplitz
  determinant for `P^(0)_0`).
- `references/tex/.../hp_arXiv_final.tex:1322-1359` — Remark 2.2, the
  block-Toeplitz determinant `Δ^(i)` and the normalising constant `NP^(0)`.
- `docs/v0p2_pillarA_hermite_pade_findings.md:116-134` — §2 "Determinantal
  Representation", the same Prop. 2.3 formula.
- `src/SharedPade.jl` — the SVD+QR path being validated; `_toeplitz_block`
  / `_upper_block` are mirrored verbatim in `oracle.jl`.
- `external/probes/calgo766-oracle/` — the sibling V1b oracle (Cabay–Jones–
  Labahn FORTRAN); V1c is the determinantal counterpart.
