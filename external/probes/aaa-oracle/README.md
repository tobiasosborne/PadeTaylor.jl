# AAA per-component pole cross-check oracle (bead padetaylor-0ln.5 / V1d)

A **third, structurally independent oracle** for the shared-denominator
keystone `PadeTaylor.SharedPade.shared_denominator_pade`
(`src/SharedPade.jl`).

## What this probe is

`shared_denominator_pade` recovers a single shared denominator `Q` for a
`d`-component vector of meromorphic functions from their **Taylor jets**,
via the null space of a block-Toeplitz matrix (Mano–Tsuda 2017 §2.2). Two
oracles already cross-check it:

| Oracle | Probe dir | Mechanism |
|---|---|---|
| Calgo 766 | `external/probes/calgo766-oracle/` | FORTRAN Beckermann–Labahn recurrence, Taylor jets |
| Block-Toeplitz determinant | `external/probes/shared-pade-determinant-oracle/` | Closed-form determinant, Taylor jets |
| **AAA (this probe)** | `external/probes/aaa-oracle/` | Barycentric fit + Loewner SVD, **sampled values** |

Both prior oracles, like `SharedPade` itself, work from Taylor
coefficients. **AAA is different at every stage of the pipeline**, which is
the point of adding it.

## The independence argument: samples vs jet

AAA — the adaptive Antoulas–Anderson algorithm of Nakatsukasa, Sète &
Trefethen (SISC 40(3), 2018) — fits a rational approximant from **sampled
function values on a point set** `Z ⊂ ℂ`, never touching a Taylor jet.
Concretely, every component of the pipeline is replaced relative to the
Taylor-jet oracles:

- **Input:** sampled values `f(Z)`, not Taylor coefficients `[c₀, c₁, …]`.
- **Representation:** a **barycentric** partial-fraction quotient
  `r = n/d` over greedily chosen support points (NST §2,
  `p.tex:155–207`), not a monomial polynomial `P(z)/Q(z)`.
- **Linear-algebra core:** the SVD of a **Loewner** matrix
  `Aₘ[i,k] = (Fᵢ − fₖ)/(Zᵢ − zₖ)` (NST eq. 3.6, `p.tex:401–415`), not the
  SVD/null-space of a block-Toeplitz matrix.
- **Pole extraction:** the finite eigenvalues of an `(m+1)×(m+1)`
  **arrowhead generalised eigenvalue problem** (NST eq. 3.8,
  `p.tex:455–474`), not the roots of a coefficient vector.

AAA *does* use an SVD internally — that is the algorithm. The independence
claimed here is at the level of **"samples vs Taylor jet"** and
**"Loewner matrix vs block-Toeplitz matrix"**, not "no SVD". When AAA's
poles land on the same points as the roots of `SharedPade`'s shared `Q`,
that agreement comes from a wholly different computational route and is
therefore genuine cross-validation evidence.

## The per-component-independence limitation (diagnostic, not ground truth)

Per `docs/v0p2_pillarA_hermite_pade_findings.md:354–373` (§6 "AAA as a
benchmark"): **AAA fits each component independently and cannot enforce a
shared denominator.** For a `d`-component system it produces `d` separate
pole sets. So AAA cannot, on its own, certify the *coefficients* of the
shared `Q` — it is a **diagnostic pole cross-check, not a coefficient-level
ground truth** (`findings.md:368–373`).

The cross-validation evidence is the **coincidence** of those `d`
independent pole sets — with each other *and* with the shared-`Q` roots.
Because AAA fits each component on its own, with no knowledge that the
components share poles, the fact that all `d` per-component pole sets
nonetheless agree is independent confirmation that the shared denominator
recovered by `SharedPade` is the correct common singular structure
(`findings.md:364–368`). A discrepancy would instead flag a Froissart
doublet or a genuine component-specific singularity (`findings.md:369–372`).

## Files

- **`aaa.jl`** — a self-contained, literate Julia port of the AAA
  algorithm. Ported from the local LaTeX ground truth
  `references/tex/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40/`:
  `p.tex` §2–§3, the reference `aaa_alg.m`, and `cleanup.m`. Public entry
  point:

  ```julia
  aaa(Z, F; tol=1e-13, mmax=100, cleanup=true) -> AAAResult
  ```

  with `AAAResult` carrying `support_points`, `values`, `weights`,
  `poles`, `residues`, `zeros`, and the per-step `errvec`. The greedy
  support-point selection, the Loewner-matrix SVD weight update, the
  arrowhead pole extraction, and the Froissart-doublet cleanup pass are all
  implemented from the `.tex`. Throws (Rule 1 — fail loud) on length
  mismatch, fewer than two samples, or a non-finite sample value (a sample
  sat on a pole).

  No new dependency is added to `Project.toml`: `aaa.jl` uses only
  `LinearAlgebra`, which the project already depends on. AAA is a test-only
  oracle.

- **`verify.jl`** — the verification record. Samples the SP.1.1 and SP.1.2
  test functions from `test/shared_pade_test.jl`, fits each component
  independently with AAA, and cross-checks the AAA poles against the roots
  of `SharedPade.shared_denominator_pade`'s shared `Q`. Does **not** add
  asserts to the project test suite — bead V1e wires `aaa.jl` into `test/`
  separately.

## Sampling contour (why the unit circle)

AAA fits from values, so a sample point must not coincide with a pole — a
sample *at* a pole is `Inf`, and `aaa.jl` throws on a non-finite value
(Rule 1). Both test functions have their poles **outside the unit disk**:

| Case | Denominator `Q` | `|root|` |
|---|---|---|
| SP.1.1 | `1 − 0.5z + 0.3z²` | `1/√0.3 ≈ 1.826` |
| SP.1.2 | `1 − 0.4z + 0.2z²` | `1/√0.2 ≈ 2.236` |

The unit circle `|z| = 1` therefore lies safely inside the analyticity
region, where every component is smooth and bounded. `verify.jl` samples
**200 equispaced points** on `|z| = 1`; AAA's greedy selection then picks
its own 3 support points from that set. (NST §3 requires `m ≤ M/2`; 200
samples admit types up to `(99, 99)`, far above the type-`(2,2)` rationals
here.)

The failure-mode check in `verify.jl` deliberately samples on a contour
through a pole (radius `≈ 1.826` for SP.1.1) and confirms `aaa.jl` throws
on the resulting `Inf` value.

## Spurious-pole filtering

NST identify numerical Froissart doublets as poles with residue
`|res| < 1e-13` (`p.tex:652–661`, `cleanup.m:3`) and remove them by
dropping the support point nearest each tiny-residue pole and re-solving
the least-squares problem once. `aaa.jl` exposes this as the `cleanup`
keyword (default `true`). For the SP.1.1/SP.1.2 inputs — exact type-`(≤2,2)`
rationals fit by a type-`(2,2)` approximant — AAA converges in 3 support
points with no Froissart doublets, so the cleanup pass is a no-op here; it
is implemented and exercised so the oracle is honest for noisier inputs
(Rule 9).

## Verification results (pinned tolerances)

Run (one Julia process — CLAUDE.md Rule 7):

```
julia --project=../../.. external/probes/aaa-oracle/verify.jl
```

Observed agreement (commit of this probe):

| Quantity | SP.1.1 (d=1) | SP.1.2 (d=2) | Pinned tol |
|---|---|---|---|
| AAA support points used | 3 | 3 (each component) | — |
| AAA poles recovered | 2 | 2 (each component) | — |
| AAA poles vs exact known-`Q` roots | `1.1e-15` | `≤ 9.9e-15` | `1e-11` |
| per-component AAA poles agree | n/a (d=1) | `9.3e-15` | `1e-11` |
| AAA poles vs `SharedPade` shared-`Q` roots | `1.2e-15` | `9.3e-15` | `1e-9` |

For SP.1.2 the **two independent per-component AAA pole sets coincide with
each other to `9.3e-15`** and with the shared-`Q` roots to the same order —
agreement at the floating-point round-off floor, far inside AAA's `1e-13`
tolerance. AAA, fitting each component separately from sampled values on
the unit circle with no shared-denominator constraint, nonetheless places
every component's poles on the same pair of points as the roots of
`SharedPade`'s shared `Q`. That coincidence — reached by a route that
shares no computational machinery with the Taylor-jet oracles — is the
independent confirmation that the shared denominator is correct.

## References

- `references/tex/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40/p.tex`
  — §2 (barycentric representations, lines 155–230), §3 (core AAA
  algorithm, lines 315–474, including the eq. 3.8 arrowhead pencil at
  lines 459–474), §7 (Froissart-doublet removal, lines 652–661).
- `references/tex/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40/aaa_alg.m`,
  `cleanup.m` — the reference MATLAB implementation ported here.
- `docs/v0p2_pillarA_hermite_pade_findings.md:354–373` — §6 "AAA as a
  benchmark": the oracle role, the per-component-independence limitation,
  the diagnostic-not-ground-truth caveat.
- `src/SharedPade.jl` — the code being cross-checked.
