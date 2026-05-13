# Worklog 020 — Classical Padé via Toeplitz `\` closes the FW reproduction gap

**Date**: 2026-05-13 (late, follow-on to worklog 019)
**Author**: Claude Opus
**Scope**: Investigation of "why does our `Float64` lose 600–25,000× of
relative accuracy vs FW 2011 Table 5.1?" — triggered by user
skepticism that 2010-era MATLAB could outperform a 2026 Julia +
LAPACK + GenericLinearAlgebra impl.  **No code shipped in this
worklog**; the empirical finding is recorded, three beads are
filed, and the implementation work is the next session.

> **Take-home (the answer to user's question)**:
>
> FW 2011 did NOT use arbitrary precision.  Their MATLAB Table 5.1
> wall times (26.5 s for z=10⁴ Padé) rule that out — `vpa` arithmetic
> would have been ~100× slower.  Their numbers are pure `double`.
>
> The gap is **algorithmic**, not precision-of-arithmetic.  FW's
> §5.1.4 (lines 346-350) explicitly describes their Padé conversion
> method: **build the Toeplitz system (5.4), solve with MATLAB's `\`
> operator (LU/QR), fall back to underdetermined min-norm `\` on
> singular systems**.  No SVD anywhere.  We use **GGT 2013
> Algorithm 2** — robust Padé via SVD — which was published *two
> years after* FW 2011 (Trefethen co-authored both).
>
> GGT's SVD adds robustness against Froissart doublets and
> near-singular Toeplitz blocks.  For well-conditioned cases like
> ℘ on the equianharmonic trajectory, that robustness is wasted —
> the extra bidiagonalization + iterative diagonalization +
> tolerance-based degree reduction cost both precision *and* wall
> time.  Empirical probe:
>
> | target | method            | wall (s) | rel-err   | vs FW Table 5.1     |
> |--------|-------------------|----------|-----------|---------------------|
> | z=30   | SVD (current)     | 5.84     | 6.6e-12   | 87× worse           |
> | z=30   | classical (FW)    | 0.01     | 1.54e-13  | 2× worse            |
> | z=10⁴  | SVD (current)     | 17.76    | 6.05e-6   | 25,800× worse       |
> | z=10⁴  | classical (FW)    | 3.45     | 6.15e-11  | **3.8× BETTER**     |
>
> Classical-Padé via `Toeplitz \` is both **more accurate**
> (`580–100,000×` per-step improvement) AND **5–580× faster** at
> `Float64`.  At z=10⁴ we now beat FW's published number.

## The skepticism that triggered this

Worklog 019 had concluded "we are 600–25,000× worse than FW's
`Float64`" and filed `padetaylor-u7o` listing four speculative
attack vectors (route F64 SVD through Jacobi; Kahan-summed Horner;
compensated `h^k`; "port FW's original MATLAB Padé impl").  The
user pushed back:

> What I am highly skeptical about is that Fornberg and W had 2010-era
> computers, and *matlab*.  They were simply worse.  I seriously
> doubt they used arb prec (or did they?)  How could we be doing worse
> with our advanced impl?  The 4 theories are OK, but I am really
> wondering.

The "advanced impl" framing prompted re-reading FW 2011's §5 with
the precision-question in mind, rather than the algorithm-question.

## Re-reading FW 2011 §5 — what did they actually do?

`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`,
matching `vpa|symbolic|matlab|long.double|extended precision|chebfun|padeapprox|symbolic toolbox|digit|precision`:

  - **Line 224**: extended precision via **Maple** for the tritronquée
    BVP boundary values (`u(±20i)` of FW eq. 4.1, computed at 32
    decimal digits).  Not for the IVP path-network.

  - **Line 200**: "Computer algebra systems, like Maple and
    Mathematica, are particularly effective for solving nonlinear
    BVPs when extended precision arithmetic is required."  Confirmed:
    BVP used arb prec (Maple), IVP did not.

  - **Line 348**: they considered the ε-algorithm (Wynn / continued-
    fraction Padé), rejected as "about 10 times slower" in MATLAB.

  - **Line 350** (the smoking gun): "Since we use MATLAB, our method
    of choice for computing Padé approximations has been the
    **Toeplitz approach and solving the system (5.4) with the
    backslash operator**."

  - **Line 346**: "In rare cases, it can happen that (5.4) becomes
    outright singular.  In the present context, this has only been
    observed in cases where the Taylor coefficients are very rapidly
    decreasing in size (or a number of them become exactly zero).
    When this occurs, we remove the last row of (5.4).  Matlab's `\`
    solver will then give the solution to the underdetermined system
    that has the smallest norm, thereby restoring numerical stability."

The FW algorithm is **classical Padé via direct linear solve**, with
a row-removal fallback for the rare singular case.  No SVD.

We adopted GGT 2013 Algorithm 2 (Gonnet, Güttel & Trefethen, *Robust
Padé Approximation via SVD*, SIAM Review 55, 2013, pp. 101-117).
GGT 2013 was published **two years after FW 2011**; FW could not
have used it.  GGT's SVD path is genuinely more robust at degenerate
blocks — but the robustness machinery costs accuracy on the smooth
generic case.

## Empirical probes

Three probes (all in `/tmp/`, ephemeral; classical-Padé impl is
~50 LOC, identical to FW's eq. (5.4) recipe):

  - **Sanity check** — `(15, 15)` Padé of `exp(t)` at `t = 0.5`:
    classical 1.3e-16 rel-err vs analytical, SVD 0.0e-16 (exact in
    F64).  Both produce the same Padé to within roundoff — we are
    not breaking the math.

  - **10-step per-step trace** along a 22.5° wedge at h=0.5, F64 vs
    BF-256 (truth):

    | step | method=:svd     | method=:classical | gap factor |
    |------|------------------|---------------------|------------|
    | 1    | 3.6e-15          | 1.8e-16             | 20×        |
    | 2    | 1.2e-14          | 1.3e-15             | 9×         |
    | 3    | **6.2e-10**      | **7.1e-13**         | **870×**   |
    | 4    | 4.2e-8           | 4.8e-11             | 880×       |
    | 5    | 4.2e-7           | 4.8e-10             | 870×       |
    | 6-10 | ~10⁻⁷ – 10⁻⁶     | ~10⁻¹⁰ – 10⁻⁹       | ~1000×     |

    Step 3 is the one closest to the `z = 1` lattice pole; both methods
    take a hit there, but classical is 870× more accurate per step.

  - **Full wedge walker** at z=30 and z=10⁴ via 5-direction `:min_u`
    wedge selection (mirrors PathNetwork's Stage 1):

    | target | method     | wall (s) | rel-err   | FW Table 5.1 |
    |--------|------------|----------|-----------|--------------|
    | z=30   | :svd       | 5.84     | 6.6e-12   | 7.62e-14     |
    | z=30   | :classical | 0.01     | 1.54e-13  | 7.62e-14     |
    | z=10⁴  | :svd       | 17.76    | 6.05e-6   | 2.34e-10     |
    | z=10⁴  | :classical | 3.45     | **6.15e-11** | 2.34e-10  |

    At z=10⁴ classical beats FW by 3.8×.  At z=30 we're 2× off FW,
    but in the same ballpark (and we have BF-256 for the headline
    case anyway — worklog 008's 2.13e-14 still beats FW's
    8.34e-14).

## Why classical wins the smooth-conditioned case

GGT 2013 Algorithm 2 (`src/RobustPade.jl`, the production path):

```
1. Build (n+1)×(n+1) Toeplitz `C̃` from c_0, …, c_n.
2. SVD `C̃ = U Σ Vᵀ` (Demmel–Kahan via LAPACK at F64,
   Demmel–Veselić Jacobi via GLA at BigFloat).
3. Classify singular values: σᵢ < tol·‖c‖₂ ⇒ "noise".
4. Solve the homogeneous system using the right singular vector
   corresponding to the smallest signal singular value.
5. Trim trailing near-zero coefficients (Chebfun reweighting,
   lines 278–280 of padeapprox.m, line ports lines 134 of GGT).
```

For a `(15, 15)` Padé of a well-conditioned analytic function with
no Froissart doublets, the singular values are clean: σ_max …
σ_15 ≈ O(1), σ_16 ≈ machine zero.  The SVD finds them fine; but the
**bidiagonalization sweep + iterative `n²` × `n_sweeps` Givens
rotations** at order 30 (matrix size 16 × 16) costs ~tens of
thousands of flops, with each rotation contributing ~1 ulp of
roundoff.  Per-Padé conversion roundoff is ~`n² · eps` ≈ `5·10⁻¹³`
absolute error in the Padé coefficients.

FW's classical method:

```
1. Build the (m × m) sub-Toeplitz from c_m+1 .. c_2m
   (or equivalently c_m, c_m-1, …, c_1 in row 1).
2. Solve T · b_tail = -[c_m+1; c_m+2; …; c_2m] with LU (or QR if
   ill-conditioned, automatically by `\`).
3. Compute numerator a_k = Σ_j b_j · c_{k-j} via convolution.
```

LU on a 15×15 matrix: ~`n³/3` ≈ 1100 flops, single sweep, ~`n · eps`
≈ `3·10⁻¹⁵` absolute error per coefficient.  **No bidiagonalization,
no iterative diagonalization, no tolerance-based degree reduction.**
Faster *and* more accurate for the case where SVD's robustness is
unnecessary.

The 580–1000× per-step accuracy improvement compounds the
non-linearly-amplifying long-range walk: at z=10⁴ over 24,000
visited nodes, that compound improvement is the difference between
`6.05·10⁻⁶` rel-err (SVD) and `6.15·10⁻¹¹` rel-err (classical).

The 5–580× wall-time improvement is mostly because `robust_pade`'s
SVD path also incurs Julia overhead (the `Vector{Tuple}` allocations
in `pade_svd`, tolerance comparisons in `_trim_and_normalise`, etc.)
that LU/QR doesn't.

## When SVD is still load-bearing

Three regimes where classical `\` fails and SVD remains correct:

  1. **Froissart doublets** — spurious pole-zero pairs in the Padé
     table that GGT 2013 §3 was designed to suppress.  Classical `\`
     can produce a Padé with a numerator-zero almost on top of a
     denominator-pole, leading to large evaluation error in the
     vicinity.  Rare in our path-network walks but possible.

  2. **Singular Toeplitz** — when the Padé block boundary lines up
     with the requested `(m, n)`, the Toeplitz becomes rank-deficient.
     FW's fix (line 346): remove a row and use `\`'s underdetermined
     min-norm solution.  Equivalent to GGT's "drop the smallest
     singular vector" only for exact rank deficiency; the SVD path
     handles graded near-singular more cleanly.

  3. **Arbitrary precision** — at `BigFloat`-256, the GenericLinearAlgebra
     Jacobi SVD has *relative-accuracy guarantees* on small singular
     values (`c · 2⁻ᵖ · σᵢ`).  This is load-bearing for ADR-0002's
     precision argument — the GGT 2013 threshold `σᵢ < tol · ‖c‖₂`
     can have small-but-genuine singular values near the threshold,
     and LU would lose them.  At `Float64` the dynamic range is
     narrower and Demmel-Kahan suffices, but for BF this matters.

The right dispatch is therefore:

| element type             | default method | fallback     |
|--------------------------|----------------|--------------|
| `Float64` / `Float32`    | `:classical`   | `:svd` (LAPACK Demmel-Kahan) |
| `Complex{F64}` / `F32`   | `:classical`   | `:svd`       |
| `BigFloat`               | `:svd`         | (none — Jacobi is already robust) |
| `Arblib.Arb`             | `:svd`         | (extension; routes through BigFloat) |

Fallback trigger: `\` returns `Inf`/`NaN` or solution satisfies
`‖T · b_tail - rhs‖ > tol`.  Cheap to check.

## Beads filed

  - **`padetaylor-txg`** P1 — *Ship classical-Padé as F64 default in
    RobustPade*.  Adds `classical_pade_diagonal` to
    `src/RobustPade.jl` (~50 LOC).  Dispatch + fallback per the
    table above.  Tightens PN.2.2 z=30 F64 rtol `1e-9 → ~1e-12` and
    PN.2.3 z=10⁴ F64 rtol `5e-5 → ~1e-10`.  Mutation-prove the
    fallback path (force-singular Toeplitz).  Worklog +
    `docs/adr/0005-classical-padé-default-at-float64.md`.

  - **`padetaylor-7zw`** P2 — *BF-256 tritronquée pin: rerun FW
    Fig 4.1 step-(i) BVP at BigFloat-256*.  BVP is already generic
    in `T <: AbstractFloat` (BV.5.1 confirms).  Adds parallel BF-256
    testset to `test/fw_fig_41_test.jl`; pins tighter `u(0)` than
    the current F64 `≤3.5·10⁻¹³`.

  - **`padetaylor-u7o`** *closed (superseded by padetaylor-txg)* —
    the four speculative attack vectors are now obsolete; the actual
    diagnosis is classical-Padé.

## Pointers

  - **The smoking-gun lines in FW 2011**:
    - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:350`
      — "method of choice ... Toeplitz approach and ... backslash operator".
    - `:346` — singular fallback (remove a row, `\` min-norm).
    - `:200`, `:224` — Maple's role (extended-precision BVP only).

  - **GGT 2013 SVD reference**:
    - `references/markdown/GGT2013_robust_pade_via_SVD_SIAMRev55/GGT2013_robust_pade_via_SVD_SIAMRev55.md`
      §3 (Froissart doublets) — why GGT introduced SVD.

  - **In-tree**:
    - `src/RobustPade.jl` — current SVD-based `robust_pade`.
    - `src/LinAlg.jl` — current `pade_svd` dispatch.
    - `docs/worklog/019-fw-table-51-z10000.md` — the predecessor
      worklog that posed the F64 gap as "open question".

## Hard-won lessons (for HANDOFF.md §"Hard-won")

**35. When porting from a paper, port the actual algorithm, not the
"modern equivalent"**.  We adopted GGT 2013 Algorithm 2 from
Chebfun's `padeapprox.m` because it's the modern recipe; FW used
classical Toeplitz `\` because it predates GGT 2013.  Two years of
algorithmic progress moved the recipe in a direction (robustness via
SVD) that costs accuracy on the well-conditioned case.  When you're
trying to reproduce a paper's published numbers, port the paper's
algorithm.

**36. "Modern equivalent" sometimes means "added work the paper
didn't need"**.  GGT 2013's robustness is real (Froissart doublets,
near-singular blocks).  But on the Painlevé pole-field problem, the
robustness machinery is unused — the singular values are clean.
The bidiagonalization + iterative diagonalization + tolerance-based
degree reduction cost ~580× per-step accuracy AND ~5-580× wall time
relative to the classical `\` for the smooth case.  Robustness has
a price; spend it only where you need it.

**37. Read the paper's wall times to constrain its arithmetic
precision**.  FW Table 5.1's 26.5 s for z=10⁴ Padé rules out `vpa` —
MATLAB's arbitrary-precision arithmetic is ~100× slower than
`double`.  Their numbers MUST be `Float64`.  When a paper claims a
specific rel-err number, the only-arithmetic-precision-that-fits-the-
wall-time is often diagnostic of the underlying impl choices.

**38. The right precision dispatch is element-type-driven, not
"always use the modern algorithm"**.  At `Float64` the dynamic
range is narrow and Demmel-Kahan SVD (or LU) suffices.  At
`BigFloat` the dynamic range opens up and relative-accuracy Jacobi
SVD is the load-bearing tool.  Dispatch by `T`, default to the
right method for each tier; expose the alternative as an opt-in
for advanced users.
