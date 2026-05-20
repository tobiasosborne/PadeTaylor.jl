# Pillar A — Hermite–Padé / Shared-Denominator Robust Padé: Literature Findings

**Session date:** 2026-05-20  
**Task:** Read-only literature survey for PadeTaylor.jl v0.2 vector-ODE extension.  
**Sources read:** GGT2013, ManoTsuda2017, NakatsukasaSeteTrefethen2018 (AAA),
GonnetPachonTrefethen2011, BeckermannLabahn1992, BeckermannLabahn1994,
LopezLagomasino2019, FidalgoLagomasino2013, and three metadata stubs.

---

## 1. GGT 2013 Baseline Recap

**Source:** `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf`, pp. 102–109.

The scalar SVD–Padé algorithm starts from the linear matching condition
`p(z) = f(z)q(z) + O(z^{m+n+1})` (GGT eq. 2.4, p. 103). Writing the
denominator coefficient vector **b** = (b₀,…,bₙ)ᵀ and the Taylor
coefficients c₀,…,c_{m+n}, the matching conditions below degree m+1 form
a Toeplitz (or generalised Toeplitz) system. Specifically, GGT eq. (2.10)
defines the n×(n+1) Toeplitz matrix C̃ whose (i,j)-entry is c_{m+i-j},
with the condition **0** = C̃**b** (p. 104). An SVD is taken:
C̃ = UΣV* (eq. 2.12, p. 104). The denominator coefficients **b** are
taken as the right singular vector corresponding to the smallest singular
value σₙ; the **b**‖₂ = 1 normalisation (eq. 2.5, p. 103) replaces the
classical b₀ = 1 normalisation, avoiding ill-conditioning when b₀ = 0
(defect-δ case). The numerator **a** is recovered from the upper rows of
the Toeplitz system (eq. 2.6 or 2.8, p. 103–104). Froissart doublets
(spurious pole–zero pairs) are suppressed because: if σₙ > 0 but is
small (noise-dominated), C̃ has effective rank ρ < n, and Algorithm 2
(p. 108) re-runs with smaller n and m until ρ = n, moving to a
lower-degree approximant in the Padé table that lies in the upper-left
of the square block (GGT §3). The relative tolerance `τ = tol·‖**c**‖₂`
(with tol = 10⁻¹⁴ by default) is the threshold for treating a singular
value as zero (p. 107–108). The QR-reweighting step appears in the
`padeapprox.m` code (GGT Figure 1, p. 109) as three lines *beyond*
Algorithm 2: `D = diag(abs(b)+sqrt(eps)); [Q,R] = qr((C̃*D).')` with
`b = D*Q(:,n+1); b = b/norm(b)`. This column-reweighted QR replaces the
final SVD null-vector extraction, improving numerical behaviour near the
tol boundary where Froissart blocks would be nearly square (GGT p. 108,
Figure 1 commentary). The reweighting is absent from the paper's
Algorithm 2 text but present in every deployed copy of `padeapprox.m`
(Chebfun `external/chebfun/padeapprox.m`, lines 278–280).

---

## 2. Shared-Denominator Algorithm: Block Matrix Construction

**Primary source:** `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`, §1.1–§2.2, pp. 4–14.  
**Supporting source:** `references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf`, §2, pp. 1–4.

### Terminology

Mano–Tsuda (following Hermite's original presentation) define two dual
problems for an L-tuple of formal power series **f** = (f₀,…,f_{L-1}):

- **Type I (Hermite–Padé):** find L polynomials Q^(i)(w) of degrees
  ≤ n−1+δᵢⱼ (one degree higher for the diagonal entry) such that
  `∑ᵢ Q^(i) fᵢ = O(w^{nL})` (ManoTsuda eq. 1.1, p. 4).
- **Type II (simultaneous Padé / shared-denominator):** find L polynomials
  P^(j)(w) of degrees ≤ n(L−1)−1 such that
  `fᵢ P^(j) − fⱼ P^(i) = O(w^{nL})` (ManoTsuda eq. 1.2, p. 4).

Type II is the "shared-denominator" structure: P^(0) plays the role of Q
(the denominator) and P^(j)/P^(0) = fⱼ/f₀ + O(w^{nL}) approximates the
ratio fⱼ/f₀ with all numerators sharing the same denominator P^(0).
For L = 2 (d = 1 component), this reduces to the ordinary scalar Padé
problem.

### The Block-Toeplitz Linear System

ManoTsuda §2.2 (pp. 12–14) gives the explicit construction for the Type
II polynomials P^(j). Setting m = n(L−1) and writing the Taylor
coefficients of fᵢ as {aⁱⱼ}_{j≥0}, ManoTsuda introduce the **rectangular
Toeplitz matrix** A^i_j(k, l) = [a^i_{j+m−n}]_{1≤m≤k, 1≤n≤l} (eq. 2.2,
p. 10). For the case j = 0 (constructing P^(0)), the matching condition
(ManoTsuda eq. 1.2 specialised) becomes a homogeneous linear system
(ManoTsuda eq. 2.6, p. 12):

```
[ A^1_m(n, m+1) ]   [b₀]   [0]
[ A^2_m(n, m+1) ] × [b₁] = [0]   (m equations per component, m = n(L−1))
[ ⋮             ]   [⋮ ]   [⋮]
[A^{L-1}_m(n,m+1)]  [bₘ]
```

This is a **block-column-stacked Toeplitz system** of size
m(L−1) × (m+1), i.e. n(L−1)² × (n(L−1)+1), for the m+1 = n(L−1)+1
unknown coefficients b₀,…,bₘ of P^(0). Each block A^i_m(n, m+1) is an
n × (m+1) Toeplitz matrix built from the coefficients of fᵢ. The full
coefficient matrix is:

```
A_full = [ A^1_m(n, m+1) ]        ← n rows, m+1 cols
         [ A^2_m(n, m+1) ]        ← n rows, m+1 cols
         [      ⋮        ]
         [A^{L-1}_m(n,m+1)]       ← n rows, m+1 cols
```

Total dimensions: **n(L−1) rows × (n(L−1)+1) columns**, so it always
has at least one nontrivial null vector (ManoTsuda p. 12, "unique up to
constants if and only if rank = m"). For general L components and
denominator degree m = n(L−1), numerator degree is also ≤ m−1, the
number of Taylor coefficients consumed is nL (m equations from each of
the L−1 non-reference components), and the total size of the linear
system is n(L−1) × (n(L−1)+1).

In the d-component notation natural for PadeTaylor (d = L−1 = number of
solution components, denominator degree m, each numerator degree ≤ m):
- **Matrix dimensions:** dm × (m+1) (before stacking)
- **d blocks,** each n×(m+1) Toeplitz, stacked vertically
- **Total system:** dm × (m+1), where n = m/d (or nearest integer)
- **Right null vector** gives denominator Q coefficients; numerators
  Pᵢ are recovered by truncated multiplication (ManoTsuda p. 13,
  formulae before Prop. 2.3)

### Determinantal Representation

ManoTsuda Remark 2.2 (p. 12) and Proposition 2.3 (p. 13) give the
**block-Toeplitz determinant** formula for P^(0):

```
P^(0)_0(w) = (1/NP^(0)) · det[ 1, w, w², …, wᵐ  ]
                                [ A^1_m(n, m+1)   ]
                                [ A^2_m(n, m+1)   ]
                                [       ⋮         ]
                                [A^{L-1}_m(n,m+1) ]
```

where NP^(0) is a normalising constant (eq. 2.7, p. 13). This is the
*block-Toeplitz determinant* Δ^(i) introduced in Remark 2.2 (p. 12):
`Δ^(i) = det[A^1_n(m,n), …, A^{i-1}_n(m,n), A^i_{n-1}(m,n), A^{i+1}_n(m,n), …, A^{L-1}_n(m,n)]`
of size m = n(L−1). The determinantal formula provides a closed-form
for the denominator but is impractical to evaluate directly; the
linear-system / SVD route is the computational path.

### Reference verdict

Mano–Tsuda is the most usable statement of the simultaneous Padé
construction *with explicit block-Toeplitz structure*. Beckermann–Labahn
1992 (`BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf`,
§2–§3, pp. 1–4) covers the same construction with the "Power Hermite
Padé" unified framework but does not write out the block matrix explicitly
in the form above; it defines PHPA via the scalar reduction trick
**F**(z) := **G**(z^s)·(1,z,…,z^{s−1})^T (eq. 5, p. 2). The
Beckermann–Labahn 1994 paper (`BeckermannLabahn1994_fast_matrix_pade_SIMAX15.pdf`)
gives the fastest *algorithm* (O(σ²) recurrence) but the structure is
the same.

---

## 3. How GGT Generalises: Scalar to Shared-Q

**Source:** ManoTsuda §2.2, pp. 12–14; BeckermannLabahn1992 §2–3, pp. 2–4.

The generalisation is direct. In the scalar case (d=1), the GGT Toeplitz
matrix C̃ is n×(n+1); its right null vector gives **b** (denominator
coefficients). In the shared-Q case (d components), the stacked
block-Toeplitz matrix is dn×(m+1) with m = n(d), i.e., denominator degree
m and consuming m+dn = m(1+1) = 2m Taylor coefficients per component.
The algorithm is:

1. Form the dn×(m+1) stacked Toeplitz matrix A_full from the Taylor
   coefficients of all d components (each component contributes an n×(m+1)
   block, as in GGT eq. 2.10 but stacked).
2. Compute SVD: A_full = UΣV*. The **right singular vector corresponding
   to the smallest singular value** gives the denominator coefficients
   **b** (= coefficients of Q). This is still one SVD of one matrix.
3. Normalise **b** so ‖**b**‖₂ = 1 (the vector analogue of GGT's
   ‖**b**‖ = 1; there is no additional per-component normalisation
   because there is only one denominator).
4. For each component i = 1,…,d, recover the numerator coefficients
   **aᵢ** by reading off the upper part of the block system (the
   "above the line" part of ManoTsuda eq. 2.6).

The key change is that the Toeplitz null-space problem scales from
n×(n+1) to dm×(m+1) (with m = denominator degree ≈ d × half-jet-length).
The SVD cost goes from O(n³) to O((dm)² × m) = O(d²m³), so the
computational cost grows with the square of the number of components.
The system remains a **single SVD of a single rectangular matrix**; no
coupled or iterated solves are needed.

### Normalisation in the vector case

The `‖**b**‖₂ = 1` condition normalises the shared denominator Q and
there is no ambiguity in the cross-component normalisation because Q is
shared. The individual numerator polynomials Pᵢ are determined up to the
same scalar multiple as Q. In exact arithmetic the scale cancels in the
ratio Pᵢ/Q; numerically, ‖**b**‖₂ = 1 keeps the denominator from
collapsing to zero (exactly as in the scalar GGT setting).

---

## 4. Froissart Suppression in the Multi-Numerator Case

**Source:** GGT2013 pp. 107–108 (Algorithm 2 §4); ManoTsuda pp. 12–14
(Remark 2.2 and uniqueness discussion); BeckermannLabahn1992 p. 4
(singular structure).

The SVD null-space tolerance argument survives in the multi-numerator
case with one important subtlety. In the scalar case a Froissart doublet
arises when the Toeplitz matrix has a small singular value not because
the function is genuinely rational at that degree, but because of
rounding noise. The robustness mechanism is: if σₘᵢₙ < τ (= tol·‖**c**‖₂),
reduce m (and n) and re-run. This works because the Padé table has square
blocks and the algorithm traverses them toward the upper-left corner
(GGT §3, Theorem 3.1).

In the shared-Q case:

- The stacked block-Toeplitz matrix A_full is dn×(m+1). It always has at
  least one right null vector (more columns than rows), so σₘᵢₙ can be
  zero even when the function is not rational. The relevant comparison
  is whether **σₘᵢₙ is small relative to σ₂** (the second-smallest
  singular value): a large gap σ₂/σₘᵢₙ ≫ 1 indicates a genuine
  one-dimensional null space (unique shared denominator up to scale).

- If multiple singular values are near zero, the null space is
  multi-dimensional. This is the vector analogue of the Padé table
  being in a square block interior (defect δ > 0). The GGT Strategy
  applies: reduce m (the denominator degree) until σₘᵢₙ is isolated.

- **Subtlety not present in scalar case:** for d ≥ 2, the block structure
  means that rank deficiency can arise even when no genuine rational
  approximant exists, because the system is underdetermined whenever
  dm ≤ m (i.e. d ≤ 1; for d ≥ 2 the system is overdetermined at
  n = m and consistent only generically). Mano–Tsuda p. 12 notes that
  the rank condition "unique up to constants if and only if rank = m"
  requires checking that the m×m square sub-determinants are nonzero;
  Beckermann–Labahn 1992 p. 4 (§4, singular-block discussion) confirms
  that the singular structure of the power Hermite–Padé table is made
  of triangular rather than square blocks in general.

- The QR-reweighting from `padeapprox.m` (column-reweighted QR
  following SVD, GGT Figure 1 p. 109) extends directly to the larger
  stacked system A_full: the reweighting is applied to the last column
  of V, which is the denominator estimate **b**, and the column
  structure of A_full is the same as in the scalar case (each column
  j of A_full corresponds to the coefficient bⱼ of Q(z) = Σ bⱼzʲ).

**Conclusion:** The tolerance–rank–reduction mechanism of GGT Algorithm 2
survives intact in the vector case; the implementation must check the
ratio σ₂/σ₁ (smallest two singular values) rather than just σ₁ alone,
and must implement the block-degree reduction analogously. The QR
reweighting step ports without change to A_full.

---

## 5. Convergence Guarantees

**Primary sources:**
- `references/hermite_pade/FidalgoLagomasino_Medina_2013_HP_meromorphic_1310.7010.pdf`,
  Theorem 1.3 and §2 (pp. 3–9).
- `references/hermite_pade/LopezLagomasino2019_intro_multiple_orthogonal_hermite_pade_1910.08548.pdf`,
  Theorem 2.1 (Nikishin perfectness) and §3 (pp. 4–8).
- Metadata stub `references/hermite_pade/AptekarevStahl1992_asymptotics_HP_polynomials_metadata.md`.

### What the theorems promise

**Fidalgo–López Lagomasino–Medina 2013, Theorem 1.3** (p. 3–4):
Let **f** = **ŝ** + **r** where **ŝ** = (ŝ_{1,1},…,ŝ_{1,m}) is a
Nikishin system of Cauchy transforms and **r** = (r₁,…,rₘ) is a vector
of rational functions whose poles lie in C \ (Δ₁ ∪ Δₘ). Let {**R_n**}
be the sequence of type II Hermite–Padé approximants (= shared-Q
approximants). Then, for j = 1,…,m:

```
lim_{n→∞, n∈Λ} P_{n,j}/Q_n = fⱼ = ŝ_{1,j} + rⱼ    inside (C \ Δ₁)'
```

where (C \ Δ₁)' denotes C \ Δ₁ minus the poles of all the rⱼ. Moreover,
each pole ζ of rⱼ of order κ *attracts exactly κ zeros of Q_n*, and
Q_n has exactly Σ(nₖ − dₖ) simple zeros inside Δ₁ (eq. 1.8, p. 4).
This is a **Stieltjes-type theorem** for the shared denominator Q: its
zeros converge to the poles of the rational perturbation **r**, and the
remaining zeros distribute along the cut Δ₁ according to the equilibrium
measure of the Nikishin system.

**López Lagomasino 2019, Theorem 2.1** (p. 4): Nikishin systems
N(σ₁,…,σₘ) with pairwise disjoint supports on the real line are *type I
and type II perfect*, meaning deg Q_**n** = |**n**| for all multi-indices
**n**. Perfectness is the key hypothesis that guarantees unique existence
(up to scale) of the shared denominator and generic non-degeneracy.

**Aptekarev–Stahl 1992** (metadata stub): The Stahl-school convergence
theory proves, for Markov/Nikishin/Angelesco systems, that Hermite–Padé
polynomials distribute their zeros according to a vector equilibrium
measure on a system of contours, and the type II approximants converge
uniformly on compact subsets of the complement of those contours. This
is the theoretical foundation for the claim "shared-Q zeros identify
the system's poles."

### Hypotheses and caveats for v0.2

1. **Nikishin structure is not guaranteed for arbitrary ODE flows.**
   The convergence theorems require the component functions (or their
   meromorphic modifications) to form a Nikishin or Angelesco system
   (supports on disjoint real intervals), which Painlevé solutions
   generically do not satisfy. The theorems still establish that the
   approximants *can* converge to the correct pole structure when the
   Taylor jet is taken from a function with poles, but for complex-plane
   ODE flows the setting is closer to the "meromorphic perturbation of
   Nikishin" class in Fidalgo–López Lagomasino–Medina, where convergence
   is in Hausdorff content (h-lim) rather than in capacity.

2. **Practical implication:** the convergence results justify the design
   choice of using the shared-Q zeros to detect movable poles, but they
   do *not* guarantee that a single SVD step with finite Taylor jet
   length will resolve all poles. The step-size heuristic (from FW 2011
   §5.2 or GGT Figure 3) must remain the primary stopping criterion; the
   convergence theorems provide the theoretical underpinning but not an
   algorithmic oracle.

3. **Conditions sufficient for PadeTaylor v0.2:**
   - The ODE components fⱼ(z) are analytic at the expansion point z₀.
   - The nearest singularity of the system is a movable pole (not an
     essential singularity or branch point), so the "rational
     modification of Nikishin" setting applies locally.
   - The denominator degree m is chosen ≥ number of expected poles in
     the step region.

---

## 6. Independent Oracle: Cabay–Jones–Labahn Calgo 766 and AAA

**Sources:**
- Metadata stub `references/hermite_pade/CabayJonesLabahn1997_algorithm_766_VECTOR_PADE_ACMTOMS23_metadata.md`.
- `references/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf`, §2–§3, pp. 1–6.
- Metadata stub `references/hermite_pade/VanIseghem1987_vector_QD_JCAM19_metadata.md`.

### Calgo 766 (Cabay–Jones–Labahn 1997)

ACM Calgo Algorithm 766, FORTRAN 90, published in ACMTOMS 23(1) 1997
(DOI 10.1145/244768.244790). The code computes **both** Padé–Hermite
(type I) and simultaneous Padé (type II = shared-Q) approximants along
the diagonal of the Padé–Hermite table, using a weakly stable algorithm
based on the Beckermann–Labahn recurrence. Its role in v0.2 is exactly
the role `padeapprox.m` plays in v0.1: a mutation-proof oracle for the
rational-approximant output. The v0.2 implementation should reproduce
Calgo 766's VECTOR_PADE output (or the Beckermann–Labahn-1994 version
at higher orders) to match accuracy on common test inputs.

**Where to obtain the FORTRAN:** The source is openly redistributable via
the ACM Calgo archive at `https://calgo.acm.org/` (directory `0766/`).
Third-party mirrors: `ACM-TOMS/CALGO` on GitHub;
`Beliavsky/Alan-Miller-Fortran/VEC_PADE.F90`. The journal paper is
paywalled (institutional access required); the code itself is open.

**Calgo 766 is not an ODE solver.** It computes approximants from a given
power series but does not implement step control, residual-based stepping,
or pole-bridging. The combination with ODE stepping is the v0.2
contribution.

### AAA as a benchmark

The AAA algorithm (Nakatsukasa–Sète–Trefethen 2018,
`NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf`) computes rational
approximations from sampled values on a discrete set Z using a
barycentric representation r(z) = n(z)/d(z), with support points and
weights selected greedily to minimise the nonlinear residual. Its
tolerance defaults to 10⁻¹³ (relative, p. 5). AAA works on *sampled
data*, not Taylor coefficients, so it is a complementary benchmark:

- **Use case as v0.2 oracle:** given that v0.2 has computed a vector
  solution trajectory along a path (evaluation data at many z points),
  AAA can independently fit a rational function to each component
  *separately* (scalar mode) and identify pole locations. If the
  shared-Q zeros from v0.2 agree with the AAA pole locations, this is
  cross-validation evidence that the shared denominator is correct.
- **Limitation:** AAA fits each component independently and therefore
  cannot enforce shared poles. Discrepancies between AAA pole locations
  for different components and the shared-Q zeros indicate either
  Froissart doublets or genuine component-specific singularities. This
  provides a diagnostic rather than a ground-truth comparison.
- AAA is freely available in Chebfun (`external/chebfun/`) and as a
  Julia port. It is already used in PadeTaylor v0.1 diagnostics.

---

## 7. Concrete Recommendation for `robust_pade` Extension

### Proposed signature

```julia
"""
shared_denominator_pade(
    jets::AbstractVector{<:AbstractVector{T}},   # d vectors of Taylor coefficients
    m::Int;                                       # denominator degree
    tol::Real = 1e-14
) -> (numerators::Vector{Vector{T}}, denominator::Vector{T})

Compute a shared-denominator (simultaneous / type-II Hermite–Padé)
approximant for a d-component vector of formal power series, each given
as a Taylor coefficient vector jets[i] = [c₀, c₁, …, c_{m+m/d}].
Returns numerator coefficient vectors P₁,…,P_d and denominator Q, all
normalised so ‖Q_coeffs‖₂ = 1.

References: ManoTsuda2017 §2.2 (eq. 2.6, p. 12); GGT2013 Algorithm 2
(p. 108); padeapprox.m lines 278–280.
"""
```

### Algorithm steps

1. **Input check.** Require `d ≥ 1`, `m ≥ 1`, `length(jets[i]) ≥ m+1`
   for all i. Throw with suggestion if any jet is too short.

2. **Form the stacked Toeplitz matrix A_full (dm × (m+1)).**
   For each component i = 1,…,d, build an n×(m+1) Toeplitz block
   (n = m rows per component, so dn = dm rows total) from jets[i].
   The (r,c) entry of block i is `jets[i][m + r - c + 1]` (zero-padded
   for indices < 0), matching ManoTsuda eq. (2.2) with k = n, l = m+1.
   Stack the d blocks vertically: `A_full = vcat(block_1, …, block_d)`.

3. **SVD of A_full.**
   ```
   U, S, V = svd(A_full)   # BigFloat SVD via GenericLinearAlgebra
   ```
   The last column of V gives the denominator coefficients **b**.
   (ADR-0002 is in force: route Matrix{Arb} through BigFloat for SVD;
   Arblib has no SVD.)

4. **Rank check and degree reduction.**
   Compute `τ = tol * norm(vcat(jets...))`. Let ρ = count(S .> τ).
   If ρ < m: reduce m to ρ, rebuild A_full, re-run. This is the
   vector analogue of GGT Algorithm 2 step 5 (p. 108).
   If σ₂/σ₁ < 10 (null space poorly isolated): warn but proceed;
   this is a heuristic threshold analogous to GGT's tol discussion.

5. **QR reweighting (port from padeapprox.m lines 278–280).**
   ```
   b = V[:, end]
   D = Diagonal(abs.(b) .+ sqrt(eps(T)))
   Q_mat, _ = qr((A_full * D)')
   b = D * Q_mat[:, end]
   b = b / norm(b)
   ```
   (Column-reweighting preserves sparsity and improves near-tol blocks;
   GGT Figure 1 commentary, p. 108.)

6. **Recover numerators.** For each component i, multiply out the
   upper rows of the block system: `a_i = upper_i * b`, where upper_i
   is the (m+1)×(m+1) upper Toeplitz block (the "above the line" part
   of ManoTsuda eq. 2.6 adapted to component i).

7. **Trailing-zero removal and normalisation.** Remove trailing small
   coefficients from **b** and each **aᵢ** (GGT Algorithm 2 steps 7–9).
   Divide all by b[1] so that Q(0) = 1 (or flag Q(0) ≈ 0 and throw with
   suggestion to use a different expansion point).

8. **Return.** Return `(numerators=[a₁,…,a_d], denominator=b)`.

### Failure modes that must throw (Rule 1)

- `length(jets[i]) < m + 1` for any i: throw with
  `detail = "jet too short for denominator degree m"`.
- After degree reduction, `ρ = 0` (all singular values below τ):
  throw with `detail = "Taylor jet is indistinguishable from zero at tolerance tol"`.
- `Q(0) = b[1] ≈ 0` (all denominator roots at 0): throw with
  `suggestion = "expand at a different base point; Q(0)=0 means a pole at the expansion center"`.
- Null-space dimension > 1 after reweighting and degree reduction:
  throw with `detail = "shared denominator is not unique; the system may not be perfect at this multi-index"`.

### BigFloat / Arb routing

The project uses `BigFloat` SVD via `GenericLinearAlgebra` (ADR-0002).
For Arb-coefficient jets:
```
A_full_big = BigFloat.(A_full)   # precision-loss caveat documented in ADR-0002
U, S, V = svd(A_full_big)
b = V[:, end]
```
Do not search for `Arblib.svd`; it does not exist (RESEARCH.md §5.1).

---

## Cross-reference summary (file paths and page numbers)

| Claim | Source | Pages |
|---|---|---|
| GGT Toeplitz matrix C̃ | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | 103–104 |
| GGT ‖b‖₂=1 normalisation | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | 103 |
| GGT Algorithm 2 (noisy/float) | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | 108 |
| padeapprox.m QR reweighting | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | 109 (Fig. 1) |
| GGT Froissart suppression mechanism | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | 107–108 |
| Mano–Tsuda Type II (simultaneous Padé) definition | `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf` | 4 (eq. 1.2) |
| Mano–Tsuda block Toeplitz construction for P^(0) | `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf` | 12–13 (eq. 2.6, Prop 2.3) |
| Mano–Tsuda block-Toeplitz determinant Δ^(i) | `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf` | 12 (Remark 2.2) |
| Mano–Tsuda rectangular Toeplitz A^i_j(k,l) | `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf` | 10 (eq. 2.2) |
| Beckermann–Labahn uniform PHPA framework | `references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf` | 1–4 |
| Beckermann–Labahn 1994 fast O(σ²) algorithm | `references/hermite_pade/BeckermannLabahn1994_fast_matrix_pade_SIMAX15.pdf` | 1–3 |
| Fidalgo–LLM Theorem 1.3 (convergence to poles) | `references/hermite_pade/FidalgoLagomasino_Medina_2013_HP_meromorphic_1310.7010.pdf` | 3–4 |
| López Lagomasino Nikishin perfectness | `references/hermite_pade/LopezLagomasino2019_intro_multiple_orthogonal_hermite_pade_1910.08548.pdf` | 4–6 |
| Aptekarev–Stahl 1992 convergence (metadata) | `references/hermite_pade/AptekarevStahl1992_asymptotics_HP_polynomials_metadata.md` | — |
| Calgo 766 oracle (metadata) | `references/hermite_pade/CabayJonesLabahn1997_algorithm_766_VECTOR_PADE_ACMTOMS23_metadata.md` | — |
| van Iseghem 1987 shared-Q algebraic object | `references/hermite_pade/VanIseghem1987_vector_QD_JCAM19_metadata.md` | — |
| AAA algorithm (scalar oracle) | `references/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf` | 1–6 |
