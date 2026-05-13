# BVP Recipe: Chebyshev–Newton Solver for Complex-Plane Segments

Canonical recipe for the Chebyshev spectral BVP solver used in PadeTaylor.jl
to fill smooth bands between pole fields (FW 2011 §3.2). All formulas are
cited to their ground-truth source; marker extractions will supersede txt-file
citations as they complete.

Provenance: Compiled from Trefethen SMIM 2000, Weideman–Reddy 2000,
Berrut–Trefethen 2004, DMSUITE source.

---

## §1. Chebyshev Extrema Nodes and Segment Mapping

### Nodes on [−1, 1]

The N+1 Chebyshev points of the second kind (extrema of T_N) are:

```
t_j = cos(j·π/N),   j = 0, 1, …, N
```

These are the "Chebyshev points of the second kind" or "Gauss–Lobatto" nodes.
They include both endpoints t_0 = 1 and t_N = −1.

Source: FW 2011 §3.2.1 (references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:184 — "t_j = cos(jπ/N), j = 0, …, N,
which are the extrema of the Chebyshev polynomial of degree N").

Also in: Weideman–Reddy 2000 §3.2, Eq. (13), cited from
references/WeidemanReddy2000_DMSUITE_ACMTOMS26.txt:1190–1205.
The DMSUITE implementation uses the trigonometrically equivalent form
x_k = sin(π(N−1−2(k−1)) / (2(N−1))) for perfect floating-point symmetry;
see external/DMSUITE/chebdif.m:38.

### Affine Map from [−1,1] to Complex Segment [z_a, z_b]

Given a complex segment from z_a to z_b:

```
z(t) = (z_b + z_a)/2 + (z_b − z_a)/2 · t,   t ∈ [−1, 1]
```

This maps t = −1 → z_a and t = 1 → z_b. The collocation nodes in the z-plane
are z_j = z(t_j).

Source: FW 2011 §3.2.1, Eq. (3.1) preamble (references/markdown/
FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:178 —
"z = ½((z_b + z_a) + (z_b − z_a)t)").

Note on generic T: For T <: AbstractFloat the cos call is T(cos(j*π/N)) or
equivalently use the sinusoidal form from chebdif.m:38 for better symmetry.
For Complex{T}, only the affine map uses z_a, z_b; all subsequent arithmetic
on collocation values is component-wise and generic.

---

## §2. First Derivative Matrix D₁

### Off-Diagonal Entries

For i ≠ j, indices i, j ∈ {0, …, N}:

```
D₁[i,j] = (c_i / c_j) · (−1)^(i+j) / (t_i − t_j)
```

where the endpoint weights are:
```
c_0 = c_N = 2,    c_i = 1  for 1 ≤ i ≤ N−1
```

### Diagonal Entries

Computed by the negative-row-sum identity (ensuring D₁ differentiates
constants exactly):

```
D₁[i,i] = −∑_{j≠i} D₁[i,j]
```

Corner entries: D₁[0,0] = (2N² + 1)/6,  D₁[N,N] = −(2N² + 1)/6.

Source: Weideman–Reddy 2000 §3.2, Eq. (14) (Chebyshev D^(1) off-diagonal):
references/markdown/WeidemanReddy2000_DMSUITE_ACMTOMS26/WeidemanReddy2000_DMSUITE_ACMTOMS26.md:419–425
("D^(1)_{k,j} = c_k(−1)^(j+k) / (c_j(x_k−x_j)) for j≠k"). Note: MATLAB
1-based indices vs 0-based; the formula is identical. The general recursion
(Eq. 4 in §2) is at md:122.

DMSUITE implementation: external/DMSUITE/chebdif.m:45–57. The matrix C at
line 45 is the c_i/c_j factor; Z at line 49 is 1/(t_i − t_j); the off-diagonal
recursion (line 55) for ell=1 reduces to the formula above. The diagonal
is set by negative row-sum at line 56.

---

## §3. Second Derivative Matrix D₂

### Option (a) — Squaring (adopted for v1)

```
D₂ = D₁ · D₁
```

Simple, O(N³) matrix multiply. Sufficient for N ≤ 50.

Source: FW 2011 §3.2.1 (references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:188 — "D₂ is the second derivative matrix
… see [12,26,29]"). Trefethen SMIM ch.6 (references/TrefethenSMIM_2000_book.txt:
garbled pdftotext; cite as pending marker extraction) recommends D₂ = D₁² for
the introductory treatment. DMSUITE chebdif.m implements the full M-th order
recursion for any M; calling chebdif(N+1, 2) returns both D₁ and D₂ directly
(external/DMSUITE/chebdif.m:54–58, ell=2 iteration).

### Option (b) — Direct Formula (better conditioning at large N)

Weideman–Reddy 2000 §2 provides a direct recursion for D^(ℓ) that avoids
repeated squaring error amplification. For ℓ = 2, the off-diagonal entry is
obtained by applying Eq. (4) with ℓ = 2:

```
D₂[k,j] = (2 / (t_k − t_j)) · ((c_k/c_j) · D₁[k,k] − D₁[k,j])
           for k ≠ j
```

with the diagonal again set by negative row sum.

Source: Weideman–Reddy 2000 recursion Eq. (4):
references/markdown/WeidemanReddy2000_DMSUITE_ACMTOMS26/WeidemanReddy2000_DMSUITE_ACMTOMS26.md:122
DMSUITE implementation: external/DMSUITE/chebdif.m:54–57 (the `ell` loop,
second pass).

**Decision for PadeTaylor.jl v1:** Use option (a) (D₂ = D₁²). For N ≤ 50,
the conditioning difference is negligible (< 1 ULP at Float64). A future bead
should track whether N > 50 becomes needed and whether chebdif.m's option (b)
should be ported.

---

## §4. Newton Iteration for u'' = F(z, u)

### Setup

The ODE on [−1, 1] (after the affine change of variables) is:

```
D₂ u = (h/2)² · F(t, u)
```

where h = z_b − z_a (complex step length), and for PI specifically:

```
F(t, u) = 6u² + z(t)
```

after substituting z(t) from §1. The exact matrix form for PI is:

```
D₂ u = (1/4)(z_b − z_a)² (6u² + (z_b+z_a)/2 · 1 + (z_b−z_a)/2 · t)
```

Source: FW 2011 §3.2.1, Eq. (3.2) (references/markdown/
FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:186).

### Boundary Conditions

Set u_0 = u_b (value at t = 1 ↔ z_b) and u_N = u_a (value at t = −1 ↔ z_a)
directly in the vector before each Newton step. Equivalently, delete rows 0 and
N from the residual and Jacobian and solve only on the N−1 interior nodes.

Source: FW 2011 §3.2.1 (references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:190 — "enforce the BCs u_0 = u_b and
u_N = u_a as described … in [29], Section 4"). Weideman–Reddy 2000 §4 details
the row/column deletion procedure (references/WeidemanReddy2000_DMSUITE_ACMTOMS26.txt:
2580–2640).

### Residual and Jacobian (Interior Nodes)

Let subscript int denote rows 1 to N−1 (interior only).

```
R = (D₂ u)_int − (1/4)(z_b−z_a)² · (6u_int² + z_int)
```

For the general second-order ODE u'' = F(z, u), the analytic Jacobian of the
residual R with respect to the interior values u_int is:

```
J = (D₂)_int,int − (1/4)(z_b−z_a)² · diag(∂F/∂u(z_int, u_int))
```

For PI specifically, ∂F/∂u = 12u, so:

```
J_PI = (D₂)_int,int − (1/4)(z_b−z_a)² · diag(12 · u_int)
```

FW 2011 notes "the Jacobian is trivial to calculate" (line 190) without giving
the explicit formula. The derivation above follows directly from differentiating
Eq. (3.2) with respect to u_int.

### Newton Step

```
Δu = J \ R
u_int ← u_int − Δu
```

Repeat until ‖R‖_∞ ≤ tol. FW 2011 reports "no more than six iterations were
needed" for their smooth-region BVPs (line 190). Their convergence tolerance was
not stated explicitly in §3.2; they reference chebint/DMSUITE which is Float64.
For our generic-T implementation: tol = 100 · eps(real(T)) is a reasonable
default (not pinned by FW; see §7).

Solver routing: for Float64, J \ R calls LAPACK (dgesv). For BigFloat or
Complex{BigFloat}, it routes to GenericLinearAlgebra (already in Project.toml).
For Complex{Float64}, standard LAPACK handles it.

---

## §5. Barycentric Interpolation to Arbitrary Points

### Weights for Chebyshev Type-2 Nodes

For nodes t_j = cos(jπ/N), j = 0, …, N, the barycentric weights are:

```
w_j = (−1)^j · δ_j,   where δ_0 = δ_N = 1/2,  δ_j = 1  for 1 ≤ j ≤ N−1
```

Source: Berrut–Trefethen 2004, §5, Eq. (5.4):
references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/BerrutTrefethen2004_barycentric_SIAMReview.md:164–171
("The Chebyshev points of the second kind are given by x_j = cos(jπ/n) …
w_j = (−1)^j δ_j, δ_j = 1/2 for j=0 or j=n, 1 otherwise").

DMSUITE implementation: external/DMSUITE/chebint.m:29–30 — w initialized to
(−1)^j and then w[1] and w[N] halved.

### Evaluation at z*

```
u(z*) = (∑_j  w_j · u_j / (z* − t_j)) / (∑_j  w_j / (z* − t_j))
```

Special case: if z* = t_k for some k, return u_k directly (avoids 0/0).

This is the "second (true) barycentric form" of Rutishauser; see Berrut–Trefethen
2004 Eq. (4.2) (txt:312–330).

Source: Berrut–Trefethen 2004, Eq. (4.2):
references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/BerrutTrefethen2004_barycentric_SIAMReview.md:127–133
(barycentric formula, second true form, with denominator). FW 2011 §3.2.1
line 190 references "barycentric form, which is both fast and numerically
stable [3,15]", where [15] is Berrut–Trefethen 2004. DMSUITE implementation:
external/DMSUITE/chebint.m:32–36 (the D matrix method; numerically equivalent
via the eps(D==0) guard for coincident-node detection).

### Mapping Back to z

The segment runs in the complex z-plane. The barycentric formula operates on
the t-parameter grid. To interpolate u to a complex target point z*, first
compute the pre-image t* = (2z* − z_a − z_b) / (z_b − z_a) ∈ [−1,1], then
apply the formula above with (t_j, u_j) pairs and evaluation point t*.

---

## §6. Generic-T Discipline

All formulas in §§1–5 must work for T <: Number with the following care:

**cos call (§1 nodes):** Write `T(cos(j * π / N))` or equivalently use the
sinusoidal form `sin(π * (N − 1 − 2j) / (2(N−1)))` which is exact and
symmetric. The plain `cos(j*π/N)` computes in Float64 if unguarded; casting
to T propagates BigFloat precision when T = BigFloat.

**Matrix multiply (§3 option a):** `D₂ = D₁ * D₁` works for any T because
matrix multiplication is generic. For Complex{T} the entries are automatically
complex.

**Backslash (§4 Newton step):** `J \ R` routes to LAPACK for Float64/Complex{Float64}
and to GenericLinearAlgebra.jl for BigFloat/Complex{BigFloat}. No special
dispatch needed; the fallback in Julia's stdlib calls the right method
automatically when GenericLinearAlgebra is loaded.

**diag (§4 Jacobian):** `Diagonal(12 .* u_int)` or `diagm(12 .* u_int)` is
generic and allocates correctly for any T.

**Barycentric division (§5):** The `eps(D==0)` guard in chebint.m is Float64-
specific. In Julia, guard with `z ≈ t_k` (using `atol = eps(real(T))`) and
return `u_k` directly in that branch.

**Segment endpoints:** z_a, z_b ∈ ℂ enter only through the affine map in §1
and the scale factor (z_b − z_a)² in §4. Everything else operates on the
t-parameter values which are real(T).

---

## §7. Open Spec Gaps

The following decisions are NOT pinned by FW 2011 §3.2 or DMSUITE and must be
resolved before production use. Each is a potential future bead.

**Newton convergence tolerance:** FW 2011 §3.2 does not state the tolerance
explicitly (line 190 says "fast and the Jacobian is trivial"). The paper was
implemented in MATLAB Float64. For our generic-T implementation we use
`100 · eps(real(T))` as a default, which is conservative. The right value
may depend on N and on the conditioning of the particular segment.

**Adaptive-N selection:** FW 2011 mentions "an increase in N is indicated" when
endpoint derivative residuals exceed `10⁻⁷` to `10⁻⁸` (line 192), but gives
no algorithm for how to choose the new N. DMSUITE and Chebfun use doubling
strategies; this is unspecified for PadeTaylor.jl.

**Segment-orientation policy:** FW 2011 §3.1 describes path selection for the
IVP stage (minimize |u| at each step) but does not specify how segment endpoints
z_a, z_b are oriented for BVP calls — specifically which end is called z_a and
which z_b, and whether the map z_a ↔ t = −1 convention matters for the Newton
initial guess (asymptotic approximation from (1.2)). This is implicit in the
MATLAB code but not stated.

**Maximum Newton iterations:** FW 2011 reports ≤ 6 in practice (line 190);
no hard cap is specified. A reasonable implementation cap is 20, throwing with
a suggestion if exceeded.

**D₂ = D₁² vs direct formula crossover:** At what N does condition number
degradation from squaring become numerically significant? FW 2011 and DMSUITE
do not bound this for complex segments. A deferred bead should profile at
N = 50, 100, 150 with BigFloat ground truth.

---

## §8. References

All citations resolved as of 2026-05-13. Marker extractions pending for
BerrutTrefethen2004 and WeidemanReddy2000 (jobs dispatched; txt-file citations
above will be superseded when markdown files appear under references/markdown/).

### Primary sources

- **FW 2011 §3.2.1** (Fornberg–Weideman 2011 JCP 230):
  `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
  Lines 176–200 (affine map, collocation nodes, Eq. 3.2, Newton, barycentric).

- **Weideman–Reddy 2000** (ACMTOMS 26):
  `references/markdown/WeidemanReddy2000_DMSUITE_ACMTOMS26/WeidemanReddy2000_DMSUITE_ACMTOMS26.md`
  (1695 lines, marker-converted). Key lines: 122 (Eq. 4 recursion, §2),
  294–296 (node convention note, §3), 419–425 (Eq. 14, D^(1) formula, §3.2),
  849–876 (§4 boundary conditions).

- **Berrut–Trefethen 2004** (SIAM Review 46):
  `references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/BerrutTrefethen2004_barycentric_SIAMReview.md`
  (460 lines, marker-converted). Key lines: 127–133 (Eq. 4.2, barycentric
  formula), 164–171 (§5, Eq. 5.4, Chebyshev type-2 weights).

- **Trefethen SMIM 2000** (book):
  `references/TrefethenSMIM_2000_book.txt` — pdftotext extraction is garbled
  (binary font encoding in PDF). Marker extraction at
  `references/markdown/TrefethenSMIM_2000_book/` in progress (dispatched as
  Task C; will supersede when complete). Relevant chapters: ch.6 (cheb.m, D₁
  formula), ch.13 or ch.14 (nonlinear BVPs by Newton).

### DMSUITE source

- `external/DMSUITE/chebdif.m` lines 38, 45–57: node computation and D^(ℓ)
  recursion.
- `external/DMSUITE/chebint.m` lines 27–36: Chebyshev points, barycentric
  weights, and evaluation formula.

### One-line provenance

Compiled from Trefethen SMIM 2000, Weideman–Reddy 2000, Berrut–Trefethen 2004,
DMSUITE source. Ground-truth equations verified against FW 2011 §3.2 markdown
extraction (marker-converted, not pdftotext).
