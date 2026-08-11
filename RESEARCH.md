# RESEARCH.md — Stage 0 deep research for PadeTaylor

> **Status: Stage 0 deliverable complete.** §1.1 and §2 are deeply
> read with cited claims; §3 is package-audit-derived; §4 is verified-
> arithmetic survey; §5 is bigfloat-SVD landscape; §6 is source-code
> hunt summary; §7 resolves all PRD open questions; §8 lists new
> questions surfaced; §9 is the go-to-Stage-1 recommendation.
>
> Per PRD's "ground truth via local PDFs" discipline, every algorithmic
> claim carries a `references/<file>.pdf` citation (line numbers refer
> to the marker-converted markdown files under `references/markdown/`).
>
> Outstanding: §1.2, §1.4 deep reads (deferred to verification phase —
> they confirm the FW 2011 algorithm on PII/PIV but add no new
> algorithmic content); §1.5 Willers 1974 (not yet acquired; ACM DL
> Cloudflare-blocked).

## 0. Scope of this document

Stage 0 of the project plan in `PRD.md`. Output of this document
informs `DESIGN.md` (Stage 1) and gates the move to implementation.
Resolves or sharpens every provisional claim in `PRD.md` against the
literature and the ecosystem. Ends with a **go/no-go recommendation**
for Stage 1.

## 1. Algorithmic core — the Fornberg–Weideman line

### 1.1 FW 2011 (J. Comput. Phys. 230, 5957–5973)

Reference: `references/FW2011_painleve_methodology_JCP230.pdf` →
markdown at `references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md`. Line citations below are to
that markdown file.

#### 1.1.1 Algorithmic spec — four independently designable layers

The "Pole Field Solver" (PFS, the FFW 2017 name) decomposes into:

1. **Taylor coefficient generation.** At each step, generate
   `c_0, c_1, …, c_n` of the local Taylor expansion `y(t+h) = Σ c_k h^k`
   `[FW2011: §2.1, lines 79–107]`. Two methods presented: (a) substitute
   the truncated expansion into `y'(t+h) = f(t+h, y(t+h))` and equate
   coefficients (easy "by hand" for PI/PII); (b) bootstrap — start
   from one term, substitute, integrate, gain one correct coefficient
   per pass; "always implementable entirely numerically"
   `[FW2011: lines 105–107; ref [2] = Barton-Willers-Zahar 1971]`.
   Method (a) is symbolic-AD; (b) is operator-overloading / recurrence.
   PadeTaylor.jl needs (b) for general analytic `f` — and §3.3 of this
   document confirms `TaylorSeries.jl::Taylor1{T}` provides exactly
   that, including for `T = Arb`.
2. **Padé conversion.** Convert truncated Taylor expansion to **diagonal
   `(v, v)` Padé** with `n = 2v` (equal numerator and denominator
   degrees) `[FW2011: §2.2, lines 109–116, eq. 2.5]`. *"A pole in
   `y(t+h)` located near `t` will just introduce a zero in the
   denominator of (2.5)"* — this is **the entire mechanism by which
   the integrator survives pole-crossings**.
3. **Step direction selection.** At each step, evaluate the Padé
   approximation in **five candidate directions** (straight at goal,
   ±22.5°, ±45°) and pick the one minimising `|y(t+h direction)|`
   `[FW2011: §3.1, lines 158–164]`. The Padé is built once and
   re-evaluated cheaply. Path is a **branching tree** — new target
   points start from the nearest previously-visited point, reusing its
   Padé data `[FW2011: lines 162–164]`.
4. **Smooth-region BVP solver.** Where RHS terms balance and the IVP
   becomes ill-conditioned, switch to Chebyshev spectral collocation
   BVP `[FW2011: §3.2.1, lines 174–194]` with Newton iteration (≤6
   iters; trivial Jacobian). Pole-field IVP supplies BVP boundary `u'`;
   match-tolerance `10⁻⁷–10⁻⁸` is the error diagnostic
   `[FW2011: lines 191–192]`.

#### 1.1.2 Standard parameters — `(order, h) = (30, 0.5)`

`(15, 15)` diagonal Padé at step length `h = 0.5`
`[FW2011: §5.1, lines 277–280, 326]`. Justified empirically against
**Weierstrass ℘** as a closed-form companion to PI (drop the `z` term
from `u'' = 6u² + z`; ℘ has the same pole structure but is elliptic, so
reference values are known to arbitrary precision). Contour plots of
(accuracy × CPU time) over `(order, h)` show a *generally favourable*
basin around `(30, 0.5)`. Worth re-running this experiment for
PadeTaylor.jl on a wider class of `f` to find whether the optimum is
PI-specific or universal.

#### 1.1.3 Pole-field-edge detection — out of scope for v1

Away from poles, `u(x, y)` satisfies Laplace's equation, so the
discrete-Laplacian residual is `O(h²)` with a tiny prefactor in smooth
regions and *much* larger in pole fields `[FW2011: §3.2.2, lines
202–208, eq. (3.3)]`. Contour at level `0.001` separates them.
Out-of-scope for v1 of PadeTaylor.jl per PRD ("real-axis only with
`step_complex` API" for v1); revisit at v2 alongside the BVP coupling.

#### 1.1.4 Step-direction selection — algorithmic options

Two equivalent routes `[FW2011: §5.3 + §5.4.1, lines 352–368]`:

- **(a) Five-direction sampling.** Default; cost is essentially free
  because the Padé is built once.
- **(b) Steepest-descent on `|u|`.** Analytic direction
  `θ = arg(−u(z₀) / u'(z₀))`. If inside the 90°-wedge towards the
  target, accept; otherwise snap to the nearer wedge edge. *"Probably
  preferred for Fortran/C implementations."* Marked `RKN12†` in their
  Table 5.1 — equally accurate but runs marginally faster.

For PadeTaylor.jl, **(b) is cleaner** — no quintupled work, cleaner
derivative math, one geometric primitive. Adopt (b) as default; expose
(a) as a fallback flag.

#### 1.1.5 The Padé routine FW actually use — not GGT 2013!

**FW use direct Toeplitz solve, not SVD-based robust Padé**
`[FW2011: §5.2, lines 328–350]`:

- Solve `T b = −c` where `T` is `v × v` Toeplitz built from
  `c_v, …, c_1`, then back-substitute for numerator coefficients via
  triangular Toeplitz `[FW2011: eqs. (5.4)–(5.5)]`.
- "It is well known that the condition number of [the Toeplitz
  system] can grow rapidly with increasing `v`" `[ref [1] = Baker &
  Graves-Morris]` — **but** "one typically obtains a rational function
  that is far more accurate than the coefficients of the denominator
  polynomial" `[FW2011: lines 343–345]`. The output `(P, Q)` is
  well-conditioned even when `T` itself is not.
- Outright-singular `T` — observed only when Taylor coefficients
  rapidly decrease or several are exactly zero — handled by removing
  the last row of (5.4); MATLAB `\` returns the min-norm solution.
- ε-algorithm / η-algorithm rejected: 10× slower in MATLAB;
  quotient-difference instabilities not bypassed.

**Implication for our design.** GGT 2013 (PRD's named canonical Padé
routine) **departs from FW's recipe**. GGT eliminates Froissart
doublets and yields *pointwise* convergence as `(L, M) → ∞`; FW's
direct Toeplitz is fragile in pathological cases but FW's tests don't
encounter those cases catastrophically. Open question for §2: **does
GGT actually do better than direct Toeplitz on FW's specific
near-singular cases** (rapidly-decaying coefs, several exact zeros)?
A/B test once GGT is read and Chebfun's `padeapprox.m` is local
(both now done — see §2 + Appendix A).

#### 1.1.6 Boundaries the paper documents

- **Cancellation post-pole.** "Some accuracy is lost because, following
  a pole, small values have been generated by cancellation of larger
  ones — a well known recipe for losing significant digits"
  `[FW2011: line 135]`. Mitigated by **path choice**: stay in low-`|u|`
  passages between poles rather than going through them.
- **Smooth regions are an IVP failure mode**, not just "less interesting"
  `[FW2011: lines 170–171, 397]`. RHS terms balance; tiny `u` changes
  cause big `u''` changes; IVP ill-conditioned. Handled by BVP switch,
  not by integrator parameter tweaks.
- **Real axis is "different"** in the right half-plane — RHS terms cannot
  approximately cancel along the positive real axis, so smooth bands
  cannot be placed there `[FW2011: line 251]`. Specific to PI;
  structurally the same observation generalises to any analytic ODE
  whose RHS does not balance there.

#### 1.1.7 Performance numbers for cross-validation

- 161×161 grid (~26 000 evaluation points) computed in **~0.75 s**
  MATLAB notebook `[FW2011: line 139]`. Stage 1 (~990 path steps at
  `h = 0.5`) ≈ 0.41 s, Stage 2 (vectorised final steps) ≈ 0.34 s
  `[FW2011: lines 164, 166]`.
- **Weierstrass ℘ at `z = 30`**: relative error `7.62 × 10⁻¹⁴` (Padé)
  in 0.020 s. **At `z = 10⁴`**: `2.34 × 10⁻¹⁰` in 26.5 s. **At
  `z = 28.261`** (high up on a pole wall): `7.92 × 10⁻¹⁰` in 0.018 s
  `[FW2011: Table 5.1, lines 385–391]`. Reference values to 16 digits
  in `[FW2011: line 295, line 372]`. **These are direct
  cross-validation targets for PadeTaylor.jl** — Float64 should match
  to within a small constant; arb-prec at, say, 100 dps should hit
  precisely the true value.

#### 1.1.8 Path-network self-consistency check — answers PRD open Q4

FW propose a **stochastic accuracy diagnostic** `[FW2011: lines
305–311]`:

- Re-run the same pole field with the same ICs but a different random
  path order; if the two runs disagree at a node, that's a numerical-
  error indicator.
- For real ICs, the conjugation symmetry upper-half ↔ lower-half plane
  is another diagnostic.

**Both are free** — no extra runs of the algorithm with elevated
precision needed. Adopt both as the path-network consistency metric
(PRD open question 4).

### 1.2 FW 2014 (FoCM 14, 985–1016) — PII survey

Reference: `references/FW2014_second_PII_exploration_FoCM14.pdf`. DOI
`10.1007/s10208-013-9156-x` (the PRD's `s10208-013-9182-8` was wrong;
correct DOI confirmed via CrossRef during acquisition).

**Status: deep read deferred to verification phase.** This paper applies
the FW 2011 PFS to PII; the *algorithm* is FW 2011's (already digested
in §1.1). What this paper adds is target *figures* of PII solution
space the v1 of PadeTaylor.jl must reproduce qualitatively (PRD
"Definition of done" item 4). Re-read when verifying — at that point
extract and pin specific reference values for cross-validation.

### 1.3 Fasondini–Fornberg–Weideman 2017 (J. Comput. Phys. 344, 36–50)

Reference: `references/FFW2017_painleve_riemann_surfaces_preprint.pdf`
(preprint version; the published JCP version requires Cloudflare-bypass
and we don't need the differences for Stage 0).

**Riemann-surface extension.** PIII / PV / PVI have *multi-valued*
transcendents on Riemann surfaces; FFW 2017 extends the PFS to compute
these on multiple sheets. Mechanism (from Page 1 of the preprint, our
§1.1.1 baseline plus a sheet-tracking layer): record path-traversal
data so that monodromy from looping around a branch point is detected
and the integration continues on the correct sheet.

**v1 vs v2 implications.** The PRD scopes Riemann-surface support to
v2: "v1 likely real-axis only" + "step_complex API". v1 of PadeTaylor.jl
should design `step_complex` so the path data structure can later be
extended to track sheets — the same `Path` abstraction Mezzarobba uses
(see §4.1) inherits cleanly here.

### 1.4 Reeger–Fornberg 2014 (Physica D 280–281, 1–13) — PIV survey

Reference: `references/ReegerFornberg2014_PIV_fundamental_domain_
PhysicaD280.pdf`. PIV applies FW 2011's PFS to the four-parameter PIV
solution space.

**Status: deep read deferred to verification phase.** Same role as
FW 2014 (§1.2): provides target figures for cross-validation when
PadeTaylor.jl reaches the qualitative-reproduction milestone. PIV is
algorithmically identical to PI/PII at the PFS level — no new
algorithmic content beyond FW 2011.

### 1.5 Willers 1974 — original IVP-with-poles

**Not yet acquired.** ACM DL Cloudflare-blocked; user-fetch link in
Appendix A. Historical-baseline reading; the algorithm has been
superseded by FW 2011 on every dimension. **Low priority.** Re-evaluate
during writeup phase to ensure we credit appropriately and to note
the continued-fraction lineage that is the precursor to Padé in this
context.

## 2. Robust Padé approximation

### 2.1 Gonnet–Güttel–Trefethen 2013 (SIAM Review 55)

Reference: `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` →
markdown at `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
GGT2013_robust_pade_via_SVD_SIREV55.md`. Line citations are to the
markdown file.

#### 2.1.1 The mathematical setup — three normalisation choices

For Padé type `(m, n)` the linearised condition is `p(z) = f(z) q(z) +
O(z^{m+n+1})` `[GGT2013: eq. (2.4)]`. With `q ≢ 0`, write
coefficients as vectors `a ∈ ℝ^{m+1}, b ∈ ℝ^{n+1}` and the condition
becomes a linear system whose lower block requires `b` to be a (right)
null vector of an **`n × (n+1)` Toeplitz matrix `C̃`**
`[GGT2013: eq. (2.10), lines 110–117]`. Three standard choices for
making this system non-vacuous:

- **`b_0 = 1`** — used by FW 2011 §5.2; reduces to a square `n × n`
  linear system on `b_1, …, b_n`; "highly ill-conditioned or singular"
  in pathological cases `[GGT2013: line 65]`.
- **`||b||_2 = 1`** — GGT 2013's choice (and Gonnet-Pachón-Trefethen
  2011's choice for rational interpolation on the unit circle). The
  null vector of the rectangular `C̃` always exists (it has more
  columns than rows). **This is the regularisation that buys
  robustness.**
- The third route — work with the *square* matrix `C` (delete the
  first column of `C̃`) and bypass it when singular — is "many
  treatments of Padé approximation" `[GGT2013: line 124]`.

The `||b||_2 = 1` normalisation is what makes Algorithm 1's
upper-left-corner-of-block hopping well-defined.

#### 2.1.2 Algorithm 1 (exact arithmetic) — the corner-hop loop

`[GGT2013: lines 140–155]`:

1. If `c_0 = … = c_m = 0`: `p = 0, q = 1`, stop.
2. If `n = 0`: trivial.
3. Compute SVD `C̃ = U Σ V*`. Let `ρ` = number of nonzero singular
   values.
4. If `ρ < n`: reduce `n → ρ`, `m → m − (n − ρ)`, restart.
5. `b` = last column of `V` (the null right singular vector); `a` =
   from the upper-block multiplication.
6. If `b_0 = … = b_{λ−1} = 0`: cancel common `z^λ` factor in `p, q`.
7. Divide by `b_0` so `b_0 = 1`.
8. Trim trailing zeros.

**Theorem 2.1**: terminates in `≤ 2 + log₂(δ + 1)` SVD computations,
where `δ` is the defect of `r_{mn}` `[GGT2013: line 157]`. In practice,
on the diagonal `(n, n)` Padé table: 99% of `(m, n)` pairs hit `ρ = n`
on the first try and stop at step 5. The hopping happens only inside
the lower-right halves of degeneracy blocks.

#### 2.1.3 Algorithm 2 (Algorithm 1 + tol thresholding, for noisy data)

`[GGT2013: lines 219–235]`. Three modifications to Algorithm 1:

1. **Singular-value threshold.** Treat `σ_i` as zero if
   `σ_i < τ = tol · ||c||_2`. **Default `tol = 10⁻¹⁴`** for Float64
   rounding-error problems.
2. **Coefficient threshold** in steps 1, 6, 8 — same `τ` (or `tol`)
   determines whether a `c_k` or `b_k` counts as zero.
3. **Optional rescaling** `f̂(z) = f(z/γ)` for some `γ > 0` to keep
   `c_0, c_1, …` of comparable magnitudes. Algorithm 1's 2-norm is
   "implicitly tied to unit-disk scaling" `[GGT2013: lines 215–217]`;
   without rescaling, geometrically-decaying or geometrically-growing
   coefficient sequences fight the algorithm. Auto-`γ` selection is
   in Fornberg 1981 `[GGT2013: ref [9]]` — same Bengt Fornberg as
   FW 2011. Fits cleanly into PadeTaylor's pipeline.

For PadeTaylor.jl at arb-prec: **set `tol ≈ 2^{−p + ε}`** where `p` is
the working precision in bits, `ε` is small slack (~10 bits). This
inherits GGT's regularisation philosophy at the precision the user
asked for, rather than hard-coding a Float64 constant.

#### 2.1.4 Failure modes GGT eliminates

- **Froissart doublets** — spurious pole-zero pairs at arbitrary
  locations, present even in exact arithmetic for some `(m, n)`,
  ubiquitous under rounding. The SVD-rank-then-reduce loop detects
  them as small `σ_n` and removes them by reducing `m, n` along the
  diagonal. Figures 5–9 of GGT 2013 give visually striking
  before/after comparisons for `tan(z⁴)`, `log(1.2 − z)`,
  `exp((z + 1.5)^{−2})`. **This is the load-bearing benefit of
  switching from FW's direct Toeplitz solve to GGT.**
- **Block-degeneracy mishandling.** Standard Padé treatments produce
  identical entries in degenerate blocks of the Padé table; some
  match `f` to less than the expected `m + n + 1` order. GGT's
  corner-hopping consistently lands at the upper-left corner of the
  block.
- **Near-singular `C` from rapidly-decaying Taylor coefficients.** The
  exact case FW 2011 §5.2 patches with min-norm-on-underdetermined-
  system; GGT detects it via threshold-`tol` and reduces `m, n` —
  more principled.

#### 2.1.5 Ill-posedness as a feature, not a bug

`[GGT2013: §7]` — GGT freely admits Padé approximation is ill-posed
*precisely when* the approximant has defect `δ > 0`. Their position is
that an algorithm should be **stable** (find the exact solution of a
slightly perturbed problem) rather than chase numerical accuracy on a
truly ill-posed instance. **Implication for PadeTaylor.jl's API**:
we should report when `m, n` were reduced (a degeneracy detected) and
expose the reduced `(μ, ν)` to the caller alongside `p, q`. The
caller can then decide whether to trust the approximant or refuse.
This matches the scientist-workbench "honest scope" discipline.

### 2.2 Chebfun `padeapprox.m` — reference implementation + reweighting trick

Local source: `external/chebfun/padeapprox.m`. Three lines marked by
`% reweighting` comments **go beyond GGT 2013 Algorithm 2**
`[GGT2013: line 236; padeapprox.m lines 111–117]`:

```matlab
[U,S,V] = svd(C, 0);
b = V(:, n+1);                          % null vector gives initial b
D = diag(abs(b) + sqrt(eps));            % reweighting preserves zeros better
[Q, R] = qr((C*D)');                     % so does final computation via QR
b = D*Q(:, n+1);  b = b / norm(b);       % compensate for reweighting
```

> **Transpose warning (corrected 2026-08-11).** This snippet previously
> reproduced `qr((C*D).')` — MATLAB's **non-conjugating** transpose. That is
> the pre-July-2018 Chebfun bug; the live source carries the inline comment
> "until July 2018 there was an erroneous `.'` here"
> (`external/chebfun/padeapprox.m:113`). For the **complex** Taylor jets that
> dominate the path-network walker, `.'` vs `'` is a silent-wrong-answer
> difference. The implementation is already correct — `src/RobustPade.jl:443`
> and `src/SharedPade.jl:290` both use `adjoint(...)`. **Do not "fix" the code
> to match an older copy of this document.**

The motivation: **for blocks corresponding to approximation accuracies
close to `tol`**, the SVD's null vector picks up small noise in
positions where `b_k` "ought to be" zero; the QR-of-reweighted-matrix
re-finds a null vector that is **more often exactly square** in the
block-structure sense. **Adopt this** in PadeTaylor.jl — the cost is
one extra QR (small) and the robustness benefit is documented.

Other notes on `padeapprox.m`:

- It accepts either Taylor coefficients directly OR a function handle;
  in the latter case, samples the function on `2048` roots of unity
  and FFTs to coefficients `[padeapprox.m lines 254–261]`. Convenience
  wrapper; PadeTaylor.jl probably wants both APIs too.
- `tol = 1e-14` default. `tol = 0` "to turn off robustness" — nice
  escape hatch.
- Output: `(r, a, b, mu, nu, poles, residues)` with `poles, residues`
  as optional outputs from `roots(b)` and a finite-difference residue
  estimate `[padeapprox.m lines 292–296]`.
- Underlying SVD: MATLAB's `svd(C, 0)` = LAPACK Demmel-Kahan. **For
  arb-prec we substitute one-sided Jacobi (Demmel-Veselić)** per §5.3.

### 2.3 Hermite–Padé (Baker & Graves-Morris) — deferred to v2

The PRD names Hermite–Padé as candidate v1 vs v2. **v1 = scalar Padé
per component** is the conservative choice — single `(P, Q)` per
component of `y`. **Hermite–Padé** would couple all components through
a *shared* denominator `Q(z)`: `y_i(z) ≈ P_i(z) / Q(z)` for each
component. The shared denominator captures pole locations as
properties of the *system* rather than per-component, which:

- **Is mathematically more honest** for a system with a global
  meromorphic structure (the Painlevé case): every pole of `y` is a
  pole of *every* component, so a shared `Q` should match.
- **Doubles or triples the matrix size** in the Padé conversion step
  for systems of dimension `d`: instead of `d` independent
  `(n+1)×(n+1)` SVDs we have one `(dn) × (dn+1)` SVD per step.
- **Has a richer literature for branch-point detection** [Baker &
  Graves-Morris]. v2 territory once we have sample target problems
  with branch-point structure.

**Recommendation: v1 = scalar Padé per component; v2 = Hermite-Padé
when target problems demand it.** Same as PRD's tentative choice;
literature does not dispute it.

## 3. Adjacent Taylor-IVP packages — design lessons

### 3.1 TIDES (Barrio et al.)

Coefficient generation is a **two-stage pipeline**: a Mathematica
preprocessor (`MathTIDES`) parses the ODE symbolically, decomposes the
RHS into a "three-variable code" (DAG of elementary ops), and emits
ANSI C / Fortran source. The emitted code evaluates Taylor coefficients
by executing AD-derived recurrences for sums, products, quotients,
powers, exp, log, sin/cos `[Jorba-Zou 2005: §2, Proposition 2.1; TIDES:
Abad et al. ACM TOMS 39(1) 2012, Algorithm 924]`. AD recurrences are
**O(p²) per step** `[Jorba-Zou 2005: §2, Cor 2.3]`, scaling acceptably
to orders 50–100. **TIDES uses variable order as well as variable step**;
high accuracy is achieved by raising p rather than shrinking h —
amortising the O(p²) cost.

**No Padé layer; no pole-handling.** The integrator is real-axis only
and explicitly described as unsuitable for stiff systems or solutions
near poles `[Jorba-Zou 2005: §1]`.

### 3.2 TaylorIntegration.jl (Pérez-Hernández, Benet)

**Gets right** — design patterns to inherit:

- The **`jetcoeffs!` / `@taylorize` dual-path API** is exemplary.
  Plain-function mode uses operator overloading on `Taylor1{U}`,
  working coefficient-by-coefficient up to `order-1` in a loop
  `[TI.jl: src/integrator/jetcoeffs.jl lines 23–43]`. The
  `@taylorize` macro performs AST transformation via Espresso.jl to
  emit allocation-free mutating versions, yielding significant
  speedups with zero change to user code
  `[TI.jl: src/parse_eqs.jl; docs/taylorize/]`.
- **Step control implements both primary and secondary Jorba-Zou
  (2005) strategies directly** `[TI.jl: src/integrator/stepsize.jl
  lines 26, 74]`, with mixed abs/rel tolerance.
- The **`TaylorSolution` struct optionally stores the full Taylor
  polynomial at each accepted step**, enabling dense output for free.
  Useful pattern for PadeTaylor's `Solution` type.

**Scope limitation** (the gap PadeTaylor fills): the independent variable
`t` is constrained to `T <: Real` throughout the public API
`[TI.jl: src/integrator/taylorinteg.jl lines 476, 509, 542, 570]`. No
complex-plane stepping, no Padé layer, no pole-detection. **f is
provided as a plain Julia function with `f(x, p, t)` or mutating
`f!(dx, x, p, t)`**, dispatched over `Taylor1` arguments via operator
overloading. Same calling convention PadeTaylor.jl should adopt.

### 3.3 TaylorSeries.jl (Benet, Sanders) — empirically validated for `Arb`

**Structurally generic.** `Taylor1{T <: Number}` — type declaration, all
arithmetic, all recurrence-based transcendentals (exp, log, sin/cos,
sqrt, power, inverse trig) dispatch only on `T <: Number` with no
hardcoded `Float64` paths in algorithmic code `[TaylorSeries.jl:
src/arithmetic.jl, src/functions.jl, src/power.jl — no Float64
type-assertions found]`. Sole `Float64` hardcoding is in the
convenience constructor `Taylor1(order::Int)` (no type arg) which
defaults to `Float64` `[src/constructors.jl lines 92–98]`. Calling
`Taylor1(Arb, order)` uses the fully generic path.

**Supported operations on `Taylor1`:** four arithmetic operators,
integer and rational powers (via recurrence `[src/power.jl]`), `sqrt`,
`exp`/`expm1`, `log`/`log1p`, `sin`/`cos`/`tan`,
`asin`/`acos`/`atan`, `sinh`/`cosh`/`tanh`/`asinh`/`acosh`/`atanh`,
`sinpi`/`cospi`, `abs`, composition via `evaluate`, and
integration/differentiation `[src/calculus.jl]`. **No polynomial
root-finding or Padé construction included** — these we add ourselves.

**Empirical validation at `T = Arb`.** A targeted Julia probe
(`external/probes/taylorseries-arb/`, see `SYNTHESIS.md` therein) ran
five tests at orders 60 and 80 with `Arblib 1.7.0` + `TaylorSeries
0.21.7` on Julia 1.12.3:

| probe | result |
|---|---|
| Construction `Taylor1(Arb, 60)` | PASS |
| Arithmetic at order 60 (`+`, `−`, `×`, `÷`) | PASS |
| Transcendentals at order 60 (`exp`, `sin`, `cos`, `log(1+x)`) | PASS — coefficients overlap exact rationals as Arb balls |
| High-order stress at order 80 (`exp(x)²` vs `exp(2x)`) | PASS |
| Precision propagation at 256 bits | PASS — radii < `2⁻²⁰⁰`, max observed `2.16e-78` |

**Bottom-line: PadeTaylor.jl can rely on `TaylorSeries.jl::Taylor1{Arb}`
directly. No hand-rolled coefficient layer is needed.** The generic
`T <: Number` design is empirically validated for both `Float64` and
`Arb` as first-class element types. Caveat: `Arblib.radius(::Arb)`
returns `Mag`, not `Float64`; wrap with `Float64(...)` for numeric
comparisons.

### 3.4 ATOMFT (Chang & Corliss) — historical reference

ATOMFT (1994, *Computers & Mathematics with Applications* 28, 209–233)
was the **first general-purpose Taylor-series IVP system** to automate
the full pipeline: a Fortran-77 source-to-source translator reads an
ODE in a natural mathematical dialect, walks the parse tree (Lex/Yacc),
eliminates common subexpressions, and emits a compiled Fortran program
that evaluates Taylor coefficients via AD-derived recurrences — the
template later adopted by Jorba-Zou's `taylor` and TIDES's `MathTIDES`.
ATOMFT additionally pioneered **interval-arithmetic step control**,
**DAE support**, **event detection**, and **stiff-system handling via
modified Taylor methods** — a broader scope than any successor.

The legacy: every subsequent Taylor-IVP package inherits the
coefficient-via-recurrence / operator-overloading dichotomy ATOMFT
introduced; the simplified "use last one or two coefficients with
asymptotic error estimate" step heuristic in current packages is a
conscious simplification of ATOMFT's rigorous interval-arithmetic
approach, chosen for speed.

### 3.5 Jorba & Zou 2005 — the canonical step-size formula

**The formula PadeTaylor.jl should adopt as its primary step-size
control** (PRD open question 2 partially resolved):

The estimated radius of convergence is `ρ_m = min(ρ_m^{(p-1)},
ρ_m^{(p)})` where for absolute tolerance `ε`:
```
ρ_m^{(j)} = (ε / ‖x_m^{[j]}‖)^{1/j}
```
or for relative tolerance `ε_r`:
```
ρ_m^{(j)} = (ε_r · ‖x_m‖ / ‖x_m^{[j]}‖)^{1/j}
```
using the two **highest** computed coefficients `j ∈ {p-1, p}`
`[Jorba-Zou 2005: §3.2, eqs. 3-5, 3-6, 3-7]`. The primary step size:
```
h_m = (ρ_m / e²) · exp(−0.7 / (p_m − 1))
```
where the safety factor `exp(−0.7 / (p−1))` evaluates to ≈ 0.90 at
`p = 8` and ≈ 0.95 at `p = 16` `[Jorba-Zou 2005: §3.2.1, eq. 3-8]`. A
secondary control guards against primary overshoot when intermediate
coefficients are anomalously small: largest `h` such that all
`‖x_m^{[j]}‖ · h^j ≤ z` for `j = 1, …, p` `[§3.2.2]`.

**TaylorIntegration.jl implements both verbatim** `[TI.jl:
src/integrator/stepsize.jl line 26]`. PadeTaylor.jl should adopt the
same, *plus* the FW Padé-denominator-root-distance heuristic as the
alternative — see §7.1 question 2.

No complex independent variable and no Padé / pole-handling appear
anywhere in Jorba-Zou; method explicitly restricted to real
`t ∈ [a, b]` and described as unsuitable for stiff or pole-proximate
problems.

## 4. Adjacent verified / arb-prec ODE work

### 4.1 Mezzarobba — D-finite / linear ODEs in Arb (`ore_algebra`)

**Papers**: Mezzarobba, *Rigorous Multiple-Precision Evaluation of D-Finite Functions in SageMath*,
arXiv:1607.01967 (extended abstract for ICMS 2016, withdrawn from proceedings over a copyright
dispute); Mezzarobba, *Truncation Bounds for Differentially Finite Series*, Annales Henri Lebesgue
**2** (2019) 99–148, DOI 10.5802/ahl.17. Earlier related work: Mezzarobba & Salvy, *Effective
Bounds for P-Recursive Sequences*, arXiv:0904.2452. Code: `ore_algebra.analytic`, a Python/SageMath
subpackage (several thousand lines, production-quality, used by Fredrik Johansson's `fungrim`
project).

**Verified enclosure technique.** The core method is a symbolic-numeric a priori tail bound on the
Taylor series of a D-finite function. Given a linear ODE with polynomial coefficients, Mezzarobba
derives a companion majorant ODE (or majorant recurrence) whose solution dominates the truncation
remainder term-by-term; the bound on the tail is then a closed-form expression evaluated in Arb
ball arithmetic, giving a rigorous enclosure of the truncated sum [Mezzarobba 2019: abstract, §1].
Coefficient arithmetic throughout uses Arb, so rounding errors in intermediate calculations are
also enclosed [Mezzarobba 2016: §2].

**Applicability to nonlinear ODEs.** The truncation-bound machinery is **fundamentally tied to
linearity**. The majorant argument exploits the fact that a linear recurrence for the Taylor
coefficients has a dominant solution controlled by the operator's Newton polygon, yielding a
Cauchy-Hadamard-style estimate for the convergence radius and a closed-form tail sum. For a
nonlinear `f(z, y)` the coefficient recurrence is nonlinear; no analogous majorant recurrence is
available without further structural assumptions. Mezzarobba's bounds therefore do **not** transfer
directly to PadeTaylor's nonlinear setting and would need to be replaced by a separate a priori
estimate (see §4.3 Synthesis below).

**Step-size strategy.** Steps are sized so that the Taylor series converges rapidly within the
disk: in practice the step is chosen to be a fixed fraction of the distance to the nearest
singularity, estimated from the leading coefficient of the differential operator — comparable in
spirit to FW's denominator-root distance heuristic [Mezzarobba 2016: §2; exact fraction is
conjecture from context].

**Architecture (design-pattern value).** `ore_algebra.analytic` separates concerns cleanly: a
`DifferentialOperator` class owns the ODE; a `Path` abstraction handles analytic continuation
along a piecewise-linear path in ℂ; and a `RigidStep` primitive encapsulates one Taylor step with
certified bounds. This operator / path / step separation is a pattern worth inheriting in
PadeTaylor's module layout regardless of whether the bounds machinery is ported.

### 4.2 COSY Infinity / Taylor models (Makino & Berz)

**Primary references**: Makino & Berz, *Efficient Control of the Dependency Problem Based on Taylor
Model Methods*, Reliable Computing **5** (1999) 3–12; Berz & Makino, *Verified Integration of ODEs
and Flows Using Differential Algebraic Methods on High-Order Taylor Models*, Reliable Computing
**4** (1998) 361–369. Wrapping-effect correction: Bünger, *Shrink Wrapping for Taylor Models
Revisited*, Numerical Algorithms **78** (2018) 1001–1017 (this paper corrects a flaw in the
original shrink-wrapping theorem). Julia implementation: `TaylorModels.jl` (JuliaIntervals
organisation), v0.10.x, actively maintained as of 2025–2026.

**What is a Taylor model.** A Taylor model (TM) of order *n* for a function *f* over domain *D*
is a pair `(P, Δ)` where `P` is a degree-*n* polynomial and `Δ` is a real interval such that
`f(x) − P(x) ∈ Δ` for all `x ∈ D` [TaylorModels.jl docs]. This is strictly more informative than
a plain interval enclosure: the polynomial part tracks smooth variation across the domain while
`Δ` bounds only the residual. Plain interval arithmetic collapses the polynomial part to a
constant interval, discarding shape information and producing wider enclosures ("dependency
problem").

**ODE integration via Taylor models.** For an IVP `y' = f(t, y)`, a TM ODE solver evaluates the
right-hand side as a TM over a box containing both the time interval and the current state
interval, builds *n* Taylor coefficient polynomials from the TM recurrence, then certifies the
step by applying a Picard fixed-point argument: one seeks a TM `(P, Δ)` such that the Picard
operator maps `Δ` into itself (Schauder theorem guarantees a fixed point once the set is compact
and convex), proving the true solution lies inside the Taylor model [nextjournal.com/UzielLinares].
The wrapping effect — accumulation of excess interval width across steps, caused by the mismatch
between the solution's natural geometry and the axis-aligned box representation — is mitigated
because the polynomial part absorbs smooth nonlinear mappings, leaving only the small residual
`Δ` to be reboxed. *Shrink wrapping* [Bünger 2018] additionally rotates the model's coordinate
frame to realign with the flow direction, further compressing `Δ`.

**TaylorModels.jl — Julia status.** The package explicitly combines `IntervalArithmetic.jl` and
`TaylorSeries.jl`, and uses `TaylorIntegration.jl` for the non-validated polynomial step before
certifying the remainder via Picard iteration. It therefore shares the same `TaylorSeries.jl`
dependency that PadeTaylor will use — a direct integration point. Validated ODE functionality has
been demonstrated only for **real-valued** problems; complex-domain ODE support is absent. The
package is research-grade (v0.10.1, ~64 stars), not production; its Picard loop is incompatible
with the Padé re-expansion step without non-trivial adaptation.

### 4.3 Synthesis — lightest-touch path to honest error bars for PadeTaylor v1

Three options exist when using Arblib.jl-backed coefficient arithmetic:

**(a) Arb mid+radius from coefficient arithmetic alone.** Every arithmetic operation in the Padé
routine is tracked by Arb's ball arithmetic, so outputs naturally carry a radius. This is free to
implement — just use `Arb` element type and read off the radius. Limitation: the radius reflects
only rounding and cancellation error in the coefficient computation; it says nothing about
**truncation error** (discarded high-order Taylor terms) or **global integration error** (drift
over many steps). The output looks like an error bar but does not bound the true ODE solution.

**(b) Option (a) plus a Cauchy–Hadamard a priori truncation bound at each step.** The truncation
error of an order-*N* Taylor expansion at `z_0` with convergence radius `ρ` is bounded by
`C · (|dz|/ρ)^{N+1} / (1 − |dz|/ρ)`, where `C` is estimated from the observed geometric decay
of computed coefficients. Adding this as an extra interval term to the Arb output gives a bound
covering truncation error. This is the pattern Mezzarobba uses for the linear case; for nonlinear
`f` the bound is heuristic rather than rigorous (the majorant argument fails), but it is
conservative in practice for well-separated singularities and requires only ~30 lines of
additional code once the convergence-radius estimate is already available for step control.

**(c) Full Taylor-model machinery.** The TaylorModels.jl Picard fixed-point loop gives a rigorous
enclosure of the true ODE solution within the step interval. This is the gold standard but
carries substantial overhead (Picard iteration per step), is currently real-only, and does not
compose with the Padé re-expansion step without a separate theory.

**Recommendation.** Option (a) for v1 — it is zero-cost and more honest than Float64 output —
with an explicit in-code caveat that the Arb radius does not bound truncation or global error.
Option (b) should be implemented as an optional flag in v1 or early v2: the incremental code is
small and the step-control logic (which estimates `ρ` already) makes it nearly free. Option (c)
is v2 or never: the real-only limitation and Padé incompatibility are genuine blockers for
pole-field traversal on ℂ. For architecture study, the Mezzarobba operator/path/step pattern
is the more useful reference than TaylorModels.jl's Picard loop, even though his truncation
bounds do not transfer to the nonlinear case.

## 5. Ecosystem inventory

### 5.1 Julia ecosystem — bigfloat / Arb SVD landscape

#### `LinearAlgebra.svd` (stdlib)

Dispatches to `LAPACK.gesdd!` for `Float32`/`Float64` only. For other
types, no generic implementation exists in stdlib: `svd` on a
`Matrix{BigFloat}` throws `MethodError` `[bigfloat-svd-synthesis: §5.1
LinearAlgebra]`.

#### `GenericLinearAlgebra.jl`

Provides generic `svd` / `svd!` via **one-sided Jacobi** (Demmel-Veselić
style) for any `T <: AbstractFloat` including `BigFloat`. CI-maintained,
moderate maturity; convergence difficulties documented for `κ > 10¹⁸`.
For GGT 2013 matrices at `n ≤ 60` (the typical range) this is the only
production-grade Julia path. **Works for `BigFloat`; `Arb`-element-type
should work given the required scalar interface (`abs`, `sqrt`, `/`)
but has not been formally tested by maintainers**
`[bigfloat-svd-synthesis: §5.1 GenericLinearAlgebra]`.

#### `Arblib.jl`

Confirmed by full source inspection (`src/matrix.jl`, `src/eigen.jl`,
`src/Arblib.jl`): **no SVD on `ArbMatrix` or `AcbMatrix`**. The matrix
layer ships LU, `ldiv!`, `norm`, `mul!`, and `approx_eig_qr!`
(wrapping Flint's C routine) — but no singular-value factorisation
`[bigfloat-svd-synthesis: §5.1 Arblib]`. Intended path: convert
`ArbMatrix` element-by-element to plain `Matrix{Arb}` and pass to
`GenericLinearAlgebra.svd`.

#### `GenericSchur.jl`

Generic Schur / eigenvalues for non-symmetric matrices, not SVD. Not
relevant to GGT 2013.

#### `TaylorSeries.jl` with `Arb` element type — empirically OK

See §3.3 — empirically validated (5/5 probes pass at order 80).

#### `Polynomials.jl` for roots of `Q` (step control)

Outstanding: needs hands-on confirmation of `Arb`-coefficient root-
finding accuracy. Likely requires `roots()` to dispatch to a generic
companion-matrix eigenvalue path. `GenericSchur.jl` provides this for
`BigFloat` and would be the natural backend.

### 5.2 TypeScript / Scientist-Workbench side

Confirmed by directory listing and grep:

- **`@workbench/bigfloat`** (per scientist-workbench ADR-0020) provides
  scalar arithmetic and transcendentals (`BigFloat`, `BigComplex`) on
  `BigInt`. **No matrix type, no SVD.**
- **`@workbench/linalg-core`** (per ADR-0014) provides `Float64`-only
  matrix operations. `Matrix` hardwired to `Float64Array`; no generic
  type parameter.
- **`tools/linalg-svd/`** (per ADR-0015/0016) — `Float64`-only.

**No bigfloat SVD substrate exists in the workbench.** Adding one is a
first-class sub-task of the TS implementation, comparable in scope to
the Padé routine itself
`[bigfloat-svd-synthesis: §5.2]`.

### 5.3 Algorithm choice — one-sided Jacobi for both implementations

**What GGT / Chebfun use** (confirmed in `external/chebfun/padeapprox.m`,
lines 93 and 106): MATLAB's `svd(C, 0)` — LAPACK `DGESVD`,
Golub-Kahan bidiagonalisation + Demmel-Kahan implicit-shift QR. Rank
decision: `sum(svd(C) > tol * norm(c))` — absolute threshold
dominated by large singular values, so Demmel-Kahan's weakness on
*small* singular values does not destabilise the count at `Float64`
where `κ · ε ≪ 1`
`[bigfloat-svd-synthesis: §5.3 ¶1]`.

**At arbitrary precision this changes.** Near a Padé-table block
boundary the gap between the last "signal" SV and the first "noise"
SV may span only a few orders of magnitude. Demmel-Kahan guarantees
absolute error `c · 2⁻ᵖ · σ_max` per SV; **one-sided Jacobi
(Demmel-Veselić 1992) guarantees relative error `c · 2⁻ᵖ · σᵢ` per
SV**. For rank-counting at arb-prec, the relative guarantee is
load-bearing — small-but-genuine singular values stay reliably above
the threshold even as precision increases
`[bigfloat-svd-synthesis: §5.3 ¶2]`.

**Cost**: GGT matrices are at most `(n+1) × (n+1)` with `n ≤ 60` in
typical use. The workbench's own SVD dispatch sends `n ≤ 500` to
Jacobi precisely because small-SV accuracy outweighs Golub-Reinsch
speed. At `n ≤ 60` Jacobi is trivially fast even with `BigFloat`
multiplication.

**Recommendation** `[bigfloat-svd-synthesis: §5.3 final block]`:

- **Julia v1**: use `GenericLinearAlgebra.svd` over `Matrix{BigFloat}`
  (convert from `ArbMatrix` element-by-element). Verify the `Arb`
  element-type path compiles; a small `Arb`-sign-behaviour shim may be
  needed. If `GenericLinearAlgebra` proves flaky on near-deficient GGT
  matrices, port the Demmel-Veselić loop from the workbench as a
  fallback — **do not reach for Golub-Reinsch**.
- **TypeScript mirror v1**: add `svdJacobiBigFloat` in a new
  `@workbench/linalg-bigfloat` package as a direct port of
  `linalg-core/src/svd.ts:svdJacobi`, replacing `Float64Array` with
  `BigFloat[]` and routing all arithmetic through `@workbench/bigfloat`.
  ~200 LOC. Skip Golub-Reinsch entirely for v1.
- **Do not** implement Demmel-Kahan bidiagonalisation from scratch at
  arb-prec. Zero-shift correctness is subtle, and Chebfun's Float64
  evidence does not transfer.

### 5.4 `@workbench/cas-core` Taylor / FPS class

`@workbench/cas-core` ships `Poly<T>` over `Field<T>` (per
scientist-workbench README) — a polynomial type, not a formal power
series type. **No FPS class exists.** Implications for the TS mirror:
the coefficient layer is built from `Poly<BigFloat>` plus an explicit
truncation order. The "Taylor1{T}" abstraction from
`TaylorSeries.jl` will need a TS analogue; rough scope ~150 LOC,
fits in a `@workbench/taylor-series` package or directly inside
`@workbench/pade-taylor`.

## 6. Source-code hunt

### 6.1 Paper appendices

- **FW 2011** has no code appendix (verified by grep over the
  markdown). The text of §5.2 describes the Toeplitz solve in detail
  (eqs. (5.4), (5.5)); these are sufficient to reconstruct the
  algorithm. References paths: `chebif`/`chebint` from the DMSUITE
  package `[FW2011: line 195; ref [29] = Weideman-Reddy 2000 ACM TOMS
  26]` for the BVP solver.
- **GGT 2013** has the `padeapprox.m` MATLAB code printed verbatim
  in Figure 1 `[GGT2013: lines 242–297]`. **This is the canonical
  reference impl** and we have the source via `external/chebfun/
  padeapprox.m` (commit 7574c77).

### 6.2 Author pages, arXiv ancillaries, and external repos

- **Weideman Stellenbosch page** (`appliedmaths.sun.ac.za/~weideman/`)
  — Apache 2.0 OK, page reachable; **no PFS code is posted directly**.
  Has links to MATLAB Differentiation Matrix Suite (DMSUITE), useful
  for the BVP solver but not the IVP integrator.
- **Fornberg CU Boulder** — page reachable; no posted PFS code.
- **arXiv ancillary files** for FW papers — not searched extensively
  during Stage 0 (acquisition agent prioritised PDF acquisition);
  worth a follow-up.
- **External repos cloned** under `external/`:
  - `chebfun/` — `padeapprox.m` verified present.
  - `TaylorSeries.jl/`, `TaylorIntegration.jl/`, `Polynomials.jl/`,
    `Arblib.jl/` — for ecosystem audit (§3, §5).
  - `probes/taylorseries-arb/` — empirical validation (§3.3).

### 6.3 Direct outreach — **NOT pursued**

The original PRD listed email outreach to Weideman, Fornberg,
Fasondini, Reeger as a Stage 0 source-code-hunt option. **Per user
direction (this session), no outreach will be attempted** — authors
are retired or no longer active.

**Mitigation: not blocking.** We have FW 2011 §5.2 (the Toeplitz solve
algorithmic detail), GGT 2013 §2 + Figure 1's `padeapprox.m` source,
and the four canonical PDFs in `references/`. That is enough to
reconstruct the algorithm end-to-end without external code.

## 7. Resolutions of `PRD.md` open questions

For each provisional claim or open question in `PRD.md`, this section
states the literature-backed resolution.

### 7.1 PRD §"Stage 1 — Design decisions to lock"

#### 1. Padé form: scalar per component (v1) vs Hermite–Padé (v2)

**Recommendation: scalar per component for v1.** Confirms PRD's
provisional choice. Justification: §2.3 — Hermite-Padé doubles or
triples matrix sizes per step in `d`-dimensional systems, has a
richer literature for branch-point detection, and is genuinely
needed only when the target problem has multi-component shared
meromorphic structure. The PI/PII/PIV target corpus is scalar
(second-order ODE rewritten as a 2-vector), so per-component scalar
Padé gives identical results to Hermite-Padé in this regime.
Re-examine when the target corpus expands to systems with
*genuinely* shared meromorphy.

#### 2. Step control: Jorba-Zou primary as default, FW Padé-denominator-root distance as flag

**Recommendation: Jorba-Zou (§3.5) as default; expose FW-style
Padé-denominator-root distance as `step_control = :pade_root` flag.**

Justification:
- Jorba-Zou is the **canonical step-size formula** for Taylor IVP
  (TaylorIntegration.jl uses it verbatim). Coefficient-decay-based,
  works on any analytic ODE, no extra computation beyond the
  coefficients we already need.
- The FW Padé-denominator-root-distance heuristic is **specifically
  designed for pole-field traversal**: the closest root of the
  denominator polynomial `Q(z)` from the current point gives a
  direct estimate of the distance to the nearest pole, which is the
  step-size constraint that matters most when integrating *through*
  a pole field. This is the case where Jorba-Zou under-performs
  because the analyticity radius is *not* set by the smooth
  convergence of the Taylor series — it is set by the next pole.
- **Both heuristics are cheap.** Default to Jorba-Zou to honour the
  Taylor-IVP design lineage; expose the FW heuristic as the
  pole-field-aware alternative.

#### 3. Order strategy: fixed `n = 30` (`(15, 15)` Padé) for v1

**Recommendation: fixed order = 30 for v1.** Direct adoption of FW
2011's empirical sweet spot `(order, h) = (30, 0.5)`. Adaptive order
is a mature pattern (Jorba-Zou, TaylorIntegration.jl, TIDES) but
adds significant code complexity and step-control coupling. At v1
scope, fixed-order is sufficient and provides a clean baseline for
v2 adaptivity.

Exposed as a constructor parameter: `PadeTaylorProblem(...; order=30)`.

#### 4. Path strategy in v1: `step_complex` primitive only

**Recommendation: real-axis only with a typed `step_complex(z, dz)`
primitive.** Confirms PRD's provisional choice. Justification: a
2-D path network from day one couples the integrator to the
path-tree data structure and the BVP coupling; doing it in v1
multiplies design surface area at exactly the moment we want to
prove out the **algorithmic core**. v1 demonstrates that single
complex-stepping works (cross-validates against FW 2011's Weierstrass
℘ values to ≤ 1 ULP at Float64); v2 builds the path network on top
of an already-trusted primitive.

The `step_complex(z, dz)` primitive should accept *any* `dz ∈ ℂ`
including pure real, pure imaginary, and arbitrary direction.

#### 5. API surface: clean core + SciML `CommonSolve` thin wrapper

**Recommendation: standalone core + optional `CommonSolve` adapter.**

The clean core API:
```julia
prob = PadeTaylorProblem(f, y0, zspan; order=30)
sol = solve_pade(prob; step_size=0.5, step_control=:jorba_zou)
sol(z)  # evaluate at any z along the path
```

The `CommonSolve` wrapper lets PadeTaylor.jl participate in
SciML's `solve(prob, alg; kwargs...)` ecosystem. Keep the wrapper
thin (~50 LOC) — translate `CommonSolve.solve` → `solve_pade`,
pass through. This way SciML users get a drop-in algorithm
choice; we don't pay maintenance cost on SciML's API surface.

#### 6. Verified arithmetic via Arb balls — v1 (option a), v1-flag/v2 (option b), never-v1 (option c)

**Recommendation: §4.3 synthesis — option (a) by default, option (b)
behind a flag.** Verbatim from §4.3:

- **(a) Mid+radius from coefficient arithmetic**: free with
  `T = Arb`, gives honest *rounding-and-cancellation* error bars.
  Default for v1.
- **(b) Cauchy-Hadamard truncation bound at each step**: optional
  flag in v1 or early v2; ~30 extra LOC; uses convergence-radius
  estimate already computed for step control.
- **(c) Full Taylor-model Picard enclosures**: incompatible with
  Padé re-expansion and the complex domain without substantial new
  theory. **Not v1.** Not before someone publishes the integrating
  theory.

#### 7. Module boundaries under ~200 LOC budget

**Recommendation:** Six modules in v1 (each `≤ 200 LOC`):

| Module | Purpose | Approximate LOC |
|---|---|---|
| `Coefficients` | wrap `Taylor1{T}` from `TaylorSeries.jl` for our consumption (compute coefficients of `f(z, y)` for given `(z, y)`) | 80 |
| `RobustPade` | port of GGT 2013 Algorithm 2 + reweighting (over `T <: Number`); calls SVD primitive injected from `LinAlg` | 180 |
| `LinAlg` | thin abstraction: `svd(A::Matrix{T})` dispatching to `LinearAlgebra.svd` for `Float64` and `GenericLinearAlgebra.svd` for `BigFloat`/`Arb` | 60 |
| `StepControl` | Jorba-Zou primary + secondary; FW Padé-root distance | 150 |
| `PadeStepper` | combine `Coefficients`, `RobustPade`, `StepControl` into one Taylor → Padé → step pipeline | 200 |
| `Problems` | IVP problem and solution types; `step_complex` API; trajectory storage | 150 |

Plus an optional `CommonSolveAdapter` (~50 LOC) and a `Probes/`
test/probe directory.

Total in scope for v1: **~1100 LOC**. Plus `~600` LOC of tests +
docstrings = **~1700 LOC total v1**. Tractable in 2-4 weeks of
focused work after Stage 1 design lock.

### 7.2 PRD §"Open questions"

#### 1. Hermite–Padé in v1 or v2

**v2.** See §7.1 question 1.

#### 2. Step control default

**Jorba-Zou primary as default; FW Padé-root distance as flag.** See
§7.1 question 2.

#### 3. Verified enclosures

**v1: option (a) Arb mid+radius (default); option (b) Cauchy-Hadamard
truncation bound (flag).** See §7.1 question 6 and §4.3.

#### 4. Path-network consistency metric

**Random-path-reordering and conjugation-symmetry checks** (FW 2011
§1.1.8). Both are free — no extra precision-elevated runs needed.
Adopt both. Not active in v1 (single path only); design consistency
hooks now so v2 can flip them on.

#### 5. Naming (`PadeTaylor.jl` vs alternatives)

**Keep `PadeTaylor.jl`.** Rationale:
- It's discoverable: searching for "Padé Taylor IVP Julia" lands it.
- It's accurate: this **is** a Taylor-then-Padé method.
- Alternatives like `MeromorphicIVP.jl` and `PoleFieldSolver.jl`
  describe the *use case* but not the *method*; users searching for
  algorithms would miss them.
- The PRD ecosystem has a TS sister (`@workbench/pade-taylor`) — name
  symmetry is valuable.

#### 6. Licence

**MIT.** Justification:
- PRD's "ecosystem fit" criterion: SciML and most of the Julia
  numerical scientific stack is MIT.
- BSD is fine but offers nothing MIT doesn't.
- AGPL (the scientist-workbench's choice) blocks adoption by
  commercial Julia users — same reason MathWorks doesn't ship LGPL
  code with MATLAB. PadeTaylor.jl wants to be a Julia ecosystem
  citizen.
- Compatible with depending on `Arblib.jl` (LGPL) and
  `GenericLinearAlgebra.jl` (MIT).

The TS mirror in scientist-workbench inherits AGPL from the
workbench (per existing tools' `LICENSE`). The two licences are
deliberately different to match each ecosystem's norms; the algorithm
description is the same.

## 8. New questions surfaced during research

Questions that did not exist when the PRD was written but that Stage 0
research has surfaced:

1. **GGT 2013 vs FW 2011 Toeplitz-direct: empirical comparison
   needed.** The two methods are not equivalent in pathological cases.
   §2.1.4 lists three failure modes GGT eliminates; §1.1.5 documents
   FW's response to the same failures. Open: does GGT actually do
   better on the *specific* near-singular cases FW encounters
   (rapidly-decaying coefficients post-pole, several exact zeros)?
   **Test plan**: when v1 is implemented, run both methods on FW's
   Weierstrass ℘ test problem (FW 2011 Table 5.1), compare relative
   errors at `z = 30`, `z = 10⁴`, `z = 28.261`. The three known
   reference values are in `references/markdown/.../FW2011_*.md`
   lines 295, 372 + Table 5.1.

2. **Auto-rescaling parameter `γ` (GGT 2013 step 1)**. The reference
   `[Fornberg 1981 ACM TOMS 7]` for automatic selection has not been
   acquired; need to either acquire it or replicate the heuristic.
   **Likely v1.5 / v2** — Algorithm 2 works with `γ = 1` for
   well-scaled problems; FW 2011's Weierstrass ℘ test runs at the
   default unit-disk scaling and works.

3. **`Polynomials.jl::roots` accuracy with `Arb` coefficients.** §5.1
   notes this is outstanding. Likely needs a generic companion-matrix
   eigenvalue path via `GenericSchur.jl` for `BigFloat`. For `Arb`
   the support story is unclear and worth a probe — either before
   v1's `step_control = :pade_root` mode is implemented, or as part
   of that mode's test suite.

4. **`Arb` element type in `GenericLinearAlgebra.svd`.** §5.1 — works
   for `BigFloat`, untested for `Arb`. Probe needed; falls back to
   `convert(Matrix{BigFloat}, ArbMatrix)` if the direct path doesn't
   compile.

5. **Symbolic AD path for coefficient generation.** §1.1.1 method (a)
   — analytically differentiate `f` (e.g. via `Symbolics.jl` or a
   from-scratch AST pass) and emit a recurrence at compile time;
   compare to operator-overloading method (b). v2 question once
   benchmarks are in.

6. **Step-direction selection — five-way sampling vs steepest-descent
   on `|u|`.** §1.1.4 — recommended (b) as default. Worth a probe on
   FW 2011 Table 5.1 to confirm the FW finding that they are
   essentially equivalent in accuracy.

7. **Path-network self-consistency at the API level.** §1.1.8 — the
   FW diagnostic is to run twice with different random path orders
   and diff. Question: should the v1 API expose a `consistency_check`
   keyword that runs the integration twice and returns the diff
   alongside the trajectory? Or leave it as a user-side test pattern?
   API question to settle in DESIGN.md.

## 9. Recommendation

**Recommendation: proceed to Stage 1 (`DESIGN.md`).**

### Why

The four load-bearing technical pieces are answered with literature-
backed confidence:

| Question | Status |
|---|---|
| Algorithm core (Taylor → Padé → step → optional path) | **Spec'd** (§1.1, §3.5) |
| Padé routine choice | **Spec'd** — GGT 2013 + reweighting (§2.1, §2.2) |
| Coefficient layer for arb-prec | **Empirically validated** — `TaylorSeries.jl::Taylor1{Arb}` works at order 80, all transcendentals (§3.3) |
| Bigfloat SVD path | **Spec'd** — `GenericLinearAlgebra.svd` for Julia; port one-sided Jacobi for TS (§5.1, §5.3) |

All seven of the PRD's "Stage 1 — Design decisions to lock" have
literature-backed recommendations in §7.1. All six of the PRD's
"Open questions" have answers in §7.2.

### What's still open — but does not block Stage 1

The seven new questions in §8 are **either v2-deferred or empirical
probes that belong inside the v1 implementation**, not Stage 0
research. Each carries a concrete probe plan; none is a design
decision that has to be made before code starts.

### What we don't have, and how we mitigate

- **No author outreach** (per user direction, §6.3). We have FW 2011
  §5.2 + GGT 2013 Figure 1 + the canonical PDFs — sufficient to
  reconstruct end-to-end.
- **Willers 1974 not yet in `references/`.** Historical reference;
  superseded by FW 2011. Acquire when convenient; do not gate Stage 1.
- **No second implementation we can defer to.** PadeTaylor.jl + the
  TS mirror in `@workbench/pade-taylor` are the reference impls;
  there's no SciPy/MATLAB equivalent of GGT-Padé-IVP that we can A/B
  against. **Mitigation**: cross-validation between Julia and TS is
  load-bearing for correctness — same algorithm, same target values,
  diff to N digits. The PRD's §"Definition of done" already requires
  this. Carry through.

### Cross-validation targets we have, ready for v1

- **Weierstrass ℘ test problem** at `(z, ref value)` ∈ {(30,
  1.095098255959744), (10⁴, 21.02530339471055), (28.261, 9.876953517025014×10⁶)}
  with target relative errors 7.62e-14, 2.34e-10, 7.92e-10
  respectively at FW's `(order, h) = (30, 0.5)`, all taken directly
  from FW 2011 Table 5.1.
- **PI tritronquée IC**: `u(0) ≈ -0.1875543083404949,
  u'(0) ≈ 0.3049055602612289` from FW 2011 eq. (4.1).
- **PI near-tronquée**: `u(0) = 0, u'(0) = 1.8518` (FW 2011 Fig. 4.3).
- **PI sectorial structure**: 5 asymptotic sectors of angular width
  `2π/5` with `u(z) ≈ ±√(-z/6) + o(1)` (FW 2011 eq. (1.2)). Visual
  reproduction.
- **Pole density**: each PI solution has `O(R^{5/2})` poles within
  radius `R` (FW 2011 line 271, citing Hille). Statistical check.

These are concrete enough to start building toward.

### Stage 1 next step

Write `DESIGN.md` covering:

1. Module decomposition (the §7.1 Q7 table).
2. Public API surface for both Julia and TS, with a deliberate-1:1
   mapping (the cross-validation tax — every public function in one
   has a counterpart in the other).
3. The shared test corpus — concrete IVP problems + reference values
   (above), with separate "happy path" and "near-singular" sub-corpora.
4. The internal protocol for the four-layer pipeline (§1.1.1)
   crossing the Julia ↔ TS boundary: what is the "Taylor coefficient
   array" type in each language, what's the "Padé approximant" type,
   etc.
5. Which `T <: Number` instantiations are first-class for the
   `bd close` / "v1 done" gate: at minimum `{Float64, BigFloat,
   Arb}` for Julia; `{number, BigFloat, BigComplex}` for TS.

After `DESIGN.md` is reviewed, Stage 2 = implementation. Estimated
v1 scope: ~1700 LOC Julia + ~1700 LOC TS (counting tests/docstrings),
2-4 weeks of focused work.

**Recommendation: proceed.**

## Appendix A — references inventory

| short name | full citation | local path | source |
|---|---|---|---|
| FW 2011 | Fornberg & Weideman, *A numerical methodology for the Painlevé equations*, J. Comput. Phys. 230 (2011) 5957–5973. DOI 10.1016/j.jcp.2011.04.007 | `references/FW2011_painleve_methodology_JCP230.pdf` | user-fetch (Cloudflare-blocked direct) |
| GGT 2013 | Gonnet, Güttel & Trefethen, *Robust Padé Approximation via SVD*, SIAM Review 55(1) (2013) 101–117. DOI 10.1137/110853236 | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | user-fetch |
| RF 2014 | Reeger & Fornberg, *Painlevé IV: A numerical study of the fundamental domain and beyond*, Physica D 280–281 (2014) 1–13. DOI 10.1016/j.physd.2014.04.006 | `references/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.pdf` | user-fetch |
| FW 2015 | Fornberg & Weideman, *A computational overview of the solution space of the imaginary Painlevé II equation*, Physica D 309 (2015) 108–118. DOI 10.1016/j.physd.2015.07.008 | `references/FW2015_imaginary_PII_PhysicaD309.pdf` | user-fetch (bonus; complements PRD's FW 2014 FoCM) |
| FFW 2017 | Fasondini, Fornberg & Weideman, *Methods for the computation of the multivalued Painlevé transcendents on their Riemann surfaces*, J. Comput. Phys. 344 (2017) 36–50. Preprint version. | `references/FFW2017_painleve_riemann_surfaces_preprint.pdf` | user-fetch |
| FW 2014 (FoCM) | Fornberg & Weideman, *A computational exploration of the second Painlevé equation*, Foundations of Computational Mathematics 14 (2014) 985–1016. DOI 10.1007/s10208-013-9156-x | `references/FW2014_second_PII_exploration_FoCM14.pdf` | Springer via TIB VPN (browser UA download) |
| Jorba & Zou 2005 | Jorba & Zou, *A software package for the numerical integration of ODE by means of high-order Taylor methods*, Experimental Mathematics 14:1 (2005) 99–117. DOI 10.1080/10586458.2005.10128904 | `references/JorbaZou2005_taylor_IVP_package_ExpMath14.pdf` | Author PS.GZ converted to PDF (`ps2pdf`); source: `https://www.maia.ub.edu/dsg/2001/0103jorba.ps.gz` |
| Mezzarobba 2019 | Mezzarobba, *Truncation Bounds for Differentially Finite Series*, arXiv:1904.01919 (2019) | `references/Mezzarobba2019_truncation_bounds_dfinite_arXiv1904.pdf` | arXiv direct download |

### External source repositories

| repo | purpose | local path | commit |
|---|---|---|---|
| chebfun | Reference impl of GGT robust Padé (`padeapprox.m`) | `external/chebfun/` | 7574c77 |
| TaylorSeries.jl | Ecosystem audit — Taylor series over generic types | `external/TaylorSeries.jl/` | b411d18 |
| TaylorIntegration.jl | Ecosystem audit — Taylor IVP integrator (Julia) | `external/TaylorIntegration.jl/` | 117bde7 |
| Polynomials.jl | Ecosystem audit — polynomial roots for step control | `external/Polynomials.jl/` | f1ee586 |
| Arblib.jl | Ecosystem audit — Arb ball arithmetic for Julia | `external/Arblib.jl/` | 83406c7 |

### Outstanding (still to acquire — user-fetch required)

| short name | rationale | direct URL for browser fetch |
|---|---|---|
| Willers 1974 | Originating IVP-with-poles algorithm. Historical baseline. ACM DL Cloudflare-blocked. | `https://dl.acm.org/doi/10.1145/361089.361100` (save as `references/Willers1974_ivp_continued_fractions_CACM17.pdf`) |
| TIDES papers (Barrio et al.) | Mature Taylor-IVP design lessons. | `https://www.unizar.es/acz/05Publicaciones/Revistas/Revista61/p_077.pdf` (main TIDES paper; save as `references/BarrioEtAl2005_TIDES_taylor_IVP_RealAcadZar.pdf`) |
| Baker & Graves-Morris (Hermite-Padé) | v2 / branch-point detection. | Book; low priority for Stage 0. |
| COSY Infinity / Berz Taylor models | Verified-bound thinking for v2. | Low priority for Stage 0. |

