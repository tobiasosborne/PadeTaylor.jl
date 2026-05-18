# Heun Functions: Research Synthesis for PadeTaylor.jl

**Date:** 2026-05-18  
**Purpose:** Ground-truth reference for the forthcoming Heun-function module.  
This document is in RESEARCH.md style — cite it by path + line number in ADRs and commits.

---

## 1. Parameter Conventions: DLMF vs Mathematica vs Maple

The general Heun ODE has six free parameters. All three authorities use the
same canonical form but differ in naming and what they treat as derived vs
independent.

### 1.1 The ODE (universal)

```
d²w/dz² + (γ/z + δ/(z-1) + ε/(z-a)) dw/dz + (αβz - q) / (z(z-1)(z-a)) w = 0
```

Four regular singular points: z = 0, 1, a, ∞.

**Fuchs relation** (constraint linking all exponent parameters):

```
α + β + 1 = γ + δ + ε
```

Sources: DLMF §31.2, eq. 31.2.1 (dlmf.nist.gov/31.2); Hortacsu 2018 eq. (6)
(`references/heun/Hortacsu2018_Heun_applications_physics_1101.0471.pdf`);
Motygin 2015 §2 eq. (1)
(`references/heun/Motygin2015_Heun_evaluation_1506.03848.pdf`).

### 1.2 Singular-point exponents

| Singular point | Exponent 1 | Exponent 2 |
|----------------|-----------|-----------|
| z = 0          | 0         | 1 − γ     |
| z = 1          | 0         | 1 − δ     |
| z = a          | 0         | 1 − ε     |
| z = ∞          | α         | β         |

Source: DLMF §31.2.

### 1.3 Function signature conventions

**DLMF** (dlmf.nist.gov/31.3):

```
Hl(a, q; α, β, γ, δ; z)
```

"Hl" = "local Heun function," holomorphic at z = 0, normalized Hl(0) = 1.
The parameter ε is always the derived quantity ε = α + β + 1 − γ − δ; it
does not appear explicitly in the DLMF signature.

**Maple** (Maple HeunG help page):

```
HeunG(a, q, α, β, γ, δ, z)
```

Same order as DLMF. ε implicit. Initial conditions:  
  HeunG(a,q,α,β,γ,δ,0) = 1,  
  HeunG'(a,q,α,β,γ,δ,0) = q/(γ·a).

**Wolfram Mathematica** (reference.wolfram.com/language/ref/HeunG.html,
introduced in Wolfram Language 12.1, 2020):

```
HeunG[a, q, α, β, γ, δ, z]
```

Same order and normalization as Maple. The ODE coefficient at z = a reads
explicitly as (1 + α + β − γ − δ)/(z − a) so ε is written out but equals
the Fuchs-constrained value. Branch cuts: [1, +∞) and [a, e^(i·arg(a))·∞)
(Wolfram blog, 2020-05-06).

**Critical difference — HeunC (confluent):**

DLMF §31.12, eq. 31.12.1:
```
d²w/dz² + (γ/z + δ/(z-1) + ε) dw/dz + (αz - q)/(z(z-1)) w = 0
```
Wolfram: `HeunC[q, α, γ, δ, ε, z]`  
Maple: uses a *completely different* parameterization with parameters (β, γ, δ, η)
and a restructured equation. Direct parameter mapping requires a conversion
table (see heun.xyz CAS comparison table).

**Summary table:**

| Symbol  | Role                     | Free? | DLMF | Maple | Mathematica |
|---------|--------------------------|-------|------|-------|-------------|
| a       | 3rd singular point       | Yes   | a    | a     | a           |
| q       | accessory parameter      | Yes   | q    | q     | q           |
| α, β    | exponents at ∞           | Yes   | α,β  | α,β   | α,β         |
| γ       | exponent at z=0          | Yes   | γ    | γ     | γ           |
| δ       | exponent at z=1          | Yes   | δ    | δ     | δ           |
| ε       | exponent at z=a          | **No** | implicit | implicit | explicit in ODE |

Source for Wolfram branch cuts: Wolfram blog 2020-05-06
("From Sine to Heun: 5 New Functions").

---

## 2. The Three-Term Recurrence for HeunG at z = 0

The local Heun function solution

```
Hl(a, q; α, β, γ, δ; z) = Σ_{j=0}^∞  c_j  z^j
```

with c₀ = 1 satisfies (DLMF §31.3, eqs. 31.3.2–31.3.4):

**Initialisation (eq. 31.3.2):**

```
a·γ·c₁ - q·c₀ = 0       →    c₁ = q / (a·γ)
```

**Three-term recurrence for j ≥ 1 (eq. 31.3.3):**

```
R_j · c_{j+1}  -  (Q_j + q) · c_j  +  P_j · c_{j-1}  =  0
```

**Recurrence coefficients (eq. 31.3.4):**

```
P_j = (j - 1 + α)(j - 1 + β)
Q_j = j [ (j - 1 + γ)(1 + a) + a·δ + ε ]
R_j = a(j + 1)(j + γ)
```

where ε = α + β + 1 − γ − δ throughout.

**Rearranged for forward computation** (Motygin 2015 §3, eq. 3):

```
P_j · b_j  =  Q_j · b_{j-1}  +  R_j · b_{j-2}
```

with coefficients (Motygin 2015 eq. 4, matching Motygin's sign/indexing convention):

```
P_j = a·j·(γ - 1 + j)
Q_j = q + (j-1)[(a+1)(γ+j-2) + ε + a·δ]
R_j = -(j - 2 + α)(j - 2 + β)
```

(The two formulations are algebraically equivalent after reindexing j.)

**Radius of convergence** of the Taylor series: |z| < min(1, |a|), since the
nearest singular points other than z = 0 are z = 1 and z = a.

Sources:
- DLMF §31.3, eqs. 31.3.2–31.3.4 (dlmf.nist.gov/31.3)
- Motygin 2015, §3, eqs. (3)–(4) (`references/heun/Motygin2015_Heun_evaluation_1506.03848.pdf`)

**Critical remark:** Unlike the hypergeometric 2F1, whose Taylor coefficients
satisfy a *two*-term recurrence (ratio of consecutive terms), the Heun
coefficients satisfy a *three*-term recurrence. This is the fundamental
obstruction to expressing HeunG in terms of simpler functions:
no contour-integral representation in terms of simpler functions is
known (Hortacsu 2018, §2 introduction; Arscott, cited therein).

---

## 3. The Accessory Parameter q and Apparent Singularities

### 3.1 What makes q "accessory"

The ODE has four regular singular points. The exponents at each singular point
are determined by γ, δ, ε (with ε Fuchs-constrained). These five parameters
together fix the *local* behaviour near each singularity. In the three-singularity
(hypergeometric) case that is the full story. In the four-singularity case one
additional free parameter survives: q appears in the *numerator* of the
connection coefficient αβz − q and does not shift any local exponent. It
controls the *global* structure — which linear combination of local Frobenius
solutions at z = 0 extends analytically to z = 1, and what monodromy that
extension picks up.

In the Heun polynomial case (α = −n, integer): q must take one of n+1
discrete eigenvalues q_{n,m} (m = 0,…,n) for the Frobenius solution at z = 0
to be simultaneously analytic at z = 1 (and then at z = a). These eigenvalues
are the eigenvalues of the (n+1) × (n+1) tridiagonal matrix built from
P_j, Q_j, R_j (DLMF §31.5). For general α (non-integer), q is an arbitrary
complex number, and the function is generically multi-valued upon analytic
continuation around z = 1 or z = a.

### 3.2 Apparent singularities

A point z₀ is an "apparent singularity" of an ODE if it appears as a singular
point of the ODE in one (transformed/gauge) form but is actually a *regular*
point of every solution. The Heun ODE on the standard four-singularity Riemann
sphere has singularities at 0, 1, a, ∞ that are *not* apparent. However, when
the Heun ODE is rewritten in a gauged form (after absorbing e.g. one logarithmic
solution), or when it arises from a higher-genus Lax problem, an additional
singularity at z = q/(αβ) can appear in intermediate computations. This point
has characteristic exponents (0, 2), so the monodromy is trivial; it is
"apparent." The accessory parameter q is, in the language of isomonodromic
deformation theory, the position of this apparent singularity on the Riemann
sphere. Changing q while keeping α, β, γ, δ, a fixed corresponds to moving
the apparent singularity; the Painlevé VI equation governs isomonodromic
deformations of the Heun equation with q as the dynamic variable.

Source: heun.xyz (Applications section: "Isomonodromy and Painlevé equations");
Lisovyy–Naidiuk 2022 §1 (arXiv:2208.01604).

---

## 4. Branch-Cut Structure of HeunG: Multiple Sheets

### 4.1 Branch points and cuts on the principal sheet

The Heun equation has four regular singular points. At z = 0, the local
exponents are {0, 1−γ}; when 1−γ ∉ ℤ, the two Frobenius solutions differ
by a non-integer power of z, making z = 0 a branch point of the *second*
local solution. Similarly z = 1, z = a, z = ∞ are branch points unless the
corresponding exponent differences are integers.

The **single-valued local Heun function** Hl (the one analytic at z = 0 with
Hl(0) = 1) is defined on a cut plane. The standard branch cuts for the
*principal* single-valued Heun function, as adopted by Wolfram Mathematica
(Wolfram blog, 2020-05-06), are:

```
Cut 1:  [1, +∞)     (running from z=1 along the positive real axis)
Cut 2:  [a, e^(i arg a)·∞)   (from z=a in the direction of arg(a))
```

The function is analytic in ℂ \ (Cut 1 ∪ Cut 2 ∪ {0}).

For real a > 1, the two cuts lie on the positive real axis: [1, a] and [a, ∞),
effectively giving a single cut [1, ∞) (but with different sheets on [1,a] vs
[a,∞) relative to z = a).

For real a ∈ (0,1), the cuts are [1,∞) and [0,a], and the region (a,1) is
analytic on the principal sheet.

### 4.2 Multi-valuedness and other sheets

Analytic continuation around any branch point lifts Hl to a new Riemann sheet
where the returned function is a different linear combination of the two local
Frobenius solutions. The Heun equation has a **192-element symmetry group**
(isomorphic to the Coxeter group D₄, acting by Möbius transformations on z
and the four singular points), yielding 192 local solutions that are all locally
equivalent (Hortacsu 2018 §2; Wikipedia "Heun function" §Symmetry). Of these,
8 are "essentially distinct" (the Maier group of order 8 acting on the local
Heun function alone).

The **path-multiplicative solutions** (DLMF §31.6) are defined by choosing
eigenvalues ν of the monodromy about a pair (s₁, s₂) of singularities:
analytic continuation around a loop encircling s₁ and s₂ multiplies the
solution by e^(2νπi). The associated accessory-parameter eigenvalues qₘ are
discrete (m = 0, 1, 2, …) and determined by a continued-fraction equation
(DLMF §31.4).

### 4.3 Convergence on the boundary circle

A key negative result: the Taylor series Σ c_j z^j does **not** converge on
the boundary of its disc of convergence |z| = min(1, |a|). This was proved
rigorously in Choun (2020), arXiv:2002.01971. The same result holds for the
confluent Heun function (Choun 2020, arXiv:2002.08170). Practical consequence:
analytic continuation past the nearest singularity must use overlapping discs
(as in Motygin's algorithms), never a direct Taylor evaluation at boundary
points.

### 4.4 Connection formulas (status)

Explicit closed-form connection formulas relating local Frobenius solutions at
different singular points are not available for generic parameters. Lisovyy and
Naidiuk (2022, arXiv:2208.01604) give *perturbative* connection formulas (in the
accessory parameter q) that can be computed to arbitrary order, and connect
these to quasiclassical Virasoro conformal blocks (AGT correspondence).

---

## 5. The Confluent Heun Equation and HeunC

### 5.1 The ODE

The confluent Heun equation (DLMF §31.12, eq. 31.12.1; Motygin 2018 eq. 1):

```
d²w/dz² + (γ/z + δ/(z-1) + ε) dw/dz + (αz - q)/(z(z-1)) w = 0
```

Arises from the general Heun equation by merging the singular point at z = a
with z = ∞, creating an irregular singularity of rank 1 at infinity. The
remaining two regular singular points are at z = 0 (exponents 0 and 1−γ)
and z = 1 (exponents 0 and 1−δ).

### 5.2 The three-term recurrence for HeunC at z = 0

The local solution cHl(z) = Σ b_n z^n with b₀ = 1 satisfies
(Motygin 2018 §2, eqs. 3–4):

```
P_n · b_n  =  Q_n · b_{n-1}  +  R_n · b_{n-2}
```

with

```
P_n = n(γ - 1 + n)
Q_n = -q + (n-1)(γ + δ - ε + n - 2)
R_n = (n - 2)ε + α
```

Initial conditions: b₋₁ = 0, b₀ = 1.

Source: Motygin 2018 §2, eqs. (3)–(4)
(`references/heun/Motygin2018_confluent_Heun_numerical_1804.01007.pdf`).

**Note:** Motygin's paper states ε and α are *independent* parameters in this
formulation; the Fuchs constraint no longer applies (there are only 3 singular
points: two regular and one irregular). This is a critical distinction from
the HeunG case.

### 5.3 Asymptotic expansion at ∞

For ε ≠ 0 (Motygin 2018 eq. 15):

```
cHA,∞(z) = (-z)^(-α/ε) Σ_{n=0}^∞ β_n · (n!) / (εz)^n
```

A second independent solution (eq. 16):

```
cHB,∞(z) = exp(-εz) · cHA,∞(q - εγ, α - ε(γ+δ), γ, δ, -ε; z)
```

The two solutions correspond physically to ingoing and outgoing waves in black
hole perturbation problems.

### 5.4 Useful symmetry transformation

```
cH(q, α, γ, δ, ε; z) = exp(-εz) · cH(q - εγ, α - ε(γ+δ), γ, δ, -ε; z)
```

Source: Motygin 2018 eq. (11). This is used in numerical evaluation when |εz|
is large (the direct Taylor series is poorly scaled; the exp-transformed series
converges faster on the complementary half-plane).

### 5.5 Branch-cut structure of HeunC

HeunC has branch points at z = 0 (if γ ∉ ℤ) and z = 1 (if δ ∉ ℤ). The
irregular singular point at ∞ is a Stokes phenomenon point, not a branch
point per se, but different sectors at ∞ connect to different linear
combinations of the two Frobenius solutions at z = 0. Stokes lines are
defined by Re(εz) = 0 (Motygin 2018 §3); Anti-Stokes lines by Im(εz) = 0.

### 5.6 Wolfram vs Maple signature

Wolfram `HeunC[q, α, γ, δ, ε, z]` matches DLMF §31.12.1 directly.  
Maple's HeunC uses a *different* parameterization: parameters (β, γ, δ, η)
in a restructured equation. The heun.xyz website explicitly flags this as the
primary source of confusion in cross-platform work.

---

## 6. Four Confluent Forms (DLMF §31.12)

| Name               | DLMF eq.    | Regular sing. | Irregular sing.      |
|--------------------|-------------|---------------|----------------------|
| Confluent (HeunC)  | 31.12.1     | z=0, z=1      | z=∞ (rank 1)        |
| Doubly confluent   | 31.12.2     | (none)        | z=0, z=∞ (rank 1,1) |
| Biconfluent (HeunB)| 31.12.3     | z=0           | z=∞ (rank 2)        |
| Triconfluent (HeunT)| 31.12.4    | (none)        | z=∞ (rank 3)        |

Source: DLMF §31.12 (dlmf.nist.gov/31.12).

In physical applications: Kerr–(Newman–)de Sitter radial and angular Teukolsky
equations reduce to the *singly confluent* form (HeunC); in the pure Kerr
limit the equation is again singly confluent; in the Schwarzschild limit the
radial Regge–Wheeler equation is also singly confluent (Fiziev 2009 §2;
Suzuki–Takasugi–Umetsu 1998 §3).

---

## 7. Reductions to Simpler Functions

**Gauss hypergeometric 2F1** (DLMF §31.7, eq. 31.7.1):

```
Hl(1, αβ; α, β, γ, δ; z) = 2F1(α, β; γ; z)     [when a = 1]
```

Setting a = 1 merges z = 1 and z = a; the Heun equation degenerates to
the hypergeometric equation. The accessory parameter q is fixed to q = αβ.

**Lamé equation** (DLMF §31.7): substitution z = sn²(ζ, k) with
α = −ν/2, β = (ν+1)/2, γ = δ = ε = 1/2, a = k² converts Heun to Lamé.

**Mathieu equation**: limit of doubly confluent Heun (DLMF §31.12 context).

**Spheroidal wave equation**: limit of confluent Heun (Hortacsu 2018 eqs. 9–10).

**Complete classification**: Reducing HeunG to 2F1 with at least one free
parameter requires (a, q) to lie in a finite set of exceptional pairs
(DLMF §31.7 remark after eq. 31.7.4). Generic HeunG cannot be expressed in
terms of hypergeometric functions.

---

## 8. Applications Inventory

### 8.1 Black hole quasi-normal modes (primary use case)

**Kerr–Newman–de Sitter** (Suzuki–Takasugi–Umetsu 1998,
`references/heun/Suzuki1998_Teukolsky_Heun_grqc9805064.pdf`):

The Teukolsky equation for spin-s massless perturbations of Kerr–Newman–de
Sitter separates into angular and radial equations, each mapping exactly to
the **general** Heun equation (four regular singular points: the two event
horizons, the cosmological horizon, and a Cauchy horizon). The parameter
mapping is given in §4.1, eq. (4.18) of that paper:

```
Angular equation:
  α = ρ₊, β = ρ₋, γ = 2A₁+1, δ = 2A₂+1, ε = 2A₃+1
  q = -u,  a_H = z_s   (location of Cauchy horizon in z-variable)

Radial equation:
  α = σ₊, β = σ₋, γ = 2B₁+s+1, δ = 2B₂+s+1, ε = 2B₃+s+1
  q = -v,  a_H = z_r
```

Eigenvalues (angular separation constants λ) are determined by
three-term recurrence plus continued-fraction eigenvalue condition (Suzuki
eq. 4.13).

**Schwarzschild/Kerr** (Fiziev 2009, `references/heun/Fiziev2009_Heun_BH_review_0902.1277.pdf`):  
Regge–Wheeler and Teukolsky equations reduce to *confluent* Heun (HeunC)
after the substitution R = r^(s+1)(r-1)^(iω) e^(iωr) H(r):

```
HeunC parameters: α=2iω, β=2s, γ=2iω, δ=2ω², η=s²-l(l+1)
```

The paper classifies all 16 local Frobenius solutions and identifies
polynomial ("algebraically special") QNMs with spectrum:

```
ω²_{l,σ_α} = σ_α (i/6)(l-1)l(l+1)(l+2),  l = 2,3,4,...
```

### 8.2 Hydrogen molecule ion H₂⁺ and Stark effect

Born–Oppenheimer approximation in prolate spheroidal coordinates:
two singly confluent Heun equations (one per coordinate), with the
internuclear distance as a parameter. Hortacsu 2018 §2, eqs. (26)–(27).

Stark effect (hydrogen in electric field): parabolic coordinates
separate Schrödinger into two *biconfluent* Heun equations.
Hortacsu 2018 §2, eqs. (23)–(24).

### 8.3 AdS/CFT and holography

Near-horizon limits of extremal Kerr–AdS₅ metrics give Heun equations;
in non-extremal limits, conformal symmetry is broken, preventing reduction
to hypergeometric. HeunC appears in 5D Eguchi–Hanson geometry (Hortacsu
2018 §4, eqs. 36–37).

### 8.4 Bethe ansatz and integrable systems

HeunC arises in the Gaudin model and spin-chain Bethe ansatz via the
"ODE/IM correspondence"; accessory parameter values correspond to
Bethe ansatz solutions. (Mentioned in heun.xyz; detailed in
Lisovyy–Naidiuk 2022 background §1.)

### 8.5 Isomonodromic deformation and Painlevé VI

The isomonodromic deformation of the Heun equation (moving the apparent
singularity) is governed by Painlevé VI. This connection makes HeunG
central to the PadeTaylor.jl project's Painlevé infrastructure
(see `docs/adr/` for Painlevé I–V roadmap).

---

## 9. Known Regimes Where Mathematica HeunG / HeunC Misbehave

### 9.1 Integer γ: logarithmic-term singularities (Motygin 2020)

Source: Motygin 2020, arXiv:2010.09053
(`references/heun/Motygin2020_regularization_Heun_2010.09053.pdf`).

When γ ∈ ℤ≤0, the coefficient P_{n★+1} in the recurrence vanishes (where
n★ = −γ), causing the (n★+1)-th derivative to blow up as O(1/(γ+n★)).
The true local solution acquires a logarithmic term:

```
Hl(z) = Σ c_n z^n + (C_{n★} log(z) + A) Σ s_n z^n
```

Neither Maple nor Mathematica handle this smoothly: both implement the
"single-valued" definition and return spurious values near γ ∈ ℤ≤0
without warning. Motygin proposes a regularized function Hl☼ that is
C∞ in γ.

Similarly, at γ = 1, the two local Frobenius solutions coincide (Hs = Hl),
producing a second degeneration (Motygin 2020 §4).

### 9.2 Boundary of the disc of convergence

The Taylor series does **not** converge on |z| = min(1, |a|) (Choun 2020,
arXiv:2002.01971). Maple's HeunG evaluates by direct Taylor summation and
returns incorrect or non-converging results for |z| close to (but outside)
the disc boundary. Mathematica (since v12.1) uses analytic continuation via
overlapping circles but may still show accuracy loss near z = 1 or z = a.

### 9.3 Large accessory parameter q or large |ε|

Motygin 2018 §5: "numerical problems are expected" when q is large or
when singular points are nearly coincident. Near the critical point
z* = q/α, there is essential loss of significance in the residual-based
convergence test (Motygin 2018 eq. 21); Motygin switches to a norm-based
criterion (eq. 22) there. Neither Maple nor Mathematica document this regime
explicitly.

### 9.4 Maple's HeunG: unreliable code, no source access

Motygin 2015 (abstract): "the only software package able to evaluate
[general] Heun functions numerically is Maple™... [its] implementation is
imperfect" and "code is not available." Motygin's own Octave/Matlab
implementation (supplied with arXiv:1506.03848 and arXiv:1804.01007) is
the only publicly available alternative as of 2020.

### 9.5 Mathematica HeunG: no closed forms, branch-cut placement

Wolfram blog (2020-05-06): "We do not know their explicit forms and are
forced to work with their generating equations." This means Mathematica
evaluates HeunG via local power series + analytic continuation, with cuts
at [1,+∞) and [a, e^(i arg a)·∞). Any path of evaluation that crosses a
cut before the destination point will land on a different sheet.

---

## 10. Recurrence Summary (Both Equations, Both Sign Conventions)

### 10.1 HeunG local solution (z = 0, normalization Hl(0) = 1)

```
c_0 = 1
c_1 = q / (a·γ)                                          [DLMF 31.3.2]
c_{j+1} = [(Q_j + q)·c_j - P_j·c_{j-1}] / R_j           [DLMF 31.3.3]
```

with

```
P_j = (j-1+α)(j-1+β)
Q_j = j[(j-1+γ)(1+a) + a·δ + ε],   ε = α+β+1-γ-δ
R_j = a(j+1)(j+γ)
```

Radius of convergence: min(1, |a|).

### 10.2 HeunC local solution (z = 0, normalization cHl(0) = 1)

```
b_0 = 1
b_{-1} = 0
b_n = [Q_n·b_{n-1} + R_n·b_{n-2}] / P_n               [Motygin 2018, eqs 3-4]
```

with

```
P_n = n(n + γ - 1)
Q_n = -q + (n-1)(γ + δ - ε + n - 2)
R_n = α + (n-2)ε
```

Radius of convergence: 1 (nearest singular point from origin is z = 1).

**Warning:** These recurrences are *numerically forward-stable* for small |z|
but may lose precision at high index j when |c_j| grows before cancellation.
Motygin 2015 §5 recommends switching to a fresh Taylor expansion centred at an
intermediate point once the series diverges relative to a norm estimate.

---

## 11. Practical Guidance for the PadeTaylor.jl Heun Module

1. **Adopt DLMF §31.2 / §31.12 conventions** as canonical. Document any
   mapping to Wolfram/Maple conventions in the module docstring.

2. **Implement HeunG and HeunC as separate types**, inheriting from a shared
   `AbstractHeunSolution` abstract type. The three-term recurrence coefficients
   differ; conflating them is a common error.

3. **Analytic continuation**: implement the Motygin overlapping-circle method
   (Motygin 2015 §5 for HeunG; Motygin 2018 §6 for HeunC). Step size ϰ·R₀
   with ϰ ≈ 0.5 and R₀ = distance to nearest singular point.

4. **Branch-cut placement**: adopt the Wolfram convention (cuts [1,∞) and
   [a·∞)) as default; document it explicitly and add a kwarg `sheet` for
   continuation on other sheets.

5. **Guard against integer γ**: detect γ ∈ ℤ≤0 at construction time and
   either raise an informative error or implement Motygin's regularization
   (Motygin 2020). Do not silently return an incorrect series.

6. **Do not evaluate on |z| = min(1, |a|)** — the series diverges (Choun 2020).
   Analytic continuation must step *past* the boundary before evaluating.

7. **For black hole applications**: the physical boundary conditions (ingoing
   wave at horizon, outgoing at spatial infinity) select specific linear
   combinations of the two local Frobenius solutions. These are the
   "minimal solution" in Leaver's continued-fraction terminology; the
   three-term recurrence then has a unique minimal solution determined by a
   continued-fraction eigenvalue equation (cf. Suzuki §4.1 eq. 4.13).

---

## 12. Citations and Source Files

All PDFs acquired in this session live in
`references/heun/`. Marker-converted markdown (where available) is in
`references/markdown/heun/`.

| Short name         | File                                                                     | arXiv / DOI                  |
|--------------------|--------------------------------------------------------------------------|------------------------------|
| DLMF §31           | dlmf.nist.gov/31 (online, no PDF; sections 31.2–31.7, 31.12 fetched)   | —                            |
| Hortacsu 2018      | `references/heun/Hortacsu2018_Heun_applications_physics_1101.0471.pdf`  | arXiv:1101.0471              |
| Motygin 2015       | `references/heun/Motygin2015_Heun_evaluation_1506.03848.pdf`            | arXiv:1506.03848             |
| Motygin 2018       | `references/heun/Motygin2018_confluent_Heun_numerical_1804.01007.pdf`   | arXiv:1804.01007             |
| Motygin 2020       | `references/heun/Motygin2020_regularization_Heun_2010.09053.pdf`        | arXiv:2010.09053             |
| Suzuki 1998        | `references/heun/Suzuki1998_Teukolsky_Heun_grqc9805064.pdf`             | arXiv:gr-qc/9805064          |
| Fiziev 2009        | `references/heun/Fiziev2009_Heun_BH_review_0902.1277.pdf`               | arXiv:0902.1277              |
| Lisovyy–Naidiuk 22 | (not downloaded; cite arXiv:2208.01604)                                  | arXiv:2208.01604             |
| Choun 2020a        | (not downloaded; cite arXiv:2002.01971)                                  | arXiv:2002.01971             |
| Choun 2020b        | (cite arXiv:2002.08170 for confluent case)                               | arXiv:2002.08170             |

Marker-converted markdown confirmed for:
- `references/markdown/heun/teukolsky/Suzuki1998_Teukolsky_Heun_grqc9805064/`
  (`Suzuki1998_Teukolsky_Heun_grqc9805064.md`)

Marker output for remaining PDFs may be in subdirectories created during this
session; check `references/markdown/heun/` after marker jobs complete.

**Not acquired (paywall):**
- Ronveaux 1995 (Oxford University Press) — TIB VPN attempt: Hindawi returned
  HTTP 402 for hindawi.com/journals/ahep/2018/8621573/; Oxford monograph not
  accessible without institutional proxy beyond TIB VPN.
- Slavyanov–Lay 2000 — same situation; OUP paywall.

**Alternative coverage**: The mathematical content from Ronveaux (parameter
conventions, transformation theory) is fully covered by DLMF §31 + Hortacsu
2018 §2 for the purposes of this module. The Slavyanov–Lay content on
confluent forms is covered by DLMF §31.12 + Motygin 2018.
