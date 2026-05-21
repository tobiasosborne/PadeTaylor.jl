# P_I^(2) Local Pole Structure: Derivation and VC-4 Certificate

**Phase A3 of bead `padetaylor-0ln.37`.**  
**Author**: research subagent (Sonnet 4.6), 2026-05-21.  
**Context**: ADR-0025 `docs/adr/0025-headline-figure-re-resolution.md`.

---

## 1. The ODE

The second member of the Painlevé-I hierarchy at t=0 (KKG 2015 eq. (1.1),
`references/tex/.../tritronquee_coeff.tex:128`):

    u'''' + 10(u')^2 + 20 u u'' + 40(u^3 + 6x) = 0

All solutions are meromorphic (Painlevé property, KKG `tritronquee_coeff.tex:3124-3125`,
citing S. Shimomura, "Painlevé property of a degenerate Garnier system of
(9/2)-type and of a certain fourth order non-linear ODE", Ann. Scuola Norm.
Sup. Pisa 29 (2000) 1–17; `tritronquee_coeff.tex:3620-3624`).

---

## 2. Dominant Balance — Step-by-Step

Near a movable singularity x_0, try u ~ A * ξ^p, ξ = x − x_0.

Leading behaviors at small ξ:

    u''''    ~  A·p(p−1)(p−2)(p−3) · ξ^{p−4}    [power p−4]
    10(u')^2 ~  10·A²·p²           · ξ^{2p−2}   [power 2p−2]
    20·u·u'' ~  20·A²·p(p−1)       · ξ^{2p−2}   [power 2p−2]
    40·u³    ~  40·A³               · ξ^{3p}      [power 3p]

For all four to balance: p−4 = 2p−2 = 3p.
From p−4 = 2p−2: **p = −2** (double pole). Check: all powers = −6. ✓

At order ξ^{−6} with p = −2:

    u''''    coefficient:  A·(−2)(−3)(−4)(−5) = 120A
    10(u')^2 coefficient:  10·A²·4             = 40A²
    20·u·u'' coefficient:  20·A²·6             = 120A²
    40·u³    coefficient:  40A³

Balance equation:

    120A + 40A² + 120A² + 40A³ = 0
    40A(3 + 4A + A²) = 0
    40A(A + 1)(A + 3) = 0

Therefore **A ∈ {0, −1, −3}**. The trivial root A = 0 gives no pole.
The two non-trivial families are A = −1 and A = −3.

Compare with P_I (FW 2011, `references/markdown/FW2011.../FW2011...md:59-61`):
P_I has a single family A = +1, also a double pole, also zero residue.
P_I^(2) differs in having negative leading coefficients and two families.

---

## 3. Painlevé Test — Laurent Recursion

### 3.1 Ansatz

    u(x) = Σ_{k≥0} a_k · ξ^{k−2},   ξ = x − x_0

with a_0 = A. Note x = x_0 + ξ so 40·6·x = 240·x_0 + 240·ξ, contributing
source terms at specific orders.

### 3.2 Resonance Polynomial

Collecting ξ^{n−6} and isolating the term linear in a_n:

    R_n · a_n = −(nonlinear terms in a_0,...,a_{n−1}) − S_n

where S_n = 240·x_0 if n=6, 240 if n=7, 0 otherwise, and

    R_n = (n−2)(n−3)(n−4)(n−5) + A·[−40(n−2) + 20(n²−5n+12)] + 120A²

The linear-in-a_n contributions assembled:

| Source term | Linear part at order n |
|-------------|----------------------|
| u'''' | (n−2)(n−3)(n−4)(n−5)·a_n |
| 10(u')^2, j=0 and j=n | −40·A·(n−2)·a_n |
| 20·u·u'', j=0 and j=n | 20·A·[(n−2)(n−3)+6]·a_n |
| 40·u³, one index = n | 120·A²·a_n |

### 3.3 Factored Forms

Substituting A = −1:

    R_n = (n+1)(n−2)(n−5)(n−8)

Resonances for A = −1: **n ∈ {−1, 2, 5, 8}**.

Substituting A = −3:

    R_n = (n+3)(n+1)(n−8)(n−10)

Resonances for A = −3: **n ∈ {−3, −1, 8, 10}**.

### 3.4 Interpretation of Resonances

For A = −1: four resonances {−1, 2, 5, 8}, all ≥ −1. The n = −1 resonance
is universal (free pole location x_0). The remaining three positive resonances
give free parameters a_2, a_5, a_8 — a full **4-parameter** Painlevé expansion,
matching a 4th-order ODE. This is the **generic family**.

For A = −3: the resonance n = −3 is a **negative resonance** (corresponds to
a ξ^{−5} term that precedes the leading ξ^{−2}). Since the expansion starts at
ξ^{−2}, there is no ξ^{−5} term and this resonance is vacuously satisfied.
The physically relevant free parameters from n ≥ −1 are only **three**: {x_0, a_8, a_10}.
This makes A = −3 a **codimension-1 subfamily** (3-parameter).

### 3.5 Compatibility Conditions (verified symbolically)

At each resonance n = r, R_r = 0 and the equation 0·a_r = RHS must be
satisfied (otherwise log ξ terms appear). Verified by explicit symbolic
computation:

| Family | Resonance | Compatibility RHS | Status |
|--------|-----------|------------------|--------|
| A = −1 | n = 2 | 0 | ✓ satisfied |
| A = −1 | n = 5 | 0 | ✓ satisfied |
| A = −1 | n = 8 | 0 | ✓ satisfied |
| A = −3 | n = 8 | 0 | ✓ satisfied |
| A = −3 | n = 10 | 0 | ✓ satisfied |

Both families are genuine Laurent series (no log terms). Consistent with
Shimomura's theorem.

---

## 4. Explicit Coefficients

### 4.1 A = −1 family (generic)

Free parameters: β ≡ a_2, γ ≡ a_5, δ ≡ a_8.

    u(x) = −ξ^{−2}                              a_0 = −1
          + 0·ξ^{−1}                             a_1 = 0  [ZERO RESIDUE]
          + β                                     a_2 = free
          + 0·ξ                                  a_3 = 0
          + 3β²·ξ²                               a_4 = 3β²
          + γ·ξ³                                 a_5 = free
          + (−10β³ + 30x_0/7)·ξ⁴               a_6
          + (−3βγ/2 + 3)·ξ⁵                     a_7 (source contributes "+3")
          + δ·ξ⁶                                 a_8 = free
          + (−12β/7)·ξ⁷                          a_9
          + O(ξ⁸)

### 4.2 A = −3 family (special, codimension-1)

Free parameters: δ ≡ a_8, ε ≡ a_10.

    u(x) = −3·ξ^{−2}                            a_0 = −3
          + 0·ξ^{−1}                             a_1 = 0  [ZERO RESIDUE]
          + 0                                    a_2 = 0
          + 0·ξ                                  a_3 = 0
          + 0·ξ²                                 a_4 = 0
          + 0·ξ³                                 a_5 = 0
          + (−10x_0/21)·ξ⁴                       a_6 = −10x_0/21  [depends on x_0!]
          + (−1)·ξ⁵                              a_7 = −1  [universal constant]
          + δ·ξ⁶                                 a_8 = free
          + 0·ξ⁷                                 a_9 = 0
          + ε·ξ⁸                                 a_10 = free
          + O(ξ⁹)

Key distinction: for A = −3, the coefficient a_6 = −10x_0/21 explicitly
depends on the pole location x_0. This encodes the codimension-1 constraint.
The coefficient a_7 = −1 is a universal constant (from the 40·6·x source
term: −240/R_7 = −240/240 = −1).

---

## 5. Which Family for V_0?

### 5.1 Genericity

- A = −1 is the **generic** Painlevé family (4 free parameters, full-dimensional
  in the solution space). Generic initial data will produce A = −1 poles.
- A = −3 is a **codimension-1** family. For an A = −3 pole to occur, the
  solution trajectory approaching x_0 must satisfy a codimension-1 constraint.
  This is a measure-zero event in the space of initial conditions.

The tritronquée V_0(x, 0) is a specific 0-parameter solution (unique in its
sector). Its poles are governed by the generic dynamics of a 4th-order ODE near
a singularity — generically, these will be A = −1.

**Expected verdict: V_0 has A = −1 poles.**

### 5.2 Literature Status

No reference in the repository explicitly states which family V_0 poles belong
to. KKG 2015 does not carry out the Painlevé test; Shimomura (2000) is
unavailable. The genericity argument above is the only available justification.

**Conservative VC-4 design: accept A ∈ {−1, −3}.**

If empirically all V_0 poles are found to be A = −1, the criterion can be
tightened in a later revision citing the empirical evidence.

---

## 6. Literature Cross-Check

| Source | Content | File:line |
|--------|---------|-----------|
| KKG 2015 TeX | Painlevé property: "all solutions meromorphic" (cites Shimomura) | `tritronquee_coeff.tex:3124-3125` |
| KKG 2015 TeX | Shimomura reference: Ann. Scuola Norm. Sup. Pisa 29 (2000) 1–17 | `tritronquee_coeff.tex:3620-3624` |
| KKG 2015 TeX | "With Shimomura's results...U_0, V_0...extend to globally-defined meromorphic functions" | `tritronquee_coeff.tex:255-260` |
| FW 2011 | P_I Laurent: u ~ +(z−z_0)^{−2}, zero residue, second free param at O(z−z_0)^4 | `FW2011...md:59-61` |
| Cosgrove 2000 | F-V = P_I^(2); generic solution has Painlevé property | PDF p. 3 |
| `v0p2_pillarC_painleve_hierarchy_findings.md` | Full ODE, vector form, sector structure, IC | §1-4 |

The KKG TeX file does not contain the words "Laurent", "resonance", "residue",
"Fuchs", or "movable singularity" — it treats the Painlevé property as an
established theorem (Shimomura) and does not derive the pole structure.
The full Laurent recursion, resonance analysis, and compatibility check in
this report are original calculations derived from first principles.

---

## 7. VC-4 Acceptance Form (Recommended)

### 7.1 Inputs

- Extracted pole location x_0 from the path-network pole field.
- The Padé-Taylor solver (same as main computation — no new integration).

### 7.2 Procedure

1. **Ring radius**: r = min(0.05, 0.1 · d_nearest) where d_nearest is the
   distance to the nearest other extracted pole. Floor: r ≥ 0.01.
   Rationale: for A = −1, the first non-constant sub-leading term appears at
   ξ^2, so the model error is O(r^2) ≈ 0.25% at r = 0.05. For A = −3, the
   first sub-leading term appears at ξ^4, so r = 0.05 gives negligible error.

2. **Sample points**: N = 32 points on the ring,
   x_j = x_0 + r · exp(2πi·j/32), j = 0,...,31.

3. **Least-squares fit**: fit the model
       u(x_j) ≈ A·(x_j − x_0)^{−2} + B·(x_j − x_0)^{−1} + C
   by complex least squares (32 equations, 3 complex unknowns = 6 real unknowns).

### 7.3 Pass/Fail Criteria

    VC-4a (leading coefficient):
      PASS if min(|A − (−1)|, |A − (−3)|) < 0.10
      [TOL_A = 0.10 is 5% of the gap 2 between the two values]

    VC-4b (zero residue):
      PASS if |B| < 0.10 · |A|
      [B = 0 exactly from Laurent theory; tolerance accounts for fit noise]

### 7.4 Mutation Proof

- **Mutation 1**: Change ODE coefficient 40 → 50. Balance gives different
  roots, fitted A ∉ {−1, −3}. → VC-4a RED.
- **Mutation 2**: Remove 10(u')^2 term (set coefficient to 0). Different A. 
  → VC-4a RED.
- **Mutation 3**: Pole order error (wrong p). Fitted A ≈ 0 or non-integer. 
  → VC-4a RED.
- **Mutation 4**: Introduce a residue by perturbing the ODE (add A·u to RHS).
  Fitted B ≠ 0. → VC-4b RED.

### 7.5 Optional Sub-leading Checks

**Zero-residue verification (stronger)**: For the fit with a 5-term model
  u ≈ A·ξ^{−2} + B·ξ^{−1} + C + D·ξ + E·ξ²,
verify B ≈ 0 and (for A = −1) D ≈ 0 (since a_3 = 0).

**A = −3 universal constant check**: if any pole has A ≈ −3, extend the
fit to 6 terms and verify the coefficient of ξ^5 ≈ −1 (the universal
a_7 = −1). This is a strong sub-leading certificate for the A = −3 branch.

---

## 8. Summary Table

| Property | A = −1 | A = −3 |
|----------|--------|--------|
| Leading coefficient | −1 | −3 |
| Pole order p | −2 | −2 |
| Residue a_1 | 0 | 0 |
| Resonances | {−1, 2, 5, 8} | {−3, −1, 8, 10} |
| Free parameters | 4 (full generic) | 3 (codimension-1) |
| First sub-leading | a_2 free (ξ^0) | a_6 = −10x_0/21 (ξ^4) |
| Universal constants | a_7 = 3 (from source) | a_7 = −1 (from source) |
| Compatibility | All satisfied | All satisfied |
| Expected for V_0 | YES (generic) | Possible but non-generic |

**VC-4 one-liner**:
> For each extracted pole x_0, fit u on a ring of radius r ≈ 0.05 to the
> model u ≈ A·(x−x_0)^{−2} + B·(x−x_0)^{−1} + C by least squares over 32
> points. PASS iff min(|A+1|, |A+3|) < 0.10 AND |B| < 0.10·|A|.

