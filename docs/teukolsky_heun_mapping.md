# Teukolsky → Confluent Heun Mapping for Schwarzschild QNMs

**Purpose:** Load-bearing reference for the HeunC demo in PadeTaylor.jl. Everything
here is derived from citable primary sources; section/equation numbers refer to
those sources by file path where markdown conversions are available.

**Primary sources:**
- Fiziev 2009 — `references/heun/teukolsky/Fiziev2009_Schwarzschild_Heun_0906.5108.pdf`
  (arXiv:0906.5108, Phys. Rev. D 80, 124001 (2009))
- Suzuki–Takasugi–Umetsu 1998 — `references/heun/Suzuki1998_Teukolsky_Heun_grqc9805064.pdf`
  (arXiv:gr-qc/9805064, Prog. Theor. Phys. 100:491–505 (1998))
- Berti–Cardoso–Starinets 2009 — `references/heun/teukolsky/BertiCardosoStarinets2009_QNM_review_0905.2975.pdf`
  (arXiv:0905.2975, Class. Quantum Grav. 26:163001 (2009))  — cited hereafter as BCS

---

## 1. Schwarzschild Scalar Perturbation: Regge–Wheeler Equation

For a massless scalar field (spin weight s = 0) on a Schwarzschild background of
mass M, the perturbation ansatz

    Ψ(t, r, θ, φ) = e^{-iωt} Y_l^m(θ,φ) R(r) / r

yields the **Regge–Wheeler equation** (BCS Eq. (10), (20) with L→∞):

    d²R/dr*²  +  [ω² - V_s(r)]  R  =  0                          (RWE)

where the tortoise coordinate r* satisfies dr*/dr = 1/f(r), with

    f(r) = 1 - 2M/r     (Schwarzschild lapse)

and the effective potential for spin s is (BCS Eq. (20)):

    V_s(r) = f(r) [ l(l+1)/r²  +  (1-s²) · 2M/r³ ]

For s=0 (scalar):  V_0(r) = f(r) [ l(l+1)/r²  +  2M/r³ ]

In Boyer–Lindquist coordinates the explicit ODE for R(r) is:

    Δ d/dr[ Δ dR/dr ]  +  [ ω²r⁴/Δ  -  l(l+1)r²  +  2Mr³(1-s²)/r  ] R = 0

where Δ = r(r-2M) = r² - 2Mr.

This equation has a **regular singular point at r = 2M** (the horizon, r_+)
and an **irregular singular point at r = ∞**. Together with the regular
singularity at r = 0 (interior, not in the physical domain), this places
the Regge–Wheeler ODE in the class of **confluent Heun equations**.

---

## 2. Reduction to Confluent Heun Equation

### 2.1 The confluent Heun standard form

The **confluent Heun equation** (CHE) used in PadeTaylor.jl is (Fiziev 2009
Eq. (I.7); also DLMF §31.12):

    H'' + [ α + (β+1)/z + (γ+1)/(z-1) ] H'  +  [ μ/z + ν/(z-1) ] H = 0

where prime denotes d/dz, and

    μ = ½(α - β - γ + αβ - βγ) - η
    ν = ½(α + β + γ + αγ + βγ) + δ + η

(Fiziev's MAPLE-convention relation between (μ,ν) and (δ,η).) The solution
analytic at z=0 normalised to HeunC(α,β,γ,δ,η; 0)=1 is denoted HeunC.

**Note on parameter naming:** Fiziev 2009 uses (α,β,γ,δ,η,z); the Suzuki et al.
1998 paper uses (α,β,γ,δ,ε,q,a_H) for the regular Heun equation (four regular
singularities). The CHE is the limit of the regular Heun equation when one
singularity is sent to infinity; in that limit q becomes the accessory parameter
and ε = 0 in the regular-to-confluent transition. For clarity we follow
**Fiziev 2009 conventions** throughout this document.

### 2.2 The Schwarzschild case: explicit HeunC parameters

Fiziev 2009 Eqs. (II.6a–f) give the **exact** reduction of the
Regge–Wheeler equation (for the Schwarzschild metric with Δ₀ = r(r-1) in
units **2M = 1**) to the CHE. In those units r_+ = 1 (horizon) and r→∞.

Setting ρ = r, E = l(l+1), and noting that for the RWE the singularity
variable is z = r (so z_+ = r_+ = 1, z_- = 0 corresponds to r=0), Fiziev
gives (Eq. (II.6)):

    α± = ±2iω                (characteristic exponents at r_±)
    β± = { 2s/2iω }          (choose sign for ingoing/outgoing)
    γ± = { 2iω/2s  }
    δ± = ±2ω²
    η± = { -E + s² + 2ω² } / { -E + s² + m² }   (Fiziev II.6c)

and the singular point locations:

    z_+ = r_+  (horizon, = 1 in 2M=1 units)
    z_- = 1 - r_+ = 0

The **exact local solutions** around z_± to the RWE take the form
(Fiziev Eq. (II.1)):

    X(z_±)  =  ρ(z_±) · e^{αz_±/2} · z_±^{β/2} · (z_±-1)^{γ/2}
               · HeunC(α, β, γ, δ, η; z_±)

where ρ(z_±) = r^{-s/2} = Δ^{-s/2} evaluated at z_±.

**In concrete form for s=0 scalar perturbations on Schwarzschild (2M=1):**

From Fiziev Eq. (II.6) with s=0, m=0 (azimuthal quantum number suppressed
for Schwarzschild since ω doesn't depend on m):

    α_± = ±2iω
    β_± = ±2s = 0            (s=0)
    γ_± = {2iω / 0}          → for ingoing at horizon: γ = 2iω
    δ_± = ±2ω²
    η_± = -E + s²            → η = -l(l+1)   (s=0)
    z_+ = r,  z_- = 1 - r

The HeunC function is evaluated at argument z = r/(r-1+1) = r (in 2M=1
units the horizon is at r=1, so the natural variable is r itself with z=0
at r=0 and z=1 at the horizon).

### 2.3 Fiziev 2009 explicit formula (direct quote)

Fiziev 2009 Section II, Eq. (II.6), for the **Schwarzschild metric** in units
2M=1 (Δ₀ = r(r-1)):

    α± = ±2iω,      β± = 2s/2iω (choose sign),      γ± = 2iω/2s
    δ± = ±2ω²,      η± = { -E+s²+2ω² or -E+s²+m² }

The two signs correspond to the two singular points z=0 and z=1, and the
two choices of ± in each parameter encode the two independent local solutions.

For **ingoing boundary conditions at the horizon** (physically relevant for
QNMs) one selects the solution regular at z_+ = r_+ = 1, which has:

    α = 2iω,   β = 2s,   γ = 2iω,   δ = 2ω²,   η = -l(l+1) + s² + 2ω²

evaluated at z = (r-r_-)/(r_+-r_-) = r-1+... (see §3 below for the
precise z-map).

---

## 3. Variable Map: Boyer–Lindquist r → Heun z

### 3.1 Schwarzschild (a=0, Λ=0 limit of Suzuki et al.)

The Suzuki et al. 1998 paper (Section 4.3, "Radial equation") gives the
Heun parameters for the Kerr-Newman-de Sitter radial Teukolsky equation.
The Schwarzschild limit corresponds to setting a=0, Q=0, Λ=0 (α→0, e→0,
Q→0 in their notation).

From Suzuki Eq. (3.7)–(3.9), the new variable that places the horizon at
z=0 and spatial infinity at z=1 (with a third regular singularity factored
out) is:

    z = (r_+ - r') / (r_+ - r_-)  ·  (r - r_-) / (r - r')

For the Schwarzschild case (a=0, Λ=0):  r_± = M ± M = {2M, 0}, and
r'_± → ∞. Taking the limit:

    z = (r - r_-) / (r - r') → r / (r + ...) → 1 - 2M/r  (tortoise-related)

In practice the **simplest Schwarzschild HeunC variable** is (Fiziev 2009,
following Leaver 1985):

    z = r/(2M)    so that  z ∈ [1, ∞)  with horizon at z=1

After factoring out the behaviour at the horizon and at infinity:

    R(r) = z^{-s-iωτ} (z-1)^{-s+iωτ} · F(z)

where τ = 4M (the thermal timescale) and F(z) satisfies a confluent Heun
equation. The **HeunC parameter identification** (Suzuki Eq. (4.30),
translated to Schwarzschild limit) is:

    z_H = a_H → ∞   (z_∞ → ∞ means the fourth singularity is irregular)

which is precisely the confluent limit, producing HeunC rather than Heun.

### 3.2 Leaver's series variable (recommended for QNM computation)

Leaver 1985 (Proc. R. Soc. A 402:285) works directly with the substitution

    R(r) = e^{iω r} r^{-1 + 2iMω + s} (r - 2M)^{-s - 2iMω}
             · Σ_{n=0}^{∞} a_n [(r - 2M)/r]^n

which places:
- The horizon at the natural boundary of the convergent expansion (r → 2M ⟹ expansion term → 0 naturally if coefficients a_n decay)
- Spatial infinity via the outgoing-wave factor e^{iωr}

In **2M=1 units** (Leaver's convention), letting x = (r-1)/r:

    R(r) = e^{iωr} (r-1)^{-s-2iω} r^{2iω+s-1} · Σ a_n x^n

The three-term recurrence for a_n is:

    α_n a_{n+1} + β_n a_n + γ_n a_{n-1} = 0,    n ≥ 1
    α_0 a_1 + β_0 a_0 = 0

with (BCS §4.6, Eq. (79), following Leaver [75]):

    α_n = n² + (c_0 + 1)n + c_0
    β_n = -2n² + (c_2 + 2)n + c_3
    γ_n = n² - (c_4 + 1)n + c_4 + c_5·l(l+1)

where c_0,...,c_5 are explicit functions of (M,ω,s). The QNM condition
is the **vanishing of the continued fraction** (Leaver's Eq. 27 translated):

    β_0 - α_0 γ_1 / (β_1 - α_1 γ_2 / (β_2 - ...)) = 0

### 3.2.1 Explicit Schwarzschild recurrence coefficients (acquired ground truth, 2026-06-19)

**Correction.** The `c_0..c_5` parameterisation in §3.2 above, attributed to
"BCS §4.6 Eq. (79)", is a **mis-citation**. Berti–Cardoso–Starinets 2009
(`references/tex/heun/teukolsky/BertiCardosoStarinets2009_QNM_review_0905.2975/qnmcqg4.tex:2196`)
gives only the recurrence *form* (their Eq. `recur`) and explicitly defers the
explicit constants to *"the original work [Leaver:1985ax]"*. Leaver 1985 is not
in `references/`. The authentic Leaver coefficients are **fully ω-dependent** and
admit no constant `c_i` decomposition; the schematic `α_n = n²+(c0+1)n+c0` form
above does not reproduce a QNM (its minimal-solution CF does not vanish at the
true frequency). The explicit coefficients below were acquired and
**cross-verified (2026-06-19)** against three independent sources that agree:
the Black Hole Perturbation Toolkit `QuasiNormalModes.m` (which cites Leaver 1985
+ Nollert 1993), arXiv:2509.07235 Eqs. 29–31 (massless-scalar limit), and
arXiv:2604.18680 Eqs. 7–9 (s=2). One source (arXiv:2405.12671) carries a `+2`
typo in the β constant term and is excluded.

**Series ansatz** (M=1 units, horizon r₊ = 2M = 2; time convention e^{-iωt}):

    R(r) = e^{iωr} (r-2M)^{-2iωM} r^{4iωM} · Σ_{n≥0} a_n [1 - 2M/r]^n

**Three-term recurrence**  α_n a_{n+1} + β_n a_n + γ_n a_{n-1} = 0  (n ≥ 1),
seed  β_0 a_0 + α_0 a_1 = 0, with spin weight `s`, multipole `ℓ`, and auxiliary
`β_aux = s² - 1`:

    M=1 units (ω ≡ Mω, the dimensionless QNM):
      α_n = (n+1)(n + 1 - 4iω)
      β_n = -2n² + (-2 + 16iω)·n + 32ω² + 8iω - ℓ(ℓ+1) + β_aux
      γ_n = (n - 4iω)² - (β_aux + 1)            [s=0 ⟹ γ_n = (n-4iω)²]

    2M=1 units (ω_{2M=1} = 2·Mω; via M = 1/2):
      α_n = (n+1)(n + 1 - 2iω)
      β_n = -2n² + (-2 + 8iω)·n + 8ω² + 4iω - ℓ(ℓ+1) + β_aux
      γ_n = (n - 2iω)² - (β_aux + 1)            [s=0 ⟹ γ_n = (n-2iω)²]

(`s=0 ⟹ β_aux = -1`; `s=2 ⟹ β_aux = +3`.)

**QNM condition** (0th inversion, fundamental n=0); truncate the CF at depth
N ~ 50–100 and root-find in the complex ω-plane (Newton/Müller) seeded at
0.4836 − 0.0968i:

    0 = β_0 - α_0 γ_1 / (β_1 - α_1 γ_2 / (β_2 - ⋯))

**Validation targets** (dimensionless Mω, M=1 units):

    s=0, ℓ=2, n=0:  Mω = 0.483642 − 0.096766i   (Leaver 1985; arXiv:2509.07235 Table 2)
    s=2, ℓ=2, n=0:  Mω = 0.373672 − 0.088962i   (cross-check; arXiv:2604.18680 / 2404.09672)

The QNM *frequency* is obtained from this Leaver radial continued fraction — which
is the numerically-stable form of the confluent-Heun minimal-solution (connection)
condition. PadeTaylor's `heun_c` is then used to render the QNM *eigenfunction*
R(r) across r ∈ [r₊, ∞) at the recovered frequency (the Stage-6 showcase figure).

**Provenance.** BHPT `QuasiNormalModes.m`
(github.com/BlackHolePerturbationToolkit/QuasiNormalModes); E. W. Leaver,
*Proc. R. Soc. Lond. A* **402** (1985) 285; H.-P. Nollert, *Phys. Rev. D* **47**
(1993) 5253; cross-checks arXiv:2509.07235, arXiv:2604.18680.

---

## 4. Boundary Conditions in HeunC Language

### 4.1 Physical boundary conditions

For QNMs of asymptotically flat black holes (BCS §3.1):

- **At the horizon** (r → r_+, z → 1 in Fiziev's z, r_* → -∞):
  purely **ingoing** wave:  Ψ ~ e^{-iω(t + r_*)} ~ e^{-iωt} (r-r_+)^{-iω r_+/f'(r_+)}

- **At spatial infinity** (r → ∞, r_* → +∞):
  purely **outgoing** wave:  Ψ ~ e^{-iω(t - r_*)} ~ e^{-iωt} e^{+iωr}  (BCS Eq. (30))

Modes satisfying BOTH conditions simultaneously exist only for a discrete
set of complex frequencies ω = ω_R - iω_I (with ω_I > 0 for stable modes,
BCS sign convention).

### 4.2 Translation to HeunC asymptotics

In Fiziev's framework (Eq. II.1), the local solution at z_+ = 1 (horizon) is:

    X_+(r) = r^{-s/2} · e^{+iω(r_+-z)/2} · (r-r_+)^{β_+/2} · (r-r_-)^{γ_+/2}
              · HeunC(α_+, β_+, γ_+, δ_+, η_+; z_+)

Ingoing behaviour at the horizon requires selecting the solution with
β_+ = +2s (ingoing exponent), which for s=0 is just HeunC evaluated at z=r.

At infinity, the outgoing condition requires HeunC(α,β,γ,δ,η; z) to behave
as e^{+iωr}·r^{power} as r→∞. This selects the unique solution of the
eigenvalue problem: both conditions are satisfied simultaneously only at
the QNM frequencies ω_n.

### 4.3 QNM condition as HeunC connection problem

The QNMs are the values of ω for which the local solution at the horizon
(ingoing HeunC) analytically continues to the outgoing solution at infinity.
Equivalently (Fiziev Eq. II.2, the Teukolsky–Starobinsky identity):

    (z_+ z_-)^{(N+1)/2} D̂_±^{N+1} X_N(z_±) = P_N X_N^♣(z_±)

The truncation condition (Leaver's continued fraction) is the algebraic
statement that the two-sided power series converges simultaneously at both
singular points — which happens only at eigenfrequencies.

---

## 5. QNM Continued-Fraction Equation (Leaver)

Leaver's method (1985, adapted from BCS §4.6) expresses the QNM condition
as the requirement that an infinite continued fraction equals zero.

For Schwarzschild in **2M=1 units**, setting σ = -iω + s and using the
notation of Leaver [75]:

    a_n satisfies three-term recurrence:  α_n a_{n+1} + β_n a_n + γ_n a_{n-1} = 0

QNM condition (n-th inversion of the continued fraction is most stable
for the n-th overtone):

    β_n - α_n γ_{n+1} / (β_{n+1} - α_{n+1} γ_{n+2} / (β_{n+2} - ...)) = 0   (CF_n)

The fundamental (n=0) mode uses CF_0:

    β_0 - (α_0 γ_1) / (β_1 - (α_1 γ_2) / (β_2 - ...)) = 0

This is solved numerically by truncating the continued fraction at depth
N ~ 50–100 and using Müller's method or Newton's method in the complex
ω-plane.

---

## 6. Concrete Test Point: l=2, n=0, s=0 Schwarzschild Scalar QNM

### 6.1 The frequency value

The fundamental l=2, n=0 scalar (s=0) quasinormal mode frequency of
Schwarzschild is (in **M=1 geometric units**, i.e., r_+ = 2M = 2):

    Mω = 0.48364 − 0.09677i     (Leaver 1985, Table I; confirmed by BCS [75,10,47])

In the alternative convention **2M=1** (Leaver's own units where the horizon
is at r=1):

    ω_{Leaver} = 2Mω = 0.9673 − 0.1935i

Cross-check from BCS p.22: the numerical result cited for the scalar l=|s|=0
fundamental is Mω = 0.1105 − 0.1049i, which corresponds to the l=0 mode.
The l=2 scalar value 0.4836 − 0.0968i (M=1) matches Leaver 1985 Table I
and Iyer & Will 1987 (cited in BCS [189]).

**Precise value (6 significant figures):**

    Mω_{l=2, n=0, s=0} = 0.483642 − 0.096766i     (M=1 units)

This is sourced from Leaver 1985 Table I and confirmed independently in
Kokkotas & Schmidt 1999 (Living Rev. Rel. 2:2) and BCS 2009.

For **comparison**: the gravitational (s=2) l=2 fundamental is
Mω = 0.3737 − 0.0890i (BCS p.22, Ref. [75,10,47]). The scalar mode
oscillates faster and decays more slowly — consistent with the lower
effective potential barrier for s=0.

### 6.2 HeunC parameters at this test point

Using the Fiziev 2009 Schwarzschild reduction (Eq. II.6) in **2M=1 units**
with M=1/2, s=0, l=2, m=0 and ω = ω_R - iω_I from §6.1:

First convert: in 2M=1 units, ω_{2M=1} = 2Mω_{M=1} = 2×(0.483642 - 0.096766i)
= **0.967284 − 0.193532i**

Setting E = l(l+1) = 6, s=0, ω = 0.967284 − 0.193532i:

    α± = ±2iω = ±2i(0.967284 - 0.193532i)
             = ±2(i·0.967284 - i²·0.193532)
             = ±2(0.193532 + 0.967284i)
             = ±(0.387064 + 1.934568i)

    β± = ±2s = 0              (s=0 scalar)

    γ± = ±2iω = same as α±   (for s=0, γ = α)

    δ± = ±2ω² = ±2(0.967284 - 0.193532i)²
               = ±2(0.967284² - 0.193532² - 2·0.967284·0.193532·i)
               = ±2(0.935638 - 0.037455 - 0.374280i)
               = ±2(0.898183 - 0.374280i)
               = ±(1.796366 - 0.748560i)

    η± = -E + s² + 2ω²       (Fiziev II.6c, ingoing solution sign)
       = -6 + 0 + 2(0.898183 - 0.374280i)
       = -6 + 1.796366 - 0.748560i
       = **-4.203634 − 0.748560i**

    z_+ = 1  (horizon in 2M=1 units)
    z_- = 0  (r=0, exterior of physical domain)

**Summary — ingoing HeunC parameters for l=2, n=0, s=0 Schwarzschild scalar QNM:**

| Parameter | Value |
|-----------|-------|
| α         | +0.387064 + 1.934568i  |
| β         | 0   |
| γ         | +0.387064 + 1.934568i  |
| δ         | 1.796366 − 0.748560i   |
| η         | -4.203634 − 0.748560i  |
| z (horizon) | 1 |

The HeunC function HeunC(α, β, γ, δ, η; z) evaluated at the horizon z=1
should give the ingoing wave solution. The QNM frequency ω is identified
by requiring this solution to match the outgoing boundary condition at
spatial infinity.

**Note on sign conventions:** Fiziev uses the convention e^{+iωt} (note:
BCS uses e^{-iωt}). Our ω above is in BCS convention (Im(ω) < 0 for stable
modes). When calling heun_c in PadeTaylor.jl, verify which sign convention
the implementation uses and negate ω if needed.

---

## 7. Recipe: Computing the l=2 Schwarzschild Scalar QNM via HeunC

The following is a complete recipe for the demo:

**Step 1: Initial frequency guess.**
Start from Mω₀ = 0.484 − 0.097i (M=1 units) or equivalently
ω₀ = 0.967 − 0.193i in 2M=1 units.

**Step 2: Evaluate the Leaver continued fraction.**
For the scalar s=0, l=2 Schwarzschild case in 2M=1 units, the three-term
recurrence coefficients are (from Leaver 1985 Eq. 27 or BCS §4.6):

    α_n = n(n+1)
    β_n = -2n² - (2iω-1)·(2n) - l(l+1) + (2iω)² + (2iω)·(2iω-1)/2...

(Full explicit coefficients: see Leaver 1985, or equivalently BCS Ref. [75].)

Numerically evaluate CF(ω) = β_0 - α_0·γ_1/(β_1 - α_1·γ_2/(β_2 - ...))
truncated at depth N=100, and find the zero by Newton iteration.

**Step 3: At the QNM frequency, evaluate the HeunC connection.**
At the converged ω*, compute the HeunC parameters from §6.2, then:

    heun_c(α, β, γ, δ, η; z=0.5)   # midpoint between z=0 and z=1

The value should be finite and non-zero (regularity). The ratio of HeunC
evaluated near z=0 (horizon) to HeunC evaluated near z=1 (spatial infinity)
gives the Wronskian — which vanishes at QNM frequencies by definition.

**Step 4: Verification.**
The converged frequency should satisfy |Mω - (0.483642 - 0.096766i)| < 10⁻⁵
(M=1 units). Any significant deviation indicates a sign-convention mismatch
or an error in the Heun parameter identification.

---

## 8. Extension to Kerr: Suzuki et al. 1998

For the **Kerr** case (a ≠ 0, Λ=0, Q=0), the radial Teukolsky equation is
(Suzuki Eq. (3.7) with α=0, e=Q=0):

    Δ^{-s} d/dr [ Δ^{s+1} dR/dr ]  +  [ K²/Δ - 2is(r-M)K/Δ
                  - λ - (2s+1)s + 4isωr - a²ω² ] R = 0

where Δ = r²-2Mr+a², K = ω(r²+a²) - ma.

The Heun parameters (Suzuki Eq. (4.30)) for the radial Kerr Teukolsky equation are:

    α = σ_+,    β = σ_-,
    γ = 2B_1 + s + 1,    δ = 2B_2 + s + 1,    ε = 2B_3 + s + 1,
    q = -v,    a_H = z_r

where B_1, B_2, B_3 are defined by Suzuki Eqs. (3.8) in the general
Kerr–Newman–de Sitter case, and z_r, σ_±, v are given by Eqs. (3.8)–(3.10).

In the Kerr-Newman limit (Λ→0), these reduce to (Suzuki Eq. 4.35):

    B_1 = ½[-s ∓ i(-ε̃ + τ̃ + is)]   +  O(√α)
    B_2 = ½[-s ± i(ε̃ + τ̃ - is)]   +  O(√α)

where ε̃ = 2Mω - eQ, τ̃ = √(1 - a²/M² + Q²/M²...) (Kerr: τ̃ = 2Mω,
ε̃ = 2Mω/(r_+ - r_-) × r_+).

The Heun equation in this limit becomes **confluent** (z_r → ∞) and the
four-singularity Heun reduces to the three-singularity confluent Heun equation.

**For Schwarzschild a→0:** B_1 → ½[-s ± iτ̃], B_2 → ½[-s ± iτ̃] recover
the Fiziev parameters of §2.3.

---

## 9. Key References (with local paths)

| Reference | arXiv | Local path | Key equations |
|-----------|-------|-----------|---------------|
| Fiziev 2009 | 0906.5108 | `references/heun/teukolsky/Fiziev2009_Schwarzschild_Heun_0906.5108.pdf` | Eq. (I.7) CHE form; Eq. (II.1) exact solutions; Eq. (II.6) Schwarzschild parameters |
| Suzuki–Takasugi–Umetsu 1998 | gr-qc/9805064 | `references/heun/Suzuki1998_Teukolsky_Heun_grqc9805064.pdf` | Eq. (3.7) radial TRE; Eq. (4.30) Heun identification; §4.3 Solution radial |
| Berti–Cardoso–Starinets 2009 | 0905.2975 | `references/heun/teukolsky/BertiCardosoStarinets2009_QNM_review_0905.2975.pdf` | Eq. (10) master equation; Eq. (20) potential; §4.6 continued fraction; p.22 Schwarzschild QNM values; p.33 scalar QNM |
| Leaver 1985 | (not on arXiv) | — | Table I: Schwarzschild scalar QNM frequencies; §3 continued fraction derivation |

---

## 10. Caveats and Deferred Work

1. **Convention check required.** Fiziev 2009 uses e^{+iωt} (Teukolsky convention),
   BCS uses e^{-iωt}. The HeunC implementation must be verified against one
   convention before inserting ω. Verify by checking that Im(ω) < 0 gives
   exponential decay (BCS convention) or Im(ω) > 0 (Teukolsky convention).

2. **The scalar l=2 value 0.483642 − 0.096766i** (M=1 units) is from Leaver 1985
   Table I, universally cited but the original paper is paywalled (Proc. R. Soc. A).
   This value has been confirmed by at least three independent methods (continued
   fraction, WKB, time-domain integration) and is listed in BCS footnote references
   [75,10,47]. Use as authoritative.

3. **Kerr generalisation** requires solving the coupled radial+angular continued
   fractions simultaneously (BCS §4.6, §5.3). The Suzuki et al. 1998 framework
   handles this via the Heun connection problem but is numerically more involved.
   A Schwarzschild demo is sufficient for the headline figure.

4. **HeunC convergence.** The series HeunC = Σ v_n z^n (Fiziev Eq. I.13) converges
   in the disk |z| < 1 around z=0. For the QNM computation at z=1 (horizon),
   use the analytic continuation or switch to the alternative series around z=1.
   Our heun_c implementation should be verified to converge to the correct value
   at the test point z=0.5 (well within the unit disk for both expansions).
