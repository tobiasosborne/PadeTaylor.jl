# Bug sweep A3 — StepControl + VectorStepControl

Auditor: read-only numerical-methods audit, 2026-06-01.

## Area

`src/StepControl.jl` (`step_jorba_zou`, `step_pade_root`) and
`src/VectorStepControl.jl` (`vector_step_jorba_zou`), audited against:

- Jorba-Zou 2005 §3.3.1 eq. 11 (first step-size control) and §3.3.2
  eq. 12 (second step-size control).
- FW 2011 §3.1 (selection of integration paths) and §5.1 (order=30,
  h=0.5 motivation).
- The canonical Julia port the docstrings cite verbatim:
  `external/TaylorIntegration.jl/src/integrator/stepsize.jl`.

Two live consumers were also traced (out of strict scope but where the
intermittency mechanism would actually fire): the `:jorba_zou` branch
of `VectorProblems.vector_solve_pade` and `VectorPathNetworkStage2.
_validity_radius`.

## References checked

- Jorba-Zou eq. 11 (the §3.3.1 first step-size control):
  `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-645`.
  - `ρ̄_j = (|x_n^[j]|)^(-1/j)` (line 622-625) — exponent `-1/j`,
    coefficient `x_n^[j]` of order `j`.
  - `ρ_n = min{ρ̄_{p_n-1}, ρ̄_{p_n}, ρ̂_{p_n-1}, ρ̂_{p_n}}` (line 629).
  - `h_n = ± ρ_n / e²` (line 633) — the `/e²` factor, sign for
    backward integration (line 635).
  - The `≈ ε_n` algebra justifying the reduction: lines 637-645.
- Jorba-Zou eq. 12 (the §3.3.2 second step-size control):
  `.../JorbaZou2005_...md:664-672` — `|x^[0]| + h_0|x^[1]| ≥ |x^[j]| h^j`,
  `j = 2,…,p`.
- TI.jl scalar `stepsize`:
  `external/TaylorIntegration.jl/src/integrator/stepsize.jl:17-35`
  (loop `k in (ord-1, ord)`; `_stepsize` = `(epsilon/aux1)^(1/k)` at
  lines 65-69; ε-resolution at lines 27-31).
- TI.jl `_second_stepsize`:
  `.../stepsize.jl:77-90` (loop `k = 1:ord-2`, `epsilon = u = one(R)`).
- TI.jl vector `stepsize` (per-component min idiom):
  `.../stepsize.jl:37-57`.
- FW 2011 §3.1 path selection (fixed `h`, min-modulus direction):
  `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:151-168`
  (`h = 0.3` / `h = 0.5`; five-direction min-modulus selection at
  lines 158-164).
- FW 2011 order=30, h=0.5 motivation:
  `.../FW2011_...md:279`.
- Tests + oracle:
  `test/stepcontrol_test.jl:1-122`, `test/_oracle_stepcontrol.jl:1-34`,
  `test/vector_step_control_test.jl:1-203`.
- Live consumers:
  `src/VectorProblems.jl:276-291`,
  `src/VectorPathNetworkStage2.jl:310-329`,
  `src/VectorPathNetworkStage2.jl:344-352` (`_rescale_by_powers`).

## Findings

### [LOW] `step_pade_root` is mis-attributed to "FW 2011 §3.1" — the paper has no such pole-distance step rule

- **Location:** `src/StepControl.jl:17-22, 84-92, 211-239` (docstrings).
- **Ground truth (cited):** FW 2011 §3.1
  (`.../FW2011_...md:151-168`) describes a *fixed* step length
  (`h = 0.3`, then `h = 0.5` "the case used throughout the present
  work", line 164) with the step **direction** chosen as the
  minimum-modulus of five candidate rays (lines 158-164). FW §5.1
  reaffirms the fixed `h = 0.5` standard (`.../FW2011_...md:279`).
  There is no "distance to the nearest forward root of the Padé
  denominator" step-size rule anywhere in §3.1 or §5.
- **Code behavior:** `step_pade_root` is the project's own safety
  heuristic (forward-projected nearest-denominator-root distance).
  Its geometry is internally correct (see verified-correct list) and
  mutation-proofed (`test/stepcontrol_test.jl:109-117`). The function
  is **not wired into any live scalar stepper** — grep shows
  `StepControl.step_pade_root` appears only in docstrings of
  `VectorWedgeStep.jl` / `VectorPathNetwork.jl`; the live vector
  pole-distance step is the *direction-agnostic* `min|t*|` in
  `VectorWedgeStep._adaptive_h` and the `_validity_radius` clamp in
  Stage 2, neither of which calls `step_pade_root`.
- **Mechanism (intermittent?):** No runtime mechanism — the function
  is dead code in the live pipeline, and where the same idea *is* used
  (the `min|t*|` clamps) it is deliberately direction-agnostic, which
  is documented. The only hazard is an agent trusting the docstring's
  "FW 2011 §3.1 idea" attribution and porting a directional cap into a
  path solver believing it is paper-faithful. This is a Law-1
  documentation-vs-ground-truth discrepancy, not a numeric defect.
- **Intermittent?** No.
- **Confidence:** 0.85 (the mis-attribution is demonstrable against
  the cited FW lines; "no live caller" verified by grep over `src/`).

## Areas verified correct

These were checked against the cited ground truth and found faithful;
listing them so a later pass need not re-derive them.

- **Jorba-Zou exponent `1/k` (not `1/(k+1)`), correctly paired with the
  order-`k` coefficient.** `src/StepControl.jl:180-184` loops
  `k in (p-1, p)`, reads `coefs[k+1]` (the order-`k` coefficient, since
  `coefs[1]` is order 0), and raises `(ε/aux)^(R(1)/k)`. This matches
  the paper's `ρ̄_j = (|x_n^[j]|)^(-1/j)`
  (`.../JorbaZou2005_...md:622-625`, exponent `-1/j` on the order-`j`
  coefficient) and TI.jl's `_stepsize` = `(epsilon/aux1)^(1/k)` over
  `k in (ord-1, ord)` (`stepsize.jl:23, 65-69`). The index/exponent
  pairing (`coefs[k+1]` ↔ `1/k`) is consistent. Mutation A in the test
  file flips `1/k → 1/(k+1)` and RED's 4.1.1/4.1.2
  (`test/stepcontrol_test.jl:103-107`).

- **Missing `/e²` is correct, not a dropped safety factor.** The paper
  writes `h_n = ρ_n/e²` (`.../JorbaZou2005_...md:633`) at the
  *adaptively chosen* order `p_n = -½ln ε + 1`. At that `p_n`,
  `1/e^{2(p_n-1)} ≈ ε_n` (the paper's own algebra,
  `.../JorbaZou2005_...md:637-645`), so the `ρ/e²` form and the direct
  "make the last term equal ε" form `(ε/|c_k|)^(1/k)` target the same
  thing. At the project's fixed `p = 30` the TI.jl direct form is the
  correct reduction. Cross-validated against TI.jl, mpmath, and
  Mathematica to 47 digits: `4.50120637033898607690…`
  (`test/stepcontrol_test.jl:48-59`, `test/_oracle_stepcontrol.jl:19-24`).
  The legend "Jorba-Zou §3.2.1 eq. 3-8 / `(ρ/e²)·exp(-0.7/(p-1))`" in
  the old DESIGN.md is correctly flagged as a hallucination in the
  module docstring (`src/StepControl.jl:24-64`) and the test header
  (`test/stepcontrol_test.jl:11-26`).

- **Float exponent, no integer-division truncation.** `R = float(real(T))`,
  so `R(1)/k` and `R(1)/j` are float divisions (`StepControl.jl:170,
  183, 197`; `VectorStepControl.jl:181, 193, 207`). No `1//k` integer
  path.

- **No branch-cut hazard on complex coefficients.** The base of every
  fractional power is real non-negative (`abs(coefs[k+1])` /
  `vnorm(...)` divided into a real ε), so `^(1/k)` never hits the
  principal-branch cut of complex `^`
  (`StepControl.jl:181-183`; `VectorStepControl.jl:186-193`).

- **`_second_stepsize` fallback bounds and `ε = 1` target match TI.jl.**
  `StepControl.jl:194-198` loops `j in 1:(p-2)` (excluding the already-
  tried `p-1, p`), solving `|c_j| h^j = 1` via `(R(1)/aux)^(R(1)/j)` —
  exactly TI.jl `_second_stepsize` `for k = 1:ord-2` with
  `epsilon = u = one(R)` (`stepsize.jl:83-87`). The fallback is
  guarded by `isinf(h)` after the primary loop
  (`StepControl.jl:186-188`), mirroring TI.jl's vector `isinf(h)` gate
  (`stepsize.jl:49`). PadeTaylor adds a Rule-1 throw when even the
  fallback finds no nonzero coefficient (`StepControl.jl:199-203`),
  where scalar TI.jl would return `Inf` — a deliberate fail-loud
  enhancement, not a divergence.

- **ε-resolution matches TI.jl exactly.** `eps_abs ≥ eps_rel·|c₀|`
  selects abs vs rel mode (`StepControl.jl:177`), identical to TI.jl
  `absepsilon ≥ relepsilon * norm(x0, Inf)` (`stepsize.jl:27`).
  Factoring it out of the loop is sound because the test is
  coefficient-independent.

- **`step_pade_root` projection geometry, strict `t > 0`, and the cap
  are correct.** `t = real((r - z_current)·conj(unit))`
  (`StepControl.jl:254`) is the standard real scalar projection of
  `r - z_current` onto the unit step direction; `t > 0` keeps only
  forward roots (`:255`); `h = min(min forward t, |Δ|)` caps at the
  nearest forward pole but never exceeds the requested distance
  (`:251, 257, 261`). Degenerate `Δ = 0` returns 0 (`:243`); constant
  denominator returns full distance (`:246`); no-forward-root falls
  back to full distance (`:261`). Verified by tests 4.1.3 (pole at
  z=2 ⇒ 2.0) and 4.1.4 (poles 3±i, real-axis projection ⇒ 3.0),
  `test/stepcontrol_test.jl:72-94`, and mutation C (real→imag)
  RED's both (`:109-117`). `Polynomial(P.b)` uses `b[1]` = constant
  term, consistent with the `PadeApproximant` convention
  (`test/stepcontrol_test.jl:73-76`). NB the *attribution* of this
  function to FW §3.1 is wrong — see the LOW finding above — but the
  math is sound on its own terms.

- **VectorStepControl is a faithful norm-lift of the scalar selector
  and reduces exactly at d=1.** `vector_step_jorba_zou`
  (`VectorStepControl.jl:158-215`) replaces `|c_k|` with
  `vnorm(_coef_vector(jets,k,d))` where `_coef_vector` assembles the
  cross-component slice `[jets[i][k+1] for i in 1:d]`
  (`:225-226`). For `d=1` this is `[jets[1][k+1]]` and
  `norm([x]) = abs(x)`, so the formula, ε-resolution, zero-skip, and
  fallback collapse onto the scalar `step_jorba_zou` — asserted
  bit-equal by test VSC.1.1 (`test/vector_step_control_test.jl:40-63`).
  Default 2-norm is the more conservative (smaller-step) choice vs the
  ∞-norm / per-component-min idiom (the documented ordering, verified
  by VSC.1.3, `:89-107`). The `_coef_vector` slice, the exponent, and
  the h-ceiling were each mutation-proofed (M1/M2/M3,
  `test/vector_step_control_test.jl:171-200`).

- **Live consumer scalings are internally consistent.**
  `VectorProblems.vector_solve_pade` (`src/VectorProblems.jl:278-281`)
  feeds *raw* z-variable jets and treats `h_jz` as a z-step capped by
  the user ceiling — correct, since `vector_taylor_coefficients`
  returns z-coefficients. `VectorPathNetworkStage2._validity_radius`
  (`src/VectorPathNetworkStage2.jl:315-317`) feeds *rescaled* jets
  `c̃_k = c_k·h_v^k` (`_rescale_by_powers`, `:344-352`, verified to
  produce `c[k+1]·h^k`) and then multiplies the returned `h_JZ` by
  `h_v` to return to z-units — also correct. Both consumers handle the
  t-vs-z variable bookkeeping consistently; no off-by-one in the power
  of `h` in the rescaling.

## Summary

No numeric/algorithmic transcription bug was found in either assigned
file. The Jorba-Zou exponent (`1/k`), the dropped `/e²` (correct
reduction, three-source oracle to 47 digits), the second-stepsize
fallback bounds, the ε-resolution, the complex-power branch safety, the
`step_pade_root` projection geometry, and the d=1 vector reduction are
all faithful to the cited ground truth and are mutation-proofed. The
only defect is a LOW-severity Law-1 documentation issue: `step_pade_root`
is attributed to "the FW 2011 §3.1 idea", but FW §3.1 uses a *fixed*
step with min-modulus direction selection and contains no
pole-distance step-size rule; the function is the project's own
(unwired) heuristic. This is not a source of intermittent
discontinuities. If the maintainer's intermittency is real, it more
plausibly lives in the live direction-agnostic `min|t*|` step laws of
`VectorWedgeStep._adaptive_h` / `VectorPathNetworkStage2` (outside this
assignment's two files) or in the core Padé/Taylor coefficient
arithmetic.
