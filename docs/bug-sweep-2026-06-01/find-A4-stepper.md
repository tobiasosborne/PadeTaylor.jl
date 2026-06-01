# Bug sweep A4 — PadeStepper + VectorStepper

## Area

`src/PadeStepper.jl` and `src/VectorStepper.jl`: one Padé–Taylor step
(Taylor jet → `h^k` rescale → Padé → evaluate at `t = 1` → mutate state),
the coefficient rescaling, the value/derivative handoff, the FFW 2017
adaptive truncation-error controller, and the shared-`Q` vector step.

Special focus per assignment: (1) the `c̃_k = h^k c_k` rescaling power and
its consistency with evaluation; (2) the value/`u'` handoff and `state.z`
advance; (3) continuity at step joins; (4) the `d = 1` reduction of the
vector stepper to the scalar path.

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:279`
  — order = 30, `(15,15)` Padé, step `h = 0.5` (the `(order÷2, order÷2)`
  default).
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:330-350`
  — FW §5.2 eqs (5.4)/(5.5) Toeplitz Padé; line 348 explicitly: "separate
  evaluations are needed for u(z+h) and u'(z+h), whereas the Toeplitz
  approach can compute the derivative by the quotient rule of calculus" —
  ground truth for the analytic-derivative handoff in `_evaluate_pade_deriv`.
- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-93`
  — adaptive step: md:76 `ε_k = c_k + Σ_{r=1..ν} b_r c_{k-r}`; md:80
  `T(h) = |ε_{n+1} h^{n+1} / p(h)|`; md:89 eq. (2)
  `q = (k·Tol/T(h))^{1/(n+1)}`; md:93 the seed-from-accepted-`|qh|` memory.
- `docs/adr/0021-vector-step-control.md:98-109` — the `d = 1` reduction
  oracle (norm-based selector collapses onto the scalar `step_jorba_zou`).
- `src/Coefficients.jl:171-203` — `taylor_coefficients_2nd` returns the
  **value** coefficients `[c_0,…,c_order]` of `u(z0+h)` (NOT `u'`).
- `src/RobustPade.jl:110-140, 367-453` — `PadeApproximant{T}` conventions:
  `a` numerator low-to-high (`a[1]=a_0`), `b` denominator low-to-high with
  `b[1]=1`; the `r(z) = Σa[k+1]z^k / Σb[k+1]z^k` convention.
- `src/SharedPade.jl:114-123, 163-268` — shared-`Q` block-Toeplitz; the
  `d=1` reduction note; QR-reweighting `adjoint(A_full*D)`.
- `external/chebfun/padeapprox.m:113` — the documented historical
  adjoint-vs-transpose bug ("until July 2018 there was an erroneous `.'`
  here"); confirmed the Julia ports use `adjoint`, not `transpose`.
- `src/PathNetwork.jl:531-649, 856-955` — consumer of the stepper handoff;
  wedge-step value read-off and the re-seeding of `h`.
- `test/_oracle_padestepper.jl` — pinned three-source oracle values.

## Findings

### [MEDIUM] Scalar `_evaluate_pade` / `_evaluate_pade_deriv` use bare `iszero` while the vector stepper uses a relative `√ε·‖Q‖` pole guard — near-pole steps divide by a tiny-but-nonzero denominator

**Location**: `src/PadeStepper.jl:376` and `src/PadeStepper.jl:404`
(`iszero(den)` / `iszero(D)`), contrasted with
`src/VectorStepper.jl:257-258` (`abs(q_at_one) > eps_T * norm(denominator)`).

**Ground truth (cited)**: The maintainer's own VectorStepper rationale,
`src/VectorStepper.jl:74-78` and `:248-254`: "a step that lands on a pole
gives `Q(1)` of the order of floating-point roundoff on `‖Q‖`, not an exact
zero; dividing by such a `Q(1)` would yield a meaningless `~10¹⁵` value
rather than a true large solution value." This is precisely why
VectorStepper rejects on `|Q(1)| ≤ √ε·‖Q‖` instead of bare `iszero`.

**Code behavior**: The scalar second-order stepper never adopted the
relative guard. When a scalar step lands very close to (but not exactly on)
a local Padé pole, `den = D(1)` is a tiny nonzero float of order
`roundoff · ‖b‖`, so `iszero(den)` is `false`, the throw is skipped, and
`new_u = num/den` (and `new_up` via `(Nt·D − N·Dt)/D²`, which divides by
`D²` — even more sensitive) returns a spurious very-large value instead of
failing loud per Rule 1.

**Mechanism (intermittent discontinuity)**: Whether `D(1)` is "exactly
representable as zero" versus "a tiny nonzero residual" is data-dependent —
it depends on the exact float bit-pattern of the rescaled coefficients,
which in turn depends on `z`, `u`, `u'`, and `h` at that step. Two steps
that land at almost the same distance from a pole can fall on opposite
sides of the `iszero` predicate: one throws (caught upstream and retried/
shrunk), the other silently returns a `~10¹⁵`-magnitude value that the
trajectory then carries forward — a visible discontinuity at that join.
In PathNetwork the consequence is softer (a wedge candidate with huge `|u|`
loses the min-`|u|` selection, `PathNetwork.jl:559-561`), but a direct
`pade_step!` caller, or a `u'`-dependent downstream, gets the spurious
value unfiltered. The `D²` in the derivative makes `new_up` the more
exposed of the two.

**Intermittent?**: Yes — data/step-dependent; only fires when a step lands
within roundoff of a pole, which the step controller usually but not always
prevents.

**Confidence**: 0.55 — the asymmetry between the two steppers is
demonstrated and the maintainer's own comments name the mechanism; the
"medium" rating reflects that step-size control normally keeps scalar steps
off poles, so this is a latent hazard rather than an always-on bug.

### [LOW] `_rescale_by_powers` repeated-multiply accumulator can desync the rescaled coefficient power from the matching evaluation only if `T` arithmetic is non-associative — verified safe, noted for completeness

**Location**: `src/PadeStepper.jl:342-350`, `src/VectorStepper.jl:286-294`.

**Ground truth (cited)**: The intended map is `c̃[k+1] = h^k · c[k+1]`
(`src/PadeStepper.jl:336`, FFW/FW substitution `h' → h·t`). Evaluation at
`t = 1` (`src/PadeStepper.jl:319`, `:480`) must use the SAME convention.

**Code behavior**: `h_pow` starts at `one(T)` and is multiplied by `h`
AFTER each store, so `out[1] = c[1]·1 = c_0` (power 0), `out[2] = c[2]·h`
(power 1), …, `out[k] = c[k]·h^{k-1}` i.e. `c̃_{k-1} = h^{k-1} c_{k-1}`.
This is exactly the intended `c̃_k = h^k c_k` in 1-based indexing.
Evaluation reads `t = 1` so the powers cancel correctly:
`P_u(1) = Σ c̃_k = Σ h^k c_k = u(z+h)`.

**Mechanism**: None found. I traced the identity end-to-end:
`u(z+h') = P_u(h'/h)` ⇒ value at `t=1` is `u(z+h)`, derivative
`u'(z+h) = P_u'(1)/h` (`src/PadeStepper.jl:320`), matching the docstring
and FW §5.2's quotient-rule derivative (FW md:348). The order-0 / order-1
seeding in `Coefficients` (`u = Taylor1(T[y0,y1])`) is consistent.

**Intermittent?**: No.

**Confidence**: 0.9 that this is correct (listed as LOW to document the
verification, not to assert a bug).

### [LOW] FFW truncation-error citation/algebra — verified self-consistent; one citation-precision caveat

**Location**: `src/PadeStepper.jl:447-487` (`ffw_truncation_error`).

**Ground truth (cited)**: FFW md:76 `ε_{n+1} = c_{n+1} + Σ_{r=1..ν} b_r
c_{n+1-r}` with **unscaled** `c_k` and the **unscaled** Padé denominator
`b_r`; md:80 `T(h) = |ε_{n+1} h^{n+1} / p(h)|`.

**Code behavior**: The code works entirely in the rescaled variable. It
forms `ε̃_{n+1} = c̃_{n+1} + Σ_{r} P_u.b[r+1]·c̃_{n+1-r}` (lines 473-477)
and returns `|ε̃_{n+1} / a_rescaled(1)|` (line 486). I verified the
algebraic identity: with `c̃_k = h^k c_k` and the rescaled-Padé denominator
`b̃_r = h^r b_r`, the cross terms collapse —
`ε̃_{n+1} = h^{n+1}c_{n+1} + Σ_r (h^r b_r)(h^{n+1-r} c_{n+1-r}) =
h^{n+1}(c_{n+1}+Σ_r b_r c_{n+1-r}) = h^{n+1} ε_{n+1}`, and
`a_rescaled(1) = a(h) = p(h)`, so `|ε̃_{n+1}/a_rescaled(1)| = T(h)` exactly.
The 1-based index `c̃_{n+1-r}` at `coefs_u_scaled[n_int+2-r]` is correct,
and the loop lower bound (`n_int+2-νeff ≥ 2`) stays in range. The `q` law
(`ffw_rescale_q`, line 515) is `(k·Tol/T_h)^{1/(order+1)}`, matching md:89
eq. (2) verbatim.

**Mechanism**: None — the controller math is faithful.

Citation caveat (Law-1 doc hygiene, not numerics): `PadeStepper.jl:80-84`
cites "FW 2011 §3.2 line 396" for the `h' = h·t` substitution, but FW's §3.2
contains no explicit `h'=h·t` formula at that line (the rescaling is
standard Padé-Taylor practice; FW's own derivation in §5.2 eqs 5.4/5.5
works directly in `h`). The math is correct; only the line-cite is loose.

**Intermittent?**: No.

**Confidence**: 0.88 that the controller is correct.

## Areas verified correct

- **Rescaling power `c̃_k = h^k c_k`** (`PadeStepper.jl:342-350`,
  `VectorStepper.jl:286-294`) — the accumulator produces exactly
  `out[k] = c[k]·h^{k-1}`, i.e. `c̃_{k-1}=h^{k-1}c_{k-1}` in 1-based form,
  and evaluation at `t=1` recovers `Σ c_k h^k`. Both stepper copies are
  identical; no off-by-one in the power.

- **Value handoff `new_u = P_u(1)`** (`PadeStepper.jl:319`) — `t=1`
  corresponds to `h'=h`, so `P_u(1)=u(z+h)`. Cross-checked against the
  pinned oracle `u_5_1_1_at_05 = 4.0044646690030875`
  (`test/_oracle_padestepper.jl:34`).

- **Derivative handoff `new_up = P_u'(1)/h`** (`PadeStepper.jl:320`,
  `_evaluate_pade_deriv` 397-409) — the chain rule on `u(z+h')=P_u(h'/h)`
  gives `u'(z+h)=P_u'(1)/h`; the `/h_T` is present and correct. The
  quotient-rule implementation `(N'D − ND')/D²` with derivative coefficient
  `N'[k]=(k-1)·a[k]` (indices `k:-1:2`) is the correct derivative of a
  low-to-high polynomial whose `a[k]` multiplies `z^{k-1}`. FW md:348
  confirms "compute the derivative by the quotient rule of calculus." No
  off-by-one in reading the derivative.

- **No `up`-coefficient off-by-one risk** — `Coefficients.taylor_coefficients_2nd`
  returns only the **value** coefficients (`Coefficients.jl:202`); the
  stepper never reads `up` coefficients, deriving `u'` analytically from the
  value-Padé. The historically-tempting "re-Padé the differentiated
  coefficients `(k+1)c_{k+1}`" path (which would invite a `(k+1)` factorial
  off-by-one) is explicitly NOT taken (`PadeStepper.jl:96-122`).

- **`state.z` advances by exactly `h`** — `state.z = state.z + h_T`
  (`PadeStepper.jl:323`, `VectorStepper.jl:267`); `h_T = T(h)` is the same
  value used in the rescale and the `/h_T` derivative divide. No drift.

- **Continuity at joins** — `pade_step_with_pade!` reads the new state off
  the Padé at the step endpoint and writes `(z+h, u(z+h), u'(z+h))`; the
  next step seeds from those exact fields. In PathNetwork the chaining
  `z_cur,u_cur,up_cur = z_new,u_new,up_new` (`PathNetwork.jl:636`) and
  `parent_idx` linkage (`:634`) carry the endpoint forward verbatim — the
  start of step n+1 IS the endpoint of step n. The canonical real-`h` Padé
  stored for dense output (`PathNetwork.jl:605`, `_local_pade`) is a
  *separate* object used only for Stage-2 interpolation; it does not feed
  the value handoff, so the wedge-step (complex-`h`) vs canonical (real-`h`)
  split cannot introduce a value discontinuity at the join.

- **`d = 1` vector reduction** — `VectorStepper` feeds `shared_denominator_pade`,
  whose `d=1` case is bit-identical to the scalar `robust_pade(...,:svd)`
  (`SharedPade.jl:46-53`); the stepper's `Pᵢ(1)/Q(1)` at `d=1` is the same
  value the scalar pipeline produces (test VS.1.1, `test/vector_stepper_test.jl:66-88`).

- **Adjoint-vs-transpose (padeapprox.m:113 historical bug class)** — both
  the scalar `RobustPade.jl:443` and vector `SharedPade.jl:234` use
  `adjoint(C*D)` / `adjoint(A_full*D)` (conjugate-transpose), matching the
  *corrected* `padeapprox.m` `'`, NOT the erroneous `.'` plain transpose.
  The stepper itself contains no transpose/adjoint, and feeds complex data
  through the correctly-adjointed primitives.

- **Diagonal `(order÷2, order÷2)`** matches FW md:279 `(15,15)` at
  `order=30` (`PadeStepper.jl:313`, `VectorStepper.jl:241`).

- **`ffw_rescale_q` exponent** `(k·Tol/T_h)^{1/(order+1)}` matches FFW md:89
  eq. (2) (`PadeStepper.jl:515`); the `T_h==0 ⇒ Inf` guard is defensive,
  not reached on the accept path.

- **Complex-`h` direction preservation** in `adaptive_pade_step!` — `h *= q`
  with real `q ∈ (0,1]` (`PadeStepper.jl:576`) shrinks magnitude only,
  preserving the wedge direction; `meta.h_used = abs(h)` (`:583`) is the
  real magnitude used to re-seed (FFW md:93).
