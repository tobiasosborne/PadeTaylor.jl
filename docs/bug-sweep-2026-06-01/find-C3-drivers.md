# Bug sweep C3 — Problems + VectorProblems: solve loop & dense-output segment selection

Auditor pass 2026-06-01. Read-only. Hunt: an intermittent discontinuity
introduced by a mis-transcribed equation / off-by-one / boundary-comparison
bug in the driver loop or the dense-output segment-ownership logic.

## Area

`src/Problems.jl` (scalar `solve_pade` + `PadeTaylorSolution` dense callable)
and `src/VectorProblems.jl` (`vector_solve_pade` + `VectorPadeTaylorSolution`
dense callable). Specifically:

- the `while state.z < z_end` solve loop and its final-step clamp;
- the dense-output segment-selection scan (which stored segment Padé owns a
  query `z`), and the strict-vs-nonstrict comparison at segment joins;
- per-segment Padé store indexing vs the breakpoints array (off-by-one);
- buffer aliasing / order dependence that could produce intermittency.

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-166`
  — Stage 1/Stage 2 path strategy. Line 164: "Once we reach **within a
  distance of h** from the target point, this path is halted. Along such a
  path, we store for each discrete point not just u and u' but also the Padé
  coefficients." This is the ground truth for the per-discrete-point Padé
  store and for the fact that FW stops *within h* (does NOT clamp to land
  exactly on a target).
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:216,220-222`
  — dense evaluation at fine-grid nodes via the local stored Padé; confirms
  "one local Padé per discrete point, evaluate at nearby nodes" model.
- `docs/adr/0001-four-layer-architecture.md:13-24,71-76` — the layer
  contract: `Problems` owns the trajectory and the per-segment store; the
  cross-layer interface is typed structs; the dense callable is "a small
  Horner loop".
- Cross-checked the steppers to ground-truth the join-continuity claim and
  buffer freshness:
  `src/PadeStepper.jl:316-326` (`new_u = _evaluate_pade(P_u, 1)`,
  `state.z = state.z + h_T`),
  `src/VectorStepper.jl:255-269` (`new_y[i] = P_i(1)/Q(1)`,
  `state.z = state.z + h_T`, returns freshly-allocated `numerators,
  denominator`),
  `src/PadeStepper.jl:367-403` (`_evaluate_pade` / `_evaluate_pade_deriv`
  Horner low-to-high, derivative coeff `(k-1)*a[k]`).

## Findings

### [LOW] Final-step clamp relies on `a + (b - a) == b`; sub-ulp residue can spawn a degenerate final micro-segment

- **Location**: `src/Problems.jl:195-208` (loop + `h_step = min(h_max_T,
  z_end - state.z)` at line 202; `state.z = state.z + h_T` at
  `src/PadeStepper.jl:323`). Mirror in `src/VectorProblems.jl:262-291`
  (`h_ceiling = min(h_T, z_end - state.z)` at line 276; `state.z = state.z +
  h_T` at `src/VectorStepper.jl:267`).
- **Ground truth (cited)**: FW2011 §3.1 line 164 stops *within h* of the
  target, never claiming exact landing; the "land exactly on z_end" behaviour
  is a v1 design choice (documented `src/VectorProblems.jl:218-219`), not a
  transcription. There is no reference equation requiring `state.z == z_end`
  bit-exactly.
- **Code behavior**: The last step sets `h_step = z_end - state.z` and then
  `state.z := state.z + h_step`. In floating point `a + (b - a)` is not
  guaranteed to equal `b`, so `state.z` may end at `z_end ± 1 ulp`. If it lands
  fractionally *below* `z_end`, the `while state.z < z_end` guard fires once
  more and produces a sub-ulp final segment with `h[end] ≈ 0`.
- **Mechanism (intermittent)**: A near-zero `h[end]` makes the dense callable's
  `t = (z - z[k])/h_k` (Problems.jl:232 / VectorProblems.jl:326) divide by an
  almost-zero denominator for queries that fall in that micro-segment,
  amplifying rounding. Whether the extra micro-segment appears depends on the
  exact FP residue of `state.z_old + (z_end - state.z_old)`, which varies with
  `z_start`, `h_max`, and accumulated step rounding — i.e. it is
  data-dependent / intermittent. In practice the micro-segment is sub-ulp wide,
  so essentially no query z (itself rounded to T) lands strictly inside it, and
  it does not produce a *visible* discontinuity. Reported for completeness; the
  clamp intent is correct and there is NO overshoot-then-extrapolate bug (the
  special-focus hypothesis does not hold — see Areas verified correct).
- **Intermittent?**: Yes (FP-residue / accumulation dependent), but negligible
  magnitude.
- **Confidence**: 0.2 that this ever manifests as an observable artifact; 0.9
  that the FP-residue reasoning itself is correct.

### [LOW] Vector dense scan and guards compare only `real(z)`; latent ambiguity if segments ever share a real part

- **Location**: `src/VectorProblems.jl:312-323` (DomainError guards on
  `real(...)` and the scan `while k < length(sol.h) && real(z_T) >
  real(sol.z[k+1])`). Contrast the scalar callable `src/Problems.jl:219-229`
  which compares the full (real) value directly.
- **Ground truth (cited)**: ADR-0001 (`docs/adr/0001-four-layer-architecture.md:13-24`)
  describes `Problems` as owning a 1-D trajectory parameterised by step length;
  the vector driver enforces a real positive `h` (`src/VectorProblems.jl:244`).
  There is no reference prescribing how to order complex breakpoints.
- **Code behavior**: Segment ownership is decided purely by `real(z)`. For the
  supported trajectories — real positive `h`, breakpoints `z[k] = z_start +
  Σh`, all sharing `imag(z_start)` with strictly increasing real parts — the
  ordering is exact and unambiguous, so this is correct *today*.
- **Mechanism (intermittent)**: It would only become a discontinuity source if
  a future complex path-network produced two breakpoints with equal real part
  but different imaginary parts; then `real(z_T) > real(z[k+1])` would pick the
  wrong segment for some queries, jumping at the (mis-ordered) join. Not
  reachable in the current pipeline (no complex stepping reaches this driver).
- **Intermittent?**: Not in the current code path (latent only).
- **Confidence**: 0.15 (latent design limitation, not an active bug). Worth a
  comment/bead when a complex vector path-network ships.

## Areas verified correct

- **Solve-loop clamp does NOT overshoot z_end** (the headline suspicion is
  refuted). `solve_pade` passes `h_step = min(h_max_T, z_end - state.z)`
  (`src/Problems.jl:202`); `vector_solve_pade` passes
  `min(h_T, z_end - state.z)` (`src/VectorProblems.jl:276`). The stepper
  advances `state.z` by exactly that `h` (`src/PadeStepper.jl:323`,
  `src/VectorStepper.jl:267`). The last segment therefore ends at `z_end`
  (modulo the sub-ulp caveat above), never beyond it, so the dense callable
  never extrapolates the last segment past `t = 1`. The special-focus
  "overshoot → extrapolate → blow up near z_end" mechanism is not present.

- **Segment-ownership comparison at joins is consistent (no jump).** Scalar
  scan `src/Problems.jl:226-229` and vector scan `src/VectorProblems.jl:320-323`
  both advance `k` only while `z > z[k+1]` (strict). At a join `z == z[k+1]`
  the strict test is false, so the *left* segment `k` owns the join and is
  evaluated at `t = (z[k+1]-z[k])/h[k] = 1`. The stepper *defines* the next
  breakpoint as that same `t=1` value (`new_u = _evaluate_pade(P_u, 1)`,
  `src/PadeStepper.jl:319`; `new_y[i] = P_i(1)/Q(1)`,
  `src/VectorStepper.jl:264`), and the next segment's Padé reproduces it at
  `t=0` (its `c_0` is the carried-in state). Left-at-`t=1` equals
  right-at-`t=0` by construction ⇒ continuous joins. The convention is uniform
  (every interior join to the left segment; final endpoint owned by the last
  segment), with no fallthrough.

- **No off-by-one between the per-segment Padé/`h` store and the breakpoint
  array.** Driver invariant: `z` and `y` start with one element
  (`src/Problems.jl:189-190`, `src/VectorProblems.jl:256-257`) and each
  iteration pushes one `z`, one `y`, one `h`, one Padé
  (`src/Problems.jl:204-207`, `src/VectorProblems.jl:287-290`). After `n`
  steps: `length(z)=length(y)=n+1`, `length(h)=length(pade)=n`. The dense
  callable reads `h[k]`, `z[k]`, `z[k+1]`, `pade[k]` with `1 ≤ k ≤ n`; the max
  index is `z[n+1]=z[end]`, in-bounds. Matches the documented structure
  (`src/Problems.jl:143-147`, `src/VectorProblems.jl:182-204`) exactly. The
  `@inbounds` annotations do not mask any out-of-range access.

- **Padé local variable matches the dense rescaling.** The stepper builds the
  Padé around `z_old = z[k]` in `t = h'/h` (`src/PadeStepper.jl:280-281`);
  the dense callable evaluates at `t = (z - z[k])/h[k]`
  (`src/Problems.jl:232`, `src/VectorProblems.jl:326`) — the same variable.
  The `u'` chain-rule factor `1/h[k]` (`src/Problems.jl:235`) matches the
  stepper's `P_u'(1)/h` recovery (`src/PadeStepper.jl:320`), consistent with
  FW2011's `h^k` rescaling rationale (module docstring `src/Problems.jl:51-65`).

- **No buffer aliasing / order dependence in the store.** Vector state is
  copied on push (`copy(state.y)`, `src/VectorProblems.jl:288`; initial
  `copy(prob.y0)`, line 257). The stepper returns freshly-allocated
  `numerators, denominator` each call (`shared_denominator_pade(...)`,
  `src/VectorStepper.jl:242`, returned at line 269), wrapped in a new
  `SharedPadeApproximant` per segment (line 290) — no shared mutable buffer
  across segments, so stored segments cannot alias the last step's
  coefficients.

- **Out-of-window guards are correct and fail-loud.** `z < z[1]` / `z > z[end]`
  throw `DomainError` with a suggestion (`src/Problems.jl:219-222`,
  `src/VectorProblems.jl:312-315`), satisfying Rule 1; queries exactly at
  `z_start` (`t=0`, segment 1) and `z_end` (`t=1`, last segment) are accepted
  and land on the correct terminal segment because the scan caps `k` at
  `length(h)` via the `k < length(sol.h)` guard.

- **Horner coefficient order and derivative factor are consistent.**
  `_evaluate_pade` / `_eval_poly` Horner low-to-high
  (`src/PadeStepper.jl:367-375`, `src/VectorProblems.jl:339-345`); derivative
  coefficient is `(k-1)*a[k]` for the 1-indexed `z^(k-1)` coefficient
  (`src/PadeStepper.jl:402-403`) — the correct derivative relation, no
  factorial/coefficient confusion.
