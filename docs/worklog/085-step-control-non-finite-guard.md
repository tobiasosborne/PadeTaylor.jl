# 085 — Step control returned NaN / 0 on non-finite coefficients (bead xbxm)

2026-09-29. Found by a read-only scout, confirmed by reading, fixed here.

## The defect

`step_jorba_zou` (`src/StepControl.jl`) and `vector_step_jorba_zou`
(`src/VectorStepControl.jl`) compute

    h = min over k ∈ {p-1, p} of (ε / |c_k|)^(1/k)

and leave through `if !isinf(h) return h end`. Two inputs defeat that exit:

- a `NaN` coefficient: `iszero(NaN)` is false, `min(h, NaN)` is `NaN`, and
  `!isinf(NaN)` is true, so `NaN` is returned as the step;
- an `Inf` coefficient: `ε / Inf = 0`, so `h = 0` is returned.

`src/VectorProblems.jl` then takes `min(h_jz, h_T, gap_mag)` with no
finiteness check, so the value reaches the stepper. Both docstrings promised
the opposite ("never returning a meaningless `Inf` or `0`").

## The fix

Both functions now throw `ArgumentError` naming the first offending index
(and, in the vector case, the component), with a suggestion. Every
coefficient is checked, not only the trailing two, because the fallback scan
reads indices `1 … p-2`. Zero coefficients remain legal.

## Tests

- `test/stepcontrol_test.jl` 4.1.7: `NaN`, `Inf`, `-Inf` in each of four
  slots; a complex `NaN`; a `BigFloat` `NaN`; and a jet with a zero entry that
  must still return the closed-form step.
- `test/vector_step_control_test.jl` VSC.1.5: `NaN` and `Inf` in each of three
  slots of the second component.

Observed after the fix: StepControl 39/39, VectorStepControl 81/81.

## Mutation proof

Mutation E: both guards replaced by `true || throw(...)`.

- StepControl: 4.1.7 RED, failures at `test/stepcontrol_test.jl:120` and
  `:121` (the captured output was cut at its first ten failure lines, so the
  total was not recorded).
- VectorStepControl: 75 passed, 6 failed, all in VSC.1.5 at
  `test/vector_step_control_test.jl:232`.

Both files restored from backup; re-run 39/39 and 81/81.

## Not done

The tolerance arguments (`eps_abs`, `eps_rel`) are not validated: a `NaN` or
non-positive tolerance is a separate caller error and is left as is.
