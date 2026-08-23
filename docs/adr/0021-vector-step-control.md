# ADR-0021 — Vector step control: norm-based Jorba–Zou

**Status**: Accepted (2026-05-20) | **Bead**: `padetaylor-0ln.10` (V3c —
`VectorStepControl.jl` norm-based vector step selector) | **Plan**:
`docs/v0p2_plan.md` row V3c.

## Decision

`VectorStepControl.vector_step_jorba_zou(jets, eps_abs; eps_rel, vnorm)`
— the v0.2 vector step-size selector for first-order systems
`y' = f(z, y)`, `y ∈ ℂ^d` — generalises `StepControl.step_jorba_zou`
from one scalar Taylor jet to `d` jets by replacing the scalar
coefficient magnitude `|c_k|` with a **vector norm** `vnorm(c_k)` of the
`k`-th coefficient *vector* `c_k = [jets[1][k+1], …, jets[d][k+1]]`:

    h = min over k ∈ {p-1, p} of (ε / vnorm(c_k))^(1/k),

with `ε = eps_abs` if `eps_abs ≥ eps_rel·vnorm(c_0)`, else
`ε = eps_rel·vnorm(c_0)`.  The `vnorm` keyword defaults to the **2-norm**;
a caller passes `vnorm = v -> norm(v, Inf)` for the `∞`-norm.

`vector_solve_pade` gains a strictly-additive kwarg
`step_policy::Symbol = :fixed`.  `:fixed` (default) is byte-identical to
the V3b fixed-step driver; `:jorba_zou` calls `vector_step_jorba_zou` at
each step and uses the supplied `h` as a **ceiling**:
`h_step = min(h_jz, h, z_end − z)`.

## Context

v0.1's `StepControl.step_jorba_zou` picks a step size from *one* scalar
Taylor jet (the TI.jl-ported fixed-order Jorba–Zou 2005 §3.3.1 eq. 11
formula; see `src/StepControl.jl`).  A vector ODE has `d` jets, so v0.2
needs a selector that picks a *single* step honouring the truncation
tolerance for the *whole* vector — the design question this ADR settles
is **how the `d` per-component constraints are combined into one step**.

The choice propagates into the `:jorba_zou` adaptive policy on the V3b
driver and into every later v0.2 adaptive solve (`P_I⁽²⁾`, Noumi–Yamada),
so it is recorded here before those beads build on it.

## Alternatives considered

Write `s_k = [|jets[1][k+1]|, …, |jets[d][k+1]|]` for the vector of
component-magnitudes of the `k`-th Taylor coefficient.  `(ε/a)^(1/k)` is
**decreasing** in `a`, so a *larger* coefficient measure ⇒ a *smaller*
(more conservative) step.

### A — per-component minimum

Run the scalar `step_jorba_zou` on each component independently, reduce
with `min`.  This is exactly what `TaylorIntegration.jl`'s own vector
`stepsize` does (`external/TaylorIntegration.jl/src/integrator/
stepsize.jl:38-58`).  Because `min over i of (ε/|jets[i][k+1]|)^(1/k)` is
`(ε / max_i |jets[i][k+1]|)^(1/k)`, the per-component-min keys the step
off `max_i |jets[i][k+1]| = ‖s_k‖_∞` — the `∞`-norm of the
coefficient-magnitude vector.

### B — `∞`-norm (chosen as the opt-in non-default)

`vnorm = v -> norm(v, Inf)` keys the step off `max_i |jets[i][k+1]|` —
the **same quantity** as per-component-min.  So alternative A and the
`∞`-norm choice **coincide**; the `∞`-norm is the cleaner single-norm
statement of the per-component-min rule.  A caller who wants TI.jl's
idiom passes `vnorm = v -> norm(v, Inf)` and gets it exactly — A is not
lost, it is a one-keyword opt-in.

### C — 2-norm (chosen as the default)

`vnorm = norm` keys the step off `‖s_k‖_2 = sqrt(Σ_i |jets[i][k+1]|²)`.

## Why the 2-norm default

Norm equivalence gives `‖s_k‖_∞ ≤ ‖s_k‖_2 ≤ √d · ‖s_k‖_∞`.  In
**absolute mode** (`ε = eps_abs`, norm-independent) the step
`(ε / ‖s_k‖)^(1/k)` is *decreasing* in the coefficient measure, so the
2-norm — being **at least as large** as the `∞`-norm — yields a step
that is **at most as large**: there the 2-norm is the more conservative
choice (test VSC.1.3).

**This ordering does NOT hold in relative mode** (corrected 2026-08-23,
bead `padetaylor-divo`; the original text of this section claimed it
held unconditionally).  `ε` is resolved as `eps_rel · vnorm(c_0)`
(`src/VectorStepControl.jl`, the `ε = …` line after `c0_norm`), i.e.
the same norm rescales the numerator, and `‖c_0‖_2 / ‖c_0‖_∞` can be as
large as `√d`.  Hand-verified counterexample, pinned by test VSC.1.3b:
`c_0 = [1, 1]`, `c_1 = c_2 = [1, 0]`, `eps_abs = 1e-12`,
`eps_rel = 1e-10` → `h_∞ = 1e-10`, `h_2 = √2 · 1e-10` (measured ratio
`1.4142135623730954`).  So the honest statement is: the 2-norm step is
never larger than the `∞`-norm step in absolute mode, and at most `√d`
times larger in relative mode.  Jorba–Zou's own relative-error estimate
(eq. 10, `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:601-604`)
divides by the **sup norm** of `max{|x_n|, |ẋ_n|}`; a caller who wants
the paper's relative mode passes `vnorm = v -> norm(v, Inf)`.

The policy stays the 2-norm: the design rule is *one consistent norm on
both sides of `ε / ‖s_k‖`*, not "always the smaller step".  Two reasons
make it the right default:

  1. **It is the natural truncation-error norm.**  The local truncation
     error of the vector jet is the Euclidean length of the
     per-component error vector; the Jorba–Zou bound `‖c_k‖·h^k ≈ ε` is
     most honestly stated in the 2-norm.

  2. **It does not discard the non-dominant components.**  The `∞`-norm
     (and per-component-min) is blind to every component but the
     largest: a system whose components are all comparable in magnitude
     but individually below the dominant one gets a step keyed off the
     dominant component alone.  The 2-norm folds every component's
     contribution into the measure, so a vector with many moderate
     components is — correctly — treated as carrying more truncation
     error than a single moderate component would.

The `∞`-norm remains available for callers who specifically want the
TI.jl per-component-min behaviour (e.g. cross-validating against a
TI.jl reference run).

## The `d = 1` reduction — the primary correctness oracle

For `d = 1` every coefficient vector `c_k = [jets[1][k+1]]` is a
length-1 vector, and `norm([x]) == abs(x)` for **both** the 2-norm and
the `∞`-norm.  The formula, the `ε` resolution, the zero-norm skip, and
the second-stepsize fallback then collapse *exactly* onto
`StepControl.step_jorba_zou`.  Hence `vector_step_jorba_zou([jet], ε)`
is bit-identical to `step_jorba_zou(jet, ε)` — test VSC.1.1 asserts this
on several jets (including the canonical `exp` jet), and it is the
tightest available check.  This reduction property is *why* the
norm-based design is safe: at `d = 1` it inherits the scalar selector's
three-source-consensus pinned value (`test/_oracle_stepcontrol.jl`).

## No vector `step_pade_root` — the shared-`Q` reuse note

`StepControl` carries a second selector, `step_pade_root` — the FW 2011
§3.1 pole-distance heuristic that caps the step at the distance to the
closest forward root of the Padé *denominator*.  This bead deliberately
ships **no vector generalisation** of it, and that is a correctness
statement, not a deferral.

The v0.2 vector pipeline carries **one shared denominator `Q`** per
integration segment (ADR-0019; `src/SharedPade.jl`): all `d` components
of a vector meromorphic solution blow up at the *same* movable poles, so
a single consistent `Q` per segment is the honest representation.
`step_pade_root` already operates on a single denominator polynomial
`P.b` — it computes the roots of one polynomial and projects them onto
the step direction.  Given the shared `Q`, `step_pade_root` therefore
applies **unchanged** to the vector case: feed it the segment's shared
`Q` and it returns the pole-distance cap for the whole vector system at
once.  There is no per-component denominator to reconcile, hence nothing
to generalise.  The vector path-network bead (V7) will call the existing
`step_pade_root` on each segment's shared `Q` directly.

## The `step_policy` cap semantics

`:jorba_zou` treats the user-supplied `h` as a **ceiling**, not a target:
`h_step = min(h_jz, h, z_end − z)` where `h_jz` is the norm-based
Jorba–Zou step.  So `h` is the largest step the integrator is *allowed*
to take; the adaptive policy may take a smaller one when the local jet
demands it, and never a larger one.  This makes `:jorba_zou` a strict
refinement of `:fixed`: switching policies can only shrink steps, never
grow them past the user's stated bound — important for a user who chose
`h` to stay within a known pole-free radius.

## Consequences

- `VectorStepControl.jl` is a new module (61 effective LOC, well under
  the Rule-6 200-LOC cap; literate top docstring).  v0.1 `StepControl.jl`
  is **untouched and bit-identical** — the additive architecture of
  `docs/v0p2_plan.md` §"Architecture decision — additive".
- The only edit to a previously-shipped v0.2 module is the additive
  `step_policy::Symbol = :fixed` kwarg on `VectorProblems.vector_solve_pade`.
  The default `:fixed` preserves the V3b fixed-step behaviour
  byte-for-byte; the V3b test suite (`test/vector_problems_test.jl`,
  72/72) is re-confirmed GREEN unchanged.
- Generic in `T`: the step is returned at `float(real(T))`, so a
  `BigFloat` / `Arb` vector jet selects a step at the matching real
  precision.
- Fail-loud (Rule 1): an empty jet collection, jets shorter than 2,
  ragged jets of unequal length, and an all-zero vector jet all throw
  `ArgumentError` with a `Suggestion` line.  An unknown `step_policy`
  symbol throws.

## References

- `src/VectorStepControl.jl` — the implementation (literate, 61 LOC).
- `src/StepControl.jl` — the scalar `step_jorba_zou` this module
  generalises (the TI.jl-ported fixed-order Jorba–Zou formula) and the
  `step_pade_root` pole-distance heuristic reused unchanged on the
  shared `Q`; the `d = 1` correctness oracle.
- `src/VectorProblems.jl` — the V3b driver gaining the `step_policy`
  kwarg.
- `test/vector_step_control_test.jl` — the V3c suite (69 assertions,
  GREEN) + the VSC.1.6 mutation-proof record (M1/M2/M3, all bit).
- Jorba & Zou 2005, Experimental Mathematics 14, §3.3.1 eq. 11 —
  `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/
  JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-645`.
- `external/TaylorIntegration.jl/src/integrator/stepsize.jl:38-90` —
  the canonical Julia vector `stepsize` (per-component-min idiom),
  contrasted in "Alternatives considered".
- FW 2011 §3.1 — the pole-distance heuristic reused unchanged.
- ADR-0019 — the shared-denominator Padé; the keystone the
  `step_pade_root` reuse note rests on.
- ADR-0020 — the `Vector{Taylor1{T}}` / `Vector{Vector{T}}` jet shape
  this selector consumes.
- `docs/v0p2_plan.md` — §"The 10 phases", V3c row; §"Architecture
  decision — additive".
- CLAUDE.md Rules 1, 4, 6, 9, 10.
