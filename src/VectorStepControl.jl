"""
    PadeTaylor.VectorStepControl

Norm-based Jorba–Zou step-size selection for the **vector** Padé–Taylor
stepper.  v0.1's `StepControl.step_jorba_zou` chooses a step from one
scalar Taylor jet.  A vector ODE `y' = f(z, y)`, `y ∈ ℂ^d`, produces `d`
jets; this module picks **one** step size honouring the truncation
tolerance for the whole vector.

## Where this fits — the additive v0.2 vector lift

`VectorStepControl` is the vector analogue of `StepControl`'s
`step_jorba_zou`.  It is a *new* module; v0.1 `StepControl.jl` is
untouched and bit-identical — the additive architecture of
`docs/v0p2_plan.md` §"Architecture decision — additive".  The only edit
to a previously-shipped v0.2 module is one strictly-additive kwarg on
`VectorProblems.vector_solve_pade` (`step_policy`, default `:fixed`,
which preserves the V3b fixed-step behaviour byte-for-byte).

There is no vector analogue of `StepControl.step_pade_root` (the FW 2011
§3.1 pole-distance heuristic).  The v0.2 vector pipeline carries **one
shared denominator `Q`** per segment (ADR-0019): all `d` components blow
up at the *same* movable poles, so the FW 2011 pole-distance heuristic
applies unchanged to that single shared `Q` — there is nothing to
generalise.  ADR-0021 records this reuse.

## The norm-based generalisation — scalar `|c_k|` → vector `vnorm(c_k)`

The scalar `step_jorba_zou` (ported from `TaylorIntegration.jl`, see
`src/StepControl.jl`) is, at fixed order `p`,

    h = min over k ∈ {p-1, p} of (ε / |c_k|)^(1/k),

with `ε = eps_abs` if `eps_abs ≥ eps_rel·|c_0|`, else `ε = eps_rel·|c_0|`.

For a vector jet the `k`-th Taylor coefficient is no longer a scalar but
the `d`-vector `c_k = [jets[1][k+1], …, jets[d][k+1]]`.  We replace the
scalar magnitude `|c_k|` with a vector norm `vnorm(c_k)` (default the
2-norm; the `vnorm` kwarg lets a caller pick the `∞`-norm or any other):

    h = min over k ∈ {p-1, p} of (ε / vnorm(c_k))^(1/k),

with `ε = eps_abs` if `eps_abs ≥ eps_rel·vnorm(c_0)`, else
`ε = eps_rel·vnorm(c_0)`.

Everything else from the scalar module's discipline is preserved:

  - **zero-norm coefficient skip** — a `c_k` whose `vnorm` is zero
    contributes no constraint (an odd-only expansion truncated at an
    even order has `c_p = 0`);
  - **second-stepsize fallback** — when both `c_{p-1}` and `c_p` have
    zero norm the primary formula yields `Inf`; we fall back to the
    TI.jl `_second_stepsize` scan over `1 ≤ j ≤ p-2`, taking the
    maximum of `(1/vnorm(c_j))^(1/j)`;
  - **Rule-1 throws** — jets shorter than 2, an all-zero vector jet,
    ragged jets of unequal length, and an empty jet collection all
    throw `ArgumentError` with a `Suggestion` line rather than
    returning a meaningless `Inf` / `0` / `NaN`.

## Why a single vector norm — not per-component-min

`TaylorIntegration.jl`'s own vector `stepsize`
(`external/TaylorIntegration.jl/src/integrator/stepsize.jl:38-58`)
takes the **per-component minimum**: it runs the scalar selector on
each component independently and reduces with `min`.  We deliberately
do **not** mirror that.  The three candidates relate as follows.  Write
`s_k = [|jets[1][k+1]|, …, |jets[d][k+1]|]` for the vector of
component-magnitudes of the `k`-th coefficient:

  - **per-component-min** keys the step off `max_i |jets[i][k+1]|`
    — equivalently the `∞`-norm of `s_k` — taken before the `(·)^(1/k)`
    and the `min over k`.  (`(ε/a)^(1/k)` is decreasing in `a`, so the
    smallest per-component step picks the largest magnitude.)
  - **`∞`-norm** (`vnorm = v -> norm(v, Inf)`) keys off
    `max_i |jets[i][k+1]|` — the *same* quantity.  So per-component-min
    and the `∞`-norm choice coincide here; the `∞`-norm is the cleaner
    statement of the same conservative rule.
  - **2-norm** (the default) keys off `sqrt(Σ_i |jets[i][k+1]|²)`.
    Since `‖s_k‖_∞ ≤ ‖s_k‖_2 ≤ √d · ‖s_k‖_∞`, the 2-norm of each
    coefficient vector is **at least as large** as its `∞`-norm.

**When the 2-norm is the more conservative choice — and when it is
not.**  In *absolute* mode (`eps_abs ≥ eps_rel · vnorm(c_0)`, so
`ε = eps_abs` is norm-independent) the step `(ε / ‖s_k‖)^(1/k)` is
decreasing in `‖s_k‖`, hence the 2-norm step is **at most as large** as
the `∞`-norm step — the 2-norm is the more conservative default there
(test VSC.1.3).  In *relative* mode the ordering **does not hold**: `ε`
itself is rescaled by the *same* `vnorm` applied to `c_0` (the
`ε = eps_rel · vnorm(c_0)` branch below), so the norm enters both the
numerator and the denominator of `ε / ‖s_k‖`.  Counterexample (bead
`padetaylor-divo`, test VSC.1.3b): `c_0 = [1, 1]`, `c_1 = c_2 = [1, 0]`,
`eps_abs = 1e-12`, `eps_rel = 1e-10` gives `h_∞ = 1e-10` but
`h_2 = √2 · 1e-10` — the 2-norm step is *larger* by exactly
`‖c_0‖_2 / ‖c_0‖_∞ = √2`, because `c_0` is spread across components
while the later coefficients are concentrated in one.  In general the
2-norm step exceeds the `∞`-norm step by at most a factor `√d` in
relative mode, and never in absolute mode.  Jorba–Zou's own relative
error (their eq. 10, `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:601-604`)
normalises by the **sup norm** of `x_n`, so `vnorm = v -> norm(v, Inf)`
is the choice closest to the paper's relative mode.

We nonetheless default to the 2-norm because it is the truncation-error
norm the Jorba–Zou bound is naturally stated in (the local truncation
error of the vector jet is the Euclidean length of the per-component
error vector), and it does not discard the contribution of the
non-dominant components the way the `∞`-norm does.  The policy is
"one consistent norm for both `ε` and `‖s_k‖`", not "always the
smaller step".  A caller who wants the per-component-min idiom (or
the paper's sup-norm relative mode) passes `vnorm = v -> norm(v, Inf)`
and gets it exactly.

## The d = 1 reduction — the primary correctness oracle

For `d = 1` every coefficient vector `c_k` is a length-1 vector, and
`norm([x]) == abs(x)` for the 2-norm (and the `∞`-norm).  The formula,
the `ε` resolution, the zero-skip, and the fallback all then collapse
*exactly* onto `StepControl.step_jorba_zou`.  Hence
`vector_step_jorba_zou([jet], ε)` is bit-identical to
`step_jorba_zou(jet, ε)` — test VSC.1.1 asserts this on several jets,
including the canonical `exp` jet, and it is the tightest available
check.

## References

  - `src/StepControl.jl` — the scalar `step_jorba_zou` this module
    generalises (the TI.jl-ported fixed-order Jorba–Zou formula); the
    `d = 1` correctness oracle.
  - Jorba & Zou 2005, Experimental Mathematics 14, §3.3.1 eq. 11 —
    `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/
    JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-645`.
  - `external/TaylorIntegration.jl/src/integrator/stepsize.jl:38-90`
    — the canonical Julia vector `stepsize` (per-component-min idiom),
    contrasted above.
  - `docs/v0p2_plan.md` — §"The 10 phases", V3c row.
  - ADR-0021 — the norm-based vector step-control decision (the
    per-component-min vs 2-norm vs ∞-norm comparison; the
    `step_pade_root` reuse note).
"""
module VectorStepControl

using LinearAlgebra: norm

export vector_step_jorba_zou

# -----------------------------------------------------------------------------
# Norm-based vector Jorba-Zou step
# -----------------------------------------------------------------------------

"""
    vector_step_jorba_zou(jets::AbstractVector{<:AbstractVector{T}},
                          eps_abs::Real; eps_rel::Real = eps_abs,
                          vnorm = norm) where {T}
        -> real-typed step h

Norm-based generalisation of `StepControl.step_jorba_zou` to a vector
ODE.  `jets` is the `d`-element collection of per-component Taylor jets
(`jets[i] = [c⁽ⁱ⁾_0, …, c⁽ⁱ⁾_p]`, all of common length `p + 1`).  The
`k`-th Taylor coefficient is the `d`-vector
`c_k = [jets[1][k+1], …, jets[d][k+1]]`, and the returned step is

    h = min over k ∈ {p-1, p} of (ε / vnorm(c_k))^(1/k),

where `ε = eps_abs` if `eps_abs ≥ eps_rel · vnorm(c_0)`, else
`ε = eps_rel · vnorm(c_0)`.  `vnorm` defaults to the 2-norm; pass
`vnorm = v -> norm(v, Inf)` for the `∞`-norm (the per-component-min
idiom — see the module docstring).  Coefficient vectors whose `vnorm`
is zero are skipped; when both `c_{p-1}` and `c_p` have zero norm the
TI.jl second-stepsize scan over `1 ≤ j ≤ p-2` is used as a fallback.

For `d = 1` this reduces bit-identically to
`StepControl.step_jorba_zou(jets[1], eps_abs; eps_rel)`.

Throws `ArgumentError` (Rule 1 — fail loud) if `jets` is empty, if any
jet has fewer than 3 entries (`p ≥ 2`; a length-2 jet puts `k = 0` in
the candidate set and `x^(1/0)` silently yields `h = 0`), if
the jets are ragged (unequal length — a malformed vector jet), or if
every coefficient vector has zero norm (no nonzero coefficient to
estimate a step from).
"""
function vector_step_jorba_zou(jets::AbstractVector{<:AbstractVector{T}},
                               eps_abs::Real; eps_rel::Real = eps_abs,
                               vnorm = norm) where {T}
    d = length(jets)
    d ≥ 1 || throw(ArgumentError(
        "vector_step_jorba_zou: jets is empty — a zero-component system " *
        "has no Taylor jet to choose a step from. " *
        "Suggestion: pass a non-empty collection of per-component jets."))

    len = length(jets[1])
    len ≥ 3 || throw(ArgumentError(
        "vector_step_jorba_zou: each jet must have length ≥ 3 (got " *
        "$len); the formula uses the last two Taylor coefficients with " *
        "p ≥ 2 — for p = 1 the k = 0 candidate silently gives h = 0. " *
        "Suggestion: pass per-component jets of order ≥ 2."))
    for i in 2:d
        length(jets[i]) == len || throw(ArgumentError(
            "vector_step_jorba_zou: ragged jets — jet 1 has length $len " *
            "but jet $i has length $(length(jets[i])); a vector Taylor " *
            "jet must have one common order across all components. " *
            "Suggestion: truncate every component's jet to a common order."))
    end

    p = len - 1
    R = float(real(T))

    # ε resolution per the scalar StepControl (TI.jl stepsize.jl:27-31),
    # with the scalar magnitude |c_0| replaced by the vector norm of the
    # leading-coefficient vector c_0 = [jets[1][1], …, jets[d][1]].
    c0_norm = vnorm(_coef_vector(jets, 0, d))
    ε = eps_abs ≥ eps_rel * c0_norm ? R(eps_abs) : R(eps_rel * c0_norm)

    h = typemax(R)
    for k in (p - 1, p)
        nrm = vnorm(_coef_vector(jets, k, d))
        iszero(nrm) && continue
        h = min(h, (ε / nrm)^(R(1) / k))
    end

    if !isinf(h)
        return h
    end

    # Fallback: the TI.jl `_second_stepsize` scan (stepsize.jl:80-90),
    # lifted to the vector norm.  Scan j = 1 .. p-2 for the largest h
    # such that vnorm(c_j) · h^j = 1, i.e. h = (1/vnorm(c_j))^(1/j).
    h2 = zero(R)
    for j in 1:(p - 2)
        nrm = vnorm(_coef_vector(jets, j, d))
        iszero(nrm) && continue
        h2 = max(h2, (R(1) / nrm)^(R(1) / j))
    end
    iszero(h2) && throw(ArgumentError(
        "vector_step_jorba_zou: every Taylor coefficient of the vector " *
        "jet has zero norm; cannot estimate a step. " *
        "Suggestion: check the integrator's right-hand side and initial " *
        "condition — a constant zero solution needs no stepping."))
    return h2
end

"""
    _coef_vector(jets, k::Integer, d::Integer) -> Vector

Assemble the `k`-th Taylor coefficient as the `d`-vector
`[jets[1][k+1], …, jets[d][k+1]]` — the cross-component slice the vector
norm is taken over.  A fresh small vector per `(k)` call; the jet count
`d` is the vector ODE dimension, modest in the v0.2 targets.
"""
_coef_vector(jets, k::Integer, d::Integer) =
    [jets[i][k + 1] for i in 1:d]

end # module VectorStepControl
