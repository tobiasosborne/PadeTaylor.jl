"""
    PadeTaylor.StepControl

Step-size selection for the Taylor-Padé stepper.  Two strategies live
here, matching the two flavours of step control discussed in the
literature for analytic IVPs:

  - `step_jorba_zou` — the Jorba-Zou 2005 §3.3.1 first-step-size
    control, restricted to a fixed truncation order `p` (this project
    uses the FW 2011 default `p = 30`).  This is the "default"
    selector: it asks "given the last two Taylor coefficients of the
    local jet, how far can I march before the truncation error
    exceeds `ε`?"  The formula is local to the Taylor jet and ignores
    nearby singularities of the actual solution — the Padé layer is
    what handles the singularities.

  - `step_pade_root` — the FW 2011 §3.1 pole-distance heuristic.
    Given the Padé approximant of the local jet, the closest root of
    the *denominator* polynomial in the forward direction tells us
    how far we can step before crossing the singularity that the
    Padé layer flagged.  This is the "alternative" / safety selector
    you cap the Jorba-Zou step against in `PadeStepper`.

## The Jorba-Zou step formula — three-source consensus

DESIGN.md §4 / HANDOFF.md cite a step formula

```
h = (ρ / e²) · exp(-0.7 / (p - 1))
```

attributed to "Jorba-Zou §3.2.1 eq. 3-8".  That formula is wrong on
two counts and is **not** what we implement.  The triangulation:

  1. The actual paper (`references/markdown/JorbaZou2005_taylor_IVP_
     package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:
     613-635`) writes the §3.3.1 first step-size control as

         ρ̄_j = |x_n^[j]|^(-1/j),
         ρ_n = min{ρ̄_{p_n-1}, ρ̄_{p_n}, ρ̂_{p_n-1}, ρ̂_{p_n}},
         h_n = ρ_n / e²,                                  (eq. 11)

     where `p_n = -½ ln ε_n + 1` is chosen *adaptively* from the
     accuracy `ε_n` (eq. 11 same line).  The `e²` denominator only
     makes sense at that adaptive `p_n`; at any *fixed* `p` (we fix
     `p = 30` per FW 2011) the `e²` is implicitly absorbed into the
     `ε` substitution.  See the worked algebra in the test file:
     `test/stepcontrol_test.jl:13-26`.

  2. The canonical Julia implementation (`external/TaylorIntegration.
     jl/src/integrator/stepsize.jl:17-35`) drops the `/e²` and
     `p`-adaptation and writes the equivalent fixed-order form

         h = min over k ∈ {ord-1, ord} of (ε / |x[k]|)^(1/k).

     Grep both the paper markdown and TI.jl's source for `0.7` or
     `safety`: zero matches.  The DESIGN.md `0.7` constant is a
     hallucination at the time of the spec write-up.

  3. mpmath at 50 dps and Mathematica at 50 dps both reproduce the
     TI.jl number for the canonical case `c[k] = 1/k!`, `p = 30`,
     `ε = 1e-12` to 47 decimal digits:
     `4.50120637033898607690318848315848021108523516609`.  Any other
     formula for the same input gives a different leading digit.

We therefore implement the TI.jl formula verbatim (with the function
signature adapted to `Vector{T}` / `eps_abs, eps_rel` rather than
`Taylor1{U}` / `absepsilon, relepsilon`), and the test suite in
`test/stepcontrol_test.jl` pins the result against the three-source
consensus oracle in `test/_oracle_stepcontrol.jl`.

## ε resolution and the TI.jl "second stepsize" fallback

The eps comparison `eps_abs ≥ eps_rel · |c₀|` (TI.jl `stepsize.jl:27`)
selects between absolute and relative tolerance modes.  When *both*
candidate coefficients `c[p-1]`, `c[p]` are zero (e.g. odd-only
expansions of `sin` truncated at the wrong parity), the primary
formula yields `Inf`; we fall back to the TI.jl `_second_stepsize`
trick of scanning lower-order coefficients for the largest `h` such
that `|c[j]| · h^j = 1` (TI.jl `stepsize.jl:77-89`).  If no nonzero
coefficient exists at all, throw — Rule 1 (fail fast, never return
garbage like `Inf` or `NaN` to the caller).

## The FW 2011 pole-distance heuristic

`step_pade_root` implements the FW 2011 §3.1 idea: the Padé
approximant `r(z) = a(z) / b(z)` exposes its (numerical) singularities
as the roots of `b(z)`.  When stepping from `z_current` toward
`target`, the largest safe step is bounded by the distance to the
*forward-most* such root projected onto the step direction.  Roots
behind us, on us, or sideways (perpendicular to the path) carry no
constraint.  Concretely:

  1. Build the polynomial from `P.b` via `Polynomials.Polynomial`
     (same convention: `b[1]` is the constant term, `b[end]` the
     leading coefficient).
  2. Compute its roots.
  3. Project each root `r` onto the unit step direction:
     `t = Re((r - z_current) · conj(unit))`.  Discard `t ≤ 0`.
  4. Cap the step at `min(min_t, |target - z_current|)`.

A constant denominator (`length(P.b) ≤ 1`) means no poles, so the
step is trivially capped at `|target - z_current|`.

## Polynomials.jl element-type caveat

`Polynomials.Polynomial(b)` and `Polynomials.roots(p)` work for
`b::Vector{Float64}` and `b::Vector{Complex{Float64}}`.  For the
arb-precision path (`T = Arblib.Arb`) the chain
`Polynomials.roots ∘ Polynomials.Polynomial` is *not* yet validated —
the v1 element type for the Padé-root step is `Float64`/`Complex{
Float64}` (matching the SVD layer's caveat in ADR-0002).  An Arb
caller would need a separate root-finder; that is out of scope and a
deferred bead per Rule 9.

## References

  - Jorba & Zou 2005, Experimental Mathematics 14, §3.3.1 eq. 11 —
    `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/
    JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-645`.
  - `external/TaylorIntegration.jl/src/integrator/stepsize.jl:17-89`
    — canonical Julia implementation, ported verbatim.
  - FW 2011 §3.1 pole-distance heuristic —
    `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md` §3.1.
  - `test/_oracle_stepcontrol.jl` — three-source-consensus pinned
    values (TI.jl, mpmath, Mathematica, 47 decimal digits).
"""
module StepControl

using Polynomials: Polynomial, roots

export step_jorba_zou, step_pade_root

# -----------------------------------------------------------------------------
# Jorba-Zou step (TI.jl-equivalent fixed-order reduction of §3.3.1 eq. 11)
# -----------------------------------------------------------------------------

"""
    step_jorba_zou(coefs::AbstractVector{T}, eps_abs::Real;
                   eps_rel::Real = eps_abs) -> T

Jorba-Zou 2005 §3.3.1 eq. 11 first step-size control, in the
fixed-order form ported from `TaylorIntegration.jl`.

Given Taylor coefficients `coefs = [c₀, c₁, …, c_p]` (length `p + 1`)
and tolerances `eps_abs`, `eps_rel`, returns the step size

    h = min over k ∈ {p-1, p} of (ε / |coefs[k+1]|)^(1/k)

where `ε = eps_abs` if `eps_abs ≥ eps_rel · |c₀|` and `ε = eps_rel ·
|c₀|` otherwise (TI.jl `stepsize.jl:27-31`).  Zero coefficients in
the candidate set are skipped.  When both `c_{p-1}` and `c_p` are
zero, falls back to the TI.jl second-stepsize scan over indices
`1 ≤ j ≤ p-2`, taking the maximum of `(1 / |c[j+1]|)^(1/j)`.

Throws `ArgumentError` if `coefs` has fewer than 2 entries (no notion
of "last two"), or if every nonzero candidate has been exhausted —
both indicate caller-side bugs that we surface loudly per Rule 1
rather than silently returning a meaningless `Inf` or `0`.
"""
function step_jorba_zou(coefs::AbstractVector{T}, eps_abs::Real;
                        eps_rel::Real = eps_abs) where {T}
    length(coefs) ≥ 2 || throw(ArgumentError(
        "step_jorba_zou: coefs must have length ≥ 2 (got $(length(coefs))); " *
        "the formula uses the last two coefficients. " *
        "Suggestion: pass a Taylor jet of order ≥ 1."))

    p = length(coefs) - 1
    R = float(real(T))
    c0_norm = abs(coefs[1])

    # ε resolution per TI.jl stepsize.jl:27-31.  The default eps_rel =
    # eps_abs reproduces the standard "absolute tolerance" mode for our
    # `c₀ ≈ 1` cases (test 4.1.1: `c₀ = 1`, eps_abs = eps_rel = 1e-12,
    # so eps_abs ≥ eps_rel · |c₀| holds with equality and we use eps_abs).
    ε = eps_abs ≥ eps_rel * c0_norm ? R(eps_abs) : R(eps_rel * c0_norm)

    h = typemax(R)
    for k in (p - 1, p)
        aux = abs(coefs[k + 1])
        iszero(aux) && continue
        h = min(h, (ε / aux)^(R(1) / k))
    end

    if !isinf(h)
        return h
    end

    # Fallback: TI.jl `_second_stepsize` (stepsize.jl:77-89).  Scan
    # j = 1 .. p-2 for the largest h such that |c[j+1]| · h^j = 1,
    # i.e. h = (1/|c[j+1]|)^(1/j); take the max.
    h2 = zero(R)
    for j in 1:(p - 2)
        aux = abs(coefs[j + 1])
        iszero(aux) && continue
        h2 = max(h2, (R(1) / aux)^(R(1) / j))
    end
    iszero(h2) && throw(ArgumentError(
        "step_jorba_zou: all Taylor coefficients are zero; cannot " *
        "estimate a step. Suggestion: check the integrator's right-" *
        "hand side and initial condition — a constant zero solution " *
        "needs no stepping."))
    return h2
end

# -----------------------------------------------------------------------------
# Pole-distance step (FW 2011 §3.1)
# -----------------------------------------------------------------------------

"""
    step_pade_root(P, z_current::Number, target::Number) -> Real

FW 2011 §3.1 pole-distance heuristic.  Given a Padé approximant `P`
(of the local Taylor jet) and a desired step `z_current → target`,
returns the largest step length along the direction `target -
z_current` that does not cross any forward-projecting root of `P.b`
(the denominator polynomial).  The result is capped at the full
distance `|target - z_current|`, so a clean local jet (no nearby
roots) always permits the full requested step.

Algorithm:

  1. If the denominator is constant (`length(P.b) ≤ 1`), there are
     no poles to dodge — return `|target - z_current|`.
  2. If `target == z_current`, the step is degenerate — return `0`.
  3. Compute the roots of `Polynomials.Polynomial(P.b)`.
  4. For each root `r`, project onto the unit step direction:
     `t = Re((r - z_current) · conj(unit))`.
  5. Keep only `t > 0` (forward projections); discard `t ≤ 0`.
  6. Return `min(min(forward), |target - z_current|)`, falling back
     to `|target - z_current|` when no forward roots exist.

For the all-real case (real `P.b`, real `z_current`, real `target`)
the result is the real-axis distance to the closest forward pole;
`Polynomials.roots` may return `ComplexF64` even for real input
(test 4.1.4: roots `3 ± i` from `b = [1, -0.6, 0.1]`), and the
projection picks up the real part automatically.
"""
function step_pade_root(P, z_current::Number, target::Number)
    Δ = target - z_current
    dlen = abs(Δ)
    iszero(Δ) && return zero(dlen)

    # Constant denominator: no poles, step is full distance.
    length(P.b) ≤ 1 && return dlen

    unit = Δ / dlen
    rs = roots(Polynomial(P.b))

    h = dlen
    any_forward = false
    for r in rs
        t = real((r - z_current) * conj(unit))
        if t > 0
            any_forward = true
            h = min(h, t)
        end
    end

    return any_forward ? h : dlen
end

end # module StepControl
