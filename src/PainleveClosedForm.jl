# src/PainleveClosedForm.jl — included into `module Painleve` (bead
# `padetaylor-icf`, ADR-0010).
#
# ## Why this file exists
#
# `PainleveNamed.jl` ships the *point-IC* named transcendents
# (`tritronquee(:I)`, `hastings_mcleod()`) — solutions singled out by a
# 16-digit `(u, u')` literature value at a specific point.  This file
# is the second kind in the ADR-0008 taxonomy: **parametrised families
# with exact closed-form formulas**.  They differ in shape from the
# point-IC kind:
#
#   - The IC is derived from the formula, not stamped as a literal.
#   - The formula is itself the oracle for self-validating tests
#     (solver output `≈ closed-form(z)` at downstream points).
#   - A discriminating parameter (`n`, `θ`, `kind`) selects among
#     several solutions of the *same* equation.
#
# ## What v1 ships
#
# Three families.  See ADR-0010 for the scoping logic.
#
#   - **PII rational `u_n`** for `α = n`, `n ∈ {1, 2, 3}`.  Formulas:
#     `FW2014_*.md:119-120` eq. 6.  `u_3` is the well-known
#     `(0, 0)`-IC rational solution; `u_1` and `u_2` carry real-axis
#     poles whose location depends on `n`, so the constructor
#     fails-loud if the caller picks `zspan[1]` AT a pole.
#
#   - **PII Airy `u_{n+1/2}`** for `α = n + 1/2`, `n ∈ {0, 1}`.
#     Formulas: `FW2014_*.md:144-159` eqs. 8-9.  Base function
#     `φ(z; θ) = cos(θ/2)·Ai(-z/2^⅓) + sin(θ/2)·Bi(-z/2^⅓)` solves
#     `φ'' = -(z/2)φ` (an Airy-type equation); `Φ = φ'/φ` obeys
#     `Φ' = -z/2 - Φ²` (FW2014 md:154), and `u_{1/2} = -Φ`,
#     `u_{3/2} = (2Φ³ + zΦ - 1)/(2Φ² + z)`.  v1 caps at `u_{3/2}`;
#     `u_{5/2}` (FW2014 eq. 10) is in-tree but its derivative
#     formula is hand-derivation-error-prone — deferred.
#
#   - **PIV entire** — `u(z) = -2z` at `(α, β) = (0, -2)` and
#     `u(z) = -(2/3)z` at `(α, β) = (0, -2/9)`.  RF 2014 md:91 cites
#     the two solutions but does NOT state the `(α, β)` pairs in-line.
#     Derived algebraically by substituting `u = -2z` (resp.
#     `-(2/3)z`) into the PIV equation and matching coefficients —
#     derivation recorded in this file's docstring for
#     `_u_piv_entire`.
#
# ## Tests
#
# `test/painleve_closed_form_test.jl` — CF.1.*, CF.2.*, CF.3.*, CF.4.*.
# Test idiom: assert the IC stored in `pp.problem.y0` matches the
# closed form exactly, then solve and assert the solver's trajectory
# matches the closed form at downstream `z` to ~1e-10.  Mutation-proof
# procedure in the test file footer.
#
# Reference: docs/adr/0010-painleve-closed-form-families.md.

using SpecialFunctions: airyai, airyaiprime, airybi, airybiprime

# -----------------------------------------------------------------------------
# Closed-form formula helpers (oracles).
# Each returns `(u, u')` at z evaluated from the analytic formula.
# -----------------------------------------------------------------------------

"""
    _u_pii_rational(n::Integer, z) -> (u, u')

PII rational solutions `u_n(z)` for `α = n ∈ {1, 2, 3}`.  Formulas from
`references/markdown/FW2014_second_PII_exploration_FoCM14/FW2014_second_PII_exploration_FoCM14.md:119-120`
eq. 6:

    u_1(z) = -1/z
    u_2(z) = (4 - 2z³) / (4z + z⁴)
    u_3(z) = 3z²(160 + 8z³ + z⁶) / (320 - 24z⁶ - z⁹)

Derivatives via the quotient rule, all in closed form.
"""
function _u_pii_rational(n::Integer, z)
    if n == 1
        return (-1 / z, 1 / z^2)
    elseif n == 2
        # f = 4 - 2z³,  f' = -6z²
        # g = 4z + z⁴,  g' = 4 + 4z³
        f, fp = 4 - 2*z^3, -6*z^2
        g, gp = 4*z + z^4,  4 + 4*z^3
        return (f / g, (fp*g - f*gp) / g^2)
    elseif n == 3
        # f = 3z²(160 + 8z³ + z⁶) = 480z² + 24z⁵ + 3z⁸
        # f' = 960z + 120z⁴ + 24z⁷
        # g = 320 - 24z⁶ - z⁹
        # g' = -144z⁵ - 9z⁸
        if iszero(z)
            return (zero(z), zero(z))    # u_3(0) = 0, u_3'(0) = 0
        end
        f  = 480*z^2 + 24*z^5 + 3*z^8
        fp = 960*z + 120*z^4 + 24*z^7
        g  = 320 - 24*z^6 - z^9
        gp = -144*z^5 - 9*z^8
        return (f / g, (fp*g - f*gp) / g^2)
    end
    throw(ArgumentError("_u_pii_rational: n must be 1, 2, or 3 (got $n)."))
end

# Airy φ(z; θ) and φ'(z; θ).  Base function for PII Airy solutions:
#   φ(z) = cos(θ/2) · Ai(-z/2^⅓) + sin(θ/2) · Bi(-z/2^⅓)
# satisfies φ'' = -(z/2)φ (FW2014 §2.3).  Chain rule: d/dz[Ai(-z·c)] =
# -c·Ai'(-z·c) where c = 2^(-1/3).
const _AIRY_C = 2^(-1//3)

function _airy_phi(z, θ)
    ω = -z * _AIRY_C
    return cos(θ / 2) * airyai(ω) + sin(θ / 2) * airybi(ω)
end

function _airy_phi_prime(z, θ)
    ω = -z * _AIRY_C
    return -_AIRY_C * (cos(θ / 2) * airyaiprime(ω) + sin(θ / 2) * airybiprime(ω))
end

"""
    _u_pii_airy(n::Integer, θ, z) -> (u, u')

PII Airy solutions `u_{n+1/2}(z; θ)` for `α = n + 1/2`, `n ∈ {0, 1}`.
Formulas from `FW2014_*.md:144-159` eqs. 8-9:

    u_{1/2}(z) = -Φ(z)
    u_{3/2}(z) = (2Φ³ + zΦ - 1) / (2Φ² + z)

where `Φ(z) = φ'(z)/φ(z)` for the Airy φ above.  `Φ' = -z/2 - Φ²`
(FW2014 md:154); derivatives via the quotient rule with Φ' folded in.

`θ ∈ [0, 2π)` parameterises the one-parameter family per `n`.
`θ = 0` gives the pure-Ai branch (the canonical case in FW2014 Fig. 2).
"""
function _u_pii_airy(n::Integer, θ, z)
    φ  = _airy_phi(z, θ)
    φp = _airy_phi_prime(z, θ)
    Φ  = φp / φ
    Φp = -z / 2 - Φ^2            # FW2014 md:154
    if n == 0
        return (-Φ, -Φp)
    elseif n == 1
        # N = 2Φ³ + zΦ - 1,        D = 2Φ² + z
        # dN/dz = (6Φ² + z)·Φ' + Φ
        # dD/dz = 4Φ·Φ' + 1
        N  = 2*Φ^3 + z*Φ - 1
        D  = 2*Φ^2 + z
        Np = (6*Φ^2 + z) * Φp + Φ
        Dp = 4*Φ * Φp + 1
        return (N / D, (Np*D - N*Dp) / D^2)
    end
    throw(ArgumentError("_u_pii_airy: n must be 0 or 1 (got $n).  v1 caps at u_{3/2}."))
end

"""
    _u_piv_entire(kind::Symbol, z) -> (u, u')

PIV entire solutions, `kind ∈ {:minus_2z, :minus_two_thirds_z}`.
RF 2014 (`references/markdown/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:91`)
cites the two solutions; the `(α, β)` parameter pairs are not stated
in-line, derived algebraically here.

Substituting `u(z) = -2z` (so `u' = -2`, `u'' = 0`) into the PIV
equation `u'' = (u')²/(2u) + (3/2)u³ + 4zu² + 2(z²−α)u + β/u`:

  - `(u')²/(2u) = 4/(-4z) = -1/z`
  - `(3/2)u³ = (3/2)·(-8z³) = -12z³`
  - `4zu² = 4z·4z² = 16z³`
  - `2(z²−α)u = -4z³ + 4αz`
  - `β/u = -β/(2z)`

Sum to zero: `(−1 − β/2)/z + 4αz + (−12 + 16 − 4)z³ = 0`, i.e.
`α = 0` and `β = -2`.  The `−(2/3)z` case is analogous; the algebra
yields `(α, β) = (0, -2/9)`.  Both are entire (linear in `z`, no
singularity).
"""
function _u_piv_entire(kind::Symbol, z)
    if kind === :minus_2z
        return (-2 * z, -2 * one(z))
    elseif kind === :minus_two_thirds_z
        return (-2 * z / 3, -2 * one(z) / 3)
    end
    throw(ArgumentError(
        "_u_piv_entire: kind must be :minus_2z or :minus_two_thirds_z " *
        "(got $(repr(kind)))."))
end

# -----------------------------------------------------------------------------
# Constructors — build the `PainleveProblem` from the closed-form IC.
# -----------------------------------------------------------------------------

"""
    pii_rational(n::Integer; zspan, order) -> PainleveProblem

The Painlevé-II rational solution `u_n(z)` for `α = n ∈ {1, 2, 3}`.
Returns a `PainleveProblem` with `equation = :II`, `params = (; α = n)`,
`name = :pii_rational`, and the IC `(u, u') = (u_n(zspan[1]),
u_n'(zspan[1]))` evaluated from the FW 2014 closed form.

Defaults: `zspan = (1.0, 5.0)` for `n ∈ {1, 2}` (avoids the pole at
`z = 0`); `zspan = (0.0, 5.0)` for `n = 3` (regular at `z = 0`,
matching FW 2014 md:218 `(u_3(0), u_3'(0)) = (0, 0)`).  `order = 30`.

Throws `ArgumentError` if `n ∉ {1, 2, 3}` or if `zspan[1]` coincides
with a real-axis pole of `u_n` (CLAUDE.md Rule 1 — silent `Inf` in
the IVP state is unacceptable).
"""
function pii_rational(n::Integer; zspan = nothing, order::Integer = 30)
    n in (1, 2, 3) || throw(ArgumentError(
        "pii_rational: n must be in {1, 2, 3} (got $n); FW 2014 eq. 6 " *
        "gives explicit formulas for these only.  Higher n via Bäcklund " *
        "recurrence is in-tree (FW2014 eq. 4) but deferred to follow-up."))

    # Default zspan: u_1 and u_2 have a pole at z = 0; u_3 is regular
    # there (and FW2014 md:218 names z = 0 as the canonical u_3 IC point).
    zsp = if zspan === nothing
        n == 3 ? (0.0, 5.0) : (1.0, 5.0)
    else
        zspan
    end

    # Pole guard: refuse `zspan[1]` at a known real-axis pole.
    z0 = zsp[1]
    if n == 1 && z0 == 0
        throw(ArgumentError(
            "pii_rational(1): zspan[1] = 0 is a pole of u_1(z) = -1/z " *
            "(would emit ±Inf into the IVP state).  Suggestion: use any " *
            "zspan[1] ≠ 0, e.g. (1.0, 5.0)."))
    end
    if n == 2 && z0 == 0
        throw(ArgumentError(
            "pii_rational(2): zspan[1] = 0 is a pole of u_2(z) = " *
            "(4-2z³)/(4z+z⁴) (denominator vanishes).  Suggestion: use any " *
            "zspan[1] ∉ {0, -∛4}."))
    end
    if n == 2 && isapprox(z0, -cbrt(4); atol = 1e-12)
        throw(ArgumentError(
            "pii_rational(2): zspan[1] ≈ -∛4 ≈ -1.587 is a pole of " *
            "u_2(z) (denominator 4z + z⁴ = z(4 + z³) vanishes).  " *
            "Suggestion: use any zspan[1] ∉ {0, -∛4}."))
    end

    u0, up0 = _u_pii_rational(n, z0)
    pp = PainleveProblem(:II; α = n, u0 = u0, up0 = up0,
                         zspan = zsp, order = order)
    return _with_name(pp, :pii_rational)
end

"""
    pii_airy(n::Integer; θ = 0.0, zspan, order) -> PainleveProblem

The Painlevé-II Airy solution `u_{n+1/2}(z; θ)` for `α = n + 1/2`,
`n ∈ {0, 1}`.  Returns a `PainleveProblem` with `equation = :II`,
`params = (; α = n + 1//2, θ = θ)`, `name = :pii_airy`, IC derived from
the closed form at `zspan[1]`.

The free parameter `θ ∈ [0, 2π)` selects within the one-parameter
family per `n`: `θ = 0` is the pure-Ai branch (FW2014 Fig. 2), `θ = π`
is the pure-Bi branch.  `θ` is stored in `pp.params` for downstream
introspection.

Defaults: `zspan = (0.0, 2.0)` (Airy `Ai/Bi` are smooth on the real
axis; the first zero of `Ai(-z/2^⅓)` is at `z ≈ 2.946`, so the default
end-point stays well clear for `θ = 0`).  `order = 30`.

Throws `ArgumentError` if `n ∉ {0, 1}` (v1 scope).
"""
function pii_airy(n::Integer; θ::Real = 0.0,
                  zspan = (0.0, 2.0), order::Integer = 30)
    n in (0, 1) || throw(ArgumentError(
        "pii_airy: n must be in {0, 1} (got $n); v1 ships u_{1/2} and " *
        "u_{3/2}.  u_{5/2} (FW2014 eq. 10) is in-tree but its derivative " *
        "formula is hand-derivation-error-prone — deferred to follow-up."))

    α = n + 1 // 2
    θ_stored = float(θ)        # normalise Irrational/Rational input to Float
    z0 = zspan[1]
    u0, up0 = _u_pii_airy(n, θ_stored, z0)
    pp = PainleveProblem(:II; α = α, u0 = u0, up0 = up0,
                         zspan = zspan, order = order)
    # Extend params with θ (`PainleveProblem(:II; α)` sets params = (; α);
    # we need (; α, θ) for full provenance).
    return PainleveProblem(pp.equation, (; α = α, θ = θ_stored), pp.problem,
                           pp.frame, pp.to_frame, pp.from_frame, :pii_airy)
end

"""
    piv_entire(kind::Symbol; zspan, order) -> PainleveProblem

The Painlevé-IV entire solutions:

  - `kind = :minus_2z`              — `u(z) = -2z`,      `(α, β) = (0, -2)`
  - `kind = :minus_two_thirds_z`    — `u(z) = -(2/3)z`,  `(α, β) = (0, -2/9)`

(RF 2014 md:91; `(α, β)` derived algebraically — see `_u_piv_entire`'s
docstring.)

Defaults: `zspan = (1.0, 5.0)` (z = 0 makes `u = 0`, where the
PIV RHS's `β/u` term is singular — `PainleveProblem(:IV)` throws there).
`order = 30`.

Throws `ArgumentError` on an unknown `kind` or `zspan[1] = 0`.
"""
function piv_entire(kind::Symbol; zspan = (1.0, 5.0), order::Integer = 30)
    α, β = if kind === :minus_2z
        (0, -2)
    elseif kind === :minus_two_thirds_z
        (0, -2 // 9)
    else
        throw(ArgumentError(
            "piv_entire: kind must be :minus_2z or :minus_two_thirds_z " *
            "(got $(repr(kind)))."))
    end

    z0 = zspan[1]
    z0 == 0 && throw(ArgumentError(
        "piv_entire: zspan[1] = 0 makes u(0) = 0, where the PIV RHS's " *
        "β/u term is singular.  Suggestion: use any zspan[1] ≠ 0."))

    u0, up0 = _u_piv_entire(kind, z0)
    pp = PainleveProblem(:IV; α = α, β = β, u0 = u0, up0 = up0,
                         zspan = zspan, order = order)
    return _with_name(pp, :piv_entire)
end
