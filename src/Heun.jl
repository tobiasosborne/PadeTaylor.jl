"""
    Heun

PadeTaylor.jl module for the **general Heun function** `HeunG(a,q,α,β,γ,δ;z)`
and the **confluent Heun function** `HeunC(q,α,γ,δ,ε;z)`.  Two textbook
special functions (Heun 1889, DLMF Ch. 31) sitting one layer up from
hypergeometric `₂F₁`: four regular singular points instead of three
(HeunG), or three regular + one irregular at infinity (HeunC).

# Background

Heun functions arise across mathematical physics:

  - Black-hole quasinormal modes — Teukolsky / Regge-Wheeler equations
    reduce to confluent Heun (Fiziev 2009; see
    `docs/teukolsky_heun_mapping.md`).
  - Hydrogen molecule ion `H₂⁺` and two-centre Coulomb problems
    (separation in spheroidal coordinates).
  - Bethe-ansatz eigenstates in integrable spin chains.
  - Lamé and Mathieu functions as special cases (DLMF §31.7).

Despite this, no good open-source numerical implementation existed
before 2015 — Motygin's GPL-3.0 MATLAB code (Motygin 2015, 2018, 2020)
remains the lone open prior art.  Mathematica's built-in `HeunG` /
`HeunC` are documented as fragile past the second singular point
(Wolfram blog 2020-05-06: *"We do not know their explicit forms and
are forced to work with their generating equations."*).

This module is the **first open-source Julia implementation** of Heun
functions (to our knowledge), and the second open-source implementation
of any language.

# Algorithm

The strategy is the **systematic application of PadeTaylor.jl's
existing path-network machinery to a different ODE**.  Three ingredients:

1. **Frobenius series at `z = 0`** for the initial condition.  DLMF
   31.3.2–31.3.4 gives the explicit three-term recurrence for HeunG;
   Motygin 2018 eqs 3–4 gives it for HeunC (DLMF Ch. 31 doesn't publish
   the confluent recurrence — refers reader to Ronveaux 1995).  Both
   recurrences validated at machine precision against Mathematica
   wolframscript oracle on the easy-parameter regime (see
   `external/probes/heun-oracle/oracles.txt`).

2. **Path-network continuation** for `|z|` beyond the Frobenius disc
   of convergence (`min(1, |a|)` for HeunG; `1` for HeunC).  Sample
   the Frobenius series at small `z₀ = 0.05` to get `(u, u')`, then
   walk to the user's target `z` via `path_network_solve` with
   `branch_points = (0, 1, a)` (HeunG) or `(0, 1)` (HeunC).  The
   walker's wedge selector + BranchTracker (ADR-0013) detour around
   the regular singular points automatically.

3. **Single-value evaluation API** — `heun_g(a,q,α,β,γ,δ; z)` returns
   the value at a single point, hiding the path-network plumbing.
   For multi-point queries, use `path_network_solve` directly on the
   Heun RHS factories `_heun_g_rhs` / `_heun_c_rhs`.

# Parameter conventions

DLMF §31.2 / §31.12 conventions throughout (same as Mathematica;
matches Maple HeunG; diverges from Maple HeunC).

  - **HeunG**: `heun_g(a, q, α, β, γ, δ; z)`, with Fuchsian constraint
    `ε = α + β + 1 - γ - δ` derived internally; normalised `HeunG(0) = 1`.
  - **HeunC**: `heun_c(q, α, γ, δ, ε; z)`, normalised `HeunC(0) = 1`.

# References

  - DLMF Ch. 31 — `references/heun/dlmf-31/{31.2,31.3,31.4,31.12}.html`.
  - `docs/heun_dlmf_summary.md` — the recurrence transcription.
  - `docs/heun_research.md` — broader Heun-numerics survey.
  - `docs/adr/0018-heun-module.md` — design decisions.
  - Motygin 2015 (arxiv:1506.03848), 2018 (arxiv:1804.01007), 2020
    (arxiv:2010.09053) — the reference open-source impl in MATLAB.
"""
module Heun

using ..Problems:    PadeTaylorProblem
using ..PathNetwork: path_network_solve

export heun_g, heun_c

# -----------------------------------------------------------------------------
# Frobenius series at z = 0
# -----------------------------------------------------------------------------

"""
    _heun_g_frobenius_at_0(a, q, α, β, γ, δ, z, N) -> (value, derivative, coeffs)

Evaluate the HeunG Frobenius series at `z` using DLMF 31.3.2–31.3.4's
three-term recurrence, truncated at `N+1` terms.

```
c_0 = 1
c_1 = q / (a γ)
For j ≥ 1:
  c_{j+1} = ((Q_j + q) c_j - P_j c_{j-1}) / R_j
  P_j = (j-1+α)(j-1+β)
  Q_j = j · ((j-1+γ)(1+a) + a δ + ε)        ε = α+β+1-γ-δ
  R_j = a (j+1)(j+γ)
```

Returns `(u(z), u'(z), c)` where `c` is the full coefficient vector.
The derivative is computed by term-by-term differentiation.

Throws `ArgumentError` if `γ` is a non-positive integer (DLMF 31.3.1
existence condition fails — logarithmic second solution appears;
Motygin 2020's regularised `Hl☼` not implemented at v1).
"""
function _heun_g_frobenius_at_0(a::T, q::T, α::T, β::T, γ::T, δ::T,
                                z, N::Integer) where {T}
    real_γ = real(γ)
    (real_γ <= 0 && isinteger(real_γ) && abs(imag(γ)) < eps(real(T))) && throw(
        ArgumentError(
            "heun_g: γ = $γ is a non-positive integer; HeunG existence " *
            "condition (DLMF 31.3.1) fails — second Frobenius solution " *
            "acquires a logarithm. Suggestion: perturb γ slightly, or " *
            "use a problem reformulation that avoids the integer-γ pole."))
    a == zero(T) && throw(ArgumentError(
        "heun_g: parameter a = 0 — singular points 0 and a collide. " *
        "Suggestion: use the confluent-Heun limit `heun_c` instead."))
    abs(a - one(T)) < eps(real(T)) && throw(ArgumentError(
        "heun_g: parameter a = 1 — singular points 1 and a collide. " *
        "Suggestion: parameter perturbation or confluent-limit reformulation."))

    ε = α + β + one(T) - γ - δ                       # Fuchsian constraint
    CT = promote_type(typeof(z), T)
    c = Vector{CT}(undef, N + 1)
    c[1] = one(CT)
    c[2] = q / (a * γ)
    for j in 1:(N-1)
        P_j = (j - 1 + α) * (j - 1 + β)
        Q_j = j * ((j - 1 + γ) * (1 + a) + a * δ + ε)
        R_j = a * (j + 1) * (j + γ)
        c[j+2] = ((Q_j + q) * c[j+1] - P_j * c[j]) / R_j
    end

    # Horner evaluation of u(z) and u'(z) = Σ k c_k z^{k-1}
    u  = c[N+1]
    up = N * c[N+1]
    for k in (N-1):-1:0
        u  = u * z + c[k+1]
        k >= 1 && (up = up * z + k * c[k+1])
    end
    return u, up, c
end

"""
    _heun_c_frobenius_at_0(q, α, γ, δ, ε, z, N) -> (value, derivative, coeffs)

Evaluate the HeunC Frobenius series at `z` using Motygin 2018 eqs 3–4
(`docs/heun_research.md:563-579`), truncated at `N+1` terms.

```
b_0 = 1, b_{-1} = 0
b_1 = -q/γ
For n ≥ 2:
  b_n = (Q_n b_{n-1} + R_n b_{n-2}) / P_n
  P_n = n (n + γ - 1)
  Q_n = -q + (n-1)(γ + δ - ε + n - 2)
  R_n = α + (n-2) ε
```

Throws `ArgumentError` if `γ = 0` (P_1 vanishes) or if `γ` is a
non-positive integer (similar logarithmic-degeneration condition).
"""
function _heun_c_frobenius_at_0(q::T, α::T, γ::T, δ::T, ε::T,
                                z, N::Integer) where {T}
    real_γ = real(γ)
    (real_γ <= 0 && isinteger(real_γ) && abs(imag(γ)) < eps(real(T))) && throw(
        ArgumentError(
            "heun_c: γ = $γ is a non-positive integer; first recurrence " *
            "denominator vanishes (Motygin 2018 eq 3 / DLMF 31.12). " *
            "Suggestion: perturb γ slightly."))

    CT = promote_type(typeof(z), T)
    b = Vector{CT}(undef, N + 1)
    b[1] = one(CT)
    b[2] = -q / γ
    for n in 2:N
        P_n = n * (n + γ - 1)
        Q_n = -q + (n - 1) * (γ + δ - ε + n - 2)
        R_n = α + (n - 2) * ε
        b[n+1] = (Q_n * b[n] + R_n * b[n-1]) / P_n
    end

    u  = b[N+1]
    up = N * b[N+1]
    for k in (N-1):-1:0
        u  = u * z + b[k+1]
        k >= 1 && (up = up * z + k * b[k+1])
    end
    return u, up, b
end

# -----------------------------------------------------------------------------
# RHS factories for path-network continuation
# -----------------------------------------------------------------------------

"""
    _heun_g_rhs(a, q, α, β, γ, δ) -> f

Return a thread-safe closure `f(z, u, up) -> u''` implementing DLMF
31.2.1 rearranged as `u'' = -(p(z) up + r(z) u)`:

```
u'' = -((γ/z + δ/(z-1) + ε/(z-a)) up + (α β z - q)/(z(z-1)(z-a)) u)
```

with `ε = α + β + 1 - γ - δ`.  The returned closure captures no
mutable state — safe to pass to `path_network_solve`'s parallel
wedge evaluator.
"""
function _heun_g_rhs(a, q, α, β, γ, δ)
    ε = α + β + 1 - γ - δ
    f(z, u, up) = -(
        (γ / z + δ / (z - 1) + ε / (z - a)) * up +
        (α * β * z - q) / (z * (z - 1) * (z - a)) * u
    )
    return f
end

"""
    _heun_c_rhs(q, α, γ, δ, ε) -> f

Return a thread-safe closure implementing DLMF 31.12.1:

```
u'' = -((γ/z + δ/(z-1) + ε) up + (α z - q)/(z(z-1)) u)
```
"""
function _heun_c_rhs(q, α, γ, δ, ε)
    f(z, u, up) = -(
        (γ / z + δ / (z - 1) + ε) * up +
        (α * z - q) / (z * (z - 1)) * u
    )
    return f
end

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

const _DEFAULT_EPS_START = 0.05
const _DEFAULT_N_TAYLOR  = 60

"""
    heun_g(a, q, α, β, γ, δ; z, h=0.1, order=30,
           epsilon_start=0.05, n_taylor=60, kwargs...) -> Complex

Evaluate the general Heun function `HeunG(a, q, α, β, γ, δ; z)` at
the complex point `z`.  Returns a single `Complex` value.  Normalised
so that `heun_g(...; z=0) = 1`.

For `|z| ≤ epsilon_start`, returns the Frobenius series sum
(DLMF 31.3.2–31.3.4) truncated at `n_taylor + 1` terms — no
path-network walk.  For larger `|z|`, samples the Frobenius series
at `z₀ = epsilon_start` for an IC tuple, then walks
`PadeTaylorProblem` via `path_network_solve` from `z₀` to `z` with
`branch_points = (0, 1, a)`.

# Parameters (DLMF §31.2 conventions)

  - `a` — singularity location (`a ∉ {0, 1}`).
  - `q` — accessory parameter.
  - `α, β` — exponent parameters at `∞`.
  - `γ, δ` — exponents at `z = 0, 1`; the exponent at `z = a` is
    `ε = α + β + 1 - γ - δ` (Fuchsian constraint, derived internally).

# Kwargs

  - `z::Complex` — evaluation point.  Required.
  - `h` — path-network step length; default `0.1` (smaller than the
    Painlevé default `0.5` because Heun's rational coefficients vary
    more rapidly near the singular points).
  - `order` — Taylor truncation degree per path-network step; default
    `30` (FW 2011 canonical).
  - `epsilon_start` — radius below which we use the Frobenius series
    directly; default `0.05`.
  - `n_taylor` — Frobenius series truncation order; default `60`.
  - Any other kwarg is forwarded to `path_network_solve` (e.g.
    `cross_branch`, `enforce_real_axis_symmetry`).

Throws `ArgumentError` on parameter values that violate existence
conditions (γ a non-positive integer, a ∈ {0, 1}).
"""
function heun_g(a, q, α, β, γ, δ;
                z,
                h::Real = 0.05,
                order::Integer = 30,
                epsilon_start::Real = _DEFAULT_EPS_START,
                n_taylor::Integer = _DEFAULT_N_TAYLOR,
                kwargs...)
    # Promote parameters to a complex floating-point type so integer
    # inputs (`heun_g(2, 1, 1, 2, 3, 4; z=...)`) work and the inner
    # recurrence has access to `eps(real(T))` for existence checks.
    CT = promote_type(typeof(complex(float(a))), typeof(complex(float(q))),
                      typeof(complex(float(α))), typeof(complex(float(β))),
                      typeof(complex(float(γ))), typeof(complex(float(δ))),
                      typeof(z), Complex{Float64})
    T = CT
    z_c  = CT(z)
    a_c  = T(a)
    q_c  = T(q)
    α_c  = T(α); β_c = T(β); γ_c = T(γ); δ_c = T(δ)

    # Direct Frobenius for small |z|
    if abs(z_c) ≤ epsilon_start
        u, _, _ = _heun_g_frobenius_at_0(a_c, q_c, α_c, β_c, γ_c, δ_c,
                                          z_c, n_taylor)
        return u
    end

    # Detect real-valued parameters → enables Schwarz-symmetric upper-half
    # walk (enforce_real_axis_symmetry=true) that matches Wolfram's
    # principal-sheet convention.
    real_params = abs(imag(a_c)) < eps() && abs(imag(q_c)) < eps() &&
                  abs(imag(α_c)) < eps() && abs(imag(β_c)) < eps() &&
                  abs(imag(γ_c)) < eps() && abs(imag(δ_c)) < eps()

    z0 = CT(epsilon_start)
    u0, up0, _ = _heun_g_frobenius_at_0(a_c, q_c, α_c, β_c, γ_c, δ_c,
                                         z0, n_taylor)
    f = _heun_g_rhs(a_c, q_c, α_c, β_c, γ_c, δ_c)
    prob = PadeTaylorProblem(f, (u0, up0), (z0, z_c); order = order)
    # Cuts oriented downward keep them out of the upper half plane so
    # the walker can detour upward around singular points and return
    # Wolfram's upper-limit principal-sheet value
    # (heun_research.md:533-539).
    if real_params
        # enforce_real_axis_symmetry is mutually exclusive with
        # branch_points (path_network_solve.jl:401 / ADR-0013 open
        # follow-ups).  Schwarz mirror suffices: the walker stays in
        # Im z ≥ 0, and the local Padé / min-|u| wedge selector avoids
        # the singularities at z=1 and z=a automatically (FW 2011 §3.1).
        sol = path_network_solve(prob, [z_c];
                                  h = h, order = order,
                                  enforce_real_axis_symmetry = true,
                                  kwargs...)
    else
        sol = path_network_solve(prob, [z_c];
                                  h = h, order = order,
                                  branch_points = (CT(0), CT(1), CT(a_c)),
                                  branch_cut_angles = (-π/2, -π/2, -π/2),
                                  kwargs...)
    end
    return sol.grid_u[1]
end

"""
    heun_c(q, α, γ, δ, ε; z, h=0.1, order=30, ...) -> Complex

Evaluate the confluent Heun function `HeunC(q, α, γ, δ, ε; z)` at `z`.
Returns a single `Complex` value normalised `heun_c(...; z=0) = 1`.

Same Frobenius-then-path-network strategy as `heun_g`; uses
`branch_points = (0, 1)` (no singularity at finite `z = a` since
that one merged into infinity).  Recurrence per Motygin 2018 eqs 3–4
(DLMF Ch. 31 doesn't publish HeunC's recurrence; refers to Ronveaux
1995).

Throws `ArgumentError` if `γ` is a non-positive integer.
"""
function heun_c(q, α, γ, δ, ε;
                z,
                h::Real = 0.05,
                order::Integer = 30,
                epsilon_start::Real = _DEFAULT_EPS_START,
                n_taylor::Integer = _DEFAULT_N_TAYLOR,
                kwargs...)
    CT = promote_type(typeof(complex(float(q))), typeof(complex(float(α))),
                      typeof(complex(float(γ))), typeof(complex(float(δ))),
                      typeof(complex(float(ε))), typeof(z), Complex{Float64})
    T = CT
    z_c = CT(z)
    q_c = T(q); α_c = T(α); γ_c = T(γ); δ_c = T(δ); ε_c = T(ε)

    if abs(z_c) ≤ epsilon_start
        u, _, _ = _heun_c_frobenius_at_0(q_c, α_c, γ_c, δ_c, ε_c,
                                          z_c, n_taylor)
        return u
    end

    real_params = abs(imag(q_c)) < eps() && abs(imag(α_c)) < eps() &&
                  abs(imag(γ_c)) < eps() && abs(imag(δ_c)) < eps() &&
                  abs(imag(ε_c)) < eps()

    z0 = CT(epsilon_start)
    u0, up0, _ = _heun_c_frobenius_at_0(q_c, α_c, γ_c, δ_c, ε_c,
                                         z0, n_taylor)
    f = _heun_c_rhs(q_c, α_c, γ_c, δ_c, ε_c)
    prob = PadeTaylorProblem(f, (u0, up0), (z0, z_c); order = order)
    if real_params
        sol = path_network_solve(prob, [z_c];
                                  h = h, order = order,
                                  enforce_real_axis_symmetry = true,
                                  kwargs...)
    else
        sol = path_network_solve(prob, [z_c];
                                  h = h, order = order,
                                  branch_points = (CT(0), CT(1)),
                                  branch_cut_angles = (-π/2, -π/2),
                                  kwargs...)
    end
    return sol.grid_u[1]
end

end # module Heun
