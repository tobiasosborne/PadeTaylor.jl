# examples/teukolsky_qnm.jl — Schwarzschild quasinormal-mode frequency via
# Leaver's radial continued fraction.
#
# ## From Teukolsky to a continued-fraction eigenvalue problem
#
# A perturbation of a Schwarzschild black hole of mass M, decomposed into
# spin-weight-`s` spherical harmonics with multipole `ℓ`, obeys the radial
# Teukolsky (equivalently, for s=0 the Regge–Wheeler) equation.  In the
# tortoise coordinate this is a 1-D wave equation with an effective potential
# barrier (`docs/teukolsky_heun_mapping.md:17-47`).  The equation has a
# **regular** singular point at the horizon r = 2M and an **irregular**
# singular point at spatial infinity r = ∞; this triple-singularity structure
# is exactly that of the **confluent Heun equation**
# (`docs/teukolsky_heun_mapping.md:51-176`).
#
# A quasinormal mode (QNM) is a solution that is **purely ingoing at the
# horizon** and **purely outgoing at infinity** simultaneously
# (`docs/teukolsky_heun_mapping.md:271-298`).  Both conditions hold only for a
# discrete set of complex frequencies ω = ω_R − iω_I (Im ω < 0 ⟹ damped, in
# the e^{−iωt} convention).  The outgoing wave at infinity is, however,
# *exponentially subdominant* to the ingoing one there, so naive forward
# integration of the ODE drowns the mode we want in the growing partner
# solution: a direct shooting method is numerically hopeless.  This is why the
# package's HeunC Frobenius recurrence about the horizon — which encodes
# regularity at r = 2M but says nothing about the irregular point r = ∞ —
# cannot by itself locate the frequency.
#
# ## Leaver's series and its three-term recurrence
#
# Leaver (1985) sidesteps the cancellation problem by building the boundary
# behaviour at *both* singular points into the ansatz analytically and
# expanding the smooth remainder in powers of the compactified coordinate
# x = 1 − 2M/r ∈ [0, 1) (`docs/teukolsky_heun_mapping.md:178-191, 226-228`):
#
#     R(r) = e^{iωr} (r − 2M)^{−2iωM} r^{4iωM} · Σ_{n≥0} a_n [1 − 2M/r]^n .
#
# Substituting into the radial equation collapses it to a **three-term
# recurrence** for the coefficients aₙ (`docs/teukolsky_heun_mapping.md:230-244`):
#
#     α_n a_{n+1} + β_n a_n + γ_n a_{n−1} = 0   (n ≥ 1),
#     β_0 a_0 + α_0 a_1 = 0                      (seed).
#
# Generically a three-term recurrence has two linearly independent solutions,
# one growing and one decaying; an arbitrary seed excites the growing one and
# `aₙ` blows up, so the series diverges.  Convergence of R(r) at infinity is
# precisely the demand that we sit on the **minimal (recessive) solution** of
# the recurrence — and that is the outgoing boundary condition in disguise.
#
# ## The continued-fraction (connection) condition
#
# Pincherle's theorem turns "the minimal solution exists with this seed" into
# the vanishing of a continued fraction (`docs/teukolsky_heun_mapping.md:246-260,
# 314-335`).  For the fundamental overtone (0th inversion):
#
#     0 = β_0 − α_0 γ_1 / (β_1 − α_1 γ_2 / (β_2 − ⋯)) .
#
# This is the numerically-stable algebraic form of the confluent-Heun
# minimal-solution / connection condition.  The QNM frequencies are its
# complex roots.  We truncate the tail at depth N ~ 100 and root-find with a
# complex Newton step seeded near the known mode.
#
# ## Conventions and units
#
# We work in **M = 1 geometric units** (horizon r₊ = 2M = 2) and the e^{−iωt}
# time convention, so the recovered ω is the dimensionless Mω that the
# literature tabulates directly.  The explicit coefficients below are the
# cross-verified ground truth recorded in `docs/teukolsky_heun_mapping.md`
# §3.2.1 (`docs/teukolsky_heun_mapping.md:209-265`); the auxiliary constant
# β_aux = s² − 1 carries the spin dependence (s=0 ⟹ β_aux = −1, the scalar
# Regge–Wheeler case; s=2 ⟹ β_aux = +3, the gravitational case).

# `PadeTaylor` is imported for the examples-env context; the Stage-6 companion
# figure (figures/teukolsky_qnm_demo.jl) uses `heun_c` to render the QNM
# eigenfunction at the frequency found here. The frequency itself is pure
# complex arithmetic, so this file adds no other dependency.
using PadeTaylor

"""
    leaver_recurrence(n, ω; s, ℓ) -> (α_n, β_n, γ_n)

The three coefficients of Leaver's Schwarzschild radial recurrence in **M = 1
units**, verbatim from `docs/teukolsky_heun_mapping.md:234-244`:

    α_n = (n+1)(n + 1 − 4iω)
    β_n = −2n² + (−2 + 16iω)·n + 32ω² + 8iω − ℓ(ℓ+1) + β_aux
    γ_n = (n − 4iω)² − (β_aux + 1)

with auxiliary β_aux = s² − 1.  `ω` is the dimensionless Mω (complex).
"""
function leaver_recurrence(n::Integer, ω::Complex; s::Integer, ℓ::Integer)
    β_aux = s^2 - 1
    nn = float(n)
    αn = (nn + 1) * (nn + 1 - 4im * ω)
    βn = -2nn^2 + (-2 + 16im * ω) * nn + 32 * ω^2 + 8im * ω - ℓ * (ℓ + 1) + β_aux
    γn = (nn - 4im * ω)^2 - (β_aux + 1)
    return αn, βn, γn
end

"""
    qnm_continued_fraction(ω; s, ℓ, N) -> f(ω)

Value of Leaver's 0th-inversion continued fraction

    f(ω) = β_0 − α_0 γ_1 / (β_1 − α_1 γ_2 / (β_2 − ⋯)) ,

whose complex zeros are the QNM frequencies
(`docs/teukolsky_heun_mapping.md:246-260`).  The tail is summed by **downward
recursion** from depth `N`: set the deepest tail T_{N+1} = 0, then
T_k = β_k − α_k γ_{k+1} / T_{k+1} for k = N, …, 1, and finally
f = β_0 − α_0 γ_1 / T_1.  Downward evaluation is the standard stable scheme
for a minimal-solution continued fraction (it converges *to* the recessive
solution rather than away from it).
"""
function qnm_continued_fraction(ω::Complex; s::Integer, ℓ::Integer, N::Integer = 100)
    # Bottom of the ladder: tail beyond depth N is truncated to zero.
    _, βN, _ = leaver_recurrence(N, ω; s = s, ℓ = ℓ)
    tail = βN
    # Recurse upward (decreasing k): T_k = β_k − α_k γ_{k+1}/T_{k+1}.
    for k in (N - 1):-1:1
        αk, βk, _ = leaver_recurrence(k, ω; s = s, ℓ = ℓ)
        _, _, γk1 = leaver_recurrence(k + 1, ω; s = s, ℓ = ℓ)
        tail = βk - αk * γk1 / tail
    end
    α0, β0, _ = leaver_recurrence(0, ω; s = s, ℓ = ℓ)
    _, _, γ1 = leaver_recurrence(1, ω; s = s, ℓ = ℓ)
    return β0 - α0 * γ1 / tail
end

"""
    find_qnm(ω₀; s, ℓ, N, tol, maxiter) -> ω*

Locate a QNM frequency as a complex root of `qnm_continued_fraction`, using
Newton's method with a finite-difference (complex-step-sized) derivative
seeded at `ω₀`.  Fails loud (Rule 1) — throws on non-convergence or a NaN/Inf
iterate rather than returning garbage.
"""
function find_qnm(ω₀::Complex; s::Integer, ℓ::Integer, N::Integer = 100,
                  tol::Real = 1e-12, maxiter::Integer = 100)
    ω = ω₀
    h = 1e-8
    for it in 1:maxiter
        f = qnm_continued_fraction(ω; s = s, ℓ = ℓ, N = N)
        if !isfinite(f)
            throw(ErrorException(
                "find_qnm: continued fraction returned non-finite value at " *
                "iterate ω=$ω (s=$s, ℓ=$ℓ, N=$N). " *
                "Suggestion: reduce depth N (a pole of the truncated tail may " *
                "sit on the iterate) or move the seed ω₀ off the singular point."))
        end
        abs(f) < tol && return ω
        fp = (qnm_continued_fraction(ω + h; s = s, ℓ = ℓ, N = N) - f) / h
        if !isfinite(fp) || fp == 0
            throw(ErrorException(
                "find_qnm: degenerate Newton derivative (f'=$fp) at ω=$ω. " *
                "Suggestion: perturb the seed ω₀ or the step h=$h."))
        end
        ω -= f / fp
    end
    throw(ErrorException(
        "find_qnm: Newton iteration failed to reach |f| < $tol within " *
        "$maxiter iterations (last ω=$ω, s=$s, ℓ=$ℓ, N=$N). " *
        "Suggestion: improve the seed ω₀, raise maxiter, or increase depth N."))
end

# Script entry point: recover both validation targets from
# `docs/teukolsky_heun_mapping.md:252-255` and report recovered Mω + |error|.
if abspath(PROGRAM_FILE) == @__FILE__
    targets = (
        (s = 0, ℓ = 2, seed = 0.4836 - 0.0968im, ref = 0.483642 - 0.096766im,
         label = "scalar  s=0 ℓ=2 n=0"),
        (s = 2, ℓ = 2, seed = 0.3737 - 0.0890im, ref = 0.373672 - 0.088962im,
         label = "grav.   s=2 ℓ=2 n=0"),
    )
    println("Schwarzschild QNM via Leaver's radial continued fraction (M=1 units)")
    println("="^70)
    for t in targets
        ω = find_qnm(complex(t.seed); s = t.s, ℓ = t.ℓ, N = 100)
        err = abs(ω - t.ref)
        mω = complex(round(real(ω); digits = 6), round(imag(ω); digits = 6))
        println("  ", rpad(t.label, 20), "  Mω = ", mω,
                "   |err| = ", round(err; sigdigits = 3))
    end
end
