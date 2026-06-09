# test/certified_oracle_test.jl -- bead `padetaylor-krgy.5` (epic padetaylor-krgy,
# Tier-2 §2.2 of docs/test_corpus/03_hardening_methodology.md:93-108).
#
# CERTIFIED BALL ORACLES -- pilot on ONE oracle family.
#
# ----------------------------------------------------------------------------
# What "certified" means here (read this twice -- it is the whole point)
# ----------------------------------------------------------------------------
# The corpus already pins package output against high-precision *floats*
# (mpmath / Mathematica / BigFloat point values).  Those goldens are
# trustworthy but not *certified*: a golden could itself be subtly wrong and
# no test would notice.  Arb ball (midpoint-radius) arithmetic fixes this --
# but NOT by the naive route one first reaches for.
#
# THE NAIVE ROUTE FAILS.  `Arblib.contains(tight_ball, package_value)` does
# NOT work: the package computes to ~1e-15 (Float64); a 256-bit Arb ball has
# radius ~1e-74, FAR tighter than the package error.  Containing the package
# value in the tight ball therefore FAILS -- the package value sits outside.
#
# THE CORRECT TWO-STEP PATTERN.  The certification value is different:
#
#   (1) CERTIFY THE ORACLE.  Compute the closed-form oracle (here the
#       Weierstrass ℘) in Arb BALL arithmetic and assert its radius is
#       provably tiny: `radius(ball) < ε_cert` (ε_cert = 1e-40).  THIS is the
#       certification -- it proves the oracle midpoint equals the true value
#       to ε_cert, upgrading "a high-precision float we trust" to "a
#       provably-correct golden".  This removes the "is our golden itself
#       right?" risk that plain high-precision floats silently carry.
#
#   (2) PIN THE PACKAGE.  THEN pin the package's computed value to the ball's
#       midpoint at the PACKAGE's own tolerance:
#       `abs(package - midpoint(ball)) < tol_pkg`,  with tol_pkg >> ε_cert.
#
# ε_cert (oracle exactness, ~1e-74 achieved) and tol_pkg (package accuracy,
# ~1e-15 / ~1e-7 achieved) live on opposite scales by ~30-60 orders of
# magnitude; the gap between them is exactly the headroom the certification
# buys.  See ADR-0029.
#
# ----------------------------------------------------------------------------
# SCOPE (reconciled, bead NOTES + ADR-0002).  Ball arithmetic certifies the
# closed-form ORACLE, NOT the package's SVD-based Padé output: Arblib ships no
# SVD (ADR-0002), so balls do NOT flow through the BigFloat-SVD path.  We
# certify the oracle side; the package value is pinned to the certified
# midpoint as an ordinary high-precision float.  Why BALL not INTERVAL:
# interval arithmetic suffers the dependency problem (overestimation); Arb's
# midpoint-radius form is tighter and auto-propagates error.
#
# ----------------------------------------------------------------------------
# Pilot family: equianharmonic Weierstrass ℘  (PREFERRED per bead).
#
# The headline test problem is `u'' = 6u²` (problems_test.jl 6.1.*), whose
# closed form is `u(z) = ℘(z + c₁; g₂=0, g₃=c₂)` with the FW 2011 IC
# `c₁ = -1, c₂ = 2` -- a g₂=0 (equianharmonic) lattice with a real double
# pole at z = 1.  Arblib EXPOSES ℘ directly via `Arblib.elliptic_p!`
# (acb_elliptic_p), so we build a certified ℘ ball and pin the package's
# `solve_pade` ℘ values to it.  (Fallback `robust_pade`-of-exp was NOT needed.)
#
# Arblib's ℘ takes the lattice as `(z, τ)` with periods (1, τ), giving FIXED
# invariants `(G₂, G₃)` read off via `Arblib.elliptic_invariants!`.  For
# g₂ = 0 the only fundamental-domain point is the equianharmonic
# τ = exp(iπ/3), which yields G₂ = 0, G₃ ≈ 820.82.  We rescale to the
# package's `g₃ = 2` by the ℘ homogeneity relation
#   ℘(z; g₂, g₃) = λ⁻² ℘(z/λ; G₂, G₃),  with g₂ = λ⁻⁴G₂, g₃ = λ⁻⁶G₃,
# so λ = (G₃ / g₃)^(1/6) (real positive; its sixth root keeps λ real and the
# lattice real-period-bearing -- λ ≈ 2.726 is exactly the real pole spacing).
# All of this is done in Acb ball arithmetic, so the resulting ℘ ball carries
# a rigorous radius.
#
# References:
#   - docs/adr/0029-certified-ball-oracles.md   (the decision + tradeoff)
#   - docs/adr/0002-bigfloat-svd-via-genericlinalg.md  (Arblib-has-no-SVD)
#   - docs/adr/0003-extensions-pattern.md       (Arblib is a test extra)
#   - docs/test_corpus/03_hardening_methodology.md:93-108  (Tier-2 §2.2)
#   - test/problems_test.jl 6.1.*               (the package values pinned)
#   - test/_oracle_problems.jl                  (the float goldens certified)
#   - Arblib acb_elliptic.jl: elliptic_p!, elliptic_invariants!  (Arb's ℘)
#   - DLMF 23.x (Weierstrass ℘ + homogeneity 23.10.1).
#
# Mutation-proof procedure recorded in the file footer (both directions).

using Test
using PadeTaylor
using Arblib
using Arblib: Acb, Arb

# ---------------------------------------------------------------------------
# Certified equianharmonic ℘ ball: u(z) = ℘(z; g₂ = 0, g₃ = 2), as a real Arb.
# ---------------------------------------------------------------------------
"""
    _wp_equianharmonic_ball(zval; prec)

Certified `℘(zval; g₂ = 0, g₃ = 2)` as a real `Arb` ball.  Built from Arblib's
lattice-(1, τ) ℘ at the equianharmonic τ = exp(iπ/3) (where G₂ = 0), rescaled
to g₃ = 2 via the homogeneity relation.  The returned ball's radius is the
rigorous Arb enclosure of the true ℘ value.
"""
function _wp_equianharmonic_ball(zval; prec::Int = 256)
    θ   = big(π) / 3
    τ   = Acb(cos(θ), sin(θ); prec = prec)        # equianharmonic period ratio
    G2  = Acb(prec = prec); G3 = Acb(prec = prec)
    Arblib.elliptic_invariants!(G2, G3, τ, prec = prec)   # G2 = 0, G3 ≈ 820.82
    λ   = (G3 / Acb(2; prec = prec))^(Acb(1) / 6)          # g₃ = λ⁻⁶ G₃ = 2
    w   = Acb(zval; prec = prec)
    res = Acb(prec = prec)
    Arblib.elliptic_p!(res, w / λ, τ, prec = prec)         # ℘ on lattice (1, τ)
    out = λ^(-2) * res                                     # homogeneity rescale
    r   = Arb(prec = prec)
    Arblib.set!(r, Arblib.realref(out))                    # real part (℘ real here)
    return r
end

@testset "Certified ball oracles (bead padetaylor-krgy.5): equianharmonic ℘" begin

    PREC    = 256
    ε_cert  = big"1e-40"   # oracle-exactness bar; achieved radii ~1e-74

    # The package's headline pole-bridge problem, IDENTICAL to problems_test
    # 6.1.1: one Padé built at z = 0 with h_max = 1.5 spans the lattice pole
    # at z = 1.  ICs are the FW 2011 16-digit goldens (md:292-295).
    fW(z, u, up) = 6 * u^2
    u0_FW  = 1.071822516416917
    up0_FW = 1.710337353176786
    prob = PadeTaylorProblem(fW, (u0_FW, up0_FW), (0.0, 1.5); order = 30)
    sol  = solve_pade(prob; h_max = 1.5)

    # -----------------------------------------------------------------------
    # CO.1 — the certification step in isolation: the published FW IC golden
    # `u_0_FW` is itself pinned to a CERTIFIED ℘ ball.  This is the purest
    # demonstration: it proves the project's trusted 16-digit IC float is
    # correct (to 4.2e-16) against a provably-tiny (~1e-74) oracle ball.
    # u(0) = ℘(0 - 1; 0, 2).
    # -----------------------------------------------------------------------
    @testset "CO.1: certified ℘ ball validates the published FW IC golden" begin
        ball = _wp_equianharmonic_ball(big(-1.0); prec = PREC)
        rad  = Arblib.radius(Arb, ball)
        mid  = BigFloat(Arblib.midpoint(ball))
        # (1) CERTIFY: the oracle ball is provably tiny.
        @test BigFloat(rad) < ε_cert                       # achieved ~1.2e-74
        # (2) PIN: the published IC golden lands inside at golden tolerance.
        tol_golden = big"1e-15"
        @test abs(big(u0_FW) - mid) < tol_golden           # achieved ~4.2e-16
    end

    # -----------------------------------------------------------------------
    # CO.2 — pin the package solve_pade ℘ value BEFORE the pole.  At z = 0.5
    # the package (Float64, order-30 (15,15) Padé) is at its accuracy floor.
    # u(0.5) = ℘(0.5 - 1; 0, 2).
    # -----------------------------------------------------------------------
    @testset "CO.2: pin package solve_pade ℘ before the pole (z = 0.5)" begin
        ball = _wp_equianharmonic_ball(big(-0.5); prec = PREC)
        rad  = Arblib.radius(Arb, ball)
        mid  = BigFloat(Arblib.midpoint(ball))
        upkg, _ = sol(0.5)
        # (1) CERTIFY oracle.
        @test BigFloat(rad) < ε_cert                       # achieved ~2.8e-74
        # (2) PIN package at its Float64 floor.
        tol_pkg = 1e-13
        @test abs(BigFloat(upkg) - mid) < tol_pkg          # achieved ~2.1e-15
        # Equivalent containment form (inflate the certified ball by tol_pkg
        # then assert it contains the package value) -- the two-step pattern
        # restated; both must agree.
        infl = Arb(prec = PREC); Arblib.set!(infl, ball)
        Arblib.add_error!(infl, Arb(tol_pkg; prec = PREC))
        @test Arblib.contains(infl, Arb(upkg; prec = PREC))
    end

    # -----------------------------------------------------------------------
    # CO.3 — pin the package solve_pade ℘ value PAST the pole.  The same Padé,
    # evaluated past the bridged pole at z = 1, where the rational-approx
    # error is larger (~1e-7).  u(1.4) = ℘(1.4 - 1; 0, 2).  Demonstrates the
    # ε_cert / tol_pkg separation widens, NOT the oracle radius.
    # -----------------------------------------------------------------------
    @testset "CO.3: pin package solve_pade ℘ past the bridged pole (z = 1.4)" begin
        ball = _wp_equianharmonic_ball(big(0.4); prec = PREC)
        rad  = Arblib.radius(Arb, ball)
        mid  = BigFloat(Arblib.midpoint(ball))
        upkg, _ = sol(1.4)
        # (1) CERTIFY oracle -- radius STILL ~1e-74 though the package error grew.
        @test BigFloat(rad) < ε_cert                       # achieved ~4.0e-74
        # (2) PIN package at the post-pole Padé-error floor.
        tol_pkg = 1e-6
        @test abs(BigFloat(upkg) - mid) < tol_pkg          # achieved ~1.3e-7
    end

    # -----------------------------------------------------------------------
    # CO.4 — cross-check the certified oracle against the ODE invariant the
    # whole construction must satisfy: (℘')² = 4℘³ - g₂℘ - g₃ with g₂=0,g₃=2,
    # certified entirely in ball arithmetic.  Computes ℘ and ℘' as Acb balls,
    # asserts the residual ball CONTAINS zero with a provably-tiny radius.
    # This guards the homogeneity rescale itself (a mis-scaled λ breaks it).
    # -----------------------------------------------------------------------
    @testset "CO.4: certified ℘ satisfies (℘')² = 4℘³ - 2 (ball residual ∋ 0)" begin
        θ   = big(π) / 3
        τ   = Acb(cos(θ), sin(θ); prec = PREC)
        G2  = Acb(prec = PREC); G3 = Acb(prec = PREC)
        Arblib.elliptic_invariants!(G2, G3, τ, prec = PREC)
        λ   = (G3 / Acb(2; prec = PREC))^(Acb(1) / 6)
        z   = Acb(big(-0.5); prec = PREC)
        P   = Acb(prec = PREC); Pp = Acb(prec = PREC)
        Arblib.elliptic_p!(P, z / λ, τ, prec = PREC)
        Arblib.elliptic_p_prime!(Pp, z / λ, τ, prec = PREC)
        P  = λ^(-2) * P                                     # ℘   (homogeneity)
        Pp = λ^(-3) * Pp                                    # ℘'  (extra λ⁻¹ from d/dz)
        resid = Pp^2 - (Acb(4) * P^3 - Acb(2))             # = (℘')² - (4℘³ - 2)
        rr = Arb(prec = PREC); Arblib.set!(rr, Arblib.realref(resid))
        @test Arblib.contains_zero(rr)                      # residual encloses 0
        @test BigFloat(Arblib.radius(Arb, rr)) < big"1e-30" # provably tiny
    end

end # @testset Certified ball oracles

# ----------------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 4) — BOTH directions, verified
# 2026-06-09, restored before commit (`git diff` of test + src clean).
#
# Direction (a) — PIN goes RED.  Perturb the pinned package value so it leaves
#   the certified midpoint by more than tol_pkg.  Concretely, in CO.2 replace
#       upkg, _ = sol(0.5)
#   with
#       upkg, _ = sol(0.5); upkg = upkg * (1 + 1e-10)
#   The package value now sits 4.0e-10 from the certified midpoint while
#   tol_pkg = 1e-13.  VERIFIED BITE: CO.2 step-(2) `abs(pkg - mid) < tol_pkg`
#   fails (4.0e-10 > 1e-13) AND the equivalent `contains(infl, pkg)` fails
#   (inflation by 1e-13 cannot reach a 4.0e-10 displacement).  The certify
#   step (1) stays GREEN -- the oracle is untouched -- exactly proving the two
#   steps are independent.  Restored.
#
# Direction (b) — CERTIFY goes RED.  Compute the oracle ball at LOW precision
#   so its radius blows past ε_cert.  Concretely, lower the working precision
#   that every call site passes:
#       PREC = 256   →   PREC = 24
#   (CO.4 builds its ball inline with PREC too, so this one knob drives all
#   four testsets.)  At 24-bit precision Arb's ℘ enclosure has radius ~1e-5,
#   far wider than ε_cert = 1e-40.  VERIFIED BITE (7 fails): CO.1 and CO.3 go
#   fully RED, CO.2's certify step `radius < ε_cert` fails, and CO.4's
#   `radius < 1e-30` fails.  Crucially CO.2's pin step (2) STILL passes (the
#   midpoint is still ~right) and CO.4's `contains_zero` STILL passes — again
#   proving the two assertions are orthogonal: (b) kills the certification,
#   not the pin.  Restored to PREC = 256.
#
#   (NOTE: mutating only the `_wp_equianharmonic_ball` DEFAULT `prec` does NOT
#   bite, because all call sites pass `prec = PREC` explicitly — the correct
#   mutation is the `PREC` binding above.  Recorded so a future reader does
#   not mis-apply the mutation and conclude the test is weak.)
#
# Both mutations confirmed RED then reverted; no src/ change at any point.
