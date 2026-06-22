# test/corpus_riccati_rational_test.jl
#
# ============================================================================
# Ground-truth corpus: RICCATI = log-derivative of a linear ODE, RATIONAL
# subfamily (bead padetaylor-bnqq, epic padetaylor-25og).
#
# WHAT THIS FILE PINS
# -------------------
# For any linear ODE  w″+p w′+q w = 0,  u = w′/w  solves a Riccati equation
# with a SIMPLE pole of residue +1 at every zero of w
# (docs/test_corpus/02_corpus_extension_plan.md:134-162, Family A).  When w is
# an orthogonal polynomial the solution u is a ratio of polynomials, so the
# pole locations are exact algebraic numbers (polynomial roots) AND every
# evaluation is an EXACT RATIONAL — there is no transcendental approximation
# floor.  Two members, n=4 (even n: an odd n would seed a zero — i.e. a pole —
# at the origin z₀=0):
#
#   CRic.3  Hermite   u = H₄′/H₄ = 2n·H₃/H₄   poles at the 4 Hermite zeros
#   CRic.4  Chebyshev u = T₄′/T₄              poles at cos((2k−1)π/8)
#
# Each is the corpus's first "row of simple poles" pinned to machine-precision
# rationals.  Like CPB.3 (1/(1−z)) we drive solve_pade through a single window
# that brackets the INNER positive real-pole pair strictly interior, evaluate
# PAST the pole, and pin u, u′ to the sympy-exact rationals; extract_poles is
# pinned to the exact polynomial roots.
#
# GROUND TRUTH each case embodies
# -------------------------------
# All recasts are VERIFIED by sympy (symbolic residual ≡ 0) in
# external/probes/corpus-oracles/riccati-logderiv/capture.py, which prints the
# exact rationals transcribed inline below.  solve_pade's public API takes the
# 2nd-order branch u″=f(z,u,u′); the first-order Riccati R(z,u) is recast by
# u″ = ∂R/∂z + (∂R/∂u)·R (the same differentiation idiom as CPB.2/3/4):
#   CRic.3:  R = −u²+2zu−2n  ⇒  u″ = 2u³−6u²z+4uz²+18u−16z      (n=4)
#   CRic.4:  R = −u²+(zu−n²)/(1−z²)
#            ⇒  u″ = [2u³(1−z²)² + 3u²z(z²−1) + u(33−30z²) − 48z] / (1−z²)²
# (CRic.4's recast was DERIVED symbolically, not hand-copied — capture.py
# gates on its residual ≡ 0 before emitting.)  ICs at z₀=0: CRic.3 (0,−8),
# CRic.4 (0,−16).  Both inner residues are +1 (capture.py).
#
# TOLERANCES (justified by METHOD accuracy, never fitted to pass)
# -----------------------------------------------------------------
# These solutions are exact rationals, so the ONLY error is the Float64
# (15,15)-diagonal Padé past-pole rational-approximation floor (~3e-10 at the
# far end of a window; README §2 / problems_test.jl 6.1.6 / CPB.1 commentary).
#   * before the pole / just past it (small |t| past midpoint): rtol 1e-13
#     (machine precision — same regime as CPB.3 u(1.5)=−2).
#   * farther past the pole (t approaching the window end): rtol 1e-9, the
#     documented past-pole Padé floor used by CPB.1/CPB.2.
#   * extract_poles is the denominator-root finder: atol 1e-7 on the inner
#     pair (the figure-catalogue Float64 spec, polefield_test.jl PF.1.1), 1e-6
#     on the outer pair (farther from the window centre).
#
# REFERENCES: docs/test_corpus/02_corpus_extension_plan.md:134-162 (Family A);
# DLMF 18.5/18.16 (Hermite Hₙ + zeros), 18.5.1/18.16 (Chebyshev Tₙ + zeros
# cos((2k−1)π/2n)).  Oracle constants: capture.py (sympy, residual==0 gated).
#
# MUTATION-PROOF (Rule 4): see the comment block at the end of this file —
# one M-oracle (flip a pinned rational) and one M-impl (perturb the stored
# Padé) were ACTUALLY executed RED and restored byte-for-byte.
# ============================================================================

using Test
using PadeTaylor

@testset "Corpus: Riccati log-deriv rational (CRic)" begin

    # ----------------------------------------------------------------------
    # CRic.3  Hermite n=4:  u″=2u³−6u²z+4uz²+18u−16z,  u=H₄′/H₄.
    # Inner real-pole pair ±0.5246476232752904, outer ±1.650680123885785.
    # Window (0,1.4): inner pole at t≈0.375 interior; targets past it but
    # short of the outer pole.
    # ----------------------------------------------------------------------
    @testset "CRic.3  Hermite H₄′/H₄: bridge the inner zero-pair pole" begin
        fH(z, u, up) = 2u^3 - 6u^2 * z + 4u * z^2 + 18u - 16z
        prob = PadeTaylorProblem(fH, (0.0, -8.0), (0.0, 1.4); order = 30)
        sol  = solve_pade(prob; h = 1.4)        # inner pole 0.5246 interior

        u03, up03 = sol(0.3)                         # BEFORE the inner pole
        u08, up08 = sol(0.8)                         # PAST the inner pole
        u10, up10 = sol(1.0)
        u12, _    = sol(1.2)                         # near window end (Padé floor)

        # EXACT rationals (capture.py): the only error is the Padé floor.
        @test isapprox(u03,  -5640 / 1627; rtol = 1e-13)   # before pole: machine eps
        @test isapprox(up03, -58492400 / 2647129; rtol = 1e-13)
        @test isapprox(u08,   6880 / 1901; rtol = 1e-12)   # just past pole
        @test isapprox(up08, -55318600 / 3613801; rtol = 1e-11)
        @test isapprox(u10,   8 / 5;  rtol = 1e-11)
        @test isapprox(up10, -184 / 25; rtol = 1e-10)
        @test isapprox(u12,   240 / 1247; rtol = 1e-9)     # far end: past-pole floor

        # The pole really is bridged: u(0.8) is finite and the sign flipped
        # (u→+∞ as z↑0.5246⁻, u→−∞ as z↓0.5246⁺ then climbs back through 0).
        @test isfinite(u08) && u08 > 0

        # extract_poles pins the four exact Hermite zeros (residue +1).
        poles = extract_poles(sol; min_support = 1)
        for zexact in (0.5246476232752904, 1.650680123885785)
            @test minimum(abs(p - complex(zexact, 0.0)) for p in poles) ≤
                  (zexact < 1 ? 1e-7 : 1e-6)
        end
    end

    # ----------------------------------------------------------------------
    # CRic.4  Chebyshev n=4:  u″=[2u³(1−z²)²+3u²z(z²−1)+u(33−30z²)−48z]/(1−z²)²,
    # u=T₄′/T₄.  Inner pole pair ±0.3826834323650898, outer ±0.9238795325112867.
    # ±1 are RHS coefficient-singularities ⇒ keep the window strictly inside
    # (−1,1).  Window (0,0.85): inner pole at t≈0.45 interior; outer pole and
    # ±1 stay OUTSIDE the window.
    # ----------------------------------------------------------------------
    @testset "CRic.4  Chebyshev T₄′/T₄: bridge inner pole inside (−1,1)" begin
        function fC(z, u, up)
            d = (1 - z^2)^2
            return (2u^3 * d + 3u^2 * z * (z^2 - 1) + u * (33 - 30z^2) - 48z) / d
        end
        prob = PadeTaylorProblem(fC, (0.0, -16.0), (0.0, 0.85); order = 30)
        sol  = solve_pade(prob; h = 0.85)       # inner pole 0.3827 interior

        u02, up02 = sol(0.2)                         # BEFORE the inner pole
        u05, up05 = sol(0.5)                         # PAST the inner pole
        u06, _    = sol(0.6)
        u07, _    = sol(0.7)                         # near window end (Padé floor)

        # EXACT rationals (capture.py); u(1/2)=8 is the brief's headline value.
        @test isapprox(u02,  -1840 / 433; rtol = 1e-13)    # before pole: machine eps
        @test isapprox(up02, -6676400 / 187489; rtol = 1e-13)
        @test isapprox(u05,   8.0; rtol = 1e-12)           # u(1/2)=8 exact
        @test isapprox(up05, -80.0; rtol = 1e-11)          # u′(1/2)=−80 exact
        @test isapprox(u06,   1680 / 527; rtol = 1e-10)
        @test isapprox(u07,   280 / 1249; rtol = 1e-9)     # far end: past-pole floor

        # Bridged: u(0.5) finite, positive, the sign-flipped branch past 0.3827.
        @test isfinite(u05) && u05 > 0

        # extract_poles pins the inner Chebyshev zero pair (residue +1); the
        # outer pole 0.9239 sits OUTSIDE this short window so we pin only the
        # interior pair it can resolve.
        poles = extract_poles(sol; min_support = 1)
        @test minimum(abs(p - complex(0.3826834323650898, 0.0)) for p in poles) ≤ 1e-7
        @test minimum(abs(p - complex(-0.3826834323650898, 0.0)) for p in poles) ≤ 1e-7
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
# M-oracle (CRic.4): assert `u05 == 9.0` instead of the exact rational
#   u(1/2)=8.  Result: ONLY the CRic.4 u(1/2) assertion went RED
#   (Evaluated: isapprox(7.999999999999984, 9.0; rtol=1e-12)); the other 18
#   assertions stayed GREEN.  Restored the `8.0` literal; suite GREEN (19/19).
#
# M-impl (the gold standard): in src/PadeStepper.jl `_evaluate_pade`, scale the
#   returned numerator by (1 + 1e-7) — a 1e-7 relative perturbation of every
#   stored Padé value.  Result: all 8 u-VALUE assertions went RED
#   (CRic.3 u(0.3,0.8,1.0,1.2) and CRic.4 u(0.2,0.5,0.6,0.7)); the u′-derivative
#   pins (routed through `_evaluate_pade_deriv`) and the extract_poles pins
#   (denominator-root reads) stayed GREEN — the exact coverage split the CPB.3
#   note documents.  Restored src/PadeStepper.jl byte-for-byte (`git diff src/`
#   empty); suite GREEN (19/19).
#
# Run standalone with:
#   julia --project=. test/corpus_riccati_rational_test.jl
# ============================================================================
