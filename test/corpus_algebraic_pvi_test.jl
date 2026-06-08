# ============================================================================
# corpus_algebraic_pvi_test.jl  —  Corpus Family F: algebraic PVI (CPvi)
# bead padetaylor-by02, epic padetaylor-25og
#
# WHAT THIS FILE PINS
#   The FIRST NON-self-referential PVI oracle in the suite.  Every v1 PVI
#   figure value is self-pinned (the package's own output, capability-map
#   gap #9).  This file pins the package's PVI machinery
#   (`pVI_transformed_rhs(alpha,beta,gamma,delta)`, src/SheetTracker.jl:155,
#   FFW2017 eq. 3 at md:144) against TWO EXTERNAL algebraic closed-form
#   solutions of the standard PVI ODE -- closing gap #9.
#
#   CMS.4 (corpus_multisheet_test.jl) already pins the same RHS closure against
#   an independent mpmath typing of the FFW formula.  CPvi EXTENDS that: the
#   oracle here is a KNOWN ALGEBRAIC SOLUTION (residual==0), not a re-typing of
#   the same equation -- a stronger external check.
#
# GROUND TRUTH each case embodies  (sympy SYMBOLIC residual==0 gate +
#   mpmath dps=50 numeric pins; external/probes/corpus-oracles/algebraic-pvi/
#   capture.py).  The STANDARD PVI ODE used (auditable wiring):
#
#     u'' = (1/2)(1/u + 1/(u-1) + 1/(u-z)) (u')^2
#           - (1/z + 1/(z-1) + 1/(u-z)) u'
#           + u(u-1)(u-z)/(z^2 (z-1)^2) *
#             ( a + b*z/u^2 + g*(z-1)/(u-1)^2 + d*z(z-1)/(u-z)^2 )
#
#     CPvi.1  u(z) = z^(-1/2)        solves PVI (a,b,g,d) = (0,   0, 1/8, -5/8)
#     CPvi.2  u(z) = (sqrt(z)+1)/2   solves PVI (a,b,g,d) = (1/2, 0, 1/8, -5/8)
#
#   Both have a branch point ONLY at z=0 -> two sheets +-sqrt; u(4)=+1/2 for
#   CPvi.1 (other sheet -1/2), u(2+3i)=0.46432545...-0.24849944...i.
#
#   THE RHS-vs-ALGEBRAIC-u'' PIN.  `pVI_transformed_rhs` is the zeta-PLANE
#   closure `(zeta, w, wp) -> w''` under z=e^zeta, w(zeta)=u(z):
#       zeta = log z,  w = u(z),  wp = dw/dzeta = z*u'(z),
#       w'' = d^2 w/d zeta^2 = z*u'(z) + z^2*u''(z).
#   The test feeds the package the algebraic state (zeta=log z, w=u, wp=z*u')
#   and asserts its w'' equals the ANALYTIC z*u' + z^2*u'' from the closed
#   form -- on BOTH sheets (+- chooses the sqrt branch) at z=4 and z=2+3i.
#   This pins the (a,b,g,d) wiring AND the e^zeta transform Jacobian against an
#   external algebraic solution.  The minus-sheet pins are independent, NOT
#   negations of the plus-sheet (verified: CPvi.2 minus u=-0.337-0.448i).
#
# TOLERANCES (justified by METHOD)
#   1e-13: two independent FLOAT64 typings of the same analytic quantity (the
#     package closure vs the closed-form z*u'+z^2*u''); round-off ~1e-15, the
#     transform's exp/log + the nonlinear param block inflate it to ~1e-16
#     (measured max |diff| over all 8 states < 1.1e-16).  1e-13 is a generous,
#     non-band-aid Float64 floor (matches CMS.4's PVI-zeta pin).
#
# LIMITATION (verified; docs/test_corpus/02_corpus_extension_plan.md discovery
#   #6, :89-97):  the branch point z=0 of BOTH solutions COINCIDES with PVI's
#   own fixed singularity.  So these do NOT test circumambulation of a GENERIC
#   (movable, off-{0,1,inf}) branch -- that needs a Dubrovin-Mazzocco /
#   Lisovyy-Tykhyy algebraic solution with a tabulated named-sheet value at an
#   off-singular z0, which is NOT in-repo (tracked separately, gap reference-
#   acquisition bead, blocked).  This file pins the RHS WIRING, not branch
#   monodromy.
#
# REFERENCES
#   docs/test_corpus/02_corpus_extension_plan.md Family F (:274-289),
#     discovery #6 (:89-97);  capability-map gap #9 (PVI self-reference).
#   references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#     FFW2017_painleve_riemann_surfaces_preprint.md  (PVI + the z=e^zeta map).
#   src/SheetTracker.jl:155 (pVI_transformed_rhs), :30-35 (frame doc).
#   test/corpus_multisheet_test.jl  CMS.4 (the call idiom this extends).
#
# MUTATION-PROOF: footer block (Rule 4) -- M-oracle (flip a pinned WPP) and
#   M-impl (perturb the gamma term in src/SheetTracker.jl pVI_transformed_rhs)
#   both ACTUALLY EXECUTED, confirmed RED, restored byte-for-byte.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor.SheetTracker: pVI_transformed_rhs

include(joinpath(@__DIR__, "_oracle_corpus_algebraic_pvi.jl"))

# zeta-frame state from the closed form on a chosen sqrt sheet (s = +-1):
#   returns (zeta = log z, w = u, wp = z*u').
# CPvi.1  u = z^(-1/2) = 1/(s*sqrt z),  u' = -1/(2 z^(3/2)),  u'' = 3/(4 z^(5/2)).
function state_cpvi1(z, s)
    rt  = s * sqrt(z)
    u   = 1 / rt
    up  = -0.5 * rt / z^2          # -1/2 * z^(-3/2)  (rt = s*z^(1/2))
    return (log(z), u, z * up)
end
# CPvi.2  u = (s*sqrt z + 1)/2,  u' = 1/(4 sqrt z),  u'' = -1/(8 z^(3/2)).
function state_cpvi2(z, s)
    rt  = s * sqrt(z)
    u   = (rt + 1) / 2
    up  = 0.25 / rt
    return (log(z), u, z * up)
end

const TOL = 1e-13

@testset "Corpus: algebraic PVI (CPvi)" begin

    # -----------------------------------------------------------------------
    # CPvi.1  u = z^(-1/2) solves PVI (0,0,1/8,-5/8).  Pin the zeta-frame RHS
    #         closure vs the analytic w'' = z*u' + z^2*u'' on BOTH sheets at
    #         z=4 (real) and z=2+3i (complex).  Also re-pin the closed-form
    #         state (u, wp) to mpmath dps=50 so the input is independently
    #         anchored (not back-computed from the package).
    # -----------------------------------------------------------------------
    @testset "CPvi.1: u=z^(-1/2), PVI(0,0,1/8,-5/8) RHS == algebraic w''" begin
        f = pVI_transformed_rhs(0.0, 0.0, 1 / 8, -5 / 8)
        cases = (
            (complex(4.0, 0.0), +1.0, CPVI1_Z4_PLUS_U,  CPVI1_Z4_PLUS_WP,  CPVI1_Z4_PLUS_WPP),
            (complex(4.0, 0.0), -1.0, CPVI1_Z4_MINUS_U, CPVI1_Z4_MINUS_WP, CPVI1_Z4_MINUS_WPP),
            (complex(2.0, 3.0), +1.0, CPVI1_Z2P3I_PLUS_U,  CPVI1_Z2P3I_PLUS_WP,  CPVI1_Z2P3I_PLUS_WPP),
            (complex(2.0, 3.0), -1.0, CPVI1_Z2P3I_MINUS_U, CPVI1_Z2P3I_MINUS_WP, CPVI1_Z2P3I_MINUS_WPP),
        )
        for (z, s, u_ref, wp_ref, wpp_ref) in cases
            (ζ, w, wp) = state_cpvi1(z, s)
            # Input state independently anchored to the mpmath dps=50 oracle.
            @test isapprox(w,  u_ref;  atol = TOL, rtol = TOL)
            @test isapprox(wp, wp_ref; atol = TOL, rtol = TOL)
            # THE non-self-referential pin: package RHS == analytic w''.
            @test isapprox(f(ζ, w, wp), wpp_ref; atol = TOL, rtol = TOL)
        end
        # Headline brief values (sheet semantics): u(4) = +1/2 / -1/2.
        @test isapprox(state_cpvi1(complex(4.0, 0.0), +1.0)[2],  0.5; atol = TOL)
        @test isapprox(state_cpvi1(complex(4.0, 0.0), -1.0)[2], -0.5; atol = TOL)
    end

    # -----------------------------------------------------------------------
    # CPvi.2  u = (sqrt(z)+1)/2 solves PVI (1/2,0,1/8,-5/8).  Companion with a
    #         NONZERO alpha (exercises the +alpha branch of the param block,
    #         which CPvi.1's alpha=0 cannot).  Minus-sheet states are genuinely
    #         distinct (not negations of plus), so each pins a separate value.
    # -----------------------------------------------------------------------
    @testset "CPvi.2: u=(sqrt(z)+1)/2, PVI(1/2,0,1/8,-5/8) RHS == algebraic w''" begin
        f = pVI_transformed_rhs(1 / 2, 0.0, 1 / 8, -5 / 8)
        cases = (
            (complex(4.0, 0.0), +1.0, CPVI2_Z4_PLUS_U,  CPVI2_Z4_PLUS_WP,  CPVI2_Z4_PLUS_WPP),
            (complex(4.0, 0.0), -1.0, CPVI2_Z4_MINUS_U, CPVI2_Z4_MINUS_WP, CPVI2_Z4_MINUS_WPP),
            (complex(2.0, 3.0), +1.0, CPVI2_Z2P3I_PLUS_U,  CPVI2_Z2P3I_PLUS_WP,  CPVI2_Z2P3I_PLUS_WPP),
            (complex(2.0, 3.0), -1.0, CPVI2_Z2P3I_MINUS_U, CPVI2_Z2P3I_MINUS_WP, CPVI2_Z2P3I_MINUS_WPP),
        )
        for (z, s, u_ref, wp_ref, wpp_ref) in cases
            (ζ, w, wp) = state_cpvi2(z, s)
            @test isapprox(w,  u_ref;  atol = TOL, rtol = TOL)
            @test isapprox(wp, wp_ref; atol = TOL, rtol = TOL)
            @test isapprox(f(ζ, w, wp), wpp_ref; atol = TOL, rtol = TOL)
        end
        # Sheet semantics: u(4) = +3/2 (plus) vs -1/2 (minus) -- NOT a negation.
        @test isapprox(state_cpvi2(complex(4.0, 0.0), +1.0)[2],  1.5; atol = TOL)
        @test isapprox(state_cpvi2(complex(4.0, 0.0), -1.0)[2], -0.5; atol = TOL)
    end

    # -----------------------------------------------------------------------
    # CPvi.3  alpha-wiring discriminator.  CPvi.1 (alpha=0) and CPvi.2
    #         (alpha=1/2) carry the SAME (b,g,d)=(0,1/8,-5/8).  At a fixed
    #         interior state the two closures must differ by exactly the alpha
    #         term  d(w'')/d alpha = w(w-1)(w-e^zeta)/(e^zeta-1)^2  (the
    #         param-block coefficient of alpha).  Pins that alpha is wired into
    #         the right place, independent of any solution.
    # -----------------------------------------------------------------------
    @testset "CPvi.3: alpha enters via w(w-1)(w-e^zeta)/(e^zeta-1)^2" begin
        f0 = pVI_transformed_rhs(0.0, 0.0, 1 / 8, -5 / 8)
        fα = pVI_transformed_rhs(1 / 2, 0.0, 1 / 8, -5 / 8)
        ζ, w, wp = complex(0.7, 0.3), complex(0.42, 0.1), complex(-0.05, -0.02)
        eζ = exp(ζ)
        dα = w * (w - 1) * (w - eζ) / (eζ - 1)^2     # analytic d(w'')/d alpha
        @test isapprox(fα(ζ, w, wp) - f0(ζ, w, wp), 0.5 * dα; atol = TOL, rtol = TOL)
    end

end # @testset Corpus algebraic PVI (CPvi)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
#   M-oracle.  Flip the sign of CPVI1_Z2P3I_PLUS_WPP in
#     test/_oracle_corpus_algebraic_pvi.jl (real part 0.11608... -> -0.11608).
#     Re-run standalone: CPvi.1's f(zeta,w,wp)==wpp_ref pin goes RED (the
#     package RHS still returns the correct +0.11608..., no longer matching the
#     flipped pin).  Restored byte-for-byte; GREEN.
#
#   M-impl.  In src/SheetTracker.jl pVI_transformed_rhs, perturb the gamma
#     term of the param block: `γ * eζ_m1 / w_m1^2` -> `γ * eζ_m1 / w_m1^2 * 2`
#     (double the gamma contribution).  Re-run standalone: BOTH CPvi.1 and
#     CPvi.2 RHS==w'' pins go RED at z=4 and z=2+3i on both sheets (gamma=1/8
#     is nonzero for both solutions; the doubled term shifts w'' by a
#     measurable O(1e-1) amount).  CPvi.3 stays GREEN (it differs only in
#     alpha, gamma cancels), confirming CPvi.3 is an independent alpha-only
#     probe.  Restored byte-for-byte; `git diff src/` empty; GREEN.
#
# STANDALONE RUN:
#   julia --project=. test/corpus_algebraic_pvi_test.jl
# ============================================================================
