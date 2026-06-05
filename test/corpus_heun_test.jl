# test/corpus_heun_test.jl
#
# ============================================================================
# Ground-truth corpus: HEUN functions, breaking the eps=-3 monoculture
# (bead padetaylor-h5ho, epic padetaylor-p1v0).
#
# WHAT THIS FILE PINS
# -------------------
# The committed Heun suite (test/heun_test.jl) + its oracle
# (external/probes/heun-oracle/oracles.txt, 42 records) live almost entirely
# on ONE parameter slice: HeunG(2,1,1,2,3,4) and HeunC(1,1,2,3,-1).  Per the
# capability map (B18) the *coefficient* oracle is an eps=-1 / eps=-3
# MONOCULTURE.  This file breaks it: closed-form reductions (2F1, 1/(1-z),
# Heun polynomial, trivial), the off-monoculture parameter regimes
# (complex a, |a|<1, |a|>1 with large q), and HeunC across eps=-3/+2 and
# fractional gamma.  Every value is verified by an INDEPENDENT mpmath
# recurrence AND Mathematica wolframscript HeunG/HeunC
# (external/probes/corpus-oracles/heun/capture.py) — NEVER by calling the
# Julia heun_g/heun_c (no self-reference).
#
# CLOSED-FORM REDUCTIONS (verified BEFORE use — Law 1):
#   CHN.1  HeunG(2,2,1,2,3/2,1;z)   = 2F1(1/2,1;3/2;2z-z^2)   (DLMF/Maier)
#   CHN.2  HeunG(2,4,1,4,2,2;z)     = 1/(1-z), all c_n = 1     (DLMF 31.7.2)
#   CHN.3a HeunG(2,q*,-1,2,3,4;z)   = 1 + (q*/6)z              (Heun poly, 31.5)
#   CHN.3b HeunG(a,0,0,b,g,d;z)     = 1                        (trivial)
# OFF-MONOCULTURE NUMERIC (mpmath recurrence + wolframscript):
#   CHN.4a complex a = 1/2+i/2 ; CHN.4b |a|=1/2 ; CHN.4c a=3,q=10
#   CHN.6  HeunC eps=-3, eps=+2, gamma=1/2
#
# ERRATA found + corrected here (docs/test_corpus/ERRATA.md, appended):
#   * heunc-epsilon-neg3 @ z=0.9: catalogue -10.40476 is a TRUNCATION artefact;
#     converged & Wolfram-confirmed value = -13.5052005053212 (CHN.6).
#   * heung-complex-a @ z=0.3: catalogue real part wrong in the 7th figure;
#     exact-rational a gives 1.10598130983985..., Wolfram-confirmed (CHN.4a).
#
# TOLERANCES are justified by METHOD + REGIME, never fitted:
#   * Frobenius-only regime (|z| <= epsilon_start=0.05): n_taylor=60 series sum
#     reaches ~1e-14 — but our eval points are mostly |z|>=0.1, so heun_g/heun_c
#     WALK via the path-network.  The walk's local-Pade continuation past the
#     z=0/z=1/z=a singular points carries ~1e-11..1e-13 on the principal disc;
#     points crowding a singular point (z=0.9 near z=1; small-a near |a|=0.5)
#     loosen to 1e-9.  Pin tolerances at 10x the observed gap (smoke probe
#     2026-06-05), not at the edge.
#   * Exact reductions (1/(1-z), polynomial, trivial) hit the walk floor too,
#     so they are pinned at 1e-12 not machine-eps (the WALK, not the series,
#     is exercised for |z|>=0.1).
#
# MUTATION-PROOF (Rule 4): footer block — a deliberate perturbation of the
# HeunG Frobenius recurrence (_heun_g_frobenius_at_0) was confirmed to RED a
# closed-form-reduction assertion, then restored byte-clean (git diff empty).
# ============================================================================

using Test
using PadeTaylor: heun_g, heun_c

include(joinpath(@__DIR__, "_oracle_corpus_heun.jl"))

# 2F1 reduction oracle as a closed form, recomputed in-test from a series so
# the identity is asserted, not merely pinned: 2F1(1/2,1;3/2;w) = atanh(sqrt w)/sqrt w
# (DLMF 15.4.x; w = 2z - z^2).  Used only as an independent cross-check.
_f2f1_half1_3half(w) = atanh(sqrt(complex(w))) / sqrt(complex(w))

@testset "Corpus: Heun (CHN) — breaking the eps=-3 monoculture" begin

    # ----------------------------------------------------------------------
    # CHN.1  heung-2f1-reduction : the BEST oracle (fully closed-form).
    # ----------------------------------------------------------------------
    @testset "CHN.1  HeunG(2,2,1,2,3/2,1;z) = 2F1(1/2,1;3/2;2z-z^2)" begin
        # First assert the reduction IDENTITY itself (independent of heun_g):
        for z in (0.1, 0.3, 0.5)
            @test isapprox(_f2f1_half1_3half(2z - z^2),
                           z == 0.3 ? HEUNG_2F1_Z0p3 :
                           z == 0.5 ? HEUNG_2F1_Z0p5 : 1.071704836826995399;
                           rtol = 1e-13)
        end
        # Then pin heun_g to the reduction.
        @test isapprox(real(heun_g(2,2,1,2,3//2,1; z=0.3+0im)), HEUNG_2F1_Z0p3; rtol=1e-12)
        @test isapprox(real(heun_g(2,2,1,2,3//2,1; z=0.5+0im)), HEUNG_2F1_Z0p5; rtol=1e-12)
        # z=0.7 crowds the branch point at z=1 → walk floor loosens.
        @test isapprox(real(heun_g(2,2,1,2,3//2,1; z=0.7+0im)), HEUNG_2F1_Z0p7; rtol=1e-9)
    end

    # ----------------------------------------------------------------------
    # CHN.2  heung-pole-at-z1 : HeunG(2,4,1,4,2,2;z) = 1/(1-z), all c_n=1.
    # ----------------------------------------------------------------------
    @testset "CHN.2  HeunG(2,4,1,4,2,2;z) = 1/(1-z), exact-rational" begin
        # On the principal disc, heun_g reproduces 1/(1-z) to the walk floor.
        @test isapprox(heun_g(2,4,1,4,2,2; z=0.3+0im), 1/(1-0.3)+0im; rtol=1e-12)
        @test isapprox(heun_g(2,4,1,4,2,2; z=0.5+0im), 2.0+0im;       rtol=1e-12)
        # z=0.9 crowds the simple pole at z=1 (the closed-form value is 10).
        @test isapprox(heun_g(2,4,1,4,2,2; z=0.9+0im), 10.0+0im;      rtol=1e-9)
    end

    # ----------------------------------------------------------------------
    # CHN.3a heung-polynomial-degree1 : terminating series, linear closed form.
    # ----------------------------------------------------------------------
    @testset "CHN.3a  Heun polynomial Hp_{1,m}(z) = 1 + (q*/6)z" begin
        # Eigenvalues q* = -6 +- 2 sqrt(6) terminate the series at c_2 = 0.
        qp, qm = HEUNG_POLY_QP, HEUNG_POLY_QM
        @test isapprox(qp, -6 + 2*sqrt(6); rtol=1e-14)
        @test isapprox(qm, -6 - 2*sqrt(6); rtol=1e-14)
        # closed form 1 + (q*/6) z reproduced by heun_g:
        @test isapprox(real(heun_g(2,qp,-1,2,3,4; z=0.5+0im)), 1 + qp/6*0.5; rtol=1e-12)
        @test isapprox(real(heun_g(2,qp,-1,2,3,4; z=0.5+0im)), HEUNG_POLY_P_Z0p5; rtol=1e-12)
        @test isapprox(real(heun_g(2,qp,-1,2,3,4; z=0.9+0im)), HEUNG_POLY_P_Z0p9; rtol=1e-11)
        @test isapprox(real(heun_g(2,qm,-1,2,3,4; z=0.5+0im)), HEUNG_POLY_M_Z0p5; rtol=1e-12)
    end

    # ----------------------------------------------------------------------
    # CHN.3b heung-trivial-alpha0-q0 : HeunG(a,0,0,b,g,d;z) = 1 for all z.
    # ----------------------------------------------------------------------
    @testset "CHN.3b  trivial HeunG(a,0,0,b,g,d;z) = 1" begin
        # c_1 = q/(a g) = 0 and P_j = (j-1+0)(j-1+b) kills all higher c_j.
        for zt in (0.3+0im, 0.7+0im, 0.9+0im)
            @test isapprox(heun_g(2,0,0,2,3,4; z=zt), 1.0+0im; rtol=1e-12)
        end
        # different a confirms it is a-independent (a=4, the oracles.txt slice).
        @test isapprox(heun_g(4,0,0,2,2,2; z=0.6+0im), 1.0+0im; rtol=1e-12)
    end

    # ----------------------------------------------------------------------
    # CHN.4a heung-complex-a : a=1/2+i/2 off the real axis [ERRATA-corrected].
    # ----------------------------------------------------------------------
    @testset "CHN.4a  HeunG complex a = 1/2+i/2 (ERRATA-corrected oracle)" begin
        v = heun_g(complex(0.5,0.5),1,1,2,3,4; z=0.3+0im)
        @test isapprox(v, HEUNG_COMPLEX_A_Z0p3; rtol=1e-11)
        # Anti-regression guard: must NOT equal the wrong catalogue value.
        # Correct vs wrong differ by ~1.8e-7 (mostly in Im); the computed value
        # sits 2.7e-15 from correct, so atol=1e-8 cleanly separates them.
        @test !isapprox(v, 1.1059813320247579739 - 0.10592175241504076687im; atol=1e-8)
    end

    # ----------------------------------------------------------------------
    # CHN.4b heung-small-a : |a| = 1/2 < 1, singular point inside unit disc.
    # ----------------------------------------------------------------------
    @testset "CHN.4b  HeunG |a|=1/2 (singular point inside unit disc)" begin
        @test isapprox(real(heun_g(1//2,1,1,2,3,4; z=0.1+0im)),  HEUNG_SMALL_A_Z0p1;  rtol=1e-12)
        @test isapprox(real(heun_g(1//2,1,1,2,3,4; z=0.3+0im)),  HEUNG_SMALL_A_Z0p3;  rtol=1e-11)
        # z=0.45 is past the |a|=0.5 singular point → walk detours it; 1e-9.
        @test isapprox(real(heun_g(1//2,1,1,2,3,4; z=0.45+0im)), HEUNG_SMALL_A_Z0p45; rtol=1e-9)
    end

    # ----------------------------------------------------------------------
    # CHN.4c heung-large-q-a3 : a=3 (|a|>1), q=10 (large) — oscillatory coeffs.
    # ----------------------------------------------------------------------
    @testset "CHN.4c  HeunG a=3, q=10 (|a|>1, large q)" begin
        @test isapprox(real(heun_g(3,10,1,2,2,3; z=0.25+0im)), HEUNG_LARGE_Q_Z0p25; rtol=1e-12)
        @test isapprox(real(heun_g(3,10,1,2,2,3; z=0.5+0im)),  HEUNG_LARGE_Q_Z0p5;  rtol=1e-12)
    end

    # ----------------------------------------------------------------------
    # CHN.6  HeunC across eps and fractional gamma (NOT eps=-1!).
    # ----------------------------------------------------------------------
    @testset "CHN.6  HeunC eps=-3 / eps=+2 / gamma=1/2" begin
        # eps=-3 (catalogue z=0.9 ERRATA-corrected to -13.5052005053212).
        @test isapprox(real(heun_c(1,1,2,3,-3; z=0.1+0im)), HEUNC_NEG3_Z0p1; rtol=1e-12)
        @test isapprox(real(heun_c(1,1,2,3,-3; z=0.5+0im)), HEUNC_NEG3_Z0p5; rtol=1e-11)
        v09 = heun_c(1,1,2,3,-3; z=0.9+0im)
        @test isapprox(real(v09), HEUNC_NEG3_Z0p9; rtol=1e-9)   # crowds z=1
        @test !isapprox(real(v09), -10.40, atol=1.0)            # anti-ERRATA guard
        # eps=+2 (positive — sign-convention probe).
        @test isapprox(real(heun_c(1,1,2,3,2; z=0.1+0im)), HEUNC_POS2_Z0p1; rtol=1e-12)
        @test isapprox(real(heun_c(1,1,2,3,2; z=0.5+0im)), HEUNC_POS2_Z0p5; rtol=1e-11)
        # gamma=1/2 (fractional Frobenius exponent — u'(0) = -q/gamma = -2).
        @test isapprox(real(heun_c(1,1,1//2,3,-1; z=0.1+0im)), HEUNC_FRACG_Z0p1; rtol=1e-12)
        @test isapprox(real(heun_c(1,1,1//2,3,-1; z=0.5+0im)), HEUNC_FRACG_Z0p5; rtol=1e-11)
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored.
#
# M-impl (aimed at the HeunG Frobenius recurrence, per the brief):
#   In src/Heun.jl `_heun_g_frobenius_at_0`, perturb the coefficient
#   recurrence numerator `(Q_j + q)` -> `(Q_j + q) * 1.0000001` (a 1e-7
#   relative perturbation of every c_{j+1}).  This poisons the IC tuple
#   sampled at z0 = 0.05 that seeds every path-network walk.
#   RESULT: CHN.1 (2F1 reduction, z=0.3/0.5), CHN.2 (1/(1-z)), CHN.3a
#   (Heun polynomial), CHN.4b/4c all went RED — the closed-form reductions
#   are the sharpest detectors because their oracle is exact.  The CHN.3b
#   trivial case (c_1 = c_2 = ... = 0) stayed GREEN, confirming it does NOT
#   exercise the perturbed body (all c_j after c_0 are zero) — a true
#   coverage fact, recorded here, not a gap.
#   Restored src/Heun.jl byte-for-byte (verified `git diff src` empty);
#   suite GREEN.
#
# M-oracle (CHN.6, the ERRATA anti-regression): replace the corrected
#   HEUNC_NEG3_Z0p9 = -13.5052... with the WRONG catalogue -10.40476.
#   RESULT: CHN.6 went RED at the z=0.9 assertion (got -13.5052, expected
#   -10.40476) and ONLY there.  Confirms the test pins the Wolfram-confirmed
#   truth, not the catalogue artefact.  Restored.
#
# Run standalone with:
#   julia --project=. test/corpus_heun_test.jl
# ============================================================================
