# test/corpus_multisheet_test.jl
#
# ============================================================================
# Ground-truth corpus: MULTI-SHEET transformed RHS + sheet values (B10)
# (bead padetaylor-5qqr, epic padetaylor-p1v0).  Highest remaining bug-risk
# tier: the FFW2017 multivalued machinery (P~III, P~V in ζ; PVI in ζ and the
# η double-exp plane) carries subtle sign / chain-rule artefacts, and the
# existing tests are Float64-only.  This file's value-add over the committed
# coord_transforms_test.jl / eta_pvi_test.jl / sheet_tracker_test.jl is:
#   (a) a CLOSED-FORM constant-solution residual on the P~III RHS;
#   (b) the exact pole residue A=±2/√γ via Laurent balance on the RHS;
#   (c) Schwarz conjugate-symmetry residual (closed-form = 0);
#   (d) NON-self-referential RHS pins (mpmath FFW-formula vs package closure);
#   (e) a BigFloat-256 recovery test of the η double-exp cancellation (gap #3).
#
# GROUND TRUTH (Law 1): the four transformed equations are cited by line from
#   references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#   FFW2017_painleve_riemann_surfaces_preprint.md
#     :43  P~III  w'' = (1/w)(w')² + (1/4)(αw²+γw³+βe^ζ+δe^{2ζ}/w)
#     :47  P~V    (the eq. matched in src/CoordTransforms.jl)
#     :144 PVI-ζ  eq.(3)
#     :154 PVI-η  eq.(5), incl. the chain-rule "-1" in the (dv/dη) bracket
# and matched against src/CoordTransforms.jl / src/SheetTracker.jl verbatim.
#
# CATALOGUE NOTE: catalogue records pIII-tilde-trivial-constant,
# pIII-tilde-residue-exact, pIII-tilde-conjugate-symmetry,
# pIII/pV-tilde-sheet-oracle, pVI-eta-cancellation, pVI-eta-subexpr-rhs.  The
# sheet-0 numbers + RHS values were INDEPENDENTLY re-derived
# (external/probes/corpus-oracles/multisheet/capture.py) — see the oracle
# file header for why this breaks the catalogue's self-reference.
#
# TOLERANCES are justified by the precision of the comparison:
#   * Float64 RHS vs mpmath FFW-formula: 1e-13 (Float64 round-off ≈ 1e-15,
#     padded for the e^{e^η} double-exp's amplified rounding).
#   * Closed-form constant-solution residual: 1e-13 (a few ulp of exp(ζ/2)).
#   * Laurent residue coefficient at x=1e-5: 1e-4 (the O(x) tail of the
#     ansatz; the leading x^{-3} coeff converges to ±4 as x→0).
#   * BigFloat-256 cancellation real part: rel 1e-20 (mpmath dps=60 ref).
#
# MUTATION-PROOF (Rule 4): footer block; one mutant aimed at the η chain-rule
# "-1" artefact (CMS.5) confirmed RED, then restored.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor.CoordTransforms: pIII_transformed_rhs, pV_transformed_rhs
import PadeTaylor.SheetTracker: pVI_transformed_rhs, pVI_eta_transformed_rhs

include(joinpath(@__DIR__, "_oracle_corpus_multisheet.jl"))

@testset "Corpus: multi-sheet transformed RHS + sheet values (CMS)" begin

    # -----------------------------------------------------------------------
    # CMS.1  P~III exact constant solution: w(ζ)=exp(ζ/2) for (0,0,1,-1).
    #        Closed-form residual check of the RHS implementation itself.
    #        FFW2017...md:43; catalogue pIII-tilde-trivial-constant.
    #        w=exp(ζ/2) ⇒ w'=½exp(ζ/2), w''=¼exp(ζ/2); the RHS must return
    #        exactly ¼exp(ζ/2) at every ζ (the e^{2ζ}/w term cancels γw³).
    # -----------------------------------------------------------------------
    @testset "CMS.1: P~III(0,0,1,-1) constant-solution RHS residual = 0" begin
        f = pIII_transformed_rhs(0.0, 0.0, 1.0, -1.0)
        for ζ in (0.3 + 0.7im, -1.2 + 2.1im, 2.0 + 3.0im, -0.5 + 0.0im)
            w   = exp(ζ / 2)
            wp  = 0.5 * exp(ζ / 2)
            wpp = 0.25 * exp(ζ / 2)            # closed-form w''
            @test isapprox(f(ζ, w, wp), wpp; atol = 1e-13, rtol = 1e-13)
        end
    end

    # -----------------------------------------------------------------------
    # CMS.2  P~III pole residue A=±2/√γ via Laurent balance on the RHS.
    #        catalogue pIII-tilde-residue-exact; FFW2017 Table 1 (md:113).
    #        Algebraic: 2A = A + γA³/4 ⇒ A²=4/γ.  Numerically: the RHS on the
    #        ansatz w=A/x near a pole has leading coefficient rhs·x³ → 2A.
    # -----------------------------------------------------------------------
    @testset "CMS.2: P~III pole residue A = ±2/√γ (γ=1 ⇒ ±2)" begin
        # Closed-form algebraic balance (no floating point).
        for A in (2.0, -2.0)
            @test 2A ≈ A + 1.0 * A^3 / 4          # 2A = A + γA³/4, γ=1
            @test A^2 ≈ 4.0                        # A²=4/γ
        end
        # Anti-regression: a WRONG residue (±1, the z-plane value) does NOT
        # balance — guards against confusing z-plane c_{-1}=±1/√γ with the
        # ζ-plane leading coefficient ±2/√γ.
        @test !(2 * 1.0 ≈ 1.0 + 1.0^3 / 4)
        # Numerical: package RHS leading x^{-3} coeff converges to 2A=±4.
        f = pIII_transformed_rhs(-0.5, -0.5, 1.0, -1.0)   # γ=1
        ζ0 = 0.3 + 0.5im
        for A in (2.0 + 0im, -2.0 + 0im)
            x = 1e-5
            w  = A / x
            wp = -A / x^2
            @test isapprox(f(ζ0 + x, w, wp) * x^3, 2A; atol = 1e-4)
        end
    end

    # -----------------------------------------------------------------------
    # CMS.3  Schwarz conjugate symmetry w(conj ζ) = conj(w(ζ)) for real
    #        params + real ICs.  catalogue pIII-tilde-conjugate-symmetry.
    #        Oracle: mpmath-RK4 integration along the mirror path gives a
    #        residual that is exactly 0 at dps=40 (capture.py).  We assert
    #        the pinned residual is 0 — a closed-form invariant of the ODE.
    # -----------------------------------------------------------------------
    @testset "CMS.3: P~III / P~V Schwarz conjugate-symmetry residual = 0" begin
        @test PIII_CONJ_SYM_ERR == 0.0
        @test PV_CONJ_SYM_ERR   == 0.0
        # Reinforce as a property of the RHS: with real params, the RHS
        # commutes with conjugation — f(conj ζ, conj w, conj w') == conj f(…).
        for (f, ζ, w, wp) in (
                (pIII_transformed_rhs(-0.5, -0.5, 1.0, -1.0),
                 0.4 + 0.9im, 0.6 + 0.2im, 0.1 - 0.3im),
                (pV_transformed_rhs(1.0, -1.0, 1.0, -0.5),
                 0.4 + 0.6im, 1.5 + 0.3im, -0.2 + 0.1im))
            lhs = f(conj(ζ), conj(w), conj(wp))
            @test isapprox(lhs, conj(f(ζ, w, wp)); atol = 1e-13, rtol = 1e-13)
        end
    end

    # -----------------------------------------------------------------------
    # CMS.4  Transformed-RHS interior values: package closure vs independent
    #        mpmath FFW-formula evaluation (NON-self-referential, gap #9).
    #        Two independent typings of FFW eq.(III~/V~/3/5) must agree to
    #        Float64 round-off.  Oracle: capture.py *_RHS_INTERIOR.
    # -----------------------------------------------------------------------
    @testset "CMS.4: P~III/P~V/PVI-ζ/PVI-η RHS vs mpmath FFW-formula" begin
        f3 = pIII_transformed_rhs(-0.5, -0.5, 1.0, -1.0)
        @test isapprox(f3(0.3 + 0.7im, 0.5 + 0.2im, 0.1 - 0.05im),
                       PIII_RHS_INTERIOR; atol = 1e-13, rtol = 1e-13)
        f5 = pV_transformed_rhs(1.0, -1.0, 1.0, -0.5)
        @test isapprox(f5(0.4 + 0.6im, 1.5 + 0.3im, -0.2 + 0.1im),
                       PV_RHS_INTERIOR; atol = 1e-13, rtol = 1e-13)
        fz = pVI_transformed_rhs(4.0, -4.0, 8.0, -8.0)
        @test isapprox(fz(0.7 + 0.3im, 0.42 + 0.1im, -0.05 - 0.02im),
                       PVI_ZETA_RHS_INTERIOR; atol = 1e-12, rtol = 1e-13)
        fe = pVI_eta_transformed_rhs(4.0, -4.0, 8.0, -8.0)
        η  = complex(log(log(10.0)), 0.5)
        @test isapprox(fe(η, 0.429534600325223 + 0.0im,
                          -0.037235820650106481 + 0.0im),
                       PVI_ETA_RHS_INTERIOR; atol = 1e-12, rtol = 1e-13)
    end

    # -----------------------------------------------------------------------
    # CMS.5  η-plane chain-rule "-1" artefact (FFW2017...md:154, the only
    #        place eq.(5) is NOT a structural copy of eq.(3)).  We pin it by
    #        comparing the η-RHS against a hand-built copy that DROPS the -1:
    #        the difference must equal +v' exactly (the -1 multiplies -(…)v'
    #        ⇒ a +1·v' contribution).  This is the chain-rule "smoking gun".
    # -----------------------------------------------------------------------
    @testset "CMS.5: PVI-η chain-rule -1 artefact contributes +v'" begin
        fe = pVI_eta_transformed_rhs(4.0, -4.0, 8.0, -8.0)
        # Same RHS with the -1 removed from the (dv/dη) bracket.
        fe_no1 = (η, v, vp) -> begin
            eη = exp(η); E = exp(eη); E_m1 = E - 1
            v_m1 = v - 1; v_mE = v - E
            first  = (1 / v + 1 / v_m1 + 1 / v_mE) * vp^2 / 2
            second = -(eη * E / E_m1 + eη * E / v_mE) * vp        # NO -1
            param  = 4.0 - 4.0 * E / v^2 + 8.0 * E_m1 / v_m1^2 -
                     8.0 * E * E_m1 / v_mE^2
            third  = v * v_m1 * v_mE * eη^2 / E_m1^2 * param
            first + second + third
        end
        η = complex(log(log(10.0)), 0.5)
        v, vp = 0.4295 + 0.0im, -0.0372 + 0.0im
        # f_with_minus1 - f_without = (+1)·(-(-1))·vp ... net = +vp.
        @test isapprox(fe(η, v, vp) - fe_no1(η, v, vp), vp;
                       atol = 1e-12, rtol = 1e-12)
    end

    # -----------------------------------------------------------------------
    # CMS.6  GAP #3: η double-exp cancellation E-1 = e^{e^η}-1 near the k=1
    #        fixed singularity η=log(2π)+iπ/2+1e-10.  catalogue
    #        pVI-eta-cancellation.  FINDING: Float64 loses the real part by a
    #        factor ~2.25e3; BigFloat-256 RECOVERS it (the formula is sound,
    #        only precision-starved).  We assert both: the Float64 casualty
    #        and the BigFloat-256 recovery.
    # -----------------------------------------------------------------------
    @testset "CMS.6: η double-exp cancellation recovers at BigFloat-256" begin
        ref_re = parse(Float64, ETA_EM1_REAL_STR)   # -1.9739e-19
        ref_im = parse(Float64, ETA_EM1_IMAG_STR)   #  6.2832e-10
        η_re_str = "1.8378770664093454835606594728112352766558"  # log(2π)

        # Float64: the real part is WRONG (sign + magnitude); documented.
        ηf = complex(parse(Float64, η_re_str) + 1e-10, π / 2)
        Ef = exp(exp(ηf)) - 1
        @test isapprox(real(Ef), ETA_EM1_REAL_FLOAT64_WRONG; rtol = 1e-6)
        @test abs((real(Ef) - ref_re) / ref_re) > 1e3   # casualty: ~2250×

        # BigFloat-256: same subtractive form (E=exp(exp η); E-1) recovers.
        setprecision(BigFloat, 256) do
            ηb = Complex{BigFloat}(
                BigFloat(η_re_str) + BigFloat("1e-10"), BigFloat(π) / 2)
            Eb = exp(exp(ηb)) - 1
            ref_re_big = parse(BigFloat, ETA_EM1_REAL_STR)
            ref_im_big = parse(BigFloat, ETA_EM1_IMAG_STR)
            @test abs((real(Eb) - ref_re_big) / ref_re_big) < 1e-20
            @test abs((imag(Eb) - ref_im_big) / ref_im_big) < 1e-20
        end
    end

end # @testset Corpus multi-sheet (CMS)

# ============================================================================
# MUTATION-PROOF (Rule 4) — verified 2026-06-05, then restored byte-clean.
#
#   Mutation M1 — src/SheetTracker.jl pVI_eta_transformed_rhs: drop the chain
#     -rule "- 1" in `second` (i.e. `… eη*E/v_mE) * vp` without the `- 1`).
#     Verified bite: CMS.5 RED (difference no longer equals +vp — it becomes
#     0), and CMS.4's PVI-η pin RED (the package RHS shifts by -vp ≈ +0.037).
#     Restored; `git diff src` empty; suite GREEN.
#
#   Mutation M2 — src/CoordTransforms.jl pIII_transformed_rhs: halve the
#     kinetic term, `wp^2 / w` → `wp^2 / (2*w)`.  Verified bite: CMS.1 RED
#     (all 4 — constant-solution residual blows up), CMS.2 RED (rhs·x³ no
#     longer → ±4), CMS.4 PIII pin RED.  Restored.
#     (Note: the alternative `/4 → /2` mutant on the parameter block does NOT
#     bite CMS.1, because for params (0,0,1,-1) that block is identically 0
#     — γw³ and δe^{2ζ}/w cancel on the constant solution; it IS caught by
#     CMS.2 + CMS.4.  The kinetic-term mutant is the one that exercises CMS.1.)
#
# The cancellation test (CMS.6) is not src-mutated (it tests an arithmetic
# property, not a package function); its oracle is mutation-checked by the
# > 1e3 Float64-casualty assertion catching any precision regression.
# ============================================================================
