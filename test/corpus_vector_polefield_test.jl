# test/corpus_vector_polefield_test.jl
#
# ============================================================================
# Ground-truth corpus: VECTOR pole-field extraction — the ADR-0028
# Froissart filter (gap #5) + the DEFERRED CM imaginary-collision case
# (bead padetaylor-jmcf, epic padetaylor-p1v0).
#
# This is the companion file to corpus_vector_test.jl, split out for the
# CLAUDE.md Rule 6 ≤200-code-LOC cap.  It pins:
#
#   CVE.5  extract_poles_shared_q's Froissart `min_residue` filter — the
#          consumer fix ADR-0028 Amendment 3 added, NEVER directly tested
#          on the extractor before (only shared_q_residue in isolation,
#          SPD.6).  A hand-built shared-Q node carries ONE genuine shared
#          pole AND one planted near-cancelling Froissart doublet; the
#          extractor must return the genuine pole and DROP the doublet.
#
#   CVE.6  cm-n2-imaginary-collision — DEFERRED.  The catalogue's
#          "headline complex-pole-bridge for an integrable system"
#          (md:568-583, confidence=MEDIUM, "not yet independently
#          verified").  Independent verification (sympy + 50-dps mpmath
#          RK4) shows the catalogue closed form is WRONG for the ODE the
#          repo integrates — this testset documents the discrepancy and
#          marks the case deferred, rather than pinning a shaky value
#          (Rule 9 — honest about limits; the brief's explicit
#          "DO NOT pin a shaky value" instruction).
#
# GROUND TRUTH (Law 1 — independently verified)
# ---------------------------------------------
#   * CVE.5: the residue separation is the ADR-0028 calibration
#     (src/SharedPadeDefect.jl:134-154; probe
#     external/probes/adr0028-froissart-consumer/): genuine shared poles
#     have residue ≥ 2.1, cell-B doublets ≤ 5.5e-8, an ~8-order gap, so
#     `min_residue = 1e-4` cleanly separates them.  The planted doublet
#     here has residue 5.0e-8 (matched to the documented cell-B floor);
#     the genuine pole ~1.0.
#   * CVE.6: external/probes/corpus-oracles/vector/verify_cm.py — sympy
#     gives x1'' - 4/(x1-x2)³ = -2√2·i·r0³/(2r0⁴-t²)^(3/2) ≠ 0 for the
#     catalogue form (i/2)√(4-2t²); the 50-dps RK4 from the imaginary IC
#     DIVERGES from it (|err| 0.045 at t=0.3, 0.78 at t=1.2).  The
#     catalogue form solves the ATTRACTIVE system x1''=-4/(x1-x2)³
#     (sympy residual 0 there), but the repo integrates the REPULSIVE
#     system x1''=+4/(x1-x2)³ (calogero_moser_test.jl cm_rhs, the AMM
#     rational-KdV convention).  The CORRECT repulsive imaginary-IC form
#     is (i/2)√(4+2t²) — which has NO real collision (the particles
#     repel; |x1| grows monotonically), so there is no pole to bridge in
#     the repo's ODE.  Reported loudly in the final summary + ERRATA.
#
# DE-DUPLICATION vs the existing suite (CLAUDE.md Rule 3)
# ------------------------------------------------------
#   * shared_pade_dispatch_test.jl SPD.6 pins shared_q_residue in
#     ISOLATION on bare coefficient vectors.  CVE.5 pins the FILTER
#     end-to-end inside extract_poles_shared_q on a stored
#     VectorPathNetworkSolution — the consumer ADR-0028 A3.2/A3.3 fixed,
#     which SPD.6 does not exercise.
#   * vector_path_network_test.jl VPN.1.3 hand-builds shared-Q nodes and
#     tests the cross-node `min_support` filter + z-plane mapping, but
#     every numerator is the constant [1.0] (no near-cancellation), so it
#     never exercises the `min_residue` Froissart filter.  CVE.5 is the
#     missing direct test of that filter.
#
# MUTATION-PROOF (Rule 4): footer block.  The SHARED-Q-targeted mutant
# perturbs the `min_residue` threshold and confirms the doublet leaks.
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor.VectorPoleField: extract_poles_shared_q
using PadeTaylor.VectorPathNetwork: VectorPathNetworkSolution
using PadeTaylor.SharedPadeDefect: shared_q_residue

@testset "Corpus: vector pole-field Froissart filter (CVE.5/6)" begin

    # ----------------------------------------------------------------------
    # CVE.5  extract_poles_shared_q filters a planted Froissart doublet.
    #
    # One node, shared denominator Q(t) = (1-t)(1-2t) = 1 - 3t + 2t² with
    # TWO roots: t=1 (genuine) and t=0.5 (the planted doublet).  Each
    # numerator P_c(t) NEARLY vanishes at t=0.5 (a near-cancelling
    # pole/zero pair) but is O(1) at t=1.  With z_node=0, h=1 the roots
    # map to z=t, so the genuine pole sits at z=1 and the doublet at
    # z=0.5.  The Froissart `min_residue` filter must keep z=1, drop z=0.5.
    # ----------------------------------------------------------------------
    @testset "CVE.5  Froissart doublet dropped, genuine shared pole kept" begin
        Q = ComplexF64[1.0, -3.0, 2.0]                 # roots t=1, t=0.5
        # P_c(t) = (1 - 2t) + ε_c·t : zero at t=0.5 up to ε_c/2 (doublet),
        # ≈ -1 at t=1 (genuine).  ε_1=1e-7, ε_2=3e-8 ⇒ residue at t=0.5
        # ≈ 5e-8 (the documented cell-B doublet floor); ≈ 1.0 at t=1.
        P1 = ComplexF64[1.0, -2.0 + 1e-7]
        P2 = ComplexF64[0.5, -1.0 + 3e-8]
        nums = [P1, P2]

        # Sanity: the residue separation matches the ADR-0028 calibration.
        r_genuine = shared_q_residue(nums, Q, ComplexF64(1.0))
        r_doublet = shared_q_residue(nums, Q, ComplexF64(0.5))
        @test r_genuine ≥ 1e-4                          # genuine survives
        @test r_doublet < 1e-4                          # doublet below thr
        @test isapprox(r_genuine, 1.0; atol = 1e-5)     # ~1.0 (header)
        @test r_doublet < 1e-7                          # ~5e-8 (header)

        z_node = 0.0 + 0.0im
        h_node = 1.0 + 0.0im                            # t* maps to z = t*
        sol = VectorPathNetworkSolution{Float64}(
            [z_node], [ComplexF64[1.0, 1.0]], [h_node], [nums], [Q], [0])

        # Default min_residue = 1e-4: ONLY the genuine pole at z=1 returns.
        poles = extract_poles_shared_q(sol; radius_t = 5.0, min_support = 1)
        @test length(poles) == 1
        @test isapprox(poles[1], 1.0 + 0.0im; atol = 1e-10)
        @test all(abs(p - (0.5 + 0im)) > 0.1 for p in poles)  # doublet absent

        # Lower min_residue below the doublet's residue ⇒ the doublet
        # LEAKS back in (z=0.5 reappears) — proof the FILTER, not the
        # geometry, is what dropped it.  This is the CVE.5 in-test
        # mutation-proof of the shared-Q Froissart machinery.
        leaked = extract_poles_shared_q(sol; radius_t = 5.0,
                                        min_support = 1, min_residue = 1e-10)
        @test length(leaked) == 2
        @test any(abs(p - (0.5 + 0im)) < 1e-9 for p in leaked)
        @test any(abs(p - (1.0 + 0im)) < 1e-9 for p in leaked)
    end

    # ----------------------------------------------------------------------
    # CVE.6  cm-n2-imaginary-collision — DEFERRED (catalogue form is wrong).
    #
    # The catalogue (md:574) pins x1(t)=(i/2)√(4-2t²) with a collision at
    # t*=√2 for the "Same ODE as cm-n2-rational-kdv" (repulsive, EOM +4).
    # Independent verification proves that form does NOT solve the
    # repulsive ODE.  We DO NOT pin it.  Instead we assert the
    # discrepancy here (so the deferral is itself tested, not just a
    # comment), and the case is recorded in ERRATA.md.
    # ----------------------------------------------------------------------
    @testset "CVE.6  CM imaginary-collision DEFERRED (catalogue form wrong)" begin
        # The repulsive CM ODE the repo integrates (calogero_moser_test.jl).
        x1cat(t) = (im / 2) * sqrt(complex(4 - 2 * t^2))   # catalogue form
        # Its 2nd derivative must equal 4/(x1-x2)³ with x2=-x1 ⇒ 4/(2x1)³
        # if the catalogue form solved the REPULSIVE system.  It does NOT:
        # at t=0.6 the residual x1''(t) - 4/(2x1)³ is O(1), not ~0.
        # (Finite-difference x1'' to avoid trusting any analytic claim.)
        t  = 0.6
        dt = 1e-4
        x1pp = (x1cat(t + dt) - 2 * x1cat(t) + x1cat(t - dt)) / dt^2
        rhs  = 4 / (2 * x1cat(t))^3
        resid = x1pp - rhs
        # The catalogue form fails the REPULSIVE ODE by an O(1) margin —
        # this is the defect that defers the case.
        @test abs(resid) > 1e-2
        # ...but it DOES solve the ATTRACTIVE ODE x1'' = -4/(2x1)³
        # (residual ~0) — confirming the catalogue picked the wrong sign
        # convention (attractive trajectory under a repulsive ODE label).
        @test abs(x1pp - (-rhs)) < 1e-3
        # The CORRECT repulsive imaginary-IC form is (i/2)√(4+2t²): it
        # solves x1''=+4/(2x1)³ and has NO real collision (4+2t²>0 ∀ real
        # t).  Pin the correctly-solving form so the deferral records the
        # right answer, not just the wrong one.
        x1ok(t) = (im / 2) * sqrt(complex(4 + 2 * t^2))
        x1ok_pp = (x1ok(t + dt) - 2 * x1ok(t) + x1ok(t - dt)) / dt^2
        rhs_ok  = 4 / (2 * x1ok(t))^3
        @test abs(x1ok_pp - rhs_ok) < 1e-3              # solves repulsive ODE
        @test 4 + 2 * t^2 > 0                            # no real collision
        # DEFERRAL marker: the keystone pole-bridge is carried instead by
        # CVE.1/CVE.2/CVE.3 (corpus_vector_test.jl), which bridge genuine
        # shared poles with externally-verified closed forms.
        @test_broken false   # cm-n2-imaginary-collision deferred (see header)
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# src/ confirmed clean via `git diff --stat src/` after each restore.
#
# M-impl (CVE.5, the SHARED-Q-targeted mutant — the Froissart threshold):
#   in src/VectorPoleField.jl `extract_poles_shared_q`, change the
#   default kwarg `min_residue::Real = 1.0e-4` → `1.0e-10` (lower the
#   filter below the doublet's 5e-8 residue).  Result: CVE.5's
#   `length(poles) == 1` went RED — the planted doublet at z=0.5 leaked
#   into the default-call result (2 poles, not 1), exactly the
#   pole-root-consumer corruption ADR-0028 A3.2 documents.  The in-test
#   `min_residue = 1e-10` leak assertion stayed GREEN (it expects the
#   leak), confirming the threshold is the load-bearing knob.  Restored
#   src/VectorPoleField.jl byte-for-byte.
#
# M-oracle (CVE.5): change P1's near-cancel coefficient `-2.0 + 1e-7` →
#   `-2.0 + 1e-1` (a doublet that does NOT nearly cancel — residue jumps
#   to ~0.05, above the 1e-4 threshold).  Result: the `r_doublet < 1e-4`
#   and the `length(poles) == 1` assertions went RED (the now-genuine
#   "doublet" is kept, 2 poles).  Confirms CVE.5 genuinely distinguishes
#   a NEAR-cancelling doublet from a real pole.  Restored the literal.
#
# M-impl (CVE.6): the deferral is a documented-discrepancy test, not an
#   impl pin; no src/ mutation applies.  The discrepancy itself is
#   re-derived (sympy + mpmath RK4) in verify_cm.py — re-running it with
#   the catalogue form's sign restored to the attractive convention makes
#   the symbolic residual 0, the cross-check that localises the error.
#
# Run standalone with:
#   julia --project=. test/corpus_vector_polefield_test.jl
# ============================================================================
