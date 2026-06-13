# test/corpus_winding_test.jl
#
# ============================================================================
# Ground-truth corpus: WINDING NUMBER + CUT GEOMETRY (B11)
# (bead padetaylor-5qqr, epic padetaylor-p1v0).  The closed-form half of the
# multi-sheet tier: winding_delta / accumulate_winding / sheet_index
# (src/SheetTracker.jl) and segment_crosses_cut / any_cut_crossed /
# step_sheet_update (src/BranchTracker.jl).  Every oracle here is HAND-
# COMPUTED exact (winding angles, Cramer line/ray intersections) and was
# INDEPENDENTLY re-verified (Python cmath + sympy) before pinning — the
# catalogue is not infallible (it carried a wrong winding SIGN for the two-
# cuts case; see CWD.4 + the ERRATA note below).
#
# WHAT THIS FILE PINS beyond the committed sheet_tracker_test.jl /
# branch_tracker_test.jl (which exercise the same functions at OTHER points):
#   * the four-arc / off-centre CCW circumambulations end-to-end (deltas +
#     accumulation + sheet_index), including the -3π/2 → +π/2 normalisation;
#   * the oblique-cut Cramer crossing with the EXACT t,s hand-values;
#   * one step crossing TWO cuts simultaneously (both counters bump);
#   * GAP #4 (padetaylor-61um, RESOLVED 2026-06-13) — winding_delta on a near-
#     antipodal single step is CORRECT: a straight chord cannot subtend |Δθ| ≥ π
#     about an exterior branch, so the principal value IS the chord's true
#     winding; the "+1.1π" the earlier @test_broken expected is a CURVED-path
#     major arc, unrecoverable from two endpoints.  The genuine Rule-1 defect
#     was the unguarded ill-conditioned (branch-GRAZING) step in the walker's
#     sheet bookkeeping — now fixed by the step_sheet_update grazing guard
#     (CBr.3 / CPN.7.3).  CWD.5 pins winding_delta's correctness.
#
# GROUND TRUTH (Law 1): winding normalisation (-π,π] is principal-value arg;
# the Cramer 2×2 system is src/BranchTracker.jl:112-121; sheet convention is
# FFW2017...md:187-189.  Catalogue records: sheet-winding-delta-four-arcs,
# sheet-ccw-winding-number, sheet-segment-crosses-cut-oblique-yes/-no,
# sheet-two-cuts-simultaneously, sheet-winding-delta-precondition-violation.
#
# ERRATA (independent verify per Law 1 / Rule 2): the catalogue's
# sheet-two-cuts-simultaneously verify_note claims each counter "increments
# by +1".  WRONG: the segment runs rightward ABOVE both branch points, so
# winding_delta is NEGATIVE (clockwise) at both — the counters DECREMENT to
# -1.  Confirmed by cmath (angle drops 1.816→0.927 rad at branch 0) and by
# the package.  This file pins the code-correct & independently-verified
# [-1,-1] and asserts it is NOT [+1,+1].
#
# TOLERANCES: angles are exact rationals × π; atol 1e-14 (≈ a few ulp of π).
#
# MUTATION-PROOF (Rule 4): footer; mutant aimed at the winding_delta (-π,π]
# wrap confirmed RED on CWD.1, then restored.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor.SheetTracker: winding_delta, accumulate_winding, sheet_index
import PadeTaylor.BranchTracker:
    segment_crosses_cut, any_cut_crossed, step_sheet_update

@testset "Corpus: winding + cut geometry (CWD)" begin

    # -----------------------------------------------------------------------
    # CWD.1  Four 90° CCW arcs at unit radius, branch at 0.  Each delta=π/2;
    #        step 3 raw = -3π/2 normalises to +π/2; total = 2π ⇒ sheet +1.
    #        catalogue sheet-winding-delta-four-arcs.
    # -----------------------------------------------------------------------
    @testset "CWD.1: four-arc CCW unit circle ⇒ each Δ=π/2, total 2π, sheet 1" begin
        path = ComplexF64[1, im, -1, -im, 1]
        deltas = [winding_delta(path[i-1], path[i], 0.0im) for i in 2:5]
        for d in deltas
            @test isapprox(d, π / 2; atol = 1e-14)
        end
        # Step 3 specifically: raw -3π/2 → +π/2 (the load-bearing wrap).
        @test isapprox(winding_delta(-1.0 + 0im, -1.0im, 0.0im), π / 2;
                       atol = 1e-14)
        acc = accumulate_winding(path, 0.0im)
        @test acc[1] == 0.0
        @test isapprox(acc, [0.0, π/2, π, 3π/2, 2π]; atol = 1e-13)
        @test sheet_index(acc[end]) == 1
    end

    # -----------------------------------------------------------------------
    # CWD.2  Full CCW circumambulation of an OFF-CENTRE branch at 2+0j:
    #        diamond path, total 2π ⇒ sheet +1.  catalogue
    #        sheet-ccw-winding-number.  Guards the (z-branch) recentre.
    # -----------------------------------------------------------------------
    @testset "CWD.2: off-centre CCW diamond around 2+0j ⇒ total 2π, sheet 1" begin
        path = ComplexF64[4, 2 + 2im, 0, 2 - 2im, 4]
        deltas = [winding_delta(path[i-1], path[i], 2.0 + 0im) for i in 2:5]
        for d in deltas
            @test isapprox(d, π / 2; atol = 1e-14)
        end
        acc = accumulate_winding(path, 2.0 + 0im)
        @test isapprox(acc[end], 2π; atol = 1e-13)
        @test sheet_index(acc[end]) == 1
        # Anti-regression: with branch at ORIGIN (not 2), this path does NOT
        # enclose ⇒ winding ≈ 0 ⇒ sheet 0.  Confirms the recentre matters.
        @test sheet_index(accumulate_winding(path, 0.0im)[end]) == 0
    end

    # -----------------------------------------------------------------------
    # CWD.3  Oblique cut at π/4 from branch 1+1j; Cramer line/ray crossing.
    #        YES: z=0+2j→3+0j  ⇒ t=2/5, s=√2/5  ⇒ crosses.
    #        NO : z=0+0j→2+0j  ⇒ t=0 (endpoint, excluded), s=-√2  ⇒ no cross.
    #        catalogue sheet-segment-crosses-cut-oblique-yes / -no.
    # -----------------------------------------------------------------------
    @testset "CWD.3: oblique-cut Cramer crossing (yes / no), exact t,s" begin
        # YES geometry.
        @test segment_crosses_cut(0 + 2im, 3 + 0im, 1 + 1im, π / 4) == true
        # Re-derive t,s by the same Cramer formula and pin them exactly.
        let d = (3 + 0im) - (0 + 2im), δ = (1 + 1im) - (0 + 2im), u = cis(π / 4)
            det = imag(d * conj(u))
            t = imag(δ * conj(u)) / det
            s = imag(δ * conj(d)) / det
            @test isapprox(t, 2 / 5; atol = 1e-14)
            @test isapprox(s, sqrt(2) / 5; atol = 1e-14)
            @test 0 < t < 1 && s > 0
        end
        # NO geometry: endpoint t=0 excluded, and s<0.
        @test segment_crosses_cut(0 + 0im, 2 + 0im, 1 + 1im, π / 4) == false
        let d = (2 + 0im) - (0 + 0im), δ = (1 + 1im) - (0 + 0im), u = cis(π / 4)
            det = imag(d * conj(u))
            t = imag(δ * conj(u)) / det
            s = imag(δ * conj(d)) / det
            @test isapprox(t, 0.0; atol = 1e-14)       # on the start endpoint
            @test s < 0
        end
    end

    # -----------------------------------------------------------------------
    # CWD.4  ONE step crosses TWO cuts at once: horizontal segment at height 2
    #        over upward cuts from 0 and 1.  Both Cramer solutions valid
    #        (t1=1/4, t2=3/4, s1=s2=2).  catalogue sheet-two-cuts-simultaneously.
    #        ERRATA: the catalogue says counters increment +1; the step is
    #        CLOCKWISE (rightward above the branches) ⇒ counters go to -1.
    # -----------------------------------------------------------------------
    @testset "CWD.4: one step crosses two cuts ⇒ both counters bump (−1)" begin
        z_cur, z_new = -0.5 + 2im, 1.5 + 2im
        branches  = (0.0 + 0im, 1.0 + 0im)
        cut_angles = (π / 2, π / 2)
        @test segment_crosses_cut(z_cur, z_new, branches[1], cut_angles[1])
        @test segment_crosses_cut(z_cur, z_new, branches[2], cut_angles[2])
        @test any_cut_crossed(z_cur, z_new, branches, cut_angles) == true
        # Independently: both winding deltas are NEGATIVE (clockwise, above).
        @test winding_delta(z_cur, z_new, branches[1]) < 0
        @test winding_delta(z_cur, z_new, branches[2]) < 0
        updated = step_sheet_update([0, 0], z_cur, z_new, branches, cut_angles)
        @test updated == [-1, -1]                      # both bumped, CW = −1
        @test updated != [1, 1]                        # catalogue's wrong sign
    end

    # -----------------------------------------------------------------------
    # CWD.5  GAP #4 — single-step precondition VIOLATION |Δθ|≥π.
    #        branch=0; z_old=0.1, z_new=0.1·e^{i·1.1π}.  TRUE subtended angle
    #        = +1.1π (CCW), but principal arg(z_new)=−0.9π ⇒ raw Δθ=−0.9π,
    #        already inside (-π,π] ⇒ winding_delta returns −0.9π.  The code
    #        SILENTLY loses a full +2π revolution (returns −0.9π for a true
    #        +1.1π) with NO guard.
    #
    #        FINDING (Rule 1): the contract is documented (SheetTracker.jl
    #        :275-277) but UNENFORCED.  We pin the documented-wrong value AND
    #        @test_broken the Rule-1-correct behaviour so it is loud and
    #        trackable.  REPORT: file a Rule-1 bead — winding_delta should
    #        throw (or the caller subdivide) when |raw Δθ| ≥ π.
    # -----------------------------------------------------------------------
    @testset "CWD.5: GAP #4 — single straight step, winding_delta is CORRECT (61um)" begin
        z_old = 0.1 + 0.0im
        z_new = 0.1 * exp(im * 1.1 * π)
        actual = winding_delta(z_old, z_new, 0.0im)
        true_subtended = 1.1 * π                        # +1.1π CCW (geometry)

        # Pin the ACTUAL (documented-wrong) behaviour: returns −0.9π.
        @test isapprox(actual, -0.9 * π; atol = 1e-12)
        # The error is exactly one missed revolution.
        @test isapprox(actual - true_subtended, -2π; atol = 1e-12)

        # WHY winding_delta is CORRECT here (61um resolution).  The raw
        # arg-difference is angle(z_new) − angle(z_old) = −0.9π − 0 = −0.9π,
        # already in (−π, π], and that IS the true winding of the STRAIGHT
        # segment [z_old, z_new] about the origin (a straight chord can never
        # subtend |Δθ| ≥ π about an exterior point).  The +1.1π "true subtended"
        # is the MAJOR-arc angle of a CURVED path — and which way a path winds
        # between two endpoints is information winding_delta's two arguments do
        # NOT carry, so +1.1π is unrecoverable from this call (worklog
        # 074:47-74).  winding_delta is a pure, correct principal-value primitive.
        raw_dtheta = angle(z_new) - angle(z_old)
        @test isapprox(raw_dtheta, -0.9 * π; atol = 1e-12)   # < π in magnitude
        @test abs(raw_dtheta) < π
        @test isapprox(actual, raw_dtheta; atol = 1e-14)     # winding_delta == raw == truth

        # ----- 61um RESOLUTION (was @test_broken expecting +1.1π) ------------
        # The prior marker asserted winding_delta SHOULD return the unwrapped
        # +1.1π.  That is geometrically impossible for a straight chord and
        # information-theoretically unrecoverable from two endpoints; it is NOT
        # a bug in winding_delta.  The genuine Rule-1 defect was the UNGUARDED
        # ill-conditioned step in the WALKER's sheet bookkeeping, now fixed by
        # the step_sheet_update grazing guard (demonstrated on a crossed-cut
        # grazing geometry in CBr.3 and end-to-end in CPN.7.3).  This step's
        # straight chord does NOT cross the arg=π cut from the origin (it leaves
        # the +real axis directly into the lower-left), so it is a pure
        # winding_delta primitive fixture; reframed per bd decision 2026-06-13.
        @test sheet_index(actual) == 0                       # |−0.9π| < π ⇒ sheet 0
    end

end # @testset Corpus winding (CWD)

# ============================================================================
# MUTATION-PROOF (Rule 4) — verified 2026-06-05, then restored byte-clean.
#
#   Mutation W1 — src/SheetTracker.jl winding_delta: drop the (-π,π] wrap
#     (delete the `if Δθ ≤ -π … elseif Δθ > π …` block, return raw Δθ).
#     Verified bite: CWD.1 RED — step 3's raw −3π/2 is no longer normalised
#     to +π/2, so the four-arc deltas / accumulation / sheet_index all break;
#     CWD.2 RED for the off-centre diamond's step-3 wrap likewise.  Restored;
#     `git diff src` empty; suite GREEN.
#
#   Mutation W2 — src/BranchTracker.jl step_sheet_update: flip the sign rule
#     to `Δθ > 0 ? -1 : 1`.  Verified bite: CWD.4 RED (returns [1,1] instead
#     of [-1,-1]).  Restored.
#
#   GAP #4 (CWD.5) was reframed when padetaylor-61um was RESOLVED (2026-06-13):
#     winding_delta is mathematically correct for the straight chord (it never
#     loses a revolution — a chord cannot subtend |Δθ| ≥ π about an exterior
#     branch), so the former @test_broken (expecting an impossible unwrapped
#     +1.1π) was removed and replaced by live pins of the correct value.  The
#     Rule-1 enforcement landed in BranchTracker.step_sheet_update's grazing
#     guard (the WALKER's sheet-bump site), mutation-proven in
#     branch_tracker_test.jl / corpus_pathnet_winding_test.jl.  See ADR-0031.
# ============================================================================
