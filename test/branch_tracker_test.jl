# test/branch_tracker_test.jl — A4 of FFW 2017 reproduction, ADR-0013.
#
# Unit tests for the BranchTracker primitives (cut-crossing predicate,
# any-cut-crossed scan, per-branch sheet update, kwarg resolution).
# Walker-level integration is covered separately in
# test/path_network_branch_test.jl.
#
# Ground truth: ADR-0013 §"Decision" + FFW2017...md:178 verbatim.

using Test
using PadeTaylor: PadeTaylor    # accessor for the internal BranchTracker module
const BT = PadeTaylor.BranchTracker

@testset "BranchTracker (A4 of FFW 2017 / ADR-0013)" begin

    # -----------------------------------------------------------------
    # BT.1.1 — segment_crosses_cut: canonical geometric cases on the
    # standard cut (branch at origin, cut along arg = π = negative real
    # axis). Same shape as the standard `log(z)` branch cut.
    # -----------------------------------------------------------------
    @testset "BT.1.1: segment_crosses_cut on the standard log cut" begin
        b  = 0.0 + 0.0im
        α  = π                                # cut along arg = π
        # Segment crossing the negative real axis from above to below:
        # (-1, +0.5) → (-1, -0.5).  Goes from arg ≈ π - small to π + small
        # — i.e., crosses the cut.
        @test BT.segment_crosses_cut(-1.0 + 0.5im, -1.0 - 0.5im, b, α)
        # Same segment reversed: still crosses (predicate is symmetric).
        @test BT.segment_crosses_cut(-1.0 - 0.5im, -1.0 + 0.5im, b, α)
        # Segment on the POSITIVE real axis side: does not cross the cut.
        @test !BT.segment_crosses_cut(1.0 + 0.5im, 1.0 - 0.5im, b, α)
        # Segment entirely above the real axis: does not cross.
        @test !BT.segment_crosses_cut(-1.0 + 0.5im, 1.0 + 0.5im, b, α)
        # Segment touching the cut endpoint only (lands exactly on the
        # branch point): NOT counted (s = 0 excluded per ADR-0013).
        @test !BT.segment_crosses_cut(0.0 + 1.0im, 0.0 - 1.0im, b, α)
        # Segment FAR from origin crossing the cut: still crosses.
        @test BT.segment_crosses_cut(-100.0 + 0.5im, -100.0 - 0.5im, b, α)
        # Segment that would cross the cut's BACKWARD extension (positive
        # real axis): does not cross (the cut is one-sided, s ≥ 0 only).
        @test !BT.segment_crosses_cut(5.0 + 0.5im, 5.0 - 0.5im, b, α)
    end

    # -----------------------------------------------------------------
    # BT.1.2 — segment_crosses_cut: parallel and endpoint-touching cases.
    # -----------------------------------------------------------------
    @testset "BT.1.2: degenerate parallel / endpoint segments" begin
        b = 0.0 + 0.0im
        α = π
        # Segment EXACTLY along the negative real axis (parallel to cut):
        # det = 0 → predicate returns false.  Documented behaviour.
        @test !BT.segment_crosses_cut(-3.0 + 0.0im, -1.0 + 0.0im, b, α)
        # Step that lands EXACTLY on the cut at t = 1.0 (endpoint).
        # `0 < t < 1` excludes t = 1, so this returns false.
        @test !BT.segment_crosses_cut(-1.0 + 0.5im, -1.0 + 0.0im, b, α)
        # Step that STARTS exactly on the cut at t = 0.0 (start point).
        @test !BT.segment_crosses_cut(-1.0 + 0.0im, -1.0 - 0.5im, b, α)
    end

    # -----------------------------------------------------------------
    # BT.1.3 — Non-default cut angle: branch at b ≠ 0 with cut along
    # arg = 0 (positive real direction from b).  This is the configuration
    # users adopt when steering the walker AWAY from the negative-real
    # direction (e.g., for a PIII branch on the imaginary axis).
    # -----------------------------------------------------------------
    @testset "BT.1.3: branch off origin, cut along arg = 0" begin
        b = 2.0 + 0.0im                       # branch at z = 2
        α = 0.0                               # cut along positive real
        # Segment crossing the positive real axis to the RIGHT of b:
        # crosses the cut.
        @test BT.segment_crosses_cut(5.0 + 0.5im, 5.0 - 0.5im, b, α)
        # Segment crossing to the LEFT of b (between origin and b):
        # does NOT cross — only the half-line s ≥ 0 from b counts.
        @test !BT.segment_crosses_cut(1.0 + 0.5im, 1.0 - 0.5im, b, α)
        # Segment crossing the cut just to the right of b:
        @test BT.segment_crosses_cut(2.5 + 0.3im, 2.5 - 0.3im, b, α)
    end

    # -----------------------------------------------------------------
    # BT.1.4 — any_cut_crossed: multi-branch scan
    # -----------------------------------------------------------------
    @testset "BT.1.4: any_cut_crossed over multiple branches" begin
        branches   = (0.0 + 0.0im, 2π * im, -2π * im)
        cut_angles = (π, π, π)                # all standard log cuts
        # Step that crosses ONLY the upper branch's cut:
        z_old = -1.0 + (2π - 0.3) * im
        z_new = -1.0 + (2π + 0.3) * im
        @test BT.any_cut_crossed(z_old, z_new, branches, cut_angles)
        # Step that crosses NO cut (well above all of them):
        @test !BT.any_cut_crossed(1.0 + 100.0im, 2.0 + 100.0im,
                                   branches, cut_angles)
        # Step that crosses the origin branch's cut:
        @test BT.any_cut_crossed(-3.0 + 0.5im, -3.0 - 0.5im,
                                  branches, cut_angles)
        # Step that crosses the LOWER branch's cut:
        @test BT.any_cut_crossed(-1.0 + (-2π + 0.3) * im,
                                  -1.0 + (-2π - 0.3) * im,
                                  branches, cut_angles)
    end

    # -----------------------------------------------------------------
    # BT.1.5 — step_sheet_update: sign of update follows winding sign
    # -----------------------------------------------------------------
    @testset "BT.1.5: step_sheet_update sign follows winding direction" begin
        branches   = (0.0 + 0.0im,)
        cut_angles = (π,)
        sheet_in   = [3]
        # CCW crossing (going from below to above the negative-real
        # axis means CW around the origin, i.e., Δθ < 0).  Compute
        # explicitly: arg(-1 - 0.5i - 0) ≈ -2.677, arg(-1 + 0.5i - 0)
        # ≈ +2.677, normalised Δθ = +2.677 - (-2.677) = +5.355 →
        # wrap by -2π gives ≈ -0.928 < 0 → CW, bumps -1.
        sheet_out = BT.step_sheet_update(sheet_in,
                                          -1.0 - 0.5im, -1.0 + 0.5im,
                                          branches, cut_angles)
        @test sheet_out == [2]
        # Reverse direction: +1.
        sheet_back = BT.step_sheet_update(sheet_in,
                                           -1.0 + 0.5im, -1.0 - 0.5im,
                                           branches, cut_angles)
        @test sheet_back == [4]
        # Non-crossing step: pass-through unchanged.
        sheet_noop = BT.step_sheet_update(sheet_in,
                                           1.0 + 0.5im, 1.0 - 0.5im,
                                           branches, cut_angles)
        @test sheet_noop == [3]
        @test sheet_noop !== sheet_in     # fresh allocation, not aliased
    end

    # -----------------------------------------------------------------
    # BT.1.6 — step_sheet_update: only the crossed branch updates
    # -----------------------------------------------------------------
    @testset "BT.1.6: per-branch isolation of sheet updates" begin
        branches   = (0.0 + 0.0im, 2π * im)
        cut_angles = (π, π)
        sheet_in   = [0, 0]
        # Cross only the upper branch's cut:
        sheet_out = BT.step_sheet_update(sheet_in,
                                          -1.0 + (2π - 0.3)im,
                                          -1.0 + (2π + 0.3)im,
                                          branches, cut_angles)
        # Sheet for the upper branch must change; sheet for origin must
        # stay at 0 (the step does not cross the origin's cut).
        @test sheet_out[1] == 0
        @test sheet_out[2] != 0
    end

    # -----------------------------------------------------------------
    # BT.1.7 — resolve_cut_angles: scalar broadcasts, tuple validated
    # -----------------------------------------------------------------
    @testset "BT.1.7: resolve_cut_angles broadcasting + validation" begin
        branches = (0im, 2π*im, -2π*im)
        # Scalar broadcasts to all branches as Float64.
        out_scalar = BT.resolve_cut_angles(branches, π)
        @test out_scalar isa NTuple{3, Float64}
        @test all(==(Float64(π)), out_scalar)
        # Matching tuple passes through (with Float64 coercion).
        out_tuple = BT.resolve_cut_angles(branches, (0.0, π, π/2))
        @test out_tuple == (0.0, Float64(π), Float64(π/2))
        # Mismatched length throws ArgumentError per CLAUDE.md Rule 1.
        @test_throws ArgumentError BT.resolve_cut_angles(branches, (π, π))
        @test_throws ArgumentError BT.resolve_cut_angles(branches,
                                                          (π, π, π, π))
        # Empty branches: any scalar resolves to ().
        @test BT.resolve_cut_angles((), π) === ()
        @test BT.resolve_cut_angles((), ()) === ()
    end

end # @testset BranchTracker

# Mutation-proof procedure (verified before commit, 2026-05-16,
# at worklog 042; 99 GREEN total under unmutated impl):
#
#   Mutation B1  --  in `segment_crosses_cut`, drop the `s > 0` check
#     (return `0 < t < 1`).  Verified bite: 5 RED of 30 BT assertions
#     (BT.1.1 (3): the "backward extension" + endpoint-only assertions
#     fire false-positive crossings; BT.1.3 (1): "LEFT of b" segment
#     fires false-positive; BT.1.5 (1): the non-crossing pass-through
#     assertion bumps the sheet incorrectly when a near-branch step
#     looks like it intersects the FULL line through `b`).
#
#   Mutation B2  --  in `step_sheet_update`, flip the sign of the bump
#     (`-1` for `Δθ > 0`, `+1` for `Δθ < 0`).  Verified bite: 9 RED
#     across BT.1.5 (2) + the integration tests' PNB.1.3 (7) — the
#     direct sheet-direction assertions read `+1` instead of `-1` (and
#     vice versa) for the canonical CW/CCW crossings.
#
#   Mutation B3  --  in `step_sheet_update`, drop the
#     `segment_crosses_cut` gate (always apply the sign bump
#     regardless of whether the step crossed the cut).  Verified
#     bite: 52 RED, cascading widely (BT.1.5 + BT.1.6 + PNB.1.3) —
#     every step that does not actually cross the cut still bumps
#     the sheet counter, so the visited_sheet chain diverges from
#     the chain-walk's expected value at almost every node.
