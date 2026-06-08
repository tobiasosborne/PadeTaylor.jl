# test/corpus_pathnet_winding_test.jl
#
# ============================================================================
# Ground-truth corpus: PATH-NETWORK WINDING WALK (CPN.7) -- THE HEADLINE
# padetaylor-61um CATCH UNDER A REAL WALK.
# (bead padetaylor-nkw8, epic padetaylor-25og.)  Family I of the v2 corpus
# plan, docs/test_corpus/02_corpus_extension_plan.md:374-381.
#
# WHAT THIS FILE PINS.  This is the FIRST time padetaylor-61um is exercised by
# an actual `path_network_solve(cross_branch=true)` Stage-1 WEDGE STEP whose
# realised |dtheta| >= pi about a tracked branch point.  The walker-side sheet
# bookkeeping (BranchTracker.step_sheet_update, which bumps visited_sheet by
# sign(SheetTracker.winding_delta) per cut crossing) is driven through the full
# walk, not a bare unit call.
#
# HOW THIS DIFFERS FROM CWD.5 / CBr.3 (the existing 61um fixtures).  CWD.5
# (corpus_winding_test.jl) and CBr.3 (corpus_elementary_branch_test.jl) are
# BARE winding_delta(z_old, z_new, branch) unit fixtures -- they pin the wrap
# arithmetic in isolation, with hand-chosen endpoints.  CPN.7 is the REAL-WALK
# complement: the endpoints are produced by the package's OWN wedge stepper,
# the cut crossing is detected by the package's OWN segment_crosses_cut, and
# the corrupted sheet index lands in the package's OWN visited_sheet field.
# Same -2*pi lost revolution, now end-to-end through the solver.
#
# GROUND TRUTH (Law 1).  Background ODE: y = sqrt(1-z^2), C=1, recast for
# path_network_solve as f(z,y,yp) = (-1 - yp^2)/y (sympy residual==0 GATE in
# capture.py; branch points +-1 = zeros of 1-z^2; IC z0=0 => y0=1, y'0=0).
# The +-sqrt sheet oracle at z=1.5+0.5i (principal vs sheet-flipped, exact
# negation) is the SEMANTIC meaning of a sheet flip -- mpmath dps=50, NEVER
# the package's own output.  Winding normalisation (-pi,pi] is principal-value
# arg (SheetTracker.jl:283-292); sheet convention FFW2017...md:187-189.
# Transcendental pins live in test/_oracle_corpus_pathnet_winding.jl.
#
# THE TWO LEGS (verified geometry; re-measured at test time).  Branch=1,
# ring radius r=0.45; two single-wedge-step cross_branch=true walks differing
# ONLY in step coarseness:
#   COARSE -- IC arg +0.82*pi, target arg +1.18*pi, h=0.48.  The single step
#     straddles the arg=pi cut with realised subtended ~+0.555*pi (< pi, CCW).
#     winding_delta is EXACT; visited_sheet=[1] is CORRECT (matches TRUE +1).
#   FINE -- IC arg -0.47*pi, target arg +0.59*pi, h=0.70.  The single step
#     straddles the cut with realised TRUE subtended ~+1.083*pi (> pi, CCW).
#     winding_delta CLAMPS into (-pi,pi] -> ~ -2.88 rad (loses a +2*pi
#     revolution) -> visited_sheet=[-1], the 61um-CORRUPTED sheet (TRUE +1).
# Same ODE/branch/cut/CCW-direction; only the subtended angle crosses the pi
# threshold -- proving the corruption is STEP-SIZE-TRIGGERED.
#
# TOLERANCES (justified by METHOD).
#   * Closed-form +-sqrt oracle: exact analytic sqrt -> atol 1e-12.
#   * Realised wedge-step angles / winding_delta: geometry of the package's
#     own deterministic (rng_seed) step -> atol 1e-9 (path-tree FP noise).
#   * visited_sheet entries are EXACT integers (== comparison, no tolerance).
#
# REFERENCES.  Plan Family I (02_corpus_extension_plan.md:374-381); bug
# padetaylor-61um (winding_delta unenforced |dtheta|>=pi precondition);
# FFW2017 sheet tracking (md:163-189); DLMF 4.x (sqrt principal branch).
# Catalogue mirror of CWD.5/CBr.3 (the bare unit fixtures); this is the
# REAL-WALK complement.
#
# MUTATION-PROOF (Rule 4): footer -- M-oracle (flip a pinned +-sqrt value ->
# reddens) + M-impl (perturb winding_delta's (-pi,pi] wrap -> the coarse-ring
# sheet assertion reddens AND the FINE @test_broken marker flips to
# "unexpectedly passes", proving it is the LIVE 61um auto-flip).  ACTUALLY
# EXECUTED then restored byte-clean.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor.SheetTracker: winding_delta, accumulate_winding, sheet_index
import PadeTaylor.BranchTracker: segment_crosses_cut

include(joinpath(@__DIR__, "_oracle_corpus_pathnet_winding.jl"))

# CPN.7 background RHS: y = sqrt(1-z^2), recast f(z,y,yp) = (-1 - yp^2)/y.
const CPN7_F = (z, y, yp) -> (-1 - yp^2) / y
const CPN7_BRANCH = 1.0 + 0.0im
cpn7_ysqrt(z) = sqrt(1 - z^2)            # mpmath-agreeing principal +sqrt

# Drive a SINGLE-wedge-step cross_branch=true real walk from IC=P to target=Q.
# Returns the (deterministic, rng_seed) PathNetworkSolution; asserts exactly
# one step was taken (two visited nodes) so the leg is an unambiguous single
# crossing edge.
function cpn7_leg(P, Q, h)
    y0  = cpn7_ysqrt(P)
    yp0 = -P / y0
    prob = PadeTaylorProblem(CPN7_F, (y0, yp0), (P, Q); order = 20)
    sol = path_network_solve(prob, ComplexF64[Q]; h = h, order = 20,
                             rng_seed = 7, branch_points = (CPN7_BRANCH,),
                             cross_branch = true)
    @assert length(sol.visited_z) == 2 "CPN.7 leg expected 1 step (2 nodes), got $(length(sol.visited_z))"
    return sol
end

@testset "Corpus: path-network winding (CPN.7) (CPN)" begin

    # -----------------------------------------------------------------------
    # CPN.7.1  Background + the +-sqrt sheet oracle.  IC at z0=0 (y0=1, y'0=0);
    #          branch points +-1; the sheet flip is exactly the sign change.
    # -----------------------------------------------------------------------
    @testset "CPN.7.1: background recast + ±√ sheet oracle (closed form)" begin
        # IC values y0=1, y'0=0 (exact) reproduce the principal +sqrt at z0=0.
        @test cpn7_ysqrt(0.0 + 0.0im) == CPN7_Y0
        @test isapprox(CPN7_Y0, 1.0 + 0.0im; atol = 1e-12)
        @test CPN7_YP0 == 0.0 + 0.0im
        # Branch points +-1 are the zeros of 1 - z^2 (structural pin).
        @test (1 - 1.0^2) == 0.0 && (1 - (-1.0)^2) == 0.0

        # The +-sqrt sheet oracle at z=1.5+0.5i.  Principal +sqrt is the mpmath
        # dps=50 pin; the sheet-flipped value is exactly its negation (the
        # SEMANTIC meaning of sheet_index = +1 after one CCW loop about z=1).
        zt = 1.5 + 0.5im
        @test isapprox(cpn7_ysqrt(zt), CPN7_YPRIN_AT_1p5_0p5i; atol = 1e-12)
        @test isapprox(CPN7_YFLIP_AT_1p5_0p5i, -CPN7_YPRIN_AT_1p5_0p5i;
                       atol = 1e-12)
        @test isapprox(CPN7_YPRIN_AT_1p5_0p5i,
                       0.6335517491618166 - 1.1838022718621541im; atol = 1e-12)
    end

    # -----------------------------------------------------------------------
    # CPN.7.2  COARSE leg -- a REAL single-step cross_branch=true walk whose
    #          wedge step subtends < pi.  winding_delta is EXACT, so the
    #          walker-side sheet bookkeeping is CORRECT: visited_sheet=[1]
    #          agrees with the TRUE CCW crossing (+1).
    # -----------------------------------------------------------------------
    @testset "CPN.7.2: COARSE step (<π) ⇒ visited_sheet CORRECT" begin
        sol = cpn7_leg(CPN7_COARSE_P, CPN7_COARSE_Q, 0.48)
        zp, zc = sol.visited_z[1], sol.visited_z[2]

        # The realised single wedge step crosses the arg=π cut from branch=1.
        @test segment_crosses_cut(zp, zc, CPN7_BRANCH, π) == true
        # Its realised subtended angle is below π (no wrap) and CCW (+).
        wd = winding_delta(zp, zc, CPN7_BRANCH)
        @test wd > 0
        @test abs(wd) < π
        @test isapprox(abs(angle((zc - CPN7_BRANCH) / (zp - CPN7_BRANCH))),
                       abs(wd); atol = 1e-9)        # principal arg-diff == wd

        # Therefore the package's visited_sheet is the TRUE sheet: +1 (one CCW
        # crossing).  EXACT integer comparison -- the load-bearing pin.
        @test sol.visited_sheet[1] == [0]           # root on principal sheet
        @test sol.visited_sheet[2] == [1]           # CORRECT +1 bump
        # Independent cross-check: step_sheet_update on the realised edge bumps
        # +1, matching the package and the TRUE CCW direction.
        @test PadeTaylor.BranchTracker.step_sheet_update(
                  [0], zp, zc, (CPN7_BRANCH,), (Float64(π),)) == [1]
    end

    # -----------------------------------------------------------------------
    # CPN.7.3  FINE leg -- a REAL single-step cross_branch=true walk whose
    #          wedge step subtends > pi.  winding_delta WRAPS into (-π,π] and
    #          returns a NEGATIVE value, so step_sheet_update bumps -1 instead
    #          of +1: visited_sheet=[-1] is the padetaylor-61um CORRUPTION.
    #          The TRUE crossing is CCW (+1).
    # -----------------------------------------------------------------------
    @testset "CPN.7.3: FINE step (>π) ⇒ visited_sheet CORRUPTED (61um)" begin
        sol = cpn7_leg(CPN7_FINE_P, CPN7_FINE_Q, 0.70)
        zp, zc = sol.visited_z[1], sol.visited_z[2]

        # The realised single wedge step crosses the same arg=π cut.
        @test segment_crosses_cut(zp, zc, CPN7_BRANCH, π) == true

        # TRUE (continuous) subtended angle of the realised edge exceeds π and
        # is CCW (+): unwrap the principal arg-difference by +2π.
        raw  = angle((zc - CPN7_BRANCH) / (zp - CPN7_BRANCH))   # principal, < 0
        true_subtended = raw < 0 ? raw + 2π : raw
        @test true_subtended > π
        @test true_subtended < 2π                               # < a full loop

        # winding_delta returns the WRAPPED (negative) value -- the silent loss
        # of one +2π revolution (the 61um precondition violation).
        wd = winding_delta(zp, zc, CPN7_BRANCH)
        @test wd < 0
        @test isapprox(wd - true_subtended, -2π; atol = 1e-9)   # lost exactly -2π
        @test isapprox(wd, raw; atol = 1e-12)                   # = principal arg-diff

        # CURRENT (documented-wrong) sheet bookkeeping: bumps -1 (sign of the
        # wrapped wd), so visited_sheet=[-1].  GREEN today -- pins the bug.
        @test sol.visited_sheet[2] == [-1]

        # ----- 61um REAL-WALK AUTO-FLIP MARKER -------------------------------
        # The TRUE CCW crossing should bump +1, giving visited_sheet=[1].  Today
        # the wrapped winding_delta yields [-1], so this is BROKEN.  It auto-
        # flips to "unexpectedly passes" the day padetaylor-61um is fixed (a
        # winding_delta that returns the TRUE +1.083π makes step_sheet_update
        # bump +1).  DO NOT delete -- this is the LIVE real-walk 61um tracker,
        # the complement to CWD.5/CBr.3's unit-level markers.
        @test_broken sol.visited_sheet[2] == [1]
    end

end # @testset Corpus path-network winding (CPN.7)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) -- ACTUALLY EXECUTED, then restored exactly.
#
#   M-oracle -- flip a pinned +-sqrt value.  In
#     test/_oracle_corpus_pathnet_winding.jl negate the real part of
#     CPN7_YPRIN_AT_1p5_0p5i (0.6335... -> -0.6335...).  Verified bite: CPN.7.1
#     RED -- `isapprox(cpn7_ysqrt(zt), CPN7_YPRIN_...)` fails (the analytic
#     +sqrt no longer matches the corrupted pin) AND the explicit-literal pin
#     fails.  Restored the literal; file GREEN again.
#
#   M-impl -- perturb the (-π,π] wrap in src/SheetTracker.jl winding_delta
#     (delete the `if Δθ ≤ -π … elseif Δθ > π …` block; return the raw
#     Δθ = angle(z_new-branch) - angle(z_old-branch)).  Verified bite (9 RED +
#     1 Unexpected-Pass of 22):
#     (i)  CPN.7.3's @test_broken FLIPS to "Unexpected Pass" -- with the wrap
#          removed, winding_delta returns the raw +1.083π-ish POSITIVE value,
#          step_sheet_update bumps +1, visited_sheet[2]==[1], so the marker
#          passes.  AND CPN.7.3's `wd < 0` / `wd - true ≈ -2π` pins go RED (wd
#          is now positive and equals the unwrapped value).  This is the
#          DEFINITIVE proof that CPN.7.3's @test_broken is the LIVE 61um
#          auto-flip: removing the exact wrap that causes 61um turns the marker
#          green.
#     (ii) CPN.7.2 (COARSE) ALSO reddens (5 RED), a STRONGER result than first
#          predicted: the COARSE edge's two principal arg() values straddle the
#          arg=π ray, so the RAW difference is -4.538 rad (OUT of (-π,π]) and
#          the wrap is load-bearing for the coarse leg too -- removing it makes
#          wd=-4.538 (<0), so `wd>0` / `abs(wd)<π` / the +1 sheet pins all fail
#          and visited_sheet bumps -1.  Recorded honestly: the wrap is NOT
#          incidental to the coarse leg; both legs depend on it (the FINE leg's
#          dependence is the BUG, the COARSE leg's is the CORRECT use).
#     (iii) CPN.7.1 unaffected (closed-form ±√, no winding).
#     Restored byte-for-byte; `git diff src/` empty; file GREEN with exactly
#     one @test_broken.
#
# STANDALONE RUN:
#   julia --project=. test/corpus_pathnet_winding_test.jl
# ============================================================================
