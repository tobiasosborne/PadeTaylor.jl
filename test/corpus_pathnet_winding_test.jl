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
#     straddles the arg=pi cut with realised subtended ~+0.555*pi (< pi, CCW)
#     and clears the branch (grazing ratio ~0.20 > graze_tol=0.1).  winding_delta
#     is EXACT; visited_sheet=[1] is CORRECT (matches TRUE +1) and does NOT throw.
#   FINE -- IC arg -0.47*pi, target arg +0.59*pi, h=0.70 (padetaylor-61um
#     RESOLVED 2026-06-13).  The single step GRAZES branch=1 (realised grazing
#     ratio ~0.06 « graze_tol=0.1).  winding_delta on the chord is CORRECT
#     (-0.94*pi, the principal value; a straight chord never subtends |Δθ|>=pi
#     about an exterior branch).  Because the step grazes the branch the sheet
#     bump is winding-AMBIGUOUS, so step_sheet_update fails loud and the walk
#     THROWS.  (The "+1.083*pi true subtended" is the CURVED major arc, not the
#     straight chord; unrecoverable from two endpoints -- worklog 074:47-74.)
# Same ODE/branch/cut; only whether the chord CLEARS or GRAZES the branch
# differs -- proving the fail-loud is STEP-SIZE-TRIGGERED (shrink h to clear it).
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
# reddens) + M-guard-off (disable the step_sheet_update grazing guard -> the
# FINE @test_throws reddens) + M-guard-overfire (raise graze_tol above COARSE's
# ratio -> the SAFE COARSE leg reddens).  ACTUALLY EXECUTED then restored
# byte-clean.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor.SheetTracker: winding_delta, accumulate_winding, sheet_index
import PadeTaylor.BranchTracker: segment_crosses_cut, step_sheet_update

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
    # CPN.7.3  FINE leg (padetaylor-61um RESOLVED 2026-06-13) -- a REAL single-
    #          step cross_branch=true walk whose wedge step GRAZES branch=1
    #          (realised grazing ratio ~0.06 « graze_tol=0.1).  winding_delta on
    #          the chord is CORRECT (returns -0.94π, the principal value — a
    #          straight chord cannot subtend |Δθ| ≥ π about an exterior branch).
    #          Because the step grazes the branch the sheet bump is winding-
    #          AMBIGUOUS, so step_sheet_update REFUSES it fail-loud and the walk
    #          THROWS — surfacing the ill-conditioned step instead of silently
    #          corrupting the sheet index.  (The earlier framing — "true
    #          subtended +1.083π, winding_delta loses a revolution, visited_sheet
    #          corrupted to [-1]" — was geometrically false: +1.083π is the
    #          MAJOR-arc angle of a CURVED path the walker never realises,
    #          unrecoverable from two endpoints.  worklog 074:47-74; ADR-0031.)
    # -----------------------------------------------------------------------
    @testset "CPN.7.3: FINE step (grazes branch) ⇒ fail-loud, winding_delta CORRECT (61um)" begin
        # winding_delta on the IDEAL FINE chord (the documented ring geometry)
        # returns the chord's TRUE winding -0.94π = the principal arg-difference.
        # This is a pure, CORRECT primitive call (no throw — winding_delta has no
        # guard; the ambiguity is the WALKER's concern, handled below).
        @test segment_crosses_cut(CPN7_FINE_P, CPN7_FINE_Q, CPN7_BRANCH, π) == true
        wd = winding_delta(CPN7_FINE_P, CPN7_FINE_Q, CPN7_BRANCH)
        @test isapprox(wd, CPN7_FINE_IDEAL_WRAPPED; atol = 1e-9)   # = -0.94π, CORRECT
        @test wd < 0 && abs(wd) < π                                # never wraps past π
        @test isapprox(wd, angle((CPN7_FINE_Q - CPN7_BRANCH) /
                                 (CPN7_FINE_P - CPN7_BRANCH)); atol = 1e-12)
        # The +1.06π "ideal subtended" is the curved major arc, off by exactly
        # the +2π neither the chord nor two endpoints carry.
        @test isapprox(wd - CPN7_FINE_IDEAL_SUBTENDED, -2π; atol = 1e-9)

        # UNIT: the FINE chord crosses the arg=π cut AND grazes branch=1
        # (ratio ~0.047 < graze_tol), so step_sheet_update REFUSES the bump.
        @test_throws DomainError step_sheet_update([0], CPN7_FINE_P, CPN7_FINE_Q,
                                                   (CPN7_BRANCH,), (Float64(π),))
        gerr = try
            step_sheet_update([0], CPN7_FINE_P, CPN7_FINE_Q, (CPN7_BRANCH,), (Float64(π),))
        catch e; e end
        @test gerr isa DomainError
        @test occursin("graze", sprint(showerror, gerr))
        @test occursin("shorten h", sprint(showerror, gerr))

        # REAL WALK: the wedge stepper realises that same grazing step, so
        # path_network_solve now THROWS at the sheet-bump site (PathNetwork.jl:627)
        # — the 61um fix in action (was: silent visited_sheet=[-1] corruption).
        @test_throws DomainError cpn7_leg(CPN7_FINE_P, CPN7_FINE_Q, 0.70)
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
#   M-guard-off -- in src/BranchTracker.jl step_sheet_update, disable the 61um
#     grazing guard (e.g. set graze_tol default to 0.0, or delete the
#     `mind < graze_tol*|d|` throw).  Verified bite: CPN.7.3's two @test_throws
#     (the unit step_sheet_update + the cpn7_leg(FINE) walk) go RED — with no
#     guard the FINE grazing step no longer throws (it would silently bump
#     sign(winding_delta) = -1, the old corruption).  CBr.3's @test_throws
#     reddens identically.  Restored byte-clean.
#   M-guard-overfire -- raise graze_tol default to 0.3 (above CPN.7.2 COARSE's
#     realised grazing ratio 0.2045).  Verified bite: CPN.7.2 reddens — the SAFE
#     coarse step (which must return [1] WITHOUT throwing) now throws, and the
#     other in-suite cross_branch walks would too.  Proves the threshold is not
#     over-firing on well-conditioned steps.  Restored byte-clean.
#   Measured (campaign env, 2026-06-13): M-guard-off → CPN.7.3 5 RED (the unit
#   @test_throws + 3 message pins + the cpn7_leg(FINE) walk @test_throws) and
#   CBr.3 1 RED; M-guard-overfire (tol=0.3) → CPN.7.2 1 ERROR (the safe COARSE
#   walk throws).  Both reverted byte-clean (`git diff src/` empty).
#   CPN.7.1 unaffected (closed-form ±√, no winding).
#
# STANDALONE RUN:
#   julia --project=. test/corpus_pathnet_winding_test.jl
# ============================================================================
