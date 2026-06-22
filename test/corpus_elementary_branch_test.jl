# test/corpus_elementary_branch_test.jl
#
# ============================================================================
# Ground-truth corpus: ELEMENTARY BRANCH-POINT ODEs (CBr)
# (bead padetaylor-g9lg, epic padetaylor-25og).  Family E of the v2 corpus
# plan, docs/test_corpus/02_corpus_extension_plan.md:245-270.
#
# WHAT THIS FILE PINS.  The v1 corpus only ever drove the Riemann-sheet /
# winding machinery (src/SheetTracker.jl, src/BranchTracker.jl) through the
# PVI/PIII coordinate transforms or the *entire* ODE u''=u' (where branches
# are pure routing constraints, decoupled from the dynamics).  This file binds
# those primitives to genuinely MULTIVALUED closed forms whose monodromy is
# elementary and exactly known:
#
#   CBr.2  log single branch.   u' = 1/z, recast u'' = -1/z^2, u = log z.
#          Real probe (solve_pade from z0=1): u(2)=log 2, u(e)=1, u'(2)=1/2.
#          Winding (the primitive against log's monodromy): a CCW unit-circle
#          path [1, i, -1, -i, 1] about the branch at 0 accumulates total 2*pi,
#          sheet_index = +1, and log jumps exactly +2*pi*i per loop.
#
#   CBr.3  atanh two-branch (THE padetaylor-61um CATCHER — the headline).
#          u' = 2/(1-z^2), recast u'' = 4z/(1-z^2)^2, u = 2*atanh z.  Branch
#          points ±1.  Real probe (solve_pade from z0=0): u(1/2)=2*atanh(1/2)
#          = log 3; pinned via the DERIVATIVE identity u'=2/(1-z^2) and the
#          value 2*atanh — NOT the raw log-vs-atanh string (sympy shows a
#          cosmetic principal-branch relabel there; see capture.py).
#          THE 61um STEP: a single step about branch z=1 subtending |Δθ|≈1.1*pi
#          CCW (z_old = 1 + 0.1*exp(-0.45*pi*i), z_new = 1 + 0.1*exp(0.65*pi*i)).
#          The TRUE subtended angle is +1.1*pi; winding_delta clamps to (-pi,pi]
#          and returns the WRAPPED -0.9*pi — a SILENTLY LOST full revolution.
#          We @test the current (documented-wrong) -0.9*pi (GREEN today) and
#          @test_broken the Rule-1-correct +1.1*pi (the AUTO-FLIP marker for
#          padetaylor-61um: flips to "unexpectedly passes" the day 61um lands).
#
#   CBr.1  sqrt finite cut.   y' = -z/y, recast y'' = -2/y^3, y = sqrt(2-z^2).
#          Branch points ±sqrt(2), a FINITE cut between them.  Principal-sheet
#          (+sqrt) value via path_network_solve(branch_points=(-√2,√2)) +
#          eval_at_sheet([0,0]) vs +sqrt(2-z^2), plus the branch geometry.
#          The cross-sheet (-sqrt) eval is DEFERRED (see CBr.1 note).
#
# HOW CBr.3 DIFFERS FROM CWD.5 (corpus_winding_test.jl).  CWD.5 is a bare unit
# fixture: branch at z=0, BOTH endpoints ON the real axis (0.1 and 0.1*e^{i1.1π}
# whose principal arg is -0.9π).  CBr.3 hits the SAME bug under a *realistic
# two-branch walk* geometry: branch at z=1 (not 0); BOTH endpoints OFF the real
# axis at NONZERO arg relative to the branch (-0.45π and +0.65π); and it lives
# inside a real atanh ODE with two branch points ±1.  Same -2*pi lost
# revolution, genuinely distinct configuration.
#
# GROUND TRUTH (Law 1).  Every value is an INDEPENDENT sympy/mpmath pin from
# external/probes/corpus-oracles/elementary-branch/capture.py (dps=50, sympy
# residual==0 GATE on each recast + the atanh==log derivative identity + the
# ±√2 / ±1 branch-point sets before any literal is emitted).  Winding
# normalisation (-pi,pi] is principal-value arg (SheetTracker.jl:283-292);
# sheet convention is FFW2017...md:187-189.  Transcendental pins live in
# test/_oracle_corpus_elementary_branch.jl.
#
# TOLERANCES (justified by METHOD).
#   * Winding angles: exact rationals × π → atol 1e-12 (a few ulp of π).
#   * solve_pade real-axis values: Padé-Taylor at order 30 on a smooth real
#     arc bounded away from the branch → atol 1e-12.
#   * path_network_solve principal-sheet eval: Stage-2 nearest-node Padé,
#     path-tree truncation noise → atol 1e-9 (plan Family E tol).
#
# REFERENCES.  Plan Family E (02_corpus_extension_plan.md:245-270); bug
# padetaylor-61um (winding_delta unenforced |Δθ|≥π precondition); DLMF 4.x
# (log / atanh), DLMF 22.x analogue (sqrt finite cut).  Catalogue mirror of
# CWD.5 in test/corpus_winding_test.jl (the bare unit fixture).
#
# MUTATION-PROOF (Rule 4): footer — M-oracle (flip a pin → reddens) + M-impl
# (perturb SheetTracker.jl winding/sheet machinery → winding asserts redden),
# both ACTUALLY EXECUTED then restored byte-clean.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor.SheetTracker: winding_delta, accumulate_winding, sheet_index
import PadeTaylor.BranchTracker: segment_crosses_cut, step_sheet_update

include(joinpath(@__DIR__, "_oracle_corpus_elementary_branch.jl"))

@testset "Corpus: elementary branch points (CBr)" begin

    # -----------------------------------------------------------------------
    # CBr.2  log single branch — u'=1/z, recast u''=-1/z^2, u=log z.
    #   (a) real probe via solve_pade from z0=1 (u(1)=0, u'(1)=1);
    #   (b) winding primitive bound to log's monodromy: CCW unit circle about
    #       the branch at 0 ⇒ total 2π, sheet_index +1, log jump +2πi/loop.
    # -----------------------------------------------------------------------
    @testset "CBr.2: log single branch (real probe + CCW monodromy)" begin
        # (a) Real probe.  Recast u''=-1/z^2 as a 2nd-order RHS (z,u,up)->u''.
        rhs = (z, u, up) -> -1 / z^2
        prob = PadeTaylorProblem(rhs, (0.0, 1.0), (1.0, 3.0); order = 30)
        sol  = solve_pade(prob; h = 0.25)
        u2, up2 = sol(2.0)
        ue, _   = sol(ℯ)
        @test isapprox(u2, CBR2_U_AT_2; atol = 1e-12)    # log 2
        @test isapprox(up2, CBR2_UP_AT_2; atol = 1e-12)  # u'(2) = 1/2
        @test isapprox(ue, CBR2_U_AT_E; atol = 1e-12)    # log e = 1

        # (b) Winding about branch 0 — the same four-arc CCW unit circle as
        # CWD.1, here BOUND to log's monodromy.  Each Δ=π/2; total 2π.
        path = ComplexF64[1, im, -1, -im, 1]
        acc  = accumulate_winding(path, 0.0im)
        @test acc[1] == 0.0
        @test isapprox(acc[end], 2π; atol = 1e-12)
        @test sheet_index(acc[end]) == 1
        # log's monodromy: one CCW loop adds exactly 2πi (purely imaginary).
        # The accumulated winding IS the imaginary part of the log jump.
        @test isapprox(acc[end], CBR2_LOG_JUMP_PER_LOOP_IMAG; atol = 1e-12)
        # Two loops ⇒ sheet +2, jump +4πi.
        path2 = vcat(path, path[2:end])
        @test sheet_index(accumulate_winding(path2, 0.0im)[end]) == 2
    end

    # -----------------------------------------------------------------------
    # CBr.3  atanh two-branch — u'=2/(1-z^2), recast u''=4z/(1-z^2)^2,
    #        u=2*atanh z.  Branch points ±1.
    #   (a) real probe via solve_pade from z0=0 (u(0)=0, u'(0)=2): u(1/2)=log3;
    #   (b) THE padetaylor-61um STEP about branch z=1 (see header for the
    #       CWD.5 distinction).
    # -----------------------------------------------------------------------
    @testset "CBr.3: atanh two-branch (real probe + 61um step)" begin
        # (a) Real probe inside (-1,1).  Recast u''=4z/(1-z^2)^2.
        rhs = (z, u, up) -> 4z / (1 - z^2)^2
        prob = PadeTaylorProblem(rhs, (0.0, 2.0), (0.0, 0.5); order = 30)
        sol  = solve_pade(prob; h = 0.1)
        uh, uph = sol(0.5)
        @test isapprox(uh, CBR3_U_AT_HALF; atol = 1e-12)   # 2*atanh(1/2)
        @test isapprox(uh, CBR3_LOG3; atol = 1e-12)        # = log 3 (cross-check)
        @test isapprox(uph, CBR3_UP_AT_HALF; atol = 1e-12) # u'(1/2) = 8/3
        # Pin the DERIVATIVE identity u'=2/(1-z^2) directly (NOT the log string).
        @test isapprox(uph, 2 / (1 - 0.5^2); atol = 1e-12)

        # Branch points ±1 are the zeros of (1-z^2) — structural pin.
        @test (1 - 1.0^2) == 0.0 && (1 - (-1.0)^2) == 0.0

        # (b) THE 61um STEP (RESOLVED 2026-06-13).  Branch at z=1; both
        # endpoints on a radius-0.1 ring at arg -0.45π / +0.65π relative to the
        # branch — a NEARLY-ANTIPODAL pair, so the straight chord GRAZES the
        # branch (ratio ~0.079).  winding_delta returns the chord's TRUE winding
        # -0.9π (the principal value).  This is CORRECT, not a "lost revolution":
        # a straight chord can never subtend |Δθ| ≥ π about an exterior branch,
        # so the +1.1π "true subtended" is the MAJOR-arc angle of a CURVED path
        # the walker never realises — unrecoverable from two endpoints
        # (worklog 074:47-74; ADR-0031).  The genuine 61um defect was the
        # UNGUARDED grazing step in the sheet bookkeeping, now fixed: a crossed-
        # cut grazing step fails loud in step_sheet_update.
        branch = 1.0 + 0.0im
        z_old  = 1 + 0.1 * exp(im * (-0.45π))
        z_new  = 1 + 0.1 * exp(im * (+0.65π))
        actual = winding_delta(z_old, z_new, branch)

        # winding_delta returns the chord's TRUE winding -0.9π = the principal
        # arg-difference (it is a pure, CORRECT primitive).
        @test isapprox(actual, CBR3_61UM_WRAPPED; atol = 1e-12)   # -0.9π
        @test isapprox(actual, -0.9π; atol = 1e-12)
        @test abs(actual) < π                                     # never wraps past π
        @test isapprox(actual, angle((z_new - branch) / (z_old - branch));
                       atol = 1e-12)                              # = principal arg-diff
        # The +1.1π "true subtended" is the curved major arc, NOT the straight
        # chord; the chord and winding_delta agree on -0.9π (they differ from the
        # major arc by exactly the 2π neither the chord nor two endpoints carry).
        @test isapprox(CBR3_61UM_TRUE_SUBTENDED, 1.1π; atol = 1e-12)
        @test isapprox(actual - CBR3_61UM_TRUE_SUBTENDED, -2π; atol = 1e-12)

        # ----- 61um RESOLUTION (was @test_broken expecting +1.1π) ------------
        # The chord crosses the arg=π cut AND grazes branch=1 (ratio ~0.079 <
        # graze_tol=0.1), so the sheet bump is winding-AMBIGUOUS: step_sheet_update
        # REFUSES it fail-loud (Rule 1) rather than silently mis-accumulating.
        @test segment_crosses_cut(z_old, z_new, branch, π) == true
        @test_throws DomainError step_sheet_update([0], z_old, z_new,
                                                   (branch,), (Float64(π),))
    end

    # -----------------------------------------------------------------------
    # CBr.1  sqrt finite cut — y'=-z/y, recast y''=-2/y^3, y=sqrt(2-z^2), C=2.
    #        Branch points ±sqrt(2); a FINITE cut between them.  IC off the cut
    #        on the imaginary axis (z0=0.5i, y0=+sqrt(2-(0.5i)^2)=1.5).
    #        Principal-sheet (+sqrt) value via path_network_solve +
    #        eval_at_sheet([0,0]); the branch geometry; the -sqrt sheet pin.
    #
    # CROSS-SHEET DEFERRAL (Rule 9, documented).  Reaching the -sqrt sheet via
    # eval_at_sheet([0,1]) requires the wedge walker to deliberately cross the
    # FINITE two-point cut and accumulate a +1 sheet bump on a node near the
    # query point — fragile for a two-branch finite cut whose default arg=π
    # half-line cuts overlap.  Per the brief we ship the principal-sheet value
    # + branch-point geometry and DEFER the clean cross-sheet eval.  The -sqrt
    # sheet value is exactly the negation (pinned below as ground truth);
    # binding it to a walker route is the deferred follow-up.
    # -----------------------------------------------------------------------
    @testset "CBr.1: sqrt finite cut (principal sheet + branch geometry)" begin
        # Branch-point geometry: ±sqrt(2), a FINITE cut between them.
        @test isapprox(CBR1_BRANCH_POS, sqrt(2.0); atol = 1e-12)
        @test isapprox(CBR1_BRANCH_NEG, -sqrt(2.0); atol = 1e-12)
        @test isapprox(CBR1_BRANCH_POS - CBR1_BRANCH_NEG, 2sqrt(2.0); atol = 1e-12)

        # IC off the cut: z0=0.5i ⇒ y0 = +sqrt(2-(0.5i)^2) = +sqrt(2.25) = 1.5.
        z0 = 0.0 + 0.5im
        @test isapprox(sqrt(2 - z0^2), CBR1_Y0; atol = 1e-12)        # 1.5
        @test isapprox(-z0 / sqrt(2 - z0^2), CBR1_YP0; atol = 1e-12) # y'0 = -z0/y0

        # Principal-sheet (+sqrt) eval via the branch-aware path network.  The
        # query points sit on the upper imaginary axis, OFF the finite cut on
        # the real axis between ±√2, so the walk stays on sheet [0,0].
        rhs  = (z, y, yp) -> -2 / y^3                  # recast y''=-2/y^3
        prob = PadeTaylorProblem(rhs, (CBR1_Y0, CBR1_YP0),
                                 (z0, 0.0 + 2.0im); order = 20)
        grid = ComplexF64[0.0 + 0.5im, 0.0 + 1.0im, 0.3 + 0.5im]
        sol  = path_network_solve(prob, grid; h = 0.25, rng_seed = 42,
                                   branch_points = (CBR1_BRANCH_NEG + 0.0im,
                                                    CBR1_BRANCH_POS + 0.0im))
        # eval_at_sheet on the IC's principal sheet [0,0] must reproduce +sqrt.
        for (z, yref) in ((0.0 + 0.5im, CBR1_YPRIN_AT_0p5i),
                          (0.0 + 1.0im, CBR1_YPRIN_AT_1i),
                          (0.3 + 0.5im, CBR1_YPRIN_AT_0p3_0p5i))
            yv, _ = eval_at_sheet(sol, z, [0, 0])
            @test isfinite(yv)
            @test isapprox(yv, yref; atol = 1e-9)        # = +sqrt(2-z^2)
        end

        # The opposite (-sqrt) sheet value is exactly the negation — pinned as
        # GROUND TRUTH (capture.py).  Binding it to a walker route is DEFERRED.
        @test isapprox(CBR1_YMINUS_AT_0p5i, -CBR1_YPRIN_AT_0p5i; atol = 1e-12)
    end

end # @testset Corpus elementary branch (CBr)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
#   M-oracle — flip a pinned value.  In test/_oracle_corpus_elementary_branch.jl
#     set CBR3_61UM_WRAPPED from -2.8274...(-0.9π) to its negation (+0.9π).
#     Verified bite: CBr.3 RED — the `isapprox(actual, CBR3_61UM_WRAPPED)`
#     assertion fails (actual is the wrapped -0.9π; the flipped pin is +0.9π).
#     Restored the literal; file GREEN again.
#
#   M-impl — perturb src/SheetTracker.jl winding/sheet machinery.
#     (i) In winding_delta, drop the (-π,π] wrap (delete the
#         `if Δθ ≤ -π … elseif Δθ > π …` block; return raw Δθ).  Verified
#         bite: CBr.2 RED — the four-arc step 3 raw -3π/2 is no longer
#         normalised to +π/2, so the accumulated total ≠ 2π and sheet_index
#         ≠ 1; CBr.1 unaffected (no |Δθ|≥π/2 steps there).  Note this mutation
#         ALSO makes CBr.3's @test_broken UNEXPECTEDLY PASS (raw +1.1π is now
#         returned) — i.e. it demonstrates the 61um auto-flip mechanism live.
#     (ii) In sheet_index, change `round(total / 2π)` to `floor`.  Verified
#         bite: CBr.2 RED — sheet_index(2π) rounds to 1 but floors to 0 (2π/2π
#         lands on the boundary with fp noise) / the two-loop 4π case shifts.
#     Both restored byte-for-byte; `git diff src/` empty.
#
# STANDALONE RUN:
#   julia --project=. test/corpus_elementary_branch_test.jl
# ============================================================================
