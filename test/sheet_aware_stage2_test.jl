# test/sheet_aware_stage2_test.jl — A5 of FFW 2017 reproduction.
# Bead `padetaylor-hed`.
#
# Tests for `path_network_solve`'s `grid_sheet` kwarg and the
# `eval_at_sheet` post-hoc accessor.  Both are sheet-aware variants
# of the existing Stage-2 nearest-visited lookup: they restrict the
# pool to visited nodes whose `visited_sheet[k]` equals the query's
# sheet.  Composes with A4's branch_points / cross_branch.
#
# Ground truth: docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md
# §"Open follow-ups: A5 sheet-aware Stage-2".

using Test
using PadeTaylor

@testset "sheet-aware Stage-2 (A5 of FFW 2017 / bead padetaylor-hed)" begin

    # ODE u'' = u' (so u = exp(z) with the chosen IC).  Branches are
    # pure routing constraints — solution is entire.  Setup mirrors
    # path_network_branch_test.jl with a grid that REQUIRES crossings
    # under cross_branch=true (so visited_sheet has multiple distinct
    # values for the tests to exercise).
    f       = (z, u, up) -> up
    z₀      = 1.0 + 0.0im
    u₀      = 1.0 + 0.0im
    up₀     = 1.0 + 0.0im
    zend    = 4.0 + 4.0im
    # Mixed upper + lower targets force the walker to cross the
    # arg=π cut at the origin (under cross_branch=true).
    grid    = ComplexF64[1.0 + 1.0im, 0.5 + 0.5im, 0.5 + 1.0im,
                          -1.0 - 0.5im, -0.5 - 1.0im]
    prob    = PadeTaylorProblem(f, (u₀, up₀), (z₀, zend); order = 20)

    # -----------------------------------------------------------------
    # SA.1.1 — Backward compatibility: grid_sheet=nothing (default)
    # gives byte-identical results to a call without the kwarg.
    # -----------------------------------------------------------------
    @testset "SA.1.1: grid_sheet=nothing default unchanged" begin
        sol_a = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                    branch_points = (0.0+0.0im,),
                                    cross_branch  = true)
        sol_b = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                    branch_points = (0.0+0.0im,),
                                    cross_branch  = true,
                                    grid_sheet    = nothing)
        @test sol_a.grid_u  == sol_b.grid_u
        @test sol_a.grid_up == sol_b.grid_up
    end

    # -----------------------------------------------------------------
    # SA.1.2 — grid_sheet matching the IC's sheet (all zeros) returns
    # *some* finite values for grid points reached via non-crossing
    # paths, and NaN where the nearest matching-sheet visited node is
    # too far away (or there is no matching-sheet node at all).
    # -----------------------------------------------------------------
    @testset "SA.1.2: grid_sheet=[0] restricts pool to root-sheet nodes" begin
        sol = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                  branch_points = (0.0+0.0im,),
                                  cross_branch  = true,
                                  grid_sheet    = fill(Int[0], length(grid)))
        # At least the upper-half grid points should resolve to a
        # finite value (they're reachable from the IC without
        # crossing — they stay on sheet 0).
        @test isfinite(sol.grid_u[1])
        @test isfinite(sol.grid_u[2])
        @test isfinite(sol.grid_u[3])
    end

    # -----------------------------------------------------------------
    # SA.1.3 — Sheet mismatch → NaN.  Query a grid point against a
    # sheet that no visited node lies on (e.g., sheet [99]).  Every
    # output must be NaN.
    # -----------------------------------------------------------------
    @testset "SA.1.3: grid_sheet of an absent sheet returns NaN" begin
        sol = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                  branch_points = (0.0+0.0im,),
                                  cross_branch  = true,
                                  grid_sheet    = fill(Int[99], length(grid)))
        @test all(z -> isnan(real(z)), sol.grid_u)
        @test all(z -> isnan(real(z)), sol.grid_up)
    end

    # -----------------------------------------------------------------
    # SA.1.4 — Sheet of a grid point can DIFFER from the IC's sheet:
    # querying a lower-half grid point against sheet [-1] (the sheet
    # reached by CW crossing of the cut from above) should resolve
    # to a finite value when the walker actually visited a node on
    # that sheet near the grid point.  When the walker DID NOT reach
    # sheet -1 (the cross_branch happens to walk only to +1 / 0),
    # the test passes vacuously since NaN is also "fail-soft" output.
    # We assert a *consistency* property instead: querying the same
    # point against sheet [k] yields exp(z) when finite (since the
    # ODE u'=u has a single-valued solution u=exp).
    # -----------------------------------------------------------------
    @testset "SA.1.4: sheet-selected eval returns the analytic value when finite" begin
        sol = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                  branch_points = (0.0+0.0im,),
                                  cross_branch  = true)
        # The walker reached SOME sheets; collect the distinct sheet
        # values it visited.
        sheets_visited = unique(sol.visited_sheet)
        for sheet in sheets_visited
            sheet_grid = fill(sheet, length(grid))
            sol_k = path_network_solve(
                prob, grid; h = 0.5, rng_seed = 42,
                branch_points = (0.0+0.0im,),
                cross_branch  = true,
                grid_sheet    = sheet_grid)
            for (i, z) in enumerate(grid)
                if isfinite(sol_k.grid_u[i])
                    # u(z) = exp(z - z₀) by IC (z₀=1, u(z₀)=1); loose
                    # tolerance to absorb path-tree truncation noise.
                    @test isapprox(sol_k.grid_u[i], exp(z - z₀); rtol = 1e-3)
                end
            end
        end
    end

    # -----------------------------------------------------------------
    # SA.1.5 — eval_at_sheet accessor: same Stage-2 output as the
    # vectorised grid_sheet call, for a single point.
    # -----------------------------------------------------------------
    @testset "SA.1.5: eval_at_sheet matches the vectorised path" begin
        sheet_choice = Int[0]
        sol_vec = path_network_solve(
            prob, grid; h = 0.5, rng_seed = 42,
            branch_points = (0.0+0.0im,),
            cross_branch  = true,
            grid_sheet    = fill(sheet_choice, length(grid)))
        sol_for_eval = path_network_solve(
            prob, grid; h = 0.5, rng_seed = 42,
            branch_points = (0.0+0.0im,),
            cross_branch  = true)
        for (i, z) in enumerate(grid)
            u_acc, up_acc = eval_at_sheet(sol_for_eval, z, sheet_choice)
            # Both should be NaN or both finite + equal.
            if isfinite(sol_vec.grid_u[i])
                @test u_acc == sol_vec.grid_u[i]
                @test up_acc == sol_vec.grid_up[i]
            else
                @test isnan(real(u_acc))
            end
        end
    end

    # -----------------------------------------------------------------
    # SA.1.6 — eval_at_sheet shape validation throws on wrong-length
    # sheet input.  Two sub-cases: (a) solution has branches but
    # sheet is wrong length; (b) solution has no branches but sheet
    # is non-empty.
    # -----------------------------------------------------------------
    @testset "SA.1.6: eval_at_sheet input validation" begin
        # Use cross_branch=true so the walker can navigate the
        # lower-half grid points without cornering itself against
        # the cut.
        sol_branched = path_network_solve(
            prob, grid; h = 0.5, rng_seed = 42,
            branch_points = (0.0+0.0im, 1.0+1.0im),
            initial_sheet = [0, 0],
            cross_branch  = true)
        @test_throws ArgumentError eval_at_sheet(sol_branched, 1.0+0im, [0])
        @test_throws ArgumentError eval_at_sheet(sol_branched, 1.0+0im, [0, 0, 0])
        sol_unbranched = path_network_solve(prob, grid; h = 0.5, rng_seed = 42)
        @test_throws ArgumentError eval_at_sheet(sol_unbranched, 1.0+0im, [0])
    end

    # -----------------------------------------------------------------
    # SA.1.7 — grid_sheet kwarg shape validation throws before the
    # Stage-1 walk starts (so an expensive solve isn't burned on a
    # malformed input).
    # -----------------------------------------------------------------
    @testset "SA.1.7: grid_sheet kwarg shape validation" begin
        # Length mismatch: grid has 5 entries, grid_sheet has 3.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5,
            branch_points = (0.0+0.0im,),
            grid_sheet    = fill(Int[0], 3))
        # Inner-vector length mismatch: branch_points has 2 entries,
        # grid_sheet[i] has 1.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5,
            branch_points = (0.0+0.0im, 1.0+0.0im),
            grid_sheet    = fill(Int[0], length(grid)))
    end

    # -----------------------------------------------------------------
    # SA.1.8 — Direct unit test of _nearest_visited_on_sheet helper:
    # hand-built visited arrays exercise sheet-matching + tiebreak +
    # not-found cases.
    # -----------------------------------------------------------------
    @testset "SA.1.8: _nearest_visited_on_sheet helper" begin
        # Fixture: 5 visited nodes split across two sheets.
        #  index | z           | sheet
        #    1   | (0,   0)    | [1]
        #    2   | (0.5, 0)    | [0]
        #    3   | (-0.5, 0)   | [0]
        #    4   | (0.3, 0)    | [1]
        #    5   | (-0.3, 0)   | [1]
        visited_z     = ComplexF64[0.0+0.0im, 0.5+0.0im, -0.5+0.0im,
                                    0.3+0.0im, -0.3+0.0im]
        visited_sheet = [Int[1], Int[0], Int[0], Int[1], Int[1]]
        helper        = PadeTaylor.PathNetwork._nearest_visited_on_sheet
        # Query (0,0) on sheet [1]: distances 0, 0.3, 0.3 → node 1
        # wins by distance 0.
        @test helper(visited_z, visited_sheet, 0.0+0.0im, [1]) == 1
        # Query (0.4, 0) on sheet [0]: distances 0.1, 0.9 → node 2.
        @test helper(visited_z, visited_sheet, 0.4+0.0im, [0]) == 2
        # Query (0, 0) on sheet [99] (no match): returns 0.
        @test helper(visited_z, visited_sheet, 0.0+0.0im, [99]) == 0
        # Tiebreak: query (0, 0) on sheet [0] — nodes 2 and 3 are
        # both at distance 0.5.  Lex tiebreak (smaller Re first)
        # picks node 3 (Re = -0.5 < +0.5).
        @test helper(visited_z, visited_sheet, 0.0+0.0im, [0]) == 3
        # Tiebreak among sheet-[1] nodes: query at (0, +0.1) →
        # node 1 at (0, 0) wins (distance 0.1 vs nodes 4 (0.316)
        # and 5 (0.316)).
        @test helper(visited_z, visited_sheet, 0.0+0.1im, [1]) == 1
    end

end # @testset sheet-aware Stage-2

# Mutation-proof procedure (verified before commit, 2026-05-16,
# at worklog 043; 30 GREEN under unmutated impl):
#
#   Mutation S1  --  in `_nearest_visited_on_sheet`, drop the
#     `visited_sheet[i] == sheet || continue` filter (so the function
#     becomes a plain `_nearest_visited` ignoring sheets).  Verified
#     bite: 4 RED of 30 — SA.1.3 (2): outputs read finite values
#     instead of NaN because the nearest visited node by Euclidean
#     distance has some real sheet, no longer filtered out; SA.1.8 (2):
#     sheet-[0] tiebreak query at (0,0) picks node 1 (the absolute-
#     nearest at distance 0, sheet [1]) instead of node 3.
#
#   Mutation S2  --  in Stage-2 dispatch, ignore `grid_sheet` (always
#     call the unrestricted `_nearest_visited`).  Verified bite:
#     6 RED of 30 — SA.1.3 (2): outputs become finite when
#     grid_sheet specifies an absent sheet; SA.1.5 (4): the
#     vectorised path's per-point u/up differ from `eval_at_sheet`'s
#     per-point u/up because the vectorised path lost its sheet
#     restriction.  SA.1.7 cannot bite this (validation happens
#     before the dispatch).
#
#   Mutation S3  --  in `eval_at_sheet`, drop the sheet-aware lookup
#     and call the unrestricted `_nearest_visited` directly.
#     Verified bite: 2 RED of 30 — SA.1.5 — `eval_at_sheet`'s output
#     diverges from the vectorised `grid_sheet=...` path on the
#     query points where the unrestricted pool picks a different
#     node than the sheet-restricted pool.
