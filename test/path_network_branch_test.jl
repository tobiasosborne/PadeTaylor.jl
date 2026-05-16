# test/path_network_branch_test.jl — A4 of FFW 2017 reproduction, ADR-0013.
#
# Integration tests for the branch-aware `path_network_solve` kwargs:
# `branch_points`, `branch_cut_angles`, `cross_branch`, `initial_sheet`.
# The geometric primitives are unit-tested separately in
# test/branch_tracker_test.jl.
#
# Ground truth: ADR-0013 §"Decision" + FFW2017...md:178 verbatim.
#
# The test ODE is `u' = u` with `u(0) = 1`, whose exact solution is
# `u(z) = exp(z)`.  This is entire (no actual branches), so branch
# enforcement is purely a routing constraint — the analytic solution
# does not change under any choice of cut.  This isolates the path-
# routing/sheet-bookkeeping behaviour from the dynamics; the load-
# bearing tests are: "default kwargs reproduce existing tree byte-
# identically" (backward-compat invariant), "refuse mode rejects
# crossings", "cross mode bumps sheet counter", "all-forbidden raises
# loud-fail".

using Test
using PadeTaylor

@testset "PathNetwork branch-aware kwargs (A4 / ADR-0013)" begin

    # ODE: u'' = u' (second-order placeholder dynamics whose exact
    # solution u(z) = exp(z) is entire — branches are pure routing
    # constraints, decoupled from the ODE).  IC sits at z = 1 + 0im,
    # OFF the branch point at z = 0 and OFF the cut along arg = π
    # (negative real axis).  Grid points all live in the upper-half
    # plane or on the right of the branch, so refuse-mode can reach
    # them without crossing.
    f       = (z, u, up) -> up
    z₀      = 1.0 + 0.0im
    u₀      = 1.0 + 0.0im
    up₀     = 1.0 + 0.0im
    zend    = 4.0 + 4.0im
    # Strictly upper-right grid: every target lies in Re ≥ 0, Im > 0,
    # so the wedge walker from z₀ = (1, 0) can reach them without
    # straying below y = 0 (where it would corner itself against the
    # cut at arg = π).  PNB.1.2 / PNB.1.4 use this grid.
    grid    = ComplexF64[1.0 + 1.0im, 0.5 + 0.5im, 0.5 + 1.0im]
    prob    = PadeTaylorProblem(f, (u₀, up₀), (z₀, zend); order = 20)

    # -----------------------------------------------------------------
    # PNB.1.1 — Backward compatibility: with `branch_points = ()`
    # (default), the new field `visited_sheet` is populated with empty
    # vectors and the rest of the solution is byte-identical to a
    # pre-A4 call.  This is the LOAD-BEARING invariant.
    # -----------------------------------------------------------------
    @testset "PNB.1.1: default branch_points=() is byte-equivalent" begin
        sol_a = path_network_solve(prob, grid; h = 0.5, rng_seed = 42)
        sol_b = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                    branch_points = ())
        @test sol_a.visited_z      == sol_b.visited_z
        @test sol_a.visited_u      == sol_b.visited_u
        @test sol_a.visited_h      == sol_b.visited_h
        @test sol_a.visited_parent == sol_b.visited_parent
        @test sol_a.grid_u         == sol_b.grid_u
        # New field present in both with the correct empty shape.
        @test all(==(Int[]), sol_a.visited_sheet)
        @test all(==(Int[]), sol_b.visited_sheet)
        @test length(sol_a.visited_sheet) == length(sol_a.visited_z)
    end

    # -----------------------------------------------------------------
    # PNB.1.2 — Refuse mode: the walker refuses to cross a cut.  Use
    # a branch at 0 with the standard log cut (arg = π).  The visited
    # tree must contain no edges that cross the cut, and the filter
    # helper itself must replace forbidden candidates with the
    # `(z_cur, Inf, 0, nothing)` failure sentinel.
    # -----------------------------------------------------------------
    @testset "PNB.1.2: refuse mode routes around the cut" begin
        sol = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                  branch_points = (0.0 + 0.0im,))
        # No visited edge may cross the cut.
        cuts_crossed = 0
        for k in 2:length(sol.visited_z)
            z_parent = sol.visited_z[sol.visited_parent[k]]
            z_child  = sol.visited_z[k]
            if PadeTaylor.BranchTracker.segment_crosses_cut(
                   z_parent, z_child, 0.0+0.0im, π)
                cuts_crossed += 1
            end
        end
        @test cuts_crossed == 0
        # Sheet counters all stay at 0 (no crossing under refuse).
        @test all(s -> s == [0], sol.visited_sheet)

        # Direct unit test of the filter helper: forbidden candidate
        # is rewritten to the `(z_cur, Inf, 0, nothing)` sentinel that
        # the downstream `_select_candidate` will skip; allowed
        # candidates pass through bit-equally.  Hand-build evals using
        # the IC visited_pade as a valid PadeApproximant sentinel.
        CT       = ComplexF64
        good     = sol.visited_pade[1]
        evals_in = Any[
            (CT(0.5 + 0.0im), CT(2.0 + 0.0im), CT(0.1), good),   # toward + real: no cross
            (CT(-1.0 - 0.5im), CT(0.3 + 0.0im), CT(0.1), good),  # crosses cut from above
            (CT(-1.0 + 0.5im), CT(0.4 + 0.0im), CT(0.1), good),  # doesn't cross
            (CT(0.5 + 0.5im), CT(0.5 + 0.0im), CT(0.1), good),   # no cross
            (CT(2.0 + 0.0im), CT(1.0 + 0.0im), CT(0.1), good),   # no cross
        ]
        z_cur_test = CT(-0.5 + 0.5im)   # above the cut
        out = PadeTaylor.PathNetwork._filter_forbidden_candidates(
            evals_in, z_cur_test, (0.0+0.0im,), (Float64(π),), CT)
        @test out[1][2] == CT(2.0 + 0.0im)    # unchanged
        @test out[2][2] == CT(Inf, 0)         # forbidden → sentinel
        @test out[2][4] === nothing
        @test out[3][2] == CT(0.4 + 0.0im)    # unchanged
        @test out[4][2] == CT(0.5 + 0.0im)    # unchanged
        @test out[5][2] == CT(1.0 + 0.0im)    # unchanged
    end

    # -----------------------------------------------------------------
    # PNB.1.3 — Cross mode: walker is permitted to cross; sheet counter
    # bumps on each crossing.  We extend the standard grid with a
    # lower-half-plane target (`-1 - 0.5im`) that REQUIRES crossing the
    # cut to reach from the upper-half IC.  Verify by counting
    # crossings along each path-tree chain and comparing against
    # `visited_sheet`.
    # -----------------------------------------------------------------
    @testset "PNB.1.3: cross mode bumps visited_sheet on each crossing" begin
        grid_cross = vcat(grid, [-1.0 - 0.5im, -0.5 - 1.0im])
        sol = path_network_solve(prob, grid_cross; h = 0.5, rng_seed = 42,
                                  branch_points = (0.0 + 0.0im,),
                                  cross_branch  = true)
        # We need at least one node with non-zero sheet (otherwise the
        # geometry-of-grid choice failed to force any crossing).
        @test any(s -> s != [0], sol.visited_sheet)
        # For every visited node, walk up the parent chain to the root,
        # summing the per-step sheet delta (final - parent.sheet).  This
        # must equal the absolute value of crossings along the chain
        # (one per crossed cut, with sign).
        for k in 1:length(sol.visited_z)
            # Reconstruct the path from root to this node.
            path = Int[k]
            while path[end] != 1
                pidx = sol.visited_parent[path[end]]
                @assert pidx > 0
                push!(path, pidx)
            end
            reverse!(path)
            # Verify sheet at each node along the chain agrees with
            # accumulated crossing-sign computed on the spot.
            expected = 0
            for j in 2:length(path)
                z_par  = sol.visited_z[path[j-1]]
                z_chld = sol.visited_z[path[j]]
                if PadeTaylor.BranchTracker.segment_crosses_cut(
                       z_par, z_chld, 0.0+0.0im, π)
                    Δθ = PadeTaylor.SheetTracker.winding_delta(
                            z_par, z_chld, 0.0+0.0im)
                    expected += Δθ > 0 ? 1 : (Δθ < 0 ? -1 : 0)
                end
                @test sol.visited_sheet[path[j]] == [expected]
            end
        end
    end

    # -----------------------------------------------------------------
    # PNB.1.4 — Initial-sheet kwarg seeds the root node.
    # -----------------------------------------------------------------
    @testset "PNB.1.4: initial_sheet seeds the root" begin
        sol = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                  branch_points = (0.0 + 0.0im,),
                                  initial_sheet = [7])
        # Refuse mode (default): no crossings, all nodes inherit root's sheet.
        @test all(s -> s == [7], sol.visited_sheet)
    end

    # -----------------------------------------------------------------
    # PNB.1.5 — Fail-loud when all candidates are forbidden.
    # Configuration: IC at (-0.5, 0), target at (-0.5, +1); branch at
    # (0, 0.3) with cut along arg = π (horizontal half-line y = 0.3,
    # x ≤ 0).  Every wedge candidate aimed at the upward target lands
    # at y > 0.3 with x < 0 — crosses the cut.  All 5 forbidden →
    # fail-loud with helpful cross_branch hint.
    # -----------------------------------------------------------------
    @testset "PNB.1.5: fail-loud message mentions cross_branch" begin
        prob_local = PadeTaylorProblem(
            f, (u₀, up₀),
            (-0.5 + 0.0im, -0.5 + 1.0im); order = 20)
        @test_throws ErrorException path_network_solve(
            prob_local, ComplexF64[-0.5 + 1.0im];
            h = 0.5, rng_seed = 42,
            branch_points     = (0.0 + 0.3im,),
            branch_cut_angles = π)
        try
            path_network_solve(prob_local, ComplexF64[-0.5 + 1.0im];
                               h = 0.5, rng_seed = 42,
                               branch_points     = (0.0 + 0.3im,),
                               branch_cut_angles = π)
            error("should have thrown")
        catch e
            @test e isa ErrorException
            @test occursin("cross_branch", e.msg)
        end
    end

    # -----------------------------------------------------------------
    # PNB.1.6 — Argument validation: cross_branch with no branch_points
    # throws; initial_sheet length mismatch throws.
    # -----------------------------------------------------------------
    @testset "PNB.1.6: kwarg validation throws on inconsistent inputs" begin
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, cross_branch = true)
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5,
            branch_points  = (0.0+0.0im, 1.0+0.0im),
            initial_sheet  = [0])           # wrong length
        # cross_branch with both branches: valid.
        sol = path_network_solve(prob, grid; h = 0.5, rng_seed = 42,
                                  branch_points  = (0.0+0.0im, 1.0+0.0im),
                                  initial_sheet  = [5, 11],
                                  cross_branch   = true)
        @test sol.visited_sheet[1] == [5, 11]
    end

    # -----------------------------------------------------------------
    # PNB.1.7 — Schwarz-reflection incompatibility (ADR-0013 §"Decision"):
    # cannot combine enforce_real_axis_symmetry with branch_points.
    # -----------------------------------------------------------------
    @testset "PNB.1.7: enforce_real_axis_symmetry + branch_points throws" begin
        prob_real = PadeTaylorProblem(
            (z, u, up) -> u, (1.0+0.0im, 0.0+0.0im),
            (0.0+0.0im, 1.0+0.0im); order = 15)
        @test_throws ArgumentError path_network_solve(
            prob_real, ComplexF64[0.5 + 0.5im, 0.5 - 0.5im];
            h = 0.3,
            enforce_real_axis_symmetry = true,
            branch_points              = (1.0 + 0.0im,))
    end

end # @testset PathNetwork branch-aware

# Mutation-proof procedure (verified before commit, 2026-05-16,
# at worklog 042; 75 GREEN total under unmutated impl):
#
#   Mutation PB1  --  in the wedge loop, drop the
#     `_filter_forbidden_candidates` call entirely (comment out the
#     `if branched && !cross_branch ... end` block).  Verified bite:
#     2 RED of 75 — PNB.1.5's @test_throws and the message-content
#     check both fire because the walker no longer detects the
#     all-forbidden case and instead steps freely through the cut.
#     PNB.1.2's direct unit test of the filter helper does NOT bite
#     this mutation (the helper itself is unchanged); the bite is the
#     wired-up-ness of the call from inside the wedge loop.
#
#   Mutation PB2  --  in the wedge-loop sheet-update block, always
#     copy the parent's sheet regardless of `cross_branch` flag.
#     Verified bite: 8 RED of 75, all in PNB.1.3 — the chain-walk
#     assertion `sol.visited_sheet[path[j]] == [expected]` fires
#     because `expected` accumulates non-zero values for nodes
#     reached by cut-crossing edges, but the mutated impl always
#     returns `[0]`.
#
#   Mutation PB3  --  in the validation block, drop the
#     `initial_sheet` length check (comment out).  Verified bite:
#     1 RED of 75 — PNB.1.6's wrong-length call no longer throws
#     `ArgumentError` (a downstream `@assert` from
#     `step_sheet_update` fires at solve time instead), so
#     `@test_throws ArgumentError` fails.
