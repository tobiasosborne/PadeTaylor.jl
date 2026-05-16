# test/extrapolate_test.jl — ADR-0015 / bead `padetaylor-afs`.
#
# Tests for the `extrapolate::Bool = false` kwarg on:
#   - `path_network_solve` (Stage-2 cell loop)
#   - `eval_at_sheet`      (A5 sheet-aware per-point accessor)
#   - `eval_at`            (new sheet-blind per-point accessor)
#
# Ground truth: ADR-0015 §"Decision" + FFW md:62 verbatim spec
# ("Padé steps are taken from each ζ_i to the points on the fine
# grid to which it is the closest point among the ζ_i" — no disc
# check).

using Test
using PadeTaylor

@testset "extrapolate kwarg (ADR-0015 / bead padetaylor-afs)" begin

    # u'' = u' → u(z) = exp(z - z₀).  IC at z = 0 + 0im, u₀ = u'₀ = 1.
    f       = (z, u, up) -> up
    z₀      = 0.0 + 0.0im
    u₀      = 1.0 + 0.0im
    up₀     = 1.0 + 0.0im
    prob    = PadeTaylorProblem(f, (u₀, up₀), (z₀, 5.0+0im); order = 20)

    # Walk to a sparse set of targets near the IC.  With h = 0.5
    # canonical step, the visited tree's nodes have visited_h ≈ 0.5,
    # so cells more than ~0.5 from any visited node fall in gaps.
    targets = ComplexF64[0.5+0.0im, 1.0+0.0im, 1.5+0.0im, 2.0+0.0im]
    sol = path_network_solve(prob, targets; h = 0.5, rng_seed = 0)

    # ---- SX.1.1 — eval_at default returns NaN beyond disc ----------
    @testset "SX.1.1: eval_at default fail-soft beyond disc" begin
        # A query point far from any visited node.  Visited nodes are
        # at {0, 0.5, 1.0, 1.5, 2.0} (the targets + IC at distance 0
        # from the closest); query at z = 10 is far from all of them.
        u, up = eval_at(sol, 10.0 + 0.0im)
        @test isnan(real(u))
        @test isnan(imag(u))
        @test isnan(real(up))
        @test isnan(imag(up))
    end

    # ---- SX.1.2 — eval_at extrapolate=true returns finite ----------
    @testset "SX.1.2: eval_at extrapolate=true returns finite" begin
        u, up = eval_at(sol, 10.0 + 0.0im; extrapolate = true)
        @test isfinite(real(u))
        @test isfinite(imag(u))
        @test isfinite(real(up))
        @test isfinite(imag(up))
        # Note: the Padé is evaluated WAY past |t| = 1 (~17 units
        # from nearest visited node at h=0.5), so the value will
        # be very inaccurate — but it IS finite, which is what
        # `extrapolate=true` promises.  We don't pin its value
        # (degenerate extrapolation, not a meaningful target).
    end

    # ---- SX.1.3 — eval_at inside disc: extrapolate kwarg no-op ----
    @testset "SX.1.3: extrapolate is no-op inside disc" begin
        # Query at z = 0.5 sits at a visited node (distance 0);
        # both modes return the same value.
        u_off, _ = eval_at(sol, 0.5 + 0.0im)
        u_on,  _ = eval_at(sol, 0.5 + 0.0im; extrapolate = true)
        @test u_off == u_on
        @test isfinite(real(u_off))
        # And in the middle of a disc (between two visited nodes).
        u_mid_off, _ = eval_at(sol, 0.75 + 0.0im)
        u_mid_on,  _ = eval_at(sol, 0.75 + 0.0im; extrapolate = true)
        @test u_mid_off == u_mid_on
        @test isfinite(real(u_mid_off))
    end

    # ---- SX.1.4 — path_network_solve(grid) kwarg accepted ---------
    # Note: path_network_solve walks to EVERY grid cell as a target,
    # so its Stage-2 lookup always finds the cell itself as the
    # nearest visited node (distance 0).  The extrapolate kwarg
    # therefore rarely changes path_network_solve's grid_u output
    # in normal operation; both modes return finite values for every
    # cell.  This test confirms the kwarg is ACCEPTED + the outputs
    # AGREE on cells inside disc.  The figure-script use case
    # (sparse walker targets + dense render lattice via eval_at) is
    # covered by SX.1.6 below.
    @testset "SX.1.4: path_network_solve accepts extrapolate kwarg" begin
        grid = ComplexF64[0.5+0im, 1.5+0im, 3.0+0im]
        sol_off = path_network_solve(prob, grid; h = 0.5, rng_seed = 0)
        sol_on  = path_network_solve(prob, grid; h = 0.5, rng_seed = 0,
                                      extrapolate = true)
        # All cells walked-to → both modes return same finite values.
        for i in 1:3
            @test isfinite(real(sol_off.grid_u[i]))
            @test sol_off.grid_u[i] == sol_on.grid_u[i]
        end
    end

    # ---- SX.1.6 — figure-script pattern: sparse walker + dense eval -
    # Walk to a sparse set of targets, then evaluate at a DENSE
    # render lattice via eval_at.  Cells far from any visited node
    # are NaN by default; extrapolate=true fills them.
    @testset "SX.1.6: sparse walker + dense eval_at fills with extrapolate" begin
        sparse_targets = ComplexF64[1.0+0im, 2.0+0im]
        sol_sparse = path_network_solve(prob, sparse_targets;
                                         h = 0.5, rng_seed = 0)
        # Render lattice including cells far from any visited node.
        render_cells = ComplexF64[0.5+0im, 1.0+0im, 8.0+0im, -3.0+0im]
        rendered_off = [eval_at(sol_sparse, z)[1] for z in render_cells]
        rendered_on  = [eval_at(sol_sparse, z; extrapolate = true)[1]
                         for z in render_cells]
        # Cells inside disc: agreement.
        for i in 1:2
            @test isfinite(real(rendered_off[i]))
            @test rendered_off[i] == rendered_on[i]
        end
        # Far cells: NaN by default, finite under extrapolate.
        for i in 3:4
            @test isnan(real(rendered_off[i]))
            @test isfinite(real(rendered_on[i]))
        end
    end

    # ---- SX.1.5 — eval_at_sheet honours kwarg ---------------------
    @testset "SX.1.5: eval_at_sheet kwarg honoured" begin
        # Branchless setup so visited_sheet is empty for every node;
        # eval_at_sheet with empty sheet must still respect kwarg.
        # Far cell:
        u_off, _ = eval_at_sheet(sol, 10.0+0im, Int[])
        u_on,  _ = eval_at_sheet(sol, 10.0+0im, Int[]; extrapolate = true)
        @test isnan(real(u_off))
        @test isfinite(real(u_on))
        # Near cell (inside disc):
        u_off_n, _ = eval_at_sheet(sol, 0.5+0im, Int[])
        u_on_n,  _ = eval_at_sheet(sol, 0.5+0im, Int[]; extrapolate = true)
        @test u_off_n == u_on_n
        @test isfinite(real(u_off_n))
    end

end # @testset extrapolate

# Mutation-proof procedure (verified before commit, 2026-05-16,
# at worklog 045; 30 GREEN under unmutated impl):
#
#   Mutation X1  --  in path_network_solve Stage-2 loop, drop the
#     `extrapolate ||` guard (replace `!extrapolate && ...` with
#     `false && ...`).  Verified bite: 0 RED — path_network_solve
#     walks to EVERY grid cell, so the Stage-2 NaN path is dead
#     code in normal operation.  The kwarg is structurally wired
#     but its observable effect on grid_u is contingent on the
#     pathological "grid cell not reached by walker" case (which
#     normally throws).  Tests SX.1.4 cover threading without
#     biting the eval branch.  This is honest mutation-proof
#     coverage with documented gap: X1 doesn't bite because the
#     mutation has no observable effect in the test geometry.
#     The eval_at + eval_at_sheet kwargs (X2 + X3) carry the
#     load-bearing weight in test coverage.
#
#   Mutation X2  --  in eval_at, drop the `extrapolate ||` guard.
#     Verified bite: 6 RED of 30 — SX.1.1 (4): all four isnan
#     assertions for the eval_at far-cell case fail (cell returns
#     finite under mutation); SX.1.6 (2): the rendered_off far-cell
#     `isnan` assertions fail.
#
#   Mutation X3  --  in eval_at_sheet, drop the `extrapolate ||`
#     guard.  Verified bite: 1 RED of 30 — SX.1.5 — the far-cell
#     eval_at_sheet returns finite instead of NaN under default
#     extrapolate=false.
