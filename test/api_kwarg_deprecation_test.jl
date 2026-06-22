# test/api_kwarg_deprecation_test.jl -- beads `xds` + `0xn` deprecation shims.
#
# v1.0 public-API canonicalisation (docs/api-review-2026-05-16.md §3(a).1,
# §3(a).2, §9 rows 1-2) renamed three step-size kwargs to `h`
# (`solve_pade`/`PadeTaylorAlg` `h_max`, `lattice_dispatch_solve` `h_path`)
# and the BVP iteration-cap kwargs to `max_iter` (`bvp_solve` /
# `vector_bvp_solve` `maxiter`; `dispatch_solve` `bvp_maxiter` →
# `bvp_max_iter`).  Each rename ships a BACKWARD-COMPATIBLE shim: the old
# name is still accepted, maps onto the new one, and emits a
# `Base.depwarn`.  This file is the regression guard for the shims.
#
# Two invariants per shim:
#   (1) WORKS — the deprecated name produces the EXACT same result as the
#       canonical name (not merely "doesn't throw"; we compare outputs).
#   (2) WARNS — the deprecated name emits a deprecation warning.
#
# Depwarn mode caveat (verified empirically 2026-06-20): `Base.depwarn`
# only emits when `Base.JLOptions().depwarn != 0`.  `Pkg.test()` launches
# its test subprocess with `--depwarn=yes` (Pkg `Operations.jl:2147`), so
# the always-on correctness gate exercises invariant (2).  A standalone
# `julia --project=. test/api_kwarg_deprecation_test.jl` inherits flag 0
# (warnings suppressed); there invariant (2) is SKIPPED with an explicit
# `@test_skip` so the file still runs green standalone while the WORKS
# invariant (1) — the load-bearing one — always runs.

using Test
using PadeTaylor

# Whether the running process actually emits depwarns.  `Pkg.test()` → yes.
const _DEPWARN_ON = Base.JLOptions().depwarn != 0

"Assert `f()` emits a `:warn`-level log iff depwarns are on; always returns f()."
macro assert_depwarn(expr)
    quote
        if _DEPWARN_ON
            @test_logs (:warn,) match_mode = :any $(esc(expr))
        else
            # Warnings suppressed in this mode; the WORKS assertion below
            # is the meaningful guard.  Flag the skipped log-assertion so
            # the count is honest (Rule 5).
            @test_skip false
        end
    end
end

@testset "API kwarg deprecation shims (beads xds, 0xn)" begin

    # =====================================================================
    # xds.1 — solve_pade: `h_max` → `h`
    # =====================================================================
    @testset "solve_pade: h_max → h" begin
        fW(z, u, up) = u                       # u'' = u ⇒ u = cosh, finite
        prob = PadeTaylorProblem(fW, (1.0, 0.0), (0.0, 0.3); order = 30)

        sol_new = solve_pade(prob; h = 0.1)
        sol_old = solve_pade(prob; h_max = 0.1)
        @test sol_old.z == sol_new.z           # WORKS: identical trajectory
        @test sol_old.y == sol_new.y
        @test sol_old.h == sol_new.h
        @assert_depwarn solve_pade(prob; h_max = 0.1)   # WARNS
    end

    # =====================================================================
    # xds.2 — PadeTaylorAlg constructor: `h_max` → `h`
    # =====================================================================
    @testset "PadeTaylorAlg: h_max → h" begin
        alg_new = PadeTaylorAlg(; h = 0.5)
        alg_old = PadeTaylorAlg(; h_max = 0.5)
        @test alg_old.h == alg_new.h == 0.5    # WORKS: field carries through
        @test alg_old.max_steps == alg_new.max_steps
        @assert_depwarn PadeTaylorAlg(; h_max = 0.5)    # WARNS
    end

    # =====================================================================
    # xds.3 — lattice_dispatch_solve: `h_path` → `h`
    # =====================================================================
    @testset "lattice_dispatch_solve: h_path → h" begin
        # Linear u'' = u with a supplied interior mask (cosh closed form);
        # mirrors LD.1.2's minimal bridgeable setup.
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)

        N  = 11
        xs = range(-1.5, 1.5; length = N)
        ys = range(-1.5, 1.5; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(1.5 * sqrt(2)))
        prob  = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 20)
        mask  = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end

        sol_new = lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                         h = 0.5, order = 20, mask = mask)
        sol_old = lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                         h_path = 0.5, order = 20, mask = mask)
        @test sol_old.grid_u == sol_new.grid_u  # WORKS: identical field
        @test sol_old.region_tag == sol_new.region_tag
        @assert_depwarn lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                               h_path = 0.5, order = 20,
                                               mask = mask)               # WARNS
    end

    # =====================================================================
    # 0xn.1 — bvp_solve (2-arg): `maxiter` → `max_iter`
    # =====================================================================
    @testset "bvp_solve (2-arg): maxiter → max_iter" begin
        # Linear u'' = u on [-1, 1] with u(±1) = 1 (cosh-shaped; converges
        # in one Newton step, so any positive cap suffices).
        f(z, u)   = u
        ∂f(z, u)  = one(u)

        sol_new = bvp_solve(f, ∂f, -1.0, 1.0, 1.0, 1.0; N = 8, max_iter = 10)
        sol_old = bvp_solve(f, ∂f, -1.0, 1.0, 1.0, 1.0; N = 8, maxiter = 10)
        @test sol_old.u_nodes == sol_new.u_nodes      # WORKS: identical solve
        @test sol_old.iterations == sol_new.iterations
        @assert_depwarn bvp_solve(f, ∂f, -1.0, 1.0, 1.0, 1.0;
                                  N = 8, maxiter = 10)                     # WARNS
    end

    # =====================================================================
    # 0xn.2 — vector_bvp_solve: `maxiter` → `max_iter`
    # =====================================================================
    @testset "vector_bvp_solve: maxiter → max_iter" begin
        f1(z, y) = [y[1]]                              # y' = y, linear
        Ba1 = reshape([1.0], 1, 1)
        Bb1 = reshape([0.0], 1, 1)
        g1  = [1.0]

        sol_new = vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1; N = 8, max_iter = 10)
        sol_old = vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1; N = 8, maxiter = 10)
        @test sol_old.Y_nodes == sol_new.Y_nodes       # WORKS: identical solve
        @test sol_old.iterations == sol_new.iterations
        @assert_depwarn vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1;
                                         N = 8, maxiter = 10)             # WARNS
    end

    # =====================================================================
    # 0xn.3 — dispatch_solve: `bvp_maxiter` → `bvp_max_iter`
    # =====================================================================
    @testset "dispatch_solve: bvp_maxiter → bvp_max_iter" begin
        # Lone BVPSegment passthrough (mirrors DP.1.2); linear u'' = u on
        # [-1, 1], u(±1) = 1.
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)

        prob = PadeTaylorProblem(f_lin_2, (1.0, 0.0), (-1.0, 1.0); order = 30)
        seg  = BVPSegment(ComplexF64(1.0), ComplexF64(1.0),
                          ComplexF64[-0.5, 0.0, 0.5])

        sol_new = dispatch_solve(prob, f_lin_1, ∂f_lin_1, [seg];
                                 h = 0.5, N_bvp = 8, bvp_max_iter = 10)
        sol_old = dispatch_solve(prob, f_lin_1, ∂f_lin_1, [seg];
                                 h = 0.5, N_bvp = 8, bvp_maxiter = 10)
        @test sol_old.grid_u == sol_new.grid_u         # WORKS: identical solve
        @test sol_old.bvp_solutions[1].u_nodes == sol_new.bvp_solutions[1].u_nodes
        @assert_depwarn dispatch_solve(prob, f_lin_1, ∂f_lin_1, [seg];
                                       h = 0.5, N_bvp = 8,
                                       bvp_maxiter = 10)                  # WARNS
    end
end
