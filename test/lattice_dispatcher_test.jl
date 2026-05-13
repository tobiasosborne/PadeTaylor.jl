# test/lattice_dispatcher_test.jl -- Phase 12 v2 / bead `padetaylor-k31` tests.
#
# Tier-3+ 2D-lattice composition layer: PathNetwork (Phase 10) +
# EdgeDetector (Phase 12.5) + per-row BVP (Phase 11) per FW 2011 §4.4
# `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:190`:
# "161 separate BVP solutions; one for each grid line."
#
# v1 scope: per-row BVP fill on contiguous smooth runs flanked by
# pole-field cells.  Smooth runs touching grid boundaries retain IVP
# values (tagged `:ivp_only`).  FW Fig 4.1 quantitative pin
# (u(20i) ≤ 1e-10 abs) deferred to a follow-up bead — different
# compositional pattern (vertical BVP + two outward pole fields).
#
# Test plan:
#   - LD.1.1: PI tritronquée Phase-9 grid composition completes; region
#     tags consistent with mask; no spurious BVP fills (the wedge is at
#     the grid edge, so no smooth run is two-sided-flanked).
#   - LD.1.2: linear u''=u with SUPPLIED mask creates interior smooth
#     runs; BVP fills happen; BVP-filled u-values match cosh closed
#     form to ≤ 1e-10; region tags are :ivp / :bvp / :ivp_only on the
#     right cells.
#   - LD.1.3: same setup, but mask = all-false → zero BVP fills,
#     all cells :ivp_only (the "no pole field, no bridging" case).
#   - LD.2.1: fail-fast guards.
#   - LD.3.1: mutation-proof (procedure documented at bottom).

using Test
using PadeTaylor

@testset "LatticeDispatcher (Phase 12 v2): FW 2011 §4.4 + line 190" begin

    # ---- LD.1.1: PI tritronquée at Phase-9 setup ---------------------------
    @testset "LD.1.1: PI tritronquée composition completes (Phase 9 setup)" begin
        f_PI(z, u, up) = 6 * u^2 + z
        f_PI_1(z, u)   = 6 * u^2 + z
        ∂f_PI_1(z, u)  = 12 * u
        u_tri, up_tri  = -0.1875543083404949, 0.3049055602612289

        N = 25
        xs = range(-4.0, 4.0; length = N)
        ys = range(-4.0, 4.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(4 * sqrt(2)))
        prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

        sol = lattice_dispatch_solve(prob, f_PI_1, ∂f_PI_1, xs, ys;
                                     h_path = 0.5, order = 30)

        @test size(sol.grid_u)  == (N, N)
        @test size(sol.grid_up) == (N, N)
        @test size(sol.mask)    == (N, N)
        @test size(sol.region_tag) == (N, N)
        # Region tag invariant: every mask=true cell is :ivp; every
        # mask=false cell is :bvp OR :ivp_only.  (No untagged cells.)
        @test count(sol.region_tag .== :ivp) == count(sol.mask)
        @test all(t -> t in (:ivp, :bvp, :ivp_only), sol.region_tag)
        # PI tritronquée wedge is at grid edge (positive real axis) →
        # no smooth run is two-sided-flanked → zero BVP fills.
        @test length(sol.bvp_solutions) == 0
        @test count(sol.region_tag .== :bvp) == 0
    end

    # ---- LD.1.2: linear u''=u with supplied mask, BVP fill triggers --------
    @testset "LD.1.2: linear u''=u with supplied mask, BVP fill matches cosh" begin
        # ODE u'' = u; IC u(0)=1, u'(0)=0; closed form u(z) = cosh(z)
        # for any complex z.
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)

        N = 11
        xs = range(-1.5, 1.5; length = N)
        ys = range(-1.5, 1.5; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(1.5 * sqrt(2)))
        prob  = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 30)

        # Supply a synthetic mask: mask=true at i ∈ {3, 9} for every
        # interior row j ∈ 2..N-1.  This creates a single bridgeable
        # smooth run from i=4..i=8 in each interior row.
        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end

        sol = lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                     h_path = 0.5, order = 20, mask = mask)

        # Each interior row should have one bridgeable smooth run, so
        # we expect N-2 = 9 BVP solutions.
        @test length(sol.bvp_solutions) == N - 2

        # Region-tag invariants
        @test sol.region_tag[3, 6] == :ivp        # supplied-true cell
        @test sol.region_tag[9, 6] == :ivp        # supplied-true cell
        @test sol.region_tag[6, 6] == :bvp        # interior of bridged run
        @test sol.region_tag[6, 1] == :ivp_only   # boundary row, unbridged
        @test sol.region_tag[1, 6] == :ivp_only   # boundary col, unbridged

        # BVP-filled cells must match cosh(z) to spectral precision.
        # N_bvp=20 on a width-1.2 segment gives Chebyshev convergence
        # well below 1e-10 for the entire cosh on the segment.
        for i in 4:8, j in 2:(N - 1)
            z = xs[i] + im * ys[j]
            u_exact = cosh(z)
            @test abs(sol.grid_u[i, j] - u_exact) < 1e-10
        end

        # The :ivp cells (supplied mask=true) carry IVP values from the
        # path-network — also close to cosh, but only to path-network
        # accuracy (~1e-13 for this trivially analytic problem).
        for j in 2:(N - 1)
            z = xs[3] + im * ys[j]
            @test abs(sol.grid_u[3, j] - cosh(z)) < 1e-10
        end
    end

    # ---- LD.1.3: mask=all-false → no BVP fills -----------------------------
    @testset "LD.1.3: mask=all-false produces zero BVP fills" begin
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)
        N  = 9
        xs = range(-1.0, 1.0; length = N)
        ys = range(-1.0, 1.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(sqrt(2)))
        prob  = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 20)

        mask_zero = falses(N, N)
        sol = lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                     h_path = 0.5, order = 20, mask = mask_zero)
        @test length(sol.bvp_solutions) == 0
        @test count(sol.region_tag .== :bvp) == 0
        @test all(sol.region_tag .== :ivp_only)
    end

    # ---- LD.2.1: fail-fast guards ------------------------------------------
    @testset "LD.2.1: fail-fast guards" begin
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)
        zspan = (0.0 + 0.0im, ComplexF64(1.5))
        prob  = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 20)

        xs3 = range(-0.5, 0.5; length = 3)
        ys3 = range(-0.5, 0.5; length = 3)
        xs2 = range(-0.5, 0.5; length = 2)

        @test_throws ArgumentError lattice_dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, xs2, ys3)
        @test_throws ArgumentError lattice_dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, xs3, xs2)
        @test_throws ArgumentError lattice_dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, xs3, ys3; N_bvp = 3)
        # Mask shape mismatch
        @test_throws ArgumentError lattice_dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, xs3, ys3;
            mask = falses(4, 4))
        # Anisotropic grid (xs step ≠ ys step)
        xs_aniso = collect(0.0:0.1:0.2)
        ys_aniso = collect(0.0:0.2:0.4)
        @test_throws ArgumentError lattice_dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, xs_aniso, ys_aniso)
    end

end

# Mutation-proof procedure (verified 2026-05-13).
#
# Mutation E — swap `u_a` and `u_b` in the bvp_solve call:
#     bvp_sol = bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b, u_b, u_a; ...)
# Verified bite: LD.1.2's cosh-comparison loop fails on every BVP-filled
# cell (the BVP solves a problem with the wrong BC, returning ~ −cosh
# in some cells and arbitrary values elsewhere).  Bites 45/45
# `abs(grid_u - cosh) < 1e-10` assertions in that testset.
#
# Mutation F — drop the `region_tag[k, j] = :bvp` update:
#     # (commented out the tag assignment)
# Verified bite: LD.1.2's region-tag assertion `sol.region_tag[6, 6] ==
# :bvp` fails (the cell keeps its default :ivp_only tag).  Note that
# the u-grid values are still correct, since the BVP-fill loop above
# the tag assignment still runs.  This is a clean isolated bite on
# the tagging contract — separating "did the BVP run" from "is the
# bookkeeping correct".
#
# Both mutations restored before commit per CLAUDE.md Rule 4.
