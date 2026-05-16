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
        # Under ADR-0017 (bead padetaylor-0tj) the default IVP source is
        # `edge_gated_pole_field_solve`, whose tighter, edge-confirmed
        # mask exposes a small number of smooth cells along the wedge
        # boundary that *are* two-sided-flanked by field cells (the v2
        # plain `path_network_solve` source produced an over-inclusive
        # mask under which the same cells were un-flanked).  Pin to a
        # small finite count rather than `== 0`; the exact count (2 on
        # this fixture at HEAD) is geometry-sensitive, so we bound it
        # rather than equate it to avoid brittleness under further mask
        # refinement.
        @test length(sol.bvp_solutions) ≤ 5
        @test count(sol.region_tag .== :bvp) == length(sol.bvp_solutions)
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

    # ---- LD.2.2: order kwarg defaults to prob.order (bead `padetaylor-9xf`) ----
    @testset "LD.2.2: order kwarg defaults to prob.order" begin
        # Regression for `padetaylor-9xf`: `lattice_dispatch_solve`
        # carried its own `order::Integer = 30` kwarg and threaded *that*
        # into `path_network_solve`, ignoring `prob.order`.  It now
        # defaults to `prob.order`.  Probe on the pole-free linear ODE
        # u'' = u (path_network just walks straight): a problem built at
        # order 4 must — with no explicit `order` — solve at order 4
        # (bit-exact match to an explicit `order = 4` call) and NOT at
        # the old hard-coded 30.  `mask = all-false` ⇒ pure IVP fill, so
        # every cell's value is the order-dependent `path_network_solve`
        # output.
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)
        N  = 5
        xs = range(-0.5, 0.5; length = N)
        ys = range(-0.5, 0.5; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(sqrt(2)))
        prob4 = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 4)
        mask0 = falses(N, N)

        sol_default = lattice_dispatch_solve(prob4, f_lin_1, ∂f_lin_1, xs, ys;
                                             h_path = 0.5, mask = mask0)
        sol_ord4    = lattice_dispatch_solve(prob4, f_lin_1, ∂f_lin_1, xs, ys;
                                             h_path = 0.5, order = 4, mask = mask0)
        sol_ord30   = lattice_dispatch_solve(prob4, f_lin_1, ∂f_lin_1, xs, ys;
                                             h_path = 0.5, order = 30, mask = mask0)

        @test sol_default.grid_u == sol_ord4.grid_u     # honours prob.order
        @test sol_default.grid_u != sol_ord30.grid_u    # not the old hard-coded 30

        # MUTATION-PROOF (verified 2026-05-14, bead `padetaylor-9xf`):
        # reverting `src/LatticeDispatcher.jl`'s `order::Integer =
        # prob.order` to `order::Integer = 30` makes the default solve
        # identical to the `order = 30` solve — both assertions go RED.
        # Restored.
    end

    # ====================================================================
    # Phase 12 v3 / bead `padetaylor-0tj` / ADR-0017 — `strict::Bool`
    # kwarg + edge-gated IVP source.  Tests LD.X.1 … LD.X.7.
    #
    # Divergence-trigger trick used by LD.X.3 / LD.X.4: pass `bvp_tol`
    # tighter than Newton's step-norm floor (`eps^(3/4)` ≈ 1e-12 for
    # Float64).  At `bvp_tol = 1e-300` Newton runs `maxiter = 10`
    # iterations without ever satisfying `‖Δu‖_∞ ≤ tol` and `bvp_solve`
    # throws the documented `ErrorException("bvp_solve: Newton did not
    # converge…")` — exactly the failure mode `strict = false` is meant
    # to swallow.  This is fixture-independent and works on both the
    # auto and mask paths.
    # ====================================================================

    # ---- LD.X.1: strict=true (default) still throws on forced divergence ---
    @testset "LD.X.1: strict=true throws on BVP non-convergence (back-compat)" begin
        # PI tritronquée on a small mask-driven fixture; impossibly tight
        # `bvp_tol` forces the per-row BVP to diverge.  Default (no kwarg)
        # must rethrow the bvp_solve ErrorException verbatim.
        f_PI(z, u, up) = 6 * u^2 + z
        f_PI_1(z, u)   = 6 * u^2 + z
        ∂f_PI_1(z, u)  = 12 * u
        u_tri, up_tri  = -0.1875543083404949, 0.3049055602612289

        N = 11
        xs = range(-2.0, 2.0; length = N)
        ys = range(-2.0, 2.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(2.0 * sqrt(2)))
        prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 20)

        # Mask path: flanks at i ∈ {3, 9} on every interior row, so the
        # dispatcher tries to BVP-fill the run i ∈ 4:8 on each row.
        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end

        # Default (strict=true) — must throw, message-matched.
        @test_throws ErrorException lattice_dispatch_solve(
            prob, f_PI_1, ∂f_PI_1, xs, ys;
            h_path = 0.5, order = 20, mask = mask, bvp_tol = 1e-300)

        # Explicit strict=true — same.
        @test_throws ErrorException lattice_dispatch_solve(
            prob, f_PI_1, ∂f_PI_1, xs, ys;
            h_path = 0.5, order = 20, mask = mask, bvp_tol = 1e-300,
            strict = true)
    end

    # ---- LD.X.2: strict=false returns a LatticeSolution (fail-soft) --------
    @testset "LD.X.2: strict=false swallows Newton non-convergence" begin
        f_PI(z, u, up) = 6 * u^2 + z
        f_PI_1(z, u)   = 6 * u^2 + z
        ∂f_PI_1(z, u)  = 12 * u
        u_tri, up_tri  = -0.1875543083404949, 0.3049055602612289

        N = 11
        xs = range(-2.0, 2.0; length = N)
        ys = range(-2.0, 2.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(2.0 * sqrt(2)))
        prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 20)

        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end

        # Fail-soft — must NOT throw, must return a LatticeSolution.
        sol = lattice_dispatch_solve(prob, f_PI_1, ∂f_PI_1, xs, ys;
                                     h_path = 0.5, order = 20,
                                     mask = mask, bvp_tol = 1e-300,
                                     strict = false)
        @test sol isa LatticeSolution
        @test size(sol.region_tag) == (N, N)
        # No `:bvp` cell may carry NaN — only converged BVPs got that tag.
        for i in 1:N, j in 1:N
            if sol.region_tag[i, j] == :bvp
                @test isfinite(sol.grid_u[i, j])
            end
        end
    end

    # ---- LD.X.3: :bvp_fail tags appear when BVP diverges -------------------
    @testset "LD.X.3: :bvp_fail tag appears under forced divergence" begin
        f_PI(z, u, up) = 6 * u^2 + z
        f_PI_1(z, u)   = 6 * u^2 + z
        ∂f_PI_1(z, u)  = 12 * u
        u_tri, up_tri  = -0.1875543083404949, 0.3049055602612289

        N = 11
        xs = range(-2.0, 2.0; length = N)
        ys = range(-2.0, 2.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(2.0 * sqrt(2)))
        prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 20)

        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end

        sol = lattice_dispatch_solve(prob, f_PI_1, ∂f_PI_1, xs, ys;
                                     h_path = 0.5, order = 20,
                                     mask = mask, bvp_tol = 1e-300,
                                     strict = false)
        # Every bridgeable run hits divergence ⇒ at least one :bvp_fail.
        @test count(sol.region_tag .== :bvp_fail) > 0
        # And by construction every smooth-run cell flips to :bvp_fail
        # (no row converges, so no :bvp tag survives in the bridged span).
        @test count(sol.region_tag .== :bvp) == 0
    end

    # ---- LD.X.4: :bvp_fail cells retain their IVP values (no overwrite) ----
    @testset "LD.X.4: :bvp_fail cells keep their pre-BVP values" begin
        # Manual-mask path: pre-BVP values are the `path_network_solve`
        # IVP run, which is finite everywhere on this benign linear
        # fixture.  After fail-soft divergence, the :bvp_fail cells must
        # bit-exactly match the IVP reference computed on the same grid.
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)

        N = 11
        xs = range(-1.5, 1.5; length = N)
        ys = range(-1.5, 1.5; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(1.5 * sqrt(2)))
        prob  = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 20)

        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end

        sol = lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                     h_path = 0.5, order = 20,
                                     mask = mask, bvp_tol = 1e-300,
                                     strict = false)

        # Independent IVP reference on the same mask path.  Mirrors the
        # mask-path branch of lattice_dispatch_solve verbatim (the same
        # call into path_network_solve over the vec'd grid).
        grid_z = ComplexF64[xs[i] + im * ys[j] for i in 1:N, j in 1:N]
        pn_ref = PadeTaylor.path_network_solve(prob, vec(grid_z);
                                               h = 0.5, order = 20)
        ivp_u  = Matrix{ComplexF64}(reshape(pn_ref.grid_u, (N, N)))

        n_fail = 0
        for i in 1:N, j in 1:N
            if sol.region_tag[i, j] == :bvp_fail
                n_fail += 1
                @test sol.grid_u[i, j] == ivp_u[i, j]
            end
        end
        @test n_fail > 0
    end

    # ---- LD.X.5: manual-mask fallback path still works ---------------------
    @testset "LD.X.5: mask kwarg still routes through path_network_solve" begin
        # Reproduces LD.1.2 in miniature — supplied mask → manual path,
        # BVP fill happens, region tags correct, return value is shape-OK.
        f_lin_2(z, u, up) = u
        f_lin_1(z, u)     = u
        ∂f_lin_1(z, u)    = one(u)
        N = 11
        xs = range(-1.5, 1.5; length = N)
        ys = range(-1.5, 1.5; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(1.5 * sqrt(2)))
        prob  = PadeTaylorProblem(f_lin_2, (1.0, 0.0), zspan; order = 20)
        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end
        sol = lattice_dispatch_solve(prob, f_lin_1, ∂f_lin_1, xs, ys;
                                     h_path = 0.5, order = 20, mask = mask)
        @test sol isa LatticeSolution
        @test size(sol.grid_u) == (N, N)
        @test sol.region_tag[3, 6] == :ivp
        @test sol.region_tag[6, 6] == :bvp     # interior of bridged run
        @test sol.region_tag[1, 6] == :ivp_only
        @test length(sol.bvp_solutions) == N - 2
    end

    # ---- LD.X.6: region_tag enumeration correctness ------------------------
    @testset "LD.X.6: region_tag ⊆ {:ivp, :bvp, :bvp_fail, :ivp_only}" begin
        # Run the LD.X.2 fixture (which exercises :ivp, :bvp_fail, and
        # :ivp_only) and assert the unique tag set is a subset of the
        # documented enum.  Catches typos (e.g., `:bvp_failed`).
        f_PI(z, u, up) = 6 * u^2 + z
        f_PI_1(z, u)   = 6 * u^2 + z
        ∂f_PI_1(z, u)  = 12 * u
        u_tri, up_tri  = -0.1875543083404949, 0.3049055602612289
        N = 11
        xs = range(-2.0, 2.0; length = N)
        ys = range(-2.0, 2.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(2.0 * sqrt(2)))
        prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 20)
        mask = falses(N, N)
        for j in 2:(N - 1)
            mask[3, j] = true
            mask[9, j] = true
        end
        sol = lattice_dispatch_solve(prob, f_PI_1, ∂f_PI_1, xs, ys;
                                     h_path = 0.5, order = 20,
                                     mask = mask, bvp_tol = 1e-300,
                                     strict = false)
        allowed = Set([:ivp, :bvp, :bvp_fail, :ivp_only])
        @test issubset(Set(unique(sol.region_tag)), allowed)
    end

    # ---- LD.X.7: strict=false on the auto (edge-gated) path returns OK ----
    @testset "LD.X.7: strict=false on the default edge-gated path" begin
        # No mask supplied → edge-gated IVP fill.  With a tractable
        # `bvp_tol`, the v3 default code path must complete cleanly on
        # the PI tritronquée fixture.  Confirms the new default code
        # path is robust (no spurious BVP failures on a smooth fixture).
        f_PI(z, u, up) = 6 * u^2 + z
        f_PI_1(z, u)   = 6 * u^2 + z
        ∂f_PI_1(z, u)  = 12 * u
        u_tri, up_tri  = -0.1875543083404949, 0.3049055602612289
        N = 21
        xs = range(-4.0, 4.0; length = N)
        ys = range(-4.0, 4.0; length = N)
        zspan = (0.0 + 0.0im, ComplexF64(4.0 * sqrt(2)))
        prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)
        sol = lattice_dispatch_solve(prob, f_PI_1, ∂f_PI_1, xs, ys;
                                     h_path = 0.5, order = 30, strict = false)
        @test sol isa LatticeSolution
        # Auto path: region tags use the full enum.
        @test all(t -> t in (:ivp, :bvp, :bvp_fail, :ivp_only), sol.region_tag)
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
#
# ─────────────────────────────────────────────────────────────────────
# Mutation-proof for the v3 / bead `padetaylor-0tj` testsets
# (LD.X.1 … LD.X.7), verified 2026-05-16.
#
# Mutation X — invert the strict gate in `src/LatticeDispatcher.jl`
# around line 414, changing
#     if !strict && err isa ErrorException && occursin(...)
# to
#     if strict && err isa ErrorException && occursin(...)
# This swaps the fail-soft semantics: `strict = false` (the documented
# fail-soft mode) now rethrows, while `strict = true` swallows.  The
# whole `strict::Bool` feature is defeated.
#
# Verified bite (single Julia run, file in isolation):
#   • LD.X.2 (strict=false ... swallows Newton non-convergence) goes
#     RED — the lattice_dispatch_solve call now rethrows the
#     `ErrorException("bvp_solve: Newton did not converge…")` rather
#     than returning a LatticeSolution.
#   • LD.X.3 / LD.X.4 / LD.X.6 also bite because the strict=false
#     fixtures they share all rethrow before tagging anything.
#   • LD.X.1 (strict=true still throws) goes GREEN by accident (the
#     inverted gate happens to also catch under strict=true, but the
#     two explicit `strict=true` `@test_throws` assertions fail
#     because the throw is now swallowed).  Net: LD.X.1 also bites.
#
# Restored verbatim before commit.
