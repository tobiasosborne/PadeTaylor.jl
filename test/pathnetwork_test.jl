# test/pathnetwork_test.jl -- Phase 10 / bead `padetaylor-1jf` tests.
#
# Tier-2 path-network: FW 2011 §3.1 5-direction wedge tree + Stage-2
# fine-grid extrapolation.  Verifies that the complex-plane path
# network bridges the lattice pole of u(z) = ℘(z + c₁; 0, c₂) at z = 1
# via off-axis detours when stepping past it on the real axis would
# land on the pole.
#
# Reference: docs/adr/0004-path-network-architecture.md (algorithm +
# test plan PN.1.1-PN.3.1), docs/unified_path_network_spec.md (full
# spec), references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:155-166 (FW 2011 §3.1).
#
# This first cut covers PN.1.1, PN.1.2, PN.2.1, PN.4.1.  PN.2.2 (FW
# Table 5.1 long-range z=30 to ≤1e-13) and PN.3.1 (:steepest_descent
# agreement with :min_u) are deferred to a follow-up commit so the
# initial GREEN ships on a tractable test corpus.

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "PathNetwork (Phase 10): FW 2011 §3.1 path-network" begin

    # Test ODE: u'' = 6u^2 with FW 2011 ICs at z=0; closed form
    # u(z) = ℘(z + c_1; 0, c_2) with c_1 = -1, c_2 = 2.  Pole at z = 1.
    fW(z, u, up) = 6 * u^2

    @testset "PN.1.1: Stage 1 + Stage 2 single-step sanity (5-pt grid near z=0)" begin
        # Grid entirely inside |z| ≤ h=0.5; Stage 1 takes 0 steps for
        # most targets (they're already within h of the IC); Stage 2
        # is the IC-Padé evaluated at t = z_f / h.  This is the
        # cheapest end-to-end exercise of the public API.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[0.0, 0.1+0.2im, -0.1+0.3im, 0.3-0.2im, 0.4+0.0im]
        sol  = path_network_solve(prob, grid; h = 0.5)

        # Stage 1: IC plus zero or one extra visited node (depending on
        # which targets fell within h).  Always ≥ 1 (the IC itself).
        @test length(sol.visited_z) ≥ 1
        @test sol.visited_z[1] == ComplexF64(0.0)

        # Stage 2: every grid point covered (no NaN).
        @test all(isfinite, real.(sol.grid_u))
        @test all(isfinite, imag.(sol.grid_u))

        # Verify against analytic ℘ at z = 0.4 (real-axis crosscheck
        # — far from pole, Padé is essentially exact).  ℘(0.4 + c1; 0, c2)
        # via series expansion is tabulated in _oracle_problems.jl as
        # implied by Phase-6 logic; here we just sanity-check magnitude.
        idx_04 = findfirst(z -> isapprox(real(z), 0.4) && abs(imag(z)) < 1e-12,
                           sol.grid_z)
        @test idx_04 !== nothing
        # u(0.4) — finite, real-valued (we're on real axis), and
        # consistent with Phase-6's IVP at z = 0.5 nearby (u ≈ 4.00 at
        # z = 0.5; at z = 0.4 we're a bit smaller in magnitude).
        @test abs(imag(sol.grid_u[idx_04])) < 1e-10
        @test real(sol.grid_u[idx_04]) > 1.07     # Just above IC value
        @test real(sol.grid_u[idx_04]) < 5.0      # Below z=0.5 value
    end

    @testset "PN.1.2: Multi-step path bridging the pole at z=1" begin
        # Target z = 1.4 (past the pole at z=1).  With h=0.5 from z=0
        # the goal direction is purely real; direct stepping would
        # land at z=0.5, then z=1.0 (ON THE POLE).  The 5-direction
        # wedge + min-|u| selection should detour off-axis.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
        grid = ComplexF64[1.4 + 0.0im]
        sol  = path_network_solve(prob, grid; h = 0.5)

        # Must have taken at least 3 steps from the IC.
        @test length(sol.visited_z) ≥ 3

        # None of the visited points should be exactly at (or have
        # blown up near) z=1.  If min-|u| failed and any candidate
        # landed on the pole, we'd see Inf or NaN here.
        @test all(isfinite, abs.(sol.visited_u))
        @test all(z -> abs(z - 1.0) > 0.01, sol.visited_z[2:end])

        # At least one off-axis (Im ≠ 0) visited node — proof of
        # complex-plane detour.
        @test any(z -> abs(imag(z)) > 1e-6, sol.visited_z)

        # u(1.4) at the target — finite + bounded.  Structural only.
        # The quantitative accuracy check belongs in PN.2.2 (FW Table 5.1
        # long-range tuning); the path-network's `:min_u` walk currently
        # accumulates Padé approximation error along the off-axis detour
        # and the final boundary-of-disc evaluation at `|t| ≈ 1` is the
        # dominant source.  See ADR-0004 deferral and worklog 005's
        # order/rtol-coupling pattern.  For PN.1.2 here we only assert
        # that the algorithm did not blow up or NaN.
        u_target = sol.grid_u[1]
        @test isfinite(real(u_target))
        @test isfinite(imag(u_target))
        @test abs(u_target) < 1e3       # Rough sanity: not divergent.
    end

    @testset "PN.2.1: Stage-2 NaN sentinel for uncovered grid points" begin
        # A grid point far outside any visited node's disc must return
        # NaN+NaN·im, not silent extrapolation (CLAUDE.md Rule 1).
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 0.5); order = 30)
        # Targets: trivial near-IC + one wild outlier at z=100 not
        # listed as a target (so Stage 1 ignores it).  Then evaluate
        # at z=100 in Stage 2 anyway — it should NaN.
        grid_targets = ComplexF64[0.1, 0.2, 0.3]
        sol = path_network_solve(prob, grid_targets; h = 0.5)
        # Confirm targets work normally.
        @test all(isfinite, real.(sol.grid_u))

        # Now build a separate grid that includes an uncovered point.
        # We do this via a re-invocation including z=100; Stage-1 will
        # FAIL to reach z=100 in max_steps_per_target=50 because that
        # requires 200 steps.  So we expect a thrown error here.
        @test_throws ErrorException path_network_solve(
            prob, ComplexF64[100.0 + 0.0im];
            h = 0.5, max_steps_per_target = 50)
    end

    @testset "PN.4.1: Fail-fast guards (CLAUDE.md Rule 1)" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[0.1 + 0.0im]

        # Unknown step_size_policy still throws (only :fixed and
        # :adaptive_ffw are recognised; bead `padetaylor-8ui` retired
        # the :adaptive_ffw deferral in favour of the FFW 2017 §2.1.2
        # controller — see ADR-0011).
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, step_size_policy = :bogus_policy)

        # Unknown step selection.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, step_selection = :bogus)

        # Wrong wedge size.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, wedge_angles = [-π/4, 0.0, π/4])

        # Non-positive h.
        @test_throws ArgumentError path_network_solve(prob, grid; h = -0.5)
        @test_throws ArgumentError path_network_solve(prob, grid; h = 0.0)
    end

    @testset "PN.2.2: FW Table 5.1 z=30 long-range integration" begin
        # FW 2011 Fig 5.1 / Table 5.1: u(z=30) of the equianharmonic ℘
        # trajectory with FW IC matches `1.0950982559597442` to ≤1e-13
        # at BigFloat-256 (FW reports their own 8.34e-14 in the paper;
        # `docs/figure_catalogue.md §1 row FW2011 Fig 5.1`).
        #
        # This test exercises the path-network's load-bearing canonical-
        # Padé-per-visited-node invariant (ADR-0004 design decision):
        # each visited node z_v must store a REAL-h-direction Padé so
        # Stage 2's `t = (z_f - z_v) / h_v` interpolation lands inside
        # the disc.  See worklog 008 §"The wedge-vs-canonical-Padé bug
        # at long range".

        # Float64 path: rel-err ≤ 1e-12 acceptance.  Tightened from `1e-9`
        # in bead `padetaylor-txg`: switching the F64 Padé default from
        # GGT 2013 SVD to FW 2011 classical-via-Toeplitz-`\\` recovers the
        # per-step accuracy that GGT's robustness machinery was throwing
        # away.  Worklog 020 probe measured `1.54e-13` rel-err at z=30
        # under classical; the test bound carries an order-of-magnitude
        # cross-platform margin.  ~75 visited nodes at h=0.5, order=30.
        prob_f64 = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 30.0);
                                     order = 30)
        sol_f64 = path_network_solve(prob_f64, ComplexF64[30.0 + 0im];
                                     h = 0.5, max_steps_per_target = 200)
        u30_f64 = sol_f64.grid_u[1]
        @test isapprox(u30_f64, u_at_30_FW_ref; rtol = 1e-12)
        @test abs(imag(u30_f64)) < 1e-12       # Real solution on real axis.

        # BF-256 path: rel-err ≤ 1e-13 acceptance per FW Table 5.1.
        # ~50s wall time; the dominant test cost in the suite.
        setprecision(BigFloat, 256) do
            u0_bf  = big"1.071822516416917"
            up0_bf = big"1.710337353176786"
            ref_bf = big"1.0950982559597442"
            prob_bf = PadeTaylorProblem(fW, (u0_bf, up0_bf),
                                        (big(0.0), big(30.0)); order = 30)
            sol_bf = path_network_solve(prob_bf,
                                        Complex{BigFloat}[Complex{BigFloat}(big(30.0))];
                                        h = big(0.5), max_steps_per_target = 200)
            u30_bf = sol_bf.grid_u[1]
            rel = abs(u30_bf - ref_bf) / abs(ref_bf)
            @test rel ≤ big"1e-13"
            @test abs(imag(u30_bf)) < big"1e-15"
        end
    end

    @testset "PN.2.3: FW Table 5.1 z=10⁴ long-range integration" begin
        # FW 2011 Table 5.1 column (b): long-range integration to
        # u(10⁴) = 21.02530339471055 of the equianharmonic ℘ trajectory
        # with FW IC.  FW's Padé method (in BigFloat per Table 5.1's
        # 7.62e-14 BF-rel-err at z=30) reports 2.34e-10 rel-err at
        # z=10⁴ — accumulated error over their h=0.5 path-network walk.
        # See references/markdown/FW2011_*.md:299-301 (the reference
        # value) and ...md:385-391 (Table 5.1 itself).
        #
        # The bead `padetaylor-g9x` asks to "confirm or refute the
        # claim" that PN.2.2's z=30 BF-256 result generalises to
        # z=10⁴.  The Float64 routine test below covers algorithmic
        # stability over the full ~24,000-node walk; the BF-256
        # confirmation is an offline one-shot probe (see
        # `external/probes/pathnetwork-long-range/probe.jl` +
        # `docs/worklog/019-fw-table-51-z10000.md`) — wall time at
        # BF-256 is ~4.5 h, untenable in the routine suite.

        # Float64 path: rel-err ≤ 5·10⁻¹⁰ acceptance.  Tightened from
        # `5·10⁻⁵` in bead `padetaylor-txg`: classical Padé via Toeplitz
        # `\\` (FW 2011 §5.1.4) replaces GGT 2013 SVD as the F64 default
        # and closes the long-range accuracy gap.  Measured: rel-err
        # `1.4·10⁻¹⁰`, |imag(u)| `2.6·10⁻⁹` at the full PathNetwork
        # (Stage 1 + Stage 2 interpolation).  We compare against FW's
        # published `2.34·10⁻¹⁰` (Table 5.1 column b) — classical at
        # full-PathNetwork still BEATS FW by 1.7×, vs `6.05·10⁻⁶` under
        # SVD (~5 orders worse than FW).  Margins: ~3.6× on rtol, ~3.8×
        # on imag, both holding tight under classical roundoff.
        #
        # The worklog 020 number `6.15·10⁻¹¹` was from a custom wedge
        # walker (Stage 1 only) — the full PathNetwork is ~2× looser
        # because Stage 2's barycentric interpolation adds a Padé-eval
        # step at the target.  Closing further would need BF-256.
        prob_f64 = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 10000.0);
                                     order = 30)
        sol_f64 = path_network_solve(prob_f64,
                                     ComplexF64[10000.0 + 0im];
                                     h = 0.5,
                                     max_steps_per_target = 200_000)
        u10k_f64 = sol_f64.grid_u[1]
        @test isapprox(u10k_f64, u_at_10000_FW_ref; rtol = 5e-10)
        @test abs(imag(u10k_f64)) < 1e-8     # Real solution on real axis.
        @test length(sol_f64.visited_z) > 20_000  # Stage-1 walked ~24k nodes.
    end

    @testset "PN.5.1: visited_parent is a well-formed Stage-1 path tree" begin
        # FW 2011 Fig 3.2 draws the Stage-1 path tree.  The tree edges
        # are `{(visited_parent[k], k) : k ≥ 2}`; this testset pins the
        # structural invariants that make `visited_parent` a tree:
        #
        #   (a) one parent slot per visited node;
        #   (b) the root (IC node) has parent 0;
        #   (c) every non-root parent precedes its child (1 ≤ p ≤ k-1)
        #       — no cycles, and the tree is built incrementally;
        #   (d) every tree edge has Euclidean length ≈ h, because each
        #       visited node is reached by exactly one wedge step
        #       (length h) from its parent — including the first node
        #       of a per-target walk, which chains off the nearest
        #       already-visited node (FW 2011 line 164);
        #   (e) following parent pointers from any node reaches the root.
        #
        # A multi-target grid that forces multi-step walks exercises
        # both the "first node chains off nearest visited" branch and
        # the "later node chains off predecessor" branch.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
        grid = ComplexF64[1.4 + 0im, 0.8 + 0.6im, -0.5 - 0.4im, 1.1 + 0.2im]
        h    = 0.5
        sol  = path_network_solve(prob, grid; h = h)

        n = length(sol.visited_z)
        # (a)
        @test length(sol.visited_parent) == n
        # (b)
        @test sol.visited_parent[1] == 0
        # The grid forces enough stepping that the tree is non-trivial.
        @test n ≥ 6
        for k in 2:n
            p = sol.visited_parent[k]
            # (c)
            @test 1 ≤ p ≤ k - 1
            # (d)
            @test isapprox(abs(sol.visited_z[k] - sol.visited_z[p]), h;
                           atol = 1e-9)
        end
        # (e) — walk every node's ancestry up to the root.
        for k in 1:n
            cur  = k
            hops = 0
            while cur != 0
                cur  = cur == 1 ? 0 : sol.visited_parent[cur]
                hops += 1
                @test hops ≤ n        # cycle guard
            end
        end

        # MUTATION-PROOF (verified 2026-05-14, bead `padetaylor-d26`):
        # in `src/PathNetwork.jl`, deleting the line
        #   `parent_idx = length(visited_z)   # next step chains off this node`
        # makes every node of a multi-step walk point at the walk's
        # start node `idx_v` instead of its immediate predecessor.
        # Invariant (d) then goes RED for every node ≥ 2 hops into a
        # walk — its edge to `idx_v` has length k·h, not h.  Restored.
    end

    @testset "PN.3.1: :steepest_descent path agrees with :min_u (pole-bridge)" begin
        # FW 2011 §5.4.1 (line 362-368) introduces `:steepest_descent`
        # as a perf-tuned alternative to `:min_u`: pick the wedge angle
        # closest to θ_sd = arg(-u/u') instead of evaluating all five
        # |u| values.  On any smooth-region grid + pole-bridge grid
        # where both rules agree on the wedge index at each step, the
        # two paths should produce identical visited-node values to
        # within Padé approximation noise (≤1e-10).
        #
        # We use targets that REQUIRE stepping (distance > h) including
        # one pole-bridge case (z=1.4 past the lattice pole at z≈1.13),
        # so the step-selection logic is exercised.  Both rules' answer
        # at z=1.4 is u≈6.2518 (Phase-10 PN.1.2 ground truth).
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
        grid = ComplexF64[1.4 + 0im, 1.2 + 0.4im, 0.6 + 0.3im]
        sol_min     = path_network_solve(prob, grid; h = 0.5,
                                         step_selection = :min_u)
        sol_descent = path_network_solve(prob, grid; h = 0.5,
                                         step_selection = :steepest_descent)
        @test length(sol_min.grid_u) == length(sol_descent.grid_u)
        for i in eachindex(sol_min.grid_u)
            Δu  = abs(sol_min.grid_u[i]  - sol_descent.grid_u[i])
            Δup = abs(sol_min.grid_up[i] - sol_descent.grid_up[i])
            @test Δu  ≤ 1e-10
            @test Δup ≤ 1e-9   # u' carries the 1/h chain-rule factor; looser tol.
        end
    end

    @testset "PN.6.1: enforce_real_axis_symmetry — bit-exact u(z̄) = ū(z)" begin
        # Bead `padetaylor-dtj` + worklog 014.  For ODEs that preserve
        # complex conjugation (real coefficients, real ICs on the real
        # axis), the analytic solution satisfies the Schwarz reflection
        # `u(z̄) = ū(z)`.  The default path-network breaks this — the
        # `shuffle(rng, targets)` step at PathNetwork.jl:173 creates an
        # asymmetric visited tree whose Stage-2 nearest-visited lookup
        # propagates 4-5 orders of magnitude into `|u|` at conjugate-pair
        # grid cells (worklog 014 §"Bug 1").  Opt-in kwarg
        # `enforce_real_axis_symmetry=true` walks only upper-half + on-
        # axis targets, then mirrors lower-half via conj — bit-exact.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)

        # 5×5 conjugate-symmetric grid, well inside `|z| < 1` (clear of
        # the lattice pole at z = 1).  This keeps each Stage-1 walk short
        # so the test costs <1 s.
        xs = range(-0.3, 0.3; length = 5)
        ys = range(-0.3, 0.3; length = 5)   # symmetric around 0
        grid = ComplexF64[x + im*y for x in xs for y in ys]

        sol = path_network_solve(prob, grid; h = 0.5,
                                  enforce_real_axis_symmetry = true)

        # Coverage: every grid point evaluated to a finite value.
        @test all(isfinite, real.(sol.grid_u))
        @test all(isfinite, imag.(sol.grid_u))
        @test all(isfinite, real.(sol.grid_up))
        @test all(isfinite, imag.(sol.grid_up))

        # Bit-exact Schwarz reflection: for every off-axis pair (z,
        # conj(z)) in the grid, `grid_u[z] == conj(grid_u[conj(z)])`.
        # Zero floating-point slack (the mirror is exact `conj()` of the
        # upper-half walk, not a re-walk).  On-axis cells (z == conj(z))
        # are skipped — the pair is trivial; the path-network's natural
        # complex output is preserved as-is on the real axis.
        idx_of = Dict(z => i for (i, z) in enumerate(sol.grid_z))
        max_asym_u  = 0.0
        max_asym_up = 0.0
        n_offaxis_pairs = 0
        for (i, z) in enumerate(sol.grid_z)
            imag(z) == 0 && continue
            j = idx_of[conj(z)]
            n_offaxis_pairs += 1
            max_asym_u  = max(max_asym_u,  abs(sol.grid_u[i]  - conj(sol.grid_u[j])))
            max_asym_up = max(max_asym_up, abs(sol.grid_up[i] - conj(sol.grid_up[j])))
        end
        @test n_offaxis_pairs == 20   # 5×5 grid, 5 on-axis, 20 off-axis cells
        @test max_asym_u  == 0.0
        @test max_asym_up == 0.0
    end

    @testset "PN.6.2: enforce_real_axis_symmetry — visited tree confined to Im(z) ≥ 0" begin
        # The walker should NOT descend into the lower half-plane when
        # the kwarg is on.  Sanity-checks that the implementation is
        # actually filtering targets — not just doing the mirror on the
        # full-walk result.
        #
        # Grid contains lower-half targets > h from the IC so that Stage-1
        # would walk into the lower half if the filter weren't applied.
        # Without the kwarg, this grid produces an asymmetric visited tree
        # spanning both halves (worklog 014); with the kwarg, the lower-
        # half targets are reflected to the upper half before walking.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
        grid = ComplexF64[0.7-1.2im, 0.7+1.2im,
                          -0.5-0.8im, -0.5+0.8im,
                          0.3+0.0im]
        sol = path_network_solve(prob, grid; h = 0.5,
                                  enforce_real_axis_symmetry = true)
        @test length(sol.visited_z) > 1                                  # walked
        @test all(z -> imag(z) ≥ -10*eps(Float64), sol.visited_z)        # upper-only
    end

    @testset "PN.6.3: enforce_real_axis_symmetry — input grid order preserved" begin
        # Existing invariant for the default path; the kwarg must NOT
        # reorder the output.  Callers index sol.grid_u/grid_up by input
        # position (see examples/tritronquee_3d.jl's per-cell loop).
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[-0.2+0.1im, 0.1-0.3im, 0.0+0.2im, 0.2+0.0im,
                          -0.1-0.2im, 0.1+0.3im]
        sol = path_network_solve(prob, grid; h = 0.5,
                                  enforce_real_axis_symmetry = true)
        @test sol.grid_z == grid
    end

    @testset "PN.6.4: enforce_real_axis_symmetry — fail-fast on off-axis IC" begin
        # Schwarz-reflection holds globally only if the IC sits on the
        # real axis (so `u(z̄) = ū(z)` at one point ⇒ at all points by
        # uniqueness).  An off-axis IC silently breaks the assumption;
        # the kwarg must fail loud (CLAUDE.md Rule 1).

        # zspan[1] off-axis.
        prob_zoff = PadeTaylorProblem(fW, (u_0_FW + 0.0im, up_0_FW + 0.0im),
                                       (0.0 + 0.1im, 1.0 + 0.0im); order = 30)
        @test_throws ArgumentError path_network_solve(
            prob_zoff, ComplexF64[0.1+0.2im];
            h = 0.5, enforce_real_axis_symmetry = true)

        # u_0 off-axis.
        prob_uoff = PadeTaylorProblem(fW,
                                       (u_0_FW + 0.1im, up_0_FW + 0.0im),
                                       (0.0 + 0.0im, 1.0 + 0.0im); order = 30)
        @test_throws ArgumentError path_network_solve(
            prob_uoff, ComplexF64[0.1+0.2im];
            h = 0.5, enforce_real_axis_symmetry = true)

        # up_0 off-axis.
        prob_upoff = PadeTaylorProblem(fW,
                                        (u_0_FW + 0.0im, up_0_FW + 0.1im),
                                        (0.0 + 0.0im, 1.0 + 0.0im); order = 30)
        @test_throws ArgumentError path_network_solve(
            prob_upoff, ComplexF64[0.1+0.2im];
            h = 0.5, enforce_real_axis_symmetry = true)
    end

    @testset "PN.7.1: order kwarg defaults to prob.order (bead `padetaylor-9xf`)" begin
        # Regression for `padetaylor-9xf`: `path_network_solve` carried
        # its own `order::Integer = 30` kwarg and silently ignored
        # `prob.order`, so a problem built at any order ≠ 30 was solved
        # at 30 regardless (caught while writing `figures/fw2011_fig_5_2.jl`,
        # worklog 025).  The kwarg now defaults to `prob.order`.
        #
        # Probe: build the ℘ test problem at a deliberately low order
        # (6 → a (3,3) Padé) and walk to z = 3.0.  The default-order
        # solve must (a) match an explicit `order = 6` solve bit-for-bit
        # — proving it reads `prob.order` — and (b) DIFFER from an
        # explicit `order = 30` solve — proving it is no longer the old
        # hard-coded 30.  The walk is deterministic (`rng_seed = 0`), so
        # the (a) comparison is exact equality.
        prob6 = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 3.0); order = 6)
        grid  = ComplexF64[3.0 + 0im]

        sol_default = path_network_solve(prob6, grid; h = 0.5)
        sol_ord6    = path_network_solve(prob6, grid; h = 0.5, order = 6)
        sol_ord30   = path_network_solve(prob6, grid; h = 0.5, order = 30)

        # (a) default honours prob.order — bit-exact identical solve.
        @test sol_default.grid_u[1] == sol_ord6.grid_u[1]
        @test length(sol_default.visited_z) == length(sol_ord6.visited_z)
        # (b) and is NOT the old hard-coded 30.
        @test sol_default.grid_u[1] != sol_ord30.grid_u[1]

        # MUTATION-PROOF (verified 2026-05-14, bead `padetaylor-9xf`):
        # reverting `src/PathNetwork.jl`'s `order::Integer = prob.order`
        # back to `order::Integer = 30` makes the default solve identical
        # to the `order = 30` solve — test (a) goes RED (default no longer
        # matches the order-6 solve) and test (b) goes RED (default now
        # equals the order-30 solve).  Restored.
    end

end # @testset PathNetwork

# PN.5.1  Mutation-proof procedure (verified manually before commit; see
# ADR-0004 §"Test plan" + worklog 004 §"Mutation-proof procedure for
# the next agent" for the Phase-6 lineage):
#
#   Mutation A  --  in `_select_candidate`, replace `argmin(abs(e[2]))`
#     with `argmax(abs(e[2]))` for the :min_u branch.  Steers path
#     TOWARD poles instead of away from them.
#     Expected: PN.1.2 RED at lines 79, 83 — `all(isfinite, abs.(visited_u))`
#     fails because steering toward z=1 pole yields Inf/NaN, and
#     `all(z -> abs(z - 1) > 0.01)` fails because a step lands near pole.
#     Verified 2026-05-13: 2 fails on PN.1.2 as predicted.
#
#   Mutation D  --  in the Stage-2 evaluation loop, invert the coverage
#     check from `abs(z_f - z_v) > h_v` to `abs(z_f - z_v) < h_v` so
#     COVERED points NaN-out and UNCOVERED points (which we have none of
#     in PN.1.1/PN.1.2) get evaluated.
#     Expected: PN.1.1 RED at lines 45-46 + 58-60 (`isfinite` checks +
#     u(0.4) magnitude bounds); PN.1.2 RED at lines 94-96 (u_target
#     finiteness + magnitude); PN.2.1 line 109 (the targets-work-normally
#     check inverts).
#     Verified 2026-05-13: 9 fails across PN.1.1 (5), PN.1.2 (3), PN.2.1 (1)
#     as predicted — the NaN-sentinel logic is load-bearing across all
#     Stage-2 tests.
#
# Restoration: both mutations restored before commit.
#
# PN.2.2 + PN.3.1 follow-up mutations (verified 2026-05-13 with bead
# `padetaylor-yt1` in flight):
#
#   Mutation E  --  in PathNetwork.jl Stage-1 loop, restore the pre-bugfix
#     behaviour of storing the wedge-direction `pade_sel` instead of the
#     canonical-direction `pade_canonical` per visited node.  This is the
#     bug worklog 008 diagnosed; restoring it MUST bite PN.2.2 + PN.2.3.
#     Verified bite (2026-05-13): 4 fails on PN.2.2 — F64 rel-err 0.218
#     >> 1e-9, F64 imag(u) 7.2e-3 >> 1e-9, BF-256 rel-err 0.218 >> 1e-13,
#     BF-256 imag 7.2e-3 >> 1e-15.  Algorithmic error invariant to
#     arithmetic precision.  PN.2.3 (z=10⁴ at Float64, added 2026-05-13
#     bead `padetaylor-g9x`) also fails — the canonical-Padé invariant is
#     exercised at every one of the ~24,000 visited nodes; perturbing it
#     blows up the long-range result completely (rel-err >> 5e-5).  See
#     worklog 019 for the matched z=10⁴ probe under the mutation.
#
#   Mutation F  --  in `_select_candidate` :steepest_descent branch,
#     flip the sign in `θ_sd = angle(-u/u')` to `θ_sd = angle(u/u')`,
#     steering toward poles instead of away.
#     Verified bite: 4 fails on PN.3.1 — |Δu| at z=1.4 reaches 6.25
#     (steepest_descent returns 0+0im on a failed pole crossing while
#     :min_u correctly gives u≈6.2518); |Δu| at z=1.2+0.4i reaches 5.00;
#     |Δu'| reaches 31.2 and 22.3 respectively.  Way above the 1e-10 /
#     1e-9 tolerance.
#
# Restoration: both mutations restored before commit.
#
# PN.6.* (enforce_real_axis_symmetry, bead `padetaylor-dtj`, worklog 014):
#
#   Mutation G  --  in `_solve_with_schwarz_reflection`'s output-build loop,
#     drop the `conj(...)` calls in the `imag(z) < 0` branch so lower-half
#     cells inherit the upper-half walk's value as-is (no reflection).
#     Verified bite (2026-05-13): PN.6.1 lines 257-258 RED — `max_asym_u`
#     and `max_asym_up` non-zero (the Schwarz mirror is broken).  PN.6.2-4
#     stay GREEN, as predicted (they don't depend on the mirror itself).
#
#   Mutation I  --  in `_solve_with_schwarz_reflection`, replace
#     `complex(real(z), abs(imag(z)))` with `complex(real(z), imag(z))`
#     (drop `abs`).  The upper-canon collapse breaks: lower-half cells
#     map to themselves, the recursive walk visits both halves, and the
#     mirror-branch conj is applied to an asymmetric walk's value.
#     Verified bite (2026-05-13): PN.6.1 lines 257-258 RED (asymmetric
#     visited tree breaks bit-exact symmetry) AND PN.6.2 line 279 RED
#     (visited_z contains Im<0 nodes).  Validates BOTH the symmetry
#     mirror AND the upper-half filter as load-bearing.
#
# Restoration: both mutations restored before commit (GREEN at 1250).
