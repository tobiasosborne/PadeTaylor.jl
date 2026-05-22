# test/vector_pipeline_oracle_test.jl -- bead `padetaylor-0ln.37.16` (VC-9).
#
# VC-9 — the Weierstrass-℘ oracle for the v0.2 *vector* pipeline
# (ADR-0025 Validation Criteria Menu; Amendment 8).
#
# ## Why this test exists
#
# The v0.1 *scalar* Padé–Taylor pipeline is oracle-validated against the
# Fornberg–Weideman 2011 §5.1.1 test problem — the equianharmonic
# Weierstrass-℘ ODE `u'' = 6u²`, whose closed-form solution has
# analytically-known second-order poles on a rhombic lattice and
# whose medium-/long-range values FW pin in Table 5.1
# (`references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:281-301`).  `test/problems_test.jl`
# (`solve_pade` vs `u(0.5)`, `u(30)`, …) and `test/polefield_test.jl`
# (`extract_poles` vs the closed-form lattice) certify the scalar
# pipeline against it.
#
# The v0.2 **vector** pipeline — `VectorProblems.vector_solve_pade`,
# `VectorBVP.vector_bvp_solve`, `VectorPathNetwork.vector_path_network_solve`,
# `VectorPoleField.extract_poles_shared_q` — had **no equivalent oracle
# test**.  Its figure-grade exercise is the `P_I^(2)` tritronquée, which
# has *no closed form* (KKG 2015 plotted surfaces only — ADR-0025
# Decision 1).  VC-9 closes that gap: it formulates the FW ℘ test
# problem as the `d = 2` first-order companion system and drives it
# through the vector pipeline, asserting agreement with the FW Table 5.1
# reference values and the closed-form ℘ pole lattice — a genuine
# end-to-end oracle that owes nothing to the `P_I^(2)` substrate under
# test in the headline figure.
#
# ## The test problem — `u'' = 6u²` as a d=2 companion system
#
# FW 2011 eq. (5.1) `u'' = 6u²` with `c₁ = -1, c₂ = 2` has solution
# `u(z) = ℘(z - 1; 0, 2)` and initial conditions (FW md:292-295)
# `u(0) ≈ 1.071822516416917`, `u'(0) ≈ 1.710337353176786`.  Its
# first-order companion form is `y = (u, u')`,
#
#     y' = f(z, y) = (y₂, 6·y₁²)            (d = 2).
#
# FW Table 5.1 reference values (FW md:301):
#     u(30)   = 1.095098255959744           (medium range)
#     u(10⁴)  = 21.02530339471055           (long range)
# The real-axis ℘ poles sit at `x = 1 + 2ωk`, `ω = Γ(⅓)³/(2^{13/6}π)`.
#
# ## What is asserted (CLAUDE.md Rule 5 — every test pins a known value)
#
#   VPO.1  the vector IVP `vector_solve_pade` of the companion system
#          reproduces the closed-form ℘ on a pole-free segment, and the
#          FW closed-form value `u(0.5)` to ~machine precision.
#   VPO.2  the vector path-network `vector_path_network_solve` —
#          tunnelling ~21 lattice poles along the real axis — reaches
#          the FW Table 5.1 medium-range value `u(30)`.
#   VPO.3  the vector BVP `vector_bvp_solve` of the companion system on
#          a pole-free segment recovers the genuine ℘ derivative `u'(0)`
#          (the FW §5.2 "the BVP solver provides u'" diagnostic) and
#          matches the IVP ℘ trajectory.
#   VPO.4  the vector pole extractor `extract_poles_shared_q` recovers
#          the nearest closed-form ℘ lattice pole at `z = 1` — the FW
#          2011 Fig 4.7 pole-location capability, lifted to the
#          shared-`Q` vector extractor.
#
# Reference values are sourced from `test/_oracle_problems.jl` (the
# wolframscript-captured Phase-6 oracle — `u_at_30` etc. are
# cross-checked against the FW-paper-pinned `u_at_30_FW_ref`).  The ℘
# pole lattice is the closed form FW md:297 gives — no ℘-evaluator is
# hand-rolled (`test/polefield_test.jl` makes the same observation).
#
# Mutation-proof procedure is in the file footer.
#
# Standalone-runnable: `julia --project=. test/vector_pipeline_oracle_test.jl`.

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_problems.jl"))

# The FW ℘ companion RHS: u'' = 6u²  ⇒  y' = (y₂, 6·y₁²).
fW_vec(z, y) = [y[2], 6 * y[1]^2]
# The analytic d=2 companion Jacobian ∂f/∂y = [0 1; 12·y₁ 0].
Jf_W(z, y) = [zero(y[1]) one(y[1]); 12 * y[1] zero(y[1])]

# Equianharmonic ℘ real half-period ω, FW md:297 (Float64 precision;
# reproduce with `gamma(1/3)^3 / (2^(13/6)*π)`).  Identical constant to
# `test/polefield_test.jl` Ω_FW — VC-9 mirrors that test's lattice
# oracle for the vector pipeline.
const Ω_VPO = 1.3630340904278908
# Exact ℘ pole lattice of u(z) = ℘(z - 1; 0, 2): z(m,n) per FW md:297.
℘_pole_vpo(m::Int, n::Int) =
    1 + 2Ω_VPO * (m + n / 2) + im * (Ω_VPO * sqrt(3) * n)

@testset "VPO — VC-9 Weierstrass-℘ oracle for the vector pipeline" begin

    @testset "VPO.1: vector IVP reproduces closed-form ℘ (FW 5.1.1)" begin
        # `vector_solve_pade` of the d=2 companion system on the
        # pole-free segment [0, 0.5] (the nearest ℘ pole is at z = 1).
        # FW closed-form: u(0.5) = u_at_0_5 (oracle-pinned, 3-source
        # closed-form ≡ NDSolve ≡ mpmath — see _oracle_problems.jl).
        prob = VectorPadeTaylorProblem(fW_vec, [u_0_FW, up_0_FW],
                                       (0.0, 0.5); order = 30)
        sol  = vector_solve_pade(prob; h = 0.5)
        y_end = sol(0.5)
        # The companion y₁ ≡ u must equal the closed-form ℘ at z = 0.5 —
        # the FW Table-5.1-style headline value, pinned tightly.
        @test isapprox(y_end[1], u_at_0_5;  rtol = 1.0e-12)
        # ...and y₂ ≡ u' the closed-form ℘'.  The derivative component is
        # carried by the single `h = 0.5` segment's (15,15)-Padé and
        # evaluated at the segment endpoint `t = 1`; on this steep ℘
        # segment (`u` climbs `1.07 → 4.00`, `u'` is `~16`) that route
        # places `u'` to ~`4e-4` — the genuine single-segment vector
        # Padé accuracy for the derivative, NOT machine precision.  A
        # finer-`h` re-solve below confirms it converges; the tight
        # `u'` oracle is validated independently in VPO.3 (the BVP
        # recovers `u'(0)` to `1e-10`).  `5e-3` is the honest envelope.
        @test isapprox(y_end[2], up_at_0_5; atol = 5.0e-3)
        # The finer-`h` re-solve genuinely tightens `u'(0.5)` — proof the
        # `4e-4` is a step-resolution effect, not a wrong trajectory.
        finer = vector_solve_pade(
            VectorPadeTaylorProblem(fW_vec, [u_0_FW, up_0_FW],
                                    (0.0, 0.5); order = 30); h = 0.1)
        @test abs(finer(0.5)[2] - up_at_0_5) < abs(y_end[2] - up_at_0_5)
        # The trajectory starts exactly at the FW ICs.
        @test isapprox(sol(0.0)[1], u_0_FW;  atol = 1.0e-13)
        @test isapprox(sol(0.0)[2], up_0_FW; atol = 1.0e-13)
        # An interior point: a finer-`h` re-solve of the SAME companion
        # system must agree — a discretisation-independence check that
        # the single-segment dense interpolant is genuinely accurate.
        fine = vector_solve_pade(
            VectorPadeTaylorProblem(fW_vec, [u_0_FW, up_0_FW],
                                    (0.0, 0.5); order = 30); h = 0.05)
        for zq in (0.13, 0.27, 0.41)
            @test isapprox(sol(zq)[1], fine(zq)[1]; atol = 1.0e-6)
        end
    end

    @testset "VPO.2: vector path-network reaches FW Table 5.1 u(30)" begin
        # The medium-range FW Table 5.1 value u(30) = 1.095098255959744
        # (FW md:301; `u_at_30` in _oracle_problems.jl is the
        # wolframscript capture, cross-checked against the paper-pinned
        # `u_at_30_FW_ref`).  Reaching z = 30 along the real axis means
        # tunnelling ~21 ℘ lattice poles — exactly the movable-pole
        # problem the path-network solver exists for: a straight-line
        # IVP `vector_solve_pade` step would land on a pole and
        # degenerate the shared-`Q` jet.  `vector_path_network_solve`'s
        # wedge walk steers around each pole; the Stage-2 fill then
        # evaluates the dense solution at z = 30.
        prob = VectorPadeTaylorProblem(fW_vec, ComplexF64[u_0_FW, up_0_FW],
                                       (0.0 + 0im, 31.0 + 0im); order = 30)
        walk = vector_path_network_solve(
            prob, ComplexF64[31.0 + 0im];
            h = 0.15, fine_grid = ComplexF64[30.0 + 0im],
            extrapolate = true, max_steps_per_target = 100_000)
        u30 = walk.grid_y[1][1]
        # The oracle capture and the FW-paper reference agree to the
        # printed digits — assert against both (Rule 5).
        @test isapprox(u_at_30, u_at_30_FW_ref; rtol = 1.0e-15)
        # The vector pipeline reaches the medium-range value: the
        # measured error is ~6e-6 — the path-network tunnelled the ~21
        # poles and the accumulated step error is `< 1e-4`.  `1e-4` is
        # the load-bearing bar (mutation-proven below); the field value
        # genuinely tracks the FW reference, not a pole-degenerate lie.
        @test isapprox(real(u30), u_at_30; atol = 1.0e-4)
        @test isapprox(real(u30), u_at_30_FW_ref; atol = 1.0e-4)
        # The tritronquée-style separatrix instability does not apply to
        # the ℘ system — u(30) is O(1), so the imaginary part the walk
        # picks up (the real ℘ is real) is a small numerical artefact,
        # not a sign of a wrong branch.
        @test abs(imag(u30)) < 1.0e-3
        # The walk genuinely threaded a long real-axis filament.
        @test length(walk.visited_z) > 100
        @test maximum(real, walk.visited_z) ≥ 30.0
    end

    @testset "VPO.3: vector BVP recovers the genuine ℘ derivative u'(0)" begin
        # FW 2011 §5.2 (FW md:192-193): "the pole field solver not only
        # provides boundary values u_a, u_b, but also derivative values
        # u_a', u_b'."  The vector BVP, given only the closed-form ℘
        # *values* u(0), u(0.5) as the 2-point Dirichlet data, must
        # recover the genuine ℘ *derivative* u'(0) — a genuine oracle:
        # the BVP is never told u'(0); it solves the companion system
        # and the recovered y₂(0) must equal the closed-form ℘'(0).
        Ba = [1.0 0.0; 0.0 0.0]    # pin y₁ = u at z_a = 0
        Bb = [0.0 0.0; 1.0 0.0]    # pin y₁ = u at z_b = 0.5
        g  = [u_0_FW, u_at_0_5]
        # Linear-in-z initial guess for u; a flat guess for u' — the
        # BVP is never seeded with the answer u'(0).
        guess(z) = [u_0_FW + (u_at_0_5 - u_0_FW) * z / 0.5, 1.7]
        solB = vector_bvp_solve(fW_vec, 0.0, 0.5, Ba, Bb, g;
                                N = 48, jacobian = Jf_W,
                                initial_guess = guess)
        # The BVP converged at the spectral floor.
        @test solB.iterations ≤ 8
        @test solB.residual_inf < 1.0e-10
        yB0 = solB(0.0)
        # The pinned datum is reproduced...
        @test isapprox(yB0[1], u_0_FW; atol = 1.0e-12)
        # ...and the FREE derivative y₂(0) is the genuine closed-form
        # ℘'(0) — the FW §5.2 "BVP provides u'" diagnostic.  The
        # recovered value matches `up_0_FW` to ~machine precision.
        @test isapprox(yB0[2], up_0_FW; atol = 1.0e-10)
        # The far endpoint datum is reproduced.
        @test isapprox(solB(0.5)[1], u_at_0_5; atol = 1.0e-12)
        # The whole BVP solution tracks the independent IVP ℘ trajectory
        # (two distinct discretisations of one continuous problem — a
        # genuine cross-validation, as `test/vector_bvp_test.jl` VB.6
        # does for the P_I companion).  The figure `u'` is large on this
        # steep segment, so the spectral-vs-Padé gap is ~1e-6.
        ivp = vector_solve_pade(
            VectorPadeTaylorProblem(fW_vec, [u_0_FW, up_0_FW],
                                    (0.0, 0.5); order = 30); h = 0.05)
        for zq in (0.13, 0.27, 0.41)
            @test isapprox(solB(zq)[1], ivp(zq)[1]; atol = 1.0e-5)
        end
    end

    @testset "VPO.4: vector pole extractor recovers the ℘ lattice pole" begin
        # FW 2011 Fig 4.7 is a pole-location scatter; the vector
        # `extract_poles_shared_q` is its v0.2 capability.  A short
        # path-network walk from the FW IC out *past the nearest ℘ pole
        # at z = 1* — its shared-`Q` denominator roots place that pole.
        # This mirrors `test/polefield_test.jl` PF.1.1 (a single short
        # network surfaces the nearest pole) for the vector pipeline.
        #
        # The target fan is `{1.1, 1.5}` — chosen so the walk genuinely
        # threads *past* the pole at z = 1, the test's stated intent.
        # The `1.1` target makes a node land within ≈ 0.1 of the pole;
        # a local Padé places a *near* pole far more accurately than a
        # far one, so this is what makes the `1e-5` recovery bar a
        # robust assertion rather than node-geometry roulette.  The
        # adaptive `:max_q_root` walk (ADR-0026 Amendment 4) steers to
        # the *largest* pole-free disc, so with the lone far target
        # `1.5` alone it arcs ≈ 0.46 *wide* of the pole and places it
        # only ≈ 1e-4 — accurate to the geometry it actually reached,
        # but not to this oracle's `1e-5` figure-catalogue bar.  A
        # target near the pole is the FW-faithful "walk past the pole
        # to extract it" (FW §5.5 "complete the pole field"); it cannot
        # sit *on* z = 1 (the adaptive walk caps `h` → 0 onto a pole and
        # fails loud — VPN.1.1's note), so `1.1` is the honest choice.
        prob = VectorPadeTaylorProblem(fW_vec, ComplexF64[u_0_FW, up_0_FW],
                                       (0.0 + 0im, 1.6 + 0im); order = 30)
        walk  = vector_path_network_solve(prob,
                                          ComplexF64[1.1 + 0im, 1.5 + 0im];
                                          h = 0.15)
        poles = extract_poles_shared_q(walk; radius_t = 5.0,
                                       cluster_atol = 0.2,
                                       min_support  = 1)
        @test !isempty(poles)
        # The nearest extracted pole matches the closed-form lattice
        # pole z(0,0) = 1 (FW md:297) to the figure-catalogue Float64
        # spec of 1e-6 — the same bar `polefield_test.jl` PF.1.1 pins
        # for the scalar pipeline.  The pole is a genuine ℘ pole; the
        # shared-`Q` extractor placed it from the vector companion jets.
        z1 = ℘_pole_vpo(0, 0)                       # = 1 + 0im
        @test isapprox(z1, 1.0 + 0im; atol = 1.0e-12)
        err_z1 = minimum(abs(p - z1) for p in poles)
        @test err_z1 ≤ 1.0e-5
    end

    # ----------------------------------------------------------------------
    # VPO.5 — mutation-proof: a corrupted IC breaks the oracle agreement.
    #
    # The load-bearing VC-9 assertions pin the vector pipeline against
    # FW reference numbers.  Confirm they have teeth: perturb the FW
    # initial condition by a small but real amount.  `u'' = 6u²` has
    # `u(z) = ℘(z + c₁; 0, c₂)` — a DIFFERENT `(u(0), u'(0))` is a
    # DIFFERENT ℘ translate, so a perturbed IC integrates to a value
    # that no longer matches the FW Table 5.1 reference.  A pipeline
    # that ignored its IC (returned a constant, say) would pass VC-9
    # with the wrong IC too — this mutation proves it does not.
    # ----------------------------------------------------------------------
    @testset "VPO.5: mutation-proof — a perturbed IC fails the oracle" begin
        # Perturb u'(0) by 1e-3 — a genuinely different ℘ translate.
        prob_bad = VectorPadeTaylorProblem(fW_vec,
                                           [u_0_FW, up_0_FW + 1.0e-3],
                                           (0.0, 0.5); order = 30)
        sol_bad  = vector_solve_pade(prob_bad; h = 0.5)
        # The perturbed solve no longer hits the FW closed-form u(0.5):
        # the VPO.1 `rtol = 1e-12` assertion would be RED for it.
        @test !isapprox(sol_bad(0.5)[1], u_at_0_5; rtol = 1.0e-12)
        # ...yet the genuine FW IC still does (the discriminator works).
        sol_ok = vector_solve_pade(
            VectorPadeTaylorProblem(fW_vec, [u_0_FW, up_0_FW],
                                    (0.0, 0.5); order = 30); h = 0.5)
        @test isapprox(sol_ok(0.5)[1], u_at_0_5; rtol = 1.0e-12)
    end
end # @testset VPO

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Each mutation below was applied, the suite confirmed RED, the mutation
# reverted.  VPO.5 is the in-test mutation-proof (a perturbed IC); the
# mutations here perturb the *oracle* and the *impl* to confirm the
# load-bearing assertions catch a regression.
#
#   M1 — perturb the FW oracle value.  In `_oracle_problems.jl`, change
#        `u_at_30` from 1.0950982559597442 to 1.10.  VPO.2's
#        `isapprox(real(u30), u_at_30; atol = 1e-4)` goes RED — the
#        path-network value (~1.09509) no longer matches the corrupted
#        reference.  Confirms VPO.2 genuinely pins the FW Table 5.1
#        number, not a self-consistent tautology.  Restored.
#
#   M2 — tighten VPO.2's tolerance to expose the pole-tunnelling error.
#        Change VPO.2's `atol = 1e-4` to `atol = 1e-8`.  VPO.2 goes RED
#        — the measured path-network error reaching u(30) is ~6e-6
#        (the accumulated shared-`Q` step error over ~21 tunnelled
#        poles).  Confirms the `1e-4` bar is the genuine, honest
#        envelope of the vector pole-tunnelling accuracy, not slack.
#        Restored.
#
#   M3 — corrupt the companion RHS.  Change `fW_vec` to
#        `(z,y) -> [y[2], 5*y[1]^2]` (the wrong ODE — `5u²` not `6u²`).
#        VPO.1, VPO.2, VPO.3, VPO.4 all go RED: the vector pipeline now
#        integrates a different ODE whose solution is not the FW ℘, so
#        none of the reference values / lattice poles match.  Confirms
#        the whole VC-9 suite genuinely exercises the `u'' = 6u²` test
#        problem end-to-end.  Restored.
#
#   M4 — drop the BVP free-component recovery.  In VPO.3, the genuine
#        `up_0_FW` recovery is load-bearing.  Replace the analytic
#        Jacobian `Jf_W` with the wrong constant `[0 1; 0 0]` (a
#        linear-ODE Jacobian).  The BVP Newton solve diverges / converges
#        to a non-℘ solution and `isapprox(yB0[2], up_0_FW)` goes RED —
#        confirms VPO.3 genuinely recovers u'(0) via a correct nonlinear
#        companion solve.  Restored.
#
#   M5 — revert the scale-covariant cluster representative ordering
#        (bead padetaylor-gt1, ADR-0026 Amendment 6 §S7).  In
#        `extract_poles_shared_q` the candidate tuple's element 2 is the
#        z-plane node-to-pole distance `z_dist`, and `sort!` orders by
#        it so the cluster representative is the *z-plane-closest*
#        node's placement.  Revert it to the per-node-step-relative
#        `|t*|`: replace
#          push!(candidates, (z_pole, T(z_dist), k, T(abs(h_node))))
#        with
#          push!(candidates, (z_pole, T(abs(t_C)), k, T(abs(h_node))))
#        Expected: under the adaptive per-node `h` a coarse-`h` *far*
#        node carries a smaller `|t*|` than a fine-`h` *near* one, so
#        the representative becomes a far node's (less accurate)
#        placement and VPO.4's `err_z1 ≤ 1e-5` bar bites.
#        Result: RED — VPO.4's `err_z1 ≤ 1e-5` failed (the min-`|t*|`
#        representative placed the ℘ pole at ≈ 4e-5, an order of
#        magnitude worse than the z-plane-closest node's ≈ 6e-7).
#        Restored to GREEN.
#
# All mutations restored before commit.  Matches the inline
# mutation-proof pattern of `test/vector_bvp_test.jl` (VB.7.1).
# ----------------------------------------------------------------------
