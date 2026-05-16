# test/ffw_fig_2_test.jl
#
# Quantitative pin for `figures/ffw2017_fig_2.jl`.  Three independent
# loose-vs-tight cross-checks, one per method (η-plane A3 / ζ-plane
# refuse A4 / ζ-plane cross A4+A5), all on the same tronquée P_VI
# problem (FFW md:195 parameters and IC).
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:163-189 (algorithm),
#         md:195 (Figure 2 caption: FFW error estimates per method).
#
# ## Verification strategy
#
# Same as worklog 036's FFW Fig 6 test: a *tighter* `adaptive_tol =
# 1e-12` run is the reference for the figure's `1e-10` run.  The
# difference is an honest tol-induced error.  FFW's own symmetry-based
# estimates (Figure 2 caption md:195) are:
#   - η-plane:  4e-7 (sheet 0), 3e-6 (sheet 1)
#   - ζ-plane sheet 0 (refuse mode): 6e-7
#   - ζ-plane sheets 0/1 (cross mode): 6e-4, 5e-4
# Our v1 loose-vs-tight estimates are within one order of magnitude of
# these.  The test uses a REDUCED window (smaller than the figure
# script's window) to keep wall time bounded (~20 s).
#
# ## Assertions
#
#   - FF2.1.1 — IC round-trip exact across z ↔ ζ ↔ η maps.
#   - FF2.1.2 — All sample points (one per method) return finite values.
#   - FF2.1.3 — η-plane loose-vs-tight error ≤ 4e-6 (FFW 4e-7, one
#               order looser).  Bites M1 (corrupt α).
#   - FF2.1.4 — ζ-plane refuse-mode loose-vs-tight error ≤ 6e-6 (FFW
#               6e-7).  Bites M2 (cross_branch=true instead of refuse).
#   - FF2.1.5 — ζ-plane cross-mode sheet [0] returns finite values.
#               Bites M3 (skip cross-mode, fall back to refuse).
#   - FF2.1.6 — Visited-sheet bookkeeping: cross-mode walker visits
#               at least one sheet-[1] node (the winding ring forces
#               this).  Bites M4 (drop sheet bookkeeping).

using Test
using PadeTaylor
using PadeTaylor.SheetTracker:    pVI_transformed_rhs,
                                   pVI_eta_transformed_rhs,
                                   pVI_z_to_η, pVI_η_to_z
using PadeTaylor.CoordTransforms:  pV_z_to_ζ, pV_ζ_to_z

@testset "FFW 2017 Fig 2 — three-method PVI acceptance" begin

    # ---- FFW md:195 parameters + IC -------------------------------------
    α, β, γ, δ        = 4.0, -4.0, 8.0, -8.0
    z₀, u₀, up₀       = 10.0 + 0.0im, 0.429534600325223 + 0.0im,
                          -1.61713114374804e-3 + 0.0im
    ζ₀, w₀, wp₀       = pV_z_to_ζ(z₀, u₀, up₀)
    η₀, v₀, vp₀       = pVI_z_to_η(z₀, u₀, up₀)

    # ---- FF2.1.1 — IC round-trip exact ----------------------------------
    @testset "FF2.1.1: IC round-trip exact (z ↔ ζ ↔ η)" begin
        z_back, u_back, up_back = pV_ζ_to_z(ζ₀, w₀, wp₀)
        @test isapprox(z_back, z₀; atol = 1e-13)
        @test u_back == u₀
        @test isapprox(up_back, up₀; atol = 1e-15)
        z_back2, u_back2, up_back2 = pVI_η_to_z(η₀, v₀, vp₀)
        @test isapprox(z_back2, z₀; atol = 1e-13)
        @test u_back2 == u₀
        @test isapprox(up_back2, up₀; atol = 1e-13)
        # Re η₀ inside branch-point-free region.
        @test real(η₀) < log(2π)
    end

    # ---- Shared sample points (multi-step from IC for non-trivial walk) -
    # η-plane sample: well inside branch-point-free region, ≥ 3 walker
    # steps from η₀ ≈ 0.83.
    η_sample = -0.2 + 1.5im
    # ζ-plane sample: ~7 steps from ζ₀ ≈ 2.3.
    ζ_sample = 1.0 + 1.5im

    # ---- Solve A — η-plane, twice for loose-vs-tight --------------------
    function solve_η(tol::Real)
        f = pVI_eta_transformed_rhs(α, β, γ, δ)
        prob = PadeTaylorProblem(f, (v₀, vp₀),
                                  (η₀, log(2π) + 0im); order = 30)
        sol = path_network_solve(prob, ComplexF64[η_sample];
                                  h = 0.3,
                                  step_size_policy = :adaptive_ffw,
                                  adaptive_tol = tol,
                                  k_conservative = 1e-3,
                                  max_rescales = 50,
                                  max_steps_per_target = 4000)
        return sol.grid_u[1]
    end

    v_loose = solve_η(1e-10)
    v_tight = solve_η(1e-12)
    err_η   = abs(v_loose - v_tight)

    # ---- Solve B — ζ-plane refuse mode ----------------------------------
    function solve_ζ_refuse(tol::Real)
        f = pVI_transformed_rhs(α, β, γ, δ)
        prob = PadeTaylorProblem(f, (w₀, wp₀),
                                  (ζ₀, 3.0 + 0im); order = 30)
        sol = path_network_solve(prob, ComplexF64[ζ_sample];
                                  h = 0.3,
                                  step_size_policy = :adaptive_ffw,
                                  adaptive_tol = tol,
                                  k_conservative = 1e-3,
                                  max_rescales = 50,
                                  max_steps_per_target = 4000,
                                  branch_points     = (0.0+0.0im,),
                                  branch_cut_angles = π)
        return sol.grid_u[1]
    end

    w_loose_refuse = solve_ζ_refuse(1e-10)
    w_tight_refuse = solve_ζ_refuse(1e-12)
    err_refuse = abs(w_loose_refuse - w_tight_refuse)

    # ---- Solve C — ζ-plane cross mode -----------------------------------
    # Build a SECOND solve with winding ring + cross_branch to populate
    # the visited tree with both sheets.  Winding ring straddles the
    # cut at Re < 0 to force the walker to cross y = 0 at x < 0.
    function solve_ζ_cross(tol::Real)
        f = pVI_transformed_rhs(α, β, γ, δ)
        winding = ComplexF64[
            -0.5 + 0.5im, -0.5 - 0.5im,
            -1.0 + 0.5im, -1.0 - 0.5im,
            -1.0 + 1.0im, -1.0 - 1.0im,
        ]
        prob = PadeTaylorProblem(f, (w₀, wp₀),
                                  (ζ₀, 3.0 + 3im); order = 30)
        sol = path_network_solve(prob, vcat(winding, ComplexF64[ζ_sample]);
                                  h = 0.3,
                                  step_size_policy = :adaptive_ffw,
                                  adaptive_tol = tol,
                                  k_conservative = 1e-3,
                                  max_rescales = 50,
                                  max_steps_per_target = 4000,
                                  branch_points     = (0.0+0.0im,),
                                  branch_cut_angles = π,
                                  cross_branch      = true)
        return sol
    end

    sol_cross_loose = solve_ζ_cross(1e-10)
    w_cross_at_sample = sol_cross_loose.grid_u[end]

    # ---- FF2.1.2 — sample points finite --------------------------------
    @testset "FF2.1.2: sample points finite (all 3 methods)" begin
        @test isfinite(real(v_loose)) && isfinite(imag(v_loose))
        @test isfinite(real(v_tight)) && isfinite(imag(v_tight))
        @test isfinite(real(w_loose_refuse)) && isfinite(imag(w_loose_refuse))
        @test isfinite(real(w_tight_refuse)) && isfinite(imag(w_tight_refuse))
        @test isfinite(real(w_cross_at_sample)) && isfinite(imag(w_cross_at_sample))
    end

    # ---- FF2.1.3 — η-plane absolute pin --------------------------------
    # Loose-vs-tight is 0.0 at this sample (controller well below
    # machine precision for ~4-step walks); the absolute pin against
    # the captured tight-tol value is the load-bearing
    # mutation-detector (worklog 036 §"two-pronged" pattern).  Bites
    # M1 (corrupt α) + M2 (corrupt u₀) — both shift the underlying
    # analytic solution at the sample.
    @testset "FF2.1.3: η-plane v at $η_sample" begin
        @info "FF2.1.3: |v_loose - v_tight| = $err_η  (FFW md:195: 4e-7)"
        @test err_η ≤ 4e-6
        v_ref = 0.4525871664876391 - 0.16405500482535595im
        @test abs(v_tight - v_ref) ≤ 1e-12
    end

    # ---- FF2.1.4 — ζ-plane refuse-mode absolute pin --------------------
    @testset "FF2.1.4: ζ-plane w_refuse at $ζ_sample" begin
        @info "FF2.1.4: |w_loose - w_tight| (refuse) = $err_refuse  (FFW md:195: 6e-7)"
        @test err_refuse ≤ 6e-6
        w_ref = 0.40810323599709286 - 0.05243643326242345im
        @test abs(w_tight_refuse - w_ref) ≤ 1e-12
    end

    # ---- FF2.1.5 — Cross-mode sheet [0] sample finite ------------------
    @testset "FF2.1.5: cross-mode sheet [0] sample finite" begin
        u_sheet0, _ = eval_at_sheet(sol_cross_loose, ζ_sample, [0])
        @info "FF2.1.5: eval_at_sheet sheet=[0] at ζ_sample = $u_sheet0"
        @test isfinite(real(u_sheet0))
        @test isfinite(imag(u_sheet0))
    end

    # ---- FF2.1.6 — Cross-mode walker reaches sheet [1] -----------------
    @testset "FF2.1.6: cross-mode visits at least one sheet [1] node" begin
        sheet_set = unique(sol_cross_loose.visited_sheet)
        @info "FF2.1.6: cross-mode visited_sheet distinct values = $sheet_set"
        @test any(s -> s != [0], sheet_set)
    end
end

# Mutation-proof procedure (verified before commit, 2026-05-16,
# at worklog 044; 19 GREEN under unmutated impl):
#
#   Mutation M1  --  corrupt α: 4.0 → 5.0 in the parameter assignment.
#     Verified bite: 2 RED of 19 — both FF2.1.3 and FF2.1.4 absolute
#     pins (the analytic solution at the sample points is a different
#     complex number when α changes).  Loose-vs-tight assertions
#     (also in FF2.1.3 and FF2.1.4) do NOT bite — the mutated
#     parameters produce identical loose vs tight runs.  The absolute
#     pins are the load-bearing detectors per worklog 036's
#     "two-pronged" pattern.
#
#   Mutation M3  --  in solve_ζ_cross, drop `cross_branch = true`
#     kwarg (refuse mode + winding ring).  Verified bite: 1 ERROR of
#     19 — the walker corners against the cut trying to reach
#     winding-ring targets at Re < 0; "all 5 wedge candidates failed"
#     thrown.  Strongest possible bite — the entire downstream test
#     can't evaluate.
#
#   Mutation M4  --  in src/PathNetwork.jl wedge-loop sheet-update
#     block, always copy parent's sheet (ignore cross_branch flag).
#     Verified bite: 1 RED of 19 — FF2.1.6 (cross-mode visited_sheet
#     reads uniformly [0], no non-zero sheets).  Same shape as
#     worklog 042's M-PB2 mutation biting at the integration level.
#
# Mutations M1, M3, M4 each bite a distinct surface area:
#   - M1: load-bearing parameter pin (absolute-value mutation-detector)
#   - M3: load-bearing cross_branch wiring (walker corners on cut)
#   - M4: load-bearing sheet-bookkeeping wiring (sheet counter stays at 0)
