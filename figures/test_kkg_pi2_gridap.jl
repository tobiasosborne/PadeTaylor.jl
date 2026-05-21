# figures/test_kkg_pi2_gridap.jl
#
# Figures-env verification of the Gridap FEM Laplace solver — voter (3)
# of the ADR-0024 triple-method majority-vote harmonic-extension fill.
#
# Run from the repo root with:
#
#     julia --project=figures figures/test_kkg_pi2_gridap.jl
#
# This is the figures-project home for the Gridap verification.  The
# main `test/` suite stays Gridap-free (Gridap is a `figures/`-project
# dependency, NOT a package weak-dep — see ADR-0024 and the helper
# file's header for the ForwardDiff-conflict root cause).
#
# What is asserted (CLAUDE.md Rule 5 — every assertion against a
# KNOWN-CORRECT value):
#
#   G.1  The Gridap FEM solve reproduces closed-form harmonic functions
#        (x²−y², exp(x)cos(y)) to FEM tolerance — order-1 Lagrangian
#        elements converge as O(h²), so at n = 64 the interior error is
#        O(1e-4); a generous 5e-3 is pinned.
#   G.2  Algebraic O(h²) convergence on mesh refinement (NOT spectral).
#   G.3  The Gridap FEM solve agrees with the in-house 2D-Chebyshev
#        spectral solver `PadeTaylor.laplace2d_solve` on the SAME
#        rectangle to within the FEM tolerance — the core triple-method
#        cross-check value (a FEM bug and a spectral bug cannot coincide).
#   G.4  Fail-fast guards (degenerate range, mesh too coarse).
#
# Reference: figures/_kkg_pi2_gridap_helper.jl; src/Laplace2D.jl;
# docs/adr/0024-laplace-harmonic-extension.md.

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_kkg_pi2_gridap_helper.jl"))

# Closed-form harmonic test functions (∇²φ = 0 exactly).
const _G_QUAD = (x, y) -> x^2 - y^2
const _G_EXP  = (x, y) -> exp(x) * cos(y)

# Edge-trace boundary data from a closed-form g on [a,b]×[c,d], in the
# `laplace2d_solve` / `laplace2d_solve_gridap` callable convention.
_bc(g, a, b, c, d) = (y -> g(a, y), y -> g(b, y),
                      x -> g(x, c), x -> g(x, d))

# Worst interior-point error of `sol` vs the closed form `g`.
function _interior_err(sol, g; rect = (0.0, 1.0, 0.0, 1.0), m = 9)
    a, b, c, d = rect
    err = 0.0
    for i in 1:m, j in 1:m
        x = a + (b - a) * i / (m + 1)
        y = c + (d - c) * j / (m + 1)
        err = max(err, abs(sol(x, y) - g(x, y)))
    end
    return err
end

@testset "Gridap FEM Laplace voter (3) — figures-env verification" begin

    @testset "G.1: FEM reproduces closed-form harmonic functions" begin
        sol_q = laplace2d_solve_gridap(_bc(_G_QUAD, 0.0, 1.0, 0.0, 1.0)...,
                                       (0.0, 1.0), (0.0, 1.0);
                                       n_x = 64, n_y = 64)
        @test sol_q isa Laplace2DSolution
        # Order-1 FEM, n = 64: O(h²) ≈ O(2.4e-4); pin a generous 5e-3.
        eq = _interior_err(sol_q, _G_QUAD)
        @test eq < 5e-3
        @test abs(sol_q(0.5, 0.5) - _G_QUAD(0.5, 0.5)) < 5e-3
        @test sol_q(0.8, 0.2) > sol_q(0.2, 0.8)   # x²−y² monotone structure
        println("  G.1  x²−y²        interior err = ", eq)

        sol_e = laplace2d_solve_gridap(_bc(_G_EXP, 0.0, 1.0, 0.0, 1.0)...,
                                       (0.0, 1.0), (0.0, 1.0);
                                       n_x = 64, n_y = 64)
        ee = _interior_err(sol_e, _G_EXP)
        @test ee < 5e-3
        @test abs(sol_e(0.37, 0.61) - _G_EXP(0.37, 0.61)) < 5e-3
        println("  G.1  exp(x)cos(y) interior err = ", ee)
    end

    @testset "G.2: algebraic O(h²) convergence on mesh refinement" begin
        errs = Float64[]
        for n in (24, 48, 96)
            sol = laplace2d_solve_gridap(_bc(_G_EXP, 0.0, 1.0, 0.0, 1.0)...,
                                         (0.0, 1.0), (0.0, 1.0);
                                         n_x = n, n_y = n)
            push!(errs, _interior_err(sol, _G_EXP))
        end
        @test errs[2] < errs[1]
        @test errs[3] < errs[2]
        # O(h²): each refinement gives a rate ≥ 1.6 (slack below the
        # ideal 2.0 for the discrete sampling of the worst point).
        @test log2(errs[1] / errs[2]) > 1.6
        @test log2(errs[2] / errs[3]) > 1.6
        # And it is NOT machine-precision spectral.
        @test errs[1] > 1e-5
        println("  G.2  errs (n=24,48,96) = ", errs,
                "  rates = ", (log2(errs[1]/errs[2]), log2(errs[2]/errs[3])))
    end

    @testset "G.3: cross-check — FEM agrees with 2D-Chebyshev spectral" begin
        # The core triple-method verification value (ADR-0024): the
        # Gridap FEM voter (3) and the in-house spectral voter (2) solve
        # the SAME Dirichlet Laplace problem and must agree to within
        # the FEM tolerance.
        rect = (0.0, 1.0, 0.0, 1.0)
        for g in (_G_QUAD, _G_EXP)
            fem  = laplace2d_solve_gridap(_bc(g, rect...)...,
                                          (0.0, 1.0), (0.0, 1.0);
                                          n_x = 64, n_y = 64)
            spec = laplace2d_solve(_bc(g, rect...)...,
                                   (0.0, 1.0), (0.0, 1.0); Nx = 28, Ny = 28)
            # Spectral hits machine precision on a smooth harmonic field.
            @test _interior_err(spec, g) < 1e-9
            # FEM vs spectral: agreement within the O(h²) FEM tolerance.
            diff = 0.0
            for i in 1:9, j in 1:9
                x = i / 10; y = j / 10
                diff = max(diff, abs(fem(x, y) - spec(x, y)))
            end
            @test diff < 5e-3
            println("  G.3  FEM-vs-spectral max disagreement = ", diff)
        end
    end

    @testset "G.4: fail-fast guards (CLAUDE.md Rule 1)" begin
        bc = _bc(_G_QUAD, 0.0, 1.0, 0.0, 1.0)
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (1.0, 1.0), (0.0, 1.0))
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (0.0, 1.0), (1.0, 0.0))
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (0.0, 1.0), (0.0, 1.0); n_x = 1)
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (0.0, 1.0), (0.0, 1.0); n_y = 0)
    end

end # @testset
