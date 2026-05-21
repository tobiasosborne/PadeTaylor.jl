# test/ext_gridap_test.jl -- bead `padetaylor-0ln.36` tests.
#
# Validates the `PadeTaylorGridapExt` package extension per ADR-0003 +
# ADR-0024: `laplace2d_solve_gridap` solves the Dirichlet Laplace
# problem `∇²φ = 0` on a rectangle by Gridap.jl finite elements — voter
# (3) of the triple-method majority-vote harmonic-extension fill, an
# algorithmically disjoint cross-check of the in-house 2D-Chebyshev
# spectral solver `laplace2d_solve` (voter (2), src/Laplace2D.jl).
#
# The assertions are against KNOWN-CORRECT values (CLAUDE.md Rule 5):
# closed-form harmonic functions and their exact boundary traces.
#
# Test groups:
#   GRD.1.1  FEM solve reproduces φ = x²−y²       to FEM accuracy.
#   GRD.1.2  FEM solve reproduces φ = exp(x)cos(y) to FEM accuracy.
#   GRD.1.3  FEM solve reproduces φ = x³−3xy²      to FEM accuracy.
#   GRD.1.4  Algebraic O(h^p) convergence: refine the mesh, the error
#            drops at the expected polynomial rate (NOT spectral).
#   GRD.2.1  Cross-check vs `laplace2d_solve`: FEM and spectral agree
#            to within FEM tolerance on the same problem (the core
#            verification value of the triple-method discipline).
#   GRD.3.1  Fail-fast guards (degenerate range, n too small).
#
# FEM accuracy note: order-1 Lagrangian elements converge as O(h²) in
# the field error (h ~ 1/n).  At n = 64 cells/axis the interior error
# on a smooth harmonic field is O(1e-4) — algebraic, NOT the geometric
# convergence of the spectral `laplace2d_solve`.  The tolerances below
# are sized to that O(h²) law and pinned by the GRD.1.4 refinement test.
#
# Reference: ext/PadeTaylorGridapExt.jl; docs/adr/0024-laplace-harmonic-
# extension.md; docs/adr/0003-extensions-pattern.md; src/Laplace2D.jl.

using Test
using PadeTaylor
using Gridap

# Closed-form harmonic test functions (∇²φ = 0 exactly).
const _G_QUAD = (x, y) -> x^2 - y^2
const _G_EXP  = (x, y) -> exp(x) * cos(y)
const _G_CUB  = (x, y) -> x^3 - 3x * y^2

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

@testset "GridapExt (bead padetaylor-0ln.36): FEM Laplace voter (3)" begin

    @testset "GRD.1.1: FEM reproduces φ = x² − y²" begin
        sol = laplace2d_solve_gridap(_bc(_G_QUAD, 0.0, 1.0, 0.0, 1.0)...,
                                     (0.0, 1.0), (0.0, 1.0); n_x = 64, n_y = 64)
        @test sol isa Laplace2DSolution
        # Order-1 FEM, n = 64: O(h²) ≈ O(2.4e-4); pin a generous 5e-3.
        @test _interior_err(sol, _G_QUAD) < 5e-3
        # The solve is genuinely non-trivial: the field varies O(1).
        @test abs(sol(0.5, 0.5) - _G_QUAD(0.5, 0.5)) < 5e-3
        @test sol(0.8, 0.2) > sol(0.2, 0.8)        # x²−y² monotone structure
    end

    @testset "GRD.1.2: FEM reproduces φ = exp(x) cos(y)" begin
        sol = laplace2d_solve_gridap(_bc(_G_EXP, 0.0, 1.0, 0.0, 1.0)...,
                                     (0.0, 1.0), (0.0, 1.0); n_x = 64, n_y = 64)
        @test _interior_err(sol, _G_EXP) < 5e-3
        @test abs(sol(0.37, 0.61) - _G_EXP(0.37, 0.61)) < 5e-3
    end

    @testset "GRD.1.3: FEM reproduces φ = x³ − 3xy²" begin
        sol = laplace2d_solve_gridap(_bc(_G_CUB, 0.0, 1.0, 0.0, 1.0)...,
                                     (0.0, 1.0), (0.0, 1.0); n_x = 64, n_y = 64)
        @test _interior_err(sol, _G_CUB) < 5e-3
        @test abs(sol(0.43, 0.29) - _G_CUB(0.43, 0.29)) < 5e-3
    end

    @testset "GRD.1.4: algebraic O(h²) convergence on mesh refinement" begin
        # Order-1 Lagrangian FEM converges as O(h²) — NOT spectral.
        # Halving h (doubling n) must shrink the error by ≈ 4×.
        errs = Float64[]
        for n in (24, 48, 96)
            sol = laplace2d_solve_gridap(_bc(_G_EXP, 0.0, 1.0, 0.0, 1.0)...,
                                         (0.0, 1.0), (0.0, 1.0);
                                         n_x = n, n_y = n)
            push!(errs, _interior_err(sol, _G_EXP))
        end
        # Monotone decrease.
        @test errs[2] < errs[1]
        @test errs[3] < errs[2]
        # O(h²): each refinement gives a rate ≥ 1.6 (allow slack below
        # the ideal 2.0 for the discrete sampling of the worst point).
        @test log2(errs[1] / errs[2]) > 1.6
        @test log2(errs[2] / errs[3]) > 1.6
        # And it is NOT machine-precision spectral — a coarse FEM mesh
        # must carry a visible algebraic error.
        @test errs[1] > 1e-5
    end

    @testset "GRD.2.1: cross-check — FEM agrees with 2D-Chebyshev spectral" begin
        # The core triple-method verification value (ADR-0024): the
        # Gridap FEM voter (3) and the in-house spectral voter (2) solve
        # the SAME Dirichlet Laplace problem and must agree to within
        # the FEM tolerance.  The spectral solver is at machine
        # precision on these harmonic fields, so FEM-vs-spectral
        # disagreement ≈ the FEM discretisation error.
        rect = (0.0, 1.0, 0.0, 1.0)
        for g in (_G_QUAD, _G_EXP, _G_CUB)
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
        end
    end

    @testset "GRD.3.1: fail-fast guards (CLAUDE.md Rule 1)" begin
        bc = _bc(_G_QUAD, 0.0, 1.0, 0.0, 1.0)
        # Degenerate ranges.
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (1.0, 1.0), (0.0, 1.0))
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (0.0, 1.0), (1.0, 0.0))
        # Mesh too coarse.
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (0.0, 1.0), (0.0, 1.0); n_x = 1)
        @test_throws ArgumentError laplace2d_solve_gridap(
            bc..., (0.0, 1.0), (0.0, 1.0); n_y = 0)
    end

end # @testset GridapExt

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Applied to `ext/PadeTaylorGridapExt.jl`, suite confirmed RED, reverted.
#
#   M1 -- the weak form.  In `laplace2d_solve_gridap`, change the bilinear
#         form `aform(u,v) = ∫( ∇(u)⋅∇(v) )*dΩ` to a mass form
#         `∫( u*v )*dΩ`.  Bite: the operator no longer discretises the
#         Laplacian, so GRD.1.1/1.2/1.3's `_interior_err < 5e-3` go RED
#         (the FEM field is wrong by O(1)) and GRD.2.1's FEM-vs-spectral
#         `diff < 5e-3` goes RED.  Restored.
#
#   M2 -- boundary-edge dispatch.  In the Dirichlet datum `g`, swap the
#         `bc_x_lo` / `bc_x_hi` branches (left edge gets the right-edge
#         data).  Bite: the boundary trace is wrong on two edges, so
#         GRD.1.1's `_interior_err < 5e-3` goes RED and GRD.2.1's
#         FEM-vs-spectral agreement goes RED -- the FEM solve carries
#         the wrong Dirichlet data.  Restored.
#
#   M3 -- the affine-scale on output.  In the dense-output Chebyshev
#         grid, drop the `(b-a)/2` factor in `xs` (collapse the grid to
#         [-1,1] regardless of the rectangle).  Bite: GRD.2.1 goes RED --
#         the FEM solution is sampled off the rectangle and disagrees
#         with the spectral solver.  Restored.
# ----------------------------------------------------------------------
