# test/vector_bvp_test.jl -- bead `padetaylor-0ln.24` (VB1) tests.
#
# First-order vector-system Chebyshev-collocation Newton BVP solver.
# `VectorBVP.vector_bvp_solve` solves  y'(z) = f(z, y),  y ∈ ℂ^d, on a
# complex segment under a general linear two-point boundary condition
# B_a·y(z_a) + B_b·y(z_b) = g.
#
# Oracles (all closed-form — Rule 5, every assertion checks a known value):
#   VB.1.1  d=1 reduction          : scalar linear ODE y' = y, closed form e^z.
#   VB.1.x  linear vector oracle   : harmonic oscillator y'=(y₂,-y₁), the
#                                    closed form (cos z, -sin z), real AND
#                                    complex segments, Float64 AND BF-256.
#   VB.2.x  D1 matrix pins         : D1 differentiates a polynomial exactly.
#   VB.3.x  autodiff Jacobian      : Taylor1 ∂f/∂y vs hand-computed analytic.
#   VB.4.x  fail-fast guards       : each throws the documented exception.
#   VB.5.x  BigFloat-256 path.
#
# Reference: docs/adr/0023-vector-bvp-solver.md, references/bvp_recipe.md,
# src/BVP.jl (the scalar second-order solver this mirrors).
#
# Standalone-runnable: `julia --project=. test/vector_bvp_test.jl`.

using Test

# VB1 does not wire VectorBVP into the PadeTaylor umbrella (that is VB3's
# job); include the module file directly.
if !isdefined(Main, :VectorBVP)
    include(joinpath(@__DIR__, "..", "src", "VectorBVP.jl"))
end
using .VectorBVP

@testset "VectorBVP (VB1): first-order vector Chebyshev-collocation BVP" begin

    # -------------------------------------------------------------------------
    # VB.1.1 — d = 1 reduction: scalar linear ODE y' = y on [0, 1].
    # Closed form y(z) = y(0)·e^z.  BC pins y at z_a: 1·y(0) + 0·y(1) = 1,
    # so B_a = [1;;], B_b = [0;;], g = [1].  Solution y(z) = e^z.
    # -------------------------------------------------------------------------
    @testset "VB.1.1: d=1 reduction matches scalar linear ODE e^z" begin
        f1(z, y) = [y[1]]
        sol = vector_bvp_solve(f1, 0.0, 1.0,
                               reshape([1.0], 1, 1), reshape([0.0], 1, 1),
                               [1.0]; N = 24)
        @test sol.iterations ≤ 2                         # linear ⇒ 1 Newton step
        @test sol.residual_inf ≤ 1e-12
        # Node values vs closed form e^z.
        for j in 1:sol.N+1
            @test isapprox(sol.Y_nodes[j][1], exp(sol.nodes_z[j]); atol = 1e-12)
        end
        # Barycentric callable at off-node points.
        for z in (0.17, 0.5, 0.83)
            @test isapprox(sol(z)[1], exp(z); atol = 1e-12)
        end
    end

    # -------------------------------------------------------------------------
    # VB.1.2 — linear vector oracle: harmonic oscillator as a d=2 system
    # y' = (y₂, -y₁).  Closed form y₁(z) = cos z, y₂(z) = -sin z.
    # BC pins y₁ at both ends:  B_a = [1 0; 0 0], B_b = [0 0; 1 0],
    # g = [cos z_a; cos z_b].  Real segment [0, π/2].
    # -------------------------------------------------------------------------
    f_ho(z, y) = [y[2], -y[1]]
    @testset "VB.1.2: harmonic oscillator d=2 on real [0, π/2]" begin
        z_a, z_b = 0.0, π/2
        B_a = [1.0 0.0; 0.0 0.0]
        B_b = [0.0 0.0; 1.0 0.0]
        g   = [cos(z_a), cos(z_b)]
        sol = vector_bvp_solve(f_ho, z_a, z_b, B_a, B_b, g; N = 28)
        @test sol.iterations ≤ 2                          # linear
        @test sol.residual_inf ≤ 1e-12
        for j in 1:sol.N+1
            zj = sol.nodes_z[j]
            @test isapprox(sol.Y_nodes[j][1],  cos(zj); atol = 1e-12)
            @test isapprox(sol.Y_nodes[j][2], -sin(zj); atol = 1e-12)
        end
        # Barycentric callable at off-node points.
        for z in (0.31, 0.9, 1.4)
            yz = sol(z)
            @test isapprox(yz[1],  cos(z); atol = 1e-12)
            @test isapprox(yz[2], -sin(z); atol = 1e-12)
        end
    end

    @testset "VB.1.3: harmonic oscillator d=2 on a complex segment" begin
        # Same system, segment from 0 to 0.6 + 0.4im — exercises the
        # complex affine scale factor s = (z_b - z_a)/2.
        z_a = 0.0 + 0.0im
        z_b = 0.6 + 0.4im
        B_a = ComplexF64[1 0; 0 0]
        B_b = ComplexF64[0 0; 1 0]
        g   = ComplexF64[cos(z_a), cos(z_b)]
        sol = vector_bvp_solve(f_ho, z_a, z_b, B_a, B_b, g; N = 30)
        @test sol.iterations ≤ 2
        @test sol.residual_inf ≤ 1e-11
        for j in 1:sol.N+1
            zj = sol.nodes_z[j]
            @test isapprox(sol.Y_nodes[j][1],  cos(zj); atol = 1e-11)
            @test isapprox(sol.Y_nodes[j][2], -sin(zj); atol = 1e-11)
        end
        # Callable at an interior complex point.
        zq = 0.3 + 0.2im
        yq = sol(zq)
        @test isapprox(yq[1],  cos(zq); atol = 1e-11)
        @test isapprox(yq[2], -sin(zq); atol = 1e-11)
    end

    @testset "VB.1.4: caller-supplied jacobian matches autodiff path" begin
        # Same harmonic-oscillator oracle, but pass the constant analytic
        # Jacobian ∂f/∂y = [0 1; -1 0] explicitly.
        z_a, z_b = 0.0, π/2
        B_a = [1.0 0.0; 0.0 0.0]
        B_b = [0.0 0.0; 1.0 0.0]
        g   = [cos(z_a), cos(z_b)]
        Jf(z, y) = [0.0 1.0; -1.0 0.0]
        sol = vector_bvp_solve(f_ho, z_a, z_b, B_a, B_b, g; N = 28, jacobian = Jf)
        for j in 1:sol.N+1
            zj = sol.nodes_z[j]
            @test isapprox(sol.Y_nodes[j][1],  cos(zj); atol = 1e-12)
            @test isapprox(sol.Y_nodes[j][2], -sin(zj); atol = 1e-12)
        end
    end

    # -------------------------------------------------------------------------
    # VB.2.x — D1 matrix pins.  The Chebyshev D1 differentiates a polynomial
    # of degree ≤ N exactly: for p(t) = t³ - 2t + 1, D1·p_nodes = p'_nodes.
    # -------------------------------------------------------------------------
    @testset "VB.2.1: D1 differentiates a cubic exactly" begin
        N  = 8
        t  = [cos(j*π/N) for j in 0:N]
        D1 = VectorBVP._chebyshev_D1(t, Float64, N)
        p      = [tj^3 - 2tj + 1   for tj in t]
        dp     = [3tj^2 - 2        for tj in t]
        @test isapprox(D1 * p, dp; atol = 1e-11)
        # D1 annihilates a constant (negative-sum-of-row identity).
        @test maximum(abs, D1 * ones(N+1)) < 1e-12
    end

    @testset "VB.2.2: D1 pinned entries vs Trefethen cheb.m formula" begin
        # N=4 corner entries (Trefethen SMIM cheb.m).
        N  = 4
        t  = [cos(j*π/N) for j in 0:N]
        D1 = VectorBVP._chebyshev_D1(t, Float64, N)
        # Top-left corner D1[1,1] = (2N²+1)/6 for Chebyshev-2 grid.
        @test isapprox(D1[1, 1], (2*N^2 + 1)/6; atol = 1e-12)
        # Bottom-right is its negative by the grid symmetry.
        @test isapprox(D1[N+1, N+1], -(2*N^2 + 1)/6; atol = 1e-12)
        # Row sums vanish.
        for i in 1:N+1
            @test abs(sum(D1[i, :])) < 1e-12
        end
    end

    # -------------------------------------------------------------------------
    # VB.3.x — autodiff Jacobian.  For a known nonlinear f the Taylor1
    # Jacobian must equal the hand-computed analytic d×d matrix.
    # -------------------------------------------------------------------------
    @testset "VB.3.1: Taylor1 autodiff Jacobian vs analytic" begin
        # f(z, y) = [y₁²·y₂ + z·y₁,  sin(y₁) + y₂³]   (d = 2).
        # Analytic ∂f/∂y =
        #   [ 2y₁y₂ + z      y₁²    ]
        #   [ cos(y₁)        3y₂²   ]
        fnl(z, y) = [y[1]^2 * y[2] + z*y[1],  sin(y[1]) + y[2]^3]
        z0 = 0.7
        y0 = [1.3, -0.4]
        Jad = VectorBVP._autodiff_jacobian(fnl, z0, y0, 2, Float64)
        Jan = [2*y0[1]*y0[2] + z0   y0[1]^2;
               cos(y0[1])           3*y0[2]^2]
        @test isapprox(Jad, Jan; atol = 1e-13)

        # Complex argument — same f, complex z and y.
        zc = 0.2 + 0.5im
        yc = [0.6 - 0.1im, 0.3 + 0.4im]
        Jadc = VectorBVP._autodiff_jacobian(fnl, zc, yc, 2, ComplexF64)
        Janc = [2*yc[1]*yc[2] + zc   yc[1]^2;
                cos(yc[1])           3*yc[2]^2]
        @test isapprox(Jadc, Janc; atol = 1e-13)
    end

    # -------------------------------------------------------------------------
    # VB.4.x — fail-fast guards (CLAUDE.md Rule 1).
    # -------------------------------------------------------------------------
    @testset "VB.4.1: fail-fast guards throw the documented exceptions" begin
        f1(z, y) = [y[1]]
        Ba1 = reshape([1.0], 1, 1); Bb1 = reshape([0.0], 1, 1); g1 = [1.0]

        # N too small.
        @test_throws ArgumentError vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1; N = 3)
        # maxiter zero.
        @test_throws ArgumentError vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1;
                                                    N = 8, maxiter = 0)
        # Negative tol.
        @test_throws ArgumentError vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1;
                                                    N = 8, tol = -1.0)
        # B_a wrong shape (2×2 but d = length(g) = 1).
        @test_throws ArgumentError vector_bvp_solve(f1, 0.0, 1.0,
                                                    [1.0 0.0; 0.0 1.0], Bb1, g1; N = 8)
        # RHS returns a wrong-length vector.
        f_bad(z, y) = [y[1], y[1]]                # length 2, d = 1
        @test_throws ArgumentError vector_bvp_solve(f_bad, 0.0, 1.0, Ba1, Bb1, g1; N = 8)
        # Non-convergent Newton: nonlinear problem, maxiter = 1, wild guess.
        fnl(z, y) = [y[1]^2 + 1.0]
        bad_init(z) = [1.0e6]
        @test_throws ErrorException vector_bvp_solve(fnl, 0.0, 1.0, Ba1, Bb1, g1;
                                                     N = 8, maxiter = 1,
                                                     initial_guess = bad_init)
        # Out-of-segment evaluation.
        sol = vector_bvp_solve(f1, 0.0, 1.0, Ba1, Bb1, g1; N = 12)
        @test_throws DomainError sol(5.0)
    end

    # -------------------------------------------------------------------------
    # VB.5.x — BigFloat-256 precision path.
    # -------------------------------------------------------------------------
    @testset "VB.5.1: BigFloat-256 harmonic oscillator" begin
        setprecision(BigFloat, 256) do
            z_a, z_b = big(0.0), big(π)/2
            B_a = BigFloat[1 0; 0 0]
            B_b = BigFloat[0 0; 1 0]
            g   = BigFloat[cos(z_a), cos(z_b)]
            sol = vector_bvp_solve(f_ho, z_a, z_b, B_a, B_b, g; N = 32)
            @test sol.residual_inf ≤ big(1e-12)
            for j in 1:sol.N+1
                zj = sol.nodes_z[j]
                @test isapprox(sol.Y_nodes[j][1],  cos(zj); atol = big(1e-25))
                @test isapprox(sol.Y_nodes[j][2], -sin(zj); atol = big(1e-25))
            end
            # Barycentric callable at BF-256.
            zq = big(1.0)
            yq = sol(zq)
            @test isapprox(yq[1],  cos(zq); atol = big(1e-25))
            @test isapprox(yq[2], -sin(zq); atol = big(1e-25))
        end
    end

end # @testset VectorBVP

# VB.6.1  Mutation-proof procedure (verified 2026-05-21 before commit;
# CLAUDE.md Rule 4 — perturb the impl, confirm specific VB.* go RED, restore):
#
#   Mutation A  --  in `vector_bvp_solve`, drop the affine scale factor `s`
#     on the residual: replace `R[...] .-= s .* CT.(fj)` with
#     `R[...] .-= CT.(fj)`.  The collocation operator dY/dt = s·f is then
#     mis-scaled by 1/s.
#     Verified bite (2026-05-21): VB.1.1, VB.1.2, VB.1.3, VB.1.4 RED and
#     VB.4.1 + VB.5.1 errored — 6 errored, 17 passed.  Newton fails to
#     converge (the linearised problem no longer has the right fixed point)
#     so the non-convergence ErrorException fires.  The `s` factor is the
#     single most load-bearing element, exactly as the analogous Mutation A
#     bit scalar BVP.jl.  Restored.
#
#   Mutation B  --  in `vector_bvp_solve`, swap the τ-method endpoint↔node
#     pairing: write `J_bc[:, 1:d] .= Ba` and `J_bc[:, nodeN+1:nodeN+d] .= Bb`
#     (and swap `Yza`/`Yzb` in the residual the same way).  The BC then
#     multiplies B_a against y(z_b) and B_b against y(z_a) — the descending
#     DMSUITE node order (t_0=+1↦z_b, t_N=-1↦z_a) is mis-read.
#     Verified bite (2026-05-21): VB.1.1, VB.1.2, VB.1.3, VB.1.4, VB.5.1
#     RED — 281 failed, 26 passed.  Newton still converges (residual
#     ~1e-15) but to the REFLECTED solution y₁(z)=cos(z_a+z_b-z), which
#     violates the user's BC.  This is the exact bug caught during VB1
#     development; the reflected-solution symptom is in ADR-0023.  Restored.
#
#   Mutation C  --  in `_autodiff_jacobian`, read coefficient `[0]` instead
#     of `[1]`: `Jf[r, i] = ft[r][0]`.  The Jacobian becomes the RHS value,
#     not its derivative.
#     Verified bite (2026-05-21): VB.3.1 RED (both the real and complex
#     asserts) — the autodiff Jacobian no longer matches the analytic
#     matrix.  VB.1.2/1.3 + VB.5.1 also bite (3 failed, 3 errored) because
#     Newton's Jacobian is then wrong on the autodiff path.  VB.1.4 stays
#     fully GREEN — it passes an explicit `jacobian` and so never touches
#     `_autodiff_jacobian`: the suite localises the fault.  Restored.
#
# All three mutations restored before commit.  Matches the inline
# mutation-proof pattern of test/bvp_test.jl (BV.6.1).

# Standalone entry point: `julia --project=. test/vector_bvp_test.jl`.
