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

# VB.6.x cross-validates the vector solver against the scalar second-order
# `BVP.bvp_solve`, which is itself oracle-pinned (DMSUITE/mpmath) in
# test/bvp_test.jl + test/_oracle_bvp.jl.  `BVP.jl` is self-contained
# (`using LinearAlgebra` only), so it is `include`d directly here — the
# same standalone-runnable idiom VB1 uses for `VectorBVP.jl`.
if !isdefined(Main, :BVP)
    include(joinpath(@__DIR__, "..", "src", "BVP.jl"))
end
using .BVP

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
                                                    N = 8, max_iter = 0)
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
                                                     N = 8, max_iter = 1,
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

    # -------------------------------------------------------------------------
    # VB.6.x — NONLINEAR cross-validation against the scalar BVP oracle.
    #
    # The scalar Painlevé-I BVP  u'' = 6u² + z  is recast as a d=2 first-order
    # companion system:  y = (u, u'),  y' = f(z, y) = (y₂, 6·y₁² + z).
    # `vector_bvp_solve` on that system, with B_a = [1 0; 0 0], B_b = [0 0; 1 0],
    # g = [u(z_a); u(z_b)] pinning y₁ at both ends, must reproduce the scalar
    # `BVP.bvp_solve` solution: vector y₁ ≡ scalar u, vector y₂ ≡ scalar u'.
    #
    # The scalar config is the PI BVP on [-18,-14] N=20 — the FW 2011 Fig 4.6
    # case validated in test/bvp_test.jl BV.1.3 against the DMSUITE/Octave
    # oracle constants `test_bvp_pi_real_N20_u` / `..._nodes_z` in
    # test/_oracle_bvp.jl (Group 3).  Cross-checking the vector solver against
    # `bvp_solve` here therefore transitively validates it against that
    # independent oracle (CLAUDE.md Rule 4 port-and-verify).
    #
    # This run uses the AUTODIFF Jacobian (jacobian = nothing), so it also
    # exercises the Taylor1 autodiff path on a genuinely nonlinear RHS — the
    # d=1 linear reduction was already covered by VB.1.1, so it is not re-added.
    # -------------------------------------------------------------------------
    @testset "VB.6.x: nonlinear PI companion system vs scalar BVP oracle" begin
        # Scalar PI ODE  u'' = 6u² + z  and its companion 4-vector form... d=2.
        f_PI_scalar(z, u)  = 6 * u^2 + z
        ∂f_PI_scalar(z, u) = 12 * u
        f_PI_vec(z, y)     = [y[2], 6 * y[1]^2 + z]   # companion: y'=(y₂,6y₁²+z)

        z_a, z_b = -18.0, -14.0
        sqrt_guess(z)    = sqrt(-z / 6)               # u ≈ √(-z/6) asymptotic
        # u' of the asymptotic guess: d/dz √(-z/6) = -1 / (2·√6·√(-z)).
        sqrt_guess_up(z) = -1 / (2 * sqrt(6) * sqrt(-z))
        u_a, u_b = sqrt_guess(z_a), sqrt_guess(z_b)

        # Scalar oracle-validated solution (BV.1.3 config, N = 20).
        scalar_sol = bvp_solve(f_PI_scalar, ∂f_PI_scalar, z_a, z_b, u_a, u_b;
                               N = 20, initial_guess = sqrt_guess)

        # Vector companion-system solve.  B_a pins y₁ at z_a, B_b pins y₁ at
        # z_b; g = [u(z_a); u(z_b)].  Autodiff Jacobian (jacobian = nothing).
        B_a = [1.0 0.0; 0.0 0.0]
        B_b = [0.0 0.0; 1.0 0.0]
        g   = [u_a, u_b]
        vec_init(z) = [sqrt_guess(z), sqrt_guess_up(z)]
        vec_sol = vector_bvp_solve(f_PI_vec, z_a, z_b, B_a, B_b, g;
                                   N = 20, initial_guess = vec_init)

        @testset "VB.6.1: Newton convergence sanity (autodiff path)" begin
            # Nonlinear ⇒ a few Newton steps; FW 2011 line 190 reports ≤ 6.
            @test vec_sol.iterations ≤ 6
            @test vec_sol.iterations ≥ 2          # genuinely nonlinear, not 1-step
            # Discrete residual floors at cond(D1⊗I)·eps — same order as the
            # scalar BV.1.3 floor (~1e-12); assert it is tiny.
            @test vec_sol.residual_inf ≤ 1e-10
            # Same Chebyshev grid + affine map ⇒ identical mapped nodes.
            @test isapprox(vec_sol.nodes_z, scalar_sol.nodes_z; atol = 1e-12)
        end

        @testset "VB.6.2: companion y₁ node values ≡ scalar u (oracle-pinned)" begin
            # vector y₁ at every node must equal the scalar BVP's u — which
            # BV.1.3 pins to test_bvp_pi_real_N20_u to atol 1e-10.
            #
            # The scalar D₂-collocation method and the vector D₁⊗I τ-method
            # are TWO DIFFERENT discretizations of the same continuous PI BVP;
            # each converges to its own discrete fixed point, and those node
            # values differ from each other by the DIFFERENCE of their
            # Chebyshev truncation errors — empirically ~7.5e-13 at N=20 on
            # [-18,-14] (the ~1e-12 floor that _oracle_bvp.jl Group 3 records).
            # That gap is set by N, NOT by the float precision: VB.6.6 below
            # confirms the BF-256 gap is the SAME ~7.5e-13.  atol = 1e-11
            # therefore certifies the cross-validation with honest headroom
            # above the (precision-independent) truncation floor.
            for j in 1:vec_sol.N+1
                @test isapprox(vec_sol.Y_nodes[j][1], scalar_sol.u_nodes[j];
                               atol = 1e-11)
            end
        end

        @testset "VB.6.3: companion y₂ node values ≡ scalar u'" begin
            # vector y₂ ≡ u'.  The scalar callable gives (u, u') at each node;
            # compare to the vector solver's second component.
            for j in 1:vec_sol.N+1
                zj          = scalar_sol.nodes_z[j]
                (_, up_sca) = scalar_sol(zj)
                @test isapprox(vec_sol.Y_nodes[j][2], up_sca; atol = 1e-9)
            end
        end

        @testset "VB.6.4: barycentric callable ≡ scalar (u, u') off-node" begin
            # The vector barycentric callable's components must track the
            # scalar callable at interior off-node points.
            for z in (-17.3, -16.0, -14.7)
                yz          = vec_sol(z)
                (u_s, up_s) = scalar_sol(z)
                @test isapprox(yz[1], u_s;  atol = 1e-10)
                @test isapprox(yz[2], up_s; atol = 1e-9)
            end
        end

        @testset "VB.6.5: explicit analytic Jacobian agrees with autodiff" begin
            # Analytic 2×2 companion Jacobian ∂f/∂y = [0 1; 12·y₁ 0].
            Jf(z, y) = [0.0 1.0; 12 * y[1] 0.0]
            vec_sol_J = vector_bvp_solve(f_PI_vec, z_a, z_b, B_a, B_b, g;
                                         N = 20, initial_guess = vec_init,
                                         jacobian = Jf)
            for j in 1:vec_sol_J.N+1
                @test isapprox(vec_sol_J.Y_nodes[j][1], scalar_sol.u_nodes[j];
                               atol = 1e-12)
                @test isapprox(vec_sol_J.Y_nodes[j][2],
                               vec_sol.Y_nodes[j][2]; atol = 1e-11)
            end
        end

        @testset "VB.6.6: BigFloat-256 companion system vs scalar BVP" begin
            # Same cross-validation at BF-256.  CRITICAL (CLAUDE.md Rule 2):
            # the scalar↔vector node-value gap is ~7.5e-13 at N=20 — and it
            # is the SAME ~7.5e-13 at Float64 and at BigFloat-256.  That is
            # the root-cause signature of a *discretization* difference (two
            # distinct Chebyshev schemes for one continuous PI BVP), NOT a
            # rounding difference: raising the precision F64 → BF-256 drives
            # each solver's *residual* from ~1e-12 to ~1e-74 yet leaves the
            # mutual gap untouched.  So the cross-validation atol is the
            # N=20 truncation floor (~1e-12), precision-independent — a
            # tighter atol would be wrong (and was caught, not patched away).
            setprecision(BigFloat, 256) do
                za, zb = big(-18.0), big(-14.0)
                sg(z)   = sqrt(-z / 6)
                sg_up(z) = -1 / (2 * sqrt(big(6)) * sqrt(-z))
                ua, ub  = sg(za), sg(zb)
                sca = bvp_solve(f_PI_scalar, ∂f_PI_scalar, za, zb, ua, ub;
                                N = 20, initial_guess = sg)
                Ba  = BigFloat[1 0; 0 0]
                Bb  = BigFloat[0 0; 1 0]
                gv  = BigFloat[ua, ub]
                vinit(z) = [sg(z), sg_up(z)]
                vec = vector_bvp_solve(f_PI_vec, za, zb, Ba, Bb, gv;
                                       N = 20, initial_guess = vinit)
                @test vec.iterations ≤ 6
                # Newton drives the BF-256 residual far below the F64 floor —
                # proof the solver itself reaches machine-precision accuracy.
                @test vec.residual_inf ≤ big(1e-40)
                # The truncation-limited cross-validation gap.
                gapB = maximum(abs(vec.Y_nodes[j][1] - sca.u_nodes[j])
                               for j in 1:vec.N+1)
                @test gapB ≤ big(1e-11)
                for j in 1:vec.N+1
                    @test isapprox(vec.Y_nodes[j][1], sca.u_nodes[j];
                                   atol = big(1e-11))
                end
            end
        end
    end

end # @testset VectorBVP

# VB.7.1  Mutation-proof procedure (Mutations A–C verified 2026-05-21 for
# VB1; Mutation D verified 2026-05-21 for VB2; CLAUDE.md Rule 4 — perturb
# the impl, confirm specific VB.* go RED, restore).  Renumbered from the
# VB1 label `VB.6.1` so it no longer collides with the new VB.6.x
# nonlinear-cross-validation testsets.
#
#   Mutation A  --  in the `_residual` helper, drop the affine scale factor
#     `s` on the residual: replace `R[rng] .-= ctx.s .* ctx.CT.(fj)` with
#     `R[rng] .-= ctx.CT.(fj)`.  The collocation operator dY/dt = s·f is then
#     mis-scaled by 1/s.  (Post-VB2 the residual lives in `_residual`; the
#     `s` factor is the same load-bearing element.)
#     Verified bite (2026-05-21, VB1): VB.1.1, VB.1.2, VB.1.3, VB.1.4 RED and
#     VB.4.1 + VB.5.1 errored.  Newton fails to converge (the linearised
#     problem no longer has the right fixed point) so the non-convergence
#     ErrorException fires.  The `s` factor is the single most load-bearing
#     element, exactly as the analogous Mutation A bit scalar BVP.jl.
#     Restored.
#
#   Mutation B  --  swap the τ-method endpoint↔node pairing: write
#     `J_bc[:, 1:d] .= Ba` and `J_bc[:, nodeN+1:nodeN+d] .= Bb` in
#     `vector_bvp_solve`, and swap `ctx.Ba`/`ctx.Bb` (or the node-N/node-0
#     slices) in `_residual`'s BC-row line.  The BC then multiplies B_a
#     against y(z_b) and B_b against y(z_a) — the descending DMSUITE node
#     order (t_0=+1↦z_b, t_N=-1↦z_a) is mis-read.
#     Verified bite (2026-05-21, VB1): VB.1.1, VB.1.2, VB.1.3, VB.1.4,
#     VB.5.1 RED.  Newton still converges (residual ~1e-15) but to the
#     REFLECTED solution y₁(z)=cos(z_a+z_b-z), which violates the user's BC.
#     This is the exact bug caught during VB1 development; the
#     reflected-solution symptom is in ADR-0023.  Restored.
#
#   Mutation C  --  in `_autodiff_jacobian`, read coefficient `[0]` instead
#     of `[1]`: `Jf[r, i] = ft[r][0]`.  The Jacobian becomes the RHS value,
#     not its derivative.
#     Verified bite (2026-05-21, VB1): VB.3.1 RED (both the real and complex
#     asserts) — the autodiff Jacobian no longer matches the analytic
#     matrix.  VB.1.2/1.3 + VB.5.1 also bite because Newton's Jacobian is
#     then wrong on the autodiff path.  VB.1.4 stays fully GREEN — it passes
#     an explicit `jacobian` and so never touches `_autodiff_jacobian`: the
#     suite localises the fault.  Restored.
#
#   Mutation D  (VB2 — the nonlinear cross-validation)  --  the same
#     `_autodiff_jacobian` perturbation as Mutation C (`ft[r][1]` → `ft[r][0]`),
#     re-verified against the VB.6.x companion-system tests.
#     Verified bite (2026-05-21, VB2): the VB.6.x testset ERRORS — building
#     the autodiff-path `vec_sol` for the nonlinear PI companion system
#     `y' = (y₂, 6y₁²+z)` throws the Newton non-convergence ErrorException,
#     because the wrong node Jacobian `[0 0; 0 0]` (the RHS value's `[0]`
#     coefficient at the linear/quadratic terms) destroys the Newton basin.
#     VB.1.4 (explicit-jacobian harmonic oscillator) stays fully GREEN — it
#     never touches `_autodiff_jacobian`, so the suite localises the fault
#     to the autodiff path.  This confirms the VB.6.x nonlinear
#     cross-validation genuinely exercises (and depends on) the Taylor1
#     autodiff Jacobian on a nonlinear RHS.  Restored.
#
# All mutations restored before commit.  Matches the inline mutation-proof
# pattern of test/bvp_test.jl (BV.6.1).

# Standalone entry point: `julia --project=. test/vector_bvp_test.jl`.
