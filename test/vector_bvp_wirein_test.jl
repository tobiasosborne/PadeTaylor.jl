"""
VB3 tests for the `VectorBVP` umbrella wire-in + the `PainleveHierarchy`
analytic-Jacobian and BVP-convenience layer (bead `padetaylor-0ln.26`;
v0.2 plan row VB3; ADR-0023).

VB1+VB2 shipped `src/VectorBVP.jl` (`vector_bvp_solve`, `VectorBVPSolution`)
but did not wire it into the `PadeTaylor` umbrella — the VB1 test
`include`s the module file directly.  VB3 wires it in and adds:

  - `painleve_hierarchy_jacobian(:I, m; t)` — the EXACT analytic
    companion-form Jacobian `∂f/∂y` of the `P_I^(m)` RHS;
  - a `PainleveHierarchyProblem` method of `vector_bvp_solve` that
    boundary-value-solves a hierarchy member, defaulting the Jacobian to
    the analytic helper.

Tests (`VW.*`), each asserting an invariant against a known value
(CLAUDE.md Rule 5):

    VW.1   wire-in            — vector_bvp_solve / VectorBVPSolution are
                                reachable as PadeTaylor package symbols.
    VW.2   Jacobian correctness (the load-bearing test) — for m=1 and
                                m=2 the analytic Jacobian agrees with the
                                `VectorBVP` Taylor1 autodiff Jacobian AND
                                an independent finite-difference Jacobian,
                                at several random (x, y) points, to
                                machine precision.  Three independent
                                computations of the same object.
    VW.3   failure modes      — m ≥ 3 and equation ≠ :I throw.
    VW.4   PainleveHierarchyProblem BVP convenience — construct an m=1
                                problem, call vector_bvp_solve(php, …),
                                assert the returned VectorBVPSolution
                                satisfies the BC and has a small residual,
                                and reproduces a known IVP trajectory.
    VW.5   mutation-proof     — recorded at end of file.

Oracles: the autodiff Jacobian `VectorBVP._autodiff_jacobian` and a
central finite-difference Jacobian (VW.2); an m=1 IVP trajectory whose
endpoints become the BVP Dirichlet data (VW.4).

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_bvp_wirein_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor

# Central finite-difference Jacobian — an independent third oracle for
# VW.2.  Step h = ∛eps is the standard central-difference optimum.
function _fd_jacobian(f, x, y)
    d  = length(y)
    h  = cbrt(eps(Float64))
    Jf = Matrix{Float64}(undef, d, d)
    for i in 1:d
        yp = copy(y); yp[i] += h
        ym = copy(y); ym[i] -= h
        col = (f(x, yp) .- f(x, ym)) ./ (2h)
        Jf[:, i] .= col
    end
    return Jf
end

@testset "VectorBVP wire-in + PainleveHierarchy Jacobian (VB3)" begin

    # -------------------------------------------------------------------------
    # VW.1 — umbrella wire-in.
    # -------------------------------------------------------------------------
    @testset "VW.1 VectorBVP wired into PadeTaylor" begin
        @test isdefined(PadeTaylor, :VectorBVP)
        @test isa(PadeTaylor.vector_bvp_solve, Function)
        @test PadeTaylor.VectorBVPSolution isa Type
        # Reachable unqualified after `using PadeTaylor` (exported).
        @test isa(vector_bvp_solve, Function)
        @test VectorBVPSolution isa Type
        @test isa(painleve_hierarchy_jacobian, Function)
    end

    # -------------------------------------------------------------------------
    # VW.2 — Jacobian correctness: analytic == autodiff == finite-difference.
    # The load-bearing test (CLAUDE.md Rule 4 cross-validation).  The
    # `painleve_hierarchy_jacobian` analytic helper and the `VectorBVP`
    # Taylor1 autodiff Jacobian are two independent computations of the
    # SAME object ∂f/∂y; a central finite-difference is a third.  They
    # must agree — autodiff to machine precision, FD to ∛eps.
    # -------------------------------------------------------------------------
    @testset "VW.2 analytic Jacobian == autodiff == finite-difference" begin
        @testset "VW.2.m=1" begin
            f  = painleve_hierarchy(:I, 1)
            Jf = painleve_hierarchy_jacobian(:I, 1)
            pts = [(0.3,  [0.4, -0.1]),
                   (-1.7, [1.2,  0.9]),
                   (2.5,  [-0.6, 2.1]),
                   (0.0,  [3.3, -4.4])]
            for (x, y) in pts
                Jan = Jf(x, y)
                Jad = PadeTaylor.VectorBVP._autodiff_jacobian(f, x, y, 2,
                                                              Float64)
                Jfd = _fd_jacobian(f, x, y)
                @test size(Jan) == (2, 2)
                # analytic vs autodiff — both exact: machine precision.
                @test Jan ≈ Jad atol = 1e-13 rtol = 1e-13
                # analytic vs finite-difference — FD-accurate.
                @test Jan ≈ Jfd atol = 1e-6 rtol = 1e-6
                # Pin the structural form: row 1 = identity shift [0 1].
                @test Jan[1, 1] == 0.0
                @test Jan[1, 2] == 1.0
                @test Jan[2, 2] == 0.0
                @test Jan[2, 1] ≈ 12 * y[1] atol = 1e-13
            end
        end

        @testset "VW.2.m=2" begin
            for t in (0.0, 0.7, -1.3)
                f  = painleve_hierarchy(:I, 2; t = t)
                Jf = painleve_hierarchy_jacobian(:I, 2; t = t)
                pts = [(0.5,  [1.0, -2.0, 0.3,  4.0]),
                       (-3.0, [0.2,  1.1, -0.7, 2.5]),
                       (2.0,  [-1.5, 0.0, 1.0, -0.5]),
                       (1.1,  [0.7,  3.2, -2.4, 0.9])]
                for (x, y) in pts
                    Jan = Jf(x, y)
                    Jad = PadeTaylor.VectorBVP._autodiff_jacobian(f, x, y,
                                                                  4, Float64)
                    Jfd = _fd_jacobian(f, x, y)
                    @test size(Jan) == (4, 4)
                    @test Jan ≈ Jad atol = 1e-12 rtol = 1e-12
                    @test Jan ≈ Jfd atol = 1e-5 rtol = 1e-5
                    # Structural pins: rows 1-3 are identity shifts.
                    @test Jan[1, :] == [0.0, 1.0, 0.0, 0.0]
                    @test Jan[2, :] == [0.0, 0.0, 1.0, 0.0]
                    @test Jan[3, :] == [0.0, 0.0, 0.0, 1.0]
                    # Row 4 — the ODE row, pinned to the hand derivation.
                    y1, y2, y3, _ = y
                    @test Jan[4, 1] ≈ -20*y3 - 120*y1^2 + 240*t atol = 1e-12
                    @test Jan[4, 2] ≈ -20*y2 atol = 1e-12
                    @test Jan[4, 3] ≈ -20*y1 atol = 1e-12
                    @test Jan[4, 4] == 0.0
                end
            end
        end
    end

    # -------------------------------------------------------------------------
    # VW.3 — failure modes (Rule 1): mirror painleve_hierarchy.
    # -------------------------------------------------------------------------
    @testset "VW.3 Jacobian failure modes" begin
        @test_throws ArgumentError painleve_hierarchy_jacobian(:I, 3)
        @test_throws ArgumentError painleve_hierarchy_jacobian(:I, 5)
        @test_throws ArgumentError painleve_hierarchy_jacobian(:I, 0)
        @test_throws ArgumentError painleve_hierarchy_jacobian(:II, 2)
        @test_throws ArgumentError painleve_hierarchy_jacobian(:VI, 1)
        # Wrong state-vector length at evaluation time must also throw.
        @test_throws ArgumentError painleve_hierarchy_jacobian(:I, 1)(0.0,
            [1.0, 2.0, 3.0])
        @test_throws ArgumentError painleve_hierarchy_jacobian(:I, 2)(0.0,
            [1.0, 2.0])
    end

    # -------------------------------------------------------------------------
    # VW.4 — PainleveHierarchyProblem BVP convenience.
    # Build an m=1 problem on a smooth interval; take Dirichlet boundary
    # data from a v0.1-consistent IVP trajectory, then BVP-solve with the
    # convenience method.  The BVP solution must (a) satisfy the BC, (b)
    # have a small Newton residual, (c) reproduce the IVP trajectory — a
    # known-correct oracle, not merely "didn't throw" (Rule 5).
    # -------------------------------------------------------------------------
    @testset "VW.4 PainleveHierarchyProblem BVP convenience" begin
        u0, up0 = 0.1, 0.05
        xspan   = (0.0, 0.6)

        # IVP oracle — the hierarchy m=1 trajectory (PI, u_xx = 6u^2 + x).
        php_ivp = PainleveHierarchyProblem(1; y0 = [u0, up0], xspan = xspan,
                                           order = 30)
        sv      = vector_solve_pade(php_ivp; h = 0.05)
        y_a     = sv(xspan[1])          # (u, u') at x_a
        y_b     = sv(xspan[2])          # (u, u') at x_b

        # Dirichlet two-point BC: y(z_a) = y_a, y(z_b) = y_b.  With
        # B_a = I, B_b = 0 this pins the z_a endpoint; B_b = I, B_a = 0
        # pins z_b.  The companion system is 2nd order, so pinning the
        # full 2-vector at one end (here z_a) plus one component at the
        # other is over-determined; the well-posed Dirichlet form for
        # u_xx = 6u^2 + x pins u at BOTH ends.  Encode that:
        #   B_a = [1 0; 0 0], B_b = [0 0; 1 0], g = [u(x_a), u(x_b)].
        B_a = [1.0 0.0; 0.0 0.0]
        B_b = [0.0 0.0; 1.0 0.0]
        g   = [y_a[1], y_b[1]]

        php = PainleveHierarchyProblem(1; y0 = [u0, up0], xspan = xspan,
                                       order = 30)
        # Convenience call — defaults jacobian to painleve_hierarchy_jacobian.
        sol = vector_bvp_solve(php, B_a, B_b, g; N = 24)

        @test sol isa VectorBVPSolution
        @test sol.residual_inf < 1e-9
        @test sol.iterations ≥ 1

        # (a) the BC is satisfied at the converged solution.
        ya_bvp = sol(xspan[1])
        yb_bvp = sol(xspan[2])
        bc = B_a * ya_bvp .+ B_b * yb_bvp .- g
        @test maximum(abs, bc) < 1e-9

        # (c) reproduces the IVP trajectory — the known-correct oracle.
        for x in 0.0:0.1:0.6
            @test sol(x)[1] ≈ sv(x)[1] atol = 1e-7 rtol = 1e-7
            @test sol(x)[2] ≈ sv(x)[2] atol = 1e-7 rtol = 1e-7
        end

        # The `jacobian = nothing` override forces the autodiff path; it
        # must converge to the SAME solution (the two Jacobians agree).
        sol_ad = vector_bvp_solve(php, B_a, B_b, g; N = 24, jacobian = nothing)
        for x in 0.0:0.15:0.6
            @test sol_ad(x)[1] ≈ sol(x)[1] atol = 1e-10
        end

        # Span-override kwargs route through to the core solver.
        sol_sub = vector_bvp_solve(php, B_a, B_b,
                                   [sv(0.1)[1], sv(0.5)[1]];
                                   z_a = 0.1, z_b = 0.5, N = 20)
        @test sol_sub.z_a ≈ 0.1
        @test sol_sub.z_b ≈ 0.5
        @test sol_sub(0.3)[1] ≈ sv(0.3)[1] atol = 1e-7
    end

end

# -----------------------------------------------------------------------------
# VW.5 — Mutation-proof record (CLAUDE.md Rule 4).
#
# Mutation applied to `src/PainleveHierarchy.jl::painleve_hierarchy_jacobian`,
# the suite re-run, the RED count recorded, then the source restored.
#
#   M1  — m=1 branch: `12 * y1` → `6 * y1` in the Jf[2,1] entry.
#         VW.2.m=1 went RED — the analytic-vs-autodiff assert, the
#         analytic-vs-FD assert, and the `Jan[2,1] ≈ 12*y[1]` structural
#         pin all fail at every one of the 4 random points (12 failing
#         assertions).  VW.2.m=2, VW.3, VW.4 stayed GREEN, exactly as
#         expected for an m=1-only mutation.  RED confirmed; source
#         restored to GREEN.
#
# The mutation proves VW.2 catches a wrong Jacobian entry: the analytic
# helper and the autodiff path are genuinely independent — perturbing one
# is detected by their disagreement, which is the cross-validation the
# load-bearing test exists to enforce.
# -----------------------------------------------------------------------------

# Standalone entry point: `julia --project=. test/vector_bvp_wirein_test.jl`.
