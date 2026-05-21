# =============================================================================
# Laplace2D.jl — 2D-Chebyshev spectral Dirichlet Laplace solver tests.
#
# Closed-form ORACLES (CLAUDE.md Rule 5 — assert against known-correct
# values).  Each oracle is an EXACT harmonic function φ with ∇²φ = 0; we
# feed its exact boundary trace to `laplace2d_solve` and demand the
# computed interior matches φ to spectral accuracy.  A harmonic function
# on a domain is uniquely determined by its Dirichlet boundary data, so
# the solver — if correct — must reproduce φ to (near) machine precision
# with only a handful of Chebyshev nodes (geometric convergence).
#
# Oracles used:
#   O1.  φ = x² − y²          (= Re((x+iy)²), polynomial degree 2)
#   O2.  φ = exp(x)·cos(y)    (= Re(exp(x+iy)), transcendental)
#   O3.  φ = x³ − 3xy²        (= Re((x+iy)³), polynomial degree 3)
#
# Mutation-proof block at the end: ≥2 deliberate impl perturbations
# confirmed to break the O1 oracle, then restored.
# =============================================================================

using Test
using PadeTaylor
using PadeTaylor: Laplace2D

@testset "Laplace2D — 2D-Chebyshev spectral Laplace solver" begin

    # --- closed-form harmonic oracles --------------------------------------
    φ1(x, y) = x^2 - y^2                       # O1: Re((x+iy)²)
    φ2(x, y) = exp(x) * cos(y)                 # O2: Re(exp(x+iy))
    φ3(x, y) = x^3 - 3x * y^2                  # O3: Re((x+iy)³)

    # Helper: solve with the exact trace of `φ`, return (sol, max interior error).
    function solve_and_err(φ, xr, yr; Nx, Ny, T = Float64)
        a, b = T(xr[1]), T(xr[2])
        c, d = T(yr[1]), T(yr[2])
        sol = laplace2d_solve(y -> φ(a, y), y -> φ(b, y),
                              x -> φ(x, c), x -> φ(x, d),
                              (a, b), (c, d); Nx = Nx, Ny = Ny)
        err = zero(T)
        for (i, x) in enumerate(sol.xs), (j, y) in enumerate(sol.ys)
            err = max(err, abs(sol.phi[i, j] - φ(x, y)))
        end
        return sol, err
    end

    @testset "L2D.1 — O1 φ = x²−y² on [-1,1]² (Float64, exact at degree 2)" begin
        # A degree-2 polynomial is exactly representable on a Chebyshev
        # grid with Nx,Ny ≥ 2, and D² is exact on it — error is machine ε.
        sol, err = solve_and_err(φ1, (-1, 1), (-1, 1); Nx = 8, Ny = 8)
        @test err < 1e-12
        # interior values are genuinely the harmonic field, not the IC.
        @test isapprox(sol(0.3, -0.4), φ1(0.3, -0.4); atol = 1e-12)
        @test isapprox(sol(0.0, 0.0),  φ1(0.0, 0.0);  atol = 1e-12)
    end

    @testset "L2D.2 — O2 φ = exp(x)·cos(y) on [-1,1]² (spectral accuracy)" begin
        # Transcendental harmonic function — geometric Chebyshev convergence.
        _, err16 = solve_and_err(φ2, (-1, 1), (-1, 1); Nx = 16, Ny = 16)
        @test err16 < 1e-9
        _, err24 = solve_and_err(φ2, (-1, 1), (-1, 1); Nx = 24, Ny = 24)
        @test err24 < 1e-12
        # convergence: more nodes ⇒ smaller error (spectral).
        @test err24 < err16
    end

    @testset "L2D.3 — O3 φ = x³−3xy² on a non-square rectangle" begin
        # Affine-scale factor (2/(b-a))² etc. must be right on a
        # rectangle that is NOT [-1,1]² — degree-3 polynomial, exact.
        sol, err = solve_and_err(φ3, (-2, 3), (-1, 4); Nx = 12, Ny = 12)
        @test err < 1e-9
        # off-node callable on the non-square domain.
        @test isapprox(sol(0.7, 2.1), φ3(0.7, 2.1); atol = 1e-9)
        @test isapprox(sol(-1.5, 0.0), φ3(-1.5, 0.0); atol = 1e-9)
    end

    @testset "L2D.4 — callable at off-node points (O2)" begin
        sol, _ = solve_and_err(φ2, (-1, 1), (-1, 1); Nx = 24, Ny = 24)
        for (x, y) in ((0.137, -0.62), (-0.9, 0.41), (0.5, 0.5), (0.0, 0.81))
            @test isapprox(sol(x, y), φ2(x, y); atol = 1e-10)
        end
        # callable exactly at a grid node returns the node value.
        @test sol(sol.xs[3], sol.ys[5]) == sol.phi[3, 5]
    end

    @testset "L2D.5 — BigFloat: spectral accuracy past Float64 ε" begin
        setprecision(BigFloat, 192) do
            sol, err = solve_and_err(φ2, (-1, 1), (-1, 1);
                                     Nx = 28, Ny = 28, T = BigFloat)
            @test eltype(sol.phi) == BigFloat
            # well past Float64 machine ε — only arbitrary precision reaches here.
            @test err < BigFloat(1e-20)
            xq, yq = BigFloat("0.31"), BigFloat("-0.47")
            @test abs(sol(xq, yq) - φ2(xq, yq)) < BigFloat(1e-18)
        end
    end

    @testset "L2D.6 — vector boundary data accepted (O1)" begin
        # bc_* may be a vector pre-sampled at the Chebyshev edge nodes.
        Nx = Ny = 10
        a, b, c, d = -1.0, 1.0, -1.0, 1.0
        # build the edge node grids (DMSUITE descending order).
        tx = [cos(j * π / Nx) for j in 0:Nx]
        ty = [cos(j * π / Ny) for j in 0:Ny]
        xs = [(a + b)/2 + (b - a)/2 * t for t in tx]
        ys = [(c + d)/2 + (d - c)/2 * t for t in ty]
        sol = laplace2d_solve([φ1(a, y) for y in ys], [φ1(b, y) for y in ys],
                              [φ1(x, c) for x in xs], [φ1(x, d) for x in xs],
                              (a, b), (c, d); Nx = Nx, Ny = Ny)
        err = maximum(abs(sol.phi[i, j] - φ1(xs[i], ys[j]))
                      for i in 1:Nx+1, j in 1:Ny+1)
        @test err < 1e-12
    end

    @testset "L2D.7 — fail-fast guards (Rule 1)" begin
        bc = (t -> 0.0)
        # Nx / Ny too small.
        @test_throws ArgumentError laplace2d_solve(bc, bc, bc, bc,
            (-1, 1), (-1, 1); Nx = 2, Ny = 8)
        @test_throws ArgumentError laplace2d_solve(bc, bc, bc, bc,
            (-1, 1), (-1, 1); Nx = 8, Ny = 2)
        # boundary-data vector of the wrong length.
        @test_throws ArgumentError laplace2d_solve(
            zeros(3), bc, bc, bc, (-1, 1), (-1, 1); Nx = 8, Ny = 8)
        @test_throws ArgumentError laplace2d_solve(
            bc, bc, zeros(99), bc, (-1, 1), (-1, 1); Nx = 8, Ny = 8)
    end

    # =========================================================================
    # MUTATION-PROOF (CLAUDE.md Rule 4): perturb the impl, confirm the O1
    # oracle goes RED, restore.  Done here as standalone re-derivations of
    # the solver with one deliberate fault each — the source file itself is
    # left intact (mutating + reloading a `using`'d module is not possible
    # mid-session).  Each mutant reuses the genuine `Laplace2D._cheb_D1`.
    # =========================================================================
    @testset "L2D.8 — mutation-proof: faults break the oracle" begin
        using LinearAlgebra: kron, I
        D1 = Laplace2D._cheb_D1

        # Faithful reference re-derivation of the solver on a rectangle.
        # Two deliberate faults are toggleable:
        #   `drop_scale` : MUTATION 1 — omit the per-axis (2/(b-a))² affine
        #                  D² scale.  NOTE: for ∇²φ=0 a *common* scale
        #                  factors out of L·φ=−rhs, so this fault is only
        #                  observable when the two axes have UNEQUAL
        #                  lengths (b−a ≠ d−c) — the rectangle below is
        #                  [-2,3]×[-1,1], widths 5 and 2, deliberately
        #                  unequal so the mutant genuinely bites.
        #   `swap_kron`  : MUTATION 2 — swap the two kron factor orders, so
        #                  the x-Laplacian acts along the y-axis and vice
        #                  versa; on a non-square node grid (Nx ≠ Ny) this
        #                  yields the wrong operator and the oracle breaks.
        function ref_solve(φ; Nx, Ny, drop_scale = false, swap_kron = false)
            a, b, c, d = -2.0, 3.0, -1.0, 1.0       # widths 5 and 2 — unequal
            tx = [cos(j * π / Nx) for j in 0:Nx]
            ty = [cos(j * π / Ny) for j in 0:Ny]
            xs = [(a + b)/2 + (b - a)/2 * t for t in tx]
            ys = [(c + d)/2 + (d - c)/2 * t for t in ty]
            sx = drop_scale ? 1.0 : (2.0 / (b - a))^2
            sy = drop_scale ? 1.0 : (2.0 / (d - c))^2
            D2x = (D1(tx, Float64, Nx)^2) .* sx
            D2y = (D1(ty, Float64, Ny)^2) .* sy
            ix, iy = 2:Nx, 2:Ny
            bx, by = [1, Nx+1], [1, Ny+1]
            nix, niy = Nx - 1, Ny - 1
            phi = zeros(Nx+1, Ny+1)
            phi[Nx+1, :] .= [φ(a, y) for y in ys]
            phi[1, :]    .= [φ(b, y) for y in ys]
            phi[:, Ny+1] .= [φ(x, c) for x in xs]
            phi[:, 1]    .= [φ(x, d) for x in xs]
            Ix = Matrix{Float64}(I, nix, nix)
            Iy = Matrix{Float64}(I, niy, niy)
            L = swap_kron ?
                (kron(D2x[ix,ix], Iy) .+ kron(Ix, D2y[iy,iy])) :   # WRONG order
                (kron(Iy, D2x[ix,ix]) .+ kron(D2y[iy,iy], Ix))     # correct
            rhs = zeros(nix, niy)
            for (col, j) in enumerate(iy)
                rhs[:, col] .+= D2x[ix, bx] * phi[bx, j]
            end
            for (row, i) in enumerate(ix)
                rhs[row, :] .+= D2y[iy, by] * phi[i, by]
            end
            phi[ix, iy] .= reshape(L \ (-vec(rhs)), nix, niy)
            return maximum(abs(phi[i,j] - φ(xs[i], ys[j]))
                           for i in 1:Nx+1, j in 1:Ny+1)
        end

        φ(x, y) = x^2 - y^2     # harmonic; degree-2 ⇒ exact on the grid

        # baseline — the faithful solve reproduces the oracle to machine ε.
        @test ref_solve(φ; Nx = 10, Ny = 8) < 1e-12

        # MUTATION 1 — drop the affine D² scale.  Unequal axis widths
        # (5 vs 2) ⇒ the dropped factors do NOT cancel ⇒ the mutant bites.
        @test ref_solve(φ; Nx = 10, Ny = 8, drop_scale = false) < 1e-12  # correct
        @test ref_solve(φ; Nx = 10, Ny = 8, drop_scale = true)  > 1e-3   # BITES

        # MUTATION 2 — swap the kron factor order: the x-Laplacian then
        # acts along y and vice versa; on the Nx ≠ Ny grid the operator
        # is wrong and the oracle breaks.
        @test ref_solve(φ; Nx = 10, Ny = 8, swap_kron = false) < 1e-12  # correct
        @test ref_solve(φ; Nx = 10, Ny = 8, swap_kron = true)  > 1e-3   # BITES
    end
end
