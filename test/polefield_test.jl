# test/polefield_test.jl -- bead `padetaylor-xvf` tests.
#
# PoleField.extract_poles: recover pole locations in the z-plane from a
# solved path-network's per-node Padé store.  This is the capability
# behind FW 2011 Fig 4.7 / 4.8 (pole-location scatter plots).
#
# Ground-truth oracle — the equianharmonic Weierstrass-℘ test problem
# of FW 2011 §5.1.1 (references/markdown/FW2011_painleve_methodology_
# JCP230/FW2011_painleve_methodology_JCP230.md:281-318):
#
#     u''(z) = 6 u(z)^2,        u(z) = ℘(z + c₁; 0, c₂),   c₁ = -1, c₂ = 2.
#
# Its solution has *analytically known* second-order poles on a rhombic
# lattice of equilateral triangles.  FW md:297 gives the real-axis
# half-period ω = Γ(1/3)³ / (2^{13/6} π) ≈ 1.363; with c₁ = -1 the
# poles of u(z) sit at
#
#     z(m, n) = 1 + 2ω(m + n/2) + i · ω√3 · n,    m, n ∈ ℤ.
#
# So the extracted poles can be checked against an exact closed form —
# no pole-finder oracle is needed (this matches `figures/fw2011_fig_5_1.jl`,
# which makes the same observation).
#
# Tests: PF.1.1 (single-node extraction surfaces the nearest pole),
# PF.1.2 (full pole field over a 2D grid — no spurious poles, all
# in-region lattice poles recovered, conjugate symmetry), PF.2.1
# (degenerate-network edge case), PF.2.2 (the cross-node-support filter
# is load-bearing).  Mutation-proof procedure is in the file footer.

using Test
using PadeTaylor
using PadeTaylor.PathNetwork: PathNetworkSolution

include(joinpath(@__DIR__, "_oracle_problems.jl"))

# Equianharmonic ℘ real half-period, FW md:297.  Hard-coded to Float64
# precision rather than pulling in SpecialFunctions as a test dep;
# reproduce with `gamma(1/3)^3 / (2^(13/6) * π)`.
const Ω_FW = 1.3630340904278908

# Exact pole lattice of u(z) = ℘(z - 1; 0, 2): z(m,n) per FW md:297.
℘_pole(m::Int, n::Int) = 1 + 2Ω_FW * (m + n / 2) + im * (Ω_FW * sqrt(3) * n)

# All lattice poles with |m| ≤ M, |n| ≤ N — the reference set the
# spurious-pole check measures against.
℘_lattice(M::Int, N::Int) = [℘_pole(m, n) for m in -M:M for n in -N:N]

# Distance from `z` to the nearest exact lattice pole.
nearest_lattice_dist(z, lattice) = minimum(abs(z - p) for p in lattice)

# Is `z` inside the closed box [xlo,xhi] × [ylo,yhi]?  `extract_poles`
# honestly reports *every* pole ≥ min_support nodes agree on, including
# genuine lattice poles just outside the explored grid that only distant
# nodes can see — those are correctly identified as real, but placed
# less accurately (FW's Fig 4.7 likewise clips to a display window).
# Accuracy assertions are therefore made over the grid's covered box.
in_box(z, xlo, xhi, ylo, yhi) =
    xlo ≤ real(z) ≤ xhi && ylo ≤ imag(z) ≤ yhi

@testset "PoleField (bead padetaylor-xvf): FW 2011 Fig 4.7/4.8 capability" begin

    # Test ODE: u'' = 6u^2, FW 2011 §5.1.1 ICs at z = 0.
    fW(z, u, up) = 6 * u^2

    @testset "PF.1.1: single-node network surfaces the nearest pole" begin
        # A grid wholly inside |z| ≤ h = 0.5: Stage 1 takes no steps,
        # so the only visited node is the IC.  With the cross-node
        # support filter disabled (`min_support = 1`) the IC node's
        # own Padé is read directly.  Its (15,15) denominator places
        # the nearest pole — z = 1, at rescaled t = 1/h = 2 — to the
        # double-pole Padé accuracy floor (~5e-7 here), within the
        # figure-catalogue Float64 spec of 1e-6.  (A single far node
        # cannot place *distant* poles well — that is exactly what the
        # multi-node support filter exercised in PF.1.2 is for.)
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[0.0, 0.2 + 0.1im, -0.15 + 0.2im, 0.25 - 0.2im]
        sol  = path_network_solve(prob, grid; h = 0.5)

        @test length(sol.visited_z) == 1          # IC only — no Stage-1 steps

        poles = extract_poles(sol; min_support = 1)
        @test !isempty(poles)

        err_z1 = minimum(abs(p - (1.0 + 0im)) for p in poles)
        @test err_z1 ≤ 1.0e-6
    end

    @testset "PF.1.2: full pole field over a 2D grid" begin
        # A 2D grid covering several lattice poles.  Offsets chosen so
        # no grid target sits exactly on a pole.  Every in-region pole
        # is then seen by many independent visited nodes, so the
        # default cross-node support filter keeps the physical poles
        # and rejects node-local artefacts.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        xs = -2.4:0.5:4.1
        ys = -2.6:0.5:2.6
        grid = ComplexF64[x + im * y for x in xs for y in ys]
        sol  = path_network_solve(prob, grid; h = 0.5)

        poles = extract_poles(sol)
        @test !isempty(poles)

        # (a) No spurious poles: every extracted pole *inside the covered
        #     grid box* lies within 1e-6 of an exact lattice pole.
        #     (Best-placed estimates are ~1e-8 in practice; the 1e-6 pin
        #     is the figure-catalogue Float64 spec.)
        lattice = ℘_lattice(6, 3)
        in_grid = filter(p -> in_box(p, -2.4, 4.1, -2.6, 2.6), poles)
        @test !isempty(in_grid)
        worst = maximum(nearest_lattice_dist(p, lattice) for p in in_grid)
        @test worst ≤ 1.0e-6

        # (b) Completeness: the poles unambiguously inside the covered
        #     region must all be recovered, each to ≤ 1e-6.
        expected = ComplexF64[
            ℘_pole(0, 0),                       #  1            (real axis)
            ℘_pole(1, 0),                       #  3.726        (real axis)
            ℘_pole(-1, 0),                      # -1.726        (real axis)
            ℘_pole(0, 1),  ℘_pole(0, -1),       #  2.363 ± 2.36i
            ℘_pole(-1, 1), ℘_pole(-1, -1),      # -0.363 ± 2.36i
        ]
        for ze in expected
            err = minimum(abs(p - ze) for p in poles)
            @test err ≤ 1.0e-6
        end

        # (c) Conjugate symmetry: the IC is real, so the pole field is
        #     symmetric across the real axis.  Every off-axis extracted
        #     pole inside the covered box has a conjugate partner.
        for p in in_grid
            if abs(imag(p)) > 1.0e-3
                partner = minimum(abs(conj(p) - q) for q in poles)
                @test partner ≤ 1.0e-6
            end
        end
    end

    @testset "PF.2.1: degenerate network (constant denominators) → no poles" begin
        # A hand-built network whose every stored Padé has a constant
        # denominator b = [1] (a pure polynomial — no poles).  extract_poles
        # must skip such nodes cleanly and return an empty vector, not
        # throw on the empty `roots` call.
        T = Float64
        flat = PadeApproximant{Complex{T}}(
            Complex{T}[1.0, 0.5], Complex{T}[1.0], 1, 0, one(Complex{T}))
        sol = PathNetworkSolution{T}(
            Complex{T}[0.0, 0.5],            # visited_z
            Complex{T}[1.0, 1.2],            # visited_u
            Complex{T}[0.0, 0.1],            # visited_up
            [flat, flat],                    # visited_pade
            T[0.5, 0.5],                     # visited_h
            Int[0, 1],                       # visited_parent
            Complex{T}[],                    # grid_z
            Complex{T}[],                    # grid_u
            Complex{T}[],                    # grid_up
        )
        @test extract_poles(sol) == Complex{T}[]
        @test extract_poles(sol; min_support = 1) == Complex{T}[]
    end

    @testset "PF.2.2: cross-node support filter is load-bearing" begin
        # On the same 2D-grid network, disabling the support filter
        # (`min_support = 1`) lets node-local artefacts through, so it
        # must yield strictly more "poles" than the default — and those
        # extras must be the off-lattice spurious ones.  This pins the
        # filter as doing real work, not decoration.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        xs = -2.4:0.5:4.1
        ys = -2.6:0.5:2.6
        grid = ComplexF64[x + im * y for x in xs for y in ys]
        sol  = path_network_solve(prob, grid; h = 0.5)

        clean    = extract_poles(sol)                       # default filter
        unfilt   = extract_poles(sol; min_support = 1)      # filter disabled
        @test length(unfilt) > length(clean)

        lattice = ℘_lattice(6, 3)
        # The default result is on-lattice within the covered box; the
        # unfiltered result drags in node-local off-lattice artefacts.
        clean_in = filter(p -> in_box(p, -2.4, 4.1, -2.6, 2.6), clean)
        @test maximum(nearest_lattice_dist(p, lattice) for p in clean_in) ≤ 1.0e-6
        @test maximum(nearest_lattice_dist(p, lattice) for p in unfilt)   > 1.0e-6
    end
end

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Each mutation below was applied to `src/PoleField.jl`, the suite
# confirmed RED, then the mutation reverted.
#
#   M1 — root → z mapping.  Change `z_v + h * t` to `z_v + t` (drop the
#        canonical-step rescale).  PF.1.1's `err_z1 ≤ 1e-6` goes RED
#        (z = 1 maps to t = 2), and PF.1.2 / PF.2.2 go RED — every pole
#        lands off the lattice.  (Observed: 3 failed, 2 errored.)
#
#   M2 — radius_t filter.  Comment out `abs(t) ≤ radius || continue`.
#        Far extrapolation roots leak in, corrupting the greedy-|t*|
#        clustering: PF.1.2 (completeness + symmetry) and PF.2.2 go RED.
#        (Observed: 43 failed.)
#
#   M3 — cross-node support filter.  Replace the final comprehension
#        guard `length(support[j]) ≥ min_support` with `≥ 1`.  PF.1.2
#        (a) goes RED (node-local artefacts reported as in-box poles)
#        and PF.2.2 goes RED (`length(unfilt) == length(clean)`, and the
#        default result is no longer on-lattice).  (Observed: 49 failed.)
# ----------------------------------------------------------------------
