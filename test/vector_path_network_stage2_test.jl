"""
F2 tests for `PadeTaylor.VectorPathNetwork`'s **Stage-2 fine-grid fill**
(bead `padetaylor-0ln.20`; v0.2 plan row F2).

V7 (`vector_path_network_test.jl`) covers the Stage-1 visited-tree
walk.  This file covers the Stage-2 fine-grid fill: when
`vector_path_network_solve` is called with a `fine_grid` kwarg, the
solve evaluates each node's stored shared-`Q` Padé densely over that
grid, producing the `grid_z`/`grid_y` solution fields — a filled
heatmap of `y(z)`, the substrate for the `P_I⁽²⁾` tritronquée surface
figure.  Stage 2 is the direct vector lift of v0.1's
`PathNetwork.path_network_solve` Stage-2 loop (`src/PathNetwork.jl:669-693`);
the disc-radius `extrapolate` gate has the ADR-0015 semantics verbatim.

These tests (`VPN.2.*`) are written RED-first per CLAUDE.md Rule 4;
each asserts an invariant against a known-correct value (Rule 5):

    VPN.2.1  closed-form oracle  — the d=3 Riccati system y_i' = -y_i²/c_i
                                   has the EXACT solution y_i(z) = c_i/(z−p).
                                   Stage-2-fill a fine grid; assert every
                                   filled cell inside a node's disc matches
                                   the closed form to local Padé accuracy.
    VPN.2.2  harmonic oscillator — the d=2 pole-free system y' = (y₂, −y₁)
                                   with IC (0,1) has the exact solution
                                   (sin z, cos z).  Stage-2-fill a real
                                   sub-interval; assert the fill matches
                                   (sin, cos).  A pole-free cross-check —
                                   no disc-coverage gaps to reason around.
    VPN.2.3  fail-soft / extrapolate — a grid point far from every visited
                                   node returns a d-vector of NaN (no
                                   throw); with extrapolate=true the same
                                   point evaluates to a finite vector.
    VPN.2.4  node-consistency    — Stage-2 at a grid point equal to a
                                   visited node reproduces that node's
                                   `visited_y` (t = 0 ⇒ Pᵢ(0)/Q(0)).
    VPN.2.5  fail-fast / empties — empty `fine_grid` throws; no `fine_grid`
                                   leaves the Stage-2 fields empty.

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_path_network_stage2_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution

# -----------------------------------------------------------------------------
# Oracle A — the d-component Riccati pole system (the V7 oracle, reused).
#
# y_i' = -y_i²/c_i  has the exact solution  y_i(z) = c_i/(z-p).
# Given IC y_i(z0) = c_i/(z0-p), the pole is p = z0 - c_i/y_i(z0).
# This is the figure-relevant case: a pole-rich vector meromorphic
# solution whose components all blow up at the SAME z = p.
# -----------------------------------------------------------------------------
function riccati_pole_problem_s2(p, z0, cs; order = 24)
    f  = (z, y) -> [-(y[i]^2) / cs[i] for i in eachindex(y)]
    y0 = [cs[i] / (z0 - p) for i in eachindex(cs)]
    return VectorPadeTaylorProblem(f, y0, (z0, z0 + 1.0 + 0im); order = order)
end
riccati_exact(z, p, cs) = [cs[i] / (z - p) for i in eachindex(cs)]

# -----------------------------------------------------------------------------
# Oracle B — the d=2 harmonic-oscillator system y' = (y₂, -y₁).
#
# With IC y(0) = (0, 1) the exact solution is y(z) = (sin z, cos z):
#   y₁' = y₂ = cos z = (sin z)',   y₂' = -y₁ = -sin z = (cos z)'.
# Pole-free — the whole plane is in some node's disc of validity, so
# this oracle has no fail-soft gaps to reason around.
# -----------------------------------------------------------------------------
function harmonic_problem(z0, z1; order = 24)
    f  = (z, y) -> [y[2], -y[1]]
    y0 = ComplexF64[sin(z0), cos(z0)]
    return VectorPadeTaylorProblem(f, y0, (z0, z1); order = order)
end
harmonic_exact(z) = ComplexF64[sin(z), cos(z)]

@testset "VectorPathNetwork Stage-2 fill (VPN.2.*)" begin

    # -------------------------------------------------------------------------
    # VPN.2.1 — closed-form oracle.  Riccati y_i(z) = c_i/(z-p); the
    # Stage-2 fill must reproduce it to local Padé accuracy on every
    # grid cell that lands inside some visited node's disc.
    # -------------------------------------------------------------------------
    @testset "VPN.2.1 closed-form Riccati oracle" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem_s2(p, z0, cs; order = 24)

        targets = ComplexF64[x + y * im
                             for x in 0.2:0.4:1.8 for y in -0.2:0.4:1.4]
        # The fine grid: a denser lattice over the same region.
        fine = ComplexF64[x + y * im
                          for x in 0.0:0.2:1.6 for y in -0.2:0.2:1.2]

        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                        fine_grid = fine)

        # The Stage-2 fields are populated and shaped correctly.
        @test length(sol.grid_z) == length(fine)
        @test length(sol.grid_y) == length(fine)
        @test sol.grid_z == fine

        # Every covered cell (a finite, non-NaN fill) must match the
        # closed form.  The min-‖y‖ walk steers around the pole p, so
        # cells very close to p may be NaN (uncovered) or large-error;
        # we assert accuracy only on covered cells comfortably away
        # from the singularity.
        n_checked = 0
        for (zf, yf) in zip(sol.grid_z, sol.grid_y)
            all(isfinite, yf) || continue          # skip fail-soft NaN cells
            abs(zf - p) > 0.35 || continue         # skip near-pole cells
            exact = riccati_exact(zf, p, cs)
            err = maximum(abs.(yf .- exact))
            @test err < 1.0e-6
            n_checked += 1
        end
        # The fill must actually have covered a substantial number of
        # cells — otherwise "every covered cell matches" is vacuous.
        @test n_checked ≥ 20
    end

    # -------------------------------------------------------------------------
    # VPN.2.2 — harmonic-oscillator oracle.  Pole-free: (sin, cos) is the
    # exact solution; every Stage-2 cell is covered and must match it.
    # -------------------------------------------------------------------------
    @testset "VPN.2.2 harmonic-oscillator oracle" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)

        # Targets + fine grid along the real axis (the solution is
        # entire, so a real-axis walk is well-posed).
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.8]
        fine    = ComplexF64[x + 0.0im for x in 0.1:0.1:1.7]

        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                        fine_grid = fine)

        @test length(sol.grid_y) == length(fine)
        # Pole-free ⇒ every cell within h of a node ⇒ no NaN gaps for
        # this h-and-target choice; assert the fill matches (sin, cos).
        n_checked = 0
        for (zf, yf) in zip(sol.grid_z, sol.grid_y)
            all(isfinite, yf) || continue
            exact = harmonic_exact(zf)
            @test abs(yf[1] - exact[1]) < 1.0e-8       # sin component
            @test abs(yf[2] - exact[2]) < 1.0e-8       # cos component
            n_checked += 1
        end
        @test n_checked ≥ 14
    end

    # -------------------------------------------------------------------------
    # VPN.2.3 — fail-soft NaN vs extrapolate=true.
    # -------------------------------------------------------------------------
    @testset "VPN.2.3 fail-soft / extrapolate" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.5]

        # A grid point planted FAR from any plausible visited node — the
        # walk never goes near 12+12im, so this cell is uncovered.
        far  = 12.0 + 12.0im
        near = 0.6 + 0.0im                  # comfortably inside a node disc
        fine = ComplexF64[near, far]

        # Default (extrapolate=false): the far cell is a d-vector of NaN,
        # the near cell is finite — fail-soft, NOT a throw (ADR-0015 / Rule 1).
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                        fine_grid = fine)
        @test all(isfinite, sol.grid_y[1])             # near: covered
        @test all(z -> isnan(real(z)) && isnan(imag(z)), sol.grid_y[2])
        @test length(sol.grid_y[2]) == 2               # still a d-vector

        # extrapolate=true: the disc-radius check is skipped, so the far
        # cell now evaluates to a finite (if inaccurate) vector.
        sol_x = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                          fine_grid = fine, extrapolate = true)
        @test all(isfinite, sol_x.grid_y[2])
        # Inside-disc cells are identical between the two modes — the
        # gate only changes behaviour for out-of-disc cells.
        @test sol.grid_y[1] == sol_x.grid_y[1]
    end

    # -------------------------------------------------------------------------
    # VPN.2.4 — node-consistency.  A grid point placed exactly on a
    # visited node must reproduce that node's `visited_y`: t = 0, and
    # the canonical shared-Q approximant satisfies Pᵢ(0)/Q(0) = yᵢ.
    # -------------------------------------------------------------------------
    @testset "VPN.2.4 node-consistency" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.5]

        # First walk with no fine grid to learn the visited nodes.
        sol0 = vector_path_network_solve(prob, targets; order = 24, h = 0.3)
        @test isempty(sol0.grid_z)                     # Stage-2 fields empty
        @test isempty(sol0.grid_y)

        # Re-solve with the visited-node positions themselves as the
        # fine grid.  Each grid point z = visited_z[k] is its own
        # nearest node ⇒ t = 0 ⇒ the fill must equal visited_y[k].
        node_grid = copy(sol0.visited_z)
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                        fine_grid = node_grid)
        @test length(sol.grid_y) == length(sol.visited_z)
        for k in eachindex(sol.visited_z)
            yk    = sol.visited_y[k]
            fill_k = sol.grid_y[k]
            @test maximum(abs.(fill_k .- yk)) < 1.0e-10
        end
    end

    # -------------------------------------------------------------------------
    # VPN.2.5 — fail-fast on an empty fine_grid; no-fine_grid invariant.
    # -------------------------------------------------------------------------
    @testset "VPN.2.5 fail-fast / no-grid invariant" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.5]

        # An empty fine_grid is a caller mistake — pass nothing to skip.
        @test_throws ArgumentError vector_path_network_solve(
            prob, targets; order = 24, h = 0.3, fine_grid = ComplexF64[])

        # No fine_grid: the Stage-2 fields stay empty, the V7 scatter
        # behaviour is byte-unaffected.
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3)
        @test sol.grid_z == ComplexF64[]
        @test sol.grid_y == Vector{ComplexF64}[]

        # The 6-arg backward-compat constructor also yields empty Stage-2
        # fields (the hand-built-fixture path).
        sol6 = VectorPathNetworkSolution{Float64}(
            ComplexF64[1.0], [ComplexF64[1.0]], ComplexF64[0.5],
            [[ComplexF64[1.0]]], [ComplexF64[1.0]], [0])
        @test isempty(sol6.grid_z)
        @test isempty(sol6.grid_y)
    end

end

# =============================================================================
# Mutation-proof record (CLAUDE.md Rule 4).
#
# Procedure: perturb `src/VectorPathNetworkStage2.jl`, rerun this file,
# confirm RED, restore.  2 mutations bit; 1 false-trail mutation is
# recorded because the "why it doesn't bite" is itself instructive.
#
#   M1 — drop the 1/h_v rescale.  In `_stage2_fill`, change
#          t = (z_f - z_v) / h_v
#        to
#          t = (z_f - z_v)
#        (evaluate the canonical Padé at the raw z-displacement instead
#        of the rescaled variable t = (z - z_node)/h the approximant
#        actually lives in).
#        Expected: every covered off-node Stage-2 cell evaluates the
#        Padé at the wrong argument, so the closed-form-oracle
#        assertions fail.
#        Result: RED — VPN.2.1 (Riccati err < 1e-6) and VPN.2.2
#        (harmonic sin/cos err < 1e-8) bit on essentially every covered
#        cell.  VPN.2.4 node-consistency correctly STILL passed: at a
#        node z_f = z_v the displacement is 0 either way (t = 0
#        unchanged), so M1 specifically corrupts only off-node cells —
#        exactly the load-bearing 1/h_v factor's job.  Restored to GREEN.
#
#   M2-false — double every Horner coefficient.  In `_eval_poly`, change
#          s = s * t + c[k]
#        to
#          s = s * t + 2 * c[k]
#        Expected (naïvely): the polynomial value is wrong, so the
#        shared-Q evaluation Pᵢ(t)/Q(t) is wrong.
#        Result: GREEN — NO test bit.  Doubling every coefficient
#        computes 2·Pᵢ(t) for the numerator AND 2·Q(t) for the
#        denominator; the factor of 2 cancels exactly in the ratio
#        2·Pᵢ(t)/(2·Q(t)) = Pᵢ(t)/Q(t).  A scalar multiple of `_eval_poly`
#        is invisible to a shared-Q quotient.  Recorded as a false trail:
#        a Horner mutation only bites if it changes the polynomial
#        non-multiplicatively — see M2 below.
#
#   M2 — reverse the Horner loop direction.  In `_eval_poly`, change
#          for k in length(c):-1:1
#        to
#          for k in 1:length(c)
#        (Horner requires high-to-low coefficient order; iterating
#        low-to-high evaluates the *reversed* polynomial tⁿ·P(1/t)).
#        Expected: numerator and denominator are reversed polynomials of
#        different degrees, so the factor does NOT cancel in the ratio —
#        the shared-Q evaluation is genuinely wrong.
#        Result: RED — 99 of 118 assertions bit: VPN.2.1, VPN.2.2 and
#        VPN.2.4 all failed (the reversed evaluation breaks the oracle
#        match AND Pᵢ(0)/Q(0) = yᵢ at the node).  Restored to GREEN.
#
# Certified bites: M1 → VPN.2.1 + VPN.2.2 RED (VPN.2.4 correctly
# unaffected — it probes t = 0 only, where the missing 1/h_v is moot);
# M2 → VPN.2.1 + VPN.2.2 + VPN.2.4 RED (99 assertions).  Restored to
# GREEN after each mutation.
# =============================================================================
