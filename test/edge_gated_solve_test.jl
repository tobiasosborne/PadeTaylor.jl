# test/edge_gated_solve_test.jl -- bead `padetaylor-dmb` tests.
#
# EdgeGatedSolve.edge_gated_pole_field_solve: confine the IVP path-network
# to the pole field by region growing, so it never wanders into (and is
# corrupted by) the smooth regions.  FW 2011 md:401 — "smooth regions are
# unstable regions for any IVP solver ... force the path selection
# algorithm to complete pole fields before stepping into smooth regions
# ... with the pole field edge detection procedure (3.2.2) this process
# can readily be fully automated."
#
# Oracle — the PI tritronquée.  The Boutroux tritronquée of
#
#     u''(z) = 6 u(z)^2 + z,    u(0) = -0.1875543083404949,
#                               u'(0) = 0.3049055602612289   (FW eq. 4.1)
#
# is pole-free in four of its five 2π/5 sectors; poles populate only the
# sector around the positive real axis (FW §2 lines 52-54, §4.5).  A
# plain `path_network_solve` over a large grid corrupts this — accumulated
# IVP error perturbs the solution off the tritronquée manifold and
# spurious poles bloom in the four pole-free sectors.  The edge-gated
# solve must NOT: its extracted pole field stays sector-confined.
#
# Tests: EG.1.1 (the solve runs and actually grew a region), EG.1.2
# (sector confinement — the four pole-free sectors stay near-empty),
# EG.1.3 (the load-bearing contrast: the plain ungated solve over the
# same grid DOES leak poles into those sectors), EG.2.1 (input
# validation).  Mutation-proof procedure in the file footer.

using Test
using PadeTaylor
using PadeTaylor.EdgeGatedSolve: edge_gated_pole_field_solve

@testset "EdgeGatedSolve (bead padetaylor-dmb): FW md:401 edge-gated IVP" begin

    # PI equation; FW eq. (4.1) Maple-computed tritronquée ICs.
    pI(z, u, up) = 6 * u^2 + z
    u_tri  = -0.1875543083404949
    up_tri =  0.3049055602612289

    # Moderate grid over [-12,12]² — large enough that the plain solve
    # demonstrably leaks (≈ a quarter of its poles fall in the pole-free
    # sectors at this extent), small enough for a fast test.
    HALF = 12.0
    N    = 33
    xs   = range(-HALF, HALF; length = N)
    ys   = range(-HALF, HALF; length = N)

    prob = PadeTaylorProblem(pI, (u_tri, up_tri), (0.0, HALF); order = 30)

    # The PI tritronquée's poles live in the sector around the +real
    # axis; the sectors around ±144° are pole-free.  We measure poles
    # by |z| > 5 (clear of the calm origin) in two angular zones.
    in_populated(z) = abs(z) > 5 && abs(angle(z)) < π / 5            # +real wedge
    in_pole_free(z) = abs(z) > 5 && 2π / 3 < abs(angle(z)) ≤ π       # back sectors

    egs = edge_gated_pole_field_solve(prob, xs, ys; h = 0.5, grow_rings = 3)

    @testset "EG.1.1: the gated solve runs and grew a region" begin
        @test egs isa EdgeGatedSolution
        @test egs.iterations ≥ 2                       # seed + ≥1 growth pass
        @test count(egs.field_mask) > 0
        # The field grew well beyond the seed disc (seed_radius default
        # max(5, 4Δ) ≈ 5 here): a [-12,12]² tritronquée wedge is dozens
        # of cells.
        @test count(egs.field_mask) ≥ 40
    end

    @testset "EG.1.2: extracted pole field is sector-confined" begin
        poles = extract_poles(egs.pn_solution)
        npop  = count(in_populated, poles)
        nfree = count(in_pole_free, poles)
        # The +real-axis sector is genuinely populated ...
        @test npop ≥ 20
        # ... and the four pole-free sectors stay (near-)empty.  The
        # gated solve targets one thin smooth ring around the field, so
        # a handful of frontier stragglers are tolerated; the ratio is
        # the real invariant.
        @test nfree ≤ 5
        @test nfree < 0.1 * npop
    end

    @testset "EG.1.3: the plain ungated solve leaks into the pole-free sectors" begin
        # Same problem, same grid, but a plain full-grid path_network_solve
        # — the IVP is forced to target the smooth-sector cells, and
        # accumulated error there bleeds spurious poles into the four
        # nominally pole-free sectors.  This is the failure the gating
        # exists to prevent; the contrast is what makes the gating a
        # load-bearing component rather than decoration.
        grid_full = ComplexF64[x + im * y for y in ys for x in xs]
        pn_plain  = path_network_solve(prob, grid_full; h = 0.5)
        poles_plain = extract_poles(pn_plain)

        nfree_plain = count(in_pole_free, poles_plain)
        nfree_gated = count(in_pole_free, extract_poles(egs.pn_solution))

        @test nfree_plain ≥ 20                         # the leak is real
        @test nfree_gated < nfree_plain ÷ 4            # gating removes it
    end

    @testset "EG.2.1: fail-fast input validation" begin
        # grow_rings = 1 can never confirm a new cell.
        @test_throws ArgumentError edge_gated_pole_field_solve(
            prob, xs, ys; grow_rings = 1)
        # Non-uniform / mismatched steps break the isotropic stencil.
        @test_throws ArgumentError edge_gated_pole_field_solve(
            prob, range(-12.0, 12.0; length = 33), range(-12.0, 12.0; length = 41))
        # Degenerate axes.
        @test_throws ArgumentError edge_gated_pole_field_solve(
            prob, [0.0, 1.0], ys)
    end
end

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Each mutation below was applied to `src/EdgeGatedSolve.jl`, the suite
# confirmed RED, then the mutation reverted.
#
#   M1 — drop the edge-detector gate.  In the growth loop, replace
#        `new_cells = reachable .& targets .& .!field` with
#        `new_cells = targets .& .!field` (admit every solved cell, not
#        only edge-confirmed cells connected to the field).  The region
#        grows to the whole grid, the IVP traverses the smooth sectors,
#        and EG.1.2 (`nfree ≤ 5`) goes RED — the gated solve now leaks
#        just like the plain one.
#
#   M2 — disable region growing.  Make `_dilate` return its `mask`
#        argument unchanged (early `return out` before the dilation
#        loop).  `targets` never extends past the seed, no new cells
#        are ever confirmed: EG.1.1 (`count(field) ≥ 40`) and EG.1.2
#        go RED.
#
#   M3 — thicken the final frontier.  Change the final solve's
#        `_dilate(field, 1)` to `_dilate(field, 8)`.  The returned
#        `pn_solution` now carries eight rings of smooth frontier whose
#        visited Padés are corrupted; spurious poles reappear in the
#        pole-free sectors and EG.1.2 (`nfree ≤ 5`, `nfree < 0.1·npop`)
#        goes RED.
#
#   (A fourth candidate — shrinking the seed radius — does NOT bite:
#   region growing bootstraps fine from even a single seed cell, since
#   `grow_rings` reaches the first poles regardless.  That robustness
#   is a feature, so it is left untested rather than forced.)
# ----------------------------------------------------------------------
