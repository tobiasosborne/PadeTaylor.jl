# test/phase9_tritronquee_test.jl -- Phase 9 / bead `padetaylor-kvi` tests.
#
# Tier-C qualitative reproduction of PI tritronquée pole-field structure
# per DESIGN.md §4 Phase 9 + FW 2011 §2 lines 52-54 (5-fold sectorial
# structure) + §4.1 lines 214-229 (eq. 4.1 tritronquée ICs).
#
# Phase 9 is "no new code" per DESIGN.md — the deliverable is a
# composition test: PathNetwork (Phase 10, IVP path-tree over 2D grid)
# + EdgeDetector (Phase 12.5, FW eq. 3.3 5-point Laplacian) on the PI
# tritronquée problem.  Asserts the four qualitative invariants that
# define a Painlevé I tritronquée pole field:
#
#   PT.1.1 — composition runs; deterministic; pole field non-empty.
#   PT.1.2 — sectorial structure: 8 of 12 30°-bins are empty (the
#            "pole-free 4 of 5 sectors" property of Boutroux's
#            tritronquée).
#   PT.1.3 — y → −y conjugate symmetry: mask is symmetric under
#            reflection across the real axis (real ICs ⇒ u(z̄) = ū(z),
#            so poles come in conjugate pairs).
#   PT.1.4 — first real-axis pole magnitude check: |u(z ≈ 2.33)| > 100,
#            confirming the leading pole on the positive real axis is
#            captured.  Pole position consistent with FW Table 5.1
#            (the next-pole-after-2.07 is at z ≈ 2.07 for the
#            tritronquée; we hit a near-pole grid cell at z=2.33).
#   PT.2.1 — discriminative test: perturbed IC (u(0)=0, u'(0)=0) gives
#            a DIFFERENT sector pattern than tritronquée.  Ensures the
#            sector test is sensitive to the actual ODE state, not a
#            generic "any solution gives this pattern" artefact.
#
# Mutation-proof procedure documented at bottom (perturbations to the
# ODE definition; bites the sector-asymmetry assertion).

using Test
using PadeTaylor

@testset "Phase 9 — PI tritronquée pole-field qualitative" begin

    # ODE: u'' = 6u² + z (Painlevé I).
    f_PI(z, u, up) = 6 * u^2 + z

    # FW 2011 eq. (4.1) Maple-computed tritronquée ICs.  Real ICs ⇒
    # u(x) is real on the real axis and u(z̄) = ū(z) globally.
    u_tri  = -0.1875543083404949
    up_tri =  0.3049055602612289

    # 25×25 grid over [-4, 4]² — captures the first pole on the
    # positive real axis (z ≈ 2.33 on this lattice) plus enough
    # surrounding cells to see the sectorial pattern.  Larger windows
    # would catch more sectors but cost suite time; 25×25 at h=0.5
    # path-network runs in ≈ 2 s.
    N = 25
    xs = range(-4.0, 4.0; length = N)
    ys = range(-4.0, 4.0; length = N)
    # Convention: grid[i, j] = xs[i] + im·ys[j].  Vec'd column-major
    # so reshape recovers the same orientation.
    grid_mat = ComplexF64[xs[i] + im * ys[j] for i in 1:N, j in 1:N]
    grid_vec = vec(grid_mat)
    h_grid   = step(xs)

    # `zspan[2]` is unused for path_network_solve (the grid drives Stage
    # 1); set it to a far-corner placeholder to satisfy
    # PadeTaylorProblem's non-degeneracy guard.
    zspan = (0.0 + 0.0im, ComplexF64(4.0 * sqrt(2)))
    prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

    sol      = path_network_solve(prob, grid_vec; h = 0.5)
    u_grid   = reshape(sol.grid_u, (N, N))   # u_grid[i, j] = u(xs[i] + im·ys[j])
    mask     = pole_field_mask(u_grid, h_grid)

    # Helper: histogram of flagged cells over 12 × 30° angular bins,
    # bin k spans [-180° + (k-1)·30°, -180° + k·30°).
    function sector_counts(mask)
        nbins  = 12
        counts = zeros(Int, nbins)
        nrow, ncol = size(mask)
        for j in 2:(ncol - 1), i in 2:(nrow - 1)
            if mask[i, j]
                x, y = xs[i], ys[j]
                θ    = atan(y, x)
                k    = clamp(Int(floor((θ + π) / (2π / nbins))) + 1, 1, nbins)
                counts[k] += 1
            end
        end
        return counts
    end

    @testset "PT.1.1: composition runs; pole field non-empty + deterministic" begin
        # Sanity check: the composition completes without error, and the
        # output is non-trivial (poles exist in the window).  Empirically
        # 67 cells flagged out of (N-2)² = 529 interior cells; the bounds
        # are generous to absorb minor wedge-walk path variations.
        @test count(mask) ≥ 20            # at least a small pole field
        @test count(mask) < (N - 2)^2 ÷ 2 # not majority-flagged
        @test all(isfinite, sol.grid_u)   # no NaN sentinels (full coverage)

        # Determinism: a second call with the same kwargs returns the
        # same mask exactly.  (path_network_solve threads a MersenneTwister
        # with rng_seed=0 by default.)
        sol2  = path_network_solve(prob, grid_vec; h = 0.5)
        mask2 = pole_field_mask(reshape(sol2.grid_u, (N, N)), h_grid)
        @test mask == mask2
    end

    @testset "PT.1.2: 4-of-5 pole-free sectors (Boutroux's tritronquée)" begin
        # PI's 5-fold sectorial structure (FW2011...md:52) gives 5
        # asymptotic sectors of width 2π/5 = 72° each.  Boutroux's
        # tritronquée has pole-free behaviour in 4 of these — i.e., a
        # contiguous ≈ 4·72° = 288° pole-free arc.  In our 12 × 30°
        # binning, that means at LEAST 8 of 12 bins are empty
        # (288°/30° = 9.6, but binning boundaries cost us a bin or two).
        counts = sector_counts(mask)
        n_empty_bins = count(==(0), counts)
        @test n_empty_bins ≥ 8

        # The pole-bearing sector is centred on the positive real axis
        # (angle 0).  The bins [-30°, 0°) (bin 6) and [0°, 30°) (bin 7)
        # together must contain a strong majority of the flagged cells —
        # ≥ 75 % is a conservative bound for the asymmetry.
        center_two = counts[6] + counts[7]
        @test center_two ≥ 0.75 * sum(counts)

        # The opposite (negative real axis) bin pair must be empty:
        # bin 1 [-180°, -150°) and bin 12 [150°, 180°).
        @test counts[1] == 0
        @test counts[12] == 0
    end

    @testset "PT.1.3: y → −y conjugate symmetry (real ICs ⇒ u(z̄) = ū(z))" begin
        # Real ICs make u(x) real on the real axis; analytic continuation
        # off the axis yields u(z̄) = ū(z) by the reflection principle.
        # Poles come in conjugate pairs, so the discrete mask is
        # invariant under j → (N + 1 − j) (the y → −y reflection on
        # our centred [-4, 4] grid).  Tolerate up to 4 mismatching
        # cells from path-network walk asymmetry (slightly different
        # path trees on top vs bottom can produce different round-off
        # residuals in EdgeDetector at the edge of the flagged region).
        mask_reflected = falses(N, N)
        for i in 1:N, j in 1:N
            mask_reflected[i, j] = mask[i, N + 1 - j]
        end
        mismatches = count(mask .!= mask_reflected)
        @test mismatches ≤ 4
    end

    @testset "PT.1.4: leading positive-real-axis pole magnitude" begin
        # The flagged cell on the real axis closest to (but past) the
        # tritronquée's first real-axis pole at z ≈ 2.07 lands at
        # x ≈ 2.33 on our 25-pt grid (xs[18] = -4 + 17/3 = 2.333…).
        # |u| at this cell exceeds 100 (empirically ≈ 387) — confirms
        # the pole is captured, not interpolated through.
        i_real = findmin(abs.(ys))[2]              # j-col where y ≈ 0
        x_near = 2.33
        j_x    = findmin(abs.(xs .- x_near))[2]   # i-row where x ≈ 2.33
        @test abs(u_grid[j_x, i_real]) > 100.0
        # ... and the corresponding cell IS flagged in the mask.
        @test mask[j_x, i_real] == true
    end

    @testset "PT.2.1: discriminative — null-IC sector pattern differs" begin
        # Reset the IC to (u(0), u'(0)) = (0, 0).  This is NOT a
        # tronquée: per FW Fig 4.3 + 4.4, u(0)=0 gives a near-tronquée
        # with pole fields on BOTH sides of the real axis, very
        # different from tritronquée's single-sector concentration.
        # The 4-of-5 sector-free property must NOT hold for this IC.
        prob_null = PadeTaylorProblem(f_PI, (0.0, 0.0), zspan; order = 30)
        sol_null  = path_network_solve(prob_null, grid_vec; h = 0.5)
        u_null    = reshape(sol_null.grid_u, (N, N))
        mask_null = pole_field_mask(u_null, h_grid)

        counts_null = let nbins = 12
            cs = zeros(Int, nbins)
            for j in 2:(N - 1), i in 2:(N - 1)
                if mask_null[i, j]
                    x, y = xs[i], ys[j]
                    θ    = atan(y, x)
                    k    = clamp(Int(floor((θ + π) / (2π / nbins))) + 1, 1, nbins)
                    cs[k] += 1
                end
            end
            cs
        end
        # Pole-bearing centre-two bins should NOT dominate as
        # overwhelmingly as in the tritronquée case (>75 %).  We assert
        # the WEAKER property: the centre-two share is <70 % for the
        # null IC, but >75 % for tritronquée — i.e., the two patterns
        # are quantitatively distinguishable.
        share_null = (counts_null[6] + counts_null[7]) / max(1, sum(counts_null))
        @test share_null < 0.70
    end

end

# Mutation-proof procedure (verified 2026-05-13).
#
# Mutation C — flip sign of the `z` driver in the ODE:
#     f_PI(z, u, up) = 6 * u^2 - z   (was + z)
# This converts PI to a related but distinct equation.  The pole
# pattern's sectorial structure rotates and the leading pole moves to
# negative-real-axis territory.  Bites:
#   - PT.1.2: the centre-two bins (positive real) no longer dominate;
#     `center_two ≥ 0.75 * sum(counts)` fails.
#   - PT.1.4: |u(2.33, 0)| is no longer >100 (the pole has moved away
#     from this lattice cell).
#
# Mutation D — perturb the tritronquée IC by 1e-2 (1000x the FW eq. 4.1
# uncertainty):
#     u_tri = -0.1875543083404949 - 0.01    # was without the -0.01
# This is far enough from tritronquée that the solution acquires poles
# in the other 4 sectors (per FW 2011 Fig 3.1's "near-tritronquée" with
# u(0)=-0.1875: pole fields appear in all 5 sectors, just pushed
# further out in 4 of them).  Within [-4, 4]² window, the 4-of-5
# pole-free property still holds approximately, so this mutation does
# NOT bite the sector test — confirming that the test is robust to
# small IC perturbations within the tritronquée's stability radius
# (an intentional feature, not a bug; the test asserts QUALITATIVE
# tritronquée structure, not bitwise-exact ICs).
#
# Both mutations restored before commit per CLAUDE.md Rule 4.
