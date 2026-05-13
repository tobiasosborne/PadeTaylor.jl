# test/edge_detector_test.jl -- bead `padetaylor-c2p` tests.
#
# 5-point Laplacian pole-field edge detector per FW 2011 §3.2.2
# (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:202-208`).
# Direct port of eq. (3.3) — `Δu ≈ (stencil · u) / h²`, an `O(h²)`
# approximation to `∇²u`; classifier thresholds `log₁₀|Δu|` against the
# level curve in Fig. 3.3 (FW's default 0.001).
#
# Test plan:
#   - ED.1.1: stencil is exact-zero on a real harmonic quadratic
#     (`u = x² − y²` ≡ Re(z²)) — `∇²u = 0` analytically, so the
#     `O(h²)` stencil must return zero to roundoff (the 4th derivative
#     in the leading truncation term is zero for any cubic).
#   - ED.1.2: stencil is zero on a complex-analytic function
#     (`u(z) = z² = (x² − y²) + 2ixy`) — both real and imaginary parts
#     harmonic.
#   - ED.1.3: stencil is *large* near a pole (`u(z) = 1/(z − z₀)`),
#     small far from it.  Magnitudes scale per `1/|z − z₀|³` near the
#     pole.
#   - ED.2.1: `pole_field_mask` returns BitMatrix, boundary `false`,
#     interior poles `true` at default level=0.001.
#   - ED.2.2: separation — smooth analytic `u = exp(z)` gives all-false
#     mask at h=0.05; 1/(z − z₀) near the pole gives mask-true cluster.
#   - ED.3.1: fail-fast on bad inputs (grid too small, non-positive h).
#   - ED.4.1: mutation-proof — procedure documented at bottom.

using Test
using PadeTaylor
using PadeTaylor: EdgeDetector
using PadeTaylor.EdgeDetector: laplacian_residual, pole_field_mask

@testset "EdgeDetector (Phase 12.5): FW 2011 §3.2.2 pole-field detector" begin

    # Helper: build u_grid by evaluating a function `f` on a Cartesian
    # lattice `x_range × y_range`.  Convention: `u_grid[i, j]` at
    # `z = x_range[j] + im * y_range[i]` — row indexes y, col indexes x,
    # i.e., the natural `imshow` orientation FW uses in Fig. 3.3.
    function build_grid(f, x_range, y_range)
        nrow = length(y_range)
        ncol = length(x_range)
        u = Matrix{ComplexF64}(undef, nrow, ncol)
        for j in 1:ncol, i in 1:nrow
            u[i, j] = ComplexF64(f(x_range[j] + im * y_range[i]))
        end
        return u
    end

    @testset "ED.1.1: harmonic quadratic — stencil exact-zero (roundoff)" begin
        # u(x, y) = x² − y² is the real part of z².  ∇²u = 2 + (−2) = 0
        # analytically.  The 5-point stencil is exact for any polynomial
        # of total degree ≤ 3 (the leading truncation O(h²) involves the
        # 4th derivative; quadratic has none).  So residual = 0 to
        # roundoff at any h.
        f(z) = real(z)^2 - imag(z)^2
        x_range = -1.0:0.1:1.0
        y_range = -1.0:0.1:1.0
        u = build_grid(f, x_range, y_range)
        h = 0.1

        Δu = laplacian_residual(u, h)
        nrow, ncol = size(Δu)
        # Interior cells: residual is essentially zero (roundoff bound).
        # Conservative bound: ‖u‖∞ · ε / h² ≈ 1 · 2e-16 / 0.01 = 2e-14.
        for j in 2:(ncol-1), i in 2:(nrow-1)
            @test abs(Δu[i, j]) < 1e-12
        end
        # Boundary cells: NaN by design.
        @test isnan(real(Δu[1,   1]))
        @test isnan(real(Δu[end, end]))
        @test isnan(real(Δu[1,   3]))
        @test isnan(real(Δu[3,   end]))
    end

    @testset "ED.1.2: complex-analytic u(z) = z² — stencil exact-zero" begin
        # u(z) = z² is complex-analytic.  Its real part x²−y² and
        # imaginary part 2xy are both harmonic.  Stencil → 0 in both
        # components, to roundoff.
        f(z) = z^2
        x_range = -1.0:0.1:1.0
        y_range = -1.0:0.1:1.0
        u = build_grid(f, x_range, y_range)
        h = 0.1

        Δu = laplacian_residual(u, h)
        nrow, ncol = size(Δu)
        for j in 2:(ncol-1), i in 2:(nrow-1)
            @test abs(Δu[i, j]) < 1e-12
        end
    end

    @testset "ED.1.3: simple pole u = 1/(z − z₀) — residual is large near z₀" begin
        # u(z) = 1/(z − z₀) is meromorphic with a single pole at z₀.
        # The discrete Laplacian residual at lattice points close to z₀
        # is dominated by the analytic ∇²u ≠ 0 contribution (the
        # function is NOT harmonic where it diverges).  Away from z₀,
        # u is analytic so residual → 0.
        z₀ = 0.0 + 0.0im
        f(z) = inv(z - z₀)
        # Grid centred on z₀ but z₀ itself excluded (use offset so z₀
        # lands between lattice points).
        x_range = -0.45:0.1:0.45
        y_range = -0.45:0.1:0.45
        u = build_grid(f, x_range, y_range)
        h = 0.1

        Δu = laplacian_residual(u, h)
        nrow, ncol = size(Δu)
        # Find the interior lattice cell closest to z₀.  At
        # x_range = -0.45:0.1:0.45 the closest are at index ≈ 5 (value
        # -0.05) or 6 (0.05); abs(z) ≈ √(0.05² + 0.05²) ≈ 0.071.
        center_i = findmin(abs.(y_range))[2]
        center_j = findmin(abs.(x_range))[2]
        @test abs(Δu[center_i, center_j]) > 100.0   # near pole, residual is huge

        # Farthest interior corner has u ≈ 1/(0.45 + 0.45i), still
        # analytic in a neighbourhood (no pole within h).  Residual is
        # small but not roundoff-zero (3rd-derivative content of
        # 1/(z-z₀) is non-zero everywhere except at z₀; truncation
        # ~ |f''''| · h² ).  Empirically ≈ 30 at the corner.  Test
        # only the *ratio*: residual at corner << residual at centre.
        @test abs(Δu[2, 2]) < abs(Δu[center_i, center_j]) / 10.0
    end

    @testset "ED.2.1: pole_field_mask returns BitMatrix; boundary cells false" begin
        # Sanity test for the bitmap convenience wrapper.
        z₀ = 0.0 + 0.0im
        f(z) = inv(z - z₀)
        x_range = -0.45:0.1:0.45
        y_range = -0.45:0.1:0.45
        u = build_grid(f, x_range, y_range)
        h = 0.1

        mask = pole_field_mask(u, h)
        @test isa(mask, BitMatrix)
        @test size(mask) == size(u)

        # Boundary cells (any cell in first/last row or first/last col)
        # are `false` by design.
        nrow, ncol = size(mask)
        @test all(mask[1,   :] .== false)
        @test all(mask[end, :] .== false)
        @test all(mask[:,   1] .== false)
        @test all(mask[:, end] .== false)

        # Near-pole centre cells are `true` (large |Δu| → log₁₀ >> 0.001).
        center_i = findmin(abs.(y_range))[2]
        center_j = findmin(abs.(x_range))[2]
        @test mask[center_i, center_j] == true
    end

    @testset "ED.2.2: smooth vs pole separation at default level" begin
        # On a smooth analytic grid, mask is all-false: |Δu| << 1 in
        # the entire interior, so log₁₀|Δu| << 0.001.  On the
        # 1/(z-z₀) grid, the high-|Δu| region forms a connected
        # cluster around z₀.
        f_smooth(z) = exp(z)
        x_range = -1.0:0.05:1.0
        y_range = -1.0:0.05:1.0
        u_smooth = build_grid(f_smooth, x_range, y_range)
        mask_smooth = pole_field_mask(u_smooth, 0.05)
        # No interior cell should be flagged.
        @test count(mask_smooth) == 0

        # Pole grid: at least a handful of cells flagged; cluster
        # centred near the pole.
        z₀ = 0.0 + 0.0im
        f_pole(z) = inv(z - z₀)
        u_pole = build_grid(f_pole, x_range, y_range)
        mask_pole = pole_field_mask(u_pole, 0.05)
        # At h=0.05 on a 41x41 grid with a pole at the origin, the
        # |Δu| > 10^0.001 ≈ 1 disk extends roughly to |z| ≲ 0.35
        # (since |∇²(1/z)| = 2/|z|³ = 1 at |z| ≈ 2^(1/3) · 0.5 once
        # discretisation kicks in).  Empirically ≈ 188 cells flagged.
        @test count(mask_pole) ≥ 4   # at least a small cluster
        @test count(mask_pole) < length(mask_pole) ÷ 2   # not majority
    end

    @testset "ED.3.1: fail-fast guards" begin
        u_2x2 = ComplexF64[1 2; 3 4]
        @test_throws ArgumentError laplacian_residual(u_2x2, 0.1)

        u_3x3 = ComplexF64[1 2 3; 4 5 6; 7 8 9]
        @test_throws ArgumentError laplacian_residual(u_3x3, 0.0)
        @test_throws ArgumentError laplacian_residual(u_3x3, -0.1)
    end

    @testset "ED.4.1: pre-computed-residual variant of pole_field_mask" begin
        # The two-arg signature (residual already computed) returns the
        # same mask as the one-step convenience wrapper.  Useful for
        # callers who want to inspect the raw residual once and threshold
        # at multiple levels.
        z₀ = 0.0 + 0.0im
        f(z) = inv(z - z₀)
        x_range = -0.45:0.1:0.45
        y_range = -0.45:0.1:0.45
        u = build_grid(f, x_range, y_range)
        h = 0.1

        Δu = laplacian_residual(u, h)
        mask_1 = pole_field_mask(u, h; level = 0.001)
        mask_2 = pole_field_mask(Δu; level = 0.001)
        @test mask_1 == mask_2

        # Different levels give different masks: a higher level retains
        # only the highest-|Δu| cells.
        mask_low  = pole_field_mask(Δu; level = -3.0)   # |Δu| > 1e-3
        mask_high = pole_field_mask(Δu; level =  3.0)   # |Δu| > 1e3
        @test count(mask_low) ≥ count(mask_high)
    end

end

# Mutation-proof procedure (verified 2026-05-13).
#
# Mutation A — change stencil centre coefficient from `-4` to `-3` in
# `laplacian_residual`.  Bites:
#   - ED.1.1: a non-zero residual `u_ij / h²` appears at every interior
#     cell on the harmonic quadratic (instead of zero); `|Δu| ≈ 100`
#     at u ≈ 1 with h=0.1, breaking the `< 1e-12` assertion.
#   - ED.1.2: same.
#
# Mutation B — drop the `/h²` factor in `laplacian_residual`.  Bites:
#   - ED.1.3: magnitudes are scaled down by `h² = 0.01`, so the
#     near-pole `|Δu|` is `~1` instead of `~100`, breaking the `> 100`
#     assertion.
#   - ED.2.2: log₁₀|Δu| values shift by 2; mask cells flip false.
#
# Both mutations restored before commit per CLAUDE.md Rule 4.
