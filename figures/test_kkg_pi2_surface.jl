# figures/test_kkg_pi2_surface.jl
#
# Acceptance test for the F3 whole-plane P_I^(2) tritronquée compute
# kernel (`figures/_kkg_pi2_surface_helpers.jl`, bead
# `padetaylor-0ln.33`).  Runs under the `figures/` Julia environment:
#
#     julia --project=figures figures/test_kkg_pi2_surface.jl
#
# It is a standalone `@testset` script — NOT part of the package `test/`
# suite — because it `include`s the Gridap-backed figure helpers, and
# Gridap is a `figures/`-project dependency only (ADR-0024 *Correction*).
#
# ## What is asserted (CLAUDE.md Rule 5 — every test pins a known value)
#
#   1. **Negative-real-axis pin.**  On `x < 0` the `V_0` tritronquée is
#      real, positive, with the KKG algebraic asymptotic
#      `u → +(6|x|)^{1/3}` (eq. (1.3); pillar C §4).  At `x ≈ -15` the
#      assembled `Re u` must match `+(6·15)^{1/3} ≈ 4.4814` and `Im u`
#      must vanish.  Tolerance `5e-3`: the `n_terms = 2` truncated
#      asymptotic series is itself only accurate to `~1.3e-4` at
#      `|x| = 15` (the V8b probe number), and the BVP-and-vote pipeline
#      adds the FEM voter's `~5e-4` and the ray-fan interpolation's
#      `~1e-3`; `5e-3` is a generous but honest envelope around all
#      three.
#
#   2. **Three-method agreement.**  Over the bulk of the pole-free
#      sector the max pairwise disagreement among Methods 1/2/3 (the
#      `spread` map) must be small.  Tolerance `1e-2`: the FEM voter is
#      `O(h²) ≈ 5e-5` on a smooth field but the inner-arc asymptotic
#      seed (v1 corner 1) injects an `O(10⁻²)` Dirichlet error at
#      `|x| ≈ 2` that the harmonic solve damps but does not erase, so
#      the worst-case spread sits near the inner arc.  We assert the
#      *median* spread is `< 3e-3` (the bulk concurs tightly) and the
#      *max* spread is `< 1e-2` (even the inner-arc corner stays
#      bounded).  The voted surface is verified to equal the
#      per-component median of the three voters.
#
#   3. **Sector coverage.**  ≥ 80% of the pole-free-sector grid cells
#      (inside the `|x| ≤ 20` disc) must carry a non-`NaN` value.
#
#   4. **Wedge.**  ≥ 1 pole in the wedge; every extracted pole confined
#      to `|arg x| ≲ 36°` (continuity with V8b's 21 poles in the wedge).
#
#   5. **Schwarz symmetry.**  `V_0(x̄) = conj(V_0(x))` (the ODE has real
#      coefficients and the tritronquée is real on `x < 0`), so in the
#      sector `Re u` is even and `Im u` is odd under `y ↦ -y`.
#
#   6. **Reproducibility.**  Two `kkg_pi2_surface()` calls are identical.
#
# This test does NOT assert "did not throw" anywhere — every block pins
# an invariant against a known-correct value (Rule 5).  Where a method
# is unreliable (the inner arc) the `spread` map exposes it and the test
# reports the number honestly rather than masking it.

using Test
using Statistics: median        # the INDEPENDENT median oracle for PI2S.2

include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))

# The compute kernel is expensive (a fan of BVPs + two Laplace solves +
# the wedge walk); run it ONCE and share the result across testsets.
const SURF = kkg_pi2_surface()
const SURF2 = kkg_pi2_surface()      # second run, for reproducibility

@testset "F3 — P_I^(2) tritronquée whole-plane surface" begin

    n = SURF_GRID_N

    @testset "PI2S.1 — negative-real-axis pin (KKG eq. 1.3)" begin
        i15 = argmin(abs.(SURF.xs .+ 15.0))
        j0  = argmin(abs.(SURF.ys))
        # The grid is constructed to land a node on x = -15, y = 0.
        @test isapprox(SURF.xs[i15], -15.0; atol = 1e-9)
        @test isapprox(SURF.ys[j0],    0.0; atol = 1e-9)

        expected = cbrt(6.0 * 15.0)            # +(6·15)^{1/3} ≈ 4.4814
        re15 = SURF.Re_u[i15, j0]
        im15 = SURF.Im_u[i15, j0]
        @test !isnan(re15)
        @test isapprox(re15, expected; atol = 5e-3)
        @test isapprox(im15,      0.0; atol = 5e-3)
        # The tritronquée is real and positive on x < 0.
        @test re15 > 0
    end

    @testset "PI2S.2 — three-method agreement / the majority vote" begin
        # Gather the sector spread map (the agreement diagnostic).
        sp = Float64[]
        for j in 1:n, i in 1:n
            s = SURF.spread[i, j]
            isnan(s) || push!(sp, s)
        end
        @test !isempty(sp)
        sort!(sp)
        med = sp[length(sp) ÷ 2]
        mx  = sp[end]
        # The bulk of the sector concurs tightly...
        @test med < 3e-3
        # ...and even the inner-arc v1 corner stays bounded.
        @test mx < 1e-2

        # The voted surface IS the per-component median of the three
        # voters.  Re-derive the median from `sector_method` with an
        # INDEPENDENT computation — `Statistics.median` on the finite
        # voters, NOT a call to the kernel's own `surf_vote` (calling
        # `surf_vote` here would be a tautology: a bug in `surf_vote`
        # would corrupt both sides equally).  This catches a kernel that
        # silently votes by the mean, the max, the first voter, etc.
        indep_vote(a, b, c) = begin
            v = filter(isfinite, (a, b, c))
            isempty(v) ? NaN : median(collect(v))
        end
        nchecked = 0; nfull = 0
        for ((i, j), m) in SURF.sector_method
            re1, im1, re2, im2, re3, im3 = m
            ref_re = indep_vote(re1, re2, re3)
            ref_im = indep_vote(im1, im2, im3)
            @test isapprox(SURF.Re_u[i, j], ref_re; atol = 1e-12, nans = true)
            @test isapprox(SURF.Im_u[i, j], ref_im; atol = 1e-12, nans = true)
            nchecked += 1
            # Count the grid points where all three voters are finite —
            # there the median is genuinely the *middle* of three
            # distinct numbers, so this block is a real discriminator
            # (the mean would differ from the median there).
            all(isfinite, (re1, re2, re3)) && (nfull += 1)
        end
        @test nchecked > 1000          # the vote ran on the whole sector
        @test nfull    > 1000          # ...and most points had 3 finite voters
    end

    @testset "PI2S.3 — pole-free-sector coverage ≥ 80%" begin
        total = 0; filled = 0
        for j in 1:n, i in 1:n
            z = ComplexF64(SURF.xs[i], SURF.ys[j])
            (abs(z) <= SURF_XY_LIM && surf_in_sector(z)) || continue
            total += 1
            isnan(SURF.Re_u[i, j]) || (filled += 1)
        end
        @test total > 0
        @test filled / total ≥ 0.80
    end

    @testset "PI2S.4 — wedge pole field" begin
        @test length(SURF.poles) ≥ 1
        # Every pole confined to the |arg x| ≲ 36° wedge (a small
        # numerical margin past the 36° Stokes line is allowed for the
        # clustered root positions).
        for p in SURF.poles
            @test abs(rad2deg(angle(p))) ≤ 36.0 + 2.0
        end
        # Continuity with V8b: the same recipe found 21 poles.
        @test length(SURF.poles) ≥ 10
    end

    @testset "PI2S.5 — Schwarz symmetry across the real axis" begin
        # V_0(x̄) = conj(V_0(x)): in the sector Re u is even, Im u odd
        # under y ↦ -y.  Check at a handful of conjugate sector pairs.
        npairs = 0
        for i in 1:n
            x = SURF.xs[i]
            x < -3.0 || continue                # stay on the negative side
            for jt in (8.0, 12.0)
                jp = argmin(abs.(SURF.ys .- jt))
                jm = argmin(abs.(SURF.ys .+ jt))
                rp = SURF.Re_u[i, jp]; rm = SURF.Re_u[i, jm]
                ip = SURF.Im_u[i, jp]; im_ = SURF.Im_u[i, jm]
                (isnan(rp) || isnan(rm)) && continue
                @test isapprox(rp,  rm; atol = 1e-6)
                @test isapprox(ip, -im_; atol = 1e-6)
                npairs += 1
            end
        end
        @test npairs > 5
    end

    @testset "PI2S.6 — reproducibility" begin
        # Two independent runs produce byte-identical matrices (NaNs in
        # the same cells; finite cells equal).
        @test size(SURF.Re_u) == size(SURF2.Re_u)
        diffs = 0
        for idx in eachindex(SURF.Re_u)
            a, b = SURF.Re_u[idx], SURF2.Re_u[idx]
            c, d = SURF.Im_u[idx], SURF2.Im_u[idx]
            if isnan(a) != isnan(b) || isnan(c) != isnan(d)
                diffs += 1
            elseif !isnan(a) && (a != b || c != d)
                diffs += 1
            end
        end
        @test diffs == 0
        @test length(SURF.poles) == length(SURF2.poles)
    end
end
