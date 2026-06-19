# test/teukolsky_qnm_test.jl — Stage 4c of the P1 Heun epic (bead padetaylor-4lh).
#
# Pins the Schwarzschild quasinormal-mode FREQUENCY recovered by Leaver's
# radial continued fraction against the cross-verified literature values in
# `docs/teukolsky_heun_mapping.md` §3.2.1 (lines 252-255).  The continued
# fraction is the numerically-stable form of the confluent-Heun
# minimal-solution / connection condition (`docs/teukolsky_heun_mapping.md:256-260`).
#
# Per CLAUDE.md Rule 5 every assertion checks a KNOWN-CORRECT value:
#   - QNM.1 pins both recovered Mω (s=0 and s=2) to the Leaver targets, rtol 1e-4.
#   - QNM.2 pins `leaver_recurrence` to the closed forms at a sample (n,ω) for
#     BOTH spins, so the β_aux = s²−1 branch (γ_n = (n−4iω)² − (β_aux+1)) is
#     pinned independently of the root-find.
#
# We `include` the example's functions directly (the example is the
# single source of truth for the recurrence and CF); the example does
# `using PadeTaylor, Printf`, both on the Pkg.test dependency path.

using Test

# Bring the example's leaver_recurrence / qnm_continued_fraction / find_qnm
# into scope.  The `if abspath(PROGRAM_FILE)==@__FILE__` guard means the
# include does NOT run the script's root-finds.
include(joinpath(@__DIR__, "..", "examples", "teukolsky_qnm.jl"))

@testset "QNM.1 — Schwarzschild QNM frequencies (Leaver continued fraction)" begin
    # PRIMARY oracle: scalar s=0, ℓ=2, n=0  (doc §3.2.1, line 254).
    ref0 = 0.483642 - 0.096766im
    ω0 = find_qnm(0.4836 - 0.0968im; s = 0, ℓ = 2, N = 100)
    @test abs(ω0 - ref0) / abs(ref0) ≤ 1e-4

    # CROSS-CHECK: gravitational s=2, ℓ=2, n=0  (doc §3.2.1, line 255).
    # Exercises the general-s β_aux branch (β_aux = +3).
    ref2 = 0.373672 - 0.088962im
    ω2 = find_qnm(0.3737 - 0.0890im; s = 2, ℓ = 2, N = 100)
    @test abs(ω2 - ref2) / abs(ref2) ≤ 1e-4
end

@testset "QNM.2 — leaver_recurrence closed forms (β_aux branch pinned)" begin
    # Sample point: n=3, ω = 0.5 - 0.1i.  We assert the coefficients against
    # the closed forms in doc §3.2.1 (lines 234-244) recomputed independently
    # here, so a sign/typo in the impl is caught for each spin.
    n = 3
    ω = 0.5 - 0.1im

    # --- s=0: β_aux = -1, so γ_n = (n - 4iω)². ---
    α0, β0, γ0 = leaver_recurrence(n, ω; s = 0, ℓ = 2)
    βaux0 = 0^2 - 1
    @test α0 ≈ (n + 1) * (n + 1 - 4im * ω) atol = 1e-12
    @test β0 ≈ -2n^2 + (-2 + 16im * ω) * n + 32 * ω^2 + 8im * ω - 2 * (2 + 1) + βaux0 atol = 1e-12
    @test γ0 ≈ (n - 4im * ω)^2 atol = 1e-12            # β_aux+1 = 0
    @test γ0 ≈ (n - 4im * ω)^2 - (βaux0 + 1) atol = 1e-12

    # --- s=2: β_aux = +3, so γ_n = (n - 4iω)² - 4. ---
    α2, β2, γ2 = leaver_recurrence(n, ω; s = 2, ℓ = 2)
    βaux2 = 2^2 - 1
    @test α2 ≈ (n + 1) * (n + 1 - 4im * ω) atol = 1e-12     # α independent of s
    @test β2 ≈ -2n^2 + (-2 + 16im * ω) * n + 32 * ω^2 + 8im * ω - 2 * (2 + 1) + βaux2 atol = 1e-12
    @test γ2 ≈ (n - 4im * ω)^2 - 4 atol = 1e-12            # β_aux+1 = 4
    @test γ2 ≈ (n - 4im * ω)^2 - (βaux2 + 1) atol = 1e-12

    # The two spins MUST differ in β and γ (guards against β_aux being dropped).
    @test β0 != β2
    @test γ0 != γ2
end
