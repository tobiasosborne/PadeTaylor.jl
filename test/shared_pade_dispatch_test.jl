# shared_pade_dispatch_test.jl — ADR-0028 dual-construction dispatch (cell B +
# relative-ODE-defect selector), tested in ISOLATION before the stepper is wired
# (Phase 1).  Standalone @testset: `julia --project=. test/shared_pade_dispatch_test.jl`.
#
# The dispatch is the production realisation of the probe
# external/probes/adr0028-dual-construction-audition/ (cells.jl `cell_B_grow`,
# `relative_defect`; validated 18/18 in pole_crossing.jl).  These tests pin:
#   SPD.1  cell B = wide-square Mano–Tsuda (m_eff=⌈m/d⌉·d): recovers a known rational;
#          d=1 reduces to robust_pade(m-1,m;:svd) (the off-diagonal scalar oracle).
#   SPD.2  shared_pade_select reproduces the probe's per-step picks, incl. the lone
#          genuine A-win (Calogero–Moser at BigFloat-256) and B elsewhere.
#   SPD.3  the degenerate guards fall back to cell A (never throw, never lie).
#   SPD.4  determinism: an eps-perturbation of the jets does not flip the cell.
#   SPD.5  complex step h is handled (conjugate-pair pole straddling the real t-axis).

using Test
using PadeTaylor
using PadeTaylor: shared_denominator_pade, robust_pade, vector_taylor_coefficients,
                  noumi_yamada_rhs
using PadeTaylor.SharedPadeCellB: build_square_cell
using PadeTaylor.SharedPadeDispatch: shared_pade_select
using LinearAlgebra: eigvals

# --- helpers ----------------------------------------------------------------
polyval(c, z) = isempty(c) ? zero(z) : sum(c[k] * z^(k - 1) for k = 1:length(c))

function ratio_jet(P, Q, N)                # formal P/Q Taylor jet to order N
    c = zeros(ComplexF64, N + 1)
    for k = 0:N
        s = (k + 1) ≤ length(P) ? ComplexF64(P[k + 1]) : 0.0 + 0im
        for j = 1:k
            qj = (j + 1) ≤ length(Q) ? ComplexF64(Q[j + 1]) : 0.0 + 0im
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / ComplexF64(Q[1])
    end
    c
end
qdiff(x, y) = (L = max(length(x), length(y));
               maximum(abs.(vcat(x, zeros(eltype(x), L - length(x))) .-
                            vcat(y, zeros(eltype(y), L - length(y))))))
function rootset(Q)
    n = length(Q) - 1; n == 0 && return ComplexF64[]
    cc = ComplexF64.(Q) ./ Q[end]
    C = zeros(ComplexF64, n, n)
    for i = 1:n-1; C[i+1, i] = 1; end
    for i = 1:n; C[i, n] = -cc[i]; end
    sort(real.(eigvals(C)))
end
# build rescaled jets for an RHS, exactly as VectorStepper does (c̃_k = h^k c_k)
function rjets(f, z0, y0, h, order)
    raw = vector_taylor_coefficients(f, z0, y0, order)
    [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:length(y0)]
end
pick(f, z0, y0, h, m) =
    shared_pade_select(rjets(f, z0, y0, h, 2m + 2), m, f, z0, h; diagnostics = true)[3].chosen_cell

@testset "ADR-0028 SharedPade dispatch" begin

    @testset "SPD.1 cell B = wide-square Mano–Tsuda" begin
        # rational with shared poles {2,3}: (1-z/2)(1-z/3) = 1 - 5/6 z + 1/6 z²
        Qt = [1.0, -5 / 6, 1 / 6]; P1 = [1.0, 2.0]; P2 = [3.0, -1.0]; m = 2
        j1 = ratio_jet(P1, Qt, 2m + 2); j2 = ratio_jet(P2, Qt, 2m + 2)
        cb = build_square_cell([j1, j2], m)
        @test cb !== nothing
        numsB, denB = cb
        @test rootset(denB) ≈ [2.0, 3.0] atol = 1e-8        # SPD.1a recovers the poles
        @test abs(polyval(numsB[1], 0.5 + 0im) / polyval(denB, 0.5 + 0im) -
                  polyval(P1, 0.5) / polyval(Qt, 0.5)) < 1e-10    # SPD.1b eval matches

        # d=1: cell B is the off-diagonal (m-1,m) = robust_pade(jet, m-1, m; :svd)
        je = ComplexF64[inv(factorial(big(k))) for k = 0:2m+2]   # exp jet
        md = 6
        cb1 = build_square_cell([je], md)
        @test cb1 !== nothing
        rpB = robust_pade(je, md - 1, md; method = :svd)
        @test qdiff(cb1[2], rpB.b ./ rpB.b[1]) / maximum(abs.(rpB.b ./ rpB.b[1])) < 1e-8  # SPD.1c
        # cell B (degree m_eff=md, numerator m_eff-1) differs from cell A (m,m)
        ca1 = shared_denominator_pade([je], md)
        @test length(cb1[1][1]) < length(ca1[1][1])             # SPD.1d μ_B < μ_A at d=1
    end

    @testset "SPD.2 selector reproduces probe picks" begin
        # entire systems → square cell B wins
        harm = (z, y) -> [y[2], -y[1]]
        @test pick(harm, 0.0 + 0im, ComplexF64[1, 0], 0.7, 15) == :square      # SPD.2a
        expp = (z, y) -> [y[1], y[1] + y[2]]
        @test pick(expp, 0.0 + 0im, ComplexF64[1, 0], 0.2, 15) == :square      # SPD.2b
        # Calogero–Moser (d=4 meromorphic): F64 regular step → B
        cm = (z, y) -> (dd = y[1] - y[2]; a = 4 / dd^3; [y[3], y[4], a, -a])
        @test pick(cm, 0.0 + 0im, ComplexF64[1, -1, 0, 0], 0.1, 15) == :square # SPD.2c
        # tan-companion pole-CROSSING (double pole at π/2) → B
        tanc = (z, y) -> [y[2], 2y[1] + 2y[1]^3]
        @test pick(tanc, 0.0 + 0im, ComplexF64[0, 1], 1.8, 12) == :square      # SPD.2d
        # d=1 short-circuit → cell A
        @test pick((z, y) -> [y[1]], 0.0 + 0im, ComplexF64[1], 0.3, 8) == :diagonal_d1  # SPD.2e

        # THE GENUINE A-WIN: Calogero–Moser at BigFloat-256 (FW Table 5.1 regime)
        setprecision(BigFloat, 256) do
            CBF = Complex{BigFloat}
            cmb = (z, y) -> (dd = y[1] - y[2]; a = 4 / dd^3; [y[3], y[4], a, -a])
            @test pick(cmb, CBF(0), CBF[1, -1, 0, 0], BigFloat("0.1"), 15) == :diagonal  # SPD.2f
            # harmonic at BF still picks B
            harmb = (z, y) -> [y[2], -y[1]]
            @test pick(harmb, CBF(0), CBF[1, 0], BigFloat("0.7"), 15) == :square          # SPD.2g
        end
    end

    @testset "SPD.3 degenerate guards fall back to cell A" begin
        je = ComplexF64[inv(factorial(big(k))) for k = 0:30]
        # d ≥ m guard: build_square_cell returns nothing
        @test build_square_cell([je, je, je], 2) === nothing          # SPD.3a d≥m
        # jet too short for m_eff+1
        short = ComplexF64[1, 0.5, 0.25]
        @test build_square_cell([short, short], 8) === nothing        # SPD.3b short jet
        # selector with a tripped guard returns cell A exactly.  m=7 (odd, d=2) ⇒
        # m_eff=8 needs 9 coeffs; a length-8 jet (order 7) feeds cell A (needs 8) but
        # trips cell B's short-jet guard ⇒ fall back to A.
        harm = (z, y) -> [y[2], -y[1]]
        rj = rjets(harm, 0.0 + 0im, ComplexF64[1, 0], 0.5, 7)         # length 8
        @test build_square_cell(rj, 7) === nothing                    # SPD.3c cell B guard trips
        nums, den = shared_pade_select(rj, 7, harm, 0.0 + 0im, 0.5)
        numsA, denA = shared_denominator_pade(rj, 7)
        @test den == denA && nums == numsA                            # SPD.3d fallback == cell A
    end

    @testset "SPD.4 determinism under eps-perturbation" begin
        for (f, y0, h, m) in (((z, y) -> [y[2], -y[1]], ComplexF64[1, 0], 0.7, 15),
                              ((z, y) -> (dd = y[1]-y[2]; a = 4/dd^3; [y[3], y[4], a, -a]),
                               ComplexF64[1, -1, 0, 0], 0.1, 15))
            jets = rjets(f, 0.0 + 0im, y0, h, 2m + 2)
            c0 = shared_pade_select(jets, m, f, 0.0 + 0im, h; diagnostics = true)[3].chosen_cell
            pert = [[cj * (1 + 50 * eps() * (isodd(k) ? 1 : -1)) for (k, cj) in enumerate(j)] for j in jets]
            c1 = shared_pade_select(pert, m, f, 0.0 + 0im, h; diagnostics = true)[3].chosen_cell
            @test c0 == c1                                            # SPD.4 stable
        end
    end

    @testset "SPD.5 complex step h handled" begin
        harm = (z, y) -> [y[2], -y[1]]
        h = 0.5 * cis(0.4)                                            # genuine complex wedge step
        res = shared_pade_select(rjets(harm, 0.0 + 0im, ComplexF64[1, 0], h, 32), 15,
                                 harm, 0.0 + 0im, h; diagnostics = true)
        @test res[3].chosen_cell in (:square, :diagonal)              # ran, picked validly
        @test isfinite(res[3].defect_A) || isfinite(res[3].defect_B)  # at least one finite score
        @test res[3].chosen_cell == :square                          # SPD.5 entire ⇒ B even at complex h
    end
end
