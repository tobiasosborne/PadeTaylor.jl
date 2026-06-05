# test/corpus_robust_pade_test.jl
#
# ============================================================================
# Ground-truth corpus: ROBUST PADÉ core (bead padetaylor-7pav,
# epic padetaylor-p1v0).
#
# WHAT THIS FILE PINS
# -------------------
# The algebraic heart of PadeTaylor.jl is `robust_pade` (GGT 2013 Algorithm 2
# + FW 2011 classical Toeplitz, dispatched by element type — ADR-0005).  This
# file nails six independent claims about it against closed-form / mpmath
# oracles, plus one cross-module byte-identity guard:
#
#   CRP.1  exp-77-exact-rational   (7,7) Padé of exp == exact binomial form
#   CRP.2  tan-33-exact-pade       (3,3) Padé of tan: p=z-z^3/15, q=1-2z^2/5
#   CRP.3  exp-iz-complex-pade     Complex{F64} (7,7): r(1) == exp(i)
#   CRP.4  adr0002-rank-boundary-svd  σ_2 recovery, F64 LAPACK vs BF256 Jacobi
#   CRP.5  geom-noisy-1on1z        tol>noise ⇒ collapse to type (0,1), root≈1
#   CRP.6  B22 cross-module        _chebyshev_D1 byte-identity across 3 copies
#
# GROUND TRUTH each case embodies
# -------------------------------
# Recipes verbatim from docs/test_corpus/01_corpus_catalogue.md "Full
# candidate records":  exp-77 (md:1459-1474), tan-33 (md:1476-1491),
# exp-iz (md:132, exp-iz-77-complex-bigfloat md:1544-1559), adr0002
# (md:1527-1542), geom-noisy (geom-series-reduction md:1493-1508).
# B22 is the confirmed 3-copy duplication of the Chebyshev differentiation
# matrix (capability-map B22) — BVP._chebyshev_D1, VectorBVP._chebyshev_D1,
# Laplace2D._cheb_D1, each declared "copied verbatim" in its docstring.
#
# DE-DUPLICATION vs the existing suite (CLAUDE.md Rule 3, skepticism):
#   * test/robustpade_test.jl 2.1.1 already pins exp(2,2); 2.1.2 pins the
#     (20,20)→(7,7) GGT *reduction* under :svd.  CRP.1 is NOT a duplicate:
#     it pins the (7,7) request from an order-14 jet against the EXACT
#     binomial coefficients (not the Octave oracle, not a reduction), and
#     it pins the no-spurious-pole (entire-function) claim — neither tested.
#   * 2.1.6 pins the noisy-geom (μ,ν)=(0,1) outcome on a Julia-RNG sample.
#     CRP.5 strengthens it: a DETERMINISTIC mpmath noise vector (no RNG) and
#     it additionally pins the denominator ROOT ≈ 1 (the genuine pole) plus
#     the exact noise-free reduction b=[1,-1].
#   * 2.1.8 pins BigFloat-256 exp(20,20).  CRP.4 is orthogonal: the ADR-0002
#     σ_2 rank-boundary stress case, F64-vs-BF256 SV resolution.
#   * tan(z^4) Froissart (2.1.4) is ALREADY pinned from the frozen Octave
#     oracle; per the bead brief we do not re-pin an identical case here.
#   * No existing test asserts the binomial closed form, the exp-iz complex
#     path, the ADR-0002 σ_2, or the _chebyshev_D1 cross-copy identity.
#
# TOLERANCES (justified by method/precision, never fitted):
#   * CRP.1  classical F64 Toeplitz \ on the κ≈9.1e12 7×7 system: observed
#     max coeff error 1.18e-11; we assert ≤ 2e-11 (the natural Float64 floor
#     for this conditioning, cf. 2.1.2 which uses 1e-9 for the same block).
#   * CRP.2  rational coeffs recovered to ~4e-17 (machine ε) — assert 1e-14.
#   * CRP.3  |r(1)-exp(i)| observed 1.57e-16 (≈ ε) — assert 1e-14.
#   * CRP.4  BF256 Jacobi recovers σ_2 to 2.2e-22 RELATIVE (ADR-0002 claim);
#     assert rtol 1e-18.  F64 LAPACK σ_2 to 0.29% relative — assert rtol 5e-3
#     (DOCUMENTED observed behaviour, see the FINDING note in CRP.4).
#   * CRP.5  Float64 SVD reduction: denominator root within 1e-6 of z=1
#     (noise floor 1e-9 propagated through the (0,1) solve).
#   * CRP.6  exact byte-identity (==), Float64 and BigFloat.
#
# MUTATION-PROOF (Rule 4): see the footer comment block — each load-bearing
# assertion was verified RED under a deliberate oracle/impl perturbation,
# then restored byte-for-byte (src/ confirmed clean via `git diff`).
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor: robust_pade, PadeApproximant
using PadeTaylor.LinAlg: pade_svd
using LinearAlgebra: norm
using Polynomials: Polynomial, roots

include(joinpath(@__DIR__, "_oracle_corpus_robust_pade.jl"))

# Lower-triangular Toeplitz block C (n×(n+1)) GGT 2013 eq.2.6 — the SAME
# bottom-n-rows block RobustPade.jl builds as Z[(m+2):(m+n+1), :].  Replicated
# here (not reaching into the private helper) so CRP.4 pins the SVD on an
# independently-constructed matrix.
function _toep_block(c::AbstractVector{T}, m::Int, n::Int) where {T}
    Z = zeros(T, m + n + 1, n + 1)
    for j = 1:(n + 1), i = j:(m + n + 1)
        Z[i, j] = c[i - j + 1]
    end
    return Z[(m + 2):(m + n + 1), :]
end

@testset "Corpus: robust Padé core (CRP)" begin

    # ----------------------------------------------------------------------
    # CRP.1  exp-77-exact-rational : (7,7) Padé of the order-14 exp jet
    # equals the exact binomial closed form (already b[1]=1 normalised),
    # and the denominator has NO root inside the unit disk (exp is entire).
    # ----------------------------------------------------------------------
    @testset "CRP.1  exp(z) (7,7) == exact binomial Padé, no spurious poles" begin
        c = Float64[1 / factorial(big(k)) for k = 0:14]
        P = robust_pade(c, 7, 7)                       # F64 default :classical
        @test P.μ == 7
        @test P.ν == 7
        @test length(P.a) == 8 && length(P.b) == 8
        # Exact binomial coefficients, b[1]=1 convention (no rescale needed —
        # the closed form is already monic-at-zero).  κ(T_7×7)≈9.1e12 ⇒ 2e-11.
        @test maximum(abs, P.a .- CRP_EXP77_A) ≤ 2e-11
        @test maximum(abs, P.b .- CRP_EXP77_B) ≤ 2e-11
        # a_1 = 1/2, b_1 = -1/2 exactly (the analytic spot-check in md:1474).
        @test isapprox(P.a[2],  0.5; atol = 2e-11)
        @test isapprox(P.b[2], -0.5; atol = 2e-11)
        # Entire function ⇒ no Padé pole inside the working disk |z|≤5; the
        # 7 conjugate-pair roots all sit at |z|≈9.94..12.1 (md:1465).
        rts = roots(Polynomial(P.b))
        @test all(abs(z) > 5.0 for z in rts)
        @test isapprox(minimum(abs, rts), 9.9436; atol = 1e-3)
    end

    # ----------------------------------------------------------------------
    # CRP.2  tan-33-exact-pade : p=z-z^3/15, q=1-2z^2/5; root at sqrt(5/2).
    # ----------------------------------------------------------------------
    @testset "CRP.2  tan(z) (3,3) Padé: p=z-z^3/15, q=1-2z^2/5, root √(5/2)" begin
        # tan series c = [0, 1, 0, 1/3, 0, 2/15, 0] (DLMF 4.14).
        c = Float64[0, 1, 0, 1/3, 0, 2/15, 0]
        P = robust_pade(c, 3, 3)
        @test P.μ == 3
        # Numerator p = z - z^3/15 exactly (rational, ~machine ε).
        @test maximum(abs, P.a .- CRP_TAN33_A) ≤ 1e-14
        # Denominator q = 1 - 2z^2/5; the classical path leaves a trailing
        # z^3 coefficient that is exactly zero (q has no cubic term).
        @test length(P.b) ≥ 3
        @test isapprox(P.b[1], CRP_TAN33_B[1]; atol = 1e-14)
        @test isapprox(P.b[2], CRP_TAN33_B[2]; atol = 1e-14)
        @test isapprox(P.b[3], CRP_TAN33_B[3]; atol = 1e-14)
        @test abs(P.b[end]) ≤ 1e-14          # trailing cubic term is ~0
        # Denominator root at z = √(5/2) = 1.5811…, NOT the true pole π/2.
        # Evaluate q(√(5/2)) ≈ 0 directly (cleaner than root-finding the
        # zero-padded polynomial); error from π/2 is the catalogue 1.034e-2.
        qval = sum(P.b[k] * CRP_TAN33_ROOT^(k - 1) for k = 1:length(P.b))
        @test isapprox(qval, 0.0; atol = 1e-13)
        @test isapprox(CRP_TAN33_ROOT, 1.5811388300841898; atol = 1e-13)
        @test isapprox(CRP_TAN33_ROOT - π/2, 1.034e-2; atol = 1e-4)
    end

    # ----------------------------------------------------------------------
    # CRP.3  exp-iz-complex-pade : Complex{Float64} (7,7), r(1) == exp(i).
    # Exercises the complex classical/SVD path (ADR-0002 Complex{T} row).
    # ----------------------------------------------------------------------
    @testset "CRP.3  exp(iz) Complex (7,7): r(1) ≈ exp(i) to machine ε" begin
        c = Complex{Float64}[(im)^k / factorial(big(k)) for k = 0:14]
        P = robust_pade(c, 7, 7)               # Complex{F64} default :classical
        @test eltype(P.a) === Complex{Float64}
        @test P.μ == 7 && P.ν == 7
        r1 = sum(P.a) / sum(P.b)               # evaluate r at z = 1
        @test isapprox(r1, CRP_EXPI_R1; atol = 1e-14)   # |err| ~ 1.6e-16
        # Poles of exp(iz)-Padé are the real-exp poles rotated by -i; the
        # closest sits at ≈ 0 - 9.94i (md:1559).
        rts = roots(Polynomial(P.b))
        zc = rts[argmin(abs.(rts))]
        @test isapprox(zc, complex(0.0, -9.9436); atol = 1e-2)
    end

    # ----------------------------------------------------------------------
    # CRP.4  adr0002-rank-boundary-svd : 10×11 Toeplitz of
    # c_k = 1 + 1e-13·(1/2)^k.  σ_2 sits 3.9% above the F64 GGT threshold.
    # BF256 Jacobi recovers σ_2 to RELATIVE accuracy (the ADR-0002 claim);
    # we ALSO document what F64 LAPACK does on the same matrix.
    #
    # FINDING (reported in the final summary): contrary to the ADR-0002
    # *worst-case* narrative, on THIS matrix F64 LAPACK does NOT misclassify
    # σ_2 below threshold — it resolves σ_2 to 0.29% relative and keeps
    # rank=2.  The catalogue verify_note (md:1542) explicitly anticipates
    # this: "If LAPACK does not misclassify in practice, reframe the test as
    # a precision comparison."  We pin the OBSERVED invariant: F64 keeps the
    # genuine SV but at low (absolute-accuracy) relative precision, while
    # BF256 Jacobi delivers full relative accuracy — exactly the precision
    # GAP ADR-0002 is built on, even though the rank flip does not trigger.
    # ----------------------------------------------------------------------
    @testset "CRP.4  ADR-0002 σ_2 rank-boundary: F64 LAPACK vs BF256 Jacobi" begin
        # Build c_k at each precision; assert the mpmath σ_2 is recovered.
        c64 = Float64[1 + 1e-13 * (0.5)^k for k = 0:20]
        C64 = _toep_block(c64, 10, 10)
        _, S64, _ = pade_svd(C64; full = false)          # LAPACK Demmel-Kahan
        thr = 1e-14 * norm(c64[1:21])
        # OBSERVED: F64 resolves σ_1 to full precision and σ_2 ABOVE threshold
        # (rank kept = 2), but σ_2 only to ~0.29% relative (DK absolute floor).
        @test isapprox(S64[1], CRP_SVD_SIGMA1; rtol = 1e-12)
        @test count(s -> s > thr, S64) == 2              # rank NOT misclassified
        @test S64[2] > thr                               # σ_2 survives at F64
        @test isapprox(S64[2], CRP_SVD_SIGMA2; rtol = 5e-3)   # only ~0.29% rel
        # σ_3 at F64 is a LAPACK noise artefact (~7e-16), NOT the true ~5e-51:
        # this is the absolute-accuracy fingerprint DK leaves behind.
        @test S64[3] > 1e-17 && S64[3] < thr

        setprecision(BigFloat, 256) do
            cb = BigFloat[1 + BigFloat(10)^-13 * (BigFloat(1) / 2)^k for k = 0:20]
            Cb = _toep_block(cb, 10, 10)
            _, Sb, _ = pade_svd(Cb; full = false)        # GenericLinAlg Jacobi
            # ADR-0002 relative-accuracy claim: σ_2 to 2^-246·σ_2 ≈ machine ε
            # of the precision (observed 2.2e-22 rel) — compare in BigFloat
            # against the full 25-digit mpmath oracle, rtol 1e-18.
            sig2 = parse(BigFloat, CRP_SVD_SIGMA2_STR)
            @test abs(Sb[2] - sig2) / sig2 < BigFloat("1e-18")
            # σ_3 at BF256 is the TRUE exact-arithmetic zero (~7.7e-77), far
            # below σ_2 — Jacobi sees the real null space, unlike LAPACK.
            @test Float64(Sb[3]) < 1e-70
        end
    end

    # ----------------------------------------------------------------------
    # CRP.5  geom-noisy-1on1z : 1/(1-z) series + tiny noise.  tol above the
    # noise level collapses the (10,10) block to type (0,1) with the genuine
    # pole at z≈1; the exact noise-free series reduces to b=[1,-1] precisely.
    # ----------------------------------------------------------------------
    @testset "CRP.5  noisy 1/(1-z): tol>noise ⇒ type (0,1), denom root ≈ 1" begin
        # Deterministic mpmath noise (no RNG); noise scale 1e-9 ≪ tol 1e-5.
        P = robust_pade(CRP_NOISY_C, 10, 10; tol = 1e-5, method = :svd)
        @test P.μ == 0
        @test P.ν == 1
        @test length(P.b) == 2
        # Genuine simple pole of 1/(1-z): the denominator root -b[1]/b[2] ≈ 1.
        root = -P.b[1] / P.b[2]
        @test isapprox(root, 1.0; atol = 1e-6)        # noise 1e-9 → ~1e-9 here
        # Exact noise-free series (c_k = 1) reduces to b = [1, -1] precisely.
        Pexact = robust_pade(ones(Float64, 21), 10, 10; tol = 1e-14, method = :svd)
        @test Pexact.μ == 0 && Pexact.ν == 1
        @test isapprox(Pexact.a, [1.0]; atol = 1e-13)
        @test isapprox(Pexact.b, [1.0, -1.0]; atol = 1e-13)
    end

    # ----------------------------------------------------------------------
    # CRP.6  B22 cross-module byte-identity : the Chebyshev first-derivative
    # matrix is duplicated verbatim in THREE modules (capability-map B22).
    # Each docstring claims "copied verbatim from BVP._chebyshev_D1"
    # (src/Laplace2D.jl:287, src/VectorBVP.jl:369).  Guard the claim against
    # silent drift: build the same Chebyshev-2 nodes and assert all three
    # builders agree ELEMENTWISE.  A disagreement is a real finding.
    # ----------------------------------------------------------------------
    @testset "CRP.6  B22 _chebyshev_D1 byte-identity across 3 module copies" begin
        D1_bvp  = PadeTaylor.BVP._chebyshev_D1
        D1_vbvp = PadeTaylor.VectorBVP._chebyshev_D1
        D1_lap  = PadeTaylor.Laplace2D._cheb_D1
        for N in (4, 8, 12)
            # Chebyshev-2 extrema, the convention BVP.jl uses (src/BVP.jl:272).
            t = Float64[cos(Float64(j) * pi / Float64(N)) for j = 0:N]
            A = D1_bvp(t, Float64, N)
            B = D1_vbvp(t, Float64, N)
            C = D1_lap(t, Float64, N)
            @test A == B           # exact byte-identity, not isapprox
            @test A == C
            @test size(A) == (N + 1, N + 1)
            # Spectral consistency: D·1 = 0 (the row-sum identity each builder
            # imposes) — pins that all three really are differentiation matrices.
            @test maximum(abs, A * ones(Float64, N + 1)) ≤ 1e-12
        end
        # BigFloat path too (the builders are generic in T<:AbstractFloat).
        setprecision(BigFloat, 128) do
            N = 6
            t = BigFloat[cos(BigFloat(j) * BigFloat(pi) / BigFloat(N)) for j = 0:N]
            @test D1_bvp(t, BigFloat, N) == D1_vbvp(t, BigFloat, N)
            @test D1_bvp(t, BigFloat, N) == D1_lap(t, BigFloat, N)
        end
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# src/ confirmed clean via `git diff --stat src/` after each restore.
#
# M1 (CRP.1, oracle): negate CRP_EXP77_B[2] in the oracle helper (−0.5 → +0.5).
#   Result: CRP.1 went RED at the b-coefficient assertion (max|b−exact| jumped
#   to 1.0) and at the a_1/b_1 spot-check; all other testsets stayed GREEN.
#   Restored the literal.
#
# M2 (CRP.2, impl): in src/RobustPade.jl `classical_pade_diagonal`, perturb the
#   numerator recurrence `s += b_full[j+1] * cv[k-j+1]` → `* 1.0000001 *`.
#   Result: CRP.2 RED at `max|a − [0,1,0,−1/15]|` (got ~1e-7) AND CRP.1 RED;
#   the rational tan numerator and exp numerator both route through this line.
#   Restored src/RobustPade.jl byte-for-byte.
#
# M3 (CRP.4, oracle): change CRP_SVD_SIGMA2_STR mantissa …758889598953 →
#   …759000000000 (perturb past the 16th digit, where only the BigFloat
#   comparison can see it).  Result: CRP.4 RED at the BF256 rtol-1e-18
#   assertion (the F64 rtol-5e-3 assertion stayed GREEN — two orders looser),
#   confirming the BF256 branch pins σ_2 to full relative accuracy while F64
#   does not.  Restored.
#
# M4 (CRP.6, impl): in src/Laplace2D.jl `_cheb_D1`, flip the off-diagonal sign
#   `iseven(i + j) ? T(1) : T(-1)` → `iseven(i + j) ? T(-1) : T(1)`.
#   Result: CRP.6 RED at `A == C` for every N (the Laplace2D copy diverged
#   from the BVP source) — exactly the silent-drift this guard exists to
#   catch.  The D·1=0 row-sum check stayed GREEN (sign flip preserves it),
#   proving the byte-identity assertion is strictly stronger.  Restored
#   src/Laplace2D.jl byte-for-byte.
#
# Run standalone with:
#   julia --project=. test/corpus_robust_pade_test.jl
# ============================================================================
