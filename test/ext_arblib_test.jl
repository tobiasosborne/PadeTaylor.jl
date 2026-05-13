# test/ext_arblib_test.jl -- Phase 8 / bead `padetaylor-jhq` tests.
#
# Validates the `PadeTaylorArblibExt` package extension per ADR-0003 +
# ADR-0002: `pade_svd(::Matrix{Arb})` routes through `BigFloat`, lifts
# results back to `Arb` with radius zero.
#
# Reference: docs/adr/0003-extensions-pattern.md (extension design);
# docs/adr/0002-bigfloat-svd-via-genericlinalg.md (radius-discard
# caveat); ext/PadeTaylorArblibExt.jl (the adapter).
#
# Test groups:
#   AB.1.1  pade_svd(::Matrix{Arb}) ≡ pade_svd(BigFloat.(A))
#   AB.1.2  pade_svd(::Matrix{Acb}) ≡ pade_svd(Complex{BigFloat}.(A))
#   AB.2.1  SVD round-trip A ≈ U·diag(S)·Vt at high precision
#   AB.3.1  full=true gives n×n Vt (load-bearing for RobustPade's
#           null-vector path)
#   AB.4.1  Mutation-proof commentary

using Test
using PadeTaylor
using PadeTaylor: LinAlg
using Arblib

@testset "ArblibExt (Phase 8): pade_svd over Arb / Acb element types" begin

    setprecision(BigFloat, 256)

    @testset "AB.1.1: pade_svd(::Matrix{Arb}) ≡ BigFloat path" begin
        # Build a small well-conditioned matrix in Arb + BigFloat in
        # parallel.  Compare SVD outputs.  The Arb path discards radii,
        # so its results should match the BigFloat-mid-point path bit-for-bit.
        A_bf = BigFloat[1 2 3; 4 5 6; 7 8 10]
        A_arb = Arblib.Arb.(A_bf)

        U_bf,  S_bf,  Vt_bf  = LinAlg.pade_svd(A_bf)
        U_arb, S_arb, Vt_arb = LinAlg.pade_svd(A_arb)

        @test eltype(U_arb)  === Arblib.Arb
        @test eltype(S_arb)  === Arblib.Arb
        @test eltype(Vt_arb) === Arblib.Arb
        @test size(U_arb) == size(U_bf)
        @test size(S_arb) == size(S_bf)
        @test size(Vt_arb) == size(Vt_bf)

        # Mid-point comparison: lifting Arb→BigFloat→Arb is identity
        # (radii are zero by construction).
        @test BigFloat.(S_arb) == S_bf
        @test all(BigFloat.(U_arb) .== U_bf)
        @test all(BigFloat.(Vt_arb) .== Vt_bf)
    end

    @testset "AB.1.2: pade_svd(::Matrix{Acb}) ≡ Complex{BigFloat} path" begin
        A_cbf = Complex{BigFloat}[1+im 2 3-im; 4 5+2im 6; 7 8-im 10+im]
        A_acb = Arblib.Acb.(A_cbf)

        U_cbf, S_cbf, Vt_cbf = LinAlg.pade_svd(A_cbf)
        U_acb, S_acb, Vt_acb = LinAlg.pade_svd(A_acb)

        @test eltype(U_acb)  === Arblib.Acb
        @test eltype(S_acb)  === Arblib.Arb       # singular values real
        @test eltype(Vt_acb) === Arblib.Acb
        @test all(Complex{BigFloat}.(U_acb) .== U_cbf)
        @test all(BigFloat.(S_acb)         .== S_cbf)
        @test all(Complex{BigFloat}.(Vt_acb) .== Vt_cbf)
    end

    @testset "AB.2.1: SVD round-trip A ≈ U·diag(S)·Vt for Arb input" begin
        # Identity test: rebuild A from the factorisation and compare
        # mid-point to the original Arb matrix.  The reconstruction
        # error reflects the BigFloat-precision Jacobi accuracy.
        A_arb = Arblib.Arb[2 -1 0; -1 2 -1; 0 -1 2]   # symmetric, well-conditioned
        U, S, Vt = LinAlg.pade_svd(A_arb)
        # Build Diagonal matrix from S.
        D = zeros(Arblib.Arb, 3, 3)
        for i in 1:3
            D[i, i] = S[i]
        end
        A_reconstructed = U * D * Vt
        # Compare mid-points to A_arb to within BF-256 precision (~1e-70).
        max_err = maximum(abs.(BigFloat.(A_arb) .- BigFloat.(A_reconstructed)))
        @test max_err < big(1e-60)
    end

    @testset "AB.3.1: full=true returns n×n Vt (null-vector path)" begin
        # RobustPade's GGT 2013 Algorithm 2 needs the null right singular
        # vector of an n×(n+1) Toeplitz; obtained as `transpose(Vt[end, :])`
        # under full=true.  Verify the ext honours full=true.
        A_arb = Arblib.Arb[1 2 3 4; 5 6 7 8; 9 10 11 12]   # 3×4 thin
        U_thin, S_thin, Vt_thin = LinAlg.pade_svd(A_arb; full = false)
        U_full, S_full, Vt_full = LinAlg.pade_svd(A_arb; full = true)

        @test size(Vt_thin) == (3, 4)   # min(m,n) = 3 rows; n = 4 cols
        @test size(Vt_full) == (4, 4)   # full n×n
        # The first min(m,n) rows of full Vt match thin Vt up to sign.
        for j in 1:size(Vt_thin, 1)
            sign_match = BigFloat(Vt_thin[j, 1]) * BigFloat(Vt_full[j, 1]) ≥ 0 ? 1 : -1
            for k in 1:size(Vt_thin, 2)
                @test isapprox(BigFloat(Vt_thin[j, k]),
                               sign_match * BigFloat(Vt_full[j, k]);
                               atol = big(1e-60))
            end
        end
    end

end # @testset ArblibExt

# AB.4.1  Mutation-proof procedure (verified 2026-05-13):
#
#   Mutation I  --  in `pade_svd(::Matrix{Arb})`, downcast `A_bf` to
#     Matrix{Float64} (instead of BigFloat) before the Jacobi call.
#     This is THE correct precision-loss mutation per
#     docs/worklog/001-stages-Z-1-2-handoff.md §"GenericLinearAlgebra
#     svd! piracy": the naïve s/Generic/Linear swap doesn't bite
#     because GLA pirates LinearAlgebra.svd! — the bite comes from
#     forcing precision loss via Float64 downcast.
#     Verified bite: 7 fails across AB.1.1 (S/U/Vt mid-point mismatch
#     at Float64 precision ~1e-15 vs BigFloat ~1e-65), AB.2.1
#     (reconstruction error 1.03e-15 >> 1e-60 BF-256 floor), AB.3.1
#     (Vt sign/precision flip at 3 entries).  Confirms BF-256 dispatch
#     is load-bearing.
#
#   Mutation J  --  in `pade_svd(::Matrix{Arb})`, replace `F.U` with
#     `F.U'` (transpose).  Breaks the SVD identity A = U·S·Vt.
#     Verified bite: 2 fails — AB.1.1 (U mid-points don't match BF
#     path), AB.2.1 (reconstruction error 3.71, huge).
#
# Restoration: both mutations restored before commit per CLAUDE.md Rule 4.
