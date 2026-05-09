# =============================================================================
# linalg_test.jl — Phase 1 tests for PadeTaylor.LinAlg.pade_svd
#
# Test plan from DESIGN.md §4 Phase 1.  The five tests cover:
#   1.1.1  Random Float64 5×6, matches LinearAlgebra.svd exactly.
#   1.1.2  Hilbert-10 Float64 reconstruction.
#   1.1.3  Hilbert-10 BigFloat at precision=256 reconstruction.
#   1.1.4  Rank-deficient [1 2; 2 4] BigFloat — small SV is genuinely zero.
#   1.1.5  Mutation-proof discipline: documented manual procedure (no asserts).
#
# Tolerances are chosen to match the algorithm's documented guarantees, not
# eyeballed from a successful run.
# =============================================================================

using Test
using LinearAlgebra
using GenericLinearAlgebra
using Random
using PadeTaylor: LinAlg

@testset "LinAlg.pade_svd" begin
    # -------------------------------------------------------------------------
    # 1.1.1 — Float64 dispatch matches LinearAlgebra.svd exactly.
    #
    # `pade_svd` over `Matrix{Float64}` must dispatch to LAPACK Demmel-Kahan
    # via `LinearAlgebra.svd` (ADR-0002 §"Why one-sided Jacobi for arb-prec"
    # — at Float64 the GGT tol regime tolerates DK).  We assert equality of
    # the singular values themselves; U and Vt are unique only up to sign
    # per column, so we check that |U·diagm(S)·Vt - A|_F is at machine
    # precision instead of element-wise equality.
    # -------------------------------------------------------------------------
    @testset "1.1.1 Float64 random 5×6 matches LinearAlgebra" begin
        Random.seed!(42)
        A = randn(5, 6)
        U, S, Vt = LinAlg.pade_svd(A)
        @test S ≈ svdvals(A) atol = 1e-14
        # Reconstruction sanity.  `pade_svd` returns the thin form: for a
        # 5×6 input, U is 5×5, S is length 5, Vt is 5×6.
        @test size(U) == (5, 5)
        @test length(S) == 5
        @test size(Vt) == (5, 6)
        @test norm(U * Diagonal(S) * Vt - A) < 1e-13 * norm(A)
    end

    # -------------------------------------------------------------------------
    # 1.1.2 — Hilbert-10 Float64 reconstruction.
    #
    # Hilbert-10 has condition number ~1.6e13 — well-conditioned compared
    # to GGT pathologies; reconstruction at Float64 should be at most a
    # factor of κ·ε ≈ 10⁻³ off in matrix norm.  We accept anything below
    # 1e-9 to leave headroom for matrix-multiplication round-off.
    # -------------------------------------------------------------------------
    @testset "1.1.2 Hilbert-10 Float64 reconstruction" begin
        n = 10
        H = [1.0 / (i + j - 1) for i = 1:n, j = 1:n]
        U, S, Vt = LinAlg.pade_svd(H)
        @test norm(U * Diagonal(S) * Vt - H) < 1e-9 * norm(H)
    end

    # -------------------------------------------------------------------------
    # 1.1.3 — Hilbert-10 BigFloat at 256 bits reconstructs to 1e-65 rel.
    #
    # Load-bearing for arb-prec.  At 256 bits (~77 decimal digits), even a
    # condition-number-1e13 matrix should reconstruct to ~10⁻⁶⁵ — proving
    # the dispatch routes BigFloat to `GenericLinearAlgebra.svd` rather than
    # falling through to a Float64 path.  If the BigFloat path is broken,
    # this fails with `MethodError(svd, ...)` (since stdlib has no generic
    # SVD; ADR-0002 §"`LinearAlgebra.svd` (stdlib)").
    # -------------------------------------------------------------------------
    @testset "1.1.3 Hilbert-10 BigFloat-256 reconstruction" begin
        setprecision(BigFloat, 256) do
            n = 10
            H = [BigFloat(1) / (BigFloat(i) + BigFloat(j) - BigFloat(1))
                 for i = 1:n, j = 1:n]
            U, S, Vt = LinAlg.pade_svd(H)
            recon = U * Diagonal(S) * Vt - H
            rel = norm(recon) / norm(H)
            @test rel < BigFloat(1e-65)
            # Also confirm the SV's themselves are BigFloat, not Float64 —
            # the dispatch must preserve element type through the SVD.
            @test eltype(S) === BigFloat
        end
    end

    # -------------------------------------------------------------------------
    # 1.1.3b — full=true exposes the (n+1)-th column of V for n × (n+1)
    # Toeplitz matrices.
    #
    # Load-bearing for RobustPade: GGT 2013 Algorithm 2 reads
    # `b = V[:, n+1]` as the null right singular vector of an `n × (n+1)`
    # Toeplitz matrix.  In thin SVD V is `(n+1) × n` and the `(n+1)`-th
    # column does not exist; with `full=true`, Vt is `(n+1) × (n+1)`
    # and `Vt[end, :]'` is the null vector.
    # -------------------------------------------------------------------------
    @testset "1.1.3b full=true gives n+1 × n+1 Vt for wide matrix" begin
        # 3 × 4 matrix.  Build it deliberately rank-deficient: the 4th
        # column is in the span of the first three, so the null space is
        # non-trivial.
        A = Float64[1 0 0 1; 0 1 0 2; 0 0 1 3]
        U, S, Vt = LinAlg.pade_svd(A; full = true)
        @test size(U)  == (3, 3)
        @test length(S) == 3
        @test size(Vt) == (4, 4)
        # The 4th column of V (= Vt'[:, 4] = Vt[4, :]) is the null vector.
        # `A * v_null` should be zero.
        v_null = vec(Vt[end, :])
        @test norm(A * v_null) < 1e-13 * norm(A)
    end

    # -------------------------------------------------------------------------
    # 1.1.4 — Rank-deficient [1 2; 2 4] BigFloat — small SV is genuinely zero.
    #
    # The matrix [1 2; 2 4] has rank 1 exactly: column 2 = 2·column 1.  A
    # well-behaved arb-prec SVD must report S[2] far below any reasonable
    # threshold a caller would set; RobustPade defaults to
    # 2^(-precision+10), so a genuine zero SV must be << 2^(-precision+100)
    # at minimum.  We assert S[2]/S[1] < 2^-100 — i.e. at least 30 decimal
    # digits below the largest SV.
    # -------------------------------------------------------------------------
    @testset "1.1.4 BigFloat-256 rank-deficient SV is genuinely zero" begin
        setprecision(BigFloat, 256) do
            A = BigFloat[1 2; 2 4]
            U, S, Vt = LinAlg.pade_svd(A)
            @test length(S) == 2
            @test S[1] > BigFloat(1)               # σ₁ ≈ √(1²+2²+2²+4²) ≈ √25 = 5
            @test S[2] / S[1] < BigFloat(2)^(-100) # genuine zero
        end
    end

    # -------------------------------------------------------------------------
    # 1.1.5 — Mutation-proof discipline (documented manual procedure).
    #
    # CLAUDE.md Rule 4 requires every load-bearing test to be mutation-
    # proven.  For the LinAlg dispatcher the load-bearing assertion is
    # "BigFloat input produces a correct *full-precision* SVD."
    #
    # FRICTION RECORDED (bead: GenericLinearAlgebra-pirates-LinearAlgebra):
    # the obvious mutation `GenericLinearAlgebra.svd → LinearAlgebra.svd`
    # does NOT bite, because `GenericLinearAlgebra` adds methods to
    # `LinearAlgebra.svd!` (intentional piracy in their design).  Both
    # call paths route to the same Jacobi impl when GenericLinearAlgebra
    # is loaded.
    #
    # The mutation that DOES bite (verified during Phase 1):
    #
    #   1. Edit `src/LinAlg.jl` to downcast the input matrix to Float64
    #      before SVD:
    #          F = GenericLinearAlgebra.svd(Matrix{Float64}(A); full = false)
    #          return (Matrix{T}(F.U), Vector{T}(F.S), Matrix{T}(F.Vt))
    #   2. Run `Pkg.test()` and verify tests 1.1.3 AND 1.1.4 RED — the
    #      reconstruction tolerance of 1e-65 fails for 1.1.3, and the
    #      genuine-zero check S[2]/S[1] < 2^-100 fails for 1.1.4 because
    #      Float64 only buys ~16 decimal digits.
    #   3. Restore the original dispatch.
    #   4. Re-run; confirm GREEN.
    #
    # If the mutation does not trigger RED, the test is not load-bearing
    # and must be strengthened.
    # -------------------------------------------------------------------------
    # No assertion — this @testset documents the discipline.  Tests 1.1.3
    # and 1.1.4 carry the empirical weight.
end
