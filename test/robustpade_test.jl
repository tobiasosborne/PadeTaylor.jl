# =============================================================================
# robustpade_test.jl — Phase 2 tests for PadeTaylor.RobustPade.robust_pade.
#
# Test plan from DESIGN.md §4 Phase 2.  Eight tests covering the GGT 2013
# Algorithm 2 + Chebfun reweighting port.  Oracle values were captured by
# running Chebfun's `padeapprox.m` (commit 7574c77) under Octave 8.4.0
# (`external/probes/padeapprox-oracle/capture.m` → `test/_oracles.jl`).
# What each case actually pins against that capture (bead padetaylor-044m):
#
#   2.1.1  exp (2,2)       — degrees + COEFFICIENTS a, b (1e-13, closed form
#                             AND the captured oracle; :classical path).
#   2.1.2  exp (20,20)     — degrees (7,7) + COEFFICIENTS a, b (1e-9) +
#                             r(z) function values (1e-14); :svd path.
#   2.1.3  log(1.2-z)      — DEGREES ONLY (10,10) + r(z) values.  The
#                             captured a, b are deliberately NOT asserted
#                             (GGT 2013 §7 ill-posedness at the block
#                             boundary; see the note in the testset).
#   2.1.4  tan(z⁴)         — degrees (20,16) + the captured POLE SET
#                             (16 poles, bidirectional set match, 1e-12).
#   2.1.5  1+z² (1,1)      — degrees (0,0) + COEFFICIENTS a = b = [1].
#   2.1.6  noisy 1/(1-z)   — DEGREES ONLY (0,1); a different RNG
#                             realisation from the Octave capture, so the
#                             captured a, b are not comparable.
#   2.1.7  mutation-proof procedure (documentation only).
#   2.1.8  BigFloat-256 exp (20,20) — type genericity + r(z) accuracy; no
#                             Octave oracle (Octave is Float64-only).
#
# Tolerances are chosen to match `padeapprox.m`'s own behaviour, not
# eyeballed.  The coefficient tolerances admit the known floating-point
# non-uniqueness of QR sign conventions between two LAPACK pipelines.
# =============================================================================

using Test
using LinearAlgebra
using Polynomials: Polynomials
using PadeTaylor: robust_pade, PadeApproximant

# Load Octave-captured oracle data.
include("_oracles.jl")

@testset "RobustPade.robust_pade" begin
    # -------------------------------------------------------------------------
    # 2.1.1 — exp(z), (m, n) = (2, 2): closed-form Padé(2,2) of exp.
    #
    # The closed-form Padé(2,2) of exp(z) is well-known:
    #     P/Q = (1 + z/2 + z²/12) / (1 - z/2 + z²/12).
    # We assert the algorithm produces these coefficients to ~1e-13 rel.
    # -------------------------------------------------------------------------
    @testset "2.1.1 exp(z) Padé(2,2) closed form" begin
        c = Float64.([1 / factorial(big(k)) for k = 0:4])
        P = robust_pade(c, 2, 2)
        @test P.μ == 2
        @test P.ν == 2
        @test P.a ≈ [1.0, 0.5, 1/12]   atol = 1e-13
        @test P.b ≈ [1.0, -0.5, 1/12]  atol = 1e-13
        # Cross-check against the Octave-captured oracle for the same case.
        @test P.a ≈ test_2_1_1_exp_2_2_a atol = 1e-13
        @test P.b ≈ test_2_1_1_exp_2_2_b atol = 1e-13
    end

    # -------------------------------------------------------------------------
    # 2.1.2 — exp(z), (m, n) = (20, 20): diagonal-stripe reduction.
    #
    # exp(z) at order 40 is resolved to machine precision well before the
    # full (20, 20) block; the algorithm hops along the diagonal and
    # settles at (μ, ν) = (7, 7) — GGT 2013 Fig. 2 left panel verbatim.
    # -------------------------------------------------------------------------
    @testset "2.1.2 exp(z) (20,20) reduces to (7,7)" begin
        c = Float64.([1 / factorial(big(k)) for k = 0:40])
        # `method = :svd` is explicit here: bead `padetaylor-txg` changed
        # the F64 default to :classical (which does NOT diagonal-hop, so
        # it would return full (20, 20)).  This test is specifically the
        # GGT 2013 Algorithm 2 reduction invariant; we keep it under :svd.
        P = robust_pade(c, 20, 20; method = :svd)
        # Load-bearing structural assertions: the diagonal-stripe
        # reduction is the GGT 2013 Algorithm 2 invariant.
        @test P.μ == test_2_1_2_exp_20_20_mu  # 7
        @test P.ν == test_2_1_2_exp_20_20_nu  # 7
        @test length(P.a) == length(test_2_1_2_exp_20_20_a)
        @test length(P.b) == length(test_2_1_2_exp_20_20_b)

        # Coefficient match: looser tolerance (1e-9) admits the
        # documented floating-point non-uniqueness of two independent
        # SVD/QR pipelines (Octave/LAPACK vs Julia/LAPACK).  GGT 2013 §7
        # explicitly notes Padé is ill-posed when defect > 0; the (7,7)
        # block boundary is exactly such a regime, so coefficient values
        # vary at the 11th digit between equally-valid implementations
        # while the *function value* (P(z)/Q(z)) is much tighter.
        @test maximum(abs, P.a .- test_2_1_2_exp_20_20_a) < 1e-9
        @test maximum(abs, P.b .- test_2_1_2_exp_20_20_b) < 1e-9

        # Functional-equality test: r(z) = exp(z) to ~16 digits at a
        # well-conditioned point.  This is the load-bearing accuracy
        # claim — the polynomials may differ in coefficient form but
        # the rational function they define agrees to machine precision.
        function _eval_pade(P, z)
            num = sum(P.a[k] * z^(k-1) for k = 1:length(P.a))
            den = sum(P.b[k] * z^(k-1) for k = 1:length(P.b))
            num / den
        end
        for z in (0.1, 0.3, 0.5)
            r_ours = _eval_pade(P, z)
            @test abs(r_ours - exp(z)) < 1e-14
        end
    end

    # -------------------------------------------------------------------------
    # 2.1.3 — log(1.2 - z), (m, n) = (20, 20): branch-cut function.
    #
    # log(1.2 - z) has a branch cut on [1.2, ∞); padeapprox should reduce
    # (20, 20) to (10, 10).  GGT 2013 Fig. 8 documents this behaviour
    # qualitatively.
    # -------------------------------------------------------------------------
    @testset "2.1.3 log(1.2-z) (20,20) reduces to (10,10)" begin
        c = zeros(Float64, 41)
        c[1] = log(1.2)
        for k = 1:40
            c[k+1] = -1.0 / (k * 1.2^k)
        end
        P = robust_pade(c, 20, 20; method = :svd)
        # Load-bearing structural reduction.
        @test P.μ == test_2_1_3_log12_20_20_mu  # 10
        @test P.ν == test_2_1_3_log12_20_20_nu  # 10

        # NOTE — coefficient match deliberately NOT asserted here.
        # Empirically the Octave-output / Julia-output coefficients
        # differ at ~1e-3 absolute, while their rational-function
        # values agree to ~1e-15 over most of the convergence disk.
        # GGT 2013 §7 explicitly warns: Padé approximation is ill-
        # posed when defect > 0.  The (10,10) reduction sits at exactly
        # such a block boundary for log(1.2-z), and two equally-valid
        # SVD/QR pipelines (LAPACK in Octave vs Julia, with different
        # sign conventions and rounding orders) produce numerically
        # distinct coefficient sets for the *same* rational function.
        # Asserting per-coef match would chase floating-point noise;
        # asserting r(z) match is the load-bearing accuracy claim.

        # Functional accuracy: r(z) ≈ log(1.2 - z) at sample points
        # well within the convergence disk (radius 1.2).  Machine
        # precision at the centre; loosens as we approach the branch
        # cut at z = 1.2 (any (10,10) rational misses log there).
        function _eval_pade(P, z)
            num = sum(P.a[k] * z^(k-1) for k = 1:length(P.a))
            den = sum(P.b[k] * z^(k-1) for k = 1:length(P.b))
            num / den
        end
        @test abs(_eval_pade(P, 0.0) - log(1.2))     < 1e-15
        @test abs(_eval_pade(P, 0.3) - log(0.9))     < 1e-15
        @test abs(_eval_pade(P, 0.6) - log(0.6))     < 1e-14
        # At z = 0.9 we are 3/4 of the way to the branch cut at 1.2;
        # accept the documented branch-cut-induced loss of accuracy.
        @test abs(_eval_pade(P, 0.9) - log(0.3))     < 1e-9
    end

    # -------------------------------------------------------------------------
    # 2.1.4 — tan(z⁴), (m, n) = (20, 20): Froissart-doublet removal.
    #
    # GGT 2013 Fig. 6 demonstrates that the robust algorithm removes 4
    # Froissart doublets from the (20, 20) approximation of tan(z⁴),
    # yielding (μ, ν) = (20, 16).  We assert (μ, ν) match exactly when
    # the same FFT-derived coefficients padeapprox saw are passed in, AND
    # that the 16 surviving poles match the captured `padeapprox.m` pole
    # set (`poles = roots(b(end:-1:1))`, padeapprox.m:150) as a SET.
    # -------------------------------------------------------------------------
    @testset "2.1.4 tan(z⁴) (20,20) removes Froissart doublets" begin
        # Use the Octave-FFT-derived coefficients to stay byte-exact with
        # what padeapprox saw.  The function `tan(z⁴)` has only every 4th
        # coefficient nonzero, starting at z⁴.
        c = real.(test_2_1_4_tan_z4_20_20_coefs)
        # Froissart doublet removal is GGT 2013 Algorithm 2 specific.
        P = robust_pade(c, 20, 20; method = :svd)
        @test P.μ == test_2_1_4_tan_z4_20_20_mu  # 20
        @test P.ν == test_2_1_4_tan_z4_20_20_nu  # 16

        # Pole-set match (bead padetaylor-044m).  Octave's `roots` and
        # Polynomials.jl's `roots` order the roots differently, so match
        # as a set: every captured pole has one of ours within `tol` and
        # vice versa, with equal cardinality.  Measured agreement on
        # 2026-08-23 (Julia 1.12, OpenBLAS LAPACK): max nearest-neighbour
        # distance 1.35e-14 in both directions.  1e-12 is pinned — ~100×
        # headroom over the measurement for LAPACK / companion-matrix
        # eigensolver variation, and still 4 orders tighter than the 1e-8
        # that would let a Froissart-doublet residue (|z| ~ 1e-7 shifts)
        # slip through.
        ours = Polynomials.roots(Polynomials.Polynomial(P.b))
        ref  = test_2_1_4_tan_z4_20_20_poles
        @test length(ours) == length(ref) == 16
        tol  = 1e-12
        @test maximum(r -> minimum(o -> abs(o - r), ours), ref) < tol
        @test maximum(o -> minimum(r -> abs(o - r), ref), ours) < tol
        # Ground-truth anchor for the oracle itself.  These are the
        # APPROXIMANT's poles, not tan(z⁴)'s: the inner ring of 8 sits at
        # |z| = 1.1195172 vs the true nearest singularities (π/2)^(1/4) =
        # 1.1195151 (2.1e-6 off — the (20,16) Padé resolves the nearest
        # ring to ~6 digits), while the outer ring at |z| = 1.4993 is an
        # approximant artefact that does NOT sit on the true second ring
        # (3π/2)^(1/4) = 1.4734 and is therefore not pinned to closed form.
        moduli = sort(abs.(ref))
        @test all(abs.(moduli[1:8] .- (π / 2)^(1 / 4)) .< 1e-5)
        @test all(abs.(moduli[9:16] .- 1.4992939919206) .< 1e-12)   # captured
    end

    # -------------------------------------------------------------------------
    # 2.1.5 — 1 + z², (m, n) = (1, 1): defect-1 ill-posedness.
    #
    # GGT 2013 §7 first equation: f(z) = 1 + z², r₁₁ = 1 (defect 1).  The
    # exact-arithmetic Padé(1,1) is the constant 1.  Our algorithm should
    # produce μ = ν = 0, a = b = [1].
    # -------------------------------------------------------------------------
    @testset "2.1.5 (1+z²) Padé(1,1) defect-1 collapse" begin
        c = Float64[1, 0, 1]
        P = robust_pade(c, 1, 1)
        @test P.μ == 0
        @test P.ν == 0
        @test P.a ≈ [1.0]  atol = 1e-14
        @test P.b ≈ [1.0]  atol = 1e-14
    end

    # -------------------------------------------------------------------------
    # 2.1.6 — noisy 1/(1-z), tol = 1e-5: noise-thresholded recovery.
    #
    # GGT 2013 Fig. 4 right panel; with tol above the noise level, the
    # algorithm recovers (μ, ν) = (0, 1) for c_k = 1 + ε·N(0,1).
    #
    # The exact a/b values depend on the noise realisation (Octave's
    # randn vs Julia's randn); we assert only the (μ, ν) outcome on a
    # Julia-side noise sample at the same scale.
    # -------------------------------------------------------------------------
    @testset "2.1.6 noisy 1/(1-z), tol=1e-5 recovers (0,1)" begin
        # Generate noise with a fixed Julia seed; not byte-equal to
        # Octave's, but the *outcome* of (μ, ν) should match.
        using Random
        rng = MersenneTwister(42)
        n_coefs = 21
        c = ones(Float64, n_coefs) .+ 1e-6 .* randn(rng, n_coefs)
        # Noise-thresholded recovery uses the tol knob inside the SVD
        # diagonal-hop loop; classical has no tol knob, so this is an
        # SVD-only test.
        P = robust_pade(c, 10, 10; tol = 1e-5, method = :svd)
        @test P.μ == 0
        @test P.ν == 1
    end

    # -------------------------------------------------------------------------
    # 2.1.7 — Mutation-proof procedure (documented manual procedure).
    #
    # Per CLAUDE.md Rule 4: every load-bearing test gets a documented
    # mutation that RED's it.
    #
    # The QR-reweighting block (`padeapprox.m` lines 111–117) is the
    # specific addition GGT 2013 Algorithm 2 does *not* document but is
    # part of Chebfun's reference impl.  Mutation: in src/RobustPade.jl
    # replace the QR-reweighting branch
    #
    #   F = qr(adjoint(C * D))
    #   b = D * F.Q[:, n+1]
    #   b ./= norm(b)
    #
    # with the SVD null vector directly:
    #
    #   b = vec(Vt[end, :])
    #
    # and run this file standalone.  2.1.3 has NO coefficient-match
    # assertion (see the header and the note in that testset), so it is
    # NOT where this mutation bites.  Verified by hand 2026-08-23 (bead
    # padetaylor-044m), 5 RED of 66: 2.1.2's two coefficient-match
    # assertions (1e-9) and, in 2.1.4, the ν == 16 degree assertion plus
    # both pole-set assertions — the un-reweighted null vector no longer
    # hits the exact zeros of b, so the Froissart-doublet trim is lost.
    # 2.1.1/2.1.3/2.1.5/2.1.6/2.1.8 stay GREEN.  This mutation is now in
    # the periodic gate catalogue, `test/mutation/run_mutation_gate.jl`.
    # -------------------------------------------------------------------------
    # Documented only; the executable form lives in the mutation gate.

    # -------------------------------------------------------------------------
    # 2.1.8 — BigFloat at precision = 256 bits: exp(z) (20,20).
    #
    # Type-genericity test, with a precision-tier subtlety: at 256-bit
    # precision the Padé(20, 20) of exp does NOT reduce to (7, 7).
    # The Float64 reduction was a precision-driven phenomenon — the
    # high-order Taylor coefficients `1/k!` for k > 16 are smaller than
    # `1e-14 · ‖c‖₂`, so the Float64 SVD detects rank deficiency.  At
    # 256 bits, `1/40! ≈ 1.2e-48` is comfortably above the BigFloat
    # default tolerance `2^(-246) ≈ 1e-74`, so no rank deficiency is
    # detected and (μ, ν) = (20, 20) survives.
    #
    # **This is correct behaviour** — arb-prec doing what it should:
    # higher precision admits more accurate Padé approximants.
    #
    # The load-bearing assertions: type genericity (eltype(P.a) ===
    # BigFloat) and functional accuracy at BigFloat precision.
    # -------------------------------------------------------------------------
    @testset "2.1.8 BigFloat-256 exp(z) (20,20) — full block, no reduction" begin
        setprecision(BigFloat, 256) do
            c = [BigFloat(1) / factorial(big(k)) for k = 0:40]
            P = robust_pade(c, 20, 20)
            # No precision-driven reduction at 256 bits.
            @test P.μ == 20
            @test P.ν == 20
            # Type genericity: the dispatch routes BigFloat through
            # GenericLinearAlgebra.svd (per ADR-0002) and preserves
            # element type all the way out.
            @test eltype(P.a) === BigFloat
            @test eltype(P.b) === BigFloat
            # Functional accuracy at BigFloat precision: r(z) = exp(z)
            # to far better than Float64 ε.  Padé(20,20) of exp gives
            # ~80 digits at z = 0.5 by the diagonal Padé error formula
            # |exp(z) - r(z)| = (n!)² / ((2n)!(2n+1)!) z^(2n+1) e^z; for
            # n = 20, z = 0.5 this is ~1e-50.  We accept 1e-40 to admit
            # any roundoff in the SVD/QR pipeline at BigFloat precision.
            function _eval_pade(P, z)
                num = sum(P.a[k] * z^(k-1) for k = 1:length(P.a))
                den = sum(P.b[k] * z^(k-1) for k = 1:length(P.b))
                num / den
            end
            for z in (BigFloat("0.1"), BigFloat("0.3"), BigFloat("0.5"))
                r_ours = _eval_pade(P, z)
                @test abs(r_ours - exp(z)) < BigFloat(1e-40)
            end
        end
    end

    # -------------------------------------------------------------------------
    # 2.1.9 — non-finite coefficients / tol are rejected (bug padetaylor-lbqb).
    #
    # Pre-fix, `tol = Inf` or any `Inf` coefficient made the r ≡ 0 test at
    # GGT step 2 vacuously true and robust_pade SILENTLY returned a = [0],
    # b = [1] (a truthful-looking lie, Rule 1); NaN coefficients went into
    # LAPACK.  Now: ArgumentError naming the overflow and the remedy.
    # -------------------------------------------------------------------------
    # MUTATION-PROOF RECORD (CLAUDE.md Rule 4), bead padetaylor-lbqb.
    #   Edit (src/RobustPade.jl, robust_pade): replace the two guards
    #     `all(isfinite, c) || throw(ArgumentError(…))` and
    #     `(0 ≤ tol < Inf) || throw(ArgumentError(…))`
    #   with `true || throw(…)`; run `julia --project=. test/robustpade_test.jl`.
    #   Result: 22 of 61 RED, all in 2.1.9 (Inf/-Inf/NaN coefficient and
    #   Inf/negative/NaN tol, both methods): the zero approximant came back
    #   silently or LAPACK threw a non-ArgumentError.  Restored → 61/61 GREEN.
    # -------------------------------------------------------------------------
    @testset "2.1.9 non-finite input validation (lbqb)" begin
        c = Float64.([1 / factorial(big(k)) for k = 0:4])
        for method in (:classical, :svd)
            for bad in (Inf, -Inf, NaN)
                cb = copy(c); cb[2] = bad
                e = try; robust_pade(cb, 2, 2; method = method); catch e; e; end
                @test e isa ArgumentError
                msg = sprint(showerror, e)
                @test occursin("overflowed", msg) && occursin("reduce h or order", msg)
            end
            for badtol in (Inf, -1e-14, NaN)
                e = try; robust_pade(c, 2, 2; tol = badtol, method = method); catch e; e; end
                @test e isa ArgumentError
                @test occursin("tol", sprint(showerror, e))
            end
        end
        # tol = 0 is a legal (exact-arithmetic) tolerance — must not throw.
        P = robust_pade(c, 2, 2; tol = 0.0, method = :svd)
        @test P.a ≈ [1.0, 0.5, 1/12] atol = 1e-13
    end
end
