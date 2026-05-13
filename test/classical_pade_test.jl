# =============================================================================
# classical_pade_test.jl — bead `padetaylor-txg` tests (CP.*)
#
# Tests for `RobustPade.classical_pade_diagonal` (FW 2011 §5.1.4 Toeplitz `\`)
# and the `method::Symbol` dispatch on `RobustPade.robust_pade`.
#
# Ground truth: `references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:330-350`:
#
#   - Eqs. (5.4) + (5.5): Toeplitz system for `b_1..b_v`, then numerator
#     by `a_0 = c_0` and `[a_1; …; a_v] = T_low · [b_1; …; b_v] + [c_1; …; c_v]`.
#   - Line 346: singular case → "remove the last row of (5.4) … MATLAB's `\\`
#     solver will then give the solution to the underdetermined system that
#     has the smallest norm".  Our v1 routes singular cases to the existing
#     SVD path instead (which is GGT 2013's principled treatment of the
#     same ill-conditioning), per `padetaylor-txg` acceptance (b).
#   - Line 350: "method of choice … Toeplitz approach and solving the system
#     (5.4) with the backslash operator".
#
# Probe results (worklog 020): at z=10⁴ F64 the classical method achieves
# 6.15e-11 rel-err in 3.45 s, vs SVD's 6.05e-6 rel-err in 17.76 s; classical
# also beats FW's published 2.34e-10 by 3.8×.  At z=30 F64, classical gives
# 1.54e-13 vs SVD's 6.6e-12 (4× better).  Per-step probe: 870× per-step
# accuracy improvement near the pole at z=1.
#
# References:
#   - `docs/worklog/020-classical-pade-toeplitz-backslash.md` — empirical
#     diagnosis, dispatch rationale, FW 2011 ground-truth pointers.
#   - `docs/adr/0005-classical-pade-default-at-float64.md` — design decision.
#   - `references/markdown/FW2011_*.md:330-350` — FW 2011 §5.1.4.
#   - `references/markdown/GGT2013_*.md` — GGT 2013 Algorithm 2 (SVD path).
# =============================================================================

using Test
using LinearAlgebra
using PadeTaylor: robust_pade, PadeApproximant
using PadeTaylor.RobustPade: classical_pade_diagonal

# Internal helper: evaluate Padé as `r(z) = (Σ a[k+1] z^k) / (Σ b[k+1] z^k)`.
function _eval_pade(P, z)
    T = promote_type(eltype(P.a), typeof(z))
    num = zero(T)
    @inbounds for k in length(P.a):-1:1
        num = num * z + P.a[k]
    end
    den = zero(T)
    @inbounds for k in length(P.b):-1:1
        den = den * z + P.b[k]
    end
    num / den
end

@testset "Classical Padé via Toeplitz \\ (FW 2011 §5.1.4)" begin

    # -------------------------------------------------------------------------
    # CP.1.1 — exp(z) closed form (2, 2): both methods recover the textbook
    #          Padé(2,2): `(1 + z/2 + z²/12) / (1 - z/2 + z²/12)`.
    #
    # We assert coefficient match — NOT function-value-vs-exp(z) match, since
    # Padé(2,2) of exp has truncation error `~z⁵/120` (≈ 3·10⁻³ at z=0.5),
    # which is the algorithm's intrinsic accuracy ceiling at that order, not
    # floating-point roundoff.  Function-value-vs-exp tests live in CP.1.2
    # at Padé(15,15) where the truncation error drops below F64 eps.
    # -------------------------------------------------------------------------
    @testset "CP.1.1: exp(z) Padé(2,2) — closed form matches both methods" begin
        c = Float64.([1 / factorial(big(k)) for k = 0:4])

        # Classical via the new entry point.
        P_cls = classical_pade_diagonal(c, 2)
        @test P_cls.μ == 2
        @test P_cls.ν == 2
        @test P_cls.a ≈ [1.0, 0.5, 1/12]   atol = 1e-13
        @test P_cls.b ≈ [1.0, -0.5, 1/12]  atol = 1e-13

        # Same via `robust_pade(...; method=:classical)`.
        P_dispatch = robust_pade(c, 2, 2; method = :classical)
        @test P_dispatch.a ≈ P_cls.a  atol = 1e-15
        @test P_dispatch.b ≈ P_cls.b  atol = 1e-15

        # Same RATIONAL function as `robust_pade(...; method=:svd)`.  Both
        # methods produce numerically identical rationals on this well-
        # conditioned case; we sample the rational at several points to
        # confirm `r_classical(z) ≡ r_svd(z)` within roundoff.  No exp(z)
        # appears here — the algorithm's truncation error at order 2 is
        # ~3·10⁻³ at |z|=0.5, far above floating-point noise.
        P_svd = robust_pade(c, 2, 2; method = :svd)
        for z in (-0.3, -0.1, 0.1, 0.3, 0.5)
            @test _eval_pade(P_cls, z) ≈ _eval_pade(P_svd, z)  atol = 1e-14
        end
    end

    # -------------------------------------------------------------------------
    # CP.1.2 — exp(z) (15, 15) F64: well-conditioned smooth case.  Classical
    #          and SVD methods should produce rational functions that agree
    #          with exp(z) to ~eps at well-resolved points.  Classical does
    #          NOT diagonal-hop, so it returns the full (15, 15); SVD reduces
    #          (per existing 2.1.2 logic, exp(z) Padé(15, 15) at F64 hits
    #          its tol-floor early and reduces to a smaller block).
    # -------------------------------------------------------------------------
    @testset "CP.1.2: exp(z) Padé(15,15) classical retains full degree" begin
        c = Float64.([1 / factorial(big(k)) for k = 0:30])

        P_cls = classical_pade_diagonal(c, 15)
        # Classical doesn't reduce; full (15, 15) preserved.
        @test P_cls.μ == 15
        @test P_cls.ν == 15

        # Function value at a sweep of well-resolved points.
        for z in (-0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5)
            r = _eval_pade(P_cls, z)
            @test abs(r - exp(z)) < 1e-14
        end
    end

    # -------------------------------------------------------------------------
    # CP.1.3 — log(1.2 - z) (10, 10): branch-cut function.  Classical does
    #          not perform the GGT-style block reduction (the SVD-side
    #          (20→10) reduction at existing 2.1.3 is SVD-specific).  We
    #          assert classical produces a rational function value matching
    #          log(1.2 - z) at well-resolved interior points.
    # -------------------------------------------------------------------------
    @testset "CP.1.3: log(1.2 - z) Padé(10,10) classical function-value" begin
        c = zeros(Float64, 21)
        c[1] = log(1.2)
        for k = 1:20
            c[k+1] = -1.0 / (k * 1.2^k)
        end
        P_cls = classical_pade_diagonal(c, 10)
        @test P_cls.μ == 10
        @test P_cls.ν == 10
        @test abs(_eval_pade(P_cls, 0.0) - log(1.2)) < 1e-15
        @test abs(_eval_pade(P_cls, 0.3) - log(0.9)) < 1e-14
        @test abs(_eval_pade(P_cls, 0.6) - log(0.6)) < 1e-13
    end

    # -------------------------------------------------------------------------
    # CP.1.4 — Singular Toeplitz fallback.  GGT 2013 §7 first equation:
    #          f(z) = 1 + z² has Padé(1, 1) = 1 (defect 1).  The 1×1
    #          Toeplitz built from c = [1, 0, 1] is exactly [0] — singular.
    #          Classical must detect singularity and fall back to the SVD
    #          path, which produces the (0, 0) collapse with a = b = [1].
    #          This is per FW 2011 line 346 (singular fallback) + bead
    #          `padetaylor-txg` acceptance (b) (SVD route on singular).
    # -------------------------------------------------------------------------
    @testset "CP.1.4: singular Toeplitz routes to SVD fallback" begin
        c = Float64[1, 0, 1]
        # Direct call to classical_pade_diagonal must signal singular by
        # throwing.  We use `SingularException` because that's what `\`
        # would have thrown on a singular `lu` factorisation.
        @test_throws SingularException classical_pade_diagonal(c, 1)

        # robust_pade(c, 1, 1; method=:classical) must catch the singular
        # exception and fall back to SVD, producing the (0, 0) collapse.
        P_dispatch = robust_pade(c, 1, 1; method = :classical)
        @test P_dispatch.μ == 0
        @test P_dispatch.ν == 0
        @test P_dispatch.a ≈ [1.0]  atol = 1e-14
        @test P_dispatch.b ≈ [1.0]  atol = 1e-14

        # method=:svd directly also produces the same collapse.
        P_svd = robust_pade(c, 1, 1; method = :svd)
        @test P_svd.μ == 0
        @test P_svd.ν == 0
    end

    # -------------------------------------------------------------------------
    # CP.1.5 — Off-diagonal (m ≠ n) routes through SVD.  Classical's
    #          eqs. (5.4)+(5.5) generalise but FW 2011 §5.2 focuses on the
    #          diagonal (v, v) case; our v1 ships diagonal-only.  Callers
    #          with off-diagonal (m, n) automatically fall back to SVD.
    # -------------------------------------------------------------------------
    @testset "CP.1.5: off-diagonal (m≠n) routes through SVD" begin
        c = Float64.([1 / factorial(big(k)) for k = 0:10])
        # (3, 2) Padé of exp: well-defined; classical_pade_diagonal NOT
        # called (we only do (m, m)); robust_pade with method=:classical
        # must transparently route to SVD.
        P_cls = robust_pade(c, 3, 2; method = :classical)
        P_svd = robust_pade(c, 3, 2; method = :svd)
        @test P_cls.μ == P_svd.μ
        @test P_cls.ν == P_svd.ν
        # Function values agree.
        for z in (-0.3, 0.0, 0.3, 0.5)
            @test _eval_pade(P_cls, z) ≈ _eval_pade(P_svd, z) atol = 1e-14
        end
    end

    # -------------------------------------------------------------------------
    # CP.1.6 — Element-type dispatch.  Default for F64 / F32 / Complex
    #          variants is :classical; default for BigFloat / Arb /
    #          generic AbstractFloat is :svd (where GGT's relative-
    #          accuracy guarantees are load-bearing).
    #
    #          Distinguishing observable: exp(z) Padé(20, 20) at F64
    #          reduces under SVD to (7, 7) (existing 2.1.2 ground truth),
    #          but stays at full (20, 20) under classical (no diagonal
    #          hopping).  We use this to fingerprint which path ran.
    # -------------------------------------------------------------------------
    @testset "CP.1.6: element-type dispatch defaults" begin
        c_f64 = Float64.([1 / factorial(big(k)) for k = 0:40])

        # F64 default: classical → (20, 20) preserved (no diagonal hop).
        P_f64_default = robust_pade(c_f64, 20, 20)
        @test P_f64_default.μ == 20
        @test P_f64_default.ν == 20

        # F64 with method=:svd: (20, 20) → (7, 7) reduction (existing 2.1.2).
        P_f64_svd = robust_pade(c_f64, 20, 20; method = :svd)
        @test P_f64_svd.μ == 7
        @test P_f64_svd.ν == 7

        # Complex{Float64} default: classical (no reduction).
        c_cf64 = Complex{Float64}.(c_f64)
        P_cf64_default = robust_pade(c_cf64, 20, 20)
        @test P_cf64_default.μ == 20
        @test P_cf64_default.ν == 20

        # BigFloat default: :svd path.  At default tol (~2^-246), exp's
        # high-order coefficients stay above threshold, so no reduction
        # (existing 2.1.8).  Force tol=1e-14 to trigger F64-like reduction;
        # under :svd → (7, 7); under :classical override → (20, 20).
        setprecision(BigFloat, 256) do
            c_bf = [BigFloat(1) / factorial(big(k)) for k = 0:40]
            P_bf_default = robust_pade(c_bf, 20, 20; tol = 1e-14)
            @test P_bf_default.μ == 7
            @test P_bf_default.ν == 7
            P_bf_cls = robust_pade(c_bf, 20, 20; method = :classical)
            @test P_bf_cls.μ == 20
            @test P_bf_cls.ν == 20
        end
    end

    # -------------------------------------------------------------------------
    # CP.1.7 — Degenerate underlying rational forces a singular Toeplitz.
    #
    # 1/(1 - z/2) has c_k = (1/2)^k.  The m × m Toeplitz built from these
    # coefficients has T[i, j] = α^(m+i-j) (with α = 1/2): every row is
    # α^(i-1) · row 1, so rank(T) = 1 → exactly singular.  Classical must
    # detect and throw SingularException; `robust_pade(..; method =
    # :classical)` then routes to :svd which recovers the true (0, 1)
    # rational `1 / (1 - z/2)`.
    #
    # This complements CP.1.4 (defect-1 collapse at low m): here the
    # singular detection fires at high m on a smoother underlying
    # function — confirming the fallback is robust to the actual rank-
    # deficiency mechanism, not just the special case `c = [1, 0, 1]`.
    # -------------------------------------------------------------------------
    @testset "CP.1.7: geometric series triggers high-m singular fallback" begin
        c = Float64[(0.5)^k for k = 0:30]

        # Direct classical call: rank-1 Toeplitz → SingularException.
        @test_throws SingularException classical_pade_diagonal(c, 15)

        # robust_pade dispatch: catches, falls through to :svd, recovers
        # the underlying (0, 1) rational.  Function value matches 1/(1 - z/2)
        # to ~eps (SVD path's Padé(15,15) of 1/(1-z/2) IS exact: P.μ = 0,
        # P.ν = 1, a = [1], b = [1, -1/2]).
        P_dispatch = robust_pade(c, 15, 15; method = :classical)
        @test P_dispatch.μ == 0
        @test P_dispatch.ν == 1
        for z in (-0.4, -0.2, 0.0, 0.2, 0.4)
            @test abs(_eval_pade(P_dispatch, z) - 1 / (1 - z/2)) < 1e-14
        end
    end

    # -------------------------------------------------------------------------
    # CP.1.8 — Argument validation.  Classical entry point fails fast on
    #          negative m (CLAUDE.md Rule 1).
    # -------------------------------------------------------------------------
    @testset "CP.1.8: fail-fast on bad arguments" begin
        @test_throws ArgumentError classical_pade_diagonal(Float64[1, 0, 1], -1)
        @test_throws ArgumentError robust_pade(Float64[1, 0, 1], 1, 1;
                                                method = :bogus)
    end

    # -------------------------------------------------------------------------
    # CP.1.9 — Complex{Float64} sanity: classical works for complex
    #          coefficients (the path-network use case is Complex{F64}).
    # -------------------------------------------------------------------------
    @testset "CP.1.9: Complex{Float64} classical Padé" begin
        # exp(iz) = cos(z) + i sin(z).  Padé(5, 5) of exp at imaginary
        # coefficients.  Use real coefficients with complex eltype.
        c = Complex{Float64}[Complex{Float64}(1 / factorial(big(k))) for k = 0:10]
        P_cls = classical_pade_diagonal(c, 5)
        @test eltype(P_cls.a) === Complex{Float64}
        @test eltype(P_cls.b) === Complex{Float64}
        @test P_cls.μ == 5
        @test P_cls.ν == 5
        # exp(0.3) closed-form.
        @test abs(_eval_pade(P_cls, Complex{Float64}(0.3)) - exp(0.3)) < 1e-14
    end

end # @testset Classical Padé

# CP.1.M — Mutation-proof procedure.  Verified 2026-05-13 in bead
# `padetaylor-txg`; all four mutations bit then restored cleanly to
# 1377/1377 GREEN.  Full log in `docs/worklog/021-classical-pade-default.md`.
#
#   Mutation P1  --  in `classical_pade_diagonal`, flip the sign of the
#     RHS in the Toeplitz solve (`rhs[i] = cv[m + i + 1]` instead of `-`).
#     Verified bite: CP.1.1 (coefficient mismatch, dispatch-equiv check,
#     5 `r_cls ≈ r_svd` checks), CP.1.3 (2 log function-value asserts),
#     CP.1.9 (Complex{F64} exp check).  PROPAGATES through the entire
#     IVP stack: Phase 5 PadeStepper (~8 fails), Phase 6 Problems (9 fails
#     in 6.1.1-6.1.5), PathNetwork (PN.1.2 + PN.2.2 + PN.2.3 errors),
#     Phase 9 tritronquée (2 fails), LatticeDispatcher (51 fails).
#     ~60+ total RED — classical's correctness is load-bearing across
#     every downstream consumer.
#
#   Mutation P2  --  in `robust_pade(... ; method = :classical)`, remove
#     the singular-fallback try/catch around `classical_pade_diagonal`.
#     Verified bite: CP.1.4 (1f + 1e), CP.1.7 (1f + 1e), RobustPade 2.1.5
#     (1+z² defect-1 collapse — 1 error: F64 default now hits singular
#     classical instead of SVD reducing), Dispatcher DP.2.1 + DP.3.1
#     (1f + 2e in junction tests), LatticeDispatcher LD.2.1 (1 fail-fast
#     guard returns SingularException instead of ArgumentError).
#     Singular fallback is load-bearing.
#
#   Mutation P3  --  in `_default_pade_method`, swap F64 / Complex{F64}'s
#     :classical default for :svd.
#     Verified bite: CP.1.6 (4 of 10 — F64 default and Complex{F64} default
#     P.μ == 20 assertions fail because SVD reduces to (7, 7)); PN.2.2
#     (2 fails — rel-err 6.6e-12 vs new rtol 1e-12, |imag| above bound);
#     PN.2.3 (2 fails — rel-err 6.05e-6 vs new rtol 5e-10).  Total 8 RED.
#     Proves BOTH the dispatch-default change AND the PN.2.2 / PN.2.3
#     rtol tightening are load-bearing on classical being the default.
#
#   Mutation P4  --  in `robust_pade(... ; method = :classical)` dispatcher,
#     remove the `m == n` guard so off-diagonal also routes to
#     classical_pade_diagonal.
#     Verified bite: CP.1.5 (4 of 6 — `P_cls.μ == P_svd.μ` fails because
#     classical_pade_diagonal returns a (3, 3) approximant instead of the
#     (3, 2) shape the (m, n) request specified; function values at the
#     four z samples disagree).  Off-diagonal routing is load-bearing.
#
# Restoration: all mutations restored before commit; 1377/1377 GREEN
# post-restore.
