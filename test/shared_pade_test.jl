# =============================================================================
# shared_pade_test.jl — Phase V1a tests for PadeTaylor.SharedPade.
#
# Bead padetaylor-0ln.2.  TDD shape = port-and-verify (CLAUDE.md Rule 4):
# the d=1 case is the verbatim scalar GGT 2013 `:svd` path, so the primary
# oracle is `robust_pade(jet, m, m; method=:svd)` itself.  The d=2 case is
# cross-validated against a hand-constructed shared denominator whose roots
# are known in closed form.
#
# Each test asserts an invariant against a known-correct value (Rule 5):
#   - SP.1.1  d=1 reduction:  shared_denominator_pade([jet], m) == scalar :svd
#   - SP.1.2  genuine d=2:    recovered Q roots == known poles; r_i values
#   - SP.1.3/4/5             the failure-mode throws (Rule 1)
#   - SP.1.6  mutation-proof: documented impl perturbations turn SP.1.1/1.2 RED
#
# Self-contained: `using Test, PadeTaylor` so it runs standalone and under
# runtests.jl.
# =============================================================================

using Test
using LinearAlgebra
using PadeTaylor: shared_denominator_pade, robust_pade

# --- helpers -----------------------------------------------------------------

# Taylor coefficients [c_0 … c_N] of P(z)/Q(z) given numerator/denominator
# coefficient vectors (low-to-high) by formal power-series long division.
# c_k = (p_k - Σ_{j=1..k} q_j c_{k-j}) / q_0.
function _ratio_jet(P::Vector{Float64}, Q::Vector{Float64}, N::Int)
    c = zeros(Float64, N + 1)
    for k = 0:N
        s = k + 1 ≤ length(P) ? P[k + 1] : 0.0
        for j = 1:k
            qj = j + 1 ≤ length(Q) ? Q[j + 1] : 0.0
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / Q[1]
    end
    return c
end

# Evaluate a polynomial (low-to-high coeffs) at z.
_polyval(coef, z) = sum(coef[k] * z^(k - 1) for k = 1:length(coef))
_ratval(num, den, z) = _polyval(num, z) / _polyval(den, z)

@testset "SharedPade.shared_denominator_pade" begin

    # -------------------------------------------------------------------------
    # SP.1.1 — d=1 reduction to the scalar :svd path.
    #
    # f(z) = (1 + 2z) / (1 - 0.5z + 0.3z²) is an exact rational of type (1,2);
    # request denominator degree m=2.  shared_denominator_pade with a single
    # jet must reproduce robust_pade(jet, 2, 2; method=:svd) — same algorithm,
    # d=1.  Hard-won lesson (bead text): compare rational-function VALUES at
    # several z, not just raw coefficients, because near-tol blocks have
    # sign/scale-ambiguous coefficient vectors.
    # -------------------------------------------------------------------------
    @testset "SP.1.1 d=1 reduces to scalar :svd" begin
        Pnum = [1.0, 2.0]
        Qden = [1.0, -0.5, 0.3]
        m = 2
        jet = _ratio_jet(Pnum, Qden, 2m)          # length 2m+1
        nums, den = shared_denominator_pade([jet], m)
        scalar = robust_pade(jet, m, m; method = :svd)

        # One numerator only for d=1.
        @test length(nums) == 1
        # Denominator normalised to Q(0)=1 (b[1]=1), matching scalar path.
        @test den[1] ≈ 1.0 atol = 1e-13
        @test scalar.b[1] ≈ 1.0 atol = 1e-13
        # Coefficient agreement (both b[1]=1-normalised, so sign is fixed).
        @test den ≈ scalar.b atol = 1e-12
        @test nums[1] ≈ scalar.a atol = 1e-12
        # Function-value agreement: the load-bearing invariant.
        for z in (0.1, -0.2, 0.37, -0.41im, 0.25 + 0.1im)
            r_shared = _ratval(nums[1], den, z)
            r_scalar = _ratval(scalar.a, scalar.b, z)
            r_true   = _ratval(Pnum, Qden, z)
            @test r_shared ≈ r_scalar atol = 1e-12
            @test r_shared ≈ r_true   atol = 1e-12
        end
    end

    # -------------------------------------------------------------------------
    # SP.1.2 — genuine d=2 shared denominator.
    #
    # Known Q(z) = 1 - 0.4z + 0.2z² with two distinct (complex) roots; two
    # numerators P_1 = 1 + 0.3z, P_2 = 2 - 0.7z + 0.5z².  f_i = P_i/Q.
    # Taylor jets fed in; recovered shared denominator's ROOTS must equal the
    # roots of Q (the shared poles) and each f_i value must be recovered.
    # -------------------------------------------------------------------------
    @testset "SP.1.2 genuine d=2 shared Q" begin
        Qden = [1.0, -0.4, 0.2]
        P1   = [1.0, 0.3]
        P2   = [2.0, -0.7, 0.5]
        m = 2
        jet1 = _ratio_jet(P1, Qden, 2m)
        jet2 = _ratio_jet(P2, Qden, 2m)
        nums, den = shared_denominator_pade([jet1, jet2], m)

        @test length(nums) == 2
        @test den[1] ≈ 1.0 atol = 1e-13

        # Shared poles: roots of recovered den == roots of Qden.
        # roots of a 0+1+2 coeff poly via the companion / quadratic formula.
        roots_of(q) = begin
            a, b, c = q[3], q[2], q[1]
            disc = sqrt(complex(b^2 - 4a * c))
            sort([(-b + disc) / (2a), (-b - disc) / (2a)], by = z -> (real(z), imag(z)))
        end
        @test roots_of(den) ≈ roots_of(Qden) atol = 1e-10

        # Function values of each component recovered.
        for z in (0.1, -0.2, 0.33, 0.2 + 0.15im)
            @test _ratval(nums[1], den, z) ≈ _ratval(P1, Qden, z) atol = 1e-12
            @test _ratval(nums[2], den, z) ≈ _ratval(P2, Qden, z) atol = 1e-12
        end
    end

    # -------------------------------------------------------------------------
    # SP.1.3 — jet too short throws (pillar A §7 failure mode 1).
    # length(jets[i]) < m+1 must throw with a suggestion/detail message.
    # -------------------------------------------------------------------------
    @testset "SP.1.3 jet too short throws" begin
        short = [1.0, 2.0]            # length 2, need ≥ m+1 = 4
        @test_throws Exception shared_denominator_pade([short], 3)
        # The error message must be informative, not a bare bounds error.
        err = try
            shared_denominator_pade([short], 3)
            nothing
        catch e
            e
        end
        @test err !== nothing
        @test occursin("jet", lowercase(sprint(showerror, err)))
    end

    # -------------------------------------------------------------------------
    # SP.1.4 — all-zero jets throw (pillar A §7 failure mode 2).
    # A jet of all zeros is indistinguishable from zero at tolerance: every
    # singular value falls below τ = tol·‖c‖, so ρ = 0.  Must throw.
    # -------------------------------------------------------------------------
    @testset "SP.1.4 all-zero jets throw" begin
        zjet = zeros(Float64, 7)
        @test_throws Exception shared_denominator_pade([zjet, zjet], 3)
        err = try
            shared_denominator_pade([zjet, zjet], 3)
            nothing
        catch e
            e
        end
        @test err !== nothing
        msg = lowercase(sprint(showerror, err))
        @test occursin("zero", msg) || occursin("tol", msg)
    end

    # -------------------------------------------------------------------------
    # SP.1.5 — Q(0) ≈ 0 throws (pillar A §7 failure mode 3).
    #
    # Q(0)=0 means a pole at the expansion centre.  We construct this directly:
    # a jet whose shared denominator has a vanishing constant term is produced
    # when f_i = P_i / (z·something), i.e. the function itself blows up at 0.
    # Take f(z) = 1/z + analytic: not a power series, so instead we force the
    # condition by constructing a rational with Q(0)=0 in the implementation's
    # null-vector frame.  Cleanest honest construction: f(z) = (1+z)/(z - z²),
    # i.e. Q = [0, 1, -1].  The Taylor expansion of f about 0 does NOT exist
    # (pole at 0), so we cannot feed a finite jet.  Per Rule 9 we instead test
    # the throw via the post-normalisation guard: a denominator whose recovered
    # b[1] is below tol.  We trigger it with a jet of a function regular at 0
    # but where m is chosen so large the null vector's leading coeff underflows
    # is NOT reliable — so this is an honest non-test of the exact path.
    #
    # Honest non-test (Rule 9): the Q(0)≈0 branch is exercised structurally —
    # we assert the guard EXISTS and fires when handed a constructed null
    # vector with b[1]=0 via the d=1 jet of (z)/(1) which has c_0=0… that
    # still yields b[1]=1.  A clean Q(0)=0 trigger needs a pole at 0, which
    # has no Taylor jet.  We therefore assert the guard's presence indirectly:
    # the implementation must `throw` rather than divide-by-zero.  This is
    # covered by code review of the `_finalise` guard; documented here so a
    # future reader knows SP.1.5 is deliberately a structural check.
    # -------------------------------------------------------------------------
    @testset "SP.1.5 Q(0)≈0 guard (structural, Rule 9 honest non-test)" begin
        # A function regular at 0 cannot produce Q(0)=0 from a finite jet;
        # see the comment block above.  We verify the well-posed path does
        # NOT spuriously trip the guard for a normal input.
        Qden = [1.0, -0.4, 0.2]
        jet  = _ratio_jet([1.0, 0.3], Qden, 4)
        nums, den = shared_denominator_pade([jet], 2)
        @test abs(den[1]) > 1e-8         # guard correctly did NOT fire
    end

    # -------------------------------------------------------------------------
    # SP.1.6 — mutation-proof (Rule 4 / Rule 5).
    #
    # Three impl mutations were applied by hand; the record below states for
    # each whether it bit, with the genuine finding documented (Rule 2/3):
    #
    #   Mutation A (NON-biting) — "wrong SVD null vector": seed the QR
    #   reweighting with `b_init = Vt[end-1, :]` (second-smallest right
    #   singular vector) instead of `Vt[end, :]`.  Result: all 34 asserts
    #   stayed GREEN.  FINDING: `b_init` only sets the reweighting diagonal
    #   `D = diag(|b_init| + √ε)`; the QR step then recomputes the true null
    #   vector of `(A_full·D)ᵀ` independently.  With the `+√ε` floor `D` is
    #   nearly uniform, so the QR result is insensitive to a wrong `b_init`.
    #   This is not a test gap — `b_init` is genuinely not load-bearing for
    #   the result, only for conditioning.
    #
    #   Mutation A′ (BITING — load-bearing) — "wrong null-space column":
    #   select `b = D * F.Q[:, m_cur]` instead of `F.Q[:, m_cur+1]`, i.e.
    #   pick a non-null column of the QR factor.  Result: 22 of 34 asserts
    #   went RED — SP.1.1 (12 fail: coefficient AND value asserts), SP.1.2
    #   (9 fail: pole-root AND value asserts), SP.1.6 canary (1 fail).  The
    #   failure-mode throw tests SP.1.3/1.4/1.5 correctly stayed GREEN (they
    #   never reach the QR step).  Reverted after confirmation.
    #
    #   Mutation B (NON-biting, by design) — "negated denominator": set
    #   `b = -b` before numerator recovery.  The final `b[1]=1` normalisation
    #   divides BOTH a and b by b[1], so the ratio P/Q is invariant under a
    #   global sign.  Documented as a deliberately non-biting mutation: it
    #   proves the normalisation is sign-robust, not that a test regresses.
    #
    # The biting mutation (A′) is the load-bearing one.  This testset re-runs
    # the SP.1.1 value invariant as a standalone canary so the mutation-proof
    # record is self-contained.
    # -------------------------------------------------------------------------
    @testset "SP.1.6 mutation-proof canary" begin
        Pnum = [1.0, 2.0]
        Qden = [1.0, -0.5, 0.3]
        jet  = _ratio_jet(Pnum, Qden, 4)
        nums, den = shared_denominator_pade([jet], 2)
        # Under mutation A this value invariant goes RED (recovered ratio is
        # the second null vector, not P/Q).
        @test _ratval(nums[1], den, 0.3) ≈ _ratval(Pnum, Qden, 0.3) atol = 1e-12
    end

end
