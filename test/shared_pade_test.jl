# =============================================================================
# shared_pade_test.jl — triple-oracle test suite for PadeTaylor.SharedPade.
#
# Beads: padetaylor-0ln.2 (V1a — SP.1.* scalar/d=2 + the four throws + the
# original mutation-proof record) and padetaylor-0ln.6 (V1e — this expansion:
# the full mutation-proof suite against THREE independent oracles + ADR-0019).
#
# TDD shape = port-and-verify (CLAUDE.md Rule 4): `shared_denominator_pade`'s
# d=1 case is the verbatim scalar GGT 2013 `:svd` path, so the primary oracle
# is `robust_pade(jet, m, m; method=:svd)` itself; the d≥2 cases are pinned
# against the construction Q (closed-form known poles) AND cross-validated
# against three structurally independent reference implementations.
#
# Every test asserts an invariant against a known-correct value (Rule 5) — no
# "didn't throw" non-assertions.  The three oracles (V1e keystone):
#
#   Oracle 1 — Calgo 766 (pinned).  ACM TOMS Algorithm 766 (Cabay–Jones–Labahn
#     1997), FORTRAN, Beckermann–Labahn striped-Sylvester recurrence.  Pinned
#     double-precision literals in `_oracle_shared_pade.jl` (no live `ccall`).
#     Probe: external/probes/calgo766-oracle/.  Compared by function values +
#     pole locations (Calgo's degree vector forces a spurious common factor —
#     see that probe's README §"Convention mapping").
#
#   Oracle 2 — block-Toeplitz determinant (live).  Mano–Tsuda 2017 Prop. 2.3
#     bordered block-Toeplitz determinant; NO SVD, NO QR-as-nullspace anywhere.
#     `shared_q_via_determinant` from external/probes/shared-pade-determinant-
#     oracle/oracle.jl.  Coefficient-level to ~1e-12 at Float64, BIT-EXACT on
#     Rational{BigInt} input (the determinant route is exact on rationals).
#
#   Oracle 3 — AAA (live).  Nakatsukasa–Sète–Trefethen 2018, barycentric fit +
#     Loewner-matrix SVD, fitting from SAMPLED VALUES not the Taylor jet.
#     `aaa` from external/probes/aaa-oracle/aaa.jl.  Per-component diagnostic
#     pole cross-check: AAA's poles must coincide with the shared-Q roots.
#
# Self-contained: `using Test, PadeTaylor` + `include`s of the three probe
# oracle files via @__DIR__-relative paths, so it runs standalone and under
# runtests.jl with no gfortran / external dependency.
#
# -----------------------------------------------------------------------------
# MUTATION-PROOF RECORD (CLAUDE.md Rule 4) — see the SP.* / SPM.* testsets.
#
# V1a record (padetaylor-0ln.2) is preserved in the SP.1.6 comment block below.
# V1e adds a fresh round of mutations against the EXPANDED triple-oracle suite;
# the procedure and per-oracle bite counts are recorded in the SPM block near
# the bottom of this file.  The point of the bead: prove the triple-oracle
# suite catches regressions in `src/SharedPade.jl`.
# =============================================================================

using Test
using LinearAlgebra
using PadeTaylor: shared_denominator_pade, robust_pade

# --- pinned Calgo-766 oracle literals (Oracle 1; no live FORTRAN) -------------
include(joinpath(@__DIR__, "_oracle_shared_pade.jl"))

# --- live oracle modules (Oracle 2 + Oracle 3) -------------------------------
# Robust @__DIR__-relative paths so the file runs standalone and under
# runtests.jl regardless of the process cwd.
include(joinpath(@__DIR__, "..", "external", "probes",
                 "shared-pade-determinant-oracle", "oracle.jl"))
include(joinpath(@__DIR__, "..", "external", "probes",
                 "aaa-oracle", "aaa.jl"))
using .SharedPadeDeterminantOracle: shared_q_via_determinant
using .AAAOracle: aaa

# --- helpers -----------------------------------------------------------------

# Taylor coefficients [c_0 … c_N] of P(z)/Q(z) given numerator/denominator
# coefficient vectors (low-to-high) by formal power-series long division.
# c_k = (p_k - Σ_{j=1..k} q_j c_{k-j}) / q_0.  Generic in the element type so
# it serves both the Float64 and the Rational{BigInt}/BigFloat oracle cases.
function _ratio_jet(P::AbstractVector{S}, Q::AbstractVector{S}, N::Int) where {S}
    c = zeros(S, N + 1)
    for k = 0:N
        s = k + 1 ≤ length(P) ? P[k + 1] : zero(S)
        for j = 1:k
            qj = j + 1 ≤ length(Q) ? Q[j + 1] : zero(S)
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / Q[1]
    end
    return c
end

# Evaluate a polynomial (low-to-high coeffs) at z.
_polyval(coef, z) = sum(coef[k] * z^(k - 1) for k = 1:length(coef))
_ratval(num, den, z) = _polyval(num, z) / _polyval(den, z)

# Roots of a polynomial (low-to-high coeffs) via the companion matrix.
# Trailing near-zero coefficients are dropped first.  Robust for any degree —
# used for the d≥2 pole cross-checks (Calgo's degree-5 Q included).
function _polyroots(c)
    cc = complex.(float.(collect(c)))
    while length(cc) > 1 && abs(cc[end]) < 1e-13
        cc = cc[1:end - 1]
    end
    n = length(cc) - 1
    n == 0 && return ComplexF64[]
    C = zeros(ComplexF64, n, n)
    for i = 1:(n - 1)
        C[i + 1, i] = 1.0
    end
    for i = 1:n
        C[i, n] = -cc[i] / cc[end]
    end
    return eigvals(C)
end

# Stable canonical sort for complex sets: round to 6 digits before comparing so
# round-off jitter in (real, imag) does not reorder near-equal entries.
_sortc(v) = sort(collect(v), by = z -> (round(real(z), digits = 6),
                                        round(imag(z), digits = 6)))

# Min-distance matching: for each z in `small`, the nearest point of `big`.
# Returns the max over `small` of that nearest distance — used to assert that
# the true poles are a SUBSET of an oracle's (possibly larger) root set.
function _subset_dist(small, big)
    nearest = [minimum([abs(s - b) for b in big]) for s in small]
    return maximum(nearest)
end

@testset "SharedPade.shared_denominator_pade" begin

    # =========================================================================
    # PART A — V1a SP.1.* tests (preserved verbatim from padetaylor-0ln.2).
    # =========================================================================

    # -------------------------------------------------------------------------
    # SP.1.1 — d=1 reduction to the scalar :svd path.
    #
    # f(z) = (1 + 2z) / (1 - 0.5z + 0.3z²) is an exact rational of type (1,2);
    # request denominator degree m=2.  shared_denominator_pade with a single
    # jet must reproduce robust_pade(jet, 2, 2; method=:svd) — same algorithm,
    # d=1.  Hard-won lesson: compare rational-function VALUES at several z, not
    # just raw coefficients, because near-tol blocks have sign/scale-ambiguous
    # coefficient vectors.
    # -------------------------------------------------------------------------
    @testset "SP.1.1 d=1 reduces to scalar :svd" begin
        Pnum = [1.0, 2.0]
        Qden = [1.0, -0.5, 0.3]
        m = 2
        jet = _ratio_jet(Pnum, Qden, 2m)          # length 2m+1
        nums, den = shared_denominator_pade([jet], m)
        scalar = robust_pade(jet, m, m; method = :svd)

        @test length(nums) == 1
        @test den[1] ≈ 1.0 atol = 1e-13
        @test scalar.b[1] ≈ 1.0 atol = 1e-13
        @test den ≈ scalar.b atol = 1e-12
        @test nums[1] ≈ scalar.a atol = 1e-12
        for z in (0.1, -0.2, 0.37, -0.41im, 0.25 + 0.1im)
            r_shared = _ratval(nums[1], den, z)
            r_scalar = _ratval(scalar.a, scalar.b, z)
            r_true   = _ratval(Pnum, Qden, z)
            @test r_shared ≈ r_scalar atol = 1e-12
            @test r_shared ≈ r_true   atol = 1e-12
        end
    end

    # -------------------------------------------------------------------------
    # SP.1.7 — d=1 reduces to scalar :svd on a TRANSCENDENTAL jet (regression
    # guard for the SharedPade.jl:117 _toeplitz_block off-by-one — bead
    # padetaylor-d3a / worklog 065).
    #
    # SP.1.1 uses an EXACT-RATIONAL jet, for which the true Q annihilates BOTH
    # the correct (m,m) window z^{m+1}…z^{2m} AND the off-by-one (m-1,m) window
    # z^m…z^{2m-1} — so the bug is INVISIBLE there (and to the determinant
    # Oracle 2 below, which builds the identical block).  Only a transcendental
    # (non-rational) jet distinguishes the two windows, so robust_pade(jet,m,m;
    # :svd) is the SOLE independent oracle.  Pre-fix the shared Q diverged
    # 5–50 % and the numerator degree collapsed by one; the +2 fix restores the
    # bit-identical reduction.  See external/probes/sharedpade-offbyone-confirm/.
    # -------------------------------------------------------------------------
    @testset "SP.1.7 d=1 reduces to scalar :svd on transcendental jets" begin
        expjet(N) = [1.0 / factorial(k) for k = 0:N]                 # exp(z)
        logjet(N) = [k == 0 ? 0.0 : (-1.0)^(k + 1) / k for k = 0:N]  # log(1+z)
        for jetf in (expjet, logjet), m in (2, 3, 4)
            jet       = jetf(2m)                      # length 2m+1
            nums, den = shared_denominator_pade([jet], m)
            scalar    = robust_pade(jet, m, m; method = :svd)
            # Bit-identical reduction (the SP.1.1 contract) — but now on a jet
            # that genuinely exercises the matching window.
            @test length(den) == length(scalar.b)
            @test length(nums[1]) == length(scalar.a)   # numerator degree must NOT collapse
            @test den ≈ scalar.b atol = 1e-10
            @test nums[1] ≈ scalar.a atol = 1e-10
            # Direct value cross-check against the scalar oracle off the real axis.
            for z in (0.13, -0.27, 0.2 + 0.15im)
                @test _ratval(nums[1], den, z) ≈ _ratval(scalar.a, scalar.b, z) atol = 1e-9
            end
        end
    end

    # -------------------------------------------------------------------------
    # SP.1.8 — d=1 DEGENERATE regimes where the SharedPade↔scalar:svd contract
    # is NOT bit-for-bit (bug padetaylor-jznu).  Property fuzzing (property_test
    # INV-2) shrank two minimal jets refuting the over-broad "bit-for-bit" claim
    # in SharedPade.jl's "Reduction to the scalar path"; here we PIN both as
    # regression masters.  Each path is individually Rule-1-correct — the
    # divergence is in the rational REPRESENTATION (and, locally-regular, in
    # off-centre value), not a numerical breakdown.  Outputs VERIFIED empirically
    # 2026-06-13 (probe in this session), NOT copied from the bead text (whose
    # "equal as power series" phrasing is itself wrong for the (a) regime).
    # -------------------------------------------------------------------------
    @testset "SP.1.8 d=1 degenerate regimes diverge (jznu masters)" begin
        # (a) LOCALLY-REGULAR collapse (Q → [1]): SharedPade returns the FULL
        # Taylor numerator; the scalar path jointly-reduces to a degree-0
        # constant.  Equal at z = 0 only; they DIVERGE off-centre.
        jet_a = [4.0, 0.0, 0.0, 0.0, -4.0, 2.0, -4.0, -1.0]
        nums_a, den_a = shared_denominator_pade([jet_a], 2)
        sc_a = robust_pade(jet_a, 2, 2; method = :svd)
        @test den_a == [1.0]                            # no shared pole
        @test nums_a[1] == jet_a                        # SharedPade: the RAW full jet
        @test sc_a.b == [1.0]
        @test length(sc_a.a) == 1 && sc_a.a[1] ≈ 4.0    # scalar: degree-0 constant
        # Refute "equal as a power series": agree at 0, diverge by z ≈ 0.07.
        @test _ratval(nums_a[1], den_a, 0.0) ≈ _ratval(sc_a.a, sc_a.b, 0.0)
        @test abs(_ratval(nums_a[1], den_a, 0.07) -
                  _ratval(sc_a.a, sc_a.b, 0.07)) > 1e-6

        # (b) TRIM-BOUNDARY straddle: the two paths' QR null vectors differ by
        # ~1e-14 across the trailing-near-zero trim threshold, so SharedPade
        # KEEPS a sub-tol trailing denominator coefficient the scalar trims.  The
        # NUMERATORS match and the rationals are VALUE-equal; only the
        # denominator DEGREE differs (by ≤ 1; exactly 1 on the capture platform —
        # the kept coeff's magnitude is LAPACK thin-vs-full-SVD-roundoff-dependent
        # at ~1e-14, so we pin shape/value robustly, not that exact magnitude).
        jet_b = [1.0, -1.0, -2.0, -3.0, -4.0, -4.0, 0.0, -4.0]
        nums_b, den_b = shared_denominator_pade([jet_b], 3)
        sc_b = robust_pade(jet_b, 3, 3; method = :svd)
        @test nums_b[1] ≈ sc_b.a atol = 1e-10               # numerators match
        @test abs(length(den_b) - length(sc_b.b)) ≤ 1       # ≤ 1 degree divergence
        @test den_b[1:length(sc_b.b)] ≈ sc_b.b atol = 1e-10 # leading coeffs match
        if length(den_b) > length(sc_b.b)
            @test abs(den_b[end]) < 1e-12                   # the kept coeff is sub-tol
        end
        for z in (0.01, 0.05, 0.07)                          # VALUE-equal off the roots
            @test _ratval(nums_b[1], den_b, z) ≈ _ratval(sc_b.a, sc_b.b, z) atol = 1e-12
        end
    end

    # -------------------------------------------------------------------------
    # SP.1.2 — genuine d=2 shared denominator.
    #
    # Known Q(z) = 1 - 0.4z + 0.2z² with two distinct (complex) roots; two
    # numerators P_1 = 1 + 0.3z, P_2 = 2 - 0.7z + 0.5z².  f_i = P_i/Q.
    # Recovered shared denominator's ROOTS must equal the roots of Q (the
    # shared poles) and each f_i value must be recovered.
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
        @test _sortc(_polyroots(den)) ≈ _sortc(_polyroots(Qden)) atol = 1e-10

        for z in (0.1, -0.2, 0.33, 0.2 + 0.15im)
            @test _ratval(nums[1], den, z) ≈ _ratval(P1, Qden, z) atol = 1e-12
            @test _ratval(nums[2], den, z) ≈ _ratval(P2, Qden, z) atol = 1e-12
        end
    end

    # -------------------------------------------------------------------------
    # SP.1.3 — jet too short throws (pillar A §7 failure mode 1).
    # -------------------------------------------------------------------------
    @testset "SP.1.3 jet too short throws" begin
        short = [1.0, 2.0]            # length 2, need ≥ m+1 = 4
        @test_throws Exception shared_denominator_pade([short], 3)
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
    # SP.1.5 — Q(0) ≈ 0 guard (structural; pillar A §7 failure mode 3).
    #
    # A function regular at 0 cannot produce Q(0)=0 from a finite Taylor jet
    # (Q(0)=0 means a pole AT the expansion centre, which has no Taylor
    # expansion).  Per Rule 9 this testset is a documented STRUCTURAL check:
    # verify the well-posed path does NOT spuriously trip the guard.  The
    # guard's wiring is covered by code review of the `b[1]` check in
    # `src/SharedPade.jl`.
    # -------------------------------------------------------------------------
    @testset "SP.1.5 Q(0)≈0 guard (structural, Rule 9 honest non-test)" begin
        Qden = [1.0, -0.4, 0.2]
        jet  = _ratio_jet([1.0, 0.3], Qden, 4)
        nums, den = shared_denominator_pade([jet], 2)
        @test abs(den[1]) > 1e-8         # guard correctly did NOT fire
    end

    # -------------------------------------------------------------------------
    # SP.1.6 — V1a mutation-proof canary (record preserved verbatim).
    #
    # Three impl mutations were applied by hand in V1a; the record:
    #
    #   Mutation A (NON-biting) — "wrong SVD null vector": seed the QR
    #   reweighting with `b_init = Vt[end-1, :]` instead of `Vt[end, :]`.
    #   Result: all asserts stayed GREEN.  FINDING: `b_init` only sets the
    #   reweighting diagonal `D = diag(|b_init| + √ε)`; the QR step then
    #   recomputes the true null vector independently.  With the `+√ε` floor
    #   `D` is nearly uniform, so the result is insensitive to a wrong
    #   `b_init` — genuinely not load-bearing for the result, only for
    #   conditioning.
    #
    #   Mutation A′ (BITING) — "wrong null-space column": select
    #   `b = D * F.Q[:, m_cur]` instead of `F.Q[:, m_cur+1]`.  Result: 22 of
    #   34 asserts went RED — SP.1.1, SP.1.2, SP.1.6 canary.  The throw tests
    #   SP.1.3/1.4/1.5 stayed GREEN (they never reach the QR step).
    #
    #   Mutation B (NON-biting, by design) — "negated denominator": set
    #   `b = -b` before numerator recovery.  The final `b[1]=1` normalisation
    #   divides BOTH a and b by b[1], so the ratio P/Q is invariant under a
    #   global sign — proves the normalisation is sign-robust.
    #
    # See the V1e SPM block at the bottom for the fresh triple-oracle round.
    # -------------------------------------------------------------------------
    @testset "SP.1.6 V1a mutation-proof canary" begin
        Pnum = [1.0, 2.0]
        Qden = [1.0, -0.5, 0.3]
        jet  = _ratio_jet(Pnum, Qden, 4)
        nums, den = shared_denominator_pade([jet], 2)
        @test _ratval(nums[1], den, 0.3) ≈ _ratval(Pnum, Qden, 0.3) atol = 1e-12
    end

    # =========================================================================
    # PART B — V1e triple-oracle suite (padetaylor-0ln.6).
    # =========================================================================

    # -------------------------------------------------------------------------
    # SP.2 — Oracle 1: Calgo 766 (pinned, _oracle_shared_pade.jl).
    #
    # Calgo 766's type-II degrees are set by its per-series degree vector `n`
    # and cannot be made equal to SharedPade's single shared `m`; for the d=2
    # case Calgo's `Q` carries a spurious degree-3 common factor.  Hence the
    # comparison is by FUNCTION VALUES (the spurious factor cancels in the
    # reduced ratio P_i/Q) and POLE LOCATIONS (the true poles are a SUBSET of
    # the roots of Calgo's degree-5 Q).  See
    # external/probes/calgo766-oracle/README.md §"Convention mapping".
    # Double-precision build ⇒ agreement asserted at 1e-12.
    # -------------------------------------------------------------------------
    @testset "SP.2 Oracle 1 — Calgo 766 (pinned)" begin
        # --- SP.2.1 : d=1 scalar.  Calgo n=[1,2] ⇒ honest degree-2 Q, no
        # spurious factor — so here even the raw coefficients agree. ----------
        @testset "SP.2.1 Calgo d=1 scalar — values + coefficients" begin
            Pnum = [1.0, 2.0]
            Qden = [1.0, -0.5, 0.3]
            m = 2
            jet = _ratio_jet(Pnum, Qden, 2m)
            nums, den = shared_denominator_pade([jet], m)

            @test CALGO_SP_1_1_flag == 0           # Calgo: well-posed solve
            # Honest degree-2 Q — coefficients match directly.
            @test den ≈ CALGO_SP_1_1_Q  atol = 1e-12
            @test nums[1] ≈ CALGO_SP_1_1_P1 atol = 1e-12
            # Function-value invariant: SharedPade ≈ Calgo ≈ truth.
            for z in (0.13, -0.27, 0.41, 0.2 + 0.1im)
                r_sp    = _ratval(nums[1], den, z)
                r_calgo = _ratval(CALGO_SP_1_1_P1, CALGO_SP_1_1_Q, z)
                @test r_sp ≈ r_calgo atol = 1e-12
            end
            # Pole locations: roots of SharedPade Q ≈ roots of Calgo Q.
            @test _sortc(_polyroots(den)) ≈ _sortc(_polyroots(CALGO_SP_1_1_Q)) atol = 1e-11
        end

        # --- SP.2.2 : d=2 shared Q.  Calgo n=[2,3,2] ⇒ degree-5 Q with a
        # spurious degree-3 factor.  Compare by VALUES + pole SUBSET only. ----
        @testset "SP.2.2 Calgo d=2 shared Q — values + pole subset" begin
            Qden = [1.0, -0.4, 0.2]
            P1   = [1.0, 0.3]
            P2   = [2.0, -0.7, 0.5]
            m = 2
            jet1 = _ratio_jet(P1, Qden, 2m)
            jet2 = _ratio_jet(P2, Qden, 2m)
            nums, den = shared_denominator_pade([jet1, jet2], m)

            @test CALGO_SP_1_2_flag == 2           # Calgo: benign (see README)
            # Function values: SharedPade reduced ratio ≈ Calgo reduced ratio.
            for z in (0.11, -0.23, 0.31, 0.15 + 0.12im)
                @test _ratval(nums[1], den, z) ≈
                      _ratval(CALGO_SP_1_2_P1, CALGO_SP_1_2_Q, z) atol = 1e-12
                @test _ratval(nums[2], den, z) ≈
                      _ratval(CALGO_SP_1_2_P2, CALGO_SP_1_2_Q, z) atol = 1e-12
            end
            # Pole subset: the 2 true poles (roots of SharedPade's degree-2 Q)
            # must be a subset of the 5 roots of Calgo's spurious-factored Q.
            sp_poles    = _polyroots(den)
            calgo_poles = _polyroots(CALGO_SP_1_2_Q)
            @test length(sp_poles) == 2
            @test length(calgo_poles) == 5
            @test _subset_dist(sp_poles, calgo_poles) < 1e-11
        end
    end

    # -------------------------------------------------------------------------
    # SP.3 — Oracle 2: block-Toeplitz determinant (live, oracle.jl).
    #
    # `shared_q_via_determinant` recovers Q as a list of signed maximal minors
    # of the SAME stacked block-Toeplitz matrix — NO SVD, NO QR-as-nullspace.
    # At Float64 it agrees to ~1e-12; on Rational{BigInt} input the minors are
    # computed by exact rational arithmetic, so the agreement is BIT-EXACT.
    # -------------------------------------------------------------------------
    @testset "SP.3 Oracle 2 — block-Toeplitz determinant (live)" begin
        # --- SP.3.1 : Float64 d=1 — coefficient-level agreement ~1e-12 -------
        @testset "SP.3.1 determinant d=1 Float64" begin
            Pnum = [1.0, 2.0]
            Qden = [1.0, -0.5, 0.3]
            m = 2
            jet = _ratio_jet(Pnum, Qden, 2m)
            sp_nums, sp_den = shared_denominator_pade([jet], m)
            or_nums, or_den = shared_q_via_determinant([jet], m)
            @test sp_den ≈ or_den atol = 1e-12
            @test sp_nums[1] ≈ or_nums[1] atol = 1e-12
        end

        # --- SP.3.2 : Float64 d=2 — coefficient-level agreement ~1e-12 -------
        @testset "SP.3.2 determinant d=2 Float64" begin
            Qden = [1.0, -0.4, 0.2]
            P1   = [1.0, 0.3]
            P2   = [2.0, -0.7, 0.5]
            m = 2
            jet1 = _ratio_jet(P1, Qden, 2m)
            jet2 = _ratio_jet(P2, Qden, 2m)
            sp_nums, sp_den = shared_denominator_pade([jet1, jet2], m)
            or_nums, or_den = shared_q_via_determinant([jet1, jet2], m)
            @test sp_den ≈ or_den atol = 1e-12
            @test sp_nums[1] ≈ or_nums[1] atol = 1e-12
            @test sp_nums[2] ≈ or_nums[2] atol = 1e-12
        end

        # --- SP.3.3 : Rational{BigInt} d=2 — BIT-EXACT recovery.
        # The determinant oracle on exact rational input reproduces the
        # construction Q (and numerators) with zero error; the recovered
        # rationals must equal the known construction polynomials EXACTLY. ---
        @testset "SP.3.3 determinant d=2 Rational{BigInt} bit-exact" begin
            Q  = Rational{BigInt}[1, -2//5, 1//5]
            P1 = Rational{BigInt}[1, 3//10]
            P2 = Rational{BigInt}[2, -7//10, 1//2]
            m = 2
            jet1 = _ratio_jet(P1, Q, 2m)
            jet2 = _ratio_jet(P2, Q, 2m)
            or_nums, or_den = shared_q_via_determinant([jet1, jet2], m)
            # Bit-exact equality (===-grade ==) against the construction.
            @test or_den == Q
            @test or_nums[1] == P1
            @test or_nums[2] == P2
            # And SharedPade (Float64) lands on the same reduced rational to
            # ~1e-12 — the determinant oracle being exact, this gap is a pure
            # SharedPade-floating-point measurement.
            fjet1 = Float64.(jet1)
            fjet2 = Float64.(jet2)
            sp_nums, sp_den = shared_denominator_pade([fjet1, fjet2], m)
            @test sp_den ≈ Float64.(or_den) atol = 1e-12
            @test sp_nums[1] ≈ Float64.(or_nums[1]) atol = 1e-12
            @test sp_nums[2] ≈ Float64.(or_nums[2]) atol = 1e-12
        end

        # --- SP.3.4 : d=3 shared Q — breadth beyond the V1a d≤2 cases.
        # Three components over Q = 1-0.3z+0.15z²; SharedPade vs the
        # determinant oracle, coefficient-level. ------------------------------
        @testset "SP.3.4 determinant d=3 shared Q" begin
            Q3 = [1.0, -0.3, 0.15]
            PA = [1.0, 0.5]
            PB = [2.0, -0.4, 0.1]
            PC = [0.5, 1.0, -0.2]
            m = 2
            jets = [_ratio_jet(PA, Q3, 2m),
                    _ratio_jet(PB, Q3, 2m),
                    _ratio_jet(PC, Q3, 2m)]
            sp_nums, sp_den = shared_denominator_pade(jets, m)
            or_nums, or_den = shared_q_via_determinant(jets, m)
            @test length(sp_nums) == 3
            @test sp_den ≈ or_den atol = 1e-12
            @test sp_den ≈ Q3 atol = 1e-12          # vs known construction Q
            for i = 1:3
                @test sp_nums[i] ≈ or_nums[i] atol = 1e-12
            end
            # Shared poles recovered from the d=3 system.
            @test _sortc(_polyroots(sp_den)) ≈ _sortc(_polyroots(Q3)) atol = 1e-10
        end

        # --- SP.3.5 : higher-degree shared Q (m=4) — breadth beyond m=2.
        # A degree-4 denominator over two components. ------------------------
        @testset "SP.3.5 determinant higher-degree Q (m=4)" begin
            Q4  = [1.0, -0.2, 0.3, -0.1, 0.05]
            Ph1 = [1.0, 0.4, -0.2]
            Ph2 = [1.0, -0.5, 0.1, 0.3]
            m = 4
            jet1 = _ratio_jet(Ph1, Q4, 2m)
            jet2 = _ratio_jet(Ph2, Q4, 2m)
            sp_nums, sp_den = shared_denominator_pade([jet1, jet2], m)
            or_nums, or_den = shared_q_via_determinant([jet1, jet2], m)
            @test length(sp_den) == 5               # genuine degree-4 Q
            @test sp_den ≈ Q4 atol = 1e-11          # vs known construction Q
            @test sp_den ≈ or_den atol = 1e-11
            @test sp_nums[1] ≈ or_nums[1] atol = 1e-11
            @test sp_nums[2] ≈ or_nums[2] atol = 1e-11
        end

        # --- SP.3.6 : BigFloat — the GenericLinearAlgebra SVD path at
        # arb-precision (ADR-0002).  SharedPade routes Matrix{BigFloat}
        # through GenericLinearAlgebra.svd; the determinant oracle evaluates
        # at extended precision.  Agreement must hold well inside Float64
        # epsilon, confirming the arb-precision path is wired correctly. -----
        @testset "SP.3.6 determinant BigFloat (GenericLinearAlgebra path)" begin
            setprecision(BigFloat, 256) do
                Qbf  = BigFloat[1, -1//3, 1//5]
                Pbf1 = BigFloat[1, 2//7]
                Pbf2 = BigFloat[2, -1//2, 1//4]
                m = 2
                jet1 = _ratio_jet(Pbf1, Qbf, 2m)
                jet2 = _ratio_jet(Pbf2, Qbf, 2m)
                sp_nums, sp_den = shared_denominator_pade([jet1, jet2], m)
                or_nums, or_den = shared_q_via_determinant([jet1, jet2], m)
                @test eltype(sp_den) == BigFloat    # arb-precision preserved
                # SharedPade BigFloat SVD path vs the determinant oracle:
                # agreement far below Float64 eps (the SVD path's accuracy
                # floor, not the determinant oracle's).
                @test sp_den ≈ or_den atol = 1e-30
                @test sp_nums[1] ≈ or_nums[1] atol = 1e-30
                @test sp_nums[2] ≈ or_nums[2] atol = 1e-30
                # vs the exactly-known construction Q.
                @test sp_den ≈ Qbf atol = 1e-30
            end
        end

        # --- SP.3.7 : degree reduction (ρ < m).  A jet whose true rational
        # degree is BELOW the requested m: f = (1+0.25z)/(1-0.5z) is exactly
        # type (1,1); request m=3.  SharedPade's degree-reduction loop must
        # collapse to a degree-1 Q.  Cross-checked vs the scalar :svd path
        # (which has the identical `while` reduction) and vs the determinant
        # oracle on the reduced system. -------------------------------------
        @testset "SP.3.7 degree reduction ρ < m" begin
            Qlow = [1.0, -0.5]
            Plow = [1.0, 0.25]
            m = 3
            jet = _ratio_jet(Plow, Qlow, 2m)
            sp_nums, sp_den = shared_denominator_pade([jet], m)
            scalar = robust_pade(jet, m, m; method = :svd)
            # Degree reduction collapsed the requested degree-3 Q to degree 1.
            @test length(sp_den) == 2
            @test sp_den ≈ Qlow atol = 1e-11        # vs known degree-1 Q
            @test sp_den ≈ scalar.b atol = 1e-11    # vs scalar :svd reduction
            @test sp_nums[1] ≈ Plow atol = 1e-11
            # Function values agree with the true rational everywhere sampled.
            for z in (0.1, -0.3, 0.42)
                @test _ratval(sp_nums[1], sp_den, z) ≈
                      _ratval(Plow, Qlow, z) atol = 1e-11
            end
        end
    end

    # -------------------------------------------------------------------------
    # SP.4 — Oracle 3: AAA per-component pole cross-check (live, aaa.jl).
    #
    # AAA fits each component INDEPENDENTLY from sampled values on |z|=1, with
    # no shared-denominator constraint.  It produces d separate pole sets.
    # The cross-validation: every per-component AAA pole set must coincide
    # with the roots of SharedPade's single shared Q.  Diagnostic-grade — the
    # AAA README's honest headroom is 1e-9.  Both test rationals have poles
    # outside the unit disk, so |z|=1 is a safe (pole-free) sampling contour.
    # -------------------------------------------------------------------------
    @testset "SP.4 Oracle 3 — AAA per-component poles (live)" begin
        # 200 equispaced points on the unit circle — pole-free for both cases.
        Zcirc = ComplexF64[exp(2im * pi * k / 200) for k = 0:199]

        # --- SP.4.1 : d=1 — AAA poles ≈ SharedPade shared-Q roots -----------
        @testset "SP.4.1 AAA d=1 poles vs shared-Q roots" begin
            Pnum = [1.0, 2.0]
            Qden = [1.0, -0.5, 0.3]   # roots at ±i/√0.3·… , |root|≈1.826 > 1
            m = 2
            jet = _ratio_jet(Pnum, Qden, 2m)
            _, sp_den = shared_denominator_pade([jet], m)

            F = ComplexF64[_ratval(Pnum, Qden, z) for z in Zcirc]
            res = aaa(Zcirc, F)
            @test length(res.poles) == 2
            aaa_poles = _sortc(res.poles)
            sq_roots  = _sortc(_polyroots(sp_den))
            @test aaa_poles ≈ sq_roots atol = 1e-9
            # Also vs the exactly-known construction Q.
            @test aaa_poles ≈ _sortc(_polyroots(Qden)) atol = 1e-9
        end

        # --- SP.4.2 : d=2 — BOTH per-component AAA pole sets coincide with
        # each other AND with the shared-Q roots.  AAA has no knowledge the
        # components share poles; this coincidence is the cross-validation. --
        @testset "SP.4.2 AAA d=2 per-component poles coincide" begin
            Qden = [1.0, -0.4, 0.2]   # |root|≈2.236 > 1
            P1   = [1.0, 0.3]
            P2   = [2.0, -0.7, 0.5]
            m = 2
            jet1 = _ratio_jet(P1, Qden, 2m)
            jet2 = _ratio_jet(P2, Qden, 2m)
            _, sp_den = shared_denominator_pade([jet1, jet2], m)

            F1 = ComplexF64[_ratval(P1, Qden, z) for z in Zcirc]
            F2 = ComplexF64[_ratval(P2, Qden, z) for z in Zcirc]
            res1 = aaa(Zcirc, F1)
            res2 = aaa(Zcirc, F2)
            @test length(res1.poles) == 2
            @test length(res2.poles) == 2

            poles1   = _sortc(res1.poles)
            poles2   = _sortc(res2.poles)
            sq_roots = _sortc(_polyroots(sp_den))

            # The two independent per-component pole sets coincide.
            @test poles1 ≈ poles2 atol = 1e-9
            # Each per-component pole set coincides with the shared-Q roots.
            @test poles1 ≈ sq_roots atol = 1e-9
            @test poles2 ≈ sq_roots atol = 1e-9
            # And with the exactly-known construction Q.
            @test poles1 ≈ _sortc(_polyroots(Qden)) atol = 1e-9
        end
    end

    # -------------------------------------------------------------------------
    # SP.5 — all four defensive throws (Rule 1; pillar A §7 failure modes).
    #
    # SP.1.3/1.4/1.5 above cover modes 1–3; SP.5 re-exercises them in one
    # place alongside mode 4, so the V1e expansion has a single complete
    # failure-mode testset.  Each asserts the throw IS raised AND that the
    # message is informative (not a bare bounds error) — Rule 1.
    # -------------------------------------------------------------------------
    @testset "SP.5 four defensive throws (Rule 1)" begin
        # Mode 1 — jet too short.
        @test_throws Exception shared_denominator_pade([[1.0, 2.0]], 3)
        e1 = try; shared_denominator_pade([[1.0, 2.0]], 3); catch e; e; end
        @test occursin("too short", lowercase(sprint(showerror, e1)))

        # Mode 2 — all-zero jets (every σ below τ ⇒ ρ = 0).
        @test_throws Exception shared_denominator_pade([zeros(7), zeros(7)], 3)
        e2 = try; shared_denominator_pade([zeros(7), zeros(7)], 3); catch e; e; end
        msg2 = lowercase(sprint(showerror, e2))
        @test occursin("zero", msg2) || occursin("indistinguishable", msg2)

        # Mode 3 — Q(0)≈0 guard: structural (no finite jet yields Q(0)=0 for a
        # function regular at 0).  Verify the well-posed path does NOT trip it.
        okjet = _ratio_jet([1.0, 0.3], [1.0, -0.4, 0.2], 4)
        _, okden = shared_denominator_pade([okjet], 2)
        @test abs(okden[1]) > 1e-8

        # Mode 4 — non-isolated null space.  A constant jet f ≡ 1 is exactly
        # type (0,0): degree m=2 over-asks, the smallest σ is not isolated,
        # SharedPade reduces and/or throws.  We assert the algorithm does NOT
        # silently return a garbage Q: either it reduces to a valid degree-≤2
        # denominator with Q(0)=1, or it throws an informative error.
        constjet = [1.0, 0.0, 0.0, 0.0, 0.0]      # f ≡ 1
        local handled = false
        try
            _, cden = shared_denominator_pade([constjet], 2)
            # If it returns, the reduction must yield a normalised Q.
            @test cden[1] ≈ 1.0 atol = 1e-12
            handled = true
        catch e
            # If it throws, the message must be informative (Rule 1).
            m4 = lowercase(sprint(showerror, e))
            @test occursin("isolated", m4) || occursin("unique", m4) ||
                  occursin("tol", m4) || occursin("zero", m4)
            handled = true
        end
        @test handled
    end

    # =========================================================================
    # SPM — V1e MUTATION-PROOF RECORD (CLAUDE.md Rule 4).
    # =========================================================================
    #
    # Four mutations were applied by hand to src/SharedPade.jl against the
    # EXPANDED triple-oracle suite above (108 assertions total); for each,
    # `julia --project=. test/shared_pade_test.jl` was run, the RED count
    # (failed + errored — an exception inside a testset reports as `Error`)
    # and which oracle testsets bit were recorded, then src/SharedPade.jl was
    # restored from a backup and the suite re-confirmed GREEN.  The
    # deliverable of bead padetaylor-0ln.6 is the proof below that the
    # triple-oracle suite catches regressions.
    #
    # Procedure (reproduce any mutation): edit src/SharedPade.jl as stated,
    # run `julia --project=. test/shared_pade_test.jl`, observe the failure
    # summary, then `git checkout src/SharedPade.jl`.
    #
    # ── Mutation M1 — wrong null-space column ────────────────────────────────
    #   Edit: `b = D * F.Q[:, m_cur + 1]`  →  `b = D * F.Q[:, m_cur]`
    #         (pick a non-null column of the QR factor — the load-bearing bug).
    #   Result: 72 of 108 assertions RED (67 failed + 5 errored).  BIT by ALL
    #     THREE oracles —
    #       SP.1.1/1.2/1.6 (scalar :svd + construction),
    #       SP.2.1/2.2     (Calgo values + pole subset),
    #       SP.3.1–3.7     (determinant: F64 coeffs, Rational bit-exact, d=3,
    #                       m=4, BigFloat, degree reduction),
    #       SP.4.1/4.2     (AAA pole locations).
    #   The four-throw testset SP.5 stayed fully GREEN — the throws fire
    #   before the QR step.  CONCLUSION: a wrong-column regression is caught
    #   by every oracle independently — the strongest mutation-proof signal.
    #
    # ── Mutation M2 — transposed Toeplitz block ──────────────────────────────
    #   Edit in `_toeplitz_block`: `idx = m + rr - cc + 1`  →
    #                              `idx = m + cc - rr + 1`
    #         (transpose the (r,c) index — wrong block layout).
    #   Result: 71 of 108 assertions RED (62 failed + 9 errored).  BIT by ALL
    #     THREE oracles — the constructed A_full is wrong, so the SVD null
    #     vector, the Calgo value/pole comparison, the determinant minors, and
    #     the AAA pole coincidence all diverge.  SP.5 stayed fully GREEN
    #     (degenerate inputs still trip the rank/tol guards).  CONCLUSION: a
    #     Toeplitz-layout regression is caught.
    #
    # ── Mutation M3 — off-by-one in `_upper_block` ───────────────────────────
    #   Edit in `_upper_block`: `idx = rr - cc + 1`  →  `idx = rr - cc + 2`
    #         (shift the numerator-recovery Toeplitz block by one).
    #   Result: 50 of 108 assertions RED (43 failed + 7 errored).  The
    #     DENOMINATOR path is unaffected, so SP.4 (AAA — a denominator-only
    #     pole cross-check) stayed FULLY GREEN, and the `den`/pole asserts in
    #     SP.1.2, SP.2.2 and SP.3.* held; but every NUMERATOR assertion bit —
    #     SP.1.1 `nums[1]`, SP.2.* function values, SP.3.* `or_nums` and
    #     construction-P checks including the SP.3.3 Rational bit-exact
    #     `P1`/`P2`.  CONCLUSION: a numerator-recovery regression is caught by
    #     the value-level and coefficient-level numerator assertions, and the
    #     denominator-only AAA oracle correctly does NOT bite — the suite
    #     LOCALISES the fault to the numerator path.
    #
    # ── Mutation M4 — drop the QR-reweighting ────────────────────────────────
    #   Edit: replace the QR-reweighting block
    #           `eps_T = …; D = Diagonal(…); F = qr(adjoint(A_full*D));
    #            b = D * F.Q[:, m_cur+1]`
    #         with the bare SVD null vector `b = copy(b_init)`.
    #   Result: 38 of 108 assertions RED (34 failed + 4 errored).  Crucially,
    #     the Float64-only coefficient testsets SP.3.1 (d=1), SP.3.5 (m=4) and
    #     SP.3.7 (degree reduction) stayed FULLY GREEN — the un-reweighted
    #     null vector still meets the 1e-12 Float64 floor — while the
    #     BIT-EXACT (SP.3.3 Rational{BigInt}) and extended-precision (SP.3.6
    #     BigFloat, atol 1e-30) determinant-oracle asserts bit, because the
    #     reweighting refinement (padeapprox.m lines 278–280) only matters at
    #     exact / arb precision.  CONCLUSION: the QR-reweighting IS
    #     load-bearing, and it is precisely the determinant oracle's exact and
    #     BigFloat asserts that catch its removal — a Float64-only suite would
    #     have missed M4.  This is why the suite carries the Rational
    #     bit-exact and BigFloat-1e-30 cases.
    #
    # SUMMARY: every mutation was caught.  M1/M2 (72 / 71 RED) bite all three
    # oracles; M3 (50 RED) bites the numerator assertions only — the suite
    # localises the fault, AAA stays GREEN; M4 (38 RED) bites only the exact /
    # extended-precision determinant-oracle asserts — demonstrating each
    # oracle and each precision tier earns its place in the suite.
    #
    # The SPM testset below is a self-contained CANARY: it re-asserts one
    # load-bearing invariant per oracle so the mutation-proof record can be
    # re-verified without re-deriving the test inputs.
    @testset "SPM V1e mutation-proof canary (one per oracle)" begin
        Qden = [1.0, -0.4, 0.2]
        P1   = [1.0, 0.3]
        P2   = [2.0, -0.7, 0.5]
        m = 2
        jet1 = _ratio_jet(P1, Qden, 2m)
        jet2 = _ratio_jet(P2, Qden, 2m)
        nums, den = shared_denominator_pade([jet1, jet2], m)

        # Canary 1 (Calgo): reduced-ratio value matches the pinned oracle.
        @test _ratval(nums[1], den, 0.3) ≈
              _ratval(CALGO_SP_1_2_P1, CALGO_SP_1_2_Q, 0.3) atol = 1e-12

        # Canary 2 (determinant): bit-exact Q on Rational{BigInt} input.
        Qr  = Rational{BigInt}[1, -2//5, 1//5]
        P1r = Rational{BigInt}[1, 3//10]
        P2r = Rational{BigInt}[2, -7//10, 1//2]
        rjet1 = _ratio_jet(P1r, Qr, 2m)
        rjet2 = _ratio_jet(P2r, Qr, 2m)
        _, rden = shared_q_via_determinant([rjet1, rjet2], m)
        @test rden == Qr

        # Canary 3 (AAA): a per-component AAA pole coincides with a shared-Q
        # root (diagnostic-grade, 1e-9).
        Zc = ComplexF64[exp(2im * pi * k / 200) for k = 0:199]
        F1 = ComplexF64[_ratval(P1, Qden, z) for z in Zc]
        aaa_poles = _sortc(aaa(Zc, F1).poles)
        sq_roots  = _sortc(_polyroots(den))
        @test aaa_poles ≈ sq_roots atol = 1e-9
    end

    # =========================================================================
    # SP.6 — scale covariance of the zero-input test (bug padetaylor-ked0).
    #
    # GGT 2013 md:213 / Algorithm 2 step 2 (md:225): every threshold is
    # RELATIVE, τ = tol·‖c‖₂, so the algorithm is homogeneous under c ↦ αc.
    # The pre-fix guard compared ‖c‖₂ against the ABSOLUTE `tol`, so a valid
    # jet scaled by α = 1e-20 threw "indistinguishable from zero".  Pinned:
    # numerators scale by α, the denominator and all degrees are unchanged,
    # for α ∈ {1e-20, 1e20}; the genuinely-zero jet still throws.
    # =========================================================================
    @testset "SP.6 rescaling invariance of the zero-input guard (ked0)" begin
        Qden = [1.0, -0.4, 0.2]
        P1   = [1.0, 0.3]
        P2   = [2.0, -0.7, 0.5]
        m    = 2
        jet1 = _ratio_jet(P1, Qden, 2m)
        jet2 = _ratio_jet(P2, Qden, 2m)
        nums0, den0 = shared_denominator_pade([jet1, jet2], m)
        @test den0 ≈ Qden atol = 1e-12
        for α in (1e-20, 1e20)
            nums, den = shared_denominator_pade([α .* jet1, α .* jet2], m)
            @test length(den) == length(den0)
            @test den ≈ den0 atol = 1e-12
            for i in 1:2
                @test length(nums[i]) == length(nums0[i])
                @test nums[i] ≈ α .* nums0[i] rtol = 1e-12
            end
        end
        # Genuinely-zero input (‖c‖ = 0) must still fail loud with the
        # Rule-1 suggestion text, independent of tol.
        for tol in (1e-14, 0.0)
            e = try; shared_denominator_pade([zeros(7), zeros(7)], 3; tol = tol); catch e; e; end
            @test e isa ErrorException
            msg = lowercase(sprint(showerror, e))
            @test occursin("zero", msg) && occursin("suggestion", msg)
        end
    end

    # =========================================================================
    # SP.7 — non-finite coefficients / tol are rejected (bug padetaylor-lbqb).
    # An upstream Taylor-jet overflow must surface as an ArgumentError with a
    # concrete suggestion, never as a silent zero / garbage approximant.
    # =========================================================================
    # =========================================================================
    # MUTATION-PROOF RECORD for SP.6 / SP.7 (beads padetaylor-ked0 / -lbqb).
    # Each mutation was applied by hand to src/SharedPade.jl, the suite run
    # with `julia --project=. test/shared_pade_test.jl`, then the file
    # restored and the suite re-confirmed GREEN (192/192).
    #
    # ── Mutation M5 — absolute zero-input test (the ked0 bug, restored) ─────
    #   Edit: `iszero(cnorm) && throw(ErrorException(`
    #         →  `cnorm ≤ tol_t && throw(ErrorException(`
    #   Result: 1 of 177 RED — SP.6 errored at α = 1e-20 (the valid jet was
    #     rejected as "indistinguishable from zero").  Every other testset
    #     GREEN: the bug is invisible to unit-amplitude jets, which is why it
    #     survived the triple-oracle suite.
    #
    # ── Mutation M6 — drop both lbqb guards ─────────────────────────────────
    #   Edit: `all(isfinite, jet) || throw(…)` → `true || throw(…)` and
    #         `(0 ≤ tol < Inf) || throw(…)`   → `true || throw(…)`.
    #   Result: 11 of 192 RED, all in SP.7 (Inf/-Inf/NaN coefficient and
    #     Inf/negative/NaN tol each either returned silently or threw the
    #     wrong exception type).
    # =========================================================================
    @testset "SP.7 non-finite input validation (lbqb)" begin
        Qden = [1.0, -0.4, 0.2]
        jet  = _ratio_jet([1.0, 0.3], Qden, 4)
        for bad in (Inf, -Inf, NaN)
            j2 = copy(jet); j2[3] = bad
            e = try; shared_denominator_pade([jet, j2], 2); catch e; e; end
            @test e isa ArgumentError
            msg = sprint(showerror, e)
            @test occursin("overflowed", msg) && occursin("reduce h or order", msg)
        end
        for badtol in (Inf, -1e-14, NaN)
            e = try; shared_denominator_pade([jet, jet], 2; tol = badtol); catch e; e; end
            @test e isa ArgumentError
            @test occursin("tol", sprint(showerror, e))
        end
    end

end
