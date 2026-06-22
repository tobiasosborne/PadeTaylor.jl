# test/corpus_higher_order_pole_test.jl
#
# ============================================================================
# Ground-truth corpus: HIGHER-ORDER POLES (orders 3,4,5) — the corpus tops
# out at DOUBLE poles before this file (bead padetaylor-gf8i, epic
# padetaylor-25og; plan Family B, docs/test_corpus/02_corpus_extension_plan.md:177-200).
#
# WHAT THIS FILE PINS
# -------------------
# The v1 corpus has dense coverage of simple + double poles but NO order-3+
# poles (plan blind-spot matrix, md:35).  This file fills orders 3,4,5 with
# three independently-derived cases, each a different pole *morphology*:
#
#   CHop.1  order-k FIXED pole   u''=k(k+1)u/(1-z)², u=(1-z)^-k, k=2..5.
#           A fixed (non-movable) pure order-k pole, residue 0.  VERIFIED
#           PATH FACT (Rule 3, against the brief's premise): the UN-rescaled
#           [k/0] Padé IS the function ⇒ the classical Toeplitz solve of the
#           raw jet is SINGULAR (probe: SingularException) and robust_pade
#           routes it to the SVD path — but solve_pade's STEPPER rescales the
#           jet by h^k, which REGULARISES the Toeplitz into a well-conditioned
#           (15,15) block, so the Float64 stepper VALUE bridges via the
#           :classical path (no crash; no SVD fallback at the stepper).  The
#           SVD path is exercised directly by the multiplicity check below
#           (un-rescaled robust_pade(c,30,30;method=:svd)) and by the
#           BigFloat-256 case (its default method IS :svd end-to-end).
#           Multiplicity = pre-collapse Q-root cluster count.
#   CHop.3  mixed-order shared pole vector  f=[2y₂, 3y₁²], IC (1,1),
#           y=((1-z)^-2,(1-z)^-3): orders 2 AND 3 sharing z=1 over ONE cubic
#           shared Q.  Extends CVE.1 (orders 1,2) to max order 3.
#   CHop.2  ℘′ TRIPLE pole  f=[y₂, 6y₁²-g₂/2], (y₁,y₂)=(℘,℘′), lemniscatic
#           g₂=1,g₃=0.  ℘′ has a triple pole (Laurent leading −2/w³),
#           residue 0; bridged via vector_path_network_solve around the real
#           lattice pole at z=2ω≈3.708.
#
# GROUND TRUTH each case embodies (all INDEPENDENTLY verified — Law 1)
# -------------------------------------------------------------------
# external/probes/corpus-oracles/higher-order-poles/capture.py (sympy 1.14 +
# mpmath 1.3) gates EVERY value before emitting it:
#   * CHop.1: sympy residual u''−k(k+1)u/(1-z)² ≡ 0 (all k); exact rationals
#     u(1.05)=(−0.05)^-k ∈ {400,−8000,160000,−3200000}, u(0.5)=2^k,
#     u'(1.05)=k(−0.05)^-(k+1).  Written inline (exact integers).
#   * CHop.3: sympy residual y₁'−2y₂ ≡ 0 AND y₂'−3y₁² ≡ 0; exact rationals
#     y(2.0)=(1,−1), y(2.5)=(4/9,−8/27).  Written inline.
#   * CHop.2: mpmath dps=50, GATED on the invariant (℘')²=4℘³−g₂℘−g₃ (≤1e-40)
#     and the −2/w³ Laurent lead; ℘ via the validated Jacobi-sn recipe, ℘'
#     via mp.diff.  Transcendental pins in _oracle_corpus_higher_order_pole.jl.
#
# TOLERANCES (justified by the PER-ORDER Padé split width / floor, never fitted)
# -----------------------------------------------------------------------------
# A pure order-k pole plants k near-coincident Q-roots; the past-pole Padé
# floor and the root-split width BOTH grow with k (measured below, matching
# plan discovery #3/#4):
#   * CHop.1 past-pole value relerr (measured): k=2 ~1e-14, k=3 ~1e-13,
#     k=4 ~3e-11, k=5 ~2.8e-9 ⇒ per-k tol {1e-13,1e-12,1e-10,1e-8}, each ~1
#     order above the measured floor.  BigFloat-256 k=2: ~3e-74 ⇒ 1e-30.
#   * CHop.1 extract_poles z-location dist (measured): k=2 ~1.8e-8, k=3
#     ~5.7e-6, k=4 ~1.1e-4, k=5 ~8e-4 ⇒ per-k tol {1e-7,1e-5,1e-3,5e-3}
#     (the Padé root-split width, NOT fitted).
#   * CHop.1 multiplicity: the centered (30,30) SVD Padé denominator has
#     EXACTLY k roots clustered at z=1 (pre-collapse count; extract_poles
#     returns LOCATIONS only — plan discovery #3).
#   * CHop.3 exact-rational shared-Q bridge: relerr ~1e-13 ⇒ atol 1e-11.
#   * CHop.2 ℘/℘′ past-pole via path-network + extrapolate: relerr ~1e-13
#     (well inside the brief's 1e-7); pole-location dist ~1.2e-4 (triple-pole
#     split width) ⇒ atol 1e-3.
#
# REFERENCES: plan Family B (md:177-200); DLMF 23.x (Weierstrass ℘);
# discovery #3 (extract_poles has no order field, md:70-76), #4 (pure-power
# Toeplitz singular ⇒ SVD, md:78-80).  Oracle: capture.py (residual==0 /
# invariant gated).
#
# MUTATION-PROOF (Rule 4): footer block — one M-oracle (flip a pinned value)
# and one M-impl (perturb the load-bearing SVD/shared-Q path) were ACTUALLY
# executed RED and restored byte-for-byte.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor: taylor_coefficients_2nd
using PadeTaylor.RobustPade: robust_pade
using PadeTaylor.VectorProblems: VectorPadeTaylorProblem, vector_solve_pade
using PadeTaylor.VectorPathNetwork: vector_path_network_solve
using PadeTaylor.VectorPoleField: extract_poles_shared_q
using Polynomials: Polynomial, roots

include(joinpath(@__DIR__, "_oracle_corpus_higher_order_pole.jl"))

# Per-k past-pole exact rationals u(1.05)=(-1/20)^-k, u(0.5)=2^k,
# u'(1.05)=k(-1/20)^-(k+1).  Exact integers (capture.py, sympy-verified).
const CHOP1_U_1p05  = Dict(2 => 400, 3 => -8000, 4 => 160000, 5 => -3200000)
const CHOP1_U_0p5   = Dict(2 => 4,   3 => 8,     4 => 16,     5 => 32)
const CHOP1_UP_1p05 = Dict(2 => -16000, 3 => 480000, 4 => -12800000, 5 => 320000000)
# Per-k tolerances justified by the measured floor / split width (header).
const CHOP1_VAL_TOL = Dict(2 => 1e-13, 3 => 1e-12, 4 => 1e-10, 5 => 1e-8)
const CHOP1_LOC_TOL = Dict(2 => 1e-7,  3 => 1e-5,  4 => 1e-3,  5 => 5e-3)

@testset "Corpus: higher-order poles (CHop)" begin

    # ----------------------------------------------------------------------
    # CHop.1  order-k FIXED pole: u''=k(k+1)u/(1-z)², u=(1-z)^-k, sweep k=2..5.
    # solve_pade bridges the pole NO-CRASH (the raw [k/0] jet is classical-
    # singular, but the stepper's h-rescaled (15,15) Toeplitz is regular —
    # see header path fact).  Pole at z=1 (t≈0.667 inside [0,1.5]), residue
    # 0, multiplicity k (the SVD path is exercised by the multiplicity check).
    # ----------------------------------------------------------------------
    @testset "CHop.1  u=(1-z)^-$k: bridge the order-$k pole at z=1" for k in 2:5
        f(z, u, up) = (k * (k + 1)) * u / (1 - z)^2
        prob = PadeTaylorProblem(f, (1.0, float(k)), (0.0, 1.5); order = 30)
        sol  = solve_pade(prob; h = 1.5)          # pole z=1 at t≈0.667, interior
        u05,  _      = sol(0.5)
        u105, up105  = sol(1.05)                       # PAST the order-k pole
        # Past-pole exact rationals (per-k tol justified by the measured floor).
        @test isapprox(u05,   CHOP1_U_0p5[k];   rtol = CHOP1_VAL_TOL[k])
        @test isapprox(u105,  CHOP1_U_1p05[k];  rtol = CHOP1_VAL_TOL[k])
        @test isapprox(up105, CHOP1_UP_1p05[k]; rtol = CHOP1_VAL_TOL[k])
        # The pole is really bridged: finite + correct sign (alternates with k).
        @test isfinite(u105) && sign(u105) == (-1)^k

        # MULTIPLICITY (plan discovery #3): the centered SVD Padé denominator
        # has EXACTLY k roots clustered at z=1 (pre-collapse count).  Build
        # the full order-(m+n)=60 jet centred at z₀=0 (robust_pade(30,30)
        # consumes m+n+1=61 coefficients; a short jet zero-pads and the SVD
        # collapses the denominator to degree 0) and root :svd's b.
        c  = taylor_coefficients_2nd(f, 0.0, 1.0, float(k), 60)
        pa = robust_pade(c, 30, 30; method = :svd)
        rts = roots(Polynomial(pa.b))
        near1 = filter(r -> abs(r - 1.0) < 0.2, rts)
        @test length(near1) == k                       # order-k ⇒ k clustered roots

        # extract_poles returns LOCATIONS only: pin the z=1 location to the
        # per-k split-width tol.  (min_support=1 — single-trajectory chain.)
        poles = extract_poles(sol; radius_t = 5.0, min_residue = 1e-8,
                              cluster_atol = 1e-1, min_support = 1)
        @test !isempty(poles)
        @test minimum(abs(p - 1.0) for p in poles) < CHOP1_LOC_TOL[k]
    end

    @testset "CHop.1.bf  k=2 at BigFloat-256: tighter past-pole accuracy" begin
        setprecision(BigFloat, 256) do
            k = 2
            f(z, u, up) = BigFloat(k * (k + 1)) * u / (1 - z)^2
            prob = PadeTaylorProblem(f, (BigFloat(1), BigFloat(k)),
                                     (BigFloat(0), BigFloat(3) / 2); order = 40)
            sol  = solve_pade(prob; h = BigFloat(3) / 2)
            u105, _ = sol(BigFloat(105) / 100)
            uex = (BigFloat(-5) / 100)^(-k)            # exact = 400
            @test isapprox(u105, uex; rtol = BigFloat("1e-30"))
        end
    end

    # ----------------------------------------------------------------------
    # CHop.3  mixed-order shared pole vector: f=[2y₂, 3y₁²], IC (1,1),
    # y=((1-z)^-2,(1-z)^-3).  y₁ order-2, y₂ order-3, SAME pole z=1, over ONE
    # cubic shared Q.  Extends CVE.1 (orders 1,2) to max order 3.
    # ----------------------------------------------------------------------
    @testset "CHop.3  y=((1-z)^-2,(1-z)^-3): mixed-order shared pole z=1" begin
        f  = (z, y) -> [2 * y[2], 3 * y[1]^2]
        y0 = [1.0, 1.0]
        prob = VectorPadeTaylorProblem(f, y0, (0.0, 3.0); order = 30)
        sol  = vector_solve_pade(prob; h = 3.0)        # pole z=1 at t≈0.333, interior
        @test length(sol.h) == 1                       # ONE shared-Q step brackets it
        # Before the pole: exact rationals (1/(1-0.5)²=4, 1/(1-0.5)³=8).
        yb = sol(0.5)
        @test isapprox(yb[1], 4.0; atol = 1e-11)
        @test isapprox(yb[2], 8.0; atol = 1e-11)
        # PAST the pole — the keystone: both DIFFERENT orders past z=1.
        y20 = sol(2.0)
        @test isapprox(y20[1],  1.0;          atol = 1e-11)   # (1-2)^-2 = 1
        @test isapprox(y20[2], -1.0;          atol = 1e-11)   # (1-2)^-3 = -1
        y25 = sol(2.5)
        @test isapprox(y25[1],  4 / 9;        atol = 1e-11)   # (1-2.5)^-2
        @test isapprox(y25[2], -8 / 27;       atol = 1e-11)   # (1-2.5)^-3
        # The two components blow up at the SAME pole with DIFFERENT orders:
        # the order-2/order-3 ratio y₂/y₁ = (1-z)^-1 holds past the pole.
        @test isapprox(y20[2] / y20[1], 1 / (1 - 2.0); atol = 1e-11)
    end

    # ----------------------------------------------------------------------
    # CHop.2  ℘′ TRIPLE pole: f=[y₂, 6y₁²-g₂/2], (y₁,y₂)=(℘,℘′), lemniscatic
    # g₂=1,g₃=0.  ℘′ has a triple pole (Laurent leading −2/w³), residue 0.
    # Bridge the real lattice pole at z=2ω≈3.708 via the path-network walk
    # (curved complex route around the real-axis pole), then extrapolate the
    # past-pole value at z_tgt=2ω+0.5.
    # ----------------------------------------------------------------------
    @testset "CHop.2  ℘′ triple pole: path-network bridge of the z=2ω pole" begin
        g2 = 1.0
        f  = (z, y) -> [y[2], 6 * y[1]^2 - g2 / 2]
        y0 = ComplexF64[CHOP2_P_Z0, CHOP2_PP_Z0]       # (℘(z₀), ℘′(z₀)) at z₀=0.9
        ztgt = complex(CHOP2_ZTGT)
        prob = VectorPadeTaylorProblem(f, y0,
                                       (complex(CHOP2_Z0), ztgt + 1.0); order = 30)
        sol = vector_path_network_solve(prob, ComplexF64[ztgt]; order = 30,
                  h = 0.5, step_policy = :max_q_root, adaptive = true,
                  on_target_failure = :throw, fine_grid = ComplexF64[ztgt],
                  extrapolate = true)
        @test length(sol.visited_z) ≥ 1
        # Past-pole value: ℘ and ℘′ on the far side of the z=2ω triple pole.
        yq = sol.grid_y[1]
        @test all(isfinite, yq)
        @test isapprox(yq[1], CHOP2_P_ZTGT;  rtol = 1e-7)   # ℘  past pole
        @test isapprox(yq[2], CHOP2_PP_ZTGT; rtol = 1e-7)   # ℘′ past pole (triple)
        # Pole LOCATION: extract_poles_shared_q places z=2ω to the triple-pole
        # split-width tol (~1.2e-4 measured).  min_support=2 (cross-node).
        poles = extract_poles_shared_q(sol; radius_t = 5.0, min_support = 2)
        @test !isempty(poles)
        @test minimum(abs(p - CHOP2_POLE_REAL) for p in poles) < 1e-3
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# `git diff src/` confirmed EMPTY after each restore.
#
# M-oracle (CHop.1, k=3): flip the pinned past-pole rational
#   `CHOP1_U_1p05[3]` from `-8000` to `+8000`.  Result: the CHop.1 k=3
#   `isapprox(u105, CHOP1_U_1p05[k]; …)` assertion went RED (got −8000,
#   expected +8000) and ONLY at k=3; k=2,4,5 + CHop.2/CHop.3 stayed GREEN.
#   Restored the `-8000` literal; suite GREEN.
#
# M-impl-A (CHop.1 + BF, the scalar stored-Padé evaluation): in
#   src/PadeStepper.jl `_evaluate_pade`, append `num *= 1.0000001` after the
#   numerator Horner loop (a 1e-7 corruption of every evaluated scalar Padé
#   value).  Result: 9 RED — every CHop.1 VALUE assertion (k=2..5 u05 + u105)
#   plus the BigFloat-256 case.  GREEN throughout: the `up105` derivative
#   pins (computed by the separate `_evaluate_pade_deriv` chain — the same
#   coverage fact CPB.* records), the multiplicity Q-root counts, and the
#   extract_poles LOCATION pins (the denominator/roots are untouched by a
#   numerator rescale).  CHop.2/CHop.3 stayed GREEN: the VECTOR stepper /
#   path-network evaluate through their own routine, not this scalar one —
#   so this mutant cleanly isolates the scalar value path.  Restored
#   src/PadeStepper.jl byte-for-byte; suite GREEN.
#
# M-impl-B (CHop.3 + CHop.2, the shared-Q vector jet convolution): in
#   src/VectorCoefficients.jl, change the integration relation
#   `y[i][k] = f_t[i][k-1] / k` → `/ (k+1)` (the same off-by-one CVE.4
#   mutation-proves).  Result: 10 RED — CHop.3 (both components past z=1 at
#   2.0 and 2.5, plus the order-2/order-3 ratio invariant) AND CHop.2 (℘ and
#   ℘′ past the triple pole, plus the pole LOCATION — a corrupted jet moves
#   the shared-Q roots).  CHop.1 (scalar `taylor_coefficients_2nd`, a
#   different jet builder) stayed GREEN.  Confirms the vector shared-Q
#   bridges are pinned by the INDEPENDENT oracle (exact rationals / mpmath
#   ℘′), not self-consistent.  Restored src/VectorCoefficients.jl
#   byte-for-byte.
#
# Supplementary M-impl-C (the SVD path the multiplicity check + BF case ride):
#   in src/RobustPade.jl `_trim_and_normalise`, `a ./= b1` → `a ./=
#   (b1 * 1.0000001)`.  Result: the BigFloat-256 case (default method :svd
#   end-to-end) went RED at 1e-30; the Float64 CHop.1 stepper VALUES stayed
#   GREEN — the verified path fact (header): the h-rescaled Float64 Toeplitz
#   is REGULAR, so the stepper uses :classical, never reaching this SVD
#   normaliser.  Restored byte-for-byte.
#
# Run standalone with:
#   julia --project=. test/corpus_higher_order_pole_test.jl
# ============================================================================
