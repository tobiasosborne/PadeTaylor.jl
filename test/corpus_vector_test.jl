# test/corpus_vector_test.jl
#
# ============================================================================
# Ground-truth corpus: VECTOR shared-denominator pole-bridging — the v0.2
# KEYSTONE (bead padetaylor-jmcf, epic padetaylor-p1v0).
#
# WHAT THIS FILE PINS
# -------------------
# The v0.2 keystone claim (src/VectorStepper.jl, "Why one shared Q"): the
# components of a vector meromorphic ODE solution all blow up at the SAME
# movable pole, so a SINGLE shared denominator Q per step can step the
# whole vector THROUGH that pole.  Until this file, NO test walked the
# vector solver across a shared pole and checked the value on the far side
# against an external closed form (catalogue gap §7.549, md:549).  This
# file nails it against three independent closed-form oracles, all
# bridged inside a SINGLE shared-Q step:
#
#   CVE.1  shared-denom-canonical  y=(1/(z-z0),-1/(z-z0)²), z0=1.5
#          — the cleanest shared-Q oracle: y1 order-1 pole, y2 order-2
#            pole, BOTH at z0.  Exact rational; machine-precision bridge.
#   CVE.2  tan-sec²                y=(tan z, sec² z), shared pole π/2
#          — y1 simple pole, y2 double pole; transcendental (mpmath).
#   CVE.3  PII u₁ companion        y=(-1/t, 1/t²), shared pole t=0
#          — PII(α=1) rational; clean cell-A/B dispatch probe.
#   CVE.4  tan-sec²-vector-jet     order-30 d=2 jet vs the independent
#            mpmath recurrence (catalogue tan-sec2-vector-jet-b13).
#
# GROUND TRUTH each case embodies (all INDEPENDENTLY verified — Law 1)
# -------------------------------------------------------------------
# Every ODE was sympy-verified residual-0 and every value re-derived from
# the closed form / mpmath, NOT from the Julia routine under test
# (external/probes/corpus-oracles/vector/{capture.py,verify_cm.py}):
#   * CVE.1: y1'=y2, y2'=2y2²/y1; sympy residual 0; jet c1[k]=(-1)^k/z0^{k+1}
#     (catalogue coupled-riccati-shared-denom-canonical, md:769-784).
#   * CVE.2: y1'=y2, y2'=2y1y2 with (tan, sec²); sympy residual 0
#     (catalogue coupled-riccati-tan-sec2, md:71/634-649).
#   * CVE.3: y1'=y2, y2'=2y1³+t·y1+1 with (-1/t, 1/t²); sympy residual 0
#     (catalogue pii-u1-companion-2vector, md:94).
#   * CVE.4: recurrence c1[n+1]=c2[n]/(n+1), c2[n+1]=2Σc1c2/(n+1) AND
#     mpmath.diff agree to 1e-40 (capture.py); catalogue md:1616.
#
# DE-DUPLICATION vs the existing suite (CLAUDE.md Rule 3, skepticism)
# ------------------------------------------------------------------
#   * vector_problems_test.jl VP.1.3 integrates a meromorphic vector
#     system y=(1/(1-z),2/(1-z)) but ONLY on [0,0.8] — short of the pole,
#     never crossing it.  CVE.1/CVE.2/CVE.3 BRACKET the pole strictly
#     inside one step and assert the value PAST it — the keystone that
#     VP.1.3 explicitly defers ("Pole-crossing … later beads").
#   * calogero_moser_test.jl (CM.*) integrates a SMOOTH repulsive CM
#     stretch [0,0.4] with NO pole — not a pole bridge.  (The catalogue's
#     cm-n2-imaginary-collision pole-crossing variant is DEFERRED — see
#     corpus_vector_polefield_test.jl CVE.6, where its closed form is
#     shown to be WRONG for the repo's repulsive ODE.)
#   * shared_pade_dispatch_test.jl SPD.2 pins the cell PICK on a
#     pole-crossing jet, but never the end-to-end stepped VALUE past the
#     pole; CVE.* pins the value (the load-bearing keystone assertion).
#   * vector_coefficients_test.jl pins the d=1 reduction + a harmonic
#     jet; no existing test pins the (tan, sec²) d=2 meromorphic jet
#     against an independent transcendental-recurrence oracle (CVE.4).
#   * VPO.2 already pins vector_path_network_solve reaching the FW ℘
#     u(30); the catalogue pi1-okamoto-2vector reduces (t=0 slice) to
#     exactly that — REFERENCED, not re-pinned (a true duplicate).
#
# TOLERANCES (justified by method/precision, never fitted)
# --------------------------------------------------------
#   * Exact-rational bridges (CVE.1, CVE.3): the shared-Q (15,15) Padé of
#     an exactly-rational jet recovers the rational to machine ε — the
#     stepper measured ≤4e-14 past the pole (probe).  atol 1e-11 carries
#     ~3 orders of headroom.
#   * CVE.2 tan/sec²: a SINGLE h=3 step spans 0→3 with the pole at the
#     interior t=0.524.  y1=tan to ~1e-11; y2=sec² (the steeper DOUBLE
#     pole) to ~2.7e-6 — the genuine single-big-step shared-Q accuracy
#     for the double-pole component (probe).  rtol 1e-9 (tan) / 1e-5
#     (sec²), each ~2 orders above the measured floor.
#   * CVE.4 jet: 50-dps mpmath strings; Float64 parse hits ~1e-13 rel
#     (machine ε on O(1)…O(7e4) coeffs); BigFloat-256 hits ~1e-70.
#
# MUTATION-PROOF (Rule 4): footer block — each load-bearing assertion was
# perturbed to confirm RED (incl. a SHARED-Q-targeted mutant: cell-B
# window) then restored byte-clean.  See corpus_vector_polefield_test.jl
# for the Froissart-filter cases (CVE.5) and the deferred CM (CVE.6).
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor.VectorProblems: VectorPadeTaylorProblem, vector_solve_pade
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients
using PadeTaylor.SharedPadeDispatch: shared_pade_select

include(joinpath(@__DIR__, "_oracle_corpus_vector.jl"))

# The rescale-by-h^k the stepper applies before the shared-Q solve — used
# only to ask WHICH cell the dispatch would pick at a bracketing step
# (the value assertions below own the keystone; the cell is diagnostic).
_resc(jet, h) = [jet[k] * h^(k - 1) for k in 1:length(jet)]
function _cell_at(f, z0, y0, h, order)
    jets = vector_taylor_coefficients(f, z0, y0, order)
    rj   = [_resc(j, h) for j in jets]
    shared_pade_select(rj, order ÷ 2, f, z0, h; diagnostics = true)[3].chosen_cell
end

@testset "Corpus: vector shared-Q pole-bridge (CVE)" begin

    # ----------------------------------------------------------------------
    # CVE.1  shared-denom-canonical: bridge z0=1.5 in ONE shared-Q step.
    # y1'=y2, y2'=2y2²/y1; y=(1/(z-1.5), -1/(z-1.5)²).  IC at z=0.
    # ----------------------------------------------------------------------
    @testset "CVE.1  y=(1/(z-1.5),-1/(z-1.5)²): bridge the shared pole z=1.5" begin
        z0 = 1.5
        f  = (z, y) -> [y[2], 2 * y[2]^2 / y[1]]
        y0 = [1 / (0 - z0), -1 / (0 - z0)^2]          # (-2/3, -4/9), exact
        prob = VectorPadeTaylorProblem(f, y0, (0.0, 3.0); order = 30)
        sol  = vector_solve_pade(prob; h = 3.0)        # pole z0 at t=0.5, interior
        @test length(sol.h) == 1                       # ONE step brackets the pole
        # Before the pole, exact-rational machine precision.
        yb = sol(0.5)
        @test isapprox(yb[1], 1 / (0.5 - z0); atol = 1e-11)   # = -1
        @test isapprox(yb[2], -1 / (0.5 - z0)^2; atol = 1e-11) # = -1
        # PAST the pole — the keystone assertion: both components match the
        # exact closed form on the far side of z0 = 1.5.
        for zq in (2.0, 2.5)
            yq = sol(zq)
            @test isapprox(yq[1], 1 / (zq - z0); atol = 1e-11)
            @test isapprox(yq[2], -1 / (zq - z0)^2; atol = 1e-11)
        end
        # ADR-0028: cell B (Mano–Tsuda) wins the pole-crossing step
        # (Amendment 2 — diagnostic; the VALUE above is load-bearing).
        @test _cell_at(f, 0.0, y0, 3.0, 30) === :square
    end

    # ----------------------------------------------------------------------
    # CVE.2  tan-sec²: bridge the shared pole at π/2 in ONE step on [0,3].
    # y1'=y2, y2'=2y1y2; y=(tan z, sec² z); IC y(0)=(0,1).
    # ----------------------------------------------------------------------
    @testset "CVE.2  y=(tan z, sec² z): bridge the shared pole at π/2" begin
        f = (z, y) -> [y[2], 2 * y[1] * y[2]]
        prob = VectorPadeTaylorProblem(f, [0.0, 1.0], (0.0, 3.0); order = 30)
        sol  = vector_solve_pade(prob; h = 3.0)        # π/2≈1.571 at t≈0.524
        @test length(sol.h) == 1
        # PAST the pole: y1=tan to the simple-pole floor; y2=sec² (DOUBLE
        # pole) to the steeper double-pole floor (tolerances in header).
        y2v = sol(2.0)
        @test isapprox(y2v[1], CRTS_TAN_2p0;  rtol = 1e-9)
        @test isapprox(y2v[2], CRTS_SEC2_2p0; rtol = 1e-5)
        y25 = sol(2.5)
        @test isapprox(y25[1], CRTS_TAN_2p5;  rtol = 1e-9)
        @test isapprox(y25[2], CRTS_SEC2_2p5; rtol = 1e-5)
        # The pole really is bridged: tan(2) is finite, negative, and the
        # naïve sec²=1+tan² identity holds on the recovered values.
        @test isfinite(y2v[1]) && y2v[1] < 0
        @test isapprox(y2v[2], 1 + y2v[1]^2; rtol = 1e-5)   # sec²=1+tan²
        @test _cell_at(f, 0.0, [0.0, 1.0], 3.0, 30) === :square
    end

    # ----------------------------------------------------------------------
    # CVE.3  PII u₁ companion: bridge the shared pole t=0 on [-1,1].
    # PII(α=1) u=-1/t; companion y=(-1/t, 1/t²); y1'=y2, y2'=2y1³+t·y1+1.
    # ----------------------------------------------------------------------
    @testset "CVE.3  PII u₁=(-1/t,1/t²): bridge the shared pole at t=0" begin
        f = (z, y) -> [y[2], 2 * y[1]^3 + z * y[1] + 1]
        y0 = [-1 / (-1.0), 1 / (-1.0)^2]               # at t=-1: (1, 1) exact
        prob = VectorPadeTaylorProblem(f, y0, (-1.0, 1.0); order = 30)
        sol  = vector_solve_pade(prob; h = 2.0)        # pole t=0 at t-resc 0.5
        @test length(sol.h) == 1
        # Before the pole (t=-0.5): exact rational.
        yb = sol(-0.5)
        @test isapprox(yb[1], 2.0; atol = 1e-11)        # -1/(-0.5)=2
        @test isapprox(yb[2], 4.0; atol = 1e-11)        # 1/(-0.5)²=4
        # PAST the pole — the keystone: both components past t=0.
        for tq in (0.5, 1.0)
            yq = sol(tq)
            @test isapprox(yq[1], -1 / tq;  atol = 1e-11)
            @test isapprox(yq[2],  1 / tq^2; atol = 1e-11)
        end
        @test _cell_at(f, -1.0, y0, 2.0, 30) === :square
    end

    # ----------------------------------------------------------------------
    # CVE.4  tan-sec²-vector-jet-b13: the d=2 order-30 jet vs the
    # INDEPENDENT mpmath recurrence (NOT the Julia routine vs itself).
    # y1'=y2, y2'=2y1y2; y1=tan(z+π/4), y2=sec²(z+π/4); IC (1,2).
    # ----------------------------------------------------------------------
    @testset "CVE.4  vector_taylor_coefficients d=2 jet vs mpmath recurrence" begin
        f = (z, y) -> [y[2], 2 * y[1] * y[2]]

        # Float64: machine-ε match against the 50-dps recurrence strings.
        cf = vector_taylor_coefficients(f, 0.0, [1.0, 2.0], 30)
        @test length(cf) == 2
        @test length(cf[1]) == 31 && length(cf[2]) == 31
        @test eltype(cf[1]) == Float64
        for k in 0:30
            @test isapprox(cf[1][k + 1], parse(Float64, CTSJ_C1[k + 1]); rtol = 1e-13)
            @test isapprox(cf[2][k + 1], parse(Float64, CTSJ_C2[k + 1]); rtol = 1e-13)
        end
        # The catalogue spot values (exact): c1[1]=2, c1[2]=2, c2[1]=4.
        @test cf[1][2] == 2.0 && cf[1][3] == 2.0 && cf[2][2] == 4.0

        # BigFloat-256: the same jet at extended precision vs the same
        # mpmath strings (the IC is the SAME on both sides — honest).
        setprecision(BigFloat, 256) do
            cb = vector_taylor_coefficients(f, BigFloat(0),
                                            BigFloat[1, 2], 30)
            @test eltype(cb[1]) == BigFloat
            for k in 0:30
                e1 = parse(BigFloat, CTSJ_C1[k + 1])
                e2 = parse(BigFloat, CTSJ_C2[k + 1])
                @test abs(cb[1][k + 1] - e1) / abs(e1) < BigFloat("1e-44")
                @test abs(cb[2][k + 1] - e2) / abs(e2) < BigFloat("1e-44")
            end
        end
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# src/ confirmed clean via `git diff --stat src/` after each restore.
#
# M-oracle (CVE.1): in CVE.1's far-side loop, flip the expected sign:
#   `isapprox(yq[1], 1/(zq-z0); …)` → `isapprox(yq[1], -1/(zq-z0); …)`.
#   Result: CVE.1 went RED at the past-pole y1 line for BOTH z=2.0 and
#   z=2.5 (the stepped value 2.0/1.0 no longer matches the negated
#   expected −2.0/−1.0); every other CVE assertion stayed GREEN.  Restored
#   the literal byte-for-byte (verified via `git diff`).
#
# M-impl (the SHARED-Q-targeted mutant — cell-B window): in
#   src/SharedPadeCellB.jl `_block_b`, change the window index
#   `idx = m + rr - cc + 1` → `m + rr - cc + 2` (cell B's +1 window → cell
#   A's +2 window — the very off-by-one ADR-0028 A1.2 corrects).  Result:
#   CVE.2 went RED — the mutated cell-B mis-builds the (tan, sec²) bridge,
#   degrading the past-pole tan(2.5) from ~1e-11 to ~2e-8 error (over the
#   1e-9 bar).  CVE.1/CVE.3 stayed GREEN: their jets are EXACTLY rational,
#   so when the broken cell B scores a worse relative-ODE defect the
#   A-default dispatch (SharedPadeDispatch tie-break) FALLS BACK to cell A,
#   which still bridges the exact-rational pole to machine precision.  This
#   is an honest coverage fact — cell B is load-bearing precisely where its
#   degree advantage matters (the transcendental double pole, CVE.2), and
#   the dispatch's A-fallback protects the rational cases (exactly the
#   ADR-0028 design).  Restored src/SharedPadeCellB.jl byte-for-byte.
#
# M-impl2 (CVE.4 jet — the 2nd-order vector convolution): in
#   src/VectorCoefficients.jl, change `y[i][k] = f_t[i][k-1] / k` →
#   `/ (k+1)` (the off-by-one in the integration relation).  Result:
#   CVE.4 RED from c1[1] onward (rel err O(1)) at both Float64 and
#   BigFloat — confirms the d=2 jet is pinned by the INDEPENDENT mpmath
#   recurrence, not self-consistent.  Restored byte-for-byte.
#
# Run standalone with:
#   julia --project=. test/corpus_vector_test.jl
# ============================================================================
