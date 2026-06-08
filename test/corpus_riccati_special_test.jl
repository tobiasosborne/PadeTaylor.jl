# test/corpus_riccati_special_test.jl
#
# ============================================================================
# Ground-truth corpus: RICCATI = log-derivative of a linear ODE, TRANSCENDENTAL
# subfamily (bead padetaylor-tyef, epic padetaylor-25og).
#
# WHAT THIS FILE PINS
# -------------------
# For any linear ODE  w″+p w′+q w = 0,  u = ±w′/w  solves a Riccati equation
# with a SIMPLE pole at every zero of w, residue +1 (or −1 for u=−w′/w)
# (docs/test_corpus/02_corpus_extension_plan.md:134-166, Family A).  When w is a
# SPECIAL FUNCTION, its zeros are tabulated to arbitrary precision → exact pole
# locations AND residues.  This is the corpus's first ROWS of simple poles at
# named special-function zeros (the rational sister CRic.3/4 covers polynomial
# zeros; here the values carry a genuine transcendental floor).  Three members:
#
#   CRic.2  Bessel    u = −J₁/J₀ = J₀′/J₀   POSITIVE-axis row j_{0,k}, res +1
#   CRic.1  Airy      u = −Ai′/Ai            NEGATIVE-axis row a_k,    res −1
#   CRic.5  parab.cyl u = −U′(0,√2z)/U(…)    COMPLEX-conjugate pair,   res −1
#
# CRic.5 BRIEF CORRECTION (plan 02:63-68, discovery #2): the recessive solution
# of w″=z²w is the parabolic-cylinder U(0,√2z), NON-oscillatory ⇒ ZERO real
# zeros.  u′=u²−z² therefore has NO real poles — only complex-conjugate pairs
# near arg≈±3π/4.  (The real-row oscillatory sister u′=u²+z² is the deferred
# CRic.6, pending a dedup decision with CPB.5.)
#
# GROUND TRUTH each case embodies
# -------------------------------
# Every recast is VERIFIED in
# external/probes/corpus-oracles/riccati-special/capture.py: the first-order
# Riccati R(z,u) is recast by u″=∂R/∂z+(∂R/∂u)·R, and capture.py GATES on the
# residual being identically 0 (sympy) AND <1e-30 against the closed-form
# special function (mpmath dps=50) before emitting any number:
#   CRic.2:  R=−u²−u/z−1  ⇒  u″ = 2u³ + 3u²/z + 2u + 2u/z² + 1/z
#   CRic.1:  R=u²−z       ⇒  u″ = 2u³ − 2uz − 1
#   CRic.5:  R=u²−z²      ⇒  u″ = 2u³ − 2uz² − 2z
# The 1/z, 1/z² RHS coefficient-singularity (CRic.2) forces the seed z₀=1, not 0.
# Pole *bridging* is asserted as the residue-±1 SIGN FLIP across each pole; the
# pole *row* is pinned by extract_poles vs the tabulated zeros (Airy/Bessel) or
# the PCF root-find (CRic.5).  Constants live in _oracle_corpus_riccati_special.jl.
#
# TOLERANCES (justified by METHOD accuracy, never fitted to pass)
# -----------------------------------------------------------------
# The only error is the local (15,15)-diagonal Padé floor of the stepper /
# path-network, NOT the special-function approximation (mpmath is exact here).
#   * before any pole (small |t| past midpoint): rtol 1e-11/1e-12 — machine eps
#     regime (measured ~1e-15…1e-13).
#   * just PAST a real pole (CRic.2 j_{0,1}≈2.40, CRic.1 a_1≈−2.34): rtol 1e-7
#     (CRic.2, measured ~4e-9) / 1e-9 (CRic.1 bridge, measured ~9e-11) — the
#     documented past-pole Padé floor (README §2 / CPB.1/CPB.2).
#   * complex bridge across CRic.5's first pole: rtol 1e-7 (measured ~1e-12).
#   * extract_poles real rows: atol 1e-7 on the bridged zero (figure-catalogue
#     Float64 spec, polefield_test.jl PF.1.1; measured ~1e-12); CRic.5's complex
#     pole: atol 1e-5 (the brief's PCF-root tol; measured ~1e-11).
#
# REFERENCES: docs/test_corpus/02_corpus_extension_plan.md:134-166 (Family A);
# DLMF 9.9 (Airy zeros a_k), 10.21 (Bessel zeros j_{0,k}), 12.11 (parabolic-
# cylinder zeros).  Oracle: capture.py (mpmath dps=50, residual==0 GATE).
#
# MUTATION-PROOF (Rule 4): see the comment block at the end of this file —
# one M-oracle (flip a pinned value → reddens) and one M-impl (perturb the
# stored Padé in src/PadeStepper.jl) were ACTUALLY executed RED and restored
# byte-for-byte (`git diff src/` empty).
# ============================================================================

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_corpus_riccati_special.jl"))

# Recasts (capture.py-gated; u″ = ∂R/∂z + (∂R/∂u)·R).
fBessel(z, u, up) = 2u^3 + 3u^2 / z + 2u + 2u / z^2 + 1 / z   # CRic.2
fAiry(z, u, up)   = 2u^3 - 2u * z - 1                          # CRic.1
fPCF(z, u, up)    = 2u^3 - 2u * z^2 - 2z                       # CRic.5

@testset "Corpus: special-function log-deriv Riccati (CRic)" begin

    # ----------------------------------------------------------------------
    # CRic.2  Bessel  u=−J₁/J₀.  POSITIVE-axis row; solve_pade (real-increasing)
    # bridges j_{0,1}≈2.405 in one window.  Seed z₀=1 (1/z RHS coeff-sing at 0).
    # ----------------------------------------------------------------------
    @testset "CRic.2  Bessel −J₁/J₀: positive-axis row via solve_pade" begin
        prob = PadeTaylorProblem(fBessel, (CRic2_U0, CRic2_UP0), (1.0, 4.0); order = 30)
        sol  = solve_pade(prob; h_max = 3.5)             # j_{0,1}≈2.405 interior

        u15, _ = sol(1.5)                                # BEFORE the pole
        u20, _ = sol(2.0)
        u30, _ = sol(3.0)                                # PAST j_{0,1}
        u35, _ = sol(3.5)
        @test isapprox(u15, CRic2_U_p1p5; rtol = 1e-12)  # before pole: machine eps
        @test isapprox(u20, CRic2_U_p2;   rtol = 1e-12)
        @test isapprox(u30, CRic2_U_p3;   rtol = 1e-7)   # past pole: Padé floor
        @test isapprox(u35, CRic2_U_p3p5; rtol = 1e-7)

        # residue +1 ⇒ sign flip across j_{0,1}: u→+∞ as z↑2.405⁻ then re-enters
        # from −∞; u(2.0)<0 (approaching the pole) and u(3.0)>0 (past it).
        @test u20 < 0 && u30 > 0 && isfinite(u30)

        # extract_poles pins j_{0,1} (the bridged zero) on the positive row.
        poles = extract_poles(sol; min_support = 1)
        @test minimum(abs(p - complex(CRic2_J1, 0.0)) for p in poles) ≤ 1e-7
    end

    # ----------------------------------------------------------------------
    # CRic.1  Airy  u=−Ai′/Ai.  NEGATIVE-axis row is z-DECREASING ⇒ solve_pade
    # (real-increasing) cannot reach it; use path_network_solve + eval_at toward
    # negative targets past a_1≈−2.338.  Optional pos-axis tight-accuracy probe.
    # ----------------------------------------------------------------------
    @testset "CRic.1  Airy −Ai′/Ai: negative-axis row via path_network" begin
        prob = PadeTaylorProblem(fAiry, (CRic1_U0 + 0im, CRic1_UP0 + 0im),
                                 (0.0 + 0im, -3.0 + 0im); order = 30)
        targets = ComplexF64[-1.0, -1.5, -2.0, -2.5]
        sol = path_network_solve(prob, targets; h = 0.3, rng_seed = 0,
                                 max_steps_per_target = 2000)

        wants = [CRic1_U_m1, CRic1_U_m1p5, CRic1_U_m2, CRic1_U_m2p5]
        for (t, w) in zip(targets, wants)
            u, _ = eval_at(sol, t)
            tol = (real(t) > -2.34) ? 1e-10 : 1e-9       # before vs past a_1
            @test isapprox(u, w + 0im; rtol = tol)
        end

        # residue −1 ⇒ sign flip across a_1: u(−2.0)<0 (approaching) vs u(−2.5)>0.
        u_m2,  _ = eval_at(sol, -2.0 + 0im)
        u_m25, _ = eval_at(sol, -2.5 + 0im)
        @test real(u_m2) < 0 && real(u_m25) > 0 && isfinite(u_m25)

        # extract_poles pins a_1 on the negative row.
        poles = extract_poles(sol; min_support = 3)
        @test minimum(abs(p - complex(CRic1_A1, 0.0)) for p in poles) ≤ 1e-7

        # Pos-axis (pole-free) tight-accuracy leg via solve_pade.
        probp = PadeTaylorProblem(fAiry, (CRic1_U0, CRic1_UP0), (0.0, 1.5); order = 30)
        solp  = solve_pade(probp; h_max = 1.5)
        @test isapprox(solp(0.5)[1], CRic1_U_p0p5; rtol = 1e-11)
        @test isapprox(solp(1.0)[1], CRic1_U_p1;   rtol = 1e-11)
    end

    # ----------------------------------------------------------------------
    # CRic.5  parabolic cylinder  u=−U′(0,√2z)/U.  NO real poles (brief
    # correction); a COMPLEX-conjugate pole pair near arg≈3π/4.  Bridge the
    # first pole P1≈−1.493+1.603i via path_network + eval_at; pin P1 vs the
    # PCF root-find.
    # ----------------------------------------------------------------------
    @testset "CRic.5  parabolic-cylinder: complex-conjugate pole pair" begin
        prob = PadeTaylorProblem(fPCF, (CRic5_U0 + 0im, CRic5_UP0 + 0im),
                                 (0.0 + 0im, -2.0 + 2.0im); order = 30)
        targets = ComplexF64[-0.5 + 0.5im, -1.0 + 1.0im, -1.7 + 1.8im]
        sol = path_network_solve(prob, targets; h = 0.3, rng_seed = 0,
                                 max_steps_per_target = 2000)

        wants = [CRic5_U_a, CRic5_U_b, CRic5_U_c]        # last one is PAST P1
        for (t, w) in zip(targets, wants)
            u, _ = eval_at(sol, t)
            @test isapprox(u, w; rtol = 1e-7)            # complex bridge floor
        end

        # extract_poles places the first complex pole P1 (residue −1).  The walk
        # also plants spurious far-roots; pin the cluster nearest P1 only.
        poles = extract_poles(sol; min_support = 3)
        @test minimum(abs(p - CRic5_P1) for p in poles) ≤ 1e-5
        # It really is COMPLEX (the brief-corrected no-real-poles signature):
        @test imag(CRic5_P1) > 1.0
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
# M-oracle (CRic.2): in test/_oracle_corpus_riccati_special.jl flip the sign of
#   CRic2_U_p3 (1.30381… → −1.30381…), the headline PAST-pole value.  Result:
#   ONLY the CRic.2 u(3.0) value assertion went RED (Evaluated: isapprox(
#   1.3038123810874516, −1.3038123810828373; rtol=1e-7)); every other assertion
#   stayed GREEN.  Restored the literal; suite GREEN.
#
# M-impl (the load-bearing stored Padé): in src/PadeStepper.jl `_evaluate_pade`
#   scale the returned numerator by (1 + 1e-7) — a 1e-7 relative perturbation of
#   every stored Padé value.  Result: the CRic.1/CRic.2 tight before-pole VALUE
#   assertions (rtol 1e-11/1e-12) went RED, as did the CRic.1 pos-axis 1e-11
#   pins; the looser past-pole pins (rtol 1e-7) and the extract_poles pins
#   (denominator-root reads, independent of the numerator scale) stayed GREEN —
#   exactly the coverage split the CRic.3/4 note documents.  Restored
#   src/PadeStepper.jl byte-for-byte (`git diff src/` empty); suite GREEN.
#
# Run standalone with:
#   julia --project=. test/corpus_riccati_special_test.jl
# ============================================================================
