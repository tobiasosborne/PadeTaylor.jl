# test/corpus_elliptic_lattice_test.jl
#
# ============================================================================
# Ground-truth corpus: NEW ELLIPTIC LATTICES + the Jacobi-TRIPLE shared-Q
# oracle (bead padetaylor-dctk, epic padetaylor-25og; plan Family C,
# docs/test_corpus/02_corpus_extension_plan.md:204-228).
#
# WHAT THIS FILE PINS
# -------------------
# The v1 corpus covers only ONE 2D lattice — the equianharmonic ℘ (60°
# rhombus, CPB.1).  This file adds the second canonical lattice (the SQUARE
# lemniscatic ℘) and the cleanest external oracle for the v0.2 shared-Q
# keystone: a d=3 vector whose three components share ONE identical pole
# lattice.  Three independently-derived cases:
#
#   CEl.1  lemniscatic ℘ (SQUARE lattice)   u''=6u²−1/2, u=℘(z−1;1,0), double
#          pole z=1.  Same scalar idiom as the CPB.1 keystone, but the SQUARE
#          lattice (ω'/ω=1) vs equianharmonic's rhombus.  The discriminating
#          signal: the lattice correction above the bare 1/(z−1)²=400 at
#          z=1.05 is ~1.25e-4 (square) vs ~4.46e-7 (rhombus) — an ~280×
#          contrast at FIXED double-pole order, isolating lattice GEOMETRY.
#   CEl.2  Jacobi cn, dn (scalar 2nd-order)  cn: u''=(2m−1)u−2m u³, IC (1,0);
#          dn: u''=(2−m)u−2u³, IC (1,0).  m=0.36 (k=0.6).  Simple poles at
#          iK'+lattice (OFF the real axis), residues −i/k (cn), −i (dn).
#          Complex bridge UP the imaginary axis through iK'.
#   CEl.3  Jacobi TRIPLE (first-order d=3 vector — the shared-Q oracle)
#          sn'=cn dn, cn'=−sn dn, dn'=−m sn cn, IC (0,1,1).  m=1/2 (self-dual
#          K=K', square lattice — links to CEl.1).  Three components, ONE
#          identical pole lattice (simple poles at iK') — the textbook
#          external oracle for the shared-Q keystone (vs VPN.1.1's DISTINCT
#          per-component poles).  Far-side EXACT algebraic at the half-period
#          past iK': (sn,cn,dn)=(−i/√k, −√((1+k)/k), −√(1+k)), k=√m.
#
# GROUND TRUTH each case embodies (all INDEPENDENTLY verified — Law 1)
# -------------------------------------------------------------------
# external/probes/corpus-oracles/elliptic-lattices/capture.py (sympy 1.14 +
# mpmath 1.3, dps=50) GATES every value before emitting it:
#   * CEl.1: recast residual u''−(6u²−1/2)≡0 (≤1e-38) AND (℘')²=4℘³−g₂℘−g₃;
#     the ℘ is the validated Jacobi-sn recipe (real e={1/2,0,−1/2}, dps-clean).
#   * CEl.2: cn''/dn'' recast residuals (≤1e-38) AND the first-order
#     identities (cn')²=(1−cn²)(1−m+m cn²), (dn')²=(1−dn²)(dn²−1+m).
#   * CEl.3: the d=3 SYSTEM residuals (≤1e-38), the invariants sn²+cn²=1,
#     dn²+m sn²=1 (≤1e-40), AND the far-side algebraic identity (≤1e-40).
# Transcendental pins live in _oracle_corpus_elliptic_lattice.jl; the exact
# algebraic targets (−i/√k etc., the bare 400) are computed in-test.
#
# TOLERANCES (justified by METHOD, never fitted to pass)
# ------------------------------------------------------
#   * CEl.1 past-pole ℘ via the Float64 (15,15)-diagonal Padé over a single
#     h=1.5 window saturates at the same ~3e-10 rational-approx floor as
#     CPB.1 (README §2; problems_test.jl 6.1.6): rtol 1e-9 near the pole
#     (z=1.05), 1e-7 far (z=1.4).  The lattice-correction SIGNALS are pinned
#     to a generous rtol 1e-4 (they differ by 280× — the test only needs to
#     resolve the ORDER of each, not its last digit).
#   * CEl.2 complex single-step bridge up the imaginary axis (h=K'≈2.0 spans
#     a half-period across the pole): ~1e-5, matching CPB.7's h=Kₚ sn bridge.
#   * CEl.3 shared-Q vector bridge through iK' via vector_path_network_solve
#     + extrapolate: ~1e-5 (brief Family C).  The exact-algebraic far-side
#     and the post-pole mpmath pin agree to ~1e-14, well inside 1e-5.
#
# REFERENCES: plan Family C (md:204-228); DLMF 22.x (Jacobi sn/cn/dn),
# 23.x (Weierstrass ℘).  Lattice contrast plan md:210.
#
# MUTATION-PROOF (Rule 4): footer block — one M-oracle (flip a pinned value)
# and the M-impl perturbing each case's load-bearing component (scalar
# stepper/Padé for CEl.1/CEl.2, the vector shared-Q jet for CEl.3) were
# ACTUALLY executed RED and restored byte-for-byte (`git diff src/` empty).
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor.VectorProblems: VectorPadeTaylorProblem
using PadeTaylor.VectorPathNetwork: vector_path_network_solve
using PadeTaylor.VectorPoleField: extract_poles_shared_q

include(joinpath(@__DIR__, "_oracle_corpus_elliptic_lattice.jl"))

@testset "Corpus: elliptic lattices + Jacobi triple (CEl)" begin

    # ----------------------------------------------------------------------
    # CEl.1  lemniscatic ℘ (SQUARE lattice): u''=6u²−1/2, double pole at z=1.
    # Mirrors the CPB.1 keystone idiom — single window z₀=0→1.5, h_max=1.5
    # (pole z=1 lands at t≈0.667, interior) — at the SQUARE lattice instead
    # of the equianharmonic rhombus.
    # ----------------------------------------------------------------------
    @testset "CEl.1  ℘(z−1;1,0): Padé bridges the SQUARE-lattice pole at z=1" begin
        f(z, u, up) = 6 * u^2 - 1 // 2
        prob = PadeTaylorProblem(f, (CEL1_U0, CEL1_UP0), (0.0, 1.5); order = 30)
        sol  = solve_pade(prob; h_max = 1.5)         # single seg, pole at t≈0.667
        u_before, _      = sol(0.5)
        u_past, up_past  = sol(1.05)                 # PAST the double pole
        u_far, _         = sol(1.4)
        @test isapprox(u_before, CEL1_U_0p5;   rtol = 1e-10)
        @test isapprox(u_past,   CEL1_U_1p05;  rtol = 1e-9)
        @test isapprox(up_past,  CEL1_UP_1p05; rtol = 1e-9)
        @test isapprox(u_far,    CEL1_U_1p4;   rtol = 1e-7)

        # THE LATTICE-DISCRIMINATING SIGNAL (plan md:210).  The bare double
        # pole gives 1/(z−1)² = 1/0.05² = 400 exactly; the SQUARE lattice adds
        # a ~1.25e-4 correction (the bridged value minus 400) — ~280× the
        # equianharmonic rhombus's ~4.46e-7.  The reconstructed correction
        # must match the SQUARE-lattice value, NOT the rhombus one — this is
        # the geometric fingerprint the test pins.
        bare = 1 / 0.05^2                            # = 400 exactly
        corr = u_past - bare
        @test isapprox(corr, CEL1_LATTICE_CORR_1p05; rtol = 1e-4)   # square
        # It is DISTINCT from the equianharmonic correction: the lemniscatic
        # signal is ~280× larger, so the square-lattice value is nowhere near
        # the rhombus one — the lattice geometry is resolved, not just the pole.
        @test corr / CEL1_EQUI_CORR_1p05 > 100
        @test !isapprox(corr, CEL1_EQUI_CORR_1p05; rtol = 0.5)
    end

    # ----------------------------------------------------------------------
    # CEl.2  Jacobi cn, dn: simple poles at iK' OFF the real axis ⇒ bridge UP
    # the imaginary axis via path_network_solve + eval_at (the CPB.7 idiom).
    # cn and dn are REAL on the imaginary axis between lattice points and FLIP
    # SIGN across the pole at iK' (PRE>0 at iK'/2, POST<0 at i·3K'/2).
    # ----------------------------------------------------------------------
    @testset "CEl.2  cn(z,0.36): imaginary bridge through the pole at iK'" begin
        m = 0.36
        fcn(z, u, up) = (2m - 1) * u - 2m * u^3
        prob = PadeTaylorProblem(fcn,
                                 (complex(1.0, 0.0), complex(0.0, 0.0)),   # cn(0)=1, cn'(0)=0
                                 (complex(0.0, 0.0), complex(0.0, 2 * CEL2_KP));
                                 order = 30)
        ztgt = complex(0.0, 3 * CEL2_KP / 2)         # PAST pole at iK'
        zpre = complex(0.0, CEL2_KP / 2)             # BEFORE the pole
        sol = path_network_solve(prob, [ztgt]; h = CEL2_KP, order = 30)
        u_pre,  _ = eval_at(sol, zpre)
        u_past, _ = eval_at(sol, ztgt)
        @test isapprox(real(u_pre),  CEL2_CN_PRE;  rtol = 1e-6)
        @test isapprox(imag(u_pre),  0.0;          atol = 1e-6)
        @test isapprox(real(u_past), CEL2_CN_POST; rtol = 1e-5)   # PAST pole, sign-flipped
        @test isapprox(imag(u_past), 0.0;          atol = 1e-5)
        # The bridge really crossed the pole: cn flipped sign (PRE>0, POST<0).
        @test real(u_pre) > 0 && real(u_past) < 0
        @test isapprox(real(u_past), -CEL2_CN_PRE; rtol = 1e-5)    # odd across iK'
    end

    @testset "CEl.2  dn(z,0.36): imaginary bridge through the pole at iK'" begin
        m = 0.36
        fdn(z, u, up) = (2 - m) * u - 2 * u^3
        prob = PadeTaylorProblem(fdn,
                                 (complex(1.0, 0.0), complex(0.0, 0.0)),   # dn(0)=1, dn'(0)=0
                                 (complex(0.0, 0.0), complex(0.0, 2 * CEL2_KP));
                                 order = 30)
        ztgt = complex(0.0, 3 * CEL2_KP / 2)
        zpre = complex(0.0, CEL2_KP / 2)
        sol = path_network_solve(prob, [ztgt]; h = CEL2_KP, order = 30)
        u_pre,  _ = eval_at(sol, zpre)
        u_past, _ = eval_at(sol, ztgt)
        @test isapprox(real(u_pre),  CEL2_DN_PRE;  rtol = 1e-6)
        @test isapprox(real(u_past), CEL2_DN_POST; rtol = 1e-5)   # PAST pole
        @test real(u_pre) > 0 && real(u_past) < 0
        @test isapprox(real(u_past), -CEL2_DN_PRE; rtol = 1e-5)
    end

    # ----------------------------------------------------------------------
    # CEl.3  Jacobi TRIPLE — the cleanest shared-Q oracle.  d=3 first-order
    # vector sn'=cn dn, cn'=−sn dn, dn'=−m sn cn, IC (0,1,1), m=1/2.  All
    # THREE components share ONE identical pole lattice (simple poles at iK');
    # bridge UP the imaginary axis through iK' via vector_path_network_solve +
    # extrapolate (the CHop.2 idiom), then pin the far-side EXACT algebraic
    # values AND confirm extract_poles_shared_q places the shared pole at iK'.
    # ----------------------------------------------------------------------
    @testset "CEl.3  Jacobi triple: ONE shared-Q lattice, far-side algebraic" begin
        m = 0.5
        f = (z, y) -> [y[2] * y[3], -y[1] * y[3], -m * y[1] * y[2]]
        y0 = ComplexF64[0.0, 1.0, 1.0]               # (sn,cn,dn)(0) = (0,1,1)
        ztgt = complex(0.0, 3 * CEL3_KP / 2)         # PAST pole at iK'
        prob = VectorPadeTaylorProblem(f, y0,
                                       (complex(0.0, 0.0), ztgt + 1.0); order = 30)
        sol = vector_path_network_solve(prob, ComplexF64[ztgt]; order = 30,
                  h = 0.6, step_policy = :max_q_root, adaptive = true,
                  on_target_failure = :throw, fine_grid = ComplexF64[ztgt],
                  extrapolate = true)
        @test length(sol.visited_z) ≥ 1

        # Far-side values: (sn,cn,dn) PAST the shared iK' pole.  The EXACT
        # algebraic targets (k=√m): sn=−i/√k, cn=−√((1+k)/k), dn=−√(1+k).
        k = sqrt(m)
        sn_alg = complex(0.0, -1 / sqrt(k))
        cn_alg = complex(-sqrt((1 + k) / k), 0.0)
        dn_alg = complex(-sqrt(1 + k), 0.0)
        yq = sol.grid_y[1]
        @test all(isfinite, yq)
        @test length(yq) == 3                        # the d=3 shared-Q state
        # Pin each component to its EXACT algebraic far-side value.
        @test isapprox(yq[1], sn_alg; atol = 1e-5)
        @test isapprox(yq[2], cn_alg; atol = 1e-5)
        @test isapprox(yq[3], dn_alg; atol = 1e-5)
        # Cross-check against the mpmath pins (same value, independent route).
        @test isapprox(imag(yq[1]), CEL3_SN_POST_IM; atol = 1e-5)
        @test isapprox(real(yq[2]), CEL3_CN_POST;    atol = 1e-5)
        @test isapprox(real(yq[3]), CEL3_DN_POST;    atol = 1e-5)
        # The two invariants hold PAST the pole (the shared-lattice signature):
        # sn²+cn²=1 and dn²+m·sn²=1 even on the far side of iK'.
        @test isapprox(yq[1]^2 + yq[2]^2, complex(1.0, 0.0); atol = 1e-4)
        @test isapprox(yq[3]^2 + m * yq[1]^2, complex(1.0, 0.0); atol = 1e-4)

        # The SHARED pole: extract_poles_shared_q roots the ONE shared Q and
        # places the simple pole at z=iK' (the half-period the bridge crossed).
        # min_support=2 (cross-node) ⇒ a genuine shared pole, not a per-node
        # artefact.  This is the keystone the file delivers: a single Q whose
        # root IS the lattice pole all three components share.
        poles = extract_poles_shared_q(sol; radius_t = 5.0, min_support = 2)
        @test !isempty(poles)
        @test minimum(abs(p - complex(0.0, CEL3_KP)) for p in poles) < 1e-2
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# `git diff src/` confirmed EMPTY after each restore.
#
# M-oracle (CEl.1): flip the pinned past-pole value `CEL1_U_1p05` in
#   _oracle_corpus_elliptic_lattice.jl from `400.000125…` to `399.999875…`
#   (subtract 2× the lattice correction — i.e. assert the WRONG-SIGN lattice
#   geometry).  Result: CEl.1's `isapprox(u_past, CEL1_U_1p05; rtol=1e-9)`
#   AND the `isapprox(corr, CEL1_LATTICE_CORR_1p05; …)` signal assertion went
#   RED, and ONLY in CEl.1; CEl.2/CEl.3 stayed GREEN.  This proves the test
#   resolves the SQUARE-lattice correction to its sign — the discriminating
#   signal is genuinely pinned.  Restored the literal; suite GREEN.
#
# M-impl-A (CEl.1 + CEl.2, the scalar stored-Padé evaluation): in
#   src/PadeStepper.jl `_evaluate_pade`, multiply the numerator constant term
#   by 1.0000001 (a 1e-7 corruption of every evaluated scalar Padé value).
#   Result: RED at CEl.1 (u_before, u_past, u_far, and the corr signal) and
#   CEl.2 cn/dn (all the real-axis value pins) — every assertion routed
#   through the scalar stepper Padé caught it.  CEl.3 stayed GREEN: the VECTOR
#   path-network evaluates through its own routine, not this scalar one — so
#   the mutant cleanly isolates the scalar value path.  Restored byte-for-byte.
#
# M-impl-B (CEl.3, the shared-Q vector jet convolution — the KEYSTONE): in
#   src/VectorCoefficients.jl, change the integration relation
#   `y[i][k] = f_t[i][k-1] / k` → `/ (k+1)` (the off-by-one CVE.4 / CHop.2
#   mutate).  Result: RED at CEl.3 — every far-side algebraic component pin
#   (sn/cn/dn), both invariants past the pole, AND the shared-pole LOCATION
#   (a corrupted jet moves the shared-Q root off iK').  CEl.1/CEl.2 (scalar
#   `taylor_coefficients_2nd`, a different jet builder) stayed GREEN.  This
#   confirms the d=3 shared-Q keystone is pinned by the INDEPENDENT oracle
#   (exact algebraic / mpmath), not self-consistently.  Restored byte-for-byte.
#
# Run standalone with:
#   julia --project=. test/corpus_elliptic_lattice_test.jl
# ============================================================================
