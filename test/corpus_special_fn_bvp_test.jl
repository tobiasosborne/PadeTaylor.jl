# test/corpus_special_fn_bvp_test.jl
#
# ============================================================================
# WHAT THIS FILE PINS
# -------------------
# Special-function BVP oracles for the Chebyshev-Newton spectral solver
# (src/BVP.jl) — a class absent from the corpus before this file.  Three named
# special functions, each solving a 2nd-order linear ODE on a complex line
# segment, pinned against an INDEPENDENT oracle (sympy / scipy / mpmath).
# Bead `padetaylor-evlo`, epic `padetaylor-25og`.  Tagged CBvx.2/CBvx.3.
#
#   CBvx.3  Kummer M(-2, 5/2, z)  (3-arg, EXACT degree-2 POLYNOMIAL)
#           z u'' + (b-z) u' - a u = 0, a=-2, b=5/2.  a a negative integer ⇒ the
#           Kummer series TERMINATES: M = 4z²/35 - 4z/5 + 1.  Singular at z=0 ⇒
#           segment [1/2, 2] (away from origin).  3-arg recast
#           u'' = ((z-5/2)u' - 2u)/z exercises the D1-coupled, NON-diagonal
#           Newton-Jacobian term (gap #6) with a z-DEPENDENT ∂f/∂u' = 1 - 5/(2z).
#           The solution lies in the collocation space (deg 2 ≪ N=24) ⇒ NO
#           truncation error ⇒ exact-rational interior pins at machine precision.
#
#   CBvx.2  Mathieu ce_2  (2-arg, externally-PINNED eigenvalue — a new class)
#           u'' + (a - 2q cos2z) u = 0, q=1, a = a_2(1) = 4.371300982735086.
#           The eigenvalue a (scipy.special.mathieu_a) is what makes the even
#           pi-periodic Mathieu function ce_2 the solution.  2-arg recast
#           u'' = -(a - 2q cos2z) u on [0, pi/2].  Periodic variable-coefficient
#           BVP — the first eigenvalue-pinned member of the corpus.
#
#   CBvx.3-Whittaker  M_{2,3/2}(z)  (2-arg, mpmath dps=50)
#           u'' + (-1/4 + k/z + (1/4-m²)/z²) u = 0, k=2, m=3/2, on [1,4].
#           Clean transcendental cross-check via mpmath.whitm.
#
# GROUND TRUTH (Rule 5 — independent of the routine under test)
# -------------------------------------------------------------
# external/probes/corpus-oracles/special-fn-bvp/capture.py (python3):
#   * sympy verifies each ODE recast residual u''-f IDENTICALLY ZERO as a
#     symbolic rational function (a GATE refusing to emit otherwise) — DLMF 13.x
#     Kummer/Whittaker, DLMF 28.x Mathieu — and cross-checks the Kummer
#     polynomial against sympy.hyperexpand(1F1(-2;5/2;z));
#   * scipy.special supplies the Mathieu eigenvalue a_2(1) AND ce_2 values;
#   * mpmath dps=50 supplies the Whittaker M values.
# The Kummer values are EXACT RATIONALS pinned inline; the Mathieu / Whittaker
# transcendental pins live in test/_oracle_corpus_special_fn_bvp.jl.
#
# TOLERANCES (justified, never fitted — observed errors in the FINDINGS footer)
# ----------------------------------------------------------------------------
#   * CBvx.3 Kummer  : exact-polynomial spectral floor ⇒ 1e-13 (observed ~1e-15).
#   * CBvx.2 Mathieu : scipy is DOUBLE precision ⇒ 1e-10 (the oracle's own
#                      ~1e-16 floor + the BVP truncation error at N=24).
#   * CBvx.3-Whittaker: mpmath 45-digit oracle, Float64 solver ⇒ 1e-12.
#
# REFERENCES
# ----------
#   * docs/test_corpus/02_corpus_extension_plan.md:309-327 (Family H).
#   * DLMF 13.2 (Kummer), 13.14 (Whittaker), 28.2/28.4 (Mathieu).
#   * src/BVP.jl:370-467 (3-arg overload, D1-coupled Jacobian), :239-345 (2-arg).
#
# MUTATION-PROOF (Rule 4): see footer — an oracle-pin flip and an impl-side
# Jacobian/residual perturbation in src/BVP.jl were each confirmed RED, then
# restored byte-for-byte.
# ============================================================================

using Test
using PadeTaylor
include(joinpath(@__DIR__, "_oracle_corpus_special_fn_bvp.jl"))

@testset "Corpus: special-function BVPs (CBvx)" begin

    # =====================================================================
    # CBvx.3 — Kummer M(-2, 5/2, z) = 4z²/35 - 4z/5 + 1  (3-arg, EXACT poly).
    # z u'' + (b-z)u' - a u = 0 ⇒ u'' = ((z-b)u' + a u)/z = ((z-5/2)u' - 2u)/z.
    # Segment [1/2, 2] (away from the z=0 singularity).  Exact-rational pins.
    # =====================================================================
    @testset "CBvx.3: Kummer M(-2,5/2,z) (3-arg, exact degree-2 polynomial)" begin
        # 3-arg RHS u''=f(z,u,u') and its analytic partials (capture.py).
        f(z, u, up)    = ((z - 5 / 2) * up - 2 * u) / z   # = ((z-b)u' + a u)/z
        ∂fu(z, u, up)  = -2 / z                           # a/z      = -2/z
        ∂fup(z, u, up) = 1 - 5 / (2z)                     # (z-b)/z  = 1 - 5/(2z)
        Mexact(z) = 4z^2 / 35 - 4z / 5 + 1
        # BC: u(1/2) = 22/35, u(2) = -1/7 (exact rationals).
        sol = bvp_solve(f, ∂fu, ∂fup, 1 / 2, 2.0, 22 / 35, -1 / 7; N = 24)
        # Exact polynomial ⇒ collocation reproduces it to the FP floor at nodes.
        @test maximum(abs, sol.u_nodes .- [Mexact(z) for z in sol.nodes_z]) < 1e-13
        # EXACT rational interior pins — the cleanest possible oracle.
        @test isapprox(sol(1.0)[1],   11 / 35; atol = 1e-13)   # u(1)   = 11/35
        @test isapprox(sol(3 / 2)[1],  2 / 35; atol = 1e-13)   # u(3/2) = 2/35
        # Derivative cross-check: M'(z) = 8z/35 - 4/5 ⇒ M'(1) = 8/35 - 4/5 = -20/35.
        @test isapprox(sol(1.0)[2], 8 / 35 - 4 / 5; atol = 1e-12)
    end

    # =====================================================================
    # CBvx.2 — Mathieu ce_2 (2-arg, externally-pinned eigenvalue a_2(1)).
    # u'' = -(a - 2q cos2z) u, a=4.371300982735086, q=1, on [0, pi/2].
    # =====================================================================
    @testset "CBvx.2: Mathieu ce_2 (2-arg, scipy eigenvalue a_2(1))" begin
        a, q = CBvx2_A, CBvx2_Q
        f(z, u)  = -(a - 2q * cos(2z)) * u
        ∂f(z, u) = -(a - 2q * cos(2z))
        sol = bvp_solve(f, ∂f, 0.0, π / 2, CBvx2_BC_A, CBvx2_BC_B; N = 24)
        # Interior barycentric values vs the independent scipy ce_2 oracle.
        @test isapprox(sol(π / 4)[1], CBvx2_U_pi4; atol = 1e-10)
        @test isapprox(sol(π / 3)[1], CBvx2_U_pi3; atol = 1e-10)
        # Anti-regression on the eigenvalue itself: ce_2(0)=BC_A must hold the BC.
        @test isapprox(sol(0.0)[1], CBvx2_BC_A; atol = 1e-12)
    end

    # =====================================================================
    # CBvx.3-Whittaker — M_{2,3/2}(z) (2-arg, mpmath dps=50 oracle) on [1,4].
    # u'' = -(-1/4 + k/z + (1/4-m²)/z²) u, k=2, m=3/2.
    # =====================================================================
    @testset "CBvx.3-Whittaker: M_{2,3/2}(z) (2-arg, mpmath dps=50)" begin
        k, m = 2.0, 3 / 2
        coeff(z) = -1 / 4 + k / z + (1 / 4 - m^2) / z^2
        f(z, u)  = -coeff(z) * u
        ∂f(z, u) = -coeff(z)
        sol = bvp_solve(f, ∂f, 1.0, 4.0, CBvxW_BC_A, CBvxW_BC_B; N = 24)
        @test isapprox(sol(2.0)[1], CBvxW_U_2;   atol = 1e-12)
        @test isapprox(sol(2.5)[1], CBvxW_U_2p5; atol = 1e-12)
        @test isapprox(sol(3.0)[1], CBvxW_U_3;   atol = 1e-12)
    end

end # @testset Corpus: special-function BVPs (CBvx)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# ----------------------------------------------------------------------------
# M-ORACLE (pin flip).  Replace the Kummer interior pin u(1) = 11/35 with the
#   wrong rational 12/35.  CBvx.3 went RED (1 fail: got 11/35, expected 12/35).
#   Restored to 11/35 → 10/10 GREEN.
#
# M-IMPL (src/BVP.jl perturbation).  Two independent mutants, each restored
#   byte-for-byte (`git status --porcelain src/` empty after each):
#   (1) 3-arg Jacobian D1-coupling, src/BVP.jl:439-440 — perturb the coupling
#       factor `(scale * inv_hd)` → `(2 * scale * inv_hd)` (double the O(1)
#       D1 coupling).  CBvx.3's z-DEPENDENT ∂f/∂u' = 1-5/(2z) makes Newton
#       diverge: it ERRORS (non-convergence in maxiter) on CBvx.3.  The 2-arg
#       cases (CBvx.2, Whittaker) — which do NOT use the 3-arg path — stayed
#       fully GREEN, isolating the failure to the D1-coupled term.  Restored.
#   (2) 2-arg residual scale, src/BVP.jl:310 — perturb `scale` → `1.01*scale`
#       in the residual assembly.  CBvx.2 (Mathieu) ERRORED (Newton converged
#       to a wrong eigenfunction / non-convergence) and CBvx.3-Whittaker FAILED
#       all 3 interior pins; the 3-arg Kummer case (separate residual path)
#       stayed GREEN.  Confirms the 2-arg path guards the Mathieu/Whittaker
#       pins.  Restored.
#
# FINDINGS (reported to the caller):
#   Observed Float64 max-node / interior errors at N=24 (measured):
#     CBvx.3 Kummer    node 1.5e-15, interior 1.3e-15 (exact poly, machine floor)
#     CBvx.2 Mathieu   interior 3.7e-14 vs scipy (well inside the 1e-10 budget)
#     CBvx.3-Whittaker interior 4.0e-15 vs mpmath (well inside 1e-12)
#   The 3-arg D1-coupled Jacobian reaches the EXACT-polynomial floor on the
#   z-dependent Kummer ∂f/∂u' — no gap-#6 mis-scaling.  Whittaker WAS included
#   (it added a clean mpmath-pinned case at no convergence risk).
#
# STANDALONE RUN:
#   julia --project=. test/corpus_special_fn_bvp_test.jl
# ============================================================================
