# test/corpus_orthopoly_bvp_test.jl
#
# ============================================================================
# WHAT THIS FILE PINS
# -------------------
# Ground-truth corpus: orthogonal-polynomial BVPs whose EXACT solution is a
# degree-5 polynomial.  Because the solution lies in the Chebyshev collocation
# space (N = 24 ≫ deg 5), spectral collocation reproduces it to the machine-
# precision floor — there is NO truncation error.  All three are 3-arg
# (u'-dependent) linear BVPs, chosen to isolate the D1-coupled, NON-diagonal
# Newton-Jacobian term  (scale·inv_hd)·diag(∂f/∂u')·D₁  in src/BVP.jl:439-440
# (capability-map gap #6).  Bead `padetaylor-40q1`, epic `padetaylor-25og`.
# Tagged CBvx.1.<case>.
#
# The node-VARYING ∂f/∂u' of the two singular cases —
#     Legendre  ∂f/∂u' = 2z/(1-z²)        Chebyshev  ∂f/∂u' = z/(1-z²)
# — is a STRONGER probe of that coupling than the v1 CBV.1-3 oracles (which
# used a constant or 1/z coefficient): a mis-scaled (scale·inv_hd) factor
# perturbs every interior node differently and cannot be absorbed into a
# uniform rescale, so the exact-rational interior pins below catch it.
#
# GROUND TRUTH (Rule 5 — independent of the routine under test)
# -------------------------------------------------------------
# external/probes/corpus-oracles/orthopoly-bvp/capture.py (python3 + sympy):
#   * for each named polynomial it verifies the ODE residual u''-f is
#     IDENTICALLY ZERO as a symbolic rational function (a gate that refuses to
#     emit otherwise) — DLMF 18.8 Hermite/Legendre/Chebyshev ODEs;
#   * it prints every BC + interior value as an EXACT RATIONAL (no
#     transcendental constants anywhere).  The literals below are transcribed
#     from that run and pinned inline (per BUILD CONTRACT — no _oracle include).
#
# NOTE (erratum vs the bead prose): the Legendre/Chebyshev boundary signs in
# the bead text are transposed.  P₅ and T₅ are ODD, and on the interior
# segments [-4/5,4/5], [-9/10,9/10] one has P₅(4/5) = -2497/6250 < 0 and
# T₅(9/10) = -3951/6250 < 0 (the negative middle terms dominate the leading
# +z⁵ inside ±1).  sympy is authoritative; the pins use the sympy signs.
#
# TOLERANCES (the exact-polynomial spectral floor — never fitted)
# ---------------------------------------------------------------
# The exact solution is a polynomial of degree 5 < N = 24, so the only error
# source is floating-point arithmetic in building D₂, the Newton solve and the
# barycentric evaluation.  All three cases hit the SAME RELATIVE floor
#   ‖u_node − P(z)‖∞ / ‖u‖∞  ≈ 1.6e-14 (Hermite), 8.5e-15 (Legendre),
#   8.1e-15 (Chebyshev) at N = 24 (measured — see FINDINGS footer).
# The metric is RELATIVE, not absolute, because the absolute floor scales with
# the solution magnitude: a spectral-collocation solve is backward-stable, so
# the node error is ~cond(D₂)·eps·‖u‖∞.  Hermite's H₅ reaches ‖u‖∞ ≈ 41 on
# [-1,1], so its ABSOLUTE node error is ~6.7e-13 = 41·1.6e-14 — fully expected,
# NOT a Jacobian bug (the relative error matches the bounded Legendre/Chebyshev
# cases, |P₅|,|T₅|<1).  Asserting an absolute 1e-13 would wrongly flag the
# large-magnitude case.  We therefore assert a RELATIVE node floor < 1e-13
# (comfortably above the observed ≤1.6e-14, far below any O(1) Jacobian mis-
# scaling) and pin interior values with rtol = 1e-12 (magnitude-aware).  These
# are LINEAR ODEs ⇒ Newton converges in 2 steps from the linear-ramp guess.
#
# REFERENCES
#   docs/test_corpus/02_corpus_extension_plan.md  Family H (CBvx.1, gap #6)
#   DLMF 18.8 (orthogonal-polynomial second-order ODEs)
#   src/BVP.jl:430-453 (3-arg Newton + the D1-coupled Jacobian term)
#
# MUTATION-PROOF (Rule 4): see footer — an M-oracle (flipped rational) and the
# load-bearing M-impl (perturbed D1-coupling factor) were each ACTUALLY RUN
# RED, then restored byte-for-byte.  `git diff src/` is empty.
# ============================================================================

using Test
using PadeTaylor

# ---- pinned EXACT-RATIONAL oracles (capture.py, sympy residual==0 verified) --
# Hermite H₅ = 32z⁵-160z³+120z on [-1,1]
const HERM_BC_A   =  8 // 1          # u(-1)
const HERM_BC_B   = -8 // 1          # u(1)
const HERM_U_mhalf = -41 // 1        # u(-1/2)
const HERM_U_qtr   = 881 // 32       # u(1/4)
const HERM_U_half  =  41 // 1        # u(1/2)
# Legendre P₅ = (63z⁵-70z³+15z)/8 on [-4/5,4/5]  (singular at ±1)
const LEG_BC_A   =  2497 // 6250     # u(-4/5)
const LEG_BC_B   = -2497 // 6250     # u(4/5)
const LEG_U_m2_5 = -3383 // 12500    # u(-2/5)
const LEG_U_0    =  0 // 1           # u(0)
const LEG_U_2_5  =  3383 // 12500    # u(2/5)
# Chebyshev T₅ = 16z⁵-20z³+5z on [-9/10,9/10]  (singular at ±1)
const CHEB_BC_A   =  3951 // 6250    # u(-9/10)
const CHEB_BC_B   = -3951 // 6250    # u(9/10)
const CHEB_U_mhalf = -1 // 2         # u(-1/2)
const CHEB_U_qtr   = 61 // 64        # u(1/4)
const CHEB_U_half  =  1 // 2         # u(1/2)

@testset "Corpus: orthogonal-polynomial BVPs (CBvx)" begin

    # =====================================================================
    # CBvx.1.hermite — u'' = 2z·u' - 2n·u, n=5, u = H₅ on [-1,1].
    # REGULAR everywhere.  ∂f/∂u = -2n = -10 (const); ∂f/∂u' = 2z (linear,
    # node-varying) → exercises the D1-coupled Jacobian term.
    # =====================================================================
    @testset "CBvx.1.hermite: u''=2z·u'-10u, u=H₅ (∂f/∂u'=2z)" begin
        n = 5
        f(z, u, up)    = 2z * up - 2n * u
        ∂fu(z, u, up)  = -2.0 * n
        ∂fup(z, u, up) = 2.0 * z
        sol = bvp_solve(f, ∂fu, ∂fup, -1.0, 1.0,
                        float(HERM_BC_A), float(HERM_BC_B); N = 24)
        H5(z) = 32z^5 - 160z^3 + 120z
        unorm = maximum(abs, sol.u_nodes)
        @test maximum(abs, sol.u_nodes .- [H5(z) for z in sol.nodes_z]) / unorm < 1e-13
        @test isapprox(sol(-0.5)[1], float(HERM_U_mhalf); rtol = 1e-12)
        @test isapprox(sol(0.25)[1], float(HERM_U_qtr);   rtol = 1e-12)
        @test isapprox(sol(0.5)[1],  float(HERM_U_half);  rtol = 1e-12)
    end

    # =====================================================================
    # CBvx.1.legendre — u'' = (2z·u' - n(n+1)·u)/(1-z²), n=5, u = P₅.
    # SINGULAR at ±1 → segment strictly inside, [-4/5,4/5].
    # ∂f/∂u = -30/(1-z²);  ∂f/∂u' = 2z/(1-z²) (NODE-VARYING coupling).
    # =====================================================================
    @testset "CBvx.1.legendre: u''=(2z·u'-30u)/(1-z²), u=P₅ (∂f/∂u'=2z/(1-z²))" begin
        n = 5
        f(z, u, up)    = (2z * up - n*(n+1) * u) / (1 - z^2)
        ∂fu(z, u, up)  = -n*(n+1) / (1 - z^2)
        ∂fup(z, u, up) = 2z / (1 - z^2)
        sol = bvp_solve(f, ∂fu, ∂fup, -0.8, 0.8,
                        float(LEG_BC_A), float(LEG_BC_B); N = 24)
        P5(z) = (63z^5 - 70z^3 + 15z) / 8
        unorm = maximum(abs, sol.u_nodes)
        @test maximum(abs, sol.u_nodes .- [P5(z) for z in sol.nodes_z]) / unorm < 1e-13
        @test isapprox(sol(-0.4)[1], float(LEG_U_m2_5); rtol = 1e-12)
        @test isapprox(sol(0.0)[1],  float(LEG_U_0);    atol = 1e-13)  # P₅(0)=0 exactly
        @test isapprox(sol(0.4)[1],  float(LEG_U_2_5);  rtol = 1e-12)
    end

    # =====================================================================
    # CBvx.1.chebyshev — u'' = (z·u' - n²·u)/(1-z²), n=5, u = T₅.
    # SINGULAR at ±1 → segment [-9/10,9/10].  ∂f/∂u = -25/(1-z²);
    # ∂f/∂u' = z/(1-z²) (NODE-VARYING coupling, distinct from Legendre's 2z).
    # =====================================================================
    @testset "CBvx.1.chebyshev: u''=(z·u'-25u)/(1-z²), u=T₅ (∂f/∂u'=z/(1-z²))" begin
        n = 5
        f(z, u, up)    = (z * up - n^2 * u) / (1 - z^2)
        ∂fu(z, u, up)  = -n^2 / (1 - z^2)
        ∂fup(z, u, up) = z / (1 - z^2)
        sol = bvp_solve(f, ∂fu, ∂fup, -0.9, 0.9,
                        float(CHEB_BC_A), float(CHEB_BC_B); N = 24)
        T5(z) = 16z^5 - 20z^3 + 5z
        unorm = maximum(abs, sol.u_nodes)
        @test maximum(abs, sol.u_nodes .- [T5(z) for z in sol.nodes_z]) / unorm < 1e-13
        @test isapprox(sol(-0.5)[1], float(CHEB_U_mhalf); rtol = 1e-12)
        @test isapprox(sol(0.25)[1], float(CHEB_U_qtr);   rtol = 1e-12)
        @test isapprox(sol(0.5)[1],  float(CHEB_U_half);  rtol = 1e-12)
    end

end # @testset Corpus orthogonal-polynomial BVPs (CBvx)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# ----------------------------------------------------------------------------
# M-oracle (pins are load-bearing): flip CHEB_U_qtr 61//64 → 63//64.  The
#   CBvx.1.chebyshev interior assertion `isapprox(sol(0.25)[1], 61/64)` reddens
#   (the routine returns the true 61/64 = 0.953125, not 63/64).  Restored.
#
# M-impl (THE BEAD'S LOAD-BEARING TARGET — the D1-coupled Jacobian factor in
#   src/BVP.jl:439-440):
#       J = D2_ii .- scale .* Diagonal(diag_∂fu) .-
#           (scale * inv_hd) .* (Diagonal(diag_∂fup) * D1_ii)
#   Perturb the coupling factor `(scale * inv_hd)` → `(2 * scale * inv_hd)`
#   (an O(1) doubling of the off-diagonal ∂f/∂u' term).  Because the residual R
#   still uses the CORRECT f, a mis-scaled Jacobian changes only Newton's STEP,
#   not its fixed point — so an O(1) error wrecks the convergence RATE: Newton
#   no longer reaches tol within maxiter and the fail-fast guard (Rule 1) throws
#   `bvp_solve(3-arg): Newton did not converge …` for ALL THREE cases — Hermite
#   (‖Δu‖∞→0.058), and the node-VARYING-∂f/∂u' Legendre + Chebyshev (‖R‖∞→1e-5,
#   3.5e-5).  All three ERROR → RED.  A 1.5× factor likewise diverges (Legendre
#   the closest, ‖Δu‖∞=4.5e-11 just past tol) — confirming the right root but a
#   broken rate.  A near-identity 1.0000001× factor still converges (GREEN):
#   only an O(1) mis-scale is fatal, exactly the gap-#6 failure mode.  Verdict:
#   the (scale·inv_hd) factor is CORRECTLY scaled — NO gap-#6 bug; the mutant is
#   caught (via non-convergence) on the node-varying cases the bead targeted.
#   Restored byte-for-byte; `git diff src/` confirmed EMPTY.
#
# FINDINGS
#   Observed RELATIVE max-node errors at N=24 (capture-cross-checked):
#     Hermite 1.6e-14   Legendre 8.5e-15   Chebyshev 8.1e-15  → all at the same
#   ~1e-14 spectral-collocation floor.  Hermite's ABSOLUTE error is ~6.7e-13
#   because ‖H₅‖∞≈41 on [-1,1]; the bead's 1e-13 ABSOLUTE target was the wrong
#   metric for a large-magnitude solution (this is the only deviation from the
#   bead — a relative floor is the correct, magnitude-aware pin; NOT a bug).
#   Erratum: the bead's Legendre/Chebyshev BC signs are transposed — sympy
#   (authoritative) gives u(4/5)=-2497/6250, u(9/10)=-3951/6250 (odd P₅,T₅).
#
# Standalone run:
#   julia --project=. test/corpus_orthopoly_bvp_test.jl
# ============================================================================
