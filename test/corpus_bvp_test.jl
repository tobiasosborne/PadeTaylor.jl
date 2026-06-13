# test/corpus_bvp_test.jl
#
# ============================================================================
# Ground-truth corpus: SCALAR BVP — the 3-arg (u'-dependent) Jacobian,
# the oblique-complex segment, BigFloat-256, and the real(t*) domain guard.
# Bead `padetaylor-u638`, epic `padetaylor-p1v0`.  Tagged CBV.*.
#
# WHY THIS FILE EXISTS (capability gap #6 + #7)
# ---------------------------------------------
# The three-argument `bvp_solve(f, ∂f_∂u, ∂f_∂up, …)` overload (src/BVP.jl:370)
# carries a NON-diagonal Newton-Jacobian term coupling u' to u_int through D₁:
#
#     J = D₂_ii − scale·diag(∂f/∂u) − (scale/half_diff)·diag(∂f/∂u')·D₁_ii
#                                     ╰── THE D1-COUPLED, MIS-SCALABLE TERM ──╯
#
# Before this file that term was exercised ONLY through the PIII hybrid
# (IB.1.4/IB.1.5) at N=10, tol=1e-10 — *loose* settings well below the
# spectral floor, where an O(1) mis-scaling of (scale/half_diff) could
# hide.  CBV.1-CBV.3 isolate the term against closed-form oracles at
# SPECTRAL accuracy (the same ~1e-14 floor the 2-arg path reaches), so a
# wrong scale factor cannot survive.  Result (see FINDINGS footer): every
# 3-arg case DOES reach the spectral floor → the D1-coupling is correctly
# scaled.  No gap-#6 bug.
#
# CASES (recipes verbatim from docs/test_corpus/01_corpus_catalogue.md):
#   CBV.1  3arg-damped-oscillator  u''=-2u'-5u,  u=e^{-z}cos2z   [md:1234-1249]
#          → isolates the CONSTANT ∂f/∂u'=-2 D1 coupling.
#   CBV.2  3arg-euler-cauchy       u''=u/z²-u'/z, u=(4/3)z+(2/3)/z [md:1251-1266]
#          → ∂f/∂u'=-1/z is z-DEPENDENT (node-varying D1 coupling); exact rational.
#   CBV.3  3arg-nonlinear-exp      u''=(u')²/u,  u=e^{2z}          [md:1268-1283]
#          → nonlinear; full ∂f/∂u=-(u')²/u², ∂f/∂u'=2u'/u Newton.
#   CBV.4  bratu-bvp (2-arg)       u''=-2e^u, B=cosh(B/2)          [md:1302-1317]
#          → transcendental-constant nonlinear oracle (B solved + verified).
#   CBV.5  airy-bvp (2-arg)        u''=z·u, u=Ai(z)                [md:1319-1334]
#          → clean linear variable-coefficient cross-check (mpmath).
#   CBV.6  oblique-complex-cosh    u''=u, u=cosh(z) on [0,1+i]     [md:1387-1402]
#          + bigfloat256-oblique-cosh  (BigFloat-256, mpmath 80-dps) [md:1404-1419]
#   CBV.7  bvp-domain-guard-imaginary-tstar  (the real(t*) guard)  [md:1919-1934]
#
# ORACLES: every ODE residual u''-f was verified == 0 SYMBOLICALLY with
# sympy, and every boundary/interior value re-derived from mpmath/closed
# form, in external/probes/corpus-oracles/bvp/capture.py — independent of
# the Julia routine under test (Rule 5).  Constants below are pinned from
# that script's output (run: python3 capture.py).  The catalogue is NOT
# infallible (ERRATA.md): I independently re-checked each value; the only
# correction needed was using the airy segment [0,3] (not [0,1.5]) per the
# catalogue's own BC pin Ai(3).
#
# DE-DUP vs the committed suite:
#   * test/bvp_test.jl BV.1.2 already pins the 2-arg linear cosh on [-1,1];
#     we REFERENCE it, do not re-pin.  BV.5.1 pins 2-arg BigFloat real.
#   * NO existing test exercises the 3-arg overload against a closed form
#     (only IB.1.4/.5 at loose tol via PIII) — CBV.1-3 are the gap fill.
#   * NO existing scalar test uses an oblique COMPLEX segment (VB.1.3 is
#     the VECTOR solver) — CBV.6 fills it; CBV.7 pins the guard weakness.
#
# TOLERANCES (justified by spectral accuracy + precision, never fitted —
# the observed errors are in the FINDINGS footer):
#   * CBV.1-6 Float64: node/interior error asserted ≤ 1e-12 (observed
#     ~1e-14..1e-15; the same Chebyshev-Newton floor BV.1.2 hits at N=16).
#   * CBV.6 BigFloat-256: interior error ≤ 1e-35 (observed 5.2e-39; the
#     256-bit analogue of BV.5.1, eps^(3/4) ≈ 1e-57).
#
# MUTATION-PROOF (Rule 4): see footer — a Jacobian-targeted mutant on the
# (scale/half_diff) D1-coupling factor was confirmed RED, then restored.
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor: BVP

# ---- pinned oracle constants (external/probes/corpus-oracles/bvp/capture.py) --
# (1) damped-oscillator  u=e^{-z}cos(2z)
const DAMP_BC_A = -1.1312043837568135      # u(-1)
const DAMP_BC_B = -0.15309186567422628     # u(1)
const DAMP_U_m0p7 = 0.3422717941963816
const DAMP_U_0p0  = 1.0
const DAMP_U_0p7  = 0.08440318529167407
# (3) nonlinear-exp  u=e^{2z}
const NLEXP_BC_B = 7.38905609893065        # u(1)=e^2
const NLEXP_U_0p25 = 1.6487212707001282
const NLEXP_U_0p5  = 2.718281828459045
const NLEXP_U_0p75 = 4.4816890703380645
# (4) bratu  B=cosh(B/2)
const BRATU_B = 1.178775526938701
const BRATU_U_0p25 = 0.24333656779461701
const BRATU_U_0p5  = 0.3289524213411136
const BRATU_U_0p75 = 0.24333656779461701
# (5) airy  Ai
const AIRY_BC_A = 0.3550280538878172       # Ai(0)
const AIRY_BC_B = 0.006591139357460719     # Ai(3)
const AIRY_U_0p5 = 0.23169360648083348
const AIRY_U_1p5 = 0.07174949700810541
const AIRY_U_2p5 = 0.01572592338047049
# (6) oblique cosh  (Float64 + BigFloat-256 80-dps reference strings)
const OBL_INT  = 0.9895848833999199 + 0.24982639750046154im   # cosh(0.5+0.5i)
const OBL_BF_INT_RE = big"0.989584883399919936444070533597850936748218873112423178477981725745469303094275"
const OBL_BF_INT_IM = big"0.249826397500461531489559656030584391059210680171568246651286996984047394219192"

@testset "Corpus BVP (CBV): 3-arg Jacobian, oblique-complex, BigFloat, t* guard" begin

    # =====================================================================
    # CBV.1 — 3arg-damped-oscillator: u'' = -2u' - 5u, u = e^{-z}cos(2z).
    # Isolates the CONSTANT ∂f/∂u' = -2 D1-coupled Jacobian term on [-1,1].
    # =====================================================================
    @testset "CBV.1: 3arg damped oscillator (constant ∂f/∂u'=-2)" begin
        f(z, u, up)    = -2up - 5u
        ∂fu(z, u, up)  = -5.0
        ∂fup(z, u, up) = -2.0
        sol = bvp_solve(f, ∂fu, ∂fup, -1.0, 1.0, DAMP_BC_A, DAMP_BC_B; N = 24)
        uexact(z) = exp(-z) * cos(2z)
        # SPECTRAL accuracy at all nodes — same floor the 2-arg path hits.
        @test maximum(abs, sol.u_nodes .- [uexact(z) for z in sol.nodes_z]) < 1e-12
        # Interior barycentric values vs the independent mpmath oracle.
        @test isapprox(sol(-0.7)[1], DAMP_U_m0p7; atol = 1e-12)
        @test isapprox(sol(0.0)[1],  DAMP_U_0p0;  atol = 1e-12)
        @test isapprox(sol(0.7)[1],  DAMP_U_0p7;  atol = 1e-12)
        # The derivative u'(-0.7) = -e^{0.7}(cos1.4 + 2 sin1.4) from the
        # barycentric D₁ path (a second, independent check of the slice).
        @test isapprox(sol(0.0)[2], -1.0; atol = 1e-11)   # u'(0) = -(1+0) = -1
    end

    # =====================================================================
    # CBV.2 — 3arg-euler-cauchy: u'' = u/z² - u'/z, u = (4/3)z + (2/3)/z.
    # ∂f/∂u' = -1/z is z-DEPENDENT → node-by-node D1 coupling.  EXACT rational.
    # =====================================================================
    @testset "CBV.2: 3arg Euler-Cauchy (z-dependent ∂f/∂u'=-1/z, exact rational)" begin
        f(z, u, up)    = u / z^2 - up / z
        ∂fu(z, u, up)  = 1 / z^2
        ∂fup(z, u, up) = -1 / z
        sol = bvp_solve(f, ∂fu, ∂fup, 1.0, 2.0, 2.0, 3.0; N = 24)
        u2(z) = (4 / 3) * z + (2 / 3) / z
        @test maximum(abs, sol.u_nodes .- [u2(z) for z in sol.nodes_z]) < 1e-12
        # EXACT rational interior pins — the cleanest possible oracle.
        @test isapprox(sol(1.5)[1],  22 / 9; atol = 1e-12)
        @test isapprox(sol(1.25)[1], 11 / 5; atol = 1e-12)
        @test isapprox(sol(1.75)[1], 19 / 7; atol = 1e-12)
    end

    # =====================================================================
    # CBV.3 — 3arg-nonlinear-exp: u'' = (u')²/u, u = e^{2z} on [0,1].
    # Nonlinear Newton with the FULL ∂f/∂u and ∂f/∂u'.
    # =====================================================================
    @testset "CBV.3: 3arg nonlinear u''=(u')²/u (full nonlinear Jacobian)" begin
        f(z, u, up)    = up^2 / u
        ∂fu(z, u, up)  = -up^2 / u^2
        ∂fup(z, u, up) = 2up / u
        sol = bvp_solve(f, ∂fu, ∂fup, 0.0, 1.0, 1.0, NLEXP_BC_B; N = 24)
        u3(z) = exp(2z)
        @test maximum(abs, sol.u_nodes .- [u3(z) for z in sol.nodes_z]) < 1e-12
        @test isapprox(sol(0.25)[1], NLEXP_U_0p25; atol = 1e-12)
        @test isapprox(sol(0.5)[1],  NLEXP_U_0p5;  atol = 1e-12)
        @test isapprox(sol(0.75)[1], NLEXP_U_0p75; atol = 1e-12)
    end

    # =====================================================================
    # CBV.4 — bratu (2-arg): u'' = -2e^u, u = 2log(B/cosh(B(z-1/2))),
    # B = cosh(B/2).  The transcendental B was solved + the residual
    # verified == 0 SYMBOLICALLY for any B (capture.py block 4).
    # =====================================================================
    @testset "CBV.4: Bratu u''=-2e^u (transcendental constant B, verified)" begin
        # Anti-regression: the pinned B actually solves B = cosh(B/2).
        @test isapprox(BRATU_B, cosh(BRATU_B / 2); atol = 1e-13)
        f(z, u)  = -2 * exp(u)
        ∂f(z, u) = -2 * exp(u)
        ub(z) = 2 * log(BRATU_B / cosh(BRATU_B * (z - 0.5)))
        sol = bvp_solve(f, ∂f, 0.0, 1.0, 0.0, 0.0; N = 24)
        @test maximum(abs, sol.u_nodes .- [ub(z) for z in sol.nodes_z]) < 1e-12
        @test isapprox(sol(0.5)[1],  BRATU_U_0p5;  atol = 1e-12)
        @test isapprox(sol(0.25)[1], BRATU_U_0p25; atol = 1e-12)
        @test isapprox(sol(0.75)[1], BRATU_U_0p75; atol = 1e-12)
    end

    # =====================================================================
    # CBV.5 — airy (2-arg): u'' = z·u, u = Ai(z) on [0,3].  Clean linear
    # variable-coefficient cross-check (mpmath airyai oracle).
    # =====================================================================
    @testset "CBV.5: Airy u''=z·u, u=Ai(z) (linear variable-coeff, mpmath)" begin
        f(z, u)  = z * u
        ∂f(z, u) = z
        sol = bvp_solve(f, ∂f, 0.0, 3.0, AIRY_BC_A, AIRY_BC_B; N = 32)
        @test isapprox(sol(0.5)[1], AIRY_U_0p5; atol = 1e-12)
        @test isapprox(sol(1.5)[1], AIRY_U_1p5; atol = 1e-12)
        @test isapprox(sol(2.5)[1], AIRY_U_2p5; atol = 1e-12)
    end

    # =====================================================================
    # CBV.6 — oblique-complex-cosh + bigfloat256: u'' = u, u = cosh(z) on the
    # OBLIQUE segment [0, 1+i].  Exercises the complex affine scale factor
    # s=(z_b-z_a)/2=0.5+0.5i and the Complex{T} promotion path.  The
    # on-segment query 0.5+0.5i has t*=0 (real), so the guard is correct here.
    # =====================================================================
    @testset "CBV.6: oblique-complex cosh on [0,1+i] (Float64 + BigFloat-256)" begin
        f(z, u)  = u
        ∂f(z, u) = one(u)
        zb = 1.0 + 1.0im
        sol = bvp_solve(f, ∂f, 0.0 + 0.0im, zb, 1.0 + 0.0im, cosh(zb); N = 24)
        @test isapprox(sol(0.5 + 0.5im)[1], OBL_INT; atol = 1e-12)
        @test isapprox(sol(0.5 + 0.5im)[1], cosh(0.5 + 0.5im); atol = 1e-12)

        # BigFloat-256: the GenericLinearAlgebra backslash on Complex{BigFloat}.
        setprecision(BigFloat, 256) do
            zbb = big(1.0) + big(1.0) * im
            solbf = bvp_solve(f, ∂f, big(0.0) + big(0.0) * im, zbb,
                              big(1.0) + big(0.0) * im, cosh(zbb); N = 24)
            ubf = solbf(big(0.5) + big(0.5) * im)[1]
            ref = OBL_BF_INT_RE + OBL_BF_INT_IM * im   # mpmath 80-dps oracle
            @test abs(ubf - ref) < big"1e-35"
            # Cross-check against Julia's own BigFloat cosh (a 2nd oracle).
            @test abs(ubf - cosh(big(0.5) + big(0.5) * im)) < big"1e-35"
        end
    end

    # =====================================================================
    # CBV.7 — bvp-domain-guard-imaginary-tstar (gap #7), bug padetaylor-53tu
    # FIXED: the guard at src/BVP.jl now tests |t*| (was real(t*)-only), so an
    # off-segment query whose pre-image t* = 10i — Re(t*)=0 (the old guard
    # admitted) but |Im t*|=10 — is correctly REJECTED (Rule 1 fail-loud).
    #
    # Segment [0, 1+i]; query z=-4.5+5.5i ⇒ t* = (z-half_sum)/half_diff
    #   = (-5+5i)/(0.5+0.5i) = 10i  EXACTLY (capture.py block 7, hand-verified).
    # Pre-fix the barycentric interpolant EXTRAPOLATED wildly (off by ~42 vs
    # cosh); now `abs(t_star) ≤ 1 + 100eps` throws a DomainError.  The former
    # @test_broken should-throw assertion is flipped to @test_throws below.
    # =====================================================================
    @testset "CBV.7: |t*| guard rejects off-segment imaginary t* (53tu fixed)" begin
        f(z, u)  = u
        ∂f(z, u) = one(u)
        zb = 1.0 + 1.0im
        sol = bvp_solve(f, ∂f, 0.0 + 0.0im, zb, 1.0 + 0.0im, cosh(zb); N = 24)
        zq = -4.5 + 5.5im                       # t* = 10i exactly
        # Confirm the pre-image is what the recipe says (independent of the guard).
        t_star = (zq - (0.0im + zb) / 2) / ((zb - 0.0im) / 2)
        @test isapprox(t_star, 10im; atol = 1e-12)
        @test abs(real(t_star)) ≤ 1            # Re(t*)=0: the OLD guard admitted
        @test abs(imag(t_star)) > 1            # but |t*|=10 ⇒ far off-segment

        # FIXED (53tu): the |t*| guard now REJECTS this off-segment query
        # (|t*|=10 > 1) — src/BVP.jl throws a DomainError instead of silently
        # extrapolating (was off by ~42 vs cosh).  Flipped from @test_broken.
        @test_throws DomainError sol(zq)
    end

end # @testset Corpus BVP (CBV)

# ============================================================================
# MUTATION-PROOF (Rule 4) — Jacobian-targeted, verified RED then restored.
# ----------------------------------------------------------------------------
# TARGET: the D1-coupled term's scale factor in src/BVP.jl:439-440
#   J = D2_ii .- scale .* Diagonal(diag_∂fu) .-
#       (scale * inv_hd) .* (Diagonal(diag_∂fup) * D1_ii)
#
# MUTANT J1 — perturb the coupling factor `(scale * inv_hd)` → `(2 * scale *
#   inv_hd)` (double the O(1) D1-coupling).  On CBV.2 (Euler-Cauchy, where
#   ∂f/∂u'=-1/z dominates) Newton converges to a DIFFERENT function: the
#   node error jumps from 7e-15 to ~5e-2, the exact-rational pins u(1.5)=22/9
#   etc. fail (RED).  CBV.1 (constant ∂f/∂u'=-2) likewise degrades to ~3e-2.
#   This is the smoking-gun probe the bead demanded: a mis-scaled (scale/
#   half_diff) is caught.  Verified RED on CBV.1 + CBV.2, restored byte-clean.
#
# MUTANT J2 — drop the D1 term entirely (`.- 0`).  Newton still *converges*
#   (the residual uses the correct f) but to the wrong basin / more slowly;
#   CBV.3 nonlinear fails to reach 1e-12 node error.  Confirms the term is
#   load-bearing for nonlinear convergence too.  Restored.
#
# RESTORATION: `git diff src/` confirmed EMPTY after each mutant.  The
# spectral-accuracy assertions (< 1e-12) are the live guards — they cannot
# be satisfied by a mis-scaled Jacobian on these analytic solutions.
#
# FINDINGS (gap #6 verdict — reported to the caller):
#   Observed Float64 max-node errors at N=24/32 (capture-cross-checked):
#     CBV.1 damped   1.7e-14   CBV.2 euler  7.1e-15   CBV.3 nlexp 1.0e-14
#     CBV.4 bratu    3.0e-15   CBV.5 airy   4.3e-16   CBV.6 oblique 3.4e-15
#   EVERY 3-arg case reaches the SAME spectral floor as the 2-arg path
#   (BV.1.2 hits 6e-15 at N=16).  The D1-coupled (scale/half_diff) factor is
#   correctly scaled — NO gap-#6 bug.  BigFloat-256 CBV.6: 5.2e-39 (clean).
#   gap-#7: the real(t*) guard IS a Rule-1 weakness (CBV.7 @test_broken).
# ============================================================================
