# test/corpus_oblique_guard_test.jl
#
# ============================================================================
# WHAT THIS FILE PINS — CBvx.4 oblique-complex BVP/VectorBVP domain guard.
# Plan docs/test_corpus/02_corpus_extension_plan.md Family H (CBvx.4).
# Bead padetaylor-p0yv (epic padetaylor-25og); CATCHES bug padetaylor-53tu.
#
# THE BUG (padetaylor-53tu).  The barycentric-callable domain guard checks
# only `real(t_star)` — src/BVP.jl:491 (BVPSolution) and the identical guard
# src/VectorBVP.jl:318 (VectorBVPSolution):
#
#     real(t_star) ≤ 1 + 100eps && real(t_star) ≥ -1 - 100eps || throw(...)
#                   ╰──────────── Im(t*) is IGNORED ────────────╯
#
# On an OBLIQUE-complex segment (both Re and Im of z_b-z_a nonzero) a query
# point whose pre-image t* has |Re t*| ≤ 1 but LARGE |Im t*| PASSES the guard,
# and the barycentric interpolant then EXTRAPOLATES silently off-segment — a
# Rule-1 violation (it should fail loud).
#
# DISTINCTION vs CBV.7 (test/corpus_bvp_test.jl).  CBV.7 pins only the
# PERPENDICULAR case: query z=-4.5+5.5i ⇒ t*=10i, Re(t*)=0 (exactly on the
# guard's "centre").  That is the easy-to-spot case.  CBvx.4 is the insidious
# INTERIOR-Re complement: t* = 0.5 + 3i has Re(t*)=0.5 STRICTLY INTERIOR to
# (-1,1) — it looks like a bona-fide in-segment query — yet Im(t*)=3 puts it
# far off the segment.  A `Re`-only guard cannot tell the two apart; an
# honest `|t*|`-based guard rejects both.
#
# GROUND TRUTH.  ODE u''=u, exact u=cosh(z), on z_a=0, z_b=1+i, N=24, 2-arg
# bvp_solve (f(z,u)=u, ∂f/∂u=1).  Vector mirror y'=(y2,y1) ⇒ y=(cosh,sinh)
# on the same oblique segment (B_a=I, B_b=0, g=(cosh0,sinh0)=(1,0); rank 2).
# The exact t* mapping (READ from src/BVP.jl:489):
#     half_sum = half_diff = 0.5+0.5i;  t*(z*) = (z* - half_sum)/half_diff.
# off-segment z* solved from t*=0.5+3i:
#     z_off = half_sum + half_diff·(0.5+3i) = -0.75 + 2.25i  (exact).
# Oracle: mpmath cosh/sinh dps=30, sympy residual==0 gate
#     (external/probes/corpus-oracles/oblique-guard/capture.py;
#      pinned in test/_oracle_corpus_oblique_guard.jl).
#
# TOLERANCES (spectral floor, never fitted):
#   * on-segment value vs cosh/sinh ≤ 1e-12 (the same Chebyshev-Newton floor
#     CBV.6 hits at N=24, observed ~3e-15 — proves the solver is CORRECT on
#     this oblique segment, so the off-segment failure is purely the guard).
#   * off-segment extrapolation error > 0.1 (documents 53tu's wrong answer;
#     observed ≫ 0.1 — a spectral interpolant diverges fast off-segment).
#
# MUTATION-PROOF: see footer — M-oracle (flip a pinned cosh value → reddens)
# and M-impl (perturb the barycentric eval / the guard in src/BVP.jl) were
# ACTUALLY EXECUTED then restored byte-for-byte.
#
# REFERENCES: plan Family H (CBvx.4); bug padetaylor-53tu; CBV.6/CBV.7 in
# test/corpus_bvp_test.jl; src/BVP.jl:491, src/VectorBVP.jl:318.
# ============================================================================

using Test
using PadeTaylor
using LinearAlgebra: rank
include(joinpath(@__DIR__, "_oracle_corpus_oblique_guard.jl"))

# `rejects(f)` — does calling `f()` throw?  The 53tu auto-flip pattern: today
# the off-segment query does NOT throw (silent extrapolation) ⇒ `false` ⇒ the
# @test_broken is "broken"; once 53tu replaces the Re-only guard with a |t*|
# guard the query WILL throw ⇒ `true` ⇒ the marker flips to "unexpectedly
# passes", signalling the fix landed.
rejects(f) = try; f(); false; catch; true; end

@testset "Corpus: oblique-complex BVP guard (CBvx)" begin

    # ====================================================================
    # CBvx.4.1 — SCALAR: u''=u, u=cosh(z) on the oblique segment [0,1+i].
    # ON-segment correctness proves the solver itself is sound here, so the
    # off-segment failure isolates purely to the real(t*)-only guard.
    # ====================================================================
    @testset "CBvx.4.1  scalar oblique cosh: on-segment correct" begin
        f(z, u)  = u
        ∂f(z, u) = one(u)
        sol = bvp_solve(f, ∂f, OG_Z_A, OG_Z_B, OG_BC_A, OG_BC_B; N = 24)

        # On-segment midpoint z_mid=0.5+0.5i has t*=0 (real) — guard correct.
        # (For even N, t*=0 is a node ⇒ this hits the coincident-node guard.)
        u_mid = sol(OG_Z_MID)[1]
        @test isapprox(u_mid, OG_COSH_MID; atol = 1e-12)      # mpmath oracle
        @test isapprox(u_mid, cosh(OG_Z_MID); atol = 1e-12)   # 2nd oracle
        # On-segment OFF-NODE z_qrt=0.65+0.65i has t*=0.3 — this EXERCISES the
        # barycentric interpolant (not the coincident-node fast path), so it is
        # the assertion the M-impl barycentric mutation reddens.
        u_qrt = sol(OG_Z_QRT)[1]
        @test isapprox(u_qrt, OG_COSH_QRT; atol = 1e-12)      # mpmath oracle
        @test isapprox(u_qrt, cosh(OG_Z_QRT); atol = 1e-12)   # 2nd oracle
        # All nodes reach the spectral floor (the CBV.6 floor, ~3e-15).
        @test maximum(abs, sol.u_nodes .- [cosh(z) for z in sol.nodes_z]) < 1e-12
    end

    # ====================================================================
    # CBvx.4.2 — SCALAR 53tu: the INTERIOR-Re(t*) silent extrapolation.
    # z_off=-0.75+2.25i ⇒ t*=0.5+3i: Re(t*)=0.5 INTERIOR (≠ CBV.7's Re=0),
    # Im(t*)=3 large.  The guard passes it; barycentric extrapolates.
    # ====================================================================
    @testset "CBvx.4.2  scalar 53tu: interior-Re t* extrapolates (no throw)" begin
        f(z, u)  = u
        ∂f(z, u) = one(u)
        sol = bvp_solve(f, ∂f, OG_Z_A, OG_Z_B, OG_BC_A, OG_BC_B; N = 24)

        # Confirm the pre-image is the INTERIOR-Re case (independent of guard).
        half_diff = (OG_Z_B - OG_Z_A) / 2
        t_star = (OG_Z_OFF - (OG_Z_A + OG_Z_B) / 2) / half_diff
        @test isapprox(t_star, OG_T_OFF; atol = 1e-12)     # t* = 0.5 + 3i
        @test 0 < real(t_star) < 1                          # STRICTLY interior
        @test abs(imag(t_star)) > 1                         # far off-segment

        # CURRENT (53tu) behaviour: NO throw, and the returned value is a
        # meaningless extrapolation wrong by ≫ 0.1 vs the true cosh(z_off).
        u_off = sol(OG_Z_OFF)[1]
        @test isfinite(real(u_off)) && isfinite(imag(u_off))   # returns finite
        @test abs(u_off - OG_COSH_OFF) > 0.1                   # mpmath true cosh
        @test abs(u_off - cosh(OG_Z_OFF)) > 0.1                # 2nd oracle

        # 53tu AUTO-FLIP MARKER (scalar).  A Rule-1-faithful |t*| guard SHOULD
        # reject the off-segment query.  Today it does NOT (false) ⇒ broken;
        # flips to "unexpectedly passes" when src/BVP.jl:491 uses |t*|.
        @test_broken rejects(() -> sol(OG_Z_OFF))
        # Anti-regression: an ON-segment query must NOT be rejected.
        @test !rejects(() -> sol(OG_Z_MID))
    end

    # ====================================================================
    # CBvx.4.3 — VECTOR MIRROR: y'=(y2,y1), y=(cosh z, sinh z) on [0,1+i].
    # Same guard (src/VectorBVP.jl:318).  B_a=I, B_b=0, g=(1,0) pins both
    # components at z_a (full rank 2); the linear system recovers (cosh,sinh).
    # ====================================================================
    @testset "CBvx.4.3  vector oblique (cosh,sinh): on-segment correct" begin
        fv(z, y) = [y[2], y[1]]                      # y1'=y2=sinh, y2'=y1=cosh
        B_a = [1.0+0im 0; 0 1.0+0im]                 # picks y(z_a)
        B_b = [0.0+0im 0; 0 0.0+0im]
        g   = [OG_BC_A, OG_SINH_A]                   # (cosh0, sinh0) = (1, 0)
        @test rank([B_a B_b]) == 2                    # well-posed
        sol = vector_bvp_solve(fv, OG_Z_A, OG_Z_B, B_a, B_b, g; N = 24)

        y_mid = sol(OG_Z_MID)                                 # node (t*=0)
        @test isapprox(y_mid[1], OG_COSH_MID; atol = 1e-12)   # mpmath cosh
        @test isapprox(y_mid[2], OG_SINH_MID; atol = 1e-12)   # mpmath sinh
        @test isapprox(y_mid, [cosh(OG_Z_MID), sinh(OG_Z_MID)]; atol = 1e-12)
        y_qrt = sol(OG_Z_QRT)                                 # off-node (t*=0.3)
        @test isapprox(y_qrt[1], OG_COSH_QRT; atol = 1e-12)   # barycentric cosh
        @test isapprox(y_qrt[2], OG_SINH_QRT; atol = 1e-12)   # barycentric sinh
        for j in 1:sol.N+1
            z = sol.nodes_z[j]
            @test isapprox(sol.Y_nodes[j], [cosh(z), sinh(z)]; atol = 1e-12)
        end
    end

    # ====================================================================
    # CBvx.4.4 — VECTOR 53tu: same interior-Re(t*) silent extrapolation at
    # the VectorBVP.jl:318 guard.
    # ====================================================================
    @testset "CBvx.4.4  vector 53tu: interior-Re t* extrapolates (no throw)" begin
        fv(z, y) = [y[2], y[1]]
        B_a = [1.0+0im 0; 0 1.0+0im]
        B_b = [0.0+0im 0; 0 0.0+0im]
        g   = [OG_BC_A, OG_SINH_A]
        sol = vector_bvp_solve(fv, OG_Z_A, OG_Z_B, B_a, B_b, g; N = 24)

        # CURRENT (53tu) behaviour: no throw, both components extrapolate.
        y_off = sol(OG_Z_OFF)
        @test all(isfinite, real.(y_off)) && all(isfinite, imag.(y_off))
        @test abs(y_off[1] - OG_COSH_OFF) > 0.1      # cosh component, mpmath
        @test abs(y_off[2] - OG_SINH_OFF) > 0.1      # sinh component, mpmath

        # 53tu AUTO-FLIP MARKER (vector).  Mirrors CBvx.4.2: flips when the
        # src/VectorBVP.jl:318 guard switches to |t*|.
        @test_broken rejects(() -> sol(OG_Z_OFF))
        @test !rejects(() -> sol(OG_Z_MID))          # on-segment NOT rejected
    end

end # @testset Corpus oblique-complex BVP guard (CBvx)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# ----------------------------------------------------------------------------
# M-ORACLE.  Flip OG_COSH_MID's real part 0.98958… → 0.88958… in
#   test/_oracle_corpus_oblique_guard.jl.  Result: CBvx.4.1's mpmath pin
#   `isapprox(u_mid, OG_COSH_MID; atol=1e-12)` AND CBvx.4.3's vector cosh pin
#   both go RED (|Δ|≈0.1 ≫ 1e-12); the cross-checks against Julia's own cosh
#   stay GREEN, proving the corrupted *mpmath* pin is what failed (not the
#   solver).  Restored byte-for-byte (untracked file → backup byte-identical).
#
#   NOTE: z_mid (t*=0) is a Chebyshev NODE for even N, so sol(z_mid) returns
#   via the coincident-node fast path — it does NOT exercise the barycentric
#   interpolant.  That is exactly why CBvx.4.1/4.3 ALSO pin the OFF-NODE point
#   z_qrt (t*=0.3): the barycentric path is the M-impl(a) target below.
#
# M-IMPL.  Two independent perturbations of src/BVP.jl, each restored
#   byte-for-byte (`git status --porcelain src/` empty after each):
#   (a) BARYCENTRIC EVAL — in `_barycentric_eval` change `num += w*u/delta`
#       to `num += 1.001*w*u/delta`.  Result: the OFF-NODE assertions in
#       CBvx.4.1 (`isapprox(u_qrt, OG_COSH_QRT)` + the Julia cross-check) go
#       RED — confirming the barycentric path is load-bearing, not vacuous.
#       CBvx.4.3 (VECTOR) stays GREEN: its interpolant lives in VectorBVP.jl,
#       so the BVP.jl mutant cannot reach it — proof the two modules' eval
#       paths are independent.  Restored.
#   (b) GUARD — in src/BVP.jl:491 replace the two `real(t_star)` checks with a
#       single `abs(t_star) ≤ 1 + 100eps` (the 53tu FIX).  Result: the scalar
#       case CBvx.4.2 FLIPS — the now-correct guard throws a DomainError at the
#       off-segment query, so `rejects(() -> sol(OG_Z_OFF))` would become true
#       (the @test_broken auto-flips to "unexpectedly passes") — exactly the
#       fix-landed signal the marker advertises.  CBvx.4.1/4.3 on-segment stay
#       fully GREEN (the fix never rejects valid in-segment queries).  CBvx.4.4
#       (VECTOR) stays BROKEN because src/VectorBVP.jl:318 was NOT touched —
#       confirming the two guards are independent loci, both needing the fix.
#       Restored: the guard is back to `real(t_star)` and both markers broken.
#
# RESTORATION: `git status --porcelain src/` confirmed EMPTY after the run.
#
# STANDALONE RUN:
#   julia --project=. test/corpus_oblique_guard_test.jl
# ============================================================================
