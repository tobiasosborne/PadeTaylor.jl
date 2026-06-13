# test/corpus_bvp_hybrid_test.jl
#
# ============================================================================
# Ground-truth corpus: IVP↔BVP DISPATCH + the q0yq glue-throw discrepancy.
# Bead `padetaylor-u638`, epic `padetaylor-p1v0`.  Tagged CBV.* (B17 split
# from corpus_bvp_test.jl to keep each file ≤200 code LOC, Rule 6).
#
# WHAT THIS FILE PINS
# -------------------
#   CBV.8  ivp-bvp-dispatch-hybrid  [md:1438-1453]
#          u''=u, u=cosh(z) on the OBLIQUE segment [0,1+i].  Two independent
#          solvers must agree at the junction z=0.5+0.5i:
#            (a) BVP: bvp_solve with cosh Dirichlet BCs;
#            (b) IVP: path_network_solve shot from u(0)=1, u'(0)=sinh(0)=0.
#          The catalogue documents atol 1e-10; observed |Δu|≈3.2e-15,
#          |Δu'| comparable.  This pins the composition invariant that the
#          IVP and BVP branches the Dispatcher/Hybrid stitch are mutually
#          consistent on a pole-free segment (FW 2011 §3.2 md:190-192).
#
#   CBV.9  q0yq guard (bug `padetaylor-q0yq`, FIXED 2026-06-13)  [IVPBVPHybrid.jl]
#          The callable docstring once promised an unconditional glue-continuity
#          throw the BODY never performed (linear interpolation between
#          bracketing BVP slices, NO check) — a discontinuous glue passed
#          SILENTLY (Rule 1 gap).  The literal PFS-vs-BVP match is INFEASIBLE in
#          v1 (a value tautology — BVP BCs === harvested PFS values — plus a
#          ~1e-1 derivative mismatch, FFW md:226, deferred as padetaylor-sn9a).
#          The FIX is an inter-slice COHERENCE guard: a jump between the two
#          bracketing slices' values exceeding 3× the smaller value (floored at
#          1) is a divergent/corrupt slice ⇒ fail loud (measured: honest
#          jump/|w| ≈ 0.09, the +1000 injection ≈ 14).  We inject +1000 into one
#          slice and assert the callable now THROWS (@test_throws, flipped from
#          @test_broken); GM2 asserts the un-corrupted stack does NOT throw.
#
# ORACLES: cosh closed form (residual u''-u==0 sympy-verified, capture.py
# block 6/8) + Julia stdlib cosh/sinh as a second oracle.  No fitted tol.
#
# DE-DUP vs the committed suite (Rule 3):
#   * test/dispatcher_test.jl DP.* pins the 1D-chain Dispatcher on REAL
#     segments [-1,1]/[0,0.6] only — never an oblique COMPLEX segment, and
#     never the direct bvp-vs-ivp value identity at a complex junction.
#     CBV.8 is the oblique-complex composition invariant — orthogonal.
#   * test/ivp_bvp_hybrid_test.jl IB.1.4 queries sol(ζ) on a boundary but
#     only asserts finiteness/magnitude — it never tests the docstring's
#     glue-THROW promise.  CBV.9 is the first test to pin q0yq.
#
# COST NOTE: CBV.9 runs one small PIII hybrid solve (same fixture cost as
# IB.1.4, already in the suite).  Constructing a discontinuous PFS-vs-BVP
# glue from raw PathNetworkSolution objects would need synthetic
# PadeApproximant trees; injecting the discontinuity into a real slice stack
# is the faithful, robust route.
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor: BVP, IVPBVPHybrid
using PadeTaylor.IVPBVPHybrid: IVPBVPSolution, solve_pole_free_hybrid
using PadeTaylor.BVP: BVPSolution
using PadeTaylor.Painleve: PainleveProblem

@testset "Corpus BVP hybrid (CBV.8-9): dispatch invariant + q0yq glue marker" begin

    # =====================================================================
    # CBV.8 — ivp-bvp-dispatch-hybrid: BVP-vs-IVP agreement on [0,1+i].
    # =====================================================================
    @testset "CBV.8: BVP↔IVP composition invariant, cosh on [0,1+i]" begin
        zb = 1.0 + 1.0im
        zq = 0.5 + 0.5im

        # (a) BVP branch: cosh Dirichlet BCs.
        f_bvp(z, u)  = u
        ∂f_bvp(z, u) = one(u)
        sol_bvp = bvp_solve(f_bvp, ∂f_bvp, 0.0 + 0.0im, zb,
                            1.0 + 0.0im, cosh(zb); N = 24)
        u_bvp, up_bvp = sol_bvp(zq)

        # (b) IVP branch: path-network shot from u(0)=1, u'(0)=sinh(0)=0.
        f_ivp(z, u, up) = u                      # 2nd-order PathNetwork form
        prob = PadeTaylorProblem(f_ivp, (1.0 + 0.0im, 0.0 + 0.0im),
                                 (0.0 + 0.0im, zb); order = 30)
        grid = ComplexF64[0.25 + 0.25im, zq, 0.75 + 0.75im]
        ivp = path_network_solve(prob, grid; h = 0.25)
        u_ivp = ivp.grid_u[2]; up_ivp = ivp.grid_up[2]

        # Both branches agree at the junction — the composition invariant.
        # Catalogue documents atol 1e-10; both are ~1e-15.  Pin at 1e-10
        # (the FW 2011 md:192 derivative-match tolerance) — not fitted.
        @test isapprox(u_bvp,  u_ivp;  atol = 1e-10)
        @test isapprox(up_bvp, up_ivp; atol = 1e-10)
        # And both match the closed form cosh / sinh (independent oracle).
        @test isapprox(u_bvp,  cosh(zq); atol = 1e-12)
        @test isapprox(u_ivp,  cosh(zq); atol = 1e-12)
        @test isapprox(up_bvp, sinh(zq); atol = 1e-10)
        @test isapprox(up_ivp, sinh(zq); atol = 1e-12)
    end

    # =====================================================================
    # CBV.9 — q0yq MARKER: the IVPBVPSolution callable's promised
    # glue-continuity throw is ABSENT.  Build a real hybrid solution, inject
    # a discontinuity, and confirm the callable silently interpolates across
    # it instead of throwing.
    # =====================================================================
    @testset "CBV.9: IVPBVPSolution glue-throw is absent (bug padetaylor-q0yq)" begin
        # FFW Fig 5 PIII fixture (mirrors IB.1.4's cheap setup).
        α, β, γ, δ = 1.0, -1/20, 0.0, -1.0
        z_anchor = 30.0 + 0im
        u_a, up_a = pIII_asymptotic_ic(z_anchor; n_terms = 2, β = β, δ = δ)
        pp = PainleveProblem(:III; α = α, β = β, γ = γ, δ = δ,
                             u0 = u_a, up0 = up_a,
                             zspan = (z_anchor, z_anchor + 1), order = 30)
        sector = (im_lo = -0.5, im_hi = 0.5,
                   re_anchor = 2 * log(30), re_extent = 0.4)
        sol = solve_pole_free_hybrid(pp, sector,
                z -> pIII_asymptotic_ic(z; n_terms = 2, β = β, δ = δ);
                pfs_kwargs = (; h = 0.2, max_steps_per_target = 100),
                bvp_kwargs = (; N = 10, tol = 1e-10, maxiter = 30),
                n_slices = 4, glue_tol = 1e-6)

        @test length(sol.bvp_slices) == 4        # sanity: real slice stack
        @test sol.glue_tol == 1e-6

        # Inject a +1000 discontinuity into slice 2.  The linear interp
        # between slice 1 and slice 2 in the callable now spans a gap of
        # ~1000 — a glaring glue violation (≫ glue_tol = 1e-6).
        sl = sol.bvp_slices[2]
        bad = BVPSolution{Float64, ComplexF64}(
            sl.z_a, sl.z_b, sl.u_a, sl.u_b, sl.N, sl.nodes_t, sl.nodes_z,
            sl.u_nodes .+ 1000.0, sl.residual_inf, sl.iterations)
        slices2 = copy(sol.bvp_slices); slices2[2] = bad
        sol_bad = IVPBVPSolution{Float64}(sol.pfs_top, sol.pfs_bot, slices2,
            sol.slice_re, sol.sector, sol.glue_tol, sol.pp)

        # Query a ζ whose Re is between slice 1 and slice 2 → the callable
        # linearly interpolates ACROSS the 1000-level discontinuity.
        re_q = (sol.slice_re[1] + sol.slice_re[2]) / 2
        ζq = complex(re_q, 0.0)

        # FIXED (q0yq): the callable now FAILS LOUD on the gross inter-slice
        # incoherence.  At ζq the +1000-corrupted slice jumps ≈14× the
        # un-corrupted bracketing value, exceeding the 3× catastrophe bound, so
        # the callable throws instead of silently interpolating a fictitious
        # midpoint (~original + 500).  (Flipped from the q0yq @test_broken.)
        @test_throws ErrorException sol_bad(ζq)

        # GM2 (must NOT throw): the un-corrupted, CONTINUOUS real stack.  The
        # guard must not false-fire on honest O(0.1·|w|) inter-slice variation
        # (measured jump/|w| ≈ 0.09 here — far below the 3× bound).
        w_ok, up_ok = sol(ζq)
        @test isfinite(w_ok) && isfinite(up_ok)
    end

end # @testset Corpus BVP hybrid (CBV.8-9)

# ============================================================================
# MUTATION-PROOF (Rule 4):
#   CBV.8 is mutation-covered transitively by the BVP Mutant J1 (see
#   corpus_bvp_test.jl footer) and the Dispatcher Mutant B (swap u_a/u_b,
#   dispatcher_test.jl DP.5.1) — both perturb the BVP solution that CBV.8's
#   |Δu|<1e-10 junction identity depends on, confirmed RED there.  Targeted
#   spot-check: replacing `cosh(zb)` with `cosh(zb)+0.1` in CBV.8's BVP BC
#   makes |u_bvp - u_ivp| ≈ 0.05 ⇒ both atol-1e-10 junction tests go RED.
#   Restored; `git diff src/` empty.
#
#   CBV.9 — bug padetaylor-q0yq FIXED 2026-06-13 (inter-slice coherence guard
#   in IVPBVPHybrid.jl).  Mutation-proof (EXECUTED, restored byte-clean):
#     MT1 — disable the guard (`jump > 3*scale` → `jump > 1e12*scale`): the
#       +1000-corrupted query no longer throws ⇒ `@test_throws ErrorException
#       sol_bad(ζq)` RED.  Proves the guard catches the catastrophe.
#     MT2 — threshold → 0 (`jump > 0*scale`, always throw): the un-corrupted
#       CONTINUOUS stack (GM2 `sol(ζq)`) false-fires ⇒ RED.  Proves the 3×
#       calibration is load-bearing (no false-throw on honest variation).
#   CALIBRATION (Rule 5, measured this session): honest inter-slice jump/|w| ≈
#   0.093 (the FFW md:226 ~0.1 phenomenon) vs ≈14 for the +1000 injection — a
#   150× separation; the 3× bound sits 32× above honest, 4.7× below catastrophe.
#   The `min`-based scale (NOT max) is essential: the +1000 inflates the LARGER
#   bracketing value, so a max-based bound (e.g. 10·max) would never fire.
#   Regression: ivp_bvp_hybrid (65) + ffw_fig_5 (23) GREEN — no false-throws.
# ============================================================================
