# test/corpus_pathnet_walls_rows_test.jl
#
# ============================================================================
# Ground-truth corpus: PATH-NETWORK WALLS & ROWS — routing around structured
# pole sets (bead padetaylor-vt02, epic padetaylor-25og).
#
# WHAT THIS FILE PINS
# -------------------
# Four FW 2011 §3.1 path-network routing targets (Family I,
# docs/test_corpus/02_corpus_extension_plan.md:345-385).  Each background ODE
# has an EXACT closed form, so both the EVALUATED value (`eval_at`) and the
# CHOSEN ROUTE / extracted poles are pinnable independently.  The walk grows a
# tree from the IC, steers each 5-direction wedge step toward smallest |u| (away
# from poles), then dense-fills.
#
#   CPN.1  logistic VERTICAL WALL.  u'=u(1-u), recast u''=u(1-u)(1-2u),
#          u=1/(1+2 e^{-z}) (c=2) ⇒ pole wall at Re z=log 2, Im z=π(2k+1).
#          IC z0=0: u(0)=1/3, u'(0)=2/9 (exact).  Targets {2±2i, -1±2i} straddle
#          the wall.  ROUTING: NO visited node lands within ε=0.05 of any wall
#          pole — the tree threads the corridor.  VALUES: eval_at(target) pinned
#          to mpmath dps=50.
#
#   CPN.4  near-coalescent pair, δ-SWEEP.  Exact u=1/(z-a)+1/(z-(a+δ)), a=1, with
#          the explicit z-forcing f(z,u,up)=2/(z-a)^3+2/(z-(a+δ))^3 (= u'';
#          residual ≡ 0, sympy GATE).  IC z0=0.2.  δ∈{0.4,0.1,0.02,0.005}.
#          Couples routing to the FROISSART filter: extract_poles resolves TWO
#          distinct poles for δ≥0.1 and MERGES to one for δ≤0.02 (the default
#          cluster_atol=0.1 floor); NO spurious doublet appears between them.
#          Past-pole values are EXACT rationals (δ=0.1: u(0.5)=-11/3, u(1.5)=9/2).
#
#   CPN.8  off-node DENSE-EVAL (reuse logistic c=2).  Build the tree on a coarse
#          real-axis target set; eval_at STRICTLY BETWEEN visited nodes pinned to
#          the closed form.  NaN sentinel fires beyond every disc
#          (extrapolate=false); extrapolate=true returns a closed-form-consistent
#          value.
#
#   CPN.2  Airy-Riccati SEMI-INFINITE REAL ROW (the :min_u adversarial case).
#          f=2·u·up-1, u=-Ai'(z)/Ai(z), IC z0=0.5.  Poles at the Airy zeros a_k
#          on the NEGATIVE axis; between zeros u has its OWN zero (small |u|) next
#          to a pole spike — adversarial for |u|-smallest steering.  Targets march
#          past the row {-1,-1.5,-3,-5} plus oblique -2+0.5i.  ROUTING: no visited
#          node within ε of any Airy zero.  VALUES pinned to mpmath dps=50.
#
# GROUND TRUTH each case embodies
# -------------------------------
# external/probes/corpus-oracles/pathnet-routing/capture.py is the residual==0
# GATE (sympy: CPN.1 recast, CPN.4 forcing, CPN.2 Riccati all ≡ 0) AND the pin
# source (mpmath dps=50).  Routing assertions pin the walk's OUTPUT against the
# KNOWN pole set — "the tree stayed > ε from every known pole" — an honest,
# independent statement.  CPN.4 exact rationals and CPN.1/4 ICs are inline.
#
# TOLERANCES (justified by METHOD, never fitted to pass)
# ------------------------------------------------------
# The closed forms are entire-minus-poles meromorphic; the only error is the
# Float64 (15,15)-Padé / path-network floor.  Short-range walks here land far
# tighter than FW Table 5.1's ~1e-13 long-range floor; we pin the METHOD floor,
# NOT the lucky machine-epsilon residual.
#   * CPN.1 straddle values (no pole on the leg): rtol 1e-10 (empirically ~1e-16).
#   * CPN.1 wall-clearance ROUTING: no node within ε=0.05 (empirically ≥1.7).
#   * CPN.4 past-pole EXACT rationals: rtol 1e-10.
#   * CPN.4 pole locations (extract_poles): atol 1e-6 (figure-catalogue F64 spec,
#     polefield_test.jl PF.1.1); the 2-vs-1 merge flip is a discrete count.
#   * CPN.8 interior off-node values: 1e-8 (interior-of-disc Padé, brief);
#     extrapolate-recovered far value: 1e-6 (honest |t|>1 extrapolation floor).
#   * CPN.2 march values past the row: rtol 1e-9 (empirically ~3.6e-14);
#     Airy-zero ROUTING clearance ε=0.2 (empirically ≥0.34).
#
# REFERENCES: docs/test_corpus/02_corpus_extension_plan.md:345-385 (Family I);
# references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:155-166 (FW 2011 §3.1 path network).
# Oracle: capture.py (mpmath dps=50, residual==0 GATE).
#
# MUTATION-PROOF (Rule 4): see the comment block at the end of this file —
# one M-oracle (flip a pinned value → reddens) and one M-impl (perturb the
# wedge :min_u selector in src/PathNetwork.jl) were ACTUALLY executed RED and
# restored byte-for-byte (`git diff src/` empty).
# ============================================================================

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_corpus_pathnet_walls_rows.jl"))

# CPN.1 / CPN.8 logistic c=2 recast; CPN.2 Airy-Riccati recast.
flog(z, u, up) = u * (1 - u) * (1 - 2u)
fairy(z, u, up) = 2 * u * up - 1

@testset "Corpus: path-network walls/rows (CPN)" begin

    # ----------------------------------------------------------------------
    # CPN.1  logistic vertical wall (c=2).  Targets {2±2i, -1±2i} straddle the
    # wall at Re z=log 2, Im z=π(2k+1).  Pin values + wall-clearance routing.
    # ----------------------------------------------------------------------
    @testset "CPN.1  logistic vertical-wall threading" begin
        targets = ComplexF64[2 + 2im, 2 - 2im, -1 + 2im, -1 - 2im]
        prob = PadeTaylorProblem(flog, (1/3 + 0im, 2/9 + 0im),
                                 (0.0 + 0im, 3.0 + 0im); order = 30)
        sol = path_network_solve(prob, targets; h = 0.5, rng_seed = 0)

        # VALUES — pinned to mpmath dps=50.
        wants = [CPN1_U_2p2i, CPN1_U_2m2i, CPN1_U_m1p2i, CPN1_U_m1m2i]
        for (t, w) in zip(targets, wants)
            u, _ = eval_at(sol, t)
            @test isapprox(u, w; rtol = 1e-10)
        end

        # ROUTING — the tree threads the corridor: no visited node within
        # ε=0.05 of any wall pole log2 + iπ(2k+1), k in -3..3.
        walls = ComplexF64[CPN1_LOG2 + im * pi * (2k + 1) for k in -3:3]
        clearance = minimum(minimum(abs(zv - w) for w in walls)
                            for zv in sol.visited_z)
        @test clearance > 0.05
    end

    # ----------------------------------------------------------------------
    # CPN.4  near-coalescent pair δ-sweep.  Exact rationals + Froissart filter:
    # extract_poles resolves 2 poles for δ≥0.1, merges to 1 for δ≤0.02; no
    # spurious doublet.  Past-pole values are exact rationals.
    # ----------------------------------------------------------------------
    @testset "CPN.4  near-coalescent pair, δ-sweep (Froissart)" begin
        a = 1.0
        # (δ, expected distinct-pole count near the pair under the default
        #  cluster_atol=0.1).  2 above the floor, 1 (merged) below.
        sweep = [(0.4, 2), (0.1, 2), (0.02, 1), (0.005, 1)]
        for (δ, n_expected) in sweep
            b = a + δ
            frat(z, u, up) = 2 / (z - a)^3 + 2 / (z - b)^3
            ucf(z)  = 1 / (z - a) + 1 / (z - b)
            upcf(z) = -(1 / (z - a)^2 + 1 / (z - b)^2)
            z0 = 0.2
            u0  = ucf(z0 + 0im)
            up0 = upcf(z0 + 0im)
            targets = ComplexF64[-1.0, 0.5, 1.5, 2.5]
            prob = PadeTaylorProblem(frat, (u0, up0), (z0 + 0im, 3.0 + 0im);
                                     order = 30)
            sol = path_network_solve(prob, targets; h = 0.3, rng_seed = 0)

            # FROISSART / resolution: count distinct poles near the pair.
            poles = extract_poles(sol)          # DEFAULT cluster_atol=0.1
            near = filter(p -> abs(real(p) - 1.0) < 0.5 && abs(imag(p)) < 0.3,
                          poles)
            @test length(near) == n_expected     # 2 above floor, 1 merged below
            # No spurious off-axis doublet between the two poles.
            for p in near
                @test abs(imag(p)) ≤ 1e-3
            end
            # Resolved members land on the exact poles a and a+δ (δ≥0.1).
            if n_expected == 2
                @test minimum(abs(p - (a + 0im)) for p in near) ≤ 1e-6
                @test minimum(abs(p - (b + 0im)) for p in near) ≤ 1e-6
            end
        end

        # Past-pole EXACT rationals at δ=0.1: u(0.5)=-11/3, u(1.5)=9/2.
        let δ = 0.1, b = a + δ
            frat(z, u, up) = 2 / (z - a)^3 + 2 / (z - b)^3
            ucf(z) = 1 / (z - a) + 1 / (z - b)
            z0 = 0.2
            prob = PadeTaylorProblem(frat,
                                     (ucf(z0 + 0im), -(1/(z0-a)^2 + 1/(z0-b)^2)),
                                     (z0 + 0im, 3.0 + 0im); order = 30)
            sol = path_network_solve(prob, ComplexF64[0.5, 1.5];
                                     h = 0.3, rng_seed = 0)
            u05, _ = eval_at(sol, 0.5 + 0im)
            u15, _ = eval_at(sol, 1.5 + 0im)
            @test isapprox(u05, -11 / 3 + 0im; rtol = 1e-10)   # exact rational
            @test isapprox(u15, 9 / 2 + 0im;   rtol = 1e-10)   # exact rational
        end
    end

    # ----------------------------------------------------------------------
    # CPN.8  off-node dense-eval (reuse logistic c=2).  eval_at between nodes
    # pinned to the closed form; NaN sentinel beyond every disc; extrapolate.
    # ----------------------------------------------------------------------
    @testset "CPN.8  off-node dense-eval + NaN sentinel + extrapolate" begin
        prob = PadeTaylorProblem(flog, (1/3 + 0im, 2/9 + 0im),
                                 (0.0 + 0im, 2.0 + 0im); order = 30)
        coarse = ComplexF64[0.5, 1.0, 1.5, 2.0]
        h = 0.5
        sol = path_network_solve(prob, coarse; h = h, rng_seed = 0)

        # Off-node interior probes (real axis, pole-free): pinned to closed form.
        for (zp, w) in [(1.17, CPN8_U_1p17), (0.83, CPN8_U_0p83)]
            u, _ = eval_at(sol, zp + 0im)
            @test isfinite(u)
            @test isapprox(u, w + 0im; rtol = 1e-8)
        end

        # NaN sentinel beyond every visited disc (extrapolate=false, Rule 1).
        zfar = ComplexF64(5.0 + 5im)
        @test all(zv -> abs(zfar - zv) > h, sol.visited_z)   # truly uncovered
        uf, upf = eval_at(sol, zfar; extrapolate = false)
        @test isnan(real(uf)) && isnan(imag(uf))
        @test isnan(real(upf)) && isnan(imag(upf))

        # extrapolate=true recovers a closed-form-consistent value (|t|>1).
        ue, _ = eval_at(sol, zfar; extrapolate = true)
        ucf(z) = 1 / (1 + 2 * exp(-z))
        @test isfinite(ue)
        @test isapprox(ue, ucf(zfar); rtol = 1e-6)
    end

    # ----------------------------------------------------------------------
    # CPN.2  Airy-Riccati semi-infinite real row (the :min_u adversarial case).
    # March past the Airy-zero row {-1,-1.5,-3,-5} + oblique -2+0.5i; pin values
    # and assert no visited node within ε of any Airy zero.
    # ----------------------------------------------------------------------
    @testset "CPN.2  Airy-Riccati negative-axis row (:min_u stress)" begin
        u0  = CPN2_U0
        up0 = CPN2_UP0
        zeros = [CPN2_AIRY_ZERO_1, CPN2_AIRY_ZERO_2, CPN2_AIRY_ZERO_3,
                 CPN2_AIRY_ZERO_4, CPN2_AIRY_ZERO_5]
        targets = ComplexF64[-1.0, -1.5, -3.0, -5.0, -2 + 0.5im]
        wants = [CPN2_U_m1, CPN2_U_m1p5, CPN2_U_m3, CPN2_U_m5, CPN2_U_m2p0p5i]
        prob = PadeTaylorProblem(fairy, (u0 + 0im, up0 + 0im),
                                 (0.5 + 0im, -6.0 + 0im); order = 30)
        sol = path_network_solve(prob, targets; h = 0.3, rng_seed = 0,
                                 max_steps_per_target = 2000)

        # VALUES — march past a_1,a_2,a_3 on the negative axis + oblique target.
        for (t, w) in zip(targets, wants)
            u, _ = eval_at(sol, t)
            @test isapprox(u, w + 0im; rtol = 1e-9)
        end

        # ROUTING — the :min_u walk threads the row: no visited node within
        # ε=0.2 of any Airy zero, despite the small-|u| valleys next to spikes.
        clearance = minimum(minimum(abs(zv - z) for z in zeros)
                            for zv in sol.visited_z)
        @test clearance > 0.2
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
# M-oracle (CPN.1): in test/_oracle_corpus_pathnet_walls_rows.jl flip the sign
#   of the imaginary part of CPN1_U_2p2i (0.29024… → -0.29024…).  Result: ONLY
#   the CPN.1 `u(2+2i)` value assertion went RED; every routing/Froissart/other
#   assertion stayed GREEN.  Restored the literal; suite GREEN.
#
# M-impl (the load-bearing wedge selector): in src/PathNetwork.jl
#   `_select_candidate`, change the :min_u branch from
#       return argmin(abs(e[2]) for e in evals)
#   to
#       return argmax(abs(e[2]) for e in evals)
#   — steering the walk TOWARD the largest |u| (toward poles) instead of away.
#   Result (verified 2026-06-08): 23 passed, 2 FAILED, 1 ERRORED.  CPN.4 ERRORS
#   (a δ-sweep walk steered into a pole throws on a pole-incident wedge step,
#   and the Froissart pole count breaks); CPN.2 has 2 FAILS (a march value no
#   longer matches the closed form, and the Airy-zero clearance collapses below
#   ε=0.2).  CPN.8 and CPN.1 happened to survive for their particular grids
#   (CPN.1's targets sit far enough off-axis that even the perverse route reached
#   them clear of the wall) — the mutation still clearly bites, proving the
#   :min_u wedge selector is load-bearing for the routing/Froissart cases.
#   Restored src/PathNetwork.jl byte-for-byte (`git diff src/` empty); suite GREEN.
#
# Run standalone with:
#   julia --project=. test/corpus_pathnet_walls_rows_test.jl
# ============================================================================
