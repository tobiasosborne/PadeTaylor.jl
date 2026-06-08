# test/corpus_vector_bvp_test.jl
#
# ============================================================================
# Ground-truth corpus: GENERAL two-point boundary-condition assembly for the
# first-order vector BVP solver (bead padetaylor-g65t, epic padetaylor-p1v0).
#
# WHAT THIS FILE PINS — and why it is new (catalogue gap #10)
# -----------------------------------------------------------
# `VectorBVP.vector_bvp_solve` solves  y'(z) = f(z, y),  y ∈ ℂ^d, under the
# GENERAL linear two-point boundary condition
#
#     B_a·y(z_a) + B_b·y(z_b) = g          (d scalar equations).
#
# This general-B assembly is ADR-0023's entire rationale (the P_I⁽²⁾
# tritronquée companion BVP is NOT Dirichlet) — yet the committed suite
# tests ONLY the Dirichlet SELECTOR B (vector_bvp_test.jl VB.* and
# vector_bvp_wirein_test.jl all use B_a=[1 0;0 0], B_b=[0 0;1 0]: each row
# picks ONE component at ONE endpoint, no cross-endpoint coupling).  The
# τ-method endpoint↔node pairing — `B_a` hits the node-N block `y(z_a)`,
# `B_b` hits the node-0 block `y(z_b)`, because the DMSUITE node order is
# DESCENDING (t_0=+1↦z_b, t_N=-1↦z_a) — is exercised here for the first
# time with a B that genuinely couples both endpoints (VectorBVP.jl:229-239).
#
#   CVB.1  vector-antiperi-sincos-bvp  full-rank COUPLING across both ends:
#          B_a=I, B_b=-I (antiperiodic), g=(-1,1), [0, π/2].  The cleanest
#          non-Dirichlet oracle — every entry of both B blocks is active.
#   CVB.2  vector-endpoint-mixing-bvp  OFF-DIAGONAL / cross-endpoint mix:
#          B_a=[[0,0],[0,1]] (picks y2 at z_a), B_b=[[1,0],[0,0]] (picks y1
#          at z_b), g=(0,1), [0, π].  Exercises the off-diagonal assembly a
#          selector B never touches.
#   CVB.3  vector-rank-deficient-bc    Rule-1 fail-loud probe: a genuinely
#          RANK-DEFICIENT combined boundary operator → singular Jacobian →
#          the solver MUST throw (never silently return garbage).
#   CVB.4  vector-offdiag-coupling-bvp GENUINE off-diagonal COUPLING B (gap-#10
#          residual, 02_corpus_extension_plan.md Family H CVB.4 note): neither a
#          selector (CVB.2) nor a scalar ±I (CVB.1) but B_a=[[1,1],[0,1]],
#          B_b=[[1,0],[1,1]] — each block COUPLES y1 with y2 via a non-zero
#          OFF-DIAGONAL entry, so the B_a·Y_N + B_b·Y_0 τ-assembly is
#          load-bearing in its off-diagonal entries (a diagonal-only mutation of
#          the assembly changes the recovered function).
#
# THE ODE AND CLOSED FORM (one oracle, three BCs)
# -----------------------------------------------
# y'=(y2,-y1) with the closed form y=(sin z, cos z):  y1'=cos z=y2 ✓,
# y2'=-sin z=-y1 ✓ (sympy residual [0,0]).  NOTE this is the (sin,cos)
# orientation, NOT the (cos,-sin) of the committed VB.* harmonic oracle —
# a distinct closed form, so CVB.* are not value-duplicates of VB.*.
#
# ORACLE VERIFICATION (Law 1 — every g re-derived from the CLOSED FORM,
# never from the Julia routine; see capture.py in
# external/probes/corpus-oracles/vector-bvp/):
#   * CVB.1: B_a·y(0)+B_b·y(π/2) = I·(0,1) + (-I)·(1,0) = (-1,1) = g. sympy
#     confirms; the BC pins the UNIQUE solution A=1,B=0 of the general
#     family y=(A sin+B cos, A cos−B sin) — i.e. exactly (sin,cos).
#     [B_a|B_b] has full rank 2.  y(π/4)=(√2/2,√2/2) exact.
#   * CVB.2: row1 = y1(π) = sin π = 0; row2 = y2(0) = cos 0 = 1; g=(0,1).
#     sympy confirms UNIQUE A=1,B=0; [B_a|B_b] full rank 2.  y(π/2)=(1,0).
#   * CVB.3: B_a=[[1,0],[2,0]], B_b=0 ⇒ row2 = 2·row1 ⇒ [B_a|B_b] rank 1
#     (verified by sympy .rank()).  g=(1/2,1) is consistent (g2=2·g1) yet
#     NON-trivial (≠0 — dodging the catalogue-flagged g=(0,0) collapse to a
#     valid trivial solution).  Solving the BC leaves a FREE parameter
#     (A = 1−√3·B): the continuous problem is genuinely underdetermined, so
#     the discrete BC rows are linearly dependent and the Jacobian singular.
#   * CVB.4: B_a=[[1,1],[0,1]], B_b=[[1,0],[1,1]] on [0, π/2].
#     B_a·y(0)   = [[1,1],[0,1]]·(0,1) = (1,1);
#     B_b·y(π/2) = [[1,0],[1,1]]·(1,0) = (1,1);  g = (1,1)+(1,1) = (2,2).
#     sympy confirms [B_a|B_b] full rank 2 ⇒ UNIQUE A=1,B=0 = (sin,cos), and
#     that a DIAGONAL-ONLY mutant (zeroing the off-diagonal B_a[1,2]) instead
#     pins A=4/3,B=2/3 ≠ (sin,cos) — i.e. the off-diagonal coupling is
#     load-bearing.  y(π/4)=(√2/2,√2/2) exact.
#
# DE-DUPLICATION (CLAUDE.md Rule 3 — skepticism)
# ----------------------------------------------
#   * vector_bvp_test.jl VB.* / vector_bvp_wirein_test.jl VW.* pin ONLY the
#     Dirichlet selector B (each B row = one component at one endpoint).
#     CVB.1/CVB.2 add the COUPLING (antiperiodic full-rank) and OFF-DIAGONAL
#     (cross-endpoint mixing) B assembly — the untested general-B path.
#   * harmonic-oscillator-bvp-2vec (catalogue P2, B16) is per the catalogue
#     "already tested in vector_bvp_test.jl" (VB.1.2, selector B, closed
#     form (cos,-sin)).  REFERENCED, not re-pinned — a true duplicate.
#   * VB.4.1 already pins the fail-fast guards (N<4, dim mismatch, …) and a
#     non-convergence throw, but NOT a rank-deficient-BC singular-Jacobian
#     throw.  CVB.3 pins that distinct Rule-1 path (the ADR-0023/_solve
#     "rank-deficient" suggestion branch, VectorBVP.jl:463-474).
#
# TOLERANCES (justified by spectral accuracy, never fitted)
# ---------------------------------------------------------
# The linear-ODE collocation Newton converges in 1 step; the analytic
# solution is entire, so the Chebyshev interpolant is spectrally exact.
# Measured floors (probe, N as below, Float64): CVB.1 max node err 7.8e-16,
# y(π/4) err 5.6e-16; CVB.2 max node err 5.1e-16, y(π/2) err 2.3e-16.
# atol = 1e-12 carries ~3-4 orders of headroom above the ~1e-15 floor.
#
# MUTATION-PROOF (Rule 4): footer block — each load-bearing assertion was
# perturbed to confirm RED (incl. the BOUNDARY-ROW-ASSEMBLY mutant: swap the
# B_a/B_b endpoint↔node pairing) then restored byte-clean; `git diff src/`
# confirmed empty.
# ============================================================================

using Test
using LinearAlgebra: rank

# Standalone-runnable: include the module file directly (the same idiom as
# vector_bvp_test.jl), so `julia --project=. test/corpus_vector_bvp_test.jl`
# runs without a full package load.
if !isdefined(Main, :VectorBVP)
    include(joinpath(@__DIR__, "..", "src", "VectorBVP.jl"))
end
using .VectorBVP

# The shared ODE  y'=(y2,-y1)  and its closed form  y=(sin z, cos z).
f_sc(z, y) = [y[2], -y[1]]
y_exact(z) = [sin(z), cos(z)]

@testset "Corpus: general two-point BC assembly (CVB)" begin

    # ----------------------------------------------------------------------
    # CVB.1  vector-antiperi-sincos-bvp — full-rank COUPLING across both
    # endpoints.  B_a=I, B_b=-I (antiperiodic), g=(-1,1) on [0, π/2].
    # Oracle: g = I·y(0) + (-I)·y(π/2) = (-1,1); unique soln (sin,cos).
    # ----------------------------------------------------------------------
    @testset "CVB.1  antiperiodic B_a=I, B_b=-I: recover (sin z, cos z)" begin
        z_a, z_b = 0.0, π/2
        B_a = [1.0 0.0; 0.0 1.0]
        B_b = [-1.0 0.0; 0.0 -1.0]
        g   = [-1.0, 1.0]                       # = B_a·y(0)+B_b·y(π/2), sympy
        # The combined boundary operator is full rank ⇒ well-posed.
        @test rank([B_a B_b]) == 2
        sol = vector_bvp_solve(f_sc, z_a, z_b, B_a, B_b, g; N = 28)
        @test sol.iterations ≤ 2                # linear ⇒ 1 Newton step
        @test sol.residual_inf ≤ 1e-12
        # Interior node values to spectral accuracy — the gap-#10 assertion.
        for j in 1:sol.N+1
            @test isapprox(sol.Y_nodes[j], y_exact(sol.nodes_z[j]); atol = 1e-12)
        end
        # Off-node barycentric callable, incl. the catalogue spot value.
        @test isapprox(sol(π/4), [sqrt(2)/2, sqrt(2)/2]; atol = 1e-12)
        for z in (0.3, 0.9, 1.4)
            @test isapprox(sol(z), y_exact(z); atol = 1e-12)
        end
    end

    # ----------------------------------------------------------------------
    # CVB.2  vector-endpoint-mixing-bvp — OFF-DIAGONAL / cross-endpoint mix.
    # B_a=[[0,0],[0,1]] picks y2 at z_a; B_b=[[1,0],[0,0]] picks y1 at z_b.
    # g=(0,1) on [0, π].  Oracle: row1=y1(π)=sin π=0, row2=y2(0)=cos 0=1.
    # ----------------------------------------------------------------------
    @testset "CVB.2  mixing B (y2@z_a, y1@z_b): recover (sin z, cos z)" begin
        z_a, z_b = 0.0, Float64(π)
        B_a = [0.0 0.0; 0.0 1.0]               # picks y2 at z_a
        B_b = [1.0 0.0; 0.0 0.0]               # picks y1 at z_b
        g   = [0.0, 1.0]                        # = (y1(π), y2(0)) = (0, 1)
        @test rank([B_a B_b]) == 2
        sol = vector_bvp_solve(f_sc, z_a, z_b, B_a, B_b, g; N = 32)
        @test sol.iterations ≤ 2
        @test sol.residual_inf ≤ 1e-12
        for j in 1:sol.N+1
            @test isapprox(sol.Y_nodes[j], y_exact(sol.nodes_z[j]); atol = 1e-12)
        end
        # The interior midpoint the catalogue calls out, and the two pinned
        # boundary scalars themselves.
        @test isapprox(sol(π/2), [1.0, 0.0]; atol = 1e-12)
        @test isapprox(sol(0.0)[2], 1.0;        atol = 1e-12)   # y2(0)=1
        @test isapprox(sol(Float64(π))[1], 0.0; atol = 1e-12)   # y1(π)=0
        for z in (0.7, 1.9, 2.8)
            @test isapprox(sol(z), y_exact(z); atol = 1e-12)
        end
    end

    # ----------------------------------------------------------------------
    # CVB.3  vector-rank-deficient-bc — Rule-1 fail-loud probe.
    # B_a=[[1,0],[2,0]], B_b=0 ⇒ row2 = 2·row1 ⇒ [B_a|B_b] rank 1 (of 2).
    # g=(1/2,1) is consistent (g2=2·g1) but NON-trivial — so this is a
    # genuine underdetermined operator, NOT the catalogue's g=(0,0) collapse.
    # Per Rule 1 the singular Jacobian MUST throw (never a silent NaN lie).
    # ----------------------------------------------------------------------
    @testset "CVB.3  rank-deficient BC throws (Rule 1, no silent garbage)" begin
        z_a, z_b = π/6, π/2
        B_a = [1.0 0.0; 2.0 0.0]               # row2 = 2·row1
        B_b = [0.0 0.0; 0.0 0.0]
        g   = [0.5, 1.0]                        # consistent (g2=2g1), ≠ 0
        # Independently confirm the operator really is rank-deficient.
        @test rank([B_a B_b]) == 1             # underdetermined: rank 1 < d=2
        # The singular boundary block makes the Newton Jacobian singular; the
        # solver re-throws it as an ErrorException with a rank-deficient
        # suggestion (VectorBVP.jl _solve, lines 463-474).  Fail-loud, Rule 1.
        @test_throws ErrorException vector_bvp_solve(f_sc, z_a, z_b,
                                                     B_a, B_b, g; N = 16)
    end

    # ----------------------------------------------------------------------
    # CVB.4  vector-offdiag-coupling-bvp — GENUINE off-diagonal COUPLING B
    # (gap-#10 residual).  B_a=[[1,1],[0,1]] couples y1+y2 at z_a; B_b=
    # [[1,0],[1,1]] couples y1+y2 at z_b; g=(2,2) on [0, π/2].  Oracle (sympy,
    # from the closed form): B_a·y(0)=(1,1), B_b·y(π/2)=(1,1), g=(2,2); full
    # rank ⇒ UNIQUE (sin,cos).  Each off-diagonal B entry is load-bearing —
    # a diagonal-only assembly would pin A=4/3,B=2/3 ≠ (sin,cos).
    # ----------------------------------------------------------------------
    @testset "CVB.4  off-diagonal coupling B: recover (sin z, cos z)" begin
        z_a, z_b = 0.0, π/2
        B_a = [1.0 1.0; 0.0 1.0]               # row1 couples y1+y2 at z_a
        B_b = [1.0 0.0; 1.0 1.0]               # row2 couples y1+y2 at z_b
        g   = [2.0, 2.0]                        # = B_a·y(0)+B_b·y(π/2), sympy
        @test rank([B_a B_b]) == 2             # full-rank coupling ⇒ well-posed
        sol = vector_bvp_solve(f_sc, z_a, z_b, B_a, B_b, g; N = 24, tol = 1e-12)
        @test sol.iterations ≤ 2               # linear ⇒ 1 Newton step
        @test sol.residual_inf ≤ 1e-12
        # Interior node values to spectral accuracy — the gap-#10 assertion.
        for j in 1:sol.N+1
            @test isapprox(sol.Y_nodes[j], y_exact(sol.nodes_z[j]); atol = 1e-12)
        end
        # Off-node barycentric callable, incl. the catalogue spot value.
        @test isapprox(sol(π/4), [sqrt(2)/2, sqrt(2)/2]; atol = 1e-12)
        for z in (0.2, 0.8, 1.3)
            @test isapprox(sol(z), y_exact(z); atol = 1e-12)
        end
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
# `git diff --stat src/` confirmed empty after each restore.
#
# M-oracle (CVB.1): negate the expected closed form in CVB.1's node loop —
#   `isapprox(sol.Y_nodes[j], y_exact(...))` → `… -y_exact(...)`.  Result:
#   CVB.1 RED at the first node (the recovered (sin,cos) ≠ its negation);
#   CVB.2/CVB.3 stayed GREEN.  Restored byte-for-byte.
#
# M-impl (the BOUNDARY-ROW-ASSEMBLY mutant — the gap-#10 target): in
#   src/VectorBVP.jl, swap the τ-method endpoint↔node pairing — write
#     `J_bc[:, nodeN+1:nodeN+d] .= Bb`   (was `.= Ba`)
#     `J_bc[:, 1:d]               .= Ba`   (was `.= Bb`)
#   and the matching swap in `_residual`'s BC line
#     `R[ctx.bc_rows] .= ctx.Bb * Y[n+1:n+d] .+ ctx.Ba * Y[1:d] .- ctx.g_CT`.
#   This mis-reads the descending DMSUITE node order (B_a must hit y(z_a) =
#   node N, B_b must hit y(z_b) = node 0).  Result: CVB.1 AND CVB.2 went RED
#   — Newton still converges (residual ~1e-14) but to a DIFFERENT function
#   that satisfies the SWAPPED BC, not the user's; the node values no longer
#   match (sin,cos).  (For CVB.1's symmetric B_a=I/B_b=-I the swap flips the
#   sign of the coupling, breaking uniqueness; for CVB.2's asymmetric mixing
#   B the swap pins the WRONG component at the WRONG endpoint.)  CVB.1/CVB.2
#   are therefore the load-bearing assembly assertions.  (CVB.3 incidentally
#   also went RED under this mutant: the swapped rank-1 block lands on the
#   node-0 columns and at N=16 the resulting Jacobian was no longer singular,
#   so the expected throw did NOT fire — an extra catch, not the targeted
#   bite.)  This is exactly the τ-method trap ADR-0023 / VB.* Mutation B warn
#   about, here localised to the general/coupling/mixing-B path.  Restored
#   byte-for-byte.
#
# M-oracle (CVB.4): flip the CVB.4 boundary right-hand side g — `g=[2.0,2.0]`
#   → `g=[3.0,2.0]`.  Result: CVB.4 RED (29 sub-asserts — the node loop, the
#   off-node loop, and the y(π/4) spot value all detect the recovered function
#   ≠ (sin,cos) once g is wrong); CVB.1/CVB.2/CVB.3 stayed GREEN (36/42/2).
#   Restored byte-for-byte.
#
# M-impl (the OFF-DIAGONAL τ-ASSEMBLY mutant — CVB.4's gap-#10 target): in
#   src/VectorBVP.jl, drop the off-diagonal contribution of B_a in the
#   `B_a·Y_N` assembly — `J_bc[:, nodeN+1:nodeN+d] .= Diagonal(Ba)` (was `.=
#   Ba`) and the matching `_residual` line `R[ctx.bc_rows] .= Diagonal(ctx.Ba)
#   * Y[n+1:n+d] .+ ctx.Bb*Y[1:d] .- ctx.g_CT` (with `Diagonal` added to the
#   LinearAlgebra import).  This keeps the τ-pairing correct but silently zeros
#   the cross-component coupling B_a[1,2].  Result: CVB.4 RED (29 sub-asserts)
#   — Newton converges but to the diagonal-only function A=4/3,B=2/3 ≠
#   (sin,cos) (sympy-predicted in capture.py CASE 4); CVB.1 (B_a=I, already
#   diagonal) and CVB.2 (B_a=[[0,0],[0,1]], already diagonal) stayed GREEN
#   (36/42 — the mutation is a no-op on a diagonal B_a), and CVB.3 still threw
#   (2 GREEN).  So CVB.4 is the UNIQUE detector of an off-diagonal assembly bug
#   — its off-diagonal coupling is provably load-bearing.  Restored
#   byte-for-byte.
#
# Run standalone with:
#   julia --project=. test/corpus_vector_bvp_test.jl
# ============================================================================
