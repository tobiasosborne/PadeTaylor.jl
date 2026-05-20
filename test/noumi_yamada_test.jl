"""
V5a tests for `PadeTaylor.NoumiYamada` (bead `padetaylor-0ln.12`;
v0.2 plan row V5a; ADR-0022).

`NoumiYamada` is the problem-builder layer for the even-parity
Noumi–Yamada `A_{2n}^{(1)}` higher-order Painlevé systems
(NY1998 `\\eqref{A2n}`):

    f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1}  −  Σ_{1≤r≤n} f_{j+2r} ) + α_j ,
    indices mod (2n+1),

a `(2n+1)`-component first-order vector ODE.  `noumi_yamada_rhs(n; α)`
is the RHS closure factory; `NoumiYamadaProblem(n; α, f0, tspan)` is the
builder that assembles a `VectorPadeTaylorProblem` solvable by
`vector_solve_pade`.

These tests (`NY.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    NY.1.1  A_2^(1) RHS         — matches hand-written f_0'=f_0(f_1−f_2)+α_0
                                   and its cyclic rotations (pillar B §1.4)
    NY.1.2  A_4^(1) RHS         — matches f_0'=f_0(f_1−f_2+f_3−f_4)+α_0
                                   and cyclic rotations
    NY.1.3  cyclic structure    — component j is component 0 with indices
                                   rotated +j mod (2n+1)
    NY.1.4  constraint invariant — Σ_j f_j(t) = t along the whole solved
                                   trajectory (Σf_j' = Σα_j = 1)
    NY.1.5  builder             — NoumiYamadaProblem builds + solves
                                   end-to-end for A_2 and A_4, finite traj
    NY.1.6  failure modes       — α not summing to 1; wrong-length α/f0;
                                   n < 1 — all throw informatively
    NY.1.7  mutation-proof      — recorded at end of file

Oracle: the A_4^(1) Type C exact rational solution `f_j(t) = t/5` for
`α = (1/5,…,1/5)` (pillar B §5.4 / §7.2; Matsuda 2012 lines 322–338).

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/noumi_yamada_test.jl`) and under `runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.NoumiYamada: noumi_yamada_rhs, NoumiYamadaProblem
using PadeTaylor.VectorProblems: VectorPadeTaylorProblem, vector_solve_pade

# -----------------------------------------------------------------------------
# Test-local oracles — hand-written RHS, kept independent of the impl.
# -----------------------------------------------------------------------------

# A_2^(1) RHS written out by hand (NY1998 Appendix line 1198 + cyclic):
#   f_0' = f_0(f_1 − f_2) + α_0,  indices mod 3.
function _a2_rhs_oracle(f, α)
    m(j) = mod(j, 3) + 1                       # 0-based j → 1-based index
    return [f[m(j)] * (f[m(j+1)] - f[m(j+2)]) + α[m(j)] for j in 0:2]
end

# A_4^(1) RHS written out by hand (NY1998 Appendix line 1199 + cyclic):
#   f_0' = f_0(f_1 − f_2 + f_3 − f_4) + α_0,  indices mod 5.
function _a4_rhs_oracle(f, α)
    m(j) = mod(j, 5) + 1
    return [f[m(j)] * (f[m(j+1)] - f[m(j+2)] + f[m(j+3)] - f[m(j+4)]) + α[m(j)]
            for j in 0:4]
end

@testset "NoumiYamada (V5a)" begin

    # -------------------------------------------------------------------------
    # NY.1.1 — A_2^(1) RHS matches the hand-written formula + cyclic rotations.
    # -------------------------------------------------------------------------
    @testset "NY.1.1 A_2 RHS" begin
        α = [0.2, 0.3, 0.5]                    # Σα = 1
        f = noumi_yamada_rhs(1; α = α)
        fvec = [1.7, -0.4, 2.1]
        got = f(0.0, fvec)
        want = _a2_rhs_oracle(fvec, α)
        @test length(got) == 3
        @test got ≈ want
        # Spell out component 0 explicitly: f_0(f_1 − f_2) + α_0.
        @test got[1] ≈ fvec[1] * (fvec[2] - fvec[3]) + α[1]
        # And component 2 (cyclic): f_2(f_0 − f_1) + α_2.
        @test got[3] ≈ fvec[3] * (fvec[1] - fvec[2]) + α[3]
    end

    # -------------------------------------------------------------------------
    # NY.1.2 — A_4^(1) RHS matches the hand-written formula + cyclic rotations.
    # -------------------------------------------------------------------------
    @testset "NY.1.2 A_4 RHS" begin
        α = [0.1, 0.15, 0.2, 0.25, 0.3]        # Σα = 1
        f = noumi_yamada_rhs(2; α = α)
        fvec = [0.7, 1.3, -0.5, 2.2, 0.9]
        got = f(0.0, fvec)
        want = _a4_rhs_oracle(fvec, α)
        @test length(got) == 5
        @test got ≈ want
        # Component 0 explicitly: f_0(f_1 − f_2 + f_3 − f_4) + α_0.
        @test got[1] ≈ fvec[1] *
            (fvec[2] - fvec[3] + fvec[4] - fvec[5]) + α[1]
        # Component 3 (cyclic): f_3(f_4 − f_0 + f_1 − f_2) + α_3.
        @test got[4] ≈ fvec[4] *
            (fvec[5] - fvec[1] + fvec[2] - fvec[3]) + α[4]
    end

    # -------------------------------------------------------------------------
    # NY.1.3 — cyclic structure: component j is component 0 rotated +j.
    # -------------------------------------------------------------------------
    @testset "NY.1.3 cyclic structure" begin
        for n in (1, 2, 3)
            d = 2n + 1
            α = fill(1.0 / d, d)
            f = noumi_yamada_rhs(n; α = α)
            # Pick an asymmetric state so a wrong rotation would be caught.
            fvec = [Float64(j)^2 - 1.3j + 0.7 for j in 0:(d - 1)]
            got = f(0.0, fvec)
            for j in 0:(d - 1)
                # Component j on `fvec` must equal component 0 evaluated on
                # the index-rotated state (f shifted by +j, α shifted by +j).
                fvec_rot = [fvec[mod(j + i, d) + 1] for i in 0:(d - 1)]
                α_rot    = [α[mod(j + i, d) + 1]    for i in 0:(d - 1)]
                f_rot    = noumi_yamada_rhs(n; α = α_rot)
                got_rot  = f_rot(0.0, fvec_rot)
                @test got[j + 1] ≈ got_rot[1]
            end
        end
    end

    # -------------------------------------------------------------------------
    # NY.1.4 — constraint invariant: Σ_j f_j(t) = t along the trajectory.
    # The flow preserves Σf_j − t because Σf_j' = Σα_j = 1.
    # Oracle: A_4 Type C, f_j = t/5, α = (1/5,…,1/5) (pillar B §5.4).
    # -------------------------------------------------------------------------
    @testset "NY.1.4 constraint invariant Σf_j = t" begin
        # A_2: generic parameters, IC satisfying Σf0 = t0.
        α2 = [0.2, 0.3, 0.5]
        f0_2 = [0.5, 0.3, 0.2]                 # Σ = 1.0 = t0
        prob2 = NoumiYamadaProblem(1; α = α2, f0 = f0_2, tspan = (1.0, 1.6))
        sol2 = vector_solve_pade(prob2; h = 0.1)
        for (zk, yk) in zip(sol2.z, sol2.y)
            @test sum(yk) ≈ zk atol = 1e-9
        end
        # Dense interior point.
        @test sum(sol2(1.27)) ≈ 1.27 atol = 1e-9

        # A_4 Type C exact rational: f_j(t) = t/5.  Σf_j = t holds exactly,
        # and each f_j stays at t/5 (a genuine cross-check of the RHS).
        α4 = fill(0.2, 5)
        f0_4 = fill(0.2, 5)                    # t0 = 1.0 → f_j = 1/5
        prob4 = NoumiYamadaProblem(2; α = α4, f0 = f0_4, tspan = (1.0, 1.8))
        sol4 = vector_solve_pade(prob4; h = 0.1)
        for (zk, yk) in zip(sol4.z, sol4.y)
            @test sum(yk) ≈ zk atol = 1e-9
            for fj in yk
                @test fj ≈ zk / 5 atol = 1e-9   # Type C exact solution
            end
        end
    end

    # -------------------------------------------------------------------------
    # NY.1.5 — builder constructs + solves end-to-end, finite trajectory.
    # -------------------------------------------------------------------------
    @testset "NY.1.5 builder end-to-end" begin
        # A_2.
        prob2 = NoumiYamadaProblem(1; α = [0.2, 0.3, 0.5],
                                   f0 = [0.5, 0.3, 0.2], tspan = (1.0, 1.5))
        @test prob2 isa NoumiYamadaProblem
        @test prob2.n == 1
        @test prob2.α == [0.2, 0.3, 0.5]
        @test prob2.problem isa VectorPadeTaylorProblem
        sol2 = vector_solve_pade(prob2; h = 0.1)
        @test all(yk -> all(isfinite, yk), sol2.y)
        @test length(sol2.y[1]) == 3

        # A_4.
        prob4 = NoumiYamadaProblem(2; α = fill(0.2, 5),
                                   f0 = fill(0.2, 5), tspan = (1.0, 1.5))
        @test prob4.n == 2
        sol4 = vector_solve_pade(prob4; h = 0.1)
        @test all(yk -> all(isfinite, yk), sol4.y)
        @test length(sol4.y[1]) == 5
    end

    # -------------------------------------------------------------------------
    # NY.1.6 — failure modes all throw informatively (Rule 1).
    # -------------------------------------------------------------------------
    @testset "NY.1.6 failure modes" begin
        # α not summing to 1.
        @test_throws ArgumentError noumi_yamada_rhs(1; α = [0.1, 0.1, 0.1])
        @test_throws ArgumentError NoumiYamadaProblem(
            1; α = [0.1, 0.1, 0.1], f0 = [0.5, 0.3, 0.2], tspan = (1.0, 1.5))
        # α of wrong length (length must be 2n+1).
        @test_throws ArgumentError noumi_yamada_rhs(1; α = [0.5, 0.5])
        @test_throws ArgumentError NoumiYamadaProblem(
            1; α = [0.5, 0.5], f0 = [0.5, 0.3, 0.2], tspan = (1.0, 1.5))
        # f0 of wrong length.
        @test_throws ArgumentError NoumiYamadaProblem(
            1; α = [0.2, 0.3, 0.5], f0 = [0.5, 0.5], tspan = (1.0, 1.5))
        # n < 1.
        @test_throws ArgumentError noumi_yamada_rhs(0; α = [1.0])
        @test_throws ArgumentError NoumiYamadaProblem(
            0; α = [1.0], f0 = [1.0], tspan = (1.0, 1.5))
        # IC violating Σf0 = tspan[1].
        @test_throws ArgumentError NoumiYamadaProblem(
            1; α = [0.2, 0.3, 0.5], f0 = [0.5, 0.5, 0.5], tspan = (1.0, 1.5))
    end

    # -------------------------------------------------------------------------
    # NY.1.7 — mutation-proof.  Recorded below; the suite must be GREEN as
    # shipped, so the mutations are documented, not left active.  Each
    # mutation was applied to `src/NoumiYamada.jl`, the suite re-run, then
    # the source restored — all three bit.
    #
    # M1 — flip the even/odd sign in the bracket
    #      `(Σf_odd − Σf_even)` → `(Σf_odd + Σf_even)`:
    #      61 failures (NY.1.1: 3, NY.1.2: 3, NY.1.4: 55).  NY.1.3 is
    #      invariant to a uniform sign flip (it compares the impl against
    #      its own cyclic rotation) and correctly does not bite — the
    #      RHS-value and trajectory tests catch it instead.
    # M2 — drop the `+ α_j` term:
    #      61 failures (NY.1.1: 3, NY.1.2: 3, NY.1.4: 55).
    # M3 — wrong index-wrap modulus `mod(j, 2n)` instead of `mod(j, 2n+1)`:
    #      25 failures (NY.1.1: 3, NY.1.2: 3, NY.1.3: 12, NY.1.4: 7).
    # -------------------------------------------------------------------------

end
