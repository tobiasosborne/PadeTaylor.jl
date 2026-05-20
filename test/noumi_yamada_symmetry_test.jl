"""
V5c tests for `PadeTaylor.NoumiYamadaSymmetry` (bead `padetaylor-0ln.14`;
v0.2 plan row V5c).

`NoumiYamadaSymmetry` ships two related capabilities for the even-parity
Noumi–Yamada `A_{2n}^{(1)}` systems:

  (a) `noumi_yamada_rational(n, kind)` — exact rational-solution oracles
      (`kind ∈ {:A, :B, :C}`), each a closed-form solution that doubles
      as its own oracle;
  (b) the affine-Weyl-group `W(A_{2n}^{(1)})` Bäcklund symmetry:
      `noumi_yamada_backlund(i, f, α)` (reflection `s_i`) and
      `noumi_yamada_rotation(f, α)` (diagram rotation `π`), plus the
      algebraic group-relation self-test.

These tests (`NYS.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    NYS.1.1  rational solutions satisfy the ODE   — f_j′ = noumi_yamada_rhs
                                                    exactly (Rational algebra)
    NYS.1.2  Type C vector-solve cross-check      — A_4 f_j(t)=t/5 to ~1e-10
    NYS.1.3  s_i² = id                            — reflection is involutive
    NYS.1.4  commutativity + braid                — W(A_{2n}^{(1)}) relations
    NYS.1.5  rotation relations                   — π^{2n+1}=id, πs_i=s_{i+1}π
    NYS.1.6  Bäcklund preserves solutions         — s_i maps the solution
                                                    variety to itself
    NYS.1.7  failure modes + mutation-proof       — bad kind/i; recorded below

Representation: `f` and `α` are length-`2n+1` vectors of `Rational`
(or any number) — the Bäcklund action is exact rational arithmetic, so
`s_i²=id`, braid, and `π` relations are checked by exact `==`.  The
rational-solution oracle `noumi_yamada_rational` returns `(α, f)` where
`f` is a closure `t -> Vector` and `α` is the parameter vector.

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/noumi_yamada_symmetry_test.jl`) and under
`runtests.jl`.

Sources: `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`
lines 136–166 (`\\eqref{WAl1}`, `\\eqref{WAl2}`, group relations);
`docs/v0p2_pillarB_noumi_yamada_findings.md` §4 (generators) + §5
(rational solutions).
"""

using Test
using PadeTaylor
using PadeTaylor.NoumiYamadaSymmetry:
    noumi_yamada_rational, noumi_yamada_backlund, noumi_yamada_rotation
using PadeTaylor.NoumiYamada: noumi_yamada_rhs
using PadeTaylor.VectorProblems: vector_solve_pade

@testset "NoumiYamadaSymmetry (V5c)" begin

    # -------------------------------------------------------------------------
    # NYS.1.1 — rational solutions satisfy the ODE exactly.
    #
    # The closed forms f_j(t) are linear in t, so f_j′ is a constant.  We
    # substitute the closed form into `noumi_yamada_rhs` at a symbolic-exact
    # rational t and assert f_j′ == RHS_j with exact `==` (Rational algebra,
    # no floating-point tolerance).
    # -------------------------------------------------------------------------
    @testset "NYS.1.1 rational solutions satisfy the ODE" begin
        for n in (1, 2)
            d = 2n + 1
            # :B is the A_4-only intermediate solution (pillar B §5.2) —
            # exercise it only for n = 2; :A and :C are valid for all n.
            kinds = n == 2 ? (:A, :B, :C) : (:A, :C)
            for kind in kinds
                α, f = noumi_yamada_rational(n, kind)
                # Closed-form derivative: f_j(t) is c_j·t, so f_j′ = c_j;
                # extract c_j as f(1) - f(0) (exact for linear f).
                f0 = f(0 // 1)
                f1 = f(1 // 1)
                fprime = f1 .- f0
                # Evaluate the RHS at a generic rational t with the exact
                # closed-form state.  Σα = 1 must hold for noumi_yamada_rhs.
                @test sum(α) == 1
                rhs = noumi_yamada_rhs(n; α = α)
                t = 7 // 3
                ft = f(t)
                @test length(ft) == d
                rhs_val = rhs(t, ft)
                @test rhs_val == fprime          # exact equality (Rational)
            end
        end
    end

    # -------------------------------------------------------------------------
    # NYS.1.2 — Type C vector-solve cross-check (all-nonzero, no zero-component
    # degeneracy).  A_4 Type C: f_j(t) = t/5, α = (1/5,…,1/5).  The vector
    # Padé solver must reproduce f_j(t) = t/5 along the trajectory.
    # -------------------------------------------------------------------------
    @testset "NYS.1.2 Type C vector-solve cross-check" begin
        n = 2
        α, f = noumi_yamada_rational(n, :C)
        αf = float.(α)
        f0 = float.(f(1.0))                       # IC at t0 = 1.0
        prob = PadeTaylor.NoumiYamada.NoumiYamadaProblem(
            n; α = αf, f0 = f0, tspan = (1.0, 1.8))
        sol = vector_solve_pade(prob; h = 0.1)
        for (zk, yk) in zip(sol.z, sol.y)
            for fj in yk
                @test isapprox(fj, zk / 5; atol = 1e-10)
            end
        end
        # Dense interior point.
        for fj in sol(1.37)
            @test isapprox(fj, 1.37 / 5; atol = 1e-10)
        end
    end

    # -------------------------------------------------------------------------
    # NYS.1.3 — s_i² = id.  Applying the reflection s_i twice returns the
    # original (f, α) exactly, on a generic non-degenerate tuple.
    # -------------------------------------------------------------------------
    @testset "NYS.1.3 s_i² = id" begin
        for n in (1, 2)
            d = 2n + 1
            # Generic non-degenerate tuple: all f_j ≠ 0, α generic Rational.
            f = [Rational(j + 2, 3) for j in 0:(d - 1)]   # 2/3, 1, 4/3, …
            α = [Rational(j + 1, 7) for j in 0:(d - 1)]
            for i in 0:(d - 1)
                f1, α1 = noumi_yamada_backlund(i, f, α)
                f2, α2 = noumi_yamada_backlund(i, f1, α1)
                @test f2 == f
                @test α2 == α
                # s_i is not the identity itself (sanity: it actually moves).
                @test (f1 != f) || (α1 != α)
            end
        end
    end

    # -------------------------------------------------------------------------
    # NYS.1.4 — commutativity (s_i s_j = s_j s_i, non-adjacent i,j) and braid
    # (s_i s_{i±1} s_i = s_{i±1} s_i s_{i±1}, adjacent).
    # -------------------------------------------------------------------------
    @testset "NYS.1.4 commutativity + braid" begin
        # Use n = 2 (A_4, d = 5): indices 0..4, plenty of non-adjacent pairs.
        n = 2
        d = 2n + 1
        f = [Rational(j + 3, 5) for j in 0:(d - 1)]
        α = [Rational(j + 2, 11) for j in 0:(d - 1)]

        sij(i, j, f, α) = noumi_yamada_backlund(j, noumi_yamada_backlund(i, f, α)...)

        # Commutativity: i, j non-adjacent (and i ≠ j) mod d.
        for i in 0:(d - 1), j in 0:(d - 1)
            adjacent = (mod(i - j, d) == 1) || (mod(j - i, d) == 1) || i == j
            adjacent && continue
            fa, αa = sij(i, j, f, α)
            fb, αb = sij(j, i, f, α)
            @test fa == fb
            @test αa == αb
        end

        # Braid: adjacent i, j.  s_i s_j s_i = s_j s_i s_j.
        braid(a, b, c, f, α) = noumi_yamada_backlund(
            c, noumi_yamada_backlund(b, noumi_yamada_backlund(a, f, α)...)...)
        for i in 0:(d - 1)
            for j in (mod(i + 1, d), mod(i - 1, d))
                fL, αL = braid(i, j, i, f, α)
                fR, αR = braid(j, i, j, f, α)
                @test fL == fR
                @test αL == αR
            end
        end
    end

    # -------------------------------------------------------------------------
    # NYS.1.5 — rotation relations: π^{2n+1} = id and π s_i = s_{i+1} π.
    # -------------------------------------------------------------------------
    @testset "NYS.1.5 rotation relations" begin
        for n in (1, 2)
            d = 2n + 1
            f = [Rational(j + 2, 3) for j in 0:(d - 1)]
            α = [Rational(j + 1, 7) for j in 0:(d - 1)]

            # π^{2n+1} = id.
            fp, αp = f, α
            for _ in 1:d
                fp, αp = noumi_yamada_rotation(fp, αp)
            end
            @test fp == f
            @test αp == α
            # One rotation actually moves a non-constant tuple.
            f1, α1 = noumi_yamada_rotation(f, α)
            @test f1 != f

            # π s_i = s_{i+1} π  for every i (NY1998 main.tex:164).
            #
            # `noumi_yamada_backlund` / `noumi_yamada_rotation` each
            # implement a SINGLE generator's *value-action* on a tuple of
            # functions: φ maps the function f_j ↦ φ(f_j), and the helper
            # returns the value of φ(f_j) given the originals' values.  An
            # automorphism acts on an *expression* by substituting into its
            # variables, so the value-action of an operator COMPOSITE runs
            # in the reverse order of the operator composition (standard
            # contravariance):  value-action(π ∘ s_i) = vs_i ∘ vπ.
            # Hence the operator identity π s_i = s_{i+1} π is checked as
            #   vs_i ∘ vπ   ==   vπ ∘ vs_{i+1}.
            for i in 0:(d - 1)
                # LHS operator π ∘ s_i  →  value-action vs_i ∘ vπ.
                fL, αL = noumi_yamada_backlund(
                    i, noumi_yamada_rotation(f, α)...)
                # RHS operator s_{i+1} ∘ π  →  value-action vπ ∘ vs_{i+1}.
                fR, αR = noumi_yamada_rotation(
                    noumi_yamada_backlund(mod(i + 1, d), f, α)...)
                @test fL == fR
                @test αL == αR
            end
        end
    end

    # -------------------------------------------------------------------------
    # NYS.1.6 — Bäcklund preserves solutions.  Apply s_i to a rational
    # solution (f, α); the transformed (f′, α′) must STILL satisfy the
    # A_{2n}^{(1)} ODE (g_j′ = noumi_yamada_rhs(α′) at g).  This connects the
    # two halves: the affine-Weyl group genuinely acts on the solution
    # variety (NY1998 Theorem `\\ref{thmA}`, main.tex:174–178).
    #
    # We use Type C (all f_j = t/(2n+1) ≠ 0 for t ≠ 0): the reflection's
    # α_i/f_i term is well-defined.  The transformed solution g_j(t) is a
    # *rational* function of t (the term α_i/f_i(t) is not linear), so we
    # must use its EXACT analytic derivative, not a finite difference.
    # Since s_i commutes with d/dt:
    #   g_i      = f_i             ⇒ g_i′      = f_i′
    #   g_{i+1}  = f_{i+1}+α_i/f_i ⇒ g_{i+1}′  = f_{i+1}′ − α_i·f_i′/f_i²
    #   g_{i−1}  = f_{i−1}−α_i/f_i ⇒ g_{i−1}′  = f_{i−1}′ + α_i·f_i′/f_i²
    #   g_j      = f_j  (else)     ⇒ g_j′      = f_j′
    # All arithmetic is exact over Rational.
    # -------------------------------------------------------------------------
    @testset "NYS.1.6 Bäcklund preserves solutions" begin
        for n in (1, 2)
            d = 2n + 1
            α, f = noumi_yamada_rational(n, :C)
            # Closed-form derivative of the ORIGINAL solution (constant per j,
            # since each f_j is linear in t).
            fprime0 = f(1 // 1) .- f(0 // 1)
            slot(j) = mod(j, d) + 1
            for i in 0:(d - 1)
                t = 2 // 1
                ft = f(t)
                f_new, α_new = noumi_yamada_backlund(i, ft, collect(α))
                @test sum(α_new) == 1
                # Exact analytic derivative of the transformed solution.
                fi  = ft[slot(i)]
                fip = fprime0[slot(i)]
                drift = α[slot(i)] * fip / fi^2          # α_i·f_i′/f_i²
                g_prime = copy(fprime0)
                g_prime[slot(i + 1)] = fprime0[slot(i + 1)] - drift
                g_prime[slot(i - 1)] = fprime0[slot(i - 1)] + drift
                rhs_new = noumi_yamada_rhs(n; α = α_new)
                got = rhs_new(t, f_new)
                @test got == g_prime                     # exact (Rational)
            end
        end
    end

    # -------------------------------------------------------------------------
    # NYS.1.7 — failure modes (Rule 1) + mutation-proof record.
    # -------------------------------------------------------------------------
    @testset "NYS.1.7 failure modes" begin
        # Bad kind.
        @test_throws ArgumentError noumi_yamada_rational(2, :Z)
        @test_throws ArgumentError noumi_yamada_rational(2, :rational)
        # n < 1.
        @test_throws ArgumentError noumi_yamada_rational(0, :A)
        # Bad generator index for s_i (out of 0..2n).
        f = [Rational(j + 2, 3) for j in 0:4]
        α = [Rational(j + 1, 7) for j in 0:4]
        @test_throws ArgumentError noumi_yamada_backlund(5, f, α)
        @test_throws ArgumentError noumi_yamada_backlund(-1, f, α)
        # Mismatched f / α lengths.
        @test_throws ArgumentError noumi_yamada_backlund(0, f, α[1:4])
        @test_throws ArgumentError noumi_yamada_rotation(f, α[1:4])
        # Reflection at a zero component (α_i/f_i ⇒ division by zero) fails
        # loud rather than emitting ±Inf.
        f_zero = [Rational(0), Rational(1), Rational(1), Rational(1), Rational(1)]
        @test_throws ArgumentError noumi_yamada_backlund(0, f_zero,
            [Rational(1, 5) for _ in 0:4])
    end

    # -------------------------------------------------------------------------
    # NYS.1.7 — mutation-proof.  Each mutation was applied to
    # `src/NoumiYamadaSymmetry.jl`, the suite re-run, the source restored.
    #
    # M1 — wrong neighbour sign in the reflection: swap the +/− in
    #      `s_i(f_{i±1}) = f_{i±1} ± α_i/f_i`.  Bit: 8 failures, all in
    #      NYS.1.6 (the ODE-residual check).  s_i²=id, braid, commutativity
    #      and the rotation relations are invariant to a *consistent*
    #      neighbour-sign swap (the second application of s_i undoes it), so
    #      NYS.1.3–1.5 correctly do not bite — only the solution-variety
    #      check ties the reflection to the actual ODE and catches it.
    # M2 — wrong α reflection: `s_i(α_{i±1}) = α_{i±1} − α_i` instead of
    #      `+ α_i`.  Bit: 11 failures + 1 error — 10 in NYS.1.4 (the braid
    #      relation breaks), and NYS.1.6 (the transformed α no longer sums
    #      to 1, so `noumi_yamada_rhs`'s Σα = 1 guard throws).
    # M3 — wrong rotation direction: `π(f_j) = f_{j-1}` instead of `f_{j+1}`.
    #      Bit: 16 failures, all in NYS.1.5 — π^{2n+1}=id survives (π⁻¹ to
    #      the d-th power is still the identity) but the `π s_i = s_{i+1} π`
    #      relation and the "one rotation moves the tuple" sanity break.
    # -------------------------------------------------------------------------

end
