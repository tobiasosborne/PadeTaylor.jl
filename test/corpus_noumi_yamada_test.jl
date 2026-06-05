# test/corpus_noumi_yamada_test.jl
#
# ============================================================================
# Ground-truth corpus: NOUMI–YAMADA A_{2n}^{(1)} systems (B21/B14) — the
# Type-C self-validating residuals, the affine-Weyl Bäcklund orbit, and the
# FIRST EXTERNAL (non-rational) validation of the vector NY walk
# (bead padetaylor-lgpc, epic padetaylor-p1v0).
#
# Split from corpus_painleve_rational_test.jl for the Rule 6 ≤200-LOC cap.
# Corpus goal: ROOT OUT CORRECTNESS BUGS — verify every oracle (sympy /
# exact Rational{BigInt} / mpmath) before pinning.
#
#   CNY.6  Type-C exact residual via Rational{BigInt}: A_2 f_j=t/3, A_4
#          f_j=t/5; noumi_yamada_rhs residual is EXACTLY 0 (no float tol),
#          self-validating.  Constraints Σf_j=t and Σα=1 checked.
#   CNY.7  affine-Weyl Bäcklund ORBIT: s_0 on Type A reproduces the
#          catalogue pole-at-t=0 rational (A_2 and A_4); the transformed
#          tuple STILL solves the ODE (Rational{BigInt}); s_i²=1 and the
#          braid hold on these SPECIFIC solution tuples.  *** ERRATA: the
#          A_2 orbit α is (-1,1,1) not the catalogue's (-1,1,0). ***
#   CNY.8  ny-a4-generic-mpmath — the first EXTERNAL non-rational oracle:
#          mpmath odefun integration of A_4^(1) from a fixed IC, pinned via
#          capture.py; the vector NY walk (vector_solve_pade) must match.
#
# DE-DUPLICATION vs the existing suite (Rule 3)
# ---------------------------------------------
#   * noumi_yamada_test.jl NY.1.1/1.2 check the A_2/A_4 RHS vs a hand
#     oracle on FLOAT states; NY.1.4 checks Σf=t with FLOAT Type C.  CNY.6
#     pins the Type-C residual to EXACTLY 0 in Rational{BigInt} (the
#     self-validating closed form), which the float tests cannot.
#   * noumi_yamada_symmetry_test.jl NYS.1.1 checks rational-ODE residuals
#     and NYS.1.3-1.6 the group relations on GENERIC (non-solution) tuples.
#     CNY.7 instead pins the SPECIFIC catalogue Bäcklund-ORBIT tuples
#     (the s_0-image of Type A, pole at t=0) — the named B14/B21 oracle,
#     never tested by id before — and exposes the catalogue α error.
#   * NYS.1.2 vector-solves Type C (rational).  CNY.8 vector-solves a
#     GENERIC NON-rational transcendent against an mpmath odefun oracle —
#     the first external validation of the NY walk against a value with no
#     closed form.
#
# TOLERANCES: Type-C / Bäcklund residuals are EXACT (Rational{BigInt}, no
# tol).  The mpmath odefun oracle is dps=40; the Float64 vector walk is
# pinned at ~1e-12 (the walk is pole-free on [1, 1.5], so full-precision).
#
# MUTATION-PROOF (Rule 4): footer.  The NY-cyclic-RHS mutant (wrong cyclic
# index) → CNY.6/CNY.8 RED.
# ============================================================================

using Test
using PadeTaylor
using PadeTaylor.NoumiYamada: noumi_yamada_rhs, NoumiYamadaProblem
using PadeTaylor.NoumiYamadaSymmetry:
    noumi_yamada_rational, noumi_yamada_backlund
using PadeTaylor.VectorProblems: vector_solve_pade

@testset "Corpus: Noumi–Yamada (CNY)" begin

    # -----------------------------------------------------------------------
    # CNY.6  Type-C exact residual (Rational{BigInt}) — self-validating.
    # f_j = t/(2n+1) for all j, α_j = 1/(2n+1): each bracket vanishes (n
    # odd-offset and n even-offset equal terms), so f_j′ = 1/(2n+1) = α_j.
    # We assert the RESIDUAL noumi_yamada_rhs(t,f) - f′ is EXACTLY the zero
    # vector in Rational{BigInt} arithmetic — no floating-point tolerance.
    # -----------------------------------------------------------------------
    @testset "CNY.6  Type-C residual is exactly 0 (Rational{BigInt})" begin
        for n in (1, 2)
            d  = 2n + 1
            α  = [big(1)//d for _ in 1:d]          # 1/(2n+1) each
            @test sum(α) == 1                       # Σα = 1 (k=1)
            t  = big(13)//7                         # generic exact rational t
            f  = [t//d for _ in 1:d]                # Type C: f_j = t/d
            @test sum(f) == t                       # Σf_j = t (constraint)
            fprime = [big(1)//d for _ in 1:d]       # f_j′ = 1/d (constant)
            rhs = noumi_yamada_rhs(n; α = α)
            @test rhs(t, f) == fprime               # EXACT (Rational{BigInt})
            # And via the impl's own oracle factory (consistency of the two
            # entry points), still exact.
            αr, fr = noumi_yamada_rational(n, :C)
            @test rhs(t, fr(t)) == fprime
        end
    end

    # -----------------------------------------------------------------------
    # CNY.7  affine-Weyl Bäcklund ORBIT: s_0(Type A) = the pole-at-t=0
    # rational.  For A_2 (mod 3) index 0 has neighbours {1,2}, BOTH getting
    # +α_0, so α: (1,0,0) ↦ (-1,1,1) — Σα stays 1.  *** The catalogue's
    # α=(-1,1,0) (Σ=0) is WRONG (ERRATA): it fails the f_2 ODE.  We pin the
    # code-correct, sympy-verified (-1,1,1). ***  For A_4 (mod 5) the
    # neighbours {1,4} are distinct ⇒ α=(-1,1,0,0,1), Σ=1 (catalogue OK).
    # -----------------------------------------------------------------------
    @testset "CNY.7  Bäcklund orbit s_0(Type A): pole at t=0, ODE-exact" begin
        # --- A_2 ---
        f_A2 = [big(1)//1, big(0)//1, big(0)//1]      # (t,0,0) value at t=1
        α_A2 = [big(1)//1, big(0)//1, big(0)//1]
        f2, α2 = noumi_yamada_backlund(0, f_A2, α_A2)
        @test f2 == [1, 1, -1]                          # catalogue IC f(1)
        @test α2 == [-1, 1, 1]                          # ERRATA-corrected α
        @test sum(α2) == 1                              # k=1 preserved
        @test α2 != [-1, 1, 0]                          # anti-regression: NOT catalogue
        # The orbit rational f(t)=(t,1/t,-1/t) solves A_2 with α=(-1,1,1):
        rhs2 = noumi_yamada_rhs(1; α = α2)
        t = big(7)//3
        f_orbit2 = [t, 1//t, -1//t]
        fprime2  = [big(1)//1, -1//t^2, 1//t^2]
        @test rhs2(t, f_orbit2) == fprime2              # EXACT residual 0
        # s_0² = identity on this specific solution tuple.
        f2b, α2b = noumi_yamada_backlund(0, f2, α2)
        @test f2b == f_A2 && α2b == α_A2

        # --- A_4 (catalogue α=(-1,1,0,0,1) IS correct) ---
        f_A4 = [big(j == 1 ? 1 : 0)//1 for j in 1:5]    # (t,0,0,0,0) at t=1
        α_A4 = [big(j == 1 ? 1 : 0)//1 for j in 1:5]
        f4, α4 = noumi_yamada_backlund(0, f_A4, α_A4)
        @test f4 == [1, 1, 0, 0, -1]                    # catalogue IC f(1)
        @test α4 == [-1, 1, 0, 0, 1]                    # catalogue α (correct)
        @test sum(α4) == 1
        rhs4 = noumi_yamada_rhs(2; α = α4)
        s = big(5)//2
        f_orbit4 = [s, 1//s, big(0)//1, big(0)//1, -1//s]
        fprime4  = [big(1)//1, -1//s^2, big(0)//1, big(0)//1, 1//s^2]
        @test rhs4(s, f_orbit4) == fprime4              # EXACT residual 0
        # Braid relation on a generic (nonzero) tuple — exercises s_i on the
        # A_4 system the orbit lives in.  s_0 s_1 s_0 = s_1 s_0 s_1.
        fg = [big(j + 2)//3 for j in 1:5]
        αg = [big(j + 1)//11 for j in 1:5]
        b(a, c, e) = noumi_yamada_backlund(e,
                       noumi_yamada_backlund(c,
                         noumi_yamada_backlund(a, fg, αg)...)...)
        fL, αL = b(0, 1, 0)
        fR, αR = b(1, 0, 1)
        @test fL == fR && αL == αR                      # braid holds
    end

    # -----------------------------------------------------------------------
    # CNY.8  ny-a4-generic-mpmath — FIRST EXTERNAL non-rational validation.
    # mpmath odefun (dps=40, capture.py) integrates A_4^(1) from f(1)=
    # (0.3,0.25,0.15,0.2,0.1), α=(0.1,0.15,0.2,0.25,0.3) (Σ=1).  The vector
    # NY walk must match the pinned downstream values, and preserve Σf=t.
    # This is the first time the NY vector path-network is checked against a
    # value with NO closed form (only an external numerical oracle).
    # -----------------------------------------------------------------------
    @testset "CNY.8  A_4 generic transcendent vs mpmath odefun oracle" begin
        α  = [0.1, 0.15, 0.2, 0.25, 0.3]
        f0 = [0.3, 0.25, 0.15, 0.2, 0.1]
        @test sum(α) ≈ 1                                # Σα = 1
        @test sum(f0) ≈ 1.0                             # Σf0 = t0 = 1
        prob = NoumiYamadaProblem(2; α = α, f0 = f0, tspan = (1.0, 1.5),
                                  order = 30)
        sol = vector_solve_pade(prob; h = 0.05)
        # mpmath odefun oracle, dps=40 (capture.py):
        oracle12 = [0.331293362787775535, 0.267785942517965550,
                    0.195164178133955150, 0.245572975077188483,
                    0.160183541483115281]
        oracle15 = [0.374497248563133786, 0.296185710762128707,
                    0.265151309115081824, 0.312561288338753043,
                    0.251604443220902640]
        f12 = sol(1.2); f15 = sol(1.5)
        # The NY vector walk matches the EXTERNAL non-rational oracle.
        @test maximum(abs.(f12 .- oracle12)) < 1e-11
        @test maximum(abs.(f15 .- oracle15)) < 1e-11
        # Constraint invariant Σf_j = t along the whole (non-rational) walk.
        for (zk, yk) in zip(sol.z, sol.y)
            @test sum(yk) ≈ zk atol = 1e-12
        end
        @test sum(f15) ≈ 1.5 atol = 1e-12
        # Genuine transcendent, NOT Type C: components are unequal (a wrong
        # RHS that collapsed to the trivial f_j=t/5 would fail this).
        @test !all(fj -> isapprox(fj, 1.5/5; atol = 1e-3), f15)
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — EXECUTED, then restored byte-clean.
# `git diff --stat src/` confirmed empty after each restore.
#
# M-impl (NY-CYCLIC-targeted — the brief's required NY cyclic-RHS mutant):
#   in src/NoumiYamada.jl `noumi_yamada_rhs`, corrupt the cyclic index in
#   the bracket sum from `slot(j + 2r - 1) - slot(j + 2r)` to
#   `slot(j + 2r) - slot(j + 2r - 1)` (swap the odd/even offsets — a wrong
#   cyclic pairing).  EXECUTED — result: 4 failures — CNY.7's exact orbit
#   residuals (2, the unequal-component tuple no longer solves the ODE) and
#   CNY.8's oracle-match (2, the walk integrates a different system) both
#   went RED.  CNY.6 (Type C) correctly does NOT bite: with every f_j equal
#   to t/d the bracket is the sum of equal-minus-equal terms = 0 under
#   EITHER offset ordering, so a sign-swap is invisible to the symmetric
#   Type-C point — exactly why CNY.7/CNY.8 carry asymmetric tuples to
#   catch it.  Restored src/NoumiYamada.jl byte-for-byte (git diff empty).
#
# M-oracle (CNY.7, the ERRATA guard): pin the A_2 orbit α to the
#   catalogue's wrong `[-1, 1, 0]`.  Result: `α2 == [-1, 1, 1]` went RED
#   (the impl correctly returns (-1,1,1)) and `sum(α2) == 1` would fail for
#   the catalogue value — confirming the test enforces the corrected,
#   ODE-satisfying α and would catch any regression to the bad recipe.
#   Restored the literal.
#
# M-oracle (CNY.8): perturb oracle15[1] by +1e-9 (below the asserted
#   1e-11 band).  Result: CNY.8's `maximum(abs(f15 - oracle15)) < 1e-11`
#   went RED — confirming the assertion pins the mpmath value tightly, not
#   a loose band.  Restored the literal.
#
# Run standalone with:
#   julia --project=. test/corpus_noumi_yamada_test.jl
# ============================================================================
