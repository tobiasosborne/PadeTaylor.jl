"""
V5b test — `A_2^(1)` ⇒ scalar PIV backward-compatibility self-validation
(bead `padetaylor-0ln.13`; v0.2 PRD acceptance item "v0.1 PIV
self-validation").

The Noumi–Yamada system `A_2^(1)` (the `n = 1` even-parity system, three
components `f_0, f_1, f_2`) is, after the documented change of
variables, the scalar fourth Painlevé equation PIV.  This file proves
that v0.2's *vector* solve of `A_2^(1)` (`NoumiYamadaProblem(1; …)` →
`vector_solve_pade`) reproduces v0.1's *scalar* PIV solve
(`PainleveProblem(:IV)` → `solve_pade`) under that change of variables —
the v0.2 ↔ v0.1 backward-compatibility guarantee.

## The source-verified A_2^(1) ↔ PIV map (k = 1 convention)

NY1999 derives the reduction in its own `k = 3` normalisation
(`Σf_j = 3x`, `Σα_j = 3`).  `NoumiYamada.jl` uses the project `k = 1`
normalisation (`Σf_j = t`, `Σα_j = 1`).  The k=1 map is **not** the
NY1999 §2 formulas read literally — the findings doc §2.2 explicitly
warns to rescale, and the rescaled map is *different in shape* (the
NY1999 `c = √(-3/2)` and `y = -f_1/c` are k=3 statements; in the k=1
convention the imaginary `c` cancels and the map is real).

Derivation.  Differentiate the project k=1 `A_2^(1)` system

    f_0' = f_0(f_1−f_2)+α_0,  f_1' = f_1(f_2−f_0)+α_1,
    f_2' = f_2(f_0−f_1)+α_2,   Σf_j = τ,  Σα_j = 1

(NY1998 main.tex:1198 + cyclic; `NoumiYamada.jl`).  Eliminating
`f_0, f_2` via `f_2−f_0 = (f_1'−α_1)/f_1`, `f_2+f_0 = τ−f_1`,
`4f_0f_2 = (f_2+f_0)² − (f_2−f_0)²` gives the scalar ODE for `f_1(τ)`:

    f_1'' = f_1'²/(2f_1) + (3/2)f_1³ − 2τ f_1²
            + ((τ²/2)+α_2−α_0) f_1 − α_1²/(2f_1).

Substituting `f_1 = y/√2`, `τ = -√2 t` (so `y' = dy/dt = -2 f_1'`,
`y'' = 2√2 f_1''`) turns this term-for-term into the canonical PIV of
`src/Painleve.jl` / RF 2014,

    y'' = (y')²/(2y) + (3/2)y³ + 4t y² + 2(t²−a)y + b/y,

with the **k = 1 parameter map**

    y = √2 · f_1           f_1 = y/√2
    τ = -√2 · t            t   = -τ/√2
    a = α_0 − α_2          b   = -2 α_1²
    α_0 + α_1 + α_2 = 1

(equivalently `a = 1 − α_1 − 2α_2`, since `α_0 = 1−α_1−α_2`).  Every
term of the substitution was checked by hand and cross-checked against
the NY1999 k=3 chain (`α_j = A_j/3`, `F_1 = -i√3 f_1`, `c = i√(3/2)` ⇒
`y = -F_1/c = √2 f_1`).  Both routes agree.

Sources, verified from the LaTeX (Law 1):
  - `references/tex/noumi_yamada/NoumiYamada1999_PIV_symmetries_okamoto_NagoyaMathJ153/main.tex`
    lines 420–576 — PIV eq. `\eqref{p4y}`, the symmetric form
    `\eqref{p4f}`, the reduction `f_1 = -cy`, `x = -t/c`,
    `c = √(-3/2)`, and the k=3 parameter relations
    `a = 1+3v_3`, `b = -2(v_1−v_2)²`, `α_j` ↔ `v` (lines 450–457,
    501–519, 568–575).
  - `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`
    line 1198 — `A_2^(1): f_0' = f_0(f_1−f_2)+α_0` (the k=1 system
    `NoumiYamada.jl` integrates).
  - `docs/v0p2_pillarB_noumi_yamada_findings.md` §2 — the survey of the
    reduction; its §2.2 map is correct in the k=3 convention but is
    NOT the k=1 map (see the report; no doc error, an expected
    rescaling).

## The checks (NYPIV.*)

  NYPIV.1  trajectory agreement — v0.2 `√2·f_1(τ)` ≡ v0.1 PIV `y(t)`
           and `-2·f_1'(τ)` ≡ v0.1 PIV `y'(t)` along the trajectory,
           under `t = -τ/√2`, tight tolerance.  The IC translation is
           documented inline.
  NYPIV.2  ODE residual — `y = √2 f_1` from the v0.2 solve satisfies
           the scalar PIV ODE residual to a tight tolerance.  Uses
           `y', y''` derived *analytically* from the A_2^(1) system
           (not finite differences), so this check depends on neither
           v0.1 nor numerical differentiation — it isolates whether the
           reduction itself holds.
  NYPIV.3  closed-form anchor — for `(a,b) = (0,-2)` the PIV closed
           form `u(z) = -2z` (`piv_entire(:minus_2z)`) is an
           independent third oracle: v0.1 PIV reproduces it, and the
           A_2^(1) parameters mapping to `(0,-2)` are `α = (0,1,0)`
           (`a = α_0−α_2 = 0`, `b = -2α_1² = -2`), whose A_2^(1) seed
           solution is `f_1 = τ` ⇒ `y = √2·τ = √2·(-√2 t) = -2t`.
           (The exact seed `f = (0,τ,0)` has zero components and is a
           shared-`Q` Padé degeneracy for the vector solver — a
           documented v1 limit — so the closed-form leg is verified on
           the scalar side; combined with NYPIV.1 it pins the chain
           v0.2 ≡ v0.1 ≡ closed form.)
  NYPIV.4  mutation-proof — recorded in the file footer.

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/noumi_yamada_piv_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.NoumiYamada: NoumiYamadaProblem, vector_solve_pade
using PadeTaylor.Painleve:    PainleveProblem, piv_entire
import PadeTaylor: solve_pade

# -----------------------------------------------------------------------------
# The source-verified A_2^(1) ↔ PIV parameter map (k = 1 convention).
# A small test-local helper — no new src/ module (bead scope).  Kept
# beside the test so a wrong map is a one-line, mutation-provable change.
# -----------------------------------------------------------------------------

# PIV parameters (a, b) for an A_2^(1) parameter triple α = (α_0,α_1,α_2):
#   a = α_0 − α_2,   b = -2 α_1²    (k = 1; see file docstring).
_a2_to_piv_params(α) = (a = α[1] - α[3], b = -2 * α[2]^2)

# Independent-variable map:  the NY system runs in τ, PIV runs in t,
# with  t = -τ/√2  (equivalently τ = -√2 t).
_tau_to_t(τ) = -τ / sqrt(2)

# Dependent-variable map:  y(t) = √2 · f_1(τ).
_f1_to_y(f1) = sqrt(2) * f1

# Derivative map:  y'(t) = dy/dt = -2 · f_1'(τ).  (Chain rule with
# dτ/dt = -√2 and y = √2 f_1: y' = √2·f_1'·(dτ/dt)·… ⇒ y' = -2 f_1'.)
_f1prime_to_yprime(f1p) = -2 * f1p

# The A_2^(1) RHS, written out by hand (NY1998 main.tex:1198 + cyclic),
# kept independent of the impl so the derivatives below are an oracle.
function _a2_rhs(f, α)
    f0, f1, f2 = f
    return (f0 * (f1 - f2) + α[1],     # f_0'
            f1 * (f2 - f0) + α[2],     # f_1'
            f2 * (f0 - f1) + α[3])     # f_2'
end

@testset "NoumiYamada A_2^(1) ⇒ PIV self-validation (V5b)" begin

    # -------------------------------------------------------------------------
    # A concrete, fully generic test case.  α chosen so all three of
    # (α_0, α_1, α_2) are distinct and nonzero ⇒ no component of the
    # solution is forced to vanish (the vector solver's shared-Q Padé
    # needs every component nonzero), and (a, b) is a generic PIV point.
    # -------------------------------------------------------------------------
    α = [0.5, 0.3, 0.2]                       # Σα = 1 (k = 1)
    pv = _a2_to_piv_params(α)
    a, b = pv.a, pv.b                         # a = 0.3, b = -0.18

    @testset "parameter map sanity" begin
        @test sum(α) ≈ 1                       # k = 1 normalisation
        @test a ≈ α[1] - α[3]
        @test a ≈ 1 - α[2] - 2 * α[3]          # the equivalent form
        @test b ≈ -2 * α[2]^2
    end

    # v0.2 vector solve of A_2^(1).  τ-window [1.0, 1.4]; IC chosen with
    # Σf0 = τ0 = 1.0 (the NoumiYamadaProblem constraint) and all
    # components distinct so the trajectory is genuinely transcendent.
    τspan = (1.0, 1.4)
    f0    = [0.6, 0.5, -0.1]                  # Σ = 1.0 = τspan[1]
    @assert sum(f0) ≈ τspan[1]
    nyprob = NoumiYamadaProblem(1; α = α, f0 = f0, tspan = τspan)
    nysol  = vector_solve_pade(nyprob; h = 0.05)

    # τ sample grid for trajectory comparisons (interior, away from the
    # exact endpoints so dense evaluation is unambiguous).
    τgrid = 1.0:0.1:1.4

    # -------------------------------------------------------------------------
    # NYPIV.1 — trajectory agreement: v0.2 A_2^(1) ≡ v0.1 scalar PIV.
    #
    # IC translation across the change of variables.  PIV runs in
    # t = -τ/√2, so as τ increases over [1.0, 1.4], t *decreases* over
    # [-1/√2, -1.4/√2].  v0.1 `solve_pade` steps z_start → z_end with
    # h_max > 0, i.e. it needs an increasing span.  We therefore seed
    # the v0.1 PIV solve from the τ = 1.4 end (the smaller t) and
    # integrate up to the τ = 1.0 end (the larger t):
    #
    #     PIV IC at t1 = -1.4/√2 :  y1  = √2 · f_1(τ=1.4)
    #                               y'1 = -2  · f_1'(τ=1.4)
    #     PIV span: (t1, t0) with t0 = -1.0/√2  (t0 > t1).
    #
    # f_1'(τ) is taken from the A_2^(1) RHS evaluated at the v0.2 state
    # — exact, not differenced.
    # -------------------------------------------------------------------------
    @testset "NYPIV.1 trajectory agreement v0.2 ≡ v0.1" begin
        τ1   = τspan[2]
        f_τ1 = nysol(τ1)
        _, f1p_τ1, _ = _a2_rhs(f_τ1, α)
        y1   = _f1_to_y(f_τ1[2])
        yp1  = _f1prime_to_yprime(f1p_τ1)
        t1   = _tau_to_t(τ1)
        t0   = _tau_to_t(τspan[1])
        @test t0 > t1                          # PIV span must be increasing

        pivprob = PainleveProblem(:IV; α = a, β = b,
                                  u0 = y1, up0 = yp1, zspan = (t1, t0))
        pivsol  = solve_pade(pivprob; h = 0.05)

        # Endpoints first: at τ = 1.4 the PIV IC IS the translated v0.2
        # state, so agreement there is exact-by-construction; the real
        # test is the *interior* and the far endpoint τ = 1.0.
        for τ in τgrid
            f_τ = nysol(τ)
            _, f1p_τ, _ = _a2_rhs(f_τ, α)
            y_v02  = _f1_to_y(f_τ[2])
            yp_v02 = _f1prime_to_yprime(f1p_τ)
            t = _tau_to_t(τ)
            y_v01, yp_v01 = pivsol(t)
            @test y_v02  ≈ y_v01  atol = 1e-9   # trajectory agreement
            @test yp_v02 ≈ yp_v01 atol = 1e-9   # derivative agreement
        end

        # Spell out the far-endpoint agreement explicitly (τ = 1.0, the
        # end furthest from the shared IC — a genuine cross-check, not a
        # tautology): √2·f_1(1.0) equals the v0.1 PIV value there.
        f_far = nysol(1.0)
        @test _f1_to_y(f_far[2]) ≈ pivsol(_tau_to_t(1.0))[1] atol = 1e-9
    end

    # -------------------------------------------------------------------------
    # NYPIV.2 — ODE residual: y = √2 f_1 satisfies the scalar PIV ODE.
    #
    # Independent of v0.1 entirely.  y', y'' are derived *analytically*
    # from the A_2^(1) system (no finite differences):
    #     y'  = -2 f_1',           f_1' from the RHS;
    #     y'' = 2√2 f_1'',         f_1'' = f_1'(f_2−f_0) + f_1(f_2'−f_0'),
    #                              every term RHS-evaluable.
    # The residual of  y'' − [(y')²/(2y) + (3/2)y³ + 4ty² + 2(t²−a)y
    #                          + b/y]  must vanish — this isolates the
    # reduction itself, not the comparison to v0.1.
    # -------------------------------------------------------------------------
    @testset "NYPIV.2 scalar-PIV ODE residual" begin
        maxres = 0.0
        for τ in 1.05:0.05:1.35               # interior sample points
            f = nysol(τ)
            f0v, f1, f2 = f
            f0p, f1p, f2p = _a2_rhs(f, α)
            # f_1'' from differentiating f_1' = f_1(f_2−f_0)+α_1.
            f1pp = f1p * (f2 - f0v) + f1 * (f2p - f0p)
            y   = _f1_to_y(f1)
            yp  = _f1prime_to_yprime(f1p)
            ypp = 2 * sqrt(2) * f1pp           # y'' = 2√2 f_1''
            t   = _tau_to_t(τ)
            res = ypp - (yp^2 / (2y) + (3 // 2) * y^3 + 4 * t * y^2 +
                         2 * (t^2 - a) * y + b / y)
            maxres = max(maxres, abs(res))
            @test abs(res) < 1e-9              # PIV ODE satisfied
        end
        # The reduction holds to integrator accuracy (≈ 4e-12 observed).
        @test maxres < 1e-9
    end

    # -------------------------------------------------------------------------
    # NYPIV.3 — closed-form anchor: PIV `u(z) = -2z` at (a,b) = (0,-2).
    #
    # Third independent oracle.  The A_2^(1) parameters mapping to the
    # closed-form PIV point (a,b) = (0,-2) are α = (0,1,0):
    #     a = α_0 − α_2 = 0,   b = -2 α_1² = -2.
    # The corresponding A_2^(1) seed solution is f_1 = τ (and f_0 =
    # f_2 = 0), so y = √2 f_1 = √2 τ = √2·(-√2 t) = -2t — exactly the
    # PIV closed form `u = -2z`.
    #
    # The exact seed f = (0,τ,0) has zero components, a shared-Q Padé
    # degeneracy for the vector solver (documented v1 limit), so the
    # closed-form ≡ vector leg is not run directly; instead this anchor
    # confirms v0.1 PIV reproduces the closed form, which — together
    # with NYPIV.1 (v0.2 ≡ v0.1) — pins the full chain
    # v0.2 A_2^(1) ≡ v0.1 PIV ≡ closed form.
    # -------------------------------------------------------------------------
    @testset "NYPIV.3 closed-form anchor (a,b)=(0,-2)" begin
        # The A_2^(1) parameter map sends α = (0,1,0) to (a,b) = (0,-2).
        α_seed = [0.0, 1.0, 0.0]
        pv_seed = _a2_to_piv_params(α_seed)
        @test pv_seed.a == 0
        @test pv_seed.b == -2

        # v0.1 PIV solve of (a,b)=(0,-2) reproduces the closed form
        # u(z) = -2z (piv_entire's documented (α,β); RF 2014 md:91).
        cfprob = piv_entire(:minus_2z; zspan = (1.0, 2.0))
        @test cfprob.params == (; α = 0, β = -2)   # matches the A_2^(1) map
        cfsol = solve_pade(cfprob; h = 0.1)
        for z in 1.0:0.25:2.0
            y, yp = cfsol(z)
            @test y  ≈ -2 * z atol = 1e-9          # closed form u = -2z
            @test yp ≈ -2     atol = 1e-9          # u' = -2
        end
    end

    # -------------------------------------------------------------------------
    # NYPIV.4 — mutation-proof.  Each mutation below was applied to the
    # parameter map / reduction *in this file*, the suite re-run
    # standalone, then restored.  All bit; the file ships GREEN at
    # 37/37.  Bite counts are total failing `@test`s observed.
    #
    # M1 — wrong dependent-variable scale `_f1_to_y(f1) = sqrt(2)*f1`
    #      → `f1` (i.e. drop the √2, the "y = -f_1/c" k=3 misreading):
    #      17 failures (NYPIV.1 trajectory + far-endpoint @tests and
    #      NYPIV.2 residual @tests all RED).
    # M2 — wrong b relation `b = -2*α[2]^2` → `b = -2*α[2]` in
    #      `_a2_to_piv_params`:  18 failures — NYPIV.1, NYPIV.2, and the
    #      "parameter map sanity" @test for `b` all RED.  (a is
    #      unchanged, so NYPIV.3's α_seed = (0,1,0) still gives
    #      b = -2·1 = -2 — M2 correctly does not bite NYPIV.3: at
    #      α_1 = 1 the two formulas coincide.  The generic
    #      α = [.5,.3,.2] case is what discriminates, and it bites.)
    # M3 — wrong independent-variable map `_tau_to_t(τ) = -τ/sqrt(2)`
    #      → `-τ` (the NY1999 "t → -x/c" misreading, dropping the √2):
    #      17 failures — NYPIV.1 (the v0.1 PIV span and sample t's are
    #      wrong) and NYPIV.2 (the 4ty² and 2(t²−a)y terms see the
    #      wrong t) both RED.
    # M4 — wrong derivative map `_f1prime_to_yprime(f1p) = -2*f1p`
    #      → `-f1p`:  17 failures — NYPIV.1 derivative-agreement @tests
    #      and NYPIV.2 residual @tests (y' enters via (y')²/(2y)) RED.
    #      Confirms the derivative leg is load-bearing.
    # -------------------------------------------------------------------------

end
