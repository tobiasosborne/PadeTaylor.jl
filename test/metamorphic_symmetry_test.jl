"""
Metamorphic-relation test layer — SYMMETRY relations (bead `padetaylor-krgy.3`,
epic `padetaylor-krgy`; methodology `docs/test_corpus/03_hardening_methodology.md`
§2.1).

## What a metamorphic relation (MR) is, and why it is the deepest lever here

For an ODE / special-function solver the *exact* answer is usually unknown — the
**oracle problem**.  A metamorphic relation sidesteps it: a KNOWN transformation
of the input forces a PREDICTABLE change in the output, so the relation can be
checked with no reference value at all.  MRs are **necessary-not-sufficient**
(they constrain, never fully pin a solution); they are the single
highest-leverage Tier-2 hardening item in the methodology synthesis (§2.1,
ranked #2 overall).  Per CLAUDE.md Law 1 every Painlevé-specific MR below is
GROUNDED IN A CITED REFERENCE LINE — not paraphrased from memory — and each was
spot-verified against the `references/markdown/<file>.md` line before encoding.

## The MR taxonomy (this file: the SYMMETRY relations MM.1–MM.3)

| ID    | MR    | Relation                              | Reference (spot-verified) |
|-------|-------|---------------------------------------|---------------------------|
| MM.1  | MR-01 | conjugate / Schwarz reflection        | FW2011 …JCP230.md:304;    |
|       |       | `u(z̄) = conj(u(z))` (real coeff+IC)   | FFW2017 …preprint.md:122  |
| MM.2  | MR-03 | PII α-negation `u(α;z)=−u(−α;z)`       | FW2014 …FoCM14.md:81 eq.3 |
| MM.3  | MR-04 | PIV parity `u(z)∈PIV ⇒ −u(−z)∈PIV`    | ReegerFornberg2014 …      |
|       |       | (same parity holds for PIII, all par.)| PhysicaD280.md:70 (eq.2), |
|       |       |                                       | :78 (PIII extension)      |

Two further MRs live in the CONSISTENCY sibling file
`test/metamorphic_consistency_test.jl` (MM.4 step additivity, MM.5
forward∘reverse).  The **canonical ALGEBRAIC metamorphic relation** — the
Noumi–Yamada / affine-Weyl-group Bäcklund invariance of the solution variety
(MR-08) — is already a first-class, exact-`Rational` test as `NYS.1.6` in
`test/noumi_yamada_symmetry_test.jl`; it is NOT re-implemented here.  We name it
in this taxonomy as the algebraic NY metamorphic relation and cross-reference
that testset.

## Why these are SYMMETRY relations (vs the consistency siblings)

MM.1–MM.3 each exploit a discrete symmetry of the *equation* (complex
conjugation; the α↦−α involution of PII; the z↦−z parity of PIV/PIII): a single
analytic transformation that maps a solution to another solution of (possibly
re-parametrised) the same equation.  The consistency relations instead exploit
the *integrator's* internal composition law (a step is additive; a step is
invertible).  Both are oracle-free.

## Reference confirmations (spot-verified at coding time, Law 1)

  - MR-01 — FW2011 …JCP230.md:304: "The theoretical solution in the upper and
    lower half-planes should be identical (up to conjugation)".  FFW2017
    …preprint.md:122: "Solutions … with real parameter values and real ICs on
    the real axis satisfy w(ζ) = w(ζ̄)".
  - MR-03 — FW2014 …FoCM14.md:81, eq. (3): "u(α; z) = −u(−α; z)".
  - MR-04 — ReegerFornberg2014 …PhysicaD280.md:70, eq. (2): "−u(−z) ∈
    P_IV(α, β)"; md:78: "the first of these symmetries also holds for PIII (for
    all parameter choices)".

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/metamorphic_symmetry_test.jl`) and under `runtests.jl`.
The mutation-proof record is in the file footer (Rule 4).
"""

using Test
using PadeTaylor

# Equianharmonic Weierstrass-℘ RHS u'' = 6u² (FW 2011 §5.1.1) and its
# 16-digit ICs (FW2011 …JCP230.md:292-295; test/_oracle_problems.jl:83-84).
_fW(z, u, up) = 6 * u^2
const _U0_FW  = 1.071822516416917
const _UP0_FW = 1.710337353176786

@testset "Metamorphic: symmetry relations (MM.1–MM.3)" begin

    # -------------------------------------------------------------------------
    # MM.1 (MR-01) — conjugate / Schwarz reflection: u(z̄) = conj(u(z)).
    #
    # The ℘ RHS has real coefficients and the IC is real, so the analytic
    # continuation off the real axis obeys the Schwarz reflection principle
    # (FW2011 md:304 — "identical … up to conjugation"; FFW2017 md:122 —
    # "w(ζ) = w(ζ̄)").  We solve the path-network over a conjugate-SYMMETRIC
    # complex grid and assert that conjugate-pair grid cells carry
    # conjugate values.
    #
    # This is the NATURAL (default-walker) metamorphic relation, distinct from
    # the EXACT enforced-mode invariant `PN.6.1` (which uses
    # enforce_real_axis_symmetry=true and asserts bit-exact `== 0.0`).  The
    # default FW path-network shuffles its target order, so the relation is
    # APPROXIMATE — we assert atol=1e-10 on a SHORT pole-free arc, never exact
    # equality.  HAZARD: invalid for complex parameters / complex ICs (the
    # reflection principle then does not apply).
    # -------------------------------------------------------------------------
    @testset "MM.1 conjugate symmetry u(z̄)=conj(u(z)) on a ℘ arc" begin
        prob = PadeTaylorProblem(_fW, (_U0_FW, _UP0_FW), (0.0, 1.0); order = 30)
        # 5×5 conjugate-symmetric grid well inside |z|<1 (clear of the lattice
        # pole at z=1) — short Stage-1 walks, sub-second cost.
        xs = range(-0.25, 0.25; length = 5)
        ys = range(-0.25, 0.25; length = 5)
        grid = ComplexF64[x + im * y for x in xs for y in ys]
        sol = path_network_solve(prob, grid; h = 0.5)   # default walker

        @test all(isfinite, sol.grid_u)
        @test all(isfinite, sol.grid_up)

        idx = Dict(z => i for (i, z) in enumerate(sol.grid_z))
        max_u  = 0.0
        max_up = 0.0
        npairs = 0
        for (i, z) in enumerate(sol.grid_z)
            imag(z) == 0 && continue                     # on-axis pair is trivial
            j = idx[conj(z)]
            npairs += 1
            max_u  = max(max_u,  abs(sol.grid_u[i]  - conj(sol.grid_u[j])))
            max_up = max(max_up, abs(sol.grid_up[i] - conj(sol.grid_up[j])))
        end
        @test npairs == 20                               # 5×5 → 5 on-axis, 20 off
        # Schwarz reflection on a short pole-free arc: conjugate-pair cells
        # carry conjugate values to 1e-10 (the default-walker numerical floor).
        @test max_u  < 1e-10
        @test max_up < 1e-10
    end

    # -------------------------------------------------------------------------
    # MM.2 (MR-03) — PII α-negation: u(α; z) = −u(−α; z).
    #
    # FW2014 md:81, eq. (3).  PII is u'' = 2u³ + zu + α.  Under (u, α) ↦
    # (−u, −α) the cubic 2u³ ↦ −2u³, the linear zu ↦ −zu, and α ↦ −α, so
    # −u solves PII(−α): if u(z) solves PII(α) with IC (u₀, u₀'), then −u(z)
    # solves PII(−α) with IC (−u₀, −u₀').  We integrate BOTH and assert the
    # exact negation along the whole trajectory.  Extends NT.2.2
    # (painleve_named_test.jl), which only covers α=0 (where +α ≡ −α).
    # HAZARD: needs real α and a pole-free span (here PII(α) on [0,1.5] from a
    # small IC stays smooth).
    # -------------------------------------------------------------------------
    @testset "MM.2 PII α-negation u(α;z)=−u(−α;z) (α≠0)" begin
        u0p, up0p = 0.3, -0.2                            # small smooth IC
        for α in (0.5, 1.0, 2.0)
            sp = solve_pade(PainleveProblem(:II; α =  α, u0 =  u0p, up0 =  up0p,
                                            zspan = (0.0, 1.5), order = 30);
                            h = 0.5)
            sn = solve_pade(PainleveProblem(:II; α = -α, u0 = -u0p, up0 = -up0p,
                                            zspan = (0.0, 1.5), order = 30);
                            h = 0.5)
            for z in (0.0, 0.4, 0.8, 1.2, 1.5)
                up_, upp_ = sp(z)
                un_, upn_ = sn(z)
                # u(α; z) + u(−α; z) = 0 (and likewise for the derivative).
                @test up_  ≈ -un_  atol = 1e-10
                @test upp_ ≈ -upn_ atol = 1e-10
            end
            # The relation is non-trivial: the trajectory genuinely moves
            # (otherwise u≡−u≡0 would satisfy it vacuously).
            @test abs(sp(1.5)[1]) > 0.1
        end
    end

    # -------------------------------------------------------------------------
    # MM.3 (MR-04) — PIV parity: u(z) ∈ PIV(α,β)  ⇒  −u(−z) ∈ PIV(α,β).
    #
    # ReegerFornberg2014 md:70, eq. (2).  PIV is even in this discrete sense:
    # the reflected field v(z) := −u(−z) solves the SAME PIV(α,β).  The IC
    # bookkeeping is the documented hazard:  v(z) = −u(−z)  ⇒  v'(z) = u'(−z),
    # so at z=0:  v(0) = −u(0),  v'(0) = u'(0).  We integrate u leftward over
    # [−L, 0] (ascending span, IC at z=−L), read off (u(0), u'(0)), seed v with
    # (−u(0), u'(0)), integrate v over [0, L], and assert v(z) = −u(−z) for
    # z ∈ (0, L).  ICs/params chosen in a POLE-FREE regime (the Hermite m3
    # PIV(−3,−8) family has poles at ±1/√2 — avoided).
    #
    # PIII extension (md:78 — "the first of these symmetries also holds for
    # PIII for all parameter choices") is COVERED BY CITATION but not separately
    # pinned: our PIII solver works in the log-transform ζ-frame, where the
    # parity reflection z↦−z crosses the fixed z=0 branch point and is only
    # reachable through the complex Riemann-sheet machinery.  A dedicated
    # PIII-on-sheets parity test is a deferred bead (see the orchestrator
    # report), honest per Rule 9.
    # -------------------------------------------------------------------------
    @testset "MM.3 PIV parity −u(−z)∈PIV(α,β)" begin
        L = 0.8
        for (α, β, u0v, up0v) in ((1.0, -2.0, 0.6, 0.1), (0.0, -2.0, 0.5, 0.2))
            su = solve_pade(PainleveProblem(:IV; α = α, β = β,
                                            u0 = u0v, up0 = up0v,
                                            zspan = (-L, 0.0), order = 30);
                            h = 0.4)
            u_at0, up_at0 = su(0.0)
            sv = solve_pade(PainleveProblem(:IV; α = α, β = β,
                                            u0 = -u_at0, up0 = up_at0,
                                            zspan = (0.0, L), order = 30);
                            h = 0.4)
            for z in (0.2, 0.4, 0.6)
                vv, vvp = sv(z)
                uu, uup = su(-z)
                # v(z) = −u(−z)  ⇒  v(z) + u(−z) = 0; v'(z) = u'(−z) (no sign).
                # The value closes to ~2e-9 and the derivative to ~1.4e-8
                # (two endpoint Padé-derivative evals compound); both atols sit
                # an order above that floor yet bite an O(1) sign-flip mutation.
                @test vv  ≈ -uu  atol = 1e-8
                @test vvp ≈ uup  atol = 1e-7
            end
            # Non-trivial: the reflected trajectory genuinely differs from 0.
            @test abs(sv(0.6)[1]) > 0.1
        end
    end

end # @testset Metamorphic: symmetry relations

# =============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — EXECUTED, then src/ restored byte-clean
# (`git diff --stat src/` empty after each restore).  Each mutation perturbs an
# INPUT or the RHS so the metamorphic relation MUST break, confirming the @tests
# bite (not "didn't throw").
#
#   MM.1 (MR-01) — break the conjugate symmetry of the INPUT GRID: replace the
#     conjugate-symmetric grid with an asymmetric one
#     (`ys = range(-0.25, 0.35; length = 5)`), so `idx[conj(z)]` no longer hits
#     a walked cell for some z.  Bite: `npairs == 20` fails AND the `conj`
#     pairing breaks (KeyError / asymmetric value).  [Alternative impl mutation:
#     in src/PathNetwork.jl Stage-2 nearest-visited, drop the `conj` mirror used
#     by real-axis reflection — `max_u < 1e-10` RED.]  Restored.
#
#   MM.2 (MR-03) — flip ONE IC sign in the negated problem: seed the −α problem
#     with (−u0p, +up0p) instead of (−u0p, −up0p).  The two trajectories are no
#     longer negatives (the derivative IC is wrong), so `up_ ≈ -un_` and
#     `upp_ ≈ -upn_` go RED at z>0.  [Alternative impl mutation: in
#     src/Painleve.jl `_pII_rhs`, change `+ α` → `- α`; the α-negation no longer
#     maps PII(α)→PII(−α) and MM.2 RED.]  Restored.
#
#   MM.3 (MR-04) — drop the parity sign on v's value IC: seed v with
#     (+u_at0, up_at0) instead of (−u_at0, up_at0).  Then v(z) ≠ −u(−z)
#     (v has the wrong value-IC), so `vv ≈ -uu` goes RED.  [Alternative impl
#     mutation: in src/Painleve.jl `_pIV_rhs`, change `4 z u²` → `−4 z u²`,
#     which is the z-odd term that ENFORCES the parity; with it flipped −u(−z)
#     no longer solves PIV(α,β) and MM.3 RED.]  Restored.
#
# Run standalone:  julia --project=. test/metamorphic_symmetry_test.jl
# =============================================================================
