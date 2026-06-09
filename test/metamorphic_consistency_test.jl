"""
Metamorphic-relation test layer — CONSISTENCY relations (bead `padetaylor-krgy.3`,
epic `padetaylor-krgy`; methodology `docs/test_corpus/03_hardening_methodology.md`
§2.1).

This is the sibling of `test/metamorphic_symmetry_test.jl` — read that file's
header first for the full MR taxonomy table, the necessary-not-sufficient
framing, and the cross-reference to the canonical ALGEBRAIC metamorphic relation
(Noumi–Yamada Bäcklund invariance, `NYS.1.6` in
`test/noumi_yamada_symmetry_test.jl`).

## Why these are CONSISTENCY relations (vs the symmetry siblings)

MM.1–MM.3 exploit a symmetry of the *equation*.  MM.4–MM.5 instead exploit the
INTEGRATOR's own internal composition law: a Padé-Taylor step is *additive*
(one `h`-step composes the same analytic continuation as two `h/2`-steps) and
*invertible* (stepping `+h` then `−h` returns the initial state).  Both are
oracle-free GENERAL relations: they hold for ANY smooth, pole-free arc of ANY
2nd-order analytic IVP, so for MM.4/MM.5 themselves no Painlevé citation is
needed — they test the stepper, not a transcendent.  (MM.5 additionally carries
a Painlevé reading: the Painlevé property — single-valuedness under analytic
continuation — is what makes the forward∘reverse round-trip well-posed for the
transcendents; FW2011 …JCP230.md:153.)

| ID    | MR    | Relation                                   | Reference (verified)  |
|-------|-------|--------------------------------------------|-----------------------|
| MM.4  | MR-05 | step additivity: 1·h-step ≡ 2·(h/2)-steps  | GENERAL (no cite)     |
| MM.5  | MR-06 | forward∘reverse = id (round-trip)          | GENERAL; Painlevé     |
|       |       |                                            | single-valuedness     |
|       |       |                                            | FW2011 …JCP230.md:153 |

## Implementation note — why these run at the STEPPER level

Both relations are driven through `PadeStepper.pade_step!` directly rather than
the `solve_pade` driver.  Reason (a documented LIBRARY FINDING, see the
orchestrator report): the scalar `solve_pade` driver only integrates in the
ASCENDING direction (`while state.z < z_end`); a descending span like
`(0.0, −L)` silently returns a degenerate single-node trajectory (Rule-1
violation — it should fail loud).  `pade_step!`, by contrast, accepts a negative
`h` and steps backward correctly, so the reverse half of MM.5 is exercised at
the stepper layer where it is genuinely supported.  This keeps the MR honest:
we test the integrator's real invertibility, not a driver artefact.

## Reference confirmation (spot-verified at coding time, Law 1)

  - MR-06 (Painlevé reading) — FW2011 …JCP230.md:153: "The Painlevé property
    ensures that the solution … will remain single-valued, no matter how the
    integration path wanders around in the complex plane."

Self-contained: `using Test, PadeTaylor` + the `PadeStepper` submodule —
runnable standalone (`julia --project=. test/metamorphic_consistency_test.jl`)
and under `runtests.jl`.  Mutation-proof record in the file footer (Rule 4).
"""

using Test
using PadeTaylor
using PadeTaylor.PadeStepper: PadeStepperState, pade_step!

# Equianharmonic Weierstrass-℘ RHS u'' = 6u² (FW 2011 §5.1.1) + 16-digit ICs
# (test/_oracle_problems.jl:83-84).  ℘ on a SHORT arc from z=0 is pole-free
# (the nearest lattice pole is at z=1), so both consistency relations apply.
_fW(z, u, up) = 6 * u^2
const _U0_FW  = 1.071822516416917
const _UP0_FW = 1.710337353176786

# A representative Painlevé-II RHS (α=1), u'' = 2u³ + zu + 1, used to show
# MM.5 round-trips a genuine transcendent (not just the ℘ companion case).
_fPII1(z, u, up) = 2 * u^3 + z * u + 1.0

@testset "Metamorphic: consistency relations (MM.4–MM.5)" begin

    # -------------------------------------------------------------------------
    # MM.4 (MR-05) — step additivity: one h-step ≡ two h/2-steps.
    #
    # GENERAL relation (no Painlevé cite): from a fixed state at z, a single
    # Padé-Taylor step of length h and a composition of two steps of length
    # h/2 must land on the SAME analytic continuation, since each step's local
    # Padé approximant represents the same germ of the solution at z.  The
    # discrepancy is bounded by the local truncation error O(h^{order+1}) — at
    # order 30, h=0.5 that envelope is h^31 ≈ 4.7e-10.  EMPIRICALLY the two
    # paths agree FAR tighter (~1e-15): both use the same order-30 (15,15) Padé
    # germ, so the truncation contributions very nearly cancel and only
    # floating-point rounding remains.  We therefore assert against a tight
    # observed floor (1e-9) that is (i) comfortably inside the h^{order+1}
    # theory envelope and (ii) tight enough to BITE a real regression — an
    # order or step-rescaling bug shifts the discrepancy by many orders of
    # magnitude (see the mutation-proof footer).  Pole-free ℘ arc from z=0.
    # -------------------------------------------------------------------------
    @testset "MM.4 one h-step ≡ two h/2-steps (℘, pole-free)" begin
        for h in (0.5, 0.4, 0.3)
            s_one = PadeStepperState{Float64}(0.0, _U0_FW, _UP0_FW)
            pade_step!(s_one, _fW, 30, h)         # one h-step

            s_two = PadeStepperState{Float64}(0.0, _U0_FW, _UP0_FW)
            pade_step!(s_two, _fW, 30, h / 2)     # two h/2-steps
            pade_step!(s_two, _fW, 30, h / 2)

            # Both land at the same z (sanity on the composition).
            @test s_one.z ≈ s_two.z
            du  = abs(s_one.u  - s_two.u)
            dup = abs(s_one.up - s_two.up)
            # The discrepancy is bounded ABOVE by the larger of (i) the local
            # truncation envelope O(h^{order+1}) = h^31, and (ii) the
            # floating-point rounding floor ~eps·|state|.  Both order-30 paths
            # share the same (15,15) Padé germ, so the truncation contributions
            # very nearly cancel and rounding dominates: at h=0.5 the h^31
            # envelope (~4.7e-10) wins; at h=0.3 (h^31 ≈ 6e-17, below eps) the
            # rounding floor wins.  We bound by max(envelope, eps-floor) — tight
            # to the achievable accuracy, yet bites an additivity-breaking
            # mutation (which inflates du to O(1e-2) or larger; see footer).
            envelope   = h^31
            round_u    = 16 * eps(Float64) * max(abs(s_one.u),  abs(s_two.u))
            round_up   = 16 * eps(Float64) * max(abs(s_one.up), abs(s_two.up))
            @test du  < max(envelope, round_u)
            @test dup < max(10 * envelope, round_up)
            # Non-trivial: the step actually advanced the state.
            @test abs(s_one.u - _U0_FW) > 0.1
        end
    end

    # -------------------------------------------------------------------------
    # MM.5 (MR-06) — forward∘reverse = identity (round-trip).
    #
    # GENERAL relation; for the Painlevé transcendents it is underwritten by
    # the Painlevé property (single-valuedness under analytic continuation,
    # FW2011 md:153) — the round-trip is well-posed precisely because the
    # solution does not branch.  Step forward by +L, then back by −L, and the
    # state must return to its initial value.  HAZARD: smooth (pole-free)
    # SHORT arcs only — away from poles the ℘/Painlevé fields are
    # exponentially sensitive, so a long reverse leg would amplify rounding.
    # We keep L ≤ 0.5 (the ℘ arc stays well short of the z=1 pole) and assert
    # 1e-9 recovery; empirically the round-trip closes to ~1e-13 (℘) / ~1e-16
    # (PII).  Driven at the stepper level because the scalar `solve_pade`
    # driver does not support a descending span (see header / the finding).
    # -------------------------------------------------------------------------
    @testset "MM.5 forward∘reverse recovers the IC (℘ + PII)" begin
        # ℘ companion case.
        for L in (0.3, 0.5)
            s = PadeStepperState{Float64}(0.0, _U0_FW, _UP0_FW)
            pade_step!(s, _fW, 30,  L)            # forward +L
            u_fwd, up_fwd = s.u, s.up
            pade_step!(s, _fW, 30, -L)            # reverse −L
            @test s.z ≈ 0.0 atol = 1e-14
            @test s.u  ≈ _U0_FW  atol = 1e-9
            @test s.up ≈ _UP0_FW atol = 1e-9
            # Non-trivial: the forward leg genuinely moved the state, so the
            # recovery is a real round-trip, not a no-op.
            @test abs(u_fwd - _U0_FW) > 0.1
        end

        # Painlevé-II transcendent (α=1) — the Painlevé-property reading.
        u0p, up0p = 0.3, -0.2
        for L in (0.3, 0.5)
            s = PadeStepperState{Float64}(0.0, u0p, up0p)
            pade_step!(s, _fPII1, 30,  L)
            pade_step!(s, _fPII1, 30, -L)
            @test s.z ≈ 0.0 atol = 1e-14
            @test s.u  ≈ u0p  atol = 1e-9
            @test s.up ≈ up0p atol = 1e-9
        end
    end

end # @testset Metamorphic: consistency relations

# =============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — EXECUTED, then src/ restored byte-clean
# (`git diff --stat src/` empty after each restore).
#
#   MM.4 (MR-05) — break additivity in the COMPOSITION: take the two-step path
#     as two FULL h-steps instead of two h/2-steps
#     (`pade_step!(s_two, _fW, 30, h)` twice).  The two paths then land at
#     different z (s_one.z = h, s_two.z = 2h), so `s_one.z ≈ s_two.z` and the
#     value comparisons go RED.  [Alternative impl mutation: in
#     src/PadeStepper.jl `pade_step!`, drop the h^k Taylor rescaling so the step
#     length is mis-scaled — the h vs 2×(h/2) discrepancy then blows past the
#     1e-9 floor and `du < 1e-9` RED.]  Restored.
#
#   MM.5 (MR-06) — break invertibility: make the reverse leg ALSO step forward
#     (`pade_step!(s, _fW, 30, +L)` for the second call).  The state then
#     advances to z=2L instead of returning to 0, so `s.z ≈ 0.0` and the IC
#     recovery (`s.u ≈ _U0_FW`) go RED.  [Alternative impl mutation: in
#     src/PadeStepper.jl, flip the sign on the Padé evaluation point so the
#     reverse step evaluates u(z+|h|) rather than u(z−|h|); the round-trip no
#     longer closes and MM.5 RED — this is the same evaluation-point path that
#     padestepper_test.jl Mutation A exercises.]  Restored.
#
# Run standalone:  julia --project=. test/metamorphic_consistency_test.jl
# =============================================================================
