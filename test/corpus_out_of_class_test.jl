# test/corpus_out_of_class_test.jl
#
# ============================================================================
# Ground-truth corpus: OUT-OF-CLASS / FAIL-LOUD — essential singularity
# (bead padetaylor-nzcj, epic padetaylor-25og; bug padetaylor-v1ub).
#
# WHAT THIS FILE PINS
# -------------------
# The package is MEROMORPHIC-ONLY.  Rule 1 (fail loud, never silently lie)
# requires that an out-of-class input — one whose solution has a NON-pole
# singularity — be refused with a throw/NaN rather than answered with a
# finite, plausible-looking lie.  The cleanest out-of-class probe (Family G,
# docs/test_corpus/02_corpus_extension_plan.md:292-306) is the ESSENTIAL
# singularity:
#
#     u'' = u·(1 + 2z) / z^4        exact solution   u(z) = e^{1/z}
#
# z = 0 is an essential singularity of e^{1/z} (Laurent series 1 + 1/z +
# 1/(2z²) + … has infinitely many negative-power terms — NOT a pole).  The
# local Taylor jet's coefficients grow super-geometrically (|c_k|^{1/k} drifts
# upward without bound); nothing in the solver inspects that drift, so the
# meromorphic-only contract has NO guard against this input.
#
# THE EMPIRICAL v1ub FINDING  (CONFIRMED — then RESOLVED by ADR-0033)
# -------------------------------------------------------------------
# Driving the REAL solve_pade (NOT an mpmath emulation) settled the open
# question for bead padetaylor-v1ub.  VERDICT: **v1ub was CONFIRMED.** The
# solver SILENTLY returned finite, plausible-looking values with NO throw and
# NO NaN as it integrated toward z=0.  Measured relerr-vs-distance-to-0 curve,
# real-increasing march z0=-1 → z<0, fixed h=0.1 (oracle e^{1/z}, dps=50):
#
#     |z|=0.5  relerr 2.1e-16   (machine precision — far from 0, accurate)
#     |z|=0.2  relerr 3.1e-14
#     |z|=0.1  relerr 4.0e-10   (still excellent)
#     |z|=0.05 relerr 4.9e-02   (degrading — the jet is going out-of-class)
#     |z|=0.02 relerr 1.2e+17   (TOTAL GARBAGE, yet u was FINITE and even had
#                                the WRONG SIGN: e^{1/z}>0 always, yet the
#                                solver returned u≈-2.4e-5 with no complaint)
#
# Accuracy collapsed near the essential singularity while the OUTPUT stayed
# finite — a pure Rule-1 silent lie.
#
# RESOLUTION (bug padetaylor-v1ub, ADR-0033).  solve_pade now ENFORCES its
# meromorphic-only contract by default (`check_in_class = true`): the driver
# watches the per-step two-order Padé convergence defect δ and throws an
# `OutOfClassError` once δ exceeds the calibrated threshold τ=1e-3 AND has
# grown monotonically over the last K=2 steps — the de-Montessus /
# Nuttall–Pommerenke signature of a jet leaving the meromorphic class (GGT
# 2013 §8; src/OutOfClass.jl).  Consequences pinned below:
#
#   - The APPROACH now THROWS once δ runs away.  Measured: solve_pade to
#     z=-0.5/-0.2/-0.1/-0.05 still returns finite accurate/degraded values
#     (δ at/below the floor — the guard is silent); solve_pade to z=-0.02
#     THROWS OutOfClassError (δ climbs 2.7e-14 → 3.2e-10 → 0.99, the gate
#     fires).  So CFail.1a/1b keep their finite pins; CFail.1c flips from
#     pinning the silent lie to asserting the throw.
#   - The across-0 single big step (CFail.1e) CANNOT trip the guard — one
#     step leaves no ≥2-step monotone history — so it still bridges to ~6
#     digits and returns finite, exactly as before (the gate is the safety).
#   - CFail.1d, the auto-flip marker (`refuses_out_of_class`), now returns
#     `true`, so its former `@test_broken` is flipped to a passing `@test`.
#
# This is a legitimate behaviour change (silent-lie → fail-loud), NOT a
# tolerance relaxation: every reframed assertion is documented inline with
# its v1ub/ADR-0033 citation, and the mutation-proof footer was re-run.
#
# CFail.2 Chazy: DEFERRED (no far-side oracle past a movable natural boundary,
# plan Family G); not shipped here.  See report / follow-on bead.
#
# TOLERANCES (justified by the EMPIRICAL regime, never fitted to pass)
# -------------------------------------------------------------------
#   * far-from-0 value pins (z=-0.5, -0.2, -0.1): rtol 1e-13 / 1e-9 — these are
#     the regime where the (15,15) Float64 Padé floor still holds because the
#     jet has not yet sensed the essential singularity (relerr measured 2e-16…
#     4e-10).  The z=-0.1 pin is loosened to 1e-8 (measured 4e-10) — the last
#     point before the out-of-class collapse.
#   * the DEGRADATION pin (z=-0.05): asserted as an INEQUALITY on relerr
#     (>1e-3) — we pin the *degradation*, not a value, plus `isfinite` (the
#     guard is still silent here: δ has not yet run away).  Structural, not
#     fitted.
#   * the COLLAPSE point (z=-0.02): post-ADR-0033 this no longer returns a
#     finite lie — it THROWS OutOfClassError.  CFail.1c asserts the throw
#     (`@test_throws`), the fail-loud replacement for the old finite-lie pin.
#   * across-0 bridge (h=0.4): relerr < 1e-4 AND isfinite — pins the silent
#     bridge to ≥4 digits (measured 5.8e-6).
#
# REFERENCES: docs/test_corpus/02_corpus_extension_plan.md:292-306 (Family G);
# bug padetaylor-v1ub.  Oracle constants:
# external/probes/corpus-oracles/out-of-class/capture.py (mpmath dps=50,
# sympy residual==0 + essential-singularity GATE).
#
# MUTATION-PROOF (Rule 4): see the comment block at the end of this file —
# one M-oracle (flip a pinned e^{1/z} value) and one M-impl (perturb the stored
# Padé in src/PadeStepper.jl) were ACTUALLY executed RED and restored
# byte-for-byte (`git diff src/` empty).
# ============================================================================

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_corpus_out_of_class.jl"))

# Essential-singularity 2nd-order recast, u'' = u(1+2z)/z^4, exact u = e^{1/z}.
fess(z, u, up) = u * (1 + 2z) / z^4

# Exact closed-form oracle helpers (Float64 — values are e^{1/z}; the high-dps
# pins above are the independent ground truth, these mirror them for the curve).
uexact(z)  = exp(1 / z)

# IC at z0 = -1 (real): u(-1)=e^{-1}, u'(-1)=-e^{-1}/(-1)^2 = -e^{-1}.
const Z0   = -1.0
const U0   = exp(-1.0)
const UP0  = -exp(-1.0)

# v1ub AUTO-FLIP HELPER.  Boolean: does the solver REFUSE this out-of-class
# input (throw, or return a non-finite sentinel) rather than silently lie?
# Avoids the "@test throws → errors the expression" trap: we return a clean
# Bool.  Post-ADR-0033 the out-of-class guard makes solve_pade throw on the
# approach, so this returns `true` — CFail.1d below is now a passing @test
# (it was @test_broken while v1ub was open).
function refuses_out_of_class(zend; h = 0.1)
    try
        prob = PadeTaylorProblem(fess, (U0, UP0), (Z0, zend); order = 30)
        sol  = solve_pade(prob; h_max = h)
        u, _ = sol(zend)
        return !isfinite(u)          # NaN/Inf would also be "did not silently lie"
    catch
        return true                  # any throw = fail-loud = refused
    end
end

@testset "Corpus: out-of-class / fail-loud (CFail)" begin

    # ----------------------------------------------------------------------
    # CFail.1a  FAR from the essential singularity the solver is ACCURATE.
    # The local jet has not yet sensed z=0; the (15,15) Float64 Padé floor
    # holds.  This is the upper, well-behaved arm of the error-inflation
    # curve — pinned to the exact e^{1/z} (mpmath dps=50).
    # ----------------------------------------------------------------------
    @testset "CFail.1a  accurate far from z=0 (z=-0.5, -0.2)" begin
        prob = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.5); order = 30)
        sol  = solve_pade(prob; h_max = 0.1)
        u_h, _ = sol(-0.5)
        @test isfinite(u_h)
        @test isapprox(u_h, CFAIL_U_AT_NEG0p5; rtol = 1e-13)

        prob2 = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.2); order = 30)
        sol2  = solve_pade(prob2; h_max = 0.1)
        u_h2, _ = sol2(-0.2)
        @test isapprox(u_h2, CFAIL_U_AT_NEG0p2; rtol = 1e-9)
    end

    # ----------------------------------------------------------------------
    # CFail.1b  ERROR-INFLATION onset.  At z=-0.1 the solver is still good
    # (~4e-10); by z=-0.05 it has DEGRADED past 1e-3 — the jet is going
    # out-of-class.  We pin the last-good value AND the onset of degradation.
    # ----------------------------------------------------------------------
    @testset "CFail.1b  error inflation onset (z=-0.1 good, z=-0.05 degraded)" begin
        prob = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.1); order = 30)
        sol  = solve_pade(prob; h_max = 0.1)
        u_h, _ = sol(-0.1)
        @test isfinite(u_h)
        @test isapprox(u_h, CFAIL_U_AT_NEG0p1; rtol = 1e-8)   # measured 4.0e-10

        # z=-0.05: STILL finite and accurate-enough that the out-of-class
        # guard (ADR-0033) does NOT fire — δ has not yet run away monotonically
        # past τ at this distance — so solve_pade returns the degraded-but-
        # finite value (the guard is correctly silent here; it only throws
        # once the δ-runaway is sustained, which happens between -0.05 and
        # -0.02, see CFail.1c).
        prob2 = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.05); order = 30)
        sol2  = solve_pade(prob2; h_max = 0.1)   # guard silent: no throw at -0.05
        u_h2, _ = sol2(-0.05)
        relerr = abs(u_h2 - CFAIL_U_AT_NEG0p05) / abs(CFAIL_U_AT_NEG0p05)
        @test isfinite(u_h2)         # finite — guard not yet fired
        @test relerr > 1e-3          # but already DEGRADED (measured ~4.9e-2)
    end

    # ----------------------------------------------------------------------
    # CFail.1c  FAIL-LOUD collapse (post-ADR-0033).  Within ~0.02 of z=0 the
    # jet has decisively left the meromorphic class: the two-order Padé
    # convergence defect δ has grown monotonically (2.7e-14 → 3.2e-10 → 0.99)
    # past τ, so solve_pade now THROWS OutOfClassError instead of returning
    # the old finite, wrong-signed lie (relerr ~1.2e17, u<0 for the strictly
    # positive e^{1/z}).  This is the v1ub fix: a confident finite lie is
    # replaced by a fail-loud refusal (Rule 1).  We assert the THROW.
    # (Reframed from the old "silent finite-but-wrong lie" pin — bug
    # padetaylor-v1ub, ADR-0033; a behaviour change, not a tolerance relax.)
    # ----------------------------------------------------------------------
    @testset "CFail.1c  fail-loud refusal near z=0 (z=-0.02 throws)" begin
        prob = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.02); order = 30)
        # The driver throws inside solve_pade as the δ-runaway gate fires —
        # the integration never produces the old finite-wrong value at -0.02.
        @test_throws OutOfClassError solve_pade(prob; h_max = 0.1)

        # Escape hatch (ADR-0033): a user who knowingly probes out-of-class
        # can opt out; then the LEGACY silent-lie behaviour is recovered
        # verbatim — finite, relerr ≫ 1, wrong sign.  This pins both the
        # opt-out and the documented legacy reality it restores.
        sol  = solve_pade(prob; h_max = 0.1, check_in_class = false)
        u_h, _ = sol(-0.02)
        ex = uexact(-0.02)
        relerr = abs(u_h - ex) / abs(ex)
        @test isfinite(u_h)          # opt-out: finite — no throw, no NaN
        @test relerr > 1.0           # but utterly wrong (measured ~1.2e17)
        @test u_h < 0                # WRONG SIGN — e^{1/z} is strictly positive
    end

    # ----------------------------------------------------------------------
    # CFail.1d  v1ub AUTO-FLIP MARKER — now FLIPPED to a passing @test.
    # The solver now REFUSES this out-of-class input on the approach side: it
    # throws OutOfClassError as the δ-runaway gate fires, so
    # `refuses_out_of_class` returns `true`.  This was a `@test_broken` while
    # v1ub was open (the solver silently returned a finite wrong value); the
    # out-of-class diagnostic (two-order Padé convergence defect + monotone
    # gate, src/OutOfClass.jl / ADR-0033) closes v1ub and flips the marker.
    # ----------------------------------------------------------------------
    @testset "CFail.1d  v1ub: solver REFUSES out-of-class input (resolved)" begin
        @test refuses_out_of_class(-0.02; h = 0.1)     # was @test_broken; v1ub fixed
    end

    # ----------------------------------------------------------------------
    # CFail.1e  ACROSS-ZERO silent bridge.  A single large step from z0=-0.2
    # (u=e^{-5}) to z=+0.2 (u=e^{+5}) integrates straight ACROSS the essential
    # singularity at z=0 and returns a finite value good to ~6 digits — the
    # "bridges z=0" headline.  Pin: finite AND relerr < 1e-4.
    # (Smaller across-0 steps may instead throw / NaN — step-dependent — so we
    # pin the documented single-large-step bridge, not a general claim.)
    # ----------------------------------------------------------------------
    @testset "CFail.1e  silent bridge ACROSS z=0 (single large step)" begin
        prob = PadeTaylorProblem(fess, (exp(-5.0), -exp(-5.0) / 0.04),
                                 (-0.2, 0.2); order = 30)
        sol  = solve_pade(prob; h_max = 0.4)          # one segment spans [-0.2,0.2]
        u_h, _ = sol(0.2)
        relerr = abs(u_h - CFAIL_U_AT_POS0p2) / abs(CFAIL_U_AT_POS0p2)
        @test isfinite(u_h)          # bridged z=0 with no throw/NaN (silent)
        @test relerr < 1e-4          # to ~6 digits (measured 5.8e-6)
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
# The NEW load-bearing assertions are the fail-loud guard: CFail.1c
# (`@test_throws OutOfClassError`) and CFail.1d (the flipped `@test`).  Both
# were mutation-proven against the out-of-class guard in src/OutOfClass.jl
# (bug padetaylor-v1ub, ADR-0033).  TWO independent guard mutations, each of
# which DISABLES firing, were run RED and restored byte-for-byte.
#
# M-guard-τ: in src/OutOfClass.jl bump `OUT_OF_CLASS_TAU` 1.0e-3 → 1.0e10, so
#   no measured δ ever exceeds the threshold and the guard never fires.
#   Result: solve_pade no longer throws on the approach — exactly CFail.1c
#   (`@test_throws`) and CFail.1d (`@test refuses…`) reddened (12 pass / 2
#   fail).  The far-accurate pins (1a, 1b) and the across-0 bridge (1e) stayed
#   GREEN — they never depended on the throw.  Restored τ; suite GREEN (14/14).
#
# M-guard-monotone (the gold standard): in src/OutOfClass.jl `check_in_class!`
#   flip the monotone-growth test `history[i] > history[i-1]` → `<`, so a
#   GROWING δ sequence never satisfies the gate and the guard never fires.
#   Result: identical reddening — only CFail.1c + CFail.1d went RED (12 pass /
#   2 fail).  This proves the HISTORY GATE itself (not just the threshold) is
#   load-bearing: without monotone-growth detection the e^{1/z} runaway is not
#   caught.  Restored the comparison byte-for-byte (`git status --porcelain
#   src/OutOfClass.jl` empty); suite GREEN (14/14).
#
# Why two mutations: the guard fires on `δ > τ AND monotone-growth-over-K`.
# M-guard-τ kills the first conjunct, M-guard-monotone the second; each alone
# reddens the same two assertions, confirming both halves of the AND are
# necessary and neither is dead.  The legacy M-oracle / M-impl value pins
# (CFail.1a value, CFail.1c opt-out sign) remain valid against the finite
# escape-hatch assertions and are still exercised by the suite.
#
# Run standalone with:
#   julia --project=. test/corpus_out_of_class_test.jl
# ============================================================================
