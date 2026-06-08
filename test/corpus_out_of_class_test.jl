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
# THE EMPIRICAL v1ub FINDING  (CONFIRMED — silent bridge)
# -------------------------------------------------------
# Driving the REAL solve_pade (NOT an mpmath emulation) settles the open
# question for bead padetaylor-v1ub.  VERDICT: **v1ub is CONFIRMED.** The
# solver SILENTLY returns finite, plausible-looking values with NO throw and
# NO NaN as it integrates toward (and, for a single large step, ACROSS) z=0.
# Measured relerr-vs-distance-to-0 curve, real-increasing march z0=-1 → z<0,
# fixed h=0.1 (oracle e^{1/z}, mpmath dps=50):
#
#     |z|=0.5  relerr 2.1e-16   (machine precision — far from 0, accurate)
#     |z|=0.2  relerr 3.1e-14
#     |z|=0.1  relerr 4.0e-10   (still excellent)
#     |z|=0.05 relerr 4.9e-02   (degrading — the jet is going out-of-class)
#     |z|=0.02 relerr 1.2e+17   (TOTAL GARBAGE, yet u is FINITE and even has
#                                the WRONG SIGN: e^{1/z}>0 always, but the
#                                solver returns u≈-2.4e-5 with no complaint)
#
# So accuracy collapses near the essential singularity while the OUTPUT stays
# finite — a pure Rule-1 silent lie.  ACROSS z=0: a single large step
# (z0=-0.2 → +0.2, h=0.4) "bridges" the essential singularity to ~6 digits
# (relerr 5.8e-6) and returns finite — the headline "bridges across 0" claim,
# confirmed.  (Smaller across-0 steps DO sometimes throw a TaylorSeries Inf-
# coefficient ArgumentError or return NaN, so the across-0 regime is step-
# dependent; the robust, step-independent confirmation is the APPROACH side,
# which is what the auto-flip marker below exercises.)
#
# Because v1ub is CONFIRMED, this file ships ONE `@test_broken` (CFail.1d): a
# `refuses(...)` helper that is currently `false` (the solver does NOT refuse)
# ⇒ "broken", and auto-flips to "unexpectedly passes" the day an out-of-class
# diagnostic (coefficient-growth / Padé-defect monotone drift) makes the
# solver throw.  Every other assertion @test-pins the documented REALITY.
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
#   * the COLLAPSE pins (z=-0.05, -0.02): asserted as INEQUALITIES on relerr
#     (>1e-3 at -0.05, >1.0 at -0.02) — we pin the *failure*, not a value, plus
#     `isfinite` (the silent-lie signature: finite-but-wrong).  These are
#     structural, not fitted.
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
# Avoids the "@test_broken throws → errors the expression" trap: we return a
# clean Bool.  Currently `false` (no out-of-class guard) ⇒ the @test_broken is
# "broken"; flips to "unexpectedly passes" when a future diagnostic lands.
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

        prob2 = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.05); order = 30)
        sol2  = solve_pade(prob2; h_max = 0.1)
        u_h2, _ = sol2(-0.05)
        relerr = abs(u_h2 - CFAIL_U_AT_NEG0p05) / abs(CFAIL_U_AT_NEG0p05)
        @test isfinite(u_h2)         # STILL finite — the silent-lie signature
        @test relerr > 1e-3          # but already DEGRADED (measured ~4.9e-2)
    end

    # ----------------------------------------------------------------------
    # CFail.1c  SILENT-LIE collapse.  Within ~0.02 of z=0 the returned value
    # is TOTAL GARBAGE (relerr ≫ 1) yet FINITE — and even the SIGN is wrong
    # (e^{1/z} > 0 for all real z, but the solver returns u < 0).  This is the
    # core Rule-1 violation: a confident, finite, plausible-looking lie.
    # ----------------------------------------------------------------------
    @testset "CFail.1c  silent finite-but-wrong lie near z=0 (z=-0.02)" begin
        prob = PadeTaylorProblem(fess, (U0, UP0), (Z0, -0.02); order = 30)
        sol  = solve_pade(prob; h_max = 0.1)
        u_h, _ = sol(-0.02)
        ex = uexact(-0.02)
        relerr = abs(u_h - ex) / abs(ex)
        @test isfinite(u_h)          # finite — no throw, no NaN
        @test relerr > 1.0           # but utterly wrong (measured ~1.2e17)
        @test u_h < 0                # WRONG SIGN — e^{1/z} is strictly positive
    end

    # ----------------------------------------------------------------------
    # CFail.1d  v1ub AUTO-FLIP MARKER (the only @test_broken in this file).
    # The solver does NOT refuse this out-of-class input on the approach side
    # — it silently returns a finite (wrong) value.  `refuses_out_of_class`
    # is currently `false` ⇒ @test_broken.  When a future out-of-class
    # diagnostic (coefficient-growth / Padé-defect monotone drift) makes the
    # solver throw, `refuses` flips to `true` and this @test_broken reports
    # "Unexpectedly passing" — the auto-flip that closes v1ub.
    # ----------------------------------------------------------------------
    @testset "CFail.1d  v1ub: solver SHOULD refuse out-of-class input" begin
        @test_broken refuses_out_of_class(-0.02; h = 0.1)
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
# M-oracle (CFail.1a): in test/_oracle_corpus_out_of_class.jl change the
#   leading digit of CFAIL_U_AT_NEG0p5 from 0.1353… to 0.1453… (a meaningful
#   ~7% flip, well above the 1e-13 tol).  Result: ONLY the CFail.1a `u_h`
#   assertion went RED (isapprox(0.13533…, 0.14533…; rtol=1e-13) = false);
#   all other 11 assertions stayed GREEN — the 1b/1c collapse pins use the
#   Float64 `uexact` helper, not this literal, so they are independent, and
#   the 1d @test_broken stayed broken.  Restored the literal; suite GREEN.
#
# M-impl (the gold standard): in src/PadeStepper.jl `_evaluate_pade`, change
#   the return `num / den` to `(num / den) * (1 + 1e-7)` — a 1e-7 relative
#   perturbation of every stored Padé value (the handle the CPB/CRic/CVrow
#   mutation-proofs use).  Result: exactly the FOUR load-bearing value/sign
#   assertions reddened (8 pass / 4 fail / 1 broken):
#     - CFail.1a  z=-0.5 (rtol 1e-13): 0.135335408… vs 0.135335283… RED;
#     - CFail.1a  z=-0.2 (rtol 1e-9):  0.006738187… vs 0.006737947… RED;
#     - CFail.1b  z=-0.1 (rtol 1e-8):  5.2906e-5    vs 4.5400e-5    RED;
#     - CFail.1c  sign pin `u_h < 0`:  the 1e-7 scaling tipped the garbage
#       value to +1.527, so the WRONG-SIGN assertion flipped RED.
#   The by-design non-quantitative pins stayed GREEN — the inequality relerr
#   pins (1b relerr>1e-3, 1c relerr>1.0) are already so far past threshold a
#   1e-7 nudge can't cross them; the 1e across-0 bridge (relerr<1e-4, measured
#   5.8e-6) is far enough below 1e-4 that the 1e-7 perturbation stays inside;
#   the 1d @test_broken (refuses==false) is structural.  Restored
#   src/PadeStepper.jl byte-for-byte (`git status --porcelain src/` empty);
#   suite GREEN (12 pass / 1 broken).
#
# Run standalone with:
#   julia --project=. test/corpus_out_of_class_test.jl
# ============================================================================
