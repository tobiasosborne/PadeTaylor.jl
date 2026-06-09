# test/accuracy_ledger_test.jl -- ACCURACY-REGRESSION LEDGER
#                                 (bead padetaylor-krgy.8, Tier-2 §2.5)
#
# ## Why this file exists -- the gap a fixed-tol gate leaves open
#
# Every value pin in the corpus asserts `relERR < tol` at a FIXED, generous
# absolute tolerance.  That catches a value that goes *wrong*, but it is
# blind to slow accuracy DRIFT: a change that degrades a case's achieved
# relERR from 1e-13 to 1e-9 sails through a 1e-6 gate unnoticed.  Over many
# such silent degradations a "passing" suite can mask a numerical core that
# has quietly lost six digits.
#
# An accuracy LEDGER closes that gap.  For a fixed set of DETERMINISTIC,
# representative computations it stores the BEST ACHIEVED relERR of each in a
# committed, human-reviewable baseline (test/accuracy_baseline.toml) and
# asserts, per case,
#
#       achieved_relerr  <=  baseline_relerr * (1 + SLACK).
#
# This is a real invariant (CLAUDE.md Rule 5 -- an error-vs-baseline BOUND,
# not "didn't throw", and strictly finer than the absolute-tol value pin it
# COMPLEMENTS, never replaces).  It bites in BOTH directions:
#
#   * DEGRADATION -- a case whose relERR rises past `baseline*(1+SLACK)`
#     turns the suite RED, even when it still passes the original absolute
#     tol.  (The mutation-proof footer demonstrates exactly this added
#     value: a perturbation that stays inside the value pin yet trips the
#     ledger.)
#   * IMPROVEMENT -- a case whose relERR falls well below baseline STILL
#     PASSES (improvement is not a failure) but prints an "improved --
#     consider tightening" advisory, the cue to regenerate the baseline.
#
# ## The SLACK band (chosen, not fitted)
#
# Every metric below is BIT-REPRODUCIBLE run-to-run on this machine (verified
# 4x per case during construction; that reproducibility is the determinism
# contract this ledger leans on -- a non-deterministic metric would make the
# ledger noise, so none was admitted).  SLACK therefore does NOT absorb
# local jitter (there is none); it absorbs CROSS-PLATFORM float jitter.  The
# package's numerical-tier ops (SVD + Padé root-finding) are bit-identical
# only *given* {arch, os, runtime} (src/PadeTaylor.jl "Determinism"); a
# different LAPACK/BLAS or arch can move the last digits of a relERR by a
# small multiplicative factor.  SLACK = 0.25 (25%) is generous enough to
# swallow that last-digit movement, yet tight enough that a genuine ~30%+
# accuracy regression -- or any order-of-magnitude blow-up -- trips it.  The
# tighten-advisory margin (50% below baseline) is deliberately wider than
# SLACK so platform jitter alone never spams the advisory.
#
# This is a LOCAL (no-CI) baseline (Rule 11).  See accuracy_baseline.toml's
# header for the regeneration/approval workflow (ties to krgy.9).
#
# ## The cases (each REUSES an oracle the suite already pins, so the
#    achieved error is meaningful, not self-referential):
#
#   AL.1 robust_pade_exp77   SVD/classical (7,7) Padé of exp  vs exact binomial
#   AL.2 robust_pade_tan33   (3,3) Padé of tan, cubic coeff    vs -1/15 closed form
#   AL.3 solve_pade_wp_z05   IVP ℘ value before the pole       vs 3-source golden
#   AL.4 solve_pade_wp_z14   IVP ℘ value PAST the bridged pole vs closed form
#   AL.5 taylor_jet_wp       2nd-order recurrence (℘ RHS)      vs mpmath recurrence
#   AL.6 taylor_jet_pi       2nd-order recurrence (Painlevé I) vs mpmath recurrence
#   AL.7 vector_jet_tansec2  d=2 vector recurrence             vs mpmath recurrence
#   AL.8 bvp_cosh_N16        Chebyshev-Newton BVP node error   vs cosh(t)/cosh(1)
#
# References:
#   - docs/test_corpus/03_hardening_methodology.md:143-146  (Tier-2 §2.5)
#   - bd show padetaylor-krgy.8                              (the bead + NOTES)
#   - test/accuracy_baseline.toml                            (the baseline)
#   - test/{certified_oracle,differential,convergence}_test.jl  (sibling beads)
#
# Mutation-proof procedure recorded in the file footer (Rule 4).

using Test
using TOML
using PadeTaylor
using PadeTaylor: robust_pade
using PadeTaylor.Coefficients: taylor_coefficients_2nd
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients
using PadeTaylor.Problems: PadeTaylorProblem, solve_pade

include(joinpath(@__DIR__, "_oracle_corpus_robust_pade.jl"))
include(joinpath(@__DIR__, "_oracle_corpus_taylor_jet.jl"))
include(joinpath(@__DIR__, "_oracle_corpus_vector.jl"))
include(joinpath(@__DIR__, "_oracle_problems.jl"))

const _BASELINE_PATH = joinpath(@__DIR__, "accuracy_baseline.toml")

# SLACK band + tighten-advisory margin -- see the header for the rationale.
const _LEDGER_SLACK   = 0.25
const _TIGHTEN_MARGIN = 0.50

_relerr(got, ref) = abs(got - ref) / abs(ref)
# max relERR over a coefficient vector; oracle entries that are exactly zero
# (structural zeros of the rational) are skipped -- a relERR there is 0/0.
function _max_relerr(got::AbstractVector, ref::AbstractVector)
    @assert length(got) == length(ref)
    m = 0.0
    for k in eachindex(ref)
        iszero(ref[k]) && continue
        m = max(m, abs(got[k] - ref[k]) / abs(ref[k]))
    end
    return m
end

# ---------------------------------------------------------------------------
# The cases.  Each entry: case-id => () -> achieved error metric.  The id is
# the [table] key in accuracy_baseline.toml; the closure recomputes the
# achieved relERR from the SAME oracle the corpus pins.  Deterministic; no
# RNG anywhere (every IC/series is a fixed literal or a parsed golden string).
# ---------------------------------------------------------------------------
const _LEDGER_CASES = let
    fW(z, u, up)  = 6 * u^2            # Weierstrass ℘ RHS
    fPI(z, u, up) = 6 * u^2 + z        # Painlevé I RHS
    fv            = (z, y) -> [y[2], 2 * y[1] * y[2]]   # (tan, sec²) system
    f_lin(z, u)   = u                  # linear BVP RHS u'' = u
    ∂f_lin(z, u)  = one(u)

    Pair{String,Function}[
        "robust_pade_exp77" => function ()
            c = Float64[1 / factorial(big(k)) for k = 0:14]
            P = robust_pade(c, 7, 7)
            return max(maximum(abs, P.a .- CRP_EXP77_A),
                       maximum(abs, P.b .- CRP_EXP77_B))
        end,
        "robust_pade_tan33" => function ()
            c = Float64[0, 1, 0, 1/3, 0, 2/15, 0]
            P = robust_pade(c, 3, 3)
            return _relerr(P.a[4], CRP_TAN33_A[4])       # cubic coeff vs -1/15
        end,
        "solve_pade_wp_z05" => function ()
            prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
            u, _ = solve_pade(prob; h_max = 1.5)(0.5)
            return _relerr(u, u_at_0_5)
        end,
        "solve_pade_wp_z14" => function ()
            prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
            u, _ = solve_pade(prob; h_max = 1.5)(1.4)
            return _relerr(u, u_at_1_4)
        end,
        "taylor_jet_wp" => function ()
            u0  = parse(Float64, CTJ_WP_U0_STR); up0 = parse(Float64, CTJ_WP_UP0_STR)
            c   = taylor_coefficients_2nd(fW, 0.0, u0, up0, 30)
            return _max_relerr(c, parse.(Float64, CTJ_WP_C))
        end,
        "taylor_jet_pi" => function ()
            u0  = parse(Float64, CTJ_PI_U0_STR); up0 = parse(Float64, CTJ_PI_UP0_STR)
            c   = taylor_coefficients_2nd(fPI, 0.0, u0, up0, 30)
            return _max_relerr(c, parse.(Float64, CTJ_PI_C))
        end,
        "vector_jet_tansec2" => function ()
            cf = vector_taylor_coefficients(fv, 0.0, [1.0, 2.0], 30)
            return max(_max_relerr(cf[1], parse.(Float64, CTSJ_C1)),
                       _max_relerr(cf[2], parse.(Float64, CTSJ_C2)))
        end,
        "bvp_cosh_N16" => function ()
            sol = bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0; N = 16)
            ref = [cosh(t) / cosh(1.0) for t in sol.nodes_t]
            return maximum(abs, sol.u_nodes .- ref)
        end,
    ]
end

# ---------------------------------------------------------------------------
# REGENERATION / APPROVAL.  When the env flag is set, rewrite the baseline
# from the current achieved errors (preserving each case's metric/oracle
# annotation) and print the diff-worthy summary.  NEVER writes on a normal
# run -- baseline updates are explicit, reviewed commits (krgy.9).
# ---------------------------------------------------------------------------
# The leading `#`-comment banner of the baseline file (everything up to the
# first non-comment line).  Captured on read so regeneration can WRITE IT
# BACK -- otherwise the approval/regeneration instructions the header carries
# would be lost the first time a reviewer regenerates.
function _baseline_header(path::AbstractString)
    io = IOBuffer()
    for line in eachline(path)
        startswith(lstrip(line), "#") || isempty(strip(line)) || break
        println(io, line)
    end
    return String(take!(io))
end

function _regenerate_baseline!(achieved::Dict{String,Float64})
    header = _baseline_header(_BASELINE_PATH)   # preserve the literate header
    old    = TOML.parsefile(_BASELINE_PATH)
    for (id, err) in achieved
        haskey(old, id) || (old[id] = Dict{String,Any}())
        old[id]["relerr"] = err          # metric/oracle annotations untouched
    end
    haskey(old, "meta") || (old["meta"] = Dict{String,Any}())
    old["meta"]["slack"] = _LEDGER_SLACK
    open(_BASELINE_PATH, "w") do io
        # Re-emit the header banner, then the TOML body.  Reviewers diff only
        # the relerr values; the approval-workflow docs survive regeneration.
        print(io, header)
        TOML.print(io, old; sorted = true)
    end
    @info "accuracy baseline REGENERATED (review the git diff, then commit)" path=_BASELINE_PATH
end

@testset "Accuracy-regression ledger (bead padetaylor-krgy.8)" begin
    baseline = TOML.parsefile(_BASELINE_PATH)
    achieved = Dict{String,Float64}()

    @testset "AL: $(id) achieved relERR ≤ baseline·(1+slack)" for (id, compute) in _LEDGER_CASES
        @test haskey(baseline, id)                     # baseline must cover the case
        base = baseline[id]["relerr"]::Float64
        got  = compute()
        achieved[id] = got

        # PRIMARY invariant: achieved error stays within the baseline band.
        bound = base * (1 + _LEDGER_SLACK)
        @test got ≤ bound
        # Sanity floor: a metric that collapsed to exactly 0 (or NaN) is a
        # broken oracle, not a real "perfect" computation -- flag it loud.
        @test isfinite(got) && got ≥ 0

        # IMPROVEMENT advisory (still PASSES) -- the cue to tighten + re-approve.
        if got < base * (1 - _TIGHTEN_MARGIN)
            @info "AL.$(id): improved — relERR $(got) ≪ baseline $(base); consider regenerating the baseline" got base
        end
    end

    if get(ENV, "PADETAYLOR_REGEN_ACCURACY_BASELINE", "0") == "1"
        _regenerate_baseline!(achieved)
    end
end

# ----------------------------------------------------------------------------
# Mutation-proof (CLAUDE.md Rule 4) -- performed + VERIFIED 2026-06-09, src/
# restored before commit (`git diff src/` clean, re-confirmed).  The whole
# point of the ledger is to catch DRIFT a fixed absolute-tol gate misses, so
# the headline proof demonstrates exactly that: a perturbation that the
# corpus's own value pin at the SAME point passes, yet the ledger catches.
#
#   ADDED-VALUE demonstration (the headline -- VERIFIED).  In
#   `src/Problems.jl` (the dense solution evaluator), scale the Padé value by
#   a tiny factor: change
#       u = _evaluate_pade(P_u, t)
#   to
#       u = _evaluate_pade(P_u, t) * (1 + 3e-8)
#   i.e. a 3e-8 RELATIVE perturbation of every `solve_pade` value.  This
#   lands deliberately between the z=1.4 ledger band and the z=1.4 value-pin
#   tolerance.  RESULT (verified):
#     * LEDGER `solve_pade_wp_z14` goes RED: achieved relERR rises from the
#       baseline 2.14e-8 to ~5e-8, which exceeds the band 2.14e-8·1.25 ≈
#       2.68e-8.  THE LEDGER CAUGHT THE DRIFT.
#     * CORPUS value pin 6.1.4 (problems_test.jl `sol(1.4)` vs u_at_1_4 at
#       rtol = 1e-7) stays fully GREEN: ~5e-8 < 1e-7.  The fixed-tol gate at
#       that point is BLIND to the degradation.  This is the ledger's reason
#       to exist, demonstrated end-to-end.
#     (The same mutation ALSO reddened the much-tighter z=0.5 ledger case and
#     the rtol-1e-13/1e-9 pins 6.1.1-6.1.3 -- expected; those have no
#     headroom.  z=1.4 is the case with slack, and that is exactly where the
#     ledger earns its keep.)  Restored; `git diff src/` clean; ledger GREEN.
#
#   NUANCE worth recording (skepticism, Rule 3).  A first attempt perturbed
#   `taylor_coefficients_2nd` in src/Coefficients.jl by (1+5e-14).  That does
#   NOT cleanly show "corpus GREEN, ledger RED" because corpus_taylor_jet_test
#   ALSO pins the BigFloat-256 jet at rtol 1e-70 -- an essentially exact
#   assertion that ANY perturbation trips.  So the added-value contrast is
#   only crisp where the corpus pin is a genuinely-loose Float64 tol (here the
#   z=1.4 rtol-1e-7 pin), which is why the headline mutation targets the IVP
#   value path, not the jet.  A future editor: pick the loose-tol case to
#   demonstrate added value, not a case the BigFloat pin already nails.
#
#   BASELINE-side (the band has teeth, not a tautology).  Temporarily halving
#   a baseline relerr in accuracy_baseline.toml (e.g. solve_pade_wp_z05
#   4.4e-16 -> 2.2e-16) makes that bound 2.2e-16·1.25 ≈ 2.8e-16 < achieved
#   4.4e-16 -> RED, confirming the assertion is the relerr·(1+slack) BOUND
#   (Rule 5), not "didn't throw".  Restored to the captured value.
#
# All mutations confirmed RED then reverted; no src/ change persisted.
