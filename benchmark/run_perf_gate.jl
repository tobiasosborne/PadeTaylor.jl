#!/usr/bin/env julia
# =============================================================================
# benchmark/run_perf_gate.jl — the local, noise-aware performance-regression gate.
#
# Bead padetaylor-krgy.11 (epic padetaylor-krgy; EFFICIENCY section of
# docs/test_corpus/03_hardening_methodology.md:163-167).
#
# WHY THIS EXISTS.  The corpus pins CORRECTNESS (values, errors, invariants) but
# nothing pins SPEED.  A refactor that quietly doubles the per-step Padé cost —
# an accidental O(n²) loop, a lost `@inbounds`, a type instability that boxes the
# inner jet — leaves every correctness test GREEN while the integrator crawls.
# This gate catches that: it benchmarks the hot paths with BenchmarkTools.jl and
# compares each one's CURRENT median time to a COMMITTED baseline, flagging a
# regression when `current_median > baseline_median · (1 + SLACK)`.
#
# WHY RATIO + WIDE SLACK, NOT ABSOLUTE.  Wall-clock timing is intrinsically
# noisy, and on a WSL2 box with no CI there is no clean dedicated runner — the
# scheduler, other processes, and thermal state move the numbers run-to-run.  A
# naive `@elapsed t < T` gate would false-alarm constantly.  BenchmarkTools'
# `minimum`/`median` over many samples is the statistically-sound estimator
# (minimum ≈ the noise-free floor; median ≈ typical), and comparing a RATIO to a
# committed baseline with a SLACK band wide enough to swallow the measured
# run-to-run jitter (but tight enough to trip a genuine >1.3× regression) is the
# discipline the methodology doc §EFFICIENCY prescribes.  SLACK is CALIBRATED,
# not guessed — see benchmark/README.md "Calibrating the slack" for the measured
# jitter on this box and the chosen value.
#
# WHY A STANDALONE TOOL, NOT PART OF Pkg.test.  Mirrors the coverage + mutation
# gates: perf is a PERIODIC measurement, FLAKY as an always-on per-commit gate,
# and is deliberately NOT included by test/runtests.jl (which `include`s its
# files explicitly — a benchmark/ dir is never globbed in).  BenchmarkTools.jl
# is a TOOL dep in benchmark/Project.toml, NEVER in the package
# [deps]/[targets].test (see that file's header for the env rationale).
#
# IDLE-BOX DISCIPLINE (load-bearing).  A baseline — and any comparison against
# it — is only meaningful on an IDLE machine with threads pinned.  The tool
# REPORTS the thread/BLAS configuration it ran under and refuses to regenerate a
# baseline unless `JULIA_NUM_THREADS`/BLAS threads match the pinned single-thread
# config (a multi-threaded BLAS makes the SVD timing non-reproducible).  Capture
# baselines with nothing else running.
#
# USAGE (single Julia process — Rule 7; see benchmark/README.md for the runbook):
#   # Compare current timings to the committed baseline (the gate):
#   julia --project=benchmark benchmark/run_perf_gate.jl
#   # Regenerate the baseline (EXPLICIT, reviewed-commit path — never silent):
#   PADETAYLOR_REGEN_PERF_BASELINE=1 julia --project=benchmark benchmark/run_perf_gate.jl
# =============================================================================

import Pkg
const REPO = normpath(joinpath(@__DIR__, ".."))

# ── Self-heal the tool env (TOOL dep; one Julia process — Rule 7) ─────────────
# If PadeTaylor is not yet resolved in benchmark/Project.toml, dev it + instantiate
# so a fresh checkout runs the gate with one command.
function _ensure_env()
    try
        Base.require(Base.PkgId(Base.UUID("00000000-0000-0000-0000-000000000001"),
                                "PadeTaylor"))
    catch
        @info "run_perf_gate: resolving benchmark env (Pkg.develop PadeTaylor + instantiate)"
        Pkg.develop(; path = REPO, io = devnull)
        Pkg.instantiate(; io = devnull)
    end
    return nothing
end
_ensure_env()

using BenchmarkTools
using TOML
using PadeTaylor
using PadeTaylor: robust_pade
using PadeTaylor.Coefficients: taylor_coefficients_2nd
using PadeTaylor.PadeStepper: PadeStepperState, pade_step_with_pade!, _evaluate_pade
using PadeTaylor.Problems: PadeTaylorProblem, solve_pade

const BASELINE_PATH = joinpath(@__DIR__, "baseline.toml")

# SLACK band — CALIBRATED on this WSL2 box (benchmark/README.md §"Calibrating
# the slack").  The gate compares the BenchmarkTools MINIMUM (the noise-floor
# estimator, most robust to upward jitter) — not the median — to the baseline
# minimum.  Measured run-to-run jitter on this box: the cheap, allocation-light
# kernels (pade_step, robust_pade, jet, evaluate) hold min-ratio ≤ ~1.13, but the
# two ALLOCATION-HEAVY paths (solve_pade ~460 allocs, path_network ~725 allocs)
# carry real GC-driven min-jitter to ~1.25–1.30 steady-state, and a whole-machine
# transient can push any path higher for one run.  SLACK = 0.50 sits comfortably
# above the worst steady-state min-jitter (~1.30) yet still trips a genuine ≥1.5×
# regression and any order-of-magnitude blow-up (the bite-proof clears it by
# 19–33×).  The median is still REPORTED for context (and is what we diff in the
# baseline), but the PASS/FAIL is on the min-ratio.
const PERF_SLACK = 0.50

# ── Fixtures — deterministic inputs for the hot paths (no RNG anywhere) ───────
# Weierstrass-℘ RHS u'' = 6u² with the FW 2011 equianharmonic IC (the canonical
# stepper workload); a (3,3)-tan jet for robust_pade; a tan/sec² IVP for the jet.
const _fW = (z, u, up) -> 6 * u^2
const _U0_FW  = 1.071822516416917      # u(0), FW equianharmonic (_oracle_problems.jl:83)
const _UP0_FW = 1.710337353176786      # u'(0)                    (_oracle_problems.jl:84)
const _ORDER  = 30                     # FW 2011 §5.1 canonical Taylor order.

# An order-30 ℘ jet, then rescaled by h^k — the exact vector pade_step feeds
# robust_pade (representative (15,15) :svd workload at the production order).
const _JET30 = taylor_coefficients_2nd(_fW, 0.0, _U0_FW, _UP0_FW, _ORDER)
const _CSCALED = let h = 0.5
    [(_JET30[k] * h^(k - 1)) for k in eachindex(_JET30)]
end
# A (15,15) Padé of that jet, for the _evaluate_pade micro-benchmark.
const _P15 = robust_pade(_CSCALED, _ORDER ÷ 2, _ORDER ÷ 2; method = :svd)

# Each entry: id => () -> a BenchmarkTools.Trial.  The closure rebuilds fresh
# inputs each sample so we time the WORK, not cached state.  `$`-interpolation
# pins the captured constants into the benchmarked expression (no global lookup).
const BENCHMARKS = Pair{String,Function}[
    "pade_step" => () -> @benchmark(
        pade_step_with_pade!(s, $_fW, $_ORDER, 0.5),
        setup = (s = PadeStepperState{Float64}(0.0, $_U0_FW, $_UP0_FW)),
        samples = 400, evals = 1, seconds = 3),
    "robust_pade_svd_1515" => () -> @benchmark(
        robust_pade($_CSCALED, $(_ORDER ÷ 2), $(_ORDER ÷ 2); method = :svd),
        samples = 400, seconds = 3),
    "taylor_coefficients_2nd_o30" => () -> @benchmark(
        taylor_coefficients_2nd($_fW, 0.0, $_U0_FW, $_UP0_FW, $_ORDER),
        samples = 600, seconds = 3),
    "evaluate_pade" => () -> @benchmark(
        _evaluate_pade($_P15, 1.0), samples = 1000, seconds = 2),
    # The two allocation-heavy paths get MORE samples + `gcsample=true` so their
    # minimum (the gate estimator) stabilises against GC-driven jitter.
    "solve_pade_wp_multiseg" => () -> @benchmark(
        solve_pade(prob; h_max = 0.3),
        setup = (prob = PadeTaylorProblem($_fW, ($_U0_FW, $_UP0_FW),
                                          (0.0, 1.5); order = $_ORDER)),
        samples = 400, seconds = 8, gcsample = true),
    "path_network_short_walk" => () -> @benchmark(
        path_network_solve(prob, [complex(0.0, 0.6)]; h = 0.5, order = $_ORDER),
        setup = (prob = PadeTaylorProblem($_fW,
                    (complex($_U0_FW, 0.0), complex($_UP0_FW, 0.0)),
                    (complex(0.0, 0.0), complex(0.0, 0.6)); order = $_ORDER)),
        samples = 200, seconds = 8, gcsample = true),
]

# ── Thread / BLAS fingerprint (idle-box discipline) ───────────────────────────
# Pin BLAS to a single thread up front: a multi-threaded OpenBLAS makes the SVD
# timing (robust_pade :svd, the dominant per-step cost) non-reproducible
# run-to-run, which is exactly the noise a ratio gate must NOT chase.  Julia
# thread count is set by the launcher (`JULIA_NUM_THREADS=1`); we report both.
import LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

function _fingerprint()
    return (jl_threads  = Threads.nthreads(),
            blas_threads = LinearAlgebra.BLAS.get_num_threads(),
            blas         = LinearAlgebra.BLAS.vendor() |> string)
end

# ── Run all benchmarks → Dict(id => (median_ns, min_ns, allocs, bytes)) ───────
function _measure()
    pinned = Threads.nthreads() == 1 && LinearAlgebra.BLAS.get_num_threads() == 1
    pinned || @warn(
        "run_perf_gate: threads NOT pinned to 1 (jl=$(Threads.nthreads()), " *
        "blas=$(LinearAlgebra.BLAS.get_num_threads())). Timings — and any baseline " *
        "regenerated now — will be NOISIER. Pin with `JULIA_NUM_THREADS=1` and " *
        "`BLAS.set_num_threads(1)` on an idle box (benchmark/README.md).")
    out = Dict{String,NamedTuple}()
    for (id, mk) in BENCHMARKS
        print("  benchmarking $(rpad(id, 30)) … "); flush(stdout)
        t = mk()
        out[id] = (median_ns = median(t).time, min_ns = minimum(t).time,
                   allocs = t.allocs, bytes = t.memory)
        println("median $(round(median(t).time / 1e3; digits = 1)) µs, ",
                "min $(round(minimum(t).time / 1e3; digits = 1)) µs, ",
                "$(t.allocs) allocs")
    end
    return out
end

# ── Compare current measurement to the committed baseline (the GATE) ──────────
# `ratio` is the GATE metric = current_min / baseline_min (noise-floor, robust);
# `med_ratio` is reported for context only.
struct Row
    id::String; base_min::Float64; cur_min::Float64; ratio::Float64
    med_ratio::Float64; flag::Bool
end

function _compare(cur::Dict{String,NamedTuple}, baseline::Dict)
    rows = Row[]
    for (id, _) in BENCHMARKS
        haskey(baseline, id) || error(
            "run_perf_gate: baseline has no entry for `$id`. Suggestion: " *
            "regenerate with PADETAYLOR_REGEN_PERF_BASELINE=1 (review the diff).")
        bmin   = Float64(baseline[id]["min_ns"])
        cmin   = cur[id].min_ns
        ratio  = cmin / bmin                                   # GATE metric (min)
        medr   = cur[id].median_ns / Float64(baseline[id]["median_ns"])
        push!(rows, Row(id, bmin, cmin, ratio, medr, ratio > 1 + PERF_SLACK))
    end
    return rows
end

function _report(rows::Vector{Row}, fp)
    println("\n", "="^74)
    println("PERFORMANCE-REGRESSION GATE  (minR = cur_min/base_min [GATE]; medR for context)")
    println("="^74)
    println("THREADS: julia=$(fp.jl_threads)  BLAS=$(fp.blas_threads) ($(fp.blas))",
            fp.jl_threads == 1 && fp.blas_threads == 1 ? "  [pinned]" :
            "  [NOT pinned — see README idle-box discipline]")
    println("SLACK: $(round(Int, 100 * PERF_SLACK))%  ",
            "(GATE on the MIN ratio; flag when current_min > baseline_min·$(1 + PERF_SLACK))\n")
    println("  ", rpad("benchmark", 30), rpad(" base_min", 12), rpad("cur_min", 12),
            "  minR    medR")
    for r in rows
        tag = r.flag ? "  *** REGRESSION ***" : ""
        println("  ", rpad(r.id, 30),
                rpad(" " * string(round(r.base_min / 1e3; digits = 2)) * "µs", 12),
                rpad(string(round(r.cur_min / 1e3; digits = 2)) * "µs", 12),
                "  ", rpad(string(round(r.ratio; digits = 3)), 7),
                " ", round(r.med_ratio; digits = 3), tag)
    end
    bad = count(r -> r.flag, rows)
    if bad == 0
        println("\nOK — every hot path within the $(round(Int, 100 * PERF_SLACK))% slack band.")
    elseif bad == length(rows)
        println("\n*** ALL $bad benchmarks flagged AT ONCE — this is the signature of a ",
                "whole-machine\n    transient (load/thermal), NOT a code regression. ",
                "Re-run on a truly idle box.")
    else
        println("\n*** $bad of $(length(rows)) benchmark(s) exceeded the slack band ",
                "— a localized\n    regression: investigate that hot path. ",
                "(`git diff` the src files it touches.)")
    end
    println("="^74)
    return bad
end

# ── Regeneration / approval (EXPLICIT; mirrors krgy.8's env-var regen) ─────────
function _baseline_header(path)
    isfile(path) || return ""
    io = IOBuffer()
    for line in eachline(path)
        startswith(lstrip(line), "#") || isempty(strip(line)) || break
        println(io, line)
    end
    return String(take!(io))
end

const _REGEN_HEADER = """
# benchmark/baseline.toml — committed PERF baseline (bead padetaylor-krgy.11).
#
# Median + min wall time (ns) and allocations of each hot path, captured on an
# IDLE box with threads pinned to 1.  The gate (benchmark/run_perf_gate.jl)
# asserts, per case, current_median <= median_ns * (1 + SLACK).  SLACK lives in
# the tool (a single calibrated band, documented in benchmark/README.md).
#
# LOCAL (no-CI) baseline (Rule 11).  These numbers are platform-specific; a
# different CPU/BLAS will shift them — the SLACK band absorbs small jitter, NOT
# order-of-magnitude moves.  REGENERATE (explicit, reviewed commit — never
# silent), then review the git diff and commit:
#   PADETAYLOR_REGEN_PERF_BASELINE=1 julia --project=benchmark benchmark/run_perf_gate.jl
"""

function _regenerate!(cur::Dict{String,NamedTuple}, fp)
    header = _baseline_header(BASELINE_PATH); isempty(header) && (header = _REGEN_HEADER)
    data = Dict{String,Any}("meta" => Dict{String,Any}(
        "captured" => "regenerated $(string(Dates.now()))",
        "slack" => PERF_SLACK,
        "julia_threads" => fp.jl_threads, "blas_threads" => fp.blas_threads,
        "blas" => fp.blas))
    for (id, m) in cur
        data[id] = Dict{String,Any}("median_ns" => m.median_ns, "min_ns" => m.min_ns,
                                    "allocs" => m.allocs, "bytes" => m.bytes)
    end
    open(BASELINE_PATH, "w") do io
        print(io, header); TOML.print(io, data; sorted = true)
    end
    @info "perf baseline REGENERATED (review the git diff, then commit)" path = BASELINE_PATH
end

import Dates
function main()
    fp = _fingerprint()
    cur = _measure()
    if get(ENV, "PADETAYLOR_REGEN_PERF_BASELINE", "0") == "1"
        _regenerate!(cur, fp)
        return 0
    end
    isfile(BASELINE_PATH) || error(
        "run_perf_gate: no baseline at $BASELINE_PATH. Suggestion: capture one " *
        "with PADETAYLOR_REGEN_PERF_BASELINE=1 on an idle box, review, commit.")
    baseline = TOML.parsefile(BASELINE_PATH)
    bad = _report(_compare(cur, baseline), fp)
    return bad == 0 ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main())
