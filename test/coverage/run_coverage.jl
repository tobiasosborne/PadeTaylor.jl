#!/usr/bin/env julia
# =============================================================================
# test/coverage/run_coverage.jl — the local code-coverage measurement tool.
#
# Bead padetaylor-krgy.7 (epic padetaylor-krgy, Tier-1 §1.4 of
# docs/test_corpus/03_hardening_methodology.md).
#
# WHY THIS EXISTS.  Coverage tells you which `src/` lines NO test exercises —
# a *necessary, not sufficient* floor.  An uncovered numerical branch is a
# guaranteed blind spot (no test ran it, so no test can have caught a
# regression in it).  The converse does NOT hold: a COVERED line is NOT a
# TESTED line — a test may execute a line incidentally and assert nothing
# about it.  That second axis (test STRENGTH) is what the mutation gate
# (`test/mutation/run_mutation_gate.jl`) measures.  Read the two together:
# coverage finds lines no test touches; mutation finds lines tests touch but
# do not pin.  Do NOT chase 100% coverage — defensive fail-loud `throw`/`error`
# guards (CLAUDE.md Rule 1), `verbose`-print blocks, display `show` methods,
# and arb-precision-only paths are legitimately cold at Float64 and are
# triaged as intentionally-uncovered in `docs/test_corpus/coverage_snapshot.md`.
#
# WHY A STANDALONE TOOL, NOT PART OF `Pkg.test`.  Mirrors the mutation gate:
# this is a PERIODIC measurement, not a per-commit gate, and is deliberately
# NOT included by `test/runtests.jl`.  Coverage.jl is a TOOL dependency, not a
# SUITE dependency — it is installed into a throwaway temp env here (never
# added to the package `[deps]`/`[extras]`/`[targets]`), exactly as the
# mutation gate avoids adding a suite dep.  The package itself stays lean.
#
# TWO WAYS TO RUN (see test/coverage/README.md for the runbook):
#
#   1. FULL (generate + aggregate) — runs the suite under `--code-coverage`,
#      which is the canonical ~22-min path:
#        julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
#      then aggregate the `src/*.cov` it drops:
#        julia test/coverage/run_coverage.jl
#
#   2. AGGREGATE-ONLY (fast) — if `src/*.cov` already exist on disk (a recent
#      `Pkg.test(coverage=true)` left them), just summarise them:
#        julia test/coverage/run_coverage.jl
#      This is the same command — the tool NEVER regenerates; it only reads the
#      `.cov` already present.  Pass `--run` to have the tool drive the full
#      `Pkg.test(coverage=true)` itself before aggregating (slow; one Julia
#      process — Rule 7).
#
# SAFETY.  `.cov` files are gitignored (`.gitignore:65 *.cov`).  This tool does
# NOT delete them (so a caller can inspect per-line `.cov` after a run); the
# snapshot workflow deletes them once the numbers are captured.  The tool never
# touches `src/`.
# =============================================================================

import Pkg

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const SRC  = joinpath(REPO, "src")

# ── Install Coverage.jl into a throwaway env (TOOL dep, never a suite dep) ────
"""Activate a temp env, add Coverage, and `import` it.  Leaves the active
project pointing at the temp dir for the rest of this process (we never need
the package env here — we only read `.cov` byte files)."""
function _load_coverage()
    dir = mktempdir()
    Pkg.activate(dir; io = devnull)
    Pkg.add("Coverage"; io = devnull)
    @eval import Coverage
    return nothing
end

# ── Optionally drive the full suite under --code-coverage (the slow path) ─────
"""Run `Pkg.test(coverage=true)` in the PACKAGE env as a child process (so the
temp env we use for Coverage.jl does not shadow the package).  One Julia child,
serial — Rule 7."""
function _generate()
    @info "run_coverage: regenerating coverage — Pkg.test(coverage=true) (~22 min)"
    jl  = Base.julia_cmd()
    cmd = setenv(`$jl --project=$REPO -e "using Pkg; Pkg.test(coverage=true)"`; dir = REPO)
    run(cmd)
    return nothing
end

# ── Aggregate the per-file `.cov` already on disk ─────────────────────────────
struct FileCov
    name    :: String   # basename
    covered :: Int      # executable lines hit ≥ once
    total   :: Int      # executable lines
    uncov   :: Vector{Int}   # 1-based line numbers with hit-count == 0
end

pct(f::FileCov) = f.total == 0 ? 100.0 : 100 * f.covered / f.total

"""Process `src/*.cov` via Coverage.jl → sorted-ascending `Vector{FileCov}`
plus the overall (covered, total).  Errors loudly when NO `.cov` hit-data is
present (Rule 1).  Note: `Coverage.process_folder` does NOT signal a missing
`.cov` by returning an empty vector — it returns one entry per `src/*.jl` with
all-`nothing` coverage ("Assuming file has no coverage"), which would otherwise
masquerade as an honest-looking 0/N.  We catch that by asserting that SOME line
was actually hit (`tc > 0`): a real `Pkg.test(coverage=true)` run hits
thousands of lines, so `tc == 0` means the `.cov` are simply absent."""
function aggregate()
    # `Base.invokelatest` is load-bearing: Coverage.jl is imported at RUNTIME
    # (`@eval import Coverage` in `_load_coverage`), so its methods are newer
    # than this enclosing world age — a direct call hits a world-age barrier.
    covs = Base.invokelatest(Coverage.process_folder, SRC)
    rows = FileCov[]
    tc = 0; tt = 0
    for fc in covs
        covered = count(x -> x !== nothing && x > 0, fc.coverage)
        total   = count(x -> x !== nothing, fc.coverage)
        uncov   = [i for (i, v) in enumerate(fc.coverage) if v !== nothing && v == 0]
        push!(rows, FileCov(basename(fc.filename), covered, total, uncov))
        tc += covered; tt += total
    end
    tc == 0 && error(
        "run_coverage: no `src/*.cov` hit-data found under $SRC (0 lines covered).  " *
        "Suggestion: generate the `.cov` first with " *
        "`julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'`, " *
        "or pass `--run` to have this tool drive it (~22 min).")
    sort!(rows; by = pct)
    return rows, tc, tt
end

# ── Summarise + flag the worst-covered files ──────────────────────────────────
function _report(rows::Vector{FileCov}, tc::Int, tt::Int; worst::Int = 12)
    println("\n", "="^72)
    println("COVERAGE REPORT  (necessary-not-sufficient; pair with the mutation gate)")
    println("="^72)
    println("OVERALL  $tc / $tt executable lines  ",
            "($(round(100 * tc / tt; digits = 2))%)")

    below = filter(r -> pct(r) < 100.0, rows)
    println("\nFILES BELOW 100% (worst first):")
    if isempty(below)
        println("  none — every file fully covered (coverage floor met; ",
                "now lean on the mutation gate for STRENGTH).")
    else
        for r in below
            println("  ", rpad(r.name, 32), lpad(string(r.covered), 5), "/",
                    lpad(string(r.total), 5), "  ",
                    rpad(string(round(pct(r); digits = 1)), 6), "%")
        end
    end

    println("\nWORST $(min(worst, length(below))) — UNCOVERED LINE NUMBERS (triage these):")
    for r in first(below, worst)
        runs = _compress(r.uncov)
        println("  src/", r.name, "  →  ", join(runs, ", "))
    end
    println("="^72)
    println("Triage verdict per region lives in docs/test_corpus/coverage_snapshot.md.")
    return nothing
end

"""Compress a sorted line-number vector into `a-b` runs for compact reporting."""
function _compress(v::Vector{Int})
    isempty(v) && return String[]
    runs = String[]
    s = v[1]; p = v[1]
    for x in v[2:end]
        if x == p + 1
            p = x
        else
            push!(runs, s == p ? string(s) : "$s-$p"); s = x; p = x
        end
    end
    push!(runs, s == p ? string(s) : "$s-$p")
    return runs
end

# ── Driver ────────────────────────────────────────────────────────────────────
function main(argv)
    if "--run" in argv
        _generate()                          # slow: regenerate .cov first
    end
    _load_coverage()                         # tool dep into throwaway env
    rows, tc, tt = aggregate()
    _report(rows, tc, tt)
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
