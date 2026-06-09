#!/usr/bin/env julia
# =============================================================================
# test/mutation/run_mutation_gate.jl — the automated mutation-testing gate.
#
# Bead padetaylor-krgy.4 (epic padetaylor-krgy, Tier-1 §1.3 of
# docs/test_corpus/03_hardening_methodology.md).
#
# WHY THIS EXISTS.  The corpus mutation-proofs its numerical-core hotspots BY
# HAND: each load-bearing `_test.jl` carries a narrative "MUTATION-PROOF" footer
# claiming "I perturbed line X, the suite went RED, I restored it."  That is
# excellent discipline (CLAUDE.md Rule 4) but it is *unmeasured and
# unrepeatable* — a refactor can silently weaken a test and the stale footer
# would never notice.  This harness turns those hand-claims into a MEASURED,
# RE-RUNNABLE mutation SCORE: for each catalogued mutant it perturbs `src/`,
# runs the ONE targeted test file that footer points at, and classifies
# KILLED (test went RED ⇒ the test genuinely asserts the invariant) vs
# SURVIVED (test stayed GREEN ⇒ a test-strength GAP → file a bead).
# Mutation score = killed / total (confirmed framing, methodology §1.3).
#
# WHY A BESPOKE HARNESS AND NOT Gremlins.jl.  Gremlins.jl (the maintained
# successor to the dead Vimes.jl; registered, julia 1.10–1.x compat) was
# assessed (external/Gremlins.jl).  It is operator-table-driven over ALL of
# `src/` and applies ONE `test_file` (default `runtests.jl`, the ~22-min full
# suite) to EVERY mutant.  Two structural mismatches make it the wrong tool
# here:  (1) it cannot map a mutant to a *per-hotspot* targeted test file — our
# four hotspots each need a DIFFERENT targeted `_test.jl`;  (2) it generates
# hundreds of mechanical mutants rather than re-checking the CURATED, already-
# hand-validated catalogue below.  Running its default model would cost ~22 min
# PER mutant (the full suite) — infeasible.  This harness instead pairs each
# curated mutant with its targeted test file: a single targeted file runs in
# seconds, and although mutating `src/` invalidates the precompile cache (one
# lazy in-process package recompile per mutant, ~10 s on this box), a tight
# ~14-mutant catalogue is a ~3–5-min PERIODIC gate — NOT a per-commit gate, and
# deliberately NOT included by `test/runtests.jl`.
#
# SAFETY (load-bearing).  `src/` must NEVER be left mutated.  For every mutant
# we read the original bytes UP FRONT, apply the perturbation, and restore in a
# `try/finally` so the original bytes are written back on EVERY path including
# an exception or Ctrl-C.  After each restore we assert via `git diff --quiet`
# that `src/<file>` is byte-clean; a non-clean tree is a hard error that aborts
# the run.  The net is self-tested at startup (`_selftest_restore`): we simulate
# a mid-mutant crash and confirm the file is restored.
#
# USAGE.  `julia --project=. test/mutation/run_mutation_gate.jl`         (full)
#         `julia --project=. test/mutation/run_mutation_gate.jl 1 2 9`   (subset
#           — 1-based catalogue indices, for fast verification).
# See test/mutation/README.md for the runbook (how to read the score, what a
# survivor means).
# =============================================================================

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const SRC  = joinpath(REPO, "src")
const JL   = Base.julia_cmd()

# ── Mutant catalogue ─────────────────────────────────────────────────────────
# Each entry: a SINGLE, unambiguous perturbation of a numerical-core hotspot,
# paired with the ONE test file whose MUTATION-PROOF footer claims to catch it.
# `old` must occur EXACTLY ONCE in `file` (asserted before applying).  The
# `expect` field records what the hand footer (or our analysis) claims:
# `:killed` for a load-bearing mutation a test must catch.
struct Mutant
    hotspot :: String   # human label for per-hotspot scoring
    file    :: String   # src/ file (relative to src/)
    test    :: String   # targeted test/ file (relative to test/)
    old     :: String
    new     :: String
    desc    :: String
    expect  :: Symbol
end

const CATALOGUE = Mutant[
    # ── Coefficients.jl  →  coefficients_test.jl ─────────────────────────────
    Mutant("Coefficients", "Coefficients.jl", "coefficients_test.jl",
        "y[k] = f_t[k-1] / k", "y[k] = f_t[k-1] / (k + 1)",
        "1st-order integration off-by-one (footer Mutation A)", :killed),
    Mutant("Coefficients", "Coefficients.jl", "coefficients_test.jl",
        "u[j] = f_t[j-2] / (j * (j-1))", "u[j] = f_t[j-2] / (j * (j+1))",
        "2nd-order recurrence off-by-one (footer Mutation C)", :killed),
    Mutant("Coefficients", "Coefficients.jl", "coefficients_test.jl",
        "up[j-1] = j * u[j]", "up[j-1] = j * u[j-1]",
        "2nd-order up-resync index slip — KNOWN GAP: coefficients_test asserts " *
        "u not the up jet (bead padetaylor-98pe)", :survived),

    # ── RobustPade.jl  →  robustpade_test.jl ─────────────────────────────────
    Mutant("RobustPade", "RobustPade.jl", "robustpade_test.jl",
        "rhs[i] = -cv[m + i + 1]", "rhs[i] = cv[m + i + 1]",
        "classical Toeplitz RHS sign flip (FW eq. 5.4) (new)", :killed),
    Mutant("RobustPade", "RobustPade.jl", "robustpade_test.jl",
        "s += b_full[j + 1] * cv[k - j + 1]", "s -= b_full[j + 1] * cv[k - j + 1]",
        "classical numerator recurrence sign (FW eq. 5.5) (new)", :killed),
    Mutant("RobustPade", "RobustPade.jl", "robustpade_test.jl",
        "ρ = count(s -> s > ts_typed, S)", "ρ = count(s -> s >= ts_typed, S)",
        "SVD rank-count threshold <→<= — LIKELY EQUIVALENT MUTANT: differs only " *
        "at an exact-threshold σ (bead padetaylor-98pe)", :survived),

    # ── SharedPadeDefect.jl  →  shared_pade_dispatch_test.jl ─────────────────
    Mutant("SharedPadeDefect", "SharedPadeDefect.jl", "shared_pade_dispatch_test.jl",
        "(k - 1) * c[k] for k = 2:length(c)", "(k - 1) * c[k+1] for k = 2:length(c)-1",
        "polyder coefficient index slip (new)", :killed),
    Mutant("SharedPadeDefect", "SharedPadeDefect.jl", "shared_pade_dispatch_test.jl",
        "_horner(numder[i], t) * qv - _horner(numerators[i], t) * qpv",
        "_horner(numder[i], t) * qv + _horner(numerators[i], t) * qpv",
        "quotient-rule derivative sign (new)", :killed),

    # ── SheetTracker.jl (winding)  →  corpus_winding_test.jl ─────────────────
    Mutant("SheetTracker", "SheetTracker.jl", "corpus_winding_test.jl",
        "        Δθ += 2 * π", "        Δθ += 0 * π",
        "winding_delta: drop the (-π,π] wrap (footer Mutation W1)", :killed),
    Mutant("SheetTracker", "SheetTracker.jl", "corpus_winding_test.jl",
        "Δθ = angle(z_new - branch) - angle(z_old - branch)",
        "Δθ = angle(z_old - branch) - angle(z_new - branch)",
        "winding_delta: swap old/new (sign flip) (new)", :killed),
    Mutant("SheetTracker", "SheetTracker.jl", "corpus_winding_test.jl",
        "round(Int, total_winding / (2 * π))", "round(Int, total_winding / (4 * π))",
        "sheet_index: wrong 2π normalisation (new)", :killed),
    Mutant("SheetTracker", "SheetTracker.jl", "corpus_winding_test.jl",
        "out[i] = out[i-1] + winding_delta(path[i-1], path[i], branch)",
        "out[i] = out[i-1] - winding_delta(path[i-1], path[i], branch)",
        "accumulate_winding: sign of increment (new)", :killed),
]

# ── Apply / restore primitives ───────────────────────────────────────────────
"""Apply `m`, return the original file bytes.  Asserts `old` is unique."""
function _apply(m::Mutant)
    path = joinpath(SRC, m.file)
    orig = read(path, String)
    n = count(_ -> true, eachmatch(Regex(_esc(m.old)), orig))
    n == 1 || error("catalogue: `old` occurs $n× (need 1) in $(m.file): $(m.desc)")
    write(path, replace(orig, m.old => m.new))
    return orig, path
end

_esc(s) = replace(s, r"([.*+?^${}()|\[\]\\])" => s"\\\1")

"""Assert `git diff` reports `relfile` byte-clean (restore succeeded)."""
function _assert_clean(relfile::String)
    rel = joinpath("src", relfile)
    ok = success(setenv(`git diff --quiet -- $rel`; dir = REPO))
    ok || error("SAFETY VIOLATION: src/$relfile not byte-clean after restore — " *
                "manual `git checkout src/$relfile` required.")
    return nothing
end

# ── One mutant: apply, run targeted test, classify, ALWAYS restore ───────────
function _run_one(m::Mutant)
    orig, path = _apply(m)
    try
        testpath = joinpath(REPO, "test", m.test)
        cmd = setenv(`$JL --project=$REPO --startup-file=no $testpath`; dir = REPO)
        out = IOBuffer()
        # `success` runs the pipeline and returns a Bool WITHOUT throwing on a
        # non-zero exit — exactly what we want, since a KILLED mutant *is* a
        # non-zero (RED) test exit.  A green test ⇒ `true` ⇒ the mutant survived.
        green = success(pipeline(cmd; stdout = out, stderr = out))
        return green ? :survived : :killed
    finally
        write(path, orig)            # restore on EVERY path (incl. exception)
        _assert_clean(m.file)        # ...and prove it
    end
end

# ── Self-test the restore net (simulate a mid-mutant crash) ──────────────────
function _selftest_restore()
    m = CATALOGUE[1]
    path = joinpath(SRC, m.file)
    before = read(path, String)
    try
        orig, _ = _apply(m)
        try
            read(path, String) == before && error("mutation did not change the file")
            throw(ErrorException("simulated mid-mutant crash"))
        finally
            write(path, orig)
            _assert_clean(m.file)
        end
    catch e
        e isa ErrorException && occursin("simulated", e.msg) || rethrow()
    end
    read(path, String) == before || error("restore self-test FAILED: file not restored")
    @info "restore self-test PASSED (src/$(m.file) byte-clean after simulated crash)"
end

# ── Driver ───────────────────────────────────────────────────────────────────
function main(argv)
    # Pre-flight: working tree must be clean for our hotspot files.
    for m in unique(c -> c.file, CATALOGUE)
        _assert_clean(m.file)
    end
    _selftest_restore()

    idx = isempty(argv) ? eachindex(CATALOGUE) : parse.(Int, argv)
    @info "running $(length(idx)) of $(length(CATALOGUE)) catalogued mutants" subset = !isempty(argv)

    results = Tuple{Mutant,Symbol}[]
    for i in idx
        m = CATALOGUE[i]
        print("[$i/$(length(CATALOGUE))] $(m.hotspot): $(m.desc) … ")
        flush(stdout)
        t = @elapsed (outcome = _run_one(m))
        emoji = outcome == :killed ? "KILLED " : "SURVIVED"
        unexpected = outcome != m.expect ? "  *** UNEXPECTED (footer claimed $(m.expect)) ***" : ""
        println("$emoji  ($(round(t; digits = 1))s)$unexpected")
        push!(results, (m, outcome))
    end

    _report(results)
    # Exit non-zero only on an UNEXPECTED outcome (result contradicts the
    # catalogued `expect`): a previously-:killed mutant that now SURVIVES (a
    # weakened/regressed test), OR a tracked :survived mutant that now gets
    # KILLED (a known gap was closed ⇒ flip its catalogue annotation to :killed).
    # Tracked :survived gaps (e.g. bead padetaylor-98pe) are NOT failures.
    bad = count(((m, o),) -> o != m.expect, results)
    return bad == 0 ? 0 : 1
end

function _report(results)
    println("\n", "="^70)
    println("MUTATION GATE REPORT")
    println("="^70)
    total  = length(results)
    killed = count(((_, o),) -> o == :killed, results)
    println("OVERALL  score = $killed / $total killed " *
            "($(round(100 * killed / total; digits = 1))%)")
    println("\nPER HOTSPOT:")
    for hs in unique(((m, _),) -> m.hotspot, results) .|> (r -> r[1].hotspot)
        sub = filter(((m, _),) -> m.hotspot == hs, results)
        k = count(((_, o),) -> o == :killed, sub)
        println("  $hs: $k / $(length(sub)) killed")
    end
    survivors = filter(((_, o),) -> o == :survived, results)
    if isempty(survivors)
        println("\nSURVIVORS: none — every catalogued mutant met its expectation.")
    else
        println("\nSURVIVORS:")
        for (m, _) in survivors
            ln  = _lineno(m)
            tag = m.expect == :survived ? "KNOWN/tracked" : "NEW GAP → file a bead"
            println("  src/$(m.file):$ln  [$tag]  — $(m.desc)")
            println("      targeted test: test/$(m.test)")
        end
    end
    unexpected = count(((m, o),) -> o != m.expect, results)
    unexpected == 0 || println("\n*** $unexpected UNEXPECTED outcome(s) — GATE FAILS ***")
    println("="^70)
end

"""1-based line of `m.old` in `src/<file>` (for survivor `file:line` report)."""
function _lineno(m::Mutant)
    src = read(joinpath(SRC, m.file), String)
    pos = findfirst(m.old, src)
    pos === nothing && return 0
    return count(==('\n'), src[1:first(pos)]) + 1
end

exit(main(ARGS))
