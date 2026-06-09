# test/quality_test.jl -- bead `padetaylor-krgy.2`: the static-analysis +
# lint quality gate (Tier-1, highest-ROI of the test-hardening sweep).
#
# ## Why this chapter exists
#
# The 90-file corpus pins VALUES (FW Table 5.1 numbers, closed-form
# Weierstrass/Painlevé, pre-captured `padeapprox.m` output) and exercises
# fail-loud error paths.  What it cannot assert is a whole class of
# *package-shape* defects that never show up as a wrong number on any one
# input: a method ambiguity that only bites a caller we never wrote; an
# unbound type parameter; a `[compat]` bound that has silently rotted; a
# stale `using` that imports a name nobody calls; a type instability that
# slows the inner step loop without changing its result.  Those are exactly
# what abstract-interpretation + package-hygiene tooling catches, and they
# are the cheapest correctness wins available (methodology doc
# `docs/test_corpus/03_hardening_methodology.md` §1.1, "static analysis +
# linting — krgy.2, P1, highest ROI").
#
# Three independent tools, each catching a disjoint defect class:
#
#   * **Aqua.jl** — package-quality auditor.  `test_all` runs eight checks:
#     method ambiguities, unbound type parameters, undefined exports,
#     Project.toml/test-extras agreement, stale (declared-but-unused) deps,
#     `[compat]` completeness, type piracy, and persistent-task leaks.  None
#     of these is an arithmetic property; all are invisible to value pins.
#
#   * **JET.jl** — abstract interpreter over Julia's own type inference.
#     `report_package` walks every method PadeTaylor defines and flags calls
#     that inference proves cannot dispatch (MethodError), undefined globals,
#     and field-access errors — statically, with zero test inputs.  This is
#     the same engine that bead krgy.12 will reuse for `@report_opt`
#     type-stability gating (shared setup, per the krgy.10 reconciliation).
#
#   * **ExplicitImports.jl** — import hygiene.  `check_no_implicit_imports`
#     proves every name a module uses is explicitly imported (no reliance on
#     `using Foo` dumping its whole namespace, which makes upstream renames
#     silent breakages); `check_no_stale_explicit_imports` proves every
#     explicit import is actually used.
#
# ## Triage record (Law 1 / Rule 9 — every finding decided, none blanket-muted)
#
# First run surfaced four real findings; each was triaged senior-style into
# fix-now / allow-list-with-justification / defer-to-bead:
#
#   1. Aqua `deps_compat` — `Printf`, `Random`, `Test` lacked `[compat]`
#      entries.  FIXED in Project.toml (added `Printf = "1"`, `Random = "1"`,
#      `Test = "1"`; stdlibs, trivially safe).
#   2. Aqua `stale_deps` — `Statistics` is declared in `[deps]` but the core
#      module never uses it; it is consumed ONLY by
#      `ext/PadeTaylorDiagnosticsExt.jl:142` (median/quantile/mean over
#      loop-closure edge stats).  Per ADR-0003's minimal-core-dep-graph
#      principle it belongs in `[weakdeps]` — but moving it there forces it
#      to become a *co-trigger* of the extension, which (verified empirically
#      this session) regresses ADR-0016's documented single-trigger
#      `using PadeTaylor, DelaunayTriangulation` activation.  That is a real
#      but non-trivial fix touching the public extension contract, so it is
#      ALLOW-LISTED here (`stale_deps = (ignore = [:Statistics],)`) and
#      DEFERRED to bead `padetaylor-fmf8`.  This is the only allow-list; the
#      other seven Aqua checks run unmodified.
#   3. ExplicitImports stale — `src/SharedPadeDefect.jl` imported `norm` from
#      LinearAlgebra but only ever calls `eigvals`.  FIXED in src (dropped
#      `norm` from the import list).
#   4. JET `report_package` — one method-error report inside PadeTaylor's own
#      modules: `quality_diagnose(::PathNetworkSolution)` at
#      `src/PathNetwork.jl:711`.  This is NOT a bug: the `quality_diagnose`
#      generic is declared empty in `Diagnostics.jl` and its
#      `PathNetworkSolution` method lives in the DelaunayTriangulation
#      extension (ADR-0016), which is not loaded during static
#      `report_package`.  The call site at `PathNetwork.jl:709-723` already
#      wraps it in a `try/catch` that rethrows the MethodError as a
#      `using DelaunayTriangulation` suggestion (Rule 1 fail-loud).  It is the
#      single documented dynamic spot the bead description anticipates
#      ("allow-list known dynamic spots") and is filtered below.  It is present
#      (count 1) when this file runs STANDALONE, but ABSENT (count 0) under the
#      full suite once an earlier test has loaded the extension (e.g.
#      `using Makie` in ext_makie_test.jl transitively loads
#      DelaunayTriangulation, defining the method) — hence the gate asserts
#      `count ≤ 1`, never `== 1`, to be deterministic across test-load orders.
#
# JET runs Julia's internal inference machinery, which on a fresh major like
# 1.12 can emit spurious reports when it descends into Base/SparseArrays
# (reinterpret/bitcast/CartesianIndices chains).  We suppress that noise the
# supported way — `target_modules = (PadeTaylor,)`, which restricts reporting
# to methods PadeTaylor itself defines — rather than by a blanket ignore
# list.  With that scope the only remaining report is the documented
# extension-dispatch spot above.
#
# Reference: docs/test_corpus/03_hardening_methodology.md §1.1; ADR-0003
# (extension / weakdep pattern); ADR-0016 (Diagnostics extension); Aqua.jl,
# JET.jl, ExplicitImports.jl upstream docs.

using Test
using PadeTaylor
using Aqua
using JET
using ExplicitImports

@testset "Quality — static analysis" begin

    @testset "Aqua.test_all (package hygiene)" begin
        # All eight Aqua checks.  The sole allow-list is `Statistics` under
        # stale_deps (extension-only stdlib dep; bead padetaylor-fmf8) — every
        # other check (ambiguities, unbound_args, undefined_exports,
        # project_extras, deps_compat, piracies, persistent_tasks) runs with
        # Aqua's defaults so a real regression in any of them fails the gate.
        Aqua.test_all(PadeTaylor; stale_deps = (ignore = [:Statistics],))
    end

    @testset "JET.report_package (abstract interpretation)" begin
        # Scope reporting to PadeTaylor's own methods (drops Base/stdlib
        # inference noise that JET surfaces on Julia 1.12); then allow-list
        # exactly the one documented extension-dispatch spot.
        report = JET.report_package(PadeTaylor; target_modules = (PadeTaylor,))
        reports = JET.get_reports(report)

        # `quality_diagnose(::PathNetworkSolution)` resolves only when the
        # DelaunayTriangulation extension is loaded (ADR-0016); under static
        # `report_package` that method is absent by design and the call site
        # (src/PathNetwork.jl:711) catches the MethodError and rethrows a
        # `using DelaunayTriangulation` suggestion.  Match it by signature
        # text and exclude it; ANY other report fails the gate.
        is_known_extension_dispatch(r) =
            occursin("quality_diagnose", sprint(show, r))
        unexpected = filter(!is_known_extension_dispatch, reports)

        # The known spot is present ONLY when the Diagnostics extension is NOT
        # loaded.  Under the FULL suite an earlier test activates it — `using
        # Makie` (ext_makie_test.jl) transitively loads DelaunayTriangulation,
        # the `PadeTaylorDiagnosticsExt` trigger, which DEFINES
        # `quality_diagnose(::PathNetworkSolution)` — so the call dispatches and
        # the report vanishes (count 0).  Run standalone, the method is absent
        # and JET reports it (count 1).  Assert ≤ 1 (never MORE, so the
        # allow-list still cannot mask a brand-new report) to stay deterministic
        # regardless of test-load order; `isempty(unexpected)` below is the
        # load-bearing invariant either way.
        @test count(is_known_extension_dispatch, reports) <= 1
        # And that nothing else inference-broke inside our modules.
        if !isempty(unexpected)
            for r in unexpected
                @info "Unexpected JET report" report = sprint(show, r)
            end
        end
        @test isempty(unexpected)
    end

    @testset "ExplicitImports (import hygiene)" begin
        # `check_no_*` return `nothing` on success and throw a descriptive
        # exception on failure; asserting `=== nothing` makes the pass an
        # invariant rather than a "didn't throw".
        @test check_no_implicit_imports(PadeTaylor) === nothing
        @test check_no_stale_explicit_imports(PadeTaylor) === nothing
    end

end # @testset "Quality — static analysis"

# ── Mutation-proof procedure (verified 2026-06-09, bead padetaylor-krgy.2) ──
#
# A static-analysis gate is only worth its runtime if it actually goes RED on
# the defects it claims to catch.  Each tool was proven to bite by injecting a
# genuine defect of its class, running THIS file, confirming the failure, then
# restoring the source exactly (`git diff src/` clean of the injection after).
#
#   Mutation A (Aqua / ExplicitImports — stale import)
#     In `src/SharedPadeDefect.jl:53`, re-add the unused name:
#       using LinearAlgebra: norm, eigvals          # `norm` never called
#     Verified bite: the "ExplicitImports" testset goes RED —
#       check_no_stale_explicit_imports(PadeTaylor) === nothing  → Test Failed
#       (StaleImportsException: module SharedPadeDefect has stale explicit
#        imports for: norm).
#     Restored to `using LinearAlgebra: eigvals`.
#
#   Mutation B (JET — undefined-call / dispatch error)
#     In `src/SharedPadeDefect.jl`, inside `relative_defect`, insert a call to
#     a function that does not exist, e.g. `__nonexistent_jet_probe__(C)`.
#     Verified bite: the "JET.report_package" testset goes RED — the new
#     report (`UndefVarReport`/`MethodErrorReport` for the bogus name) is NOT
#     matched by the quality_diagnose allow-list, so `unexpected` is non-empty
#       @test isempty(unexpected)   → Test Failed (and the @info dumps it).
#     Restored by deleting the injected line.
#
#   Mutation C (Aqua — method ambiguity)
#     Append two ambiguous methods to a src module, e.g.
#       __amb(::Int, ::Any) = 1
#       __amb(::Any, ::Int) = 2     # __amb(1, 1) is ambiguous
#     Verified bite: `Aqua.test_all`'s "Method ambiguity" sub-testset goes RED
#     (Aqua lists the unresolved (Int, Int) pair).  Restored by deleting both.
#
# All three mutations restored before hand-off; `git diff src/` shows only the
# real `norm`-removal fix, none of the injections.
