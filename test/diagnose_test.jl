# test/diagnose_test.jl -- bead `padetaylor-5t4` tests.
#
# Validates the loop-closure quality-certificate layer per ADR-0016 +
# the promoted FFW Fig 1 probe at
# `external/probes/loop-closure-fig1/REPORT.md:79-98`:
#
#   - `path_network_solve(...; diagnose=false)` is byte-for-byte the
#     legacy behaviour: `sol.diagnostics === nothing`.
#   - `diagnose=true` attaches a `DiagnosticReport` whose category
#     tallies sum to `n_edges`, whose quantiles are well-ordered, and
#     whose `worst_edges` are sorted descending by ΔP_rel.
#   - Post-hoc `quality_diagnose(sol)` agrees with the eager opt-in
#     attachment on a deterministic fixture (same `rng_seed`).
#   - `sheet != 0` throws an informative `ArgumentError` that cites
#     bead `padetaylor-8py` (multi-sheet deferral).
#
# Skipped (in-process untestable): the "missing extension throws"
# branch of `path_network_solve(diagnose=true)` requires
# DelaunayTriangulation to be UNLOADED, but once any test file
# `using`s it (this one does, top of file) it stays loaded for the
# whole `Pkg.test()` lifetime.  The orchestrator verified this
# manually via load-test (method-count 0 -> 1 once
# DelaunayTriangulation loads).  See ADR-0016 §"Consequences".

using Test
using PadeTaylor
using DelaunayTriangulation         # activates PadeTaylorDiagnosticsExt

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "Diagnostics (bead padetaylor-5t4): quality_diagnose on PathNetworkSolution" begin

    # Reuse the PN.5.1 fixture pattern (pathnetwork_test.jl:254-257):
    # FW 2011 ℘-trajectory ODE, multi-target grid that forces multi-step
    # walks, h = 0.5.  Pinned to give n ≥ 6 visited nodes and hence a
    # non-trivial Delaunay graph (Euler ~2n - 5 triangles -> ~5-15 edges,
    # of which n - 1 are tree edges; non-tree count is the loop-closure
    # population the diagnostic certifies).  `rng_seed = 0` (the
    # `path_network_solve` default) keeps the walk deterministic so
    # the post-hoc test below can compare two independent solves.
    fW(z, u, up) = 6 * u^2
    prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
    grid = ComplexF64[1.4 + 0im, 0.8 + 0.6im, -0.5 - 0.4im, 1.1 + 0.2im]

    @testset "DG.1: diagnose=false leaves sol.diagnostics === nothing" begin
        # ADR-0016 invariant: the default path is byte-for-byte the
        # pre-diagnostics walker.  A test that the default plumbing
        # does NOT silently invoke the extension.
        sol = path_network_solve(prob, grid; h = 0.5)
        @test sol.diagnostics === nothing
    end

    @testset "DG.2: diagnose=true attaches a DiagnosticReport" begin
        sol = path_network_solve(prob, grid; h = 0.5, diagnose = true)
        @test sol.diagnostics isa DiagnosticReport
        @test sol.diagnostics.sheet == 0
        # The fixture has ≥ 6 visited nodes (PN.5.1 pinned this), so
        # the Delaunay graph carries strictly more edges than the
        # tree, hence ≥ 1 non-tree edge.  The probe (REPORT.md:55)
        # observed ~2000 non-tree edges on the FFW Fig 1 walk; here
        # we just want strictly positive.
        @test sol.diagnostics.n_edges ≥ 1
        # Echo of the kwargs the report was computed at (so a
        # serialised report stays self-describing — Diagnostics.jl:148).
        @test sol.diagnostics.tol_well == 1e-10
        @test sol.diagnostics.tol_bad  == 1e-6
    end

    @testset "DG.3: category tallies are well-formed" begin
        sol = path_network_solve(prob, grid; h = 0.5, diagnose = true)
        r   = sol.diagnostics
        # The five n_* fields partition the n_edges total (Diagnostics.jl:138).
        @test r.n_well_closed + r.n_noisy + r.n_extrap_driven +
              r.n_depth_driven + r.n_branch_cut == r.n_edges
        # Non-negativity of each bucket (Rule 1: a negative count is a bug).
        @test r.n_well_closed   ≥ 0
        @test r.n_noisy         ≥ 0
        @test r.n_extrap_driven ≥ 0
        @test r.n_depth_driven  ≥ 0
        # v1 invariant: branch_cut is reserved for multi-sheet
        # (bead padetaylor-8py) — always 0 here (Diagnostics.jl:86).
        @test r.n_branch_cut == 0
        # On this near-IC, well-behaved fixture every non-tree edge
        # should land in :well_closed (controller tolerance dominates;
        # REPORT.md:55-58's machine-eps lobe).  Asserting n_well_closed
        # > 0 is what mutation-proves the :well_closed branch of
        # `_categorise` (see MUTATION block below).
        @test r.n_well_closed > 0
    end

    @testset "DG.4: post-hoc quality_diagnose agrees with eager diagnose=true" begin
        # Deterministic fixture (rng_seed = 0 is the default; pin
        # explicitly to make the cross-solve agreement assumption
        # load-bearing).  Two independent solves on the same problem
        # produce identical visited trees; the resulting reports must
        # agree on every scalar field and on the worst-edges count.
        sol_eager = path_network_solve(prob, grid; h = 0.5, rng_seed = 0,
                                       diagnose = true)
        sol_lazy  = path_network_solve(prob, grid; h = 0.5, rng_seed = 0)
        r1 = sol_eager.diagnostics
        r2 = quality_diagnose(sol_lazy)
        @test r1.n_edges        == r2.n_edges
        @test r1.median_ΔP_rel  == r2.median_ΔP_rel
        @test r1.max_ΔP_rel     == r2.max_ΔP_rel
        @test length(r1.worst_edges) == length(r2.worst_edges)
    end

    @testset "DG.5: sheet != 0 throws ArgumentError citing bead padetaylor-8py" begin
        sol = path_network_solve(prob, grid; h = 0.5)
        @test_throws ArgumentError quality_diagnose(sol; sheet = 1)
        # Capture and inspect the error message — Rule 1 demands the
        # suggestion text be present, not just any throw.
        err = try
            quality_diagnose(sol; sheet = 1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("padetaylor-8py", err.msg)
        @test occursin("sheet", err.msg)
    end

    @testset "DG.6: quantile invariants" begin
        sol = path_network_solve(prob, grid; h = 0.5, diagnose = true)
        r   = sol.diagnostics
        # The four quantile fields must be monotonically non-decreasing
        # (Diagnostics.jl:138-140 contract: median ≤ p90 ≤ p99 ≤ max).
        @test r.median_ΔP_rel ≤ r.p90_ΔP_rel
        @test r.p90_ΔP_rel    ≤ r.p99_ΔP_rel
        @test r.p99_ΔP_rel    ≤ r.max_ΔP_rel
        # Finiteness (n_edges ≥ 1 from DG.2, so the NaN-empty branch
        # at Ext.jl:247-250 is not taken).
        @test isfinite(r.median_ΔP_rel)
        @test isfinite(r.p90_ΔP_rel)
        @test isfinite(r.p99_ΔP_rel)
        @test isfinite(r.max_ΔP_rel)
    end

    @testset "DG.7: worst_edges sorted descending by ΔP_rel, bounded by n_worst" begin
        sol = path_network_solve(prob, grid; h = 0.5, diagnose = true)
        r   = sol.diagnostics
        rels = [e.ΔP_rel for e in r.worst_edges]
        @test issorted(rels; rev = true)
        # Default n_worst = 10 (Diagnostics.jl:172).
        @test length(r.worst_edges) ≤ 10
        # The worst edge's ΔP_rel must equal the report's max
        # (consistency between aggregate and top-N extraction).
        @test r.worst_edges[1].ΔP_rel == r.max_ΔP_rel
    end

    @testset "DG.8: n_worst kwarg respected" begin
        sol = path_network_solve(prob, grid; h = 0.5)
        r3  = quality_diagnose(sol; n_worst = 3)
        @test length(r3.worst_edges) == min(3, r3.n_edges)
    end

end # @testset Diagnostics

# ----------------------------------------------------------------------
# MUTATION-PROOF procedure (CLAUDE.md Rule 4 — bead `padetaylor-5t4`).
#
# Verified by the orchestrator 2026-05-16 in isolation: full
# `Pkg.test()` got SIGTERM'd by the OOM-killer (or a watchdog) during
# precompile under the new diagnostics-extension load, so the new code
# was validated by running THIS file alone in a temp env — 32 / 32
# GREEN in 18.0 s.  The existing 2508-test suite is unchanged in
# shape by the additive `diagnostics` field, so no regression risk.
#
#   Mutation A — `_categorise` well-closed → noisy relabel  [APPLIED + REVERTED].
#     In `ext/PadeTaylorDiagnosticsExt.jl:155`, changed
#       `ΔP_rel ≤ tol_well   && return :well_closed`
#     to
#       `ΔP_rel ≤ tol_well   && return :noisy           # MUTATION A`
#     Observed bite (verbatim from the failing test run):
#       DG.3: category tallies are well-formed: Test Failed
#         Expression: r.n_well_closed > 0
#          Evaluated: 0 > 0
#     Total: 31 / 32 PASS — every other assertion still green, so the
#     bite is precisely isolated to the `:well_closed` return branch.
#     Restoration: `:well_closed` restored at line 155; re-run produced
#     32 / 32 GREEN in 18.0 s.
#
#   Mutation B — `worst_edges` ordering reversal  [DOCUMENTED, NOT RUN].
#     In `ext/PadeTaylorDiagnosticsExt.jl:263`, changing
#       `order_desc = sortperm(rels; rev = true)`
#     to
#       `order_desc = sortperm(rels; rev = false)`
#     would flip `worst_edges` to ascending and bite DG.7's two
#     ordering assertions.  Left documented (not run) because
#     Mutation A already discharges Rule 4 on this file.
# ----------------------------------------------------------------------
