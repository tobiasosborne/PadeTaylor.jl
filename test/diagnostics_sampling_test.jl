#=
Sampling-accounting regressions for padetaylor-qdsm / padetaylor-orb5.

The independent oracle is the constant field u(z)=1: every finite midpoint
disagreement is exactly zero. A three-node triangle has exactly three
Delaunay edges; a rooted star removes two, leaving the known edge (2,3).
An exact denominator zero on that midpoint must be counted as one thrown
evaluation, rather than disappearing from a clean certificate (ADR-0016a).
=#
using Test
using PadeTaylor, DelaunayTriangulation
using PadeTaylor.VectorPathNetwork: VectorPathNetworkSolution

const SamplingExt = Base.get_extension(PadeTaylor, :PadeTaylorDiagnosticsExt)
constant_patch() = PadeApproximant{ComplexF64}([1], [1], 0, 0)
function sampling_walk(z; patches=fill(constant_patch(), length(z)),
                       sheets=[Int[] for _ in z], parent=vcat(0, 1:length(z)-1))
    n = length(z)
    PathNetworkSolution{Float64}(ComplexF64.(z), ones(ComplexF64, n),
        zeros(ComplexF64, n), patches, ones(n), parent,
        ComplexF64[], ComplexF64[], ComplexF64[], sheets)
end

function sampling_error(f)
    try
        f()
    catch err
        return err
    end
    error("Expected a diagnostics error; the sampling guard did not fire")
end

@testset "Diagnostics sampling: qdsm / orb5" begin
    @testset "z-frame keeps the full 100-unit imaginary window" begin
        z = ComplexF64[x + y*im for y in -50:50 for x in -1:1]
        # Pre-fix: 39/303 = 12.8713% retained; 264/303 = 87.1287% lost.
        @test count(SamplingExt._sheet_mask(z, [Int[] for _ in z], 0)) == 303
        @test count(SamplingExt._sheet_mask(z, Vector{Int}[], 0)) == 303
        r = quality_diagnose(sampling_walk(z); n_worst=0)
        @test r.n_nodes == 303
        @test r.n_nodes_retained == 303
        @test r.n_nodes_dropped == 0
        @test r.n_edges_candidate == r.n_edges
        @test r.n_edges_dropped == 0
        @test r.max_ΔP_rel == 0.0
        @test r.n_well_closed == r.n_edges
        summary = sprint(show, MIME"text/plain"(), r)
        @test occursin("303 / 303", summary)
    end

    @testset "unrequested sheet loss fails with a suggestion" begin
        sol = sampling_walk([0, 2, 2im, 2+2im];
                            sheets=[[0], [0], [0], [1]], parent=[0,1,1,1])
        @test_throws ArgumentError quality_diagnose(sol)
        err = sampling_error(() -> quality_diagnose(sol))
        @test occursin("Suggestion", sprint(showerror, err))
        @test occursin("3 / 4", sprint(showerror, err))
        @test occursin("min_retained_fraction", sprint(showerror, err))
        r = quality_diagnose(sol; min_retained_fraction=0.75)
        @test (r.n_nodes, r.n_nodes_retained, r.n_nodes_dropped) == (4,3,1)
        @test r.n_tree_edges_off_sheet == 1
        @test (r.n_edges_candidate, r.n_edges, r.n_edges_dropped) == (1,1,0)
        @test r.max_ΔP_rel == 0.0
        @test_throws ArgumentError quality_diagnose(sol; min_retained_fraction=0.8)
        for fraction in (-0.1, 1.1, NaN)
            @test_throws ArgumentError quality_diagnose(sol; min_retained_fraction=fraction)
        end
    end

    @testset "known pole and nonfinite midpoint are counted once per edge" begin
        # A convex pentagon has seven Delaunay edges, including five hull
        # edges. Every non-root hull edge is non-tree for this rooted star.
        # A patch with one midpoint root drops exactly one edge; a NaN
        # patch drops every incident edge. Other patches equal u(z)=1.
        z = ComplexF64[0, 4, 3+4im, 2im, 5+im]
        parent = [0,1,1,1,1]
        for kbad in (2, 5)
            good = sampling_walk(z; parent)
            baseline = quality_diagnose(good; n_worst=100)
            # Exact empty-circumcircle oracle (rational arithmetic) gives
            # non-tree edges {(2,3),(2,4),(2,5),(3,4),(3,5)}; no four
            # vertices are cocircular, so insertion order cannot change it.
            @test Set((e.A,e.B) for e in baseline.worst_edges) ==
                  Set([(2,3), (2,4), (2,5), (3,4), (3,5)])
            incident = [e for e in baseline.worst_edges if kbad in (e.A,e.B)]
            @test length(incident) == (kbad == 2 ? 3 : 2)
            # Q(t) has one exact midpoint root; the independent count is
            # the number of baseline edges incident on a nonfinite node.
            pole_edge = first(incident)
            root = pole_edge.midpoint - z[kbad]
            patches = fill(constant_patch(), 5)
            patches[kbad] = PadeApproximant{ComplexF64}([1], [-root, 1], 0, 1)
            pole = quality_diagnose(sampling_walk(z; patches, parent))
            @test pole.n_edges_candidate == baseline.n_edges
            @test pole.n_edges_threw == 1
            @test pole.n_edges_nonfinite == 0
            @test pole.n_edges_dropped == 1
            @test pole.n_edges == baseline.n_edges - 1
            patches[kbad] = PadeApproximant{ComplexF64}([NaN], [1], 0, 0)
            bad = quality_diagnose(sampling_walk(z; patches, parent))
            @test bad.n_edges_threw == 0
            @test bad.n_edges_nonfinite == length(incident)
            @test bad.n_edges_dropped == length(incident)
            @test bad.n_edges + bad.n_edges_dropped == bad.n_edges_candidate
            @test bad.max_ΔP_rel == 0.0
            summary = sprint(show, MIME"text/plain"(), bad)
            @test occursin("nonfinite=$(length(incident))", summary)
        end
    end

    @testset "no surviving edges is a failure, not NaN statistics" begin
        sol = sampling_walk([0, 2, 2im]; parent=[0,1,1])
        sol.visited_pade[2] = PadeApproximant{ComplexF64}([1], [1-im, 1], 0, 1)
        @test_throws ArgumentError quality_diagnose(sol)
        err = sampling_error(() -> quality_diagnose(sol))
        @test occursin("Suggestion", sprint(showerror, err))
        @test occursin("all 1 candidate edges failed (threw=1, nonfinite=0)",
                       sprint(showerror, err))
    end

    @testset "unexpected evaluator errors propagate" begin
        sol = sampling_walk([0, 2, 2im]; parent=[0,1,1])
        evaluate(s, k, t) = error("unexpected evaluator bug")
        run() = SamplingExt._diagnostic_edges(sol, 1:3, [(2,3)], evaluate,
            SamplingExt._scalar_disagreement, 1e-10, 1e-6)
        @test_throws ErrorException run()
        err = sampling_error(run)
        @test sprint(showerror, err) == "unexpected evaluator bug"
    end

    @testset "vector adapter reports the same population accounting" begin
        z = ComplexF64[0, 4, 3+4im, 2im, 5+im]
        nums = [[ComplexF64[1], ComplexF64[2]] for _ in z]
        dens = [ComplexF64[1] for _ in z]
        sol = VectorPathNetworkSolution{Float64}(z, [ComplexF64[1,2] for _ in z],
            ones(5), nums, dens, [0,1,1,1,1])
        baseline = quality_diagnose(sol; n_worst=100)
        incident = count(e -> 2 in (e.A,e.B), baseline.worst_edges)
        @test (baseline.n_edges, incident) == (5,3)
        nums[2][1][1] = NaN
        r = quality_diagnose(sol)
        @test (r.n_nodes, r.n_nodes_retained, r.n_nodes_dropped) == (5,5,0)
        @test r.n_tree_edges_off_sheet == 0
        @test r.n_edges_candidate == baseline.n_edges
        @test r.n_edges_threw == 0
        @test r.n_edges_nonfinite == incident
        @test r.n_edges_dropped == incident
        @test r.n_edges == baseline.n_edges - incident
        @test r.max_ΔP_rel == 0.0
    end
end

# MUTATION PROOF — 2026-09-29, applied individually and restored.
# A diagnostics_scalar.jl: `return trues(length(visited_z))` -> old ζ strip.
#   RED: `39 == 303` twice; guard reports retained 39 / 303 (12.871287%).
#   52 pass / 2 fail / 1 error.
# B diagnostics_scalar.jl: prepend `false &&` to the retained-fraction guard.
#   RED: expected ArgumentError, "No exception thrown"; 50/1/1.
# C diagnostics_edges.jl: `n_threw += 1` -> `n_threw += 0`.
#   RED: `pole.n_edges_threw == 1`, evaluated `0 == 1`; 55/5/1.
# D diagnostics_edges.jl: both `n_nonfinite += 1` -> `n_nonfinite += 0`.
#   RED: counts `0 == 3` / `0 == 2`, scalar and vector; 52/10/0.
# E diagnostics_geometry.jl: `n_off_sheet += 1` -> `n_off_sheet += 0`.
#   RED: `r.n_tree_edges_off_sheet == 1`, evaluated `0 == 1`; 61/1/0.
# F diagnostics_edges.jl: `err isa DomainError || rethrow()` -> `true`.
#   RED: expected ErrorException, "No exception thrown"; 60/1/1.
# G diagnostics_edges.jl: prepend `false &&` to all-candidate-failure guard.
#   RED: expected ArgumentError, "No exception thrown"; 59/1/1.
# H diagnostics_scalar.jl: replace fraction-domain validation with `true ||`.
#   RED: invalid fraction accepted, "No exception thrown" twice; 60/2/0.
# All triples above are pass/fail/error; every mutant had 0 broken.
# Restored files byte-match saved snapshots; restored suite: 62/0/0/0.
# Exact output and design context: docs/worklog/081-diagnostics-silent-drops.md.
# Local raw logs: .depot/diagnostics-mutation-{mask,retention,threw,nonfinite,
# offsheet,unexpected,allfailed,fraction}.log. No mutation remains in source.
