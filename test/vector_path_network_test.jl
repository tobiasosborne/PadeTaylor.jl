"""
V7 tests for `PadeTaylor.VectorPathNetwork` and
`PadeTaylor.VectorPoleField` (bead `padetaylor-0ln.16`; v0.2 plan row
V7).

`vector_path_network_solve(prob, targets; …)` is the *minimal*
Stage-1 vector analogue of v0.1's `path_network_solve`: seed at the
problem's IC, and for each target walk from the nearest visited node
toward it via a 5-direction wedge, picking the wedge direction by the
pole-avoidance criterion min-‖y‖.  Each landed node stores the
shared-`Q` approximant `(numerators, denominator)` produced by
`VectorStepper.vector_pade_step_with_pade!`.

`extract_poles_shared_q(sol)` roots every visited node's *single*
shared denominator `Q`, maps each root `t` to the z-plane via
`z = z_node + t·h_node`, and clusters across nodes — the shared-`Q`
analogue of v0.1's `PoleField.extract_poles`, but cleaner: one `Q` per
node rather than `d` per-component denominators.

These tests (`VPN.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    VPN.1.1  known-pole oracle   — d=3 Riccati system y_i' = -y_i²/c_i
                                   has the exact solution
                                   y_i(z) = c_i/(z-p); all components
                                   blow up at the same z=p.  Walk a
                                   network over a region containing p
                                   and assert extract_poles_shared_q
                                   recovers p to ~1e-8.
    VPN.1.2  A_4^(1) walk        — walk the Noumi–Yamada A_4^(1)
                                   system; the visited tree is finite,
                                   non-empty, visited_parent encodes a
                                   valid tree, extracted poles finite.
    VPN.1.3  pole-extraction     — roots of a hand-built shared Q map
                                   correctly via z = z_node + root·h;
                                   cross-node clustering merges dupes.
    VPN.1.4  wedge avoidance     — on VPN.1.1's system the walk steers
                                   around the pole (no visited node
                                   lands inside a small disc of p).
    VPN.1.5  failure modes       — empty targets, degenerate prob;
                                   mutation-proof recorded at end.

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_path_network_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution
using PadeTaylor.VectorPoleField: extract_poles_shared_q

# -----------------------------------------------------------------------------
# Oracle 1 — the d-component Riccati pole system.
#
# y_i' = -y_i²/c_i  has the exact solution  y_i(z) = c_i/(z-p),  since
#   y_i' = -c_i/(z-p)²  and  -y_i²/c_i = -(c_i/(z-p))²/c_i = -c_i/(z-p)².
# Given the IC y_i(z_0) = c_i/(z_0-p), the pole is p = z_0 - c_i/y_i(z_0).
# Every component blows up at the SAME z=p — the shared-pole property
# the shared-Q machinery exploits.
# -----------------------------------------------------------------------------

# Build a Riccati-pole VectorPadeTaylorProblem with a chosen pole `p`,
# IC point `z0`, and per-component constants `cs`.  IC is y_i(z0).
function riccati_pole_problem(p, z0, cs; order = 24)
    f = (z, y) -> [-(y[i]^2) / cs[i] for i in eachindex(y)]
    y0 = [cs[i] / (z0 - p) for i in eachindex(cs)]
    return VectorPadeTaylorProblem(f, y0, (z0, z0 + 1.0 + 0im); order = order)
end

@testset "VectorPathNetwork V7 (VPN.*)" begin

    # -------------------------------------------------------------------------
    # VPN.1.1 — known-pole oracle.  d=3 Riccati system, exact pole p.
    # -------------------------------------------------------------------------
    @testset "VPN.1.1 known-pole oracle" begin
        p   = 1.0 + 0.6im                       # the exact movable pole
        z0  = -0.4 + 0.0im                      # IC far enough away
        cs  = ComplexF64[1.0, 0.7, -1.3]        # per-component constants
        prob = riccati_pole_problem(p, z0, cs; order = 24)

        # Targets: a small grid in the complex plane containing p.
        targets = ComplexF64[x + y * im
                             for x in 0.2:0.4:1.8 for y in -0.2:0.4:1.4]

        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.25)
        @test sol isa VectorPathNetworkSolution
        @test length(sol.visited_z) ≥ 1

        poles = extract_poles_shared_q(sol; radius_t = 6.0,
                                       cluster_atol = 0.2, min_support = 2)
        @test !isempty(poles)
        # The recovered pole nearest p must be close to the exact p.
        d = minimum(abs(q - p) for q in poles)
        @test d < 1.0e-8
    end

    # -------------------------------------------------------------------------
    # VPN.1.2 — A_4^(1) Noumi–Yamada walk.
    # -------------------------------------------------------------------------
    @testset "VPN.1.2 A_4^(1) walk" begin
        n = 2                                   # A_4^(1): d = 2n+1 = 5
        d = 2n + 1
        α = fill(1.0 / d + 0im, d)              # Σα = 1, k=1 normalisation
        t0 = 0.3 + 0.0im
        f0 = fill(t0 / d, d)                    # Σf0 = t0 (Type C point)
        ny = NoumiYamadaProblem(n; α = α, f0 = f0,
                                tspan = (t0, t0 + 1.0 + 0im), order = 20)

        targets = ComplexF64[x + y * im
                             for x in 0.5:0.5:1.5 for y in -0.5:0.5:0.5]
        sol = vector_path_network_solve(ny.problem, targets;
                                        order = 20, h = 0.2)

        @test sol isa VectorPathNetworkSolution
        @test 1 ≤ length(sol.visited_z) < 10_000     # finite, non-empty
        # visited_parent encodes a valid tree: root has parent 0, every
        # other node points to an earlier (already-visited) index.
        @test sol.visited_parent[1] == 0
        for k in 2:length(sol.visited_parent)
            @test 1 ≤ sol.visited_parent[k] < k
        end
        # Every node carries d numerators + one shared denominator.
        for num in sol.visited_numerators
            @test length(num) == d
        end
        poles = extract_poles_shared_q(sol)
        @test all(isfinite, poles)               # finite or empty — sensible
    end

    # -------------------------------------------------------------------------
    # VPN.1.3 — pole-extraction mechanics.
    # -------------------------------------------------------------------------
    @testset "VPN.1.3 pole-extraction mechanics" begin
        # Hand-build a one-node solution whose shared Q = (1 - t) has a
        # single root t = 1.  With z_node = 2+0im, h = 0.5, the pole
        # must map to z = z_node + 1·h = 2.5.
        z_node = 2.0 + 0.0im
        h_node = 0.5 + 0.0im
        Q = ComplexF64[1.0, -1.0]                # Q(t) = 1 - t, root t=1
        N = [ComplexF64[1.0]]                    # numerator (d=1)
        sol1 = VectorPathNetworkSolution{Float64}(
            [z_node], [ComplexF64[1.0]], [h_node], [N], [Q], [0])
        poles = extract_poles_shared_q(sol1; radius_t = 5.0, min_support = 1)
        @test length(poles) == 1
        @test abs(poles[1] - (z_node + 1.0 * h_node)) < 1.0e-12

        # Two nodes see the SAME physical pole at z = 2.5; from node 1 the
        # root is t=1 (h=0.5 → 2+0.5), from node 2 (z=2.25, h=0.25) the
        # root is t=1 (2.25+0.25 = 2.5).  Clustering must merge them.
        z2  = 2.25 + 0.0im
        h2  = 0.25 + 0.0im
        sol2 = VectorPathNetworkSolution{Float64}(
            [z_node, z2], [ComplexF64[1.0], ComplexF64[1.0]],
            [h_node, h2], [N, N], [Q, Q], [0, 1])
        merged = extract_poles_shared_q(sol2; radius_t = 5.0,
                                        cluster_atol = 0.1, min_support = 2)
        @test length(merged) == 1               # two estimates → one cluster
        @test abs(merged[1] - 2.5) < 1.0e-12

        # Cross-node support filter: a third node carries a SPURIOUS root
        # — Q3(t) = 1 - t/3 has root t = 3, which maps to a z-plane
        # location no other node confirms.  With min_support = 2 that
        # lone-node artefact must be dropped; only the physical pole at
        # z = 2.5 (seen by nodes 1 & 2) survives.
        z3  = 0.0 + 0.0im
        h3  = 1.0 + 0.0im
        Q3  = ComplexF64[1.0, -1.0 / 3.0]        # root t = 3 → z = 3.0
        sol3 = VectorPathNetworkSolution{Float64}(
            [z_node, z2, z3],
            [ComplexF64[1.0], ComplexF64[1.0], ComplexF64[1.0]],
            [h_node, h2, h3], [N, N, N], [Q, Q, Q3], [0, 1, 0])
        filtered = extract_poles_shared_q(sol3; radius_t = 5.0,
                                          cluster_atol = 0.1, min_support = 2)
        @test length(filtered) == 1              # the lone-node root dropped
        @test abs(filtered[1] - 2.5) < 1.0e-12
        # With the filter disabled (min_support = 1) the spurious root
        # surfaces — proving it is genuinely present and the filter is
        # what removes it.
        unfiltered = extract_poles_shared_q(sol3; radius_t = 5.0,
                                            cluster_atol = 0.1,
                                            min_support = 1)
        @test length(unfiltered) == 2
    end

    # -------------------------------------------------------------------------
    # VPN.1.4 — wedge pole-avoidance.  The min-‖y‖ wedge selection steers
    # the walk around the pole; no visited node lands ON p.
    # -------------------------------------------------------------------------
    @testset "VPN.1.4 wedge pole-avoidance" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.2:0.4:1.8 for y in -0.2:0.4:1.4]
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.25)
        # No visited node lands inside a small disc of the pole — the
        # min-‖y‖ wedge selection avoids the singularity.
        min_dist = minimum(abs(z - p) for z in sol.visited_z)
        @test min_dist > 1.0e-3
    end

    # -------------------------------------------------------------------------
    # VPN.1.5 — failure modes.
    # -------------------------------------------------------------------------
    @testset "VPN.1.5 failure modes" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)

        # Empty targets — nothing to walk to.
        @test_throws ArgumentError vector_path_network_solve(
            prob, ComplexF64[]; order = 24, h = 0.25)

        # Non-positive h — degenerate step.
        @test_throws ArgumentError vector_path_network_solve(
            prob, ComplexF64[0.5 + 0.0im]; order = 24, h = 0.0)

        # Wrong wedge_angles count — FW 2011 fixes the 5-direction wedge.
        @test_throws ArgumentError vector_path_network_solve(
            prob, ComplexF64[0.5 + 0.0im]; order = 24, h = 0.25,
            wedge_angles = [0.0, 0.1])

        # extract_poles_shared_q on an empty solution — no nodes.
        empty_sol = VectorPathNetworkSolution{Float64}(
            ComplexF64[], Vector{ComplexF64}[], ComplexF64[],
            Vector{Vector{ComplexF64}}[], Vector{ComplexF64}[], Int[])
        @test isempty(extract_poles_shared_q(empty_sol))
    end

end

# =============================================================================
# Mutation-proof record (CLAUDE.md Rule 4).
#
# Procedure: perturb `src/VectorPathNetwork.jl` / `src/VectorPoleField.jl`,
# rerun this file, confirm RED, restore.  3 meaningful mutations applied.
#
#   M1 — VectorPoleField: wrong z-plane root mapping.  Change
#        `C(z_node + h_node * t_C)`  →  `C(z_node + t_C)`  (drop the
#        h factor that converts the rescaled root t* to a z-plane
#        displacement).
#        Expected: VPN.1.1 (pole off by the missing scaling) and
#        VPN.1.3 (mapped pole ≠ 2.5) go RED.
#        Result: RED — 6 assertions bit (4 fail + 2 error): VPN.1.1
#        pole-recovery, VPN.1.3 single-node mapping + two-node merged
#        location + the three-node filtered/unfiltered locations.
#        Restored to GREEN.
#
#   M2 — VectorPathNetwork: ignore the min-‖y‖ wedge criterion.  Replace
#        the `_wedge_step` selection loop with a fixed pick of wedge
#        index 3 (straight at the goal direction — no pole avoidance).
#        Expected: the walk drives straight into the shared movable
#        pole instead of steering around it.
#        Result: RED — 2 errors (VPN.1.1 and VPN.1.4): the straight-line
#        walk steps onto the pole, `vector_pade_step_with_pade!` throws
#        `DomainError` (Q(1) ≈ 0), the walk fails loud (Rule 1).  This
#        is exactly the failure the min-‖y‖ criterion exists to prevent.
#        Restored to GREEN.
#
#   M3 — VectorPoleField: drop the cross-node cluster support filter.
#        Change the final comprehension guard
#        `length(support[j]) ≥ min_support`  →  `true`.
#        Expected: VPN.1.3's two-node `min_support = 2` case no longer
#        merges-and-filters correctly — a single-node estimate that
#        should be discarded now leaks through, so the asserted pole
#        count `length(merged) == 1` can break when stray candidates
#        appear; the load-bearing role of the cross-node filter.
#        Result: RED — VPN.1.3's `length(filtered) == 1` assertion bit
#        (with the filter disabled the lone-node spurious root at z = 3
#        leaks through, so the pole count is 2 not 1).  Restored to
#        GREEN.
#
# Certified bites: M1 → 6 RED, M2 → 2 RED, M3 → 1 RED.  Restored to
# GREEN after each mutation.
# =============================================================================
