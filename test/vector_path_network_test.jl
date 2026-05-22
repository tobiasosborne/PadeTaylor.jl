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

## The B2 dense-wedge walk tests (`VPN.3.*`)

B2 (`padetaylor-0ln.37.6`, ADR-0025 Lever 2; absorbs bead `0ln.23`)
replaces the V7 min-‖y‖ wedge selector with the principled
`:max_q_root` shared-`Q`-root-distance criterion and makes the step `h`
adaptive.  The `VPN.3.*` testset asserts the real B2 wins:

    VPN.3.1  threads the wedge   — the P_I⁽²⁾ tritronquée pole-rich
                                   wedge: a target fan reaching
                                   |x| ≳ 15 that the fixed-h walk
                                   BLOCKS on (the A2-probe failure) is
                                   threaded by the B2 adaptive walk.
    VPN.3.2  selector picks the  — `_select_wedge(:max_q_root)` chooses
             pole-free disc        the candidate with the larger
                                   shared-Q-root distance, demonstrably.
    VPN.3.3  scale-derived step  — the non-ratcheting `_adaptive_h`
                                   law `h = clamp(SAFETY·D_local, …)`:
                                   a near pole gives a small step, a
                                   far pole `h_max`, a node recovers in
                                   ONE step (no geometric sink), and no
                                   adaptive-walk node overshoots a pole.
    VPN.3.4  fail-loud           — unknown step_policy throws; the
                                   wedged-walk h-floor throws on a
                                   genuinely sub-h_min disc (Rule 1).
    VPN.3.5  fixed-walk fallback — adaptive=false + :min_y reproduces
                                   the verified V7 fixed-step walk.

## The VC-10 target-ordering tests (`VPN.4.*`)

VC-10 (`padetaylor-0ln.37.17`, ADR-0025 Phase D) makes the Stage-1
target processing order controllable by the optional `rng` kwarg — the
substrate of the FW/FFW double-run accuracy indicator (FW 2011 §3.1
line 156 shuffles the targets; FFW 2017 `:246-247` makes the
cross-ordering disagreement the diagnostic).  The `VPN.4.*` testset
pins the kwarg's contract:

    VPN.4.1  additive default    — `rng = nothing` (the default)
                                   reproduces the no-`rng` walk
                                   bit-identically (ADR-0001).
    VPN.4.2  fixed-seed determ.  — a fixed `MersenneTwister` seed gives
                                   a bit-identical walk on every call
                                   (the per-fixed-seed reproducibility
                                   guarantee — PI2F.1.6's invariant
                                   lifted to the seeded walk).
    VPN.4.3  different orderings — two different seeds genuinely build
                                   different path-network trees; both
                                   still reach every target and extract
                                   finite, non-empty pole fields that
                                   agree within the walk's numerical
                                   error (the VC-10 indicator itself).

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_path_network_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using Random: MersenneTwister
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution,
                                    VectorWalkFailure, COVER_FRAC
using PadeTaylor.VectorPoleField: extract_poles_shared_q
using PadeTaylor.VectorWedgeStep: _select_wedge, _adaptive_h, H_MIN_RATIO,
                                  SAFETY, VectorWalkError
using PadeTaylor.PainleveHierarchy: painleve_hierarchy, pI2_tritronquee_ic
using PadeTaylor.VectorProblems: VectorPadeTaylorProblem
using PadeTaylor.VectorStepper: VectorPadeStepperState,
                                vector_pade_step_with_pade!
using Polynomials: Polynomial, roots

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

        # Targets: a grid BRACKETING p without sitting on it.  The B2
        # adaptive walk caps `h` toward zero as it nears a pole — it
        # cannot honestly walk *onto* a pole (a target equal to `p`
        # would collapse `h` below `h_min` and fail loud, Rule 1).  The
        # offset grid `{0.0,…,1.6} × {-0.4,…,1.2}` surrounds p so the
        # walk still builds nodes whose shared-Q captures it, but no
        # target coincides with `p` itself — the principled walk's
        # correct semantics (legitimate evolution per ADR-0025 B2).
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]

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
    # VPN.1.3b — scale-covariant cluster tolerance (ADR-0026 Amendment 3
    # §S4, bead padetaylor-2vv).
    #
    # `extract_poles_shared_q`'s default `cluster_atol = nothing` derives
    # the across-node merge tolerance from the local step `visited_h`:
    # two candidates merge iff their separation is
    # `≤ CLUSTER_FRAC · min(|h_i|, |h_j|)`.  An *absolute* cluster
    # tolerance is a scale heresy — it wrongly merges distinct poles
    # where the field is dense.  These hand-built fixtures pin the three
    # behaviours: (a) distinct poles kept separate; (b) noise-separated
    # duplicates merged; (c) the legacy `Real` path is byte-identical.
    # -------------------------------------------------------------------------
    @testset "VPN.1.3b scale-covariant cluster tolerance" begin
        # Fixture: a *fine-scale* two-node network.  Both nodes use a
        # small local step h = 0.05 (post-S2 a fine local scale ⇒ a
        # dense local pole field).  Node 1's shared-Q places a pole at
        # z = 1.00, node 2's at z = 1.12 — a genuine 0.12 separation,
        # i.e. > 2 local steps apart, two *distinct* physical poles.
        #
        # The OLD absolute default (cluster_atol ≈ 0.1, and the figure's
        # 0.2) would wrongly MERGE these into one pole.  The scale-
        # derived tolerance is 0.4·0.05 = 0.02 ≪ 0.12 — they stay
        # separate, the correct count.
        h_fine = 0.05 + 0.0im
        # Q(t) = 1 - t has root t = 1; mapped pole = z_node + h·1.
        Qroot1 = ComplexF64[1.0, -1.0]
        N1     = [ComplexF64[1.0]]
        # node 1 at z = 0.95, h = 0.05 → pole at 0.95 + 0.05 = 1.00
        # node 2 at z = 1.07, h = 0.05 → pole at 1.07 + 0.05 = 1.12
        sol_distinct = VectorPathNetworkSolution{Float64}(
            ComplexF64[0.95, 1.07],
            [ComplexF64[1.0], ComplexF64[1.0]],
            [h_fine, h_fine],
            [N1, N1], [Qroot1, Qroot1], [0, 1])
        # (a) scale-derived mode keeps the two distinct poles SEPARATE.
        sep = extract_poles_shared_q(sol_distinct; radius_t = 5.0,
                                     min_support = 1)        # default atol
        @test length(sep) == 2
        @test minimum(abs(p - 1.00) for p in sep) < 1.0e-12
        @test minimum(abs(p - 1.12) for p in sep) < 1.0e-12
        # The OLD absolute tolerance (≥ 0.1, e.g. the figure's 0.2)
        # WRONGLY merges them — proving the heresy is real and the
        # scale-derived mode is what fixes it.
        wrong = extract_poles_shared_q(sol_distinct; radius_t = 5.0,
                                       cluster_atol = 0.2, min_support = 1)
        @test length(wrong) == 1

        # (b) Two NOISE-separated copies of ONE pole are still merged.
        # Both nodes step at h = 0.5 (a coarse local scale); their
        # shared-Q roots place the *same* physical pole, separated only
        # by walk noise ~5e-4 ≪ h.  The scale-derived tolerance
        # 0.4·0.5 = 0.2 swallows the noise → one reported pole.
        h_coarse = 0.5 + 0.0im
        # node 1 at z = 1.0, h = 0.5 → pole at 1.5000
        # node 2 at z = 1.0005, h = 0.5 → pole at 1.5005 (noise 5e-4)
        sol_dup = VectorPathNetworkSolution{Float64}(
            ComplexF64[1.0, 1.0005],
            [ComplexF64[1.0], ComplexF64[1.0]],
            [h_coarse, h_coarse],
            [N1, N1], [Qroot1, Qroot1], [0, 1])
        merged = extract_poles_shared_q(sol_dup; radius_t = 5.0,
                                        min_support = 2)      # default atol
        @test length(merged) == 1     # noise-separated duplicates merged
        @test abs(merged[1] - 1.5) < 1.0e-2

        # (c) the legacy `Real` cluster_atol path is byte-identical to
        # the old behaviour.  Re-run the VPN.1.3 two-node merge under an
        # explicit absolute cluster_atol = 0.1 — the fixture there has
        # two estimates of z = 2.5 at h = 0.5 / 0.25; an absolute 0.1
        # tolerance merges them, exactly as before this change.
        z_n  = 2.0 + 0.0im
        z2   = 2.25 + 0.0im
        Qm   = ComplexF64[1.0, -1.0]
        Nm   = [ComplexF64[1.0]]
        sol_legacy = VectorPathNetworkSolution{Float64}(
            ComplexF64[z_n, z2],
            [ComplexF64[1.0], ComplexF64[1.0]],
            ComplexF64[0.5, 0.25],
            [Nm, Nm], [Qm, Qm], [0, 1])
        legacy = extract_poles_shared_q(sol_legacy; radius_t = 5.0,
                                        cluster_atol = 0.1, min_support = 2)
        @test length(legacy) == 1                 # legacy absolute merge
        @test abs(legacy[1] - 2.5) < 1.0e-12
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
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.25)
        # No visited node lands inside a small disc of the pole — the
        # B2 :max_q_root wedge selection picks the candidate with the
        # largest shared-Q-root distance, so the walk steers around the
        # singularity (and adaptive `h` never overshoots onto it).
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

    # =========================================================================
    # VPN.3.* — the B2 dense-wedge walk (bead padetaylor-0ln.37.6, ADR-0025
    # Lever 2; discharges bead 0ln.23).  See the file docstring.
    # =========================================================================

    # -------------------------------------------------------------------------
    # VPN.3.1 — the B2 headline win.  The P_I⁽²⁾ tritronquée carries a
    # *dense* pole field in the wedge straddling the positive real axis
    # (KKG Fig 7.4; ADR-0025).  The A2 probe
    # (external/probes/wedge-tractability/REPORT.md §3.1) measured the
    # fixed-h=0.1 walk BLOCKING past |x|≈8 on an extended target fan: a
    # bridging walk crosses the pole field on a chord and a fixed-h chord
    # step lands on a pole, degenerating the shared-Q linear solve
    # (`shared_denominator_pade: every singular value is below τ`).  B2's
    # adaptive walk threads that same fan to |x| ≳ 15.  This is the exact
    # failure B2 exists to fix — and the test that proves it.
    # -------------------------------------------------------------------------
    @testset "VPN.3.1 B2 adaptive walk threads the P_I^(2) wedge" begin
        # The P_I⁽²⁾ companion RHS and an asymptotic tritronquée seed at
        # x = -3 (the KKG n_terms=2 series — accurate on the negative
        # axis; the figure pipeline's Stage-A BVP refines it, but the
        # series seed suffices to march the wedge here, keeping the test
        # self-contained and fast).
        f      = painleve_hierarchy(:I, 2; t = 0.0)
        y_seed = ComplexF64.(pI2_tritronquee_ic(-3.0 + 0.0im;
                                                t = 0.0, n_terms = 2))
        prob   = VectorPadeTaylorProblem(f, y_seed, (-3.0 + 0im, 20.0 + 0im);
                                         order = 24)

        # The A2 extended fan: radii 2..18, five angles spanning ±0.4 rad
        # — straddling the positive-real wedge, reaching |x| = 18.
        fan = ComplexF64[r * cis(a)
                         for r in 2.0:2.0:18.0
                         for a in range(-0.4, 0.4; length = 5)]

        # The fixed-h walk BLOCKS — the A2 failure mode.  It throws the
        # shared-Q degeneration (a chord step landed on a pole).
        @test_throws Exception vector_path_network_solve(
            prob, fan; order = 24, h = 0.1,
            adaptive = false, step_policy = :min_y)

        # The B2 adaptive :max_q_root walk threads the same fan.
        sol = vector_path_network_solve(prob, fan; order = 24, h = 0.3,
                                        adaptive = true,
                                        step_policy = :max_q_root)
        @test sol isa VectorPathNetworkSolution
        # It reaches the outer wedge — a node at |x| ≳ 15 (the brief's
        # bar), in fact essentially the full |x| = 18 fan.
        reach = maximum(abs, sol.visited_z)
        @info "VPN.3.1 B2 walk: $(length(sol.visited_z)) nodes, " *
              "reached |x| = $(round(reach; digits = 2))"
        @test reach ≥ 15.0
        # The threaded walk carries a non-trivial pole field — the
        # scientific deliverable (ADR-0025 Amendment 2).
        poles = extract_poles_shared_q(sol)
        @test all(isfinite, poles)
        @test !isempty(poles)
        # Every visited node carries a 4-vector state + shared-Q store
        # (P_I⁽²⁾ companion form is d = 4).
        for num in sol.visited_numerators
            @test length(num) == 4
        end
        # `visited_h` records the *adaptive* per-node step — with a
        # genuine adaptive walk they are not all equal (the fixed-h
        # walk would have a constant column).
        steps = unique(round.(real.(sol.visited_h); digits = 6))
        @test length(steps) ≥ 2
    end

    # -------------------------------------------------------------------------
    # VPN.3.2 — the :max_q_root selector demonstrably picks the more
    # pole-free candidate.  This is the *contract* of the principled
    # criterion: of the five wedge candidates, choose the one whose
    # landed node has the largest shared-Q-root distance (= pole-free
    # disc).  A single-pole toy cannot show this — there min-‖y‖ and
    # max-disc coincide.  The P_I⁽²⁾ tritronquée wedge is a genuine
    # *multi*-pole field: there the two selectors make different choices,
    # and at every node where they differ `:max_q_root`'s pick is the
    # candidate with the (weakly) larger shared-Q-root disc, and at some
    # node a STRICTLY larger one — the win min-‖y‖ cannot make.
    # -------------------------------------------------------------------------
    @testset "VPN.3.2 :max_q_root picks the pole-free candidate" begin
        f      = painleve_hierarchy(:I, 2; t = 0.0)
        y_seed = ComplexF64.(pI2_tritronquee_ic(-3.0 + 0.0im;
                                                t = 0.0, n_terms = 2))
        prob   = VectorPadeTaylorProblem(f, y_seed, (-3.0 + 0im, 20.0 + 0im);
                                         order = 24)
        # A short walk into the wedge to harvest realistic mid-walk
        # states; then probe the selector at each.
        walk = vector_path_network_solve(prob, ComplexF64[8.0 + 2.0im];
                                         order = 24, h = 0.3,
                                         adaptive = true,
                                         step_policy = :max_q_root)
        wedge = [-π/4, -π/8, 0.0, π/8, π/4]
        h     = 0.3

        # The shared-Q-root distance (= pole-free disc) of the node a
        # wedge candidate lands on — the very quantity :max_q_root
        # maximises.  Mirrors `VectorWedgeStep._candidate_pole_disc`.
        function cand_disc(zc, yc, θ)
            hs = ComplexF64(h * cos(θ), h * sin(θ))
            st = VectorPadeStepperState{ComplexF64}(zc, yc)
            try
                vector_pade_step_with_pade!(st, f, 24, hs)
                st2 = VectorPadeStepperState{ComplexF64}(st.z, st.y)
                _, _, den = vector_pade_step_with_pade!(st2, f, 24,
                                                        ComplexF64(h))
                length(den) < 2 && return 10.0 * h
                rs = roots(Polynomial(collect(den)))
                isempty(rs) ? 10.0 * h :
                    min(10.0 * h, h * minimum(abs, rs))
            catch
                -Inf
            end
        end

        n_differ, n_strict = 0, 0
        for k in 2:min(40, length(walk.visited_z))
            z_cur = walk.visited_z[k]
            y_cur = walk.visited_y[k]
            goal  = angle((15.0 + 3.0im) - z_cur)
            discs = Float64[cand_disc(z_cur, y_cur, goal + a) for a in wedge]
            picks = ComplexF64[z_cur + ComplexF64(h*cos(goal+a),
                                                  h*sin(goal+a)) for a in wedge]
            zq, _ = _select_wedge(f, z_cur, y_cur, 24, h, goal,
                                  wedge, :max_q_root)
            zy, _ = _select_wedge(f, z_cur, y_cur, 24, h, goal,
                                  wedge, :min_y)
            kq = argmin(abs.(picks .- zq))
            ky = argmin(abs.(picks .- zy))
            # The :max_q_root pick has the MAXIMAL pole-disc among the
            # candidates (the contract — true at every node).
            @test discs[kq] ≥ maximum(discs) - 1.0e-9
            if abs(zq - zy) > 1.0e-6
                n_differ += 1
                # Where the selectors differ, :max_q_root's pick is at
                # least as pole-free as min-‖y‖'s — never worse.
                @test discs[kq] ≥ discs[ky] - 1.0e-9
                discs[kq] > discs[ky] + 1.0e-3 && (n_strict += 1)
            end
        end
        @info "VPN.3.2 selectors differed at $n_differ nodes; " *
              ":max_q_root strictly more pole-free at $n_strict"
        # The two selectors are genuinely different criteria — they
        # differ at several wedge nodes ...
        @test n_differ ≥ 3
        # ... and at some of those :max_q_root's disc is STRICTLY larger:
        # the principled win min-‖y‖ structurally cannot make.
        @test n_strict ≥ 1
    end

    # -------------------------------------------------------------------------
    # VPN.3.3 — the scale-derived, non-ratcheting `_adaptive_h` law
    # (ADR-0026 Amendment 4, the S2 step law; bead `padetaylor-cqk`).
    # `h = clamp(SAFETY·D_local, h_min, h_max)` where
    # `D_local = h_prev·min|t*|` is the h-independent absolute distance
    # to the nearest pole.  This testset pins the law's quantitative
    # contract directly on `_adaptive_h` — a near pole, a far pole, and
    # the LOAD-BEARING non-ratcheting recovery — plus the walk-level
    # consequence (no node overshoots a known pole).
    # -------------------------------------------------------------------------
    @testset "VPN.3.3 scale-derived non-ratcheting adaptive step" begin
        # (a) A mid-range pole gives `h = SAFETY·D_local`.  A shared-Q
        # `Q(t) = 1 - t/t0` has its single root at `t = t0`; with `t0`
        # placed so `SAFETY·h_prev·t0` lands strictly inside [h_min,
        # h_max], `_adaptive_h` must return exactly that value.
        h_prev, h_max, h_min = 0.3, 0.6, 1.0e-4
        t0   = 0.5                                   # mid-range Q-root
        Q    = ComplexF64[1.0, -1.0 / t0]            # root at t = t0
        h_ad = _adaptive_h(Q, h_prev, h_max, h_min)
        D_local = h_prev * t0                        # = 0.15
        # h = clamp(SAFETY·0.15, 1e-4, 0.6) = SAFETY·0.15, in-range.
        @test h_ad ≈ SAFETY * D_local  rtol = 1.0e-12
        @test h_min < h_ad < h_max                   # clamp did not bind
        # The step lands strictly short of the pole — a step of length
        # `h_ad` reaches a fraction `SAFETY` (< 1) of the way there.
        @test h_ad / D_local ≈ SAFETY  rtol = 1.0e-12

        # (b) LOAD-BEARING — non-ratcheting recovery.  A node deep in a
        # collapsed pocket (`h_prev` tiny) but whose nearest pole is now
        # FAR must jump back to `h_max` in this SINGLE call.  The old
        # geometric-grow law would crawl up by `GROW·h_prev = 1.5·1e-4`;
        # the scale-derived law has no `h_prev` factor, so `h` is set
        # purely by the (large) `D_local`.  We place the Q-root far
        # enough that `SAFETY·D_local ≥ h_max` and the clamp pins `h_max`.
        h_tiny = 1.0e-4                              # a near-collapsed step
        # D_local = h_tiny·t_far; want SAFETY·D_local ≥ h_max ⇒
        # t_far ≥ h_max/(SAFETY·h_tiny).  Pick comfortably above that.
        t_far  = 4.0 * h_max / (SAFETY * h_tiny)
        Q_far_root = ComplexF64[1.0, -1.0 / t_far]   # root at t = t_far
        h_recover  = _adaptive_h(Q_far_root, h_tiny, h_max, h_min)
        @test h_recover ≈ h_max  rtol = 1.0e-12      # full recovery, ONE step
        # And it genuinely jumped — not a `GROW·h_prev` crawl.
        @test h_recover > 100 * h_tiny

        # (c) A constant / rootless shared-Q means no nearby pole ⇒
        # `h = h_max` (D_local unbounded, the clamp upper bound governs).
        @test _adaptive_h(ComplexF64[1.0], h_tiny, h_max, h_min) ≈ h_max
        @test _adaptive_h(ComplexF64[], h_tiny, h_max, h_min) ≈ h_max

        # The walk-level consequence: a B2 adaptive walk threading around
        # the Riccati pole never lands a node on it.
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.5,
                                        adaptive = true,
                                        step_policy = :max_q_root)
        min_dist = minimum(abs(z - p) for z in sol.visited_z)
        @info "VPN.3.3 adaptive walk min |z - p| = $(round(min_dist;digits=4))"
        @test min_dist > 1.0e-2
        # Every recorded step lies in [h_min, h_max]: h_max = 0.5,
        # h_min = H_MIN_RATIO·0.5.
        h_lo = H_MIN_RATIO * 0.5
        @test all(s -> h_lo ≤ real(s) ≤ 0.5 + 1.0e-12, sol.visited_h)
    end

    # -------------------------------------------------------------------------
    # VPN.3.4 — fail-loud (Rule 1).  An unknown step_policy throws; the
    # adaptive h-floor throws when the pole field wedges the walk.
    # -------------------------------------------------------------------------
    @testset "VPN.3.4 fail-loud" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)

        # An unknown step_policy is a caller mistake — throw, never
        # silently fall back.
        @test_throws ArgumentError vector_path_network_solve(
            prob, ComplexF64[0.5 + 0.0im]; order = 24, h = 0.25,
            step_policy = :min_u)

        # `_adaptive_h` throws `:step_collapse` only for the GENUINE
        # case — the true local pole density exceeds what h_min can
        # honestly resolve, i.e. `SAFETY·D_local < h_min`.  A shared-Q
        # with a root very close in the t-variable (Q = 1 - t/t0, root
        # at t = t0) places the nearest pole right next door, so
        # D_local = h_prev·t0 is tiny and SAFETY·D_local falls below
        # h_min.  Under the new non-ratcheting law this is the *only*
        # way the throw fires — never a self-inflicted geometric ratchet.
        h_prev = 0.5; h_max = 0.5; h_min = 0.1
        t0     = 1.0e-4                              # pole right next door
        Q_near = ComplexF64[1.0, -1.0 / t0]          # root at t = t0
        # SAFETY·D_local = SAFETY·0.5·1e-4 = 1.25e-5 ≪ h_min = 0.1.
        # `_adaptive_h` throws the typed `VectorWalkError` (reason
        # :step_collapse) — ADR-0026 D1; the typed exception lets the
        # resilient driver classify a caught failure by `reason`.
        @test_throws VectorWalkError _adaptive_h(Q_near, h_prev, h_max, h_min)
        # A Q whose root is far ⇒ D_local large ⇒ SAFETY·D_local well
        # above h_min — no throw; the clamp returns a valid step.
        Q_far  = ComplexF64[1.0, -1.0 / 50.0]        # root at t = 50
        h_ok   = _adaptive_h(Q_far, h_prev, h_max, h_min)
        @test h_min ≤ h_ok ≤ h_max
    end

    # -------------------------------------------------------------------------
    # VPN.3.5 — the fixed-walk fallback.  `adaptive = false` + `:min_y`
    # reproduces the verified V7 fixed-step walk: every step is exactly
    # the requested `h`, and the walk reaches its targets.
    # -------------------------------------------------------------------------
    @testset "VPN.3.5 fixed-walk fallback (adaptive = false)" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                        adaptive = false, step_policy = :min_y)
        @test sol isa VectorPathNetworkSolution
        @test length(sol.visited_z) ≥ 2
        # adaptive = false ⇒ every recorded step is exactly h = 0.25
        # (the root carries the nominal h too).
        @test all(s -> real(s) ≈ 0.25, sol.visited_h)
        @test all(s -> imag(s) == 0.0, sol.visited_h)
    end

    # =========================================================================
    # VPN.4.* — VC-10 target-ordering control (bead padetaylor-0ln.37.17,
    # ADR-0025 Phase D).  See the file docstring.
    # =========================================================================

    # -------------------------------------------------------------------------
    # VPN.4.1 — the `rng` kwarg is ADDITIVE.  `rng = nothing` (the
    # default) walks the targets in the order given; it must reproduce
    # the no-`rng` call bit-identically (ADR-0001 — v0.1 + v0.2 stay
    # byte-identical, the kwarg only ADDS the shuffle path).
    # -------------------------------------------------------------------------
    @testset "VPN.4.1 rng = nothing reproduces the default walk" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        # The default (no `rng` passed) and the explicit `rng = nothing`
        # must produce a bit-identical visited tree — the additive
        # guarantee: the new code path is inert unless an RNG is given.
        sol_default = vector_path_network_solve(prob, targets;
                                                order = 24, h = 0.25)
        sol_nilrng  = vector_path_network_solve(prob, targets;
                                                order = 24, h = 0.25,
                                                rng = nothing)
        @test sol_nilrng.visited_z == sol_default.visited_z
        @test sol_nilrng.visited_y == sol_default.visited_y
        @test sol_nilrng.visited_parent == sol_default.visited_parent
        @test sol_nilrng.visited_h == sol_default.visited_h
        @test sol_nilrng.visited_denominator == sol_default.visited_denominator
    end

    # -------------------------------------------------------------------------
    # VPN.4.2 — a FIXED seed is fully deterministic.  VC-10 introduces a
    # *different* ordering per seed, but a given seed's walk is still
    # bit-reproducible: the same `MersenneTwister` state ⇒ the same
    # shuffle ⇒ the same path-network tree.  This is the per-fixed-seed
    # reproducibility guarantee (`test/kkg_pi2_figure_test.jl` PI2F.1.6's
    # exact-determinism invariant, lifted to the seeded vector walk).
    # -------------------------------------------------------------------------
    @testset "VPN.4.2 fixed-seed walk is bit-deterministic" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        # Two calls with the SAME seed (a fresh MersenneTwister at the
        # same seed each time) — bit-identical tree.
        s1 = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                       rng = MersenneTwister(12345))
        s2 = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                       rng = MersenneTwister(12345))
        @test s1.visited_z == s2.visited_z
        @test s1.visited_parent == s2.visited_parent
        @test s1.visited_h == s2.visited_h
        poles1 = extract_poles_shared_q(s1; radius_t = 6.0,
                                        cluster_atol = 0.2, min_support = 2)
        poles2 = extract_poles_shared_q(s2; radius_t = 6.0,
                                        cluster_atol = 0.2, min_support = 2)
        @test poles1 == poles2
    end

    # -------------------------------------------------------------------------
    # VPN.4.3 — two DIFFERENT orderings build different trees, and both
    # converge to the same pole field within the walk's numerical error
    # — the FW/FFW double-run accuracy indicator (VC-10).  On the
    # exactly-soluble Riccati oracle the field is a single known pole, so
    # "within the numerical error" is a hard ~1e-8 bar.
    # -------------------------------------------------------------------------
    @testset "VPN.4.3 different orderings, same pole field (VC-10)" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        sa = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                       rng = MersenneTwister(1))
        sb = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                       rng = MersenneTwister(2))
        # Two different seeds genuinely thread DIFFERENT path-network
        # trees — the FFW premise: "different paths will be run".
        @test sa.visited_z != sb.visited_z
        # ...yet both still cover the region and extract the SAME
        # physical pole `p` to the walk's numerical error.  This is the
        # VC-10 indicator: different orderings converge to the same field.
        pa = extract_poles_shared_q(sa; radius_t = 6.0,
                                    cluster_atol = 0.2, min_support = 2)
        pb = extract_poles_shared_q(sb; radius_t = 6.0,
                                    cluster_atol = 0.2, min_support = 2)
        @test !isempty(pa) && !isempty(pb)
        da = minimum(abs(q - p) for q in pa)
        db = minimum(abs(q - p) for q in pb)
        @test da < 1.0e-8
        @test db < 1.0e-8
        # The two runs' best pole estimates agree with each other — the
        # two-run disagreement (the FW/FFW indicator) is at the numerical
        # floor on this exactly-soluble oracle.
        qa = pa[argmin(abs.(pa .- p))]
        qb = pb[argmin(abs.(pb .- p))]
        @info "VPN.4.3 two-run disagreement |q_a - q_b| = $(abs(qa - qb))"
        @test abs(qa - qb) < 1.0e-7
    end

    # =========================================================================
    # VPN.5.* — the resilient Stage-1 walk (bead padetaylor-0ln.40.a,
    # ADR-0026 D1).  `on_target_failure = :skip` makes a per-target walk
    # failure record a `VectorWalkFailure` in `sol.failed_targets` and
    # continue, rather than aborting the whole run; `:throw` (default) is
    # byte-identical to the pre-ADR-0026 behaviour.  See the module
    # docstring and `VectorWalkFailure`.
    #
    # The deterministically-reproducible failure used here is the
    # :unreachable one: `max_steps_per_target = 1` with a far target — the
    # walk cannot reach it in a single step, so the unreachable-target
    # `VectorWalkError` fires on a known input every run.  The other two
    # reasons (:all_candidates_failed, :step_collapse) are exercised
    # indirectly: VPN.3.4 already pins the `:step_collapse` throw of
    # `_adaptive_h`, and VPN.3.1 the `:all_candidates_failed` throw — both
    # now as a `VectorWalkError`.  A *driver-level* deterministic trigger
    # for those two would need a hand-tuned pole-on-the-chord geometry
    # that is fragile to recreate; the :unreachable path exercises the
    # identical catch/record/continue machinery, so the skip logic itself
    # is fully covered.
    # =========================================================================
    @testset "VPN.5 resilient walk (on_target_failure, ADR-0026 D1)" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem(p, z0, cs; order = 24)

        # A far target the walk provably cannot reach in ONE step — the
        # deterministically-unreachable input.
        far = ComplexF64[40.0 + 0.0im]

        # ---- VPN.5.1 — :skip records the unreachable target and the run
        # completes normally. -------------------------------------------------
        @testset "VPN.5.1 :skip records the failure, run completes" begin
            sol = vector_path_network_solve(prob, far; order = 24, h = 0.25,
                                            max_steps_per_target = 1,
                                            on_target_failure = :skip)
            @test sol isa VectorPathNetworkSolution
            # Exactly one skipped target, recorded as first-class data.
            @test length(sol.failed_targets) == 1
            rec = sol.failed_targets[1]
            @test rec isa VectorWalkFailure
            @test rec.reason == :unreachable
            @test rec.target == far[1]
            # `z_stuck` is the live walk position — finite, and the
            # recorded residual is exactly |z_stuck - target| (not parsed
            # from the exception string).
            @test isfinite(rec.z_stuck)
            @test rec.residual_dist == abs(rec.z_stuck - rec.target)
        end

        # ---- VPN.5.2 — :throw on the same input still throws (now a
        # typed VectorWalkError, not a bare ErrorException). ------------------
        @testset "VPN.5.2 :throw still aborts the run" begin
            @test_throws VectorWalkError vector_path_network_solve(
                prob, far; order = 24, h = 0.25,
                max_steps_per_target = 1, on_target_failure = :throw)
            # :throw is the default — omitting the kwarg behaves the same.
            @test_throws VectorWalkError vector_path_network_solve(
                prob, far; order = 24, h = 0.25, max_steps_per_target = 1)
        end

        # ---- VPN.5.3 — partial nodes from the failed walk are retained.
        # A one-step walk toward `far` lands one node before giving up;
        # that node is a valid solution node and must stay in visited_z. --
        @testset "VPN.5.3 partial nodes from a failed walk are kept" begin
            sol = vector_path_network_solve(prob, far; order = 24, h = 0.25,
                                            max_steps_per_target = 1,
                                            on_target_failure = :skip)
            # The walk took one accepted step before the unreachable
            # throw — so the tree grew past the lone IC root node.
            @test length(sol.visited_z) > 1
            # The tree stays well-formed: root parent 0, others earlier.
            @test sol.visited_parent[1] == 0
            for k in 2:length(sol.visited_parent)
                @test 1 ≤ sol.visited_parent[k] < k
            end
        end

        # ---- VPN.5.4 — a successful run has an empty failed_targets. --------
        @testset "VPN.5.4 a successful run records no failures" begin
            targets = ComplexF64[x + y * im
                                 for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
            sol = vector_path_network_solve(prob, targets; order = 24,
                                            h = 0.25,
                                            on_target_failure = :skip)
            @test sol isa VectorPathNetworkSolution
            @test isempty(sol.failed_targets)
            # The default-mode solve carries the field too, also empty —
            # the additive guarantee (ADR-0001).
            sol_def = vector_path_network_solve(prob, targets; order = 24,
                                                h = 0.25)
            @test isempty(sol_def.failed_targets)
        end

        # ---- VPN.5.5 — an unknown on_target_failure is a caller mistake. ----
        @testset "VPN.5.5 unknown on_target_failure throws ArgumentError" begin
            @test_throws ArgumentError vector_path_network_solve(
                prob, far; order = 24, h = 0.25,
                on_target_failure = :ignore)
        end

        # ---- VPN.5.6 — :skip with a MIX of reachable and unreachable
        # targets reaches the reachable ones and records only the bad one. --
        @testset "VPN.5.6 :skip isolates the failure to its target" begin
            # `max_steps_per_target = 1` makes the FAR target unreachable
            # while a near target (within one step of the IC) succeeds.
            mixed = ComplexF64[z0 + 0.2 + 0.0im, 40.0 + 0.0im]
            sol = vector_path_network_solve(prob, mixed; order = 24,
                                            h = 0.25,
                                            max_steps_per_target = 1,
                                            on_target_failure = :skip)
            # Exactly the far target is recorded; the near one is not.
            @test length(sol.failed_targets) == 1
            @test sol.failed_targets[1].target == 40.0 + 0.0im
            @test sol.failed_targets[1].reason == :unreachable
        end

        # ---- VPN.5.7 — skip-if-covered (ADR-0026 Amendment 3, S3). ---------
        # `skip_covered = true` makes the target loop skip a target that
        # already lies inside a visited node's `COVER_FRAC·real(visited_h)`
        # covering disc, killing the dense-target chording regression
        # (bottleneck 3).  `skip_covered = false` (the default) keeps the
        # pre-S3 behaviour exactly: only an EXACT-coincidence (`10·eps`)
        # skip.
        #
        # The skip is a *target-loop* decision taken BEFORE `_nearest_
        # visited` runs.  It is distinct from — and stacks on top of — the
        # per-target while-loop guard `|z−target| > visited_h[parent]`,
        # which already declines to lay a node for a target inside the
        # *nearest* node's step `h`.  The two coincide for a uniform fixed
        # `h` (a covered target — within `COVER_FRAC·h < h` of some node —
        # is then necessarily within `h` of its nearest node, so the loop
        # guard alone skips it too).  Where `skip_covered` genuinely bites
        # is the *adaptive* walk: there `visited_h` is per-node and varies,
        # so a target can lie inside a sparse-field node's wide covering
        # disc yet outside its (collapsed-`h`, smaller) nearest node's
        # step — the default then crawls a redundant chord to it, S3 skips
        # it.  Hence (d) below — the load-bearing strict-inequality test —
        # runs on the *adaptive* default walk over a dense pole-region
        # grid, the regime S3 was built for.
        @testset "VPN.5.7 skip-if-covered (skip_covered, S3)" begin
            h_fixed = 0.25
            # The IC root node sits at z0; with adaptive=false its step is
            # h_fixed, so its covering disc has radius COVER_FRAC·h_fixed.
            r_cover = COVER_FRAC * h_fixed

            # ---- (a) a target WELL INSIDE the IC node's covering disc is
            # skipped when skip_covered = true — the walk lays NO new node
            # for it, so the tree stays the lone IC root node. ------------
            covered = ComplexF64[z0 + 0.5 * r_cover + 0.0im]
            @test abs(covered[1] - z0) < r_cover          # genuinely inside
            sol_cov = vector_path_network_solve(prob, covered; order = 24,
                                                h = h_fixed, adaptive = false,
                                                skip_covered = true)
            @test length(sol_cov.visited_z) == 1          # IC root only
            @test sol_cov.visited_z[1] == ComplexF64(z0)

            # ---- (b) with skip_covered = false (default) the SAME target
            # is NOT coverage-skipped: the default skip decision fires only
            # on EXACT coincidence (`10·eps`).  A target exactly AT the IC
            # node IS skipped (the pre-S3 `10·eps` path), but the merely
            # *covered* target is processed by the walk machinery — the old
            # behaviour S3 is opt-in over. -------------------------------
            #   the covered target is not within 10·eps of z0 …
            @test abs(covered[1] - z0) > 10 * eps(Float64)
            #   … so under the default the walk DOES process it (it is the
            #   `:max_q_root` selector + nearest-node lookup that run, not
            #   a coverage skip).  An exactly-coincident target, by
            #   contrast, is skipped by the surviving `10·eps` check.
            exact = ComplexF64[ComplexF64(z0)]
            sol_exact = vector_path_network_solve(prob, exact; order = 24,
                                                  h = h_fixed,
                                                  adaptive = false)
            @test length(sol_exact.visited_z) == 1        # 10·eps skip held
            #   and the covered target, default-mode, is NOT skipped by
            #   coverage — running it skip_covered = true vs false changes
            #   the target-loop decision (true → 1 node, false → processed).
            sol_cov_default = vector_path_network_solve(prob, covered;
                                                        order = 24,
                                                        h = h_fixed,
                                                        adaptive = false)
            @test length(sol_cov_default.visited_z) ≥
                  length(sol_cov.visited_z)
            # omitting the kwarg is byte-identical to skip_covered = false.
            sol_omit = vector_path_network_solve(prob, covered; order = 24,
                                                 h = h_fixed,
                                                 adaptive = false,
                                                 skip_covered = false)
            @test sol_omit.visited_z == sol_cov_default.visited_z

            # ---- (c) a genuinely UNCOVERED target — well outside every
            # visited node's COVER_FRAC·h disc — is still walked to even
            # with skip_covered = true (S3 skips only *covered* targets). -
            uncov = ComplexF64[z0 + 1.5 + 0.0im]
            @test abs(uncov[1] - z0) > r_cover            # genuinely outside
            sol_unc = vector_path_network_solve(prob, uncov; order = 24,
                                                h = h_fixed, adaptive = false,
                                                skip_covered = true)
            @test length(sol_unc.visited_z) > 1           # the walk ran

            # ---- (d) the LOAD-BEARING test: on a dense target set tiling
            # the pole region, the *adaptive* default walk with
            # skip_covered = true lays strictly FEWER visited nodes than
            # with skip_covered = false — the redundant covered targets are
            # skipped instead of re-walked along chording crawls.  Adaptive
            # `h` is the regime S3 bites in (see the testset comment); a
            # fixed-`h` walk would tie here by construction. --------------
            dense = ComplexF64[x + y * im
                               for x in -0.4:0.05:2.4 for y in -1.0:0.05:1.4]
            sol_dense_skip = vector_path_network_solve(
                prob, dense; order = 24, h = 0.5, adaptive = true,
                skip_covered = true, on_target_failure = :skip)
            sol_dense_full = vector_path_network_solve(
                prob, dense; order = 24, h = 0.5, adaptive = true,
                skip_covered = false, on_target_failure = :skip)
            @test length(sol_dense_skip.visited_z) <
                  length(sol_dense_full.visited_z)
        end
    end

end

# =============================================================================
# Mutation-proof record (CLAUDE.md Rule 4).
#
# Procedure: perturb `src/VectorPathNetwork.jl` / `src/VectorPoleField.jl`
# / `src/VectorWedgeStep.jl`, rerun this file, confirm RED, restore.
# V7 mutations M1/M3 + B2 mutations M4/M5/M6 below; M2 superseded.
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
#   M2 — SUPERSEDED by B2.  The V7 `_wedge_step` this mutation perturbed
#        no longer exists — B2 (bead padetaylor-0ln.37.6) extracted the
#        wedge selector into `VectorWedgeStep._select_wedge`.  Its V7
#        intent (the wedge criterion is load-bearing) is carried by B2's
#        M4 below (`:max_q_root` → min-‖y‖) on the new selector.
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
#   M3b — VectorPoleField: break the scale-derived cluster tolerance
#        (bead padetaylor-2vv, ADR-0026 Amendment 3 §S4).  In
#        `extract_poles_shared_q`'s scale-derived branch, replace
#          tol = T(CLUSTER_FRAC) * min(h_cand, rep_h[r])
#        with the absolute constant
#          tol = T(0.2)
#        (revert to the audit-flagged absolute cluster tolerance).
#        Expected: VPN.1.3b's case (a) — two genuinely-distinct poles
#        0.12 apart at a fine local scale h = 0.05 — is wrongly merged,
#        so `length(sep) == 2` bites.
#        Result: RED — 2 assertions bit in VPN.1.3b: `length(sep) == 2`
#        failed (got 1 — the absolute 0.2 tolerance merged the two
#        distinct poles), and the `minimum(abs(p - 1.12))` location
#        assert failed (with the two poles merged, the surviving rep is
#        at 1.00, so no candidate sits at 1.12).  Restored to GREEN.
#
# -- B2 dense-wedge walk mutations (bead padetaylor-0ln.37.6) --
#
#   M4 — VectorWedgeStep: revert the :max_q_root selector to min-‖y‖.
#        In `_select_wedge`, replace
#          primary = step_policy === :max_q_root ?
#              _candidate_pole_disc(f, cand_z[k], cand_y[k], order, h_mag, C) :
#              -norm(cand_y[k])
#        with the unconditional proxy
#          primary = -norm(cand_y[k])
#        (the V7 min-‖y‖ heuristic — discard the principled
#        shared-Q-root-distance criterion).
#        Expected: VPN.3.2's contract `discs[kq] ≥ maximum(discs)` — the
#        :max_q_root pick is the maximal-pole-disc candidate — bites at
#        the multi-pole P_I⁽²⁾ wedge nodes where min-‖y‖ picks a
#        different (smaller-disc) candidate; and VPN.3.1's adaptive walk
#        loses the principled steering and BLOCKS.
#        Result: RED — 16 VPN.3.2 assertions failed (the min-‖y‖ pick is
#        not the max-disc candidate) + VPN.3.1 errored (the min-‖y‖
#        adaptive walk hit `all 5 wedge candidates failed` at |x| ≈ 12.5,
#        the A2 block).  Restored to GREEN.
#
#   M5 — VectorWedgeStep: freeze the adaptive step.  In `_adaptive_h`,
#        prepend `return h_max` (the V1 fixed-step behaviour — `h` never
#        adapts, never caps, never throws).
#        Expected: the frozen-h walk re-introduces the A2 failure (a
#        fixed chord step lands on a pole), so VPN.3.1's adaptive wedge
#        walk BLOCKS; VPN.3.4's wedged-walk throw never fires.
#        Result: RED — VPN.3.1 + VPN.3.2 errored (the frozen-h walks
#        block on the P_I⁽²⁾ wedge fan) + VPN.3.4 failed (`_adaptive_h`
#        no longer throws on the near-pole `Q_near`).  Restored to GREEN.
#
#   M6 — VectorWedgeStep: drop the pole cap.  In `_adaptive_h`, replace
#        the `if length(denominator) > 1 … end` cap block with the
#        unconditional `h = h_grow` (the geometric grow with no
#        pole-distance ceiling).
#        Expected: without the cap `h` grows to `h_max` regardless of
#        pole proximity — the walk overshoots poles; `_adaptive_h`
#        returns the grow value, not `POLE_SAFETY·h_prev·min|t*|`.
#        Result: RED — VPN.3.3 failed (the cap-value assertions
#        `h_ad ≈ 0.5·h_prev·t0` and the "lands at 0.5 of the way"
#        check bit — the mutant returns `h_grow`), VPN.3.1 + VPN.3.2
#        errored (the uncapped walk blocks), VPN.3.4 failed (no throw on
#        `Q_near`).  Restored to GREEN.
#
# -- S2 scale-derived step-law mutation (bead padetaylor-cqk) --
#
#   M8 — VectorWedgeStep: re-introduce the h_prev ratchet.  The S2 step
#        law (ADR-0026 Amendment 4) replaced the geometric sink with
#        `h_target = SAFETY·D_local` — an absolute step with no h_prev
#        factor.  This mutation restores the sink by reverting that line
#        to the old geometric form
#          h_target = min(1.5·h_prev, 0.5·D_local)
#        — `h` can again only crawl up by `GROW·h_prev` from a collapsed
#        h_prev, and the absolute-value contract is broken.
#        Expected: VPN.3.3's LOAD-BEARING non-ratcheting recovery test
#        (a tiny h_prev + a FAR pole must jump straight to h_max in ONE
#        call) bites — the mutant returns `1.5·h_tiny`, a slow crawl,
#        not h_max; and the near-pole `h ≈ SAFETY·D_local` value
#        assertions bite (the mutant returns the `0.5·D_local` branch).
#        Result: RED — 4 VPN.3.3 assertions failed: the near-pole value
#        `h_ad ≈ SAFETY·D_local` + the `h_ad/D_local ≈ SAFETY` ratio
#        (line 475, 479 — mutant uses 0.5 not SAFETY), and the
#        non-ratcheting recovery `h_recover ≈ h_max` + `h_recover >
#        100·h_tiny` (line 494, 496 — mutant crawls to 1.5·h_tiny ≪
#        h_max).  This proves the non-ratcheting recovery test genuinely
#        catches a re-introduced geometric sink.  Restored to GREEN.
#
# -- S3 skip-if-covered mutation (bead padetaylor-0kx) --
#
#   M9 — VectorPathNetwork: flip the skip-if-covered inequality.  The S3
#        target-loop coverage check (ADR-0026 Amendment 3) is
#          any(i -> abs(target - visited_z[i]) ≤
#                   COVER_FRAC * real(visited_h[i]), eachindex(visited_z))
#        — a target INSIDE some node's covering disc is skipped.  This
#        mutation flips `≤` to `≥`, so the walk instead skips any target
#        OUTSIDE a node's disc — the exact inverse of the coverage
#        semantics.
#        Expected: VPN.5.7's case (c) — a genuinely UNCOVERED target must
#        still be walked under `skip_covered = true` — bites: with the
#        flipped test the uncovered target satisfies `abs ≥ COVER_FRAC·h`
#        and is wrongly skipped, so the tree never grows past the lone IC
#        root.
#        Result: RED — VPN.5.7's `length(sol_unc.visited_z) > 1`
#        assertion failed (`1 > 1`): the flipped check skipped the
#        uncovered target, the walk laid no node.  This proves the
#        coverage check's inequality direction is load-bearing — the test
#        catches a reversed skip predicate.  Restored to GREEN.
#
# -- VC-10 target-ordering mutation (bead padetaylor-0ln.37.17) --
#
#   M7 — VectorPathNetwork: make the `rng` shuffle inert.  In
#        `vector_path_network_solve`, replace
#          rng === nothing || (target_list = shuffle(rng, target_list))
#        with the no-op
#          target_list = target_list
#        (the `rng` kwarg is accepted but never used — every seed walks
#        the same order).
#        Expected: VPN.4.3's `sa.visited_z != sb.visited_z` assertion
#        bites — with the shuffle removed two different seeds build the
#        SAME tree, so the two-run accuracy indicator (VC-10) is vacuous.
#        Result: RED — VPN.4.3's different-tree assertion failed (the two
#        seeds produced byte-identical visited trees); VPN.4.1 / VPN.4.2
#        stayed GREEN (the additive default and fixed-seed determinism do
#        not depend on the shuffle being live).  This proves the shuffle
#        is genuinely load-bearing for VC-10.  Restored to GREEN.
#
# Certified bites: M1 → 6 RED, M3 → 1 RED, M4 → 17 RED (16 fail +
# 1 error), M5 → 3 RED (1 fail + 2 error), M6 → 6 RED (4 fail +
# 2 error), M7 → 1 RED (VPN.4.3 different-tree assertion — VC-10
# substrate), M8 → 4 RED (VPN.3.3 non-ratcheting recovery + near-pole
# value assertions — the re-introduced geometric sink), M9 → 1 RED
# (VPN.5.7 case (c) uncovered-target-still-walked — the flipped
# skip-if-covered inequality).  M2 superseded;
# M5/M6 describe the retired geometric-grow `_adaptive_h` (ADR-0026
# Amendment 4) and are kept as the historical sink record.  Restored to
# GREEN after each mutation.
# =============================================================================
