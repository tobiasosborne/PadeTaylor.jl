# test/polefield_test.jl -- bead `padetaylor-xvf` tests.
#
# PoleField.extract_poles: recover pole locations in the z-plane from a
# solved path-network's per-node Padé store.  This is the capability
# behind FW 2011 Fig 4.7 / 4.8 (pole-location scatter plots).
#
# Ground-truth oracle — the equianharmonic Weierstrass-℘ test problem
# of FW 2011 §5.1.1 (references/markdown/FW2011_painleve_methodology_
# JCP230/FW2011_painleve_methodology_JCP230.md:281-318):
#
#     u''(z) = 6 u(z)^2,        u(z) = ℘(z + c₁; 0, c₂),   c₁ = -1, c₂ = 2.
#
# Its solution has *analytically known* second-order poles on a rhombic
# lattice of equilateral triangles.  FW md:297 gives the real-axis
# half-period ω = Γ(1/3)³ / (2^{13/6} π) ≈ 1.363; with c₁ = -1 the
# poles of u(z) sit at
#
#     z(m, n) = 1 + 2ω(m + n/2) + i · ω√3 · n,    m, n ∈ ℤ.
#
# So the extracted poles can be checked against an exact closed form —
# no pole-finder oracle is needed (this matches `figures/fw2011_fig_5_1.jl`,
# which makes the same observation).
#
# Tests: PF.1.1 (single-node extraction surfaces the nearest pole),
# PF.1.2 (full pole field over a 2D grid — no spurious poles, all
# in-region lattice poles recovered, conjugate symmetry), PF.2.1
# (degenerate-network edge case), PF.2.2 (the cross-node-support filter
# is load-bearing), PF.3.* (bead `padetaylor-26r` — extraction from a
# single-trajectory `solve_pade` result, `PadeTaylorSolution`).
# Mutation-proof procedure is in the file footer.

using Test
using PadeTaylor
using PadeTaylor.PathNetwork: PathNetworkSolution

include(joinpath(@__DIR__, "_oracle_problems.jl"))

# Equianharmonic ℘ real half-period, FW md:297.  Hard-coded to Float64
# precision rather than pulling in SpecialFunctions as a test dep;
# reproduce with `gamma(1/3)^3 / (2^(13/6) * π)`.
const Ω_FW = 1.3630340904278908

# Exact pole lattice of u(z) = ℘(z - 1; 0, 2): z(m,n) per FW md:297.
℘_pole(m::Int, n::Int) = 1 + 2Ω_FW * (m + n / 2) + im * (Ω_FW * sqrt(3) * n)

# All lattice poles with |m| ≤ M, |n| ≤ N — the reference set the
# spurious-pole check measures against.
℘_lattice(M::Int, N::Int) = [℘_pole(m, n) for m in -M:M for n in -N:N]

# Distance from `z` to the nearest exact lattice pole.
nearest_lattice_dist(z, lattice) = minimum(abs(z - p) for p in lattice)

# Is `z` inside the closed box [xlo,xhi] × [ylo,yhi]?  `extract_poles`
# honestly reports *every* pole ≥ min_support nodes agree on, including
# genuine lattice poles just outside the explored grid that only distant
# nodes can see — those are correctly identified as real, but placed
# less accurately (FW's Fig 4.7 likewise clips to a display window).
# Accuracy assertions are therefore made over the grid's covered box.
in_box(z, xlo, xhi, ylo, yhi) =
    xlo ≤ real(z) ≤ xhi && ylo ≤ imag(z) ≤ yhi

@testset "PoleField (bead padetaylor-xvf): FW 2011 Fig 4.7/4.8 capability" begin

    # Test ODE: u'' = 6u^2, FW 2011 §5.1.1 ICs at z = 0.
    fW(z, u, up) = 6 * u^2

    @testset "PF.1.1: single-node network surfaces the nearest pole" begin
        # A grid wholly inside |z| ≤ h = 0.5: Stage 1 takes no steps,
        # so the only visited node is the IC.  With the cross-node
        # support filter disabled (`min_support = 1`) the IC node's
        # own Padé is read directly.  Its (15,15) denominator places
        # the nearest pole — z = 1, at rescaled t = 1/h = 2 — to the
        # double-pole Padé accuracy floor (~5e-7 here), within the
        # figure-catalogue Float64 spec of 1e-6.  (A single far node
        # cannot place *distant* poles well — that is exactly what the
        # multi-node support filter exercised in PF.1.2 is for.)
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[0.0, 0.2 + 0.1im, -0.15 + 0.2im, 0.25 - 0.2im]
        sol  = path_network_solve(prob, grid; h = 0.5)

        @test length(sol.visited_z) == 1          # IC only — no Stage-1 steps

        poles = extract_poles(sol; min_support = 1)
        @test !isempty(poles)

        err_z1 = minimum(abs(p - (1.0 + 0im)) for p in poles)
        @test err_z1 ≤ 1.0e-6
    end

    @testset "PF.1.2: full pole field over a 2D grid" begin
        # A 2D grid covering several lattice poles.  Offsets chosen so
        # no grid target sits exactly on a pole.  Every in-region pole
        # is then seen by many independent visited nodes, so the
        # default cross-node support filter keeps the physical poles
        # and rejects node-local artefacts.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        xs = -2.4:0.5:4.1
        ys = -2.6:0.5:2.6
        grid = ComplexF64[x + im * y for x in xs for y in ys]
        sol  = path_network_solve(prob, grid; h = 0.5)

        poles = extract_poles(sol)
        @test !isempty(poles)

        # (a) No spurious poles: every extracted pole *inside the covered
        #     grid box* lies within 1e-6 of an exact lattice pole.
        #     (Best-placed estimates are ~1e-8 in practice; the 1e-6 pin
        #     is the figure-catalogue Float64 spec.)
        lattice = ℘_lattice(6, 3)
        in_grid = filter(p -> in_box(p, -2.4, 4.1, -2.6, 2.6), poles)
        @test !isempty(in_grid)
        worst = maximum(nearest_lattice_dist(p, lattice) for p in in_grid)
        @test worst ≤ 1.0e-6

        # (b) Completeness: the poles unambiguously inside the covered
        #     region must all be recovered, each to ≤ 1e-6.
        expected = ComplexF64[
            ℘_pole(0, 0),                       #  1            (real axis)
            ℘_pole(1, 0),                       #  3.726        (real axis)
            ℘_pole(-1, 0),                      # -1.726        (real axis)
            ℘_pole(0, 1),  ℘_pole(0, -1),       #  2.363 ± 2.36i
            ℘_pole(-1, 1), ℘_pole(-1, -1),      # -0.363 ± 2.36i
        ]
        for ze in expected
            err = minimum(abs(p - ze) for p in poles)
            @test err ≤ 1.0e-6
        end

        # (c) Conjugate symmetry: the IC is real, so the pole field is
        #     symmetric across the real axis.  Every off-axis extracted
        #     pole inside the covered box has a conjugate partner.
        for p in in_grid
            if abs(imag(p)) > 1.0e-3
                partner = minimum(abs(conj(p) - q) for q in poles)
                @test partner ≤ 1.0e-6
            end
        end
    end

    @testset "PF.2.1: degenerate network (constant denominators) → no poles" begin
        # A hand-built network whose every stored Padé has a constant
        # denominator b = [1] (a pure polynomial — no poles).  extract_poles
        # must skip such nodes cleanly and return an empty vector, not
        # throw on the empty `roots` call.
        T = Float64
        flat = PadeApproximant{Complex{T}}(
            Complex{T}[1.0, 0.5], Complex{T}[1.0], 1, 0, one(Complex{T}))
        sol = PathNetworkSolution{T}(
            Complex{T}[0.0, 0.5],            # visited_z
            Complex{T}[1.0, 1.2],            # visited_u
            Complex{T}[0.0, 0.1],            # visited_up
            [flat, flat],                    # visited_pade
            T[0.5, 0.5],                     # visited_h
            Int[0, 1],                       # visited_parent
            Complex{T}[],                    # grid_z
            Complex{T}[],                    # grid_u
            Complex{T}[],                    # grid_up
        )
        @test extract_poles(sol) == Complex{T}[]
        @test extract_poles(sol; min_support = 1) == Complex{T}[]
    end

    @testset "PF.2.2: cross-node support filter is load-bearing" begin
        # On the same 2D-grid network, disabling the support filter
        # (`min_support = 1`) lets node-local artefacts through, so it
        # must yield strictly more "poles" than the default — and those
        # extras must be the off-lattice spurious ones.  This pins the
        # filter as doing real work, not decoration.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        xs = -2.4:0.5:4.1
        ys = -2.6:0.5:2.6
        grid = ComplexF64[x + im * y for x in xs for y in ys]
        sol  = path_network_solve(prob, grid; h = 0.5)

        clean    = extract_poles(sol)                       # default filter
        unfilt   = extract_poles(sol; min_support = 1)      # filter disabled
        @test length(unfilt) > length(clean)

        lattice = ℘_lattice(6, 3)
        # The default result is on-lattice within the covered box; the
        # unfiltered result drags in node-local off-lattice artefacts.
        clean_in = filter(p -> in_box(p, -2.4, 4.1, -2.6, 2.6), clean)
        @test maximum(nearest_lattice_dist(p, lattice) for p in clean_in) ≤ 1.0e-6
        @test maximum(nearest_lattice_dist(p, lattice) for p in unfilt)   > 1.0e-6
    end

    @testset "PF.3.1: extraction from a single solve_pade trajectory" begin
        # The Phase-6 pole-bridge problem (Problems.jl module docstring):
        # u'' = 6u² with the FW ICs has a lattice pole at z = 1.  A single
        # solve_pade segment with h_max = 1.5 brackets it (rescaled
        # t = 1/1.5 ≈ 0.667, well inside the unit interval), so the
        # segment's stored Padé denominator places the pole — and
        # extract_poles must read it back off a `PadeTaylorSolution`
        # exactly as it does off a path-network's per-node store.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol  = solve_pade(prob; h = 1.5)
        @test sol isa PadeTaylorSolution
        @test length(sol.h) == 1                  # single segment

        poles = extract_poles(sol)                # default min_support = 1
        @test !isempty(poles)
        @test eltype(poles) <: Complex
        err_z1 = minimum(abs(p - (1.0 + 0im)) for p in poles)
        @test err_z1 ≤ 1.0e-6
    end

    @testset "PF.3.2: PadeTaylorSolution default min_support = 1 is load-bearing" begin
        # A single solve_pade trajectory is a chain, not a fan: the pole
        # at z = 1 is seen by exactly one segment here.  The default
        # `min_support = 1` for the `PadeTaylorSolution` method keeps it;
        # the path-network default of 3 would (correctly, for that
        # method) demand three independent sightings and discard it.
        # This pins the changed default as doing real work.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol  = solve_pade(prob; h = 1.5)

        @test !isempty(extract_poles(sol))                    # default = 1
        @test isempty(extract_poles(sol; min_support = 3))    # one segment < 3
    end

    @testset "PF.4.1: varying-h walk does not drop poles (scale-covariance, ADR-0026 §S7)" begin
        # The scalar twin of the vector §S7 fix (`VectorPoleField`,
        # ADR-0026 Amendment 6/7).  A trajectory whose per-segment step
        # `h` VARIES across nodes must not silently lose a genuine pole.
        #
        # Repro: stitch two `solve_pade` runs at different fixed `h` into
        # one `PadeTaylorSolution` (the legitimate adaptive/non-uniform-h
        # case).  A COARSE run over [0, 0.5] at h = 0.5 sets the walk's
        # step ceiling h_max = 0.5; a FINE run over [0.5, 0.85] at
        # h = 0.05 approaches the equianharmonic-℘ lattice pole at
        # z = ℘_pole(0,0) = 1 with small steps.  The five fine nodes at
        # z_ctr ∈ {0.5, 0.55, 0.6, 0.65, 0.7} each see the pole at z=1 at
        # z-distance {0.5, 0.45, 0.4, 0.35, 0.3} — all GREATER than
        # 5·h_fine = 5·0.05 = 0.25, so the OLD per-|t*| filter
        # (`|t*| ≤ radius_t`) discards every one of them; only the coarse
        # z=0 node and the two fine nodes at z_ctr ∈ {0.75, 0.8}
        # (z-distance ≤ 0.25) survive it — exactly THREE distinct nodes.
        # The NEW z-plane-distance filter (`|h·t*| ≤ radius_t·h_max`,
        # h_max = 0.5 ⇒ radius_z = 2.5) keeps all eight nodes' sightings.
        #
        # Measured per-node accuracy of the z=1 estimate (see probe):
        #   z_ctr=0.0  (h=0.5,  |t*|=2.0, z-dist=1.0): err 4.57e-7  (worst)
        #   z_ctr=0.8  (h=0.05, |t*|=4.0, z-dist=0.2): err 2.94e-9  (best)
        # The OLD greedy clustering orders by |t*|, so the min-|t*| node
        # (the FAR coarse z=0 node, |t*|=2.0) becomes the representative —
        # its placement is the *worst* of the surviving set (4.57e-7).
        # The NEW clustering orders by z-plane distance, so the z-closest
        # node (z_ctr=0.8, z-dist=0.2) wins — placement 2.94e-9.
        prob_c = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 0.5); order = 30)
        sol_c  = solve_pade(prob_c; h = 0.5)
        u05, up05 = sol_c(0.5)
        prob_f = PadeTaylorProblem(fW, (u05, up05), (0.5, 0.85); order = 30)
        sol_f  = solve_pade(prob_f; h = 0.05)

        Tt = Float64
        Yt = Tuple{Float64, Float64}
        Pt = eltype(sol_c.pade)
        z_all    = vcat(sol_c.z,    sol_f.z[2:end])     # drop duplicated z=0.5
        y_all    = vcat(sol_c.y,    sol_f.y[2:end])
        h_all    = vcat(sol_c.h,    sol_f.h)
        pade_all = vcat(sol_c.pade, sol_f.pade)
        sol = PadeTaylorSolution{Tt, Yt, Pt}(z_all, y_all, h_all, pade_all)

        # Precondition the repro actually has varying h with a high ceiling.
        @test maximum(h_all) ≈ 0.5
        @test minimum(h_all) ≈ 0.05

        # `min_support = 4` is the discriminator: the OLD per-|t*| filter
        # leaves the z=1 cluster with only 3 surviving nodes (< 4), so it
        # drops the pole entirely.  The NEW z-distance filter keeps all 8,
        # so the cluster clears the bar.  RED-before / GREEN-after.
        pole_z1 = ℘_pole(0, 0)                       # exact lattice pole: z = 1
        poles   = extract_poles(sol; min_support = 4)

        near = filter(p -> abs(p - pole_z1) < 0.2, poles)
        # (a) Filter half — the genuine pole is recovered at all.  RED on
        #     the OLD per-|t*| filter (cluster support 3 < 4).
        @test length(near) == 1

        # (b) Sort-key half — the representative is the z-plane-CLOSEST
        #     node (err 2.94e-9), not the min-|t*| far node (err 4.57e-7).
        #     1e-8 sits cleanly between the two: RED on the OLD min-|t*|
        #     ordering, GREEN on the NEW z-distance ordering.
        @test abs(near[1] - pole_z1) ≤ 1.0e-8

        # (c) No spurious far poles *inside the covered region*.  The
        #     trajectory runs along the real axis z ∈ [0, 0.85]; the only
        #     lattice pole it can place accurately is z = 1 (just past the
        #     end).  As the file header and PF.1.2/2.2 document, far nodes
        #     also "see" genuine off-region poles but place them poorly
        #     (here `-0.22 ± 2.06i`-class candidates), so accuracy is
        #     asserted only over the covered box — exactly the `in_box`
        #     discipline of PF.1.2.  Widening the far-root window must not
        #     invent a spurious pole IN-region: every in-box pole is on
        #     the exact ℘ lattice to the Float64 catalogue spec.
        lattice = ℘_lattice(2, 2)
        in_cov  = filter(p -> in_box(p, -0.5, 1.2, -0.5, 0.5), poles)
        @test !isempty(in_cov)
        @test maximum(nearest_lattice_dist(p, lattice) for p in in_cov) ≤ 1.0e-6
    end

    # ----------------------------------------------------------------------
    # PF.5.1 — dense 2D ELLIPTIC-LATTICE pole field: no over-split duplicates
    # (bug padetaylor-fzse).  THE regression the self-merge was built for.
    #
    # On the dense equianharmonic-℘ field of figures/demo_lattice_singularities,
    # one physical lattice pole is seen by MANY nodes at slightly different t*,
    # and those z-estimates spread by MORE than cluster_atol (0.1) — so the
    # greedy first-fit fragments each into 2-4 duplicate reps (each still
    # clearing min_support).  The post-pass single-linkage self-merge collapses
    # the duplicates; this test pins it by comparing the DEFAULT (self-merge)
    # against the LEGACY clustering (`merge_atol = 0`) on the SAME solve.
    #
    # CALIBRATED on external/probes/fzse-calibration/probe.jl (vs the closed-form
    # lattice, in + out of window): legacy → 32 reps with 4 over-split DUP-poles;
    # the DISJOINT-support self-merge → 24 reps, **2** dup-poles, PRECISION 1.0
    # (every rep within 0.5·spacing of a real lattice point — there are NO
    # interstitial false positives at the default; the bead's "FPs" were genuine
    # OUT-OF-WINDOW poles + an artifact of tuning cluster_atol to 0.4).
    #
    # The 2 residual dups are NOT cross-node over-split — they are the ℘ DOUBLE
    # pole resolved by a SINGLE node into two roots ~0.15-0.2 apart (shared
    # support), which the disjoint-support predicate deliberately does NOT merge
    # (it is indistinguishable BY DISTANCE from a genuine near-coalescent PAIR,
    # which CPN.4 / CRic.3 require kept distinct).  Collapsing a double pole to
    # one location-with-multiplicity-2 is the job of the multiplicity API
    # (bead padetaylor-90oh); here the load-bearing claim is that the CROSS-NODE
    # froth is removed (4 → 2 dup-poles) with no spurious reps.
    #
    # Recall is bounded by the WALK's node coverage (bottom-edge rows are sparsely
    # visited), NOT by the clustering, so it is unchanged by the self-merge — and
    # deliberately NOT asserted = 17 (that would conflate extraction with walk
    # density).
    # ----------------------------------------------------------------------
    @testset "PF.5.1: dense ℘-lattice field — self-merge kills cross-node over-split (fzse)" begin
        SP = 2 * Ω_FW                                          # lattice spacing 2.726
        full_lat = ℘_lattice(8, 8)                             # in + out of window
        win_lat  = filter(p -> in_box(p, -4.0, 6.0, -5.0, 5.0), full_lat)
        ndup(reps) = count(t -> count(r -> abs(r - t) ≤ 0.5SP, reps) > 1, win_lat)
        ncov(reps) = count(t -> any(r -> abs(r - t) ≤ 0.5SP, reps), win_lat)

        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        xs = range(-4.0, 6.0; length = 121); ys = range(-5.0, 5.0; length = 121)
        grid = ComplexF64[complex(x, y) for y in ys for x in xs]
        sol  = path_network_solve(prob, grid; h = 0.5,
                                  max_steps_per_target = 4000,
                                  enforce_real_axis_symmetry = true)

        legacy = extract_poles(sol; merge_atol = 0)            # greedy first-fit only
        merged = extract_poles(sol)                            # default self-merge

        # (a) The bug EXISTS without the fix: the dense field over-splits.
        @test ndup(legacy) ≥ 3
        # (b) The fix REDUCES the over-split (the cross-node froth is collapsed).
        #     The residual (≤ the ℘ double-pole doublet-splits) is a multiplicity
        #     matter deferred to padetaylor-90oh; the load-bearing claim is that
        #     a DISJOINT-support physical-pole duplicate no longer survives.
        @test ndup(merged) < ndup(legacy)
        @test ndup(merged) ≤ 2
        # (c) The self-merge only REMOVES reps (collapses duplicates).
        @test length(merged) < length(legacy)
        # (d) PRECISION = 1.0: every reported pole is a real lattice point
        #     (in OR out of window) — no interstitial spurious rep.  0.5·spacing
        #     is generous, but the over-split froth this test targets sits ON the
        #     lattice; an interstitial FP would be ≥ ~0.9 from any pole.
        @test all(r -> nearest_lattice_dist(r, full_lat) ≤ 0.5SP, merged)
        @test maximum(nearest_lattice_dist(r, full_lat) for r in merged) < 0.5
        # (e) Recall is NOT reduced by the merge (coverage is a walk property).
        @test ncov(merged) == ncov(legacy)
        @test ncov(merged) ≥ 12                                # walk-coverage floor
    end
end

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Each mutation below was applied to `src/PoleField.jl`, the suite
# confirmed RED, then the mutation reverted.
#
#   M1 — root → z mapping.  Change `z_v + h * t` to `z_v + t` (drop the
#        canonical-step rescale).  PF.1.1's `err_z1 ≤ 1e-6` goes RED
#        (z = 1 maps to t = 2), and PF.1.2 / PF.2.2 go RED — every pole
#        lands off the lattice.  (Observed: 3 failed, 2 errored.)
#
#   M2 — radius_t filter.  Comment out `abs(t) ≤ radius || continue`.
#        Far extrapolation roots leak in, corrupting the greedy-|t*|
#        clustering: PF.1.2 (completeness + symmetry) and PF.2.2 go RED.
#        (Observed: 43 failed.)
#
#   M3 — cross-node support filter.  Replace the final comprehension
#        guard `length(support[j]) ≥ min_support` with `≥ 1`.  PF.1.2
#        (a) goes RED (node-local artefacts reported as in-box poles)
#        and PF.2.2 goes RED (`length(unfilt) == length(clean)`, and the
#        default result is no longer on-lattice).  (Observed: 49 failed.)
#
#   M4 — PadeTaylorSolution method default (bead `padetaylor-26r`).
#        Change the `extract_poles(sol::PadeTaylorSolution; …)` default
#        `min_support = 1` to `min_support = 3`.  PF.3.1 goes RED
#        (`!isempty(poles)` fails / the `minimum` over an empty result
#        errors — the single segment cannot reach support 3) and PF.3.2's
#        first assertion goes RED.  Restored.
#
#   Shared-core re-verification.  M1–M3 act on `_extract_poles_core`,
#        which both `extract_poles` methods now call.  M1 was directly
#        re-applied: changing `z_ctr + h * t` → `z_ctr + t` takes down
#        PF.1.1 (path-network) *and* PF.3.1 (single trajectory),
#        confirming the shared core is exercised by both paths.  M2 / M3
#        live in the same shared core and were already proven against the
#        path-network path (PF.1.2 / PF.2.2).
#
#   M5 — §S7 z-plane-distance FILTER (bead padetaylor-bez, 2026-06-02).
#        In `_extract_poles_core`, revert the far-root filter from the
#        z-plane test `z_dist ≤ radius_z` back to the legacy per-node
#        `abs(t) ≤ radius`.  PF.4.1 (a) goes RED — `length(near) == 1`
#        evaluates `0 == 1`: under the per-|t*| filter only 3 distinct
#        nodes survive the z=1 cluster (the coarse z=0 node + the two
#        fine nodes at z-distance ≤ 0.25), below the `min_support = 4`
#        bar, so the genuine pole is dropped; (b)/(c) then fail on the
#        empty `near`/`in_cov`.  All of PF.1.*/2.*/3.* stay GREEN —
#        for the uniform-`h` walks they exercise, `h ≈ h_max`, so the
#        z-plane test reduces exactly to the legacy `|t*| ≤ radius_t`.
#        (Observed: 31 passed, 2 failed, 2 errored.)  Restored.
#
#   M6 — §S7 z-plane-distance SORT KEY (bead padetaylor-bez, 2026-06-02).
#        Keeping the new filter, revert ONLY the cluster sort key:
#        change the candidate tuple's 2nd field from `RT(z_dist)` back
#        to `RT(abs(t))`.  PF.4.1 (a) STILL PASSES (the new filter keeps
#        all 8 nodes' sightings, support 8 ≥ 4, so z=1 is recovered) —
#        but (b) `abs(near[1] - pole_z1) ≤ 1e-8` goes RED, evaluating
#        `4.56910677626432e-7 ≤ 1e-8`: ordering by `|t*|` crowns the
#        FAR coarse z=0 node (|t*|=2.0, the smallest |t*| of the
#        surviving set) the representative, and its placement (4.57e-7)
#        is two orders worse than the z-closest fine node's (2.94e-9)
#        that the z-distance ordering picks.  The representative drifts
#        exactly as documented.  (Observed: 34 passed, 1 failed.)
#        Restored.  The two halves are ONE atomic §S7 fix (ADR-0026
#        Amendment 7): M5 proves the filter half, M6 the sort-key half.
#
#   M7 — self-merge OFF (bug padetaylor-fzse, 2026-06-14).  In
#        `_extract_poles_core` change the default merge radius
#        `mtol = merge_atol === nothing ? h_max : …` to `… ? zero(h_max) :`
#        (no post-pass self-merge by default).  PF.5.1 (b) `ndup(merged) <
#        ndup(legacy)` and (c) `length(merged) < length(legacy)` go RED — the
#        default no longer collapses the cross-node froth (merged == legacy).
#        Restored.
#   M8 — DISJOINT-support condition dropped.  In `_single_linkage_merge`
#        delete `&& isdisjoint(supports[a], supports[b])` (distance-only
#        merge).  The near-COALESCENT-pair fixtures go RED — CPN.4
#        (corpus_pathnet_walls_rows, the δ-sweep coalescent pair) and CRic.3
#        (corpus_riccati_rational, the Hermite zero-pair) — because two
#        distinct poles a single node resolves (shared support, separation
#        < h_max) would now wrongly merge.  This proves the disjoint-support
#        predicate — not distance — is what protects coalescent pairs while
#        still collapsing cross-node froth.  Restored.
# ----------------------------------------------------------------------
