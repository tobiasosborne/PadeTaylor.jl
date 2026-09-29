# test/windowed_pole_ownership_test.jl
#
# ============================================================================
# WINDOWED POLE OWNERSHIP: the composite keeps every physical pole EXACTLY
# ONCE (bead padetaylor-580u, worklog 084).  Guards `WindowedTiling._own_poles`
# as called by `WindowedComposite.windowed_extract_poles`.
#
# THE BUG.  The pre-580u rule kept window wi's pole estimate p iff
# `_nearest_center(p) == wi` — decided per ESTIMATE.  Two overlapping windows
# estimate one physical pole at two positions ~1e-8 apart, so a pole lying ON
# a core line had its two estimates land on opposite sides of it and was kept
# TWICE (each in its own core) or ZERO times (each in the other's core).  On
# a real-symmetric problem every real-axis pole sits exactly on the
# `Im z = 0` core line of an even tiling — this fixture's case.
#
# THE ORACLE.  Weierstrass ℘ on the equianharmonic lattice, u'' = 6u² from the
# same ICs as external/probes/fig47-seam-diagnosis/p6_ingn_pole_truth.jl: the
# true poles are p(m,n) = 1 + 2ω m + 2ω e^{iπ/3} n EXACTLY, so "kept once" is
# checked against ground truth, not against another run.  [-20,20]² with
# window_extent 20 ⇒ a 2×2 tiling whose core lines are Re z = 0 and Im z = 0.
#
# WHY THE GATE IS DUPLICATE-SENSITIVE (unlike FSEAM.2's match_frac, which is
# blind to a pole duplicated in both runs): WPO.2 asserts per-exact-pole
# MULTIPLICITY == 1 — zero drops AND zero duplicates — plus zero cross-window
# kept pairs.  WPO.1 proves the fixture genuinely exercises the hazard (the
# legacy per-estimate rule, re-derived inline, DOES duplicate and drop here),
# so a passing WPO.2 certifies the fix, not an easy fixture.
#
# MEASURED (seed 0, order 30; worklog 084): legacy rule 4 dropped + 3 exact
# poles duplicated (4 cross-window kept pairs) of 246; fixed rule 0 / 0.
# Mutation proof in worklog 084 §4.
# ============================================================================

using PadeTaylor, Test, Printf
using PadeTaylor.WindowedTiling: _nearest_center, _own_poles

@testset "Windowed pole ownership: each physical pole kept once (WPO)" begin
    f(z, u, up) = 6u^2
    u0, up0 = 1.071822516416917, 1.710337353176786        # ℘ ICs (p6 probe)
    HALF = 20.0
    xs = ys = range(-HALF, HALF; length = 41)
    prob = PadeTaylorProblem(f, (u0, up0), (0.0, HALF); order = 30)
    wsol = windowed_path_network_solve(prob, xs, ys; window_extent = 20.0,
               overlap = 6.0, h = 0.5, order = 30, rng_seed = 0)
    C = wsol.centers

    ω = 1.3629079683730673; g1 = 2ω + 0im; g2 = 2ω * cis(pi / 3)
    lat   = ComplexF64[1 + m * g1 + n * g2 for m in -30:30 for n in -30:30]
    exact = [p for p in lat if abs(real(p)) ≤ HALF && abs(imag(p)) ≤ HALF]
    tol = 0.3                                   # ≪ lattice spacing 2ω ≈ 2.73
    mult(ps) = [count(p -> abs(p - e) < tol, ps) for e in exact]
    xpairs(ps, w) = count(((i, j),) -> w[i] != w[j] && abs(ps[i] - ps[j]) ≤ 0.5,
                          [(i, j) for i in eachindex(ps) for j in i+1:length(ps)])

    per_win = [extract_poles(ws) for ws in wsol.window_sols]

    # ------------------------------------------------------------------------
    # WPO.1  NON-TRIVIALITY — the legacy per-estimate rule fails HERE.
    # ------------------------------------------------------------------------
    @testset "WPO.1: legacy per-estimate ownership duplicates and drops" begin
        @test length(C) == 4                              # genuine 2×2 tiling
        legacy = ComplexF64[]; lw = Int[]
        for wi in eachindex(per_win), p in per_win[wi]
            _nearest_center(p, C) == wi && (push!(legacy, p); push!(lw, wi))
        end
        m = mult(legacy)
        @printf("[WPO.1] legacy: kept=%d exact=%d drops=%d dup=%d xpairs=%d\n",
                length(legacy), length(exact), count(==(0), m), count(≥(2), m),
                xpairs(legacy, lw))
        @test count(==(0), m) ≥ 1          # legacy DROPS a real pole ...
        @test count(≥(2), m) ≥ 1           # ... and DUPLICATES one
    end

    # ------------------------------------------------------------------------
    # WPO.2  THE FIX — every exact pole kept exactly once, nothing spurious.
    # ------------------------------------------------------------------------
    @testset "WPO.2: composite multiplicity is exactly 1 per exact pole" begin
        kept = windowed_extract_poles(wsol)
        # Recover each kept pole's window to count cross-window pairs.
        kw = [findfirst(wi -> any(==(p), per_win[wi]), eachindex(per_win)) for p in kept]
        m = mult(kept)
        spur = count(p -> minimum(abs(p - e) for e in lat) ≥ tol, kept)
        @printf("[WPO.2] fixed:  kept=%d drops=%d dup=%d xpairs=%d spurious=%d\n",
                length(kept), count(==(0), m), count(≥(2), m), xpairs(kept, kw), spur)
        @test all(==(1), m)                # no drop, no duplicate
        @test xpairs(kept, kw) == 0        # no cross-window pair within 0.5
        @test spur == 0
        # Same answer as calling the ownership kernel directly (the driver
        # adds nothing but the default boundary_atol = cluster_atol = 0.1).
        @test kept == _own_poles(per_win, C, 0.1)
    end

    # ------------------------------------------------------------------------
    # WPO.3  KWARG ROUTING — `merge_atol` now reaches PoleField.extract_poles
    #        (it used to be captured as a composite-level dedup), and a
    #        non-positive boundary_atol fails loud.
    # ------------------------------------------------------------------------
    @testset "WPO.3: merge_atol forwards; boundary_atol validated" begin
        fwd = [extract_poles(ws; merge_atol = 3.0) for ws in wsol.window_sols]
        @test fwd != per_win     # merge_atol = 3 (> spacing 2ω) genuinely changes extraction
        @test windowed_extract_poles(wsol; merge_atol = 3.0) == _own_poles(fwd, C, 0.1)
        @test_throws ArgumentError windowed_extract_poles(wsol; boundary_atol = 0)
    end
end
