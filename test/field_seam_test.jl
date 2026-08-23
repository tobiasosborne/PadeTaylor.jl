# test/field_seam_test.jl
#
# ============================================================================
# FIELD SEAM: SEED-INVARIANCE GATE for the FW 2011 Fig 4.7 "seam" (bead
# padetaylor-sny7, P0 epic padetaylor-vwgl).  Guards the cure shipped in
# src/WindowedComposite.jl (`windowed_path_network_solve`).
#
# THE BUG (worklog 077 root cause, 078 measurement).  The monolithic
# `path_network_solve` walks ONE tree from the IC across the whole grid.  Its
# branch structure is `rng_seed`-dependent (Stage-1 target order is shuffled),
# and two seed-different subtrees accumulate DIFFERENT IVP error; where they
# meet, the reconstructed pole lattice has a grain boundary — a path-dependence
# seam.  FW never solved Fig 4.7 monolithically: they tiled it into 25
# independent ≤20×20 windows, each from the shared IC (FW2011 md:147).
# `windowed_path_network_solve` productionises that tiling; this file is the
# regression gate that the tiling actually removes the seed-dependence.
#
# WHY THIS GATE IS NON-GAMEABLE (the crux — worklog 078 §1 cross-cutting #4).
# A naive "two seeds agree" test is gameable: an impl that IGNORED the global
# seed would make both runs bit-identical and pass trivially while curing
# nothing.  So this file asserts BOTH halves —
#   FSEAM.1 (non-triviality): the two global seeds GENUINELY re-randomised the
#     walk — `window_seeds` differ AND at least one window's internal walk
#     (`visited_z`) differs.  This proves the pipeline is live, so a passing
#     FSEAM.2 certifies real path-independence, not a frozen pipeline.
#   FSEAM.2 (the cure): the EXTRACTED POLE FIELD is seed-invariant — bidirectional
#     pole-set match ≥ 90% and soft count agreement.
#
# MEASURED RED/GREEN (this exact config; ~18-38s for the two composite solves):
#   monolithic two-seed pole-set match @0.5 = 77.6%  (RED — fails 0.90 by 12pt)
#   windowed composite                       = 97.2%  (GREEN — passes by 7pt)
# Threshold 0.90 sits cleanly between the two regimes.  The complementary
# check that the composite does not DROP REAL poles (only spurious seed-
# dependent ones) is bead padetaylor-ingn (Weierstrass-℘ recall/precision),
# OUT OF SCOPE for this seed-invariance gate.
#
# REFERENCES.  FW2011 md:147 (the 5×5 / 18-run composite construction);
# docs/worklog/078-fig47-seam-cure-scoping.md §4-5 (P5 confirmation: |Δ|
# 151→12, match 77→99.4%, recommended plan); src/WindowedComposite.jl module
# docstring "Why per-window seeding is load-bearing".  Mutation-proof in the
# file footer.
# ============================================================================

using PadeTaylor, Test, Printf

# match_frac(A, B; tol) — directional pole-set agreement: the fraction of A's
# poles that lie within `tol` of SOME pole of B (a one-sided Hausdorff-style
# recall).  Asserting it in BOTH directions (A→B and B→A) catches both a pole
# that moved AND a pole that appeared/vanished between seeds.  Empty A returns
# 0.0 so the gate fails loud (Rule 1) rather than passing on a degenerate
# 0/0; the testset also pins the pole sets non-empty.
function match_frac(A::AbstractVector, B::AbstractVector; tol::Real = 0.5)
    isempty(A) && return 0.0
    return count(a -> any(b -> abs(a - b) ≤ tol, B), A) / length(A)
end

@testset "Field seam: path-independence via windowed composite (FSEAM)" begin

    # PI tritronquée; FW eq. (4.1) ICs (the Fig 4.7 problem).
    pI(z, u, up) = 6 * u^2 + z
    u0, up0 = -0.1875, 0.3049
    prob = PadeTaylorProblem(pI, (u0, up0), (0.0, 30.0); order = 20)

    # 61×61 lattice over [-30,30]² (spacing 1.0) — large enough that the
    # monolithic seam is visible, small enough for a fast gate.
    xs = ys = range(-30.0, 30.0; length = 61)

    # The two composite solves at two GLOBAL seeds.  Built ONCE here and
    # reused by both inner testsets (do NOT re-solve — exactly two windowed
    # solves total keeps the gate to ~18-38s).
    wsol0  = windowed_path_network_solve(prob, xs, ys; window_extent = 20.0,
                 overlap = 6.0, h = 0.5, order = 20, rng_seed = 0)
    wsol42 = windowed_path_network_solve(prob, xs, ys; window_extent = 20.0,
                 overlap = 6.0, h = 0.5, order = 20, rng_seed = 42)

    # ------------------------------------------------------------------------
    # FSEAM.1  NON-TRIVIALITY — the global seed is genuinely threaded into
    #          every window, so the two runs really did re-randomise the walk.
    #          Without this, a seed-ignoring pipeline would pass FSEAM.2 for
    #          free; this is the anti-gaming guard that gives FSEAM.2 teeth.
    # ------------------------------------------------------------------------
    @testset "FSEAM.1: the two global seeds genuinely re-randomised the walk" begin
        # Same tiling ⇒ same window count (precondition for the eachindex scan).
        @test length(wsol0.window_sols) == length(wsol42.window_sols)

        # The global seed IS threaded into every window: per-window seeds differ.
        @test wsol0.window_seeds != wsol42.window_seeds

        # ... and that re-randomisation reaches the COMPUTATION: at least one
        # window's internal Stage-1 walk (its visited-node sequence) differs.
        n_diff = count(wi -> wsol0.window_sols[wi].visited_z !=
                             wsol42.window_sols[wi].visited_z,
                       eachindex(wsol0.window_sols))
        @printf("[FSEAM.1] windows=%d  seed-differing walks=%d  window_seeds differ=%s\n",
                length(wsol0.window_sols), n_diff,
                wsol0.window_seeds != wsol42.window_seeds)
        @test n_diff ≥ 1
    end

    # ------------------------------------------------------------------------
    # FSEAM.2  THE CURE — the composited POLE FIELD is seed-invariant.  Pole-
    #          set agreement ≥ 0.90 in BOTH directions plus soft count
    #          agreement.  Monolithic scores 77.6% here (RED); the windowed
    #          composite scores 97.2% (GREEN).  Threshold pinned at the
    #          measured floor, not relaxed (Rule 2).
    # ------------------------------------------------------------------------
    @testset "FSEAM.2: composited pole field is seed-invariant (≥0.90 both ways)" begin
        c0  = windowed_extract_poles(wsol0)
        c42 = windowed_extract_poles(wsol42)

        # Both seeds must actually find a pole field, else the fractions are
        # meaningless (fail loud rather than pass on empty sets).
        @test !isempty(c0)
        @test !isempty(c42)

        m_0_42 = match_frac(c0, c42; tol = 0.5)   # recall of seed-0 poles in seed-42
        m_42_0 = match_frac(c42, c0; tol = 0.5)   # recall of seed-42 poles in seed-0
        ncap   = max(length(c0), length(c42))
        @printf("[FSEAM.2] |poles| seed0=%d seed42=%d  match(0→42)=%.1f%%  match(42→0)=%.1f%%  Δcount=%d\n",
                length(c0), length(c42), 100 * m_0_42, 100 * m_42_0,
                abs(length(c0) - length(c42)))

        # Seed-invariance, asserted both directions (catches moved poles AND
        # appeared/vanished poles).
        @test m_0_42 ≥ 0.90
        @test m_42_0 ≥ 0.90
        # Soft count agreement: pole counts within 10% of the larger set.
        @test abs(length(c0) - length(c42)) ≤ 0.10 * ncap
    end
end

# ----------------------------------------------------------------------------
# FSEAM.3  SINGLE-WINDOW DEGENERACY WARNING (bead padetaylor-fqwf).
#          `_tile_centers` yields ONE window on an axis when window_extent ≥
#          that axis's span, i.e. the "windowed" composite silently becomes
#          the monolithic path-dependent solve the module exists to replace.
#          The impl must @warn (naming the bead and the fix) on that axis,
#          and must stay SILENT when the tiling is genuine.  Cheap problem:
#          u'' = u' on a 5×5 lattice, so this adds well under a second.
# ----------------------------------------------------------------------------
@testset "FSEAM.3: warn when window_extent ≥ span collapses to one window" begin
    f    = (z, u, up) -> up
    prob = PadeTaylorProblem(f, (1.0 + 0im, 1.0 + 0im), (0.0 + 0im, 3.0 + 0im);
                             order = 20)
    xs = ys = range(-1.0, 1.0; length = 5)       # span 2.0 on both axes

    # window_extent 20 ≥ span 2 on BOTH axes ⇒ exactly two warnings (x, y),
    # each naming the bead and the fix.
    rx = r"window_extent=20.0 ≥ x-span=2.0.*padetaylor-fqwf.*reduce window_extent"
    ry = r"window_extent=20.0 ≥ y-span=2.0.*padetaylor-fqwf.*reduce window_extent"
    wsol1 = @test_logs (:warn, rx) (:warn, ry) windowed_path_network_solve(
        prob, xs, ys; window_extent = 20.0, overlap = 0.0, h = 0.5)
    @test length(wsol1.centers) == 1                 # it really is 1×1
    @test all(isfinite, wsol1.grid_u)                # and still a valid solve

    # Mixed: collapse on y only (xs span 4 > extent 3; ys span 2 < 3).
    xs4 = range(-2.0, 2.0; length = 9)
    wsol_y = @test_logs (:warn, r"ONE window on y") windowed_path_network_solve(
        prob, xs4, ys; window_extent = 3.0, overlap = 0.0, h = 0.5)
    @test length(wsol_y.centers) == 2                # 2×1 tile grid

    # Genuine tiling (extent 1 on span 2 ⇒ 2×2): NO warning may be emitted.
    wsol4 = @test_logs min_level = Base.CoreLogging.Warn windowed_path_network_solve(
        prob, xs, ys; window_extent = 1.0, overlap = 0.5, h = 0.5)
    @test length(wsol4.centers) == 4
end


# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — the orchestrator EXECUTES each mutation,
# confirms the predicted RED, then restores the source byte-clean.
#
#   M-monolithic — swap the CURE for the BUG it replaces.  In FSEAM.2, replace
#     the two `windowed_extract_poles(wsol{0,42})` pole sets with the monolithic
#     baseline at the SAME config over the FULL grid:
#         grid = ComplexF64[x + im * y for y in ys for x in xs]
#         c0  = extract_poles(path_network_solve(prob, grid; h = 0.5, order = 20, rng_seed = 0))
#         c42 = extract_poles(path_network_solve(prob, grid; h = 0.5, order = 20, rng_seed = 42))
#     Predicted bite: FSEAM.2 goes RED — the bidirectional match collapses to
#     ~77.6% (measured), failing the 0.90 threshold.  Proves the metric
#     DISCRIMINATES the path-dependence bug (it is not passing on noise).
#     FSEAM.1 is unaffected (still a live seed-varying pipeline).  Restore.
#
#   M-frozen-seed — give the anti-gaming guard something to catch.  In
#     src/WindowedComposite.jl, make `_window_seed` IGNORE the global seed:
#         _window_seed(rng_seed::Integer, wi::Integer) = wi
#     Predicted bite: FSEAM.1 goes RED — `window_seeds` are now identical across
#     the two global seeds (`wsol0.window_seeds == wsol42.window_seeds`) and
#     every window's walk is identical (`n_diff == 0`).  Proves FSEAM.1 has
#     teeth: a seed-ignoring pipeline is caught here, so a green FSEAM.2 cannot
#     be a frozen-pipeline artefact.  FSEAM.2 would still pass (trivially seed-
#     invariant) — which is EXACTLY the gaming FSEAM.1 exists to expose.  Restore.
#
# Note: the complementary guarantee that the composite suppresses only SPURIOUS
# seed-dependent poles (and does not DROP REAL poles) is bead padetaylor-ingn
# (Weierstrass-℘ recall/precision truth anchor), out of scope for this seed-
# invariance gate.
#
# Mutation F3 (bead padetaylor-fqwf, 2026-08-23; FSEAM.3 only, GREEN 7/7
#   unmutated) — in `src/WindowedComposite.jl` `_warn_single_window`, change
#   the guard `n == 1 && span > 0 || return` to `return` (never warn).
#   Measured bite: 2 RED of 7 — the two `@test_logs (:warn, …)` blocks
#   (x+y, and y-only) each fail with "Log Test Failed" (0 warnings seen);
#   the silent-tiling check stays GREEN as it must.  Restored.
#
# STANDALONE RUN:
#   julia --project=. test/field_seam_test.jl
# ============================================================================
