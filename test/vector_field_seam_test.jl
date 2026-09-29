# test/vector_field_seam_test.jl
#
# ============================================================================
# VECTOR FIELD SEAM: two-seed FIELD-LEVEL path-independence gate for
# `vector_path_network_solve` (bead padetaylor-lg4y).  The vector port of the
# scalar gate `test/field_seam_test.jl` (FSEAM, bead sny7 / ADR-0034).
#
# WHY THIS FILE EXISTS.  Worklog 077 root-caused the scalar Fig 4.7 seam as
# WALK path-dependence: two seed-different Stage-1 trees accumulate different
# IVP error and the field disagrees where their subtrees meet.  The vector walk
# (src/VectorPathNetwork.jl, `rng` kwarg :461/:655, shuffle :723) had only
# VPN.4.3 (test/vector_path_network_test.jl:853-885): one Riccati pole value
# agrees to 1e-7.  Nothing compared the vector FIELD across two walks.
#
# THE FIXTURE.  The shipped P_I^(2) tritronquée wedge march — the exact Stage-B
# computation of figures/_kkg_pi2_helpers.jl (`kkg_pole_field`: BVP-derived
# seed at z = -3, order 24, h = 0.1, the 20-point r·e^{iθ} target fan).  P_I^(2)
# has the Painlevé property, so the solution is single-valued meromorphic: two
# walks that reach the same z MUST agree to the walk's numerical error.
#
# THE METRIC (why not the scalar pole-set match).  The scalar gate compares
# extracted POLE SETS.  On this fixture the shared-Q pole extraction is itself
# jittery at the 0.05-0.5 scale (PI2S.10 measured median two-run pole
# disagreement 0.35; this worklog measured 81-91 % pole-set match @0.5 between
# seeds whose FIELDS agree to 4e-7), so a pole-set metric cannot tell a walk
# seam from extraction noise.  This gate therefore compares the FIELD itself:
# run A is evaluated (Stage-2 `fine_grid`, true-radius gated, extrapolate=false)
# at every visited node of run B, and compared with B's stored node state
# `visited_y` — all four components of y = (u, u', u'', u''').  Only nodes
# inside A's honest validity discs are compared (the rest come back NaN by
# the B1 gate's contract); the covered count is pinned so the gate cannot pass
# on an empty overlap.
#
# MEASURED (docs/worklog/082-vector-two-seed-field-gate.md):
#   seed 0 at seed-42 nodes : 158/824 covered, max rel diff 5.8e-6 (u: 3.8e-7)
#   seed 42 at seed-0 nodes : 164/848 covered, max rel diff 3.4e-6 (u: 3.7e-7)
#     (both overlaps reach the far wedge |z| ≥ 5.5: 96 / 100 covered nodes)
#   DEFAULT order (rng = nothing) at seed-0 nodes: 215/848 covered, 65 nodes
#     > 1e-4, max rel diff 2.18 — a REAL field disagreement, confined to the
#     FAR wedge (every bad node has |z| ∈ [5.82, 8.0], arg ∈ [-0.06, 0.52]);
#     the inner wedge (|z| < 5.5, 77 nodes) agrees to 9.5e-8.  Independent
#     cross-check (worklog 082 §3): the default order re-run at h = 0.05, and
#     seed 7 at h = 0.05 and 0.1, all agree with seed 0 in the far wedge
#     (≤ 7.7e-5) and all disagree with the default h = 0.1 walk (0.24-0.29).
#     So the outlier is the DEFAULT-order h = 0.1 walk — the exact walk the
#     shipped P_I^(2) figure (figures/_kkg_pi2_helpers.jl:339) runs.  Pinned
#     below as VFSEAM.3 (@test_broken) — NOT papered over.
#
# Mutation proof in the file footer.
# ============================================================================

using Test, PadeTaylor, Printf, Random
include(joinpath(@__DIR__, "..", "figures", "_kkg_pi2_helpers.jl"))

# cross_field(wa, wb) — per-component relative disagreement of walk `wa`'s
# Stage-2 field (evaluated at `wb.visited_z`) against `wb`'s stored node
# states.  Relative to max(1, |y_B|) so pole-adjacent large values are not
# over-weighted.  Returns (rel::Matrix, covered::BitVector): rel[k, c] for node
# k, component c; covered[k] iff every component is finite (inside A's disc).
function cross_field(wa_at_b, wb)
    n, d = length(wb.visited_z), length(first(wb.visited_y))
    rel = [abs(wa_at_b.grid_y[k][c] - wb.visited_y[k][c]) /
           max(1.0, abs(wb.visited_y[k][c])) for k in 1:n, c in 1:d]
    covered = vec(all(isfinite, rel; dims = 2))
    return rel, covered
end

@testset "Vector field seam: two-seed path-independence (VFSEAM)" begin
    bvp    = kkg_bvp_solve()
    f      = painleve_hierarchy(:I, 2; t = KKG_T)
    y_seed = ComplexF64.(bvp(ComplexF64(KKG_Z_SEED)))
    prob   = VectorPadeTaylorProblem(f, y_seed,
                 (ComplexF64(KKG_Z_SEED), 8.0 + 0.0im); order = KKG_PN_ORDER)
    targets = kkg_target_wedge()
    walk(rng; fine = nothing) = vector_path_network_solve(prob, targets;
        order = KKG_PN_ORDER, h = KKG_PN_H, rng = rng, fine_grid = fine)

    w0  = walk(MersenneTwister(0))
    w42 = walk(MersenneTwister(42))

    # ------------------------------------------------------------------------
    # VFSEAM.1 NON-TRIVIALITY (anti-gaming, as FSEAM.1): the two seeds really
    # built different trees, so a green VFSEAM.2 is path-independence, not a
    # frozen pipeline that ignores `rng`.
    # ------------------------------------------------------------------------
    @testset "VFSEAM.1: the two seeds genuinely re-randomised the walk" begin
        @test w0.visited_z != w42.visited_z
        @test length(w0.visited_z) > 100 && length(w42.visited_z) > 100
    end

    # ------------------------------------------------------------------------
    # VFSEAM.2 THE GATE: each seeded walk's field reproduces the other seed's
    # node states to ≤ 1e-4 relative, all four components, both directions.
    # Measured max 5.8e-6; the seam regime (VFSEAM.3) is 1e-1..2.  Threshold
    # sits between the regimes by > 1 decade on each side.
    # ------------------------------------------------------------------------
    @testset "VFSEAM.2: seeded walks agree field-wide (both directions)" begin
        for (wa_rng, wb, tag) in ((MersenneTwister(0), w42, "0@42"),
                                  (MersenneTwister(42), w0, "42@0"))
            wa = walk(wa_rng; fine = wb.visited_z)
            rel, cov = cross_field(wa, wb)
            @printf("[VFSEAM.2] %s covered=%d/%d  max rel per comp = %s\n",
                    tag, count(cov), length(cov),
                    join((@sprintf("%.2e", maximum(rel[cov, c]))
                          for c in axes(rel, 2)), ", "))
            nfar = count(k -> cov[k] && abs(wb.visited_z[k]) ≥ 5.5, eachindex(cov))
            @printf("[VFSEAM.2] %s covered nodes with |z| ≥ 5.5: %d\n", tag, nfar)
            @test count(cov) ≥ 100           # non-empty overlap (measured 158/164)
            @test nfar ≥ 50                  # overlap reaches the far wedge (96/100)
            @test maximum(rel[cov, :]) ≤ 1.0e-4
        end
    end

    # ------------------------------------------------------------------------
    # VFSEAM.3 THE SEAM FOUND (bead lg4y measurement).  The DEFAULT target
    # order (rng = nothing — the order the shipped figure uses) disagrees with
    # the seeded consensus by O(1) relative on 65 FAR-wedge nodes (|z| ≥ 5.82).
    # The invariant "every ordering agrees with every other" is therefore FALSE
    # today: marked @test_broken so a fix flips it to an unexpected pass (Julia
    # reports that as an Error — delete the marker then).  The localisation
    # (no bad node at |z| < 5.5; the inner wedge agrees ≤ 1e-4) and the size of
    # the defect (≤ 100 bad nodes) are pinned as TRUE facts so it cannot grow
    # silently.
    # ------------------------------------------------------------------------
    @testset "VFSEAM.3: default ordering vs seeded consensus (known seam)" begin
        wd = walk(nothing; fine = w0.visited_z)
        rel, cov = cross_field(wd, w0)
        bad = findall(k -> cov[k] && maximum(rel[k, :]) > 1.0e-4, eachindex(cov))
        rmin = isempty(bad) ? Inf : minimum(k -> abs(w0.visited_z[k]), bad)
        inner = findall(k -> cov[k] && abs(w0.visited_z[k]) < 5.5, eachindex(cov))
        @printf("[VFSEAM.3] default@0 covered=%d/%d  bad(>1e-4)=%d  max rel=%.2e  min |z_bad|=%.3f  inner(|z|<5.5) covered=%d max rel=%.2e\n",
                count(cov), length(cov), length(bad), maximum(rel[cov, :]), rmin,
                length(inner), maximum(rel[inner, :]))
        @test count(cov) ≥ 100
        @test_broken maximum(rel[cov, :]) ≤ 1.0e-4      # the seam (bead lg4y)
        @test rmin ≥ 5.5                                  # far wedge only
        @test maximum(rel[inner, :]) ≤ 1.0e-4             # inner wedge agrees
        @test length(bad) ≤ 100                           # measured 63; not growing
    end
end

# ============================================================================
# MUTATION-PROOF RECORD (Rule 4) — executed 2026-09-29, source restored
# byte-clean after each (`git diff src/` empty).  Unmutated: 12 pass, 1 broken (~44 s wall).
#
#   M-frozen-seed — src/VectorPathNetwork.jl:723, replace the shuffle line
#     `rng === nothing || (target_list = shuffle(rng, target_list))` by a
#     comment.  Measured bite: VFSEAM.1 RED (`w0.visited_z != w42.visited_z`
#     fails — identical 795-node trees); VFSEAM.2 then passes trivially at
#     4.9e-16 on 795/795 nodes — exactly the gaming VFSEAM.1 exists to expose;
#     VFSEAM.3's @test_broken unexpectedly passes (Error).  Restored.
#   M-drift — src/VectorPathNetwork.jl, directly after the `_select_wedge`
#     call (:819-821) insert `y_new = y_new .* (1 + 1e-5)`: a per-step relative
#     error, so accumulated error depends on path LENGTH (a manufactured walk
#     path-dependence).  Measured bite: VFSEAM.2 RED in both directions (max
#     rel 1.75e+02 for 0@42, 8.11e+01 for 42@0 vs threshold 1e-4); VFSEAM.3
#     localisation also RED (min |z_bad| = 1.82, inner max rel 20.0).  Restored.
#
# STANDALONE RUN:
#   julia --project=. test/vector_field_seam_test.jl
# ============================================================================
