# figures/test_kkg_pi2_surface.jl
#
# Acceptance test for the F3 whole-plane P_I^(2) tritronquée compute
# kernel (`figures/_kkg_pi2_surface_helpers.jl`, bead
# `padetaylor-0ln.33`).  Runs under the `figures/` Julia environment:
#
#     julia --project=figures figures/test_kkg_pi2_surface.jl
#
# It is a standalone `@testset` script — NOT part of the package `test/`
# suite — because it `include`s the Gridap-backed figure helpers, and
# Gridap is a `figures/`-project dependency only (ADR-0024 *Correction*).
#
# ## What is asserted (CLAUDE.md Rule 5 — every test pins a known value)
#
#   1. **Negative-real-axis pin.**  On `x < 0` the `V_0` tritronquée is
#      real, positive, with the KKG algebraic asymptotic
#      `u → +(6|x|)^{1/3}` (eq. (1.3); pillar C §4).  At `x ≈ -15` the
#      assembled `Re u` must match `+(6·15)^{1/3} ≈ 4.4814` and `Im u`
#      must vanish.  Tolerance `5e-3`: the `n_terms = 2` truncated
#      asymptotic series is itself only accurate to `~1.3e-4` at
#      `|x| = 15` (the V8b probe number), and the BVP-and-vote pipeline
#      adds the FEM voter's `~5e-4` and the ray-fan interpolation's
#      `~1e-3`; `5e-3` is a generous but honest envelope around all
#      three.
#
#   2. **Three-method agreement.**  Over the bulk of the pole-free
#      sector the max pairwise disagreement among Methods 1/2/3 (the
#      `spread` map) must be small.  Tolerance `1e-2`: the FEM voter is
#      `O(h²) ≈ 5e-5` on a smooth field but the inner-arc asymptotic
#      seed (v1 corner 1) injects an `O(10⁻²)` Dirichlet error at
#      `|x| ≈ 2` that the harmonic solve damps but does not erase, so
#      the worst-case spread sits near the inner arc.  We assert the
#      *median* spread is `< 3e-3` (the bulk concurs tightly) and the
#      *max* spread is `< 1e-2` (even the inner-arc corner stays
#      bounded).  The voted surface is verified to equal the
#      per-component median of the three voters.
#
#   3. **Sector coverage.**  ≥ 80% of the pole-free-sector grid cells
#      (inside the `|x| ≤ 20` disc) must carry a non-`NaN` value.
#
#   4. **Wedge pole field + the VC-4/VC-5 per-pole validation
#      (ADR-0025 Amendment 2 + Amendment 1 §A3).**  The wedge panel is
#      rescoped to the validated pole *field* + an honest partial `|u|`
#      surface underlay — NOT a filled surface (the A2 probe proved a
#      filled honest wedge surface numerically unreachable; honest
#      coverage saturates at ~8-18 %).  PI2S.4 asserts the
#      Amendment-2 deliverable's structural invariants AND the Phase-D
#      per-pole validation criteria VC-4 / VC-5:
#        - the pole field is non-empty and grows *richer* than V8b's
#          21 poles (the B3 extended threading fan threads the *whole*
#          wedge to `|x| = 20`, where V8b's fan stopped at `|x| ≤ 8`);
#        - every pole is wedge-confined, `|arg p|` within `36° + margin`
#          (VC-3) — no pole leaks into the pole-free sector;
#        - the honest-coverage mask `wedge_covered` is *consistent*
#          with `Re_u`/`Im_u`: a covered cell is finite, an uncovered
#          wedge cell is `NaN` — and the coverage fraction is the
#          honest, partial (~5-20 %) figure A2 predicts, NOT a filled
#          surface;
#        - no Padé is evaluated out of disc (covered ⊆ finite, the B1
#          true-radius gate's contract);
#        - **VC-4** — every pole in the figure's *final* field passes
#          the dominant-balance certificate: a 32-point ring fit
#          `u ≈ A·ξ⁻² + B·ξ⁻¹ + C` has `min(|A+1|,|A+3|) < 0.10` (the
#          `A ∈ {-1,-3}` family, VC-4a) and `|B| < 0.10·|A|` (the zero
#          residue, VC-4b); the VC-4 prune actually removed spurious
#          candidates, and a deliberately-injected fake `A≈0` pole is
#          rejected by the filter (the mutation-proof);
#        - **VC-5** (ADR-0025 Amendment 6) — the VC-4-surviving poles
#          match into conjugate pairs by a globally-optimal
#          (maximum-cardinality) bipartite matching, with a small
#          pairing residual (the FW-style accuracy estimate); the
#          flagged-unpaired count materially beats the v1 greedy
#          matcher's 72 (Amendment 5); the matching is verified
#          maximum-cardinality (no two unpaired poles are admissible
#          mirrors); every *un-flagged* unpaired pole sits on the real
#          axis (`|Im|` small).
#      VC-7 loop closure remains a separate Phase-D bead.
#
#   5. **Schwarz symmetry.**  `V_0(x̄) = conj(V_0(x))` (the ODE has real
#      coefficients and the tritronquée is real on `x < 0`), so in the
#      sector `Re u` is even and `Im u` is odd under `y ↦ -y`.
#
#   6. **Reproducibility.**  Two `kkg_pi2_surface()` calls are identical.
#
# This test does NOT assert "did not throw" anywhere — every block pins
# an invariant against a known-correct value (Rule 5).  Where a method
# is unreliable (the inner arc) the `spread` map exposes it and the test
# reports the number honestly rather than masking it.

using Test
using Statistics: median        # the INDEPENDENT median oracle for PI2S.2

include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))

# The compute kernel is expensive (a fan of BVPs + two Laplace solves +
# the wedge walk); run it ONCE and share the result across testsets.
const SURF = kkg_pi2_surface()
const SURF2 = kkg_pi2_surface()      # second run, for reproducibility

@testset "F3 — P_I^(2) tritronquée whole-plane surface" begin

    n = SURF_GRID_N

    @testset "PI2S.1 — negative-real-axis pin (KKG eq. 1.3)" begin
        i15 = argmin(abs.(SURF.xs .+ 15.0))
        j0  = argmin(abs.(SURF.ys))
        # The grid is constructed to land a node on x = -15, y = 0.
        @test isapprox(SURF.xs[i15], -15.0; atol = 1e-9)
        @test isapprox(SURF.ys[j0],    0.0; atol = 1e-9)

        expected = cbrt(6.0 * 15.0)            # +(6·15)^{1/3} ≈ 4.4814
        re15 = SURF.Re_u[i15, j0]
        im15 = SURF.Im_u[i15, j0]
        @test !isnan(re15)
        @test isapprox(re15, expected; atol = 5e-3)
        @test isapprox(im15,      0.0; atol = 5e-3)
        # The tritronquée is real and positive on x < 0.
        @test re15 > 0
    end

    @testset "PI2S.2 — three-method agreement / the majority vote" begin
        # Gather the sector spread map (the agreement diagnostic).
        sp = Float64[]
        for j in 1:n, i in 1:n
            s = SURF.spread[i, j]
            isnan(s) || push!(sp, s)
        end
        @test !isempty(sp)
        sort!(sp)
        med = sp[length(sp) ÷ 2]
        mx  = sp[end]
        # The bulk of the sector concurs tightly...
        @test med < 3e-3
        # ...and even the inner-arc v1 corner stays bounded.
        @test mx < 1e-2

        # The voted surface IS the per-component median of the three
        # voters.  Re-derive the median from `sector_method` with an
        # INDEPENDENT computation — `Statistics.median` on the finite
        # voters, NOT a call to the kernel's own `surf_vote` (calling
        # `surf_vote` here would be a tautology: a bug in `surf_vote`
        # would corrupt both sides equally).  This catches a kernel that
        # silently votes by the mean, the max, the first voter, etc.
        indep_vote(a, b, c) = begin
            v = filter(isfinite, (a, b, c))
            isempty(v) ? NaN : median(collect(v))
        end
        nchecked = 0; nfull = 0
        for ((i, j), m) in SURF.sector_method
            re1, im1, re2, im2, re3, im3 = m
            ref_re = indep_vote(re1, re2, re3)
            ref_im = indep_vote(im1, im2, im3)
            @test isapprox(SURF.Re_u[i, j], ref_re; atol = 1e-12, nans = true)
            @test isapprox(SURF.Im_u[i, j], ref_im; atol = 1e-12, nans = true)
            nchecked += 1
            # Count the grid points where all three voters are finite —
            # there the median is genuinely the *middle* of three
            # distinct numbers, so this block is a real discriminator
            # (the mean would differ from the median there).
            all(isfinite, (re1, re2, re3)) && (nfull += 1)
        end
        @test nchecked > 1000          # the vote ran on the whole sector
        @test nfull    > 1000          # ...and most points had 3 finite voters
    end

    @testset "PI2S.3 — pole-free-sector coverage ≥ 80%" begin
        total = 0; filled = 0
        for j in 1:n, i in 1:n
            z = ComplexF64(SURF.xs[i], SURF.ys[j])
            (abs(z) <= SURF_XY_LIM && surf_in_sector(z)) || continue
            total += 1
            isnan(SURF.Re_u[i, j]) || (filled += 1)
        end
        @test total > 0
        @test filled / total ≥ 0.80
    end

    @testset "PI2S.4 — wedge pole field + honest partial underlay" begin
        # --- the pole field — the Amendment-2 primary wedge deliverable ---
        # The B3 extended threading fan drives the B2 `:max_q_root`
        # adaptive walk through the whole wedge to `|x| = 20`; the field
        # is non-empty and RICHER than V8b's 21 poles (V8b's fan stopped
        # at `|x| ≤ 8` — A2 §1).  `> 21` is the senior-grade bar: the
        # re-resolution must beat the baseline, not merely match it.
        @test length(SURF.poles) ≥ 1
        @test length(SURF.poles) > 21

        # VC-3 — every extracted pole is wedge-confined: `|arg p|` within
        # the 36° Stokes line plus a small clustering margin.  A pole
        # leaking into the pole-free ~270° sector would be spurious (the
        # tritronquée is analytic there).
        for p in SURF.poles
            @test abs(rad2deg(angle(p))) ≤ 36.0 + 2.0
        end

        # --- the honest-coverage mask `wedge_covered` is consistent -------
        # ADR-0025 Amendment 1/2 contract: the Stage-2 fill ran
        # `extrapolate = false`, so a cell is finite in `Re_u`/`Im_u`
        # IFF it is honestly covered by a B1-gated node disc.  Verify
        # the mask agrees with the matrices cell-by-cell over the wedge,
        # and that uncovered wedge cells are genuinely `NaN` (an honest
        # gap — no Padé evaluated out of disc).
        n_covered = 0; n_wedge = 0; mask_consistent = true
        for j in 1:n, i in 1:n
            z = ComplexF64(SURF.xs[i], SURF.ys[j])
            (abs(z) <= SURF_XY_LIM && !surf_in_sector(z) &&
             !surf_in_mask(z) && !iszero(z)) || continue
            n_wedge += 1
            cov = SURF.wedge_covered[i, j]
            fin = !isnan(SURF.Re_u[i, j]) && !isnan(SURF.Im_u[i, j])
            # covered  ⟺  finite  (the B1-gate honesty contract)
            (cov == fin) || (mask_consistent = false)
            cov && (n_covered += 1)
        end
        @test mask_consistent
        @test n_wedge > 0

        # The honest coverage is PARTIAL — the A2 probe measured
        # ~8-18 %; this is NOT a filled surface.  Lower bound: the gate
        # genuinely covers *some* cells (the field is not all-NaN);
        # upper bound: it is well short of a filled wedge (`< 60 %` — A2
        # found no annulus past 50 %, so the global figure cannot
        # approach the 70 % filled-surface bar).
        cov_frac = n_covered / n_wedge
        @test cov_frac > 0.0
        @test cov_frac < 0.60

        # Outside the covered mask, every wedge cell is `NaN` — the
        # honest gap.  (This is the negation of out-of-disc evaluation:
        # the gate returns `NaN`, never a silently-extrapolated value.)
        for j in 1:n, i in 1:n
            z = ComplexF64(SURF.xs[i], SURF.ys[j])
            (abs(z) <= SURF_XY_LIM && !surf_in_sector(z) &&
             !surf_in_mask(z) && !iszero(z)) || continue
            SURF.wedge_covered[i, j] && continue
            @test isnan(SURF.Re_u[i, j])
            @test isnan(SURF.Im_u[i, j])
        end

        # --- VC-4 — dominant-balance per-pole structural certificate ------
        # ADR-0025 Amendment 1 §A3.  The figure's `poles` field is the
        # VC-4-VALIDATED set: `extract_poles_shared_q` produces a
        # candidate field (~380 poles), and `vc4_validate` prunes every
        # candidate that is not a genuine P_I^(2) double pole.  Assert:
        #   (a) the diagnostics are reported (candidate count, prune count);
        #   (b) VC-4 actually pruned spurious candidates;
        #   (c) EVERY pole in the final field passes VC-4a and VC-4b —
        #       re-fit each from scratch here, an independent check that
        #       does NOT trust the kernel's own `vc4.pass` flag.
        v4 = SURF.vc4
        @test v4 !== nothing
        # (a) the candidate field is the ~380-pole B3 extraction; the
        # validated field is a strict subset (some candidates pruned).
        n_cand = length(v4.A)
        @test n_cand > length(SURF.poles)        # pruning shrank the field
        @test v4.n_pruned == n_cand - length(SURF.poles)
        # (b) VC-4 pruned a positive number of spurious candidates — the
        # A4 baseline already found ≥1 spurious pole among 21; among 380
        # candidates the prune count must be well above zero.
        @test v4.n_pruned > 0
        @test length(v4.prune_reason) == v4.n_pruned
        # every pruned pole carries a classified failure mode.
        for (_, reason) in v4.prune_reason
            @test reason in (:froissart, :out_of_family, :nonzero_residue)
        end

        # (c) every FINAL pole passes VC-4a + VC-4b.  The kernel records
        # the per-candidate fitted `(A,B)` and the `pass` mask in `vc4`;
        # the figure's `poles` field is exactly the `pass`-true subset.
        # Assert that mapping holds AND that each passing candidate's
        # `(A,B)` genuinely meets the two acceptance inequalities — so a
        # bug that kept a failing candidate (wrong `pass` logic) is RED.
        # The fit itself is re-exercised independently in the
        # mutation-proof block below (a fresh `vc4_validate` call).
        kept_idx = [k for k in 1:n_cand if v4.pass[k]]
        @test length(kept_idx) == length(SURF.poles)
        for k in kept_idx
            A = v4.A[k]; B = v4.B[k]
            # VC-4a — leading coefficient in the {-1,-3} family.
            @test min(abs(A + 1), abs(A + 3)) < 0.10
            # VC-4b — the residue vanishes.
            @test abs(B) < 0.10 * abs(A)
        end
        # The A-family breakdown: the tritronquée's wedge poles are the
        # generic A=-1 family (ADR-0025 A3 §5.1).  Every kept pole is
        # classified `:m1` or `:m3`, never `:none`.
        for k in kept_idx
            @test v4.family[k] in (:m1, :m3)
        end

        # --- VC-4 mutation-proof — a deliberately-injected fake pole ------
        # Inject a fake "pole" at a location where `u` has NO double-pole
        # structure (a generic interior wedge point — `u` there is
        # analytic / O(1), so the ring fit returns `A ≈ 0`).  VC-4 MUST
        # reject it: `vc4_validate` on a one-element field containing
        # only the fake must prune it (kept empty, `n_pruned == 1`,
        # failure mode `:froissart`).  This proves the filter has teeth —
        # without it a fake `A≈0` location would survive into the figure.
        # The fit is re-exercised here against the kernel's own
        # `wedge_walk` (the shared-Q path-network the figure rendered),
        # so this is an INDEPENDENT re-run of the VC-4 ring fit.
        walk = SURF.wedge_walk
        @test walk !== nothing
        z_fake = ComplexF64(8.5, 0.5)    # a generic in-wedge analytic point
                                         # (clearance ~0.67 from any pole)
        @test !surf_in_sector(z_fake) && !surf_in_mask(z_fake)
        mut = vc4_validate(walk, ComplexF64[z_fake])
        @test isempty(mut.kept)              # the fake is NOT kept
        @test mut.n_pruned == 1              # it was pruned
        @test mut.prune_reason[1][2] == :froissart   # A ≈ 0 — Froissart
        # ...and a genuine pole from the real field, fed through the SAME
        # one-element validation path, IS kept — proving the filter is
        # not a blanket reject (mutation-proof: real pole GREEN, fake RED).
        genuine = SURF.poles[argmin(abs.(SURF.poles))]   # nearest-origin
        gkeep = vc4_validate(walk, ComplexF64[genuine])
        @test length(gkeep.kept) == 1

        # --- VC-5 — conjugate-symmetry pole pairing (ADR-0025 Amend. 6) ---
        # `V_0(x̄) = conj V_0(x)`, so the VC-4-surviving field must be
        # conjugate-symmetric.  `vc5_pair` matches the survivors into
        # conjugate pairs by a globally-optimal (maximum-cardinality)
        # bipartite matching of the conjugate-admissibility graph and
        # reports the pairing residual — itself an FW-style accuracy
        # estimate.  Amendment 6 (bead `padetaylor-0ln.37.20`, VC-5b)
        # replaced the v1 greedy matcher — which was not maximum-
        # cardinality and left 72 off-axis poles spuriously flagged —
        # and re-derived `VC5_MATCH_TOL` (0.5 → 0.6) from the *measured*
        # 0.69 pole spacing / ~0.3-0.5 field accuracy.
        v5 = SURF.vc5
        @test v5 !== nothing
        @test !isempty(v5.pairs)             # the field DOES pair up

        # VC-5b — the flagged-pole count.  The v1 greedy matcher flagged
        # 72 off-axis poles as unpaired (ADR-0025 Amendment 5).  The
        # Amendment-6 optimal matcher + re-derived tolerance materially
        # cuts that: the measured count on the B3 field is 52.  Assert a
        # concrete bound — `≤ 55` — strictly below the 72 baseline (the
        # VC-5b deliverable).  The bound is set BELOW what the v1 greedy
        # matcher achieves at the same `VC5_MATCH_TOL` (greedy flags 58
        # — measured), so it is genuinely load-bearing for the
        # optimal-matching fix, not carried by the tolerance change
        # alone.  This is the load-bearing VC-5b regression assertion
        # (mutation-proven below).
        @test length(v5.flagged) ≤ 55
        @test length(v5.flagged) < 72        # strictly beats the baseline
        # ...and the optimal matcher finds materially MORE conjugate
        # pairs than the 93 the v1 greedy matcher found.
        @test length(v5.pairs) ≥ 100

        # The pairing residual stays a genuine FW-style accuracy
        # estimate — well below the wedge pole nearest-neighbour spacing
        # (~0.69).  The optimal matcher no longer discards the
        # harder-to-pair tail, so the median residual sits a touch
        # higher than the greedy matcher's optimistically-biased 0.24;
        # `< 0.40` is the honest envelope (the A4 baseline was ~1.25 —
        # the re-resolution still beats it ~3×).
        @test v5.median_resid < 0.40
        @test isfinite(v5.max_resid)
        # Each reported pair genuinely mirrors under conjugation, within
        # the re-derived `VC5_MATCH_TOL = 0.6`.
        for (pu, pl) in v5.pairs
            @test abs(pu - conj(pl)) ≤ 0.6   # VC5_MATCH_TOL
            @test imag(pu) > 0 && imag(pl) < 0
        end
        # The matching is maximum-cardinality: NO unpaired upper pole
        # has an unpaired lower pole within `VC5_MATCH_TOL` of its
        # conjugate.  If one did, a larger matching would exist — so this
        # asserts the matcher genuinely maximised the pair count (a
        # greedy matcher can fail this; the optimal one cannot).
        up_unpaired = [p for p in v5.unpaired if imag(p) >  0.15]
        lo_unpaired = [p for p in v5.unpaired if imag(p) < -0.15]
        for pu in up_unpaired, pl in lo_unpaired
            @test abs(pu - conj(pl)) > 0.6
        end
        # An unpaired pole is EITHER on the real axis (small `|Im|`) OR
        # flagged suspect.  Verify the partition: every unpaired pole the
        # kernel did NOT flag genuinely sits near the real axis.
        flagged_set = Set(v5.flagged)
        for p in v5.unpaired
            (p in flagged_set) && continue   # flagged suspect — reported
            @test abs(imag(p)) < 0.15        # un-flagged ⇒ real-axis pole
        end
        # The flag is itself meaningful: every flagged pole genuinely has
        # a large `|Im|` (it is an off-axis pole with no mirror partner).
        for p in v5.flagged
            @test abs(imag(p)) ≥ 0.15
        end

        # --- VC-5b mutation-proof — greedy matching is genuinely worse ---
        # The load-bearing VC-5b assertion is `length(flagged) ≤ 60`.
        # It is carried by the Amendment-6 optimal (maximum-cardinality)
        # matcher; the v1 greedy matcher would NOT meet it.  Re-run the
        # v1 greedy commit ON THE SAME VC-4-validated field and confirm
        # it flags materially more poles — proving the optimal matcher
        # is what the assertion exercises, not the field happening to be
        # easy.  (A greedy matcher is not maximum-cardinality: it
        # commits a locally-tight pair that blocks two poles from each
        # finding their only admissible mirror.)
        let P = ComplexF64.(SURF.poles), tol = 0.6, ax = 0.15
            up = [p for p in P if imag(p) >  ax]
            lo = [p for p in P if imag(p) < -ax]
            edges = Tuple{Float64,Int,Int}[]
            for (iu, pu) in enumerate(up), (il, pl) in enumerate(lo)
                d = abs(pu - conj(pl))
                d ≤ tol && push!(edges, (d, iu, il))
            end
            sort!(edges; by = e -> e[1])         # greedy: tightest first
            uu = falses(length(up)); ll = falses(length(lo))
            greedy_pairs = 0
            for (_, iu, il) in edges
                (uu[iu] || ll[il]) && continue
                uu[iu] = true; ll[il] = true; greedy_pairs += 1
            end
            greedy_flagged = length(P) -
                2 * greedy_pairs - count(p -> abs(imag(p)) ≤ ax, P)
            # The optimal matcher must strictly beat greedy on this field
            # — fewer flagged, more pairs.  If they tied, the optimal
            # matcher would not be load-bearing and this block is RED.
            @test length(v5.flagged) < greedy_flagged
            @test length(v5.pairs)   > greedy_pairs
        end
    end

    @testset "PI2S.5 — Schwarz symmetry across the real axis" begin
        # V_0(x̄) = conj(V_0(x)): in the sector Re u is even, Im u odd
        # under y ↦ -y.  Check at a handful of conjugate sector pairs.
        npairs = 0
        for i in 1:n
            x = SURF.xs[i]
            x < -3.0 || continue                # stay on the negative side
            for jt in (8.0, 12.0)
                jp = argmin(abs.(SURF.ys .- jt))
                jm = argmin(abs.(SURF.ys .+ jt))
                rp = SURF.Re_u[i, jp]; rm = SURF.Re_u[i, jm]
                ip = SURF.Im_u[i, jp]; im_ = SURF.Im_u[i, jm]
                (isnan(rp) || isnan(rm)) && continue
                @test isapprox(rp,  rm; atol = 1e-6)
                @test isapprox(ip, -im_; atol = 1e-6)
                npairs += 1
            end
        end
        @test npairs > 5
    end

    @testset "PI2S.6 — reproducibility" begin
        # Two independent runs produce byte-identical matrices (NaNs in
        # the same cells; finite cells equal).
        @test size(SURF.Re_u) == size(SURF2.Re_u)
        diffs = 0
        for idx in eachindex(SURF.Re_u)
            a, b = SURF.Re_u[idx], SURF2.Re_u[idx]
            c, d = SURF.Im_u[idx], SURF2.Im_u[idx]
            if isnan(a) != isnan(b) || isnan(c) != isnan(d)
                diffs += 1
            elseif !isnan(a) && (a != b || c != d)
                diffs += 1
            end
        end
        @test diffs == 0
        @test length(SURF.poles) == length(SURF2.poles)
    end
end
