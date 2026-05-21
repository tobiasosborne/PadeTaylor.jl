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
#   2. **Three-method agreement + the C1 lifted resolution.**  Over the
#      bulk of the pole-free sector the max pairwise disagreement among
#      Methods 1/2/3 (the `spread` map) must be small.  Tolerance `1e-2`:
#      the FEM voter is `O(h²) ≈ 5e-5` on a smooth field but the
#      asymptotic-seed truncation floor (`n_terms = 2`,
#      `O(|x|^{-7/3})`) leaves an `O(10⁻³)` boundary-layer error near
#      the inner arc `|x| ≈ 2` that no harmonic-solve `N` can erase, so
#      the worst-case spread sits near the inner arc.  We assert the
#      *median* spread is `< 3e-3` (the bulk concurs tightly) and the
#      *max* spread is `< 1e-2` (even the inner-arc corner stays
#      bounded).  The voted surface is verified to equal the
#      per-component median of the three voters.  PI2S.2 also pins the
#      ADR-0025 **Amendment 11** C1 sector re-resolution (bead
#      `padetaylor-0ln.37.9`, retiring v1 corner C6): the `401²` render
#      grid + `2°` ray fan, and the inner-arc *max* spread is asserted
#      against the measured asymptotic-seed truncation floor (~6.7·10⁻³
#      — Amendment 11 found the inner-arc spread is NOT a Laplace
#      under-resolution, contra Amendment 10; the C1 convergence probe
#      falsified that, see `_kkg_pi2_surface_helpers.jl` `SURF_BVP_N` /
#      `SURF_LAP_NX`) and verified r-localised to the inner annulus.
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
#
#   7. **VC-7 loop-closure quality certificate (ADR-0025 Amendment 3 +
#      Amendment 7).**  `quality_diagnose` (ADR-0016) — the loop-closure
#      certificate — is brought to the figure's vector wedge walk via
#      the D-VC7 vector adapter (`ext/PadeTaylorDiagnosticsExt.jl`).
#      PI2S.7 is the VC-7 figure acceptance gate: it asserts the
#      certificate is non-degenerate and structurally sound (the
#      category tallies partition `n_edges`, the quantiles are monotone
#      — the same shape `test/diagnose_test.jl` pins for the scalar
#      method), the well-closed lobe genuinely closes to machine
#      precision (proving the shared-`Q` `Pᵢ(t)/Q(t)` evaluation
#      pipeline is correct end-to-end), and a meaningful fraction of the
#      genuine in-disc loop closures close tightly.  The load-bearing
#      bars are mutation-proven by corrupting a visited node's shared-`Q`
#      store.  Amendment 7 records the measured numbers and why the raw
#      ΔP_rel distribution is NOT a before/after re-resolution proxy
#      (the metric is independent of the B1/B2 levers).
#
#   5. **Schwarz symmetry.**  `V_0(x̄) = conj(V_0(x))` (the ODE has real
#      coefficients and the tritronquée is real on `x < 0`), so in the
#      sector `Re u` is even and `Im u` is odd under `y ↦ -y`.
#
#   6. **Reproducibility.**  Two `kkg_pi2_surface()` calls are identical.
#
#   10. **VC-10 — the two-run pole-field accuracy indicator (ADR-0025
#      Amendment 9).**  FW 2011 §3.1 shuffles the Stage-1 targets into a
#      random order; FFW 2017 (`:246-247`) makes the cross-ordering
#      disagreement a diagnostic — two runs with different target
#      orderings "differ by approximately the numerical error".  PI2S.10
#      runs the wedge walk twice with two `MersenneTwister` seeds,
#      VC-4-validates both pole fields, matches them nearest-neighbour
#      (the VC-5 maximum-cardinality matcher, reused), and asserts the
#      matched-pole disagreement is small (well below the ~0.69 pole
#      spacing — the manufactured FW-style accuracy certificate of a
#      field with no external oracle).  Each fixed seed's walk stays
#      fully deterministic (PI2F.1.6 reproducibility, lifted to the
#      seeded vector walk).
#
# This test does NOT assert "did not throw" anywhere — every block pins
# an invariant against a known-correct value (Rule 5).  Where a method
# is unreliable (the inner arc) the `spread` map exposes it and the test
# reports the number honestly rather than masking it.

using Test
using Statistics: median        # the INDEPENDENT median oracle for PI2S.2
# `import` (NOT `using`): loading `DelaunayTriangulation` activates the
# `PadeTaylorDiagnosticsExt` package extension (so the VC-7 vector
# `quality_diagnose` method becomes callable), but a bare `using` would
# also splat `DelaunayTriangulation`'s exports — including `Triangulation`
# — into `Main`, colliding with the unqualified `Triangulation` the
# Gridap-backed `_kkg_pi2_gridap_helper.jl` resolves from `Gridap`.
import DelaunayTriangulation     # activates PadeTaylorDiagnosticsExt (VC-7)
using PadeTaylor.Diagnostics: quality_diagnose   # the VC-7 certificate

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
        # --- C1 — the lifted render resolution (ADR-0025 Amendment 11) ----
        # Bead `padetaylor-0ln.37.9` retires v1 corner C6 — the
        # resolution-limited `121²` grid / `6°` ray fan.  Assert the
        # figure is rendered at the lifted resolution: a `401²` Cartesian
        # raster and a `2°` ray fan.  A regression to the old coarse
        # parameters (re-introducing C6) is RED here.
        @test SURF_GRID_N == 401
        @test SURF_RAY_DPHI_DEG == 2.0
        @test n == 401
        @test length(SURF.xs) == 401 && length(SURF.ys) == 401
        @test size(SURF.Re_u) == (401, 401)
        # The `401` raster still lands a node exactly on `x = 0` and on
        # `x = -15` (the PI2S.1 pin) — the C6 retirement keeps the
        # odd-grid invariant.
        @test minimum(abs.(SURF.xs)) < 1e-9
        @test minimum(abs.(SURF.xs .+ 15.0)) < 1e-9

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

        # --- C1 — the inner-arc vote spread (ADR-0025 Amendment 11) -------
        # v1 corner C1: voters 2/3's inner-arc (`|x| ≈ 2`) Dirichlet datum
        # was the KKG `n_terms = 2` asymptotic series, and worklog 057
        # measured the inner-arc triple-method spread at a coarse
        # ~2.9·10⁻³ median.  C2 (Amendment 10) retired C1 — each ray BVP
        # extends inward to `SURF_R_INNER_BC < SURF_R_MIN`, so the inner
        # arc is a BVP-*interior* datum and voters 2/3 source the arc
        # from the ray fan.  The inner-arc spread MEDIAN sits materially
        # below the worklog ~2.9·10⁻³ baseline (`< 5·10⁻⁴`).
        #
        # THE INNER-ARC *MAX* SPREAD — ADR-0025 Amendment 11, bead
        # `padetaylor-0ln.37.9` (the C1 sector re-resolution).  Amendment
        # 10 hypothesised the inner-arc *max* spread (~3.9·10⁻³ on the
        # old `121²` grid) was a voter-2/3 Laplace *under-resolution* and
        # assigned a higher-Laplace-`N` fix to bead `0ln.37.9`.  The C1
        # convergence probe **falsified that hypothesis** (CLAUDE.md
        # Rule 2): both Laplace voters converge geometrically and a
        # higher `N` does NOT move the inner-arc error — solving the
        # spectral problem even with *perfect* dedicated-BVP edge traces
        # leaves the same ~4.8·10⁻³ interior error, `N`-independent.  The
        # genuine root cause is the **asymptotic-seed truncation floor**:
        # each ray BVP pins `(u,u')` at `R_INNER_BC` to the KKG
        # `n_terms = 2` series, whose `O(|x|^{-7/3})` error propagates
        # inward; solving the same ray at several `r_in` gives interior
        # values disagreeing by ~5·10⁻³ at `|x| = 2`.  This floor is
        # irreducible by any `N` / fan density.  Moving `R_INNER_BC`
        # *does* shrink the v1–v2 spread — but only by making the two
        # voters agree on the SAME seed-floored value (the gameable-
        # spread anti-pattern this very testset warns of, below), NOT by
        # improving accuracy: `R_INNER_BC` is therefore left at the
        # Amendment-10 / PI2S.11-certified `1.05`.  The `401²` grid
        # genuinely SHARPENS the figure (retiring v1 corner C6) and, in
        # doing so, samples grid points closer to the worst seed-floor
        # zone at `|x| = 2` — so the honest inner-arc *max* spread is
        # ~6.7·10⁻³ (the seed floor fully exposed; the old `121²`
        # 3.9·10⁻³ was a coarse-sampling artefact, never a real lower
        # floor).  The bound below (`< 8·10⁻³`) is the measured seed
        # floor plus headroom — an HONEST exposure, not a regression.
        # Closing it needs a `n_terms ≥ 3` seed (KKG eq. 7.2 `c₇…`),
        # unimplemented — a deferred bead in Amendment 11.
        inner_sp = Float64[]
        for j in 1:n, i in 1:n
            s = SURF.spread[i, j]
            isnan(s) && continue
            abs(ComplexF64(SURF.xs[i], SURF.ys[j])) ≤ 3.0 &&
                push!(inner_sp, s)
        end
        @test !isempty(inner_sp)
        sort!(inner_sp)
        inner_med = inner_sp[cld(length(inner_sp), 2)]
        inner_max = inner_sp[end]
        @test inner_med < 5e-4               # C2 — well below worklog ~2.9e-3
        # The inner-arc MAX is the asymptotic-seed truncation floor
        # (Amendment 11) — honestly bounded by the measured ~6.7·10⁻³
        # plus headroom.  This is NOT the < 3.9·10⁻³ Amendment 10
        # anticipated: the C1 probe proved that target rested on a
        # falsified premise (a Laplace under-resolution); the genuine
        # floor is the seed and the `401²` grid exposes it honestly.
        @test inner_max < 8e-3

        # The seed floor is r-LOCALISED — a boundary layer at the inner
        # arc that decays steeply outward.  Verify the spread genuinely
        # belongs to the inner annulus, not a global defect: the spread
        # in the outer half of the inner band (`2.6 ≤ |x| ≤ 3`) is
        # materially smaller than at the arc itself.  If the inner-arc
        # max were a global Laplace bug it would not decay with radius —
        # this assertion is what distinguishes the seed-floor diagnosis
        # from the (falsified) under-resolution one.
        outer_band = Float64[]
        for j in 1:n, i in 1:n
            s = SURF.spread[i, j]
            isnan(s) && continue
            r = abs(ComplexF64(SURF.xs[i], SURF.ys[j]))
            (2.6 ≤ r ≤ 3.0) && push!(outer_band, s)
        end
        @test !isempty(outer_band)
        @test maximum(outer_band) < inner_max   # the floor decays outward
        @info "PI2S.2 — C1 inner-arc spread (|x| ≤ 3, ADR-0025 Amend. 11)" *
              "\n  inner-arc spread median = $(round(inner_med; sigdigits = 4))" *
              " (worklog-057 baseline ≈ 2.9e-3)" *
              "\n  inner-arc spread MAX    = $(round(inner_max; sigdigits = 4))" *
              " (asymptotic-seed truncation floor — Amendment 11)" *
              "\n  |x|∈[2.6,3.0] max       = $(round(maximum(outer_band); sigdigits = 4))" *
              " (the floor decays outward — it is a seed boundary layer)" *
              "\n  whole-sector spread median = $(round(med; sigdigits = 4))"

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

    @testset "PI2S.7 — VC-7 loop-closure quality certificate" begin
        # ADR-0025 VC-7 (Amendment 3 §"D-VC7 scope" + Amendment 7; bead
        # `padetaylor-0ln.37.14`).  `quality_diagnose` (ADR-0016) is the
        # loop-closure quality certificate: on every non-tree Delaunay
        # edge of the visited-node cloud it records the midpoint
        # disagreement `ΔP_rel` of the two endpoints' canonical Padé
        # approximants, classifies each edge, and aggregates.  D-VC7
        # adds the **vector adapter** — `quality_diagnose` on the
        # `VectorPathNetworkSolution` the headline figure's wedge walk
        # produces — and this testset is the VC-7 figure acceptance gate.
        #
        # WHAT THE METRIC ACTUALLY MEASURES (ADR-0025 Amendment 7).
        # `quality_diagnose` evaluates each node's CANONICAL shared-`Q`
        # Padé at Delaunay-edge midpoints.  That metric is governed by
        # the geometry `Delaunay-edge-length vs canonical step h` — it
        # is INDEPENDENT of the Stage-2 `extrapolate` flag and the B1
        # true-radius gate (those gate only the separate Stage-2 grid
        # fill).  So VC-7 is a structural quality CERTIFICATE of the
        # walk's loop-closure consistency; it is *not* a before/after
        # proxy for the B1/B2 re-resolution (Amendment 7 records the
        # measured numbers and corrects the A4 §2b over-attribution).
        #
        # The gate therefore asserts the certificate's genuine,
        # load-bearing invariants: the certificate is non-degenerate and
        # structurally sound, the well-closed lobe genuinely closes to
        # machine precision (proving the shared-`Q` `Pᵢ(t)/Q(t)`
        # evaluation pipeline is correct end-to-end), and a meaningful
        # fraction of the genuine in-disc loop closures close tightly.
        # The load-bearing assertions are mutation-proven below.
        walk = SURF.wedge_walk
        @test walk !== nothing

        rep = quality_diagnose(walk)

        # --- (a) the certificate is non-degenerate -------------------------
        # The B3 threading fan builds an ~1800-node walk; its Delaunay
        # graph has thousands of non-tree edges — a rich loop-closure
        # population.  `> 100` is a generous floor (the measured count
        # is ~3700).
        @test rep.n_edges > 100
        # Single-sheeted companion system ⇒ the report's `sheet` is 0.
        @test rep.sheet == 0

        # --- (b) structural soundness — the scalar-method invariants -------
        # The five category tallies partition `n_edges` (Diagnostics.jl
        # contract); the quantiles are monotone non-decreasing.  This is
        # the same shape `test/diagnose_test.jl` asserts for the scalar
        # method — the vector adapter must honour it identically.
        @test rep.n_well_closed + rep.n_noisy + rep.n_extrap_driven +
              rep.n_depth_driven + rep.n_branch_cut == rep.n_edges
        @test rep.n_branch_cut == 0          # reserved — v1
        @test rep.median_ΔP_rel ≤ rep.p90_ΔP_rel
        @test rep.p90_ΔP_rel    ≤ rep.p99_ΔP_rel
        @test rep.p99_ΔP_rel    ≤ rep.max_ΔP_rel
        @test issorted([e.ΔP_rel for e in rep.worst_edges]; rev = true)

        # --- (c) the well-closed lobe genuinely closes ---------------------
        # The LOAD-BEARING positive invariant.  Among all non-tree
        # edges, the best-closed one must close to machine precision
        # (`min ΔP_rel < 1e-8`): two independently-walked shared-`Q`
        # patches agree at their shared midpoint to ~1e-10.  This proves
        # the vector adapter's `Pᵢ(t)/Q(t)` Horner evaluation, the
        # midpoint rescaling `t = (M − z_X)/h_X`, and the vector-norm
        # `ΔP_rel` are all correct end-to-end — a corrupted node's
        # numerators or denominator would destroy this lobe (the
        # mutation-proof below corrupts exactly that and confirms RED).
        all_rels = Float64[e.ΔP_rel for e in
                           quality_diagnose(walk; n_worst = 10^9).worst_edges]
        @test minimum(all_rels) < 1e-8

        # --- (d) the in-disc loop closures — the honest consensus facet ----
        # An edge is a genuine loop-closure test only when its midpoint
        # is INSIDE the canonical disc of both endpoints (`extrap_max ≤
        # 1`); the rest extrapolate and the disagreement is dominated by
        # Padé extrapolation, not loop-closure error (Diagnostics.jl
        # module docstring).  On the figure walk a meaningful fraction
        # of the in-disc edges close tightly (`ΔP_rel ≤ tol_bad`); the
        # measured fraction is ~0.38.  `> 0.25` is the load-bearing bar
        # — comfortably passed, and mutation-sensitive (corrupting a
        # node collapses both this fraction and the lobe in (c)).
        all_edges = quality_diagnose(walk; n_worst = 10^9).worst_edges
        in_disc   = [e for e in all_edges if e.extrap_max ≤ 1.0]
        @test !isempty(in_disc)
        n_tight   = count(e -> e.ΔP_rel ≤ rep.tol_bad, in_disc)
        frac_tight = n_tight / length(in_disc)
        @test frac_tight > 0.25

        # --- measured distribution (Amendment 7 — reported, not masked) ----
        # Honest reporting: VC-7's value is the certificate, and the
        # certificate is only meaningful if its numbers are surfaced.
        @info "PI2S.7 — VC-7 loop-closure certificate (figure wedge walk)" *
              "\n  visited nodes        = $(length(walk.visited_z))" *
              "\n  non-tree edges       = $(rep.n_edges)" *
              "\n  ΔP_rel  median       = $(round(rep.median_ΔP_rel; sigdigits = 4))" *
              "\n  ΔP_rel  p90 / p99    = $(round(rep.p90_ΔP_rel; sigdigits = 4)) / " *
                                          "$(round(rep.p99_ΔP_rel; sigdigits = 4))" *
              "\n  well/noisy/extr/depth= $(rep.n_well_closed) / $(rep.n_noisy) / " *
                  "$(rep.n_extrap_driven) / $(rep.n_depth_driven)" *
              "\n  in-disc edges        = $(length(in_disc)) " *
                  "($(round(100 * length(in_disc) / rep.n_edges; digits = 1))%)" *
              "\n  in-disc tight frac   = $(round(frac_tight; digits = 3))" *
              "\n  best-closed ΔP_rel   = $(round(minimum(all_rels); sigdigits = 3))"

        # --- VC-7 mutation-proof — a corrupted node breaks the certificate -
        # The load-bearing positive invariant is (c) `min ΔP_rel < 1e-8`
        # — the well-closed lobe.  Corrupt **one endpoint** of the
        # best-closed edge: add a fixed real offset to that node's
        # constant numerator coefficient.  The corruption is
        # deliberately ASYMMETRIC and ADDITIVE — `ΔP_rel` is a *relative*
        # norm ratio, so scaling *both* endpoints (or scaling one
        # multiplicatively while the other stays) would leave a
        # both-poisoned edge invariant; an additive shift to a single
        # endpoint genuinely moves `y_A` away from `y_B`.  With the
        # best-closed edge's endpoint A poisoned, that edge can no
        # longer close to machine precision: its `ΔP_rel` in the
        # corrupted walk jumps above `tol_well`.  This proves invariant
        # (c) has teeth — a bug corrupting a node's shared-`Q` numerator
        # is caught by VC-7, not silently passed.
        be_i      = argmin([e.ΔP_rel for e in all_edges])
        best_edge = all_edges[be_i]
        @test best_edge.ΔP_rel ≤ rep.tol_well        # the lobe edge itself
        let vnum = [[copy(p) for p in node]
                    for node in walk.visited_numerators]
            # Poison endpoint A's first numerator's constant term.
            vnum[best_edge.A][1][1] += 1.0
            WT  = typeof(walk)
            bad = WT(walk.visited_z, walk.visited_y, walk.visited_h,
                     vnum, walk.visited_denominator, walk.visited_parent,
                     walk.grid_z, walk.grid_y, walk.visited_jets)
            bad_edges = quality_diagnose(bad; n_worst = 10^9).worst_edges
            # Re-find the SAME (A,B) edge in the corrupted walk's report.
            ai, bi = best_edge.A, best_edge.B
            poisoned_edge = nothing
            for e in bad_edges
                if (e.A == ai && e.B == bi) || (e.A == bi && e.B == ai)
                    poisoned_edge = e; break
                end
            end
            @test poisoned_edge !== nothing
            # The corrupted edge no longer closes to machine precision —
            # invariant (c) `min ΔP_rel ≤ tol_well` genuinely fails for
            # this edge on the corrupted walk (RED), so it is genuinely
            # load-bearing.
            @test poisoned_edge.ΔP_rel > rep.tol_well
            @test poisoned_edge.ΔP_rel > best_edge.ΔP_rel
        end
    end

    @testset "PI2S.8 — VC-8 BVP endpoint higher-derivative match" begin
        # ADR-0025 VC-8 (Validation Criteria Menu; bead
        # `padetaylor-0ln.37.15`; Amendment 8).  FW 2011 §5.2's
        # endpoint-derivative diagnostic
        # (`references/markdown/FW2011_painleve_methodology_JCP230/
        # FW2011_painleve_methodology_JCP230.md:192-193`) — estimate the
        # endpoint derivatives of a converged Chebyshev BVP solution via
        # the `D₁` matrix and match them to known values; a mismatch
        # beyond `~10⁻⁷·10⁻⁸` says the collocation `N` is too small —
        # adapted to the figure's `d = 4` first-order vector **companion**
        # BVP.  The companion form asserts the exact chain identities
        # `y[k]' = y[k+1]` (`k = 1,2,3`); `surf_vc8_companion_check`
        # certifies that a converged collocation solution satisfies them:
        # the spectral derivative `D₁·y[k]/s` must reproduce `y[k+1]` at
        # the endpoints (and throughout) to FW tolerance.  This is a
        # genuine self-consistency invariant of the spectral solution —
        # it owes nothing to the asymptotic seed and catches an
        # under-resolved `N`.
        #
        # The figure solves BVPs in `surf_ray_bvp` (every ray of the
        # Region-1 fan) and `surf_anchor_bvp` (the Region-2 wedge
        # anchor).  VC-8 is asserted on a representative ray BVP and on
        # the anchor BVP — the figure's two BVP code paths.
        f8, Jf8 = surf_companion()

        # The FW §5.2 tolerance band: the diagnostic's "increase N"
        # threshold is `10⁻⁷·10⁻⁸`.  A converged figure BVP must sit
        # comfortably below it — `1e-7` is the load-bearing bar.
        FW_VC8_TOL = 1.0e-7

        for (label, sol8) in (
                ("ray φ=180° (negative real axis)",
                 surf_ray_bvp(f8, Jf8, deg2rad(180.0))),
                ("ray φ=120°", surf_ray_bvp(f8, Jf8, deg2rad(120.0))),
                ("anchor BVP [-20,-2]", surf_anchor_bvp()))
            d8 = surf_vc8_companion_check(sol8)

            # (a) the companion-consistency identity holds to FW
            # tolerance — both at the endpoints (FW §5.2's check) and
            # over the whole node set.  A collocation `N` too small
            # would leave `D₁·y[k]/s − y[k+1]` above this band.
            @test d8.consistency_endpoint < FW_VC8_TOL
            @test d8.consistency_max      < FW_VC8_TOL
            # The endpoint mismatch cannot exceed the global max (the
            # endpoints are a subset of the nodes) — a structural sanity
            # bound that pins the two quantities are consistently
            # computed.
            @test d8.consistency_endpoint ≤ d8.consistency_max + 1e-15

            # (b) the BVP's own Newton residual is at the spectral floor
            # — VC-8 reads a *converged* solution.
            @test d8.residual_inf < 1.0e-8

            # (c) the asymptotic-IC cross-check.  `VectorBVP` pins only
            # `y[1],y[2]`; the free `y[3]=u''`, `y[4]=u'''` at the
            # deep-asymptotic endpoint `|z_a|=20` must agree with the
            # `pI2_tritronquee_ic` seed to the seed's own `n_terms=2`
            # truncation accuracy — `O(10⁻⁵)` at `|x|=20` (the V8b probe
            # number).  This certifies the BVP converged to the
            # *tritronquée branch*, not merely to a companion solution.
            @test d8.ic_y3_err < 1.0e-4
            @test d8.ic_y4_err < 1.0e-4
            # ...and to FAR better than O(1) — the seed and the solve
            # genuinely agree, they are not unrelated.
            @test d8.ic_y3_err < 1.0e-3
            @test d8.ic_y4_err < 1.0e-3

            @info "PI2S.8 — VC-8 BVP companion-consistency ($label)" *
                  "\n  N                       = $(d8.N)" *
                  "\n  companion-consistency max = " *
                      "$(round(d8.consistency_max; sigdigits = 4))" *
                  "\n  endpoint consistency      = " *
                      "$(round(d8.consistency_endpoint; sigdigits = 4))" *
                  "\n  y₃ vs asymptotic seed     = " *
                      "$(round(d8.ic_y3_err; sigdigits = 4))" *
                  "\n  y₄ vs asymptotic seed     = " *
                      "$(round(d8.ic_y4_err; sigdigits = 4))" *
                  "\n  BVP Newton residual       = " *
                      "$(round(d8.residual_inf; sigdigits = 4))"
        end

        # --- VC-8 mutation-proof — an under-resolved BVP fails VC-8 -------
        # The load-bearing invariant is `consistency_max < FW_VC8_TOL`.
        # Solve the SAME ray BVP at a deliberately starved collocation
        # `N` — small enough that the Chebyshev interpolant cannot
        # resolve the tritronquée's curvature on the `[2,20]` ray — and
        # confirm the companion-consistency residual blows past the FW
        # band.  This proves VC-8 is the genuine "increase N" detector
        # FW §5.2 describes, not a check the figure passes for free: at
        # the figure `N = SURF_BVP_N = 96` it is comfortably GREEN, at a
        # starved `N` it is RED.  (`N = 6` is the smallest the τ-method
        # admits a 4-vector BC with — `VectorBVP` floors at `N ≥ 4`.)
        starved = vector_bvp_solve(
            f8,
            SURF_R_MAX * cis(deg2rad(180.0)),
            SURF_R_MIN * cis(deg2rad(180.0)),
            ComplexF64[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0],
            ComplexF64[0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0],
            let sa = pI2_tritronquee_ic(SURF_R_MAX * cis(deg2rad(180.0));
                                        t = SURF_T, n_terms = 2),
                sb = pI2_tritronquee_ic(SURF_R_MIN * cis(deg2rad(180.0));
                                        t = SURF_T, n_terms = 2)
                ComplexF64[sa[1], sa[2], sb[1], sb[2]]
            end;
            N = 6, tol = SURF_BVP_TOL, maxiter = SURF_BVP_MAXITER,
            jacobian = Jf8,
            initial_guess = z -> pI2_tritronquee_ic(z; t = SURF_T,
                                                    n_terms = 2))
        d_starved = surf_vc8_companion_check(starved)
        # The starved solve's companion-consistency residual is far
        # above the FW band — VC-8 RED, the "increase N" signal fires.
        @test d_starved.consistency_max > FW_VC8_TOL
        # ...and a converged figure-`N` solve on the same ray is GREEN —
        # the discriminator is genuine.
        d_full = surf_vc8_companion_check(
            surf_ray_bvp(f8, Jf8, deg2rad(180.0)))
        @test d_full.consistency_max < FW_VC8_TOL
        @test d_starved.consistency_max > d_full.consistency_max
        @info "PI2S.8 — VC-8 mutation-proof (starved N=6 vs figure N)" *
              "\n  starved N=6  companion-consistency max = " *
                  "$(round(d_starved.consistency_max; sigdigits = 4))" *
              "\n  figure  N=$(d_full.N) companion-consistency max = " *
                  "$(round(d_full.consistency_max; sigdigits = 4))"
    end

    @testset "PI2S.10 — VC-10 two-run pole-field accuracy indicator" begin
        # ADR-0025 VC-10 (Validation Criteria Menu; Amendment 9).  The
        # FW/FFW double-run accuracy indicator (FW 2011 §3.1 line 156
        # shuffles the Stage-1 targets; FFW 2017
        # `references/markdown/FFW2017_painleve_riemann_surfaces_
        # preprint/FFW2017_painleve_riemann_surfaces_preprint.md:246-247`
        # — "if we compute the same solution twice, different paths will
        # be run, resulting in solutions that should differ by
        # approximately the numerical error").  The wedge walk is run
        # twice with two different `MersenneTwister` target orderings;
        # each run independently extracts + VC-4-validates a pole field;
        # the two fields are matched nearest-neighbour and the
        # matched-pole disagreement is the practical, oracle-free
        # accuracy estimate of the wedge pole field.
        #
        # The wedge pole field has NO external oracle (KKG published no
        # pole table — ADR-0025 Decision 1), so this manufactured
        # FW-style criterion is the figure's accuracy certificate.
        anchor8 = surf_anchor_bvp()

        # The wedge grid points — exactly the set `kkg_pi2_surface`
        # builds for the Stage-2 fill (inside the disc, in the wedge,
        # off the masked strip, not the origin).
        wedge_pts10 = ComplexF64[]
        for j in 1:n, i in 1:n
            z = ComplexF64(SURF.xs[i], SURF.ys[j])
            abs(z) > SURF_XY_LIM && continue
            surf_in_mask(z)      && continue
            surf_in_sector(z)    && continue
            iszero(z)            && continue
            push!(wedge_pts10, z)
        end
        @test !isempty(wedge_pts10)

        vc10 = surf_vc10_two_run(anchor8, wedge_pts10)

        # --- (a) each run produced a genuine, VC-4-validated field -------
        # Both differently-ordered runs thread the whole wedge and
        # extract a rich field (richer than V8b's 21 poles — the B3 fan
        # bar).  A run that collapsed (empty / tiny field) would make the
        # indicator meaningless.
        @test length(vc10.poles_a) > 21
        @test length(vc10.poles_b) > 21
        @test all(isfinite, vc10.poles_a)
        @test all(isfinite, vc10.poles_b)

        m = vc10.match
        # The two fields DO pair up — most poles are found by both runs.
        @test !isempty(m.pairs)
        @test length(m.pairs) ≥ 100

        # --- (b) the load-bearing accuracy bar ---------------------------
        # The matched-pole disagreement IS the FW/FFW accuracy estimate.
        # The measured value (ADR-0025 Amendment 9) is median ≈ 0.35,
        # max ≈ 0.59.  Assert the median is genuinely SMALL — well below
        # the wedge pole nearest-neighbour spacing (~0.69, the VC-5
        # measurement).  `< 0.45` is the load-bearing envelope: below the
        # spacing (so the disagreement is genuine same-pole jitter, not
        # confusion between distinct lattice poles) and a touch above the
        # measured 0.35.  This is the criterion's deliverable — the wedge
        # pole field is accurate to ~half the pole spacing.
        @test m.median_disagree < 0.45
        # The disagreement is genuinely better than the pole spacing —
        # the accuracy indicator is meaningful, not vacuous.
        @test m.median_disagree < 0.69
        @test isfinite(m.max_disagree)
        # Every matched pair is within the match tolerance — a genuine
        # same-physical-pole pairing, not a cross-lattice mis-match.
        for (pa, pb) in m.pairs
            @test abs(pa - pb) ≤ SURF_VC10_MATCH_TOL
        end
        # The disagreements list aligns with the pairs.
        @test length(m.disagreements) == length(m.pairs)
        # The reported median / max equal the list's median / max — the
        # NamedTuple is internally consistent (independent re-derivation).
        ds = sort(m.disagreements)
        @test m.median_disagree ≈ ds[cld(length(ds), 2)]
        @test m.max_disagree ≈ ds[end]

        # --- (c) the single-run-only poles -------------------------------
        # Per FFW, the poles found by only ONE run are a second facet of
        # the indicator: the far-wedge A2 tractability ceiling means each
        # ordering resolves a slightly different subset of the dense pole
        # lattice.  The count is reported (not masked) and bounded — a
        # MAJORITY of each field must still be matched (the two runs
        # genuinely agree, they do not produce disjoint fields).
        @test length(m.only_a) < length(vc10.poles_a) / 2
        @test length(m.only_b) < length(vc10.poles_b) / 2
        # The matched + only-this-run poles partition each field.
        @test length(m.pairs) + length(m.only_a) == length(vc10.poles_a)
        @test length(m.pairs) + length(m.only_b) == length(vc10.poles_b)

        # --- (d) per-fixed-seed determinism — VC-10 ≠ a broken PI2S.6 ----
        # VC-10 deliberately introduces a DIFFERENT ordering per seed,
        # but each individual seed's walk stays fully deterministic (the
        # reproducibility guarantee `test/kkg_pi2_figure_test.jl`
        # PI2F.1.6 pins, here lifted to the seeded vector walk).  Re-run
        # one seed and confirm bit-identical visited tree + pole field.
        run_a1 = surf_wedge_fill(anchor8, wedge_pts10;
                                 rng = MersenneTwister(SURF_VC10_SEED_A))
        run_a2 = surf_wedge_fill(anchor8, wedge_pts10;
                                 rng = MersenneTwister(SURF_VC10_SEED_A))
        @test run_a1.walk.visited_z == run_a2.walk.visited_z
        @test run_a1.walk.visited_parent == run_a2.walk.visited_parent
        @test run_a1.poles == run_a2.poles
        # ...and a DIFFERENT seed genuinely builds a different tree — the
        # two-run indicator is exercising real path divergence, not the
        # same walk twice.
        run_b1 = surf_wedge_fill(anchor8, wedge_pts10;
                                 rng = MersenneTwister(SURF_VC10_SEED_B))
        @test run_a1.walk.visited_z != run_b1.walk.visited_z

        @info "PI2S.10 — VC-10 two-run pole-field accuracy indicator" *
              "\n  $(vc10.message)" *
              "\n  matched pairs        = $(length(m.pairs))" *
              "\n  disagreement median  = " *
                  "$(round(m.median_disagree; sigdigits = 4))" *
              "\n  disagreement max     = " *
                  "$(round(m.max_disagree; sigdigits = 4))" *
              "\n  single-run-only      = $(length(m.only_a)) (A) + " *
                  "$(length(m.only_b)) (B)"

        # --- VC-10 mutation-proof — a same-ordering run has zero error ---
        # The load-bearing fact is that the disagreement is a genuine
        # error indicator — driven by the DIFFERENT orderings, not
        # spurious.  Run the wedge walk twice with the SAME seed and
        # confirm the two pole fields are IDENTICAL (zero disagreement,
        # every pole matched, no single-run-only pole).  If
        # `surf_vc10_two_run` were measuring something other than
        # ordering-induced path divergence, a same-seed pair would still
        # show disagreement — it does not.  Conversely the genuine
        # two-seed run above shows a non-zero median: the indicator has
        # teeth precisely because different orderings diverge.
        same = surf_vc10_match(run_a1.poles, run_a1.poles)
        @test length(same.pairs) == length(run_a1.poles)
        @test isempty(same.only_a) && isempty(same.only_b)
        @test same.max_disagree == 0.0
        # The genuine two-ordering indicator is strictly non-zero — the
        # different paths really did diverge (the FFW premise).
        @test m.median_disagree > 0.0
    end

    @testset "PI2S.11 — C2/C3 sector voter-1 reconstruction" begin
        # ADR-0025 Amendment 10 (beads `padetaylor-0ln.37.10` C2 +
        # `0ln.37.11` C3).  C2 retires v1 corner C1 — the inner-arc
        # Dirichlet datum is now a BVP-interior sample, not the KKG
        # asymptotic series.  C3 retires v1 corner C5 — voter 1 is
        # reconstructed exact-in-radius (the stored ray `VectorBVP
        # Solution`s, barycentric) + cubic-in-angle (Catmull-Rom across
        # four bracketing rays), retiring the lossy bilinear polar
        # raster.  A *harmonic* voter 1 was rejected by the C3 audition:
        # it would be a Laplace solve `≡ voter 2` and collapse the
        # ADR-0024 triple-method independence.  This testset pins both
        # decisions' load-bearing invariants and mutation-proves them.
        fan11 = surf_ray_fan()
        f11, Jf11 = surf_companion()

        # --- C2 — the ray BVP segment extends inward past the inner arc --
        # The genuine C2 fix is the *extended-inward* segment: each ray
        # BVP runs `[SURF_R_MAX, SURF_R_INNER_BC]` with
        # `SURF_R_INNER_BC < SURF_R_MIN`, so the inner arc `|x| = R_MIN`
        # is an INTERIOR collocation point — not the pinned endpoint.
        @test SURF_R_INNER_BC < SURF_R_MIN
        ray11 = surf_ray_bvp(f11, Jf11, deg2rad(196.0))
        @test isapprox(abs(ray11.z_b), SURF_R_INNER_BC; atol = 1e-9)
        @test isapprox(abs(ray11.z_a), SURF_R_MAX;      atol = 1e-9)
        # The inner arc is strictly inside the BVP segment — the C2
        # precondition that makes the inner-arc datum an ODE-constrained
        # interior value rather than the asymptotic seed.
        @test SURF_R_INNER_BC < SURF_R_MIN < SURF_R_MAX

        # --- C2 — the inner-arc datum beats the asymptotic seed ----------
        # Independent oracle: a high-`N`, deep-inward dedicated ray BVP,
        # converged to ~1e-14.  The C2 fan inner-arc datum must be
        # materially closer to it than the `n_terms = 2` series is.
        function oracle_ray(φ)
            CT = ComplexF64
            z_a = SURF_R_MAX * cis(φ); z_b = 1.02 * cis(φ)
            Ba = zeros(CT, 4, 4); Bb = zeros(CT, 4, 4)
            Ba[1,1] = 1; Ba[2,2] = 1; Bb[3,1] = 1; Bb[4,2] = 1
            sa = pI2_tritronquee_ic(z_a; t = SURF_T, n_terms = 2)
            sb = pI2_tritronquee_ic(z_b; t = SURF_T, n_terms = 2)
            vector_bvp_solve(f11, z_a, z_b, Ba, Bb,
                             CT[sa[1], sa[2], sb[1], sb[2]];
                             N = 220, tol = 1e-9, maxiter = 40,
                             jacobian = Jf11,
                             initial_guess = z -> pI2_tritronquee_ic(z;
                                 t = SURF_T, n_terms = 2))
        end
        n_c2_better = 0; n_c2_tested = 0
        for degφ in (60.0, 120.0, 196.0, 240.0, 300.0)
            φ = deg2rad(degφ)
            oref = oracle_ray(φ)(ComplexF64(SURF_R_MIN * cis(φ)))[1]
            fan_arc = surf_ray_arc_eval(fan11, SURF_R_MIN, φ)
            seed_arc = pI2_tritronquee_ic(ComplexF64(SURF_R_MIN * cis(φ));
                                          t = SURF_T, n_terms = 2)[1]
            err_fan  = abs(fan_arc  - oref)
            err_seed = abs(seed_arc - oref)
            # The C2 BVP-sourced inner-arc datum is closer to the oracle
            # than the asymptotic series — at every tested angle.
            @test err_fan < err_seed
            # ...and it is genuinely accurate (the C2 probe: ~1e-4).
            @test err_fan < 1e-3
            err_fan < err_seed && (n_c2_better += 1)
            n_c2_tested += 1
        end
        @test n_c2_better == n_c2_tested

        # --- C2 mutation-proof — the legacy `[R_MAX, R_MIN]` segment ------
        # Solve the SAME ray with the inner BC pinned at `SURF_R_MIN`
        # (the pre-C2 behaviour): the inner endpoint is then the pinned
        # asymptotic seed, and the value at `|x| = R_MIN` is the seed
        # *exactly* — measurably less accurate than the C2 interior
        # datum.  This proves the extended-inward segment is genuinely
        # load-bearing for C2, not cosmetic.
        legacy = surf_ray_bvp(f11, Jf11, deg2rad(196.0);
                              r_in = SURF_R_MIN)
        legacy_arc = legacy(ComplexF64(SURF_R_MIN * cis(deg2rad(196.0))))[1]
        seed196 = pI2_tritronquee_ic(
            ComplexF64(SURF_R_MIN * cis(deg2rad(196.0)));
            t = SURF_T, n_terms = 2)[1]
        # The legacy inner endpoint IS the asymptotic seed (the 2+2 split
        # BC pins it) — the exact equality the C2-exploration probe found.
        @test isapprox(legacy_arc, seed196; atol = 1e-12)
        oref196 = oracle_ray(deg2rad(196.0))(
            ComplexF64(SURF_R_MIN * cis(deg2rad(196.0))))[1]
        # C2's interior datum is strictly closer to the oracle than the
        # legacy pinned-seed endpoint — the mutation-proof discriminator.
        @test abs(surf_ray_arc_eval(fan11, SURF_R_MIN, deg2rad(196.0)) -
                  oref196) < abs(legacy_arc - oref196)

        # --- C3 — voter 1 is exact in the radius -------------------------
        # The C3 reconstruction evaluates the stored ray `VectorBVP
        # Solution` directly — a barycentric spectral interpolant, exact
        # at ANY radius.  On a fan ray, `surf_ray_eval` at a non-raster
        # radius must equal that ray's BVP solution to ~1e-9.  The
        # retired bilinear scheme interpolated a 40-row raster and could
        # not — it carried ~3.6e-4 radial error (the C3 audition).
        j_mid = cld(length(fan11.phis), 2)
        φ_on  = fan11.phis[j_mid]
        sol_on = fan11.sols[j_mid]
        c3_radial_max = 0.0
        for r in (3.3, 7.7, 11.1, 16.6)        # deliberately off-raster
            z = ComplexF64(r * cis(φ_on))
            c3_radial_max = max(c3_radial_max,
                                abs(surf_ray_eval(fan11, z) - sol_on(z)[1]))
        end
        # Exact in the radius — far below the retired bilinear scheme's
        # ~3.6e-4 on-ray radial error.
        @test c3_radial_max < 1e-8

        # --- C3 — voter 1 is accurate between rays (cubic-angular) -------
        # At an inter-ray midpoint a DEDICATED BVP straight through that
        # angle is the independent oracle.  The C3 cubic reconstruction's
        # inter-ray error must be materially below the retired bilinear
        # scheme's ~1.3e-3 max (the C3 audition measured ~6.9e-5).
        c3_inter_max = 0.0
        for jr in (10, 20, 30)
            φ_mid = (fan11.phis[jr] + fan11.phis[jr + 1]) / 2
            sol_mid = surf_ray_bvp(f11, Jf11, φ_mid)
            for r in (4.0, 9.0, 14.0)
                z = ComplexF64(r * cis(φ_mid))
                c3_inter_max = max(c3_inter_max,
                                   abs(surf_ray_eval(fan11, z) - sol_mid(z)[1]))
            end
        end
        # Materially below the retired bilinear scheme's ~1.3e-3 max —
        # the load-bearing C3 bar (`< 3e-4` is comfortably above the
        # audition's ~6.9e-5 yet far below bilinear's 1.3e-3).
        @test c3_inter_max < 3e-4

        # --- C3 mutation-proof — the retired bilinear polar raster -------
        # Reconstruct voter 1 the OLD way: sample the fan onto a 40-row
        # polar raster and bilinearly interpolate.  On a fan ray at an
        # off-raster radius it carries a real radial interpolation error
        # — it CANNOT meet the C3 exact-radius bar.  This proves the
        # exact-radius reconstruction is genuinely load-bearing.
        let radii_b = collect(range(SURF_R_MIN, SURF_R_MAX; length = 40)),
            phis_b  = fan11.phis
            Ub = Matrix{ComplexF64}(undef, length(radii_b), length(phis_b))
            for (jb, sol) in enumerate(fan11.sols), (ib, r) in
                    enumerate(radii_b)
                Ub[ib, jb] = sol(ComplexF64(r * cis(phis_b[jb])))[1]
            end
            function bilin_v1(z)
                r = abs(z); φ = mod2pi(angle(z))
                ib = clamp(searchsortedlast(radii_b, r), 1,
                           length(radii_b) - 1)
                jb = clamp(searchsortedlast(phis_b, φ), 1,
                           length(phis_b) - 1)
                tr = (r - radii_b[ib]) / (radii_b[ib+1] - radii_b[ib])
                tφ = (φ - phis_b[jb]) / (phis_b[jb+1] - phis_b[jb])
                return (1-tr)*(1-tφ)*Ub[ib,jb] + tr*(1-tφ)*Ub[ib+1,jb] +
                       (1-tr)*tφ*Ub[ib,jb+1]   + tr*tφ*Ub[ib+1,jb+1]
            end
            bilin_radial_max = 0.0
            for r in (3.3, 7.7, 11.1, 16.6)
                z = ComplexF64(r * cis(φ_on))
                bilin_radial_max = max(bilin_radial_max,
                                       abs(bilin_v1(z) - sol_on(z)[1]))
            end
            # The retired bilinear scheme carries a real on-ray radial
            # error (~1e-4) — far above the C3 exact-radius `1e-8` bar:
            # the reconstruction change is genuinely load-bearing.
            @test bilin_radial_max > 1e-6
            @test bilin_radial_max > c3_radial_max
        end

        # --- C3 — the three voters stay algorithmically independent ------
        # The audition's binding constraint: voter 1 must NOT become a
        # Laplace solve (≡ voter 2).  Verify the three voters genuinely
        # disagree at a mid-sector point — a FEM voter, a spectral voter
        # and a direct-ODE voter, three disjoint computations.  If voter
        # 1 had been made harmonic, voter 1 ≡ voter 2 to ~1e-12 here.
        lap11 = surf_laplace_voters(fan11)
        z_ind = ComplexF64(-8.0, 3.0)
        v1_ind = surf_ray_eval(fan11, z_ind)
        l2_ind = surf_laplace_eval_one(lap11.rem2, lap11.imm2,
                                       lap11.rect, z_ind)
        l3_ind = surf_laplace_eval_one(lap11.rem3, lap11.imm3,
                                       lap11.rect, z_ind)
        @test l2_ind !== nothing && l3_ind !== nothing
        # Voter 1 (direct ODE solve) is NOT bit-identical to voter 2
        # (spectral Laplace) — they are genuinely different computations.
        @test abs(real(v1_ind) - l2_ind.re) > 1e-9
        # ...yet all three agree to figure tolerance (the vote is sound).
        @test abs(real(v1_ind) - l2_ind.re) < 1e-2
        @test abs(l2_ind.re - l3_ind.re)    < 1e-2

        @info "PI2S.11 — C2/C3 sector voter-1 reconstruction" *
              "\n  C2 inner-arc datum vs oracle (5 angles): all beat the seed" *
              "\n  C3 on-ray radial error  = $(round(c3_radial_max; sigdigits = 3))" *
              " (retired bilinear: ~3.6e-4)" *
              "\n  C3 inter-ray error      = $(round(c3_inter_max; sigdigits = 3))" *
              " (retired bilinear: ~1.3e-3)"
    end
end
