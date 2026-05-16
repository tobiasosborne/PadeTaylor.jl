# test/ffw_fig_3_test.jl
#
# Quantitative pin for `figures/ffw2017_fig_3.jl` — the PVI phase-portrait
# reveal on the ζ + z Riemann surfaces.  Same PVI problem as Fig 2
# (FFW md:195 verbatim parameters + IC), same A4 cross-mode + A5
# sheet-aware Stage-2 stack, but a different load-bearing claim:
# Fig 3 reveals the **multi-sheet structure** via phase portraits,
# so the pinning is about the sheet bookkeeping reaching the
# expected per-sheet values, not just visiting non-[0] sheets.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:178-199
#         (caption + sheet parametrisation).
#
# ## Verification strategy
#
# Sheet [0]: a tight-tol absolute pin at a well-walked sample.  Both
# loose (`adaptive_tol = 1e-10`, the figure-script setting) and tight
# (`1e-12`) runs converge to the same value at this sample (tronquée
# pole-free region) — the absolute pin against the tight value is the
# load-bearing mutation-detector (worklog 044 §"two-pronged").
#
# Sheet [1]: the load-bearing Fig-3 claim is that sheet bookkeeping
# is active and the alternate-sheet evaluator returns a DIFFERENT
# value than the principal-sheet evaluator at the same ζ point.
# Without A4+A5's sheet bookkeeping the alternate-sheet query would
# either NaN (no sheet-[1] nodes visited) or return the principal-
# sheet value — either way the test bites.
#
# ## Assertions (FF3.1.1 - FF3.1.6)
#
#   - FF3.1.1 IC round-trip exact via `pV_z_to_ζ`/`pV_ζ_to_z` (FFW md:195
#             numerics).
#   - FF3.1.2 Sample points finite on both sheets [0] and [1].
#   - FF3.1.3 ζ-plane absolute pin (sheet 0) — tight-tol reference value;
#             loose-vs-tight degenerate to 0.0 in the well-walked region.
#   - FF3.1.4 Conjugate-symmetry pin: `arg w(ζ̄) + arg w(ζ) ≈ 0` for
#             real PVI parameters + real IC, on sheet 0 (Schwarz
#             reflection of analytic w(ζ) across the real axis).
#   - FF3.1.5 Sheet-population pin: `sol.visited_sheet` contains a non-
#             [0] sheet with > 5 visited nodes (cross-mode populates
#             both sides reproducibly under our seed=2 + winding ring).
#   - FF3.1.6 Multi-sheet phase pin: at a fixed ζ point near the
#             walker's cut crossing, `arg w` queried on sheets [0] vs
#             [1] differs by a measurable amount.  Without sheet
#             bookkeeping (M4) both queries return the same value
#             (or NaN); with sheet bookkeeping they're distinct.
#
# Target wall ≤ 12 s on a developer laptop.

using Test
using PadeTaylor
using PadeTaylor.SheetTracker:   pVI_transformed_rhs
using PadeTaylor.CoordTransforms: pV_z_to_ζ, pV_ζ_to_z

@testset "FFW 2017 Fig 3 — PVI phase-portrait Riemann surface acceptance" begin

    # ---- FFW md:195 parameters + IC (verbatim, same as Fig 2) ----------
    α, β, γ, δ        = 4.0, -4.0, 8.0, -8.0
    z₀, u₀, up₀       = 10.0 + 0.0im, 0.429534600325223 + 0.0im,
                          -1.61713114374804e-3 + 0.0im
    ζ₀, w₀, wp₀       = pV_z_to_ζ(z₀, u₀, up₀)

    # ---- FF3.1.1 IC round-trip exact -----------------------------------
    @testset "FF3.1.1: IC round-trip exact (z ↔ ζ)" begin
        z_back, u_back, up_back = pV_ζ_to_z(ζ₀, w₀, wp₀)
        @test isapprox(z_back, z₀; atol = 1e-13)
        @test u_back == u₀
        @test isapprox(up_back, up₀; atol = 1e-15)
        # FFW md:195 numerics literal pin (catches typo in IC).
        @test real(u₀) == 0.429534600325223
        @test real(up₀) == -1.61713114374804e-3
    end

    # ---- Build the figure-script Stage-1 target set (sparse base +
    #      winding ring) and solve at both tolerances ------------------
    function rect_targets(re_lo, re_hi, im_lo, im_hi, dx, dy)
        tgs = ComplexF64[]
        re = re_lo
        while re ≤ re_hi + 1e-12
            im_ = im_lo
            while im_ ≤ im_hi + 1e-12
                push!(tgs, complex(re, im_))
                im_ += dy
            end
            re += dx
        end
        return tgs
    end

    ζ_RE_LO, ζ_RE_HI = -1.5, 2.5
    ζ_IM_LO, ζ_IM_HI = -π + 0.05, 3π - 0.05
    targets_base = rect_targets(0.05, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI, 0.35, 0.7)
    winding_ring = ComplexF64[
        -0.5 + 0.5im, -0.5 - 0.5im,
        -1.0 + 0.5im, -1.0 - 0.5im,
        -1.0 + 1.0im, -1.0 - 1.0im,
        -0.5 + 1.0im, -0.5 - 1.0im,
    ]
    stage1_targets = vcat(targets_base, winding_ring)

    f = pVI_transformed_rhs(α, β, γ, δ)
    prob = PadeTaylorProblem(f, (w₀, wp₀), (ζ₀, 3.0 + 3im); order = 30)

    function solve_at(tol::Real)
        return path_network_solve(prob, stage1_targets;
                                   h = 0.3,
                                   step_size_policy = :adaptive_ffw,
                                   adaptive_tol = tol,
                                   k_conservative = 1e-3,
                                   max_rescales = 50,
                                   max_steps_per_target = 4000,
                                   branch_points     = (0.0+0.0im,),
                                   branch_cut_angles = π,
                                   cross_branch      = true,
                                   rng_seed          = 2)
    end

    sol_loose = solve_at(1e-10)
    sol_tight = solve_at(1e-12)

    # ---- Sample points -------------------------------------------------
    # FF3.1.3 sheet-0 sample: well inside the walked region, upper-half,
    # ~3 walker steps from the IC root.
    ζ_sample = 1.0 + 1.5im
    w_loose_s0, _ = eval_at_sheet(sol_loose, ζ_sample, [0]; extrapolate=true)
    w_tight_s0, _ = eval_at_sheet(sol_tight, ζ_sample, [0]; extrapolate=true)
    err_s0 = abs(w_loose_s0 - w_tight_s0)

    # FF3.1.6 cut-straddling samples: ζ just below Im ζ = π (sheet [0]
    # side) and ζ just above (sheet [1] side, reached by winding-ring
    # cut crossing).  Different sheets at the same ζ point should give
    # different `arg w`.
    r_cut = 2.0; ε_cut = 0.05
    ζ_below = complex(log(r_cut), π - ε_cut)
    w_below_s0, _ = eval_at_sheet(sol_tight, ζ_below, [0]; extrapolate=true)
    w_below_s1, _ = eval_at_sheet(sol_tight, ζ_below, [1]; extrapolate=true)

    # ---- FF3.1.2 sample points finite on both sheets --------------------
    @testset "FF3.1.2: sample points finite, sheets [0] and [1]" begin
        @test isfinite(real(w_loose_s0)) && isfinite(imag(w_loose_s0))
        @test isfinite(real(w_tight_s0)) && isfinite(imag(w_tight_s0))
        @test isfinite(real(w_below_s0)) && isfinite(imag(w_below_s0))
        @test isfinite(real(w_below_s1)) && isfinite(imag(w_below_s1))
    end

    # ---- FF3.1.3 sheet-0 absolute pin ----------------------------------
    # Captured 2026-05-16 from a tight-tol run with the figure-script
    # target set + seed=2.  Loose-vs-tight is 0.0 at this sample (the
    # walker converges to the same path; FF2.1.3 worklog 044 §"two-
    # pronged" pattern).  The absolute pin is the load-bearing
    # mutation-detector — flips under M1 (corrupt α).
    @testset "FF3.1.3: ζ-plane w_sheet0 at $ζ_sample" begin
        @info "FF3.1.3: |w_loose - w_tight| (sheet 0) = $err_s0"
        @test err_s0 ≤ 6e-6
        # Reference captured 2026-05-16 from solve_at(1e-12) with
        # rng_seed = 2 and the figure-script target set.
        w_ref = 0.40810323599649645 - 0.052436433262320184im
        @test abs(w_tight_s0 - w_ref) ≤ 1e-12
    end

    # ---- FF3.1.4 conjugate symmetry on sheet 0 --------------------------
    # Schwarz reflection: for real PVI parameters + real IC,
    # `w(ζ̄) = w(ζ)` complex-conjugated, so `arg w(ζ̄) = -arg w(ζ)`
    # mod 2π.  Sample a z-plane point on sheet 0 and its conjugate.
    @testset "FF3.1.4: conjugate symmetry (sheet 0)" begin
        z_test = 5.0 + 1.0im
        ζ_test = complex(log(abs(z_test)), angle(z_test))
        ζ_conj = complex(log(abs(z_test)), -angle(z_test))
        w_test, _ = eval_at_sheet(sol_tight, ζ_test, [0]; extrapolate=true)
        w_conj, _ = eval_at_sheet(sol_tight, ζ_conj, [0]; extrapolate=true)
        @info "FF3.1.4: arg w(ζ_test) + arg w(ζ_conj) = " *
              "$(angle(w_test) + angle(w_conj))"
        @test isfinite(real(w_test))
        @test isfinite(real(w_conj))
        # arg sum ≈ 0 to ≤ 1e-8 rad.  Schwarz reflection holds on the
        # well-walked principal sheet because the visited tree's lookup
        # is conjugate-symmetric here (cross-mode doesn't break this on
        # sheet [0] alone).
        @test abs(angle(w_test) + angle(w_conj)) ≤ 1e-8
    end

    # ---- FF3.1.5 visited-sheet population --------------------------------
    @testset "FF3.1.5: cross-mode walker populates non-[0] sheet" begin
        counts = Dict{Vector{Int},Int}()
        for s in sol_tight.visited_sheet
            counts[s] = get(counts, s, 0) + 1
        end
        @info "FF3.1.5: tight-tol visited_sheet distribution = $counts"
        nonzero_counts = [n for (s, n) in counts if s != [0]]
        @test !isempty(nonzero_counts)
        @test maximum(nonzero_counts) > 5
    end

    # ---- FF3.1.6 multi-sheet phase distinction --------------------------
    # The load-bearing Fig-3 claim: at the same ζ point, sheet [0] vs
    # sheet [1] queries return distinct values because the visited
    # trees on the two sheets have distinct nearest nodes (the sheet
    # bookkeeping is wired through to `eval_at_sheet`).  Under M4
    # (drop sheet bookkeeping in the walker's wedge loop) either:
    #   - sheet [1] would have NO visited nodes ⇒ NaN return ⇒ RED
    #     (we assert finiteness), OR
    #   - sheet [1] queries would resolve to the same nearest visited
    #     as sheet [0] ⇒ identical return ⇒ |Δarg| = 0 ⇒ RED.
    # Empirically the unmutated tight-tol run gives |Δarg| ≈ 0.018 rad
    # at our sample; we assert > 1e-3 rad (50× safety margin).
    @testset "FF3.1.6: sheet [0] vs sheet [1] arg distinct" begin
        @test isfinite(real(w_below_s0))
        @test isfinite(real(w_below_s1))
        Δarg = abs(angle(w_below_s0) - angle(w_below_s1))
        @info "FF3.1.6: |arg w(ζ_below,s0) - arg w(ζ_below,s1)| = $Δarg"
        @test Δarg > 1e-3
        # Pin the captured tight-tol value (Δarg ≈ 0.0176 at our seed +
        # target set).  Tolerance picked to allow walker-internal
        # numerical drift across Julia versions while still catching
        # bookkeeping regressions.
        Δarg_ref = 0.017594460150588437
        @test abs(Δarg - Δarg_ref) ≤ 1e-6
    end
end

# Mutation-proof procedure (verified 2026-05-16):
#
#   Mutation M1  --  corrupt α: 4.0 → 5.0 in the parameter assignment
#     at the top of the @testset.  Expected bite: FF3.1.3 absolute pin
#     RED (w_tight at sample becomes a different complex number when
#     α changes the underlying analytic solution).
#
#   Mutation M3  --  in `solve_at`, drop `cross_branch = true` kwarg
#     (refuse mode + winding ring).  Expected bite: 1 ERROR — the
#     walker corners against the cut trying to reach winding-ring
#     targets at Re < 0; `path_network_solve` throws
#     "all 5 wedge candidates failed".  All downstream tests error.
#
#   Mutation M4  --  in src/PathNetwork.jl wedge-loop sheet-update
#     block, always copy parent's sheet (ignore cross_branch flag).
#     Expected bite: FF3.1.5 RED (no non-[0] sheets visited) + FF3.1.6
#     RED (either NaN or Δarg = 0).
#
# Bite counts captured in `docs/worklog/046-ffw2017-fig-3.md`.
