# external/probes/stokes-strip-audition/probe.jl
#
# Phase-E1 audition (bead `padetaylor-0ln.37.18`, ADR-0025 v1-corner C4).
#
# THE CORNER.  The headline figure NaN-masks a ±3°-wide strip straddling
# each ±36° Stokes line (`surf_in_mask`: `|arg x| ∈ [33°,39°]`).  The
# sector ray fan is inset `SURF_SECTOR_MARGIN_DEG = 4°` from the Stokes
# line (the sector solver covers `|arg x| > 40°`); the wedge pole field
# reaches `|arg| ≈ 34°`.  The band `[34°,40°]` is covered by neither —
# rendered grey.
#
# THE AUDITION QUESTION.  The pole-free *sector* solution is analytic
# (harmonic) right up to the Stokes line itself — the 4° ray-fan inset
# and the 3° mask were a *conservative cushion*, not a hard numerical
# limit.  How close to the 36° Stokes line can the sector ray fan and
# the harmonic voters be pushed while staying honest (triple-method vote
# spread controlled)?  Can the 4° margin / 3° mask be reduced?
#
# METHOD.  The sector ray fan's angular extent is governed by
# `SURF_SECTOR_MARGIN_DEG`: the fan runs `arg x ∈ (36°+margin,
# 324°-margin)` and the two Laplace voters solve the conformal
# `w = log x` rectangle on that same angular window.  We rebuild the
# sector pipeline (`surf_ray_fan` + `surf_laplace_voters`) at a sweep of
# margins and, for each, measure the triple-method vote spread in
# annular angular bands marching toward the Stokes line — the same
# `surf_vote` median spread the figure's `spread` map records.  A margin
# is "honest" if its near-Stokes vote spread stays within the bulk
# triple-method envelope (the figure asserts bulk median `< 3e-3`, max
# `< 1e-2`; PI2S.2).
#
# We also probe the *genuine* near-Stokes accuracy with an independent
# oracle: a dedicated through-the-point BVP straight along the query ray
# (the same oracle idiom PI2S.11's C2/C3 audition used) — so we separate
# "the voters agree" (cheap, gameable) from "the voters are accurate".
#
# Run:  julia --project=figures external/probes/stokes-strip-audition/probe.jl

using Printf
using Statistics: median

const FIGDIR = joinpath(@__DIR__, "..", "..", "..", "figures")
include(joinpath(FIGDIR, "_kkg_pi2_surface_helpers.jl"))

using PadeTaylor.VectorBVP: vector_bvp_solve

# ----------------------------------------------------------------------
# Re-implement `surf_ray_fan` / `surf_laplace_voters` parametrised on the
# sector margin (the production helpers hard-read `SURF_SECTOR_MARGIN_DEG`).
# Everything else — the BVP recipe, the Laplace solvers, the C3
# reconstruction — is the production code, unchanged.
# ----------------------------------------------------------------------

"Ray fan at an explicit angular margin (degrees inside the Stokes line)."
function fan_at_margin(margin_deg::Real)
    f, Jf = surf_companion()
    φ_lo = deg2rad(SURF_WEDGE_HALF_DEG + margin_deg)
    φ_hi = deg2rad(360.0 - SURF_WEDGE_HALF_DEG - margin_deg)
    nφ   = max(2, round(Int, (φ_hi - φ_lo) / deg2rad(SURF_RAY_DPHI_DEG)) + 1)
    phis = collect(range(φ_lo, φ_hi; length = nφ))
    nr    = 40
    radii = collect(range(SURF_R_MIN, SURF_R_MAX; length = nr))
    sols  = Vector{VectorBVPSolution}(undef, nφ)
    U     = Matrix{ComplexF64}(undef, nr, nφ)
    for (j, φ) in enumerate(phis)
        sol = surf_ray_bvp(f, Jf, φ)
        sols[j] = sol
        for (i, r) in enumerate(radii)
            U[i, j] = sol(ComplexF64(r * cis(φ)))[1]
        end
    end
    return (phis = phis, radii = radii, sols = sols, U = U)
end

"The conformal rectangle for an explicit margin."
rect_at_margin(margin_deg) =
    (s_lo = log(SURF_R_MIN), s_hi = log(SURF_R_MAX),
     θ_lo = deg2rad(SURF_WEDGE_HALF_DEG + margin_deg),
     θ_hi = deg2rad(360.0 - SURF_WEDGE_HALF_DEG - margin_deg))

"Voters 2 & 3 on an explicit-margin rectangle (production solvers)."
function voters_at_margin(fan, rect)
    s_lo, s_hi = rect.s_lo, rect.s_hi
    θ_lo, θ_hi = rect.θ_lo, rect.θ_hi
    re_inner(θ) = real(surf_ray_arc_eval(fan, SURF_R_MIN, θ))
    re_outer(θ) = real(surf_ray_arc_eval(fan, SURF_R_MAX, θ))
    re_elo(s)   = real(surf_ray_eval(fan, ComplexF64(exp(s) * cis(θ_lo))))
    re_ehi(s)   = real(surf_ray_eval(fan, ComplexF64(exp(s) * cis(θ_hi))))
    im_inner(θ) = imag(surf_ray_arc_eval(fan, SURF_R_MIN, θ))
    im_outer(θ) = imag(surf_ray_arc_eval(fan, SURF_R_MAX, θ))
    im_elo(s)   = imag(surf_ray_eval(fan, ComplexF64(exp(s) * cis(θ_lo))))
    im_ehi(s)   = imag(surf_ray_eval(fan, ComplexF64(exp(s) * cis(θ_hi))))
    rem2 = laplace2d_solve(re_inner, re_outer, re_elo, re_ehi,
                           (s_lo, s_hi), (θ_lo, θ_hi);
                           Nx = SURF_LAP_NX, Ny = SURF_LAP_NY)
    imm2 = laplace2d_solve(im_inner, im_outer, im_elo, im_ehi,
                           (s_lo, s_hi), (θ_lo, θ_hi);
                           Nx = SURF_LAP_NX, Ny = SURF_LAP_NY)
    rem3 = laplace2d_solve_gridap(re_inner, re_outer, re_elo, re_ehi,
                                  (s_lo, s_hi), (θ_lo, θ_hi);
                                  n_x = SURF_FEM_NX, n_y = SURF_FEM_NY)
    imm3 = laplace2d_solve_gridap(im_inner, im_outer, im_elo, im_ehi,
                                  (s_lo, s_hi), (θ_lo, θ_hi);
                                  n_x = SURF_FEM_NX, n_y = SURF_FEM_NY)
    return (rem2 = rem2, imm2 = imm2, rem3 = rem3, imm3 = imm3, rect = rect)
end

"Independent oracle: a dedicated through-the-point BVP straight along φ."
function oracle_at(f, Jf, z)
    φ = angle(z)
    CT = ComplexF64
    z_a = SURF_R_MAX * cis(φ); z_b = 1.02 * cis(φ)
    Ba = zeros(CT, 4, 4); Bb = zeros(CT, 4, 4)
    Ba[1,1] = 1; Ba[2,2] = 1; Bb[3,1] = 1; Bb[4,2] = 1
    sa = pI2_tritronquee_ic(z_a; t = SURF_T, n_terms = 2)
    sb = pI2_tritronquee_ic(z_b; t = SURF_T, n_terms = 2)
    sol = vector_bvp_solve(f, z_a, z_b, Ba, Bb,
                           CT[sa[1], sa[2], sb[1], sb[2]];
                           N = 220, tol = 1e-9, maxiter = 40,
                           jacobian = Jf,
                           initial_guess = z -> pI2_tritronquee_ic(z;
                               t = SURF_T, n_terms = 2))
    return sol(ComplexF64(z))[1]
end

# ----------------------------------------------------------------------
# The sweep.
# ----------------------------------------------------------------------

# Angular bands (degrees off the +real axis) on the SECTOR side of the
# 36° Stokes line, marching DOWN toward 36°.  The sector solver lives at
# `arg x > 36°`; the masked strip's sector half is `[36°,39°]`.  The
# audition asks how far into `[36°,40°]` the sector pipeline can honestly
# reach.  Each band is a 1°-wide (or finer) annular wedge sampled on a
# polar grid; we report the per-band triple-method vote spread.
const BANDS = [(44.0, 45.0), (41.0, 42.0), (40.0, 41.0), (39.0, 40.0),
               (38.0, 39.0), (37.0, 38.0), (36.5, 37.0), (36.1, 36.5)]
# Radii to sample within each band (mid-sector — away from the inner-arc
# seed floor that PI2S.2 already characterises).
const PROBE_RADII = collect(5.0:1.5:18.0)

println("="^72)
println("Phase-E1 Stokes-strip audition — bead padetaylor-0ln.37.18")
println("="^72)
println()
println("Sweep: sector ray-fan margin (° inside the 36° Stokes line).")
println("Each row = a margin; columns = triple-method vote spread")
println("(median over a polar grid) in 1°-wide angular bands toward 36°.")
println()

f0, Jf0 = surf_companion()

# margin -> per-band (spread_med, spread_max, oracle_err_med, oracle_err_max)
results = Dict{Float64,Any}()

for margin in (4.0, 3.0, 2.0, 1.0, 0.5)
    fan  = fan_at_margin(margin)
    rect = rect_at_margin(margin)
    lap  = voters_at_margin(fan, rect)
    band_stats = NamedTuple[]
    for (b_lo, b_hi) in BANDS
        # The band is only covered when the fan reaches it: the fan's
        # angular edge is at 36°+margin.  A band entirely inside the
        # margin cushion has no voter datum.
        spreads = Float64[]; orerrs = Float64[]
        for degθ in range(b_lo + 0.1, b_hi - 0.1; length = 4)
            for sgn in (+1.0,)        # upper Stokes line; lower is mirror
                φ = deg2rad(sgn * degθ)
                for r in PROBE_RADII
                    z = ComplexF64(r * cis(φ))
                    # Voter 1 — direct ODE solve (C3 reconstruction).
                    m1 = surf_ray_eval(fan, z)
                    # Voters 2/3 — the two Laplace solves.
                    l2 = surf_laplace_eval_one(lap.rem2, lap.imm2, lap.rect, z)
                    l3 = surf_laplace_eval_one(lap.rem3, lap.imm3, lap.rect, z)
                    (isnan(real(m1)) || l2 === nothing || l3 === nothing) &&
                        continue
                    re1, im1 = real(m1), imag(m1)
                    _, sre = surf_vote(re1, l2.re, l3.re)
                    _, sim = surf_vote(im1, l2.im, l3.im)
                    push!(spreads, max(sre, sim))
                    # Oracle accuracy of the *voted* value.
                    vre, _ = surf_vote(re1, l2.re, l3.re)
                    vim, _ = surf_vote(im1, l2.im, l3.im)
                    oref = oracle_at(f0, Jf0, z)
                    push!(orerrs,
                          abs(ComplexF64(vre, vim) - oref))
                end
            end
        end
        push!(band_stats, (band = (b_lo, b_hi),
                           n = length(spreads),
                           spread_med = isempty(spreads) ? NaN :
                               median(spreads),
                           spread_max = isempty(spreads) ? NaN :
                               maximum(spreads),
                           orerr_med = isempty(orerrs) ? NaN :
                               median(orerrs),
                           orerr_max = isempty(orerrs) ? NaN :
                               maximum(orerrs)))
    end
    results[margin] = band_stats
    @printf("margin = %.1f°  (fan edge at arg x = %.1f°)\n",
            margin, 36.0 + margin)
    @printf("  %-12s %4s %12s %12s %12s %12s\n",
            "band(°)", "n", "spread_med", "spread_max",
            "orerr_med", "orerr_max")
    for bs in band_stats
        @printf("  [%4.1f,%4.1f]  %4d %12.3e %12.3e %12.3e %12.3e\n",
                bs.band[1], bs.band[2], bs.n,
                bs.spread_med, bs.spread_max,
                bs.orerr_med, bs.orerr_max)
    end
    println()
end

# ----------------------------------------------------------------------
# Does the solution genuinely degrade AT the Stokes line?  Walk a single
# ray's oracle BVP angle from 30° toward 36° and report the oracle error
# of the n_terms=2 seed-pinned solve — if the BVP itself blows up near
# the Stokes line, that is a real reason to mask; if it stays bounded,
# the mask is purely the fan-inset artefact.
# ----------------------------------------------------------------------
println("-"^72)
println("Near-Stokes BVP health — does the sector ODE solve degrade AT 36°?")
println("(dedicated through-the-point BVP, Newton residual + companion")
println(" consistency; a blow-up near 36° would be a real masking reason)")
println("-"^72)
for degθ in (30.0, 33.0, 35.0, 35.5, 35.9, 36.0)
    φ = deg2rad(degθ)
    CT = ComplexF64
    z_a = SURF_R_MAX * cis(φ); z_b = 1.02 * cis(φ)
    Ba = zeros(CT, 4, 4); Bb = zeros(CT, 4, 4)
    Ba[1,1] = 1; Ba[2,2] = 1; Bb[3,1] = 1; Bb[4,2] = 1
    sa = pI2_tritronquee_ic(z_a; t = SURF_T, n_terms = 2)
    sb = pI2_tritronquee_ic(z_b; t = SURF_T, n_terms = 2)
    local sol
    ok = true
    try
        sol = vector_bvp_solve(f0, z_a, z_b, Ba, Bb,
                               CT[sa[1], sa[2], sb[1], sb[2]];
                               N = 220, tol = 1e-9, maxiter = 40,
                               jacobian = Jf0,
                               initial_guess = z -> pI2_tritronquee_ic(z;
                                   t = SURF_T, n_terms = 2))
    catch err
        ok = false
        @printf("  arg x = %5.2f°  BVP FAILED: %s\n", degθ,
                sprint(showerror, err))
        continue
    end
    d8 = surf_vc8_companion_check(sol)
    @printf("  arg x = %5.2f°  Newton resid = %.2e   companion-consist = %.2e\n",
            degθ, d8.residual_inf, d8.consistency_max)
end

println()
println("="^72)
println("Audition complete.  See REPORT.md for the decision.")
println("="^72)
