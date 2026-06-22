# figures/_kkg_pi2_sector.jl
#
# Region 1 — the pole-free-sector triple-method pipeline of the KKG 2015
# P_I^(2) tritronquée headline figure.  This is the Rule-6 file-size
# split of the F3 compute kernel `_kkg_pi2_surface_helpers.jl`, sibling
# to `_kkg_pi2_vc45.jl` / `_kkg_pi2_vc8.jl` / `_kkg_pi2_vc10.jl`; the
# kernel `include`s it.  It holds the three sector voters and the
# per-grid-point majority vote (ADR-0024) — everything Region-1.
#
# ## The triple-method majority vote (ADR-0024 — unchanged here)
#
# On the pole-free sector `arg x ∈ (36°, 324°)` the tritronquée `V_0` is
# analytic, so `Re u` and `Im u` are harmonic.  ADR-0024 computes the
# surface THREE algorithmically-disjoint ways and majority-votes per
# grid point:
#
#   * **Voter 1 — the ray-fan BVP** (`surf_ray_fan` / `surf_ray_eval`):
#     a fan of radial `vector_bvp_solve` rays through the sector — a
#     direct ODE solve, genuinely independent of the two harmonic
#     solvers.
#   * **Voter 2 — the in-house 2D-Chebyshev Laplace solve**
#     (`laplace2d_solve`): the sector is conformally a rectangle under
#     `w = log x`; a tensor-product Chebyshev spectral Laplacian solves
#     the Dirichlet problem there.
#   * **Voter 3 — the Gridap FEM Laplace solve** (`laplace2d_solve_
#     gridap`): the same rectangle Dirichlet problem, an independent
#     finite-element discretisation.
#
# `surf_vote` reduces the three per component by the median — the
# 3-sample median is the majority vote (an outlier voter is discarded,
# the two that concur win) — and records the max pairwise spread as the
# agreement-map diagnostic.  **The vote and the two Laplace solvers are
# untouched by C2/C3** (ADR-0025 Amendment 10): C2 changes only the
# *source* of voters 2/3's arc Dirichlet data, C3 only voter 1's
# raster→Cartesian reconstruction.
#
# ## C2 — the inner arc is now a BVP-interior datum (retires v1 corner C1)
#
# v1 corner C1 (worklog 057, ADR-0025 Amendment 10): the inner-arc
# (`|x| = R_MIN ≈ 2`) Dirichlet datum of voters 2/3 was the KKG
# `n_terms = 2` asymptotic series — an *asymptotic* (large-`|x|`)
# expansion, only `O(10⁻³)`–`O(10⁻²)` accurate at `|x| ≈ 2`.  The first
# instinct — "source the arc data from the ray fan instead" — does *not*
# by itself fix the accuracy: the C2-exploration probe
# (`external/probes/c2-arc-data/probe.jl`) measured that the production
# ray BVP's 2+2 split boundary condition **pins** `(u, u')` at the inner
# endpoint `z_b = R_MIN·e^{iφ}` to the asymptotic seed *exactly*
# (`|prod(R_MIN) − seed(R_MIN)| = 0.0`).  Sourcing `fan.U[1,:]` from
# such a fan would just re-serve the same seed.
#
# The genuine fix the probe found and this file implements: solve each
# ray BVP over an **extended inward segment** `[R_MAX, R_INNER_BC]` with
# `R_INNER_BC = SURF_R_INNER_BC < R_MIN`, so the inner arc `|x| = R_MIN`
# is an *interior* collocation point.  The BVP value there is then the
# genuine ODE-constrained solution — the asymptotic-seed error is
# confined to the `[R_INNER_BC, R_MIN]` boundary layer below the arc and
# does not reach it.  The probe measured the mean inner-arc error drop
# from `~5·10⁻⁴` (seed) to `~1·10⁻⁴` (BVP-interior, `R_INNER_BC = 1.05`)
# — a genuine ~5× accuracy gain, not a spread-metric cosmetic.  Voters
# 2/3 then read **all four** `w`-rectangle edges from the fan: the side
# edges as before (the bracketing rays), and now the inner/outer arcs
# from the fan's innermost/outermost interior samples too (`surf_ray_
# arc_data`).  The outer arc was already `~2·10⁻⁸` accurate from the
# seed at `|x| = 20` (the probe confirmed extending outward moves it
# only `~2·10⁻⁸`), so the outer endpoint stays at `R_MAX`; uniform
# fan-sourcing of every edge just removes the asymptotic series from the
# voter-2/3 boundary path entirely (C1 fully retired).
#
# The wider `[R_MAX, R_INNER_BC]` segment has more curvature near its
# inner end; `SURF_BVP_N = 128` keeps it spectrally resolved (the
# C2-exploration probe measured companion-consistency `~1·10⁻¹¹` at
# `N = 128` for `[20, 1.05]`, vs `~1·10⁻⁷` — marginal — at `N = 96`).
#
# ## C3 — voter-1 reconstruction: exact-radius + cubic-angular (v1 corner C5)
#
# v1 corner C5: voter 1 sampled the ray fan onto a `40`-row polar raster
# and **bilinearly** interpolated that raster onto the Cartesian grid.
# Bilinear interpolation is not a harmonic reconstruction, and the C3
# audition (`external/probes/c2-arc-data/c3_audition.jl`) measured it is
# materially lossy at the production `6°` fan resolution: inter-ray
# bilinear error median `7.9·10⁻⁴` / max `1.3·10⁻³` — *larger* than the
# bulk triple-method spread (`~5.5·10⁻⁴`), so the bilinear voter was the
# least-accurate of the three and its interpolation error fed straight
# into the vote.
#
# The audition's verdict (ADR-0025 Amendment 10): a *harmonic* voter 1
# — a Laplace solve on the `w`-rectangle with the ray fan's own edge
# data — would be `≡ voter 2` (in-house 2D-Chebyshev `laplace2d_solve`)
# and collapse the triple-method independence ADR-0024 requires; it is
# **rejected**.  Plain bilinear is **also rejected** — it is genuinely
# too lossy.  The chosen reconstruction keeps voter 1 a *direct ODE
# solve* (independent of both Laplace voters) and only fixes the
# raster→Cartesian sampling:
#
#   * **exact in the radius.**  A converged `VectorBVPSolution` *is* a
#     barycentric spectral interpolant — exact (`~10⁻⁹`) at any radius.
#     `surf_ray_fan` stores the ray solutions; `surf_ray_eval` evaluates
#     them directly at the query radius.  The `40`-row raster and its
#     radial interpolation error (`~3.6·10⁻⁴` on-ray, audition) are
#     gone — the radial direction is now machine-exact.
#   * **cubic (Catmull-Rom) in the angle.**  For a query between rays
#     `j` and `j+1`, the four bracketing rays `j-1 … j+2` are each
#     evaluated exactly at the query radius and blended by a Catmull-Rom
#     cubic.  The audition measured inter-ray error drop from `1.3·10⁻³`
#     (bilinear) to `6.9·10⁻⁵` max — a ~19× cut, now well below the
#     bulk spread and below the FEM voter.
#
# Voter 1 stays the ray-fan BVP voter; voters 2 and 3 stay the spectral
# and FEM Laplace solvers — three genuinely disjoint computations
# (ADR-0024 preserved).
#
# References:
#   * `docs/adr/0025-headline-figure-re-resolution.md` — Amendment 10
#     (C2 + C3); the v1-Corner Retirement Ledger rows C1, C5.
#   * `docs/adr/0024-laplace-harmonic-extension.md` — the triple-method
#     majority-vote decision and the conformal `w = log x` rectangle.
#   * `external/probes/c2-arc-data/probe.jl` — the C2 exploration:
#     the pinned-endpoint finding + the extend-inward accuracy gain.
#   * `external/probes/c2-arc-data/c3_audition.jl` — the C3 audition:
#     bilinear-vs-exact-r-cubic-φ measured interpolation error.
#   * `src/VectorBVP.jl` — `vector_bvp_solve`; the barycentric-callable
#     `VectorBVPSolution` voter 1 now evaluates exactly.
#   * `src/Laplace2D.jl`, `figures/_kkg_pi2_gridap_helper.jl` — voters
#     2 and 3, untouched by C2/C3.

# ======================================================================
# Region 1, Voter 1 — the ray-fan BVP (C2-extended, C3-reconstructed)
# ======================================================================

"""
    surf_companion() -> (f, Jf)

The `P_I^(2)` companion RHS and analytic Jacobian at `t = SURF_T`.
"""
surf_companion() = (painleve_hierarchy(:I, 2; t = SURF_T),
                    painleve_hierarchy_jacobian(:I, 2; t = SURF_T))

"""
    surf_ray_bvp(f, Jf, φ; r_in = SURF_R_INNER_BC) -> VectorBVPSolution

Solve the `P_I^(2)` tritronquée BVP on the single radial ray
`x = ξ·e^{iφ}`, `ξ ∈ [r_in, SURF_R_MAX]`.  The 2+2 split boundary
condition pins `(u, u')` at both ray endpoints from the F1-generalised
`pI2_tritronquee_ic` (valid on any `V_0`-sheet ray); the interior guess
is the same seed sampled at every collocation node.

**C2 (ADR-0025 Amendment 10).**  `r_in` defaults to
`SURF_R_INNER_BC = 1.05`, *below* the inner arc `SURF_R_MIN = 2`.  The
inner endpoint's `(u, u')` are still pinned to the asymptotic seed —
but the seed's truncation error is now confined to the
`[SURF_R_INNER_BC, SURF_R_MIN]` boundary layer, and the arc `R_MIN`
itself is an *interior* collocation point whose value the ODE
determines (the C2-exploration probe: `~1·10⁻⁴` accurate vs the seed's
`~5·10⁻⁴` directly at the arc).  A 3-argument call keeps the legacy
`[R_MAX, R_MIN]` behaviour available only by passing `r_in = SURF_R_MIN`
explicitly; the figure and VC-8 use the default extended segment.
"""
function surf_ray_bvp(f, Jf, φ::Real; r_in::Real = SURF_R_INNER_BC)
    CT = ComplexF64
    z_a = SURF_R_MAX * cis(φ)        # outer end (smaller seed error)
    z_b = r_in * cis(φ)              # inner end — below the inner arc (C2)
    Ba = zeros(CT, SURF_D, SURF_D)
    Bb = zeros(CT, SURF_D, SURF_D)
    Ba[1, 1] = 1; Ba[2, 2] = 1       # rows 1-2: (u, u') at z_a
    Bb[3, 1] = 1; Bb[4, 2] = 1       # rows 3-4: (u, u') at z_b
    seed_a = pI2_tritronquee_ic(z_a; t = SURF_T, n_terms = 2)
    seed_b = pI2_tritronquee_ic(z_b; t = SURF_T, n_terms = 2)
    g = CT[seed_a[1], seed_a[2], seed_b[1], seed_b[2]]
    return vector_bvp_solve(f, z_a, z_b, Ba, Bb, g;
                            N        = SURF_BVP_N,
                            tol      = SURF_BVP_TOL,
                            max_iter  = SURF_BVP_MAXITER,
                            jacobian = Jf,
                            initial_guess =
                                z -> pI2_tritronquee_ic(z; t = SURF_T,
                                                        n_terms = 2))
end

"""
    surf_ray_fan() -> NamedTuple

Voter 1.  Solve the BVP on a fan of rays covering the pole-free sector
(angle step `SURF_RAY_DPHI_DEG`, a few degrees inside the Stokes lines).
Each ray's full `VectorBVPSolution` is retained — the C3 reconstruction
(`surf_ray_eval`) evaluates them directly, exact in the radius.

Returns `(phis, radii, sols, U)`:
  - `phis`  : the ray angles (radians, ascending);
  - `radii` : the `[SURF_R_MIN, SURF_R_MAX]` raster — now an *interior*
              radial window of every BVP (the segment runs to
              `SURF_R_INNER_BC < SURF_R_MIN`, C2);
  - `sols`  : `Vector{VectorBVPSolution}`, one per ray angle — the
              exact barycentric spectral solutions C3 evaluates and
              C2's voter-2/3 boundary callables read;
  - `U`     : `Matrix{ComplexF64}`, `U[i,j] = u` at `radii[i]·e^{im·
              phis[j]}` — the polar raster, kept for inspection only
              (the C3 reconstruction no longer interpolates it).
"""
function surf_ray_fan()
    f, Jf = surf_companion()
    # The fan runs across the sector arg x ∈ (36°+margin, 324°-margin).
    φ_lo = deg2rad(SURF_WEDGE_HALF_DEG + SURF_SECTOR_MARGIN_DEG)
    φ_hi = deg2rad(360.0 - SURF_WEDGE_HALF_DEG - SURF_SECTOR_MARGIN_DEG)
    nφ   = max(2, round(Int, (φ_hi - φ_lo) / deg2rad(SURF_RAY_DPHI_DEG)) + 1)
    phis = collect(range(φ_lo, φ_hi; length = nφ))
    # The radial raster spans the inner..outer arc; every point is now a
    # BVP-interior radius (the segment extends to SURF_R_INNER_BC).
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

"""
    _surf_catmull_rom(p0, p1, p2, p3, t) -> value

Uniform Catmull-Rom cubic at parameter `t ∈ [0,1]` between knots `p1`
and `p2`, with `p0`/`p3` the flanking knots.  Used by `surf_ray_eval`
for the C3 angular blend across four bracketing fan rays — a cubic
through exact ray values is a far better inter-ray reconstruction than
the retired bilinear scheme (C3 audition: ~19× lower error).
"""
function _surf_catmull_rom(p0, p1, p2, p3, t)
    return 0.5 * ((2 * p1) + (-p0 + p2) * t +
                  (2 * p0 - 5 * p1 + 4 * p2 - p3) * t^2 +
                  (-p0 + 3 * p1 - 3 * p2 + p3) * t^3)
end

"""
    surf_ray_eval(fan, z) -> ComplexF64

Voter 1's value at the complex point `z` — the C3 reconstruction
(ADR-0025 Amendment 10).  **Exact in the radius**: the bracketing rays'
stored `VectorBVPSolution`s are evaluated directly at `|z|` (the
barycentric spectral interpolant, `~10⁻⁹`).  **Cubic in the angle**: a
Catmull-Rom blend across the four rays bracketing `arg z` — in the two
*boundary* angular cells, where the four-ray stencil is degenerate
(no genuine flanking ray on the Stokes-line side, and the angular field
is steep / 6°-under-sampled there), the blend falls back to **linear**
(monotone, bounded — no extrapolation swing).  This retires the
bilinear polar-raster scheme (v1 corner C5) without making voter 1 a
Laplace solve (which would be `≡ voter 2`).

Returns `NaN + NaN·im` when `z` is outside the fan's `[r, φ]` coverage —
the caller treats that as "voter 1 has no datum here".

A query landing *exactly* on the `[r, φ]` window boundary — the
voter-2/3 boundary-data callables of `surf_laplace_voters` sample the
four rectangle edges — is in-window: the bounds test carries a small
relative tolerance and then clamps the query into the closed window, so
a rounding wobble at `|z| = R_MIN` / `R_MAX` or `arg z = φ_lo` / `φ_hi`
is never mistaken for an out-of-coverage point.
"""
function surf_ray_eval(fan, z)
    r = abs(z)
    φ = mod2pi(angle(z))             # arg z mapped into (0, 2π)
    phis, sols = fan.phis, fan.sols
    rtol = 1.0e-9 * SURF_R_MAX
    φtol = 1.0e-9 * (phis[end] - phis[1])
    (r < SURF_R_MIN - rtol || r > SURF_R_MAX + rtol) &&
        return ComplexF64(NaN, NaN)
    (φ < phis[1] - φtol    || φ > phis[end] + φtol) &&
        return ComplexF64(NaN, NaN)
    # In-window (within tolerance): clamp into the closed window so the
    # angular search and the radial BVP evaluation never index past it.
    r = clamp(r, SURF_R_MIN, SURF_R_MAX)
    φ = clamp(φ, phis[1], phis[end])
    nφ = length(phis)
    j  = clamp(searchsortedlast(phis, φ), 1, nφ - 1)
    tφ = (φ - phis[j]) / (phis[j + 1] - phis[j])
    # Each bracketing ray is evaluated EXACTLY at the query radius r —
    # the BVP barycentric callable, no radial interpolation error.
    ray_at(jr) = sols[jr](ComplexF64(r * cis(phis[jr])))[1]
    # The Catmull-Rom cubic needs four GENUINE knots (`j-1 … j+2`).  In
    # the two boundary angular cells the four-ray stencil is degenerate
    # — there is no genuine flanking ray on the Stokes-line side — and a
    # one-sided clamped-stencil "cubic" *extrapolates*: near the Stokes
    # line the angular profile is steep and 6°-under-sampled, and a
    # clamped cubic there was measured to swing the voter-1 datum ~8·10⁻³
    # off the two-bracketing-ray bracket the harmonic voters concur
    # within (the C3-audition test points were all mid-sector — they
    # never exercised the boundary cells).  So: full Catmull-Rom only
    # where the stencil is genuine; LINEAR (monotone, bounded — the old
    # behaviour) in the boundary cells.  The cubic's ~19× interior gain
    # is kept; the boundary cells incur no regression.
    if j ≥ 2 && j ≤ nφ - 2
        return _surf_catmull_rom(ray_at(j - 1), ray_at(j),
                                 ray_at(j + 1), ray_at(j + 2), tφ)
    else
        return (1 - tφ) * ray_at(j) + tφ * ray_at(j + 1)
    end
end

# ======================================================================
# Region 1, Voters 2 & 3 — the two Laplace voters on the w = log x
# conformal rectangle
# ======================================================================

"""
    surf_sector_rectangle() -> NamedTuple

The conformal-map rectangle for the harmonic solves.  Under `w = log x`
the pole-free sector becomes a rectangle in `(s, θ) = (log|x|, arg x)`:
`s ∈ [log r_min, log r_max]`, `θ ∈ [φ_lo, φ_hi]` (the same angular
window as the ray fan).  Returns `(s_lo, s_hi, θ_lo, θ_hi)`.
"""
function surf_sector_rectangle()
    φ_lo = deg2rad(SURF_WEDGE_HALF_DEG + SURF_SECTOR_MARGIN_DEG)
    φ_hi = deg2rad(360.0 - SURF_WEDGE_HALF_DEG - SURF_SECTOR_MARGIN_DEG)
    return (s_lo = log(SURF_R_MIN), s_hi = log(SURF_R_MAX),
            θ_lo = φ_lo, θ_hi = φ_hi)
end

"""
    surf_ray_arc_eval(fan, r, θ) -> ComplexF64

The ray fan's value at radius `r`, angle `θ` — the C3 exact-radius +
cubic-angle reconstruction, reused to read the fan's inner/outer arcs.
`r` is `SURF_R_MIN` (inner arc) or `SURF_R_MAX` (outer arc); both are
interior radii of every ray BVP (C2), so the value is the genuine ODE
solution, not the asymptotic seed.  Out-of-window `θ` clamps to the
fan's angular edge (the rectangle's θ-range equals the fan's by
construction, so a clamp only ever fires on a Chebyshev node landing
exactly on the boundary).
"""
function surf_ray_arc_eval(fan, r::Real, θ::Real)
    return surf_ray_eval(fan, ComplexF64(r * cis(clamp(θ, fan.phis[1],
                                                       fan.phis[end]))))
end

"""
    surf_laplace_voters(fan) -> NamedTuple

Voters 2 & 3.  Build the Dirichlet boundary data on the `w = log x`
rectangle and run both Laplace solvers — the in-house 2D-Chebyshev
`laplace2d_solve` (voter 2) and the Gridap FEM `laplace2d_solve_gridap`
(voter 3) — once for `Re u` and once for `Im u`.

**C2 (ADR-0025 Amendment 10) — all four edges are BVP data.**  Every
edge of the `w`-rectangle now comes from the ray fan:
  - the two **side edges** `θ = θ_lo, θ_hi` are the fan's first/last
    rays, read exactly via `surf_ray_eval` as a function of `s`;
  - the **inner / outer arcs** `s = log r_min, log r_max` are the fan's
    innermost/outermost radial samples — *not* the KKG asymptotic
    series (v1 corner C1 retired).  Because the ray BVPs extend inward
    past `R_MIN` (`surf_ray_bvp`'s `r_in = SURF_R_INNER_BC`), the inner
    arc is an interior BVP datum (`~1·10⁻⁴` accurate, the C2 probe), not
    the `O(10⁻³·10⁻²)` asymptotic seed.

The two Laplace solvers and the median vote are otherwise unchanged
(ADR-0024 preserved).  Returns `(rem2, imm2, rem3, imm3, rect)` — the
four `Laplace2DSolution`s (Re/Im × voter-2/3) and the rectangle.
"""
function surf_laplace_voters(fan)
    rect = surf_sector_rectangle()
    s_lo, s_hi = rect.s_lo, rect.s_hi
    θ_lo, θ_hi = rect.θ_lo, rect.θ_hi

    # All four w-rectangle edges read the ray fan exactly (C3
    # reconstruction); the inner/outer arcs are interior BVP samples (C2).
    #   bc_x_lo : s = s_lo edge (inner arc),  function of θ
    #   bc_x_hi : s = s_hi edge (outer arc),  function of θ
    #   bc_y_lo : θ = θ_lo edge (side ray),   function of s
    #   bc_y_hi : θ = θ_hi edge (side ray),   function of s
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

"""
    surf_laplace_eval_one(re_sol, im_sol, rect, z) -> Union{NamedTuple, Nothing}

Evaluate a Laplace voter pair (voter 2 *or* voter 3) at the complex
point `z`: map `z` through `w = log z` into the rectangle and read the
`Re u` / `Im u` solutions.  Returns `(re, im)` or `nothing` when `z`
falls outside the rectangle (`|z| ∉ [r_min, r_max]` or `arg z` outside
the angular window).
"""
function surf_laplace_eval_one(re_sol, im_sol, rect, z)
    s = log(abs(z))
    θ = mod2pi(angle(z))
    (s < rect.s_lo || s > rect.s_hi) && return nothing
    (θ < rect.θ_lo || θ > rect.θ_hi) && return nothing
    return (re = re_sol(s, θ), im = im_sol(s, θ))
end

# ======================================================================
# Region 1 — the per-grid-point majority vote (ADR-0024; unchanged)
# ======================================================================

"""
    surf_vote(a, b, c) -> (voted, spread)

The 3-sample majority vote, per component.  `a, b, c` are the three
voter estimates of one scalar (`Re u` or `Im u`).  The **median** of
the three is the vote — it discards a single outlier and keeps the two
that concur.  `spread` is the maximum pairwise disagreement, the
agreement-map diagnostic.

If a voter is `NaN` (a method has no datum at this point) it is dropped:
with two finite voters the vote is their mean and `spread` their
difference; with one, that one; with none, `NaN`.
"""
function surf_vote(a, b, c)
    vals = filter(isfinite, (a, b, c))
    if length(vals) == 3
        x, y, z = vals
        voted  = x + y + z - min(x, y, z) - max(x, y, z)   # median
        spread = max(x, y, z) - min(x, y, z)
        return voted, spread
    elseif length(vals) == 2
        return (vals[1] + vals[2]) / 2, abs(vals[1] - vals[2])
    elseif length(vals) == 1
        return vals[1], 0.0
    else
        return NaN, NaN
    end
end
