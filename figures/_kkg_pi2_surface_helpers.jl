# figures/_kkg_pi2_surface_helpers.jl
#
# F3 — the whole-plane compute kernel for the KKG 2015 P_I^(2)
# tritronquée surface figure (bead `padetaylor-0ln.33`, v0.2 plan).
# Makie-free: this file assembles `Re u` and `Im u` of the tritronquée
# `V_0(x, 0)` over the complex `x`-plane square `[-20,20]²` and returns
# the matrices ready for rendering; the figure script and the
# acceptance test `figures/test_kkg_pi2_surface.jl` both `include` it so
# they exercise the identical computation.
#
# This is the integration phase: F1 (the generalised
# `pI2_tritronquee_ic`), F2 (the Stage-2 fine-grid fill), the in-house
# 2D-Chebyshev Laplace solver (`src/Laplace2D.jl`) and the Gridap FEM
# helper (`figures/_kkg_pi2_gridap_helper.jl`) are all *already built*;
# F3 composes them. It runs in the `figures/` Julia environment because
# voter (3) uses Gridap, which lives only there (ADR-0024 *Correction*).
#
# ## The construction — two regions, one Cartesian grid
#
# The Type-I tritronquée `V_0` is a separatrix solution: pole-free on a
# wide angular sector and pole-rich elsewhere.  From KKG 2015 eq. (1.4)
# (pillar C §4, `docs/v0p2_pillarC_painleve_hierarchy_findings.md:481`)
# the pole-free sector, reduced mod 2π, is centred on the **negative**
# real axis with angular half-width `6π/7 − (3/7)·arctan(1/√5) ≈ 144°`.
# So:
#
#   * **Region 1 — the pole-free sector**, `arg x ∈ (36°, 324°)`
#     (equivalently `|arg x − 180°| < 144°`): `V_0` is analytic, hence
#     `Re u` and `Im u` are **harmonic**, and the surface is computed
#     THREE algorithmically-disjoint ways and **majority-voted per grid
#     point** (ADR-0024; the user requirement).
#
#   * **Region 2 — the pole-rich wedge**, `arg x ∈ (−36°, 36°)`: the
#     vector path-network walk + the F2 Stage-2 fine-grid fill, exactly
#     the V8b approach (`figures/_kkg_pi2_helpers.jl`).
#
# ### Region 1 — the triple-method majority vote
#
# **Method 1 — the ray-fan BVP (the direct voter).**  A fan of radial
# rays `x = ξ·e^{iφ}`, `ξ ∈ [r_min, r_max]`, with `φ` stepping across
# the sector.  Each ray is a `vector_bvp_solve` of the `d = 4` companion
# `P_I^(2)` system (probe-verified recipe, `_kkg_pi2_helpers.jl` Stage
# A), boundary-pinned at both ends by the F1-generalised
# `pI2_tritronquee_ic` — now valid on any `V_0`-sheet ray, not just the
# negative real axis.  The resulting polar grid is bilinearly
# interpolated onto the Cartesian grid.  (Bilinear interpolation across
# a fan of *exact* BVP rays is a sampling step, not a harmonic
# reconstruction — Methods 2/3 are the principled harmonic fills; the
# vote reconciles all three.)
#
# **Method 2 — the in-house 2D-Chebyshev Laplace solve.**  A pole-free
# *angular sector* `{r₀ ≤ |x| ≤ r₁, θ₀ ≤ arg x ≤ θ₁}` is **conformally
# a rectangle** under `w = log x`: `Re w = log|x| ∈ [log r₀, log r₁]`,
# `Im w = arg x ∈ [θ₀, θ₁]`.  The Laplacian is conformally invariant,
# so the sector Laplace problem IS a rectangle Laplace problem in `w`
# (ADR-0024 Context).  `PadeTaylor.laplace2d_solve` solves it on the
# rectangle with Dirichlet data: the two edge rays `θ = θ_lo, θ_hi` are
# Method-1's edge BVP rays (reused — no extra solve); the inner and
# outer arcs `s = log r_min, log r_max` are the `pI2_tritronquee_ic`
# asymptotic series.  Run once for `Re u`, once for `Im u`, then mapped
# back to the `x`-plane.
#
# **Method 3 — the Gridap FEM Laplace solve.**  `laplace2d_solve_gridap`
# (`figures/_kkg_pi2_gridap_helper.jl`), the same rectangle + the same
# Dirichlet data, via an algorithmically-disjoint finite-element
# discretisation — a FEM bug and a spectral bug cannot coincide.
#
# **The vote.**  At every Cartesian grid point inside the sector the
# three estimates are reduced per component by the **median** (the
# 3-sample median = the majority vote: an outlier voter is discarded,
# the two that concur win).  An **agreement / spread map** records the
# maximum pairwise disagreement per point — a figure-quality diagnostic
# of where the three methods concur (and, honestly, where they do not).
#
# ### Region 2 — the pole-rich wedge (ADR-0025 Amendment 2)
#
# The wedge is **not** a filled surface.  ADR-0025 Amendment 2 rescopes
# it after the A2 tractability probe (`external/probes/wedge-
# tractability/REPORT.md`) measured the *honest* coverage achievable by
# a path-network walk + per-node A1-gated disc and found it saturates
# at ~8.6 % (solver tol) — ~18 % even at display tol with order-48 jets.
# A filled honest wedge *surface* is numerically unreachable with the
# local-Padé-tiling architecture.  The wedge panel is therefore the
# **validated pole FIELD** — the extracted pole *locations* in the FW
# 2011 Fig 4.7 idiom — plus an **honest partial `|u|` surface underlay**
# only where the B1-gated node discs genuinely cover.
#
# Concretely Region 2 is:
#
#   1. **Pole field — the primary deliverable.**  An extended threading
#      fan (`surf_wedge_targets`) drives the B2 `:max_q_root` adaptive
#      walk through the *whole* wedge out to `|x| = 20`, so the walk's
#      nodes pass near every pole; `extract_poles_shared_q` with
#      cross-node `min_support ≥ 2` (VC-6) roots the per-node shared `Q`
#      and clusters the roots into the validated pole field.  The walk
#      is the B2 default (the fixed-`h = 0.1` V8b walk *blocks* past
#      `|x| ≈ 8` — A2 §3.1); the adaptive walk threads to `|x| ≈ 18-20`.
#
#   2. **Honest partial surface underlay.**  The F2 Stage-2 `fine_grid`
#      fill runs with `extrapolate = false`, so the B1 true-radius gate
#      (ADR-0025 Amendment 1, gate v-b) applies: a grid point outside
#      the nearest node's verified disc of validity gets `NaN` — no
#      Padé is ever evaluated outside its disc.  The kernel returns a
#      `covered` mask flagging exactly the (~13-18 %) cells the gated
#      discs honestly carry; everywhere else in the wedge is `NaN` and
#      renders neutral grey.  The partial surface is honest; the gaps
#      are explicit (Rule 1 fail-loud — an honest gap, never an
#      extrapolated lie).
#
# Per-pole validation **VC-4** (dominant-balance `A ∈ {-1,-3}`) and
# **VC-5** (conjugate-symmetry pairing) are wired in: the candidate
# pole field of `extract_poles_shared_q` (~380 poles, VC-6 cross-node
# filter only) is passed through `vc4_validate` — each candidate fitted
# to a local Laurent model on a 32-point ring, kept iff the leading
# coefficient is in the `{-1,-3}` dominant-balance family (VC-4a) and
# the residue vanishes (VC-4b), spurious candidates (Froissart doublets
# / out-of-family) pruned — so `kkg_pi2_surface().poles` is the
# VC-4-validated field.  `vc5_pair` then matches the survivors into
# conjugate pairs and reports the FW-style pairing-residual accuracy
# estimate.  The VC-4/VC-5 machinery lives in the sibling helper
# `figures/_kkg_pi2_vc45.jl` (CLAUDE.md Rule 6 file-size split).
# VC-7 (loop closure) remains a separate Phase-D bead.
#
# ### The stitch
#
# Region 1 (voted) + Region 2 are written onto one `[-20,20]²`
# Cartesian grid for `Re u` and `Im u`.  The thin strips straddling the
# sector/wedge boundary (`||arg x| − 36°| < STITCH_MASK_DEG`) are masked
# with `NaN`: the asymptotic seed degrades near the Stokes lines and
# neither region's solver is trustworthy there (Rule 9 v1 corner).
#
# ## v1 corners (CLAUDE.md Rule 9 — each with its forcing condition)
#
#   1. **Inner-arc asymptotic seed.**  Method 2/3's inner arc
#      (`|x| = r_min ≈ 2`) is the KKG `n_terms = 2` asymptotic series,
#      which is an *asymptotic* (large-`|x|`) expansion — at `|x| ≈ 2`
#      its truncation error is `O(10⁻²)`, far worse than at the outer
#      arc.  The harmonic solve damps this inward error (a Dirichlet
#      datum's error is not amplified in the interior), and the ray-fan
#      voter (Method 1, a genuine BVP) is *exact* on the inner arc, so
#      the vote's median is anchored by an accurate voter there.  *Forcing
#      condition for v2:* if a figure pin inside `|x| ≲ 3` must hold to
#      better than `1e-2`, replace the inner-arc datum with a small-`|x|`
#      Taylor/BVP patch.  (This is the sector's corner — Phase C.)
#
#   2. **RETIRED — wedge Stage-2 extrapolation.**  The shipped figure
#      called the Stage-2 fill with `extrapolate = true`, evaluating the
#      shared-`Q` Padé approximants far outside their verified disc (the
#      A4 baseline probe measured midpoints 50× outside the disc; ~70 %
#      of loop closures catastrophic).  ADR-0025 Amendment 1's B1
#      true-radius gate retires it: `surf_wedge_fill` now runs the fill
#      with `extrapolate = false`, so no Padé is evaluated outside its
#      B1-gated disc and the honest gaps render as `NaN`/grey.
#
#   3. **RETIRED — hand-tuned wedge step `h = 0.1`.**  ADR-0025
#      Amendment 4's B2 adaptive `:max_q_root` walk (the package
#      default) retires the hand-tune: the step direction dodges the
#      nearest shared-`Q` pole and `h` adapts to the local pole density.
#      The fixed-`h = 0.1` V8b walk *blocks* the extended fan past
#      `|x| ≈ 8` (A2 §3.1); the B2 walk threads it to `|x| ≈ 18-20`.
#
#   4. **Boundary-strip masking.**  `STITCH_MASK_DEG` of arc on each
#      side of the `±36°` Stokes lines is `NaN`-masked: the asymptotic
#      seed and both region solvers degrade at the sector boundary.
#      *Forcing condition for v2:* a uniform connection formula across
#      the Stokes line would let the two regions be stitched without a
#      gap (Phase E).
#
#   The wedge panel itself is no longer a v1 corner: ADR-0025
#   Amendment 2 rescopes it to the achievable senior-grade deliverable
#   — the validated pole field + honest partial underlay — and the A2
#   probe records the filled-surface ambition as a deferred bead with a
#   2D-fill-architecture forcing condition.
#
# References:
#   * `docs/adr/0025-headline-figure-re-resolution.md` — the
#     re-resolution ADR; **Amendment 1** the B1 true-radius gate,
#     **Amendment 2** the wedge-panel rescope (pole field + honest
#     partial underlay), **Amendment 4** the B2 adaptive `:max_q_root`
#     walk.  This kernel's Region 2 is the B3 deliverable of that ADR.
#   * `docs/adr/0024-laplace-harmonic-extension.md` — the triple-method
#     majority-vote decision; the conformal `w = log x` rectangle map;
#     the Gridap-as-figure-helper correction.
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015: eq. (1.1) the ODE, eq. (1.4) the `V_0` pole-free
#     sector, §7 Figs 7.4/7.5 the Re/Im surface this kernel assembles.
#   * `external/probes/wedge-tractability/REPORT.md` — A2: the wedge
#     honest-coverage saturation that forced the Amendment-2 rescope.
#   * `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — §4 the sector
#     geometry and the (sign-corrected) tritronquée IC.
#   * `figures/_kkg_pi2_helpers.jl` — V8b: the BVP recipe reused here
#     for Region 2's negative-axis anchor.
#   * `figures/_kkg_pi2_gridap_helper.jl` — voter (3).
#   * `src/Laplace2D.jl`, `src/VectorBVP.jl`, `src/VectorPathNetwork.jl`,
#     `src/VectorPathNetworkStage2.jl`, `src/VectorWedgeStep.jl`,
#     `src/PainleveHierarchy.jl`.

using PadeTaylor
using PadeTaylor.PainleveHierarchy: painleve_hierarchy,
                                    painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP:        vector_bvp_solve, VectorBVPSolution
using PadeTaylor.VectorProblems:   VectorPadeTaylorProblem
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution
using PadeTaylor.VectorPoleField:  extract_poles_shared_q

# `_kkg_pi2_gridap_helper.jl` brings in `laplace2d_solve_gridap` and
# (transitively) Gridap.  Include it once; guard against a double
# include when both helpers are pulled into one session.
if !isdefined(@__MODULE__, :laplace2d_solve_gridap)
    include(joinpath(@__DIR__, "_kkg_pi2_gridap_helper.jl"))
end

# VC-4 (dominant-balance per-pole certificate) and VC-5 (conjugate-
# symmetry pairing) — the Phase-D per-pole validation that certifies the
# genuine wedge poles and prunes the spurious ones (ADR-0025 Amendment 1
# §A3; beads `padetaylor-0ln.37.12` / `0ln.37.13`).
if !isdefined(@__MODULE__, :vc4_validate)
    include(joinpath(@__DIR__, "_kkg_pi2_vc45.jl"))
end

# ======================================================================
# Figure-defining constants — fixed here, the kernel is a pure function
# of them.
# ======================================================================

# The `P_I^(2)` parameter.  The KKG tritronquée surface figure is t = 0.
const SURF_T = 0.0

# Companion-system order: d = 2m = 4 for P_I^(2) (m = 2).
const SURF_D = 4

# The Cartesian render grid: [-20,20]² (KKG Figs 7.4/7.5).
const SURF_XY_LIM = 20.0
const SURF_GRID_N = 121          # 121×121 — odd so a node sits on 0.

# --- Region 1: the pole-free sector --------------------------------------
# Pole-free sector half-width 6π/7 − (3/7)arctan(1/√5) ≈ 144°, centred on
# the negative real axis (arg x = 180°).  So the sector is
# arg x ∈ (36°, 324°); the wedge is |arg x| < 36°.
const SURF_WEDGE_HALF_DEG = 36.0
# The ray fan stays a few degrees inside the Stokes lines — the seed
# degrades at the boundary (v1 corner 4).
const SURF_SECTOR_MARGIN_DEG = 4.0

# Ray-fan radii.  r_min ≈ 2 (the asymptotic seed is meaningless below
# |x| ≈ 1; 2 keeps a margin), r_max = 20 (the grid corner radius is
# 20√2 ≈ 28, but the sector pins are set at 20 and the corners outside
# the disc |x| ≤ 20 are NaN-masked — see SURF_R_MAX use).
const SURF_R_MIN = 2.0
const SURF_R_MAX = 20.0

# Angular step of the ray fan (≈ 6°).  The sector spans ≈ 280° usable.
const SURF_RAY_DPHI_DEG = 6.0

# BVP discretisation per ray.  N+1 Chebyshev nodes; N = 96 is converged
# (V8b probe: identical to ~1e-7 for N ≥ 40).
const SURF_BVP_N       = 96
const SURF_BVP_TOL     = 1.0e-9
const SURF_BVP_MAXITER = 40

# The two Laplace voters' rectangle resolution (in conformal w = log x
# coordinates).  Spectral (Method 2): a handful of nodes suffice.  FEM
# (Method 3): O(h²), needs a finer mesh.
const SURF_LAP_NX = 40           # radial (s = log|x|) Chebyshev nodes
const SURF_LAP_NY = 56           # angular (θ) Chebyshev nodes
const SURF_FEM_NX = 64           # FEM mesh subdivisions, radial
const SURF_FEM_NY = 88           # FEM mesh subdivisions, angular

# --- Region 2: the pole-rich wedge (ADR-0025 Amendment 2) ----------------
# The wedge deliverable is the validated pole FIELD + an honest partial
# `|u|` underlay (NOT a filled surface — A2 proved that unreachable).
#
# Seed: a BVP-anchored on-tritronquée state at `z = -3` (V8b Stage A).
const SURF_Z_SEED  = -3.0 + 0.0im

# Taylor order of the path-network jets.  Order 24 is the A2-verified
# value: A2 §5b.5 showed an order-48 *walk* degenerates the shared-`Q`
# Padé linear solve outright (all-below-τ singular spectrum), so order
# 24 is the senior-grade ceiling for the shared-`Q` route.
const SURF_PN_ORDER = 24

# The B2 adaptive-walk step *ceiling* `h_max`.  With `adaptive = true`
# (the package default — ADR-0025 Amendment 4) `h` is capped per step
# at `POLE_SAFETY·h·min|t*|` near a pole and grows geometrically toward
# this ceiling where the field is sparse.  `0.1` is the V8b ceiling;
# the adaptive controller shrinks it automatically in the dense field,
# so it is no longer hand-tuned — it is the *upper bound*, not the step.
const SURF_PN_H = 0.1

# The B2 walk truncation tolerance.  Passed to the B1 true-radius gate
# (`s(1e-8) = 0.36`, the production default — ADR-0025 Amendment 1) and
# to the adaptive step controller.  `1e-8` is the solver tolerance.
const SURF_PN_TOL = 1.0e-8

# --- The extended threading fan (the B3 production fan) ------------------
# ADR-0025 Amendment 2 rescopes B3 from "tile the wedge area" to "thread
# the walk through the *whole* wedge to extract the pole field densely
# and accurately to `|x| = 20`".  The fan is designed so the B2
# `:max_q_root` adaptive walk threads a node filament past every pole.
#
# Design (A2 §5.2 recommends a 9-angle `dr = 1` fan as the densest that
# threads; B2's adaptive walk threads it to `|x| ≈ 18-20` — Amendment 4):
#
#   * **Angles** — 9 rays fanned across `±0.5 rad ≈ ±28.6°`.  The wedge
#     half-width is `36°`; `±28.6°` keeps a `~7°` margin inside the
#     Stokes line so the fan stays well within the pole region and the
#     walk does not stray into the masked boundary strip.
#   * **Radii** — `dr = 1` spacing from `|x| = 2` out to `|x| = 20`,
#     i.e. 19 radial shells.  `9 × 19 = 171` targets — the A2 §5.2
#     "densest-completing" fan.  Denser fans (13-angle / `dr = 0.5`)
#     BLOCK (A2 §3.3): the extra inter-target bridging walks each get a
#     fresh chance to step onto a pole.
#
# Targets are ordered radius-major so the walk completes inner shells
# before reaching out — the threading order the B2 walk extends from.
const SURF_TARGET_RADII  = Tuple(2.0:1.0:20.0)             # 19 shells
const SURF_TARGET_ANGLES = Tuple(range(-0.5, 0.5; length = 9))  # rad

# Pole-extraction filter (`extract_poles_shared_q`).  `min_support = 2`
# is VC-6: a cluster is a physical pole only when ≥ 2 *distinct* nodes
# independently root it — keeps Froissart doublets (spurious
# single-node roots) out of the field.
const SURF_RADIUS_T     = 5.0
const SURF_CLUSTER_ATOL = 0.2
const SURF_MIN_SUPPORT  = 2

# --- The stitch ----------------------------------------------------------
# Half-width (degrees) of the NaN-masked strip on each Stokes line.
const SURF_STITCH_MASK_DEG = 3.0

# ======================================================================
# Geometry helpers
# ======================================================================

"""
    surf_in_sector(z) -> Bool

`true` when `z` lies in the `V_0` pole-free sector — `arg z` within
`SURF_WEDGE_HALF_DEG` of *neither* Stokes line, i.e. the angular
distance of `arg z` from the positive real axis exceeds the wedge
half-width.  `z = 0` is treated as in the wedge (its `arg` is ill-defined
and the origin sits inside the pole region for `t = 0`).
"""
function surf_in_sector(z)
    iszero(z) && return false
    return abs(rad2deg(angle(z))) > SURF_WEDGE_HALF_DEG
end

"""
    surf_in_mask(z) -> Bool

`true` when `z` is inside the NaN-masked boundary strip — within
`SURF_STITCH_MASK_DEG` of either `±SURF_WEDGE_HALF_DEG` Stokes line.
"""
function surf_in_mask(z)
    iszero(z) && return false
    d = abs(abs(rad2deg(angle(z))) - SURF_WEDGE_HALF_DEG)
    return d < SURF_STITCH_MASK_DEG
end

# ======================================================================
# Region 1, Method 1 — the ray-fan BVP
# ======================================================================

"""
    surf_companion() -> (f, Jf)

The `P_I^(2)` companion RHS and analytic Jacobian at `t = SURF_T`.
"""
surf_companion() = (painleve_hierarchy(:I, 2; t = SURF_T),
                    painleve_hierarchy_jacobian(:I, 2; t = SURF_T))

"""
    surf_ray_bvp(f, Jf, φ) -> VectorBVPSolution

Solve the `P_I^(2)` tritronquée BVP on the single radial ray
`x = ξ·e^{iφ}`, `ξ ∈ [SURF_R_MIN, SURF_R_MAX]`.  The 2+2 split boundary
condition pins `(u, u')` at both ray endpoints from the F1-generalised
`pI2_tritronquee_ic` (valid on any `V_0`-sheet ray); the interior guess
is the same seed sampled at every collocation node.
"""
function surf_ray_bvp(f, Jf, φ::Real)
    CT = ComplexF64
    z_a = SURF_R_MAX * cis(φ)        # outer end (smaller seed error)
    z_b = SURF_R_MIN * cis(φ)        # inner end
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
                            maxiter  = SURF_BVP_MAXITER,
                            jacobian = Jf,
                            initial_guess =
                                z -> pI2_tritronquee_ic(z; t = SURF_T,
                                                        n_terms = 2))
end

"""
    surf_ray_fan() -> NamedTuple

Method 1.  Solve the BVP on a fan of rays covering the pole-free sector
(angle step `SURF_RAY_DPHI_DEG`, a few degrees inside the Stokes lines).
Returns `(phis, radii, U)`:
  - `phis`  : the ray angles (radians, ascending);
  - `radii` : the common radial sample raster `[SURF_R_MIN, SURF_R_MAX]`;
  - `U`     : `Matrix{ComplexF64}` of `u = V_0(x)`, `U[i,j] = u` at
              `radii[i]·e^{im·phis[j]}` — the polar grid Methods 2/3 are
              voted against and the Cartesian fill bilinearly samples.
"""
function surf_ray_fan()
    f, Jf = surf_companion()
    # The fan runs across the sector arg x ∈ (36°+margin, 324°-margin).
    φ_lo = deg2rad(SURF_WEDGE_HALF_DEG + SURF_SECTOR_MARGIN_DEG)
    φ_hi = deg2rad(360.0 - SURF_WEDGE_HALF_DEG - SURF_SECTOR_MARGIN_DEG)
    nφ   = max(2, round(Int, (φ_hi - φ_lo) / deg2rad(SURF_RAY_DPHI_DEG)) + 1)
    phis = collect(range(φ_lo, φ_hi; length = nφ))
    # Radial raster: a modest number of points — each BVP solution is a
    # smooth barycentric interpolant, so the raster is just a sampling.
    nr    = 40
    radii = collect(range(SURF_R_MIN, SURF_R_MAX; length = nr))
    U = Matrix{ComplexF64}(undef, nr, nφ)
    for (j, φ) in enumerate(phis)
        sol = surf_ray_bvp(f, Jf, φ)
        for (i, r) in enumerate(radii)
            U[i, j] = sol(ComplexF64(r * cis(φ)))[1]
        end
    end
    return (phis = phis, radii = radii, U = U)
end

"""
    surf_ray_eval(fan, z) -> ComplexF64

Bilinearly interpolate the ray-fan polar grid `fan.U` at the
complex point `z` (in polar coordinates `(|z|, arg z)`).  Returns
`NaN + NaN·im` when `z` is outside the fan's `[r,φ]` coverage — the
caller treats that as "Method 1 has no datum here".
"""
function surf_ray_eval(fan, z)
    r = abs(z)
    # Map arg z into the fan's angular window [φ_lo, φ_hi] ⊂ (0, 2π).
    φ = mod2pi(angle(z))
    phis, radii, U = fan.phis, fan.radii, fan.U
    (r < radii[1] || r > radii[end]) && return ComplexF64(NaN, NaN)
    (φ < phis[1]  || φ > phis[end])  && return ComplexF64(NaN, NaN)
    i = searchsortedlast(radii, r); i = clamp(i, 1, length(radii) - 1)
    j = searchsortedlast(phis, φ);  j = clamp(j, 1, length(phis) - 1)
    tr = (r - radii[i]) / (radii[i + 1] - radii[i])
    tφ = (φ - phis[j])  / (phis[j + 1] - phis[j])
    u00 = U[i, j];     u10 = U[i + 1, j]
    u01 = U[i, j + 1]; u11 = U[i + 1, j + 1]
    return (1 - tr) * (1 - tφ) * u00 + tr * (1 - tφ) * u10 +
           (1 - tr) * tφ * u01      + tr * tφ * u11
end

# ======================================================================
# Region 1, Methods 2 & 3 — the two Laplace voters on the w = log x
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
    surf_laplace_voters(fan) -> NamedTuple

Methods 2 & 3.  Build the Dirichlet boundary data on the `w = log x`
rectangle and run both Laplace solvers — the in-house 2D-Chebyshev
`laplace2d_solve` (Method 2) and the Gridap FEM `laplace2d_solve_gridap`
(Method 3) — once for `Re u` and once for `Im u`.

Boundary data:
  - the two **edge rays** `θ = θ_lo, θ_hi` reuse `fan`'s edge BVP rays
    (the first and last fan columns) — no extra solve;
  - the inner / outer **arcs** `s = log r_min, log r_max` are the KKG
    `pI2_tritronquee_ic` asymptotic series (v1 corner 1: the inner arc
    is the weaker datum, damped by the harmonic solve and anchored by
    the exact Method-1 voter).

Returns `(rem2, imm2, rem3, imm3, rect)` — the four `Laplace2DSolution`
objects (Re/Im × Method-2/3) and the rectangle.
"""
function surf_laplace_voters(fan)
    rect = surf_sector_rectangle()
    s_lo, s_hi = rect.s_lo, rect.s_hi
    θ_lo, θ_hi = rect.θ_lo, rect.θ_hi

    # Edge-ray data: u on the first / last fan column, as a function of
    # the radial coordinate s = log|x|.  The fan's radii raster is in r;
    # interpolate it in s (its own barycentric ray solution is smooth,
    # but the fan only stores a discrete raster, so a 1-D linear
    # interpolation in s reads it back).
    radii = fan.radii
    s_grid = log.(radii)
    function ray_col_eval(col::AbstractVector{ComplexF64}, s)
        s ≤ s_grid[1]   && return col[1]
        s ≥ s_grid[end] && return col[end]
        k  = searchsortedlast(s_grid, s)
        k  = clamp(k, 1, length(s_grid) - 1)
        tk = (s - s_grid[k]) / (s_grid[k + 1] - s_grid[k])
        return (1 - tk) * col[k] + tk * col[k + 1]
    end
    col_lo = fan.U[:, 1]                      # θ = θ_lo edge
    col_hi = fan.U[:, end]                    # θ = θ_hi edge

    # Arc data: the asymptotic series at fixed |x|, as a function of θ.
    arc(r, θ) = pI2_tritronquee_ic(ComplexF64(r * cis(θ));
                                   t = SURF_T, n_terms = 2)[1]

    # Dirichlet callables.  The rectangle's x-axis is s, y-axis is θ:
    #   bc_x_lo : s = s_lo edge (inner arc),  function of θ
    #   bc_x_hi : s = s_hi edge (outer arc),  function of θ
    #   bc_y_lo : θ = θ_lo edge,              function of s
    #   bc_y_hi : θ = θ_hi edge,              function of s
    re_inner(θ) = real(arc(SURF_R_MIN, θ))
    re_outer(θ) = real(arc(SURF_R_MAX, θ))
    re_elo(s)   = real(ray_col_eval(col_lo, s))
    re_ehi(s)   = real(ray_col_eval(col_hi, s))
    im_inner(θ) = imag(arc(SURF_R_MIN, θ))
    im_outer(θ) = imag(arc(SURF_R_MAX, θ))
    im_elo(s)   = imag(ray_col_eval(col_lo, s))
    im_ehi(s)   = imag(ray_col_eval(col_hi, s))

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
    surf_laplace_eval(lap, z) -> Union{NamedTuple, Nothing}

Evaluate a Laplace voter pair (Method 2 *or* Method 3) at the complex
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
# Region 1 — the per-grid-point majority vote
# ======================================================================

"""
    surf_vote(a, b, c) -> (voted, spread)

The 3-sample majority vote, per component.  `a, b, c` are the three
method estimates of one scalar (`Re u` or `Im u`).  The **median** of
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

# ======================================================================
# Region 2 — the pole-rich wedge (V8b path-network + F2 Stage-2 fill)
# ======================================================================

"""
    surf_wedge_targets() -> Vector{ComplexF64}

The B3 extended threading fan for the wedge walk (ADR-0025 Amendment 2).
The Cartesian product `SURF_TARGET_RADII × SURF_TARGET_ANGLES` — 9 rays
fanned across `±0.5 rad` (`±28.6°`, a `~7°` margin inside the `36°`
Stokes line), `dr = 1` radial spacing from `|x| = 2` out to `|x| = 20`
— `9 × 19 = 171` targets.

Unlike V8b's 20-target `|x| ≤ 8` fan (copied from the pole-scatter
figure), this fan threads the B2 `:max_q_root` adaptive walk through the
*whole* wedge so the walk's node filament passes near every pole out to
`|x| = 20`.  Targets are emitted **radius-major** (inner shells first):
the walk completes the inner field before reaching out, the threading
order the adaptive walk extends along.
"""
surf_wedge_targets() =
    ComplexF64[r * cis(θ) for r in SURF_TARGET_RADII
                          for θ in SURF_TARGET_ANGLES]

"""
    surf_wedge_fill(bvp_sol, grid_pts) -> NamedTuple

Region 2 — the validated pole field + the honest partial `|u|` underlay
(ADR-0025 Amendment 2).

Seed a `VectorPadeTaylorProblem` from the BVP-anchored on-tritronquée
state `bvp_sol(SURF_Z_SEED)`, run the `vector_path_network_solve` wedge
walk over the B3 extended threading fan `surf_wedge_targets()`, and use
the F2 Stage-2 `fine_grid` fill to evaluate `u = V_0` on `grid_pts` (the
Cartesian grid points inside the wedge).

Two honesty-defining choices:

  * the walk runs with the **B2 adaptive `:max_q_root`** policy (the
    package default) — `h = SURF_PN_H` is the step *ceiling*, not a
    hand-tuned fixed step; the controller dodges poles and adapts `h`
    to the local density, threading the fan to `|x| ≈ 18-20` where the
    fixed-`h = 0.1` V8b walk *blocks* past `|x| ≈ 8` (A2 §3.1);

  * the Stage-2 fill runs with **`extrapolate = false`** — the B1
    true-radius gate (ADR-0025 Amendment 1) applies, so a grid point
    outside the nearest node's verified disc gets `NaN`.  **No Padé is
    ever evaluated outside its disc.**  `tol = SURF_PN_TOL` selects the
    gate's `s(tol)` safety factor.

The **pole field** — `extract_poles_shared_q` with cross-node
`min_support ≥ 2` (VC-6) — produces a *candidate* field (~380 poles
with the B3 fan).  VC-6 is a cross-node artefact filter; it does not
catch a Froissart doublet or a mis-located cluster.  So the candidate
field is then passed through the **Phase-D per-pole validation**:

  * **VC-4** (`vc4_validate`, ADR-0025 Amendment 1 §A3) — for each
    candidate fit `u ≈ A·ξ⁻² + B·ξ⁻¹ + C` on a 32-point ring; keep it
    iff `min(|A+1|,|A+3|) < 0.10` (VC-4a, the dominant-balance family)
    **and** `|B| < 0.10·|A|` (VC-4b, the zero residue).  A candidate
    failing either is spurious and **pruned** — `poles` is the
    VC-4-surviving field, the field the figure renders.
  * **VC-5** (`vc5_pair`) — match the VC-4-surviving poles into
    conjugate pairs (`V₀(x̄)=conj V₀(x)`); the pairing residual is an
    FW-style accuracy estimate, a flagged off-axis unpaired pole is
    suspect.  VC-5 is a diagnostic — it *reports*, it does not prune.

Returns `(walk, poles, u, covered, vc4, vc5, message)`:
  - `walk`     : the `VectorPathNetworkSolution`;
  - `poles`    : the **VC-4-validated** pole field — the figure's final
                 field (the FW Fig 4.7 pole-location scatter);
  - `u`        : `Vector{ComplexF64}` of `V_0` at `grid_pts` — finite
                 only where a B1-gated node disc honestly covers,
                 `NaN + NaN·im` everywhere else (no extrapolation);
  - `covered`  : `Vector{Bool}`, `true` iff that grid point is honestly
                 covered (`u[k]` finite) — the honest-coverage mask;
  - `vc4`      : the `vc4_validate` diagnostics — per-candidate `A`,
                 `B`, `dA`, `pass`, `family`, the pruned count and
                 per-pruned-pole failure reason;
  - `vc5`      : the `vc5_pair` diagnostics — conjugate pairs, the
                 pairing residuals, median/max, unpaired + flagged;
  - `message`  : a status note (node count, candidate / validated pole
                 counts, `|x|` frontier, honest-coverage fraction,
                 VC-4 prune count, VC-5 residual).
"""
function surf_wedge_fill(bvp_sol::VectorBVPSolution,
                         grid_pts::AbstractVector{ComplexF64})
    f      = painleve_hierarchy(:I, 2; t = SURF_T)
    y_seed = ComplexF64.(bvp_sol(ComplexF64(SURF_Z_SEED)))
    prob   = VectorPadeTaylorProblem(f, y_seed,
                                     (ComplexF64(SURF_Z_SEED),
                                      ComplexF64(20.0 + 0.0im));
                                     order = SURF_PN_ORDER)
    targets = surf_wedge_targets()
    # B2 adaptive `:max_q_root` walk (the package default); the Stage-2
    # fill is gated honest — `extrapolate = false`, B1 true-radius gate.
    walk = vector_path_network_solve(prob, targets;
                                     order       = SURF_PN_ORDER,
                                     h           = SURF_PN_H,
                                     fine_grid   = grid_pts,
                                     extrapolate = false,
                                     tol         = SURF_PN_TOL)
    candidates = extract_poles_shared_q(walk;
                                        radius_t     = SURF_RADIUS_T,
                                        cluster_atol = SURF_CLUSTER_ATOL,
                                        min_support  = SURF_MIN_SUPPORT)
    # Phase-D per-pole validation: VC-4 prunes the spurious candidates
    # (Froissart doublets / out-of-family / non-zero-residue), VC-5
    # pairs the survivors by conjugate symmetry.  `poles` is the
    # VC-4-validated field — the field the figure renders.
    vc4   = vc4_validate(walk, ComplexF64.(candidates))
    poles = vc4.kept
    vc5   = vc5_pair(poles)
    # The Stage-2 grid_y is ordered as grid_pts; u is its first
    # component.  A `NaN` slot is an honest B1 gate gap (a grid point
    # outside every node's verified disc) — never an extrapolated lie.
    u = ComplexF64[any(isnan, gy) ? ComplexF64(NaN, NaN) : gy[1]
                   for gy in walk.grid_y]
    covered = Bool[!isnan(real(uz)) for uz in u]
    # The |x| frontier the walk's node filament reached.
    frontier = isempty(walk.visited_z) ? 0.0 :
               maximum(abs, walk.visited_z)
    cov_frac = isempty(u) ? 0.0 : count(covered) / length(u)
    msg = "wedge: $(length(walk.visited_z)) visited nodes, " *
          "$(length(candidates)) candidate poles → " *
          "$(length(poles)) VC-4-validated " *
          "($(vc4.n_pruned) pruned), " *
          "|x|-frontier $(round(frontier; digits=1)), " *
          "$(count(covered))/$(length(u)) grid points honestly covered " *
          "($(round(100 * cov_frac; digits=1))%), " *
          "VC-5 conjugate-pair residual median " *
          "$(round(vc5.median_resid; digits=4))"
    return (walk = walk, poles = poles, u = u, covered = covered,
            vc4 = vc4, vc5 = vc5, message = msg)
end

# ======================================================================
# Region 2's negative-axis BVP anchor (V8b Stage A)
# ======================================================================

"""
    surf_anchor_bvp() -> VectorBVPSolution

The negative-real-axis tritronquée BVP that anchors the wedge
path-network seed — V8b's Stage-A recipe (`[-20,-2]` segment, 2+2
`(u,u')` boundary condition).  Reused so the wedge walk starts from a
genuine on-tritronquée state at `SURF_Z_SEED = -3`.
"""
function surf_anchor_bvp()
    f, Jf = surf_companion()
    CT = ComplexF64
    x_l = CT(-20.0); x_r = CT(-2.0)
    Ba = zeros(CT, SURF_D, SURF_D)
    Bb = zeros(CT, SURF_D, SURF_D)
    Ba[1, 1] = 1; Ba[2, 2] = 1
    Bb[3, 1] = 1; Bb[4, 2] = 1
    seed_l = pI2_tritronquee_ic(x_l; t = SURF_T, n_terms = 2)
    seed_r = pI2_tritronquee_ic(x_r; t = SURF_T, n_terms = 2)
    g = CT[seed_l[1], seed_l[2], seed_r[1], seed_r[2]]
    return vector_bvp_solve(f, x_l, x_r, Ba, Bb, g;
                            N        = 128,
                            tol      = SURF_BVP_TOL,
                            maxiter  = SURF_BVP_MAXITER,
                            jacobian = Jf,
                            initial_guess =
                                z -> pI2_tritronquee_ic(z; t = SURF_T,
                                                        n_terms = 2))
end

# ======================================================================
# Top-level driver — the whole-plane stitch
# ======================================================================

"""
    kkg_pi2_surface() -> NamedTuple

Assemble `Re u` / `Im u` of the `P_I^(2)` tritronquée `V_0(x, 0)` over
the Cartesian square `[-20,20]²`.

Construction (see the file docstring):
  1. **Region 1 — the pole-free sector** is computed three ways
     (ray-fan BVP / in-house 2D-Chebyshev Laplace / Gridap FEM Laplace)
     and **majority-voted per grid point** (median); an agreement /
     spread map records the max pairwise disagreement.
  2. **Region 2 — the pole-rich wedge** (ADR-0025 Amendment 2) is the
     validated pole *field* (`extract_poles_shared_q`, the primary
     deliverable) plus an honest partial `|u|` surface underlay — the
     B2 adaptive walk + the B1-gated Stage-2 fill, finite only where the
     verified node discs cover (~13-18 % of the wedge), `NaN` elsewhere.
  3. The two regions are stitched onto one grid; the `±36°` Stokes-line
     strips are `NaN`-masked.

Returns a `NamedTuple`:
  - `xs`, `ys`            : the `SURF_GRID_N`-point axes of `[-20,20]`;
  - `Re_u`, `Im_u`        : `SURF_GRID_N × SURF_GRID_N` matrices,
                            `Re_u[i,j] = Re V_0(xs[i] + im·ys[j])`
                            (`NaN` outside the disc `|x| ≤ 20`, in the
                            masked strips, in the wedge wherever the B1
                            gate found no honest datum, and where no
                            sector method has a datum);
  - `spread`              : the agreement map — per sector grid point,
                            the max pairwise disagreement of the three
                            methods (`max(Re-spread, Im-spread)`);
  - `poles`               : the **VC-4-validated** wedge pole field —
                            `extract_poles_shared_q` candidates
                            (`min_support ≥ 2`, VC-6) pruned by the
                            VC-4 dominant-balance per-pole certificate
                            (ADR-0025 Amendment 1 §A3) — the
                            Amendment-2 primary wedge deliverable (the
                            FW Fig 4.7 pole locations);
  - `vc4`                 : the VC-4 diagnostics — per-candidate `A`,
                            `B`, `dA`, `pass`, `family`, the pruned
                            count and per-pruned-pole failure reason
                            (`nothing` when the wedge has no points);
  - `vc5`                 : the VC-5 conjugate-pairing diagnostics —
                            the pairs, the residuals, median/max, the
                            unpaired and flagged-suspect poles
                            (`nothing` when the wedge has no points);
  - `wedge_walk`          : the Region-2 `VectorPathNetworkSolution`
                            (`nothing` when the wedge has no points) —
                            the shared-`Q` path-network the VC-4 ring
                            fit and the pole extraction both read, so a
                            test can re-exercise VC-4 against it;
  - `sector_method`       : per sector grid point, the per-method
                            `(re1, im1, re2, im2, re3, im3)` estimates —
                            so a test can inspect each voter (a
                            `Dict` keyed by `(i, j)`);
  - `wedge_covered`       : `SURF_GRID_N × SURF_GRID_N` `Bool` matrix,
                            `true` exactly where Region 2's honest
                            partial surface underlay carries a B1-gated
                            datum (the honest-coverage mask — covered
                            cells are finite in `Re_u`/`Im_u`,
                            uncovered cells are `NaN`).  No cell is
                            ever an extrapolation (ADR-0025 Amendment
                            1/2: `extrapolate = false`);
  - `ray_fan`             : the Method-1 polar grid (for inspection);
  - `message`             : the Region-2 status note (node count, pole
                            count, `|x|` frontier, coverage fraction).

Deterministic — a fixed recipe; two calls return identical matrices.
"""
function kkg_pi2_surface()
    xs = collect(range(-SURF_XY_LIM, SURF_XY_LIM; length = SURF_GRID_N))
    ys = collect(range(-SURF_XY_LIM, SURF_XY_LIM; length = SURF_GRID_N))
    n  = SURF_GRID_N

    Re_u   = fill(NaN, n, n)
    Im_u   = fill(NaN, n, n)
    spread = fill(NaN, n, n)
    wedge_covered = fill(false, n, n)
    sector_method = Dict{Tuple{Int,Int},
                         NTuple{6,Float64}}()

    # ---- Region 1: the three sector voters --------------------------------
    fan  = surf_ray_fan()
    lap  = surf_laplace_voters(fan)

    # ---- Region 2: the wedge BVP anchor + path-network fill ---------------
    anchor = surf_anchor_bvp()
    # Collect the wedge grid points (inside the disc, in the wedge, not in
    # the masked strip) into a flat list for the Stage-2 fine_grid call.
    wedge_idx = Tuple{Int,Int}[]
    wedge_pts = ComplexF64[]
    for j in 1:n, i in 1:n
        z = ComplexF64(xs[i], ys[j])
        abs(z) > SURF_XY_LIM       && continue        # outside the disc
        surf_in_mask(z)            && continue        # boundary strip
        surf_in_sector(z)          && continue        # not the wedge
        iszero(z)                  && continue        # origin → leave NaN
        push!(wedge_idx, (i, j))
        push!(wedge_pts, z)
    end
    wedge = isempty(wedge_pts) ?
        (walk = nothing, poles = ComplexF64[], u = ComplexF64[],
         covered = Bool[], vc4 = nothing, vc5 = nothing,
         message = "wedge: no grid points") :
        surf_wedge_fill(anchor, wedge_pts)

    # ---- assemble the grid ------------------------------------------------
    for j in 1:n, i in 1:n
        z = ComplexF64(xs[i], ys[j])
        # Outside the |x| ≤ 20 disc: leave NaN (the figure is the disc).
        abs(z) > SURF_XY_LIM && continue
        # The masked Stokes-line strips: leave NaN (v1 corner 4).
        surf_in_mask(z) && continue

        if surf_in_sector(z)
            # Region 1 — the majority vote.
            m1   = surf_ray_eval(fan, z)
            l2   = surf_laplace_eval_one(lap.rem2, lap.imm2, lap.rect, z)
            l3   = surf_laplace_eval_one(lap.rem3, lap.imm3, lap.rect, z)
            re1, im1 = real(m1), imag(m1)
            re2, im2 = l2 === nothing ? (NaN, NaN) : (l2.re, l2.im)
            re3, im3 = l3 === nothing ? (NaN, NaN) : (l3.re, l3.im)
            vre, sre = surf_vote(re1, re2, re3)
            vim, sim = surf_vote(im1, im2, im3)
            Re_u[i, j]   = vre
            Im_u[i, j]   = vim
            spread[i, j] = max(sre, sim)
            sector_method[(i, j)] =
                (re1, im1, re2, im2, re3, im3)
        end
    end

    # ---- write Region 2's honest partial underlay into the grid -----------
    # Only the B1-gated covered cells carry a datum; everywhere else in
    # the wedge stays `NaN` — an honest gap, never an extrapolation
    # (ADR-0025 Amendment 1/2: the Stage-2 fill ran `extrapolate=false`).
    if wedge.walk !== nothing
        for (k, (i, j)) in enumerate(wedge_idx)
            uz = wedge.u[k]
            if isnan(real(uz))
                # B1 gate found no honest datum here — leave NaN, honest.
                continue
            end
            Re_u[i, j] = real(uz)
            Im_u[i, j] = imag(uz)
            wedge_covered[i, j] = true
        end
    end

    return (xs = xs, ys = ys, Re_u = Re_u, Im_u = Im_u,
            spread = spread, poles = wedge.poles,
            vc4 = wedge.vc4, vc5 = wedge.vc5,
            wedge_walk = wedge.walk,
            sector_method = sector_method,
            wedge_covered = wedge_covered,
            ray_fan = fan, message = wedge.message)
end
