# figures/_kkg_pi2_surface_helpers.jl
#
# F3 — the whole-plane compute kernel for the KKG 2015 P_I^(2)
# tritronquée surface figure (bead `padetaylor-0ln.33`, v0.2 plan).
# Makie-free: this file assembles `Re u` and `Im u` of the tritronquée
# `V_0(x, 0)` over the complex `x`-plane square `[-4,4]²` and returns
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
# The three sector voters live in the sibling helper
# `figures/_kkg_pi2_sector.jl` (the Rule-6 file-size split); its
# top-of-file docstring is the chapter on the triple-method pipeline.
# In brief:
#
# **Voter 1 — the ray-fan BVP (the direct voter).**  A fan of radial
# rays `x = ξ·e^{iφ}`, `φ` stepping across the sector.  Each ray is a
# `vector_bvp_solve` of the `d = 4` companion `P_I^(2)` system,
# boundary-pinned at both ends by the F1-generalised
# `pI2_tritronquee_ic`.  **C2 (ADR-0025 Amendment 10)** widens each
# ray's segment inward to `[R_MAX, R_INNER_BC]`, `R_INNER_BC < R_MIN`,
# so the inner arc `|x| = R_MIN` is a BVP *interior* point — its value
# is the genuine ODE solution, not the pinned asymptotic seed.
# **C3 (ADR-0025 Amendment 10)** reconstructs voter 1 onto the Cartesian
# grid by exact-radius barycentric evaluation of the stored ray
# `VectorBVPSolution`s + a Catmull-Rom cubic angular blend across four
# bracketing rays — retiring the lossy bilinear polar-raster scheme.
#
# **Voter 2 — the in-house 2D-Chebyshev Laplace solve.**  A pole-free
# *angular sector* `{r₀ ≤ |x| ≤ r₁, θ₀ ≤ arg x ≤ θ₁}` is **conformally
# a rectangle** under `w = log x`: `Re w = log|x| ∈ [log r₀, log r₁]`,
# `Im w = arg x ∈ [θ₀, θ₁]`.  The Laplacian is conformally invariant,
# so the sector Laplace problem IS a rectangle Laplace problem in `w`
# (ADR-0024 Context).  `PadeTaylor.laplace2d_solve` solves it on the
# rectangle with Dirichlet data.  **C2** sources *all four* rectangle
# edges from the ray fan: the side rays as before, and now the inner /
# outer arcs from the fan's innermost / outermost interior radial
# samples — *not* the `pI2_tritronquee_ic` asymptotic series (v1
# corner C1 retired).  Run once for `Re u`, once for `Im u`.
#
# **Voter 3 — the Gridap FEM Laplace solve.**  `laplace2d_solve_gridap`
# (`figures/_kkg_pi2_gridap_helper.jl`), the same rectangle + the same
# (now BVP-sourced) Dirichlet data, via an algorithmically-disjoint
# finite-element discretisation — a FEM bug and a spectral bug cannot
# coincide.
#
# **The vote.**  At every Cartesian grid point inside the sector the
# three estimates are reduced per component by the **median** (the
# 3-sample median = the majority vote: an outlier voter is discarded,
# the two that concur win).  An **agreement / spread map** records the
# maximum pairwise disagreement per point — a figure-quality diagnostic
# of where the three methods concur (and, honestly, where they do not).
#
# ### Region 2 — the pole-rich wedge (ADR-0025 Amendment 14 / ADR-0026 Amendment 10)
#
# The headline figure's wedge renders **FW-faithful filled** with
# **honest provenance marking** — the ADR-0025 Amendment 14 / ADR-0026
# Amendment 10 dual-fill design (`tf9`).  The earlier "honest partial
# underlay only, ~5–22 % of the wedge `NaN` elsewhere" rendering
# (Amendment 2 → Amendment 9) shipped technically-honest but visually-
# inadequate figures — a ~95 %-blank wedge that did not visualise the
# tritronquée solution at all.  The D1–S8 investigation arc (ADR-0026)
# proved the no-extrapolation contract is *stricter than FW 2011 itself*:
# FW's published pole-field figures FILL precisely because FW's Stage-2
# evaluates each Padé over its full step with no validity gate
# (FW-md:395-397).  Amendment 14 keeps the strict contract as the
# package default for every non-headline figure and the public solver
# API, and grants the headline figure a per-figure override — render
# FW-style filled, mark provenance.
#
# Concretely Region 2 is:
#
#   1. **Pole field — the primary deliverable.**  A dense Cartesian
#      lattice (`surf_wedge_dense_targets`) drives the resilient
#      `:max_q_root` adaptive walk through the *whole* wedge out to
#      `|x| = 4`, so the walk's nodes pass near every pole;
#      `extract_poles_shared_q` with cross-node `min_support ≥ 2`
#      (VC-6) roots the per-node shared `Q` and clusters the roots into
#      the validated pole field.  Order 48 is the measured saturation
#      knee (ADR-0026 Amendment 5).
#
#   2. **Dual-fill `|u|` surface.**  One walk, two Stage-2 evaluations:
#      a **certified fill** (`extrapolate = false`, B1 true-radius
#      gate, ~29 % of the wedge at order 48 — the measured Amendment 9
#      figure) and an **FW-style filled overlay** (`extrapolate = true`,
#      every cell finite, identical to FW 2011's published rendering).
#      The walk is shared (the expensive part); the second fill is a
#      Horner sweep against the same `visited_*` tree.  The figure
#      script composites both with the `wedge_covered` mask: certified
#      cells render at full opacity, extrapolated cells at reduced
#      alpha.  A reader can tell certified from extrapolated by the
#      alpha channel — the figure SHOWS the full surface but LABELS
#      its provenance.
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
# `figures/_kkg_pi2_vc45.jl`; VC-8 (the BVP endpoint companion-
# consistency certificate, ADR-0025 Amendment 8) in
# `figures/_kkg_pi2_vc8.jl`; VC-10 (the FW/FFW two-run pole-field
# accuracy indicator, ADR-0025 Amendment 9) in `figures/_kkg_pi2_vc10.jl`
# — all CLAUDE.md Rule 6 file-size splits.  VC-7 (loop closure) remains
# a separate Phase-D bead.
#
# ### The stitch
#
# Region 1 (voted) + Region 2 are written onto one `[-4,4]²`
# Cartesian grid for `Re u` and `Im u`.  A thin `±1°` strip straddling
# the sector/wedge boundary (`||arg x| − 36°| < SURF_STITCH_MASK_DEG`)
# is masked with `NaN`: it is the band the sector fan does not reach
# (`arg x < 36° + margin`) and the wedge B1-gated underlay does not
# honestly cover.  ADR-0025 Amendment 12's Phase-E1 audition narrowed
# this strip `±3° → ±1°` (v1 corner C4 — see corner 4 below).
#
# ## v1 corners (CLAUDE.md Rule 9 — each with its forcing condition)
#
#   1. **RETIRED — inner-arc asymptotic seed (C1).**  Voters 2/3's inner
#      arc (`|x| = r_min ≈ 2`) was the KKG `n_terms = 2` asymptotic
#      series — an *asymptotic* (large-`|x|`) expansion, `O(10⁻³·10⁻²)`
#      at `|x| ≈ 2`.  ADR-0025 Amendment 10's C2 retires it: each ray
#      BVP now extends inward to `R_INNER_BC = 1.05` (`surf_ray_bvp`),
#      so the inner arc is a BVP *interior* datum (`~1·10⁻⁴` accurate —
#      the C2 probe), and voters 2/3 read the inner / outer arcs from
#      the ray fan instead of the series (`surf_laplace_voters`).  The
#      C2 instinct "just read `fan.U[1,:]`" was *probe-falsified* — the
#      production 2+2 split BC pins the BVP inner endpoint to the seed
#      exactly — so the genuine fix is the *extended-inward* segment,
#      not a re-route of the same pinned data.  Detail:
#      `figures/_kkg_pi2_sector.jl`, `external/probes/c2-arc-data/`.
#
#   2. **RETIRED (then narrowed-with-provenance) — wedge Stage-2
#      extrapolation.**  The shipped V8b figure called the Stage-2 fill
#      with `extrapolate = true`, evaluating the shared-`Q` Padé
#      approximants far outside their verified disc (the A4 baseline
#      probe measured midpoints 50× outside the disc; ~70 % of loop
#      closures catastrophic).  ADR-0025 Amendment 1's B1 true-radius
#      gate retired the blanket `extrapolate = true`.  ADR-0025
#      Amendment 14 (this figure's `tf9` closeout) then re-introduces
#      `extrapolate = true` as a **per-figure, provenance-marked
#      overlay** — `surf_wedge_fill` runs the Stage-2 fill twice on the
#      same walk (once `extrapolate = false` for the certified core,
#      once `extrapolate = true` for the FW-2011-style filled overlay),
#      and the figure script paints certified cells at full opacity and
#      extrapolated cells at reduced alpha so the reader can read
#      provenance off the alpha channel.  The strict B1 gate remains the
#      *package default* and applies to every non-headline figure.
#
#   3. **RETIRED — hand-tuned wedge step `h = 0.1`.**  ADR-0025
#      Amendment 4's B2 adaptive `:max_q_root` walk (the package
#      default) retires the hand-tune: the step direction dodges the
#      nearest shared-`Q` pole and `h` adapts to the local pole density.
#      The fixed-`h = 0.1` V8b walk *blocks* the extended fan past
#      `|x| ≈ 8` (A2 §3.1); the B2 walk threads it to `|x| ≈ 18-20`.
#
#   4. **NARROWED (C4) — boundary-strip masking.**  ADR-0025
#      Amendment 12's Phase-E1 audition (bead `padetaylor-0ln.37.18`)
#      retires the conservative cushion: `SURF_SECTOR_MARGIN_DEG` and
#      `SURF_STITCH_MASK_DEG` are narrowed `4° → 1°` and `3° → 1°`.  The
#      audition (`external/probes/stokes-strip-audition/`) measured the
#      triple-method vote spread + an independent dedicated-BVP oracle
#      error in 1°-wide bands marching toward the 36° Stokes line: the
#      sector ODE solve does NOT degrade at the line (Newton residual
#      `1.7·10⁻¹¹`, companion-consistency `1.8·10⁻¹⁰` *at* `36°`), and
#      pushed to a `1°` margin the near-Stokes band `[37°,38°]` stays
#      inside the figure's bulk honesty envelope (vote spread max
#      `3.3·10⁻³`, oracle error max `2.2·10⁻³` — both below the
#      already-accepted `~6.7·10⁻³` inner-arc seed floor).  The grey
#      Stokes strip shrinks from `6°` (`[33°,39°]`) to `2°`
#      (`[35°,37°]`) — a 3× reduction.  *Forcing condition for the
#      residual `±1°` mask:* the sector fan is a *sector* solver and the
#      wedge a disjoint region — neither can straddle the Stokes line
#      itself; a uniform connection formula across the line would close
#      the last `2°` (Amendment-12 deferred bead).
#
#   5. **RETIRED — bilinear ray-fan voter (C5).**  Voter 1 sampled the
#      ray fan onto a `40`-row polar raster and *bilinearly*
#      interpolated it onto the Cartesian grid — interpolation error
#      `~1·10⁻³` at the `6°` fan resolution (the C3 audition,
#      `external/probes/c2-arc-data/c3_audition.jl`), larger than the
#      bulk triple-method spread.  ADR-0025 Amendment 10's C3 retires
#      it: `surf_ray_eval` now evaluates the stored ray
#      `VectorBVPSolution`s exactly in the radius (barycentric) and
#      blends four bracketing rays by a Catmull-Rom cubic in the angle
#      (`~7·10⁻⁵` — a ~19× cut).  A *harmonic* voter 1 was rejected by
#      the audition — it would be `≡ voter 2` and collapse the ADR-0024
#      triple-method independence.
#
#   6. **RETIRED — the `121²` grid / `6°` ray fan (C6).**  ADR-0025
#      Amendment 11's C1 (bead `padetaylor-0ln.37.9`) lifts the
#      Cartesian raster to `SURF_GRID_N = 401` and the ray fan to
#      `SURF_RAY_DPHI_DEG = 2°` (≈ 141 rays) — a genuinely sharp
#      headline figure, still inside the Phase-F runtime gate.  The
#      coupled inner-arc-spread job assigned to the same bead found
#      (CLAUDE.md Rule 2) that the inner-arc spread is *not* a Laplace
#      under-resolution: it is the asymptotic-seed truncation floor
#      (~5·10⁻³ at `|x| = 2`), irreducible by any `N` — see the
#      `SURF_BVP_N` / `SURF_LAP_NX` comments and ADR-0025 Amendment 11.
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
#     walk, **Amendment 10** the C2 BVP-sourced arc data + C3 voter-1
#     reconstruction.  This kernel's Region 2 is the B3 deliverable.
#   * `figures/_kkg_pi2_sector.jl` — Region 1: the three sector voters
#     and the median vote (the Rule-6 split; C2/C3 live there).
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

# `Printf` for the verbose-flush kernel phase markers inside
# `kkg_pi2_surface` (memory `long-julia-agents-verbose-flush`).  This
# file is `include`d inside the `KKGSurfaceKernel` module of the render
# script, which has no Printf in scope; declare it here so the helpers
# remain self-contained.
using Printf

# VC-4 (dominant-balance per-pole certificate) and VC-5 (conjugate-
# symmetry pairing) — the Phase-D per-pole validation that certifies the
# genuine wedge poles and prunes the spurious ones (ADR-0025 Amendment 1
# §A3; beads `padetaylor-0ln.37.12` / `0ln.37.13`).
if !isdefined(@__MODULE__, :vc4_validate)
    include(joinpath(@__DIR__, "_kkg_pi2_vc45.jl"))
end

# VC-8 (BVP endpoint higher-derivative match, FW 2011 §5.2) — the
# Phase-D companion-consistency certificate of the figure's BVP solves
# (ADR-0025 Amendment 8; bead `padetaylor-0ln.37.15`).  Sibling helper,
# the Rule-6 file-size split, mirroring `_kkg_pi2_vc45.jl`.
if !isdefined(@__MODULE__, :surf_vc8_companion_check)
    include(joinpath(@__DIR__, "_kkg_pi2_vc8.jl"))
end

# Region 1 — the pole-free-sector triple-method pipeline (ADR-0024's
# three voters + the median vote).  Sibling helper, the Rule-6
# file-size split; ADR-0025 Amendment 10 (C2 BVP-sourced arc data + C3
# voter-1 reconstruction) lives here.  The `include` is placed AFTER the
# figure-defining constants below, so the sector code sees `SURF_R_MIN`,
# `SURF_BVP_N`, … — see the `include` call further down.

# ======================================================================
# Figure-defining constants — fixed here, the kernel is a pure function
# of them.
# ======================================================================

# The `P_I^(2)` parameter.  The KKG tritronquée surface figure is t = 0.
const SURF_T = 0.0

# Companion-system order: d = 2m = 4 for P_I^(2) (m = 2).
const SURF_D = 4

# The Cartesian render grid: [-4,4]² (KKG Figs 7.4/7.5, headline-figure
# wedge shrink — bead `padetaylor-apn`).
#
# The shipped headline figure rendered the wedge at `SURF_XY_LIM = 20.0`,
# which yielded a 3528-pole wedge field — visual noise, too many for a
# headline figure.  Shrinking to `R = 4.0` gives a wedge area scaling as
# `R²-4` (the inner arc is `|x| ≥ 2`), so the pole field drops by a
# factor of `(16-4)/(400-4) = 12/396 ≈ 1/33` to ~100–300 poles —
# legible for a headline figure.  This is a **figure-config change
# only**: every solver knob (Taylor order, walk step ceiling,
# scale-covariant tolerances) stays unchanged, and the kernel is
# scale-covariant by construction (`SURF_RADIUS_T` is dimensionless,
# `h_max = SURF_PN_H` independent of R).
#
# C1 (ADR-0025 Amendment 11, bead `padetaylor-0ln.37.9`) — retires v1
# corner C6.  The shipped figure rendered a `121²` raster; that was
# resolution-limited (worklog 057).  `401` lifts it to a genuinely sharp
# raster — odd so a node still sits on `x = 0` and on `x = -3` (the
# PI2S.1 ridge-check pin, originally `x = -15` at `R = 20`).  Voters
# 2/3 are spectrally callable at any point and voter 1 (post-C3) is
# exact-in-radius + cubic-angular, so the finer raster is mostly cheap
# `O(N²)` evaluation; the C1 runtime probe (`bead 0ln.37.9`) measured
# the full `kkg_pi2_surface()` comfortably inside the Phase-F runtime
# gate at `401²`.
const SURF_XY_LIM = 4.0
const SURF_GRID_N = 401          # 401×401 — odd so a node sits on 0.

# --- Region 1: the pole-free sector --------------------------------------
# Pole-free sector half-width 6π/7 − (3/7)arctan(1/√5) ≈ 144°, centred on
# the negative real axis (arg x = 180°).  So the sector is
# arg x ∈ (36°, 324°); the wedge is |arg x| < 36°.
const SURF_WEDGE_HALF_DEG = 36.0
# The ray fan stays `SURF_SECTOR_MARGIN_DEG` inside each Stokes line, so
# the sector solver covers `arg x > 36° + margin`.
#
# C4 (ADR-0025 Amendment 12, bead `padetaylor-0ln.37.18`) — the Phase-E1
# Stokes-strip audition.  The shipped figure used `margin = 4°` (a
# conservative cushion).  The audition (`external/probes/stokes-strip-
# audition/`) measured the triple-method vote spread + an independent
# dedicated-BVP oracle error in 1°-wide angular bands marching toward the
# 36° Stokes line, for a margin sweep `4° → 0.5°`.  Two findings forced
# the narrowing:
#
#   * the sector ODE solve does NOT degrade at the Stokes line — a
#     dedicated through-the-point BVP at `arg x = 36.0°` has Newton
#     residual `1.7·10⁻¹¹` and companion-consistency `1.8·10⁻¹⁰`, as
#     healthy as a mid-sector ray.  The pole-free sector solution is
#     analytic (harmonic) right up to the Stokes line; the 4° inset was
#     a cushion, not a numerical limit;
#   * pushed to `margin = 1°` (fan edge at `37°`), the near-Stokes band
#     `[37°,38°]` has vote spread median `1.8·10⁻³` / max `3.3·10⁻³` and
#     oracle error median `6.6·10⁻⁴` / max `2.2·10⁻³` — comfortably
#     inside the figure's bulk honesty envelope (PI2S.2: median `< 3·10⁻³`,
#     max `< 1·10⁻²`) and below the inner-arc seed floor (`~6.7·10⁻³`,
#     Amendment 11) the figure already accepts as honest.
#
# `margin = 1°` is the honest frontier: `margin = 0.5°` (fan edge `36.5°`)
# was probed and rejected — its `[36.5°,37°]` band oracle error climbs to
# max `4.1·10⁻³` (still inside the `< 1·10⁻²` envelope, but the C3
# voter-1 boundary-cell linear fallback is steeply 2°-under-sampled there;
# the senior-grade call keeps `1°` of headroom).  The mask narrows in
# lockstep (`SURF_STITCH_MASK_DEG`).
const SURF_SECTOR_MARGIN_DEG = 1.0

# Ray-fan radii.  r_min ≈ 2 (the inner arc of the conformal rectangle),
# r_max = 4 (the grid corner radius is 4√2 ≈ 5.66, but the sector pins
# are set at 4 and the corners outside the disc |x| ≤ 4 are
# NaN-masked — see SURF_R_MAX use).  `r_max` is the wedge outer radius
# (bead `padetaylor-apn`); shrinking from 20 → 4 cuts the wedge area
# by ~33× and so cuts the pole field from ~3528 to ~100–300 poles.
const SURF_R_MIN = 2.0
const SURF_R_MAX = 4.0

# C2 (ADR-0025 Amendment 10) — the ray BVP's *inner boundary-condition*
# radius.  Each sector ray BVP solves over `[SURF_R_MAX, SURF_R_INNER_BC]`
# with `SURF_R_INNER_BC < SURF_R_MIN`, so the inner arc `|x| = SURF_R_MIN`
# is an *interior* collocation point — its value is the genuine
# ODE-constrained solution, not the asymptotic seed that pins the
# endpoint.  The C2-exploration probe (`external/probes/c2-arc-data/`)
# measured `1.05` recovers the inner-arc datum to `~1·10⁻⁴` (vs the
# seed's `~5·10⁻⁴` directly at the arc) and that every fan ray BVP
# still converges; `pI2_tritronquee_ic` floors at `|x| ≥ 1`, so `1.05`
# keeps a margin.  The boundary-layer seed error lives in
# `[SURF_R_INNER_BC, SURF_R_MIN]`, below the rendered window.
const SURF_R_INNER_BC = 1.05

# Angular step of the ray fan.  The sector spans ≈ 280° usable.
#
# C1 (ADR-0025 Amendment 11, bead `padetaylor-0ln.37.9`) — retires v1
# corner C6's `6°` fan.  `2°` (≈ 141 rays) is the runtime knee: the C1
# boundary-data probe (`bead 0ln.37.9`) measured the cubic-angular
# `surf_ray_eval` inner-arc reconstruction error fall from `~5·10⁻⁴` at
# `6°` to `~2·10⁻⁸` at `2°` — far below the asymptotic-seed accuracy
# floor (see `SURF_BVP_N` below), so the ray fan ceases to be an
# accuracy contributor at `2°`.  A still-finer `1°` fan was probed and
# rejected: it costs ~8 s more (one BVP per ray) for *zero* measurable
# gain — the residual error is the seed floor, not the fan.
const SURF_RAY_DPHI_DEG = 2.0

# BVP discretisation per ray.  N+1 Chebyshev nodes.  C2 widened the ray
# segment to `[SURF_R_MAX, SURF_R_INNER_BC]` (ADR-0025 Amendment 10);
# the wider segment has more curvature near its inner end, so `N = 128`
# (the C2 probe: companion-consistency `~1·10⁻¹¹` at `N = 128` for
# `[20, 1.05]`, vs a marginal `~1·10⁻⁷` at `N = 96`).  At the
# headline-figure `SURF_R_MAX = 4` (bead `padetaylor-apn`) the segment
# `[4, 1.05]` is *shorter* and so even less curvature-demanding;
# `N = 128` is well past converged for that shorter span and is
# retained unchanged.
#
# C1 (ADR-0025 Amendment 11) — `N = 128` is *already converged*: the C1
# probe measured the ray BVP interior value at `N = 128` agrees with an
# `N = 320` solve to `~1·10⁻¹²`.  The residual inner-arc error is NOT
# an `N` deficit.  It is the **asymptotic-seed truncation floor**: each
# ray BVP pins `(u, u')` at the inner endpoint `R_INNER_BC = 1.05` to
# the KKG `n_terms = 2` series, whose `O(|x|^{-7/3})` truncation error
# at `|x| ≈ 1` propagates inward as a boundary layer.  Solving the same
# ray at several `r_in ∈ [1.05, 1.8]` (each a legitimate seed-pinned
# BVP) gives interior values that disagree by **~5·10⁻³ at `|x| = 2`,
# ~4·10⁻³ at 2.5, ~7·10⁻⁴ at 5, ~1.7·10⁻⁵ at 10** — this `r_in`-spread
# *is* the seed floor, and it is the irreducible source of the
# inner-arc vote spread.  `pI2_tritronquee_ic` implements only
# `n_terms ∈ {1, 2}`, so the seed cannot be sharpened here; closing the
# inner-arc floor needs `n_terms ≥ 3` (KKG 2015 eq. 7.2 `c₇…`),
# recorded as a deferred bead in Amendment 11.
const SURF_BVP_N       = 128
const SURF_BVP_TOL     = 1.0e-9
const SURF_BVP_MAXITER = 40

# The two Laplace voters' rectangle resolution (in conformal w = log x
# coordinates).  Spectral (Method 2): a handful of nodes suffice.  FEM
# (Method 3): O(h²), needs a finer mesh.
#
# C1 (ADR-0025 Amendment 11, bead `padetaylor-0ln.37.9`) — the deep
# finding.  ADR-0025 Amendment 10 hypothesised the inner-arc vote spread
# (~3.9·10⁻³) was a voter-2/3 Laplace *under-resolution* and assigned a
# higher-`N` fix to this bead.  The C1 convergence probe **falsified
# that hypothesis** (CLAUDE.md Rule 2 — root cause, not band-aid): both
# Laplace voters *already converge geometrically* — the spectral solve
# does not move from `Nx = 40` to `Nx = 120` (it stays pinned ~4.8·10⁻³
# off a dedicated-BVP oracle), and the FEM voter plateaus by
# `n_x = 64`.  Solving the spectral problem with *perfect*
# dedicated-BVP edge traces leaves the same ~4.8·10⁻³ interior error,
# N-independent.  The genuine root cause is the **asymptotic-seed
# truncation floor** (see `SURF_BVP_N`): the inner-arc vote spread is
# irreducible by *any* Laplace `N`.  So these constants are set to the
# measured spectral/FEM *convergence knee*, NOT raised speculatively —
# raising `SURF_LAP_NX` past ~48 in fact makes the spectral solve
# *noisier* (it over-resolves the seed-floored, slightly-non-harmonic
# edge data: self-error 2.2·10⁻⁵ at Nx=48, rising to 3.1·10⁻⁴ at
# Nx=64).  `48/64` is the cleanest spectral knee; `64/88` is where the
# FEM voter has plateaued.  Amendment 11 records the falsification and
# the quantified floor in full.
const SURF_LAP_NX = 48           # radial (s = log|x|) Chebyshev nodes
const SURF_LAP_NY = 64           # angular (θ) Chebyshev nodes
const SURF_FEM_NX = 64           # FEM mesh subdivisions, radial
const SURF_FEM_NY = 88           # FEM mesh subdivisions, angular

# --- Region 2: the pole-rich wedge (ADR-0025 Amendment 2) ----------------
# The wedge deliverable is the validated pole FIELD + an honest partial
# `|u|` underlay (NOT a filled surface — A2 proved that unreachable).
#
# Seed: a BVP-anchored on-tritronquée state at `z = -3` (V8b Stage A).
const SURF_Z_SEED  = -3.0 + 0.0im

# Taylor order of the path-network jets.
#
# **`tf9` (ADR-0026 Amendment 10, ADR-0025 Amendment 14) — raised 36 → 48,
# the measured order-saturation knee.**  The Amendment-9 converged finding
# (the D1–S8 investigation arc) is that the honest B1 disc is **100 %
# truncation-limited**; the *only* walk-side lever on coverage is the
# Taylor order.  ADR-0026 Amendment 5's S5-diag order sweep measured the
# certified coverage knee directly:
#
#     order 24 → 36 → 48 → 72   :   coverage 13 % → 23 % → 29 % → 29 %
#
# Order 48 is the **measured saturation knee** — the biggest honest
# certified core (~29 %) at ~2× the order-36 runtime; order 72 is wasted
# (coverage flat).  The order-48 *walk*'s pre-resilient A2 §5b.5
# regression (the fixed-`h` shared-`Q` SVD degeneration that motivated the
# earlier "order 48 not viable" claim) is **superseded by ADR-0026
# Amendment 5**: the resilient adaptive stack (`adaptive = true`,
# `on_target_failure = :skip`, `skip_covered = true`) runs order 48
# cleanly with ~0 failures.
#
# The shipped figure renders the wedge in **dual-fill provenance form**
# (ADR-0025 Amendment 14 / ADR-0026 Amendment 10):
#
#   * a **certified core** — the B1-gate-honest cells, ~29 % of the
#     wedge at order 48 — rendered at full opacity.  This is the senior-
#     grade region: no Padé evaluated outside its truncation-bounded
#     disc;
#
#   * an **FW-2011-style filled overlay** — every Padé evaluated over
#     its full step with no validity gate, exactly the rendering FW's
#     published pole-field figures use (FW-md:395-397).  Filled to the
#     wedge boundary; rendered at reduced opacity to mark it as
#     extrapolated.
#
# The certified-coverage progression across the D1–S8 arc:
# `~5 % → 8 % → 22 % → 29 %` at order 24 → corrected stack → order 36 →
# order 48.  Beyond order 48 the certified-core ceiling is the
# no-extrapolation honesty contract itself (ADR-0025), not the Taylor
# order; closing the remaining ~71 % visually is the FW-faithful
# extrapolated overlay, marked accordingly.
#
# The companion `SAFETY` lowering (Amendment 5's S6b) was *reverted* —
# it regresses the `src/` pole-extraction tests through a `radius_t`/`h`
# coupling (see the `const SAFETY` comment in `src/VectorWedgeStep.jl`);
# the wedge walk therefore runs at the verified `SAFETY = 0.25`.
const SURF_PN_ORDER = 48

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
# and accurately to `|x| = SURF_R_MAX`".  The fan is designed so the B2
# `:max_q_root` adaptive walk threads a node filament past every pole.
#
# Design (A2 §5.2 recommends a 9-angle `dr = 1` fan as the densest that
# threads; B2's adaptive walk threads it across the wedge — Amendment 4):
#
#   * **Angles** — 9 rays fanned across `±0.5 rad ≈ ±28.6°`.  The wedge
#     half-width is `36°`; `±28.6°` keeps a `~7°` margin inside the
#     Stokes line so the fan stays well within the pole region and the
#     walk does not stray into the masked boundary strip.
#   * **Radii** — `dr = 1` spacing from `|x| = 2` out to
#     `|x| = SURF_R_MAX`.  At the headline-figure `R = 4`, that is 3
#     radial shells (`2.0, 3.0, 4.0`); at the historical `R = 20`, 19
#     shells × 9 angles = 171 targets, the A2 §5.2 "densest-completing"
#     fan.  Denser fans (13-angle / `dr = 0.5`) BLOCK (A2 §3.3): the
#     extra inter-target bridging walks each get a fresh chance to step
#     onto a pole.
#
# **NOT in the live kernel path.**  `surf_wedge_targets()` (the
# constructor that uses these constants) is the *legacy* sparse polar
# fan; the live `surf_wedge_fill` calls `surf_wedge_dense_targets()`
# (the ADR-0026 D2 dense Cartesian lattice).  These constants are kept
# in lockstep with `SURF_R_MAX` for correctness so a future change of
# generator path stays consistent.
#
# Targets are ordered radius-major so the walk completes inner shells
# before reaching out — the threading order the B2 walk extends from.
const SURF_TARGET_RADII  = Tuple(2.0:1.0:4.0)              # radial shells
const SURF_TARGET_ANGLES = Tuple(range(-0.5, 0.5; length = 9))  # rad

# Pole-extraction filter (`extract_poles_shared_q`).  `min_support = 2`
# is VC-6: a cluster is a physical pole only when ≥ 2 *distinct* nodes
# independently root it — keeps Froissart doublets (spurious
# single-node roots) out of the field.
#
# `SURF_RADIUS_T = 5.0` is **dimensionless** — the ADR-0026 Amendment 6
# §S7 scale-covariance fix.  The filter keeps a root iff its z-plane
# distance from the node satisfies `h_node·|t*| ≤ SURF_RADIUS_T·h_max`,
# i.e. "within `SURF_RADIUS_T` steps at the walk's step ceiling
# `h_max = SURF_PN_H`".  Because `h_max` is independent of the render
# disc radius `SURF_R_MAX`, this constant does **not** change when we
# shrink the wedge from `R = 20` to `R = 4` (bead `padetaylor-apn`):
# the absolute z-radius `5.0·SURF_PN_H = 0.5` is the same in both
# configurations, exactly as the scale-covariance memory mandates.
const SURF_RADIUS_T     = 5.0
const SURF_CLUSTER_ATOL = 0.2
const SURF_MIN_SUPPORT  = 2

# --- The stitch ----------------------------------------------------------
# Half-width (degrees) of the NaN-masked strip on each Stokes line.
#
# C4 (ADR-0025 Amendment 12, bead `padetaylor-0ln.37.18`) — the Phase-E1
# Stokes-strip audition narrows the mask `3° → 1°` in lockstep with
# `SURF_SECTOR_MARGIN_DEG` (`4° → 1°`).  The mask must cover exactly the
# band NEITHER region honestly carries: the sector solver reaches down to
# `arg x = 36° + margin`, and the wedge B1-gated underlay is sparse below
# `|arg| ≈ 34°` (the wedge fan tops out at `±28.6°`), so the genuinely
# uncovered band straddling each Stokes line is `[36°−margin, 36°+margin]`
# — half-width `= SURF_SECTOR_MARGIN_DEG`.  Setting the mask half-width
# equal to the margin is the tightest honest choice: every in-sector cell
# the fan does NOT reach (`arg x ∈ (36°, 36°+margin)`) is caught by the
# mask, and every cell the fan DOES reach (`arg x > 36°+margin`) is
# voted, not masked.  The headline figure's grey Stokes strip shrinks
# from `6°` wide (`[33°,39°]`) to `2°` (`[35°,37°]`) — a 3× reduction —
# with every newly-filled cell backed by the audition's measured
# spread / oracle-error evidence (Amendment 12).  The residual `±1°`
# mask is retained with a forcing condition: the sector fan cannot
# straddle the Stokes line itself (the fan is a *sector* solver, the
# wedge a disjoint region), and the wedge underlay is sparse there; a
# uniform connection formula across the Stokes line would close the
# last `2°` (the deferred-bead forcing condition, Amendment 12).
const SURF_STITCH_MASK_DEG = 1.0

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
# Region 1 — the pole-free-sector triple-method pipeline
# ======================================================================
# `surf_companion`, `surf_ray_bvp`, `surf_ray_fan`, `surf_ray_eval`,
# `surf_sector_rectangle`, `surf_ray_arc_eval`, `surf_laplace_voters`,
# `surf_laplace_eval_one`, `surf_vote` — the three ADR-0024 voters and
# the median vote.  The Rule-6 file-size split (mirroring the VC-4/5/8/10
# helpers); `include`d here, after the figure-defining constants, so the
# sector code sees `SURF_R_MIN`, `SURF_R_INNER_BC`, `SURF_BVP_N`, ….
# ADR-0025 Amendment 10 (C2 BVP-sourced arc data + C3 voter-1
# reconstruction) is implemented in this sibling.
if !isdefined(@__MODULE__, :surf_ray_fan)
    include(joinpath(@__DIR__, "_kkg_pi2_sector.jl"))
end

# ======================================================================
# Region 2 — the pole-rich wedge (V8b path-network + F2 Stage-2 fill)
# ======================================================================

"""
    surf_wedge_targets() -> Vector{ComplexF64}

The B3 extended threading fan for the wedge walk (ADR-0025 Amendment 2).
The Cartesian product `SURF_TARGET_RADII × SURF_TARGET_ANGLES` — 9 rays
fanned across `±0.5 rad` (`±28.6°`, a `~7°` margin inside the `36°`
Stokes line), `dr = 1` radial spacing from `|x| = 2` out to
`|x| = SURF_R_MAX`.  At the headline-figure `R = 4` that is `9 × 3 =
27` targets; at the historical `R = 20`, `9 × 19 = 171`.

**Legacy / not in the live kernel path.**  The live `surf_wedge_fill`
calls `surf_wedge_dense_targets()` (the ADR-0026 D2 dense Cartesian
lattice).  This sparse polar fan is kept available for legacy callers
and parity tests; the docstring rationale below predates the dense
lattice and is retained for context.

Unlike V8b's 20-target `|x| ≤ 8` fan (copied from the pole-scatter
figure), this fan threads the B2 `:max_q_root` adaptive walk through the
*whole* wedge so the walk's node filament passes near every pole out to
`|x| = SURF_R_MAX`.  Targets are emitted **radius-major** (inner shells
first): the walk completes the inner field before reaching out, the
threading order the adaptive walk extends along.
"""
surf_wedge_targets() =
    ComplexF64[r * cis(θ) for r in SURF_TARGET_RADII
                          for θ in SURF_TARGET_ANGLES]

# --- The dense disc-spaced wedge target lattice (ADR-0026 D2) -------------
# Default Cartesian lattice spacing for `surf_wedge_dense_targets`.  The
# B1 honest disc radius (ADR-0025 Amendment 1, the true-radius Stage-2
# gate) is median ≈ 0.06, p90 ≈ 0.16 — so a spacing of ≈ 2× the p90
# disc radius keeps adjacent filaments' validity discs overlapping.
# This is a *starting* value, not a hard-coded fixed one: `s` is a kwarg
# of the generator and ADR-0026 D2 mandates rendering, measuring the
# achieved coverage, and tightening until it saturates.
const SURF_DENSE_SPACING = 0.25
# Inner radius of the dense lattice — the rendered window starts at
# `|x| = 2` (`SURF_R_MIN`); the lattice mirrors that floor.
const SURF_DENSE_R_INNER = 2.0
# Angular half-extent of the dense lattice, in degrees.  The wedge
# half-width is `SURF_WEDGE_HALF_DEG = 36°`; the `±1°` Stokes mask
# (`SURF_STITCH_MASK_DEG`) blanks `35°–37°`, so a `34°` lattice extent
# fills right up to the honest edge of the wedge underlay.  The legacy
# `surf_wedge_targets` fan stopped at `±28.6°`, leaving the genuine
# `28.6°–34°` band untargeted; the dense lattice closes it.
const SURF_DENSE_HALF_DEG = 34.0

"""
    surf_wedge_dense_targets(; s = SURF_DENSE_SPACING,
                               r_inner = SURF_DENSE_R_INNER,
                               half_deg = SURF_DENSE_HALF_DEG,
                               r_outer = SURF_R_MAX)
        -> Vector{ComplexF64}

A **dense, disc-spaced** wedge target set for the resilient Stage-1 walk
(ADR-0026 D2).

## Why a uniform Cartesian lattice, not a polar fan

The legacy `surf_wedge_targets` is a sparse polar fan — 19 radial
shells × 9 angular rays = 171 points.  Two structural problems make it
unable to *tile* the pole-rich wedge:

  1. **It is sparse.**  At the historical `R = 20`, 171
     filament-endpoints over the wedge area (`≈ 235` square units out
     to `|x| = 20`) leave wide gaps between adjacent path-network
     filaments — the headline figure's wedge was `~95 %` blank
     (worklog 059).
  2. **A polar fan has non-uniform areal density.**  A constant-`dθ`
     fan crowds targets near the origin and thins them outward: the
     ring at `|x| = 2` carries the same 9 rays as the ring at
     `|x| = SURF_R_MAX`, so the outer arc length per target is
     `SURF_R_MAX/2 ×` the inner.  The path-network's validity discs
     are roughly `|x|`-independent in size (the B1 honest radius is
     `median ≈ 0.06`), so a fan *over-tiles* the inner wedge and
     *under-tiles* the outer.

FW 2011's analogue tiles its domain with a **uniform `40×40 = 1600`
Cartesian grid** of Stage-1 nodes feeding a `161×161` fine grid
(`references/markdown/FW2011_painleve_methodology_JCP230/`
`FW2011_painleve_methodology_JCP230.md:153-164`) — a *uniform areal
density* so every fine-grid point is within one coarse step of a node.
This generator is the wedge-shaped analogue: a uniform Cartesian
lattice of spacing `s`, clipped to the wedge sector.

## The disc-tiling rationale

The B1 true-radius gate (ADR-0025 Amendment 1) measures each visited
node's *honest* validity disc — median radius `≈ 0.06`, p90 `≈ 0.16`.
For the Stage-2 fill to cover the wedge without gaps, adjacent
filaments' discs must overlap, so the target spacing `s` must be of
order the disc *diameter*.  The default `s = SURF_DENSE_SPACING = 0.25`
is `≈ 2×` the p90 disc radius — the ADR-0026 D2 starting point; the
coverage probe (`figures/probe_wedge_coverage.jl`) sweeps `s` and
measures the achieved coverage fraction so `s` is a *measured*
parameter, never assumed.  `s` is therefore a kwarg, not a constant.

## Clipping and ordering

A lattice point `x = (i·s, j·s)` is kept iff it lies in the wedge
annulus:

  * `r_inner ≤ |x| ≤ r_outer` — inside the rendered window
    (`r_inner ≈ SURF_R_MIN = 2`, `r_outer = SURF_R_MAX = 4`);
  * `|arg x| ≤ half_deg` — inside the wedge sector.  The default
    `half_deg = SURF_DENSE_HALF_DEG = 34°` fills up to `2°` shy of the
    `36°` Stokes line; the residual `35°–37°` band is the `±1°`
    `SURF_STITCH_MASK_DEG` strip, masked downstream.  The legacy fan
    stopped at `±28.6°`, leaving `28.6°–34°` a genuine untargeted gap;
    `34°` closes it.

The surviving lattice points are emitted **radius-major** — sorted by
`|x|` ascending.  This preserves the legacy fan's deliberate ordering
(FW 2011 §5.5, "complete the pole field before stepping into smooth
regions"): the resilient walk threads the dense inner pole field first
and only then reaches outward, so each per-target bridging walk starts
from a nearby already-visited node rather than chording across the
wedge.

## Kwargs

  - `s::Real`        — Cartesian lattice spacing.  Default
                       `SURF_DENSE_SPACING = 0.25` (`≈ 2×` the B1 p90
                       disc radius).  Tighten toward `0.18` for denser
                       coverage; loosen toward `0.35` for a faster
                       walk.  Must be positive.
  - `r_inner::Real`  — inner radius of the kept annulus.  Default
                       `SURF_DENSE_R_INNER = 2.0` (`= SURF_R_MIN`).
  - `half_deg::Real` — angular half-extent in degrees.  Default
                       `SURF_DENSE_HALF_DEG = 34.0`.
  - `r_outer::Real`  — outer radius.  Default `SURF_R_MAX = 4.0`.

Throws `ArgumentError` (CLAUDE.md Rule 1) for a non-positive `s`, a
non-positive `r_inner`, an `r_outer ≤ r_inner`, or a `half_deg` outside
`(0, 90]` — each a caller mistake that would yield an empty or
ill-defined lattice rather than a silent empty return.
"""
function surf_wedge_dense_targets(; s::Real = SURF_DENSE_SPACING,
                                    r_inner::Real = SURF_DENSE_R_INNER,
                                    half_deg::Real = SURF_DENSE_HALF_DEG,
                                    r_outer::Real = SURF_R_MAX)
    s > 0 || throw(ArgumentError(
        "surf_wedge_dense_targets: lattice spacing s must be positive " *
        "(got $s).  Suggestion: pass s ≈ 0.18–0.35 — of order the B1 " *
        "honest disc diameter (ADR-0025: median ≈ 0.06, p90 ≈ 0.16)."))
    r_inner > 0 || throw(ArgumentError(
        "surf_wedge_dense_targets: r_inner must be positive (got " *
        "$r_inner).  Suggestion: pass r_inner ≈ SURF_R_MIN = 2.0."))
    r_outer > r_inner || throw(ArgumentError(
        "surf_wedge_dense_targets: r_outer ($r_outer) must exceed " *
        "r_inner ($r_inner) — the kept annulus would be empty.  " *
        "Suggestion: pass r_outer = SURF_R_MAX = 4.0."))
    (0 < half_deg ≤ 90) || throw(ArgumentError(
        "surf_wedge_dense_targets: half_deg must lie in (0, 90] (got " *
        "$half_deg); it is the wedge angular half-extent in degrees.  " *
        "Suggestion: pass half_deg ≈ 34 — just inside the 36° Stokes " *
        "line."))
    # A uniform Cartesian lattice on `[-r_outer, r_outer]²` at spacing
    # `s`, clipped to the wedge annulus.  `nmax` is the largest integer
    # lattice index whose coordinate stays within `r_outer`.
    nmax = floor(Int, r_outer / s)
    pts = ComplexF64[]
    for i in -nmax:nmax, j in -nmax:nmax
        z = ComplexF64(i * s, j * s)
        r = abs(z)
        (r_inner ≤ r ≤ r_outer) || continue
        abs(rad2deg(angle(z))) ≤ half_deg || continue
        push!(pts, z)
    end
    # Radius-major emission — the FW §5.5 "complete the pole field
    # before stepping into smooth regions" ordering (see the docstring).
    sort!(pts; by = abs)
    return pts
end

"""
    surf_wedge_fill(bvp_sol, grid_pts) -> NamedTuple

Region 2 — the validated pole field + a **dual-provenance** `|u|` surface
(ADR-0025 Amendment 14 / ADR-0026 Amendment 10, `tf9`).

Seed a `VectorPadeTaylorProblem` from the BVP-anchored on-tritronquée
state `bvp_sol(SURF_Z_SEED)`, run the `vector_path_network_solve` wedge
walk over the **dense disc-spaced** target lattice
`surf_wedge_dense_targets()` (ADR-0026 D2) **once**, and evaluate the F2
Stage-2 `fine_grid` fill on `grid_pts` (the Cartesian grid points inside
the wedge) **twice** — both fills read the same Stage-1 walk so the
second one costs only the per-cell Padé evaluation.

The dual fill is the ADR-0026 Amendment 9 closeout decision: the
D1–S8 investigation arc proved the no-extrapolation honesty contract
(ADR-0025) is **stricter than FW 2011 itself** — FW's published pole-
field figures FILL precisely because FW's Stage-2 evaluates each Padé
over its full step with no validity gate (FW-md:395-397).  Rendering the
wedge FW-faithful (filled) without lying about provenance demands two
fills, not one.

Three provenance-defining choices:

  * the **certified** Stage-2 fill runs with **`extrapolate = false`** —
    the B1 true-radius gate (ADR-0025 Amendment 1) applies, so a grid
    point outside the nearest node's verified disc gets `NaN`.  **No
    Padé is ever evaluated outside its disc.**  `tol = SURF_PN_TOL`
    selects the gate's `s(tol)` safety factor.  This is the honest core
    of the figure — the set of grid points the B1 gate certifies — and
    it is what `covered` flags.  The certified fraction (~29 % at order
    48 — the measured knee, ADR-0026 Amendment 9) is the figure's
    *senior-grade* deliverable: a region where no Padé is evaluated
    beyond its truncation-bounded disc;

  * the **display** Stage-2 fill runs with **`extrapolate = true`** —
    the legacy ADR-0015 fill mode, identical to what FW 2011's published
    figures use.  Every Padé is evaluated over its full step, with no
    validity gate, so every wedge cell is finite.  This is the
    FW-faithful filled-wedge rendering — what KKG 2015 Figs 7.4/7.5
    visualise.  ADR-0025's strict no-extrapolation rule **remains the
    package default** and applies to every non-headline figure; this is
    a per-figure override, justified by FW precedent and honest *only
    when* the figure marks provenance (the figure script must visually
    distinguish certified from extrapolated cells — reduced alpha,
    hatch, or desaturation — so the reader can tell which is which);

  * **the walk itself is shared.**  One Stage-1 walk + two Stage-2
    re-evaluations, not two walks: the resilient adaptive `:max_q_root`
    walk is the expensive part, and both fills read the same visited
    tree.  The second fill calls
    `PadeTaylor.VectorPathNetworkStage2._stage2_fill` directly with
    `extrapolate = true` (and the walk's cached `visited_jets`,
    `visited_numerators`, `visited_denominator`).  This is the additive
    helper ADR-0026 Amendment 10 introduces; it adds no `src/` API
    surface — the figure helper reaches into the existing internal
    function rather than the package re-exporting it, since the dual
    fill is a *figure-specific* pattern (the legitimate scope of an
    `extrapolate = true` override).

Three further honesty-defining choices (unchanged from the pre-`tf9`
pipeline):

  * the walk runs the **S2 scale-derived adaptive step law**
    (`adaptive = true`) with the **S3 skip-if-covered** target gate
    (`skip_covered = true`) — ADR-0026 Amendments 3–4.  This supersedes
    the interim R1 fixed-`h` config (Amendment 2): R1 was the right
    *interim* fix for the now-cured collapse, but a fixed `h` is itself
    a scale heresy — `h = 0.1` is ~2× the honest B1 disc, so the walk
    steps past where each Padé is honestly valid and leaves gaps along
    every filament (Amendment 3, bottleneck 1).  Amendment 4's S1
    diagnostic root-caused the *old* `_adaptive_h` collapse as a genuine
    **geometric sink** — its pole cap `POLE_SAFETY·h_prev·min|t*|` was
    purely multiplicative, so `h` ratcheted to zero in a sustained pole
    field — and replaced it with the absolute, non-ratcheting law
    `h = clamp(SAFETY·D_local, h_min, h_max)`, where `D_local =
    h_prev·min|t*|` is the `h`-independent absolute pole-free disc
    radius (S1 EXP D).  `h` now jumps straight to a fraction of the
    *true* local disc each step — a node leaving a dense pocket recovers
    to `h_max` in one step, there is no sink — and tracks the
    `|x|⁻¹ᐟ⁶` pole-spacing scaling of the P_I⁽²⁾ field for free.
    `adaptive = true` is therefore now *correct*, and FW-faithful in
    spirit: the step is sized to the local pole field, not hand-tuned.
    The **S3** `skip_covered = true` gate then skips any target already
    inside a visited node's `COVER_FRAC·visited_h` covering disc, so the
    dense lattice does not make the walk re-traverse covered ground —
    that re-walking chorded across the wedge and clipped poles into the
    non-monotonic `:all_candidates_failed` regression (Amendment 3,
    bottleneck 3);

  * the walk runs with **`on_target_failure = :skip`** — the resilient
    Stage-1 walk (ADR-0026 D1).  A dense target lattice has thousands
    of targets; under the default `:throw` policy a single unreachable
    target would abort the whole walk and cost every other target's
    coverage.  With `:skip` a failed per-target walk is recorded as a
    `VectorWalkFailure` in `walk.failed_targets` and the walk continues.
    This is fail-*loud-by-accounting*, never a silent swallow (CLAUDE.md
    Rule 1): the count of skipped targets is surfaced in `message`, and
    `walk.failed_targets` carries one first-class record per skip;

  * the **certified Stage-2 fill** runs with **`extrapolate = false`** —
    the B1 true-radius gate (ADR-0025 Amendment 1) applies; the
    **display Stage-2 fill** runs with **`extrapolate = true`** —
    FW-2011-style filled (see the dual-fill rationale above).
    `tol = SURF_PN_TOL` selects the certified fill's gate `s(tol)`
    safety factor.  The two fills share the Stage-1 walk.

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
  * **VC-5** (`vc5_pair`, ADR-0025 Amendment 6) — match the
    VC-4-surviving poles into conjugate pairs (`V₀(x̄)=conj V₀(x)`) by a
    globally-optimal (maximum-cardinality) bipartite matching of the
    conjugate-admissibility graph; the pairing residual is an FW-style
    accuracy estimate, a flagged off-axis unpaired pole is suspect.
    VC-5 is a diagnostic — it *reports*, it does not prune, and it
    never *constructs* a pole from its mirror, so the residual stays a
    genuine accuracy cross-check (FFW 2017 `:120-124`).

Returns `(walk, poles, u, u_extrap, covered, vc4, vc5, message)`:
  - `walk`     : the `VectorPathNetworkSolution`;
  - `poles`    : the **VC-4-validated** pole field — the figure's final
                 field (the FW Fig 4.7 pole-location scatter);
  - `u`        : `Vector{ComplexF64}` of `V_0` at `grid_pts` — the
                 **B1-certified** fill (`extrapolate = false`): finite
                 only where the verified node disc covers, `NaN +
                 NaN·im` everywhere else.  The senior-grade honest core
                 of the figure;
  - `u_extrap` : `Vector{ComplexF64}` of `V_0` at `grid_pts` — the
                 **FW-2011-style filled** display surface
                 (`extrapolate = true`): every Padé evaluated over its
                 full step with no validity gate, every cell finite.
                 Paired with `u` so the figure script can composite the
                 two and visually mark provenance (full opacity on
                 certified, reduced alpha on extrapolated — ADR-0026
                 Amendment 10);
  - `covered`  : `Vector{Bool}`, `true` iff that grid point is B1-gate
                 certified (`u[k]` finite) — the honest-coverage mask.
                 Note that `covered ⟺ isfinite(u)` (the certified
                 contract is unchanged); `u_extrap` is finite
                 *everywhere*, so `covered` is the **provenance** mask
                 the figure paints alpha by;
  - `vc4`      : the `vc4_validate` diagnostics — per-candidate `A`,
                 `B`, `dA`, `pass`, `family`, the pruned count and
                 per-pruned-pole failure reason;
  - `vc5`      : the `vc5_pair` diagnostics — conjugate pairs, the
                 pairing residuals, median/max, unpaired + flagged;
  - `message`  : a status note (node count, dense-target count and the
                 number skipped by the resilient walk, candidate /
                 validated pole counts, `|x|` frontier, honest-coverage
                 fraction, VC-4 prune count, VC-5 residual).
"""
function surf_wedge_fill(bvp_sol::VectorBVPSolution,
                         grid_pts::AbstractVector{ComplexF64};
                         rng = nothing)
    f      = painleve_hierarchy(:I, 2; t = SURF_T)
    y_seed = ComplexF64.(bvp_sol(ComplexF64(SURF_Z_SEED)))
    # `zspan` carries (start, nominal endpoint).  Only `zspan[1]` is
    # load-bearing — the path-network walk seeds at `zspan[1]` and walks
    # to user-supplied `targets`, not to `zspan[2]`.  The endpoint
    # `ComplexF64(SURF_R_MAX + 0im)` records the outer reach of the
    # render disc; it tracks `SURF_R_MAX` in lockstep (bead
    # `padetaylor-apn`).
    prob   = VectorPadeTaylorProblem(f, y_seed,
                                     (ComplexF64(SURF_Z_SEED),
                                      ComplexF64(SURF_R_MAX + 0.0im));
                                     order = SURF_PN_ORDER)
    targets = surf_wedge_dense_targets()
    # S2 scale-derived adaptive `:max_q_root` walk — `adaptive = true`
    # (ADR-0026 Amendments 3–4).  This supersedes the interim R1
    # fixed-`h` config: Amendment 4's `_adaptive_h` is now the absolute,
    # non-ratcheting law `h = clamp(SAFETY·D_local, h_min, h_max)`, with
    # no multiplicative `h_prev` dependence — the geometric sink that
    # forced R1 is gone by construction, so `adaptive = true` is correct
    # and sizes each step to the *local* honest B1 disc rather than the
    # scale-blind fixed `h = 0.1`.  `skip_covered = true` is the S3 gate:
    # a target already inside a visited node's `COVER_FRAC·visited_h`
    # covering disc is skipped, so the dense lattice does not re-walk
    # covered ground (the chording that clipped poles — Amendment 3,
    # bottleneck 3).  `h = SURF_PN_H` now seeds `h_max` (the step
    # *ceiling*, not a fixed step).  This first solve produces the
    # **certified** Stage-2 fill (`extrapolate = false`, B1 true-radius
    # gate); the **display** fill is then re-evaluated from the same
    # walk with `extrapolate = true` (the dual-fill design — see the
    # docstring above).  `on_target_failure = :skip` (ADR-0026 D1) makes
    # the dense walk resilient: a per-target walk failure is recorded in
    # `walk.failed_targets` rather than aborting the whole solve — a
    # dense lattice has thousands of targets and `:throw` would cost
    # every target's coverage on the first failure.  `rng` controls the
    # Stage-1 target processing order — `nothing` (the figure default)
    # walks the radius-major lattice in order, an `AbstractRNG` shuffles
    # it for the VC-10 double-run accuracy indicator
    # (`surf_vc10_two_run`).
    walk = vector_path_network_solve(prob, targets;
                                     order             = SURF_PN_ORDER,
                                     h                 = SURF_PN_H,
                                     adaptive          = true,
                                     skip_covered      = true,
                                     fine_grid         = grid_pts,
                                     extrapolate       = false,
                                     tol               = SURF_PN_TOL,
                                     rng               = rng,
                                     on_target_failure = :skip,
                                     verbose           = true)
    # ----------------------------------------------------------------------
    # Dual-fill provenance evaluation (ADR-0025 Amendment 14 / ADR-0026
    # Amendment 10, `tf9`).  The first `vector_path_network_solve` call
    # ran the Stage-2 fill once with `extrapolate = false` — the
    # B1-certified core stored in `walk.grid_y`.  We now re-evaluate the
    # same fine grid against the same visited tree with
    # `extrapolate = true` to produce the FW-faithful filled display
    # surface.  The walk is the expensive part (every visited node, every
    # shared-Q Padé construction); the second fill is a Horner sweep over
    # the cached `visited_numerators` / `visited_denominator` and costs
    # only the per-cell Padé evaluation — no new walk, no rebuild.
    #
    # Reaching into `PadeTaylor.VectorPathNetworkStage2._stage2_fill`
    # (and the sibling private `_nearest_visited`) is the legitimate
    # figure-side override of the no-extrapolation rule: ADR-0025
    # Amendment 14 keeps `extrapolate = false` the package default for
    # every non-headline figure and the public solver API, and ADR-0026
    # Amendment 10 explicitly scopes `extrapolate = true` to the
    # headline figure's display overlay — a per-figure provenance-marked
    # exception, justified by FW 2011 precedent (FW's published pole-
    # field figures evaluate each Stage-2 Padé over its full step with
    # no validity gate, FW-md:395-397).
    grid_y_extrap = PadeTaylor.VectorPathNetworkStage2._stage2_fill(
        walk.grid_z,
        walk.visited_z, walk.visited_h,
        walk.visited_numerators, walk.visited_denominator,
        walk.visited_jets,
        true,                   # extrapolate = true — the FW-style overlay
        Float64(SURF_PN_TOL),
        prob.f, walk.visited_y, Int(SURF_PN_ORDER),
        ComplexF64,
        PadeTaylor.VectorPathNetwork._nearest_visited)
    # S4 (ADR-0026 Amendment 3) — `cluster_atol = nothing` puts
    # `extract_poles_shared_q` in scale-derived clustering mode: the
    # cluster tolerance is a fraction of the *local* pole spacing rather
    # than the absolute `SURF_CLUSTER_ATOL = 0.2`, which wrongly
    # merged/split poles and corrupted the pole count.
    #
    # The `[wedge ...]` phase markers around the post-walk validation are
    # the same "verbose mode + eager flush" discipline as the existing
    # `[kernel ...]` markers in `kkg_pi2_surface()` (memory
    # `long-julia-agents-verbose-flush`): each validation phase emits a
    # stdout line and `flush(stdout)` so the log file streams progress
    # live under background-redirection.  `t_validate` is local to the
    # validation block — the surface kernel's `tphase` is not in scope
    # inside `surf_wedge_fill`.  Bead `padetaylor-dqu`.
    t_validate = time()
    @printf("  [wedge %5.1fs] extract_poles_shared_q (radius_t=%.1f, min_support=%d)...\n",
            time() - t_validate, SURF_RADIUS_T, SURF_MIN_SUPPORT); flush(stdout)
    candidates = extract_poles_shared_q(walk;
                                        radius_t     = SURF_RADIUS_T,
                                        cluster_atol = nothing,
                                        min_support  = SURF_MIN_SUPPORT)
    @printf("  [wedge %5.1fs] extract_poles_shared_q → %d candidates\n",
            time() - t_validate, length(candidates)); flush(stdout)
    # Phase-D per-pole validation: VC-4 prunes the spurious candidates
    # (Froissart doublets / out-of-family / non-zero-residue), VC-5
    # pairs the survivors by conjugate symmetry.  `poles` is the
    # VC-4-validated field — the field the figure renders.
    @printf("  [wedge %5.1fs] vc4_validate (%d candidates)...\n",
            time() - t_validate, length(candidates)); flush(stdout)
    vc4   = vc4_validate(walk, ComplexF64.(candidates))
    @printf("  [wedge %5.1fs] vc4_validate → %d kept (%d pruned)\n",
            time() - t_validate, length(vc4.kept), vc4.n_pruned); flush(stdout)
    poles = vc4.kept
    @printf("  [wedge %5.1fs] vc5_pair (%d poles)...\n",
            time() - t_validate, length(poles)); flush(stdout)
    vc5   = vc5_pair(poles)
    @printf("  [wedge %5.1fs] vc5_pair done\n",
            time() - t_validate); flush(stdout)
    # The B1-certified Stage-2 grid_y is ordered as grid_pts; `u` is its
    # first component (the `V_0` value at each wedge grid point).  A
    # `NaN` slot is an honest B1 gate gap (a grid point outside every
    # node's verified disc) — never an extrapolated lie.  This is the
    # certified core: `covered ⟺ isfinite(u)` is the figure's honesty
    # contract (ADR-0025 Amendment 1/14).
    u = ComplexF64[any(isnan, gy) ? ComplexF64(NaN, NaN) : gy[1]
                   for gy in walk.grid_y]
    covered = Bool[!isnan(real(uz)) for uz in u]
    # The FW-style display fill, evaluated over the same grid against the
    # same walk with no validity gate — every cell finite, the wedge
    # filled.  `u_extrap[k]` is the FW-2011-style rendering of `V_0` at
    # `grid_pts[k]`; the figure script paints `covered[k]` cells at full
    # opacity (using `u[k]`) and `!covered[k]` cells at reduced opacity
    # (using `u_extrap[k]`), so the reader can read provenance off the
    # alpha channel.
    u_extrap = ComplexF64[any(isnan, gy) ? ComplexF64(NaN, NaN) : gy[1]
                          for gy in grid_y_extrap]
    # The |x| frontier the walk's node filament reached.
    frontier = isempty(walk.visited_z) ? 0.0 :
               maximum(abs, walk.visited_z)
    cov_frac = isempty(u) ? 0.0 : count(covered) / length(u)
    n_extrap_finite = count(uz -> !isnan(real(uz)), u_extrap)
    extrap_frac = isempty(u_extrap) ? 0.0 :
                  n_extrap_finite / length(u_extrap)
    # The resilient-walk accounting (ADR-0026 D1): with
    # `on_target_failure = :skip` the walk records one VectorWalkFailure
    # per target it could not reach.  A non-empty `failed_targets` is the
    # fail-loud-by-accounting signal — surface its count in `message`.
    msg = "wedge: $(length(walk.visited_z)) visited nodes, " *
          "$(length(targets)) dense targets " *
          "($(length(walk.failed_targets)) skipped), " *
          "$(length(candidates)) candidate poles → " *
          "$(length(poles)) VC-4-validated " *
          "($(vc4.n_pruned) pruned), " *
          "|x|-frontier $(round(frontier; digits=1)), " *
          "$(count(covered))/$(length(u)) grid points B1-certified " *
          "($(round(100 * cov_frac; digits=1))%) — full opacity; " *
          "FW-style overlay $(n_extrap_finite)/$(length(u_extrap)) " *
          "filled ($(round(100 * extrap_frac; digits=1))%) — reduced " *
          "alpha; " *
          "VC-5 conjugate-pair residual median " *
          "$(round(vc5.median_resid; digits=4))"
    return (walk = walk, poles = poles, u = u, u_extrap = u_extrap,
            covered = covered, vc4 = vc4, vc5 = vc5, message = msg)
end

# VC-10 — the FW/FFW two-run pole-field accuracy indicator (ADR-0025
# Validation Criteria Menu VC-10; bead `padetaylor-0ln.37.17`,
# Amendment 9).  Sibling helper — the Rule-6 file-size split, mirroring
# `_kkg_pi2_vc45.jl` (VC-4/VC-5) and `_kkg_pi2_vc8.jl` (VC-8).  Pulled
# in here, after `surf_wedge_fill` (which it runs twice) and after
# `_kkg_pi2_vc45.jl` (whose `_vc5_augment!` matcher core it reuses).
if !isdefined(@__MODULE__, :surf_vc10_two_run)
    include(joinpath(@__DIR__, "_kkg_pi2_vc10.jl"))
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
the Cartesian square `[-SURF_XY_LIM, SURF_XY_LIM]²` (the headline
figure renders `SURF_XY_LIM = 4.0`; bead `padetaylor-apn`).

Construction (see the file docstring):
  1. **Region 1 — the pole-free sector** is computed three ways
     (ray-fan BVP / in-house 2D-Chebyshev Laplace / Gridap FEM Laplace)
     and **majority-voted per grid point** (median); an agreement /
     spread map records the max pairwise disagreement.
  2. **Region 2 — the pole-rich wedge** (ADR-0025 Amendment 14 /
     ADR-0026 Amendment 10) is the validated pole *field* + a **dual-
     fill** Re/Im surface: the certified fill (`extrapolate = false`,
     B1-honest, ~29 % of the wedge) in `Re_u` / `Im_u` and the FW-
     2011-style filled overlay (`extrapolate = true`, every cell finite,
     identical to FW's published rendering — FW-md:395-397) in
     `Re_u_extrap` / `Im_u_extrap`.  The figure script composites the
     two with `wedge_covered` to mark provenance via alpha.
  3. The two regions are stitched onto one grid; the `±36°` Stokes-line
     strips are `NaN`-masked.

Returns a `NamedTuple`:
  - `xs`, `ys`            : the `SURF_GRID_N`-point axes of
                            `[-SURF_XY_LIM, SURF_XY_LIM]`;
  - `Re_u`, `Im_u`        : `SURF_GRID_N × SURF_GRID_N` matrices, the
                            **B1-certified** rendering —
                            `Re_u[i,j] = Re V_0(xs[i] + im·ys[j])`
                            (`NaN` outside the disc `|x| ≤ SURF_XY_LIM`,
                            in the masked strips, in the wedge wherever
                            the B1
                            gate found no honest datum, and where no
                            sector method has a datum).  In the sector
                            this is the triple-method vote; in the
                            wedge the no-extrapolation honesty contract
                            applies — every finite cell is B1-gate-
                            certified, every uncovered cell is `NaN`;
  - `Re_u_extrap`,
    `Im_u_extrap`         : `SURF_GRID_N × SURF_GRID_N` matrices, the
                            **FW-2011-style filled overlay** — in the
                            wedge every cell finite (every Padé
                            evaluated over its full step with no
                            validity gate, the rendering FW 2011 itself
                            uses for its published pole-field figures);
                            in the sector identical to `Re_u`/`Im_u`
                            (the sector is voted, not Padé-evaluated —
                            provenance marking applies only to the
                            wedge).  Paired with `Re_u`/`Im_u` so the
                            figure script can composite the two and
                            paint extrapolated cells at reduced alpha
                            (ADR-0025 Amendment 14 / ADR-0026
                            Amendment 10);
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
    # The FW-style extrapolated overlay matrices (ADR-0025 Amendment 14 /
    # ADR-0026 Amendment 10).  In the wedge they carry the FW-faithful
    # filled rendering — every Padé evaluated over its full step with no
    # validity gate, the rendering FW 2011 itself uses.  In the sector
    # and outside the wedge they mirror Re_u/Im_u (the sector is voted,
    # not Padé-evaluated; the dual-fill provenance applies only to the
    # wedge underlay).  The figure script composites them with `Re_u`/
    # `Im_u` and the `wedge_covered` mask to mark provenance via alpha.
    Re_u_extrap = fill(NaN, n, n)
    Im_u_extrap = fill(NaN, n, n)
    spread = fill(NaN, n, n)
    wedge_covered = fill(false, n, n)
    sector_method = Dict{Tuple{Int,Int},
                         NTuple{6,Float64}}()

    # The phase markers below use the "verbose mode + eager flush"
    # discipline (memory `long-julia-agents-verbose-flush`): each
    # major boundary inside this kernel emits a stdout line and
    # `flush(stdout)` so the log file streams progress live under
    # background-redirection.  Without these, the order-48 render's
    # ~hour-long single call would log nothing until completion.
    tphase = time()
    @printf("  [kernel %5.1fs] Region 1.1 — surf_ray_fan...\n",
            time() - tphase); flush(stdout)

    # ---- Region 1: the three sector voters --------------------------------
    fan  = surf_ray_fan()
    @printf("  [kernel %5.1fs] Region 1.2 — surf_laplace_voters (2D Cheb + Gridap FEM)...\n",
            time() - tphase); flush(stdout)
    lap  = surf_laplace_voters(fan)

    # ---- Region 2: the wedge BVP anchor + path-network fill ---------------
    @printf("  [kernel %5.1fs] Region 2.1 — surf_anchor_bvp...\n",
            time() - tphase); flush(stdout)
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
    @printf("  [kernel %5.1fs] Region 2.2 — surf_wedge_fill (%d wedge grid points, order=%d, the big one)...\n",
            time() - tphase, length(wedge_pts), SURF_PN_ORDER); flush(stdout)
    wedge = isempty(wedge_pts) ?
        (walk = nothing, poles = ComplexF64[], u = ComplexF64[],
         u_extrap = ComplexF64[], covered = Bool[],
         vc4 = nothing, vc5 = nothing,
         message = "wedge: no grid points") :
        surf_wedge_fill(anchor, wedge_pts)
    @printf("  [kernel %5.1fs] surf_wedge_fill done — %s\n",
            time() - tphase,
            wedge === nothing ? "nothing" :
            (wedge.walk === nothing ? "empty" : "$(length(wedge.walk.visited_z)) nodes, $(length(wedge.poles)) poles")); flush(stdout)

    # ---- assemble the grid ------------------------------------------------
    @printf("  [kernel %5.1fs] Region 1.3 — per-cell sector voting loop (%d×%d grid)...\n",
            time() - tphase, n, n); flush(stdout)
    for j in 1:n, i in 1:n
        z = ComplexF64(xs[i], ys[j])
        # Outside the |x| ≤ SURF_XY_LIM disc: leave NaN (the figure is the disc).
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
            # The sector is voted, not Padé-evaluated — the extrapolated
            # overlay simply mirrors the certified vote here.  Provenance
            # marking (the alpha override) applies only to the wedge
            # underlay; the sector is always full-opacity.
            Re_u_extrap[i, j] = vre
            Im_u_extrap[i, j] = vim
            spread[i, j] = max(sre, sim)
            sector_method[(i, j)] =
                (re1, im1, re2, im2, re3, im3)
        end
    end

    # ---- write Region 2's dual-fill underlay into the grid ----------------
    # ADR-0025 Amendment 14 / ADR-0026 Amendment 10 — the headline figure
    # renders the wedge with TWO provenance-distinguished fills:
    #
    #   * the **certified** fill (`extrapolate = false`, B1 true-radius
    #     gate) — written into `Re_u` / `Im_u`.  A cell finite in these
    #     matrices is B1-gate-honest: no Padé evaluated outside its
    #     truncation-bounded disc.  The figure script paints these cells
    #     at full opacity.  `wedge_covered[i,j] == true` ⟺ both matrices
    #     finite here ⟺ certified.
    #
    #   * the **FW-style extrapolated** fill (`extrapolate = true`,
    #     identical to FW 2011's published pole-field figures — FW-
    #     md:395-397) — written into `Re_u_extrap` / `Im_u_extrap`.
    #     Every wedge cell is finite (no validity gate); cells the
    #     certified fill leaves `NaN` get a value here too.  The figure
    #     script paints `!wedge_covered[i,j]` cells from these matrices
    #     at reduced opacity — the alpha drop marks them as evaluated
    #     past the verified disc.
    #
    # This dual provenance is FW-faithful (the wedge fills as FW renders
    # it) AND honest (the certified core is visually distinguishable
    # from the extrapolated overlay).  The certified rule remains the
    # package default and applies to every non-headline figure
    # (ADR-0025); this is a per-figure override justified by FW
    # precedent.
    if wedge.walk !== nothing
        for (k, (i, j)) in enumerate(wedge_idx)
            uz       = wedge.u[k]
            uz_extra = wedge.u_extrap[k]
            # FW-style overlay first — every wedge cell where the Padé
            # evaluation produced a finite number.  This is the always-
            # filled, FW-faithful surface; reduced alpha at render time
            # marks it as extrapolated.
            if !isnan(real(uz_extra))
                Re_u_extrap[i, j] = real(uz_extra)
                Im_u_extrap[i, j] = imag(uz_extra)
            end
            # Certified core: only B1-gated covered cells.  An uncovered
            # cell stays `NaN` in `Re_u` / `Im_u` — an honest gap; the
            # display surface above carries the extrapolated rendering.
            if isnan(real(uz))
                continue
            end
            Re_u[i, j] = real(uz)
            Im_u[i, j] = imag(uz)
            wedge_covered[i, j] = true
        end
    end

    return (xs = xs, ys = ys, Re_u = Re_u, Im_u = Im_u,
            Re_u_extrap = Re_u_extrap, Im_u_extrap = Im_u_extrap,
            spread = spread, poles = wedge.poles,
            vc4 = wedge.vc4, vc5 = wedge.vc5,
            wedge_walk = wedge.walk,
            sector_method = sector_method,
            wedge_covered = wedge_covered,
            ray_fan = fan, message = wedge.message)
end
