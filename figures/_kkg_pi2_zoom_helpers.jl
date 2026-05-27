# figures/_kkg_pi2_zoom_helpers.jl
#
# THE ZOOM-VARIANT KERNEL for the P_I^(2) tritronquée `V₀(x, 0)` on
# the **inner** complex-plane square `[-2, 2]² ⊂ |x| ≤ 2√2`.  Built as
# the companion to `_kkg_pi2_surface_helpers.jl`: same physics (the same
# `painleve_hierarchy(:I, 2; t = SURF_T)` companion system, the same
# `surf_anchor_bvp` BVP-anchored seed at `z = SURF_Z_SEED = -3 + 0im`),
# but a deliberately simpler integration *topology* — one straight Padé
# walk from the seed to every grid point, no wedge-walk + sector-vote +
# Stokes-strip mask multi-region scaffolding.
#
# ## Why a separate kernel and not just a re-render of the R = 8 cache
#
# The R = 8 figures (`kkg_pi2_abs_heatmap.png`,
# `kkg_pi2_abs_phase_surface.png`) are produced by the full
# `kkg_pi2_surface()` kernel, which stitches a pole-free sector (BVP
# ray-fan + Laplace voters) onto a pole-rich wedge (path-network walk)
# across the `±36°` Stokes lines.  Each region is rendered correct on
# its own; the join is visually rough in a narrow band along the
# wedge/sector boundary at `|x| ≈ 2`–3.  The user reports the
# concrete symptom: "several false poles around `Im x = -2`, `Re x = 2.5`."
# That band is exactly the inner part of the wedge where the FW-style
# extrapolated fill from the wedge meets the median-vote sector field
# at different effective scales.
#
# The zoom region `[-2, 2]²` is a square of side 4; the circumscribed
# circle has radius `2√2 ≈ 2.83`, so the entire zoom is comfortably
# inside the disc `|x| < 2.83`.  More importantly, the smooth pole-free
# sector covers `|arg x| > 36°`; only the narrow wedge slice
# `Re x ∈ [√(4 − Im²), 2]` with `|arg x| < 36°` could intersect the
# multi-region join, and even that slice's tip sits at `(2, 0)`.  The
# bulk of `[-2, 2]²` is squarely in the smooth, pole-free sector — a
# single Padé walk from the BVP-anchored seed at `z = -3 + 0im` (which
# is **outside** the zoom square, `Re = -3 < -2`) integrates *inward*
# through clean smooth-sector physics with no wedge contact at all.
#
# So the zoom kernel is **strictly simpler** than the surface kernel:
# one VectorPadeTaylorProblem, one `vector_path_network_solve` walk
# with the same `(order, h, tol)` knobs as the wedge of the R = 8
# kernel, and a 8-target perimeter seed-set that guarantees Stage-1
# reaches every corner of the zoom square before Stage-2 fills the
# 801² fine lattice by Padé evaluation at the nearest visited node.
#
# ## Resolution choice
#
# 801 × 801 = 641 601 grid points over the 4-unit-wide square gives a
# pixel size of `dx = 4 / 800 = 0.005`, i.e. **8 × finer than the R = 8
# cache's `dx = 0.04`**.  That is the resolution lift the user
# requested.  Stage-2 is `O(N_grid)` Padé evaluations at the
# nearest-visited node, so total runtime is dominated by the Stage-2
# fill (~30–120 s on a modern CPU); the kernel returns in under
# ~3 min.
#
# ## What this kernel intentionally does NOT do
#
#   * no sector ray-fan, no Laplace voters, no median vote — there are
#     no wedge poles to vote around inside `|x| ≤ 2`;
#   * no Stokes-strip mask — the strips at `±36°` start at `|x| = 2`
#     (outside the zoom for `|y| < 2`), so the zoom has no strip to
#     mask;
#   * no provenance gate — Stage-2 runs with `extrapolate = true`,
#     matching the worklog-062 "FW-2011-style filled everywhere"
#     convention the zoom inherits from its parent figures;
#   * no pole-field extraction — the smooth-sector zoom has no
#     poles in `|x| ≤ 2` to extract.
#
# Failure semantics follow Rule 1 (fail loud): `on_target_failure =
# :skip` records corner-target failures in `walk.failed_targets` so
# the render script can ASSERT zero failures (a non-empty
# `failed_targets` means the seed walk hit something unexpected on
# its way to a perimeter corner — investigate, do not silently fill).
#
# References:
#   * `figures/_kkg_pi2_surface_helpers.jl` — the parent kernel; this
#     file re-uses `surf_anchor_bvp`, `SURF_Z_SEED`, `SURF_T` from
#     there.
#   * `docs/worklog/062-complex-function-figures.md` — the R = 8
#     headline figures the zoom complements.
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015 §7 (the V₀ tritronquée); the inner `|x| ≤ 2` region
#     is the smooth-sector core of the tritronquée surface.

# Pull in the parent kernel's helpers — `surf_anchor_bvp`,
# `SURF_Z_SEED`, `SURF_T`, the `painleve_hierarchy(:I, 2)` companion
# system constants — by `include`-ing the surface helpers verbatim.
# This is the Law-2-correct path: re-use the existing, documented
# anchor BVP rather than duplicating its 20-line definition.
include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))

# `Printf` for the verbose-flush phase markers and the message
# formatting.  The parent helpers already `using Printf`, but be
# explicit so the zoom kernel module remains readable as a unit.
using Printf

# ======================================================================
# Zoom-defining constants
# ======================================================================

# The zoom window — a square of side `2·ZOOM_XY_LIM` centred at the
# origin.  `2.0` is the user-requested zoom radius; the entire square
# `[-2, 2]²` sits inside the disc `|x| ≤ 2√2 ≈ 2.83`.
const ZOOM_XY_LIM = 2.0

# 801 × 801 = 641 601 grid points.  `dx = 4 / 800 = 0.005`, i.e. 8 ×
# finer than the R = 8 cache's `dx = 0.04`.  Odd so a node sits on
# `x = 0` and `y = 0` (the negative-real ridge passes through both).
const ZOOM_GRID_N = 801

# Walk knobs — match the R = 8 wedge walk so the zoom is identical in
# *physics* to the surface kernel's wedge solve, only the *topology* is
# simpler (no wedge fan, just one seed → perimeter sweep).
const ZOOM_PN_ORDER = 48
const ZOOM_PN_H     = 0.1
const ZOOM_PN_TOL   = 1.0e-8

# The `P_I^(2)` time parameter.  Match the surface kernel's `SURF_T`
# (= 0.0) so zoom and surface are slices of the same flow.
const ZOOM_T = 0.0

# ======================================================================
# The zoom kernel
# ======================================================================

"""
    kkg_pi2_zoom() -> NamedTuple

Compute `V₀(x, 0)` of the `P_I^(2)` tritronquée on the dense Cartesian
lattice over `[-ZOOM_XY_LIM, ZOOM_XY_LIM]²` via a single
`vector_path_network_solve` walk seeded at `surf_anchor_bvp()` evaluated
at `SURF_Z_SEED = -3 + 0im`.  Stage-1 walks to 8 perimeter targets
(corners + midpoints) to span the zoom square; Stage-2 fills the 801²
fine grid by Padé evaluation at the nearest visited node.

Returns a `NamedTuple`:
  - `xs`, `ys`    — the `ZOOM_GRID_N` axes of `[-ZOOM_XY_LIM, ZOOM_XY_LIM]`;
  - `Re_u`, `Im_u` — `ZOOM_GRID_N × ZOOM_GRID_N` matrices with
                    `Re_u[i, j] = Re V₀(xs[i] + im·ys[j])`.  Indexing
                    matches Makie's heatmap convention;
  - `walk`        — the `VectorPathNetworkSolution` (kept for
                    diagnostics — not used by the rendering scripts);
  - `message`     — status note: kernel runtime, visited-node count,
                    failed-target count.

Throws via Rule 1 (fail loud) when `surf_anchor_bvp` fails or when the
`VectorPadeTaylorProblem` constructor rejects the seed.  Stage-2
failures DO NOT throw — they are recorded as NaN cells in `Re_u`,
`Im_u`; the render script's invariants flag them.
"""
function kkg_pi2_zoom()
    @printf("KKG zoom — V₀(x, 0) over [-%.1f, %.1f]² at %d×%d\n",
            ZOOM_XY_LIM, ZOOM_XY_LIM, ZOOM_GRID_N, ZOOM_GRID_N); flush(stdout)

    t0 = time()
    @printf("  [zoom %5.1fs] surf_anchor_bvp...\n", time() - t0); flush(stdout)
    anchor = surf_anchor_bvp()

    @printf("  [zoom %5.1fs] building VectorPadeTaylorProblem...\n",
            time() - t0); flush(stdout)
    f      = painleve_hierarchy(:I, 2; t = ZOOM_T)
    y_seed = ComplexF64.(anchor(ComplexF64(SURF_Z_SEED)))
    prob   = VectorPadeTaylorProblem(f, y_seed,
                                     (ComplexF64(SURF_Z_SEED),
                                      ComplexF64(ZOOM_XY_LIM + 0.0im));
                                     order = ZOOM_PN_ORDER)

    # The Cartesian grid axes.
    xs = collect(range(-ZOOM_XY_LIM, ZOOM_XY_LIM; length = ZOOM_GRID_N))
    ys = collect(range(-ZOOM_XY_LIM, ZOOM_XY_LIM; length = ZOOM_GRID_N))

    # 8 perimeter targets — the four corners plus the four edge midpoints
    # of the zoom square.  Cheap to walk to, and they guarantee Stage-1
    # spans every corner of the zoom region before Stage-2 fills the
    # interior by extrapolation from the visited node tree.
    L = ZOOM_XY_LIM
    targets = ComplexF64[
        -L - L*im,   L - L*im,   -L + L*im,   L + L*im,    # corners
        -L + 0.0im,  L + 0.0im,  0.0 - L*im,  0.0 + L*im   # edge midpoints
    ]

    # Build the fine grid column-major so `walk.grid_y[k]` corresponds
    # to `(xs[i], ys[j])` with `k = (j - 1) · N + i`.  Then
    # `reshape(u_flat, N, N)` gives `u[i, j] = V₀(xs[i], ys[j])` — the
    # Makie heatmap convention.  Built explicitly rather than via a
    # comprehension to make the index mapping unambiguous.
    fine_grid = Vector{ComplexF64}(undef, ZOOM_GRID_N * ZOOM_GRID_N)
    k = 0
    for j in 1:ZOOM_GRID_N, i in 1:ZOOM_GRID_N
        k += 1
        fine_grid[k] = ComplexF64(xs[i], ys[j])
    end

    @printf("  [zoom %5.1fs] vector_path_network_solve (8 targets, %d fine grid)...\n",
            time() - t0, length(fine_grid)); flush(stdout)
    walk = vector_path_network_solve(prob, targets;
        order             = ZOOM_PN_ORDER,
        h                 = ZOOM_PN_H,
        adaptive          = true,
        skip_covered      = true,
        fine_grid         = fine_grid,
        extrapolate       = true,
        tol               = ZOOM_PN_TOL,
        on_target_failure = :skip,
        verbose           = true)

    @printf("  [zoom %5.1fs] walk done — %d visited nodes, %d failed targets\n",
            time() - t0, length(walk.visited_z),
            length(walk.failed_targets)); flush(stdout)

    # `walk.grid_y[k]` is a length-4 `Vector{ComplexF64}` (the
    # `P_I^(2)` companion-system state at fine_grid[k]).  Component 1
    # is `V₀`.  Reshape directly: column-major fine_grid → `u[i, j] =
    # V₀(xs[i], ys[j])`.  No permutedims needed.
    u_flat = [g[1] for g in walk.grid_y]
    u      = reshape(u_flat, ZOOM_GRID_N, ZOOM_GRID_N)

    msg = @sprintf("zoom kernel: %.1fs — %d visited nodes, %d/%d targets failed; grid %d×%d over [-%.1f, %.1f]²",
                   time() - t0, length(walk.visited_z),
                   length(walk.failed_targets), length(targets),
                   ZOOM_GRID_N, ZOOM_GRID_N, ZOOM_XY_LIM, ZOOM_XY_LIM)
    @printf("  %s\n", msg); flush(stdout)

    return (xs = xs, ys = ys,
            Re_u = real.(u), Im_u = imag.(u),
            walk = walk, message = msg)
end
