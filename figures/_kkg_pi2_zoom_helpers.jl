# figures/_kkg_pi2_zoom_helpers.jl
#
# THE ZOOM-VARIANT KERNEL for the P_I^(2) tritronquée `V₀(x, 0)` on
# the **inner** complex-plane square `[-3, 3]² ⊂ |x| ≤ 3√2`.  Built as
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
# The zoom region `[-3, 3]²` is a square of side 6; the circumscribed
# circle has radius `3√2 ≈ 4.24`, so the entire zoom extends well
# beyond the inner disc `|x| < 3`.  More importantly, the smooth pole-free
# sector covers `|arg x| > 36°`; the wedge slice
# `Re x ∈ [√(9 − Im²), 3]` with `|arg x| < 36°` now fully intersects
# the zoom, and the tips of poles at `|x| ∈ [2, 3]` are visible in
# the wedge interior.  The bulk of `[-3, 3]²` is still in the smooth,
# pole-free sector — a single Padé walk from the BVP-anchored seed at
# `z = -3 + 0im` (which is **on** the left edge of the zoom square,
# `Re = -3`) integrates *inward* through smooth-sector physics, though
# walk steps crossing the wedge pole field may encounter difficult
# extrapolation near the poles.
#
# So the zoom kernel is **strictly simpler** than the surface kernel:
# one VectorPadeTaylorProblem, one `vector_path_network_solve` walk
# with the same `(order, h, tol)` knobs as the wedge of the R = 8
# kernel, and an 11 × 11 dense target lattice that guarantees Stage-1
# threads visited nodes through the entire zoom square before Stage-2
# fills the 1001² fine lattice by Padé evaluation at the nearest
# visited node — every fine-grid cell lands within ~0.05 units of a
# visited node, well inside the order-48 Padé's convergence disc.
#
# ## Resolution choice
#
# 1001 × 1001 = 1 002 001 grid points over the 6-unit-wide square gives a
# pixel size of `dx = 6 / 1000 = 0.006`, i.e. **~6.7 × finer than the R = 8
# cache's `dx = 0.04`**.  That is the resolution lift the user
# requested.  Stage-2 is `O(N_grid)` Padé evaluations at the
# nearest-visited node, so total runtime is dominated by the Stage-2
# fill (~30–120 s on a modern CPU); the kernel returns in under
# ~3 min.
#
# ## What this kernel intentionally does NOT do
#
#   * no sector ray-fan, no Laplace voters, no median vote — the wedge
#     poles at `|x| ∈ [2, 3]` appear in the [-3, 3]² zoom but are
#     rendered via the single Padé walk rather than a separate sector vote;
#   * no Stokes-strip mask — the strips at `±36°` start at `|x| = 2`
#     (they cross the zoom interior now), but the single-walk topology
#     does not apply a mask; poles appear as spikes in the rendered surface;
#   * `extrapolate = false` — at the [-3, 3]² zoom the seed-anchored
#     walks travel up to √45 ≈ 6.7 units, well beyond the per-node
#     Padé h = 0.1 convergence disc.  Raw `extrapolate = true` Padé
#     output past the disc produced SPURIOUS Q-roots (ghost poles)
#     in the analytic smooth sector — user-flagged at the [-3, 3]²
#     top/bottom edges (e.g. `(0, ±2.5)` at angle 90°, way outside
#     the wedge `|arg x| < 36°`).  Each Padé's Q has random zeros
#     in its analyticity region; these zeros do NOT represent real
#     poles of the underlying solution, and VC-4 cross-node
#     validation (in the parent surface kernel) is what filters them
#     out.  The zoom kernel skips VC-4 for simplicity, so we instead
#     pay the cost of `extrapolate = false` (cells past the nearest
#     Padé's disc become `NaN`) and recover a filled figure via
#     nearest-non-NaN BFS in the rendering scripts — interpolated
#     values from valid neighbours are MORE faithful than raw
#     past-disc Padé output.  Pole locations the walk actually
#     resolved (within disc) survive intact;
#   * no pole-field extraction — poles at `|x| ∈ [2, 3]` in the wedge
#     are visible but not separately catalogued in this kernel.
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
# origin.  `3.0` is the user-requested zoom radius (second iteration);
# the square `[-3, 3]²` extends beyond the inner disc `|x| ≤ 2` and
# crosses the wedge inner boundary `|x| = 2`, so the [-3, 3]² zoom
# shows the smooth-sector solution PLUS the inner-wedge poles (a few
# at `|x| ∈ [2, 3]`).
const ZOOM_XY_LIM = 3.0

# 1001 × 1001 = 1 002 001 grid points.  `dx = 6 / 1000 = 0.006`, i.e.
# ~6.7 × finer than the R = 8 cache's `dx = 0.04`.  Odd so a node sits
# on `x = 0` and `y = 0` (the negative-real ridge passes through both)
# AND on `x = -3` (the seed location — the IC pin).
const ZOOM_GRID_N = 1001

# Walk knobs — match the R = 8 wedge walk so the zoom is identical in
# *physics* to the surface kernel's wedge solve, only the *topology* is
# simpler (no wedge fan, just one seed → perimeter sweep).
const ZOOM_PN_ORDER = 48
const ZOOM_PN_H     = 0.05  # halved from 0.1 (bead `padetaylor-(zoom-h-half)`) to tighten each Padé convergence disc — cleaner in-disc values at the cost of more past-disc cells (the BFS fill rate may rise from 34% to ~60%, but the in-disc Padé evaluations are closer to convergence so the past-disc-Q-root artifacts the user flagged should weaken).
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
at `SURF_Z_SEED = -3 + 0im`.  Stage-1 walks to 121 dense lattice targets
(11×11 over [-3,3]², spacing 0.6) so the Stage-2 nearest-visited-node
fill stays within the order-48 Padé's convergence disc everywhere;
Stage-2 fills the 1001² fine grid by Padé evaluation at the nearest
visited node.

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

    # Dense target lattice — 21 × 21 = 441 targets evenly spaced over
    # [-ZOOM_XY_LIM, ZOOM_XY_LIM]² at spacing dr = 2·ZOOM_XY_LIM / 20 = 0.3.
    # Walks from the seed at z = -3 + 0im to all 121 targets thread the walk
    # through the entire zoom square; the walk's `skip_covered = true` +
    # `:max_q_root` adaptive stepping mean overlapping paths share nodes, so
    # the total visited-node count is a few thousand (not 121 × <steps/target>).
    # This guarantees the Stage-2 nearest-visited-node assignment has every
    # fine-grid cell within ~0.05 units of a visited node — well inside the
    # order-48 Padé's convergence disc, eliminating the Voronoi-cell rays
    # and ghost poles visible in the 8-target [-3, 3]² version (commit
    # d350ad5 + iteration above).  Bead `padetaylor-(zoom-densify)`.
    target_lattice_n = 21
    xs_t = range(-ZOOM_XY_LIM, ZOOM_XY_LIM; length = target_lattice_n)
    ys_t = range(-ZOOM_XY_LIM, ZOOM_XY_LIM; length = target_lattice_n)
    targets = ComplexF64[Complex{Float64}(x, y) for y in ys_t for x in xs_t]
    # Skip the seed-coincident target (z = -3 + 0im) — `vector_path_network_solve`
    # auto-skips it via the `10·eps` exact-coincidence check, but being explicit
    # documents intent.

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

    @printf("  [zoom %5.1fs] vector_path_network_solve (%d targets, %d fine grid)...\n",
            time() - t0, length(targets), length(fine_grid)); flush(stdout)
    walk = vector_path_network_solve(prob, targets;
        order             = ZOOM_PN_ORDER,
        h                 = ZOOM_PN_H,
        adaptive          = true,
        skip_covered      = true,
        fine_grid         = fine_grid,
        extrapolate       = false,
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
