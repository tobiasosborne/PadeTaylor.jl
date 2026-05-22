# figures/probe_wedge_coverage.jl
#
# ADR-0026 D2/D3 — the dense-wedge coverage PROBE.
#
# This is a *measurement* script, not a figure renderer.  ADR-0026 D2
# replaces the headline figure's sparse 171-point polar fan with a
# dense, disc-spaced Cartesian target lattice (`surf_wedge_dense_targets`)
# and D1's resilient walk (`on_target_failure = :skip`).  D2 explicitly
# mandates that the target spacing `s` is a *measured* parameter: render
# at a starting `s`, measure the achieved honest coverage fraction, and
# tighten `s` until coverage saturates or runtime forces a stop.
#
# Rendering the full `kkg_pi2_surface()` (401² Cartesian grid + three
# Region-1 voters + Gridap FEM + Makie) is far too slow to use as the
# measurement loop.  This probe isolates the Region-2 wedge walk: it
# builds the *exact* `VectorPadeTaylorProblem` `surf_wedge_fill` builds
# (the BVP-anchored seed at `SURF_Z_SEED = -3`), runs the resilient
# `vector_path_network_solve` over `surf_wedge_dense_targets()` at each
# of several spacings, and reports — per spacing — the honest coverage
# fraction over a representative (coarse) wedge grid plus the
# resilient-walk failure accounting.
#
# Coverage is an *area fraction*, so a coarser-than-401² fine grid is
# adequate: the probe uses a ~200-point-per-side Cartesian grid clipped
# to the wedge.  The coverage figure is the fraction of those wedge grid
# cells the B1 true-radius gate honestly carries (`extrapolate = false`,
# so no Padé is evaluated outside its verified disc — an honest gap,
# never an extrapolated lie).
#
# ## What the probe answers (ADR-0026 D3)
#
# D3: "a vector walk that skips a non-trivial fraction is exhibiting a
# bug; FW's near-100 % target success is the benchmark."  So the probe
# reports, per spacing, the number of *skipped* targets and the
# breakdown of `VectorWalkFailure.reason`.  If even the COARSEST spacing
# (fewest targets, fastest) skips a large fraction of targets, that is a
# deeper bug — the D3 FW-lens investigation — and there is no point
# grinding finer spacings; the probe prints a STOP banner and the
# operator reads the failure-reason breakdown.
#
# Run (the figures project dev-links the main package, so it sees the
# `on_target_failure` kwarg):
#
#   cd /home/tobias/Projects/PadeTaylor.jl
#   julia --project=figures figures/probe_wedge_coverage.jl
#
# The probe runs the coarsest spacing first (fastest) for an early
# signal, then tightens.  It may take many minutes — expected.

using Printf

include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))

using PadeTaylor.VectorPathNetwork: vector_path_network_solve

# --- Probe parameters --------------------------------------------------

# Spacings to sweep, COARSEST FIRST (fewest targets ⇒ fastest ⇒ the
# earliest signal — ADR-0026 D3's "stop early if the coarsest fails").
const PROBE_SPACINGS = (0.35, 0.25, 0.18)

# The coverage fine grid: a Cartesian grid clipped to the wedge.  ~200
# points per side over `[-20, 20]` — coarse vs the figure's 401², but
# coverage is an area fraction, so this resolves it well enough.
const PROBE_GRID_N = 200

# The "large fraction skipped" threshold (ADR-0026 D3).  If the coarsest
# spacing skips more than this fraction of its targets, the probe stops:
# a non-trivial skip fraction is a bug, not a coverage knob.
const PROBE_STOP_SKIP_FRAC = 0.25

"""
    probe_wedge_grid() -> Vector{ComplexF64}

The coverage measurement grid: a `PROBE_GRID_N × PROBE_GRID_N` Cartesian
lattice on `[-20, 20]²`, clipped to exactly the wedge cells
`kkg_pi2_surface` would fill — inside the `|x| ≤ 20` disc, in the wedge
(`surf_in_sector` false), off the `±1°` Stokes mask, not the origin.
The achieved coverage fraction is `count(covered) / length(grid)`.
"""
function probe_wedge_grid()
    xs = range(-SURF_XY_LIM, SURF_XY_LIM; length = PROBE_GRID_N)
    ys = range(-SURF_XY_LIM, SURF_XY_LIM; length = PROBE_GRID_N)
    pts = ComplexF64[]
    for y in ys, x in xs
        z = ComplexF64(x, y)
        abs(z) > SURF_XY_LIM && continue
        surf_in_mask(z)      && continue
        surf_in_sector(z)    && continue
        iszero(z)            && continue
        push!(pts, z)
    end
    return pts
end

"""
    probe_one(bvp_sol, grid_pts, s) -> NamedTuple

Run the Region-2 wedge walk once at lattice spacing `s` and measure
coverage.  Builds the `VectorPadeTaylorProblem` exactly as
`surf_wedge_fill` does (the BVP-anchored seed at `SURF_Z_SEED`), runs
the resilient `vector_path_network_solve` (`on_target_failure = :skip`)
over `surf_wedge_dense_targets(; s = s)` with `grid_pts` as the Stage-2
`fine_grid`, and returns the row data: target / node / skip counts, the
`VectorWalkFailure.reason` breakdown, the honest coverage fraction, and
the wall-time.
"""
function probe_one(bvp_sol, grid_pts, s)
    f      = painleve_hierarchy(:I, 2; t = SURF_T)
    y_seed = ComplexF64.(bvp_sol(ComplexF64(SURF_Z_SEED)))
    prob   = VectorPadeTaylorProblem(f, y_seed,
                                     (ComplexF64(SURF_Z_SEED),
                                      ComplexF64(20.0 + 0.0im));
                                     order = SURF_PN_ORDER)
    targets = surf_wedge_dense_targets(; s = s)

    t0 = time()
    walk = vector_path_network_solve(prob, targets;
                                     order             = SURF_PN_ORDER,
                                     h                 = SURF_PN_H,
                                     fine_grid          = grid_pts,
                                     extrapolate        = false,
                                     tol                = SURF_PN_TOL,
                                     on_target_failure  = :skip)
    wall = time() - t0

    # Honest coverage: a grid slot is covered iff its Stage-2 `grid_y`
    # carries a finite first component (no `NaN` — the B1 gate let it
    # through).  Mirrors `surf_wedge_fill`'s `covered` mask exactly.
    covered  = Bool[!any(isnan, gy) for gy in walk.grid_y]
    cov_frac = isempty(covered) ? 0.0 : count(covered) / length(covered)

    # Failure-reason breakdown over `walk.failed_targets`.
    reasons = Dict{Symbol,Int}()
    for fail in walk.failed_targets
        reasons[fail.reason] = get(reasons, fail.reason, 0) + 1
    end

    return (s = s, n_targets = length(targets),
            n_nodes = length(walk.visited_z),
            n_failed = length(walk.failed_targets),
            reasons = reasons, cov_frac = cov_frac,
            n_grid = length(covered), wall = wall)
end

reasons_str(reasons) =
    isempty(reasons) ? "none" :
    join(["$k=$v" for (k, v) in sort(collect(reasons); by = first)], " ")

function probe_print_row(r)
    @printf("  s=%.2f | targets=%5d | nodes=%6d | failed=%5d (%.1f%%) | %s | cov=%.4f (%d/%d) | %.1fs\n",
            r.s, r.n_targets, r.n_nodes, r.n_failed,
            100 * r.n_failed / max(r.n_targets, 1),
            reasons_str(r.reasons), r.cov_frac,
            round(Int, r.cov_frac * r.n_grid), r.n_grid, r.wall)
end

function main()
    println("ADR-0026 D2/D3 — dense-wedge coverage probe")
    println("=" ^ 72)
    println("Building the Region-2 BVP anchor (surf_anchor_bvp) ...")
    t_anchor = @elapsed bvp_sol = surf_anchor_bvp()
    @printf("  anchor BVP solved in %.1fs\n", t_anchor)

    grid_pts = probe_wedge_grid()
    @printf("  coverage grid: %d wedge points (%d-per-side Cartesian)\n",
            length(grid_pts), PROBE_GRID_N)
    println("-" ^ 72)
    println("Sweeping spacings (coarsest first) ...")

    rows = NamedTuple[]
    for (k, s) in enumerate(PROBE_SPACINGS)
        @printf("[%d/%d] spacing s = %.2f ...\n", k, length(PROBE_SPACINGS), s)
        r = probe_one(bvp_sol, grid_pts, s)
        push!(rows, r)
        probe_print_row(r)

        # ADR-0026 D3 early-stop: if the COARSEST spacing skips a large
        # fraction of targets, that is a deeper bug — stop, do not grind
        # finer spacings.
        if k == 1
            skip_frac = r.n_failed / max(r.n_targets, 1)
            if skip_frac > PROBE_STOP_SKIP_FRAC
                println("-" ^ 72)
                @printf("STOP (ADR-0026 D3): the coarsest spacing skipped %.1f%% of targets\n",
                        100 * skip_frac)
                println("  (> $(round(Int, 100*PROBE_STOP_SKIP_FRAC))% threshold).  A non-trivial")
                println("  skip fraction is a bug, not a coverage knob — flag for the")
                println("  D3 FW-lens walk-failure investigation.  Reason breakdown:")
                println("    ", reasons_str(r.reasons))
                println("  Not grinding finer spacings.")
                return rows
            end
        end
    end

    println("-" ^ 72)
    println("Summary — spacing / targets / nodes / failed / reasons / cov / wall:")
    for r in rows
        probe_print_row(r)
    end
    return rows
end

main()
