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
# ## The S5a re-probe (ADR-0026 Amendments 3–4 — the full corrected stack)
#
# ADR-0026 Amendment 2's Decision R1 (fixed `h`, `adaptive = false`) was
# the *interim* fix for the dense-walk collapse.  Amendment 2 root-caused
# that collapse: the old B2 `_adaptive_h` controller was a *geometric
# sink* — its pole cap `POLE_SAFETY·h_prev·min|t*|` was purely
# multiplicative, so `h` ratcheted to zero in a sustained pole field.
# Amendment 4's S1 diagnostic confirmed the sink and replaced the law
# with the absolute, non-ratcheting `h = clamp(SAFETY·D_local, h_min,
# h_max)` (no multiplicative `h_prev` term).  `adaptive = true` is now
# *correct*: the step is sized to the local honest B1 disc, tracking the
# `|x|⁻¹ᐟ⁶` pole-spacing scaling, instead of the scale-blind fixed
# `h = 0.1` (which is ~2× the honest disc — Amendment 3 bottleneck 1).
#
# This probe runs the **S5a re-probe**: the *full corrected stack* —
# `adaptive = true` (S2), `skip_covered = true` (S3, kills the
# covered-ground re-walking that chorded poles), scale-derived clustering
# (S4), dense targets over the *full* `|x| ≤ 20` wedge, swept across
# spacings.  The decisive question (per ADR-0026 Amendment 3): does the
# corrected stack climb coverage toward "filled" — the KKG Figs 7.4/7.5
# surface — or does it plateau at the ~12 % R1 ceiling?  The probe
# reports, per spacing, the honest coverage fraction **broken out by
# radial band** — inner `|x| ∈ [2,8]`, mid `[8,14]`, outer `[14,20]` —
# so a coverage that collapses with `|x|` is visible.
#
# ## What the probe answers (ADR-0026 D3)
#
# D3: "a vector walk that skips a non-trivial fraction is exhibiting a
# bug; FW's near-100 % target success is the benchmark."  So the probe
# reports, per spacing, the number of *skipped* targets and the
# breakdown of `VectorWalkFailure.reason`.
#
# Run (the figures project dev-links the main package, so it sees the
# `on_target_failure` and `adaptive` kwargs):
#
#   cd /home/tobias/Projects/PadeTaylor.jl
#   julia --project=figures figures/probe_wedge_coverage.jl
#
# The probe runs the coarsest spacing first (fastest) for an early
# signal, then tightens.  It may take 10–30 min per spacing — expected.

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

# Radial bands for the banded coverage breakdown (ADR-0026 Amendment 2,
# the R1 re-probe).  The open question: does a fixed `h` thread the
# outer wedge?  `[2,8]` inner — the dense inner pole field; `[8,14]` mid;
# `[14,20]` outer — where pole density grows toward the grid edge.
const PROBE_BANDS = ((2.0, 8.0, "inner [2,8]"),
                     (8.0, 14.0, "mid   [8,14]"),
                     (14.0, 20.0, "outer [14,20]"))

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
    banded_coverage(grid_pts, covered) -> Vector{NamedTuple}

Break the honest coverage fraction out by radial band (ADR-0026
Amendment 2, the R1 re-probe).  For each `(r_lo, r_hi, label)` in
`PROBE_BANDS`, count the grid points with `r_lo ≤ |x| < r_hi` and the
fraction of them the B1 gate honestly covers.  This is the measurement
that answers "does a fixed `h` thread the outer wedge, or does coverage
collapse with `|x|`?"
"""
function banded_coverage(grid_pts, covered)
    rows = NamedTuple[]
    for (r_lo, r_hi, label) in PROBE_BANDS
        in_band = [r_lo ≤ abs(grid_pts[k]) < r_hi
                   for k in eachindex(grid_pts)]
        n_band  = count(in_band)
        n_cov   = count(k -> in_band[k] && covered[k],
                        eachindex(grid_pts))
        frac    = n_band == 0 ? 0.0 : n_cov / n_band
        push!(rows, (label = label, n_band = n_band,
                     n_cov = n_cov, frac = frac))
    end
    return rows
end

"""
    probe_one(bvp_sol, grid_pts, s; adaptive, skip_covered, h) -> NamedTuple

Run the Region-2 wedge walk once at lattice spacing `s` and measure
coverage.  Builds the `VectorPadeTaylorProblem` exactly as
`surf_wedge_fill` does (the BVP-anchored seed at `SURF_Z_SEED`), runs
the resilient `vector_path_network_solve` (`on_target_failure = :skip`)
over `surf_wedge_dense_targets(; s = s)` with `grid_pts` as the Stage-2
`fine_grid`, and returns the row data: target / node / skip counts, the
`VectorWalkFailure.reason` breakdown, the honest coverage fraction (total
AND banded by radius — the S5a re-probe), and the wall-time.

`adaptive` / `skip_covered` / `h` select the walk's step strategy.  The
S5a re-probe runs the *full corrected stack* — `adaptive = true` (S2),
`skip_covered = true` (S3) — with `h = SURF_PN_H` seeding the adaptive
step *ceiling* `h_max` (ADR-0026 Amendments 3–4).
"""
function probe_one(bvp_sol, grid_pts, s;
                   adaptive::Bool = true, skip_covered::Bool = true,
                   h::Real = SURF_PN_H)
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
                                     h                 = h,
                                     adaptive           = adaptive,
                                     skip_covered       = skip_covered,
                                     step_policy        = :max_q_root,
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
    bands    = banded_coverage(grid_pts, covered)

    # Failure-reason breakdown over `walk.failed_targets`.
    reasons = Dict{Symbol,Int}()
    for fail in walk.failed_targets
        reasons[fail.reason] = get(reasons, fail.reason, 0) + 1
    end

    return (s = s, adaptive = adaptive, skip_covered = skip_covered,
            h = h,
            n_targets = length(targets),
            n_nodes = length(walk.visited_z),
            n_failed = length(walk.failed_targets),
            reasons = reasons, cov_frac = cov_frac,
            bands = bands,
            n_grid = length(covered), wall = wall)
end

reasons_str(reasons) =
    isempty(reasons) ? "none" :
    join(["$k=$v" for (k, v) in sort(collect(reasons); by = first)], " ")

function probe_print_row(r)
    mode = (r.adaptive ? "adapt" : "fixed") *
           (r.skip_covered ? "+skipcov" : "")
    @printf("  s=%.2f h=%.2f %s | targets=%5d | nodes=%6d | failed=%5d (%.1f%%) | %s\n",
            r.s, r.h, mode,
            r.n_targets, r.n_nodes, r.n_failed,
            100 * r.n_failed / max(r.n_targets, 1),
            reasons_str(r.reasons))
    @printf("           cov=%.4f (%d/%d total) | %.0fs (%.1f min)\n",
            r.cov_frac, round(Int, r.cov_frac * r.n_grid), r.n_grid,
            r.wall, r.wall / 60)
    for b in r.bands
        @printf("           band %-13s cov=%.4f (%d/%d)\n",
                b.label, b.frac, b.n_cov, b.n_band)
    end
end

function main()
    println("ADR-0026 D2/D3 + Amendments 3-4 — S5a full-corrected-stack coverage probe")
    println("=" ^ 72)
    println("Building the Region-2 BVP anchor (surf_anchor_bvp) ...")
    t_anchor = @elapsed bvp_sol = surf_anchor_bvp()
    @printf("  anchor BVP solved in %.1fs\n", t_anchor)

    grid_pts = probe_wedge_grid()
    @printf("  coverage grid: %d wedge points (%d-per-side Cartesian)\n",
            length(grid_pts), PROBE_GRID_N)
    # Report the per-band grid populations once up front.
    for b in banded_coverage(grid_pts, falses(length(grid_pts)))
        @printf("    band %-13s : %d grid points\n", b.label, b.n_band)
    end
    println("-" ^ 72)
    println("S5a re-probe — adaptive = true (S2), skip_covered = true (S3),")
    println("scale-derived clustering (S4), h_max = $(SURF_PN_H), FULL wedge.")
    println("Sweeping spacings (coarsest first) ...")

    rows = NamedTuple[]
    for (k, s) in enumerate(PROBE_SPACINGS)
        @printf("[%d/%d] spacing s = %.2f (adaptive, skip_covered, h_max = %.2f) ...\n",
                k, length(PROBE_SPACINGS), s, SURF_PN_H)
        flush(stdout)
        r = probe_one(bvp_sol, grid_pts, s;
                      adaptive = true, skip_covered = true, h = SURF_PN_H)
        push!(rows, r)
        probe_print_row(r)
        flush(stdout)
    end

    println("-" ^ 72)
    println("Summary — S5a re-probe (full corrected stack, banded coverage):")
    for r in rows
        probe_print_row(r)
    end

    return rows
end

main()
