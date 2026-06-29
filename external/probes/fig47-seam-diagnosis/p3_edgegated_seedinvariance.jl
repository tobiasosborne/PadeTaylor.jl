# P3 — THE DECISIVE MEASUREMENT (bead vwgl / sny7, worklog 077 follow-up).
#
# Every prior probe (p2_jumpmap_and_seeddiff.jl, seed_poles.jl) measured the
# PLAIN path_network_solve over the FULL grid. But the FIGURE (fw2011_fig_4_7.jl)
# renders panel (e) via `edge_gated_pole_field_solve` — the IVP CONFINED to the
# pole field (FW md:401) — at the deterministic default seed 0. So the figure's
# seam is a fixed seed-0 artifact, and nobody has tested whether it is
# seed-DEPENDENT through the edge-gated confinement.
#
# DECISIVE QUESTION: does confining the walk to the pole field (what the figure
# already does) COLLAPSE the seed-variance (=> the seam is a SMOOTH-SECTOR
# leakage artifact, regime A, already cured for the pole scatter), or does the
# grain boundary SURVIVE inside the pole field (=> a genuine pole-field-INTERIOR
# grain boundary, regime B, the hard residual edge-gating does NOT fix)?
#
# edge_gated_pole_field_solve does NOT expose rng_seed (its internal
# _solve_targets calls path_network_solve with no seed => seed 0), so a naive
# two-seed call through it is VACUOUS. Correct experiment: take the seed-0
# confinement mask, then re-solve THAT confined target set at seeds 0 and 42.

using PadeTaylor, Printf
const EGS = PadeTaylor.EdgeGatedSolve

pI(z, u, up) = 6u^2 + z
const HALF = 50.0; const N = 101
xs = range(-HALF, HALF; length = N)
ys = range(-HALF, HALF; length = N)
# Panel (e), exactly as fw2011_fig_4_7.jl PANELS[5].
prob = PadeTaylorProblem(pI, (-0.1875, 0.3049), (0.0, HALF); order = 30)

# Pole-set agreement metric: for each pole in A, nearest in B; report matched
# fraction within `tol`, median & max nearest-neighbour distance.
function poleset_agreement(A, B; tol = 0.5)
    isempty(A) && return (0.0, NaN, NaN)
    nn = Float64[]
    for p in A
        d = Inf
        for q in B
            dd = abs(p - q); dd < d && (d = dd)
        end
        push!(nn, d)
    end
    sort!(nn)
    matched = count(<=(tol), nn) / length(nn)
    med = nn[(length(nn) + 1) ÷ 2]
    (matched, med, maximum(nn))
end

# Two-seed field diff over a shared cell list (same z order both seeds).
function field_seeddiff(u0, u1, zs, poles; poleguard = 1.2)
    dmax = 0.0; nd = 0; nt = 0; zmx = 0.0 + 0im
    for k in eachindex(zs)
        a = u0[k]; b = u1[k]
        (isfinite(a) && isfinite(b)) || continue
        any(p -> abs(p - zs[k]) < poleguard, poles) && continue  # off-pole only
        nt += 1
        rel = abs(a - b) / (abs(a) + abs(b) + 1e-30)
        rel > dmax && (dmax = rel; zmx = zs[k])
        rel > 1e-2 && (nd += 1)
    end
    (dmax, nd, nt, zmx)
end

println("="^72)
println("ARM B — EDGE-GATED CONFINED solve, two-seed (THE DECISIVE NEW ARM)")
println("="^72); flush(stdout)

tb = time()
egs = edge_gated_pole_field_solve(prob, xs, ys; h = 0.5, order = 30, grow_rings = 4)
@printf("edge-gated: %d passes, %d pole-field cells, %d visited nodes; %.1f s\n",
        egs.iterations, count(egs.field_mask), length(egs.pn_solution.visited_z), time() - tb)

# Confined target set = field_mask dilated by one ring (== the driver's final solve set).
final_targets = EGS._dilate(egs.field_mask, 1)
zs_conf = ComplexF64[]
for j in 1:N, i in 1:N
    final_targets[i, j] && push!(zs_conf, complex(xs[i], ys[j]))
end
@printf("confined target cells: %d (of %d)\n", length(zs_conf), N * N); flush(stdout)

solB = Dict{Int,Any}()
for seed in (0, 42)
    t0 = time()
    pn = path_network_solve(prob, zs_conf; h = 0.5, order = 30, rng_seed = seed, extrapolate = true)
    ps = extract_poles(pn)
    solB[seed] = (pn, ps)
    @printf("  [confined] seed=%d nodes=%d poles=%d  %.1f s\n",
            seed, length(pn.visited_z), length(ps), time() - t0); flush(stdout)
end
(pn0B, ps0B) = solB[0]; (pn42B, ps42B) = solB[42]
dmaxB, ndB, ntB, zmxB = field_seeddiff(pn0B.grid_u, pn42B.grid_u, zs_conf, ps0B)
matchB, medB, maxnnB = poleset_agreement(ps0B, ps42B)
println("\n-- ARM B RESULT (edge-gated / pole-field-confined) --")
@printf("pole COUNT  seed0=%d  seed42=%d   |Δ|=%d\n", length(ps0B), length(ps42B), abs(length(ps0B) - length(ps42B)))
@printf("pole-set agreement (tol 0.5): matched=%.1f%%  median NN=%.3f  max NN=%.3f\n", 100matchB, medB, maxnnB)
@printf("two-seed FIELD diff (off-pole): max rel=%.3e at z=%+.1f%+.1fi ; cells>1e-2 = %d/%d (%.1f%%)\n",
        dmaxB, real(zmxB), imag(zmxB), ndB, ntB, 100ndB / max(ntB, 1)); flush(stdout)

println("\n" * "="^72)
println("ARM A — PLAIN FULL-GRID solve, two-seed (BASELINE, re-confirm 8159/7987)")
println("="^72); flush(stdout)

grid = ComplexF64[complex(x, y) for y in ys for x in xs]
solA = Dict{Int,Any}()
for seed in (0, 42)
    t0 = time()
    pn = path_network_solve(prob, grid; h = 0.5, order = 30, rng_seed = seed, extrapolate = true)
    ps = extract_poles(pn)
    solA[seed] = (pn, ps)
    @printf("  [full] seed=%d nodes=%d poles=%d  %.1f s\n",
            seed, length(pn.visited_z), length(ps), time() - t0); flush(stdout)
end
(pn0A, ps0A) = solA[0]; (pn42A, ps42A) = solA[42]
dmaxA, ndA, ntA, zmxA = field_seeddiff(pn0A.grid_u, pn42A.grid_u, grid, ps0A)
matchA, medA, maxnnA = poleset_agreement(ps0A, ps42A)
println("\n-- ARM A RESULT (plain full-grid baseline) --")
@printf("pole COUNT  seed0=%d  seed42=%d   |Δ|=%d\n", length(ps0A), length(ps42A), abs(length(ps0A) - length(ps42A)))
@printf("pole-set agreement (tol 0.5): matched=%.1f%%  median NN=%.3f  max NN=%.3f\n", 100matchA, medA, maxnnA)
@printf("two-seed FIELD diff (off-pole): max rel=%.3e at z=%+.1f%+.1fi ; cells>1e-2 = %d/%d (%.1f%%)\n",
        dmaxA, real(zmxA), imag(zmxA), ndA, ntA, 100ndA / max(ntA, 1)); flush(stdout)

println("\n" * "="^72)
println("VERDICT  (A = plain full grid, B = edge-gated pole-field-confined)")
println("="^72)
@printf("pole-count seed-variance:  A |Δ|=%d   B |Δ|=%d\n", abs(length(ps0A) - length(ps42A)), abs(length(ps0B) - length(ps42B)))
@printf("two-seed field >1e-2:      A %.1f%%    B %.1f%%\n", 100ndA / max(ntA, 1), 100ndB / max(ntB, 1))
@printf("two-seed field max rel:    A %.2e  B %.2e\n", dmaxA, dmaxB)
println("""
INTERPRETATION:
  B ~ A   => grain boundary SURVIVES confinement = pole-field-INTERIOR (regime B,
            the hard residual; edge-gating does NOT cure the visible Fig 4.7 seam).
  B << A  => confinement COLLAPSES seed-variance = SMOOTH-SECTOR leakage (regime A;
            edge-gating already cures the pole scatter; figure seam is small/elsewhere).
""")
