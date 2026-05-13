# external/probes/pathnetwork-long-range/capture.jl
#
# Offline probe — confirm `u(z=10⁴) ≈ 21.02530339471055` (Fornberg &
# Weideman 2011, J. Comput. Phys. 230, Table 5.1 column (b),
# references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:299-301 and :385-391) via the
# Phase-10 path-network solver at `Float64` and `BigFloat`-256.
#
# Multi-z sweep — z ∈ {30, 100, 500, 1000, 10000}.  Lower-z runs land
# fast and let us calibrate scaling before committing to the z=10⁴ tail.
# Each row prints {wall, visited_nodes, |imag|, rel-err vs FW reference}.
#
# Routine test coverage at `Float64` lives in `test/pathnetwork_test.jl`
# PN.2.3 (~20 s wall, rtol = 5e-5).  The BF-256 confirmation is
# expensive — `~0.5 s/node × n_nodes` at 8 threads.  Wall envelope:
#   z=30    →    75 nodes ≈     40 s
#   z=100   →   250 nodes ≈    130 s
#   z=500   →  1200 nodes ≈    600 s   (~10 min)
#   z=1000  →  2400 nodes ≈   1200 s   (~20 min)
#   z=10000 → 24000 nodes ≈  12000 s   (~3.5 h)
# Total sweep wall: ~4 h.  Re-runnable on demand.
#
# Re-run:
#
#     JULIA_NUM_THREADS=8 julia --project=. \
#         external/probes/pathnetwork-long-range/capture.jl \
#         > external/probes/pathnetwork-long-range/oracles.txt 2>&1
#
# See `docs/worklog/019-fw-table-51-z10000.md` for the captured result
# and the bead `padetaylor-g9x` close note.

using PadeTaylor
using Printf
using Base.Threads: nthreads

fW(z, u, up) = 6u^2

println("=" ^ 78)
println("PathNetwork long-range probe — u(z) of ℘(z-1; 0, 2) trajectory")
println("Reference: Fornberg & Weideman 2011 J. Comput. Phys. 230, Table 5.1.")
println("FW2011_*.md:299-301 (reference values) + :385-391 (table itself).")
println("Parameters: h = 0.5, order = 30 (FW canonical recipe, §5.1.2).")
@printf("Julia threads: %d\n", nthreads())
println("=" ^ 78)
println()

# FW reference values for u(z) on the equianharmonic-℘ trajectory.
# z=30, z=10000 come directly from FW 2011 line 301 (the (5.3) numbers).
# z=100, z=500, z=1000 are NOT in FW Table 5.1; we cross-check the
# path-network result at those points against the FW values pinned at
# z=30 and z=10000 by inspecting the rel-err scaling (should grow
# linearly in step count for Float64; should stay ≤ 1e-13 for BF-256).
const TARGETS = [
    (30.0,    big"1.0950982559597442"),          # FW (5.3) line 301
    (100.0,   nothing),                           # no FW pin; report only
    (500.0,   nothing),
    (1000.0,  nothing),
    (10000.0, big"21.02530339471055"),           # FW (5.3) line 301 (b)
]

# ──────────────────────────────────────────────────────────────────────
# Float64 sweep — algorithmic-stability check + visited-node-count
# baseline.  Bounded by Float64 roundoff at long range.
# ──────────────────────────────────────────────────────────────────────
println("Float64 sweep")
println("─" ^ 78)
@printf("%-10s %12s %10s %15s %15s %15s\n",
        "target z", "wall (s)", "n_nodes", "u(z)", "|imag|", "rel-err vs FW")
println("─" ^ 78)

let u0 = 1.071822516416917, up0 = 1.710337353176786
    for (z_target, ref_big) in TARGETS
        prob = PadeTaylorProblem(fW, (u0, up0), (0.0, z_target); order = 30)
        t0   = time()
        sol  = path_network_solve(prob,
                                   ComplexF64[ComplexF64(z_target)];
                                   h = 0.5,
                                   max_steps_per_target = 200_000)
        wall = time() - t0
        u    = sol.grid_u[1]
        rel  = ref_big === nothing ? NaN :
               Float64(abs(u - Float64(ref_big)) / abs(Float64(ref_big)))
        n    = length(sol.visited_z)
        @printf("%-10.0f %12.2f %10d %+15.6e %15.4e %15.4e\n",
                z_target, wall, n, real(u), abs(imag(u)), rel)
        flush(stdout)
    end
end
println()

# ──────────────────────────────────────────────────────────────────────
# BigFloat-256 sweep — headline confirmation.  Bound: rel-err ≤ 1e-13
# (PN.2.2 z=30 standard).  Wall dominated by the z=10⁴ tail.
# ──────────────────────────────────────────────────────────────────────
println("BigFloat-256 sweep   (verbose=true, progress_every=500)")
println("─" ^ 78)
setprecision(BigFloat, 256) do
    u0  = big"1.071822516416917"
    up0 = big"1.710337353176786"
    for (z_target, ref_big) in TARGETS
        zT_big = big(z_target)
        prob = PadeTaylorProblem(fW, (u0, up0),
                                  (big(0.0), zT_big); order = 30)
        target = Complex{BigFloat}[Complex{BigFloat}(zT_big)]
        @printf("\n>>> z=%.0f BF-256\n", z_target)
        flush(stdout)
        t0  = time()
        sol = path_network_solve(prob, target;
                                  h = big(0.5),
                                  max_steps_per_target = 200_000,
                                  verbose = true,
                                  progress_every = 500)
        wall = time() - t0
        u    = sol.grid_u[1]
        n    = length(sol.visited_z)
        if ref_big === nothing
            @printf("    z=%-8.0f wall=%.1fs  n=%d  u=%s  |imag|=%.4e  (no FW pin)\n",
                    z_target, wall, n, repr(Float64(real(u))),
                    Float64(abs(imag(u))))
        else
            rel = abs(u - ref_big) / abs(ref_big)
            @printf("    z=%-8.0f wall=%.1fs  n=%d  u=%s\n",
                    z_target, wall, n, repr(Float64(real(u))))
            @printf("                       |imag|=%.4e  rel-err vs FW=%.4e  (≤1e-13 ? %s)\n",
                    Float64(abs(imag(u))), Float64(rel),
                    rel ≤ big"1e-13" ? "YES" : "NO")
        end
        flush(stdout)
    end
end
println()
println("=" ^ 78)
println("Sweep complete.")
