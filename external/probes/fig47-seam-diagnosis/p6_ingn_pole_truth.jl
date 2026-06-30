# P6 (bead padetaylor-ingn) — does the windowed composite DROP REAL poles, or only
# suppress SPURIOUS ones? The composite finds ~220 fewer poles than monolithic on
# panel e; this verifies (Rule 9) those are spurious, not real, using the
# Weierstrass-℘ exact-lattice oracle (u''=6u², known equianharmonic lattice).
# ℘ is doubly-periodic (no smooth sectors → no bloom), so this cleanly isolates
# windowing's effect on REAL-pole preservation. Scale HALF=30 so windowing engages
# (3×3 windows). Acceptance: composite recall >= monolithic recall (no real poles
# dropped) AND composite precision >= monolithic precision (spurious suppressed).
using PadeTaylor, Printf
f(z, u, up) = 6u^2
u0, up0 = 1.071822516416917, 1.710337353176786
const HALF = 30.0; const N = 61
xs = range(-HALF, HALF; length = N); ys = range(-HALF, HALF; length = N)
prob = PadeTaylorProblem(f, (u0, up0), (0.0, HALF); order = 30)

# exact equianharmonic lattice p(m,n) = 1 + m·gen1 + n·gen2
ω = 1.3629079683730673; gen1 = 2ω + 0im; gen2 = 2ω*cis(pi/3)
exact = ComplexF64[]
for m in -20:20, n in -20:20
    p = 1 + m*gen1 + n*gen2
    -HALF ≤ real(p) ≤ HALF && -HALF ≤ imag(p) ≤ HALF && push!(exact, p)
end
fulllat = ComplexF64[1 + m*gen1 + n*gen2 for m in -24:24 for n in -24:24]
const tol = 0.3
recall(ps)    = count(e -> any(p -> abs(p-e) < tol, ps), exact) / length(exact)
precision(ps) = count(p -> any(e -> abs(p-e) < tol, fulllat), ps) / max(length(ps), 1)
spurious(ps)  = count(p -> minimum(abs(p-e) for e in fulllat) ≥ tol, ps)

@printf("exact ℘ poles in [-%g,%g]²: %d\n", HALF, HALF, length(exact))

t = time()
pm = extract_poles(path_network_solve(prob, ComplexF64[complex(x,y) for y in ys for x in xs];
                                      h = 0.5, order = 30, rng_seed = 0))
@printf("[monolithic %.0fs]\n", time()-t); flush(stdout)
t = time()
pc = windowed_extract_poles(windowed_path_network_solve(prob, xs, ys;
                            window_extent = 20.0, overlap = 6.0, h = 0.5, order = 30, rng_seed = 0))
@printf("[windowed   %.0fs]\n", time()-t); flush(stdout)

@printf("\n%-12s | poles | recall | precision | spurious\n", "method")
@printf("%-12s | %5d | %.3f  | %.3f     | %d\n", "MONOLITHIC", length(pm), recall(pm), precision(pm), spurious(pm))
@printf("%-12s | %5d | %.3f  | %.3f     | %d\n", "WINDOWED",   length(pc), recall(pc), precision(pc), spurious(pc))

drop_ok = recall(pc) ≥ recall(pm) - 0.005          # composite keeps the real poles
prec_ok = precision(pc) ≥ precision(pm) - 0.005     # and is at least as clean
println(drop_ok && prec_ok ?
        "\nINGN PASS: windowed keeps real poles (recall>=mono) AND is >=as clean (precision>=mono) => the fewer poles are SPURIOUS, not real" :
        "\nINGN FAIL: recall_ok=$drop_ok precision_ok=$prec_ok — investigate dropped real poles")
