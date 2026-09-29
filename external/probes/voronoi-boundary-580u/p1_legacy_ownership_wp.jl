# Bead 580u measurement (worklog 084). Usage: julia --project=. <this> HALF ORDER [pI]
# The `kp` list below re-derives the LEGACY per-estimate ownership rule.
# Bead 580u measurement: does windowed_extract_poles double-count / drop a
# pole straddling a Voronoi core line?  ℘ equianharmonic lattice oracle
# (same fixture as external/probes/fig47-seam-diagnosis/p6_ingn_pole_truth.jl).
using PadeTaylor, Printf
const WT = PadeTaylor.WindowedTiling
f(z, u, up) = 6u^2
u0, up0 = 1.071822516416917, 1.710337353176786
HALF = parse(Float64, get(ARGS, 1, "30")); N = Int(2HALF) + 1
ORD = parse(Int, get(ARGS, 2, "30")); PI = get(ARGS, 3, "wp") == "pI"
xs = range(-HALF, HALF; length = N); ys = xs
prob = PI ? PadeTaylorProblem((z,u,up)->6u^2+z, (-0.1875, 0.3049), (0.0, HALF); order = ORD) : PadeTaylorProblem(f, (u0, up0), (0.0, HALF); order = ORD)
ω = 1.3629079683730673; g1 = 2ω + 0im; g2 = 2ω * cis(pi / 3)
fulllat = ComplexF64[1 + m * g1 + n * g2 for m in -40:40 for n in -40:40]
exact = [p for p in fulllat if abs(real(p)) ≤ HALF && abs(imag(p)) ≤ HALF]

t = time()
for seed in (0, 42)
wsol = windowed_path_network_solve(prob, xs, ys; window_extent = 20.0,
         overlap = 6.0, h = 0.5, order = ORD, rng_seed = seed)
@printf("seed=%d solve %.0fs  windows=%d\n", seed, time() - t, length(wsol.centers))
C = wsol.centers
kp = ComplexF64[]; kw = Int[]
allp = ComplexF64[]; allw = Int[]
for wi in eachindex(wsol.window_sols)
    for p in extract_poles(wsol.window_sols[wi])
        push!(allp, p); push!(allw, wi)
        WT._nearest_center(p, C) == wi && (push!(kp, p); push!(kw, wi))
    end
end
# (pre-580u: @assert kp == windowed_extract_poles(wsol) — the legacy rule WAS the driver)
function coredist(p)
    d = sort([(abs(p - c), k) for (k, c) in enumerate(C)])
    c1 = C[d[1][2]]; c2 = C[d[2][2]]
    (d[2][1]^2 - d[1][1]^2) / (2abs(c1 - c2))
end
for r in (0.1, 0.3, 0.5, 1.0)
    npair = 0; nsame = 0
    for i in eachindex(kp), j in i+1:length(kp)
        abs(kp[i] - kp[j]) ≤ r || continue
        kw[i] == kw[j] ? (nsame += 1) : (npair += 1)
        r == 0.3 && @printf("    pair w%d:%s  w%d:%s\n", kw[i], string(round(kp[i];digits=4)), kw[j], string(round(kp[j];digits=4)))
    end
    @printf("  r=%.1f  cross-window kept pairs=%d  same-window kept pairs=%d\n", r, npair, nsame)
end
# oracle: multiplicity of each exact pole among kept poles (tol 0.3)
tol = 0.3
mult = [count(p -> abs(p - e) < tol, kp) for e in exact]
spur = count(p -> minimum(abs(p - e) for e in fulllat) ≥ tol, kp)
maxerr = maximum(minimum(abs(p - e) for e in fulllat) for p in kp)
@printf("  kept=%d exact=%d  mult0=%d mult1=%d mult≥2=%d  spurious=%d  max|kept-exact|=%.2e\n",
        length(kp), length(exact), count(==(0), mult), count(==(1), mult),
        count(≥(2), mult), spur, maxerr)
# drops: exact poles with 0 kept; of those, were they seen by some window (dropped by ownership)?
for (k, e) in enumerate(exact)
    mult[k] == 0 || continue
    seen = [(allw[i], allp[i]) for i in eachindex(allp) if abs(allp[i] - e) < tol]
    @printf("  DROP exact=%s coredist=%.3f seen_by=%s\n", string(round(e; digits=3)), coredist(e), string(seen))
end
# the boundary band: exact poles within 0.5 of a core line, and how they were resolved
band = [k for (k, e) in enumerate(exact) if coredist(e) < 0.5]
bmult = mult[band]
@printf("  core-line band (<0.5): exact=%d  mult0=%d mult1=%d mult≥2=%d\n",
        length(band), count(==(0), bmult), count(==(1), bmult), count(≥(2), bmult))
# per-window estimate spread for band poles resolved by ≥2 windows
spreads = Float64[]; nmultiwin = 0
for k in band
    e = exact[k]
    ests = [allp[i] for i in eachindex(allp) if abs(allp[i] - e) < tol]
    wins = unique([allw[i] for i in eachindex(allp) if abs(allp[i] - e) < tol])
    length(wins) ≥ 2 && (nmultiwin += 1; push!(spreads, maximum(abs(a - b) for a in ests for b in ests)))
end
@printf("  band poles seen by ≥2 windows=%d  max inter-window estimate spread=%.2e  min coredist of band=%.2e\n",
        nmultiwin, isempty(spreads) ? NaN : maximum(spreads), minimum(coredist.(exact[band])))
fx = windowed_extract_poles(wsol)   # post-580u cluster-then-assign
fm = [count(p -> abs(p - e) < tol, fx) for e in exact]
@printf("  FIXED rule: kept=%d mult0=%d mult1=%d mult≥2=%d\n",
        length(fx), count(==(0), fm), count(==(1), fm), count(≥(2), fm))
end
