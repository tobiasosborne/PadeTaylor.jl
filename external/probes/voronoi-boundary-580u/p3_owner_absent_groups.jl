# 580u review follow-up (worklog 084 §6): count _own_poles groups whose
# owner (_nearest_center of the group mean) has NO member in the group;
# split singletons vs multi-window; check ℘ discards against the exact lattice.
using PadeTaylor, Printf
using PadeTaylor.WindowedTiling: _nearest_center
ω = 1.3629079683730673; g1 = 2ω + 0im; g2 = 2ω * cis(pi / 3)
lat = ComplexF64[1 + m * g1 + n * g2 for m in -40:40 for n in -40:40]
function groups(ps, atol)
    z = reduce(vcat, ps); win = reduce(vcat, [fill(i, length(ps[i])) for i in eachindex(ps)])
    n = length(z); par = collect(1:n)
    root(a) = (while par[a] != a; a = par[a]; end; a)
    for a in 1:n, b in a+1:n
        (win[a] != win[b] && abs(z[a]-z[b]) ≤ atol) || continue
        ra, rb = root(a), root(b); ra == rb || (par[rb] = ra)
    end
    G = Dict{Int,Vector{Int}}(); for a in 1:n; push!(get!(G, root(a), Int[]), a); end
    z, win, collect(values(G))
end
function run(name, prob, H, ord, oracle)
    xs = range(-H, H; length = Int(2H)+1)
    w = windowed_path_network_solve(prob, xs, xs; window_extent=20.0, overlap=6.0, h=0.5, order=ord, rng_seed=0)
    ps = [extract_poles(s) for s in w.window_sols]
    z, win, G = groups(ps, 0.1)
    sing = 0; multi = 0; singtrue = 0; multitrue = 0; maxdiam = 0.0
    for g in G
        maxdiam = max(maxdiam, maximum(abs(z[a]-z[b]) for a in g for b in g))
        o = _nearest_center(sum(z[g]) / length(g), w.centers)
        any(win[a] == o for a in g) && continue
        istrue = oracle && minimum(abs(z[g[1]] - e) for e in lat) < 0.3
        if length(g) == 1; sing += 1; singtrue += istrue
        else; multi += 1; multitrue += istrue
            @printf("   multi owner-absent: owner=%d wins=%s z=%s\n", o, string(win[g]), string(z[g[1]]))
        end
    end
    @printf("%s: groups=%d max group diameter=%.2e  owner-absent singletons=%d (true-lattice %s)  multi=%d (true-lattice %s)\n",
            name, length(G), maxdiam, sing, oracle ? string(singtrue) : "n/a", multi, oracle ? string(multitrue) : "n/a")
end
wp = H -> PadeTaylorProblem((z,u,up)->6u^2, (1.071822516416917, 1.710337353176786), (0.0, H); order = 30)
run("wp[-20,20]", wp(20.0), 20.0, 30, true)
run("wp[-30,30]", wp(30.0), 30.0, 30, true)
run("PI[-30,30]", PadeTaylorProblem((z,u,up)->6u^2+z, (-0.1875, 0.3049), (0.0, 30.0); order = 20), 30.0, 20, false)
