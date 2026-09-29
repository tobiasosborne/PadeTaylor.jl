using PadeTaylor, Printf
using PadeTaylor.WindowedTiling: _nearest_center, _own_poles
xs = ys = range(-30.0, 30.0; length = 61)
prob = PadeTaylorProblem((z,u,up)->6u^2+z, (-0.1875, 0.3049), (0.0, 30.0); order = 20)
for seed in (0, 42)
    w = windowed_path_network_solve(prob, xs, ys; window_extent=20.0, overlap=6.0, h=0.5, order=20, rng_seed=seed)
    pw = [extract_poles(s) for s in w.window_sols]
    leg = ComplexF64[]; lw = Int[]
    for wi in eachindex(pw), p in pw[wi]; _nearest_center(p, w.centers) == wi && (push!(leg,p); push!(lw,wi)); end
    fix = windowed_extract_poles(w)
    x01 = count(((i,j),)-> lw[i]!=lw[j] && abs(leg[i]-leg[j]) ≤ 0.1, [(i,j) for i in eachindex(leg) for j in i+1:length(leg)])
    onlyleg = count(p -> !any(==(p), fix), leg); onlyfix = count(p -> !any(==(p), leg), fix)
    @printf("PI seed=%d legacy=%d fixed=%d  legacy xwin pairs≤0.1=%d  legacy-only=%d fixed-only=%d\n", seed, length(leg), length(fix), x01, onlyleg, onlyfix)
end
