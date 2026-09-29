# Probe 7 (bead rxbo): does a near-pole passage on the ancestor chain separate
# the default walk's bad far nodes from its good ones?  Chain signature:
# max ‖y‖ and min adaptive h over all ancestors.
using PadeTaylor, Printf, Serialization, LinearAlgebra, Statistics
(; wd, r0, r5) = deserialize(joinpath(@__DIR__, "walks.jls"))
relv(a, b) = maximum(abs(a[c]-b[c])/max(1,abs(b[c])) for c in 1:4)
refy(k) = (a = r0.grid_y[k]; all(isfinite, a) ? a : r5.grid_y[k])
function sig(w, k)
    my, mh = 0.0, Inf; j = k
    while j > 0; my = max(my, norm(w.visited_y[j])); mh = min(mh, real(w.visited_h[j])); j = w.visited_parent[j]; end
    my, mh
end
far = [k for k in eachindex(wd.visited_z) if abs(wd.visited_z[k]) ≥ 5.5 && all(isfinite, refy(k))]
rel = [relv(refy(k), wd.visited_y[k]) for k in far]
for (lab, sel) in (("bad (>1e-4)", rel .> 1e-4), ("good (≤1e-4)", rel .≤ 1e-4))
    s = [sig(wd, k) for k in far[sel]]
    isempty(s) && continue
    @printf("default far %-13s n=%3d  chain max‖y‖: min %.2e median %.2e max %.2e | chain min h: min %.4f median %.4f max %.4f\n",
        lab, count(sel), minimum(first.(s)), median(first.(s)), maximum(first.(s)), minimum(last.(s)), median(last.(s)), maximum(last.(s)))
end
s0 = [sig(r0, k) for k in eachindex(r0.visited_z) if abs(r0.visited_z[k]) ≥ 5.5]
@printf("seed0 far all            n=%3d  chain max‖y‖: min %.2e median %.2e max %.2e | chain min h: min %.4f median %.4f max %.4f\n",
    length(s0), minimum(first.(s0)), median(first.(s0)), maximum(first.(s0)), minimum(last.(s0)), median(last.(s0)), maximum(last.(s0)))
# the near-pole node(s) shared by the bad chains
cnt = Dict{Int,Int}()
for k in far[rel .> 1e-4]; j = k; while j > 0; norm(wd.visited_y[j]) > 1e5 && (cnt[j] = get(cnt, j, 0) + 1); j = wd.visited_parent[j]; end; end
for (j, c) in sort(collect(cnt), by = last, rev = true)[1:min(5, end)]
    @printf("  node %d z=%s ‖y‖=%.2e h=%.4f on %d bad chains\n", j, string(round(wd.visited_z[j], digits=4)), norm(wd.visited_y[j]), real(wd.visited_h[j]), c)
end
