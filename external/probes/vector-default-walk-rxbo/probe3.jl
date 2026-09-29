# Probe 3 (bead rxbo): (a) per-edge LOCAL error of the default chain, stepping
# from the REFERENCE state at the parent; (b) tree depth / path length of
# far-wedge nodes, default vs seed 0.
using PadeTaylor, Printf, Serialization, LinearAlgebra, Statistics
(; wd, r0, r5) = deserialize(joinpath(@__DIR__, "walks.jls"))
include(joinpath(@__DIR__, "..", "..", "..", "figures", "_kkg_pi2_helpers.jl"))
f = painleve_hierarchy(:I, 2; t = KKG_T)
n = length(wd.visited_z)
refy(k) = (a = r0.grid_y[k]; all(isfinite, a) ? a : r5.grid_y[k])
relv(a, b) = maximum(abs(a[c]-b[c])/max(1,abs(b[c])) for c in 1:4)
rel = [relv(refy(k), wd.visited_y[k]) for k in 1:n]
kw = 493
chain = let c = Int[], j = kw; while j > 0; pushfirst!(c, j); j = wd.visited_parent[j]; end; c end
println("edge-local error along chain of node $kw (reference-state restart):")
for i in 2:length(chain)
    p, k = chain[i-1], chain[i]
    yp = refy(p); yk = refy(k)
    (all(isfinite, yp) && all(isfinite, yk)) || continue
    st = PadeTaylor.VectorPadeStepperState{ComplexF64}(wd.visited_z[p], Vector{ComplexF64}(yp))
    PadeTaylor.vector_pade_step_with_pade!(st, f, KKG_PN_ORDER, wd.visited_z[k]-wd.visited_z[p])
    loc = relv(st.y, yk)
    # also the default's own step from its own parent state vs its stored child (sanity: 0)
    (i % 10 == 0 || loc > 1e-6 || i > 240) && @printf("%4d %4d->%4d |z|=%.3f step=%.4f  local=%.2e  accumulated=%.2e\n", i, p, k, abs(wd.visited_z[k]), abs(wd.visited_z[k]-wd.visited_z[p]), loc, rel[k])
end
depth(w, k) = (d = 0; j = k; while j > 0; d += 1; j = w.visited_parent[j]; end; d)
plen(w, k) = (s = 0.0; j = k; while w.visited_parent[j] > 0; s += abs(w.visited_z[j]-w.visited_z[w.visited_parent[j]]); j = w.visited_parent[j]; end; s)
for (tag, w) in (("default", wd), ("seed0", r0))
    far = [k for k in eachindex(w.visited_z) if abs(w.visited_z[k]) ≥ 5.5]
    ds = [depth(w, k) for k in far]; ls = [plen(w, k) for k in far]
    @printf("%-8s far nodes %d  depth median %.0f max %d   path length median %.2f max %.2f  (|z| range for far: straight-line ≈ %.1f)\n",
        tag, length(far), median(ds), maximum(ds), median(ls), maximum(ls), 3 + 6.0)
end
