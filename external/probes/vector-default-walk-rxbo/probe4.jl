# Probe 4 (bead rxbo): restart-bisection along the default chain of node 493.
# Replace the state at chain position i by the REFERENCE state, re-step the
# SAME z-sequence to node 493, compare with the reference there.  Also push a
# 1e-10 relative perturbation of the reference state from position i to 493
# (the amplification factor of the IVP along the chain).
using PadeTaylor, Printf, Serialization, LinearAlgebra, Random
(; wd, r0, r5) = deserialize(joinpath(@__DIR__, "walks.jls"))
include(joinpath(@__DIR__, "..", "..", "..", "figures", "_kkg_pi2_helpers.jl"))
f = painleve_hierarchy(:I, 2; t = KKG_T)
refy(k) = (a = r0.grid_y[k]; all(isfinite, a) ? a : r5.grid_y[k])
relv(a, b) = maximum(abs(a[c]-b[c])/max(1,abs(b[c])) for c in 1:4)
kw = 493
chain = let c = Int[], j = kw; while j > 0; pushfirst!(c, j); j = wd.visited_parent[j]; end; c end
function run_from(i, y0)
    st = PadeTaylor.VectorPadeStepperState{ComplexF64}(wd.visited_z[chain[i]], Vector{ComplexF64}(y0))
    for j in i+1:length(chain)
        PadeTaylor.vector_pade_step_with_pade!(st, f, KKG_PN_ORDER, wd.visited_z[chain[j]] - st.z)
    end
    st.y
end
yend = refy(kw)
@printf("stored default end vs ref: %.3e\n", relv(wd.visited_y[kw], yend))
@printf("replay from node 1 (seed state) vs stored default end: %.3e\n", relv(run_from(1, wd.visited_y[1]), wd.visited_y[kw]))
rng = MersenneTwister(1)
for i in (1, 30, 60, 100, 130, 150, 170, 200, 230, 250, 280, 300, 320)
    y0 = refy(chain[i]); all(isfinite, y0) || (println("pos $i: ref uncovered"); continue)
    e_ref = relv(run_from(i, y0), yend)
    δ = 1e-10 .* (randn(rng, ComplexF64, 4)) .* max.(1, abs.(y0))
    amp = relv(run_from(i, y0 .+ δ), run_from(i, y0)) / 1e-10
    @printf("restart pos %3d (node %4d, |z|=%.3f): end err vs ref = %.3e   amplification(1e-10 kick) = %.3e\n",
            i, chain[i], abs(wd.visited_z[chain[i]]), e_ref, amp)
end
