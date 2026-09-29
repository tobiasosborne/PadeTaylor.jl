# Probe 6 (bead rxbo): is the amplification across the near-pole passage
# (default chain of node 493, positions 161-201) a property of the ODE (then
# it is step-size independent) or of the discrete Padé step (then it changes
# when each edge is split into s sub-steps)?  Linear-regime kick ε = 1e-12.
using PadeTaylor, Printf, Serialization, LinearAlgebra, Random
(; wd, r0, r5) = deserialize(joinpath(@__DIR__, "walks.jls"))
include(joinpath(@__DIR__, "..", "..", "..", "figures", "_kkg_pi2_helpers.jl"))
f = painleve_hierarchy(:I, 2; t = KKG_T)
relv(a, b) = maximum(abs(a[c]-b[c])/max(1,abs(b[c])) for c in 1:4)
chain = let c = Int[], j = 493; while j > 0; pushfirst!(c, j); j = wd.visited_parent[j]; end; c end
function run_seg(i1, i2, y0, s)
    st = PadeTaylor.VectorPadeStepperState{ComplexF64}(wd.visited_z[chain[i1]], Vector{ComplexF64}(y0))
    for j in i1+1:i2, _ in 1:s
        dz = (wd.visited_z[chain[j]] - wd.visited_z[chain[j-1]]) / s
        PadeTaylor.vector_pade_step_with_pade!(st, f, KKG_PN_ORDER, dz)
    end
    st.y
end
for (i1, i2) in ((141, 161), (161, 181), (171, 191), (181, 201), (201, 241), (241, 326))
    y0 = wd.visited_y[chain[i1]]
    for s in (1, 2, 4)
        base = run_seg(i1, i2, y0, s)
        a = maximum(1:3) do r
            δ = 1e-12 .* randn(MersenneTwister(r), ComplexF64, 4) .* max.(1, abs.(y0))
            relv(run_seg(i1, i2, y0 .+ δ, s), base) / 1e-12
        end
        @printf("seg pos %3d→%3d (|z| %.3f→%.3f) substeps=%d  amplification=%.3e  end‖y‖=%.2e  vs stored=%.2e\n",
            i1, i2, abs(wd.visited_z[chain[i1]]), abs(wd.visited_z[chain[i2]]), s, a, norm(base),
            relv(base, wd.visited_y[chain[i2]]))
    end
end
