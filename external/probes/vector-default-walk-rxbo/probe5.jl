# Probe 5 (bead rxbo): (a) chain positions 140-205 of default node 493 with
# per-position amplification to the end; (b) the same kick-amplification
# along the SEED-0 walk's own chain to its node nearest z_493.
using PadeTaylor, Printf, Serialization, LinearAlgebra, Random
(; wd, r0, r5) = deserialize(joinpath(@__DIR__, "walks.jls"))
include(joinpath(@__DIR__, "..", "..", "..", "figures", "_kkg_pi2_helpers.jl"))
f = painleve_hierarchy(:I, 2; t = KKG_T)
const VPF = PadeTaylor.VectorPoleField
relv(a, b) = maximum(abs(a[c]-b[c])/max(1,abs(b[c])) for c in 1:4)
chainof(w, k) = let c = Int[], j = k; while j > 0; pushfirst!(c, j); j = w.visited_parent[j]; end; c end
function run_from(w, chain, i, y0)
    st = PadeTaylor.VectorPadeStepperState{ComplexF64}(w.visited_z[chain[i]], Vector{ComplexF64}(y0))
    for j in i+1:length(chain)
        PadeTaylor.vector_pade_step_with_pade!(st, f, KKG_PN_ORDER, w.visited_z[chain[j]] - st.z)
    end
    st.y
end
npole(w, k) = (den = w.visited_denominator[k]; length(den) ≤ 1 ? Inf :
    minimum(abs.(real(w.visited_h[k]) .* VPF.roots(VPF.Polynomial(den)))))
function amp(w, chain, i; ε = 1e-13)
    rng = MersenneTwister(i); y0 = w.visited_y[chain[i]]
    δ = ε .* randn(rng, ComplexF64, 4) .* max.(1, abs.(y0))
    relv(run_from(w, chain, i, y0 .+ δ), run_from(w, chain, i, y0)) / ε
end
for (tag, w, kend) in (("default", wd, 493),
                       ("seed0", r0, argmin(k -> abs(r0.visited_z[k] - wd.visited_z[493]), eachindex(r0.visited_z))))
    ch = chainof(w, kend)
    @printf("\n== %s: end node %d z=%s, chain length %d\n", tag, kend, string(round(w.visited_z[kend], digits=4)), length(ch))
    for i in vcat(1:10:length(ch)-1)
        k = ch[i]
        @printf("  pos %3d node %4d z=%8.4f%+8.4fi |z|=%.3f h=%.4f nearestQpole=%.4f ‖y‖=%.2e  amp(1e-13 kick→end)=%.3e\n",
            i, k, real(w.visited_z[k]), imag(w.visited_z[k]), abs(w.visited_z[k]), real(w.visited_h[k]), npole(w, k), norm(w.visited_y[k]), amp(w, ch, i))
    end
end
