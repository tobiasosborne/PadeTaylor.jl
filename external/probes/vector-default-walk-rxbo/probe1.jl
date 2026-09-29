# Probe 1 (bead rxbo): locate the FIRST edge of the default h=0.1 walk whose
# child departs from the seed-0 reference field while its parent agrees.
using PadeTaylor, Printf, Random, LinearAlgebra, Serialization
include(joinpath(@__DIR__, "..", "..", "..", "figures", "_kkg_pi2_helpers.jl"))
bvp = kkg_bvp_solve()
f = painleve_hierarchy(:I, 2; t = KKG_T)
y_seed = ComplexF64.(bvp(ComplexF64(KKG_Z_SEED)))
prob = VectorPadeTaylorProblem(f, y_seed, (ComplexF64(KKG_Z_SEED), 8.0+0.0im); order = KKG_PN_ORDER)
targets = kkg_target_wedge()
wd = vector_path_network_solve(prob, targets; order=KKG_PN_ORDER, h=KKG_PN_H)
# references: seed 0 h=0.1 and default h=0.05 evaluated at the default's nodes
r0 = vector_path_network_solve(prob, targets; order=KKG_PN_ORDER, h=KKG_PN_H, rng=MersenneTwister(0), fine_grid=wd.visited_z)
r5 = vector_path_network_solve(prob, targets; order=KKG_PN_ORDER, h=0.05, fine_grid=wd.visited_z)
serialize(joinpath(@__DIR__, "walks.jls"), (wd=wd, r0=r0, r5=r5))
n = length(wd.visited_z)
relk(ref, k) = maximum(abs(ref.grid_y[k][c]-wd.visited_y[k][c])/max(1,abs(wd.visited_y[k][c])) for c in 1:4)
rel = [begin a=relk(r0,k); isfinite(a) ? a : relk(r5,k) end for k in 1:n]
cov = isfinite.(rel)
bad = [k for k in 1:n if cov[k] && rel[k] > 1e-4]
@printf("default nodes %d, covered %d, bad %d\n", n, count(cov), length(bad))
# first departing edges: child bad, parent covered & good
first_edges = [k for k in bad if (p = wd.visited_parent[k]; p>0 && cov[p] && rel[p] ≤ 1e-6)]
@printf("clean->bad edges: %d\n", length(first_edges))
# the earliest (by index) departure along each bad node's chain
roots = Set{Int}()
for k in bad
    j = k; last_bad = k
    while j > 0
        if cov[j] && rel[j] > 1e-4; last_bad = j; end
        cov[j] && rel[j] ≤ 1e-6 && break
        j = wd.visited_parent[j]
    end
    push!(roots, last_bad)
end
println("earliest bad ancestors: ", sort(collect(roots)))
for k in sort(collect(roots))
    p = wd.visited_parent[k]
    h = real(wd.visited_h[p])
    den = wd.visited_denominator[p]
    rts = length(den) > 1 ? PadeTaylor.VectorPoleField.roots(PadeTaylor.VectorPoleField.Polynomial(den)) : ComplexF64[]
    zp = [wd.visited_z[p] + h*t for t in rts if abs(t) < 4]
    @printf("\nnode %d z=%s |z|=%.3f rel=%.3e ; parent %d z=%s rel(parent)=%s\n", k, string(round(wd.visited_z[k],digits=4)), abs(wd.visited_z[k]), rel[k], p, string(round(wd.visited_z[p],digits=4)), cov[p] ? @sprintf("%.2e",rel[p]) : "uncov")
    @printf("  step = %s |step|=%.4f h_parent=%.4f h_child=%.4f\n", string(round(wd.visited_z[k]-wd.visited_z[p],digits=4)), abs(wd.visited_z[k]-wd.visited_z[p]), h, real(wd.visited_h[k]))
    seg(a,b,x) = (t = clamp(real((x-a)*conj(b-a))/abs2(b-a),0,1); abs(a+t*(b-a)-x))
    for q in sort(zp, by=x->seg(wd.visited_z[p],wd.visited_z[k],x))[1:min(3,end)]
        @printf("  parent-Q pole z=%s  dist to segment=%.4f  |t|=%.3f\n", string(round(q,digits=4)), seg(wd.visited_z[p],wd.visited_z[k],q), abs((q-wd.visited_z[p])/h))
    end
end
