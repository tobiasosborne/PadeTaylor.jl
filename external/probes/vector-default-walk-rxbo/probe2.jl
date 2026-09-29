# Probe 2 (bead rxbo): error growth along the parent chain of the worst
# default-walk node, with per-edge step length, adaptive h, and nearest
# parent-Q pole (z-plane) — from the walks serialised by probe1.jl.
using PadeTaylor, Printf, Serialization, LinearAlgebra
const VPF = PadeTaylor.VectorPoleField
(; wd, r0, r5) = deserialize(joinpath(@__DIR__, "walks.jls"))
n = length(wd.visited_z)
relk(ref, k) = maximum(abs(ref.grid_y[k][c]-wd.visited_y[k][c])/max(1,abs(wd.visited_y[k][c])) for c in 1:4)
rel = [begin a=relk(r0,k); isfinite(a) ? a : relk(r5,k) end for k in 1:n]
function poles(k)
    den = wd.visited_denominator[k]; h = real(wd.visited_h[k])
    length(den) ≤ 1 && return ComplexF64[]
    [wd.visited_z[k] + h*t for t in VPF.roots(VPF.Polynomial(den))]
end
kw = argmax(k -> isfinite(rel[k]) ? rel[k] : -1.0, 1:n)
chain = let c = Int[], j = kw; while j > 0; pushfirst!(c, j); j = wd.visited_parent[j]; end; c end
@printf("worst node %d z=%s rel=%.3e ; chain length %d\n", kw, string(wd.visited_z[kw]), rel[kw], length(chain))
for (i,k) in enumerate(chain)
    p = wd.visited_parent[k]
    pl = poles(k); dmin = isempty(pl) ? Inf : minimum(abs.(pl .- wd.visited_z[k]))
    st = p > 0 ? abs(wd.visited_z[k]-wd.visited_z[p]) : 0.0
    (i ≤ 3 || rel[k] > 1e-9 || !isfinite(rel[k])) &&
    @printf("%4d node %4d z=%8.4f%+8.4fi |z|=%.3f  rel=%9.2e  step=%.4f h=%.4f nearestQpole=%.4f ‖y‖=%.3e\n",
        i, k, real(wd.visited_z[k]), imag(wd.visited_z[k]), abs(wd.visited_z[k]), rel[k], st, real(wd.visited_h[k]), dmin, norm(wd.visited_y[k]))
end
