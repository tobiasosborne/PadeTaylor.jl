# =============================================================================
# d2_mechanism.jl — why does the corrected (m,m) window (+2) give Q(0)≈0 on
# the ENTIRE harmonic d=2 system [cos, -sin] while the buggy (m-1,m) (+1) gave
# Q(0)≈1?  And does degree-reduction recover a usable shared Q?
#
# This informs the C1<->C2 coupling decision (worklog 065 / beads d3a+3p9c).
# Read-only: replicates the block build + raw SVD, no src edits.
# =============================================================================
using LinearAlgebra, Printf

# harmonic Taylor jets about z=0, y=[1,0]:  y1=cos, y2=-sin
cosjet(N) = Float64[ iseven(k) ? Float64((-1)^(k ÷ 2) / factorial(big(k))) : 0.0 for k = 0:N ]
nsinjet(N) = Float64[ isodd(k) ? Float64(-(-1)^((k - 1) ÷ 2) / factorial(big(k))) : 0.0 for k = 0:N ]
rescale(c, h) = [h^(k - 1) * c[k] for k = 1:length(c)]

# block builders for the two windows
block(c, m, shift) = [ (idx = m + r - cc + shift; 1 ≤ idx ≤ length(c) ? c[idx] : 0.0)
                       for r = 1:m, cc = 1:(m + 1) ]   # shift=1 -> buggy (m-1,m); shift=2 -> (m,m)

function probe(m, h)
    jets = [rescale(cosjet(2m + 1), h), rescale(nsinjet(2m + 1), h)]
    for (name, shift) in (("(m-1,m) +1 buggy", 1), ("(m,m) +2 fixed", 2))
        A = vcat([block(c, m, shift) for c in jets]...)   # 2m x (m+1)
        S = svdvals(A)
        F = svd(A; full = true)
        b = F.Vt[end, :]
        cnorm = norm(vcat(jets...))
        τ = 1e-14 * cnorm
        ρ = count(>(τ), S)
        @printf("  m=%2d %-18s: σ_min=%.2e σ_min/σ_max=%.2e  ρ=%d/%d  b[1]=%+.3e  |b[1]|/‖b‖=%.2e\n",
                m, name, S[end], S[end] / S[1], ρ, m + 1, b[1], abs(b[1]) / norm(b))
    end
end

println("Harmonic d=2 [cos,-sin], h=0.7 (VS.1.2 uses m=order÷2=15).")
println("b = smallest-right-singular-vector; b[1]=Q(0) (pre-QR-reweight, pre-normalise).")
println("ρ = rank above τ; ρ=m ⇒ isolated null space, ρ=m+1 ⇒ full column rank.\n")
for m in (2, 3, 5, 8, 12, 15)
    probe(m, 0.7)
end
