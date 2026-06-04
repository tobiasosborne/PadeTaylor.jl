# addendum.jl — skeptical stress-tests of the audition's "square cell B wins
# everywhere" finding, in the two regimes most likely to favour cell A:
#   (1) BigFloat-256 extended precision (the ADR's stated concern);
#   (2) multi-pole crossings (where A's higher matching order should help).
# Run from repo root:
#   julia --project=. external/probes/adr0028-dual-construction-audition/addendum.jl

include(joinpath(@__DIR__, "cells.jl"))
using Printf
using PadeTaylor: vector_taylor_coefficients

setprecision(BigFloat, 256)
const CBF = Complex{BigFloat}

eval1(num, Q) = sum(num) / sum(Q)

function bf_step(name, kind, rhs, y0, h, m, oracle)
    order = 2m + 4
    raw = vector_taylor_coefficients(rhs, CBF(0), CBF.(y0), order)
    jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:length(y0)]
    cA = cell_A(jets, m); cBf = cell_B_floor(jets, m)
    yA = [eval1(cA.numerators[i], cA.denominator) for i = 1:length(y0)]
    yB = [eval1(cBf.numerators[i], cBf.denominator) for i = 1:length(y0)]
    eA = maximum(abs(yA[i] - oracle[i]) for i = 1:length(y0))
    eB = maximum(abs(yB[i] - oracle[i]) for i = 1:length(y0))
    @printf("  %-16s %-12s | err_A=%.3e  err_Bfloor=%.3e | %s\n",
            name, kind, Float64(eA), Float64(eB),
            eB < eA / 3 ? "B wins ×$(round(Float64(eA/eB),sigdigits=2))" :
            eA < eB / 3 ? "A WINS ×$(round(Float64(eB/eA),sigdigits=2))" : "≈ tie")
end

println("="^88)
println("ADDENDUM 1 — BigFloat-256: does square cell B still dominate at extended precision?")
println("="^88)
hbig = BigFloat("0.7")
bf_step("harmonic", "entire", (z, y) -> [y[2], -y[1]], CBF[1, 0], hbig, 15,
        CBF[cos(hbig), -sin(hbig)])
hb2 = BigFloat("0.1")
bf_step("calogero-moser", "meromorphic",
        (z, y) -> (dd = y[1] - y[2]; a = 4 / dd^3; [y[3], y[4], a, -a]),
        CBF[1, -1, 0, 0], hb2, 15,
        (x1 = sqrt(1 + hb2^2 / 2); v1 = hb2 / (2x1); CBF[x1, -x1, v1, -v1]))
# ℘ companion at BigFloat: oracle only pinned to F64 digits, so compare A vs B
# by their MUTUAL agreement to the F64 oracle (both should land near 4.0044…)
hb3 = BigFloat("0.5")
bf_step("weierstrass-℘", "meromorphic", (z, y) -> [y[2], 6 * y[1]^2],
        CBF[big"1.071822516416917", big"1.710337353176786"], hb3, 15,
        CBF[big"4.0044646690030875", big"15.964278048239492"])

println()
println("="^88)
println("ADDENDUM 2 — multi-pole crossings (tan, poles at π/2, 3π/2): A=(m,m) vs B=(m-1,m)")
println("  Does the higher-order diagonal capture MULTIPLE in-step poles better?")
println("="^88)
function tan_jet(N)
    c = zeros(Float64, N + 1)
    for k = 0:N-1
        conv = sum(c[j+1] * c[k-j+1] for j = 0:k)
        c[k+2] = ((k == 0 ? 1.0 : 0.0) + conv) / (k + 1)
    end
    ComplexF64.(c)
end
for h in (3.0, 5.0, 6.0)          # 3.0: 1 pole; 5.0: 2 poles; 6.0: 2 poles, near 3π/2
    npole = count(p -> 0 < p < h, (π/2, 3π/2))
    for m in (12, 18, 24)
        jet = tan_jet(2m + 4)
        jr = [h^(k - 1) * jet[k] for k = 1:length(jet)]
        pA = robust_pade(jr, m, m; method = :svd)
        pB = robust_pade(jr, m - 1, m; method = :svd)
        ev(p) = sum(p.a) / sum(p.b)
        tru = tan(h)
        eA = abs(ev(pA) - tru); eB = abs(ev(pB) - tru)
        @printf("  h=%.1f (%d poles in-step) m=%2d | err_A(m,m)=%.2e  err_B(m-1,m)=%.2e | %s\n",
                h, npole, m, eA, eB,
                eA < eB / 3 ? "A WINS" : eB < eA / 3 ? "B wins" : "≈ tie")
    end
end
