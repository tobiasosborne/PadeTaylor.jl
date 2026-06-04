# addendum2.jl — which SQUARE cell-B realization is right when d∤m?
#   B_floor = ⌊m/d⌋·d   (wide square, LOWER degree)
#   B_grow  = ⌈m/d⌉·d   (wide square, HIGHER degree ≥ m)   ← candidate best-of-both
# Tested vs A on F64 battery and on the BigFloat-256 case where A beat B_floor.
# Run:  julia --project=. external/probes/adr0028-dual-construction-audition/addendum2.jl

include(joinpath(@__DIR__, "cells.jl"))
include(joinpath(@__DIR__, "battery.jl"))
using Printf
using PadeTaylor: vector_taylor_coefficients

eval1(num, Q) = sum(num) / sum(Q)
function cellerr(s, cell)
    yv = [eval1(cell.numerators[i], cell.denominator) for i = 1:s.d]
    s.oracle === :conservation ? abs(sum(yv) - (s.z0 + s.h)) :
        (tru = s.oracle(); maximum(abs(yv[i] - tru[i]) for i = 1:s.d))
end

println("="^96)
println("F64 battery — A vs B_floor (⌊m/d⌋) vs B_grow (⌈m/d⌉, wide square higher degree)")
println("="^96)
@printf("%-16s %-6s %2s %3s | %-10s %-10s %-10s | degrees A/Bf/Bg\n",
        "system", "kind", "d", "m", "err_A", "err_Bfloor", "err_Bgrow")
println("-"^96)
for s in build_battery()
    jets = rescaled_jets(s)
    cA = cell_A(jets, s.m); cBf = cell_B_floor(jets, s.m); cBg = cell_B_grow(jets, s.m)
    @printf("%-16s %-6s %2d %3d | %-10.2e %-10.2e %-10.2e | %d/%d/%d\n",
            s.name, String(s.kind), s.d, s.m, cellerr(s, cA), cellerr(s, cBf), cellerr(s, cBg),
            length(cA.denominator) - 1, length(cBf.denominator) - 1, length(cBg.denominator) - 1)
end

println()
println("="^96)
println("BigFloat-256 — the case A beat B_floor (CM ×3.5e7): does B_grow close the gap?")
println("="^96)
setprecision(BigFloat, 256)
const CBF = Complex{BigFloat}
function bf3(name, rhs, y0, h, m, oracle)
    raw = vector_taylor_coefficients(rhs, CBF(0), CBF.(y0), 2m + 4)
    jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:length(y0)]
    cA = cell_A(jets, m); cBf = cell_B_floor(jets, m); cBg = cell_B_grow(jets, m)
    e(c) = (yv = [eval1(c.numerators[i], c.denominator) for i = 1:length(y0)];
            Float64(maximum(abs(yv[i] - oracle[i]) for i = 1:length(y0))))
    @printf("  %-16s | err_A=%.2e (ν=%d)  err_Bfloor=%.2e (ν=%d)  err_Bgrow=%.2e (ν=%d)\n",
            name, e(cA), length(cA.denominator) - 1, e(cBf), length(cBf.denominator) - 1,
            e(cBg), length(cBg.denominator) - 1)
end
hb = BigFloat("0.1")
bf3("calogero-moser", (z, y) -> (dd = y[1] - y[2]; a = 4 / dd^3; [y[3], y[4], a, -a]),
    CBF[1, -1, 0, 0], hb, 15, (x1 = sqrt(1 + hb^2 / 2); v1 = hb / (2x1); CBF[x1, -x1, v1, -v1]))
hh = BigFloat("0.7")
bf3("harmonic", (z, y) -> [y[2], -y[1]], CBF[1, 0], hh, 15, CBF[cos(hh), -sin(hh)])
nyα = CBF[0.30, 0.10, 0.25, 0.15, 0.20]
# NY conservation oracle handled separately (Σf=t): print drift for each
let
    rhs = noumi_yamada_rhs(2; α = nyα)
    y0 = CBF[0.7, -0.3, 0.5, -0.55, -0.35]; h = BigFloat("0.3"); m = 12
    raw = vector_taylor_coefficients(rhs, CBF(0), y0, 2m + 4)
    jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:5]
    drift(c) = Float64(abs(sum(eval1(c.numerators[i], c.denominator) for i = 1:5) - h))
    cA = cell_A(jets, m); cBf = cell_B_floor(jets, m); cBg = cell_B_grow(jets, m)
    @printf("  %-16s | drift_A=%.2e (ν=%d)  drift_Bfloor=%.2e (ν=%d)  drift_Bgrow=%.2e (ν=%d)\n",
            "noumi-yamada-A₄", drift(cA), length(cA.denominator) - 1, drift(cBf),
            length(cBf.denominator) - 1, drift(cBg), length(cBg.denominator) - 1)
end
