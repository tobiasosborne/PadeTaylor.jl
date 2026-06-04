# pole_crossing.jl — close the last ADR-0028 open item: a selector that works
# when a vector step CROSSES a pole (where the held-out-point discriminator is
# invalid).  Scores four discriminators against ground truth in BOTH regimes:
#
#   point — held-out-POINT surrogate truth (the 8/8 winner on pole-free steps;
#           expected to FAIL across a pole — that's the gap we're closing).
#   D1    — relative ODE defect ‖ỹ'−h·f‖/‖h·f‖ (research-backed UNIFORM candidate).
#   D2    — Q²-weighted defect (rigorous pole-crossing form).
#   D3    — PRODUCTION PAIRING: held-out-point where no in-step pole, else D1.
#
# Pole-crossing test case: tan-companion u''=2u+2u³, exact oracle (tan z, sec²z),
# double pole at π/2 — a genuine d=2 shared-pole system with a closed form.
# Run:  julia --project=. external/probes/adr0028-dual-construction-audition/pole_crossing.jl

include(joinpath(@__DIR__, "cells.jl"))
include(joinpath(@__DIR__, "battery.jl"))
using Printf
using PadeTaylor: vector_taylor_coefficients, noumi_yamada_rhs

# --- discriminator wrappers (return :A or :B; A=diagonal, B=B_grow) ----------
disc_point(cA, cBg, jets) = held_out_point(cA, jets, 0.9) ≤ held_out_point(cBg, jets, 0.9) ? :A : :B
disc_d1(cA, cBg, rhs, z0, h) = relative_defect(cA, rhs, z0, h) ≤ relative_defect(cBg, rhs, z0, h) ? :A : :B
disc_d2(cA, cBg, rhs, z0, h) = q_defect(cA, rhs, z0, h) ≤ q_defect(cBg, rhs, z0, h) ? :A : :B
function disc_d3(cA, cBg, jets, rhs, z0, h)
    pm = min(nearest_pole_modulus(cA.denominator), nearest_pole_modulus(cBg.denominator))
    pm < 0.95 ? disc_d1(cA, cBg, rhs, z0, h) : disc_point(cA, cBg, jets)
end

const SC = Dict(:point => [0, 0], :D1 => [0, 0], :D2 => [0, 0], :D3 => [0, 0])
const SCREG = Dict{Tuple{Symbol,Symbol},Vector{Int}}()
reg!(d, reg, ok) = (SC[d][1] += ok; SC[d][2] += 1;
                    v = get!(SCREG, (d, reg), [0, 0]); v[1] += ok; v[2] += 1)

function judge(name, regime, eA, eBg, picks)
    truewin = eA ≤ eBg ? :A : :B
    tie = abs(log10(max(eA, 1e-300)) - log10(max(eBg, 1e-300))) < 0.5
    ok(p) = tie || p == truewin
    line = @sprintf("  %-22s %-9s eA=%.1e eB=%.1e true=%s |", name, String(regime), eA, eBg, String(truewin))
    for d in (:point, :D1, :D2, :D3)
        o = ok(picks[d]); reg!(d, regime, o)
        line *= @sprintf(" %s=%s%s", String(d), String(picks[d]), o ? "✓" : "✗")
    end
    println(line)
end

ev1(c) = [sum(c.numerators[i]) / sum(c.denominator) for i = 1:length(c.numerators)]

println("="^118)
println("POLE-CROSSING DISCRIMINATOR AUDITION — selectors vs ground truth (cell A vs cell B_grow), BOTH regimes")
println("="^118)

# --- regime 1: pole-FREE (the existing battery; held-out-point is valid) -----
println("\n[regime: POLE-FREE]  (held-out-point should work here)")
for s in build_battery()
    jets = rescaled_jets(s)
    cA = cell_A(jets, s.m); cBg = cell_B_grow(jets, s.m)
    eA = step_error(s, cA); eBg = step_error(s, cBg)
    picks = Dict(:point => disc_point(cA, cBg, jets),
                 :D1 => disc_d1(cA, cBg, s.rhs, s.z0, s.h),
                 :D2 => disc_d2(cA, cBg, s.rhs, s.z0, s.h),
                 :D3 => disc_d3(cA, cBg, jets, s.rhs, s.z0, s.h))
    judge(s.name, :polefree, eA, eBg, picks)
end

# --- regime 2: pole-CROSSING (tan-companion; held-out-point INVALID) ---------
println("\n[regime: POLE-CROSSING]  tan-companion u''=2u+2u³, pole at π/2≈1.571 (held-out-point should FAIL)")
tan_rhs = (z, y) -> [y[2], 2y[1] + 2y[1]^3]
for h in (1.8, 2.0, 2.3, 2.6, 3.0)
    m = 12
    raw = vector_taylor_coefficients(tan_rhs, 0.0 + 0im, ComplexF64[0, 1], 2m + 4)
    jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:2]
    cA = cell_A(jets, m); cBg = cell_B_grow(jets, m)
    tru = ComplexF64[tan(h), 1 + tan(h)^2]
    eA = maximum(abs.(ev1(cA) .- tru)); eBg = maximum(abs.(ev1(cBg) .- tru))
    picks = Dict(:point => disc_point(cA, cBg, jets),
                 :D1 => disc_d1(cA, cBg, tan_rhs, 0.0 + 0im, h),
                 :D2 => disc_d2(cA, cBg, tan_rhs, 0.0 + 0im, h),
                 :D3 => disc_d3(cA, cBg, jets, tan_rhs, 0.0 + 0im, h))
    judge(@sprintf("tan-companion h=%.1f", h), :crossing, eA, eBg, picks)
end

# --- regime 3: BigFloat-256 (the precision-dependent A-win must be caught) ----
println("\n[regime: BigFloat-256]  (the genuine-pole×precision A-win — does the defect catch it?)")
setprecision(BigFloat, 256)
const CBF = Complex{BigFloat}
function bf_judge(name, rhs, y0, h, m, oracle; conservation = false)
    raw = vector_taylor_coefficients(rhs, CBF(0), CBF.(y0), 2m + 4)
    jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:length(y0)]
    cA = cell_A(jets, m); cBg = cell_B_grow(jets, m)
    err(c) = (yv = ev1(c); conservation ? Float64(abs(sum(yv) - h)) :
              Float64(maximum(abs(yv[i] - oracle[i]) for i = 1:length(y0))))
    picks = Dict(:point => disc_point(cA, cBg, jets),
                 :D1 => disc_d1(cA, cBg, rhs, CBF(0), h),
                 :D2 => disc_d2(cA, cBg, rhs, CBF(0), h),
                 :D3 => disc_d3(cA, cBg, jets, rhs, CBF(0), h))
    judge(name, :bigfloat, err(cA), err(cBg), picks)
end
hb = BigFloat("0.1")
bf_judge("calogero-moser-BF", (z, y) -> (dd = y[1]-y[2]; a = 4/dd^3; [y[3], y[4], a, -a]),
         CBF[1, -1, 0, 0], hb, 15, (x1 = sqrt(1+hb^2/2); v1 = hb/(2x1); CBF[x1, -x1, v1, -v1]))
hh = BigFloat("0.7")
bf_judge("harmonic-BF", (z, y) -> [y[2], -y[1]], CBF[1, 0], hh, 15, CBF[cos(hh), -sin(hh)])

# --- regime 4: BigFloat-256 pole-CROSSING (the hard case: does A ever win a
#     crossing at extended precision, and does each rule catch it?) -----------
println("\n[regime: BF-256 POLE-CROSSING]  tan-companion across π/2 at 256-bit, varying degree m")
tan_rhs_bf = (z, y) -> [y[2], 2y[1] + 2y[1]^3]
for h in (BigFloat("1.8"), BigFloat("2.2"))
    for m in (12, 20, 28)
        raw = vector_taylor_coefficients(tan_rhs_bf, CBF(0), CBF[0, 1], 2m + 4)
        jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:2]
        cA = cell_A(jets, m); cBg = cell_B_grow(jets, m)
        tru = CBF[tan(h), 1 + tan(h)^2]
        eA = Float64(maximum(abs.(ev1(cA) .- tru))); eBg = Float64(maximum(abs.(ev1(cBg) .- tru)))
        picks = Dict(:point => disc_point(cA, cBg, jets),
                     :D1 => disc_d1(cA, cBg, tan_rhs_bf, CBF(0), h),
                     :D2 => disc_d2(cA, cBg, tan_rhs_bf, CBF(0), h),
                     :D3 => disc_d3(cA, cBg, jets, tan_rhs_bf, CBF(0), h))
        judge(@sprintf("tan-BF h=%.1f m=%d", Float64(h), m), :bfcross, eA, eBg, picks)
    end
end

# --- scoreboard --------------------------------------------------------------
println("\n" * "="^118)
println("SCOREBOARD")
println("="^118)
@printf("  %-8s %-11s %-13s %-11s %-13s %-8s\n", "rule", "polefree", "crossing", "bigfloat", "bf-crossing", "TOTAL")
for d in (:point, :D1, :D2, :D3)
    pf = get(SCREG, (d, :polefree), [0, 0]); cr = get(SCREG, (d, :crossing), [0, 0])
    bf = get(SCREG, (d, :bigfloat), [0, 0]); bc = get(SCREG, (d, :bfcross), [0, 0])
    @printf("  %-8s %-11s %-13s %-11s %-13s %d/%d\n", String(d),
        "$(pf[1])/$(pf[2])", "$(cr[1])/$(cr[2])", "$(bf[1])/$(bf[2])", "$(bc[1])/$(bc[2])", SC[d][1], SC[d][2])
end
println("\n  point = held-out-point | D1 = relative defect | D2 = Q²-defect | D3 = point+defect-fallback")
