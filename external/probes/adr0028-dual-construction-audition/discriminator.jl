# discriminator.jl — prototype the REDESIGNED selector for ADR-0028.
#
# The audition (FINDINGS.md item 1) showed the ADR's (R,g,K) Pareto mis-ranks
# 2/5 — but the failures were caused by the conditioning-g fallback + Rfloor
# deferral, NOT by the held-out residual R itself.  Here we score candidate
# selection rules against GROUND TRUTH (the oracle endpoint error) across the
# F64 battery + the BigFloat-256 regime where cell A genuinely wins (CM):
#
#   rule_pureR   — argmin held-out residual R (multi-coeff), tiebreak K↑ then A;
#                  NO conditioning g, NO Rfloor.
#   rule_ADR     — the ADR-as-written (R→g→K→A, ε-bands).
#   rule_point   — argmin held-out-POINT surrogate-truth error at tstar=0.9.
#   rule_hybrid  — held-out-point primary, pureR as tiebreak.
#
# Cell pair = A (GGT diagonal) vs B_grow (the corrected wide-square cell B).
# Run:  julia --project=. external/probes/adr0028-dual-construction-audition/discriminator.jl

include(joinpath(@__DIR__, "cells.jl"))
include(joinpath(@__DIR__, "battery.jl"))
using Printf
using PadeTaylor: vector_taylor_coefficients, noumi_yamada_rhs

const TSTAR = 0.9

function rules(cA, cBg, jets)
    Kc = min(cA.K, cBg.K)
    RA = held_out_multi(cA, jets, Kc; span = 6); RB = held_out_multi(cBg, jets, Kc; span = 6)
    PA = held_out_point(cA, jets, TSTAR);        PB = held_out_point(cBg, jets, TSTAR)
    pureR = RA ≤ RB ? :A : :B
    point = PA ≤ PB ? :A : :B
    # hybrid: point decides unless within 1 order, then pureR
    hybrid = abs(log10(max(PA, 1e-300)) - log10(max(PB, 1e-300))) > 1 ? point : pureR
    adr = select([cA, cBg], jets)[1].label == :diagonal ? :A : :B
    return (; pureR, point, hybrid, adr, RA, RB, PA, PB)
end

const SCORE = Dict(:pureR => 0, :point => 0, :hybrid => 0, :adr => 0)
const NTOT = Ref(0)

function judge(name, kind, eA, eBg, r)
    truewin = eA ≤ eBg ? :A : :B
    tie = abs(log10(max(eA, 1e-300)) - log10(max(eBg, 1e-300))) < 0.5
    ok(c) = tie || c == truewin
    NTOT[] += 1
    for k in (:pureR, :point, :hybrid, :adr)
        ok(getfield(r, k)) && (SCORE[k] += 1)
    end
    mark(c) = ok(c) ? "✓" : "✗"
    @printf("  %-16s %-11s eA=%.1e eB=%.1e true=%s | pureR=%s%s point=%s%s hyb=%s%s ADR=%s%s\n",
        name, kind, eA, eBg, String(truewin),
        String(r.pureR), mark(r.pureR), String(r.point), mark(r.point),
        String(r.hybrid), mark(r.hybrid), String(r.adr), mark(r.adr))
end

eval1(num, Q) = sum(num) / sum(Q)

println("="^110)
println("DISCRIMINATOR PROTOTYPE — selection rules vs ground-truth winner (A vs B_grow)")
println("="^110)
println("F64 battery:")
for s in build_battery()
    jets = rescaled_jets(s)
    cA = cell_A(jets, s.m); cBg = cell_B_grow(jets, s.m)
    eA = step_error(s, cA); eBg = step_error(s, cBg)
    judge(s.name, String(s.kind), eA, eBg, rules(cA, cBg, jets))
end

println("\nBigFloat-256 (the regime where A genuinely wins on genuine poles):")
setprecision(BigFloat, 256)
const CBF = Complex{BigFloat}
function bf_judge(name, kind, rhs, y0, h, m, oracle; conservation = false)
    raw = vector_taylor_coefficients(rhs, CBF(0), CBF.(y0), 2m + 6)
    jets = [[h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:length(y0)]
    cA = cell_A(jets, m); cBg = cell_B_grow(jets, m)
    function err(c)
        yv = [eval1(c.numerators[i], c.denominator) for i = 1:length(y0)]
        conservation ? Float64(abs(sum(yv) - h)) :
            Float64(maximum(abs(yv[i] - oracle[i]) for i = 1:length(y0)))
    end
    judge(name, kind, err(cA), err(cBg), rules(cA, cBg, jets))
end
hb = BigFloat("0.1")
bf_judge("calogero-moser", "mero-BF", (z, y) -> (dd = y[1]-y[2]; a = 4/dd^3; [y[3], y[4], a, -a]),
         CBF[1, -1, 0, 0], hb, 15, (x1 = sqrt(1+hb^2/2); v1 = hb/(2x1); CBF[x1, -x1, v1, -v1]))
hh = BigFloat("0.7")
bf_judge("harmonic", "ent-BF", (z, y) -> [y[2], -y[1]], CBF[1, 0], hh, 15, CBF[cos(hh), -sin(hh)])
bf_judge("noumi-yamada-A₄", "mero-BF", noumi_yamada_rhs(2; α = CBF[0.30, 0.10, 0.25, 0.15, 0.20]),
         CBF[0.7, -0.3, 0.5, -0.55, -0.35], BigFloat("0.3"), 12, nothing; conservation = true)

println("\n" * "="^110)
@printf("SCORES (out of %d):  pureR %d  |  held-out-point %d  |  hybrid %d  |  ADR-as-written %d\n",
        NTOT[], SCORE[:pureR], SCORE[:point], SCORE[:hybrid], SCORE[:adr])
println("="^110)
