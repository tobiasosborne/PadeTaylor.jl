# audition.jl — the ADR-0028 audition.
#
# Builds cells (A)=GGT diagonal (m,m), (B_ceil)=Mano–Tsuda ⌈m/d⌉ rows [the
# ADR's literal cell B], and (B_floor)=Mano–Tsuda ⌊m/d⌋+shrink-m [genuinely
# square], measures TRUE per-step error vs oracle, the held-out residual R
# (single + multi-coefficient), the conditioning gap g, and the selector's
# choice; then answers the seven sign-off questions + a pole-crossing probe.
#
# Run from repo root:
#   julia --project=. external/probes/adr0028-dual-construction-audition/audition.jl

include(joinpath(@__DIR__, "cells.jl"))
include(joinpath(@__DIR__, "battery.jl"))
using Printf
using LinearAlgebra: eigvals

function roots_of(c)
    n = length(c) - 1
    n == 0 && return ComplexF64[]
    cc = ComplexF64.(c) ./ c[end]
    C = zeros(ComplexF64, n, n)
    for i = 1:n-1; C[i+1, i] = 1; end
    for i = 1:n; C[i, n] = -cc[i]; end
    eigvals(C)
end

function pole_disagreement(QA, QB)
    ra = roots_of(QA); rb = roots_of(QB)
    (isempty(ra) || isempty(rb)) && return NaN
    maximum(minimum(abs(z - w) for w in rb) / max(abs(z), 1e-30) for z in ra)
end

lg(x) = x ≤ 0 ? -Inf : log10(x)
const REPORT = String[]
rec(s) = (println(s); push!(REPORT, s))

rec("="^104)
rec("ADR-0028 DUAL-CONSTRUCTION AUDITION")
rec("  cell A = GGT diagonal (m,m) [shipped] | B_ceil = M–T ⌈m/d⌉ [ADR's cell B] | B_floor = M–T ⌊m/d⌋+shrink [square]")
rec("="^104)
rec("")

battery = build_battery()
rows = NamedTuple[]

rec(@sprintf("%-16s %-5s %2s %3s | %-9s %-9s %-9s | %-9s %-9s | %-7s %-7s",
    "system", "kind", "d", "m", "err_A", "err_Bceil", "err_Bfloor",
    "R1_A", "R1_Bf", "g_A", "g_Bf"))
rec("-"^104)
for s in battery
    jets = rescaled_jets(s)
    cA = cell_A(jets, s.m); cBc = cell_B_ceil(jets, s.m); cBf = cell_B_floor(jets, s.m)
    eA = step_error(s, cA); eBc = step_error(s, cBc); eBf = step_error(s, cBf)
    KcF = min(cA.K, cBf.K)
    R1A = held_out_residual(cA, jets, KcF)[1]; R1Bf = held_out_residual(cBf, jets, KcF)[1]
    RmA = held_out_multi(cA, jets, KcF); RmBf = held_out_multi(cBf, jets, KcF)
    push!(rows, (; s, eA, eBc, eBf, R1A, R1Bf, RmA, RmBf,
                 gA = gap(cA), gBf = gap(cBf), gBc = gap(cBc),
                 KA = cA.K, KBf = cBf.K,
                 pd = pole_disagreement(cA.denominator, cBf.denominator)))
    rec(@sprintf("%-16s %-5s %2d %3d | %-9.2e %-9.2e %-9.2e | %-9.2e %-9.2e | %-7.1e %-7.1e",
        s.name, String(s.kind), s.d, s.m, eA, eBc, eBf, R1A, R1Bf, gap(cA), gap(cBf)))
end

# --- Q1: does the principled square cell B_floor recover accuracy? ----------
rec(""); rec("="^104)
rec("Q1 — B_floor (square M–T) vs A:  err_A/err_Bfloor  (>1 ⇒ B_floor wins)")
rec("="^104)
for r in rows
    ratio = r.eBf == 0 ? Inf : r.eA / r.eBf
    v = ratio > 3 ? "B_floor BETTER ×$(round(ratio,sigdigits=2))" :
        ratio < 1/3 ? "A BETTER ×$(round(1/ratio,sigdigits=2))" : "≈ tie"
    rec(@sprintf("  %-16s %-12s err_A=%.2e  err_Bfloor=%.2e   %s", r.s.name, String(r.s.kind), r.eA, r.eBf, v))
end

# --- Q2: does the held-out metric pick the true winner? (A vs B_floor) ------
rec(""); rec("="^104)
rec("Q2 — METRIC FIDELITY: selector(A,B_floor) vs true-error winner")
rec("    single-coeff R  vs  multi-coeff R (span 3);  Rfloor=1e-11 absolute-tie floor")
rec("="^104)
n1 = 0; nm = 0
for s in battery
    jets = rescaled_jets(s)
    cA = cell_A(jets, s.m); cBf = cell_B_floor(jets, s.m)
    eA = step_error(s, cA); eBf = step_error(s, cBf)
    truewin = eA ≤ eBf ? :diagonal : :square_floor
    isTie = abs(lg(eA) - lg(eBf)) < 0.5
    s1, _ = select([cA, cBf], jets; Rfloor = 1e-11, multi = false)
    sm, _ = select([cA, cBf], jets; Rfloor = 1e-11, multi = true)
    ok1 = s1.label == truewin || isTie; okm = sm.label == truewin || isTie
    global n1 += ok1; global nm += okm
    rec(@sprintf("  %-16s true=%-12s  R1→%-12s %s   Rmulti→%-12s %s",
        s.name, String(truewin), String(s1.label), ok1 ? "✓" : "✗",
        String(sm.label), okm ? "✓" : "✗"))
end
rec(@sprintf("  SCORE: single-coeff %d/%d   |   multi-coeff %d/%d", n1, length(battery), nm, length(battery)))

# --- Q3: is g a clean conditional-build gate? -------------------------------
rec(""); rec("="^104)
rec("Q3 — CONDITIONING GAP g_A AS A 'BUILD B' GATE")
rec("="^104)
for r in rows
    aBad = r.eA > 10 * r.eBf
    rec(@sprintf("  %-16s %-12s  g_A=%.2e   A-is-worse-than-B_floor? %s",
        r.s.name, String(r.s.kind), r.gA, aBad ? "YES (want gate→build B)" : "no"))
end
worse = [r.gA for r in rows if r.eA > 10 * r.eBf]
okA = [r.gA for r in rows if r.eA ≤ 10 * r.eBf]
rec(@sprintf("  g_A where A worse: %s", string(round.(log10.(worse), digits = 1))))
rec(@sprintf("  g_A where A ok:    %s", string(round.(log10.(okA), digits = 1))))
rec("  → if these ranges OVERLAP, g_A cannot gate the build decision.")

# --- Q5: determinism --------------------------------------------------------
rec(""); rec("="^104)
rec("Q5 — DETERMINISM under eps-perturbation (selector A vs B_floor, Rfloor=1e-11, multi)")
rec("="^104)
for s in battery
    jets = rescaled_jets(s)
    sa, _ = select([cell_A(jets, s.m), cell_B_floor(jets, s.m)], jets; Rfloor = 1e-11, multi = true)
    pert = [[cj * (1 + 50 * eps() * (isodd(k) ? 1 : -1)) for (k, cj) in enumerate(j)] for j in jets]
    sb, _ = select([cell_A(pert, s.m), cell_B_floor(pert, s.m)], pert; Rfloor = 1e-11, multi = true)
    rec(@sprintf("  %-16s chosen=%-12s perturbed=%-12s %s",
        s.name, String(sa.label), String(sb.label), sa.label == sb.label ? "STABLE" : "*** FLIPPED ***"))
end

# --- Q7: the d∤m fork (the headline) ----------------------------------------
rec(""); rec("="^104)
rec("Q7 — d∤m FORK:  B_ceil ⌈m/d⌉ (over-determined)  vs  B_floor ⌊m/d⌋+shrink (square)")
rec("="^104)
for r in rows
    dm = r.s.d * cld(r.s.m, r.s.d) == r.s.m
    diff = abs(lg(r.eBc) - lg(r.eBf))
    rec(@sprintf("  %-16s d=%d m=%d  err_Bceil=%.2e  err_Bfloor=%.2e   %s",
        r.s.name, r.s.d, r.s.m, r.eBc, r.eBf,
        dm ? "(d|m: identical)" : diff > 1 ? "DIFFER ×$(round(10.0^diff,sigdigits=2)) (ceil WORSE)" : "≈ similar"))
end

# --- Q6: disagreement diagnostic --------------------------------------------
rec(""); rec("="^104)
rec("Q6 — A/B_floor DISAGREEMENT (R_gap, pole_disagree) — does it flag ambiguous steps?")
rec("="^104)
for r in rows
    rec(@sprintf("  %-16s %-12s R_gap=R1_A/R1_Bf=%.2e  pole_disagree=%.2e  (A %s B_floor)",
        r.s.name, String(r.s.kind), r.R1Bf == 0 ? Inf : r.R1A / r.R1Bf, r.pd,
        r.eA < r.eBf ? "<" : "≥"))
end

# --- POLE-CROSSING PROBE: does diagonal (m,m) beat off-diagonal (m-1,m) across
#     a genuine transcendental pole? (the ADR's premise) -----------------------
rec(""); rec("="^104)
rec("POLE-CROSSING PROBE — tan(z), pole at π/2 ≈ 1.5708;  A=robust_pade(m,m) vs B=robust_pade(m-1,m)")
rec("  (d=1 cell-accuracy measurement; in production d=1 short-circuits to A)")
rec("="^104)
# tan jet: u'=1+u², u(0)=0 ⇒ u=tan(z).  Build the Taylor jet of tan directly.
function tan_jet(N)
    # derivative recursion via series of tan: a_{k} from u'=1+u² ; integrate
    c = zeros(Float64, N + 1)            # c[k+1] = coeff of z^k in tan
    # u' = 1 + u^2 ; (k+1) c_{k+1} = [z^k](1+u^2)
    for k = 0:N-1
        conv = sum(c[j+1] * c[k-j+1] for j = 0:k)
        rhs = (k == 0 ? 1.0 : 0.0) + conv
        c[k+2] = rhs / (k + 1)
    end
    ComplexF64.(c)
end
for h in (1.0, 1.3, 1.8, 2.0)        # 1.0,1.3 = before pole; 1.8,2.0 = past pole
    for m in (10, 15)
        jet = tan_jet(2m + 4)
        jr = [h^(k - 1) * jet[k] for k = 1:length(jet)]
        pA = robust_pade(jr, m, m; method = :svd)
        pB = robust_pade(jr, m - 1, m; method = :svd)
        evalp(p) = sum(p.a) / sum(p.b)        # P(1)/Q(1)
        tru = tan(h)
        eA = abs(evalp(pA) - tru); eB = abs(evalp(pB) - tru)
        win = eA < eB / 3 ? "A wins" : eB < eA / 3 ? "B wins" : "≈ tie"
        rec(@sprintf("  h=%.1f (z=%.2f, %s pole) m=%2d | err_A(m,m)=%.2e  err_B(m-1,m)=%.2e | %s",
            h, h, h > π / 2 ? "PAST" : "pre ", m, eA, eB, win))
    end
end

open(joinpath(@__DIR__, "REPORT.md"), "w") do io
    println(io, "# ADR-0028 dual-construction audition — raw run\n\n```")
    for l in REPORT; println(io, l); end
    println(io, "```")
end
rec(""); rec("(REPORT.md written)")
