# external/probes/pade-validity-radius/diag3.jl
#
# Gate-(v) round diagnostic.  Reads node_radii_gate_v.csv (produced by
# probe_gate_v.jl) and node_radii.csv (probe.jl) and answers three questions
# the head-to-head needs to be airtight:
#
#  A.  FINE safety sweep for gate (v-b) Jorba-Zou: the largest s with
#      over-rate <= 5% AND ~0%, with exact coverage and worst-over factors.
#  B.  The deep-wedge tail: do the 26 pole-adjacent nodes (tiny R_iv) get
#      OVER-estimated by gate (v-b)?  This is where a per-node gate could
#      betray honesty.
#  C.  A pole-adjacency clamp for (v-b): R = min(s*R_vb, 0.5*h*min|t*|),
#      does the clamp scrub the residual tail without costing coverage?
#
# Run AFTER probe_gate_v.jl:
#   julia external/probes/pade-validity-radius/diag3.jl

using Printf

const HERE = @__DIR__
const H    = 0.10

# -- load gate-v csv ---------------------------------------------------------
rows = Vector{Dict{String,Float64}}()
open(joinpath(HERE, "node_radii_gate_v.csv")) do io
    hdr = split(strip(readline(io)), ',')
    for ln in eachline(io)
        vals = parse.(Float64, split(strip(ln), ','))
        push!(rows, Dict(hdr[i] => vals[i] for i in eachindex(hdr)))
    end
end
N = length(rows)
@printf("loaded %d nodes from node_radii_gate_v.csv\n", N)

pct(v,p) = isempty(v) ? NaN : sort(v)[clamp(cld(p*length(v),100),1,length(v))]

# -- a small scorer ----------------------------------------------------------
function score(Rg::Vector{Float64}, Remp::Vector{Float64})
    over=0; n=0; cov=0.0; worst=0.0
    for k in eachindex(Rg)
        (isfinite(Rg[k]) && isfinite(Remp[k]) && Remp[k]>0) || continue
        n += 1
        ρ = Rg[k]/Remp[k]
        ρ > 1.0 && (over += 1)
        worst = max(worst, ρ)
        cov += π*Rg[k]^2
    end
    return (over=over, n=n, rate=over/max(n,1), cov=cov, worst=worst)
end

# =============================================================================
println("\n" * "="^74)
println("A.  Gate (v-b) Jorba-Zou -- FINE safety sweep (tol 1e-8)")
println("="^74)

R_iv = [rows[k]["R_iv_1e8"] for k in 1:N]
R_vb = [rows[k]["R_vb_1e8"] for k in 1:N]
R_ii = [rows[k]["R_ii"]     for k in 1:N]

# global-kappa reference
s_kap = score(fill(0.30*H, N), R_iv)
@printf("  global-kappa(0.30h) reference : over=%.1f%%  coverage=%.4f\n\n",
        100*s_kap.rate, s_kap.cov)

@printf("  %-8s %-14s %-12s %-10s %-10s\n",
        "safety","over","over-rate","coverage","worst-over")
for s in 0.50:-0.02:0.20
    Rg = Float64[s*R_vb[k] for k in 1:N]
    sc = score(Rg, R_iv)
    flag = sc.rate <= 0.05 ? (sc.rate <= 0.01 ? "  <= 1%" : "  <= 5%") : ""
    @printf("  s=%.2f   %3d/%-3d        %5.1f%%      %.4f    %6.2fx%s\n",
            s, sc.over, sc.n, 100*sc.rate, sc.cov, sc.worst, flag)
end

# =============================================================================
println("\n" * "="^74)
println("B.  The deep-wedge tail -- where do (v-b)'s over-estimates land?")
println("="^74)

# pick the safety factor that just achieves <= 5%
local s_pick = 0.40
for s in 0.50:-0.02:0.20
    sc = score(Float64[s*R_vb[k] for k in 1:N], R_iv)
    if sc.rate <= 0.05
        global s_pick = s
        break
    end
end
@printf("  using s=%.2f (the largest <=5%% safety)\n\n", s_pick)

# list every over-estimated node
println("  over-estimated nodes (R_gate > R_iv) with this gate:")
@printf("  %-5s %-9s %-9s %-9s %-9s %-9s\n",
        "node","R_iv","R_vb","R_gate","ratio","h*min|t*|")
global overs = 0
for k in 1:N
    Rg = s_pick*R_vb[k]
    (isfinite(Rg) && isfinite(R_iv[k]) && R_iv[k]>0) || continue
    if Rg > R_iv[k]
        global overs += 1
        @printf("  %-5d %-9.4f %-9.4f %-9.4f %-9.2f %-9.4f\n",
                k, R_iv[k], R_vb[k], Rg, Rg/R_iv[k], R_ii[k])
    end
end
@printf("  -> %d over-estimated nodes total\n", overs)

# how do the tiny-R_iv nodes fare?
tiny = [k for k in 1:N if isfinite(R_iv[k]) && R_iv[k] <= 0.06]
@printf("\n  deep-wedge tail: %d nodes with R_iv <= 0.06\n", length(tiny))
global tover = 0
for k in tiny
    Rg = s_pick*R_vb[k]
    isfinite(Rg) && Rg > R_iv[k] && (global tover += 1)
end
@printf("  of those, %d are over-estimated by gate (v-b) at s=%.2f\n",
        tover, s_pick)

# =============================================================================
println("\n" * "="^74)
println("C.  Pole-adjacency clamp:  R = min(s*R_vb, c*h*min|t*|)")
println("="^74)
println("  Does clamping against the nearest Q-root scrub the residual tail")
println("  WITHOUT costing coverage on the truncation-limited majority?\n")

@printf("  %-22s %-10s %-12s %-10s %-10s\n",
        "gate","over","over-rate","coverage","worst")
for (label, Rg) in (
        (@sprintf("(v-b) s=%.2f bare", s_pick),
            Float64[s_pick*R_vb[k] for k in 1:N]),
        (@sprintf("(v-b) s=%.2f clamp 0.5", s_pick),
            Float64[min(s_pick*R_vb[k], 0.5*R_ii[k]) for k in 1:N]),
        (@sprintf("(v-b) s=%.2f clamp 0.7", s_pick),
            Float64[min(s_pick*R_vb[k], 0.7*R_ii[k]) for k in 1:N]),
        # can we PUSH the safety higher if the clamp protects the tail?
        ("(v-b) s=0.46 clamp 0.5",
            Float64[min(0.46*R_vb[k], 0.5*R_ii[k]) for k in 1:N]),
        ("(v-b) s=0.50 clamp 0.5",
            Float64[min(0.50*R_vb[k], 0.5*R_ii[k]) for k in 1:N]),
    )
    sc = score(Rg, R_iv)
    @printf("  %-22s %3d/%-3d   %6.2f%%      %.4f    %6.2fx\n",
            label, sc.over, sc.n, 100*sc.rate, sc.cov, sc.worst)
end

println("\n[diag3 done]")
