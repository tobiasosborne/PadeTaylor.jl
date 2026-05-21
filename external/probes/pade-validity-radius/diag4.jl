# external/probes/pade-validity-radius/diag4.jl
#
# Gate-(v) round -- FINAL production-formula calibration.
# Reads node_radii_gate_v.csv.  Pins the clamped Jorba-Zou gate
#
#     R_gate = min( s(tol) * R_JZ ,  0.5 * h * min|t*| )
#
# at all three tolerances, sweeping s for each so the over-rate is ~<=2%
# (the same conservative tail as global-kappa's 0.8-1.3%), and reports the
# honest-coverage multiple vs global-kappa.  Also dumps the worst residual
# over-estimators per tol so the production recommendation is honest about
# the remaining tail.
#
#   julia external/probes/pade-validity-radius/diag4.jl

using Printf
const HERE = @__DIR__
const H    = 0.10
const KAPPA = Dict(1e-6=>0.40, 1e-8=>0.30, 1e-10=>0.25)

rows = Vector{Dict{String,Float64}}()
open(joinpath(HERE,"node_radii_gate_v.csv")) do io
    hdr = split(strip(readline(io)),',')
    for ln in eachline(io)
        v = parse.(Float64, split(strip(ln),','))
        push!(rows, Dict(hdr[i]=>v[i] for i in eachindex(hdr)))
    end
end
N = length(rows)
pct(v,p) = isempty(v) ? NaN : sort(v)[clamp(cld(p*length(v),100),1,length(v))]

function score(Rg, Remp)
    over=0;n=0;cov=0.0;worst=0.0;wf=Float64[]
    for k in eachindex(Rg)
        (isfinite(Rg[k])&&isfinite(Remp[k])&&Remp[k]>0) || continue
        n+=1; ρ=Rg[k]/Remp[k]
        ρ>1.0 && (over+=1; push!(wf,ρ))
        worst=max(worst,ρ); cov+=π*Rg[k]^2
    end
    (over=over,n=n,rate=over/max(n,1),cov=cov,worst=worst,
     wf90=isempty(wf) ? 0.0 : pct(sort(wf),90))
end

println("="^76)
println("FINAL gate-(v) calibration -- clamped Jorba-Zou, all tolerances")
println("  R_gate = min( s(tol)*R_JZ , 0.5*h*min|t*| )")
println("="^76)

for (tol, ivkey, vbkey) in ((1e-6,"R_iv_1e6","R_vb_1e6"),
                            (1e-8,"R_iv_1e8","R_vb_1e8"),
                            (1e-10,"R_iv_1e10","R_vb_1e10"))
    R_iv = [rows[k][ivkey] for k in 1:N]
    R_vb = [rows[k][vbkey] for k in 1:N]
    R_ii = [rows[k]["R_ii"] for k in 1:N]
    kap  = KAPPA[tol]
    s_kap = score(fill(kap*H,N), R_iv)

    @printf("\n--- tol = %.0e  (global-kappa kappa=%.2f: over=%.1f%%, cov=%.4f) ---\n",
            tol, kap, 100*s_kap.rate, s_kap.cov)
    @printf("  %-8s %-12s %-10s %-10s %-12s\n",
            "s","over-rate","coverage","worst","cov vs kappa")
    local s_at2 = NaN
    for s in 0.60:-0.02:0.20
        Rg = Float64[min(s*R_vb[k], 0.5*R_ii[k]) for k in 1:N]
        sc = score(Rg, R_iv)
        mark = ""
        if sc.rate <= 0.02 && isnan(s_at2)
            s_at2 = s; mark = "  <- <=2% pick"
        end
        @printf("  s=%.2f   %6.2f%%      %.4f    %5.2fx     x%6.2f%s\n",
                s, 100*sc.rate, sc.cov, sc.worst, sc.cov/s_kap.cov, mark)
    end
end

# -- detailed final pick at tol 1e-8 -----------------------------------------
println("\n" * "="^76)
println("PRODUCTION PICK detail -- tol 1e-8")
println("="^76)
let
    R_iv = [rows[k]["R_iv_1e8"] for k in 1:N]
    R_vb = [rows[k]["R_vb_1e8"] for k in 1:N]
    R_ii = [rows[k]["R_ii"] for k in 1:N]
    # pick s giving <=2% over
    local s_pick = 0.30
    for s in 0.60:-0.02:0.20
        sc = score(Float64[min(s*R_vb[k],0.5*R_ii[k]) for k in 1:N], R_iv)
        if sc.rate <= 0.02
            s_pick = s; break
        end
    end
    Rg = Float64[min(s_pick*R_vb[k], 0.5*R_ii[k]) for k in 1:N]
    sc = score(Rg, R_iv)
    s_kap = score(fill(0.30*H,N), R_iv)
    @printf("  picked s = %.2f\n", s_pick)
    @printf("  over-rate     : %.2f%%  (%d/%d)   [global-kappa: %.2f%%]\n",
            100*sc.rate, sc.over, sc.n, 100*s_kap.rate)
    @printf("  honest cov    : %.4f          [global-kappa: %.4f]\n",
            sc.cov, s_kap.cov)
    @printf("  coverage gain : x%.2f vs global-kappa\n", sc.cov/s_kap.cov)
    @printf("  worst-over    : %.2fx   (p90 of over-factors: %.2fx)\n",
            sc.worst, sc.wf90)
    rg = sort(filter(isfinite,Rg))
    @printf("  R_gate dist   : p10=%.3f med=%.3f p90=%.3f  (vs global 0.030)\n",
            pct(rg,10), pct(rg,50), pct(rg,90))
    # how often does the clamp bind vs the JZ term?
    local nclamp = 0
    for k in 1:N
        (isfinite(R_vb[k])&&isfinite(R_ii[k])) || continue
        0.5*R_ii[k] < s_pick*R_vb[k] && (nclamp += 1)
    end
    @printf("  clamp binds   : %d/%d nodes (%.1f%%) -- pole-adjacency term\n",
            nclamp, N, 100*nclamp/N)
end
println("\n[diag4 done]")
