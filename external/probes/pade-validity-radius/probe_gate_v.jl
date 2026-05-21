# external/probes/pade-validity-radius/probe_gate_v.jl
#
# Phase-A1 exploration spike (bead padetaylor-0ln.37.1) -- GATE (v) ROUND.
# THROWAWAY CODE.  No production file in src/ or figures/ is modified.  Run:
#
#     julia --project=figures external/probes/pade-validity-radius/probe_gate_v.jl
#
# ============================================================================
# MISSION -- the 5th-gate audition
# ============================================================================
# The first probe round (probe.jl, REPORT.md sections 1-4) rejected gates
# (ii) h*min|t*| and (iii) nearest-uncaptured, and recommended a GLOBAL gate
#
#     R_gate = kappa(tol) * h   with kappa(1e-8) = 0.30   (production default)
#
# The senior-NA objection (this round): kappa is a GLOBAL fudge factor, but
# probe.jl section 4.2 measured the honest radius R_iv/h varying PER NODE
# (p10-p90 = 0.6h-2.0h).  kappa=0.30 gates ~3.3x below the MEDIAN honest
# radius -- it wastes honest coverage to protect the worst ~5% of nodes.
# Each node's OWN jet coefficients determine its convergence radius; a per-node
# estimate should recover materially more honest coverage at the same safety.
#
# GATE (v): a PER-NODE convergence-radius estimate from the node's own
# coefficient decay.  Two computations are auditioned, the better recommended:
#
#   (v-a)  Cauchy-Hadamard / ratio estimate on the TRAILING coefficients of
#          the Pade numerators+denominator already stored in visited_*.
#          Zero refit cost -- a formula on coefficients in hand.
#
#   (v-b)  Reuse the project's existing machinery -- the Jorba-Zou step
#          formula (VectorStepControl.vector_step_jorba_zou): a per-node
#          truncation-limited radius from the decay of the node's Taylor jet.
#          The jet is rebuilt at order 24 (same primitive, same cost as
#          gate (ii)'s root-finding -- no fill-time root-finding beyond (ii)).
#
# Both are then compared, on the SAME 389 V8b nodes and the SAME order-64
# Taylor-oracle anchor R_iv, against the recommended global-kappa gate:
#   - over-estimation rate (the unacceptable failure -- must be <= 5%);
#   - total honest coverage (sum of honest disc areas);
#   - worst-case over-estimate factor.
# ============================================================================

using Printf
using LinearAlgebra: norm, eigvals

using PadeTaylor
using PadeTaylor.PainleveHierarchy: painleve_hierarchy,
                                    painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP:         vector_bvp_solve, VectorBVPSolution
using PadeTaylor.VectorProblems:    VectorPadeTaylorProblem
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients
using PadeTaylor.VectorStepControl: vector_step_jorba_zou

# -- companion-matrix root finder (verbatim from probe.jl) -------------------
function poly_roots(c::AbstractVector)
    cc = collect(ComplexF64, c)
    last = findlast(x -> abs(x) > 1.0e-14 * (1 + maximum(abs, cc)), cc)
    last === nothing && return ComplexF64[]
    cc = cc[1:last]
    n = length(cc) - 1
    n < 1 && return ComplexF64[]
    a = cc ./ cc[end]
    C = zeros(ComplexF64, n, n)
    for i in 1:n-1
        C[i+1, i] = 1.0
    end
    for i in 1:n
        C[i, n] = -a[i]
    end
    return eigvals(C)
end

# -- V8b walk parameters (verbatim from probe.jl) ----------------------------
const KKG_T          = 0.0
const KKG_X_L        = -20.0 + 0.0im
const KKG_X_R        =  -2.0 + 0.0im
const KKG_BVP_N      = 128
const KKG_BVP_TOL    = 1.0e-9
const KKG_BVP_MAXITER = 40
const KKG_D          = 4
const KKG_Z_SEED     = -3.0 + 0.0im
const KKG_PN_ORDER   = 24
const KKG_PN_H       = 0.1
const KKG_TARGET_RADII  = (2.0, 4.0, 6.0, 8.0)
const KKG_TARGET_ANGLES = (-0.5, -0.25, 0.0, 0.25, 0.5)

const ORACLE_ORDER = 64
const ORACLE_LOW   = 56
const TOLS         = (1.0e-6, 1.0e-8, 1.0e-10)
const NDIR         = 24
const RAY_DR       = 0.02
const RAY_MAX      = 3.0
const CONV_TOL     = 1.0e-12

println("="^78)
println("Phase-A1 GATE (v) round  (bead padetaylor-0ln.37.1)")
println("="^78)

# ----------------------------------------------------------------------------
# 1.  Reproduce the V8b walk (verbatim from probe.jl).
# ----------------------------------------------------------------------------
println("\n[1] Reproducing the V8b P_I^(2) wedge walk ...")

f  = painleve_hierarchy(:I, 2; t = KKG_T)
Jf = painleve_hierarchy_jacobian(:I, 2; t = KKG_T)

local bvp
let
    global bvp
    CT = ComplexF64
    Ba = zeros(CT, KKG_D, KKG_D); Bb = zeros(CT, KKG_D, KKG_D)
    Ba[1,1] = 1; Ba[2,2] = 1
    Bb[3,1] = 1; Bb[4,2] = 1
    seed_l = pI2_tritronquee_ic(KKG_X_L; t = KKG_T, n_terms = 2)
    seed_r = pI2_tritronquee_ic(KKG_X_R; t = KKG_T, n_terms = 2)
    g = CT[seed_l[1], seed_l[2], seed_r[1], seed_r[2]]
    bvp = vector_bvp_solve(f, KKG_X_L, KKG_X_R, Ba, Bb, g;
                           N = KKG_BVP_N, tol = KKG_BVP_TOL,
                           maxiter = KKG_BVP_MAXITER, jacobian = Jf,
                           initial_guess = z -> pI2_tritronquee_ic(z;
                                                t = KKG_T, n_terms = 2))
end
@printf("    BVP: %d Newton iters, residual_inf = %.3e\n",
        bvp.iterations, bvp.residual_inf)

y_seed = ComplexF64.(bvp(ComplexF64(KKG_Z_SEED)))
prob = VectorPadeTaylorProblem(f, y_seed,
                               (ComplexF64(KKG_Z_SEED), ComplexF64(8.0+0.0im));
                               order = KKG_PN_ORDER)
targets = ComplexF64[r*cis(θ) for r in KKG_TARGET_RADII
                               for θ in KKG_TARGET_ANGLES]
walk = vector_path_network_solve(prob, targets;
                                 order = KKG_PN_ORDER, h = KKG_PN_H)
N = length(walk.visited_z)
@printf("    walk: %d visited nodes (h = %.2f, order = %d)\n",
        N, KKG_PN_H, KKG_PN_ORDER)

# ----------------------------------------------------------------------------
# 2.  Evaluation helpers (verbatim from probe.jl).
# ----------------------------------------------------------------------------
function evalpoly_lh(c::AbstractVector, t)
    s = zero(promote_type(eltype(c), typeof(t)))
    @inbounds for k in length(c):-1:1
        s = s*t + c[k]
    end
    return s
end

function node_pade_eval(walk, k, z)
    z_v = walk.visited_z[k]
    h_v = real(walk.visited_h[k])
    t   = (z - z_v)/h_v
    Q   = walk.visited_denominator[k]
    q_t = evalpoly_lh(Q, t)
    return ComplexF64[evalpoly_lh(num, t)/q_t for num in walk.visited_numerators[k]]
end

function node_Q_roots(walk, k)
    Q = walk.visited_denominator[k]
    length(Q) < 2 && return (ComplexF64[], ComplexF64[])
    ts = poly_roots(Q)
    isempty(ts) && return (ComplexF64[], ComplexF64[])
    z_v = walk.visited_z[k]; h_v = real(walk.visited_h[k])
    zs = ComplexF64[z_v + h_v*t for t in ts]
    return (ts, zs)
end

function rel_disagree(yA::Vector{ComplexF64}, yB::Vector{ComplexF64})
    den = norm(yA) + norm(yB)
    den == 0 && return 0.0
    return norm(yA .- yB)/den
end

# ----------------------------------------------------------------------------
# 3.  Order-64 Taylor-series oracle + gate (iv) anchor (verbatim from probe.jl).
# ----------------------------------------------------------------------------
function node_oracle_jet(walk, k)
    z_v = ComplexF64(walk.visited_z[k])
    y_v = ComplexF64.(walk.visited_y[k])
    h_v = real(walk.visited_h[k])
    local jets
    try
        jets = vector_taylor_coefficients(f, z_v, y_v, ORACLE_ORDER)
    catch
        return nothing
    end
    return [ComplexF64[jet[q+1]*h_v^q for q in 0:length(jet)-1] for jet in jets]
end

function taylor_partial(jet_i::Vector{ComplexF64}, t, deg::Int)
    s = zero(ComplexF64)
    @inbounds for q in min(deg, length(jet_i)-1):-1:0
        s = s*t + jet_i[q+1]
    end
    return s
end

function gate_iv(walk, k, jet, tol)
    jet === nothing && return NaN
    z_v = ComplexF64(walk.visited_z[k])
    h_v = real(walk.visited_h[k])
    Rmin = Inf
    for d in 0:NDIR-1
        dir = cis(2π*d/NDIR)
        r = RAY_DR
        last_conv_r = 0.0
        ray_cross = Inf
        while r <= RAY_MAX
            z = z_v + r*dir
            t = (z - z_v)/h_v
            refy = Vector{ComplexF64}(undef, length(jet))
            lowy = Vector{ComplexF64}(undef, length(jet))
            for i in eachindex(jet)
                refy[i] = taylor_partial(jet[i], t, ORACLE_ORDER)
                lowy[i] = taylor_partial(jet[i], t, ORACLE_LOW)
            end
            conv = rel_disagree(refy, lowy) < CONV_TOL
            if !conv
                break
            end
            last_conv_r = r
            local pady, ok
            ok = true
            try
                pady = node_pade_eval(walk, k, z)
            catch
                ok = false
            end
            if !ok || any(!isfinite, pady)
                ray_cross = r; break
            end
            if rel_disagree(pady, refy) > tol
                ray_cross = r; break
            end
            r += RAY_DR
        end
        ray_R = isfinite(ray_cross) ? ray_cross : last_conv_r
        Rmin = min(Rmin, ray_R)
    end
    return isfinite(Rmin) ? Rmin : NaN
end

# ============================================================================
# 4.  GATE (v) -- per-node convergence radius from coefficient decay.
# ============================================================================
#
# Both variants estimate the truncation-limited radius (section 3.3 of
# REPORT.md: R_honest is bound by jet order + coefficient decay, NOT pole
# geometry).  Both are stated in the rescaled variable t = (z - z_v)/h_v and
# then multiplied by h_v to give a z-plane radius -- so they are directly
# comparable to R_iv and R_gate.
#
# ----------------------------------------------------------------------------
# 4a.  Gate (v-a): Cauchy-Hadamard / ratio estimate on the STORED Pade
#      coefficients.  Zero refit cost.
#
# The node stores, in the rescaled t-variable, d numerator polys P_i(t) (deg
# <= 12) and one shared denominator Q(t) (deg <= 12).  The natural
# convergence-radius proxy for the underlying *function* is the radius at
# which the order-24 Pade's TRUNCATION TAIL stops being negligible.  But the
# Pade numerator/denominator are only degree-12 each -- a [12/12] form has no
# "trailing tail" the way an order-64 Taylor jet does; its coefficients are
# the SOLUTION of a linear system, not samples of an analytic function, so a
# raw Cauchy-Hadamard limsup on 12 noisy coefficients is ill-conditioned.
#
# What IS available cheaply and IS a faithful coefficient-decay signal: the
# *rescaled Taylor jet* the Pade was built from.  REPORT.md section 2.2: the
# canonical store rebuilds a Taylor jet via `vector_taylor_coefficients` and
# then does the GGT linear solve.  The rescaled Taylor coefficients
# c_k^(i) * h_v^k decay at a rate set by 1/R_taylor where R_taylor is the
# t-plane distance to the nearest singularity.  A Cauchy-Hadamard estimate
#
#     R_CH = ( vnorm(c_K) )^(-1/K)         (K = last reliable index)
#
# on the *trailing* rescaled-Taylor coefficient is the cheap honest signal.
# We take a 3-term geometric-mean ratio for noise robustness.  This needs the
# order-24 jet -- rebuilt with the production primitive, the SAME cost as
# gate (ii)'s root-finding (REPORT.md sec.4.1 budget).  No NEW machinery.
#
# 4b.  Gate (v-b): the Jorba-Zou step formula -- the project's OWN per-node
#      truncation-radius estimator (VectorStepControl.vector_step_jorba_zou).
#      h_JZ = min over k in {p-1,p} of (eps / vnorm(c_k))^(1/k), on the SAME
#      rescaled order-24 jet.  This is literally "how far can I march before
#      the order-24 truncation error exceeds eps" -- exactly R_honest's
#      definition.  eps is the gate tolerance.  This is reuse of existing,
#      tested machinery -- the recommended route if it wins.
# ----------------------------------------------------------------------------

# -- rebuild the rescaled order-24 Taylor jet for node k (the jet the
#    canonical Pade was built from; same primitive as the oracle, order 24).
function node_taylor_jet24(walk, k)
    z_v = ComplexF64(walk.visited_z[k])
    y_v = ComplexF64.(walk.visited_y[k])
    h_v = real(walk.visited_h[k])
    local jets
    try
        jets = vector_taylor_coefficients(f, z_v, y_v, KKG_PN_ORDER)
    catch
        return nothing
    end
    # rescale: coefficient of t^q is c_q * h_v^q
    return [ComplexF64[jet[q+1]*h_v^q for q in 0:length(jet)-1] for jet in jets]
end

# coefficient vector across components: vnorm of c_k = [jet1[k+1],...]
function coef_norm(jet24, k::Int)
    return norm(ComplexF64[jet24[i][k+1] for i in eachindex(jet24)])
end

# -- Gate (v-a): Cauchy-Hadamard ratio estimate on trailing rescaled-Taylor
#    coefficients.  R_CH (t-plane) then * h_v.  tol-INDEPENDENT (it estimates
#    the convergence radius, then a safety factor maps it to a gate).
function gate_v_cauchy_hadamard(jet24, h_v)
    jet24 === nothing && return NaN
    p = length(jet24[1]) - 1            # = 24
    # trailing window: indices p-2, p-1, p (the most reliable decay signal,
    # least contaminated by the leading low-order behaviour).
    Rs = Float64[]
    for k in (p-2, p-1, p)
        nk = coef_norm(jet24, k)
        nk > 0 || continue
        push!(Rs, nk^(-1.0/k))          # Cauchy-Hadamard root estimate
    end
    isempty(Rs) && return NaN
    # geometric mean over the trailing window (noise-robust)
    R_t = exp(sum(log, Rs)/length(Rs))
    return R_t * h_v
end

# -- Gate (v-b): the project's Jorba-Zou step formula on the rescaled jet.
#    vector_step_jorba_zou expects the jet in the variable it returns h IN;
#    we pass the rescaled jet (coeffs of t^q), so it returns a radius in t,
#    which we then * h_v.  eps = gate tolerance.
function gate_v_jorba_zou(jet24, h_v, tol)
    jet24 === nothing && return NaN
    local h_t
    try
        h_t = vector_step_jorba_zou(jet24, tol)
    catch
        return NaN
    end
    return real(h_t) * h_v
end

# ----------------------------------------------------------------------------
# 5.  Compute all gates for all nodes.
# ----------------------------------------------------------------------------
println("\n[2] Computing gates (iv anchor, v-a, v-b) for $N nodes ...")

R_ii   = fill(NaN, N)                       # for the pole-adjacency clamp
R_iv   = Dict(tol => fill(NaN, N) for tol in TOLS)
R_va   = fill(NaN, N)                       # tol-independent CH estimate
R_vb   = Dict(tol => fill(NaN, N) for tol in TOLS)

for k in 1:N
    h_v      = real(walk.visited_h[k])
    ts, _    = node_Q_roots(walk, k)
    R_ii[k]  = isempty(ts) ? Inf : h_v * minimum(abs, ts)

    jet24    = node_taylor_jet24(walk, k)
    R_va[k]  = gate_v_cauchy_hadamard(jet24, h_v)

    jetO     = node_oracle_jet(walk, k)
    for tol in TOLS
        R_iv[tol][k] = gate_iv(walk, k, jetO, tol)
        R_vb[tol][k] = gate_v_jorba_zou(jet24, h_v, tol)
    end
    if k % 50 == 0
        @printf("    ... node %d/%d done\n", k, N)
    end
end
println("    all nodes processed.")

# ----------------------------------------------------------------------------
# 6.  Comparison harness.
# ----------------------------------------------------------------------------
pct(v, p) = isempty(v) ? NaN : sort(v)[clamp(cld(p*length(v),100),1,length(v))]

# Evaluate a gate-radius vector against the anchor R_emp:
#   over-rate   = fraction of nodes with R_gate > R_emp (DISHONEST)
#   coverage    = sum of pi*R_gate^2 over honest-comparable nodes
#   worst-over  = max R_gate/R_emp
function score(Rg::Vector{Float64}, Remp::Vector{Float64}, N)
    over = 0; ncomp = 0
    cov = 0.0; worst = 0.0
    overfac = Float64[]
    for k in 1:N
        (isfinite(Rg[k]) && isfinite(Remp[k]) && Remp[k] > 0 && Rg[k] >= 0) || continue
        ncomp += 1
        ρ = Rg[k]/Remp[k]
        if ρ > 1.0
            over += 1
            push!(overfac, ρ)
        end
        worst = max(worst, ρ)
        cov += π * Rg[k]^2
    end
    return (over=over, ncomp=ncomp, overrate=over/max(ncomp,1),
            coverage=cov, worst=worst,
            overfac90 = isempty(overfac) ? 0.0 : pct(sort(overfac),90))
end

function report_gate(name, Rg, Remp, N; ref=nothing)
    s = score(Rg, Remp, N)
    rg_med = pct(sort(filter(isfinite, Rg)), 50)
    line = @sprintf("  %-26s : OVER %3d/%-3d (%5.1f%%)  cov=%.4f  worst=%6.2fx  R_gate med=%.3f",
                     name, s.over, s.ncomp, 100*s.overrate, s.coverage,
                     s.worst, rg_med)
    if ref !== nothing && ref > 0
        line *= @sprintf("  cov x%.2f vs kappa", s.coverage/ref)
    end
    println(line)
    return s
end

println("\n" * "="^78)
println("RESULTS -- gate (v) vs recommended global-kappa")
println("="^78)

const KAPPA = Dict(1.0e-6 => 0.40, 1.0e-8 => 0.30, 1.0e-10 => 0.25)

for tol in TOLS
    Remp = R_iv[tol]
    @printf("\n--- tol = %.0e  (anchor R_iv: n=%d/%d finite, median=%.3f) ---\n",
            tol, count(isfinite, Remp), N, pct(sort(filter(isfinite,Remp)),50))

    # recommended global-kappa gate (probe.jl section 4.1)
    kap = KAPPA[tol]
    R_kappa = fill(kap*KKG_PN_H, N)
    s_kappa = report_gate(@sprintf("global-kappa (kappa=%.2f)", kap),
                          R_kappa, Remp, N)
    refcov = s_kappa.coverage

    # recommended HYBRID from probe.jl section 4.3:
    #   R = min(kappa*h, 0.5*h*min|t*|)
    R_hyb = Float64[min(kap*KKG_PN_H, 0.5*R_ii[k]) for k in 1:N]
    report_gate("global-kappa hybrid", R_hyb, Remp, N; ref=refcov)

    # gate (v-a): Cauchy-Hadamard, tol-independent.  Needs a safety factor s
    # mapping the convergence-radius estimate to an honest gate.  Sweep s.
    println("  -- gate (v-a) Cauchy-Hadamard, safety sweep --")
    for s in (1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2)
        Rg = Float64[s*R_va[k] for k in 1:N]
        report_gate(@sprintf("    (v-a) CH x s=%.2f", s), Rg, Remp, N; ref=refcov)
    end

    # gate (v-b): Jorba-Zou step formula at eps = tol.  Also sweep a safety
    # factor: the JZ step is itself an honest truncation bound at eps=tol, so
    # s=1.0 is the "raw" reuse; smaller s adds margin.
    println("  -- gate (v-b) Jorba-Zou step (eps=tol), safety sweep --")
    for s in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
        Rg = Float64[s*R_vb[tol][k] for k in 1:N]
        report_gate(@sprintf("    (v-b) JZ  x s=%.2f", s), Rg, Remp, N; ref=refcov)
    end
end

# ----------------------------------------------------------------------------
# 7.  Per-node correlation: does gate (v) TRACK R_iv, or just hit it on average?
#     A per-node gate is only worth it if its ratio to R_iv has SMALL spread.
# ----------------------------------------------------------------------------
println("\n" * "="^78)
println("PER-NODE TRACKING -- ratio gate/R_iv distribution (the spread test)")
println("="^78)
println("A per-node gate beats global-kappa only if ratio spread is TIGHTER")
println("than R_iv/h's own spread (p10-p90 = 0.6-2.0, a 3.3x span).")

for tol in TOLS
    Remp = R_iv[tol]
    @printf("\n--- tol = %.0e ---\n", tol)
    for (name, Rg) in (("R_iv/h (global-kappa basis)",
                         Float64[isfinite(Remp[k]) ? Remp[k]/KKG_PN_H : NaN for k in 1:N]),
                        ("R_va/R_iv (CH, raw)", Float64[R_va[k]/Remp[k] for k in 1:N]),
                        ("R_vb/R_iv (JZ, raw)", Float64[R_vb[tol][k]/Remp[k] for k in 1:N]))
        v = sort(filter(isfinite, Rg))
        isempty(v) && continue
        @printf("  %-30s : p10=%.2f p25=%.2f med=%.2f p75=%.2f p90=%.2f  span=%.1fx\n",
                name, pct(v,10), pct(v,25), pct(v,50), pct(v,75), pct(v,90),
                pct(v,90)/max(pct(v,10),1e-9))
    end
end

# ----------------------------------------------------------------------------
# 8.  FINAL head-to-head at tol 1e-8: pick the best honest variant of (v) and
#     compare honest coverage vs global-kappa, holding over-rate <= 5%.
# ----------------------------------------------------------------------------
println("\n" * "="^78)
println("HEAD-TO-HEAD  (tol = 1e-8, over-rate budget <= 5%)")
println("="^78)

let tol = 1.0e-8
    Remp = R_iv[tol]
    kap  = KAPPA[tol]
    R_kappa = fill(kap*KKG_PN_H, N)
    s_kappa = score(R_kappa, Remp, N)
    @printf("  global-kappa (0.30h)      : over=%.1f%%  coverage=%.4f\n",
            100*s_kappa.overrate, s_kappa.coverage)

    # for each (v) variant, pick the LARGEST safety factor with over-rate<=5%.
    function best_honest(label, Rraw_fn, safeties)
        best = nothing
        for s in safeties
            Rg = Float64[s*Rraw_fn(k) for k in 1:N]
            sc = score(Rg, Remp, N)
            if sc.overrate <= 0.05
                best = (s=s, sc=sc); break
            end
        end
        if best === nothing
            @printf("  %-26s : NO safety factor achieves <=5%% over-rate\n", label)
            return nothing
        end
        @printf("  %-26s : s=%.2f  over=%.1f%%  coverage=%.4f  (x%.2f vs kappa)  worst=%.2fx\n",
                label, best.s, 100*best.sc.overrate, best.sc.coverage,
                best.sc.coverage/s_kappa.coverage, best.sc.worst)
        return best
    end

    b_va = best_honest("(v-a) Cauchy-Hadamard", k->R_va[k],
                        (1.0,0.9,0.8,0.7,0.6,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1))
    b_vb = best_honest("(v-b) Jorba-Zou", k->R_vb[tol][k],
                        (1.0,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.4,0.3))
end

# ----------------------------------------------------------------------------
# 9.  Machine-readable dump.
# ----------------------------------------------------------------------------
open(joinpath(@__DIR__, "node_radii_gate_v.csv"), "w") do io
    println(io, "node,re_z,im_z,abs_z,R_ii,R_va,",
                "R_iv_1e8,R_vb_1e8,R_iv_1e6,R_vb_1e6,R_iv_1e10,R_vb_1e10")
    for k in 1:N
        z = walk.visited_z[k]
        @printf(io, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                k, real(z), imag(z), abs(z), R_ii[k], R_va[k],
                R_iv[1e-8][k], R_vb[1e-8][k], R_iv[1e-6][k], R_vb[1e-6][k],
                R_iv[1e-10][k], R_vb[1e-10][k])
    end
end
println("\n  wrote node_radii_gate_v.csv")
println("[done]")
