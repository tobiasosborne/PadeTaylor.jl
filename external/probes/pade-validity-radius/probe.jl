# external/probes/pade-validity-radius/probe.jl
#
# Phase-A1 exploration spike (bead padetaylor-0ln.37.1).  THROWAWAY CODE.
# No production file in src/ or figures/ is modified.  Run with:
#
#     julia --project=figures external/probes/pade-validity-radius/probe.jl
#
# ============================================================================
# MISSION
# ============================================================================
# The headline P_I^(2) tritronquee figure renders the wedge pole field via a
# Stage-2 fine-grid fill that, per node, gates Pade validity at the STEP SIZE
# h (`abs(z_f - z_v) > h_v` -> NaN unless extrapolate=true).  The figure
# currently sets extrapolate=true, evaluating shared-Q Pade approximants far
# outside any verified disc.  ADR-0025 Lever 1 retires that corner.
#
# This spike measures the HONEST radius of validity R_honest of each visited
# node's shared-Q Pade, auditioning four candidate gate criteria against
# measured ground truth:
#
#   Gate (i)   R_i   = h_v                       -- the current step gate.
#   Gate (ii)  R_ii  = h_v * min|t*|             -- nearest denominator root.
#   Gate (iii) R_iii = distance to nearest UNCAPTURED singularity.
#   Gate (iv)  R_iv  = EMPIRICAL ground truth.
#
# ============================================================================
# THE EMPIRICAL ANCHOR -- design, after the diag.jl investigation
# ============================================================================
# diag.jl showed an early draft pinned R_iv == RAY_DR for every node.  Root
# cause (NOT a data property): the overlap test marched a ray from node A
# toward a *distant* neighbour B and took the min first-crossing over all
# neighbours.  A far neighbour B's Pade is wrong near z_A -- but that is B's
# error, not A's -- so one bad-but-far B collapsed the anchor.  The overlap
# test only legitimately bounds node A where BOTH A and B are individually
# trustworthy; it cannot use a far node as an oracle.
#
# Corrected design.  The empirical anchor R_iv(node k) is the radius at which
# node k's OWN order-24 shared-Q Pade first loses tol accuracy, measured by
# the order-64 Taylor-series oracle:
#
#   * Build the SAME Taylor jet the canonical Pade was built from, but to
#     order 64 instead of 24.  The order-64 truncated partial sum is an
#     algorithm-DISJOINT reference inside its radius of convergence (= the
#     distance to the nearest singularity): it never saw a Pade linear solve.
#   * March radius r outward over a 24-direction fan.  At each (r,dir) compare
#     the order-24 Pade against the order-64 Taylor sum, but ONLY where the
#     Taylor sum has converged (degree-64 vs degree-56 partial sums agree).
#   * R_iv = min over directions of the first r where rel-error exceeds tol.
#
# This is a faithful LOWER bound on R_honest: the Pade legitimately extends
# accuracy beyond the Taylor radius through a pole it modelled, where the
# oracle simply goes silent -- it never over-reports.  An adjacent-pair
# overlap cross-check (genuinely neighbouring nodes only) corroborates it.
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

# Polynomial root finding via the companion matrix.  Polynomials.jl is a
# transitive dep of PadeTaylor but not a direct figures-env dep; rather than
# perturb the env we use a self-contained companion-eigenvalue root finder --
# numerically identical to `Polynomials.roots` (which itself builds the
# companion matrix and calls `eigvals`).  `c` is low-to-high coefficients.
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

# ----------------------------------------------------------------------------
# 0.  V8b walk parameters -- copied verbatim from figures/_kkg_pi2_helpers.jl.
# ----------------------------------------------------------------------------
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

# Probe knobs.
const ORACLE_ORDER = 64       # high-order Taylor jet for the empirical anchor
const ORACLE_LOW   = 56       # convergence-check companion degree
const TOLS         = (1.0e-6, 1.0e-8, 1.0e-10)
const NDIR         = 24       # directional fan for the radial march
const RAY_DR       = 0.02     # radial march resolution
const RAY_MAX      = 3.0      # how far out a ray probes
const CONV_TOL     = 1.0e-12  # Taylor-sum convergence threshold (relative)

println("="^78)
println("Phase-A1 Pade validity-radius spike  (bead padetaylor-0ln.37.1)")
println("="^78)

# ----------------------------------------------------------------------------
# 1.  Reproduce the V8b walk.
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
# 2.  Evaluation helpers.
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
# 3.  Order-64 Taylor-series oracle per node.
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

# ----------------------------------------------------------------------------
# 4.  Gate (iv) -- the empirical anchor (corrected).
#
# March radius r outward over an NDIR-direction fan.  At each (r,dir) the
# order-64 Taylor sum is the reference -- but ONLY while it is converged
# (deg-64 vs deg-56 partial sums agree to CONV_TOL).  R_iv = min over
# directions of the first r where the order-24 Pade's rel-error vs the oracle
# exceeds tol.  A ray whose oracle goes silent (Taylor diverged) before any
# crossing contributes its last-converged radius as a lower bound.
# ----------------------------------------------------------------------------
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
            # converged iff the high/low partial sums agree
            conv = rel_disagree(refy, lowy) < CONV_TOL
            if !conv
                break                       # oracle no longer trustworthy
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

# ----------------------------------------------------------------------------
# 5.  Adjacent-pair overlap cross-check.
#
# An independent corroboration of the oracle anchor that uses NO Taylor jet.
# For genuinely adjacent node pairs (|z_A - z_B| < 1.5*h, i.e. one walk step)
# both Pades are individually trustworthy near the segment; the radius from
# z_A at which the two first disagree past tol is a second R_honest estimate.
# Reported as an aggregate sanity figure, not folded into the per-node anchor.
# ----------------------------------------------------------------------------
function overlap_crosscheck(walk, tol)
    crossings = Float64[]
    for k in 2:length(walk.visited_z)
        p = walk.visited_parent[k]
        p == 0 && continue
        z_A = ComplexF64(walk.visited_z[k])
        z_B = ComplexF64(walk.visited_z[p])
        dAB = abs(z_B - z_A)
        (0.3*KKG_PN_H < dAB < 1.5*KKG_PN_H) || continue
        dir = (z_B - z_A)/dAB
        first_cross = RAY_MAX
        r = RAY_DR
        while r <= RAY_MAX
            z = z_A + r*dir
            local yA, yB, ok
            ok = true
            try
                yA = node_pade_eval(walk, k, z)
                yB = node_pade_eval(walk, p, z)
            catch
                ok = false
            end
            if !ok || any(!isfinite, yA) || any(!isfinite, yB)
                first_cross = r; break
            end
            if rel_disagree(yA, yB) > tol
                first_cross = r; break
            end
            r += RAY_DR
        end
        push!(crossings, first_cross)
    end
    return crossings
end

# ----------------------------------------------------------------------------
# 6.  Gate (iii) -- distance to nearest UNCAPTURED singularity.
#
# Operational definition.  Build a GLOBAL pole reference set: the union over
# ALL nodes of every Q-root z-image with |t*| <= 3 (a root within 3 steps of
# its own centre is a well-placed pole estimate), greedy-clustered, keeping
# clusters with cross-node support >= 2 (the spurious-root guard).  An
# [<=12/12] Pade can resolve at most m = order/2 = 12 poles; node k's own
# Q-roots with |t*| <= 3 are the poles k CAPTURED and stays accurate THROUGH.
# R_iii = z-plane distance from z_node to the nearest GLOBAL pole NOT among
# k's captured set -- the nearest singularity k's Q did not model.
# ----------------------------------------------------------------------------
const GLOBAL_POLE_TOL = 0.18

function build_global_poles(walk)
    raw = ComplexF64[]
    for k in 1:length(walk.visited_z)
        ts, zs = node_Q_roots(walk, k)
        for (t, z) in zip(ts, zs)
            abs(t) <= 3.0 || continue
            push!(raw, z)
        end
    end
    reps = ComplexF64[]; cnt = Int[]
    for z in raw
        j = findfirst(r -> abs(z-r) <= GLOBAL_POLE_TOL, reps)
        if j === nothing
            push!(reps, z); push!(cnt, 1)
        else
            cnt[j] += 1
        end
    end
    return ComplexF64[reps[j] for j in eachindex(reps) if cnt[j] >= 2]
end

function gate_iii(walk, k, global_poles)
    z_v = ComplexF64(walk.visited_z[k])
    ts, zs = node_Q_roots(walk, k)
    captured = ComplexF64[zs[i] for i in eachindex(zs) if abs(ts[i]) <= 3.0]
    best = Inf
    for gp in global_poles
        any(c -> abs(gp-c) <= GLOBAL_POLE_TOL, captured) && continue
        d = abs(gp - z_v)
        d < best && (best = d)
    end
    return best
end

# ----------------------------------------------------------------------------
# 7.  Compute all gates.
# ----------------------------------------------------------------------------
println("\n[2] Computing four validity radii for each of $N nodes ...")

global_poles = build_global_poles(walk)
@printf("    global pole reference set: %d clustered poles (support>=2)\n",
        length(global_poles))

R_i   = fill(NaN, N)
R_ii  = fill(NaN, N)
R_iii = fill(NaN, N)
R_iv  = Dict(tol => fill(NaN, N) for tol in TOLS)

for k in 1:N
    h_v = real(walk.visited_h[k])
    R_i[k]   = h_v
    ts, _    = node_Q_roots(walk, k)
    R_ii[k]  = isempty(ts) ? Inf : h_v * minimum(abs, ts)
    R_iii[k] = gate_iii(walk, k, global_poles)
    jet      = node_oracle_jet(walk, k)
    for tol in TOLS
        R_iv[tol][k] = gate_iv(walk, k, jet, tol)
    end
    if k % 50 == 0
        @printf("    ... node %d/%d done\n", k, N)
    end
end
println("    all nodes processed.")

xcheck = Dict(tol => overlap_crosscheck(walk, tol) for tol in TOLS)

# ----------------------------------------------------------------------------
# 8.  Aggregate + report.
# ----------------------------------------------------------------------------
pct(v, p) = isempty(v) ? NaN : sort(v)[clamp(cld(p*length(v),100),1,length(v))]
stats(v)  = (lo=minimum(v), p10=pct(v,10), med=pct(v,50),
             mean=sum(v)/length(v), p90=pct(v,90), hi=maximum(v))

println("\n" * "="^78)
println("RESULTS")
println("="^78)

@printf("\nGate (i)   R_i  = h_v               : constant %.3f (all %d nodes)\n",
        KKG_PN_H, N)
let s = stats(filter(isfinite, R_ii))
    @printf("Gate (ii)  R_ii = h*min|t*|         : lo=%.3f p10=%.3f med=%.3f mean=%.3f p90=%.3f hi=%.3f\n",
            s.lo, s.p10, s.med, s.mean, s.p90, s.hi)
end
let s = stats(filter(isfinite, R_iii))
    @printf("Gate (iii) R_iii= nearest-uncaptured: lo=%.3f p10=%.3f med=%.3f mean=%.3f p90=%.3f hi=%.3f\n",
            s.lo, s.p10, s.med, s.mean, s.p90, s.hi)
end

for tol in TOLS
    Remp = R_iv[tol]
    v = filter(isfinite, Remp)
    s = stats(v)
    @printf("\n--- Empirical anchor R_iv  (tol = %.0e) ---\n", tol)
    @printf("  R_iv (oracle): lo=%.3f p10=%.3f med=%.3f mean=%.3f p90=%.3f hi=%.3f  (n=%d/%d finite)\n",
            s.lo, s.p10, s.med, s.mean, s.p90, s.hi, length(v), N)
    xc = filter(isfinite, xcheck[tol])
    if !isempty(xc)
        sx = stats(xc)
        @printf("  overlap cross-check (adjacent pairs, n=%d): p10=%.3f med=%.3f p90=%.3f\n",
                length(xc), sx.p10, sx.med, sx.p90)
    end

    # Predictive quality.  ratio = R_gate / R_emp.
    #   ratio < 1  => UNDERESTIMATE (wastes coverage -- too many NaN slots).
    #   ratio > 1  => OVERESTIMATE (DISHONEST -- Pade used where it is wrong).
    for (name, Rg) in (("(i)  h_v", R_i), ("(ii) h*min|t*|", R_ii),
                        ("(iii) uncaptured", R_iii))
        ratios = Float64[]
        overcount = 0
        for k in 1:N
            (isfinite(Rg[k]) && isfinite(Remp[k]) && Remp[k] > 0) || continue
            ρ = Rg[k]/Remp[k]
            push!(ratios, ρ)
            ρ > 1.0 && (overcount += 1)
        end
        if isempty(ratios)
            @printf("    gate %-18s : no comparable nodes\n", name)
        else
            sr = sort(ratios)
            @printf("    gate %-18s : ratio p10=%.2f med=%.2f p90=%.2f  worst-over=%.2fx  OVER=%d/%d\n",
                    name, pct(sr,10), pct(sr,50), pct(sr,90),
                    maximum(sr), overcount, length(ratios))
        end
    end
end

# Honest-coverage-radius model + recommended-gate evaluation.
let tol = 1.0e-8
    Remp = filter(isfinite, R_iv[tol])
    ratios = sort(Remp ./ KKG_PN_H)
    @printf("\n--- Honest-coverage model (tol 1e-8) ---\n")
    @printf("  R_honest / h : p10=%.2f p25=%.2f med=%.2f p75=%.2f p90=%.2f\n",
            pct(ratios,10), pct(ratios,25), pct(ratios,50),
            pct(ratios,75), pct(ratios,90))

    # Candidate production gate: R_gate = SAFETY * h * min|t*|, capped, with
    # the safety factor chosen so it never over-estimates R_iv.  Find the
    # largest SAFETY in {0.1,...,1.0} with zero over-nodes.
    @printf("\n--- Recommended-gate calibration (gate (ii) family) ---\n")
    Remp_all = R_iv[tol]
    for safety in (1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2)
        over = 0; covered = Float64[]
        for k in 1:N
            (isfinite(R_ii[k]) && isfinite(Remp_all[k]) && Remp_all[k] > 0) || continue
            Rg = safety * R_ii[k]
            Rg > Remp_all[k] && (over += 1)
            push!(covered, Rg / KKG_PN_H)
        end
        @printf("  safety=%.2f : OVER-nodes=%3d/%d   median R_gate/h=%.2f\n",
                safety, over, length(covered), pct(sort(covered),50))
    end
end

# Representative nodes.
println("\n--- Representative nodes (tol 1e-8) ---")
let tol = 1.0e-8
    Remp = R_iv[tol]
    order = sortperm([abs(walk.visited_z[k]) for k in 1:N])
    picks = [order[1], order[cld(N,4)], order[cld(N,2)],
             order[cld(3N,4)], order[end]]
    @printf("  %-4s %-20s %-7s %-9s %-9s %-9s %-9s %-9s\n",
            "node","z","|z|","R_i","R_ii","R_iii","R_iv","R_iv/R_ii")
    for k in picks
        z = walk.visited_z[k]
        rr = (isfinite(R_ii[k]) && isfinite(Remp[k])) ? Remp[k]/R_ii[k] : NaN
        @printf("  %-4d (%+6.2f,%+6.2f)     %-7.2f %-9.3f %-9.3f %-9.3f %-9.3f %-9.2f\n",
                k, real(z), imag(z), abs(z), R_i[k], R_ii[k], R_iii[k],
                Remp[k], rr)
    end
end

# Machine-readable dump.
open(joinpath(@__DIR__, "node_radii.csv"), "w") do io
    println(io, "node,re_z,im_z,abs_z,R_i,R_ii,R_iii,",
                "R_iv_1e6,R_iv_1e8,R_iv_1e10")
    for k in 1:N
        z = walk.visited_z[k]
        @printf(io, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                k, real(z), imag(z), abs(z), R_i[k], R_ii[k], R_iii[k],
                R_iv[1e-6][k], R_iv[1e-8][k], R_iv[1e-10][k])
    end
end
println("\n  wrote node_radii.csv")
println("[done]")
