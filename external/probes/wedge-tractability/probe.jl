# external/probes/wedge-tractability/probe.jl
#
# Phase-A2 exploration spike (bead padetaylor-0ln.37.2).  THROWAWAY CODE.
# No production file in src/ or figures/ is modified.  Run with:
#
#     julia --project=figures external/probes/wedge-tractability/probe.jl
#
# ============================================================================
# MISSION
# ============================================================================
# The P_I^(2) tritronquee headline figure renders the pole-rich wedge
# (|arg x| < 36 deg, straddling the positive real axis) by a vector
# path-network walk + Stage-2 fill.  The current V8b walk
# (figures/_kkg_pi2_helpers.jl) uses 20 targets all within |x| <= 8, fixed
# step h = 0.1.  The figure window is |x| <= 20.
#
# This spike measures whether HONEST (non-extrapolated) coverage of the wedge
# can reach |x| = 20, or where the honest frontier lands.  "Honest" = a
# shared-Q Pade is never evaluated outside its A1-verified disc of validity.
#
# It does four things (ADR-0025 Phase A2):
#   (2) extend the target fan progressively toward |x| = 20 and run the walk;
#   (3) compute each visited node's honest disc via the A1 gate (v-b) and the
#       UNION honest coverage of the wedge as a function of |x|;
#   (4) find where the fixed-h = 0.1 walk first DEGENERATES, and how the
#       reachable frontier moves with h in {0.10, 0.05, 0.03};
#   (5) a verdict feeding bead B3.
#
# ============================================================================
# THE A1 GATE (v-b) -- the per-node honest disc
# ============================================================================
# From external/probes/pade-validity-radius/REPORT.md sec 7 (the production
# gate for B1):
#
#     R_gate = min( s(tol) * h_JZ * h_v ,  0.5 * h_v * min|t*| )
#       h_JZ    = vector_step_jorba_zou(rescaled order-24 jet, tol)
#       min|t*| = min |root of the node's shared-Q|
#       s(1e-8) = 0.36   (production default)
#
# h_JZ is computed by feeding the RESCALED order-24 jet (coefficients of t^k,
# i.e. c_k * h_v^k) to the project's own VectorStepControl.vector_step_jorba_zou
# with eps_abs = tol.  The returned step is in the rescaled t-variable; multiply
# by h_v for the z-plane radius.  This is exactly the A1 winner: a truncation-
# limited per-node honest disc, with a pole-adjacency clamp.
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
using PadeTaylor.VectorStepControl:  vector_step_jorba_zou

# Companion-eigenvalue polynomial root finder -- identical to Polynomials.roots
# (which builds the companion matrix and calls eigvals).  `c` low-to-high.
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

# A1 gate (v-b) constants.
const GATE_TOL   = 1.0e-8     # production default tolerance
const GATE_S     = 0.36       # s(1e-8) per A1 REPORT sec 7
const GATE_CLAMP = 0.5        # the 0.5 * h * min|t*| pole-adjacency clamp

# Wedge geometry.  The figure window is |x| <= 20; the pole-rich wedge is
# |arg x| < 36 deg.  A3 found the wedge poles confined to |arg x| < ~29 deg.
# We measure honest coverage of the WEDGE annuli, |arg x| <= 36 deg, the
# region the figure must fill.
const WEDGE_HALF_DEG = 36.0
const WEDGE_HALF     = WEDGE_HALF_DEG * pi / 180   # radians

println("="^78)
println("Phase-A2 wedge-tractability spike  (bead padetaylor-0ln.37.2)")
println("="^78)

# ----------------------------------------------------------------------------
# 1.  BVP + seed -- the Stage-A tritronquee, reused for every walk.
# ----------------------------------------------------------------------------
println("\n[1] Solving the P_I^(2) tritronquee BVP (Stage A) ...")

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
const Y_SEED = ComplexF64.(bvp(ComplexF64(KKG_Z_SEED)))

# ----------------------------------------------------------------------------
# 2.  Target-fan construction -- progressively wider/denser fans to |x| = 20.
#
# The V8b fan is radii (2,4,6,8) x angles (-0.5..0.5 rad), 20 targets.  We
# build fans reaching radii up to 20, staying inside the wedge.  Angular
# coverage stays at +-0.5 rad (~28.6 deg) -- the V8b convention, comfortably
# inside the 36 deg wedge and bracketing the A3 pole region (~29 deg).  The
# radial spacing matches the chosen step h (a fan denser than h is wasteful;
# coarser leaves gaps the walk must bridge).
# ----------------------------------------------------------------------------
function build_fan(r_max::Float64, dr::Float64, n_ang::Int)
    radii  = collect(2.0:dr:r_max)
    angles = n_ang == 1 ? [0.0] :
             collect(range(-0.5, 0.5; length = n_ang))
    return ComplexF64[r * cis(θ) for r in radii for θ in angles]
end

# ----------------------------------------------------------------------------
# 3.  The A1 gate (v-b) -- per-node honest disc radius (z-plane).
#
# For node k: rebuild the order-24 jet (the SAME jet the canonical Pade was
# built from), rescale it (coeffs of t^k = c_k * h_v^k), feed to the project
# Jorba-Zou estimator, multiply by h_v.  Clamp by 0.5 * h_v * min|t*|.
# Returns NaN if the jet cannot be built (a degenerate, blown-up node).
# ----------------------------------------------------------------------------
function honest_disc(walk, k; tol = GATE_TOL, s = GATE_S)
    z_v = ComplexF64(walk.visited_z[k])
    y_v = ComplexF64.(walk.visited_y[k])
    h_v = real(walk.visited_h[k])
    # Rebuild the order-24 jet; rescale to the t-variable.
    local jets
    try
        jets = vector_taylor_coefficients(f, z_v, y_v, KKG_PN_ORDER)
    catch
        return NaN
    end
    rescaled = [ComplexF64[jet[q+1]*h_v^q for q in 0:length(jet)-1]
                for jet in jets]
    local h_JZ
    try
        h_JZ = vector_step_jorba_zou(rescaled, tol)
    catch
        return NaN
    end
    R_trunc = s * h_JZ * h_v
    # Pole-adjacency clamp: 0.5 * h_v * min|t*| over the node's shared-Q roots.
    Q  = walk.visited_denominator[k]
    ts = length(Q) < 2 ? ComplexF64[] : poly_roots(Q)
    R_clamp = isempty(ts) ? Inf : GATE_CLAMP * h_v * minimum(abs, ts)
    return min(R_trunc, R_clamp)
end

# ----------------------------------------------------------------------------
# 4.  Union honest coverage of a wedge annulus.
#
# For an annulus {r0 <= |x| <= r1, |arg x| <= WEDGE_HALF}, Monte-Carlo-sample
# the annulus on a fixed deterministic grid and count the fraction of sample
# points covered by at least one node's honest disc.  A point z is covered iff
# |z - z_k| <= R_k for some node k.  Deterministic grid (no RNG) so the spike
# is reproducible.
# ----------------------------------------------------------------------------
function annulus_points(r0, r1; n_r = 18, n_a = 60)
    pts = ComplexF64[]
    for ir in 1:n_r
        r = r0 + (r1 - r0) * (ir - 0.5) / n_r
        for ia in 1:n_a
            θ = -WEDGE_HALF + 2*WEDGE_HALF * (ia - 0.5) / n_a
            push!(pts, r * cis(θ))
        end
    end
    return pts
end

function coverage_fraction(pts, znodes, Rnodes)
    isempty(pts) && return 0.0
    covered = 0
    for z in pts
        hit = false
        @inbounds for k in eachindex(znodes)
            isfinite(Rnodes[k]) || continue
            if abs(z - znodes[k]) <= Rnodes[k]
                hit = true; break
            end
        end
        hit && (covered += 1)
    end
    return covered / length(pts)
end

# ----------------------------------------------------------------------------
# 5.  Run one walk: returns (ok, walk, n_nodes, runtime, msg, max_radius).
#
# `max_radius` is the largest |z| among visited nodes whose state norm did not
# blow up (||y|| finite and < BLOWUP) -- the honest reach of the walk.
# ----------------------------------------------------------------------------
const BLOWUP = 1.0e6

function run_walk(targets, h)
    prob = VectorPadeTaylorProblem(f, Y_SEED,
                                   (ComplexF64(KKG_Z_SEED),
                                    ComplexF64(20.0 + 0.0im));
                                   order = KKG_PN_ORDER)
    t0 = time()
    local walk
    try
        walk = vector_path_network_solve(prob, targets;
                                         order = KKG_PN_ORDER, h = h)
    catch err
        rt = time() - t0
        m = first(split(sprint(showerror, err), '\n'))
        return (ok = false, walk = nothing, n = 0, rt = rt,
                msg = "$(typeof(err)): $m", max_r = NaN)
    end
    rt = time() - t0
    # Honest reach: largest |z| of a non-degenerate node.
    max_r = 0.0
    for k in 1:length(walk.visited_z)
        ny = norm(walk.visited_y[k])
        (isfinite(ny) && ny < BLOWUP) || continue
        r = abs(walk.visited_z[k])
        r > max_r && (max_r = r)
    end
    return (ok = true, walk = walk, n = length(walk.visited_z),
            rt = rt, msg = "OK", max_r = max_r)
end

# ----------------------------------------------------------------------------
# 6.  TASK 2+3 -- extended fans, honest coverage vs radius.
# ----------------------------------------------------------------------------
println("\n[2] Extended target fans -- walks toward |x| = 20 (h = 0.1) ...")

# Progressive fans: each reaches a larger r_max.  dr = 2.0 radial spacing,
# 5 angles (the V8b convention).  We test r_max in {8, 12, 16, 20}.
const FAN_RMAXS = (8.0, 12.0, 16.0, 20.0)
const FAN_DR    = 2.0
const FAN_NANG  = 5

fan_results = Dict{Float64,Any}()
for rmax in FAN_RMAXS
    fan = build_fan(rmax, FAN_DR, FAN_NANG)
    res = run_walk(fan, 0.1)
    fan_results[rmax] = res
    @printf("    fan r_max=%4.1f (%3d targets): %-7s nodes=%4d  runtime=%6.2fs  honest-reach |z|=%.2f\n",
            rmax, length(fan), res.ok ? "OK" : "BLOCKED",
            res.n, res.rt, res.max_r)
    res.ok || @printf("        -> %s\n", res.msg)
end

# Coverage-vs-radius for the largest fan that completed (or the largest fan).
println("\n[3] Honest coverage vs |x| (A1 gate v-b, tol=1e-8, s=0.36) ...")

# Pick the walk with the most honest reach.
best_rmax = nothing
for rmax in FAN_RMAXS
    global best_rmax
    r = fan_results[rmax]
    if r.ok && (best_rmax === nothing || r.n >= fan_results[best_rmax].n)
        best_rmax = rmax
    end
end

if best_rmax === nothing
    println("    NO fan completed -- see the fixed-h failure-mode probe below.")
else
    walk = fan_results[best_rmax].walk
    N = length(walk.visited_z)
    @printf("    using fan r_max=%.1f: %d visited nodes\n", best_rmax, N)
    println("    computing per-node honest discs ...")
    znodes = ComplexF64[ComplexF64(walk.visited_z[k]) for k in 1:N]
    Rnodes = Float64[honest_disc(walk, k) for k in 1:N]
    n_nan  = count(!isfinite, Rnodes)
    finite_R = filter(isfinite, Rnodes)
    @printf("    honest disc radii: n=%d finite (%d degenerate), ", length(finite_R), n_nan)
    if !isempty(finite_R)
        sr = sort(finite_R)
        @printf("med=%.4f p90=%.4f max=%.4f\n",
                sr[cld(length(sr),2)], sr[cld(9*length(sr),10)], sr[end])
    else
        println()
    end

    # Coverage of each wedge annulus [r, r+2].
    println("\n    --- Honest coverage of wedge annuli (|arg x| <= 36 deg) ---")
    @printf("    %-12s %-10s %-12s %-10s\n",
            "annulus", "coverage", "nodes-in", "nodes-touch")
    global_cov = Float64[]
    for r0 in 0.0:2.0:18.0
        r1 = r0 + 2.0
        pts = annulus_points(r0, r1)
        cov = coverage_fraction(pts, znodes, Rnodes)
        push!(global_cov, cov)
        # nodes geometrically inside the annulus wedge
        n_in = count(z -> r0 <= abs(z) <= r1 && abs(angle(z)) <= WEDGE_HALF,
                     znodes)
        # nodes whose honest disc intersects the annulus wedge at all
        n_touch = 0
        for k in 1:N
            isfinite(Rnodes[k]) || continue
            zk = znodes[k]; ar = abs(zk)
            if ar + Rnodes[k] >= r0 && ar - Rnodes[k] <= r1
                n_touch += 1
            end
        end
        @printf("    [%4.1f,%4.1f]   %6.1f%%     %-12d %-10d\n",
                r0, r1, 100*cov, n_in, n_touch)
    end

    # Honest frontier: the largest r such that the annulus [r-2,r] still has
    # >= COVER_THRESH coverage.
    cover_thresh = 0.50
    frontier = 0.0
    for (i, cov) in enumerate(global_cov)
        r1 = 2.0*i
        cov >= cover_thresh && (frontier = r1)
    end
    @printf("\n    Honest frontier (>=%.0f%% annulus coverage): |x| ~ %.1f\n",
            100*cover_thresh, frontier)
end

# ----------------------------------------------------------------------------
# 7.  TASK 4 -- fixed-h failure-mode probe (manual wedge ray-walk).
#
# `vector_path_network_solve` aborts the WHOLE solve on a single throw (it does
# not return partial results), so a fan reaching past the death radius reports
# only "BLOCKED, 0 nodes".  To find WHERE the fixed-h walk first degenerates we
# replicate the production wedge-walk step-by-step (the min-||y|| 5-direction
# wedge of VectorPathNetwork._wedge_step) and march from the seed toward a far
# wedge target, catching the throw.  The degeneration radius is the |x| of the
# last node reached before the throw / before ||y|| > BLOWUP.
#
# Wedge angles: DEFAULT_WEDGE = [-pi/4,-pi/8,0,pi/8,pi/4] (FW 2011 sec 3.1) --
# verbatim from VectorPathNetwork.
# ----------------------------------------------------------------------------
println("\n[4] Fixed-h failure-mode probe (manual wedge ray-walk) ...")

const WEDGE = [-pi/4, -pi/8, 0.0, pi/8, pi/4]

# One min-||y|| wedge step toward `goal_dir` from (z,y).  Returns
# (z_new, y_new, ok).  ok = false if all 5 candidates throw.
function wedge_step(z, y, h, goal_dir)
    best_z = ComplexF64(0); best_y = ComplexF64[]; best_val = Inf
    for θoff in WEDGE
        θ = goal_dir + θoff
        h_step = ComplexF64(h*cos(θ), h*sin(θ))
        st = VectorPadeStepperState{ComplexF64}(ComplexF64(z),
                                                ComplexF64.(y))
        try
            vector_pade_step_with_pade!(st, f, KKG_PN_ORDER, h_step)
            v = norm(st.y)
            if isfinite(v) && v < best_val
                best_val = v; best_z = st.z; best_y = st.y
            end
        catch
        end
    end
    return best_val == Inf ? (ComplexF64(0), ComplexF64[], false) :
                             (best_z, best_y, true)
end

# March from the seed toward `target` via the wedge walk at fixed h.
# Returns (reach_r, n_steps, cause): reach_r is |z| of the last good node;
# cause is :reached / :throw / :blowup.
function ray_walk(target, h; max_steps = 4000)
    z = ComplexF64(KKG_Z_SEED); y = copy(Y_SEED)
    reach_r = abs(z); n = 0
    while abs(z - target) > h && n < max_steps
        n += 1
        goal = angle(target - z)
        z_new, y_new, ok = wedge_step(z, y, h, goal)
        if !ok
            return (reach_r, n, :throw)
        end
        ny = norm(y_new)
        if !isfinite(ny) || ny >= BLOWUP
            return (reach_r, n, :blowup)
        end
        z, y = z_new, y_new
        reach_r = abs(z)
    end
    return (reach_r, n, n >= max_steps ? :maxsteps : :reached)
end

const STEP_SIZES = (0.10, 0.05, 0.03)
# March toward |x|=20 on the positive real axis (arg x = 0) -- the densest
# part of the pole field, the hardest ray; and toward 20*cis(0.4) (upper
# wedge).  The min over targets is the honest fixed-h frontier.
const PROBE_TARGETS = (ComplexF64(20.0, 0.0),
                       ComplexF64(20.0*cos(0.4), 20.0*sin(0.4)),
                       ComplexF64(20.0*cos(-0.4), 20.0*sin(-0.4)))

@printf("    %-8s %-26s %-12s %-9s %-10s\n",
        "h", "target", "reach |x|", "steps", "cause")
fixedh_frontier = Dict{Float64,Float64}()
for h in STEP_SIZES
    worst = Inf
    for tgt in PROBE_TARGETS
        t0 = time()
        reach_r, n, cause = ray_walk(tgt, h)
        rt = time() - t0
        @printf("    %-8.2f (%+5.1f,%+5.1f)  ->  20    %-12.2f %-9d %s  (%.1fs)\n",
                h, real(tgt), imag(tgt), reach_r, n, cause, rt)
        worst = min(worst, reach_r)
    end
    fixedh_frontier[h] = worst
    @printf("        => fixed-h=%.2f honest reach (min over rays): |x| = %.2f\n",
            h, worst)
end

# ----------------------------------------------------------------------------
# 8.  Incremental-fan honest frontier -- the largest fan that COMPLETES,
#     and the union honest coverage it achieves.
#
# Because the full solve aborts on a throw, we find the largest r_max for
# which the production walk completes by bisection-style scan, then report
# the union honest coverage of that completed walk.  Tested per step size.
# ----------------------------------------------------------------------------
println("\n[5] Largest completing fan + honest frontier, per step size ...")
@printf("    %-8s %-12s %-9s %-12s %-12s %-9s\n",
        "h", "max r_max", "nodes", "reach |z|", "frontier", "runtime")

cover_thresh2 = 0.50
for h in STEP_SIZES
    # Scan r_max downward from 20 to 8 in steps of 2; take the largest that
    # completes.
    best = nothing
    for rmax in (20.0, 18.0, 16.0, 14.0, 12.0, 10.0, 8.0)
        fan = build_fan(rmax, FAN_DR, FAN_NANG)
        res = run_walk(fan, h)
        if res.ok
            best = (rmax = rmax, res = res)
            break
        end
    end
    if best === nothing
        @printf("    %-8.2f no fan >= r_max 8 completed\n", h)
        continue
    end
    walk = best.res.walk
    N = length(walk.visited_z)
    znodes = ComplexF64[ComplexF64(walk.visited_z[k]) for k in 1:N]
    Rnodes = Float64[honest_disc(walk, k) for k in 1:N]
    frontier = 0.0
    for i in 1:10
        r0 = 2.0*(i-1); r1 = 2.0*i
        cov = coverage_fraction(annulus_points(r0, r1), znodes, Rnodes)
        cov >= cover_thresh2 && (frontier = r1)
    end
    @printf("    %-8.2f %-12.1f %-9d %-12.2f %-12.1f %.2fs\n",
            h, best.rmax, N, best.res.max_r, frontier, best.res.rt)
end

# ----------------------------------------------------------------------------
# 9.  TASK 5 -- dense-fan coverage scaling + tiling estimate.
#
# Section [3] showed the V8b-style fan (5 angles, dr=2) reaches only ~7%
# honest coverage: it visits a sparse RAY NETWORK, not a 2D fill.  The honest
# discs (med 0.062) are tiny, so coverage needs MANY more nodes.  Here we run
# a much DENSER fan -- angles and radial spacing matched to the honest-disc
# scale -- and measure the achievable coverage vs |x|, plus the node count.
#
# Tiling estimate: with median honest disc area pi*R_med^2 and a packing
# efficiency eta ~ 0.6 (discs of one size, hex-ish packing of a 2D region),
# the node count to tile a wedge annulus of area A is ~ A / (eta*pi*R_med^2).
# ----------------------------------------------------------------------------
println("\n[6] Dense-fan coverage scaling + tiling estimate ...")

# Coverage scaling: try fans of increasing target density and report which
# COMPLETE, and the honest coverage each achieves.  Denser fans add inter-
# target bridging walks that cross the pole field -> more pole-hit risk, so
# the densest fans may BLOCK.  We report the trade-off explicitly.
fan_specs = [
    ("5ang dr2.0",  collect(2.0:2.0:20.0), collect(range(-0.5,0.5;length=5))),
    ("9ang dr1.0",  collect(2.0:1.0:20.0), collect(range(-0.5,0.5;length=9))),
    ("13ang dr0.5", collect(2.0:0.5:20.0), collect(range(-0.5,0.5;length=13))),
]
for h in (0.05, 0.03)
    @printf("\n    --- coverage scaling at step h=%.2f ---\n", h)
    for (label, radii, angles) in fan_specs
        fan = ComplexF64[r*cis(θ) for r in radii for θ in angles]
        res = run_walk(fan, h)
        if !res.ok
            @printf("    %-14s (%3d tgt): BLOCKED\n", label, length(fan))
            continue
        end
        walk = res.walk
        N = length(walk.visited_z)
        znodes = ComplexF64[ComplexF64(walk.visited_z[k]) for k in 1:N]
        Rnodes = Float64[honest_disc(walk, k) for k in 1:N]
        # overall wedge coverage to |x|=20
        cov20 = coverage_fraction(annulus_points(0.0, 20.0; n_r=80, n_a=120),
                                  znodes, Rnodes)
        frontier = 0.0
        for i in 1:10
            r0 = 2.0*(i-1); r1 = 2.0*i
            c = coverage_fraction(annulus_points(r0, r1), znodes, Rnodes)
            c >= 0.50 && (frontier = r1)
        end
        @printf("    %-14s (%3d tgt): %5d nodes, cov(|x|<=20)=%5.1f%%, frontier=%.0f, %.1fs\n",
                label, length(fan), N, 100*cov20, frontier, res.rt)
    end
end

let h = 0.03
    dense = ComplexF64[r*cis(θ) for r in collect(2.0:1.0:20.0)
                                 for θ in collect(range(-0.5,0.5;length=9))]
    @printf("\n    detailed annulus profile -- 9ang dr1.0 fan, h=%.2f ...\n", h)
    res = run_walk(dense, h)
    if !res.ok
        @printf("    fan BLOCKED: %s\n", res.msg)
    else
        walk = res.walk
        N = length(walk.visited_z)
        znodes = ComplexF64[ComplexF64(walk.visited_z[k]) for k in 1:N]
        Rnodes = Float64[honest_disc(walk, k) for k in 1:N]
        finite_R = filter(isfinite, Rnodes)
        R_med = isempty(finite_R) ? NaN :
                sort(finite_R)[cld(length(finite_R),2)]
        @printf("    walk: %d nodes, runtime %.2fs, honest disc R_med=%.4f\n",
                N, res.rt, R_med)
        println("\n    --- Dense-fan honest coverage of wedge annuli ---")
        @printf("    %-12s %-10s %-12s %-14s\n",
                "annulus", "coverage", "nodes-in", "tile-need")
        for i in 1:10
            r0 = 2.0*(i-1); r1 = 2.0*i
            cov = coverage_fraction(annulus_points(r0, r1), znodes, Rnodes)
            n_in = count(z -> r0 <= abs(z) <= r1 &&
                              abs(angle(z)) <= WEDGE_HALF, znodes)
            # wedge-annulus area = (2*WEDGE_HALF) * (r1^2 - r0^2)/2
            A = WEDGE_HALF * (r1^2 - r0^2)
            tile_need = isnan(R_med) ? NaN : A / (0.6*pi*R_med^2)
            @printf("    [%4.1f,%4.1f]   %6.1f%%     %-12d %-14.0f\n",
                    r0, r1, 100*cov, n_in, tile_need)
        end
        # total tiling estimate for the whole wedge to |x|=20
        A_total = WEDGE_HALF * 20.0^2
        tile_total = isnan(R_med) ? NaN : A_total / (0.6*pi*R_med^2)
        @printf("\n    Whole-wedge (|x|<=20) area = %.1f ; honest disc area ", A_total)
        @printf("= %.4f\n", pi*R_med^2)
        @printf("    Tiling estimate: ~%.0f honest nodes to fill |x|<=20 wedge\n",
                tile_total)
        # frontier
        frontier = 0.0
        for i in 1:10
            r0 = 2.0*(i-1); r1 = 2.0*i
            cov = coverage_fraction(annulus_points(r0, r1), znodes, Rnodes)
            cov >= 0.50 && (frontier = r1)
        end
        @printf("    Dense-fan honest frontier (>=50%% coverage): |x| ~ %.1f\n",
                frontier)
    end
end

# ============================================================================
# 10.  PHASE-A2b -- DISPLAY-TOLERANCE recompute (the follow-up).
#
# The A1 gate (v-b) ran at SOLVER tolerances (1e-6/1e-8/1e-10); the 8.6%
# headline honest-coverage figure used tol=1e-8.  But the headline FIGURE is
# a DISPLAY artifact: a clamped Re/Im heatmap + a 3D surface.  A display
# colour channel resolves only ~1e-2..1e-3 relative accuracy -- you cannot
# SEE better than that on a clamped colormap.  So rendering the wedge SURFACE
# honest-to-1e-3 is legitimate ("honest for display"); the pole LOCATIONS --
# the scientific content -- stay separately extracted + validated to 1e-8+.
#
# The honest A1 disc grows as tol relaxes: h_JZ = vector_step_jorba_zou(jet,
# eps) grows with eps, and so does s(tol).  This section recomputes the
# honest UNION coverage at DISPLAY tolerances tol in {1e-2, 1e-3, 1e-4}.
#
# A1 calibrated s(tol) only for 1e-6/1e-8/1e-10 (s = 0.34/0.36/0.30 -- the
# ~1-2% over-rate band).  There is no calibrated s at display tol.  We hold
# s = 0.36 fixed (the production default, the most conservative interior
# value of the calibrated band) so the dominant tol-dependence is carried by
# h_JZ itself, which the Jorba-Zou estimator reads directly from the per-node
# coefficient decay.  Using a FIXED s is the conservative choice: the true
# display-tol s would, if anything, be a touch larger (the over-rate band
# loosens as tol relaxes), so the coverage reported here is a LOWER bound.
# ============================================================================
println("\n" * "="^78)
println("[A2b] DISPLAY-TOLERANCE honest-coverage recompute (the follow-up)")
println("="^78)

const DISPLAY_TOLS = (1.0e-2, 1.0e-3, 1.0e-4)

# honest_disc with an explicit Taylor order (24 production / 48 high-order)
# and an explicit tolerance.  Same A1 gate (v-b): R = min(s*h_JZ*h, 0.5*h*min|t*|).
function honest_disc_ord(walk, k, order; tol = GATE_TOL, s = GATE_S)
    z_v = ComplexF64(walk.visited_z[k])
    y_v = ComplexF64.(walk.visited_y[k])
    h_v = real(walk.visited_h[k])
    local jets
    try
        jets = vector_taylor_coefficients(f, z_v, y_v, order)
    catch
        return NaN
    end
    rescaled = [ComplexF64[jet[q+1]*h_v^q for q in 0:length(jet)-1]
                for jet in jets]
    local h_JZ
    try
        h_JZ = vector_step_jorba_zou(rescaled, tol)
    catch
        return NaN
    end
    R_trunc = s * h_JZ * h_v
    Q  = walk.visited_denominator[k]
    ts = length(Q) < 2 ? ComplexF64[] : poly_roots(Q)
    R_clamp = isempty(ts) ? Inf : GATE_CLAMP * h_v * minimum(abs, ts)
    return min(R_trunc, R_clamp)
end

# Helper: median / p10 / p90 of the finite radii.
function radii_stats(R)
    fr = sort(filter(isfinite, R))
    isempty(fr) && return (med = NaN, p10 = NaN, p90 = NaN, nfin = 0)
    return (med = fr[cld(length(fr),2)],
            p10 = fr[cld(length(fr),10)],
            p90 = fr[cld(9*length(fr),10)],
            nfin = length(fr))
end

# The "best completing fan" for the display-tol study: the A2 §3.3 winner,
# the 9-angle dr=1.0 fan at h=0.03 -- the densest fan that COMPLETES.
println("\n[A2b-1] Rebuilding the best completing fan (9ang dr1.0, h=0.03) ...")
const A2B_FAN = ComplexF64[r*cis(θ) for r in collect(2.0:1.0:20.0)
                                     for θ in collect(range(-0.5,0.5;length=9))]
a2b_res = run_walk(A2B_FAN, 0.03)
if !a2b_res.ok
    println("    fan BLOCKED: ", a2b_res.msg, "  -- cannot run A2b study")
else
    walk = a2b_res.walk
    N = length(walk.visited_z)
    znodes = ComplexF64[ComplexF64(walk.visited_z[k]) for k in 1:N]
    @printf("    walk OK: %d nodes, runtime %.1fs, honest-reach |z|=%.2f\n",
            N, a2b_res.rt, a2b_res.max_r)

    # --- TASK 1+2: honest union coverage + R_med at each display tol ----------
    println("\n[A2b-2] Honest UNION coverage of |x|<=20 wedge, ORDER 24, vs tol")
    @printf("    %-10s %-10s %-12s %-10s %-10s %-9s\n",
            "tol", "R_med", "R_p10/p90", "cov<=20", "frontier", "tile-N")
    # full-wedge coverage grid (dense)
    full_pts = annulus_points(0.0, 20.0; n_r = 80, n_a = 120)
    a2b_o24 = Dict{Float64,Any}()
    # baseline: solver tol 1e-8 for the side-by-side
    for tol in (1.0e-8, DISPLAY_TOLS...)
        R = Float64[honest_disc_ord(walk, k, 24; tol = tol, s = GATE_S)
                    for k in 1:N]
        st = radii_stats(R)
        cov = coverage_fraction(full_pts, znodes, R)
        frontier = 0.0
        for i in 1:10
            r0 = 2.0*(i-1); r1 = 2.0*i
            c = coverage_fraction(annulus_points(r0, r1), znodes, R)
            c >= 0.50 && (frontier = r1)
        end
        A_total  = WEDGE_HALF * 20.0^2
        tileN    = isnan(st.med) ? NaN : A_total / (0.6*pi*st.med^2)
        a2b_o24[tol] = (R = R, st = st, cov = cov, frontier = frontier)
        @printf("    %-10.0e %-10.4f %.4f/%.4f  %7.1f%%   %-10.1f %.0f\n",
                tol, st.med, st.p10, st.p90, 100*cov, frontier, tileN)
    end

    # --- TASK 4: coverage-vs-radius at tol=1e-3, order 24 ---------------------
    println("\n[A2b-3] Coverage-vs-radius at tol=1e-3 (ORDER 24) -- does the")
    println("        honest surface reach |x|=20 ?")
    R3 = a2b_o24[1.0e-3].R
    @printf("    %-12s %-10s %-12s\n", "annulus", "coverage", "nodes-in")
    for i in 1:10
        r0 = 2.0*(i-1); r1 = 2.0*i
        cov = coverage_fraction(annulus_points(r0, r1), znodes, R3)
        n_in = count(z -> r0 <= abs(z) <= r1 &&
                          abs(angle(z)) <= WEDGE_HALF, znodes)
        @printf("    [%4.1f,%4.1f]   %6.1f%%     %-12d\n", r0, r1, 100*cov, n_in)
    end

    # --- TASK 3: node count to tile the wedge at tol=1e-3, order 24 ----------
    let st3 = a2b_o24[1.0e-3].st
        A_total = WEDGE_HALF * 20.0^2
        @printf("\n[A2b-4] Tile-need at tol=1e-3 (order 24): R_med=%.4f, disc area",
                st3.med)
        @printf(" %.5f\n", pi*st3.med^2)
        @printf("        whole-wedge area=%.1f  ->  ~%.0f honest nodes to tile\n",
                A_total, A_total/(0.6*pi*st3.med^2))
        @printf("        (this completing fan has %d nodes -> %.1fx short)\n",
                N, (A_total/(0.6*pi*st3.med^2))/N)
    end

    # --- TASK 5: raise Taylor order 24 -> 48 (one extra walk) ----------------
    # The A1 gate disc is TRUNCATION-limited; the Jorba-Zou estimator reads
    # the coefficient decay of the supplied jet.  Feed it the order-48 jet:
    # more terms => the decay is sampled further => h_JZ (and R) grow.  We do
    # ONE extra walk at order 48 so the shared-Q denominators and visited
    # states are themselves order-48 consistent (the canonical Pade the node
    # carries is order-48), then recompute discs from order-48 jets.
    println("\n[A2b-5] Raising Taylor order 24 -> 48 (one extra walk) ...")
    const KKG_ORDER_HI = 48
    function run_walk_ord(targets, h, order)
        prob = VectorPadeTaylorProblem(f, Y_SEED,
                                       (ComplexF64(KKG_Z_SEED),
                                        ComplexF64(20.0 + 0.0im));
                                       order = order)
        t0 = time()
        local w
        try
            w = vector_path_network_solve(prob, targets; order = order, h = h)
        catch err
            return (ok = false, walk = nothing, n = 0, rt = time()-t0,
                    msg = first(split(sprint(showerror, err), '\n')))
        end
        return (ok = true, walk = w, n = length(w.visited_z),
                rt = time()-t0, msg = "OK")
    end
    res48 = run_walk_ord(A2B_FAN, 0.03, KKG_ORDER_HI)
    if !res48.ok
        @printf("    order-48 WALK BLOCKED: %s\n",
                first(split(res48.msg, ';')))
        println("    -> the order-48 shared-Q Pade linear solve degenerates:")
        println("       m = order/2 = 24 deg num/den; the order-48 jets give")
        println("       an all-below-tol singular spectrum.  An order-48 *walk*")
        println("       is NOT viable for this wedge under the shared-Q route.")
        println()
        println("    However the A1 honest DISC is a property of the JET, not")
        println("    the shared-Q Pade: vector_step_jorba_zou reads the order-N")
        println("    coefficient decay directly.  So we can still measure what")
        println("    order 48 does to R_med / coverage by rebuilding order-48")
        println("    JETS at the order-24 walk's (already honest) node states")
        println("    and re-gating -- the truncation bound, decoupled from the")
        println("    (failed) order-48 Pade construction.")
        println()
        println("    --- ORDER-48 JETS on the order-24 walk nodes: coverage vs tol ---")
        @printf("    %-10s %-10s %-12s %-10s %-10s %-9s\n",
                "tol", "R_med", "R_p10/p90", "cov<=20", "frontier", "tile-N")
        a2b_o48j = Dict{Float64,Any}()
        for tol in (1.0e-8, DISPLAY_TOLS...)
            R = Float64[honest_disc_ord(walk, k, KKG_ORDER_HI;
                                        tol = tol, s = GATE_S) for k in 1:N]
            st = radii_stats(R)
            cov = coverage_fraction(full_pts, znodes, R)
            frontier = 0.0
            for i in 1:10
                r0 = 2.0*(i-1); r1 = 2.0*i
                c = coverage_fraction(annulus_points(r0, r1), znodes, R)
                c >= 0.50 && (frontier = r1)
            end
            A_total = WEDGE_HALF * 20.0^2
            tileN   = isnan(st.med) ? NaN : A_total/(0.6*pi*st.med^2)
            a2b_o48j[tol] = (R = R, st = st, cov = cov, frontier = frontier)
            @printf("    %-10.0e %-10.4f %.4f/%.4f  %7.1f%%   %-10.1f %.0f\n",
                    tol, st.med, st.p10, st.p90, 100*cov, frontier, tileN)
        end
        println("\n    --- ORDER-48 JETS: coverage-vs-radius at tol=1e-3 ---")
        R48j_3 = a2b_o48j[1.0e-3].R
        @printf("    %-12s %-10s %-12s\n", "annulus", "coverage", "nodes-in")
        for i in 1:10
            r0 = 2.0*(i-1); r1 = 2.0*i
            cov = coverage_fraction(annulus_points(r0, r1), znodes, R48j_3)
            n_in = count(z -> r0 <= abs(z) <= r1 &&
                              abs(angle(z)) <= WEDGE_HALF, znodes)
            @printf("    [%4.1f,%4.1f]   %6.1f%%     %-12d\n",
                    r0, r1, 100*cov, n_in)
        end
        c24 = 100*a2b_o24[1.0e-3].cov
        c48 = 100*a2b_o48j[1.0e-3].cov
        @printf("\n    VERDICT INPUT  tol=1e-3:  order24 cov=%.1f%%  ", c24)
        @printf("order48-jet cov=%.1f%%  (>=70%% target)\n", c48)
        @printf("    order24 R_med=%.4f -> order48-jet R_med=%.4f  (%.2fx)\n",
                a2b_o24[1.0e-3].st.med, a2b_o48j[1.0e-3].st.med,
                a2b_o48j[1.0e-3].st.med / a2b_o24[1.0e-3].st.med)
        @printf("    NOTE: a real order-48 walk is BLOCKED -- order 48 cannot")
        @printf(" enlarge the\n          shipped figure's discs unless the")
        @printf(" shared-Q route is replaced.\n")
    end
    if res48.ok
        w48 = res48.walk
        N48 = length(w48.visited_z)
        zn48 = ComplexF64[ComplexF64(w48.visited_z[k]) for k in 1:N48]
        @printf("    order-48 walk OK: %d nodes, runtime %.1fs (order-24 was %.1fs)\n",
                N48, res48.rt, a2b_res.rt)
        println("\n    --- ORDER 48: honest union coverage vs tol ---")
        @printf("    %-10s %-10s %-12s %-10s %-10s %-9s\n",
                "tol", "R_med", "R_p10/p90", "cov<=20", "frontier", "tile-N")
        a2b_o48 = Dict{Float64,Any}()
        for tol in (1.0e-8, DISPLAY_TOLS...)
            R = Float64[honest_disc_ord(w48, k, KKG_ORDER_HI;
                                        tol = tol, s = GATE_S) for k in 1:N48]
            st = radii_stats(R)
            cov = coverage_fraction(full_pts, zn48, R)
            frontier = 0.0
            for i in 1:10
                r0 = 2.0*(i-1); r1 = 2.0*i
                c = coverage_fraction(annulus_points(r0, r1), zn48, R)
                c >= 0.50 && (frontier = r1)
            end
            A_total = WEDGE_HALF * 20.0^2
            tileN   = isnan(st.med) ? NaN : A_total/(0.6*pi*st.med^2)
            a2b_o48[tol] = (R = R, st = st, cov = cov, frontier = frontier)
            @printf("    %-10.0e %-10.4f %.4f/%.4f  %7.1f%%   %-10.1f %.0f\n",
                    tol, st.med, st.p10, st.p90, 100*cov, frontier, tileN)
        end
        # order-48 coverage-vs-radius at tol=1e-3
        println("\n    --- ORDER 48: coverage-vs-radius at tol=1e-3 ---")
        R48_3 = a2b_o48[1.0e-3].R
        @printf("    %-12s %-10s %-12s\n", "annulus", "coverage", "nodes-in")
        for i in 1:10
            r0 = 2.0*(i-1); r1 = 2.0*i
            cov = coverage_fraction(annulus_points(r0, r1), zn48, R48_3)
            n_in = count(z -> r0 <= abs(z) <= r1 &&
                              abs(angle(z)) <= WEDGE_HALF, zn48)
            @printf("    [%4.1f,%4.1f]   %6.1f%%     %-12d\n",
                    r0, r1, 100*cov, n_in)
        end
        # side-by-side verdict line
        c24 = 100*a2b_o24[1.0e-3].cov
        c48 = 100*a2b_o48[1.0e-3].cov
        @printf("\n    VERDICT INPUT  tol=1e-3:  order24 cov=%.1f%%  ", c24)
        @printf("order48 cov=%.1f%%  (>=70%% target)\n", c48)
        @printf("    order24 R_med=%.4f -> order48 R_med=%.4f  (%.2fx)\n",
                a2b_o24[1.0e-3].st.med, a2b_o48[1.0e-3].st.med,
                a2b_o48[1.0e-3].st.med / a2b_o24[1.0e-3].st.med)
    end
end

println("\n[done]  See REPORT.md sec 7 for the verdict.")
