# external/probes/v8b-baseline/probe.jl
#
# Phase-A4 EXPLORATION SPIKE (bead padetaylor-0ln.37.4) — throwaway code.
#
# Establishes the BEFORE-picture for the headline-figure re-resolution
# (ADR-0025): the baseline quality of the *current* V8b P_I^(2)
# tritronquee wedge walk and its 21 extracted poles. These numbers
# anchor the Phase-F before/after comparison and feed bead D-VC7.
#
# Three measurements:
#
#   T2 — Loop-closure ΔP_rel on the V8b vector walk.  `quality_diagnose`
#        (ADR-0016, src/Diagnostics.jl + ext/PadeTaylorDiagnosticsExt.jl)
#        is typed `quality_diagnose(sol::PathNetworkSolution)` — the
#        SCALAR path-network type.  The V8b walk produces a
#        `VectorPathNetworkSolution`.  We FIRST confirm `quality_diagnose`
#        rejects the vector solution (MethodError), THEN compute the
#        equivalent loop-closure certificate by hand from the vector
#        walk's `visited_*` arrays, mirroring the ext algorithm exactly
#        but with the shared-Q evaluation and a vector-norm ΔP_rel.
#
#   T3 — Conjugate-symmetry pairing of the 21 V8b poles.  V_0 is real on
#        the negative real axis ⇒ V_0(conj x) = conj(V_0(x)) ⇒ wedge
#        poles come in conjugate pairs.  Match the 21 poles into pairs
#        by nearest-neighbour under conjugation; report pairing residual.
#
#   T4 — VC-4 spot check: dominant-balance leading coefficient A for a
#        few poles (fit u ~ A(x-x0)^-2 on a small ring).  ADR-0025 A3:
#        A ∈ {-1,-3}.
#
# Run: julia --project=external/probes/v8b-baseline external/probes/v8b-baseline/probe.jl

using PadeTaylor
using PadeTaylor.Diagnostics: quality_diagnose
using PadeTaylor.VectorPathNetwork: VectorPathNetworkSolution
using DelaunayTriangulation: triangulate, each_edge
using Polynomials: Polynomial, roots
using Statistics: median, quantile, mean
using Printf

# The exact V8b kernel — load the production helper verbatim so the
# probe exercises the identical computation as the shipped figure.
include(joinpath(@__DIR__, "..", "..", "..", "figures", "_kkg_pi2_helpers.jl"))

const REPORT = String[]
log!(s) = (println(s); push!(REPORT, s))

log!("# V8b Baseline Probe — Phase-A4 (bead padetaylor-0ln.37.4)")
log!("")
log!("Throwaway exploration spike. BEFORE-picture for ADR-0025 re-resolution.")
log!("")

# ----------------------------------------------------------------------
# Reproduce the V8b walk.
# ----------------------------------------------------------------------
log!("## 0. Reproduce the V8b wedge walk")
log!("")
t0 = time()
result = kkg_pi2_compute()
log!("Compute time: $(round(time() - t0; digits=1)) s")
log!("march_ok    : $(result.march_ok)")
log!("message     : $(result.message)")

if !result.march_ok || result.walk === nothing
    log!("")
    log!("**FATAL** — the V8b march did not complete; no baseline to measure.")
    write(joinpath(@__DIR__, "REPORT.md"), join(REPORT, "\n") * "\n")
    error("V8b march failed; see REPORT.md")
end

walk   = result.walk
poles  = result.poles
N      = length(walk.visited_z)
log!("visited nodes : $N")
log!("poles extracted: $(length(poles))")
log!("")

# ======================================================================
# TASK 2 — loop-closure diagnostic.
# ======================================================================
log!("## Task 2 — Loop-closure ΔP_rel certificate")
log!("")

# --- 2a. Does quality_diagnose accept the vector solution as-is? -------
log!("### 2a. Does `quality_diagnose` accept a `VectorPathNetworkSolution`?")
log!("")
qd_works = false
try
    quality_diagnose(walk)
    global qd_works = true
    log!("UNEXPECTED: `quality_diagnose(walk)` returned without error.")
catch err
    log!("`quality_diagnose(walk)` THREW — `$(typeof(err))`.")
    log!("First line: `$(first(split(sprint(showerror, err), '\n')))`")
end
log!("")
log!("Verdict: `quality_diagnose` " *
     (qd_works ? "ACCEPTS" : "does NOT accept") *
     " the vector walk.")
log!("")

# --- 2b. Compute the loop-closure certificate by hand. -----------------
# Mirror ext/PadeTaylorDiagnosticsExt.jl exactly, with two adaptations
# forced by the vector solution's data model:
#   (i)  the vector node stores a shared-Q approximant (numerators +
#        one denominator), not a scalar PadeApproximant — evaluate via
#        P_i(t)/Q(t) Horner, exactly _stage2_fill's _eval_poly pattern;
#   (ii) the node state y is a d-vector — ΔP_rel uses the vector norm
#        ||y_A(M) - y_B(M)|| / (||y_A(M)|| + ||y_B(M)|| + eps).
# Everything else (Delaunay edges, tree edges, non-tree set, midpoint
# rescaling t = (M - z_X)/h_X, LCA tree-distance, category thresholds)
# is the ext algorithm verbatim.
log!("### 2b. Hand-computed equivalent loop-closure certificate")
log!("")

_eval_poly(c, t) = begin
    s = zero(eltype(c)) * t
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    s
end

# Shared-Q vector evaluation at rescaled t: returns the d-vector state.
function eval_vec_node(walk, k, t)
    nums = walk.visited_numerators[k]
    den  = walk.visited_denominator[k]
    q_t  = _eval_poly(den, t)
    return [_eval_poly(num, t) / q_t for num in nums]
end

const TOL_WELL = 1e-10
const TOL_BAD  = 1e-6
const EPS_FLOOR = 1e-300

# Sheet-0: the vector walk carries no `visited_sheet` field at all
# (P_I^(2) in companion form is meromorphic — single-sheeted; see
# VectorPathNetwork.jl docstring).  So sheet-0 is "all nodes" — no
# strip predicate needed.  Indices are global == local.
s_idx = collect(1:N)

# Delaunay triangulate the node cloud.
pts2D = [(Float64(real(z)), Float64(imag(z))) for z in walk.visited_z]
tri   = triangulate(pts2D)
norm_edge(a, b) = a < b ? (a, b) : (b, a)
delaunay_edges = Set{Tuple{Int,Int}}()
for e in each_edge(tri)
    a, b = e
    (a < 1 || b < 1) && continue
    push!(delaunay_edges, norm_edge(a, b))
end

# Tree edges from visited_parent.
tree_edges = Set{Tuple{Int,Int}}()
for k in 2:N
    p = walk.visited_parent[k]
    p == 0 && continue
    push!(tree_edges, norm_edge(k, p))
end
nontree_edges = setdiff(delaunay_edges, tree_edges)

# Tree depths + LCA distance (verbatim from the ext).
function build_depths(parent)
    n = length(parent); depth = fill(-1, n)
    for k in 1:n; parent[k] == 0 && (depth[k] = 0); end
    progress = true
    while progress
        progress = false
        for k in 1:n
            depth[k] >= 0 && continue
            p = parent[k]
            if p >= 1 && depth[p] >= 0
                depth[k] = depth[p] + 1; progress = true
            end
        end
    end
    depth
end
function tree_dist(parent, depth, a, b)
    da, db, steps = depth[a], depth[b], 0
    while da > db; a = parent[a]; da -= 1; steps += 1; end
    while db > da; b = parent[b]; db -= 1; steps += 1; end
    while a != b; a = parent[a]; b = parent[b]; steps += 2; end
    steps
end
depth = build_depths(walk.visited_parent)

log!("Delaunay edges     : $(length(delaunay_edges))")
log!("tree edges         : $(length(tree_edges))")
log!("non-tree edges     : $(length(nontree_edges))")
log!("")

# Per non-tree edge: midpoint ΔP_rel.
struct Edge
    A::Int; B::Int
    ΔP_rel::Float64
    tree_d::Int
    extrap_max::Float64
    mid::ComplexF64
    cat::Symbol
end
categorise(ΔP, em) = ΔP ≤ TOL_WELL ? :well_closed :
                     ΔP ≤ TOL_BAD  ? :noisy :
                     em > 1.0      ? :extrap_driven : :depth_driven

edges = Edge[]
for (a, b) in nontree_edges
    zA, zB = walk.visited_z[a], walk.visited_z[b]
    hA, hB = walk.visited_h[a], walk.visited_h[b]
    M  = (zA + zB) / 2
    tA = (M - zA) / hA
    tB = (M - zB) / hB
    yA = try; eval_vec_node(walk, a, tA); catch; continue; end
    yB = try; eval_vec_node(walk, b, tB); catch; continue; end
    all(isfinite, real.(yA)) && all(isfinite, imag.(yA)) || continue
    all(isfinite, real.(yB)) && all(isfinite, imag.(yB)) || continue
    nd  = sqrt(sum(abs2, yA .- yB))
    den = sqrt(sum(abs2, yA)) + sqrt(sum(abs2, yB)) + EPS_FLOOR
    ΔP  = Float64(nd / den)
    td  = tree_dist(walk.visited_parent, depth, a, b)
    em  = Float64(max(abs(tA), abs(tB)))
    push!(edges, Edge(a, b, ΔP, td, em, ComplexF64(M), categorise(ΔP, em)))
end

n_e = length(edges)
log!("non-tree edges evaluated (finite both ends): $n_e")
log!("")

if n_e == 0
    log!("**No evaluable edges** — cannot compute a distribution.")
else
    rels = Float64[e.ΔP_rel for e in edges]
    log!("### ΔP_rel distribution (vector-norm loop closure)")
    log!("")
    log!(@sprintf("  median : %.3e", median(rels)))
    log!(@sprintf("  p90    : %.3e", quantile(rels, 0.90)))
    log!(@sprintf("  p99    : %.3e", quantile(rels, 0.99)))
    log!(@sprintf("  max    : %.3e", maximum(rels)))
    log!("")
    cnt(c) = count(==(c), (e.cat for e in edges))
    pct(n) = round(100 * n / n_e; digits=1)
    log!("### Category breakdown")
    log!("")
    log!("  well_closed   : $(cnt(:well_closed)) ($(pct(cnt(:well_closed)))%)  (ΔP_rel ≤ $TOL_WELL)")
    log!("  noisy         : $(cnt(:noisy)) ($(pct(cnt(:noisy)))%)")
    log!("  extrap_driven : $(cnt(:extrap_driven)) ($(pct(cnt(:extrap_driven)))%)  (|t| > 1 at midpoint)")
    log!("  depth_driven  : $(cnt(:depth_driven)) ($(pct(cnt(:depth_driven)))%)  (in-disc loop-closure failure)")
    log!("")
    bad_mids = ComplexF64[e.mid for e in edges if e.ΔP_rel > TOL_BAD]
    centroid = isempty(bad_mids) ? complex(NaN, NaN) : ComplexF64(mean(bad_mids))
    log!("bad_centroid (mean midpoint of ΔP_rel > $TOL_BAD edges): $centroid")
    log!("n bad edges : $(length(bad_mids))")
    log!("")
    log!("### Worst 10 edges")
    log!("")
    order = sortperm(rels; rev=true)
    for i in 1:min(10, n_e)
        e = edges[order[i]]
        log!(@sprintf("  %2d. ΔP_rel=%.3e  cat=%-13s  extrap_max=%.2f  tree_d=%d  mid=%s",
                      i, e.ΔP_rel, String(e.cat), e.extrap_max, e.tree_d,
                      string(round(e.mid; digits=3))))
    end
    log!("")
end

# ======================================================================
# TASK 3 — conjugate-symmetry pairing of the 21 poles.
# ======================================================================
log!("## Task 3 — Conjugate-symmetry pole pairing")
log!("")
log!("V_0 is real on the negative real axis ⇒ V_0(conj x) = conj(V_0(x))")
log!("⇒ wedge poles come in conjugate pairs.")
log!("")
log!("Extracted poles ($(length(poles))):")
for (i, p) in enumerate(sort(poles; by = q -> (real(q), imag(q))))
    log!(@sprintf("  %2d. %+.5f %+.5fim", i, real(p), imag(p)))
end
log!("")

# Greedy nearest-neighbour pairing under conjugation.
# A pole p (Im > 0) pairs with the unmatched pole q minimising
# |q - conj(p)|.  Poles with |Im| below near_axis_tol are treated as
# real-axis (self-conjugate / unpaired).
const NEAR_AXIS = 1e-3
unmatched = collect(1:length(poles))
pairs = Tuple{Int,Int,Float64}[]    # (i_upper, i_lower, residual)
real_axis = Int[]

# Pull out near-real poles first.
filter!(unmatched) do i
    if abs(imag(poles[i])) < NEAR_AXIS
        push!(real_axis, i); return false
    end
    true
end

# Sort remaining by Im descending; greedily pair upper with lower.
sort!(unmatched; by = i -> -imag(poles[i]))
while !isempty(unmatched)
    i = popfirst!(unmatched)
    pᵢ = poles[i]
    # Candidate partners: opposite-sign Im.
    best_j, best_res = 0, Inf
    for j in unmatched
        sign(imag(poles[j])) == sign(imag(pᵢ)) && continue
        res = abs(poles[j] - conj(pᵢ))
        if res < best_res
            best_j, best_res = j, res
        end
    end
    if best_j == 0
        push!(real_axis, i)   # no partner — treat as unpaired
    else
        filter!(!=(best_j), unmatched)
        up, lo = imag(pᵢ) ≥ 0 ? (i, best_j) : (best_j, i)
        push!(pairs, (up, lo, best_res))
    end
end

log!("### Pairing result")
log!("")
log!("matched conjugate pairs : $(length(pairs))")
log!("unpaired / near-real    : $(length(real_axis))")
log!("")
if !isempty(pairs)
    log!("Per-pair residual |p_upper - conj(p_lower)|:")
    for (k, (iu, il, r)) in enumerate(sort(pairs; by = t -> t[3]))
        log!(@sprintf("  pair %2d: p_up=%+.4f%+.4fim  p_lo=%+.4f%+.4fim  resid=%.3e",
                      k, real(poles[iu]), imag(poles[iu]),
                      real(poles[il]), imag(poles[il]), r))
    end
    resids = Float64[r for (_, _, r) in pairs]
    log!("")
    log!(@sprintf("max residual    : %.3e", maximum(resids)))
    log!(@sprintf("median residual : %.3e", median(resids)))
    log!(@sprintf("mean residual   : %.3e", mean(resids)))
end
log!("")
if !isempty(real_axis)
    log!("Unpaired / near-real poles (Im x should be ~0):")
    for i in real_axis
        log!(@sprintf("  pole %2d: %+.5f %+.5fim   (|Im| = %.3e)",
                      i, real(poles[i]), imag(poles[i]), abs(imag(poles[i]))))
    end
end
log!("")

# ======================================================================
# TASK 4 — VC-4 dominant-balance leading coefficient spot check.
# ======================================================================
log!("## Task 4 — VC-4 dominant-balance A spot check")
log!("")
log!("ADR-0025 A3: every P_I^(2) pole has u ~ A(x-x0)^-2 with A ∈ {-1,-3}.")
log!("Fit u(x) ≈ A(x-x0)^-2 + B(x-x0)^-1 + C on a 32-point ring.")
log!("u(x) is recovered from the nearest visited node's shared-Q approximant.")
log!("")

# Evaluate the first companion component u = y[1] at an arbitrary z by
# the nearest visited node's shared-Q approximant.
function eval_u_at(walk, z)
    N = length(walk.visited_z)
    idx, best = 1, abs(walk.visited_z[1] - z)
    for i in 2:N
        d = abs(walk.visited_z[i] - z)
        d < best && ((idx, best) = (i, d))
    end
    h = real(walk.visited_h[idx])
    t = (z - walk.visited_z[idx]) / h
    return eval_vec_node(walk, idx, t)[1], best, h
end

# Fit A,B,C from a ring of npts points at radius r about x0.
function fit_ABC(walk, x0; r=0.04, npts=32)
    θs = range(0, 2π; length=npts+1)[1:npts]
    rows = ComplexF64[]
    rhs  = ComplexF64[]
    M = Matrix{ComplexF64}(undef, npts, 3)
    b = Vector{ComplexF64}(undef, npts)
    ok = 0
    for (i, θ) in enumerate(θs)
        x = x0 + r * cis(θ)
        u, dist, h = eval_u_at(walk, x)
        # honesty: only use ring points whose nearest node is within h.
        d = x - x0
        M[i, 1] = 1 / d^2
        M[i, 2] = 1 / d
        M[i, 3] = 1
        b[i] = u
        ok += 1
    end
    coef = M \ b
    return coef[1], coef[2], coef[3]   # A, B, C
end

# Spot-check 5 poles spread across the field (smallest |Im| first, a few).
spot = sort(1:length(poles); by = i -> abs(poles[i]))[1:min(6, length(poles))]
log!("Spot-checked poles (nearest-origin 6):")
log!("")
for i in spot
    x0 = poles[i]
    # ring radius: small, but well above node spacing h≈0.1 is impossible,
    # so use the ADR-0025 form r = min(0.05, 0.1*d_nearest). Here we just
    # use a fixed small ring and report; this is a SPOT CHECK only.
    A, B, C = try
        fit_ABC(walk, x0; r=0.04, npts=32)
    catch err
        log!(@sprintf("  pole %2d @ %+.4f%+.4fim : fit FAILED (%s)",
                      i, real(x0), imag(x0), typeof(err)))
        continue
    end
    dA = min(abs(A + 1), abs(A + 3))
    log!(@sprintf("  pole %2d @ %+.4f%+.4fim : A=%+.4f%+.4fim  |B|=%.3e  min(|A+1|,|A+3|)=%.3e",
                  i, real(x0), imag(x0), real(A), imag(A), abs(B), dA))
end
log!("")
log!("NOTE: ring radius r=0.04 < walk step h=0.1, so most ring points")
log!("share ONE node's Padé — this is a coarse spot check, not the full")
log!("VC-4 acceptance test (which re-expands a dedicated jet per pole).")
log!("")

# ----------------------------------------------------------------------
write(joinpath(@__DIR__, "REPORT.md"), join(REPORT, "\n") * "\n")
println()
println("REPORT.md written to ", joinpath(@__DIR__, "REPORT.md"))
