# external/probes/fzse-calibration/probe.jl
#
# BLOCKING calibration probe for bug padetaylor-fzse (ADR / worklog 075).
# Measures, on the ACTUAL 121x121 equianharmonic path-network field, the two
# thresholds the read-only audition could only INFER (Law 1):
#   (i)  each cluster's within-cluster member SPREAD (the agreement-spread gate);
#   (ii) each candidate root's Pade RESIDUE |N(t*)/D'(t*)| (the Froissart gate).
# Goal: confirm the 4 FALSE POSITIVES separate from the 17 genuine lattice poles
# by some (spread, residue) so the composed fix (self-merge + agreement gate) is
# calibrated, not guessed.  Independent oracle = closed-form ℘ lattice.
#
# Run: julia --project=/tmp/ptcampaign external/probes/fzse-calibration/probe.jl

using PadeTaylor
using Polynomials: Polynomial, roots, derivative
using Printf

# --- equianharmonic field (verbatim from figures/demo_lattice_singularities.jl) -
fwp(z, u, up) = 6 * u^2
prob = PadeTaylorProblem(fwp, (1.071822516416917, 1.710337353176786),
                         (0.0, 1.5); order = 30)
const XLO, XHI, YLO, YHI, NX, NY = -4.0, 6.0, -5.0, 5.0, 121, 121
xs = range(XLO, XHI; length = NX); ys = range(YLO, YHI; length = NY)
grid = ComplexF64[complex(x, y) for y in ys for x in xs]
@printf("solving %d×%d path-network (h=0.5, order=30) ...\n", NX, NY); flush(stdout)
t0 = time()
sol = path_network_solve(prob, grid; h = 0.5, max_steps_per_target = 4000,
                         enforce_real_axis_symmetry = true)
@printf("  solved in %.1f s; %d tree nodes\n", time() - t0, length(sol.visited_z))
flush(stdout)

# --- closed-form ℘(z-1;0,2) lattice (INDEPENDENT oracle) ------------------------
const Ω = 1.3630340904278908                  # equianharmonic half-period (FW md:297)
const SPACING = 2Ω
true_all = ComplexF64[1 + 2Ω*(m + n/2) + im*(Ω*sqrt(3)*n) for m in -6:6 for n in -6:6]
inwin(p) = XLO ≤ real(p) ≤ XHI && YLO ≤ imag(p) ≤ YHI
true_in = filter(inwin, true_all)
@printf("true in-window lattice poles: %d ; spacing 2Ω = %.5f ; 0.5·spacing = %.4f\n",
        length(true_in), SPACING, 0.5SPACING); flush(stdout)

# --- replicate _extract_poles_core candidate gathering, KEEPING the residue -----
function gather(sol; radius_t = 5.0, min_residue = 1e-8)
    centers, pades, hs = sol.visited_z, sol.visited_pade, sol.visited_h
    h_max = maximum(abs, hs); radius_z = radius_t * h_max
    cands = NamedTuple[]
    for k in eachindex(hs)
        P = pades[k]; z_ctr = centers[k]; h = hs[k]
        length(P.b) ≥ 2 || continue
        D = Polynomial(P.b); N = Polynomial(P.a); Dp = derivative(D)
        for t in roots(D)
            z_dist = abs(h * t)
            z_dist ≤ radius_z || continue
            res = N(t) / Dp(t)
            (isfinite(res) && abs(res) ≥ min_residue) || continue
            push!(cands, (z = ComplexF64(z_ctr + h*t), zd = Float64(z_dist),
                          k = k, res = ComplexF64(res)))
        end
    end
    return cands
end

# greedy first-fit clustering (the real code), tracking member indices ----------
function cluster(cands; cluster_atol = 0.1, min_support = 3)
    order = sortperm(cands; by = c -> c.zd)
    reps = ComplexF64[]; members = Vector{Int}[]; support = Vector{Set{Int}}()
    for i in order
        c = cands[i]
        j = findfirst(r -> abs(c.z - r) ≤ cluster_atol, reps)
        if j === nothing
            push!(reps, c.z); push!(members, [i]); push!(support, Set{Int}((c.k,)))
        else
            push!(members[j], i); push!(support[j], c.k)
        end
    end
    keep = [j for j in eachindex(reps) if length(support[j]) ≥ min_support]
    return reps, members, support, keep
end

mindist(p, set) = isempty(set) ? Inf : minimum(q -> abs(p - q), set)

# FULL lattice (in + out of window) — boundary nodes legitimately see poles
# just OUTSIDE the window, so a rep near an out-of-window lattice point is a
# GENUINE pole, NOT a false positive.  Only a rep far from EVERY lattice point
# is a true spurious FP.
const FULL_LATTICE = ComplexF64[1 + 2Ω*(m + n/2) + im*(Ω*sqrt(3)*n)
                                for m in -10:10 for n in -10:10]

# --- AUTHORITATIVE cross-check: the REAL extract_poles vs the full lattice ----
function audit_real(cluster_atol, min_support, radius_t)
    reps = extract_poles(sol; cluster_atol, min_support, radius_t)
    tp   = filter(p -> mindist(p, FULL_LATTICE) ≤ 0.5SPACING, reps)
    fp   = filter(p -> mindist(p, FULL_LATTICE) > 0.5SPACING, reps)
    covered_in = filter(p -> mindist(p, reps) ≤ 0.5SPACING, true_in)
    missed_in  = filter(p -> mindist(p, reps) > 0.5SPACING, true_in)
    @printf("\n### REAL extract_poles  atol=%.2f min_support=%d radius_t=%.1f ###\n",
            cluster_atol, min_support, radius_t)
    @printf("  reps=%d | true-positive(near full lattice)=%d | GENUINE-FP(far from ALL)=%d\n",
            length(reps), length(tp), length(fp))
    @printf("  in-window recall: %d/%d covered | %d missed | precision(vs full lattice)=%.2f\n",
            length(covered_in), length(true_in), length(missed_in), length(tp)/max(length(reps),1))
    # duplicates: in-window poles with >1 rep within 0.5·spacing
    dups = 0
    for p in true_in
        n = count(r -> abs(r - p) ≤ 0.5SPACING, reps)
        n > 1 && (dups += 1)
    end
    @printf("  in-window poles with DUPLICATE reps (>1 within 0.5·spacing): %d\n", dups)
    # Decisive: distribution of rep-to-nearest-lattice distance.  A GOOD estimate
    # is < ~0.3; the bead claimed aliases at 0.93-1.24 (hidden inside the 0.5·
    # spacing=1.363 TP radius).  If NONE land in [0.4, 1.3], there are no aliases.
    dists = sort([mindist(p, FULL_LATTICE) for p in reps])
    nclose = count(<(0.4), dists); nalias = count(d -> 0.4 ≤ d ≤ 1.3, dists)
    @printf("  rep→nearest-lattice dist: <0.4 (good)=%d | [0.4,1.3] (ALIAS zone)=%d | max=%.3f\n",
            nclose, nalias, maximum(dists))
    if !isempty(fp)
        println("  -- GENUINE false positives (far from EVERY lattice point) --")
        for p in fp
            @printf("     %+.3f%+.3fim  dist-to-nearest-lattice=%.3f\n",
                    real(p), imag(p), mindist(p, FULL_LATTICE))
        end
    end
    if !isempty(missed_in)
        println("  -- MISSED in-window poles --")
        for p in missed_in
            @printf("     %+.3f%+.3fim  (Im=%.3f)\n", real(p), imag(p), imag(p))
        end
    end
    return reps, fp, missed_in
end

function analyse(cands, cluster_atol, min_support)
    reps, members, support, keep = cluster(cands; cluster_atol, min_support)
    kept_reps = [reps[j] for j in keep]
    fp = Int[]; genuine = Int[]
    for j in keep
        (mindist(reps[j], true_in) > 0.5SPACING ? push!(fp, j) : push!(genuine, j))
    end
    missed = [p for p in true_in if mindist(p, kept_reps) > 0.5SPACING]
    @printf("\n=== cluster_atol=%.2f, min_support=%d ===\n", cluster_atol, min_support)
    @printf("  kept reps=%d  genuine=%d  FALSE-POS=%d  MISSED=%d  (precision=%.2f recall=%.2f)\n",
            length(keep), length(genuine), length(fp), length(missed),
            length(genuine)/max(length(keep),1), length(genuine)/length(true_in))
    # per-cluster spread + residue, FP vs genuine
    function clinfo(j)
        ms = members[j]; zs = [cands[i].z for i in ms]; rep = reps[j]
        spread = isempty(zs) ? 0.0 : maximum(z -> abs(z - rep), zs)
        resids = [abs(cands[i].res) for i in ms]
        return (spread = spread, nmem = length(ms), nsupp = length(support[j]),
                resmin = minimum(resids), resmax = maximum(resids),
                resmed = sort(resids)[cld(length(resids),2)])
    end
    println("  -- FALSE POSITIVES --")
    for j in fp
        ci = clinfo(j)
        @printf("     rep=%+.3f%+.3fim  d2true=%.3f  spread=%.4f  supp=%d  res[min/med/max]=%.2e/%.2e/%.2e\n",
                real(reps[j]), imag(reps[j]), mindist(reps[j], true_in),
                ci.spread, ci.nsupp, ci.resmin, ci.resmed, ci.resmax)
    end
    println("  -- GENUINE (spread / residue distribution) --")
    gspreads = Float64[]; gres = Float64[]
    for j in genuine
        ci = clinfo(j); push!(gspreads, ci.spread); push!(gres, ci.resmed)
    end
    if !isempty(gspreads)
        @printf("     spread:  min=%.4f  med=%.4f  max=%.4f\n",
                minimum(gspreads), sort(gspreads)[cld(length(gspreads),2)], maximum(gspreads))
        @printf("     res_med: min=%.2e med=%.2e max=%.2e\n",
                minimum(gres), sort(gres)[cld(length(gres),2)], maximum(gres))
    end
    if !isempty(fp)
        fspreads = [clinfo(j).spread for j in fp]
        @printf("  >>> SEPARATION: FP spread min=%.4f vs genuine spread max=%.4f  (gate viable iff FP_min > genuine_max)\n",
                minimum(fspreads), isempty(gspreads) ? NaN : maximum(gspreads))
    end
    return missed
end

# --- THE FIX: post-pass single-linkage self-merge of the survivors -------------
# Merge reps within merge_atol (keeping the best-placed = smallest-index, since
# extract_poles returns reps in z_dist discovery order).  At the TIGHT default
# atol=0.1 the survivors are all genuine (aliases dropped by min_support), so
# merging collapses the over-split duplicates with zero alias risk.
function selfmerge(reps, merge_atol)
    n = length(reps); parent = collect(1:n)
    find(x) = parent[x] == x ? x : (parent[x] = find(parent[x]))
    for a in 1:n, b in (a+1):n
        if abs(reps[a] - reps[b]) ≤ merge_atol
            ra, rb = minmax(find(a), find(b)); parent[rb] = ra
        end
    end
    return [reps[a] for a in 1:n if find(a) == a]
end

function merged_audit(cluster_atol, min_support, merge_atol)
    reps0 = extract_poles(sol; cluster_atol, min_support, radius_t = 5.0)
    reps  = selfmerge(reps0, merge_atol)
    fp = count(p -> mindist(p, FULL_LATTICE) > 0.5SPACING, reps)
    cov = count(p -> mindist(p, reps) ≤ 0.5SPACING, true_in)
    dups = count(p -> count(r -> abs(r - p) ≤ 0.5SPACING, reps) > 1, true_in)
    alias = count(p -> 0.4 ≤ mindist(p, FULL_LATTICE) ≤ 1.3, reps)
    @printf("  SELF-MERGE atol=%.2f msupp=%d merge=%.2f: %d→%d reps | dups=%d | FP=%d | alias=%d | recall=%d/%d\n",
            cluster_atol, min_support, merge_atol, length(reps0), length(reps),
            dups, fp, alias, cov, length(true_in))
end

println("\n=== SELF-MERGE calibration (at the TIGHT default atol; sweep merge_atol) ===")
for m in (0.4, 0.6, 0.8, 1.0)
    merged_audit(0.1, 3, m)
end
println("  (recall fix: min_support 3→2 recovers the seen-by-2 bottom-edge poles)")
for m in (0.4, 0.6)
    merged_audit(0.1, 2, m)
end

# AUTHORITATIVE: the real extract_poles vs the FULL lattice (in + out of window).
audit_real(0.1, 3, 5.0)      # stock default
audit_real(0.4, 3, 5.0)      # bead's "tuned" cluster_atol
audit_real(0.4, 2, 5.0)      # + the proposed recall fix (min_support 3->2)
audit_real(0.6, 2, 5.0)      # wider merge + recall fix

# Member-spread / residue diagnostics (only meaningful if GENUINE FPs exist).
cands = gather(sol; radius_t = 5.0, min_residue = 1e-8)
@printf("\ntotal candidates after far-root + Froissart(1e-8) filter: %d\n", length(cands))
println("\n(member-spread classification below uses IN-WINDOW poles only — note the")
println(" REAL-extract_poles audit above is the authoritative precision/recall.)")
analyse(cands, 0.4, 3)
