# P4 — MECHANISM SELECTION (bead vwgl): pick the in-pole-field cure on evidence.
#
# P3 proved the seam is a pole-field-INTERIOR grain boundary (regime B): edge-
# gating leaves the pole-count seed-variance unchanged (|Δ| 172 -> 171). The cure
# must act INSIDE the pole field. Two live mechanisms remain:
#   #1 bounded-window composite (containment via short shared trunks, FW md:147)
#   #4/#5 overlap reconciliation (re-derive node states from cross-branch constraints)
# This probe measures the load-bearing premise of each.
#
# (a) WINDOWING: does an independent <=20x20 window solve shrink the two-seed
#     pole-count/position variance vs the monolithic solve over the SAME region?
#     Tested for a NEAR window (short trunk) AND a FAR window (long trunk — the
#     verifier's objection: a far tile's trunk may be as long as the monolithic's).
# (b) RECONCILIATION: is the disc-overlap graph rich in CROSS-BRANCH edges (two
#     spatially-overlapping nodes whose nearest common tree-ancestor is far back)?
#     Those are the only constraints that can see the seam. If overlaps are almost
#     all parent-child (same branch), reconciliation is blind to it (= naive K-NN).

using PadeTaylor, Printf
pI(z, u, up) = 6u^2 + z
const HALF = 50.0; const N = 101
xs = range(-HALF, HALF; length = N); ys = range(-HALF, HALF; length = N)
prob = PadeTaylorProblem(pI, (-0.1875, 0.3049), (0.0, HALF); order = 30)

function poleset_agreement(A, B; tol = 0.5)
    isempty(A) && return (0.0, NaN, NaN)
    matched = 0; nn = Float64[]
    for p in A
        d = Inf
        for q in B; dd = abs(p - q); dd < d && (d = dd); end
        push!(nn, d); d <= tol && (matched += 1)
    end
    sort!(nn)
    (matched / length(A), nn[(length(nn)+1)÷2], maximum(nn))
end
restrict(ps, xlo, xhi, ylo, yhi) = filter(p -> xlo <= real(p) <= xhi && ylo <= imag(p) <= yhi, ps)

println("="^72); println("MONOLITHIC full-grid solve, seeds 0 & 42 (shared baseline)"); println("="^72); flush(stdout)
mono = Dict{Int,Any}()
for seed in (0, 42)
    t0 = time()
    pn = path_network_solve(prob, ComplexF64[complex(x,y) for y in ys for x in xs];
                            h = 0.5, order = 30, rng_seed = seed, extrapolate = true)
    mono[seed] = (pn, extract_poles(pn))
    @printf("  mono seed=%d nodes=%d poles=%d  %.1fs\n", seed, length(pn.visited_z), length(mono[seed][2]), time()-t0); flush(stdout)
end

# --- (a) WINDOWING: near (contains IC) vs far (long trunk) ------------------
println("\n" * "="^72); println("(a) WINDOWING — windowed vs monolithic two-seed pole variance, same core"); println("="^72); flush(stdout)
# (label, window bounds, core inset)
WINDOWS = [("NEAR [-10,10]^2 (short trunk, has IC)", -10.0, 10.0, -10.0, 10.0),
           ("FAR  [18,38]x[-10,10] (long trunk)",     18.0, 38.0, -10.0, 10.0)]
for (lab, xlo, xhi, ylo, yhi) in WINDOWS
    ins = 2.0; cx0,cx1,cy0,cy1 = xlo+ins, xhi-ins, ylo+ins, yhi-ins
    # monolithic restricted to this core
    m0 = restrict(mono[0][2], cx0,cx1,cy0,cy1); m42 = restrict(mono[42][2], cx0,cx1,cy0,cy1)
    mmatch,_,mmax = poleset_agreement(m0, m42)
    # independent windowed solve over just this window's cells, both seeds
    zc = ComplexF64[complex(x,y) for y in ys for x in xs if xlo<=x<=xhi && ylo<=y<=yhi]
    wp = Dict{Int,Any}()
    for seed in (0,42)
        t0=time(); pn = path_network_solve(prob, zc; h=0.5, order=30, rng_seed=seed, extrapolate=true)
        wp[seed] = restrict(extract_poles(pn), cx0,cx1,cy0,cy1)
        @printf("    [win %s] seed=%d nodes=%d core-poles=%d  %.1fs\n", lab[1:4], seed, length(pn.visited_z), length(wp[seed]), time()-t0); flush(stdout)
    end
    wmatch,_,wmax = poleset_agreement(wp[0], wp[42])
    println("  -- $lab --")
    @printf("    MONO   core poles s0=%d s42=%d |Δ|=%d  match@0.5=%.1f%%  maxNN=%.2f\n",
            length(m0), length(m42), abs(length(m0)-length(m42)), 100mmatch, mmax)
    @printf("    WINDOW core poles s0=%d s42=%d |Δ|=%d  match@0.5=%.1f%%  maxNN=%.2f\n",
            length(wp[0]), length(wp[42]), abs(length(wp[0])-length(wp[42])), 100wmatch, wmax); flush(stdout)
end

# --- (b) RECONCILIATION: cross-branch disc-overlap richness -----------------
println("\n" * "="^72); println("(b) RECONCILIATION — disc-overlap graph: cross-branch edge richness"); println("="^72); flush(stdout)
pn = mono[0][1]
VZ = pn.visited_z; VH = pn.visited_h; PAR = pn.visited_parent; Nn = length(VZ)
# ancestor set up to depth D (1-indexed; 0 = root sentinel)
function ancestors(k, D)
    s = Set{Int}(); c = k
    for _ in 1:D
        c == 0 && break; push!(s, c); c = PAR[c]
    end
    s
end
const D = 16
ANC = [ancestors(k, D) for k in 1:Nn]
# spatial bins of size = max overlap radius (h_j + h_k ~ 1.0 for fixed h=0.5)
const BIN = 1.0
binkey(z) = (floor(Int, real(z)/BIN), floor(Int, imag(z)/BIN))
bins = Dict{Tuple{Int,Int},Vector{Int}}()
for k in 1:Nn; push!(get!(bins, binkey(VZ[k]), Int[]), k); end
# annulus where the seam lives (mid-field), to report a localized number too
inann(z) = (r = abs(z); 10.0 <= r <= 40.0)
nedge=0; ncross=0; nedge_a=0; ncross_a=0; mult=zeros(Int, Nn)
for k in 1:Nn
    bk = binkey(VZ[k])
    for dbi in -1:1, dbj in -1:1
        nb = get(bins, (bk[1]+dbi, bk[2]+dbj), Int[])
        for j in nb
            j <= k && continue
            d = abs(VZ[k]-VZ[j])
            d <= (VH[k]+VH[j]) || continue          # discs intersect
            nedge += 1; mult[k]+=1; mult[j]+=1
            cross = !(j in ANC[k]) && !(k in ANC[j]) && isempty(intersect(ANC[k], ANC[j]))
            cross && (ncross += 1)
            if inann(VZ[k]) && inann(VZ[j])
                nedge_a += 1; cross && (ncross_a += 1)
            end
        end
    end
end
sort!(mult)
@printf("nodes=%d  overlap edges=%d  (avg degree=%.1f)\n", Nn, nedge, 2nedge/Nn)
@printf("disc-coverage multiplicity per node: median=%d  p90=%d  max=%d  (frac with >=2 overlaps: %.1f%%)\n",
        mult[(Nn+1)÷2], mult[clamp(round(Int,0.9Nn),1,Nn)], mult[end], 100count(>=(2),mult)/Nn)
@printf("CROSS-BRANCH overlap edges (common ancestor > %d hops back):\n", D)
@printf("   global: %d / %d = %.1f%% of all overlaps\n", ncross, nedge, 100ncross/max(nedge,1))
@printf("   mid-field annulus 10<=|z|<=40 (seam zone): %d / %d = %.1f%%\n", ncross_a, nedge_a, 100ncross_a/max(nedge_a,1))

println("\n" * "="^72); println("MECHANISM-SELECTION READOUT"); println("="^72)
println("""
(a) WINDOWING works iff WINDOW |Δ|/maxNN << MONO for BOTH windows (esp. FAR).
    If FAR window is NOT better than mono -> long-trunk objection holds; windowing
    only helps the near field and needs a window-local IC to cure far tiles.
(b) RECONCILIATION viable iff cross-branch overlap % is substantial (esp. in the
    seam annulus). If ~0% -> overlaps are all same-branch, reconciliation is blind
    to the seam (= naive K-NN). If high -> rich constraints to re-derive states.
""")
