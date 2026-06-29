# P4b — RECONCILIATION viability (part (b) of mechanism selection; p4 part (a)
# already showed windowing shrinks pole-count seed-variance ~5-13x, far incl.).
#
# Is the disc-overlap graph rich in CROSS-BRANCH edges — two spatially-overlapping
# nodes whose nearest common tree-ancestor is far back? Those are the ONLY
# constraints that can see the seam. If overlaps are almost all parent-child
# (same branch), reconciliation is blind to it (= naive K-NN, already refuted).
#
# (Analysis wrapped in a function to avoid Julia top-level soft-scope on the
# accumulators — the bug that aborted p4 part (b).)

using PadeTaylor, Printf
pI(z, u, up) = 6u^2 + z
const HALF = 50.0; const N = 101
xs = range(-HALF, HALF; length = N); ys = range(-HALF, HALF; length = N)
prob = PadeTaylorProblem(pI, (-0.1875, 0.3049), (0.0, HALF); order = 30)

t0 = time()
pn = path_network_solve(prob, ComplexF64[complex(x,y) for y in ys for x in xs];
                        h = 0.5, order = 30, rng_seed = 0, extrapolate = true)
@printf("monolithic seed-0: nodes=%d  %.1fs\n", length(pn.visited_z), time()-t0); flush(stdout)

function analyze_overlap(VZ, VH, PAR; D = 16, BIN = 1.0)
    Nn = length(VZ)
    ancestors(k) = (s = Set{Int}(); c = k; for _ in 1:D; c == 0 && break; push!(s, c); c = PAR[c]; end; s)
    ANC = [ancestors(k) for k in 1:Nn]
    binkey(z) = (floor(Int, real(z)/BIN), floor(Int, imag(z)/BIN))
    bins = Dict{Tuple{Int,Int},Vector{Int}}()
    for k in 1:Nn; push!(get!(bins, binkey(VZ[k]), Int[]), k); end
    inann(z) = (r = abs(z); 10.0 <= r <= 40.0)
    nedge = 0; ncross = 0; nedge_a = 0; ncross_a = 0; mult = zeros(Int, Nn)
    for k in 1:Nn
        bk = binkey(VZ[k])
        for dbi in -1:1, dbj in -1:1
            for j in get(bins, (bk[1]+dbi, bk[2]+dbj), Int[])
                j <= k && continue
                abs(VZ[k]-VZ[j]) <= (VH[k]+VH[j]) || continue   # discs intersect
                nedge += 1; mult[k] += 1; mult[j] += 1
                cross = !(j in ANC[k]) && !(k in ANC[j]) && isempty(intersect(ANC[k], ANC[j]))
                cross && (ncross += 1)
                if inann(VZ[k]) && inann(VZ[j]); nedge_a += 1; cross && (ncross_a += 1); end
            end
        end
    end
    sort!(mult)
    (Nn, nedge, ncross, nedge_a, ncross_a, mult, D)
end

Nn, nedge, ncross, nedge_a, ncross_a, mult, D =
    analyze_overlap(pn.visited_z, pn.visited_h, pn.visited_parent)

println("="^72); println("(b) RECONCILIATION — disc-overlap graph cross-branch richness"); println("="^72)
@printf("nodes=%d  overlap edges=%d  (avg degree=%.1f)\n", Nn, nedge, 2nedge/Nn)
@printf("disc-coverage multiplicity/node: median=%d  p90=%d  max=%d  (>=2 overlaps: %.1f%%)\n",
        mult[(Nn+1)÷2], mult[clamp(round(Int,0.9Nn),1,Nn)], mult[end], 100count(>=(2),mult)/Nn)
@printf("CROSS-BRANCH overlaps (nearest common ancestor > %d hops):\n", D)
@printf("   global:                 %d / %d = %.1f%%\n", ncross, nedge, 100ncross/max(nedge,1))
@printf("   seam annulus 10<=|z|<=40: %d / %d = %.1f%%\n", ncross_a, nedge_a, 100ncross_a/max(nedge_a,1))
println("""
READOUT: reconciliation is viable iff cross-branch overlap % is substantial in the
seam annulus. ~0% => overlaps are same-branch, reconciliation blind to the seam.
High => rich cross-branch constraints exist to re-derive node states against.
""")
