# froissart_residue.jl — calibrate the residue threshold that separates a GENUINE
# shared-Q pole from a cell-B Froissart doublet (the ADR-0028 Phase-2a finding,
# bead padetaylor-?).  The wiring is gated on a fix to the vector pole-root
# consumers; this probe measures, on the real ℘ FW Table 5.1 walk, whether
# genuine poles and cell-B spurious near-node poles separate cleanly by residue
# so a residue filter (matching scalar PoleField) is well-founded.
#
# Method: run the WORKING cell-A walk to z=31 (the VPO.2 setup), then at every
# visited node rebuild cell A (shared_denominator_pade) AND cell B
# (build_square_cell), find each one's nearest denominator root, its z-plane
# distance h·|t*|, and its genuineness residue max_c |P_c(t*)/Q'(t*)|; also note
# which cell the ODE defect would select.  We expect: cell A's nearest roots are
# genuine ℘ poles (substantial residue); cell B introduces near-node doublets on
# regular steps (residue ~ machine ε).  The gap is the threshold.
#
# Run from repo root:
#   julia --project=. external/probes/adr0028-froissart-consumer/froissart_residue.jl

using PadeTaylor
using PadeTaylor: shared_denominator_pade, vector_path_network_solve,
                  VectorPadeTaylorProblem
using PadeTaylor.SharedPadeCellB:  build_square_cell
using PadeTaylor.SharedPadeDefect: relative_defect, guard_root_estimates
using Printf

const u_0_FW  = 1.071822516416917
const up_0_FW = 1.710337353176786
fW_vec(z, y) = [y[2], 6 * y[1]^2]

horner(c, t) = (s = zero(t) * zero(eltype(c)); for k = length(c):-1:1; s = s * t + c[k]; end; s)
polyder(c) = length(c) < 2 ? eltype(c)[zero(eltype(c))] : [(k - 1) * c[k] for k = 2:length(c)]

# nearest root t* of Q, its z-distance h·|t*|, and the genuineness residue
# max_c |P_c(t*)/Q'(t*)| (a Froissart doublet ⇒ tiny; a genuine pole ⇒ O(1)+).
function probe_cell(nums, den, h)
    rts = guard_root_estimates(den)
    isempty(rts) && return (Inf, Inf, NaN)
    _, j = findmin(abs.(rts))
    tstar = rts[j]
    zdist = abs(h) * abs(tstar)
    qp = horner(ComplexF64.(polyder(den)), tstar)
    resid = abs(qp) < 1e-300 ? Inf :
            maximum(abs(horner(ComplexF64.(nums[c]), tstar) / qp) for c = 1:length(nums))
    return (abs(tstar), zdist, resid)
end

println("="^96)
println("FROISSART-CONSUMER CALIBRATION — ℘ FW Table 5.1 walk (z→31), cell A vs cell B per node")
println("="^96)

prob = VectorPadeTaylorProblem(fW_vec, ComplexF64[u_0_FW, up_0_FW],
                               (0.0 + 0im, 31.0 + 0im); order = 30)
walk = vector_path_network_solve(prob, ComplexF64[31.0 + 0im];
                                 h = 0.15, fine_grid = ComplexF64[30.0 + 0im],
                                 extrapolate = true, max_steps_per_target = 100_000,
                                 on_target_failure = :skip)
N = length(walk.visited_z)
@printf("visited nodes: %d   (max Re z = %.2f)\n\n", N, maximum(real, walk.visited_z))

m = 15
# accumulators
genA = Float64[]; genB = Float64[]            # residue at nearest root, A and B
closerB = NamedTuple[]                         # nodes where B's nearest pole is much closer than A's
nB_built = 0; nB_picked = 0; nB_picked_closer = 0
for i = 1:N
    global nB_built, nB_picked, nB_picked_closer
    jet = walk.visited_jets[i]; h = walk.visited_h[i]
    rescaled = [[h^(k - 1) * jet[c][k] for k = 1:length(jet[c])] for c = 1:length(jet)]
    numsA, denA = shared_denominator_pade(rescaled, m)
    cb = build_square_cell(rescaled, m)
    tA, zA, rA = probe_cell(numsA, denA, h)
    push!(genA, rA)
    cb === nothing && continue
    nB_built += 1
    numsB, denB = cb
    tB, zB, rB = probe_cell(numsB, denB, h)
    push!(genB, rB)
    # which cell would the defect pick?
    dA = relative_defect(numsA, denA, fW_vec, walk.visited_z[i], h)
    dB = relative_defect(numsB, denB, fW_vec, walk.visited_z[i], h)
    picksB = isfinite(dB) && dB < dA * (1 - 100 * eps())
    picksB && (nB_picked += 1)
    # cell B introduces a markedly closer pole (the collapse risk)?
    if zB < 0.3 * zA
        nB_picked_closer += picksB ? 1 : 0
        push!(closerB, (i = i, rez = real(walk.visited_z[i]), zA = zA, zB = zB,
                        rA = rA, rB = rB, picksB = picksB))
    end
end

quant(v, p) = isempty(v) ? NaN : sort(v)[clamp(ceil(Int, p * length(v)), 1, length(v))]
@printf("cell A nearest-root residue: min=%.1e  median=%.1e  max=%.1e  (n=%d)\n",
        minimum(genA), quant(genA, 0.5), maximum(genA), length(genA))
@printf("cell B nearest-root residue: min=%.1e  median=%.1e  max=%.1e  (n=%d built)\n",
        isempty(genB) ? NaN : minimum(genB), quant(genB, 0.5),
        isempty(genB) ? NaN : maximum(genB), length(genB))
@printf("\ncell B selected by defect at %d / %d built nodes\n", nB_picked, nB_built)
@printf("nodes where cell B's nearest pole is < 0.3× cell A's z-distance: %d (of which defect picks B: %d)\n",
        length(closerB), nB_picked_closer)

println("\n--- the 'cell B much closer pole' nodes (the collapse risk): residue of A (genuine) vs B (doublet?) ---")
@printf("%-5s %-8s %-11s %-11s %-11s %-11s %s\n", "node", "Re z", "zA(genuine)", "zB(cellB)", "resA", "resB", "defect→B?")
for c in first(closerB, 25)
    @printf("%-5d %-8.2f %-11.2e %-11.2e %-11.2e %-11.2e %s\n",
            c.i, c.rez, c.zA, c.zB, c.rA, c.rB, c.picksB ? "YES" : "no")
end

# separation summary: residues of B's *near-node* roots vs A's genuine poles
nearB_res = [c.rB for c in closerB]
@printf("\nSEPARATION: cell-B near-node-pole residues: min=%.1e median=%.1e max=%.1e\n",
        isempty(nearB_res) ? NaN : minimum(nearB_res), quant(nearB_res, 0.5),
        isempty(nearB_res) ? NaN : maximum(nearB_res))
@printf("            cell-A genuine-pole residues:    min=%.1e median=%.1e max=%.1e\n",
        minimum(genA), quant(genA, 0.5), maximum(genA))
println("="^96)
