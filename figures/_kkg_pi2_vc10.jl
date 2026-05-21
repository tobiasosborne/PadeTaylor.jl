# figures/_kkg_pi2_vc10.jl
#
# VC-10 — the FW/FFW two-run pole-field accuracy indicator of the KKG
# 2015 P_I^(2) tritronquée headline figure (ADR-0025 Phase D, bead
# `padetaylor-0ln.37.17`; Amendment 9).  This helper is `include`d by
# the figure kernel `_kkg_pi2_surface_helpers.jl` — the Rule-6 file-size
# split, mirroring `_kkg_pi2_vc45.jl` (VC-4/VC-5) and `_kkg_pi2_vc8.jl`
# (VC-8).
#
# ## Why a two-run indicator
#
# ADR-0025 Validation Criteria Menu VC-10.  FW 2011 §3.1 (line 156)
# shuffles the Stage-1 targets into a random order before the
# path-network walk; FFW 2017 (`references/markdown/FFW2017_painleve_
# riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.
# md:246-247`) turns that randomness into a *diagnostic*: "the PFS
# method selects the target nodes in Stage 1 in a random order ... if we
# compute the same solution twice, different paths will be run,
# resulting in solutions that should differ by approximately the
# numerical error."  The two-run pole-location disagreement is a
# practical, oracle-free accuracy estimate of the wedge pole field —
# exactly the manufactured FW-style criterion ADR-0025 Decision 1
# commits to (the tritronquée wedge field has no external oracle: KKG
# published no pole table).
#
# ## How VC-10 works
#
# `surf_vc10_two_run` runs the figure's wedge walk
# (`surf_wedge_fill`) **twice** with two different `MersenneTwister`
# target orderings (the `rng` kwarg of `vector_path_network_solve`,
# ADR-0025 Amendment 9 — a different ordering builds a different
# path-network tree, a different node filament, a different
# independently-extracted pole field).  Each run's candidate field is
# VC-4-validated (`vc4_validate`); the two validated fields are matched
# nearest-neighbour by `surf_vc10_match`, a maximum-cardinality
# bipartite matching that **reuses `_vc5_augment!`** — the Kuhn
# augmenting-path core of the VC-5 conjugate matcher (`_kkg_pi2_vc45.jl`).
# Per FFW the matched-pole disagreement IS the accuracy estimate, and
# the poles found by only one run are a second facet (each ordering's
# adaptive-`h` history resolves a slightly different subset of the dense
# far-wedge lattice — the A2 tractability ceiling).
#
# ## Determinism — VC-10 ≠ a broken reproducibility guarantee
#
# The walk for each *fixed* seed is fully deterministic (a fixed seed ⇒
# a fixed shuffle ⇒ a fixed tree, bit-identical) — the per-fixed-seed
# reproducibility guarantee `test/kkg_pi2_figure_test.jl` PI2F.1.6
# pins.  VC-10 is the *complementary* fact: *different* orderings
# converge to the *same* pole field within the reported ~0.35 accuracy.
# Both coexist; the figure itself passes no `rng` and so walks the fixed
# default order.
#
# References:
#   * `docs/adr/0025-headline-figure-re-resolution.md` — the
#     Validation Criteria Menu (VC-10), Amendment 9 (the measured
#     two-run disagreement, the `rng` kwarg).
#   * `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#     FFW2017_painleve_riemann_surfaces_preprint.md:246-247` — the
#     FW/FFW double-run accuracy criterion.
#   * `figures/_kkg_pi2_vc45.jl` — `_vc5_augment!`, the maximum-
#     cardinality augmenting-path core this helper reuses.
#   * `src/VectorPathNetwork.jl` — the `rng` kwarg (the controllable
#     Stage-1 target ordering VC-10 stands on).

using Random: MersenneTwister
using PadeTaylor.VectorBVP: VectorBVPSolution

# VC-10 two-run target-shuffle seeds.  Two fixed, distinct seeds — the
# runs stay reproducible (a fixed seed ⇒ a fixed shuffle ⇒ a fixed
# tree), but the two seeds give two different Stage-1 orderings hence
# two independently-threaded path-network trees.
const SURF_VC10_SEED_A = 20250521
const SURF_VC10_SEED_B = 19920811

# VC-10 pole-to-pole match tolerance.  Two poles from the two runs are
# "the same physical pole" when within this z-plane distance.  Set from
# the measured VC-5 field geometry (`_kkg_pi2_vc45.jl`): the pole
# nearest-neighbour spacing is ~0.69, so a match radius of 0.6 (87 % of
# the spacing) cannot mis-pair two distinct lattice poles, while being
# wide enough to admit a genuine same-pole pair separated only by the
# two-run numerical error.  Reuses the VC-5 `VC5_MATCH_TOL` rationale.
const SURF_VC10_MATCH_TOL = 0.6

"""
    surf_vc10_match(field_a, field_b; tol = SURF_VC10_MATCH_TOL)
        -> NamedTuple

Match two pole fields `field_a`, `field_b` (each a `Vector{Complex}`)
nearest-neighbour into "same physical pole" pairs.  An edge `(p_a, p_b)`
is admissible when `|p_a - p_b| ≤ tol`; among the admissible edges a
**maximum-cardinality** bipartite matching is computed by the Kuhn
augmenting-path core `_vc5_augment!` (shared with the VC-5 conjugate
matcher, `_kkg_pi2_vc45.jl`).  Each `field_a` pole's admissible list is
sorted by ascending distance so the augmenting search prefers the
tightest pairing — keeping the reported disagreement an honest
FW/FFW-style accuracy estimate, not an arbitrary one.

Returns `(pairs, disagreements, median_disagree, max_disagree,
only_a, only_b)`:
  - `pairs`           : `Vector{Tuple{ComplexF64,ComplexF64}}` of
                        matched `(p_a, p_b)`;
  - `disagreements`   : the per-pair `|p_a - p_b|`;
  - `median_disagree`, `max_disagree` : the FW/FFW accuracy estimate
                        (`NaN` when no pair matched);
  - `only_a`, `only_b`: the poles found by only one of the two runs.
"""
function surf_vc10_match(field_a::AbstractVector{<:Complex},
                         field_b::AbstractVector{<:Complex};
                         tol::Real = SURF_VC10_MATCH_TOL)
    A = ComplexF64.(field_a)
    B = ComplexF64.(field_b)
    na, nb = length(A), length(B)
    # The admissibility graph: an edge wherever two poles are within
    # `tol`.  Each A-pole's adjacency list is sorted by ascending
    # distance so the augmenting search prefers the closest partner.
    adj = [Int[] for _ in 1:na]
    for ia in 1:na
        for ib in 1:nb
            abs(A[ia] - B[ib]) ≤ tol && push!(adj[ia], ib)
        end
        sort!(adj[ia]; by = ib -> abs(A[ia] - B[ib]))
    end
    # Maximum-cardinality matching by Kuhn augmenting paths (the VC-5
    # `_vc5_augment!` core, reused — `_kkg_pi2_vc45.jl`).
    match_b = zeros(Int, nb)            # B index → matched A, 0 = free
    for ia in 1:na
        _vc5_augment!(ia, adj, match_b, fill(false, nb))
    end
    pairs         = Tuple{ComplexF64,ComplexF64}[]
    disagreements = Float64[]
    matched_a     = falses(na)
    for ib in 1:nb
        ia = match_b[ib]
        ia == 0 && continue
        matched_a[ia] = true
        push!(pairs, (A[ia], B[ib]))
        push!(disagreements, abs(A[ia] - B[ib]))
    end
    matched_b = falses(nb)
    for ib in 1:nb
        match_b[ib] != 0 && (matched_b[ib] = true)
    end
    only_a = A[.!matched_a]
    only_b = B[.!matched_b]
    med = isempty(disagreements) ? NaN :
          sort(disagreements)[cld(length(disagreements), 2)]
    mx  = isempty(disagreements) ? NaN : maximum(disagreements)
    return (pairs = pairs, disagreements = disagreements,
            median_disagree = med, max_disagree = mx,
            only_a = only_a, only_b = only_b)
end

"""
    surf_vc10_two_run(bvp_sol, grid_pts;
                      seed_a = SURF_VC10_SEED_A,
                      seed_b = SURF_VC10_SEED_B) -> NamedTuple

VC-10 — the FW/FFW double-run pole-field accuracy indicator (ADR-0025
Validation Criteria Menu VC-10, Amendment 9; FFW 2017 `:246-247`).

Run the figure's wedge walk (`surf_wedge_fill`) **twice** with two
different `MersenneTwister` target orderings (`seed_a`, `seed_b`).  Each
run is fully deterministic for its fixed seed — a fixed seed ⇒ a fixed
target shuffle ⇒ a fixed path-network tree — but the two seeds drive
two *different* Stage-1 orderings, hence two independently-threaded node
filaments and two independently-extracted, VC-4-validated pole fields.

The two fields are matched nearest-neighbour by `surf_vc10_match` (the
VC-5 maximum-cardinality augmenting-path matcher, reused).  Per FFW the
matched-pole disagreement IS a practical accuracy estimate of the pole
field — the walk's own numerical error surfaced without an external
oracle.

Returns `(walk_a, walk_b, poles_a, poles_b, match, message)`:
  - `walk_a`, `walk_b`   : the two `VectorPathNetworkSolution`s;
  - `poles_a`, `poles_b` : the two VC-4-validated pole fields;
  - `match`              : the `surf_vc10_match` result — `pairs`,
                           `disagreements`, `median_disagree`,
                           `max_disagree`, `only_a`, `only_b`;
  - `message`            : the FW/FFW accuracy-indicator status note.
"""
function surf_vc10_two_run(bvp_sol::VectorBVPSolution,
                           grid_pts::AbstractVector{ComplexF64};
                           seed_a::Integer = SURF_VC10_SEED_A,
                           seed_b::Integer = SURF_VC10_SEED_B)
    run_a = surf_wedge_fill(bvp_sol, grid_pts;
                            rng = MersenneTwister(seed_a))
    run_b = surf_wedge_fill(bvp_sol, grid_pts;
                            rng = MersenneTwister(seed_b))
    m = surf_vc10_match(run_a.poles, run_b.poles)
    msg = "VC-10 two-run: run A $(length(run_a.poles)) poles " *
          "($(length(run_a.walk.visited_z)) nodes), " *
          "run B $(length(run_b.poles)) poles " *
          "($(length(run_b.walk.visited_z)) nodes); " *
          "$(length(m.pairs)) matched, " *
          "disagreement median " *
          "$(round(m.median_disagree; digits = 4)) / max " *
          "$(round(m.max_disagree; digits = 4)); " *
          "single-run-only $(length(m.only_a)) (A) + " *
          "$(length(m.only_b)) (B)"
    return (walk_a = run_a.walk, walk_b = run_b.walk,
            poles_a = run_a.poles, poles_b = run_b.poles,
            match = m, message = msg)
end
