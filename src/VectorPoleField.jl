"""
    PadeTaylor.VectorPoleField

Pole extraction from a solved **vector** path-network — the shared-`Q`
analogue of v0.1's `PoleField.extract_poles`, and the tool the v0.2
`A_4⁽¹⁾` / `P_I⁽²⁾` pole-field scatter figures (beads V8a / V8b) read
their pole locations from.

## Why the shared `Q` makes this *cleaner* than v0.1's `PoleField`

v0.1's `PoleField` (read `src/PoleField.jl`, "Why this module exists")
extracts poles from a *scalar* path-network: each visited node stores
one `(m,n)` Padé `P(t)/D(t)`, and the node's poles are the roots of
its single denominator `D`.  Lifting that to a vector ODE the naïve
way — one independent scalar Padé per component — gives **`d`
denominators per node**, hence `d` inconsistent root sets that agree
only up to numerical noise (ADR-0019).

The v0.2 vector pipeline avoids this entirely.
`VectorStepper.vector_pade_step_with_pade!` fits all `d` components
with a **shared** denominator: `d` numerators `P₁,…,P_d` over **one**
denominator `Q` (a type-II Hermite–Padé approximant — ADR-0019;
Mano–Tsuda 2017).  A vector meromorphic solution's components all blow
up at the *same* movable poles, so the single shared `Q` is the
honest representation: **every component's poles are the roots of one
polynomial**.  Pole extraction therefore roots *one* `Q` per node, not
`d` per-component denominators — strictly cleaner than the scalar
case, with no per-component reconciliation step.

## The extraction procedure (mirrors v0.1 `PoleField._extract_poles_core`)

For every visited node `k`:

  1. Root the stored shared denominator `Q` (a polynomial, low-to-high
     coefficients) with `Polynomials.roots` — exactly v0.1's
     `PoleField` call.
  2. Map each root `t*` back to the z-plane by `z = z_node + t*·h_node`
     — the local approximant lives in the rescaled variable
     `t = (z - z_node)/h_node`, so a denominator root `t*` is a pole at
     `z_node + t*·h_node` (`src/PoleField.jl`, "z = z_v + h·t*").
  3. Discard far roots — keep only `|t*| ≤ radius_t`; a local Padé
     built from a Taylor jet of radius `~|h|` reliably places only
     singularities within a few steps of its centre.

Then **cluster across nodes** (the cross-node support filter v0.1's
`_extract_poles_core` uses): many neighbouring nodes each "see" the
same physical pole; greedy clustering in increasing `|t*|` makes the
best-placed (smallest-`|t*|`) candidate the cluster representative,
and a cluster is reported as a physical pole only when at least
`min_support` *distinct* nodes independently place a root in it.  A
node-local linear-system artefact does not recur at the same `z` from
independent nodes, so the cross-node filter is the load-bearing
spurious-pole guard — exactly as in the scalar `PoleField`.

## Fail-fast / fail-soft

An empty solution (no visited nodes) yields an empty pole list — not
an error; "no nodes ⇒ no poles" is a true answer, not a breakdown.  A
constant denominator (`Q = [1]`, no roots) contributes nothing from
that node.  `Polynomials.roots` failure on a degenerate `Q`
propagates unchanged.

## References

  - `src/PoleField.jl` — the v0.1 scalar pole extractor whose
    denominator-root / z-plane-map / cross-node-cluster procedure this
    module mirrors (`_extract_poles_core`, `:119-169`).
  - `src/VectorPathNetwork.jl` — produces the `VectorPathNetworkSolution`
    this module consumes; the per-node shared-`Q` store.
  - `docs/adr/0019-shared-denominator-pade.md` — the shared `Q` whose
    roots are *every* component's pole field.
  - `docs/v0p2_plan.md` — §"The 10 phases", V7 row.
"""
module VectorPoleField

using Polynomials:        Polynomial, roots
using ..VectorPathNetwork: VectorPathNetworkSolution

export extract_poles_shared_q

# --- Scale-covariant cluster tolerance (ADR-0026 Amendment 3 §S4) ---------
#
# `CLUSTER_FRAC` is the *one* free knob of the scale-derived clustering
# mode: two pole candidates born from nodes `i`, `j` are merged into one
# reported pole iff their z-plane separation is `≤ CLUSTER_FRAC · (the
# smaller of the two nodes' local steps)` — i.e. within a fraction of a
# local step.
#
# Why a fraction of the local step `h`, not an absolute length.  An
# absolute cluster tolerance is a *scale heresy* (project memory
# `scale-covariance-core-principle`; ADR-0026 Amendment 3 §S4): in a
# varying-scale pole field it wrongly *merges* distinct poles where the
# field is dense and *splits* one pole where it is sparse.  After step
# S2 each node's step `visited_h[i]` is *scale-adaptive* — sized near
# the local pole-free disc, hence near the local scale — so a tolerance
# tied to `visited_h` tracks the local scale for free.
#
# Why `CLUSTER_FRAC ≈ 0.4`.  Two regimes must be separated:
#   * Distinct physical poles are ~one pole-spacing apart.  ADR-0026
#     Amendment 3 measured the honest pole-free disc at ≈ 0.106 × the
#     local pole spacing, and S2 sizes `h` near that honest disc — so a
#     pole spacing is *several* `h`.  A tolerance of `0.4·h` is well
#     under one spacing and never collapses two distinct poles.
#   * Numerical duplicates of ONE pole, seen by two neighbouring nodes,
#     are separated only by walk noise — far below `h`.  A tolerance of
#     `0.4·h` (≈ half a local step) comfortably swallows that noise.
# `0.4` sits inside the reasoned 0.3–0.5 band: large enough to merge the
# noise-separated duplicates, small enough to never merge two distinct
# poles.  TUNABLE — the one knob; revisit if a pole census shows the
# extractor over- or under-merging.
const CLUSTER_FRAC = 0.4

"""
    extract_poles_shared_q(sol::VectorPathNetworkSolution{T};
                           radius_t     = 5.0,
                           cluster_atol = nothing,
                           min_support  = 2) -> Vector{Complex{T}}

Pole locations of the vector system carried by `sol`, in the z-plane.

For every visited node the roots of the stored **shared** denominator
`Q` are mapped back to the z-plane (`z = z_node + t*·h_node`), filtered
by `radius_t`, and clustered across nodes:

  - `radius_t`     — keep only roots with `|t*| ≤ radius_t`; a local
                     shared-`Q` Padé does not reliably place distant
                     singularities (default `5.0`, a few steps).  `t`
                     is the *dimensionless* rescaled variable `Δz/h`,
                     and post-S2 `h` is itself scale-derived (sized
                     near the local pole-free disc), so `radius_t` as
                     "a few steps" is **already scale-covariant** —
                     `radius_t · h` is a genuine local-scale z-distance.
                     It needs no change post-S2 (ADR-0026 Amendment 3
                     §S4): the audit-flag on `radius_t` is discharged
                     because the heresy was an *absolute* length, and
                     `radius_t` was never one.
  - `cluster_atol` — the across-node merge tolerance.  Two regimes:
      * `nothing` (default) — **scale-derived** mode.  Two candidates
        from nodes `i`, `j` merge iff their z-plane separation is
        `≤ CLUSTER_FRAC · min(|visited_h[i]|, |visited_h[j]|)` — a
        fraction of the *local* step.  Because `visited_h` is itself
        scale-adaptive after step S2, this tolerance tracks the local
        pole scale and never wrongly merges/splits (ADR-0026
        Amendment 3 §S4; see `CLUSTER_FRAC`).
      * a `Real` — **legacy absolute** mode.  Surviving roots within
        this fixed z-plane distance of each other are the same
        physical pole.  Kept for backward compatibility — existing
        callers that pin an explicit `cluster_atol` are byte-identical.
  - `min_support`  — a cluster is reported only when at least this many
                     *distinct* visited nodes independently place a
                     root in it — the cross-node spurious-pole filter
                     (see the module docstring).  Pass `1` to disable
                     it (e.g. a single-node network).

Because the denominator `Q` is *shared* across all `d` components, the
returned poles are every component's poles at once — a single
consistent estimate of the vector system's movable singularities.
Each reported pole is its cluster representative — the candidate seen
at the smallest `|t*|`, i.e. by the closest node.  Returns one
`Complex{T}` per physical pole, in order of discovery.
"""
function extract_poles_shared_q(sol::VectorPathNetworkSolution{T};
                                radius_t::Real               = 5.0,
                                cluster_atol::Union{Nothing,Real} = nothing,
                                min_support::Integer         = 2) where {T}
    C      = Complex{T}
    radius = T(radius_t)

    # Legacy absolute mode iff the caller pinned a `Real`; the
    # scale-derived mode (the default) is signalled by `nothing`.
    legacy_atol = cluster_atol === nothing ? nothing : T(cluster_atol)

    # (pole_z, |t*|, node-index, local-h) candidates gathered from every
    # node.  `local_h` — the magnitude of the birthing node's step — is
    # the per-candidate local scale the scale-derived cluster tolerance
    # is a fraction of (ADR-0026 Amendment 3 §S4).
    candidates = Tuple{C, T, Int, T}[]
    for k in eachindex(sol.visited_denominator)
        Q      = sol.visited_denominator[k]
        z_node = sol.visited_z[k]
        h_node = sol.visited_h[k]

        # A constant denominator (Q == [1]) has no roots — the local
        # approximant is a polynomial, no poles to report from this node.
        length(Q) ≥ 2 || continue

        for t in roots(Polynomial(Q))
            t_C = C(t)
            abs(t_C) ≤ radius || continue          # far-root artefact
            # Map the rescaled-variable root to the z-plane.  The h
            # factor is load-bearing: t lives in t = (z - z_node)/h.
            push!(candidates,
                  (C(z_node + h_node * t_C), T(abs(t_C)), k, T(abs(h_node))))
        end
    end

    # Greedy clustering in increasing |t*|: the best-placed (smallest
    # |t*|) candidate to land in a cluster becomes its representative;
    # later candidates only add cross-node support.  A cluster is a
    # physical pole only when ≥ min_support distinct nodes land a root
    # in it — node-local artefacts never accrue cross-node support.
    #
    # The merge test compares a candidate against each representative.
    # In scale-derived mode the tolerance for a (candidate, rep) pair is
    # `CLUSTER_FRAC · min(local-h of the two)` — the conservative
    # smaller-step choice, so a candidate born at a fine local scale is
    # never absorbed into a rep born at a coarse one (and vice versa).
    # In legacy mode it is the fixed absolute `legacy_atol`.
    sort!(candidates; by = c -> c[2])
    reps     = C[]
    support  = Vector{Set{Int}}()
    rep_h    = T[]                       # local-h of each rep's birth node
    for (p, _, k, h_cand) in candidates
        j = findfirst(eachindex(reps)) do r
            tol = legacy_atol === nothing ?
                  T(CLUSTER_FRAC) * min(h_cand, rep_h[r]) :
                  legacy_atol
            abs(p - reps[r]) ≤ tol
        end
        if j === nothing
            push!(reps, p)
            push!(support, Set{Int}((k,)))
            push!(rep_h, h_cand)
        else
            push!(support[j], k)
        end
    end
    return [reps[j] for j in eachindex(reps) if length(support[j]) ≥ min_support]
end

end # module VectorPoleField
