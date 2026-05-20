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

"""
    extract_poles_shared_q(sol::VectorPathNetworkSolution{T};
                           radius_t     = 5.0,
                           cluster_atol = 1.0e-1,
                           min_support  = 2) -> Vector{Complex{T}}

Pole locations of the vector system carried by `sol`, in the z-plane.

For every visited node the roots of the stored **shared** denominator
`Q` are mapped back to the z-plane (`z = z_node + t*·h_node`), filtered
by `radius_t`, and clustered across nodes:

  - `radius_t`     — keep only roots with `|t*| ≤ radius_t`; a local
                     shared-`Q` Padé does not reliably place distant
                     singularities (default `5.0`, a few steps).
  - `cluster_atol` — surviving roots within this z-plane distance of
                     each other are the same physical pole; the
                     cluster radius must be below the inter-pole
                     spacing.
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
                                radius_t::Real       = 5.0,
                                cluster_atol::Real   = 1.0e-1,
                                min_support::Integer = 2) where {T}
    C      = Complex{T}
    radius = T(radius_t)
    catol  = T(cluster_atol)

    # (pole_z, |t*|, node-index) candidates gathered from every node.
    candidates = Tuple{C, T, Int}[]
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
            push!(candidates, (C(z_node + h_node * t_C), T(abs(t_C)), k))
        end
    end

    # Greedy clustering in increasing |t*|: the best-placed (smallest
    # |t*|) candidate to land in a cluster becomes its representative;
    # later candidates only add cross-node support.  A cluster is a
    # physical pole only when ≥ min_support distinct nodes land a root
    # in it — node-local artefacts never accrue cross-node support.
    sort!(candidates; by = c -> c[2])
    reps    = C[]
    support = Vector{Set{Int}}()
    for (p, _, k) in candidates
        j = findfirst(r -> abs(p - r) ≤ catol, reps)
        if j === nothing
            push!(reps, p)
            push!(support, Set{Int}((k,)))
        else
            push!(support[j], k)
        end
    end
    return [reps[j] for j in eachindex(reps) if length(support[j]) ≥ min_support]
end

end # module VectorPoleField
