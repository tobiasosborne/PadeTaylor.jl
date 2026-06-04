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
  3. Discard far roots — keep only roots whose **z-plane distance**
     `h_node·|t*|` from the node is within a scale-stable z-radius; a
     local Padé built from a Taylor jet of radius `~|h|` reliably
     places only singularities within a few steps of its centre.  The
     z-radius is `radius_t · h_max` — a few steps measured at the
     walk's step *ceiling* `h_max`, NOT at the per-node adaptive `h`
     (see "The `radius_t` scale-covariance fix" below).

## The `radius_t` scale-covariance fix (ADR-0026 Amendment 6 §S7)

The far-root filter used to keep roots by `|t*| ≤ radius_t` — `|t*|`
being the root magnitude in the *rescaled* variable `t = Δz/h_node`.
That is a window of "a fixed number of *steps*."  It is a **scale-
fixing heresy**: when the adaptive step `h_node` shrinks (a smaller
`SAFETY`, a denser pole pocket), the *same physical pole* at a fixed
z-distance `D` from a node maps to a *larger* `|t*| = D/h_node`, falls
outside the `radius_t` window, and is discarded.  Fewer nodes then
root it, and the `min_support ≥ 2` cross-node filter empties the pole
field — measured: lowering `SAFETY` 0.25→0.10 emptied four `src/`
test fields (ADR-0026 Amendment 6).

The fix accepts a root by its **z-plane distance** `h_node·|t*|`
against a z-radius tied to a **scale-stable** quantity that does *not*
shrink with the per-node `h`: `h_max = maximum(|visited_h|)`, the
walk's step ceiling.  `h_max` is `h`-independent — the walk's `h`
ceiling is set once by the caller and the root node always carries
exactly `C(h_max)` (`VectorPathNetwork.jl`), so `maximum(|visited_h|)`
recovers it exactly.  The acceptance test becomes

    h_node · |t*|  ≤  radius_t · h_max

— "the pole is within `radius_t` steps *at the ceiling scale*."  This
is genuinely scale-covariant: it is a pure z-plane-distance test and
contains no absolute length (`radius_t` is dimensionless, `h_max`
carries the problem's scale, so `radius_t·h_max` rescales with the
problem — poles spaced 0.01 or 100 are treated identically).  For a
default near-uniform-`h` walk every `h_node ≈ h_max`, so the test
reduces to the old `|t*| ≤ radius_t` — **legacy behaviour is exactly
preserved** for the default-`h` case.  When `h_node` shrinks, the
`t*`-window correspondingly widens (`|t*| ≤ radius_t·h_max/h_node`),
so the same physical pole keeps being rooted and `min_support` voting
stays healthy.

Then **cluster across nodes** (the cross-node support filter v0.1's
`_extract_poles_core` uses): many neighbouring nodes each "see" the
same physical pole; greedy clustering in increasing **z-plane
node-to-pole distance** makes the best-placed candidate — the one
whose birthing node is physically closest to the pole — the cluster
representative, and a cluster is reported as a physical pole only when
at least `min_support` *distinct* nodes independently place a root in
it.  (v0.1's scalar `PoleField` ordered by `|t*|`; the vector walk's
per-node adaptive `h` makes `|t*|` an unreliable proxy for node
proximity — see "Why the representative is the z-plane-closest node"
in `extract_poles_shared_q`'s docstring.)  A node-local linear-system
artefact does not recur at the same `z` from independent nodes, so the
cross-node filter is the load-bearing spurious-pole guard — exactly as
in the scalar `PoleField`.

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
using ..SharedPadeDefect: shared_q_residue

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
by their z-plane distance, and clustered across nodes:

  - `radius_t`     — the far-root filter, in units of "steps at the
                     walk's step *ceiling* `h_max`".  A root `t*` from
                     a node with step `h_node` is kept iff its z-plane
                     distance from the node satisfies

                         h_node · |t*|  ≤  radius_t · h_max

                     where `h_max = maximum(|visited_h|)` is the
                     walk's step ceiling — a fixed, `h`-independent
                     scale.  A local shared-`Q` Padé does not reliably
                     place distant singularities (default `5.0`, a few
                     steps at the ceiling scale).

                     This is the ADR-0026 Amendment 6 §S7 scale-
                     covariance fix.  The *old* filter kept roots by
                     `|t*| ≤ radius_t` — a window of a fixed number of
                     *per-node* steps `h_node` — which discarded a
                     fixed-z-distance pole as soon as the adaptive
                     `h_node` shrank (`|t*| = D/h_node` grew past the
                     window).  Filtering on the z-plane distance
                     against the `h`-independent `h_max` cures that:
                     the test is a pure z-distance comparison with no
                     absolute length, so it is genuinely scale-
                     covariant, and for a default near-uniform-`h`
                     walk (`h_node ≈ h_max`) it reduces *exactly* to
                     the legacy `|t*| ≤ radius_t` — backward
                     compatible.  See the module docstring, "The
                     `radius_t` scale-covariance fix".
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
Each reported pole is its cluster representative — the candidate
placed by the **node closest to it in the z-plane** (the smallest
node-to-pole z-distance `h_node·|t*|`).  Returns one `Complex{T}` per
physical pole, in order of discovery.

## Why the representative is the z-plane-closest node, not min-`|t*|`

The greedy clustering picks each cluster's representative as the
*first* candidate to land in it; the candidates are visited in an
order that must put the **best-placed** one first.  v0.1's scalar
`PoleField` ordered by `|t*|` — under a *uniform* step `h` the
smallest-`|t*|` root is the one from the node physically closest to
the pole, hence the most reliably placed (a Padé places a near pole
better than a far one).  But the vector walk has a **per-node adaptive
`h`** (ADR-0026 Amendment 4): a node *closer* in the z-plane but with
a *smaller* `h` carries a *larger* `|t*| = z_dist/h`.  Ordering by
`|t*|` would then crown a coarse-`h` *far* node over a fine-`h` *near*
one — and the near node's placement is the accurate one (measured: a
℘-pole walk placed the pole at `4·10⁻⁵` from its min-`|t*|` node but
`6·10⁻⁷` from the z-plane-closest node — ADR-0026 Amendment 6 §S7).
`|t*|` is a per-node-step-relative quantity, not scale-stable — the
same heresy as the old `radius_t` filter.  The clustering therefore
orders candidates by the **z-plane node-to-pole distance**
`h_node·|t*|`, the genuinely scale-covariant "which node is closest"
measure.
"""
function extract_poles_shared_q(sol::VectorPathNetworkSolution{T};
                                radius_t::Real               = 5.0,
                                cluster_atol::Union{Nothing,Real} = nothing,
                                min_support::Integer         = 2,
                                min_residue::Real            = 1.0e-4) where {T}
    C = Complex{T}

    # The far-root filter is a z-plane-distance test against a
    # scale-STABLE z-radius (ADR-0026 Amendment 6 §S7).  The scale is
    # `h_max` — the walk's step ceiling — recovered as the largest
    # per-node step `|visited_h|`.  `h_max` is `h`-independent: the
    # caller sets the ceiling once and the root node always carries
    # exactly `C(h_max)`, so the maximum recovers it.  A root `t*` from
    # a node with step `h_node` is kept iff `h_node·|t*| ≤ radius_z`.
    # For an empty solution there are no nodes — `radius_z` is unused
    # (the loop below is empty) — so a `0` placeholder is harmless.
    h_max    = isempty(sol.visited_h) ? zero(T) :
               maximum(abs, sol.visited_h)
    radius_z = T(radius_t) * h_max

    # Legacy absolute mode iff the caller pinned a `Real`; the
    # scale-derived mode (the default) is signalled by `nothing`.
    legacy_atol = cluster_atol === nothing ? nothing : T(cluster_atol)

    # (pole_z, z_dist, node-index, local-h) candidates gathered from
    # every node.  `z_dist = |h_node·t*|` is the z-plane distance from
    # the birthing node to the pole — the scale-covariant "how close is
    # this node" measure the greedy clustering orders by (see the
    # docstring, "Why the representative is the z-plane-closest node").
    # `local_h` — the magnitude of the birthing node's step — is the
    # per-candidate local scale the scale-derived cluster tolerance is a
    # fraction of (ADR-0026 Amendment 3 §S4).
    candidates = Tuple{C, T, Int, T}[]
    for k in eachindex(sol.visited_denominator)
        Q      = sol.visited_denominator[k]
        z_node = sol.visited_z[k]
        h_node = sol.visited_h[k]

        # A constant denominator (Q == [1]) has no roots — the local
        # approximant is a polynomial, no poles to report from this node.
        length(Q) ≥ 2 || continue
        nums_k = sol.visited_numerators[k]

        for t in roots(Polynomial(Q))
            # Froissart-doublet filter (ADR-0028 dual-construction): a cell-B
            # near-cancelling pole–zero pair is not a real pole.  Skip a root
            # whose shared-Q residue `max_c |P_c(t)/Q'(t)|` is below
            # `min_residue` (calibrated 1e-4; genuine poles ≥ 2.1, doublets
            # ≤ 5.5e-8 — probe `external/probes/adr0028-froissart-consumer/`).
            shared_q_residue(nums_k, Q, t) ≥ min_residue || continue
            t_C = C(t)
            # Far-root filter — ADR-0026 Amendment 6 §S7.  Accept by
            # the root's z-plane DISTANCE `|h_node·t*|` against the
            # scale-stable z-radius `radius_z = radius_t·h_max`, NOT by
            # `|t*|` (a window of a fixed number of per-node steps,
            # which discards a fixed-z-distance pole when `h_node`
            # shrinks).  The h factor is load-bearing: t lives in
            # t = (z - z_node)/h, so the pole sits at z-distance
            # |h_node·t*|.
            z_pole   = C(z_node + h_node * t_C)
            z_dist   = abs(z_pole - C(z_node))      # = |h_node·t*|
            z_dist ≤ radius_z || continue
            push!(candidates,
                  (z_pole, T(z_dist), k, T(abs(h_node))))
        end
    end

    # Greedy clustering in increasing z-plane node-to-pole distance
    # `z_dist`: the best-placed candidate — the one whose birthing node
    # is *physically closest* to the pole — lands first and becomes the
    # cluster representative; later candidates only add cross-node
    # support.  Ordering by `z_dist` (not the per-node-step-relative
    # `|t*|`) is the ADR-0026 Amendment 6 §S7 scale-covariance fix — see
    # the docstring, "Why the representative is the z-plane-closest
    # node".  A cluster is a physical pole only when ≥ min_support
    # distinct nodes land a root in it — node-local artefacts never
    # accrue cross-node support.
    #
    # The merge test compares a candidate against each representative.
    # In scale-derived mode the tolerance for a (candidate, rep) pair is
    # `CLUSTER_FRAC · min(local-h of the two)` — the conservative
    # smaller-step choice, so a candidate born at a fine local scale is
    # never absorbed into a rep born at a coarse one (and vice versa).
    # In legacy mode it is the fixed absolute `legacy_atol`.
    sort!(candidates; by = c -> c[2])    # c[2] = z_dist — closest first
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
