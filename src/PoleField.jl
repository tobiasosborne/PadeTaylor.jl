"""
    PadeTaylor.PoleField

Pole extraction from a solved path-network.

## Why this module exists

FW 2011 Figures 4.7 and 4.8 are *pole-location plots* — scatter plots
of where the Painlevé-I solution's poles sit in the complex plane.
The `PathNetwork` driver already carries everything needed to draw
them; it just never reads it back out.

`extract_poles` accepts two solution sources, because two solvers
carry a per-node Padé store: a `PathNetworkSolution` (the FW Fig
4.7/4.8 case — a fan of overlapping nodes) and a `PadeTaylorSolution`
(a single `solve_pade` trajectory — a chain of segments).  The
extraction procedure is identical for both; only the per-node arrays
they expose differ, and the shared `_extract_poles_core` below is
where the logic actually lives.

At every visited node `z_v` the path-network stores one local Padé
approximant `P(t) = N(t) / D(t)` in the rescaled variable
`t = (z - z_v) / h` (see `PathNetwork`'s module docstring and the
`t = h'/h` convention documented in `PadeStepper`). A rational
function's poles are exactly the zeros of its denominator, so the
poles of the solution *as seen from node `k`* are

    z = z_v + h · t*,    for every root  t*  of  D(t).

`extract_poles` takes the union of these candidates over the whole
visited tree, discards the roots a local Padé cannot be trusted to
place, and clusters what survives — many neighbouring nodes each "see"
the same physical pole, and we want one estimate per pole, taken from
the node that sees it best.

## Which roots to trust

A degree-`ν` denominator has `ν` roots; most are not physical poles of
the solution:

  - **Far roots.** A local `(m, n)` Padé built from a Taylor jet of
    radius `~h` reliably locates only singularities within a few `h`
    of its centre. A root is kept iff its **z-plane distance** from the
    birthing node, `|h·t*|`, is within `radius_z = radius_t · h_max`,
    where `h_max = maximum(|hs|)` is the walk's step *ceiling* — a
    fixed, `h`-INDEPENDENT scale. Roots beyond that are extrapolation
    artefacts and are dropped. The earlier filter kept roots by
    `|t*| ≤ radius_t` (a window of a fixed number of *per-node* steps
    `h`); under a varying / adaptive `h` that silently discards a pole
    at a fixed z-distance as soon as `h` shrinks (`|t*| = D/h` grows
    past the window), emptying the pole field. Filtering on z-plane
    distance against the `h`-independent `h_max` cures that and is
    genuinely scale-covariant; for a near-uniform-`h` walk
    (`h ≈ h_max`) it reduces exactly to the legacy `|t*| ≤ radius_t`.
    This is the §S7 scale-covariance fix ported from the vector twin
    `VectorPoleField.extract_poles_shared_q` (ADR-0026 Amendment 6/7;
    `src/VectorPoleField.jl:49-97`).

  - **Froissart doublets.** A near-singular Padé linear system — the
    Float64 classical-Toeplitz path has no SVD rank guard (ADR-0005) —
    can emit spurious pole–zero pairs: a denominator root `t*` that is
    *also* very nearly a numerator root, so the rational function has
    negligible residue there. The residue `N(t*) / D'(t*)` of a
    genuine pole is `O(1)` on the solution's scale; a Froissart
    doublet's residue sits at the noise floor. Roots whose residue
    magnitude is below `min_residue` are dropped.

  - **Node-local artefacts in general.** The decisive test is
    *cross-node agreement*. A physical pole is seen — at different
    `t*` — by every node within range of it; FW's Fig 4.7 is itself a
    "composite of independent runs". A spurious root, whatever its
    origin, is local to one node's linear system and does not recur at
    the same `z` from independent nodes. So a cluster of candidate
    poles is accepted as physical only when at least `min_support`
    *distinct* nodes independently place a root there.

The node nearest a physical pole in the z-plane places it most
accurately. Clustering is therefore *greedy in increasing z-plane
node-to-pole distance* `|h·t*|`: the first (best-placed) candidate to
land in a cluster becomes its representative, and every later
candidate — from a node that sees the pole less well — only adds
cross-node support. The earlier code ordered by `|t*|`, the
*per-node-step-relative* distance; that is the same scale-fixing
heresy as the old far-root filter — under a varying `h` a node closer
in z but with a smaller `h` carries a *larger* `|t*| = z_dist/h`, so
`|t*|`-ordering crowns a coarse far node over a fine near one and the
representative drifts to the worse-placed sighting. Ordering by the
z-plane distance `|h·t*|` is the genuinely scale-covariant "which node
is closest" measure (§S7; `src/VectorPoleField.jl:223-242`).

A second-order pole (the Weierstrass-℘ test problem, and the Painlevé
transcendents themselves, have double poles) is a double root of `D`.
A finite-order Padé splits it into two near-coincident simple roots
straddling the true location; the split shrinks with `|t*|` and is
already sub-`cluster_atol` for the well-placed nodes that win the
representative slot, so the doublet collapses into its own cluster
without special handling.

## Ground truth

Verified against the equianharmonic Weierstrass-℘ test problem of
FW 2011 §5.1.1 — `u'' = 6u²`, the PI equation with the `z`-term
dropped. Its solution `u(z) = ℘(z + c₁; 0, c₂)` has analytically known
second-order poles on a rhombic lattice; `test/polefield_test.jl`
carries the lattice formula and the pinned tolerances.

References:
  - `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:147` — Fig 4.7/4.8 are
    pole-location composites of independent path-network runs;
    `:281-318` — the Weierstrass-℘ test problem and its real-axis
    pole spacing `x = 1 + 2ωk`.
  - `docs/adr/0004-path-network-architecture.md` — the per-node Padé
    store this module reads back.
"""
module PoleField

using Polynomials: Polynomial, roots, derivative
using ..PathNetwork: PathNetworkSolution
using ..Problems:    PadeTaylorSolution

export extract_poles

# -----------------------------------------------------------------------------
# Shared core
# -----------------------------------------------------------------------------
# Both the path-network and the single-trajectory `extract_poles` methods do
# the same thing: walk a sequence of node-local Padé approximants, each with
# its own centre `z_ctr` and canonical step `h`, take the denominator roots
# `t*` mapped to the `z`-plane (`z = z_ctr + h·t*`), filter the untrustworthy
# ones, and cluster what survives. The only difference between the two callers
# is *where the per-node arrays live* — `visited_*` on a `PathNetworkSolution`,
# `z` / `pade` / `h` on a `PadeTaylorSolution` (a single `solve_pade`
# trajectory carries exactly the same per-segment Padé store). So the logic
# lives once, here, and the two public methods are thin adapters that hand it
# the right three arrays. `centers` may be longer than `pades` / `hs` (a
# `PadeTaylorSolution` stores `n+1` breakpoints but `n` segments); iteration is
# over `eachindex(hs)`, so the trailing breakpoint is harmlessly ignored.
function _extract_poles_core(centers, pades, hs;
                             radius_t::Real, min_residue::Real,
                             cluster_atol::Real, min_support::Integer)
    RT     = real(float(eltype(centers)))
    CT     = Complex{RT}
    radius = RT(radius_t)
    minres = RT(min_residue)
    catol  = RT(cluster_atol)

    # Far-root filter on a SCALE-STABLE z-radius (ADR-0026 Amendment 6/7
    # §S7; ported from the vector twin `VectorPoleField.extract_poles_shared_q`,
    # `src/VectorPoleField.jl:259-261`).  The scale is `h_max` — the
    # walk's step *ceiling* — recovered as the largest per-node step
    # `|hs|`.  `h_max` is `h`-INDEPENDENT: a varying-`h` (stitched /
    # adaptive) walk shrinks individual `hs[k]`, but the ceiling is fixed
    # by the coarsest segment, so `maximum(abs, hs)` recovers it.  For an
    # empty `hs` there are no nodes — `radius_z` is unused (the loop is
    # empty) — so a `zero(RT)` placeholder is harmless (mirrors the
    # vector code's empty-`visited_h` guard).
    h_max    = isempty(hs) ? zero(RT) : RT(maximum(abs, hs))
    radius_z = radius * h_max

    # (pole_z, z_dist, node-index) candidates gathered from every node.
    # `z_dist = |h·t*|` is the z-plane distance from the birthing node to
    # the pole — the scale-covariant "how close is this node" measure the
    # greedy clustering orders by (see below, and the vector twin's
    # docstring "Why the representative is the z-plane-closest node").
    candidates = Tuple{CT, RT, Int}[]
    for k in eachindex(hs)
        P     = pades[k]
        z_ctr = centers[k]
        h     = hs[k]

        # A constant denominator (b == [1]) has no roots — the local
        # Padé is a polynomial, no poles to report from this node.
        length(P.b) ≥ 2 || continue

        D  = Polynomial(P.b)
        N  = Polynomial(P.a)
        Dp = derivative(D)

        for t in roots(D)
            # Far-root filter — ADR-0026 §S7.  Accept by the root's
            # z-plane DISTANCE `|h·t*|` against the scale-stable z-radius
            # `radius_z = radius_t·h_max`, NOT by `|t*|` (a window of a
            # fixed number of *per-node* steps `h`, which silently
            # discards a fixed-z-distance pole as soon as the adaptive
            # `h` shrinks: `|t*| = D/h` grows past the window — measured
            # to empty pole fields on the vector side, ADR-0026 Amd 6).
            # The `h` factor is load-bearing: t lives in t = (z - z_ctr)/h,
            # so the pole sits at z-distance `|h·t*|`.  For a default
            # near-uniform-`h` walk every `h ≈ h_max`, so the test reduces
            # to the legacy `|t*| ≤ radius_t` — backward compatible.
            z_dist = abs(h * t)
            z_dist ≤ radius_z || continue
            res = N(t) / Dp(t)                    # Padé residue at t*
            (isfinite(res) && abs(res) ≥ minres) || continue   # Froissart
            push!(candidates, (CT(z_ctr + h * t), RT(z_dist), k))
        end
    end

    # Greedy clustering in increasing z-plane node-to-pole distance
    # `z_dist`: the best-placed candidate — the one whose birthing node is
    # *physically closest* to the pole in the z-plane — lands first and
    # becomes the cluster representative; later candidates only contribute
    # cross-node support.  Ordering by `z_dist` (not the per-node-step-
    # relative `|t*|`) is the §S7 scale-covariance fix: with a varying
    # `h`, a node *closer* in z but with a *smaller* `h` carries a
    # *larger* `|t*| = z_dist/h`, so the old `|t*|` order crowned a coarse
    # far node over a fine near one — the near node's placement is the
    # accurate one (`src/VectorPoleField.jl:223-242`).  A cluster is a
    # physical pole only when ≥ min_support distinct nodes land a root in
    # it — node-local artefacts never accrue cross-node support.
    sort!(candidates; by = c -> c[2])
    reps    = CT[]                    # cluster representatives
    support = Vector{Set{Int}}()      # distinct source nodes per cluster
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

"""
    extract_poles(sol::PathNetworkSolution{T};
                  radius_t     = 5.0,
                  min_residue  = 1.0e-8,
                  cluster_atol = 1.0e-1,
                  min_support  = 3) -> Vector{Complex{T}}

Pole locations of the solution carried by `sol`, in the `z`-plane.

For every visited node the roots of the stored Padé denominator are
mapped back to the `z`-plane (`z = z_v + h · t*`), filtered, and then
clustered:

  - `radius_t`     — the far-root filter, in units of "steps at the
                     walk's step *ceiling* `h_max = maximum(|h|)`".  A
                     root `t*` from a node with step `h` is kept iff its
                     z-plane distance from the node satisfies
                     `|h·t*| ≤ radius_t · h_max`; a local Padé does not
                     reliably place distant singularities (default `5.0`,
                     a few canonical steps at the ceiling scale).  This is
                     the §S7 scale-covariance fix (see the module
                     docstring "Far roots"); for a near-uniform-`h` walk
                     it reduces exactly to the legacy `|t*| ≤ radius_t`.
  - `min_residue`  — drop Froissart doublets: keep a root only when the
                     Padé residue `|N(t*) / D'(t*)| ≥ min_residue`.
  - `cluster_atol` — surviving roots within this distance (in `z`) of
                     each other are the same physical pole; the cluster
                     must be smaller than the inter-pole spacing and
                     larger than a split second-order-pole doublet.
  - `min_support`  — a cluster is reported as a pole only when at least
                     this many *distinct* visited nodes independently
                     placed a root in it. This is the load-bearing
                     spurious-pole filter (see the module docstring);
                     pass `min_support = 1` to disable it, e.g. when
                     reading poles off a single-node network.

Each reported pole is the cluster representative — the candidate placed
by the node closest to it *in the z-plane* (the smallest node-to-pole
z-distance `|h·t*|`; §S7, not the per-node-step-relative `|t*|`).
Returns one `Complex{T}` per physical pole, in order of discovery.
"""
extract_poles(sol::PathNetworkSolution{T};
              radius_t::Real       = 5.0,
              min_residue::Real    = 1.0e-8,
              cluster_atol::Real   = 1.0e-1,
              min_support::Integer = 3) where {T} =
    _extract_poles_core(sol.visited_z, sol.visited_pade, sol.visited_h;
                        radius_t, min_residue, cluster_atol, min_support)

"""
    extract_poles(sol::PadeTaylorSolution;
                  radius_t     = 5.0,
                  min_residue  = 1.0e-8,
                  cluster_atol = 1.0e-1,
                  min_support  = 1) -> Vector{<:Complex}

Pole locations of a single-trajectory `solve_pade` result, in the
`z`-plane.

A `PadeTaylorSolution` carries the same per-segment Padé store a
`PathNetworkSolution` does — `sol.pade[k]` is the local approximant
covering segment `k`, centred at `sol.z[k]` with canonical step
`sol.h[k]` — so its poles are extracted by exactly the same
denominator-root / residue-filter / clustering procedure (see the
module docstring); only the array names differ.

The one changed default is `min_support = 1`.  The path-network method
defaults to `3` because its raison d'être is *cross-node* agreement —
a physical pole is seen by many overlapping nodes.  A single
`solve_pade` trajectory is a chain, not a fan: a given pole is
typically bracketed by only one or two consecutive segments, so
demanding three independent sightings would discard every real pole.
A single-trajectory solve is the "single-node network" case the
path-network docstring already flags for `min_support = 1`.  Raise it
only if the trajectory genuinely revisits a region from independent
segments.

`radius_t`, `min_residue`, `cluster_atol` carry the same meaning as in
the `PathNetworkSolution` method.
"""
extract_poles(sol::PadeTaylorSolution;
              radius_t::Real       = 5.0,
              min_residue::Real    = 1.0e-8,
              cluster_atol::Real   = 1.0e-1,
              min_support::Integer = 1) =
    _extract_poles_core(sol.z, sol.pade, sol.h;
                        radius_t, min_residue, cluster_atol, min_support)

end # module PoleField
