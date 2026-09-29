"""
    PadeTaylor.WindowedTiling

The **tiling geometry** beneath `WindowedComposite`: where the window
cores sit, which window owns a grid cell, how each window is seeded,
and the one degenerate configuration a caller must be told about.
`WindowedComposite.windowed_path_network_solve` drives these helpers;
nothing here solves an ODE.  The module is split out of
`WindowedComposite.jl` purely for the CLAUDE.md Rule 6 200-LOC cap
(bead `padetaylor-fqwf`, the commit that added the single-window
warning and pushed the composite file to 211 effective LOC) and has no
reverse dependency on it, so it loads first and `WindowedComposite`
`using`s it — the same one-way arrow as `VectorPathNetworkStage2` →
`VectorPathNetwork`.

## Ground truth — FW 2011's composite recipe

FW 2011 md:147 (`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:147`): each Fig 4.7 pole field "is
a 5×5 composite of 25 completely independent runs (each starting from
the same ICs at z = 0, and covering areas of the extent 20×20 in the
complex plane)"; Fig 4.8 "is a composite of 18 independent runs".  The
surrounding text (md:141-147) presents this as a parallelisation
convenience, but ADR-0034 (`docs/adr/0034-bounded-window-composite.md`)
shows it is load-bearing for correctness: the monolithic
`path_network_solve` walks one seed-dependent tree and the meeting
lines of its subtrees carry a path-dependence seam (worklogs 077/078).
Independent windows from the shared IC have no second subtree to
disagree with.  The three geometric primitives below are what "5×5
composite of independent runs" means concretely.

## The tiling geometry (`_tile_centers`)

Along each axis the core centres are laid down at `extent` spacing.
We use `n = ceil(span / extent)` cores (never fewer than one) and
centre the run of `n` cores on the axis midpoint, so the cores
overhang the domain *symmetrically* rather than all on one side.
The union of the cores `[c_k ± extent/2]` therefore covers `[lo, hi]`,
and because the cores are equal-sized and equally spaced, the Voronoi
partition of the plane by the centre set *is* the core tiling — every
domain point lies in exactly one core.  That identity is what lets the
composite be assembled by nearest-centre lookup instead of an explicit
rectangle membership test, and it is also why the solve extent (core
`± overlap`) can be discarded after compositing: ownership is decided
by the core, never by the overlap.

## Voronoi ownership (`_nearest_center`)

A grid cell belongs to the window whose centre is nearest.  We compare
*squared* distances — monotone in distance, so the argmin is the same,
and it saves a `sqrt` per candidate — with a strict `<` and a
lowest-index tie-break.  The tie-break is not cosmetic: a cell exactly
on a core boundary is equidistant from two centres, and a
non-deterministic choice there would make the composite depend on
iteration order.  With the strict `<`, the lower window index always
wins on the boundary, which matches the closed-on-the-left convention
of the tiling, and the assignment is reproducible bit-for-bit.

## Pole ownership: cluster first, then assign (`_own_poles`)

The same nearest-centre rule does NOT transfer verbatim to *poles*.  A
grid cell has one exact position, so `_nearest_center(cell)` is
single-valued.  A physical pole in the overlap of two windows has *two*
positions — each window's own independently-computed estimate — and
they differ in the last digits (measured 1e-8 … 6e-5 on the ℘ lattice,
worklog 084).  Asking "is this estimate in my core?" once per estimate
is a partition of the *plane*, not of the *pole set*: when the pole sits
on a core line the two estimates can land on opposite sides of it and
the pole is kept **twice** (each estimate is in its own window's core)
or **zero** times (each is in the other's core).  This is not a
hypothetical corner: a real-coefficient problem with real initial data
is symmetric about the real axis, its real-axis poles sit exactly on the
`Im z = 0` core line of any even tiling, and on the ℘ lattice over
`[-20,20]²` (2×2 windows) the per-estimate rule measured 4 duplicated
and 4 dropped poles out of 246 (bead `padetaylor-580u`, worklog 084).

`_own_poles` therefore decides ownership once per *physical pole*.  It
first links estimates from **different** windows that lie within
`atol` of each other (a sweep over real-part-sorted candidates, union-
find), then gives each linked group ONE canonical position — the mean
of its members — and ONE owner, `_nearest_center(mean)`.  Only the
owner's own estimate(s) are kept.  So a straddling pole is kept exactly
once (the deterministic tie-break settles a mean that lands exactly on
the line), and a pole seen by only one window is judged exactly as
before.  Same-window estimates are never linked: collapsing one
window's own near-coincident roots is `PoleField.extract_poles`'s job
(its `cluster_atol` / `merge_atol`), not the composite's.  `atol` must
sit far above the cross-window spread and far below the inter-pole
spacing; the caller defaults it to `extract_poles`' own `cluster_atol`
(the scale at which the extractor already calls two roots one pole).

## Per-window seeding (`_window_seed`)

Each window is solved at its own seed derived from the caller's global
`rng_seed` and the window index.  Two distinct global seeds must
re-randomise *every* window, otherwise the seed-invariance gate
(`test/field_seam_test.jl` FSEAM.1) could be passed by a pipeline that
silently ignored the seed.  We combine arithmetically,
`rng_seed * 1_000_003 + wi`, and deliberately NOT with `hash`: `hash`
is not stable across Julia versions, which would break the package's
bit-reproducibility contract without any test noticing.  The prime
multiplier is far above any realistic window count, so `(rng_seed, wi)`
pairs never collide for sane inputs.

## The single-window hazard (`_warn_single_window`)

If `extent ≥ span` on an axis, `_tile_centers` returns ONE centre on
that axis.  A 1×1 tile grid is not a composite at all — it is exactly
the monolithic, path-dependent solve the windowed API exists to
replace, and ADR-0034 measured that at `window_extent ≈ domain` the
composite is no better than monolithic (HALF=25: 91 % vs 93 %).  The
collapse is silent by construction (the code path is identical), so
the caller is told with a `@warn` naming the axis, the numbers, bead
`padetaylor-fqwf`, and the fix (reduce `window_extent`).  We warn
rather than throw because a single window is a *legal* request — a
domain genuinely smaller than one window, or a caller deliberately
reproducing the monolithic result through this API — so Rule 1's
"fail loud" is met by a log line the caller cannot miss without
refusing a well-posed solve.  `test/field_seam_test.jl` FSEAM.3 pins
both the warning and its absence on a genuine tiling.

## References

  - ADR-0034 — `docs/adr/0034-bounded-window-composite.md`.
  - FW 2011 md:141-147 — the composite-of-independent-runs recipe.
  - `src/WindowedComposite.jl` — the consumer.
"""
module WindowedTiling

export _window_seed, _nearest_center, _tile_centers, _warn_single_window,
       _own_poles

# Deterministic, version-stable per-window seed.  Two distinct global
# seeds re-randomise every window (the non-gameability property the
# seed-invariance gate relies on — see the module docstring "Why
# per-window seeding is load-bearing").  We combine arithmetically and
# NOT with `hash`: `hash` is not stable across Julia versions, which
# would silently break the package's bit-reproducibility contract.  The
# `1_000_003` multiplier (a prime well above any realistic window count)
# keeps `(rng_seed, wi)` pairs collision-free for sane inputs.
_window_seed(rng_seed::Integer, wi::Integer) = Int(rng_seed) * 1_000_003 + wi

# Index of the centre nearest to `z`, with a lowest-index tie-break.
# We compare *squared* distances (`abs2`): monotone in distance, so the
# argmin is identical, and it avoids a `sqrt` per candidate.  The strict
# `<` is what makes the tie-break deterministic — on the exact core
# boundary the lower window index wins, matching the Voronoi partition.
function _nearest_center(z::Complex, centers::AbstractVector{<:Complex})
    best = 1
    @inbounds bestd = abs2(z - centers[1])
    @inbounds for k in 2:length(centers)
        d = abs2(z - centers[k])
        if d < bestd
            bestd = d
            best  = k
        end
    end
    return best
end

# Voronoi ownership of per-window pole estimates, decided once per physical
# pole (see the module docstring "Pole ownership").  `ps[wi]` is window
# `wi`'s extracted pole list.  Returns the kept poles in window-major
# order — identical to the per-estimate rule whenever no pole straddles.
function _own_poles(ps::AbstractVector{<:AbstractVector{Complex{T}}},
                    centers::AbstractVector{<:Complex}, atol::Real) where {T}
    z   = reduce(vcat, ps; init = Complex{T}[])
    win = reduce(vcat, [fill(wi, length(ps[wi])) for wi in eachindex(ps)]; init = Int[])
    n   = length(z)
    parent = collect(1:n)
    root(a) = (while parent[a] != a; parent[a] = parent[parent[a]]; a = parent[a]; end; a)
    # Link cross-window estimates within `atol`: sweep a real-part-sorted
    # order so only candidates inside the `atol` strip are compared.
    ord = sortperm(real.(z))
    # (Two nested loops, NOT `for a in …, b in …`: a `break` in Julia's
    # comma form would exit both loops and silently skip every later `a`.)
    for (ii, a) in enumerate(ord)
        for jj in ii+1:n
            b = ord[jj]
            real(z[b]) - real(z[a]) > atol && break
            (win[a] != win[b] && abs(z[a] - z[b]) ≤ atol) || continue
            ra, rb = root(a), root(b)
            ra == rb || (parent[rb] = ra)
        end
    end
    # One canonical position (the group mean) ⇒ one owner per group.
    gsum = Dict{Int,Complex{T}}(); gcnt = Dict{Int,Int}()
    for a in 1:n
        r = root(a)
        gsum[r] = get(gsum, r, zero(Complex{T})) + z[a]
        gcnt[r] = get(gcnt, r, 0) + 1
    end
    owner = Dict(r => _nearest_center(gsum[r] / gcnt[r], centers) for r in keys(gsum))
    return Complex{T}[z[a] for a in 1:n if owner[root(a)] == win[a]]
end

# Core-centre coordinates tiling `[lo, hi]` at `extent` spacing.  We use
# `ceil(span/extent)` cores (≥1) and centre the run on the domain
# midpoint, so the union of the cores `[c ± extent/2]` covers `[lo, hi]`
# symmetrically.  Equal-size, equally-spaced cores ⇒ their Voronoi
# partition is exactly this tiling, so every domain point lands in one.
function _tile_centers(lo::T, hi::T, extent::T) where {T}
    span  = hi - lo
    n     = max(1, ceil(Int, span / extent))
    start = (lo + hi) / 2 - n * extent / 2
    return T[start + (k - T(0.5)) * extent for k in 1:n]
end

# Bead padetaylor-fqwf.  A 1×1 tile grid on an axis of positive span IS
# the monolithic path-dependent solve this module replaces (ADR-0034).
# Warn, not throw: one window is legal (tiny domain), just seam-uncured.
function _warn_single_window(axis::Symbol, n::Int, span::Real, extent::Real)
    n == 1 && span > 0 || return
    @warn "windowed_path_network_solve: window_extent=$extent ≥ $axis-span=" *
          "$span gives ONE window on $axis, so the composite degenerates to " *
          "the monolithic path-dependent solve (bead padetaylor-fqwf). " *
          "Suggestion: reduce window_extent so several windows tile each axis."
end

end # module WindowedTiling
