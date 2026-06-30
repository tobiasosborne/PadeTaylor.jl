"""
    PadeTaylor.WindowedComposite

Bounded-window composite solve — tile the plane, solve each window
*independently* from the initial condition, and composite by
Voronoi-core ownership.  The confirmed cure for the FW 2011 Fig 4.7
"seam".

## Why this module exists — the path-dependence seam

`PathNetwork.path_network_solve` walks one monolithic tree from the IC
across the *whole* grid.  Two facts about that tree combine into a
bug.  (1) The Stage-1 target order is randomised (`rng_seed` shuffles
the targets), so the tree's branch structure — which node is whose
parent, which trunk a far cell hangs off — depends on the seed.
(2) A long trunk accumulates IVP error, and that error *differs*
between two seeds because the trunks differ.  Where two seed-different
subtrees meet, their independently-accumulated errors disagree: a
**grain boundary** in the reconstructed field.  Worklog 077
root-caused it; worklog 078 measured it.

Critically, worklog 078 §2 showed the seam is **regime B** — a pole-
field *interior* grain boundary, not a smooth-sector artefact.
Edge-gating (`EdgeGatedSolve`, already shipped, already in the
figure) confines the walk to the pole field and leaves the seam
essentially untouched (pole-count seed-variance |Δ| 172 → 171).  The
cure has to act *inside* the pole field, on the path-dependence
itself.

## The cure — FW's own remedy (md:147)

FW 2011 never solved Fig 4.7 monolithically.  They tiled it:

> "each pole field shown in Fig. 4.7a–e is a 5×5 composite of 25
> completely independent runs (each starting from the same ICs at
> z = 0, and covering areas of the extent 20×20 in the complex
> plane).  Similarly, Fig. 4.8 is a composite of 18 independent runs."
> — `references/markdown/FW2011_painleve_methodology_JCP230/
>    FW2011_painleve_methodology_JCP230.md:147`

The mechanism is exactly why it works.  Every cell in one window
hangs off the *same* short trunk back to the shared IC, so the cells'
*differential* positions are seed-stable even though the trunk itself
is not (worklog 078 §3a: a far window collapses to |Δ|=1).  The
cross-branch differential that produced the seam is eliminated *by
construction* — there is no second subtree to disagree with.

## The algorithm

  1. **Tile.** Lay down core centres on a regular `window_extent`-
     spaced grid that covers `[xs[1],xs[end]] × [ys[1],ys[end]]`,
     centred so the cores tile the domain (`_tile_centers`).  The
     cores are the Voronoi sites; for equal-size, equally-spaced
     cores the Voronoi partition *is* the core tiling, so every grid
     cell lies in exactly one core.
  2. **Solve independently.** Each window solves a region that is its
     core expanded by `overlap` on every side (clamped to the
     domain).  The overlap gives interior nodes that resolve poles
     sitting on the core's own edge; it is solved, then discarded.
     Each window calls `path_network_solve` from the shared IC at its
     *own* seed (see below).
  3. **Composite by Voronoi-core.** A window writes a solved cell into
     the global field only if that cell's nearest core centre is the
     window's own (`_nearest_center == wi`).  Because the core ⊂ the
     solve extent, the owning window always solved the cell, so every
     cell is written exactly once — a hard partition, no double
     counting, no boundary seam.

## Why per-window seeding is load-bearing (non-gameable gate)

Each window is seeded `_window_seed(rng_seed, wi)` — a deterministic,
version-stable combine that is genuinely distinct per `(global seed,
window)`.  This is *not* a convenience: the seed-invariance gate
(`padetaylor-sny7`) re-runs the whole composite at two global seeds
and asserts the pole field barely moves.  If a window ignored the
global seed (e.g. always seed 0), the gate would pass trivially — the
composite would be seed-independent because nothing re-randomised.
Threading the global seed into *every* window means two global seeds
re-randomise *every* window's internal walk, so a passing gate
certifies genuine path-independence, not a frozen pipeline.  We
combine arithmetically (`rng_seed·1_000_003 + wi`), **never** `hash`,
because `hash` is not stable across Julia versions and would break the
package's bit-reproducibility contract.

## Confirmed numbers (worklog 078 §4, p5 acceptance probe)

A 5×5 composite over `[-50,50]²` against the monolithic baseline, two
global seeds, poles in `[-49,49]²`:

  - pole-count seed-variance |Δ| **151 → 12** (~12×);
  - pole-set agreement @0.5 **77.1 % → 99.4 %**;
  - median pole displacement **0.076 → 0.000**;
  - tile-boundary band match 98.9 % vs interior 99.2 % — **no
    boundary seam**; the overlap-core scheme handled it;
  - wall time **~214 s → ~102 s** (~2× faster — the windows are
    independent and each walk is short).

## Honest residuals carried forward (Rule 9)

  - **|Δ|≈12** is *containment*, not annihilation — a handful of poles
    near the noise floor still wobble.  This matches FW's own posture;
    the `sny7` gate tolerance sits at this measured floor, not a
    relaxed one.
  - The composite finds **~220 fewer poles** than the monolithic walk
    (7223 vs 7443).  These are almost certainly the spurious *seed-
    dependent* poles correctly suppressed — the count gap *is* the bug
    — but that is a *named verification* (`padetaylor-ingn`, via the
    Weierstrass-℘ truth anchor), not an assumption baked in here.

## References

  - `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:147` — the 5×5 / 18-run
    composite construction this module productionises.
  - `docs/worklog/077-*` — root cause of the seam (path dependence).
  - `docs/worklog/078-fig47-seam-cure-scoping.md` — mechanism
    selection and the p5 confirmation numbers above.
  - `src/PathNetwork.jl`, `src/PoleField.jl` — the per-window solve
    and pole extraction this module composites.
"""
module WindowedComposite

using ..Problems:    PadeTaylorProblem
using ..PathNetwork: path_network_solve, PathNetworkSolution
using ..PoleField:   extract_poles

export windowed_path_network_solve, windowed_extract_poles,
       WindowedCompositeSolution

"""
    WindowedCompositeSolution{T}

Output of `windowed_path_network_solve`.  The composited field plus
every per-window solve, so a caller can re-extract poles, inspect a
single window, or audit the Voronoi assignment.  Fields:

  - `xs`, `ys` — the lattice axes (`grid_u[i,j]` is at `xs[i] +
    im·ys[j]`).
  - `grid_u::Matrix{Complex{T}}` — the composited `u` field; cell
    `(i,j)` carries the value from the window that owns it.
  - `window_assign::Matrix{Int}` — the owning window index per cell
    (`0` means unassigned — a coverage bug that the solve refuses to
    return, Rule 1).
  - `windows::Vector{NTuple{4,T}}` — each window's *solved* extent
    `(xlo, xhi, ylo, yhi)` (core expanded by `overlap`, clamped).
  - `centers::Vector{Complex{T}}` — the window **core** centres, i.e.
    the Voronoi sites that decide cell and pole ownership.
  - `window_seeds::Vector{Int}` — the per-window seed actually passed
    to `path_network_solve` (`_window_seed(rng_seed, wi)`).
  - `window_sols::Vector{PathNetworkSolution{T}}` — the independent
    per-window solves, in window order.  Pass these (filtered to each
    window's Voronoi core) through `windowed_extract_poles`.
"""
struct WindowedCompositeSolution{T <: AbstractFloat}
    xs            :: AbstractVector{T}
    ys            :: AbstractVector{T}
    grid_u        :: Matrix{Complex{T}}
    window_assign :: Matrix{Int}
    windows       :: Vector{NTuple{4,T}}
    centers       :: Vector{Complex{T}}
    window_seeds  :: Vector{Int}
    window_sols   :: Vector{PathNetworkSolution{T}}
end

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

"""
    windowed_path_network_solve(prob::PadeTaylorProblem,
                                xs::AbstractVector{<:Real},
                                ys::AbstractVector{<:Real};
                                window_extent = 20.0, overlap = 6.0,
                                h = 0.5, order = prob.order,
                                rng_seed = 0, extrapolate = true,
                                verbose = false)
        -> WindowedCompositeSolution

Solve `prob` on the lattice `xs × ys` as a **bounded-window
composite** — FW 2011's own Fig 4.7 construction (md:147).  See the
module docstring for the seam this cures and the mechanism.

Tile the domain into cores of side `window_extent`; solve each core
expanded by `overlap` independently from the IC at a per-window seed
derived from `rng_seed`; composite by Voronoi-core ownership.  Defaults
reproduce FW's recipe: `window_extent = 20.0` ("extent 20×20"),
`overlap = 6.0` (interior nodes to resolve core-edge poles).

Kwargs:
  - `window_extent` — core side length (the Voronoi-site spacing).
  - `overlap` — margin added to every side of a core to form its solve
    extent; discarded after compositing (`overlap = 0` ⇒ bare cores).
  - `h`, `order`, `extrapolate` — passed through to each window's
    `path_network_solve`.  `extrapolate = true` (default) fills the
    panel as the figure scripts need.
  - `rng_seed` — the *global* seed; window `wi` is solved at
    `_window_seed(rng_seed, wi)`, so two global seeds re-randomise
    every window (non-gameable; see the module docstring).
  - `verbose` — per-window progress to `stdout`.

Throws `ArgumentError` on a malformed request (fewer than 3 nodes per
axis, non-increasing axes, `window_extent ≤ 0`, `overlap < 0`, or a
window that captures no grid cell because the grid is too coarse for
`window_extent`), and refuses to return a field with any unassigned
cell (Rule 1 — a coverage bug, never a silent NaN).
"""
function windowed_path_network_solve(prob::PadeTaylorProblem,
                                     xs::AbstractVector{<:Real},
                                     ys::AbstractVector{<:Real};
                                     window_extent::Real = 20.0,
                                     overlap::Real       = 6.0,
                                     h::Real             = 0.5,
                                     order::Integer      = prob.order,
                                     rng_seed::Integer   = 0,
                                     extrapolate::Bool   = true,
                                     verbose::Bool       = false)
    nx = length(xs); ny = length(ys)
    nx ≥ 3 && ny ≥ 3 || throw(ArgumentError(
        "windowed_path_network_solve: need length(xs), length(ys) ≥ 3 " *
        "(got $nx, $ny). Suggestion: a window composite needs an interior " *
        "lattice; pass at least a 3×3 grid."))
    issorted(xs; lt = ≤) && issorted(ys; lt = ≤) || throw(ArgumentError(
        "windowed_path_network_solve: xs and ys must be strictly " *
        "increasing. Suggestion: pass `range(lo, hi; length=N)` axes."))
    window_extent > 0 || throw(ArgumentError(
        "windowed_path_network_solve: window_extent must be > 0 (got " *
        "$window_extent). Suggestion: FW 2011 md:147 uses 20.0."))
    overlap ≥ 0 || throw(ArgumentError(
        "windowed_path_network_solve: overlap must be ≥ 0 (got $overlap). " *
        "Suggestion: FW 2011's composites use an overlap margin of ~6.0."))

    T   = float(typeof(real(first(xs))))
    CT  = Complex{T}
    ext = T(window_extent); ov = T(overlap); half = ext / 2
    xlo_d = T(xs[1]); xhi_d = T(xs[end]); ylo_d = T(ys[1]); yhi_d = T(ys[end])

    # Tile the domain into core centres = the Voronoi sites.
    xc = _tile_centers(xlo_d, xhi_d, ext)
    yc = _tile_centers(ylo_d, yhi_d, ext)
    centers = CT[complex(cx, cy) for cy in yc for cx in xc]
    nwin = length(centers)

    grid_u        = fill(CT(NaN, NaN), nx, ny)
    window_assign = zeros(Int, nx, ny)
    windows       = Vector{NTuple{4,T}}(undef, nwin)
    window_seeds  = Vector{Int}(undef, nwin)
    window_sols   = Vector{PathNetworkSolution{T}}(undef, nwin)

    for wi in 1:nwin
        c  = centers[wi]
        cx = real(c); cy = imag(c)
        # Solve extent = core ± overlap, clamped to the domain.
        xlo = max(xlo_d, cx - half - ov); xhi = min(xhi_d, cx + half + ov)
        ylo = max(ylo_d, cy - half - ov); yhi = min(yhi_d, cy + half + ov)
        windows[wi] = (xlo, xhi, ylo, yhi)

        # Grid cells (and their complex coordinates) inside the extent.
        ij = Tuple{Int,Int}[]; zs = CT[]
        @inbounds for j in 1:ny, i in 1:nx
            x = T(xs[i]); y = T(ys[j])
            (xlo ≤ x ≤ xhi && ylo ≤ y ≤ yhi) || continue
            push!(ij, (i, j)); push!(zs, complex(x, y))
        end
        isempty(zs) && throw(ArgumentError(
            "windowed_path_network_solve: window $wi/$nwin (core centre $c) " *
            "captured no grid cell. Suggestion: window_extent=$window_extent " *
            "is too small for this grid spacing — increase window_extent or " *
            "refine the lattice."))

        seed = _window_seed(rng_seed, wi)
        window_seeds[wi] = seed
        verbose && @info "windowed_path_network_solve" window=wi total=nwin center=c cells=length(zs)

        # Independent per-window solve from the shared IC.  Rethrow with
        # window context (Rule 1) so a per-window breakdown is locatable.
        local pn
        try
            pn = path_network_solve(prob, zs; h = h, order = order,
                                    rng_seed = seed, extrapolate = extrapolate,
                                    verbose = verbose)
        catch err
            throw(ErrorException(
                "windowed_path_network_solve: window $wi/$nwin (core centre " *
                "$c, solve extent ($xlo,$xhi)×($ylo,$yhi), seed $seed) failed: " *
                "$(sprint(showerror, err)). Suggestion: inspect this window in " *
                "isolation with path_network_solve and the same kwargs."))
        end
        window_sols[wi] = pn

        # Composite by Voronoi core: keep only cells this window owns.
        @inbounds for (k, (i, j)) in enumerate(ij)
            _nearest_center(zs[k], centers) == wi || continue
            grid_u[i, j]        = pn.grid_u[k]
            window_assign[i, j] = wi
        end
    end

    # Every domain cell lies in exactly one core, whose owning window
    # always solves it — so a remaining 0 is a coverage bug, not data.
    if any(==(0), window_assign)
        n_unassigned = count(==(0), window_assign)
        throw(ErrorException(
            "windowed_path_network_solve: $n_unassigned grid cell(s) were " *
            "never assigned to a window — a Voronoi coverage bug. Suggestion: " *
            "this should be impossible (core ⊂ solve extent); file a bead with " *
            "the (window_extent=$window_extent, overlap=$overlap) that hit it."))
    end

    return WindowedCompositeSolution{T}(xs, ys, grid_u, window_assign,
                                        windows, centers, window_seeds,
                                        window_sols)
end

"""
    windowed_extract_poles(wsol::WindowedCompositeSolution;
                           merge_atol = nothing, extract_kwargs...)
        -> Vector{Complex{T}}

Composite pole field of a `WindowedCompositeSolution`: extract poles
from each window's independent solve, keep only those in the window's
own Voronoi core, and union across windows.  This is the pole-side
counterpart of the field composite — the same hard core partition that
prevents a tile-boundary seam prevents double-counting a pole that two
overlapping windows both resolve.

`extract_kwargs...` pass straight through to `PoleField.extract_poles`
(`radius_t`, `min_residue`, `cluster_atol`, `min_support`, …).

`merge_atol` is an *optional* boundary-dedup refinement: when not
`nothing`, a greedy keep-first pass drops any pole within `merge_atol`
of an already-kept pole.  The default `nothing` is faithful to the p5
confirmation, which reached 99.4 % seed-agreement with **no** final
merge — the Voronoi-core partition already prevents almost all
double-counting, so the merge is only for the rare pole that straddles
a core line at the cluster scale.
"""
function windowed_extract_poles(wsol::WindowedCompositeSolution{T};
                                merge_atol::Union{Nothing,Real} = nothing,
                                extract_kwargs...) where {T}
    centers = wsol.centers
    kept = Complex{T}[]
    for wi in eachindex(wsol.window_sols)
        ps = extract_poles(wsol.window_sols[wi]; extract_kwargs...)
        for p in ps
            # Voronoi-core ownership: this window keeps only its own poles.
            _nearest_center(p, centers) == wi || continue
            # Optional light boundary dedup (greedy, keep first).
            if merge_atol !== nothing &&
               any(q -> abs(p - q) ≤ merge_atol, kept)
                continue
            end
            push!(kept, p)
        end
    end
    return kept
end

end # module WindowedComposite
