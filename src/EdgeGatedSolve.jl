"""
    PadeTaylor.EdgeGatedSolve

Edge-gated pole-field solve — confine the IVP path-network to the
pole field, so it never wanders into (and is corrupted by) the
smooth regions.

## Why this module exists

`PathNetwork.path_network_solve` is an *initial-value* solver. FW 2011
is blunt about its limit (md:401):

> "smooth regions are unstable regions for any IVP solver, and the
> associated loss of accuracy ought not be carried back into the pole
> field. One way to prevent this is to **force the path selection
> algorithm to complete pole fields before stepping into smooth
> regions** … with the pole field edge detection procedure described
> in Section 3.2.2, this process can readily be **fully automated**."

The failure is visible and severe. Integrate the PI *tritronquée*
(pole-free in four of its five sectors) across a large grid with a
plain `path_network_solve`, and accumulated IVP error perturbs the
solution off the tritronquée manifold — spurious poles bloom in the
sectors that should be empty. Over `[-50,50]²` the extracted pole
field is nearly angle-uniform; the four pole-free sectors fill in.

The cure is exactly FW's: **never target a smooth cell.** If the IVP's
target grid contains only pole-field cells, the wedge walk stays in
the low-`|u|` valleys between poles and never traverses the unstable
smooth terrain.

## The algorithm — region growing with a robust gate

We cannot know the pole-field extent before solving, and we must not
solve the smooth cells to find out. The resolution is to grow the
targeted region outward from the initial condition, one frontier at a
time, admitting only cells the edge detector confirms as pole-field —
*and* only cells that stay connected to the field already found:

  1. **Seed.** `field` ← grid cells within `seed_radius` of the IC.
     For a Painlevé IVP the pole field originates at the IC; that is
     where the solver is designed to start (FW §3.1).
  2. **Grow.** Repeat to a fixpoint:
       a. `targets` ← `field` dilated by `grow_rings` cells.
       b. `path_network_solve(prob, targets)` — the IVP, confined.
       c. `pole_field_mask` (FW eq. 3.3) on the solved values.
       d. **Clean the mask.** The FW edge level is calibrated for fine
          grids; on a coarse lattice the discretisation residual
          `Δu ∼ h²·u⁗/12` is large, so the raw mask carries scattered
          false positives. Morphological *opening* (erode then dilate)
          removes specks and one-cell-wide bridges; the genuine pole
          field is a thick connected blob and survives.
       e. **Flood-fill.** Keep only mask cells *connected* to the
          existing `field`. A false positive stranded in a smooth
          sector — or a thin bridge the opening missed — is never
          reached, so the field cannot leak across a smooth gap.
       f. `field` ← `field` ∪ {newly reached pole-field cells}.
       g. Stop when no new cell is admitted.
  3. **Final solve.** One last `path_network_solve` over `field`
     dilated by a single ring — the thinnest smooth frontier that
     still lets every field cell be reached. The returned
     `PathNetworkSolution`'s visited Padé store is then confined to
     the pole field (plus one harmless ring), so
     `PoleField.extract_poles` reads back a faithful, sector-confined
     pole field.

The opening + flood-fill pair is what makes the per-cell edge
classifier — noisy on its own at practical grid resolutions — into a
reliable region gate. It is the automatic stand-in for FW's "we
implemented this manually after having inspected … where the pole
field edges appeared" (md:401).

## Grid resolution matters

The §3.2.2 edge detector separates pole-field from smooth only when
the lattice is fine enough that the smooth-region discretisation
residual stays well below the pole-field signal — FW used a 161×161
grid over `[-10,10]²` (spacing `≈0.125`). At a spacing of `≈2` the two
distributions overlap and even opening + flood-fill cannot recover a
clean gate. Use a lattice spacing of roughly `≲ 1` for this driver;
`edge_gated_pole_field_solve` does not silently accept a doomed grid —
it is the caller's responsibility per the FW reference.

## What this module does *not* do — and why

It does **not** fill the smooth regions with the BVP solver. For
*bounded* smooth bands flanked by pole fields, `LatticeDispatcher.
lattice_dispatch_solve` already does the per-row Chebyshev–Newton fill
(FW md:190). For the motivating case — the tritronquée — the smooth
sectors are *unbounded*, which is not a horizontally-BVP-bridgeable
geometry (FW handles that family with the imaginary-axis BVP of
Fig 4.1, bead `padetaylor-gky`). And for a *pole-location* plot the
smooth sectors carry no poles — leaving them empty is correct.

## References

  - `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:401` — smooth regions are
    unstable for any IVP solver; gate the path selection with the
    §3.2.2 edge detector. `:202-208` — the edge detector itself.
  - `docs/adr/0004-path-network-architecture.md` — the path-network
    this module gates.
  - `src/EdgeDetector.jl`, `src/PathNetwork.jl` — the composed pieces.
"""
module EdgeGatedSolve

using ..Problems:     PadeTaylorProblem
using ..PathNetwork:  path_network_solve, PathNetworkSolution
using ..EdgeDetector: pole_field_mask

export edge_gated_pole_field_solve, EdgeGatedSolution

"""
    EdgeGatedSolution{T}

Output of `edge_gated_pole_field_solve`. Fields:

  - `xs`, `ys` — the lattice axes (`grid[i,j] = xs[i] + im·ys[j]`).
  - `field_mask::BitMatrix` — `true` ⟺ cell `(i,j)` was confirmed
    pole-field by the edge detector during region growing.
  - `pn_solution::PathNetworkSolution{T}` — the final IVP solve,
    confined to `field_mask` dilated by one ring. Pass *this* to
    `PoleField.extract_poles` for the sector-confined pole scatter.
  - `u_grid::Matrix{Complex{T}}` — the final solve's `u` values
    scattered onto the `(nx, ny)` lattice (`u_grid[i,j]` at
    `xs[i] + im·ys[j]`); cells outside `field_mask` dilated by one
    ring are `NaN + NaN·im`. The convenient 2D view of `pn_solution`.
  - `iterations::Int` — region-growing passes run before the fixpoint.
"""
struct EdgeGatedSolution{T <: AbstractFloat}
    xs          :: AbstractVector{T}
    ys          :: AbstractVector{T}
    field_mask  :: BitMatrix
    pn_solution :: PathNetworkSolution{T}
    u_grid      :: Matrix{Complex{T}}
    iterations  :: Int
end

# Chebyshev-distance (8-connected) dilation of `mask` by `r` cells.
function _dilate(mask::BitMatrix, r::Integer)
    nx, ny = size(mask)
    out = copy(mask)
    for _ in 1:r
        prev = copy(out)
        @inbounds for j in 1:ny, i in 1:nx
            prev[i, j] && continue
            for dj in -1:1, di in -1:1
                (di == 0 && dj == 0) && continue
                ii = i + di; jj = j + dj
                if 1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny && prev[ii, jj]
                    out[i, j] = true
                    break
                end
            end
        end
    end
    return out
end

# Chebyshev-distance (8-connected) erosion of `mask` by `r` cells. A
# cell survives only if every cell within Chebyshev distance `r` is
# `true`; grid-edge cells (a window hanging off the lattice) erode away.
function _erode(mask::BitMatrix, r::Integer)
    nx, ny = size(mask)
    out = copy(mask)
    for _ in 1:r
        prev = copy(out)
        @inbounds for j in 1:ny, i in 1:nx
            prev[i, j] || continue
            for dj in -1:1, di in -1:1
                (di == 0 && dj == 0) && continue
                ii = i + di; jj = j + dj
                if !(1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny) || !prev[ii, jj]
                    out[i, j] = false
                    break
                end
            end
        end
    end
    return out
end

# Morphological opening: erode then dilate by `r`. Removes features
# thinner than `2r` cells (scattered false positives, one-cell-wide
# bridges) while leaving thick blobs at their original size.
_open(mask::BitMatrix, r::Integer) = _dilate(_erode(mask, r), r)

# Cells of `mask` reachable from any `seed` cell by an 8-connected
# walk that stays inside `mask`. The connected component(s) of `mask`
# that contain the seed — a stranded false positive is excluded.
function _flood_fill(seed::BitMatrix, mask::BitMatrix)
    nx, ny = size(mask)
    reached = falses(nx, ny)
    stack = Tuple{Int,Int}[]
    @inbounds for j in 1:ny, i in 1:nx
        if seed[i, j] && mask[i, j]
            reached[i, j] = true
            push!(stack, (i, j))
        end
    end
    @inbounds while !isempty(stack)
        i, j = pop!(stack)
        for dj in -1:1, di in -1:1
            (di == 0 && dj == 0) && continue
            ii = i + di; jj = j + dj
            (1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny) || continue
            (mask[ii, jj] && !reached[ii, jj]) || continue
            reached[ii, jj] = true
            push!(stack, (ii, jj))
        end
    end
    return reached
end

# Solve the IVP over exactly the `true` cells of `targets`, and scatter
# the result back onto an `(nx, ny)` `u_grid` (`NaN` where unsolved).
function _solve_targets(prob, xs, ys, targets::BitMatrix, h, order)
    nx, ny = size(targets)
    CT = Complex{float(eltype(xs))}
    ij = Tuple{Int,Int}[]
    zs = CT[]
    for j in 1:ny, i in 1:nx
        targets[i, j] || continue
        push!(ij, (i, j))
        push!(zs, CT(xs[i] + im * ys[j]))
    end
    pn = path_network_solve(prob, zs; h = h, order = order)
    u_grid = fill(CT(NaN, NaN), nx, ny)
    @inbounds for (k, (i, j)) in enumerate(ij)
        u_grid[i, j] = pn.grid_u[k]
    end
    return pn, u_grid
end

"""
    edge_gated_pole_field_solve(prob::PadeTaylorProblem,
                                xs::AbstractVector{<:Real},
                                ys::AbstractVector{<:Real};
                                h = 0.5, order = prob.order,
                                edge_level = :auto,
                                seed_radius = nothing,
                                grow_rings = 3,
                                open_radius = 1,
                                max_iter = 500,
                                verbose = false) -> EdgeGatedSolution

Solve the 2nd-order IVP `prob` on the lattice `xs × ys`, confining the
path-network to the pole field by region growing (see the module
docstring for the algorithm and the FW 2011 md:401 rationale).

Kwargs:
  - `h`, `order` — passed through to `path_network_solve`.
  - `edge_level` — `pole_field_mask` threshold on `log₁₀|Δu|`.  Default
    `:auto` resolves to the h-aware `LEVEL0 + 2·log₁₀(min(h, H0)/H0)`
    (anchor `(0.25, 0.001)`, clamped at `H0`) — see `EdgeDetector`'s
    module docstring.  Pass `0.001` to reproduce FW Fig 3.3 verbatim.
  - `seed_radius` — radius around the IC (`prob.zspan[1]`) of the
    initial `field` seed. `nothing` ⇒ `max(3, 2.5·Δgrid)`, enough to
    reach the first poles of a Painlevé IVP without seeding a wide
    smooth patch.
  - `grow_rings` — cells the `field` is dilated by per pass; larger ⇒
    fewer (but heavier) `path_network_solve` calls. `grow_rings ≥ 2`.
  - `open_radius` — morphological-opening radius applied to each
    pass's edge mask before the flood-fill (`0` disables it). Removes
    false-positive specks and bridges thinner than `2·open_radius`.
  - `max_iter` — safety cap on growth passes.

Requires `xs`, `ys` strictly increasing with a common, uniform step
(the FW eq. 3.3 stencil is isotropic). Use a spacing `≲ 1` — see the
module docstring "Grid resolution matters". Throws `ArgumentError` on
a malformed lattice or `grow_rings < 2`.
"""
function edge_gated_pole_field_solve(prob::PadeTaylorProblem,
                                     xs::AbstractVector{<:Real},
                                     ys::AbstractVector{<:Real};
                                     h::Real          = 0.5,
                                     order::Integer   = prob.order,
                                     edge_level::Union{Real,Symbol} = :auto,
                                     seed_radius      = nothing,
                                     grow_rings::Integer  = 3,
                                     open_radius::Integer = 1,
                                     max_iter::Integer    = 500,
                                     verbose::Bool        = false)
    nx = length(xs); ny = length(ys)
    nx ≥ 3 && ny ≥ 3 || throw(ArgumentError(
        "edge_gated_pole_field_solve: need length(xs), length(ys) ≥ 3 " *
        "(got $nx, $ny); detail: the FW eq. 3.3 stencil needs an interior."))
    issorted(xs; lt = ≤) && issorted(ys; lt = ≤) || throw(ArgumentError(
        "edge_gated_pole_field_solve: xs and ys must be strictly increasing."))
    grow_rings ≥ 2 || throw(ArgumentError(
        "edge_gated_pole_field_solve: grow_rings must be ≥ 2 (got " *
        "$grow_rings); detail: a pass classifies field+grow_rings−1, " *
        "so grow_rings = 1 can never confirm a new cell."))
    open_radius ≥ 0 || throw(ArgumentError(
        "edge_gated_pole_field_solve: open_radius must be ≥ 0 (got $open_radius)."))

    T = float(eltype(xs))
    Δx = T(xs[2] - xs[1]); Δy = T(ys[2] - ys[1])
    isapprox(Δx, Δy; rtol = 1e-10) || throw(ArgumentError(
        "edge_gated_pole_field_solve: xs and ys must share one uniform " *
        "step (got Δx = $Δx, Δy = $Δy); detail: the 5-point stencil is " *
        "isotropic."))
    h_grid = Δx

    seed_r = seed_radius === nothing ? max(T(3), T(5) * h_grid / 2) : T(seed_radius)
    z_ic   = Complex{T}(prob.zspan[1])

    # --- Seed: the cells within `seed_r` of the IC ------------------------
    field = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        if abs(Complex{T}(xs[i] + im * ys[j]) - z_ic) ≤ seed_r
            field[i, j] = true
        end
    end
    any(field) || throw(ArgumentError(
        "edge_gated_pole_field_solve: seed_radius = $seed_r captured no " *
        "grid cell near the IC $z_ic; detail: increase seed_radius or " *
        "place the lattice over the IC."))

    # --- Grow: admit only edge-confirmed cells connected to the field -----
    iters = 0
    local pn
    while iters < max_iter
        iters += 1
        targets = _dilate(field, grow_rings)
        pn, u_grid = _solve_targets(prob, xs, ys, targets, h, order)

        # Raw FW eq. 3.3 classification, then clean it: the field is
        # pole-field by definition, opening kills false-positive specks
        # and thin bridges, the field is restored (opening can nibble a
        # blob's edge), and the flood-fill keeps only what is connected
        # to the field already found.
        mask = pole_field_mask(u_grid, h_grid; level = edge_level)
        mask .|= field
        open_radius > 0 && (mask = _open(mask, open_radius))
        mask .|= field
        reachable = _flood_fill(field, mask)

        new_cells = reachable .& targets .& .!field
        n_new = count(new_cells)
        verbose && @info "edge_gated_pole_field_solve" pass=iters new=n_new field=count(field)
        n_new == 0 && break
        field .|= new_cells
    end

    # --- Final solve over field + one thin frontier ring ------------------
    final_targets = _dilate(field, 1)
    pn, u_grid = _solve_targets(prob, xs, ys, final_targets, h, order)

    return EdgeGatedSolution{T}(xs, ys, field, pn, u_grid, iters)
end

end # module EdgeGatedSolve
