"""
    PadeTaylor.LatticeDispatcher

Tier-3+ 2D-lattice composition layer (Phase 12 v2, bead
`padetaylor-k31`) — stitches the IVP path-network (Phase 10), the
5-point Laplacian edge detector (Phase 12.5), and the Chebyshev-Newton
BVP solver (Phase 11) into a single 2D-grid solve that automatically
fills smooth bands per FW 2011 §4.4 (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:190`):

> "It remains to enforce the BCs … The total computer time to fill in a
> smooth band (161 separate BVP solutions; one for each grid line) is
> generally less than a second."

i.e., FW resolved 2D pole-field-plus-smooth-band figures by running ONE
BVP per horizontal grid line, with Dirichlet BCs taken from the IVP
solver's values at the two pole-field-edge cells flanking each smooth
run.  This module ships exactly that v1 algorithm.

## What's *not* in v1

v1 does **not** ship the FW Fig 4.1 quantitative acceptance pin
(`u(20i)` to `1e-10` abs).  The composition machinery and per-row BVP
fill are in place; the FW 4.1-specific BVP geometry (BVP between
`+20i` and `-20i` plus two outward-running pole fields) is a
*different* compositional pattern than "horizontal BVP rows", and is
left for a follow-up bead.  v1 ships the *horizontal-row* algorithm
(FW2011...md:190 verbatim).

v1 also does not classify cells beyond the binary `{pole-field, smooth}`
distinction.  The 5-point Laplacian classifier is faithful to FW
2011 §3.2.2 but does not separately tag "smooth-extending-to-infinity"
vs "smooth-bounded-on-both-sides".  A smooth run touching the grid
boundary cannot be BVP-bridged (one missing BC); v1 retains its IVP
values and tags it `:ivp_only` for downstream filtering.

## The algorithm

Given a 2D Cartesian lattice `xs × ys` and a Painlevé-like 2nd-order
IVP problem `prob = PadeTaylorProblem(f_ivp, (u₀, u'₀), ...)`:

  1. **IVP fill**: call `path_network_solve(prob, vec(grid); h_path,
     order)` over the flat grid; get per-cell `(u, u')`.
  2. **Partition**: reshape to a `(nx, ny)` matrix; call
     `pole_field_mask(u_grid, h_grid; level=edge_level)` to get a
     BitMatrix where `mask[i, j] = true` ⟺ cell `(xs[i], ys[j])` is
     classified as pole-field (high `|Δu|`).
  3. **Per-row BVP fill**: for each interior row `j ∈ 2..ny-1` (rows
     where the EdgeDetector mask is defined):
       a. Walk along the row, identifying maximal contiguous smooth
          runs (cells with `mask=false`) flanked by pole-field cells
          (`mask=true`) on both sides.
       b. For each such bridgeable run `[i_lo, i_hi]`, take the IVP
          values at the flanking cells (`i_lo - 1` and `i_hi + 1`) as
          Dirichlet BCs, then `bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b,
          u_a, u_b; N = N_bvp)`.
       c. Replace the IVP values at the smooth cells with the BVP's
          barycentric interpolant evaluated at each grid point in the
          run; tag those cells `:bvp`.
       d. Cells in smooth runs touching a grid boundary are retained
          as `:ivp_only` (BVP requires two-sided BCs).
  4. **Return**: a `LatticeSolution{T}` with the stitched 2D `u_grid`,
     `up_grid`, the partition `mask`, the per-cell `region_tag`, the
     underlying `PathNetworkSolution`, and the vector of BVP solutions.

The signature takes both an IVP RHS (via `prob.f`, expected
`f(z, u, u')`) and a BVP RHS (`bvp_f(z, u)` and `bvp_∂f_∂u(z, u)`,
1st-order in `u` per the FW BVP's 2nd-derivative formulation).  For
autonomous-in-`u'` ODEs (the FW target class — PI, PII, PIV etc. all
satisfy `u'' = f(z, u)` with no `u'` dependence), pass `bvp_f(z, u) =
f_ivp(z, u, _)` ignoring the third arg.

## Public API

  - `LatticeSolution{T}` — composed 2D output with per-cell region
    tag (`:ivp`, `:bvp`, `:ivp_only`), the partition mask, and the
    underlying IVP + BVP sub-solutions for diagnostic inspection.

  - `lattice_dispatch_solve(prob, bvp_f, bvp_∂f_∂u, xs, ys; kwargs...)` —
    the public driver.  Returns `LatticeSolution{T}`.

  - The `mask` kwarg accepts a pre-computed `BitMatrix` for callers who
    want to bypass the automatic EdgeDetector classification (the FW
    2011 line 401 "manual classification" workflow).

## Fail-fast contract (CLAUDE.md Rule 1)

Throws `ArgumentError` on:
  - `length(xs) < 3` or `length(ys) < 3`.
  - `xs` or `ys` not strictly increasing.
  - `prob.y0` not a 2-tuple (lattice dispatcher requires 2nd-order IVP).
  - `mask` supplied with shape ≠ `(nx, ny)`.
  - `edge_level` invalid (Inf, NaN).
  - `N_bvp < 4`.

Throws `ErrorException` on:
  - Path-network failure (delegated from `path_network_solve`).
  - BVP non-convergence on any row (delegated from `bvp_solve`).

## References

  - Fornberg & Weideman, *A numerical methodology for the Painlevé
    equations*, J. Comput. Phys. 230 (2011) 5957–5973, §4.4 +
    line 190 (the "161 BVPs, one per grid line" line).
  - `references/markdown/FW2011_*.md:190, 218-222, 249-261`.
  - `docs/adr/0004-path-network-architecture.md`.
"""
module LatticeDispatcher

using ..Problems:    PadeTaylorProblem
using ..PathNetwork: path_network_solve, PathNetworkSolution
using ..BVP:         bvp_solve, BVPSolution
using ..EdgeDetector: pole_field_mask

export LatticeSolution, lattice_dispatch_solve

"""
    LatticeSolution{T}

Composed 2D-grid output from `lattice_dispatch_solve`.  Fields:

  - `xs::AbstractVector{T}`, `ys::AbstractVector{T}` — the lattice axes.
  - `grid_z::Matrix{Complex{T}}` — `grid_z[i, j] = xs[i] + im·ys[j]`.
  - `grid_u`, `grid_up::Matrix{Complex{T}}` — stitched solution and
    derivative.  Cells tagged `:bvp` carry the BVP barycentric
    interpolant; cells tagged `:ivp` or `:ivp_only` carry the
    PathNetwork output.
  - `mask::BitMatrix` — the EdgeDetector partition.  `true` ⟺
    classified pole-field.
  - `region_tag::Matrix{Symbol}` — per-cell tag in `{:ivp, :bvp,
    :ivp_only}`.  See the module docstring for definitions.
  - `pn_solution::PathNetworkSolution{T}` — the underlying IVP solve.
  - `bvp_solutions::Vector{BVPSolution{T, Complex{T}}}` — one BVP per
    bridged smooth run.
"""
struct LatticeSolution{T <: AbstractFloat}
    xs            :: AbstractVector{T}
    ys            :: AbstractVector{T}
    grid_z        :: Matrix{Complex{T}}
    grid_u        :: Matrix{Complex{T}}
    grid_up       :: Matrix{Complex{T}}
    mask          :: BitMatrix
    region_tag    :: Matrix{Symbol}
    pn_solution   :: PathNetworkSolution{T}
    bvp_solutions :: Vector{BVPSolution{T, Complex{T}}}
end

"""
    lattice_dispatch_solve(prob, bvp_f, bvp_∂f_∂u, xs, ys;
                           h_path = 0.5, order = prob.order,
                           edge_level = :auto,
                           N_bvp = 20, bvp_tol = nothing,
                           mask = nothing) -> LatticeSolution

See module docstring for the algorithm.  `order` defaults to
`prob.order` and is threaded into `path_network_solve` for the IVP
fill; pass it explicitly only to override the problem's own order.
"""
function lattice_dispatch_solve(prob::PadeTaylorProblem,
                                bvp_f, bvp_∂f_∂u,
                                xs::AbstractVector{<:Real},
                                ys::AbstractVector{<:Real};
                                h_path::Real      = 0.5,
                                order::Integer    = prob.order,
                                edge_level::Union{Real,Symbol} = :auto,
                                N_bvp::Integer    = 20,
                                bvp_tol           = nothing,
                                mask              = nothing)

    nx = length(xs); ny = length(ys)
    nx ≥ 3 || throw(ArgumentError(
        "lattice_dispatch_solve: length(xs) must be ≥ 3 (got $nx); " *
        "detail: EdgeDetector needs a 3-cell interior."))
    ny ≥ 3 || throw(ArgumentError(
        "lattice_dispatch_solve: length(ys) must be ≥ 3 (got $ny); " *
        "detail: EdgeDetector needs a 3-cell interior."))
    issorted(xs; lt = ≤) || throw(ArgumentError(
        "lattice_dispatch_solve: xs must be strictly increasing."))
    issorted(ys; lt = ≤) || throw(ArgumentError(
        "lattice_dispatch_solve: ys must be strictly increasing."))
    length(prob.y0) == 2 || throw(ArgumentError(
        "lattice_dispatch_solve: prob.y0 must be a 2-tuple " *
        "(u₀, u'₀) for 2nd-order IVP; got $(length(prob.y0)) entries."))
    (edge_level === :auto || (edge_level isa Real && isfinite(edge_level))) ||
        throw(ArgumentError(
            "lattice_dispatch_solve: edge_level must be a finite Real " *
            "or the `:auto` sentinel (got $(repr(edge_level)))."))
    N_bvp ≥ 4 || throw(ArgumentError(
        "lattice_dispatch_solve: N_bvp must be ≥ 4 (got $N_bvp); " *
        "detail: Chebyshev-Newton needs at least 5 nodes."))

    T  = float(eltype(xs))
    CT = Complex{T}

    # Equispaced lattice requirement: the EdgeDetector stencil assumes
    # uniform spacing.  We use `step(xs)` if available, else compute
    # from xs[2]-xs[1] and trust the caller (validated above).
    h_grid_x = T(xs[2] - xs[1])
    h_grid_y = T(ys[2] - ys[1])
    isapprox(h_grid_x, h_grid_y; rtol = 1e-10) || throw(ArgumentError(
        "lattice_dispatch_solve: xs and ys must have the same step " *
        "size (got Δx = $h_grid_x, Δy = $h_grid_y); detail: 5-point " *
        "stencil is isotropic."))

    # --- Step 1: 2D-grid IVP via path-network -------------------------------
    grid_mat = CT[xs[i] + im * ys[j] for i in 1:nx, j in 1:ny]
    pn_sol = path_network_solve(prob, vec(grid_mat); h = h_path, order = order)

    u_grid  = Matrix{CT}(reshape(pn_sol.grid_u,  (nx, ny)))
    up_grid = Matrix{CT}(reshape(pn_sol.grid_up, (nx, ny)))

    # --- Step 2: Partition via EdgeDetector (or use supplied mask) ----------
    if mask === nothing
        mask_used = pole_field_mask(u_grid, h_grid_x; level = edge_level)
    else
        size(mask) == (nx, ny) || throw(ArgumentError(
            "lattice_dispatch_solve: mask must have shape ($nx, $ny); " *
            "got $(size(mask))."))
        mask_used = convert(BitMatrix, mask)
    end

    # --- Step 3: Per-row BVP fill -------------------------------------------
    # Default tag: cells with mask=true are :ivp; mask=false cells are
    # :ivp_only until/unless they get bridged to :bvp.
    region_tag = Matrix{Symbol}(undef, nx, ny)
    for i in 1:nx, j in 1:ny
        region_tag[i, j] = mask_used[i, j] ? :ivp : :ivp_only
    end

    bvp_solutions = BVPSolution{T, CT}[]
    # Iterate interior rows only — the EdgeDetector mask is undefined
    # at j = 1 and j = ny (boundary rows have mask = false everywhere).
    for j in 2:(ny - 1)
        i = 2
        while i ≤ nx - 1
            if mask_used[i, j]
                i += 1
                continue
            end
            # Found the start of a smooth run.  Walk to its end.
            run_start = i
            while i ≤ nx - 1 && !mask_used[i, j]
                i += 1
            end
            run_end = i - 1
            # Both flanks must be interior IVP cells (mask=true).
            left_flank  = run_start - 1
            right_flank = run_end + 1
            if 2 ≤ left_flank ≤ nx - 1 && 2 ≤ right_flank ≤ nx - 1 &&
               mask_used[left_flank, j] && mask_used[right_flank, j]

                z_a = CT(xs[left_flank]  + im * ys[j])
                z_b = CT(xs[right_flank] + im * ys[j])
                u_a = u_grid[left_flank,  j]
                u_b = u_grid[right_flank, j]

                bvp_sol = bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b, u_a, u_b;
                                    N = N_bvp, tol = bvp_tol)
                push!(bvp_solutions, bvp_sol)

                for k in run_start:run_end
                    z_k = CT(xs[k] + im * ys[j])
                    u_k, up_k = bvp_sol(z_k)
                    u_grid[k, j]  = u_k
                    up_grid[k, j] = up_k
                    region_tag[k, j] = :bvp
                end
            end
        end
    end

    return LatticeSolution{T}(xs, ys, grid_mat, u_grid, up_grid,
                              mask_used, region_tag, pn_sol, bvp_solutions)
end

end # module LatticeDispatcher
