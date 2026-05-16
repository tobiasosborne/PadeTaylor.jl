"""
    PadeTaylor.LatticeDispatcher

Tier-3+ 2D-lattice composition layer (Phase 12 v2, bead
`padetaylor-k31`; v3 bead `padetaylor-0tj`, ADR-0017) — stitches the
edge-gated path-network (Phase 12.7), the 5-point Laplacian edge
detector (Phase 12.5), and the Chebyshev-Newton BVP solver (Phase 11)
into a single 2D-grid solve that automatically fills smooth bands per
FW 2011 §4.4 (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:190`):

> "It remains to enforce the BCs … The total computer time to fill in a
> smooth band (161 separate BVP solutions; one for each grid line) is
> generally less than a second."

i.e., FW resolved 2D pole-field-plus-smooth-band figures by running ONE
BVP per horizontal grid line, with Dirichlet BCs taken from the IVP
solver's values at the two pole-field-edge cells flanking each smooth
run.  This module ships exactly that v1 algorithm.

## Why the IVP source is `EdgeGatedSolve` (v3)

`v1` and `v2` of this dispatcher seeded the IVP step with a plain
`path_network_solve` over the *entire* grid.  That is wrong in the
sense FW 2011 already flagged at md:401:

> "smooth regions are unstable regions for any IVP solver, and the
> associated loss of accuracy ought not be carried back into the pole
> field."

When a per-row BVP later read its Dirichlet BCs from those smooth-band
IVP cells, the BCs were quietly corrupted, Newton diverged on the
first non-convergent row, and the dispatcher threw — fail-loud per
Rule 1, but for the wrong reason (the input was the problem, not the
solver).  See bead `padetaylor-0tj` for the symptom.  v3 routes the
IVP through `EdgeGatedSolve.edge_gated_pole_field_solve`, which
confines the IVP to cells the edge detector confirms as pole-field;
smooth-band cells are left `NaN + NaN·im` by construction (Rule 1
honest: an unvisited cell carries no value).  Those cells then get
filled by the per-row BVP fill in Step 3, with BCs read from
genuine pole-field cells at the smooth band's edges.

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

  1. **Edge-gated IVP fill** (default, `mask === nothing`): call
     `edge_gated_pole_field_solve(prob, xs, ys; h, order, edge_level)`.
     The returned `u_grid` carries `NaN + NaN·im` outside the
     edge-confirmed pole field (plus its one-ring frontier);
     `up_grid` is reconstructed by scattering the underlying
     `pn_solution.grid_up` onto the same dilated-field cells.  When
     `mask` is supplied (the FW md:401 manual classification
     workflow), v1's plain `path_network_solve` path is used so the
     caller's classification rules the partition.
  2. **Partition**: when `mask === nothing`, the gated solve's
     `field_mask` (the edge-detector-confirmed pole-field cells) IS
     the partition — only these cells serve as BVP flanks.  Cells in
     the one-ring frontier carry valid IVP values from the gated
     solve but are smooth-classified and so participate as the
     leftmost / rightmost cells of a smooth band rather than as
     flanks.  When `mask` is supplied, the caller's mask rules.
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

## Fail-soft mode (`strict::Bool = true`)

A `strict::Bool = true` kwarg gates the fail-fast response to BVP
non-convergence.

  - `strict = true` (default, backward-compatible): a non-convergent
    BVP row throws `ErrorException`, the v1/v2 contract.
  - `strict = false`: a non-convergent BVP row is caught, the smooth
    cells retain their *input* values (which under the edge-gated
    default are already `NaN + NaN·im` — the honest "no value here"
    signal per Rule 1; under the manual-classification fallback they
    retain the IVP-computed values), and the row's cells are tagged
    `:bvp_fail` instead of `:bvp`.  Downstream callers can then mask
    out `:bvp_fail` cells for plotting/analysis rather than losing
    the whole solve.

The catch is intentionally narrow: only the `ErrorException` whose
message contains the substring `"bvp_solve: Newton did not converge"`
is swallowed.  Any other exception (`ArgumentError` from a malformed
BVP setup, an `InexactError` from a wrong-precision argument, etc.)
rethrows verbatim — fail-soft is for the documented Newton-divergence
case only, not for swallowing unrelated bugs.

## Public API

  - `LatticeSolution{T}` — composed 2D output with per-cell region
    tag (`:ivp`, `:bvp`, `:bvp_fail`, `:ivp_only`), the partition
    mask, and the underlying IVP + BVP sub-solutions for diagnostic
    inspection.

  - `lattice_dispatch_solve(prob, bvp_f, bvp_∂f_∂u, xs, ys; kwargs...)` —
    the public driver.  Returns `LatticeSolution{T}`.

  - The `mask` kwarg accepts a pre-computed `BitMatrix` for callers who
    want to bypass the automatic edge-gated classification (the FW
    2011 line 401 "manual classification" workflow).  Supplying `mask`
    also routes Step 1 through the v1/v2 plain-`path_network_solve`
    path so the caller's classification is the ground truth.

## Fail-fast contract (CLAUDE.md Rule 1)

Throws `ArgumentError` on:
  - `length(xs) < 3` or `length(ys) < 3`.
  - `xs` or `ys` not strictly increasing.
  - `prob.y0` not a 2-tuple (lattice dispatcher requires 2nd-order IVP).
  - `mask` supplied with shape ≠ `(nx, ny)`.
  - `edge_level` invalid (Inf, NaN).
  - `N_bvp < 4`.

Throws `ErrorException` on:
  - Path-network failure (delegated from `path_network_solve` /
    `edge_gated_pole_field_solve`).
  - BVP non-convergence on any row (delegated from `bvp_solve`), *only*
    when `strict = true` (the default).  With `strict = false`, the
    non-convergent row is tagged `:bvp_fail` and the solve continues.

## References

  - Fornberg & Weideman, *A numerical methodology for the Painlevé
    equations*, J. Comput. Phys. 230 (2011) 5957–5973, §4.4 +
    line 190 (the "161 BVPs, one per grid line" line); md:401 (the
    "smooth regions are unstable for any IVP solver" line that
    motivates the v3 edge-gated IVP source).
  - `references/markdown/FW2011_*.md:190, 218-222, 249-261, 401`.
  - `docs/adr/0004-path-network-architecture.md`.
  - `docs/adr/0017-lattice-dispatcher-strict-mode.md` (this v3 fix).
  - Bead `padetaylor-0tj` (the BC-corruption symptom).
"""
module LatticeDispatcher

using ..Problems:       PadeTaylorProblem
using ..PathNetwork:    path_network_solve, PathNetworkSolution
using ..BVP:            bvp_solve, BVPSolution
using ..EdgeDetector:   pole_field_mask
using ..EdgeGatedSolve: edge_gated_pole_field_solve, EdgeGatedSolution

export LatticeSolution, lattice_dispatch_solve

"""
    LatticeSolution{T}

Composed 2D-grid output from `lattice_dispatch_solve`.  Fields:

  - `xs::AbstractVector{T}`, `ys::AbstractVector{T}` — the lattice axes.
  - `grid_z::Matrix{Complex{T}}` — `grid_z[i, j] = xs[i] + im·ys[j]`.
  - `grid_u`, `grid_up::Matrix{Complex{T}}` — stitched solution and
    derivative.  Cells tagged `:bvp` carry the BVP barycentric
    interpolant; cells tagged `:ivp` or `:ivp_only` carry the
    PathNetwork output.  Cells tagged `:bvp_fail` (fail-soft mode
    only — see the module docstring) retain their *input* values
    — `NaN + NaN·im` under the default edge-gated IVP path, the
    IVP-computed value under the manual-classification fallback.
  - `mask::BitMatrix` — the partition.  `true` ⟺ classified
    pole-field.  Under the default edge-gated IVP path this is the
    `EdgeGatedSolve.field_mask` dilated by one ring; under the manual
    fallback it is the caller-supplied `mask` (or the
    `pole_field_mask` result if the caller passed `nothing`).
  - `region_tag::Matrix{Symbol}` — per-cell tag in `{:ivp, :bvp,
    :bvp_fail, :ivp_only}`.  See the module docstring for definitions.
    `:bvp_fail` only appears when `strict = false` and a per-row BVP
    Newton iteration failed to converge.
  - `pn_solution::PathNetworkSolution{T}` — the underlying IVP solve.
    Under the default edge-gated path this is the `pn_solution` field
    of the intermediate `EdgeGatedSolution`; under the manual fallback
    it is the plain `path_network_solve` output.
  - `bvp_solutions::Vector{BVPSolution{T, Complex{T}}}` — one BVP per
    successfully-bridged smooth run.  A row that failed under
    `strict = false` contributes nothing here — the caller diagnoses
    failure via `region_tag .== :bvp_fail`.
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
                           mask = nothing,
                           strict = true) -> LatticeSolution

See module docstring for the algorithm.  `order` defaults to
`prob.order` and is threaded into the IVP step for the path network.

Kwargs:
  - `h_path` — IVP step size for the path-network walker.  Forwarded
    as `h_path` to `path_network_solve` (manual-fallback path) and as
    `h` to `edge_gated_pole_field_solve` (default path — the kwarg
    name on edge-gated is `h`, not `h_path`).
  - `order` — Padé / Taylor order.  Forwarded to whichever IVP driver
    runs.
  - `edge_level` — edge-detector threshold.  Forwarded to either the
    edge-gated driver or `pole_field_mask` per the path.
  - `N_bvp`, `bvp_tol` — BVP collocation count and Newton tolerance.
  - `mask` — pre-computed `BitMatrix` overriding the auto-classification
    (FW md:401 manual workflow).  When supplied, the IVP step falls
    back to the v1/v2 plain `path_network_solve` over the full grid;
    when `nothing`, the edge-gated v3 path is taken.
  - `strict` — `true` (default) throws on per-row BVP non-convergence
    (back-compat); `false` swallows the *exact* "Newton did not
    converge" exception and tags the affected cells `:bvp_fail`.  See
    the module docstring "Fail-soft mode" section.
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
                                mask              = nothing,
                                strict::Bool      = true)

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

    grid_mat = CT[xs[i] + im * ys[j] for i in 1:nx, j in 1:ny]

    # --- Step 1+2: IVP fill + partition -------------------------------------
    #
    # Two paths.  When the caller supplies `mask`, they have already
    # committed to a classification (FW md:401 manual workflow), so we
    # run the v1/v2 plain `path_network_solve` over the full grid and
    # let that mask drive Step 3.  When `mask === nothing` (the common
    # case), we route through `edge_gated_pole_field_solve` so the IVP
    # never visits smooth cells — addressing the BC corruption that
    # caused bead `padetaylor-0tj` (FW2011_*.md:401 verbatim: "smooth
    # regions are unstable regions for any IVP solver, and the
    # associated loss of accuracy ought not be carried back into the
    # pole field").
    local pn_sol::PathNetworkSolution{T}
    local u_grid::Matrix{CT}
    local up_grid::Matrix{CT}
    local mask_used::BitMatrix

    if mask === nothing
        gated = edge_gated_pole_field_solve(prob, xs, ys;
                                            h          = h_path,
                                            order      = order,
                                            edge_level = edge_level)
        pn_sol  = gated.pn_solution
        u_grid  = copy(gated.u_grid)
        # Re-scatter the gated PathNetwork's `grid_up` onto the (nx, ny)
        # lattice.  The gated solve's `pn_solution.grid_z[k]` corresponds
        # to (xs[i] + im·ys[j]) for some (i, j) in the dilated field; we
        # iterate the lattice in column-major order and the field-mask
        # to recover the same indexing convention used internally by
        # `EdgeGatedSolve._solve_targets`.  Cells outside that final
        # one-ring dilation remain `NaN + NaN·im` — the honest signal
        # that the IVP did not visit them (Rule 1).
        up_grid = fill(CT(NaN, NaN), nx, ny)
        final_targets = _dilate_one(gated.field_mask)
        k = 0
        @inbounds for j in 1:ny, i in 1:nx
            final_targets[i, j] || continue
            k += 1
            up_grid[i, j] = pn_sol.grid_up[k]
        end
        # The partition for Step 3 is the edge-confirmed `field_mask`
        # itself, NOT the one-ring dilation: only cells the edge
        # detector classified as pole-field should serve as BVP flanks.
        # The one-ring dilation matters for the IVP `u_grid` / `up_grid`
        # scatter (the gated solve walks one ring beyond the field, so
        # those cells carry valid IVP values), but those frontier cells
        # are smooth-classified and the per-row scan treats them as
        # part of the smooth band to be bridged.
        mask_used = copy(gated.field_mask)
    else
        size(mask) == (nx, ny) || throw(ArgumentError(
            "lattice_dispatch_solve: mask must have shape ($nx, $ny); " *
            "got $(size(mask))."))
        pn_sol = path_network_solve(prob, vec(grid_mat);
                                    h = h_path, order = order)
        u_grid  = Matrix{CT}(reshape(pn_sol.grid_u,  (nx, ny)))
        up_grid = Matrix{CT}(reshape(pn_sol.grid_up, (nx, ny)))
        mask_used = convert(BitMatrix, mask)
    end

    # --- Step 3: Per-row BVP fill -------------------------------------------
    # Default tag: cells with mask=true are :ivp; mask=false cells are
    # :ivp_only until/unless they get bridged to :bvp (success) or
    # :bvp_fail (fail-soft mode + Newton non-convergence).
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

                # Strict (default) re-throws bvp_solve's
                # `ErrorException("bvp_solve: Newton did not converge…")`
                # verbatim — preserving the v1/v2 fail-fast contract.
                # Fail-soft (`strict = false`) catches *exactly* that
                # exception (matched on message substring — the
                # narrowest safe net) and tags the cells `:bvp_fail`;
                # any other exception rethrows so we never silently
                # swallow an unrelated bug.
                bvp_sol = try
                    bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b, u_a, u_b;
                              N = N_bvp, tol = bvp_tol)
                catch err
                    if !strict && err isa ErrorException &&
                       occursin("bvp_solve: Newton did not converge", err.msg)
                        for k in run_start:run_end
                            region_tag[k, j] = :bvp_fail
                        end
                        nothing
                    else
                        rethrow()
                    end
                end

                if bvp_sol !== nothing
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
    end

    return LatticeSolution{T}(xs, ys, grid_mat, u_grid, up_grid,
                              mask_used, region_tag, pn_sol, bvp_solutions)
end

# Chebyshev (8-connected) dilation by one ring — small helper so we
# don't reach into `EdgeGatedSolve`'s private `_dilate`.  Used to lift
# `gated.field_mask` to the final one-ring partition (matches the final
# `_solve_targets` call inside `edge_gated_pole_field_solve`).
function _dilate_one(mask::BitMatrix)
    nx, ny = size(mask)
    out = copy(mask)
    @inbounds for j in 1:ny, i in 1:nx
        mask[i, j] && continue
        for dj in -1:1, di in -1:1
            (di == 0 && dj == 0) && continue
            ii = i + di; jj = j + dj
            if 1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny && mask[ii, jj]
                out[i, j] = true
                break
            end
        end
    end
    return out
end

end # module LatticeDispatcher
