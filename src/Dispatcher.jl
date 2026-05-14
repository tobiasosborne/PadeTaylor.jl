"""
    PadeTaylor.Dispatcher

Tier-3 composition layer (Phase 12 v1, bead `padetaylor-8lk`) — stitches
IVP path-network segments and BVP Chebyshev-Newton segments into a single
ordered chain along a 1D path through the complex plane, with a
**junction derivative-match diagnostic** per FW 2011 §4.4
(`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:192, 249-261`).

This module is the **second-from-top public-API layer**: it consumes
`PathNetwork.path_network_solve` and `BVP.bvp_solve`, both already
shipped at commits `0ada60f` (Phase 10) and `cc7d8ca` (Phase 11).  No
new algorithm is invented here; the contribution is the orchestration
discipline and the diagnostic contract.

## The algorithm (FW 2011 §4.4, lines 249-261)

FW 2011 use the IVP↔BVP composition in the simplest form:
  1. The user specifies a sequence of complex-plane segments, each
     marked **IVP** (pole-field / path-network applicable) or **BVP**
     (smooth region / Chebyshev-Newton applicable).
  2. The dispatcher walks the chain segment-by-segment, propagating
     state across junctions:
       - **IVP→BVP**: the BVP's left Dirichlet BC `u_a` is taken from
         the IVP's terminal `u`; the user supplies `u_b` (the right BC).
       - **BVP→IVP**: the next IVP's IC pair `(u_0, u_0')` is taken
         from the BVP's right endpoint `(u_BVP(z_b), u'_BVP(z_b))`
         (barycentric-recovered derivative).
  3. At each interior junction (IVP→BVP or BVP→IVP), the dispatcher
     records `(|Δu|, |Δu'|)` as `junction_match[k]`:
       - `Δu` is zero by construction (BC enforced from terminus).
       - `Δu'` is the **diagnostic**: the BVP's barycentric derivative
         at its left endpoint differs from the IVP's analytic
         derivative there because the BVP is a finite-N spectral
         approximation.  FW 2011 line 192: a typical Δu' is `1e-7` to
         `1e-8`; an increase in `N` is indicated if not.

## Why this is "FW 2011 §4.4 verbatim", not an invented algorithm

FW 2011 line 192 explicitly states (verbatim): "If the agreement is
not to within some tolerance (typically set to 10⁻⁷ or 10⁻⁸), an
increase in N, the degree of the interpolant, is indicated."  The
**tolerance is a quality indicator**, not a convergence test; the
diagnostic is recorded but the dispatcher does NOT throw on mismatch
in default mode.  Users can choose to refine `N_bvp` if they see
junction mismatches.  This matches FW's intended workflow.

The strict-throw mode (`derivative_match_strict = true`) is provided
for callers who want a hard guarantee — e.g., regression tests.

## v1 scope (1D ordered chain)

Phase 12 v1 ships the **1D ordered-chain dispatcher** described above.
Each segment is an `IVPSegment` or a `BVPSegment` struct; the chain is
a `Vector` of these, processed in order.

**Deferred to v2** (bead `padetaylor-k31`): the 2D lattice dispatcher
with automatic 5-point Laplacian edge detection per FW 2011 §3.2.2
(`FW2011...md:204-208`).  v2 will partition a 2D lattice into
pole-field cells (IVP) and smooth cells (BVP), solve one BVP per row
of contiguous smooth cells (FW 2011 line 190: "161 separate BVP
solutions; one for each grid line"), and stitch the result into a
single dense 2D output with per-cell region tagging.

## Public API

  - `IVPSegment(z_end, dense_grid)` — IVP segment with terminal
    coordinate `z_end` and a user-supplied evaluation grid.
  - `BVPSegment(z_end, u_b, dense_grid; initial_guess = nothing)` —
    BVP segment.  The left BC `u_a` is taken from the previous
    segment's terminus (or `prob.y0[1]` for the first segment).
  - `DispatcherSolution{T}` — composed result.
  - `dispatch_solve(prob_ivp, f_bvp, ∂f_∂u_bvp, segments; kwargs...)` —
    the public driver.

## Fail-fast contract (CLAUDE.md Rule 1)

Throws `ArgumentError` on empty `segments`, negative `h`, negative
`derivative_match_tol`, `N_bvp < 4`.  Throws `ErrorException` on
non-convergent BVP segments (delegated from `bvp_solve`) and on
path-network failures (delegated from `path_network_solve`).
With `derivative_match_strict = true`, also throws on any junction
mismatch exceeding `derivative_match_tol`.

## References

  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-208`
    — §3.1 path-network + §3.2 BVP solver.
  - Same file, lines 249-261 — §4.4 IVP/BVP composition.
  - Same file, line 192 — the 1e-7 / 1e-8 derivative-match tolerance.
  - `docs/adr/0004-path-network-architecture.md` — Tier-3 deferral plan.
  - `docs/unified_path_network_spec.md §8` — dispatcher algorithm
    synthesis from the 5-paper FW-family literature dive.
  - `docs/figure_catalogue.md §1` — Tier-3 figures unlocked by v1.
  - `src/PathNetwork.jl`, `src/BVP.jl` — the consumed modules.
"""
module Dispatcher

using ..Problems:    PadeTaylorProblem
using ..PathNetwork: PathNetworkSolution, path_network_solve
using ..BVP:         BVPSolution, bvp_solve

export IVPSegment, BVPSegment, DispatcherSolution, dispatch_solve

# =============================================================================
# Segment specs
# =============================================================================

"""
    IVPSegment{T}(z_end, dense_grid)

Specifies one IVP segment of a Dispatcher chain.  The segment starts at
the previous segment's terminus (or `prob.y0` for the first segment)
and walks the FW 2011 §3.1 path-network until reaching within `h` of
`z_end`.  `dense_grid` is the user's evaluation grid within this
segment; values are filled via Stage-2 barycentric Padé evaluation.
"""
struct IVPSegment{T <: AbstractFloat}
    z_end      :: Complex{T}
    dense_grid :: Vector{Complex{T}}
end

# Convenience constructor for ComplexF64.
IVPSegment(z_end::Complex{T}, dense_grid::AbstractVector{<:Complex}) where T <: AbstractFloat =
    IVPSegment{T}(z_end, Vector{Complex{T}}(dense_grid))

"""
    BVPSegment{T, F}(z_end, u_b, dense_grid, initial_guess)

Specifies one BVP segment.  The left Dirichlet BC `u_a` is taken from
the previous segment's terminus (the dispatcher fills it in; the user
does NOT supply it).  `u_b` is the right Dirichlet BC at `z_end`.
`dense_grid` is the user's evaluation grid; values are filled via
barycentric interpolation from the Chebyshev nodes.  `initial_guess`
is a callable `z -> Complex{T}` (or `nothing` for the default linear
ramp from `u_a` to `u_b`).
"""
struct BVPSegment{T <: AbstractFloat, F}
    z_end         :: Complex{T}
    u_b           :: Complex{T}
    dense_grid    :: Vector{Complex{T}}
    initial_guess :: F
end

BVPSegment(z_end::Complex{T}, u_b::Complex,
           dense_grid::AbstractVector{<:Complex};
           initial_guess = nothing) where T <: AbstractFloat =
    BVPSegment{T, typeof(initial_guess)}(
        z_end, Complex{T}(u_b), Vector{Complex{T}}(dense_grid), initial_guess)

# =============================================================================
# Solution container
# =============================================================================

"""
    DispatcherSolution{T}

Composed result of `dispatch_solve`.  Carries:
  - `segment_kinds`     — `:ivp` / `:bvp` per segment, in chain order.
  - `ivp_solutions`     — vector of `PathNetworkSolution{T}`, one per `:ivp` segment.
  - `bvp_solutions`     — vector of `BVPSolution{T, Complex{T}}`, one per `:bvp` segment.
  - `junction_z`        — complex coordinates of interior junctions.
  - `junction_match`    — `(Δu, Δu')` per junction, both ≥ 0 absolute.
  - `derivative_match_tol` — the diagnostic threshold used.
  - `grid_z`, `grid_u`, `grid_up`, `grid_region` — concatenated dense
    output across all segments' `dense_grid`s, in chain order.
"""
struct DispatcherSolution{T <: AbstractFloat}
    segment_kinds         :: Vector{Symbol}
    ivp_solutions         :: Vector{PathNetworkSolution{T}}
    bvp_solutions         :: Vector{BVPSolution{T, Complex{T}}}
    junction_z            :: Vector{Complex{T}}
    junction_match        :: Vector{Tuple{T, T}}
    derivative_match_tol  :: T
    grid_z                :: Vector{Complex{T}}
    grid_u                :: Vector{Complex{T}}
    grid_up               :: Vector{Complex{T}}
    grid_region           :: Vector{Symbol}
end

# =============================================================================
# Public driver (stub — implementation follows in Phase 12 GREEN commit)
# =============================================================================

"""
    dispatch_solve(prob_ivp, f_bvp, ∂f_∂u_bvp, segments;
                   derivative_match_tol = 1e-7,
                   derivative_match_strict = false,
                   h = 0.5, order = prob_ivp.order,
                   N_bvp = 20, bvp_tol = nothing, bvp_maxiter = 10,
                   wedge_angles = nothing) -> DispatcherSolution

Compose an ordered chain of IVP and BVP segments into a single dense
output.  See module docstring for the FW 2011 §4.4 algorithm.

Arguments:
  - `prob_ivp::PadeTaylorProblem` — IVP problem (provides `f`, IC `y0`,
    `order`).  The first segment uses `prob_ivp.zspan[1]` as its
    starting coordinate and `prob_ivp.y0` as IC.
  - `f_bvp(z, u) -> Number` — 1st-order BVP rhs (caller must adapt
    from `prob_ivp.f`'s 2nd-order signature; typically equals `f_ivp`
    on its first two args if the ODE is autonomous in `u'`).
  - `∂f_∂u_bvp(z, u) -> Number` — analytic partial.
  - `segments` — `AbstractVector` of `IVPSegment` / `BVPSegment` structs.

Kwargs control IVP step parameters, BVP collocation degree, junction
tolerance + strict-mode behaviour.

Throws on bad inputs per CLAUDE.md Rule 1.  Propagates underlying
exceptions from `path_network_solve` / `bvp_solve`.
"""
function dispatch_solve(prob_ivp::PadeTaylorProblem,
                        f_bvp,
                        ∂f_∂u_bvp,
                        segments::AbstractVector;
                        derivative_match_tol::Real = 1e-7,
                        derivative_match_strict::Bool = false,
                        h::Real = 0.5,
                        order::Integer = prob_ivp.order,
                        N_bvp::Integer = 20,
                        bvp_tol = nothing,
                        bvp_maxiter::Integer = 10,
                        wedge_angles = nothing)

    # ---- validate ----------------------------------------------------------
    isempty(segments) && throw(ArgumentError(
        "dispatch_solve: segments must not be empty.  " *
        "Suggestion: provide at least one IVPSegment or BVPSegment."))
    derivative_match_tol ≥ 0 || throw(ArgumentError(
        "dispatch_solve: derivative_match_tol must be ≥ 0 (got $derivative_match_tol)."))
    h > 0 || throw(ArgumentError(
        "dispatch_solve: h must be positive (got $h)."))
    N_bvp ≥ 4 || throw(ArgumentError(
        "dispatch_solve: N_bvp must be ≥ 4 (got $N_bvp).  " *
        "Suggestion: increase to ≥ 10 for any nontrivial nonlinear problem; " *
        "Chebyshev-Newton needs at least 3 interior collocation nodes."))
    order ≥ 2 || throw(ArgumentError(
        "dispatch_solve: order must be ≥ 2 (got $order)."))

    # ---- type promotion ----------------------------------------------------
    # The chain is anchored on the IVP problem's z-coordinate type.
    z_start_real = real(prob_ivp.zspan[1])
    T = float(typeof(z_start_real))
    CT = Complex{T}
    h_T = T(h)
    tol_T = T(derivative_match_tol)

    # ---- accumulators ------------------------------------------------------
    segment_kinds  = Symbol[]
    ivp_solutions  = PathNetworkSolution{T}[]
    bvp_solutions  = BVPSolution{T, CT}[]
    junction_z     = CT[]
    junction_match = Tuple{T, T}[]
    grid_z         = CT[]
    grid_u         = CT[]
    grid_up        = CT[]
    grid_region    = Symbol[]

    # ---- chain state -------------------------------------------------------
    # `cur_*` holds the running terminus state (z, u, u') propagated from
    # one segment to the next.  Initialised from the IVP problem's IC.
    cur_z  = CT(prob_ivp.zspan[1])
    cur_u  = CT(prob_ivp.y0[1])
    cur_up = CT(prob_ivp.y0[2])

    # ---- walk the chain ----------------------------------------------------
    for (k, seg) in enumerate(segments)
        new_cur_u  :: CT = cur_u
        new_cur_up :: CT = cur_up
        new_cur_z  :: CT = cur_z

        if seg isa IVPSegment
            # Construct a fresh IVP problem for this segment: same `f`,
            # IC at the running terminus, zspan = (cur_z, seg.z_end).
            # The original prob_ivp.zspan[2] is irrelevant here — the
            # segment specifies its own endpoint.
            ivp_prob = PadeTaylorProblem(prob_ivp.f, (cur_u, cur_up),
                                         (cur_z, CT(seg.z_end));
                                         order = order)
            # Path-network walks to every target in `targets` (the user's
            # dense_grid + the terminus, so we can extract terminal state).
            targets = Vector{CT}(vcat(seg.dense_grid, [seg.z_end]))
            ivp_sol = if wedge_angles === nothing
                path_network_solve(ivp_prob, targets; h = h_T)
            else
                path_network_solve(ivp_prob, targets;
                                   h = h_T, wedge_angles = wedge_angles)
            end
            push!(ivp_solutions, ivp_sol)
            push!(segment_kinds, :ivp)

            # Append dense-grid outputs (exclude the terminus we appended).
            n_dense = length(seg.dense_grid)
            for i in 1:n_dense
                push!(grid_z, ivp_sol.grid_z[i])
                push!(grid_u, ivp_sol.grid_u[i])
                push!(grid_up, ivp_sol.grid_up[i])
                push!(grid_region, :ivp)
            end

            # Terminal state: the LAST target, which is seg.z_end by construction.
            new_cur_z  = ivp_sol.grid_z[n_dense + 1]
            new_cur_u  = ivp_sol.grid_u[n_dense + 1]
            new_cur_up = ivp_sol.grid_up[n_dense + 1]

        elseif seg isa BVPSegment
            # BVP segment with u_a = cur_u (BC enforced from previous terminus).
            ig = seg.initial_guess
            bvp_sol = if ig === nothing && bvp_tol === nothing
                bvp_solve(f_bvp, ∂f_∂u_bvp, cur_z, CT(seg.z_end),
                          cur_u, seg.u_b;
                          N = N_bvp, maxiter = bvp_maxiter)
            elseif ig === nothing
                bvp_solve(f_bvp, ∂f_∂u_bvp, cur_z, CT(seg.z_end),
                          cur_u, seg.u_b;
                          N = N_bvp, maxiter = bvp_maxiter, tol = bvp_tol)
            elseif bvp_tol === nothing
                bvp_solve(f_bvp, ∂f_∂u_bvp, cur_z, CT(seg.z_end),
                          cur_u, seg.u_b;
                          N = N_bvp, maxiter = bvp_maxiter, initial_guess = ig)
            else
                bvp_solve(f_bvp, ∂f_∂u_bvp, cur_z, CT(seg.z_end),
                          cur_u, seg.u_b;
                          N = N_bvp, maxiter = bvp_maxiter,
                          tol = bvp_tol, initial_guess = ig)
            end
            push!(bvp_solutions, bvp_sol)
            push!(segment_kinds, :bvp)

            # If this is NOT the first segment, record the junction at cur_z.
            # Δu = |cur_u - u_BVP(cur_z)|; equals 0 by BC construction but
            # we compute explicitly so a mutation breaking the BC enforcement
            # would bite.  Δu' = |cur_up - u'_BVP(cur_z)| via barycentric.
            if k > 1
                (u_bvp_left, up_bvp_left) = bvp_sol(cur_z)
                Δu  = T(abs(cur_u  - u_bvp_left))
                Δup = T(abs(cur_up - up_bvp_left))
                push!(junction_z, cur_z)
                push!(junction_match, (Δu, Δup))
                if derivative_match_strict && (Δu > tol_T || Δup > tol_T)
                    throw(ErrorException(
                        "dispatch_solve: derivative-match strict-mode " *
                        "violation at junction z=$cur_z.  |Δu| = $Δu, " *
                        "|Δu'| = $Δup > tol = $tol_T.  " *
                        "Suggestion: (a) increase N_bvp (higher spectral " *
                        "resolution), (b) refine h (denser IVP segments), " *
                        "(c) raise derivative_match_tol, or (d) disable " *
                        "strict mode and inspect junction_match diagnostics."))
                end
            end

            # Append dense-grid outputs via barycentric eval.
            for z in seg.dense_grid
                (u_z, up_z) = bvp_sol(z)
                push!(grid_z, z)
                push!(grid_u, u_z)
                push!(grid_up, up_z)
                push!(grid_region, :bvp)
            end

            # Terminal state at z_end via barycentric.
            new_cur_z = CT(seg.z_end)
            (u_at_end, up_at_end) = bvp_sol(seg.z_end)
            new_cur_u  = u_at_end
            new_cur_up = up_at_end

        else
            throw(ArgumentError(
                "dispatch_solve: unknown segment type $(typeof(seg)).  " *
                "Suggestion: use IVPSegment or BVPSegment."))
        end

        # ---- record IVP-side junction (k>1, seg isa IVPSegment) -----------
        # For X→IVP transitions, the IVP starts with IC = (cur_u, cur_up) =
        # previous terminus.  Δu = Δu' = 0 by construction.  We still record
        # the junction so the count is right and mutation-proof procedures
        # can target it.
        if k > 1 && seg isa IVPSegment
            push!(junction_z, cur_z)
            push!(junction_match, (zero(T), zero(T)))
            # No strict-mode trigger here — zero match never fails.
        end

        cur_z, cur_u, cur_up = new_cur_z, new_cur_u, new_cur_up
    end

    return DispatcherSolution{T}(segment_kinds, ivp_solutions, bvp_solutions,
                                 junction_z, junction_match, tol_T,
                                 grid_z, grid_u, grid_up, grid_region)
end

end # module Dispatcher
