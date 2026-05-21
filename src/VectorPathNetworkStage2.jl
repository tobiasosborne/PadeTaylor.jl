"""
    PadeTaylor.VectorPathNetworkStage2

The **Stage-2 fine-grid fill** of the vector path-network — the dense
evaluation step `VectorPathNetwork.vector_path_network_solve` delegates
to when the caller passes a `fine_grid` (bead `padetaylor-0ln.20`).

## Why this is a separate module

V7 (`src/VectorPathNetwork.jl`, bead `padetaylor-0ln.16`) shipped only
the Stage-1 *tree* construction — the wedge walk that builds a visited
tree of per-node shared-`Q` approximants.  The `A_4⁽¹⁾`/`P_I⁽²⁾`
pole-field *scatter* figures need only that: the pole locations are the
roots of each node's shared denominator `Q`, extracted by
`VectorPoleField`.

The headline `P_I⁽²⁾` tritronquée *surface* figure needs more — a dense
filled heatmap of `y(z)` over the pole-rich wedge.  That is Stage 2.
It is the direct vector lift of v0.1's
`PathNetwork.path_network_solve` Stage-2 step (`src/PathNetwork.jl:29-32`,
"Stage 2 — fine-grid extrapolation"; the fill loop at `:669-693`).

Stage 2 lives in its own module purely for the CLAUDE.md Rule 6 200
effective-LOC cap: `VectorPathNetwork.jl` was already ~145 effective
LOC after V7, and folding the Stage-2 fill + its Horner evaluator + the
expanded docstrings inline pushed it to the cap.  The split is along
the natural seam — Stage 1 builds the tree, Stage 2 reads it — and is
*additive*: `VectorPathNetwork`'s public API is unchanged bar the two
new `vector_path_network_solve` kwargs.

## No dependency on `VectorPathNetwork`

This module is deliberately free of any `VectorPathNetwork` import, so
it can be `include`d *before* `VectorPathNetwork` and the latter can
`using ..VectorPathNetworkStage2: _stage2_fill` at its own load time.
The one thing Stage 2 needs from Stage 1 — the nearest-visited-node
scan — is passed in as a *function argument* (`nearest_visited`) rather
than imported, so the dependency arrow points one way only.

## The fill (mirrors v0.1 `PathNetwork`'s Stage-2 loop)

For each fine-grid point `z_f`:

  1. find the nearest visited node `idx_v` via the supplied
     `nearest_visited` scan (the same lexicographic-tiebreak
     nearest-node lookup Stage 1 uses);
  2. rescale to `t = (z_f − z_v) / h_v` with the *real* canonical step
     `h_v = real(visited_h[idx_v])`.  `visited_h` is typed `Complex{T}`,
     but every node's canonical shared-`Q` store is built with a real
     step (`VectorPathNetwork._canonical_pade`), so the imaginary part
     is `0` and `real(...)` recovers the magnitude the rescaled
     variable `t` is defined against;
  3. if `!extrapolate` and `abs(z_f − z_v) > h_v`, the grid point lies
     outside the disc of validity of every nearby local Padé — the
     slot gets a `d`-vector of `NaN + NaN·im` (fail-soft per ADR-0015 /
     CLAUDE.md Rule 1; **not** a throw — the visited tree is honest,
     the grid point is simply uncovered);
  4. otherwise evaluate component `i` as
     `Pᵢ(t) / Q(t)` — every component over the *same* shared `Q(t)`,
     the v0.2 keystone applied at evaluation time, exactly
     `VectorProblems`' dense callable (`src/VectorProblems.jl:310-330`).

The `extrapolate` flag has the ADR-0015 semantics verbatim: default
`false` honours the fail-soft Rule-1 contract; `true` skips the
disc-radius check so figure scripts get a gap-free panel, at the cost
of degraded accuracy past `|t| = 1`.

## References

  - `src/VectorPathNetwork.jl` — the Stage-1 walk that produces the
    visited tree this module fills against; the `vector_path_network_solve`
    driver that delegates here.
  - `src/PathNetwork.jl` — the v0.1 scalar Stage-2 fill (`:669-693`)
    this module lifts to the shared-`Q` vector case.
  - `src/VectorProblems.jl` — the `_eval_poly` Horner evaluator
    (`:339-345`) and the `Pᵢ(t)/Q(t)` dense-evaluation pattern.
  - `docs/adr/0015-extrapolate-stage-2.md` — the `extrapolate`
    disc-radius policy.
  - `docs/adr/0019-shared-denominator-pade.md` — the shared `Q`.
"""
module VectorPathNetworkStage2

# -----------------------------------------------------------------------------
# Stage-2 fine-grid fill
# -----------------------------------------------------------------------------

"""
    _stage2_fill(grid_z, visited_z, visited_h, visited_numerators,
                 visited_denominator, extrapolate, ::Type{C},
                 nearest_visited) -> Vector{Vector{C}}

Evaluate the dense vector solution `y(z)` at every point of `grid_z`
against the nearest visited node's canonical shared-`Q` approximant.

`nearest_visited` is the Stage-1 nearest-visited-node scan
(`VectorPathNetwork._nearest_visited`), passed in so this module needs
no `VectorPathNetwork` import (see the module docstring).  It is called
as `nearest_visited(visited_z, z_f) -> Int`.

Returns one `d`-vector per grid point.  A grid point more than
`real(visited_h[idx])` from its nearest visited node gets a `d`-vector
of `NaN + NaN·im` unless `extrapolate` is `true` (ADR-0015).
"""
function _stage2_fill(grid_z::Vector{C},
                      visited_z::Vector{C},
                      visited_h::Vector{C},
                      visited_numerators::Vector{Vector{Vector{C}}},
                      visited_denominator::Vector{Vector{C}},
                      extrapolate::Bool,
                      ::Type{C},
                      nearest_visited) where {C}
    T      = real(C)
    nan_C  = C(T(NaN), T(NaN))
    grid_y = Vector{Vector{C}}(undef, length(grid_z))
    for (i, z_f) in enumerate(grid_z)
        idx_v      = nearest_visited(visited_z, z_f)
        z_v        = visited_z[idx_v]
        h_v        = real(visited_h[idx_v])      # canonical step is real
        numerators = visited_numerators[idx_v]
        if !extrapolate && abs(z_f - z_v) > h_v
            # Grid point outside every nearby disc of validity — fail-soft.
            grid_y[i] = fill(nan_C, length(numerators))
        else
            t   = (z_f - z_v) / h_v
            q_t = _eval_poly(visited_denominator[idx_v], t)
            grid_y[i] = C[_eval_poly(num, t) / q_t for num in numerators]
        end
    end
    return grid_y
end

"""
    _eval_poly(c::AbstractVector{C}, t) -> C

Horner evaluation of the polynomial `Σ c[k+1] tᵏ` (coefficients
low-to-high) — the convention the shared-`Q` numerators and denominator
are stored in.  A private copy of `VectorProblems._eval_poly`
(`src/VectorProblems.jl:339-345`): copied rather than imported so this
module does not reach into a sibling module's private symbol
(CLAUDE.md — module-boundary discipline).
"""
function _eval_poly(c::AbstractVector{C}, t) where {C}
    s = zero(C)
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

end # module VectorPathNetworkStage2
