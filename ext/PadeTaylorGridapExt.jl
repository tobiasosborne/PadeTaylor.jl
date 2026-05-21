"""
    PadeTaylorGridapExt

Package-extension providing the **Gridap.jl finite-element Laplace
backend** — voter (3) of the triple-method majority-vote
harmonic-extension fill (ADR-0024).  Loaded automatically when both
`PadeTaylor` and `Gridap` are present:

```julia
using PadeTaylor, Gridap

# φ = x² − y² is harmonic; supply its exact edge trace as boundary data.
g = (x, y) -> x^2 - y^2
sol = laplace2d_solve_gridap(y -> g(0.0, y), y -> g(1.0, y),
                             x -> g(x, 0.0), x -> g(x, 1.0),
                             (0.0, 1.0), (0.0, 1.0); n_x = 48, n_y = 48)
sol(0.37, 0.61)            # ≈ 0.37² − 0.61²
```

## Why a third, FEM-based voter (ADR-0024)

`Re u` / `Im u` of the analytic `P_I⁽²⁾` tritronquée are harmonic on the
pole-free sector, hence each is the unique solution of a Dirichlet
Laplace problem.  ADR-0024 fixes the discipline of computing that fill
**three algorithmically disjoint ways** and majority-voting per grid
point: (1) ray-fan `vector_bvp_solve`, (2) the in-house 2D-Chebyshev
spectral solve (`laplace2d_solve`, `src/Laplace2D.jl`), and (3) **this
extension's Gridap FEM solve**.  A FEM discretisation error and a
spectral discretisation error are structurally different and cannot hide
the same systematic bug — that disjointness is the whole point of the
third voter.  Gridap is, per ADR-0003, a *weak* dependency: a user who
does not need the FEM cross-check pays neither its precompile time nor
its dependency stack.

## Algorithm — Gridap Dirichlet Laplace on a rectangle

The extension solves the **same** problem as `laplace2d_solve`:
`∇²φ = 0` on `[a,b] × [c,d]` with Dirichlet data on all four edges.

  1. **Mesh.** A `CartesianDiscreteModel` of the rectangle with
     `n_x × n_y` quadrilateral cells.  Cartesian models tag the four
     boundary edges as faces `1…8` ("tag_1"…"tag_8") and the whole
     boundary collectively as `"boundary"` — a continuous boundary
     trace (the exact trace of one harmonic function is automatically
     corner-consistent) is therefore imposed in one stroke with the
     `"boundary"` tag.

  2. **FE space.** A first-order (`order = 1`) Lagrangian `TestFESpace`
     with `dirichlet_tags = "boundary"`, and the matching `TrialFESpace`
     carrying the Dirichlet datum `g(x) = φ_edge` — a single function of
     the physical point that dispatches on which edge `x` lies on.

  3. **Weak form.** Laplace's equation `∇²φ = 0` has weak form
     `∫_Ω ∇u · ∇v dΩ = 0` for all test `v` (no source term).  We
     assemble it on a degree-2 quadrature `Measure` and solve the
     resulting linear `AffineFEOperator` — one linear solve, no Newton
     (the problem is linear), mirroring `laplace2d_solve`.

  4. **Dense output.** The Gridap solution `uh` is a `CellField`, which
     is itself evaluatable at any physical point — but point evaluation
     of a raw `CellField` needs a cell-search.  We instead sample `uh`
     on a Cartesian grid and wrap it in a `Laplace2DSolution`, whose
     tensor-barycentric callable then provides `(x,y) -> φ` with the
     **same interface** as `laplace2d_solve` so F3 treats both Laplace
     voters uniformly.

## Accuracy

FEM with order-`p` Lagrangian elements converges as `O(h^{p+1})` in the
`L²` norm (`h ~ 1/n`).  With `order = 1` and `n ~ 48` cells per axis the
interior error on a smooth harmonic field is `O(10^{-3}…10^{-4})` —
algebraic, NOT spectral (contrast `laplace2d_solve`, which converges
geometrically with a handful of nodes).  This is expected and is the
reason ADR-0024 keeps Gridap only as a *cross-check* voter, never the
sole method.  The test file pins the convergence rate by mesh
refinement.

## References

  - `docs/adr/0024-laplace-harmonic-extension.md` — the triple-method
    majority-vote decision; this extension is voter (3).
  - `docs/adr/0003-extensions-pattern.md` — the weak-dep extension
    pattern this file follows.
  - `src/Laplace2D.jl` — voter (2); `laplace2d_solve`, the in-house
    2D-Chebyshev solver this FEM solve cross-checks, and the
    `Laplace2DSolution` container reused for the dense output.
  - Gridap.jl tutorial 1 ("Poisson equation") — the
    `CartesianDiscreteModel` / `TestFESpace` / `AffineFEOperator`
    idiom adapted here.
"""
module PadeTaylorGridapExt

using PadeTaylor: PadeTaylor
using PadeTaylor.Laplace2D: Laplace2DSolution
import PadeTaylor.Laplace2D: laplace2d_solve_gridap

using Gridap

# =============================================================================
# laplace2d_solve_gridap — the FEM method (voter (3))
# =============================================================================

"""
    laplace2d_solve_gridap(bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi,
                           (a, b), (c, d); n_x = 32, n_y = 32,
                           Nx_out = nothing, Ny_out = nothing)
        -> Laplace2DSolution

Gridap FEM solution of the Dirichlet Laplace problem `∇²φ = 0` on the
rectangle `[a,b] × [c,d]`.  See the `PadeTaylor.laplace2d_solve_gridap`
docstring for the contract; the four `bc_*` arguments and the rectangle
mirror `laplace2d_solve`.

`n_x` / `n_y` are the FEM mesh subdivision counts (the
`CartesianDiscreteModel` has `n_x × n_y` cells).  `Nx_out` / `Ny_out`
control the resolution of the Chebyshev output grid the FEM field is
sampled onto (default: a grid roughly matching the mesh, capped so the
returned `Laplace2DSolution`'s barycentric callable stays well
conditioned).

The four `bc_*` data are **callables** `t -> φ` (a function of the edge
coordinate, exactly as `laplace2d_solve` accepts).  They are reconciled
into a single Dirichlet datum `g(physical point)` that dispatches on
which edge the point lies on.
"""
function laplace2d_solve_gridap(bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi,
                                xrange, yrange;
                                n_x::Integer = 32, n_y::Integer = 32,
                                Nx_out::Union{Integer,Nothing} = nothing,
                                Ny_out::Union{Integer,Nothing} = nothing)

    # ---- fail-fast guards (CLAUDE.md Rule 1) -------------------------------
    n_x ≥ 2 || throw(ArgumentError(
        "laplace2d_solve_gridap: n_x must be ≥ 2 (got $n_x).  Suggestion: " *
        "use ≥ 32 — FEM is O(h^{p+1}), so a coarse mesh gives a poor " *
        "cross-check of the spectral laplace2d_solve."))
    n_y ≥ 2 || throw(ArgumentError(
        "laplace2d_solve_gridap: n_y must be ≥ 2 (got $n_y).  Suggestion: " *
        "use ≥ 32."))

    a, b = xrange
    c, d = yrange
    T = float(promote_type(typeof(a), typeof(b), typeof(c), typeof(d)))
    a, b, c, d = T(a), T(b), T(c), T(d)
    b > a || throw(ArgumentError(
        "laplace2d_solve_gridap: degenerate x-range ($a, $b) — need b > a.  " *
        "Suggestion: pass the rectangle corners as ((a,b),(c,d)) with a < b."))
    d > c || throw(ArgumentError(
        "laplace2d_solve_gridap: degenerate y-range ($c, $d) — need d > c."))

    # ---- single Dirichlet datum g(physical point) -------------------------
    # A Cartesian `"boundary"` tag carries every edge; the datum dispatches
    # on which edge the point sits on.  An interior point never reaches `g`
    # (it is a Dirichlet datum), but a robust default keeps it total.
    half = (b - a + d - c) / 4              # rough scale for the edge test
    atol = max(eps(Float64), 1e-9 * half)
    function g(p)
        x, y = p[1], p[2]
        if abs(x - a) ≤ atol
            return Float64(bc_x_lo(y))
        elseif abs(x - b) ≤ atol
            return Float64(bc_x_hi(y))
        elseif abs(y - c) ≤ atol
            return Float64(bc_y_lo(x))
        elseif abs(y - d) ≤ atol
            return Float64(bc_y_hi(x))
        else
            # Unreachable for a Dirichlet datum; bilinear blend as a guard.
            return Float64(bc_y_lo(x))
        end
    end

    # ---- mesh, FE spaces, weak form, solve --------------------------------
    domain = (Float64(a), Float64(b), Float64(c), Float64(d))
    model  = CartesianDiscreteModel(domain, (Int(n_x), Int(n_y)))

    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = TestFESpace(model, reffe; conformity = :H1,
                    dirichlet_tags = "boundary")
    U = TrialFESpace(V, g)

    Ω  = Triangulation(model)
    dΩ = Measure(Ω, 2)                       # degree-2 quadrature

    # Laplace's equation: weak form ∫ ∇u·∇v dΩ = 0 (no source term).
    aform(u, v) = ∫( ∇(u) ⋅ ∇(v) ) * dΩ
    lform(v)    = ∫( 0.0 * v ) * dΩ

    op = AffineFEOperator(aform, lform, U, V)
    uh = solve(op)                           # one linear solve, no Newton

    # ---- dense output: sample uh on a Chebyshev grid ----------------------
    # Wrap in a `Laplace2DSolution` so F3 treats the FEM voter and the
    # spectral voter (`laplace2d_solve`) through one uniform interface.
    Nx = Nx_out === nothing ? min(Int(n_x), 64) : Int(Nx_out)
    Ny = Ny_out === nothing ? min(Int(n_y), 64) : Int(Ny_out)
    Nx ≥ 3 || (Nx = 3)
    Ny ≥ 3 || (Ny = 3)

    πT = T(π)
    tx = T[cos(T(j) * πT / T(Nx)) for j in 0:Nx]
    ty = T[cos(T(j) * πT / T(Ny)) for j in 0:Ny]
    xs = T[(a + b) / 2 + (b - a) / 2 * t for t in tx]
    ys = T[(c + d) / 2 + (d - c) / 2 * t for t in ty]

    # `uh` is a CellField — evaluate it at physical points.  Searching the
    # mesh cell per point is what `Gridap` does for CellField evaluation;
    # we keep interior points strictly inside the rectangle to avoid the
    # boundary-face ambiguity in the cell search.
    sx = (b - a)
    sy = (d - c)
    nudge = 1e-12
    phi = Matrix{T}(undef, Nx + 1, Ny + 1)
    for j in 1:(Ny + 1), i in 1:(Nx + 1)
        x = xs[i]; y = ys[j]
        # Clamp evaluation points just inside the closed rectangle.
        xe = clamp(Float64(x), Float64(a) + nudge * sx, Float64(b) - nudge * sx)
        ye = clamp(Float64(y), Float64(c) + nudge * sy, Float64(d) - nudge * sy)
        # On the four edges use the exact Dirichlet datum (the FEM solution
        # equals `g` there by construction — and avoids the cell search
        # landing on a boundary face).
        if i == 1 || i == Nx + 1 || j == 1 || j == Ny + 1
            phi[i, j] = T(g(Point(Float64(x), Float64(y))))
        else
            phi[i, j] = T(uh(Point(xe, ye)))
        end
    end

    return Laplace2DSolution{T}(a, b, c, d, xs, ys, phi)
end

end # module PadeTaylorGridapExt
