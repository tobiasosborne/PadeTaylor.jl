"""
    PadeTaylor.Laplace2D

In-house **2D-Chebyshev spectral Laplace solver** on a rectangle.

`laplace2d_solve` solves the Dirichlet problem

    ∇²φ = φ_xx + φ_yy = 0    on    [a,b] × [c,d],
    φ = (given data)         on the four edges,

by a tensor-product Chebyshev spectral collocation. It is the in-house
*voter (2)* of the triple-method majority-vote harmonic-extension fill
for the `P_I⁽²⁾` tritronquée surface figure (ADR-0024): `Re u` and
`Im u` of an analytic tritronquée are harmonic, so on the pole-free
sector each is the unique solution of a Dirichlet Laplace problem with
its own boundary values. Voter (1) is the ray-fan `vector_bvp_solve`;
voter (3) is the Gridap.jl FEM extension (`PadeTaylorGridapExt`,
bead `padetaylor-0ln.36`).

## Algorithm — tensor-product Chebyshev Laplacian (Trefethen SMIM Prog. 16)

  1. **Grids.** Place `Nx+1` Chebyshev extrema in `x` and `Ny+1` in `y`
     (DMSUITE convention: `t₀ = +1`, `t_N = -1`), affine-mapped to
     `[a,b]` and `[c,d]`. The boundary nodes are the first and last of
     each axis (`{1, Nx+1}`, `{1, Ny+1}`); the rest are interior.

  2. **Differentiation matrices.** `D₁` is the standard Chebyshev
     first-derivative matrix (Trefethen `cheb.m` / Weideman–Reddy
     `chebdif.m`); the builder `_cheb_D1` is copied **verbatim** from
     `BVP._chebyshev_D1` (we do not reach into another module's private
     symbol — same discipline as `VectorBVP.jl`, ADR-0023). `D₂ = D₁·D₁`,
     scaled by the affine chain-rule factor `(2/(b-a))²` resp.
     `(2/(d-c))²` — squared because `D₂` is a *second* derivative.

  3. **Tensor-product Laplacian on interior nodes.** With the unknowns
     in lexicographic order — `φ` flattened as a length-`(Nx-1)(Ny-1)`
     vector, `x` the fast index — the interior discrete Laplacian is the
     Kronecker sum

         L = kron(D²x_ii, I_y) + kron(I_x, D²y_ii)

     where `D²x_ii = D²x[int,int]`, `I_y = I_{Ny-1}`, etc. (Trefethen
     SMIM Program 16 / Program 23,
     `references/markdown/TrefethenSMIM_2000_book/…md:753–760`).

  4. **Dirichlet BC absorbed into the RHS.** `∇²φ = 0` on every interior
     node reads `Lᵢᵢ·φ_int + (boundary terms) = 0`. The boundary terms
     are the `D²` columns that hit boundary nodes acting on the known
     edge data; moved to the right-hand side they form `rhs`. This is
     the 2D analogue of `BVP.jl`'s `bc_col` boundary-contribution column.
     Concretely: for each axis, `D²x` splits into `D²x[int,int]` and
     `D²x[int,bnd]`; the contribution of the left/right edges enters
     through `D²x[int,bnd]`, of the bottom/top edges through `D²y[int,bnd]`.

  5. **One linear solve.** `φ_int = L \\ (-rhs)`. The problem is **linear**
     — a single solve, no Newton iteration (contrast `BVP.bvp_solve`).
     Generic in `T`: Float64 routes `\\` through LAPACK, BigFloat through
     `GenericLinearAlgebra` (the `BVP.jl` / `VectorBVP.jl` idiom).

  6. **Dense output.** The solution object carries `φ` on the full
     tensor grid (boundary + interior) and is callable, `sol(x,y) -> φ`,
     by a tensor product of the Berrut–Trefethen 2004 Chebyshev-2
     barycentric formula: barycentric-interpolate each grid row in `x`,
     then interpolate the resulting column in `y`.

## Scope

`Laplace2D` solves Laplace on a **plain rectangle** and is fully
domain-agnostic. The conformal `w = log x` map that sends an angular
sector to a rectangle, and the `P_I⁽²⁾` tritronquée boundary data, are
**not** this module's concern — they belong to F3 (bead
`padetaylor-0ln.33`). Keeping `Laplace2D` general means the same solver
serves any rectangle Dirichlet Laplace problem.

## Fail-fast contract (CLAUDE.md Rule 1)

Throws `ArgumentError` with a suggestion on:
  - `Nx < 3` or `Ny < 3` (need ≥ 2 interior nodes per axis for a
    meaningful tensor-product Laplacian).
  - A boundary-data vector whose length does not match the corresponding
    Chebyshev edge node count.

## References

  - `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md`
    — SMIM ch. 6 (`D₁`/`D₂`), Program 16 / 23 (`:753–760`).
  - `src/BVP.jl:525–547` — `_chebyshev_D1`, the verbatim source of the
    `D₁` builder; `:279–297` — the `D₂ = D₁²` + affine-scale + `bc_col`
    boundary-contribution idiom.
  - `references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/…` —
    the Chebyshev-2 barycentric formula.
  - `docs/adr/0024-laplace-harmonic-extension.md` — the triple-method
    majority-vote decision this module is voter (2) of.
"""
module Laplace2D

using LinearAlgebra: kron, I
using GenericLinearAlgebra: GenericLinearAlgebra  # extends `\` for BigFloat

export laplace2d_solve, Laplace2DSolution, laplace2d_solve_gridap

# =============================================================================
# Solution container
# =============================================================================

"""
    Laplace2DSolution{T <: AbstractFloat}

The Dirichlet Laplace solution on the rectangle `[a,b] × [c,d]`.

Fields:
  - `a, b, c, d`        : rectangle corners.
  - `xs::Vector{T}`     : `Nx+1` Chebyshev `x`-nodes (descending, `b … a`).
  - `ys::Vector{T}`     : `Ny+1` Chebyshev `y`-nodes (descending, `d … c`).
  - `phi::Matrix{T}`    : `φ` on the full tensor grid, `phi[i,j] = φ(xs[i], ys[j])`.

Call as `sol(x, y) -> φ` for tensor-barycentric evaluation anywhere in
the rectangle.
"""
struct Laplace2DSolution{T <: AbstractFloat}
    a   :: T
    b   :: T
    c   :: T
    d   :: T
    xs  :: Vector{T}
    ys  :: Vector{T}
    phi :: Matrix{T}
end

# =============================================================================
# Public driver
# =============================================================================

"""
    laplace2d_solve(bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi,
                    (a, b), (c, d); Nx = 24, Ny = 24) -> Laplace2DSolution

Solve `∇²φ = 0` on the rectangle `[a,b] × [c,d]` with Dirichlet data on
the four edges.

The four `bc_*` arguments give the boundary data; each may be a
**callable** `t -> φ` evaluated at the Chebyshev edge nodes, or a
**vector** already sampled at those nodes:

  - `bc_x_lo(y)` / `bc_x_hi(y)` : data on the `x = a` / `x = b` edges,
    a function of `y` (vectors are length `Ny+1`, ordered like `ys`).
  - `bc_y_lo(x)` / `bc_y_hi(x)` : data on the `y = c` / `y = d` edges,
    a function of `x` (vectors are length `Nx+1`, ordered like `xs`).

Kwargs `Nx`, `Ny` are the per-axis Chebyshev subdivision counts; the
grids have `Nx+1` and `Ny+1` nodes. Both must be `≥ 3`.

Edge data must be consistent at the four corners (the boundary is
continuous); `laplace2d_solve` does not enforce this — supplying the
exact trace of one harmonic function automatically is consistent.

Throws `ArgumentError` on `Nx < 3` / `Ny < 3` or boundary-vector
length mismatch.
"""
function laplace2d_solve(bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi,
                         xrange, yrange; Nx::Integer = 24, Ny::Integer = 24)

    Nx ≥ 3 || throw(ArgumentError(
        "laplace2d_solve: Nx must be ≥ 3 (got $Nx).  Suggestion: use " *
        "≥ 16 for spectral accuracy on a non-trivial harmonic field — " *
        "the tensor-product Laplacian needs ≥ 2 interior nodes per axis."))
    Ny ≥ 3 || throw(ArgumentError(
        "laplace2d_solve: Ny must be ≥ 3 (got $Ny).  Suggestion: use ≥ 16."))

    a, b = xrange
    c, d = yrange
    T = float(promote_type(typeof(a), typeof(b), typeof(c), typeof(d)))
    a, b, c, d = T(a), T(b), T(c), T(d)

    # ---- Chebyshev grids (DMSUITE convention: node 1 = +1 ↦ b, last = -1 ↦ a)
    πT = T(π)
    tx = T[cos(T(j) * πT / T(Nx)) for j in 0:Nx]
    ty = T[cos(T(j) * πT / T(Ny)) for j in 0:Ny]
    xs = T[(a + b) / 2 + (b - a) / 2 * t for t in tx]      # length Nx+1
    ys = T[(c + d) / 2 + (d - c) / 2 * t for t in ty]      # length Ny+1

    # ---- D₂ on each axis, affine-scaled.  Squared chain-rule factor. --------
    D1x = _cheb_D1(tx, T, Nx)
    D1y = _cheb_D1(ty, T, Ny)
    D2x = (D1x * D1x) .* (T(2) / (b - a))^2
    D2y = (D1y * D1y) .* (T(2) / (d - c))^2

    # ---- interior / boundary index split -----------------------------------
    ix, iy = 2:Nx, 2:Ny                  # interior node indices
    bx, by = [1, Nx + 1], [1, Ny + 1]    # boundary node indices (hi, lo)
    nix, niy = Nx - 1, Ny - 1

    # ---- assemble the boundary trace on the full grid -----------------------
    # φ is known on all four edges; build the full Matrix{T} grid with the
    # interior left zero (it is overwritten by the solve).
    phi = zeros(T, Nx + 1, Ny + 1)
    xlo = _sample(bc_x_lo, ys, niy + 2)   # x = a edge, varies over y
    xhi = _sample(bc_x_hi, ys, niy + 2)   # x = b edge
    ylo = _sample(bc_y_lo, xs, nix + 2)   # y = c edge, varies over x
    yhi = _sample(bc_y_hi, xs, nix + 2)   # y = d edge
    # xs node 1 = b (hi), node Nx+1 = a (lo); likewise ys.
    phi[Nx + 1, :] .= xlo
    phi[1, :]      .= xhi
    phi[:, Ny + 1] .= ylo
    phi[:, 1]      .= yhi

    # ---- tensor-product interior Laplacian L (Trefethen SMIM Prog. 16) ------
    # Unknown vector is φ[ix, iy] flattened column-major (x fast), matching
    # Julia `vec`.  L = kron(D²y_ii, I_x) + kron(I_y, D²x_ii) acts on it.
    D2x_ii = D2x[ix, ix]
    D2y_ii = D2y[iy, iy]
    L = kron(Matrix{T}(I, niy, niy), D2x_ii) .+
        kron(D2y_ii, Matrix{T}(I, nix, nix))

    # ---- boundary contribution → RHS ---------------------------------------
    # ∇²φ = 0 ⇒ Lᵢᵢ φ_int = −(D² columns hitting boundary nodes · φ_bnd).
    # x-edges: D2x[ix, bx] · φ[bx, iy] feeds each y-row.
    # y-edges: D2y[iy, by] · φ[ix, by]ᵀ feeds each x-column.
    rhs = zeros(T, nix, niy)
    D2x_ib = D2x[ix, bx]                 # nix × 2
    D2y_ib = D2y[iy, by]                 # niy × 2
    @inbounds for (col, j) in enumerate(iy)
        rhs[:, col] .+= D2x_ib * phi[bx, j]
    end
    @inbounds for (row, i) in enumerate(ix)
        rhs[row, :] .+= D2y_ib * phi[i, by]
    end

    # ---- one linear solve (LAPACK for Float64, GenericLinearAlgebra for BigFloat)
    φ_int = L \ (-vec(rhs))
    phi[ix, iy] .= reshape(φ_int, nix, niy)

    return Laplace2DSolution{T}(a, b, c, d, xs, ys, phi)
end

"""
    _sample(bc, nodes, expected_len) -> Vector

Resolve a boundary-data argument: a callable is evaluated at each node;
a vector is length-checked and returned (promoted to the node eltype).
Fail-fast on a length mismatch (CLAUDE.md Rule 1).
"""
function _sample(bc, nodes::AbstractVector{T}, expected_len::Integer) where T
    if bc isa AbstractVector
        length(bc) == expected_len || throw(ArgumentError(
            "laplace2d_solve: boundary-data vector has length $(length(bc)) " *
            "but the matching Chebyshev edge has $expected_len nodes.  " *
            "Suggestion: pass a callable `t -> φ`, or sample the data at " *
            "the edge's Chebyshev nodes (descending order, +endpoint first)."))
        return T[T(v) for v in bc]
    end
    return T[T(bc(t)) for t in nodes]
end

# =============================================================================
# Callable: tensor-barycentric evaluation at (x, y)
# =============================================================================

"""
    (sol::Laplace2DSolution)(x, y) -> φ

Evaluate the Laplace solution at an arbitrary `(x, y)` in the rectangle
by a tensor product of the Berrut–Trefethen 2004 Chebyshev-2 barycentric
formula: interpolate each grid row of `phi` in `x` to get a column in
`y`, then interpolate that column in `y`.
"""
function (sol::Laplace2DSolution{T})(x, y) where T
    xT, yT = T(x), T(y)
    Nyp1 = length(sol.ys)
    # Interpolate phi[:, j] in x for every y-node j → a length-(Ny+1) column.
    col = T[_bary(sol.xs, view(sol.phi, :, j), xT) for j in 1:Nyp1]
    return _bary(sol.ys, col, yT)
end

# =============================================================================
# Helpers (un-exported; module-private)
# =============================================================================

"""
    _cheb_D1(t, ::Type{T}, N) -> Matrix{T}

`(N+1)×(N+1)` Chebyshev first-derivative matrix in `t ∈ [-1,1]`.

Copied **verbatim** from `BVP._chebyshev_D1` (CLAUDE.md discipline: do
not reach into another module's private symbol — duplicate the ~20-line
builder, as ADR-0023's `VectorBVP.jl` does). Off-diagonal entries
`D[i,j] = (cᵢ/cⱼ)(-1)^{i+j}/(tᵢ-tⱼ)` with `c₀ = c_N = 2`, interior
`c = 1`; the diagonal is fixed by the negative-row-sum identity so that
`D·1 = 0`.
"""
function _cheb_D1(t::AbstractVector{T}, ::Type{T}, N::Integer) where T <: AbstractFloat
    np1 = N + 1
    D   = Matrix{T}(undef, np1, np1)
    @inbounds for j in 1:np1, i in 1:np1
        if i == j
            D[i, j] = zero(T)
        else
            ci   = (i == 1 || i == np1) ? T(2) : T(1)
            cj   = (j == 1 || j == np1) ? T(2) : T(1)
            sign = iseven(i + j) ? T(1) : T(-1)
            D[i, j] = (ci / cj) * sign / (t[i] - t[j])
        end
    end
    @inbounds for i in 1:np1
        s = zero(T)
        for j in 1:np1
            i == j && continue
            s += D[i, j]
        end
        D[i, i] = -s
    end
    return D
end

"""
    _bary(x_nodes, vals, x_star) -> value

Berrut–Trefethen 2004 second-form barycentric interpolation on Chebyshev-2
nodes (mapped to `[a,b]` — the affine map preserves the barycentric
weights `wⱼ = (-1)ʲ`, halved at the endpoints). Coincident-node guard
returns the node value directly.
"""
function _bary(x_nodes::AbstractVector{T}, vals::AbstractVector, x_star::T) where T
    n    = length(x_nodes)
    atol = eps(T)
    @inbounds for j in 1:n
        abs(x_star - x_nodes[j]) ≤ atol && return T(vals[j])
    end
    num = zero(T)
    den = zero(T)
    @inbounds for j in 1:n
        w     = (iseven(j - 1) ? T(1) : T(-1)) *
                ((j == 1 || j == n) ? T(0.5) : T(1))
        δ     = x_star - x_nodes[j]
        num  += w * T(vals[j]) / δ
        den  += w / δ
    end
    return num / den
end

# =============================================================================
# Voter (3): Gridap FEM Laplace solver — package-extension stub
# =============================================================================

"""
    laplace2d_solve_gridap(bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi,
                           (a, b), (c, d); n_x = 32, n_y = 32) -> Laplace2DSolution

Solve the **same** Dirichlet Laplace problem `∇²φ = 0` on the rectangle
`[a,b] × [c,d]` as [`laplace2d_solve`](@ref), but by **Gridap.jl finite
elements** instead of Chebyshev spectral collocation.

This is *voter (3)* of the triple-method majority-vote harmonic-extension
fill (ADR-0024): an algorithmically disjoint third solver, so F3 can
majority-vote `ray-fan BVP` vs `2D-Chebyshev` vs `Gridap`. A FEM bug and
a spectral bug will not coincide.

The method is provided by the `PadeTaylorGridapExt` package extension
(ADR-0003 weak-dep pattern) — it exists only when `Gridap` is loaded
alongside `PadeTaylor`. Calling it **without** a `Gridap` load throws an
`ArgumentError` directing the caller to `using Gridap` (CLAUDE.md
Rule 1, fail loud — no silent `nothing`).

The four `bc_*` arguments and the rectangle `(a,b)`, `(c,d)` mirror
`laplace2d_solve`'s contract; the kwargs `n_x` / `n_y` are the FEM mesh
subdivision counts (the `CartesianDiscreteModel` has `n_x × n_y` cells).
The return value's interface (a `Laplace2DSolution`-shaped object with a
`φ` grid and a callable `(x,y) -> φ`) is compatible with
`laplace2d_solve`'s so F3 treats the two Laplace voters uniformly.
"""
function laplace2d_solve_gridap(args...; kwargs...)
    throw(ArgumentError(
        "laplace2d_solve_gridap: the Gridap FEM Laplace backend (voter (3) " *
        "of the ADR-0024 triple-method tritronquée sector fill) is a weak " *
        "dependency and is not currently loaded.  Suggestion: `using Gridap` " *
        "alongside `using PadeTaylor` — that activates the `PadeTaylorGridapExt` " *
        "package extension which provides this method (ADR-0003)."))
end

end # module Laplace2D
