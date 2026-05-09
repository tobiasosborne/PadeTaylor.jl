"""
    PadeTaylor.LinAlg

SVD dispatcher with relative-accuracy guarantees on small singular values
for the arb-prec tier.

## Design rationale

GGT 2013 Algorithm 2 classifies a singular value as "zero" if
`σᵢ < tol · ‖c‖₂` (where `c` is the input Taylor-coefficient vector).
At Float64 the typical regime is `tol = 10⁻¹⁴` and matrix condition
numbers of `~10¹²`, so absolute-accuracy LAPACK Demmel-Kahan
(`DGESVD`) is fine — it guarantees `c · 2⁻ᵖ · σ_max` per SV, well
below the threshold. **At arb-prec the regime changes**: near a Padé
table block boundary the gap between the smallest "signal" SV and the
largest "noise" SV may span only a few orders of magnitude, far less
than the dynamic range of the working precision. We need
*relative-accuracy* guarantees on small SVs.

`GenericLinearAlgebra.svd` provides this via one-sided Jacobi
(Demmel-Veselić 1992): error `c · 2⁻ᵖ · σᵢ` per SV. Small SVs stay
reliably above the threshold regardless of `κ(A)`. Workbench's own
SVD dispatcher routes `n ≤ 500` to Jacobi for the same reason
(`scientist-workbench/tools/linalg-svd/`, ADR-0014/0015/0016).

For GGT matrices at the typical `n ≤ 60` Jacobi is trivially fast.

The `Arb`-element-type path lives in `ext/PadeTaylorArblibExt.jl`
(per ADR-0003, loaded only when `using Arblib` is active in the user's
session): convert `Matrix{Arb}` → `Matrix{BigFloat}` (mid-point
rounding, radius discarded — see ADR-0002 caveat) before dispatch.

See `RESEARCH.md §5.1` for the empirical landscape and ADR-0002 for
the full design argument.

## Why we expose `full::Bool`

GGT 2013 Algorithm 2 needs the **null right singular vector** of an
`n × (n+1)` Toeplitz matrix `C̃`. In the thin SVD (`full=false`), Vt is
`n × (n+1)` and V is `(n+1) × n` — the `(n+1)`-th column of full V is
not produced. With `full=true`, Vt is `(n+1) × (n+1)` and the null
vector is `transpose(Vt[end, :])`. RobustPade calls `pade_svd(C̃;
full=true)` for this reason. Default `full=false` keeps the rest of
the API the cheaper thin path.

## API

  - `pade_svd(A::AbstractMatrix{T}; full::Bool = false) -> (U, S, Vt)`

For an `m × n` input with `full=false` (thin SVD): `U` is
`m × min(m,n)`, `S` is `min(m,n)`-vector, `Vt` is `min(m,n) × n`.
With `full=true`: `U` is `m × m`, `S` is `min(m,n)`-vector, `Vt` is
`n × n`. Both backends honour the `full` kwarg identically.
"""
module LinAlg

using LinearAlgebra:        LinearAlgebra, svd
using GenericLinearAlgebra: GenericLinearAlgebra

# -----------------------------------------------------------------------------
# Float64 / Float32 / Complex{F64} / Complex{F32}: stdlib `svd` (LAPACK).
# -----------------------------------------------------------------------------

# Type union for the LAPACK-supported real and complex types.
const _LAPACK_FLOAT = Union{Float32, Float64}
const _LAPACK_TYPE  = Union{_LAPACK_FLOAT, Complex{Float32}, Complex{Float64}}

"""
    pade_svd(A::AbstractMatrix{T}; full::Bool = false)
        where T <: Union{Float32, Float64, Complex{Float32}, Complex{Float64}}
        -> (U::Matrix{T}, S::Vector{<:Real}, Vt::Matrix{T})

Dispatches to LAPACK Demmel-Kahan via `LinearAlgebra.svd`.
"""
function pade_svd(A::AbstractMatrix{T}; full::Bool = false) where {T <: _LAPACK_TYPE}
    F = LinearAlgebra.svd(A; full = full)
    return (F.U, F.S, F.Vt)
end

# -----------------------------------------------------------------------------
# Generic AbstractFloat (e.g. BigFloat): GenericLinearAlgebra one-sided Jacobi.
# -----------------------------------------------------------------------------

"""
    pade_svd(A::AbstractMatrix{T}; full::Bool = false)
        where T <: AbstractFloat (non-LAPACK)
        -> (U::Matrix{T}, S::Vector{T}, Vt::Matrix{T})

Dispatches to `GenericLinearAlgebra.svd` (one-sided Jacobi
Demmel-Veselić). Provides `c · 2⁻ᵖ · σᵢ` relative-accuracy guarantees
on every singular value, load-bearing for the GGT 2013 rank-counting
threshold at arb-prec.
"""
function pade_svd(A::AbstractMatrix{T}; full::Bool = false) where {T <: AbstractFloat}
    F = GenericLinearAlgebra.svd(A; full = full)
    return (F.U, F.S, F.Vt)
end

"""
    pade_svd(A::AbstractMatrix{Complex{T}}; full::Bool = false)
        where T <: AbstractFloat (non-LAPACK)
        -> (U::Matrix{Complex{T}}, S::Vector{T}, Vt::Matrix{Complex{T}})

Same as the real `AbstractFloat` path, for complex element types.
"""
function pade_svd(A::AbstractMatrix{Complex{T}}; full::Bool = false) where {T <: AbstractFloat}
    F = GenericLinearAlgebra.svd(A; full = full)
    return (F.U, F.S, F.Vt)
end

end # module LinAlg
