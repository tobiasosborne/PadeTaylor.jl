"""
    PadeTaylor.LinAlg

SVD dispatcher for LAPACK-supported and arbitrary-precision element types.

## Design rationale

GGT 2013 Algorithm 2 classifies a singular value as "zero" if
`σᵢ < tol · ‖c‖₂` (where `c` is the input Taylor-coefficient vector).
At Float64 the typical regime is `tol = 10⁻¹⁴` and matrix condition
numbers of `~10¹²`, so the usual `O(2⁻ᵖ · σ₁)` absolute singular-value
error is below the rank threshold. The arbitrary-precision default
preserves that separation: `default_tol = 2^(-p+10)` leaves ten bits
between working precision and `tol`, while a GGT Toeplitz block of the
typical `n ≤ 60` has `σ₁ ≤ ‖C̃‖F ≤ √n · ‖c‖₂`. Thus the threshold
`tol · ‖c‖₂` remains comfortably above the absolute SVD error scale.
The CRP.4 corpus cases pin the resulting BigFloat rank decisions
empirically (`test/corpus_robust_pade_test.jl:190-201`).

`GenericLinearAlgebra.svd` does not use one-sided Jacobi. It reduces
the matrix to bidiagonal form with Householder reflectors, then applies
bidiagonal QR with the Demmel-Kahan zero-shift iteration when required
(`~/.julia/packages/GenericLinearAlgebra/X90Kh/src/svd.jl:328-345,
648-654, 81-87, 235-244`).

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
vector is the last column of `V`, i.e. `conj(Vt[end, :])` — `Vt` is
the *conjugate* transpose `V'`, so for complex input the bare row
`Vt[end, :]` is **not** a null vector (test 1.1.3c); for real input
the conjugate is a no-op. RobustPade calls `pade_svd(C̃; full=true)`
for this reason (it takes `abs.()` of the row and recomputes `b` by
adjoint-QR, so it never consumes the row unconjugated). Default
`full=false` keeps the rest of the API the cheaper thin path.

## API

  - `pade_svd(A::AbstractMatrix{T}; full::Bool = false) -> (U, S, Vt)`

For an `m × n` input with `full=false` (thin SVD): `U` is
`m × min(m,n)`, `S` is `min(m,n)`-vector, `Vt` is `min(m,n) × n`.
With `full=true`: `U` is `m × m`, `S` is `min(m,n)`-vector, `Vt` is
`n × n`. Both backends honour the `full` kwarg identically.
"""
module LinAlg

using LinearAlgebra:        LinearAlgebra
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
# Generic AbstractFloat (e.g. BigFloat): GenericLinearAlgebra SVD.
# -----------------------------------------------------------------------------

"""
    pade_svd(A::AbstractMatrix{T}; full::Bool = false)
        where T <: AbstractFloat (non-LAPACK)
        -> (U::Matrix{T}, S::Vector{T}, Vt::Matrix{T})

Dispatches to `GenericLinearAlgebra.svd`: Householder bidiagonalisation
followed by bidiagonal QR, with the Demmel-Kahan zero-shift iteration
used when required. The ten-bit gap between working precision and
`default_tol` keeps its absolute error scale below the GGT 2013
rank-counting threshold at the supported matrix sizes; see the module
docstring for the bound and empirical corpus pin.
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
