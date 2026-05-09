"""
    PadeTaylor.RobustPade

Robust Padé approximation via SVD, port of Gonnet–Güttel–Trefethen 2013
Algorithm 2 with Chebfun's QR-reweighting.

## Provenance

This module is a faithful port of `external/chebfun/padeapprox.m`
(commit `7574c77`, lines 60–143).  Three lines in `padeapprox.m`
(111–117 — the QR reweighting) go beyond GGT 2013 Algorithm 2 itself
and are part of Chebfun's reference implementation; they are
preserved here for the same reason they exist there: blocks
corresponding to approximation accuracies close to `tol` come out more
often exactly square (see `references/markdown/GGT2013_*.md:236–241`).

## Algorithm in brief

Given Taylor coefficients `c = [c_0, c_1, …]` of an analytic `f` and
target type `(m, n)`:

  1. **Special case.**  If `|c_0|, …, |c_m|` are all below a threshold,
     return `r ≡ 0` (`μ = -∞`, `ν = 0`).
  2. **Diagonal hopping.**  Build the lower-triangular Toeplitz `Z`
     of size `(m+n+1) × (n+1)` from `c`; let `C = Z[m+2:m+n+1, :]`
     (an `n × (n+1)` matrix).  Compute SVD of `C` and let `ρ` be the
     number of singular values exceeding `tol · ‖c‖₂`.  If `ρ < n`,
     replace `(m, n)` with `(m - (n - ρ), ρ)` and restart.
  3. **Null vector → denominator.**  Once `ρ = n`, the right singular
     vector for the (n+1)-th column of full-V gives `b` (denominator
     coefficients).  The QR-reweighting trick (Chebfun-only):
     refine `b` via QR factorisation of `(C·diag(|b| + √ε))ᵀ` for
     better preservation of exact zeros in the block-degenerate cases.
  4. **Numerator from upper block.**  `a = Z[1:m+1, 1:n+1] · b`.
  5. **Trim.**  Cancel any common `z^λ` factor (leading-zero `b` and
     `a`).  Discard trailing near-zero `b`.  Discard trailing
     small-relative-to-`tol·‖c‖₂` `a`.
  6. **Normalise.**  Divide both by `b[1]` so the denominator is
     monic-at-zero (`b[1] = 1`).

The output `PadeApproximant{T}` carries `a`, `b`, `μ = length(a) - 1`,
`ν = length(b) - 1`.

## Determinism contract

Default symbolic tier — bit-identical cross-platform forever for
`T <: Rational` substrate (we don't ship one in v1; mentioned for
forward-compat).  For `T <: AbstractFloat` (our v1 element types) the
behaviour is `numerical: true` — bit-identical given the platform
fingerprint, modulo the floating-point non-uniqueness in `qr` and
`svd` sign conventions.  Tests admit the latter via 1e-12 relative
match against the Octave-captured oracle.

## References

  - GGT 2013 §2 + Algorithm 2 — `references/markdown/GGT2013_*.md:31–235`.
  - `external/chebfun/padeapprox.m` — the canonical reference impl.
  - ADR-0001 (own-repo) — four-layer architecture rationale.
  - ADR-0002 (own-repo) — bigfloat-SVD via GenericLinearAlgebra.
"""
module RobustPade

using LinearAlgebra: norm, qr, Diagonal, adjoint
using ..LinAlg:      pade_svd

export PadeApproximant, robust_pade

# -----------------------------------------------------------------------------
# Type
# -----------------------------------------------------------------------------

"""
    PadeApproximant{T}

Output of `robust_pade`.  Fields:

  - `a::Vector{T}` — numerator coefficients, low-to-high (`a[1] = a_0`,
    polynomial is `a[1] + a[2] z + a[3] z² + …`).
  - `b::Vector{T}` — denominator coefficients, low-to-high; `b[1] = 1`
    by convention (output is normalised).
  - `μ::Int`       — exact numerator degree (= `length(a) - 1`).
  - `ν::Int`       — exact denominator degree (= `length(b) - 1`).
  - `γ::T`         — rescaling parameter (= `1` in v1; reserved for
    Fornberg-1981-style auto-rescaling in v1.5).

The convention `r(z) = (Σ a[k+1] z^k) / (Σ b[k+1] z^k)` matches GGT
2013 §2.

When the special case `r ≡ 0` is detected (`|c_0|, …, |c_m|` all below
threshold), we return `a = [0]`, `b = [1]`, `μ = typemin(Int)`,
`ν = 0` — `μ = -∞` per GGT 2013 convention, encoded as `typemin(Int)`.
"""
struct PadeApproximant{T}
    a::Vector{T}
    b::Vector{T}
    μ::Int
    ν::Int
    γ::T
end

PadeApproximant{T}(a, b, μ, ν) where {T} =
    PadeApproximant{T}(Vector{T}(a), Vector{T}(b), μ, ν, one(T))

# -----------------------------------------------------------------------------
# Default tolerances
# -----------------------------------------------------------------------------

"""
    default_tol(::Type{T}) -> Real

Default singular-value-threshold tolerance for `T`.  At Float64 we
match Chebfun's `padeapprox.m` default of `1e-14`.  At higher precision
we use `2.0^(-precision(T) + 10)` — 10 bits of slack above the
working-precision floor (per ADR-0002 / RESEARCH.md §2.1.3).
"""
default_tol(::Type{Float64})    = 1.0e-14
default_tol(::Type{Float32})    = 1.0f-6
default_tol(::Type{T}) where {T <: AbstractFloat} =
    T(2)^(-precision(T) + 10)
default_tol(::Type{Complex{T}}) where {T <: AbstractFloat} = default_tol(T)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

"""
    _lower_tri_toeplitz(c::AbstractVector{T}, M::Int, ncol::Int) where T

Lower-triangular Toeplitz `Z` of size `M × ncol` with `Z[i, j] = c[i - j + 1]`
for `i ≥ j`, else `0`.  GGT 2013 eq. (2.6).

(Pre-allocates the result array; we don't inherit Chebfun's MATLAB
`toeplitz(col, row)` — direct construction is simpler than dropping
`ToeplitzMatrices.jl` as a dep.)
"""
function _lower_tri_toeplitz(c::AbstractVector{T}, M::Int, ncol::Int) where {T}
    Z = zeros(T, M, ncol)
    for j = 1:ncol, i = j:M
        Z[i, j] = c[i - j + 1]
    end
    return Z
end

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

"""
    robust_pade(c::AbstractVector{T}, m::Int, n::Int; tol = default_tol(T))
        -> PadeApproximant{T}

Robust Padé approximant of type `(m, n)` to the Taylor series with
coefficients `c[1] = c₀, c[2] = c₁, …, c[m+n+1] = c_{m+n}`.

If `length(c) < m + n + 1`, the input is zero-padded to length
`m + n + 1` (matching `padeapprox.m` line 62).

Throws `ArgumentError` on negative `m` or `n`.

This is a port of Chebfun's `padeapprox.m` (commit 7574c77, lines
60–143) corresponding to GGT 2013 Algorithm 2 plus the QR-reweighting
trick at lines 111–117 (which goes beyond Algorithm 2 as published).
"""
function robust_pade(c::AbstractVector{T}, m::Integer, n::Integer;
                     tol = default_tol(T)) where {T}
    m ≥ 0 || throw(ArgumentError("robust_pade: m must be non-negative; got m=$m"))
    n ≥ 0 || throw(ArgumentError("robust_pade: n must be non-negative; got n=$n"))

    # Pad / truncate `c` to exactly `m + n + 1` entries (padeapprox.m:62).
    cv = Vector{T}(undef, m + n + 1)
    for k = 1:(m + n + 1)
        cv[k] = k ≤ length(c) ? T(c[k]) : zero(T)
    end

    # Absolute tolerance for SV thresholding (padeapprox.m:66).
    cnorm = norm(cv)
    ts = tol * cnorm

    # Tolerance type: tol may be a Real, but we need to compare against
    # T-typed magnitudes inside the loop.  Coerce `tol` and `ts` to the
    # real-typed scalar that matches T's magnitude.
    tol_typed = real(T) <: AbstractFloat ? real(T)(tol) : tol
    ts_typed  = real(T) <: AbstractFloat ? real(T)(ts)  : ts

    # Special case `r ≡ 0`: all of `c[1:m+1]` below `tol · ‖c‖∞`
    # (padeapprox.m:69).  Use ‖.‖∞ here, matching the MATLAB code; the
    # threshold is on max-magnitude, not 2-norm-scaled.
    cinf = maximum(abs, cv)
    if maximum(abs, @view(cv[1:(m + 1)])) ≤ tol_typed * cinf
        return PadeApproximant{T}([zero(T)], [one(T)], typemin(Int), 0, one(T))
    end

    # ---------------------------------------------------------------------
    # Diagonal hopping (padeapprox.m:80-102).  Build Toeplitz Z, extract
    # bottom n rows as C, SVD, count above-threshold SVs.
    # ---------------------------------------------------------------------
    local Z, C
    while true
        if n == 0
            # Trivial case: no denominator beyond `b[1] = 1`.  a = c[1:m+1].
            a = Vector{T}(cv[1:(m + 1)])
            return _trim_and_normalise(a, [one(T)], cv, ts_typed, tol_typed)
        end

        Z = _lower_tri_toeplitz(cv, m + n + 1, n + 1)
        C = Z[(m + 2):(m + n + 1), :]               # n × (n+1)

        # Rank check via thin SVD (no need for full V here).
        _, S, _ = pade_svd(C; full = false)
        ρ = count(s -> s > ts_typed, S)

        if ρ == n
            break
        end
        m -= n - ρ
        n  = ρ
    end

    # ---------------------------------------------------------------------
    # n > 0 here.  Compute b from the null right singular vector of C
    # (padeapprox.m:106-117, including the QR-reweighting refinement).
    # ---------------------------------------------------------------------

    # Full SVD so Vt is `(n+1) × (n+1)` and we can read the null vector
    # off `Vt[end, :]`.
    _, _, Vt = pade_svd(C; full = true)
    b_init = Vector{T}(Vt[end, :])

    # QR-reweighting (padeapprox.m:111-117 — beyond GGT 2013 Algorithm 2).
    # The reweighting `D = diag(|b| + √ε)` makes the null vector returned
    # by QR more often exactly hit the genuine zeros of `b` for blocks at
    # accuracy near `tol`.  See RESEARCH.md §2.2.
    eps_T = real(T) <: AbstractFloat ? sqrt(eps(real(T))) : sqrt(eps(Float64))
    D = Diagonal([abs(bk) + eps_T for bk in b_init])
    F = qr(adjoint(C * D))                        # (C*D)' in MATLAB
    # F.Q is `(n+1) × (n+1)` (full Householder Q); the (n+1)-th column
    # is in the null space of `(C*D)`.
    b = D * F.Q[:, n + 1]
    b ./= norm(b)

    # Numerator (padeapprox.m:120).
    a = Z[1:(m + 1), 1:(n + 1)] * b

    return _trim_and_normalise(Vector{T}(a), Vector{T}(b), cv, ts_typed, tol_typed)
end

# -----------------------------------------------------------------------------
# Trim leading/trailing zeros and normalise (padeapprox.m:122-138, 134-138)
# -----------------------------------------------------------------------------

function _trim_and_normalise(a::Vector{T}, b::Vector{T}, cv::Vector{T},
                              ts_typed, tol_typed) where {T}
    # Cancel common `z^λ` factor: count leading near-zeros in b
    # (padeapprox.m:122-127).  Note the MATLAB code uses `tol`, not `ts`.
    lam_idx = findfirst(x -> abs(x) > tol_typed, b)
    if lam_idx === nothing
        # Pathological: all of b is below tol.  This means the numerical
        # rank is zero — the special r ≡ 0 case should have caught it.
        # Fail loudly; this signals a bug upstream.
        error("RobustPade._trim_and_normalise: all denominator coefficients below tol; " *
              "should have been caught by the r ≡ 0 special case. " *
              "Suggestion: file a bead with the input c, m, n.")
    end
    lam = lam_idx - 1
    if lam > 0
        b = b[(lam + 1):end]
        a = a[(lam + 1):end]
    end

    # Trailing near-zeros of b (padeapprox.m:130).
    last_b = findlast(x -> abs(x) > tol_typed, b)
    last_b === nothing && error("RobustPade: post-trim b is empty; algorithm bug.")
    b = b[1:last_b]

    # Trailing near-zeros of a (padeapprox.m:134).  Note `ts`, not `tol`.
    last_a = findlast(x -> abs(x) > ts_typed, a)
    if last_a === nothing
        # All of a below ts — degenerate; force a = [0] for safety.
        a = T[zero(T)]
    else
        a = a[1:last_a]
    end

    # Normalise to b[1] = 1 (padeapprox.m:137-138).
    b1 = b[1]
    a ./= b1
    b ./= b1

    return PadeApproximant{T}(a, b, length(a) - 1, length(b) - 1, one(T))
end

end # module RobustPade
