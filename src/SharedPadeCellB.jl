"""
    PadeTaylor.SharedPadeCellB

The **square Mano–Tsuda cell** ("cell B") of the ADR-0028 dual-construction
shared-denominator Padé dispatch — the second of the two Padé-table cells the
vector stepper builds per step and selects between.

## Why a second cell exists

The shipped `SharedPade.shared_denominator_pade` is the **GGT diagonal `(m,m)`**
cell ("cell A"): a `(d·m)×(m+1)` block-Toeplitz null-space system. On a *genuine
shared pole* its over-determined stack has a sharp isolated null space and
recovers the pole crisply. But on *entire / locally-regular* jets — the smooth
stretches between poles — that same over-determination is rank-deficient, so the
smallest right-singular vector is only a *least-squares* denominator that places
spurious poles, and the accuracy cost **grows with the component count `d`**
(ADR-0028 §Context: harmonic ~5e-9, Noumi–Yamada A₄ ~6.4e-6).

Cell B cures this on the regular stretches. It is the **genuinely wide-square**
simultaneous (type-II Hermite–Padé) system of Mano–Tsuda 2017 (eq:lspa,
`references/tex/hermite_pade/ManoTsuda2017_*/hp_arXiv_final.tex:1405-1428`),
realised at degree

    m_eff = ⌈m/d⌉ · d      (≥ m),

with `n = ⌈m/d⌉` rows per block, the **`+1` window** (block top-left `c_{m_eff}`,
one lower than cell A's `c_{m_eff+1}`), and degree-`(m_eff−1)` numerators. Because
`d·n = m_eff` exactly, the stacked matrix is `m_eff × (m_eff+1)` — *wide*, so it
always has a clean ≥1-dimensional null space (no least-squares over-fit), which is
why the audition (probe `external/probes/adr0028-dual-construction-audition/`,
`cell_B_grow`) found it recovers **6–8 orders** on entire/regular/high-`d` steps.

The degree is rounded **up** (`⌈·⌉`, not `⌊·⌋`) so the cell keeps as much degree as
possible: the audition showed the higher-degree wide-square variant matches the
lower-degree one on entire jets *and* closes most of the gap to cell A on
meromorphic steps (ADR-0028 Amendment 1.2).

## Robustness — fall back to cell A, never lie (Rule 1)

Cell B deliberately omits the ADR-0027 ρ-reduction loop (the wide-square system
does not rank-collapse the way the tall cell-A stack does). Instead, every
degenerate input returns `nothing`, a sentinel that tells the dispatcher to use
**cell A** — the proven, oracle-validated, meromorphic-safe cell:

  1. `d > m` — `m_eff = d > m` with one row per block; ill-posed.
  2. jet too short for `m_eff + 1` coefficients.
  3. the denominator collapses below `tol` after the QR step.

The shared `z^λ` cancellation and trailing-trim mirror
`SharedPade.shared_denominator_pade` (`SharedPade.jl:270-311`) and the scalar
`RobustPade._trim_and_normalise`.

## References

  - `docs/adr/0028-sharedpade-dual-construction-pareto-dispatch.md` (Amendments 1 & 2).
  - `references/tex/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285/hp_arXiv_final.tex`
    eq:lspa :1405-1428 (the square `m×(m+1)` simultaneous system), :1381 (`m=n(L−1)`).
  - `src/SharedPade.jl` — cell A; the QR-reweighting / `z^λ` / normalise pattern reused here.
  - `external/probes/adr0028-dual-construction-audition/cells.jl` — `cell_B_grow`, the
    validated prototype this module ports.
"""
module SharedPadeCellB

using LinearAlgebra: norm, qr, Diagonal, adjoint
using ..LinAlg:     pade_svd
using ..RobustPade: default_tol

export build_square_cell

# `n × (m+1)` Toeplitz block for one component, `+1` window: block[r,c] =
# c_{m+r-c} (top-left c_m), one coefficient lower than cell A's `+2` window.
function _block_b(c::AbstractVector{T}, m::Int, n::Int) where {T}
    blk = zeros(T, n, m + 1)
    @inbounds for cc = 1:(m + 1), rr = 1:n
        idx = m + rr - cc + 1            # c[idx] = c_{m+r-c}; top-left c_m
        if 1 ≤ idx ≤ length(c)
            blk[rr, cc] = T(c[idx])
        end
    end
    return blk
end

# `nrows × (m+1)` upper Toeplitz recovering a degree-(nrows-1) numerator: a = up·b,
# up[r,c] = c_{r-c}.  For cell B nrows = m_eff ⇒ degree-(m_eff-1) numerators.
function _upper_b(c::AbstractVector{T}, nrows::Int, m::Int) where {T}
    up = zeros(T, nrows, m + 1)
    @inbounds for cc = 1:(m + 1), rr = cc:nrows
        idx = rr - cc + 1
        if 1 ≤ idx ≤ length(c)
            up[rr, cc] = T(c[idx])
        end
    end
    return up
end

"""
    build_square_cell(jets, m; tol = default_tol(T))
        -> (numerators, denominator)::Tuple | nothing

Build the wide-square Mano–Tsuda cell B at degree `m_eff = ⌈m/d⌉·d`.  Returns the
same `(numerators::Vector{Vector{T}}, denominator::Vector{T})` 2-tuple shape as
`shared_denominator_pade`, normalised to `Q(0) = 1`.

Returns `nothing` (⇒ the dispatcher falls back to cell A) on any degenerate input:
`d ≥ m`, a jet too short for `m_eff + 1` coefficients, or a denominator that
collapses below `tol` after the QR step.  Never throws and never returns a NaN/zero
lie — cell A is the safety net.
"""
function build_square_cell(jets::AbstractVector{<:AbstractVector{T}},
                           m::Integer;
                           tol::Real = default_tol(T)) where {T}
    d = length(jets)
    m = Int(m)
    # Guard 1 (RISK 4): d > m ⇒ ⌈m/d⌉=1, m_eff=d>m, one row per block; ill-posed.
    # (d == m is fine: m_eff = m, a genuine m×(m+1) wide-square system.)
    d > m && return nothing

    n = cld(m, d)
    m_eff = n * d
    # Guard 2 (RISK 4): jets must carry ≥ m_eff+1 coefficients.
    minjet = minimum(length, jets)
    m_eff + 1 > minjet && return nothing

    tol_t = real(T) <: AbstractFloat ? real(T)(tol) : tol
    cnorm = norm(reduce(vcat, [collect(j) for j in jets]))
    cnorm ≤ tol_t && return nothing             # genuinely-zero input ⇒ let cell A throw
    τ = tol_t * cnorm

    blocks = [_block_b(jets[i], m_eff, n) for i = 1:d]
    A_full = reduce(vcat, blocks)               # (d·n) × (m_eff+1) = m_eff × (m_eff+1)
    _, _, Vt = pade_svd(A_full; full = true)
    b_init = Vector{T}(Vt[end, :])

    # QR-reweighting (padeapprox.m 111–117; identical to SharedPade.jl:264-268).
    eps_T = real(T) <: AbstractFloat ? sqrt(eps(real(T))) : sqrt(eps(Float64))
    D = Diagonal([abs(bk) + eps_T for bk in b_init])
    F = qr(adjoint(A_full * D))
    b = D * F.Q[:, size(A_full, 2)]
    b ./= norm(b)

    raw_nums = [Vector{T}(_upper_b(jets[i], m_eff, m_eff) * b) for i = 1:d]

    # Leading-z^λ cancellation (SharedPade.jl:284-293).
    lam_idx = findfirst(x -> abs(x) > tol_t, b)
    lam_idx === nothing && return nothing       # Guard 3 (RISK 4): denominator below tol
    lam = lam_idx - 1
    if lam > 0
        b = b[(lam + 1):end]
        raw_nums = [a[(lam + 1):end] for a in raw_nums]
    end

    b1 = b[1]
    denominator = Vector{T}(b ./ b1)
    last_b = findlast(x -> abs(x) > tol_t, denominator)
    last_b === nothing && return nothing
    denominator = denominator[1:last_b]

    numerators = Vector{Vector{T}}(undef, d)
    for i = 1:d
        a_i = raw_nums[i] ./ b1
        last_a = findlast(x -> abs(x) > τ, a_i)
        numerators[i] = last_a === nothing ? T[zero(T)] : a_i[1:last_a]
    end
    return (numerators, denominator)
end

end # module SharedPadeCellB
