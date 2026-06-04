# cells.jl — the two shared-denominator Padé "cells" of ADR-0028, plus the
# held-out-residual metric and the Pareto/lexicographic selector.
#
# This is a PROBE (external/probes/), not src/.  It exists to gather wild
# evidence for the seven maintainer sign-off decisions in
# docs/adr/0028-sharedpade-dual-construction-pareto-dispatch.md BEFORE the
# gated production build.  Nothing here ships; src/ is untouched.
#
# Ground truth:
#   - cell (A) = the shipped GGT diagonal (m,m): src/SharedPade.jl (window +2,
#     idx=m+r-c+2, top-left c_{m+1}; m rows/block; degree-m numerator;
#     matching O(z^{2m+1})).
#   - cell (B) = Mano-Tsuda square (m-1,m): hp_arXiv_final.tex eq:lspa :1405-1428
#     (window +1, idx=m+r-c+1, top-left c_m; n=m/d rows/block; degree-(m-1)
#     numerator; matching O(z^{m+n+1})).  Verified faithful by the Mano-Tsuda
#     scout; the ONLY source-silent point is the d∤m row count (n=⌈m/d⌉ vs
#     ⌊m/d⌋+shrink-m), which we audition both ways (Q7).
#
# Both cells reuse the QR-reweighting + z^λ cancellation + ρ-reduction exactly
# as SharedPade.jl, so the only deliberate differences are the three knobs
# above.  cell_A is asserted bit-equal (mod sign) to the shipped
# shared_denominator_pade in smoke.jl.

using PadeTaylor: shared_denominator_pade, robust_pade
using PadeTaylor.LinAlg: pade_svd
using PadeTaylor.RobustPade: default_tol
using LinearAlgebra: norm, qr, Diagonal, adjoint

# --- block builders ---------------------------------------------------------

# One Toeplitz block, n × (m+1).  `window`=2 reproduces cell A's c_{m+1}
# top-left (idx=m+r-c+2); `window`=1 gives cell B's c_m top-left (idx=m+r-c+1).
function _block(c::AbstractVector{T}, m::Int, n::Int, window::Int) where {T}
    blk = zeros(T, n, m + 1)
    @inbounds for cc = 1:(m + 1), rr = 1:n
        idx = m + rr - cc + window
        if 1 ≤ idx ≤ length(c)
            blk[rr, cc] = T(c[idx])
        end
    end
    return blk
end

# Numerator-recovery upper Toeplitz, `nrows × (m+1)`, up[r,c]=c_{r-c}.
# nrows=m+1 ⇒ degree-m numerator (cell A); nrows=m ⇒ degree-(m-1) (cell B).
function _upper(c::AbstractVector{T}, nrows::Int, m::Int) where {T}
    up = zeros(T, nrows, m + 1)
    @inbounds for cc = 1:(m + 1), rr = cc:nrows
        idx = rr - cc + 1
        if 1 ≤ idx ≤ length(c)
            up[rr, cc] = T(c[idx])
        end
    end
    return up
end

# --- the generic cell -------------------------------------------------------

"""
    PadeCell — the result of building one shared-denominator Padé cell.

`numerators`/`denominator` are the recovered (P_i, Q).  `S` is the singular
spectrum of the accepted A_full (for the conditioning gap g = σ_m/σ_{m+1}),
`m_cur`/`n` the accepted construction degree + rows-per-block, `K` the matching
order (first unused Taylor order is K+1), `Q0` the pre-cancellation |Q(0)|.
"""
struct PadeCell{T}
    numerators  :: Vector{Vector{T}}
    denominator :: Vector{T}
    S           :: Vector{<:Real}
    m_cur       :: Int
    n           :: Int
    K           :: Int
    Q0          :: Float64
    label       :: Symbol
end

"""
    pade_cell(jets, m; window, rows_per_block, num_offset, K_of, reduce, label, tol)

Unified block-Toeplitz null-space solver.  `window`∈{1,2} sets the matching
window; `rows_per_block(m_cur,d)` the per-block row count; `num_offset`∈{0,1}
the numerator block height (m_cur+num_offset rows ⇒ degree m_cur+num_offset-1);
`K_of(m_cur,n)` the matching order.  `reduce=true` runs SharedPade's ADR-0027
ρ-reduction loop; `false` builds once at the requested degree (strict square).
"""
function pade_cell(jets::AbstractVector{<:AbstractVector{T}}, m::Integer;
                   window::Int, rows_per_block, num_offset::Int, K_of,
                   reduce::Bool, label::Symbol,
                   tol::Real = default_tol(T)) where {T}
    d = length(jets)
    tol_t = real(T) <: AbstractFloat ? real(T)(tol) : tol
    cnorm = norm(reduce_vcat(jets))
    τ = tol_t * cnorm

    m_cur = Int(m)
    local A_full, S, Vt, n
    while true
        n = rows_per_block(m_cur, d)
        blocks = [_block(jets[i], m_cur, n, window) for i = 1:d]
        A_full = reduce_vcat_mat(blocks)
        _, S, Vt = pade_svd(A_full; full = true)
        ρ = count(s -> s > τ, S)
        (!reduce || ρ ≥ m_cur) && break
        m_cur = ρ == 0 ? m_cur - 1 : ρ
        m_cur < 1 && error("pade_cell($label): reduced below degree 1")
    end

    b_init = Vector{T}(Vt[end, :])
    Q0 = Float64(abs(b_init[1]) / max(norm(b_init), eps()))
    eps_T = real(T) <: AbstractFloat ? sqrt(eps(real(T))) : sqrt(eps(Float64))
    D = Diagonal([abs(bk) + eps_T for bk in b_init])
    F = qr(adjoint(A_full * D))
    b = D * F.Q[:, m_cur + 1]
    b ./= norm(b)

    nrows = m_cur + num_offset
    raw_nums = [Vector{T}(_upper(jets[i], nrows, m_cur) * b) for i = 1:d]

    lam_idx = findfirst(x -> abs(x) > tol_t, b)
    lam_idx === nothing && error("pade_cell($label): denominator below tol")
    lam = lam_idx - 1
    if lam > 0
        b = b[(lam + 1):end]
        raw_nums = [a[(lam + 1):end] for a in raw_nums]
    end
    b1 = b[1]
    denominator = Vector{T}(b ./ b1)
    last_b = findlast(x -> abs(x) > tol_t, denominator)
    denominator = denominator[1:last_b]
    numerators = Vector{Vector{T}}(undef, d)
    for i = 1:d
        a_i = raw_nums[i] ./ b1
        last_a = findlast(x -> abs(x) > τ, a_i)
        numerators[i] = last_a === nothing ? T[zero(T)] : a_i[1:last_a]
    end
    return PadeCell{T}(numerators, denominator, S, m_cur, n, K_of(m_cur, n), Q0, label)
end

reduce_vcat(jets) = reduce(vcat, [collect(j) for j in jets])
reduce_vcat_mat(blocks) = reduce(vcat, blocks)

# --- the two (three) named cells --------------------------------------------

cell_A(jets, m; tol = default_tol(eltype(eltype(jets)))) =
    pade_cell(jets, m; window = 2, rows_per_block = (mc, d) -> mc,
              num_offset = 1, K_of = (mc, n) -> 2mc, reduce = true,
              label = :diagonal, tol = tol)

# cell B, ⌈m/d⌉ rows (the ADR's d∤m round-up; over-determined when d∤m).
cell_B_ceil(jets, m; tol = default_tol(eltype(eltype(jets)))) =
    pade_cell(jets, m; window = 1, rows_per_block = (mc, d) -> cld(mc, d),
              num_offset = 0, K_of = (mc, n) -> mc + n, reduce = true,
              label = :square_ceil, tol = tol)

# cell B, ⌊m/d⌋ rows + shrink m to n·d (the genuinely-square Mano-Tsuda system,
# no reduction loop).  When d∤m this LOWERS the denominator degree.
function cell_B_floor(jets, m; tol = default_tol(eltype(eltype(jets))))
    d = length(jets)
    n = max(fld(Int(m), d), 1)
    m_eff = n * d
    pade_cell(jets, m_eff; window = 1, rows_per_block = (mc, dd) -> n,
              num_offset = 0, K_of = (mc, nn) -> mc + n, reduce = false,
              label = :square_floor, tol = tol)
end

# cell B, ⌈m/d⌉ rows + GROW m to n·d (≥ m): a genuinely-square (wide) system
# like B_floor but at HIGHER degree — clean null AND high order.  The candidate
# best-of-both for the extended-precision genuine-pole regime where B_floor's
# reduced degree stalls.
function cell_B_grow(jets, m; tol = default_tol(eltype(eltype(jets))))
    d = length(jets)
    n = cld(Int(m), d)
    m_eff = n * d
    pade_cell(jets, m_eff; window = 1, rows_per_block = (mc, dd) -> n,
              num_offset = 0, K_of = (mc, nn) -> mc + n, reduce = false,
              label = :square_grow, tol = tol)
end

# --- held-out cross-validation residual (ADR-0028 §2) -----------------------

# coeff_k of (P_i - f_i·Q) in the rescaled (t) frame.  P_i,Q polynomials,
# f_i the jet.  For k>deg(P_i) the P_i term is 0.
function _series_coeff(p::AbstractVector{T}, q::AbstractVector{T},
                       c::AbstractVector{T}, k::Int) where {T}
    pk = (k + 1) ≤ length(p) ? p[k + 1] : zero(T)
    s = zero(T)
    @inbounds for j = 0:(length(q) - 1)
        ci = k - j
        if 0 ≤ ci ≤ length(c) - 1
            s += q[j + 1] * c[ci + 1]
        end
    end
    return pk - s
end

"""
    held_out_residual(cell, jets, Kc; tol) -> (R::Float64, per::Vector{Float64})

The ADR-0028 §2 metric at the COMMON held-out order Kc+1 (an order NEITHER cell
matched): per component, |coeff_{Kc+1}(P_i − f_i·Q)| normalised by ‖Q‖ × the
matched-window coefficient mass + ε.  Combined across components by `max`
(worst component governs the shared Q, ADR-0019).
"""
function held_out_residual(cell::PadeCell{T}, jets, Kc::Int; tol = default_tol(T)) where {T}
    q = cell.denominator
    cnorm = norm(reduce_vcat(jets))
    ε = real(T)(tol) * cnorm
    nq = norm(q)
    per = Float64[]
    for i = 1:length(jets)
        c = jets[i]
        num = abs(_series_coeff(cell.numerators[i], q, Vector{T}(c), Kc + 1))
        mass = sum(abs, @view c[1:min(Kc + 1, length(c))])
        push!(per, Float64(num / (nq * mass + ε)))
    end
    return (maximum(per), per)
end

"""
    held_out_multi(cell, jets, Kc; span, tol) -> Float64

Like `held_out_residual` but probes a RANGE of unmatched orders Kc+1 … Kc+span
(max over orders, then max over components).  Motivated by the audition finding
that a SINGLE held-out coefficient can be nailed by an over-determined least-
squares fit whose overall function is nonetheless inaccurate — one coefficient
is too weak a witness.
"""
function held_out_multi(cell::PadeCell{T}, jets, Kc::Int; span::Int = 3,
                        tol = default_tol(T)) where {T}
    q = cell.denominator
    cnorm = norm(reduce_vcat(jets))
    ε = real(T)(tol) * cnorm
    nq = norm(q)
    worst = 0.0
    for i = 1:length(jets)
        c = Vector{T}(jets[i])
        mass = sum(abs, @view c[1:min(Kc + 1, length(c))])
        for off = 1:span
            k = Kc + off
            k + 1 > length(c) && break
            num = abs(_series_coeff(cell.numerators[i], q, c, k))
            worst = max(worst, Float64(num / (nq * mass + ε)))
        end
    end
    return worst
end

# --- held-out POINT discriminator (redesign candidate) ----------------------

# Horner eval of a low-to-high coefficient vector.
function _horner(c::AbstractVector, t)
    s = zero(t) * zero(eltype(c))
    @inbounds for k = length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

"""
    held_out_point(cell, jets, tstar) -> Float64

Surrogate-truth discriminator: at a held-out point `tstar` INSIDE the Taylor
radius (no pole in (0,tstar]), the FULL jet's Taylor sum (more terms than either
cell matched) is a more accurate local value than either rational; score a cell
by how far its `P_i(tstar)/Q(tstar)` strays from that surrogate truth, worst
component.  Unlike the single-coefficient residual, this integrates ALL held-out
coefficients into a function value — robust to an LS fit that nails one
coefficient.  (Valid only when `tstar` is within radius — true for the
dispatch-relevant regimes: entire jets + regular meromorphic stretches.)
"""
function held_out_point(cell::PadeCell{T}, jets, tstar::Real) where {T}
    t = real(T)(tstar)
    q = _horner(cell.denominator, t)
    worst = 0.0
    for i = 1:length(jets)
        pred = _horner(cell.numerators[i], t) / q
        tru = _horner(Vector{T}(jets[i]), t)
        worst = max(worst, Float64(abs(pred - tru)))
    end
    return worst
end

# --- conditioning gap + Pareto/lexicographic selection (ADR-0028 §3) --------

# g = σ_m/σ_{m+1} (smallest retained σ over smallest σ).  Large ⇒ isolated
# 1-D null space (genuine shared pole); small/→1 ⇒ least-squares regime.
function gap(cell::PadeCell)
    S = cell.S
    length(S) < 2 && return Inf
    return Float64(S[end - 1] / max(S[end], eps()))
end

"""
    select(cands, jets; tol) -> (chosen::PadeCell, table)

Pareto frontier over (R↓, g↑, K↑) with the ADR-0028 §3 deterministic
lexicographic tie-break R→g→K→cell-label, ε_rel = 100·eps banded so jitter
cannot flip the choice.  `cands` is a vector of PadeCell; R is computed at the
common held-out order across all candidates.
"""
function select(cands::Vector{<:PadeCell{T}}, jets; tol = default_tol(T),
                Rfloor::Float64 = 0.0, multi::Bool = false) where {T}
    Kc = minimum(c.K for c in cands)
    Rs = multi ? [held_out_multi(c, jets, Kc; tol = tol) for c in cands] :
                 [held_out_residual(c, jets, Kc; tol = tol)[1] for c in cands]
    gs = [gap(c) for c in cands]
    Ks = [c.K for c in cands]
    εrel = 100 * eps(real(T))
    # lexicographic comparison with ε-bands; returns true if x strictly preferred to y.
    # `Rfloor` (proposed fix for the ε-band finding): when BOTH residuals are below
    # an absolute floor, both cells are effectively exact ⇒ treat R as a tie and let
    # conditioning g decide.  Rfloor=0.0 reproduces the ADR-as-written behaviour.
    function prefer(ix, iy)
        rx, ry = Rs[ix], Rs[iy]
        tie_R = max(rx, ry) ≤ Rfloor || abs(rx - ry) ≤ εrel * max(rx, ry)
        !tie_R && return rx < ry
        gx, gy = gs[ix], gs[iy]
        abs(gx - gy) > εrel * max(gx, gy) && return gx > gy
        Ks[ix] != Ks[iy] && return Ks[ix] > Ks[iy]
        return ix < iy           # final: candidate order = (A) first
    end
    best = 1
    for i = 2:length(cands)
        prefer(i, best) && (best = i)
    end
    table = (R = Rs, g = gs, K = Ks, chosen = cands[best].label)
    return (cands[best], table)
end
