# =============================================================================
# oracle.jl — block-Toeplitz-DETERMINANT shared-Q oracle (probe V1c).
#
# This file is an INDEPENDENT reference implementation of the type-II
# Hermite–Padé (simultaneous-Padé / shared-denominator) approximant, written
# to mutation-prove `PadeTaylor.SharedPade.shared_denominator_pade`.
#
# ## Why a second route at all
#
# `SharedPade.shared_denominator_pade` (`src/SharedPade.jl`) recovers the
# shared denominator `Q` as the *SVD null vector* of a stacked block-Toeplitz
# matrix `A_full`, then refines it with a column-reweighted QR step ported
# from `padeapprox.m`.  An oracle that *also* took an SVD (or used QR as a
# null-space extractor) would not be independent: a sign-flip, an off-by-one
# in the Toeplitz layout, or a wrong-column bug in the SVD+QR path could hide
# behind a matching SVD+QR oracle.  This oracle therefore recovers `Q` by a
# route with **no singular-value decomposition and no QR null-space step
# anywhere** — the determinantal / exact-linear-algebra route.
#
# ## The determinantal construction (Mano–Tsuda 2017, Prop. 2.3)
#
# Mano–Tsuda 2017 §2.2 ("Simultaneous Padé polynomials"), Proposition 2.3,
# gives a closed-form *bordered block-Toeplitz determinant* for the type-II
# denominator polynomial `P^(0)_0` (the object playing the role of `Q`).
# Reading the LaTeX source directly
# (`references/tex/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285/
#   hp_arXiv_final.tex:1490-1503`):
#
#     P^(0)_0(w) = (1/NP^(0)) · det [ 1, w, w², …, wᵐ          ]
#                                   [ A¹_m(n, m+1)              ]
#                                   [ A²_m(n, m+1)              ]
#                                   [        ⋮                 ]
#                                   [ A^{L-1}_m(n, m+1)         ]
#
# The lower part is the stacked block-Toeplitz matrix.  In Mano–Tsuda's
# minimal system each block `Aⁱ_m(n, m+1)` has `n` rows and `m = n(L-1)`, so
# the stack is exactly `m × (m+1)` — square-rectangular, and its maximal
# minors are defined.  `SharedPade`, by contrast, uses the GGT square-block
# convention `n = m` rows per block (`src/SharedPade.jl:31`), so its `A_full`
# is the *taller* `dm × (m+1)`: the extra rows are redundant matching
# equations.  To take maximal minors this oracle first reduces `A_full` to an
# `m × (m+1)` square subsystem by selecting `m` linearly independent rows
# (`_select_independent_rows` — fraction-free elimination, NOT an SVD/QR
# rank-revealing factorisation).  An exact `P/Q` satisfies every matching
# equation, so any `m` independent rows yield the same null vector; the
# selection is a bookkeeping step, not an approximation.  The top border
# `[1, w, …, wᵐ]` makes the selected matrix square, `(m+1)×(m+1)`.
#
# Expanding that determinant along its **top row** (Laplace / cofactor
# expansion) gives a polynomial in `w` whose coefficient of `wᵏ` is the
# signed cofactor of the `(1, k+1)` entry:
#
#     b_k  =  (-1)^k · det( A_full with column (k+1) deleted )
#          =  (-1)^k · M_k                          (k = 0, …, m)
#
# where `M_k` is the `k`-th *maximal minor* of the `m × (m+1)` stacked
# block-Toeplitz matrix `A_full` — the determinant of the `m × m` square
# matrix left after striking out column `k+1`.  The vector
# `b = (b_0, …, b_m)` is, up to the global scale `1/NP^(0)`, the coefficient
# vector of the shared denominator `Q`.
#
# This `b` is, by construction, a right null vector of `A_full`: for each
# row `r` of `A_full`, the dot product `∑_k A_full[r,k+1]·b_k` equals the
# Laplace expansion of an `(m+1)×(m+1)` determinant whose top row is row `r`
# of `A_full` *duplicated* below — a determinant with two equal rows, hence
# zero.  So `A_full · b = 0` exactly.  That is the same defining property
# `SharedPade` extracts via SVD; we obtain it here as a list of `m+1`
# determinants.  **No SVD, no QR-as-nullspace.**  Cramer's rule is the same
# fact in a different dress: fixing `b₀ = 1` and solving the `m × m` square
# subsystem `A_red · (b₁,…,b_m)ᵀ = -A_full[:,1]` by Cramer gives
# `b_j = det(A_red with column j replaced by -A_full[:,1]) / det(A_red)`,
# and `det(A_red) = M_0`, `det(replaced) = (-1)^j M_j` — algebraically the
# identical minor list.  We implement the minor-list form directly: it
# treats every coefficient symmetrically and never needs `b₀ ≠ 0`.
#
# ## Exact / arbitrary-precision arithmetic — why it makes the oracle TIGHT
#
# A determinant of an `m × m` matrix is a polynomial in the entries, so for
# `Rational{BigInt}` input the minors `M_k` are computed **exactly** (no
# roundoff at all) by fraction-free or rational Gaussian elimination — Julia's
# generic `LinearAlgebra.det` on a `Matrix{Rational{BigInt}}` does exactly
# this.  The oracle output is then the *exact* rational `Q` and `Pᵢ`, not a
# roundoff-limited estimate.  For `BigFloat` input the determinants are
# evaluated at extended precision, far below Float64 epsilon.  This is what
# lets V1e assert agreement with `SharedPade` (a Float64 SVD path) at ~1e-12
# and *exact* equality for rational inputs — the oracle is not itself a
# source of numerical error, so any discrepancy is a real `SharedPade` bug.
#
# ## Layout equivalence with SharedPade._toeplitz_block
#
# `SharedPade._toeplitz_block(c, m)` builds an `m × (m+1)` block with
# `block[r,c] = c_{m+r-c}` (1-based: `c[idx]`, `idx = m+r-c+1`).  Mano–Tsuda's
# rectangular Toeplitz `A^i_j(k,l)` has `(r,c)`-entry `a^i_{j+r-c}`
# (`hp_arXiv_final.tex:1061-1077`).  The denominator block in eq. (2.6) /
# Prop. 2.3 is `A^i_m(n, m+1)` — i.e. `j = m`, `k = n` rows, `l = m+1`
# columns.  `SharedPade` uses `n = m` rows per block (square-block GGT
# convention, `src/SharedPade.jl:31-44`).  With `j = m, k = n = m, l = m+1`
# the Mano–Tsuda entry is `a^i_{m+r-c}` — *identical* to `c_{m+r-c}`.  This
# oracle therefore mirrors `_toeplitz_block` verbatim (`_toeplitz_block`
# below) so the two routes operate on the *same* `A_full`; only the
# null-vector extraction differs.  That is precisely the disjointness V1c
# needs: same matrix, algorithmically independent solve.
#
# ## Output convention (matches SharedPade for direct V1e comparison)
#
# `shared_q_via_determinant` returns `(numerators, denominator)` with:
#   - coefficients low-to-high (`z⁰` first);
#   - denominator normalised so `Q(0) = b[1] = 1`;
#   - trailing near-zero coefficients trimmed;
# exactly as `shared_denominator_pade` does, so V1e can compare the two
# tuples directly (by function value / pole location — see verify.jl).
#
# ## Fail loud (CLAUDE.md Rule 1)
#
# If every maximal minor `M_k` vanishes, `A_full` has rank `< m`: there is
# no unique shared `Q` (the system is not perfect at this multi-index —
# Mano–Tsuda p. 12, "unique iff rank = m").  We throw rather than return a
# zero/NaN lie.  Likewise if the normalising minor `M_0 = Q(0)` vanishes —
# a pole sits at the expansion centre, so `Q(0)=1` normalisation is
# impossible — we throw with a suggestion.
#
# ## References
#   - `references/tex/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285/
#      hp_arXiv_final.tex:1365-1528` — §2.2 Simultaneous Padé polynomials,
#      eq. (2.6) (the `lspa` system) and Proposition 2.3 (the bordered
#      block-Toeplitz determinant for `P^(0)_0`).
#   - `references/tex/.../hp_arXiv_final.tex:1322-1359` — Remark 2.2, the
#      block-Toeplitz determinant `Δ^(i)` and `NP^(0) = ±Δ^(L)`.
#   - `docs/v0p2_pillarA_hermite_pade_findings.md:116-134` — §2
#      "Determinantal Representation": the same Prop. 2.3 formula, paraphrased.
#   - `src/SharedPade.jl:99-141` — `_toeplitz_block` / `_upper_block`, the
#      block layout this oracle mirrors.
# =============================================================================

module SharedPadeDeterminantOracle

using LinearAlgebra: det

export shared_q_via_determinant

# -----------------------------------------------------------------------------
# Toeplitz block builder — a verbatim mirror of SharedPade._toeplitz_block.
# -----------------------------------------------------------------------------

"""
    _toeplitz_block(c, m) -> Matrix

The `m × (m+1)` block-Toeplitz block for one component, `block[r,c]` equal
to the Taylor coefficient `c_{m+r-c}` (zero outside the available range).

This is a verbatim mirror of `SharedPade._toeplitz_block` (`src/SharedPade.jl`),
which is in turn `A^i_m(n, m+1)` of Mano–Tsuda 2017 eq. (2.6) at `n = m`
(`hp_arXiv_final.tex:1405-1424`).  Mirroring it exactly is deliberate: the
two routes (`SharedPade`'s SVD and this oracle's determinant) must act on the
*same* stacked matrix, so the only difference is the null-vector solve.
"""
function _toeplitz_block(c::AbstractVector{T}, m::Int) where {T}
    blk = zeros(T, m, m + 1)
    @inbounds for cc = 1:(m + 1), rr = 1:m
        idx = m + rr - cc + 1            # 1-based index into c (c[1] = c_0)
        if 1 ≤ idx ≤ length(c)
            blk[rr, cc] = T(c[idx])
        end
    end
    return blk
end

"""
    _upper_block(c, m) -> Matrix

The `(m+1) × (m+1)` upper Toeplitz block recovering a numerator from the
shared `b`: `a = upper · b`, with `upper[r,c] = c_{r-c}` (zero for `r < c`).
A verbatim mirror of `SharedPade._upper_block` — the numerator recovery
step `Pᵢ = [fᵢ · P^(0)_0]` truncated to degree `m`
(Mano–Tsuda 2017, the `[ fᵢ P₀^(0) ]^{m-1}_0` formula, `hp_arXiv_final.tex:1539`).
"""
function _upper_block(c::AbstractVector{T}, m::Int) where {T}
    up = zeros(T, m + 1, m + 1)
    @inbounds for cc = 1:(m + 1), rr = cc:(m + 1)
        idx = rr - cc + 1
        if 1 ≤ idx ≤ length(c)
            up[rr, cc] = T(c[idx])
        end
    end
    return up
end

# -----------------------------------------------------------------------------
# Row selection — pick m independent rows from the dm-row stacked matrix.
# -----------------------------------------------------------------------------

"""
    _select_independent_rows(A, want) -> Vector{Int}

Return `want` row indices of `A` that are linearly independent, found by
**fraction-free Gaussian elimination with partial (column-by-column) pivot
search** — Bareiss-style elimination, *not* an SVD or a QR rank-revealing
factorisation.

The Mano–Tsuda §2.2 simultaneous-Padé system (eq. 2.6,
`hp_arXiv_final.tex:1405-1424`) is *square-rectangular*: with block height
`n` and `m = n(L-1)` the stacked matrix `A^1_m(n,m+1);…;A^{L-1}_m(n,m+1)`
has exactly `m` rows and `m+1` columns, so its maximal minors are defined.
`SharedPade`, however, uses the GGT square-block convention `n = m` rows per
block — giving a *taller* `dm × (m+1)` stacked matrix (`src/SharedPade.jl:31`).
The extra rows are redundant matching equations (an exact `P/Q` satisfies
more than the minimal `m`).  To take maximal minors we must first reduce to
an `m × (m+1)` square subsystem by selecting `m` independent rows — exactly
the bead-text instruction "select `m` independent rows of the stacked
matrix".

The selection is exact for `Rational{BigInt}` input: a row is accepted iff,
after eliminating it against the already-chosen pivot rows, its leading
remaining entry is nonzero.  No tolerance, no SVD.  For `AbstractFloat`
input the pivot search uses the largest-magnitude entry and a relative
threshold — the floating path is a numerical convenience; the *exact*
`Rational` path is the load-bearing one for the tight-oracle claim.

Throws (Rule 1) if fewer than `want` independent rows exist — then
`rank(A) < m` and no unique shared `Q` exists.
"""
function _select_independent_rows(A::AbstractMatrix{T}, want::Int) where {T}
    nrows, ncols = size(A)
    want ≤ nrows || throw(ArgumentError(
        "_select_independent_rows: asked for $want rows but A has only $nrows."))
    # Work on a mutable copy; eliminate accepted rows out of the rest so the
    # next acceptance test sees the residual (independence) directly.
    work   = Matrix{T}(A)
    chosen = Int[]
    isfloat = T <: AbstractFloat
    # Relative zero-threshold for the floating path (exact path: literal 0).
    scale = isfloat ? maximum(abs, A; init = zero(real(T))) : zero(real(T))
    ztol  = isfloat ? sqrt(eps(real(T))) * (scale + one(real(T))) : zero(real(T))
    iszero_entry(x) = isfloat ? abs(x) ≤ ztol : iszero(x)

    used_pivot_col = Int[]
    for r = 1:nrows
        length(chosen) == want && break
        # Find a column where this (already-reduced) row has a nonzero entry,
        # skipping columns already used as pivots.
        pcol = 0
        if isfloat
            best = ztol
            for c = 1:ncols
                c in used_pivot_col && continue
                if abs(work[r, c]) > best
                    best = abs(work[r, c]); pcol = c
                end
            end
        else
            for c = 1:ncols
                c in used_pivot_col && continue
                if !iszero(work[r, c])
                    pcol = c; break
                end
            end
        end
        pcol == 0 && continue           # row is dependent on chosen rows
        push!(chosen, r)
        push!(used_pivot_col, pcol)
        piv = work[r, pcol]
        # Eliminate this pivot column from all *later* rows so future
        # independence tests are against the residual.
        for rr = (r + 1):nrows
            factor = work[rr, pcol] / piv
            iszero_entry(factor) && continue
            for c = 1:ncols
                work[rr, c] -= factor * work[r, c]
            end
        end
    end
    length(chosen) == want || throw(ErrorException(
        "_select_independent_rows: only $(length(chosen)) independent rows " *
        "found, need $want; rank(A_full) < m, so no unique shared " *
        "denominator exists (Mano–Tsuda 2017 p. 12 — unique iff rank = m). " *
        "Suggestion: reduce the denominator degree m or extend the jets."))
    return chosen
end

# -----------------------------------------------------------------------------
# Maximal minors via cofactor / Cramer route — the SVD-free core.
# -----------------------------------------------------------------------------

"""
    _maximal_minors(A) -> Vector

Given an `m × (m+1)` matrix `A`, return the `m+1` *signed maximal minors*

    b_k = (-1)^k · det( A with column (k+1) deleted )      (k = 0, …, m)

The vector `b = (b_0, …, b_m)` is the cofactor expansion, along the top
row, of the bordered `(m+1)×(m+1)` determinant

    det [ 1, w, …, wᵐ ;  A ]   (Mano–Tsuda 2017 Prop. 2.3),

so `b` is a right null vector of `A`: `A·b = 0` exactly (a determinant with
two equal rows vanishes — see the module docstring).  This is the
shared-denominator coefficient vector, up to a global scale.

Each `det(...)` is an `m × m` determinant.  For `T = Rational{BigInt}` Julia
evaluates it by exact fraction-free / rational elimination, so `b` is the
*exact* coefficient vector; for `BigFloat` it is evaluated at extended
precision.  **No SVD, no QR.**
"""
function _maximal_minors(A::AbstractMatrix{T}) where {T}
    m, mp1 = size(A)
    mp1 == m + 1 || throw(ArgumentError(
        "_maximal_minors: expected an m×(m+1) matrix; got $(m)×$(mp1)."))
    allcols = collect(1:mp1)
    b = Vector{T}(undef, mp1)
    for k = 0:m
        keep = deleteat!(copy(allcols), k + 1)        # columns ≠ k+1
        sub  = A[:, keep]                             # m × m square block
        b[k + 1] = (-1)^k * det(sub)
    end
    return b
end

# -----------------------------------------------------------------------------
# Main entry point.
# -----------------------------------------------------------------------------

"""
    shared_q_via_determinant(jets, m; tol = nothing)
        -> (numerators::Vector{Vector{T}}, denominator::Vector{T})

Compute a shared-denominator (type-II Hermite–Padé) approximant for a
`d`-component vector of formal power series, **by the block-Toeplitz
determinant route** (Mano–Tsuda 2017 Proposition 2.3) — independent of the
SVD path in `SharedPade.shared_denominator_pade`.

`jets[i] = [c₀, c₁, …]` are the Taylor coefficients of component `i`; each
must have length `≥ m+1`.  `m` is the shared denominator degree.

Returns `d` numerator coefficient vectors and one denominator `Q`, all
low-to-high and normalised so `Q(0) = 1`, with trailing near-zeros trimmed —
the *same* output convention as `shared_denominator_pade`, so V1e can
compare the two directly.

`tol` is the trailing-zero-trim threshold.  For exact (`Rational`) element
types pass `tol = 0` (the default for non-`AbstractFloat` `T`) so nothing
is trimmed spuriously; for floating types it defaults to `√eps`.

Throws (CLAUDE.md Rule 1) when:
  - any jet is shorter than `m+1` (cannot fill the block);
  - all `m+1` maximal minors vanish — `rank(A_full) < m`, so no unique
    shared `Q` exists (Mano–Tsuda p. 12, "unique iff rank = m");
  - `Q(0)` (the minor `M_0`) vanishes — a pole sits at the expansion centre.

Algorithm: build the `dm × (m+1)` stacked block-Toeplitz matrix `A_full`
(identical to `SharedPade`'s), take its `m+1` signed maximal minors as the
denominator coefficient vector, normalise to `Q(0)=1`, then recover each
numerator from the upper Toeplitz block.
"""
function shared_q_via_determinant(jets::AbstractVector{<:AbstractVector{T}},
                                  m::Integer;
                                  tol = nothing) where {T}
    d = length(jets)
    d ≥ 1 || throw(ArgumentError(
        "shared_q_via_determinant: need at least one jet; got d=0. " *
        "Suggestion: pass a non-empty vector of Taylor-coefficient vectors."))
    m ≥ 1 || throw(ArgumentError(
        "shared_q_via_determinant: denominator degree m must be ≥ 1; got m=$m."))

    mi = Int(m)
    for (i, jet) in enumerate(jets)
        length(jet) ≥ mi + 1 || throw(ArgumentError(
            "shared_q_via_determinant: jet $i too short for denominator " *
            "degree m=$mi: length(jet)=$(length(jet)), need ≥ m+1=$(mi+1). " *
            "Suggestion: extend the Taylor jet to at least m+1 coefficients."))
    end

    # Trim/normalisation threshold.  Exact types ⇒ 0 (trim only true zeros);
    # floating types ⇒ √eps unless the caller overrides.
    trimtol = if tol !== nothing
        tol
    elseif T <: AbstractFloat
        sqrt(eps(T))
    else
        zero(real(float(one(T))))
    end

    # --- Build the dm × (m+1) stacked block-Toeplitz matrix A_full ----------
    # Identical layout to SharedPade.shared_denominator_pade.
    blocks = [_toeplitz_block(jets[i], mi) for i = 1:d]
    A_full = reduce(vcat, blocks)                  # (d·m) × (m+1)

    # --- Reduce to an m × (m+1) square subsystem ----------------------------
    # Maximal minors are defined for an m×(m+1) matrix.  SharedPade's GGT
    # square-block convention makes A_full taller (dm rows); the extra rows
    # are redundant matching equations.  Select m independent rows by
    # fraction-free elimination (NOT SVD) — Mano–Tsuda §2.2's minimal
    # simultaneous-Padé system is exactly m×(m+1).  This `_select_*` call
    # also surfaces rank deficiency (rank < m ⇒ throw).
    sel    = _select_independent_rows(A_full, mi)
    A_sq   = A_full[sel, :]                        # m × (m+1)

    # --- Denominator b = signed maximal minors of A_sq (SVD-free) -----------
    b = _maximal_minors(A_sq)

    # Failure mode: every maximal minor is zero ⇒ rank(A_full) < m ⇒ the
    # shared denominator's null space is not 1-dimensional.  Fail loud.
    is_zero(x) = T <: AbstractFloat ? abs(x) ≤ trimtol : iszero(x)
    all(is_zero, b) && throw(ErrorException(
        "shared_q_via_determinant: every maximal minor of the stacked " *
        "block-Toeplitz matrix vanishes; rank(A_full) < m, so no unique " *
        "shared denominator exists. " *
        "Detail: the simultaneous-Padé system is not perfect at this " *
        "multi-index (Mano–Tsuda 2017 p. 12 — unique iff rank = m). " *
        "Suggestion: reduce the denominator degree m or extend the jets."))

    # Failure mode: Q(0) = M_0 ≈ 0 ⇒ a pole at the expansion centre.
    is_zero(b[1]) && throw(ErrorException(
        "shared_q_via_determinant: recovered Q(0) = b[1] = $(b[1]) vanishes; " *
        "the shared denominator has a root at the expansion centre, so the " *
        "Q(0)=1 normalisation is impossible. " *
        "Suggestion: expand the series at a different base point."))

    # --- Normalise to Q(0)=1, trim trailing near-zeros ----------------------
    denominator = b ./ b[1]                         # b[1] = 1 exactly now
    keep_b = T <: AbstractFloat ?
        findlast(x -> abs(x) > trimtol, denominator) :
        findlast(!iszero, denominator)
    keep_b === nothing && error(
        "shared_q_via_determinant: denominator collapsed after normalisation;" *
        " algorithm bug — file a bead with the input jets.")
    denominator = denominator[1:keep_b]

    # --- Recover numerators from the upper Toeplitz block -------------------
    # Pᵢ = [ fᵢ · P^(0)_0 ] truncated to degree m  (Mano–Tsuda Prop. 2.3,
    # the `[ fᵢ P₀^(0) ]^{m-1}_0` numerator formula).  `_upper_block · b`
    # performs that truncated polynomial multiplication; divide by b[1] to
    # carry the Q(0)=1 normalisation through to the numerators.
    numerators = Vector{Vector{T}}(undef, d)
    for i = 1:d
        a_i = (_upper_block(jets[i], mi) * b) ./ b[1]
        keep_a = T <: AbstractFloat ?
            findlast(x -> abs(x) > trimtol, a_i) :
            findlast(!iszero, a_i)
        numerators[i] = keep_a === nothing ? T[zero(T)] : a_i[1:keep_a]
    end

    return (numerators, denominator)
end

end # module SharedPadeDeterminantOracle
