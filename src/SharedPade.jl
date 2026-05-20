"""
    PadeTaylor.SharedPade

Shared-denominator robust Padé — the keystone primitive lifting Padé
approximation from scalar to vector.

## Why one shared denominator

A vector meromorphic ODE solution `f = (f₁,…,f_d)` has all `d` components
blowing up at the *same* movable poles: the singularities are a property
of the flow, not of the individual component.  The right approximant is
therefore `d` numerator polynomials `P₁,…,P_d` over **one** shared
denominator `Q` — a type-II Hermite–Padé (simultaneous Padé) approximant
in the sense of Mano–Tsuda 2017 eq. (1.2)
(`references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`,
p. 4): `P^(0)` plays the role of `Q`, and `P^(j)/P^(0) ≈ f_j/f₀`.  The
shared `Q` is what makes the recovered pole set a single, consistent
estimate of the system's singularities — fitting each component
independently (e.g. with AAA) gives `d` *different* pole sets that only
agree up to noise.

## The block-Toeplitz construction

The scalar GGT 2013 algorithm (`RobustPade.jl`, `method=:svd`) solves the
matching condition `P(z) = f(z)Q(z) + O(z^{m+n+1})` as a homogeneous
Toeplitz null-space problem: the `n×(n+1)` matrix `C̃` with `C̃[r,c] =
c_{m+r-c}` has its denominator `Q` as the smallest-singular-value right
null vector.

Mano–Tsuda §2.2 (eq. 2.6, p. 12) generalises this directly.  For `d`
components we stack one `m×(m+1)` Toeplitz block per component vertically:

```
        ┌ block₁ ┐   ← m rows from jets[1], block[r,c] = c⁽¹⁾_{m+r-c}
A_full =│ block₂ │   ← m rows from jets[2]
        │   ⋮    │
        └ block_d ┘   ← m rows from jets[d]
```

`A_full` is `dm × (m+1)`.  Its smallest-singular-value right null vector
is the shared denominator `Q` — every block contributes `m` matching
equations against the *same* unknown coefficient vector `b₀,…,bₘ`, which
is exactly the constraint "all components share `Q`".  One SVD of one
rectangular matrix; no coupled or iterated solves (pillar A §3).

## Reduction to the scalar path (the primary correctness oracle)

For `d = 1`, `A_full` is a single `m×(m+1)` block, *identical* to the
scalar GGT `C̃` at `n = m`.  The SVD, the QR-reweighting, the numerator
recovery from the upper Toeplitz block, and the `b[1]=1` normalisation
are all the verbatim scalar operations.  Hence `shared_denominator_pade`
with one jet reproduces `robust_pade(jet, m, m; method=:svd)` bit-for-bit
modulo SVD sign conventions — this is asserted in `SP.1.1`.

## The QR-reweighting (ported unchanged from padeapprox.m)

After the SVD picks the null vector `b₀ = Vt[end,:]`, Chebfun's
`padeapprox.m` (lines 278–280, *beyond* GGT 2013 Algorithm 2 itself)
refines it: `D = diag(|b₀| + √ε)`, `QR((A_full·D)ᵀ)`, `b = D·Q[:,end]`.
The column-reweighting better preserves the genuine exact zeros of `b`
for blocks at accuracy near the tolerance (pillar A §4, last bullet:
"the QR reweighting step ports without change to A_full").  We port it
unchanged from the scalar `RobustPade.jl` path.

## Defensive throws (Rule 1 — fail loud)

Four failure modes throw rather than return a NaN/zero lie, each per
pillar A §7 "Failure modes":

  1. **jet too short** — `length(jets[i]) < m+1`: cannot fill the block.
  2. **all-zero jets** — every singular value below `τ = tol·‖c‖`: the
     input is indistinguishable from zero at the requested tolerance.
  3. **Q(0) ≈ 0** — recovered `b[1]` below `tol`: a pole sits at the
     expansion centre, so `Q(0)=1` normalisation would divide by ~0.
  4. **non-isolated null space** — `ρ < m` after degree reduction
     leaves the smallest singular value not isolated: the shared
     denominator is not unique (the system is not perfect at this
     multi-index, Mano–Tsuda p. 12 "unique iff rank = m").

## References

  - `docs/v0p2_pillarA_hermite_pade_findings.md:46–127, 379–472` —
    the spec: block construction, algorithm steps 1–8, failure modes.
  - `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`
    §2.2 eq. (2.6), p. 12 — the block-Toeplitz null-space system.
  - `src/RobustPade.jl` — the scalar `d=1` special case; conventions
    (`_lower_tri_toeplitz`, QR-reweighting, `b[1]=1`) mirrored here.
  - `external/chebfun/padeapprox.m` lines 278–280 — QR-reweighting port.
"""
module SharedPade

using LinearAlgebra: norm, qr, Diagonal, adjoint
using ..LinAlg:      pade_svd
using ..RobustPade:  default_tol

export shared_denominator_pade

# -----------------------------------------------------------------------------
# Toeplitz block builder
# -----------------------------------------------------------------------------

"""
    _toeplitz_block(c, m) -> Matrix

The `m × (m+1)` Toeplitz block for one component, with `block[r,c]` equal
to the Taylor coefficient `c_{m+r-c}` (zero outside the available range).

This is the bottom `m` rows of the scalar GGT lower-triangular Toeplitz
`Z` at `(μ,ν) = (m,m)`: `RobustPade` builds `Z` of size `(2m+1)×(m+1)`
with `Z[i,j]=c[i-j+1]` and slices `C = Z[m+2:2m+1, :]`, giving
`C[r,c] = c[m+1+r-c+1] = c_{m+r-c}`.  We construct that slice directly so
the `d=1` stack is bit-identical to the scalar `C̃`.
"""
function _toeplitz_block(c::AbstractVector{T}, m::Int) where {T}
    blk = zeros(T, m, m + 1)
    @inbounds for cc = 1:(m + 1), rr = 1:m
        idx = m + rr - cc + 1            # 1-based index into c (c[1]=c_0)
        if 1 ≤ idx ≤ length(c)
            blk[rr, cc] = T(c[idx])
        end
    end
    return blk
end

"""
    _upper_block(c, m) -> Matrix

The `(m+1) × (m+1)` upper Toeplitz block recovering a numerator from the
shared `b`: `a = upper · b`, with `upper[r,c] = c_{r-c}` (zero for `r<c`).
Mirrors `RobustPade`'s `Z[1:m+1, 1:n+1]` numerator step at `n=m`.
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
# Main entry point
# -----------------------------------------------------------------------------

"""
    shared_denominator_pade(jets, m; tol = default_tol(T))
        -> (numerators::Vector{Vector{T}}, denominator::Vector{T})

Compute a shared-denominator (simultaneous / type-II Hermite–Padé)
approximant for a `d`-component vector of formal power series.

`jets[i] = [c₀, c₁, …]` are the Taylor coefficients of component `i`;
each must have length `≥ m+1`.  `m` is the (shared) denominator degree.
Returns `d` numerator coefficient vectors `P₁,…,P_d` (low-to-high) and a
single denominator `Q`, all normalised so `Q(0) = b[1] = 1`.

For `d = 1` this reduces exactly to `robust_pade(jets[1], m, m;
method=:svd)`.  Algorithm = pillar A §7 steps 1–8; see the module
docstring for the block-Toeplitz rationale and the four defensive throws.
"""
function shared_denominator_pade(jets::AbstractVector{<:AbstractVector{T}},
                                 m::Integer;
                                 tol::Real = default_tol(T)) where {T}
    d = length(jets)
    d ≥ 1 || throw(ArgumentError(
        "shared_denominator_pade: need at least one jet; got d=0. " *
        "Suggestion: pass a non-empty vector of Taylor-coefficient vectors."))
    m ≥ 1 || throw(ArgumentError(
        "shared_denominator_pade: denominator degree m must be ≥ 1; got m=$m."))

    # --- Failure mode 1: jet too short (pillar A §7) -------------------------
    for (i, jet) in enumerate(jets)
        length(jet) ≥ m + 1 || throw(ArgumentError(
            "shared_denominator_pade: jet $i too short for denominator " *
            "degree m=$m: length(jet)=$(length(jet)), need ≥ m+1=$(m+1). " *
            "Suggestion: extend the Taylor jet to at least m+1 coefficients."))
    end

    tol_t = real(T) <: AbstractFloat ? real(T)(tol) : tol
    cnorm = norm(reduce(vcat, [collect(j) for j in jets]))
    τ = tol_t * cnorm

    # --- Steps 2–4: build A_full, SVD, rank check + degree reduction --------
    # Degree reduction is the vector analogue of RobustPade's `while` loop:
    # if the smallest-σ null space is not isolated (ρ < m), shrink m and
    # rebuild (pillar A §7 step 4 / §4 "reduce m until σ_min is isolated").
    local A_full, S, Vt, m_cur
    m_cur = Int(m)
    while true
        blocks = [_toeplitz_block(jets[i], m_cur) for i = 1:d]
        A_full = reduce(vcat, blocks)            # (d·m_cur) × (m_cur+1)
        _, S, Vt = pade_svd(A_full; full = true)
        # ρ = number of singular values above the noise threshold τ.
        ρ = count(s -> s > τ, S)

        # Failure mode 2: all singular values below τ — the jets are
        # indistinguishable from zero at this tolerance.
        ρ == 0 && throw(ErrorException(
            "shared_denominator_pade: every singular value is below " *
            "τ = tol·‖c‖ = $τ; the Taylor jets are indistinguishable from " *
            "zero at tolerance tol=$tol_t. " *
            "Suggestion: check the jets are non-trivial, or loosen tol."))

        # A_full has m_cur+1 columns; a unique shared denominator (up to
        # scale) needs rank exactly m_cur, i.e. an isolated 1-D null space.
        if ρ ≥ m_cur
            break
        end
        # Degree reduction: drop to ρ and rebuild.  m_cur stays ≥ 1 because
        # ρ ≥ 1 here (ρ == 0 threw above).
        m_cur = ρ
    end

    # --- Step 5: QR-reweighting (padeapprox.m lines 278–280) ----------------
    # Ported unchanged from RobustPade's `:svd` path.  `Vt[end,:]` is the
    # smallest-singular-value right null vector = denominator estimate b₀.
    b_init = Vector{T}(Vt[end, :])

    # Failure mode 4: null space not isolated.  After degree reduction we
    # have ρ ≥ m_cur; if ρ > m_cur the smallest σ is degenerate (multiple
    # near-equal σ at the bottom) — the shared denominator is not unique.
    n_near = count(s -> s ≤ τ, S)
    n_near > 1 && throw(ErrorException(
        "shared_denominator_pade: $n_near singular values lie at/below " *
        "τ = $τ; the shared denominator's null space is not isolated. " *
        "Detail: the system may not be perfect at this multi-index " *
        "(Mano–Tsuda 2017 p. 12 — unique iff rank = m). " *
        "Suggestion: reduce the denominator degree m or extend the jets."))

    eps_T = real(T) <: AbstractFloat ? sqrt(eps(real(T))) : sqrt(eps(Float64))
    D = Diagonal([abs(bk) + eps_T for bk in b_init])
    F = qr(adjoint(A_full * D))                  # (A_full·D)ᵀ in MATLAB
    b = D * F.Q[:, m_cur + 1]
    b ./= norm(b)

    # --- Steps 6–7: recover numerators, then normalise to Q(0)=1 ------------
    # Failure mode 3: Q(0) = b[1] ≈ 0 — a pole at the expansion centre.
    abs(b[1]) > τ || abs(b[1]) > tol_t || throw(ErrorException(
        "shared_denominator_pade: recovered Q(0) = b[1] ≈ $(b[1]) is below " *
        "tolerance; the shared denominator has (effectively) a root at the " *
        "expansion centre. " *
        "Suggestion: expand at a different base point — Q(0)=0 means a pole " *
        "sits at the expansion centre."))

    # Trailing near-zero trimming (pillar A §7 step 7), mirroring
    # RobustPade._trim_and_normalise so the d=1 stack reduces *exactly* to
    # the scalar `:svd` output: the scalar path drops trailing |bₖ|≤tol from
    # Q and trailing |aₖ|≤τ from each numerator.  No leading-z^λ cancellation
    # is needed here — `b[1]=Q(0)` is normalised to 1 (the Q(0)≈0 guard above
    # already rejected a root at the expansion centre).
    denominator = Vector{T}(b ./ b[1])           # b[1] = 1 exactly now
    last_b = findlast(x -> abs(x) > tol_t, denominator)
    last_b === nothing && error(
        "shared_denominator_pade: denominator collapsed below tol after " *
        "normalisation; algorithm bug — file a bead with the input jets.")
    denominator = denominator[1:last_b]

    numerators = Vector{Vector{T}}(undef, d)
    for i = 1:d
        a_i = Vector{T}((_upper_block(jets[i], m_cur) * b) ./ b[1])
        last_a = findlast(x -> abs(x) > τ, a_i)
        numerators[i] = last_a === nothing ? T[zero(T)] : a_i[1:last_a]
    end

    return (numerators, denominator)
end

end # module SharedPade
