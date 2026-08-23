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
c_{m+r-c+1}` (top-left `C̃[1,1] = c_{m+1}`, i.e. the matching window
`z^{m+1}…z^{m+n}`) has its denominator `Q` as the smallest-singular-value
right null vector.

Mano–Tsuda §2.2 (eq. 2.6, p. 12) generalises this directly.  For `d`
components we stack one `m×(m+1)` Toeplitz block per component vertically:

```
        ┌ block₁ ┐   ← m rows from jets[1], block[r,c] = c⁽¹⁾_{m+r-c+1}
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
are all the verbatim scalar operations.  Hence, *in the full-rank,
non-degree-reduced regime*, `shared_denominator_pade` with one jet reproduces
`robust_pade(jet, m, m; method=:svd)` bit-for-bit modulo SVD sign conventions
(asserted in `SP.1.1` and `SP.1.7`).

This bit-for-bit equality is **scoped to the full-rank regime** (bug
`padetaylor-jznu`, pinned by `SP.1.8`).  In two DEGENERATE regimes the two
*separate* constructions can diverge, because the scalar path truncates its
input to `m+n+1` coefficients and jointly reduces `(m, n)`, whereas this path
keeps the full jet and reduces a single `m_cur`:

  * **locally-regular collapse** (`Q = [1]`): the scalar path returns a
    jointly-reduced numerator (e.g. `a = [4.0]` — a degree-0 constant) while
    this path returns the FULL Taylor numerator.  These agree only at `z = 0`
    and diverge off-centre (reldiff ~2e-5 by `z ≈ 0.07`) — different in
    DEGREE *and* VALUE;
  * **trim-boundary straddle**: a trailing denominator coefficient sitting
    within ~`tol` (≈1e-14, LAPACK-roundoff-dependent) is kept by one path and
    trimmed by the other, so the rationals are VALUE-equal but differ by one
    in DEGREE.

Both representations are individually Rule-1-correct (neither returns a
NaN/zero/garbage lie); the divergence is in the rational REPRESENTATION (and,
locally-regular, in off-centre value), not a numerical breakdown.  Downstream
this is benign: the stepper evaluates real ODE jets (full-rank in practice),
and V7/V8 pole extraction roots only the DENOMINATOR (a constant `Q = [1]` has
no roots).

## The QR-reweighting (ported unchanged from padeapprox.m)

After the SVD picks the null vector `b₀ = Vt[end,:]`, Chebfun's
`padeapprox.m` (lines 111–117, *beyond* GGT 2013 Algorithm 2 itself)
refines it: `D = diag(|b₀| + √ε)`, `QR((A_full·D)ᵀ)`, `b = D·Q[:,end]`.
The column-reweighting better preserves the genuine exact zeros of `b`
for blocks at accuracy near the tolerance (pillar A §4, last bullet:
"the QR reweighting step ports without change to A_full").  We port it
unchanged from the scalar `RobustPade.jl` path.

## Graceful reduction and defensive throws (Rule 1 — fail loud)

Degree reduction (ADR-0027) makes the construction robust to the common
case of a locally-regular / entire jet, which has no genuine degree-`m`
shared denominator.  The loop reduces only when the stack is genuinely
rank-deficient (`ρ = count(σ>τ) < m_cur` ⇒ drop to `ρ`, the vector analogue
of GGT's "reduce `n` to `ρ`"), and accepts when `ρ ≥ m_cur`: either an
isolated 1-D null space (`ρ == m_cur`, genuine shared poles) or full column
rank (`ρ == m_cur+1`, no exact null — the smallest-σ right vector is the
legitimate *least-squares* best shared denominator, e.g. an entire jet; see
bead `padetaylor-unk`).  A spurious pole at the expansion centre
(`Q(0) = b[1] ≈ 0`) that the least-squares `Q` may carry is handled by
cancelling the common `z^λ` factor from `b` and every numerator (mirroring
the scalar `:svd` path), **not** by throwing.

Three failure modes still throw rather than return a NaN/zero lie:

  1. **jet too short** — `length(jets[i]) < m+1`: cannot fill the block.
  2. **identically-zero jets** — `‖c‖₂ = 0` exactly: nothing to approximate.
     The predicate is deliberately `iszero`, not `‖c‖ ≤ tol`: GGT 2013's
     thresholds are all RELATIVE (`τ = tol·‖c‖₂`, md:213 and Algorithm 2
     step 2, md:225), so the whole construction is homogeneous under
     `c ↦ αc`, and an absolute zero-test was the one line that broke that
     covariance — it falsely rejected valid small-amplitude jets with
     `‖c‖ < tol` (bug `padetaylor-ked0`; pinned by `SP.6`, α ∈ {1e-20, 1e20}).
  3. **non-finite input** — any `Inf`/`NaN` coefficient, or `tol ∉ [0, ∞)`:
     an upstream Taylor-jet overflow must surface as an `ArgumentError`
     ("reduce h or order"), never as a silent zero / garbage `Q` (bug
     `padetaylor-lbqb`; pinned by `SP.7`).  Same guard as `robust_pade`.

## References

  - `docs/v0p2_pillarA_hermite_pade_findings.md:46–127, 379–472` —
    the spec: block construction, algorithm steps 1–8, failure modes.
  - `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`
    §2.2 eq. (2.6), p. 12 — the block-Toeplitz null-space system.
  - `src/RobustPade.jl` — the scalar `d=1` special case; conventions
    (`_lower_tri_toeplitz`, QR-reweighting, `b[1]=1`) mirrored here.
  - `external/chebfun/padeapprox.m` lines 111–117 — QR-reweighting port.
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
to the Taylor coefficient `c_{m+r-c+1}` (zero outside the available
range); the top-left entry is `block[1,1] = c_{m+1}`.

This is the bottom `m` rows of the scalar GGT lower-triangular Toeplitz
`Z` at `(μ,ν) = (m,m)`: `RobustPade` builds `Z` of size `(2m+1)×(m+1)`
with `Z[i,j]=c[i-j+1]` and slices `C = Z[m+2:2m+1, :]`, giving
`C[r,c] = c[m+1+r-c+1] = c_{m+r-c+1}` (1-based array index `m+r-c+2`).
We construct that slice directly so the `d=1` stack is bit-identical to
the scalar `C̃`, which is the matching window `z^{m+1}…z^{2m}` for the
`(m,m)` approximant.
"""
function _toeplitz_block(c::AbstractVector{T}, m::Int) where {T}
    blk = zeros(T, m, m + 1)
    @inbounds for cc = 1:(m + 1), rr = 1:m
        idx = m + rr - cc + 2            # c[idx]=c_{m+r-c+1}; top-left c_{m+1} (GGT (m,m) window)
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

For `d = 1`, full-rank inputs without degree reduction reproduce the scalar
`robust_pade` path. Degenerate inputs can instead reduce to the full Taylor
jet over `Q = [1]` (`docs/worklog/067-sharedpade-c1-c2-minimal-core.md:36-40`).
Algorithm = pillar A §7 steps 1–8; see the module
docstring for the block-Toeplitz rationale and the defensive throws.
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
        # Bug padetaylor-lbqb: an Inf/NaN coefficient would make every σ-vs-τ
        # comparison vacuous and could yield a silent garbage Q (Rule 1).
        all(isfinite, jet) || throw(ArgumentError(
            "shared_denominator_pade: jet $i contains Inf/NaN coefficients — " *
            "the Taylor jet overflowed; reduce h or order."))
    end
    (0 ≤ tol < Inf) || throw(ArgumentError(
        "shared_denominator_pade: tol must satisfy 0 ≤ tol < Inf; got tol=$tol. " *
        "Suggestion: pass a finite relative tolerance such as 1e-14."))

    tol_t = real(T) <: AbstractFloat ? real(T)(tol) : tol
    cnorm = norm(reduce(vcat, [collect(j) for j in jets]))
    τ = tol_t * cnorm

    # Failure mode 2 (genuinely-zero input): ‖c‖ = 0 exactly, so there is
    # nothing to approximate.  The predicate is `iszero`, NOT `cnorm ≤ tol`:
    # every other threshold in this module is relative (τ = tol·‖c‖₂, GGT 2013
    # md:213/225), so an absolute test here was the one non-scale-covariant
    # line and falsely rejected valid jets scaled by α < tol/‖c‖ (bug
    # padetaylor-ked0; SP.6 pins α ∈ {1e-20, 1e20}).  This is DISTINCT from a
    # locally-regular jet whose high-order window vanishes (that reduces to
    # Q=1 below) — here the whole jet, low-order included, is exactly zero.
    iszero(cnorm) && throw(ErrorException(
        "shared_denominator_pade: the Taylor jets are identically zero " *
        "(‖c‖ = 0), so there is nothing to approximate. " *
        "Suggestion: check the jets are non-trivial."))

    # --- Steps 2–4: build A_full, SVD, rank check + degree reduction --------
    # Vector analogue of RobustPade's reduction `while` loop (ADR-0027).
    local A_full, S, Vt, m_cur
    m_cur = Int(m)
    while true
        blocks = [_toeplitz_block(jets[i], m_cur) for i = 1:d]
        A_full = reduce(vcat, blocks)            # (d·m_cur) × (m_cur+1)
        _, S, Vt = pade_svd(A_full; full = true)
        # ρ = number of singular values above the noise threshold τ.
        ρ = count(s -> s > τ, S)

        # Accept when the matrix has rank ≥ m_cur (ADR-0027). Two accepted regimes:
        #   ρ == m_cur — exactly one σ ≤ τ: an isolated 1-D null space; the shared
        #     denominator annihilates the matching equations (genuine shared poles).
        #   ρ == m_cur+1 — full column rank (no σ ≤ τ): NO exact null exists, but
        #     the smallest-σ right vector is the *legitimate least-squares* best
        #     shared denominator (Mano–Tsuda's "unique iff rank = m" fails, yet a
        #     well-defined least-squares Q remains — the regime of an entire /
        #     locally-regular jet, e.g. bead padetaylor-unk).  Accepting it is
        #     essential: at high precision (BigFloat τ ≈ 1e-74) an entire system is
        #     full column rank at EVERY degree, so reducing-until-σ≤τ would spiral
        #     to m=0 and throw on a perfectly steppable jet.  A spurious pole at
        #     the expansion centre that this least-squares Q may carry is handled
        #     by the z^λ cancellation below, not by rejecting the Q.
        ρ ≥ m_cur && break

        # ρ < m_cur — rank-deficient → reduce toward the supported degree and
        # rebuild (the vector analogue of GGT's "reduce n to ρ", md:136).
        # ρ == 0 means the (m_cur,m_cur) matching window carries no pole
        # structure at τ: the +2 window reads the HIGH-order coefficients, which
        # vanish for a locally-regular / small-step jet — drop one degree and
        # retry rather than throw.
        m_cur = ρ == 0 ? m_cur - 1 : ρ

        if m_cur < 1
            # No shared denominator of any degree ≥ 1: the jets are locally
            # regular (polynomial-like at τ).  Return Q = 1 with the Taylor-
            # polynomial numerators — the vector analogue of RobustPade's
            # n==0 branch (RobustPade.jl:407-410); the stepper then takes the
            # plain Taylor step.  (Genuinely-zero input was rejected above.)
            numerators = Vector{Vector{T}}(undef, d)
            for i = 1:d
                a_i = Vector{T}(collect(jets[i]))
                last_a = findlast(x -> abs(x) > τ, a_i)
                numerators[i] = last_a === nothing ? T[zero(T)] : a_i[1:last_a]
            end
            return (numerators, T[one(T)])
        end
    end

    # --- Step 5: QR-reweighting (padeapprox.m lines 111–117) ----------------
    # Ported unchanged from RobustPade's `:svd` path.  `Vt[end,:]` is the
    # smallest-singular-value right null vector = denominator estimate b₀.
    b_init = Vector{T}(Vt[end, :])

    # Null-space isolation is now enforced by the `ρ == m_cur` break above
    # (exactly one σ ≤ τ ⇒ an isolated 1-D null space, Mano–Tsuda p.12 "unique
    # iff rank = m").  The old `n_near > 1` guard was dead code for d≥2 and is
    # removed (ADR-0027).

    eps_T = real(T) <: AbstractFloat ? sqrt(eps(real(T))) : sqrt(eps(Float64))
    D = Diagonal([abs(bk) + eps_T for bk in b_init])
    F = qr(adjoint(A_full * D))                  # (A_full·D)ᵀ in MATLAB
    b = D * F.Q[:, m_cur + 1]
    b ./= norm(b)

    # --- Steps 6–7: recover numerators, cancel common z^λ, normalise Q(0)=1 -
    # Recover each raw numerator a_i = upper_i · b BEFORE any cancellation, so
    # the leading z^λ factor is cancelled consistently from b and every a_i.
    raw_nums = [Vector{T}(_upper_block(jets[i], m_cur) * b) for i = 1:d]

    # Leading-z^λ cancellation, mirroring RobustPade._trim_and_normalise
    # (RobustPade.jl:461-476).  A spurious pole at the expansion centre —
    # Q(0) = b[1] ≈ 0, the generic outcome for an entire / locally-regular jet
    # whose shared-Padé denominator places a root at t=0 (ADR-0027) — shows up
    # as λ leading near-zeros SHARED by b and every numerator (a_i[1] = c₀·b[1]
    # ≈ 0).  Cancel the common z^λ factor rather than throwing: the scalar
    # `:svd` path does exactly this, and the stepper evaluates at t=1 where the
    # spurious origin pole–zero pair cancels.  Throw only if ALL of b is below
    # tol (genuine rank-0 — should have been caught by the ρ==0 mode above).
    lam_idx = findfirst(x -> abs(x) > tol_t, b)
    lam_idx === nothing && throw(ErrorException(
        "shared_denominator_pade: every denominator coefficient is below " *
        "tol after the QR step; the jets carry no recoverable shared " *
        "denominator. Suggestion: check the jets are non-trivial, or loosen tol."))
    lam = lam_idx - 1
    if lam > 0
        b        = b[(lam + 1):end]
        raw_nums = [a[(lam + 1):end] for a in raw_nums]
    end

    # Normalise to Q(0) = b[1] = 1, then trim trailing near-zeros (scalar
    # thresholds: |bₖ|≤tol off Q, |aₖ|≤τ off each numerator) so the d=1 stack
    # reduces *exactly* to robust_pade(jet,m,m;:svd).
    b1 = b[1]
    denominator = Vector{T}(b ./ b1)
    last_b = findlast(x -> abs(x) > tol_t, denominator)
    last_b === nothing && error(
        "shared_denominator_pade: denominator collapsed below tol after " *
        "normalisation; algorithm bug — file a bead with the input jets.")
    denominator = denominator[1:last_b]

    numerators = Vector{Vector{T}}(undef, d)
    for i = 1:d
        a_i = raw_nums[i] ./ b1
        last_a = findlast(x -> abs(x) > τ, a_i)
        numerators[i] = last_a === nothing ? T[zero(T)] : a_i[1:last_a]
    end

    return (numerators, denominator)
end

end # module SharedPade
