# =============================================================================
# aaa.jl — a self-contained, literate Julia port of the AAA algorithm.
#
# AAA (adaptive Antoulas–Anderson) rational approximation, from
#
#   Y. Nakatsukasa, O. Sète, L. N. Trefethen, "The AAA algorithm for
#   rational approximation", SIAM J. Sci. Comput. 40(3), 2018, A1494–A1522.
#
# Local LaTeX source (the ground truth this file is ported from, Law 1):
#   references/tex/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40/
#     p.tex        — §2 "Rational barycentric representations" (lines 155–230),
#                    §3 "Core AAA algorithm" (lines 315–474);
#     aaa_alg.m    — the reference MATLAB `aaa` and its `prz` subroutine;
#     cleanup.m    — the Froissart-doublet removal pass (§7, "removing").
#
# WHY THIS PROBE EXISTS — the orthogonal-oracle role (bead V1d)
# -----------------------------------------------------------------------------
# `PadeTaylor.SharedPade.shared_denominator_pade` recovers one shared
# denominator Q from the *Taylor jets* of d component functions, via the
# null space of a block-Toeplitz matrix.  Two oracles already cross-check
# it — Calgo 766 (FORTRAN Beckermann–Labahn recurrence) and a block-Toeplitz
# *determinant* — and both, like SharedPade itself, consume Taylor
# coefficients.  AAA is a *structurally different* third oracle:
#
#   * it fits a rational approximant from **sampled function values** on a
#     point set Z ⊂ ℂ, never touching a Taylor jet;
#   * its denominator is a **barycentric** partial fraction, not a monomial
#     polynomial;
#   * its weight update is the SVD of a **Loewner** matrix (NST eq. 3.6),
#     not of a block-Toeplitz matrix;
#   * its poles come from a **generalised arrowhead eigenvalue problem**
#     (NST eq. 3.8 / p.tex:459–472), not from rooting a coefficient vector.
#
# Every stage of the pipeline is replaced.  If AAA — sampling f_i on a
# contour in the analyticity region and fitting each component on its own —
# lands its poles on the *same* points as the roots of SharedPade's shared
# Q, that agreement is independent evidence the shared denominator is
# correct.  (AAA does use an SVD internally: that is the algorithm.  The
# independence is "samples vs jet" and "Loewner vs block-Toeplitz", not
# "no SVD".)
#
# THE BARYCENTRIC FORM (NST §2, p.tex:155–207)
# -----------------------------------------------------------------------------
# After m steps AAA represents its approximant as the quotient of two
# partial fractions over m support points z_1,…,z_m with data values
# f_j = f(z_j) and weights w_j (p.tex:158–161):
#
#       r(z) = n(z) / d(z) ,    n(z) = Σ_j w_j f_j / (z - z_j) ,
#                                d(z) = Σ_j w_j     / (z - z_j) .
#
# Multiplying through by the node polynomial ℓ(z) = Π_j (z - z_j) shows r
# is genuinely of type (m-1, m-1): the apparent type (2m-1, 2m-1) collapses
# because the ℓ(z) factor cancels (p.tex:179–204).  At each z_j with
# w_j ≠ 0 the formula is a removable ∞/∞ singularity with limit f_j, so r
# *interpolates* the data at the support points (p.tex:221–230).  This is
# why the least-squares fit below is taken over Z minus the support points.
#
# WHY GREEDY SUPPORT-POINT SELECTION (NST §3, p.tex:361–375)
# -----------------------------------------------------------------------------
# Picking support points greedily — at step m, the sample where the current
# nonlinear residual |f - r| is largest (p.tex:361–364, aaa_alg.m:24) — is
# what lets AAA avoid the exponential ill-conditioning of a fixed (e.g.
# equispaced) point set.  Each new support point is placed exactly where
# the model is worst, so the barycentric basis functions 1/(z - z_j) stay
# well separated and the Loewner matrix stays well conditioned (p.tex:
# 279–292).  The iteration stops when the residual falls below
# tol·‖f‖_∞ (default tol 1e-13, p.tex:366–369, aaa_alg.m:36).
#
# HOW POLES ARE FOUND (NST §3, p.tex:455–474; aaa_alg.m:43–50 `prz`)
# -----------------------------------------------------------------------------
# The poles of r are (generically) the zeros of the barycentric denominator
# d(z) = Σ_j w_j/(z - z_j).  NST compute them as the finite eigenvalues of
# an (m+1)×(m+1) generalised eigenvalue problem in *arrowhead* form
# (p.tex:459–472):
#
#       ⎡ 0   w_1  w_2  ⋯  w_m ⎤        ⎡ 0          ⎤
#       ⎢ 1   z_1              ⎥        ⎢   1        ⎥
#       ⎢ 1        z_2         ⎥ = λ    ⎢     1      ⎥
#       ⎢ ⋮            ⋱      ⎥        ⎢       ⋱    ⎥
#       ⎣ 1                z_m ⎦        ⎣          1 ⎦
#
# The right pencil B = diag(0,1,…,1) is singular, so two eigenvalues are
# infinite (p.tex:473–474); the remaining m-1 finite eigenvalues are the
# zeros of d, i.e. the poles of r.  Zeros of r come from the identical
# pencil with w_j replaced by w_j f_j (aaa_alg.m:49).
#
# SPURIOUS-POLE FILTERING (NST §7 "removing", cleanup.m)
# -----------------------------------------------------------------------------
# A genuine pole of f carries a residue of order unity; a numerical
# Froissart doublet is a pole whose residue is ~1e-13 or smaller — a
# pole-zero pair so close it contributes nothing (p.tex:652–661,
# cleanup.m:3).  NST's standard recipe drops such doublets by removing the
# support point nearest each tiny-residue pole and re-solving the
# least-squares problem once.  We expose that as an explicit, optional
# `cleanup` flag so the verification can report both the raw and the
# cleaned pole sets.
#
# This file deliberately adds NO dependency to the project Project.toml
# (CLAUDE.md, probe constraint): it uses only LinearAlgebra, which the
# project already depends on.  AAA is a test-only oracle.
# =============================================================================

module AAAOracle

using LinearAlgebra: svd, eigen, Diagonal, norm, I

export aaa, AAAResult

"""
    AAAResult

The output of [`aaa`](@ref): the barycentric data plus the poles, residues
and zeros of the recovered rational approximant.

Fields (all match the NST 2018 notation, p.tex §3):

  - `support_points :: Vector{Complex}` — the greedily selected z_j.
  - `values :: Vector{Complex}`         — the data f_j = f(z_j).
  - `weights :: Vector{Complex}`        — the barycentric weights w_j.
  - `poles :: Vector{Complex}`          — zeros of d(z) (eq. 3.8 arrowhead).
  - `residues :: Vector{Complex}`       — residue of r at each pole.
  - `zeros :: Vector{Complex}`          — zeros of r (zeros of n(z)).
  - `errvec :: Vector{Float64}`         — max sample error at each step.
"""
struct AAAResult{T<:Complex}
    support_points::Vector{T}
    values::Vector{T}
    weights::Vector{T}
    poles::Vector{T}
    residues::Vector{T}
    zeros::Vector{T}
    errvec::Vector{Float64}
end

# -----------------------------------------------------------------------------
# Barycentric evaluation (NST eq. 2.1, p.tex:158–161; aaa_alg.m `rhandle`)
# -----------------------------------------------------------------------------

"""
    bary_eval(zz, z, f, w) -> value(s)

Evaluate the barycentric rational `r(zz) = n(zz)/d(zz)` from NST eq. (2.1)
at one point or a vector of points.  At a support point the formula is the
removable ∞/∞ singularity whose value is the data value `f_j`
(p.tex:221–230); we detect the resulting `NaN` and substitute `f_j`
exactly, as the reference `rhandle` does (aaa_alg.m:56–59).
"""
function bary_eval(zz::Number, z::AbstractVector, f::AbstractVector,
                   w::AbstractVector)
    # Exact hit on a support point: return the interpolated data value.
    for j in eachindex(z)
        if zz == z[j]
            return complex(f[j])
        end
    end
    num = zero(complex(float(zz)))
    den = zero(complex(float(zz)))
    for j in eachindex(z)
        t = w[j] / (zz - z[j])
        num += t * f[j]
        den += t
    end
    return num / den
end

bary_eval(zz::AbstractVector, z, f, w) = [bary_eval(p, z, f, w) for p in zz]

# -----------------------------------------------------------------------------
# Pole / zero / residue extraction (NST eq. 3.8, p.tex:455–474; `prz`)
# -----------------------------------------------------------------------------

"""
    _arrowhead_roots(z, c) -> Vector{Complex}

The zeros of the barycentric partial fraction `Σ_j c_j/(ζ - z_j)`, computed
as the finite eigenvalues of the (m+1)×(m+1) arrowhead generalised
eigenvalue problem `E v = λ B v` of NST eq. (3.8) (p.tex:459–472):

  - `E` has a zero in the (1,1) corner, the coefficients `c` along the
    first row, ones down the first column, and `diag(z)` in the trailing
    m×m block;
  - `B = diag(0, 1, …, 1)` is singular, so two eigenvalues are infinite
    (p.tex:473–474).

The finite eigenvalues are returned.  With `c = w` these are the poles of
`r`; with `c = w .* f` they are the zeros of `r` (aaa_alg.m:44–50).
"""
function _arrowhead_roots(z::AbstractVector{T}, c::AbstractVector{T}) where {T}
    m = length(z)
    E = zeros(T, m + 1, m + 1)
    B = zeros(T, m + 1, m + 1)
    @inbounds for j = 1:m
        E[1, j + 1] = c[j]          # first row: the partial-fraction weights
        E[j + 1, 1] = one(T)        # first column: ones
        E[j + 1, j + 1] = z[j]      # trailing block: diag(z)
        B[j + 1, j + 1] = one(T)    # B = diag(0, 1, …, 1) — singular
    end
    λ = eigen(E, B).values
    # Drop the (generically two) infinite eigenvalues of the singular pencil.
    finite_thresh = 1e13 * (one(real(T)) + maximum(abs, z))
    return T[v for v in λ if isfinite(v) && abs(v) < finite_thresh]
end

"""
    _poles_residues_zeros(z, f, w) -> (poles, residues, zeros)

NST's `prz` (aaa_alg.m:43–50): poles = zeros of `d`, zeros = zeros of `n`,
and the residue of `r` at each pole `p` is estimated by the four-point
contour quadrature `r(p + dz)·dzᵀ/4` over `dz = 1e-5·exp(2πi k/4)`,
k = 1..4 — a tiny circle around the pole (aaa_alg.m:47–48).
"""
function _poles_residues_zeros(z::AbstractVector{T}, f::AbstractVector{T},
                               w::AbstractVector{T}) where {T}
    poles = _arrowhead_roots(z, w)
    zers  = _arrowhead_roots(z, w .* f)
    dz = T[1e-5 * exp(2im * pi * k / 4) for k = 1:4]
    residues = T[sum(bary_eval(p .+ dz, z, f, w) .* dz) / 4 for p in poles]
    return poles, residues, zers
end

# -----------------------------------------------------------------------------
# The core AAA iteration (NST §3, p.tex:315–451; aaa_alg.m:15–41)
# -----------------------------------------------------------------------------

"""
    aaa(Z, F; tol=1e-13, mmax=100, cleanup=true)
        -> AAAResult

Compute a rational approximation of the data `F` on the sample set `Z` by
the AAA algorithm of Nakatsukasa–Sète–Trefethen 2018 (§3, p.tex:315–474).

`Z` and `F` are equal-length vectors: `F[k]` is the value of the function
being fitted at sample point `Z[k]`.  The sample points must lie in the
analyticity region of the function — a sample *at* a pole would be `Inf`
(see the probe README for the contour choice used in `verify.jl`).

Keyword arguments mirror the reference `aaa` (aaa_alg.m:16–17):

  - `tol`     — relative stopping tolerance; the iteration halts once the
                max residual over the samples drops below `tol·‖F‖_∞`
                (p.tex:366–369).  Default `1e-13`.
  - `mmax`    — maximum number of support points (max type `(mmax-1,
                mmax-1)`).  Default `100`.
  - `cleanup` — if `true`, run NST's Froissart-doublet removal pass once
                (§7, cleanup.m): poles with `|residue| < 1e-13` and their
                nearest support points are dropped and the weights
                re-solved.  Default `true`.

Returns an [`AAAResult`](@ref) carrying the barycentric data
`(support_points, values, weights)` and the recovered `poles`, `residues`
and `zeros`.

Throws (CLAUDE.md Rule 1 — fail loud, no silent NaN):

  - `ArgumentError` if `Z` and `F` differ in length, if fewer than two
    samples are supplied, or if any sample value is non-finite (a sample
    sat on a pole).
  - `ErrorException` if the Loewner subset would have fewer than `m` rows
    (`m ≤ M/2` is required for the least-squares fit, p.tex:358–359).
"""
function aaa(Z::AbstractVector, F::AbstractVector;
             tol::Real = 1e-13, mmax::Integer = 100, cleanup::Bool = true)
    length(Z) == length(F) || throw(ArgumentError(
        "aaa: sample-point vector Z and value vector F must have equal " *
        "length; got length(Z)=$(length(Z)), length(F)=$(length(F)). " *
        "Suggestion: pass F = f.(Z)."))
    M = length(Z)
    M ≥ 2 || throw(ArgumentError(
        "aaa: need at least 2 sample points to fit a rational; got M=$M. " *
        "Suggestion: sample the function on a denser contour."))

    Zc = collect(complex(float.(Z)))
    Fc = collect(complex(float.(F)))
    all(isfinite, Fc) || throw(ArgumentError(
        "aaa: a sample value is non-finite (Inf/NaN) — a sample point " *
        "almost certainly sits on a pole of the function. " *
        "Suggestion: choose a sampling contour inside the analyticity " *
        "region, away from the poles."))

    Fnorm = norm(Fc, Inf)

    # `active` indexes the samples NOT yet promoted to support points; these
    # are the rows Z^(m) of the Loewner least-squares problem (p.tex:328).
    active = collect(1:M)
    sup_idx = Int[]                     # indices (into Z) of support points
    z = ComplexF64[]                    # support points  z_j
    f = ComplexF64[]                    # data values     f_j
    w = ComplexF64[]                    # barycentric weights w_j
    errvec = Float64[]

    # Step-0 model: the constant equal to the mean of the data (aaa_alg.m:22).
    R = fill(complex(sum(Fc) / M), M)

    mmax_eff = min(Int(mmax), M ÷ 2)
    mmax_eff ≥ 1 || throw(ErrorException(
        "aaa: M=$M samples admit no Loewner fit (need m ≤ M/2 with m ≥ 1). " *
        "Suggestion: supply at least 2 sample points."))

    for m = 1:mmax_eff
        # --- greedy support-point selection (p.tex:361–364, aaa_alg.m:24) ---
        # Promote the still-active sample where the residual |F - R| is worst.
        jrel = argmax(abs.(Fc[active] .- R[active]))
        j = active[jrel]
        push!(sup_idx, j)
        push!(z, Zc[j])
        push!(f, Fc[j])
        deleteat!(active, jrel)

        # --- Loewner matrix on the active subset (NST eq. 3.6, p.tex:401) --
        # A_full[i,k] = (F_i - f_k) / (Z_i - z_k); built via A = S_F C - C S_f
        # (eq. 3.7, p.tex:438) for fidelity to the reference.
        Zr = Zc[active]
        Fr = Fc[active]
        C = [1.0 / (Zr[i] - z[k]) for i in eachindex(Zr), k = 1:m]
        A = Diagonal(Fr) * C - C * Diagonal(f)

        # --- weight update: SVD of the Loewner matrix (p.tex:413–415) ------
        # w is the right singular vector of the smallest singular value.
        Vt = svd(A).V
        w = Vt[:, end]

        # --- barycentric model on the active subset (p.tex:444–451) --------
        N = C * (w .* f)                # numerator values  n(Z^(m))
        D = C * w                       # denominator values d(Z^(m))
        R = copy(Fc)                    # r interpolates f at support points
        R[active] = N ./ D

        err = norm(Fc .- R, Inf)
        push!(errvec, err)
        err ≤ tol * Fnorm && break      # converged (p.tex:366–369)
    end

    # --- pole / zero / residue extraction (NST eq. 3.8 + prz) --------------
    poles, residues, zers = _poles_residues_zeros(z, f, w)

    # --- optional Froissart-doublet removal (NST §7, cleanup.m) ------------
    if cleanup && !isempty(poles)
        spurious = findall(p -> abs(p) < 1e-13, residues)
        if !isempty(spurious)
            # Drop the support point nearest each tiny-residue pole, then
            # re-solve the least-squares problem once (cleanup.m:6–20).
            drop = Int[]
            for s in spurious
                d2 = abs.(z .- poles[s])
                push!(drop, argmin(d2))
            end
            keep = setdiff(1:length(z), unique(drop))
            z = z[keep]
            f = f[keep]
            m2 = length(z)
            if m2 ≥ 1
                # Re-form the Loewner matrix over Z minus the surviving
                # support points and re-solve for the weights.
                kept_sup = sup_idx[keep]
                active2 = setdiff(1:M, kept_sup)
                Zr = Zc[active2]
                Fr = Fc[active2]
                C = [1.0 / (Zr[i] - z[k]) for i in eachindex(Zr), k = 1:m2]
                A = Diagonal(Fr) * C - C * Diagonal(f)
                w = svd(A).V[:, end]
                poles, residues, zers = _poles_residues_zeros(z, f, w)
            end
        end
    end

    return AAAResult{ComplexF64}(z, f, w, poles, residues, zers, errvec)
end

end # module AAAOracle
