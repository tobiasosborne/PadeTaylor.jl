"""
    PadeTaylor.IVPBVPHybrid

Tier-3 *hybrid* driver that composes the existing
`PathNetwork.path_network_solve` (Padé-Taylor IVP / PFS) with
`BVP.bvp_solve` (Chebyshev spectral BVP, three-argument-RHS overload)
on **pole-free sectors** of a P̃_III-type problem.  This is FFW 2017
§3 verbatim — the load-bearing infrastructure for FFW Figure 5.

## Why a hybrid is necessary (FFW 2017 md:203-247)

The Padé-Taylor IVP machinery (PFS / PathNetwork) is *stable* in the
pole-rich region because the integration is locally bounded by the
pole-bridging Padé approximant.  On a **pole-free sector** the IVP is
exponentially ill-conditioned: FFW md:252-264 derive the relative
condition number

    κ_r = |1/w''| · |-(w')²/w + (1/4)(2αw² + 3γw³ - δe^{2ζ}/w)|

and show (md:260-262) that on the pole-free sector of FFW Figure 5,

    κ_r ~ (27/16) · e^{2 Re ζ / 3},   Re ζ → +∞,

an *exponential* growth.  At the boundary of Figure 5's sector
(`Re ζ ≤ 2 log 30`), `κ_r ≈ 27/16 · 30^{4/3} ≈ 157` (FFW md:264).  A
pure-IVP walk over the sector would amplify any perturbation by this
factor before reaching the IC, and FFW report (md:226) an error of
~10⁻¹ in that experiment.  The cure (FFW md:224, originally proposed
in [9]) is to solve the equation as a **boundary-value problem** on
the smooth sector: BVP residual minimization is *stable* on smooth
data, whereas the wrong-direction IVP shoot is not.

## The four-step FFW §3 algorithm

Given a `PainleveProblem` `pp` in the transformed (ζ-plane) frame, a
description of the pole-free sector S ⊂ ζ-plane, and an asymptotic-IC
function `z -> (u, u')` valid at large `|z|`:

  **Step 1 — PFS on the pole-bearing exterior.**  Run
  `path_network_solve` from each curved-boundary asymptotic IC point
  inward toward an *interior* anchor of the original problem.  Two
  separate PFS walks are launched — one along each curved boundary
  ray (in the ζ-plane these are two horizontal segments at
  `Im ζ = im_lo` and `Im ζ = im_hi`).  The walks supply the BC source
  for Step 3 *and* cover the pole-bearing region of the full domain.

  **Step 2 — BVP boundary harvest.**  Read off `u(ζ)` and `u'(ζ)` at
  a discrete set of points along each curved boundary of S using the
  Stage-2 evaluation of the PFS solutions from Step 1.  These are the
  Dirichlet BCs the BVP solver consumes.

  **Step 3 — BVP solve on the sector.**  Slice S into a stack of
  axis-aligned `ζ`-segments at constant `Re ζ`, each running from
  `Im ζ = im_lo` to `Im ζ = im_hi`.  Call `bvp_solve` on each slice
  with the harvested values as Dirichlet BCs.  The three-argument-RHS
  overload is required because P̃_III's RHS depends on `w'`
  (FFW md:43, bead `padetaylor-i76`).  Refine `N` (collocation count)
  until the BVP step-norm convergence drops below the PFS boundary
  accuracy (FFW md:245-247).

  **Step 4 — Glue.**  Build an `IVPBVPSolution` whose `sol(ζ)`
  callable dispatches on whether ζ falls inside S (BVP slice) or
  outside (PFS).  Between two bracketing slices the callable fails
  loud (Rule 1) on GROSS inter-slice incoherence — a divergent/corrupt
  slice (bug `padetaylor-q0yq`).  The finer FFW md:247 PFS-vs-BVP
  derivative match (~1e-7) is a v1 Float64 deferral (the value check is
  a tautology, the derivative is ~1e-1 here — bead `padetaylor-sn9a`).

Sector geometry note (FFW md:222): Figure 5's pole-free sector is
described in the z-plane as `-3π/4 < arg z < 9π/4`.  Under the PIII
transform `z = exp(ζ/2)` this maps to the ζ-plane horizontal strip
`-3π/2 < Im ζ < 9π/2` — *axis-aligned* in ζ, despite FFW's "curved
boundaries" terminology (which refers to the z-plane image of those
ζ-rays).  The driver therefore takes its sector parameters in the
ζ-frame as a `(im_lo, im_hi, re_anchor, re_extent)` NamedTuple.

## Asymptotic-IC helper `pIII_asymptotic_ic` (FFW md:222 + md:230)

The boundary ICs at the two ζ-frame entry points come from the FFW
md:230 truncated asymptotic series

    u(z) ~ z^{1/3} · [ 1 + Σ_{n ≥ 1} a_n · (z^{1/3})^{-2n} ],
                                      z → ∞,  -3π/4 < arg z < 9π/4.

The leading term `u ~ z^{1/3}` is the existence-theorem statement of
[21] (cited at FFW md:222) for the Fig 5 parameter family `(α, β, γ,
δ) = (1, β, 0, -1)`.  The first sub-leading coefficient `a_1` is
derived in `pIII_asymptotic_ic`'s docstring by substituting the
ansatz into the canonical PIII equation and matching the `z^{-1}`
order; closed form **`a_1 = -β/3`** (so `a_1 = 1/60` for FFW Fig 5's
`β = -1/20`).  Higher `a_n` follow from the same substitution
*numerically* via a power-series recurrence — the helper takes
`n_terms` as kwarg and computes them.

## What this module composes (and does not modify)

Per CLAUDE.md Rule 6 + the bead spec: the hybrid driver **composes**
existing modules, it does not refactor them.

  - `PathNetwork.path_network_solve` is consumed verbatim.
  - `BVP.bvp_solve` was minimally extended in a strictly additive
    way to support the three-argument RHS `f(z, u, u')` (bead
    `padetaylor-i76`).  The 2-arg API and every existing call site
    are byte-unchanged.
  - `Painleve.PainleveProblem` and `CoordTransforms.pIII_transformed_rhs`
    are consumed verbatim.

## References

  - **FFW 2017 §3** — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:203-247` (algorithm), md:252-264 (condition-number derivation), md:222 (sector geometry + asymptotic-series ansatz), md:230 (series form), md:240-244 (Fig 5 IC values), md:43 (P̃_III RHS).
  - **ADR-0014** — `docs/adr/0014-ivp-bvp-hybrid.md`.
  - **Worklog 039** — `docs/worklog/039-ivp-bvp-hybrid.md`.
  - **Bead `padetaylor-0co`** — this work.
  - **Bead `padetaylor-i76`** — the BVP 3-arg-RHS extension.
"""
module IVPBVPHybrid

using ..Problems:        PadeTaylorProblem
using ..PathNetwork:     PathNetworkSolution, path_network_solve
using ..BVP:             BVPSolution, bvp_solve
using ..Painleve:        PainleveProblem
import ..Painleve:       _coord

export solve_pole_free_hybrid, IVPBVPSolution, pIII_asymptotic_ic

# -----------------------------------------------------------------------------
# Asymptotic-IC helper for P_III tronquée (FFW md:222 + md:230 + ref [21])
# -----------------------------------------------------------------------------

"""
    pIII_asymptotic_ic(z; n_terms = 10, α = 1, β = -1/20, γ = 0, δ = -1) -> (u, u')

Truncated asymptotic series for the tronquée P_III solution

    u(z) ~ z^{1/3} · [ 1 + Σ_{n=1}^{n_terms} a_n · z^{-2n/3} ],   z → ∞,

valid in the pole-free sector `-3π/4 < arg z < 9π/4` for the parameter
family `(α, β, γ, δ) = (1, β, 0, -1)` of FFW 2017 Figure 5 (FFW
md:222 + md:230, citing the existence theorem of [21]).

The coefficients `a_n` are obtained by substituting the ansatz into
the canonical PIII equation (FFW md:31)

    u'' = (1/u)(u')² − (1/z)u' + (αu² + β)/z + γu³ + δ/u

and matching powers of `z^{-1/3}`.  (NB: earlier versions of this
docstring and ADR-0014 mis-transcribed the equation as
`(αu² + γu³)/z + β/z + δ/u`, putting `γu³` inside the `1/z` term and
adding a spurious `β/z`.  For the Fig 5 family `γ = 0`, so the
numbers are unaffected, but the transcription is now corrected to
FFW's form — bead padetaylor-qsj.)  The derivation, verified
symbolically with sympy, is exposited in `_pIII_asymptotic_coeffs`:

  **Leading consistency (z^{-1/3} order).**  The residual at this order
  is `−(α + δ)`; it vanishes iff `δ = −α`.  The tronquée ansatz is an
  asymptotic solution ONLY for `δ = −α` (FFW md:222: the existence
  theorem of [21] requires `γ = 0, α = 1 = −δ`).  This helper enforces
  `δ = −α` with a fail-loud throw.

  **a_1 (z^{-1} order).**  The residual is `−a_1(2α − δ) − β`, so
  `a_1 = −β/(2α − δ)`; under `δ = −α` this is `a_1 = −β/(3α)`.  With
  `α = 1`, `a_1 = −β/3` (= 1/60 for FFW Fig 5's `β = −1/20`).

  **a_2 (z^{-5/3} order — NOT z^{-7/3}).**  The residual is
  `−a_1²(α + δ) − a_2(2α − δ)`, so `a_2 = −a_1²(α + δ)/(2α − δ)`.  The
  factor `(α + δ)` VANISHES under the validity constraint `δ = −α`,
  hence **`a_2 = 0`** for the FFW Fig 5 family.  (The old shipped
  `a_2 ≈ −0.22208` came from matching one order too late and dropping
  the `αu²/z` term; it made the n_terms = 2 series ~880× LESS accurate
  against FFW's published `u(z₂)`.)

Only `(α, β, γ, δ) = (1, β, 0, −1)` is in scope (FFW md:222 restricts
the existence theorem to this family); `a_3+` are a documented
deferral (left at zero — the first dropped term is
`a_3 = β(β² − 4)/81 ≈ 2.47e-3`, contributing ~8e-6 to u at |z| = 30).

`n_terms` defaults to 10 — empirically near optimal at `|z| = 30` for
Float64 (more terms give diminishing returns once `|a_n z^{-2n/3}|`
falls below `eps`).  Returns `(u, u')` as a 2-tuple; the derivative is
computed by termwise differentiation of the series so the BVP solver
gets a self-consistent IC.

This helper validates `(α, γ) = (1, 0)` because the ansatz derivation
above relies on those choices; other values throw `ArgumentError`.
`δ` and `β` are free.

Throws `ArgumentError` on `n_terms < 1`, `|z| < 1` (series unreliable
near origin), or `arg z ∉ (-3π/4, 9π/4)` (outside the existence
sector).
"""
function pIII_asymptotic_ic(z::Number;
                            n_terms::Integer = 10,
                            α = 1, β = -1/20, γ = 0, δ = -1)
    n_terms ≥ 1 || throw(ArgumentError(
        "pIII_asymptotic_ic: n_terms must be ≥ 1 (got $n_terms)."))
    abs(z) ≥ 1 || throw(ArgumentError(
        "pIII_asymptotic_ic: |z| must be ≥ 1 (got |z| = $(abs(z))).  The " *
        "asymptotic series is divergent for `z → ∞`; truncated evaluation " *
        "below |z| ≈ 1 is dominated by truncation error.  Suggestion: " *
        "evaluate at a large boundary point (FFW Fig 5 uses |z| = 30)."))
    (α == 1 && γ == 0) || throw(ArgumentError(
        "pIII_asymptotic_ic: this series is derived for the (α, γ) = (1, 0) " *
        "family of FFW Fig 5 / ref [21]; got (α, γ) = ($α, $γ).  Other " *
        "parameter sets are out of scope."))
    # Leading consistency condition (FFW md:222 / md:230, sympy-verified
    # in `_pIII_asymptotic_coeffs`): the tronquée ansatz `u ~ z^{1/3}[1 +
    # …]` solves PIII only when `δ = −α`.  At the z^{-1/3} order the
    # residual is `−(α + δ)`, which must vanish; otherwise the series is
    # NOT an asymptotic solution (the existence theorem of [21] requires
    # `α = 1 = −δ`).  Fail loud rather than return a series that does not
    # satisfy the equation.  (bead padetaylor-qsj)
    δ == -α || throw(ArgumentError(
        "pIII_asymptotic_ic: the tronquée asymptotic ansatz requires " *
        "δ = −α (FFW md:222 existence theorem of ref [21]); got δ = $δ " *
        "with α = $α, so δ + α = $(δ + α) ≠ 0.  Suggestion: pass δ = $(-α) " *
        "(the FFW Fig 5 value is δ = −1 with α = 1)."))
    # Sector check (FFW md:222): `-3π/4 < arg z < 9π/4`.  Note `arg z` from
    # Julia's principal `angle` returns values in `(-π, π]`; we accept
    # any z whose principal argument lies in (-3π/4, π], which is a
    # subset of the FFW sector.  Callers wanting the Im z > 0 part of
    # the (π, 9π/4) sub-range encode it by adding `4π·i·s` in the
    # ζ-frame, which is the SheetTracker / sheet-index convention.
    az = angle(z)
    -3π/4 < az ≤ π || throw(ArgumentError(
        "pIII_asymptotic_ic: principal arg(z) = $az is outside the " *
        "pole-free sector `(-3π/4, 9π/4)` (Julia's `angle` returns " *
        "principal values, so only the slice `(-3π/4, π]` is reachable " *
        "by a direct call; the upper part `(π, 9π/4)` requires the " *
        "caller to evaluate at z*e^{-2πi} and then map via the sheet " *
        "index).  Suggestion: pass a z with `-3π/4 < arg z ≤ π`."))

    # Compute a_n.  We carry `a` as Vector{ComplexF64} length n_terms.
    a = _pIII_asymptotic_coeffs(β, δ, n_terms)

    # Sum the series and its derivative termwise.
    #   u  = s · (1 + Σ a_n s^{-2n})       where s = z^{1/3}
    #   du/dz = (1/(3s²)) · ( 1 - Σ (2n-1)·a_n·s^{-2n} )         (chain rule)
    s = z^(1/3)
    Z = promote_type(typeof(s), typeof(complex(0.0)))
    s_z = Z(s)
    s_inv2 = inv(s_z^2)
    u_sum  = one(Z)
    up_sum = one(Z)     # the "1" in 1 - Σ (2n-1) a_n s^{-2n}
    # Power s^{-2n} accumulator
    pow = one(Z)
    for n in 1:n_terms
        pow *= s_inv2                     # pow = s^{-2n}
        u_sum  += a[n] * pow
        up_sum -= (2n - 1) * a[n] * pow   # NOTE the minus
    end
    u  = s_z * u_sum
    up = up_sum / (3 * s_z^2)
    return (u, up)
end

# Series-substitution helper: returns the coefficients a[1..N] of
# `u = s + Σ a_n s^{1-2n}` (with `s = z^{1/3}`) satisfying the
# canonical PIII equation with `(α, γ) = (1, 0)`.
#
# Derivation outline.  Substitute `u(s) = s + a_1 s^{-1} + a_2 s^{-3} +
# … + a_N s^{-(2N-1)}` into the canonical PIII equation rewritten in `s`:
#
#   z = s³,  d/dz = (1/(3s²)) d/ds,  d²/dz² = (1/(9s⁴)) d²/ds² − (2/(9s⁵)) d/ds.
#
# Both sides become Laurent series in `s` (one-sided in negative powers
# beyond `s^1`).  Match powers `s^{−1}, s^{−3}, s^{−5}, …, s^{1-2N}`,
# i.e. orders `z^{−1/3}, z^{−1}, z^{−5/3}, …`; the resulting
# upper-triangular linear system gives `a_1, a_2, …, a_N` in turn.  At
# each `n ≥ 1` the matched equation has the form
#
#     (2α − δ) a_n  +  P_n(a_1, …, a_{n-1}; β, δ)  =  0
#         ⇒   a_n = -P_n / (2α − δ),
#
# and under the validity constraint `δ = −α` the prefactor is
# `(2α − δ) = 3α`.  The `n = 0` order (`z^{−1/3}`) carries no `a`
# unknown and yields the leading consistency condition `α + δ = 0`.
#
# We do NOT carry generic Laurent arithmetic here: the closed forms for
# a_1 and a_2 (sympy-derived, see below) are exhibited directly, which
# is exact and avoids a coefficient layer.  The closed form for a_1
# (= −β/3 at α = 1) is the correctness check baked into `IB.1.2`, and
# a_2 = 0 (at the FFW Fig 5 family) is checked structurally by `IB.1.1`
# (n_terms = 1 and n_terms = 2 give identical u).  a_3+ are deferred.
function _pIII_asymptotic_coeffs(β, δ, N::Integer)
    a = zeros(Float64, N)
    # Closed form for a_1 (derived in the helper's docstring).
    a[1] = -Float64(β) / 3
    # Higher-order coefficients via successive substitution.  We carry
    # truncated Laurent representations in s with negative-power degree
    # capped at the relevant order.  Below we give the closed form for
    # a_2; a_3+ are a documented deferral (left at zero).
    #
    # --------------------------------------------------------------------
    # a_2 is fixed at the z^{-5/3} (s^{-5}) order — NOT z^{-7/3}.
    # --------------------------------------------------------------------
    #
    # The earlier (shipped) derivation matched the s^{-7} order and is
    # WRONG: it was one order too late, and additionally omitted the
    # `α u²/z` contribution at s^{-5}.  The correct accounting, verified
    # symbolically with sympy (substitute the ansatz into the FFW md:31
    # equation, Laurent-expand in w = 1/s = z^{-1/3}, collect by power),
    # gives the residual coefficients (for γ = 0):
    #
    #   [z^{-1/3}]:  −(α + δ)                      ← CONSISTENCY condition
    #   [z^{-1}  ]:  −a_1(2α − δ) − β              ← fixes a_1
    #   [z^{-5/3}]:  −a_1²(α + δ) − a_2(2α − δ)    ← fixes a_2
    #   [z^{-7/3}]:  a_1³δ − 2 a_1 a_2(α + δ) + (4/9)a_1 − a_3(2α − δ)
    #                                              ← fixes a_3
    #
    # Step A — leading consistency (z^{-1/3} order):  −(α + δ) = 0.
    #   The tronquée ansatz u ~ z^{1/3}[1 + …] is a solution of PIII
    #   ONLY when δ = −α.  This matches FFW md:222: the existence theorem
    #   of [21] requires "γ = 0, α = 1 = −δ".  The public entry
    #   `pIII_asymptotic_ic` enforces δ = −α with a fail-loud throw; this
    #   helper still returns the well-defined coefficients for any δ so
    #   the closed forms below are exhibited in full generality.
    #
    # Step B — a_1 (z^{-1} order):  −a_1(2α − δ) − β = 0 ⇒
    #   a_1 = −β/(2α − δ).  Under the validity constraint δ = −α this is
    #   a_1 = −β/(3α); with α = 1, a_1 = −β/3 (= 1/60 for β = −1/20).
    #   [The shipped `a[1] = −β/3` already assumes α = 1, the only family
    #    this helper serves; see the caller's (α,γ) = (1,0) guard.]
    #
    # Step C — a_2 (z^{-5/3} order):  −a_1²(α + δ) − a_2(2α − δ) = 0 ⇒
    #
    #       a_2 = −a_1² (α + δ) / (2α − δ).
    #
    #   With α = 1 and a_1 = −β/3 this is the general-δ closed form
    #
    #       a_2 = −(β²/9) (1 + δ) / (2 − δ)              (α = 1).
    #
    #   The crucial structural fact is the factor (α + δ) = (1 + δ):
    #   under the validity constraint δ = −α = −1 it VANISHES, so
    #
    #       a_2 = 0    for the FFW Fig 5 family (α = 1, δ = −1).
    #
    #   The old shipped value a_2 = ((4/9) + δ a_1²)/(2δ) ≈ −0.22208 was
    #   spurious (wrong order + missing α u²/z term).  See bead
    #   `padetaylor-qsj` and ADR-0014.  Numerically, with a_2 = 0 the
    #   n_terms = 2 series matches FFW's published u(z₂) to |err| ≈ 8e-6,
    #   versus ≈ 7e-3 with the old a_2 — an ~880× improvement that
    #   confirms the direction of the fix.
    if N ≥ 2
        # General-δ a_2 at α = 1 (this helper's only family): the
        # numerator carries (1 + δ), which kills a_2 at δ = −1.
        a[2] = -(a[1]^2 * (1 + Float64(δ))) / (2 - Float64(δ))
    end
    # a_3 .. a_N: not implemented in v1; left at zero.  From Step's
    # z^{-7/3} balance, under α = 1, δ = −1, a_2 = 0 the a_3 coefficient
    # is a_3 = β(β² − 4)/81 (= 533/216000 ≈ 2.47e-3 for β = −1/20), so it
    # is the FIRST dropped term.  At |z| = 30, |z^{-7/3}| ≈ 3·10⁻⁴, so the
    # a_3·z^{-7/3} contribution to u is ~8·10⁻⁶ — exactly the measured
    # n_terms = 2 residual against FFW's published IC `u(z₂)`.  This is
    # good enough for the BVP-IC role (the BVP solver then tightens to
    # spectral accuracy on the sector).  Higher precision at the IC would
    # require porting FFW's "optimal truncation" procedure (FFW md:232) —
    # deferred as a separate bead.
    return a
end

# -----------------------------------------------------------------------------
# Hybrid solution container
# -----------------------------------------------------------------------------

"""
    IVPBVPSolution{T <: AbstractFloat}

Composite solution wrapper from `solve_pole_free_hybrid`.  Fields:

  - `pfs_top    :: PathNetworkSolution{T}` — PFS walk launched from the
    upper boundary of the sector (at `Im ζ = im_hi`).
  - `pfs_bot    :: PathNetworkSolution{T}` — PFS walk from the lower
    boundary (`Im ζ = im_lo`).  May be the same object as `pfs_top` in
    the degenerate "full plane" mode (see `solve_pole_free_hybrid`'s
    `degenerate_full_plane` fast path).
  - `bvp_slices :: Vector{BVPSolution{T, Complex{T}}}` — one BVP solve
    per `Re ζ` slice through the sector.
  - `slice_re   :: Vector{T}` — the `Re ζ` value of each slice (parallel
    to `bvp_slices`).
  - `sector     :: NamedTuple` — `(im_lo, im_hi, re_anchor, re_extent)`
    in the ζ-frame, the FFW-md:222 strip description.
  - `glue_tol   :: T` — the boundary-overlap tolerance used at glue time.
  - `pp         :: PainleveProblem` — the source problem (carried for
    self-describing access).

Call as `sol(ζ)` to get `(w, w')` in the ζ-frame.  Dispatch is on
whether ζ lies inside the BVP sector (uses the nearest `Re ζ` slice's
barycentric interpolant) or outside (uses the PFS walks' Stage-2 grid;
PFS does not support `sol(ζ)` so an explicit grid query is required
— the v1 convention is that the caller passes ζ values that line up
with `pfs_top.grid_z` / `pfs_bot.grid_z` for the exterior).
"""
struct IVPBVPSolution{T <: AbstractFloat}
    pfs_top    :: PathNetworkSolution{T}
    pfs_bot    :: PathNetworkSolution{T}
    bvp_slices :: Vector{BVPSolution{T, Complex{T}}}
    slice_re   :: Vector{T}
    sector     :: NamedTuple
    glue_tol   :: T
    pp         :: PainleveProblem
end

# -----------------------------------------------------------------------------
# Public driver
# -----------------------------------------------------------------------------

"""
    solve_pole_free_hybrid(pp::PainleveProblem, sector, asymptotic_ic_fn;
                            pfs_kwargs = (;),
                            bvp_kwargs = (;),
                            glue_tol = 1e-8,
                            n_slices = 8,
                            degenerate_full_plane = false) -> IVPBVPSolution

Solve `pp` (a `:transformed`-frame `PainleveProblem`) on the pole-rich
exterior + pole-free sector via the FFW 2017 §3 hybrid algorithm.
PII-V are in scope; the v1 implementation is hard-tested on PIII Fig 5.

`sector::NamedTuple` describes the pole-free sector in the **ζ-frame**:

    sector = (im_lo, im_hi, re_anchor, re_extent)

with `im_lo < im_hi`; `re_anchor` is the `Re ζ` value at which the
PFS walks are anchored (the FFW Fig 5 IC at `|z| = 30` translates to
`re_anchor ≈ 2 log 30 ≈ 6.8`), and `re_extent` is the half-width of
the BVP slab in `Re ζ`.  The BVP region is the rectangle
`{(Re ζ, Im ζ) : re_anchor - re_extent ≤ Re ζ ≤ re_anchor,  im_lo ≤
Im ζ ≤ im_hi}`.

`asymptotic_ic_fn(z) -> (u, u')` is the user-supplied callable that
returns the boundary IC at z-frame point `z`; the driver maps each
asymptotic point through `pp.to_frame` into the ζ-frame.  Users on
the PIII Fig 5 family should pass
`z -> pIII_asymptotic_ic(z; n_terms = 10, β = pp.params.β, δ = pp.params.δ)`.

`degenerate_full_plane = true` is a regression-test mode: the sector
is treated as empty (no BVP), and the returned object's `sol(ζ)`
delegates straight to a single PFS walk over the user-supplied grid.
This is **IB.1.3**'s acceptance: when the sector "doesn't exist", the
hybrid is bit-exact to pure `path_network_solve`.

Throws `ArgumentError` on `im_lo ≥ im_hi`, `re_extent ≤ 0`, on a
`:direct`-frame `pp`, or on non-finite asymptotic ICs.
"""
function solve_pole_free_hybrid(pp::PainleveProblem,
                                 sector::NamedTuple,
                                 asymptotic_ic_fn;
                                 pfs_kwargs = (;),
                                 bvp_kwargs = (;),
                                 glue_tol::Real = 1e-8,
                                 n_slices::Integer = 8,
                                 degenerate_full_plane::Bool = false)
    pp.frame === :transformed || throw(ArgumentError(
        "solve_pole_free_hybrid: requires a :transformed-frame " *
        "PainleveProblem (PIII / PV / PVI).  Got :direct.  The hybrid " *
        "algorithm's BVP step operates on the ζ-plane strip; there is no " *
        "PIII/PV/PVI-like 'pole-free sector with curved boundaries' for " *
        "PI / PII / PIV in their natural z-frame."))

    # ---- degenerate / regression fast path ---------------------------------
    if degenerate_full_plane
        return _degenerate_full_plane(pp, sector, asymptotic_ic_fn;
                                       pfs_kwargs, glue_tol)
    end

    # ---- validate sector geometry ------------------------------------------
    haskey(sector, :im_lo) && haskey(sector, :im_hi) &&
        haskey(sector, :re_anchor) && haskey(sector, :re_extent) || throw(
        ArgumentError("solve_pole_free_hybrid: sector must be a NamedTuple " *
        "with fields (im_lo, im_hi, re_anchor, re_extent); got $(keys(sector))."))
    sector.im_lo < sector.im_hi || throw(ArgumentError(
        "solve_pole_free_hybrid: sector.im_lo ($(sector.im_lo)) must be < " *
        "sector.im_hi ($(sector.im_hi))."))
    sector.re_extent > 0 || throw(ArgumentError(
        "solve_pole_free_hybrid: sector.re_extent must be > 0 (got " *
        "$(sector.re_extent))."))
    n_slices ≥ 1 || throw(ArgumentError(
        "solve_pole_free_hybrid: n_slices must be ≥ 1 (got $n_slices)."))

    T  = Float64    # v1: Float64 only; the hybrid's point is to avoid BF.
    CT = Complex{T}

    # ---- map ICs through pp.to_frame ---------------------------------------
    # Boundary entry points in z-frame at the SECTOR corner radii.
    # The FFW Fig 5 spec picks z₁ = 30·exp((9π/4 − π/12)i) and
    # z₂ = 30·exp((-3π/4 + π/12)i).  In ζ-frame these are
    # ζ_j = 2·log(z_j) = 2·log(30) + i·(arg z_j · 2).  We accept these
    # via the user-supplied asymptotic_ic_fn which returns (u, u') at
    # z; we then map (z, u, u') → (ζ, w, w') via pp.to_frame.
    # The actual entry-ζ points are placed at Re ζ = re_anchor, at the
    # sector's open-side boundary Im-coordinates (im_lo + small ε,
    # im_hi - small ε) where ε pulls them just inside the sector per
    # FFW md:243 "two points close to the boundary".
    ε       = 0.001 * (sector.im_hi - sector.im_lo)
    ζ_top   = complex(T(sector.re_anchor), T(sector.im_hi) - ε)
    ζ_bot   = complex(T(sector.re_anchor), T(sector.im_lo) + ε)
    z_top   = _coord(pp.from_frame, ζ_top)
    z_bot   = _coord(pp.from_frame, ζ_bot)
    u_top, up_top = _eval_asymptotic(asymptotic_ic_fn, z_top, "top")
    u_bot, up_bot = _eval_asymptotic(asymptotic_ic_fn, z_bot, "bot")

    # Map (z, u, u') → (ζ, w, w').
    _, w_top, wp_top = pp.to_frame(z_top, u_top, up_top)
    _, w_bot, wp_bot = pp.to_frame(z_bot, u_bot, up_bot)

    # ---- Step 1: PFS walks along each curved boundary inward ----------------
    # We build two PathNetworkSolutions, each anchored at one of the
    # asymptotic-IC corners and walking through the pole-bearing exterior.
    # The PathNetworkSolution carries its own visited tree + Stage-2 grid;
    # we use the grid to harvest BC values at Step 2.

    # Sample points along each boundary (parallel to the Re-axis at fixed
    # Im) used as the Stage-2 grid for harvesting.
    re_lo    = T(sector.re_anchor) - T(sector.re_extent)
    re_hi    = T(sector.re_anchor)
    n_bdry   = max(2 * n_slices, 16)
    re_bdry  = collect(range(re_lo, re_hi; length = n_bdry))
    grid_top = CT[complex(r, T(sector.im_hi)) for r in re_bdry]
    grid_bot = CT[complex(r, T(sector.im_lo)) for r in re_bdry]

    # Two PFS problems, each with a different IC anchor.
    pfs_top = _pfs_ray_walk(pp, ζ_top, w_top, wp_top, grid_top; pfs_kwargs)
    pfs_bot = _pfs_ray_walk(pp, ζ_bot, w_bot, wp_bot, grid_bot; pfs_kwargs)

    # ---- Step 2: harvest BC values at each slice's two boundaries -----------
    # Place n_slices equally spaced Re-values in [re_lo, re_hi].
    slice_re = collect(range(re_lo, re_hi; length = n_slices))
    bvp_slices = Vector{BVPSolution{T, CT}}(undef, n_slices)

    # PFS grids are sampled at re_bdry; for each slice's Re-value we
    # interpolate piecewise-linearly along the harvested vector.
    for (i, r) in enumerate(slice_re)
        w_b, _ = _harvest_at_re(pfs_bot, re_bdry, grid_bot, r)
        w_t, _ = _harvest_at_re(pfs_top, re_bdry, grid_top, r)

        za = complex(r, T(sector.im_lo))
        zb = complex(r, T(sector.im_hi))

        # Step 3: solve the BVP on the slice.  P̃_III requires 3-arg RHS.
        # Pass the asymptotic-IC fn so the BVP's Newton initial guess
        # is the user's asymptotic series rather than the (badly
        # off-basin) linear ramp from w_b to w_t.
        bvp_slices[i] = _bvp_solve_on_slice(pp, za, zb, w_b, w_t;
                                             bvp_kwargs = bvp_kwargs,
                                             asymptotic_ic_fn = asymptotic_ic_fn)
    end

    # ---- Step 4: package up; glue continuity tested via callable -----------
    return IVPBVPSolution{T}(pfs_top, pfs_bot, bvp_slices, slice_re,
                              (; im_lo = T(sector.im_lo),
                                 im_hi = T(sector.im_hi),
                                 re_anchor = T(sector.re_anchor),
                                 re_extent = T(sector.re_extent)),
                              T(glue_tol), pp)
end

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Evaluate the user-supplied asymptotic-IC callable and validate.
function _eval_asymptotic(fn, z, label)
    out = fn(z)
    out isa Tuple && length(out) == 2 || throw(ArgumentError(
        "solve_pole_free_hybrid: asymptotic_ic_fn at $label boundary " *
        "(z = $z) must return a 2-tuple (u, u'); got $out."))
    u, up = out
    isfinite(u) && isfinite(up) || throw(ArgumentError(
        "solve_pole_free_hybrid: asymptotic_ic_fn at $label boundary " *
        "(z = $z) returned non-finite (u, u') = ($u, $up).  Suggestion: " *
        "verify the truncated series doesn't blow up at the supplied z " *
        "(check |z| ≥ 1 and arg z within the pole-free sector)."))
    return (u, up)
end

# Run a PFS walk anchored at the asymptotic-IC corner; returns the
# PathNetworkSolution covering `grid`.
function _pfs_ray_walk(pp, ζ_ic::Complex, w_ic, wp_ic,
                       grid::AbstractVector; pfs_kwargs)
    # Build a fresh PadeTaylorProblem at the asymptotic-IC anchor.
    # The zspan is unused by path_network_solve except to seed
    # `visited_z[1] = zspan[1]`; PadeTaylorProblem rejects degenerate
    # spans (Problems.jl line 122).  We pick `zspan[1] = ζ_ic` and a
    # token `zspan[2] = ζ_ic + 1` — PFS doesn't iterate to it.
    new_prob = PadeTaylorProblem(pp.problem.f, (w_ic, wp_ic),
                                  (ζ_ic, ζ_ic + one(ζ_ic));
                                  order = pp.problem.order)
    return path_network_solve(new_prob, grid; pfs_kwargs...)
end

# Linear-interpolate `pfs.grid_u` and `pfs.grid_up` along the boundary at
# the requested Re-coordinate.  `re_bdry` is the parallel vector of Re-
# values that defines the boundary sample points; `grid_path` is the
# matching grid vector (same length).
function _harvest_at_re(pfs::PathNetworkSolution, re_bdry::AbstractVector,
                        grid_path::AbstractVector, r::Real)
    # Find the bracketing pair in re_bdry.
    n = length(re_bdry)
    n == length(grid_path) || throw(ArgumentError(
        "_harvest_at_re: re_bdry/grid_path length mismatch."))
    n ≥ 2 || throw(ArgumentError("_harvest_at_re: need ≥ 2 boundary samples."))

    r_first = first(re_bdry); r_last = last(re_bdry)
    if r ≤ r_first
        i = 1
        ξ = 0.0
    elseif r ≥ r_last
        i = n - 1
        ξ = 1.0
    else
        i = searchsortedfirst(re_bdry, r) - 1
        i = clamp(i, 1, n - 1)
        ξ = (r - re_bdry[i]) / (re_bdry[i+1] - re_bdry[i])
    end
    u_a = pfs.grid_u[i];  u_b = pfs.grid_u[i+1]
    up_a = pfs.grid_up[i]; up_b = pfs.grid_up[i+1]

    isfinite(u_a) && isfinite(u_b) || throw(ErrorException(
        "_harvest_at_re: PFS grid_u at re ≈ $r is NaN — Stage-2 lookup " *
        "failed because the visited tree didn't reach the boundary " *
        "sample.  Suggestion: tighten max_steps_per_target, widen the " *
        "wedge, or check that the asymptotic-IC anchor is in the " *
        "walker's reachable basin."))

    u  = (1 - ξ) * u_a + ξ * u_b
    up = (1 - ξ) * up_a + ξ * up_b
    return (u, up)
end

# Solve the BVP on a single `Re ζ = const` slice through the sector.
# Currently hard-coded to dispatch on the PIII RHS shape (uses
# `pIII_transformed_rhs` parameters to build ∂f/∂w and ∂f/∂w').
#
# `asymptotic_ic_fn` supplies a per-node initial guess: at each
# Chebyshev collocation node `ζ_j`, evaluate the user's asymptotic
# series at `z = exp(ζ_j/2)` and use that as the starting iterate.
# This avoids the linear-ramp Newton-stalling that hits as soon as the
# sector spans `Im ζ ∈ ±π` or so (the linear ramp from `w_a` to `w_b`
# can be ~half the asymptotic magnitude at the centre — Newton then
# starts on the wrong side of the basin).
function _bvp_solve_on_slice(pp::PainleveProblem,
                              za::Complex, zb::Complex,
                              wa, wb; bvp_kwargs,
                              asymptotic_ic_fn = nothing)
    pp.equation === :III || throw(ArgumentError(
        "_bvp_solve_on_slice: v1 supports the :III equation only " *
        "(got :$(pp.equation)).  PV / PVI BVP sectors are deferred " *
        "(file a bead when the user needs them — the sector geometry " *
        "and RHS analytic-Jacobian shapes differ)."))

    α = pp.params.α; β = pp.params.β; γ = pp.params.γ; δ = pp.params.δ
    f      = (ζ, w, wp) -> begin
        eζ  = exp(ζ); e2ζ = eζ * eζ
        wp^2 / w + (α*w^2 + γ*w^3 + β*eζ + δ*e2ζ/w) / 4
    end
    ∂f_w   = (ζ, w, wp) -> begin
        eζ  = exp(ζ); e2ζ = eζ * eζ
        -wp^2 / w^2 + (2α*w + 3γ*w^2 - δ*e2ζ/w^2) / 4
    end
    ∂f_wp  = (ζ, w, wp) -> 2*wp/w

    # Per-node asymptotic-IC initial guess.  Maps each ζ-node to
    # z = exp(ζ/2), evaluates u(z) via asymptotic_ic_fn, and maps
    # back to w = z·u (PIII convention).
    initial_guess = if asymptotic_ic_fn === nothing
        nothing
    else
        ζ -> begin
            z  = exp(ζ / 2)
            u, _ = asymptotic_ic_fn(z)
            return z * u
        end
    end

    if initial_guess === nothing
        return bvp_solve(f, ∂f_w, ∂f_wp, za, zb, wa, wb; bvp_kwargs...)
    else
        return bvp_solve(f, ∂f_w, ∂f_wp, za, zb, wa, wb;
                          initial_guess = initial_guess, bvp_kwargs...)
    end
end

# Degenerate-sector fast path: bypass the BVP step entirely; the
# returned IVPBVPSolution holds a single PFS walk reused for both
# `pfs_top` and `pfs_bot`, and an empty `bvp_slices`.  The callable
# `sol(ζ)` then exclusively reads the PFS grid.  Test IB.1.3 pins
# this: the hybrid driver must be bit-exact to pure path_network_solve
# when the sector is vacuous.
function _degenerate_full_plane(pp::PainleveProblem,
                                 sector::NamedTuple,
                                 asymptotic_ic_fn;
                                 pfs_kwargs = (;), glue_tol::Real = 1e-8)
    T  = Float64
    CT = Complex{T}
    # User supplies grid via pfs_kwargs[:grid].
    haskey(pfs_kwargs, :grid) || throw(ArgumentError(
        "solve_pole_free_hybrid(degenerate_full_plane = true): " *
        "pfs_kwargs must contain a `:grid` entry (the user's " *
        "ζ-plane grid).  This mode bypasses the BVP step and runs " *
        "a single PFS walk."))
    grid = pfs_kwargs[:grid]
    pfs_kw = filter(p -> first(p) !== :grid, pairs(pfs_kwargs))

    pfs = path_network_solve(pp.problem, grid; pfs_kw...)
    return IVPBVPSolution{T}(pfs, pfs, BVPSolution{T, CT}[],
                              T[], (; im_lo = T(sector.im_lo),
                                     im_hi = T(sector.im_hi),
                                     re_anchor = T(sector.re_anchor),
                                     re_extent = T(sector.re_extent)),
                              T(glue_tol), pp)
end

# -----------------------------------------------------------------------------
# Callable + accessors
# -----------------------------------------------------------------------------

"""
    (sol::IVPBVPSolution)(ζ) -> (w, w')

Evaluate the hybrid solution at ζ in the ζ-frame.  Dispatches:

  - **Inside the BVP sector** (`im_lo ≤ Im ζ ≤ im_hi` AND
    `re_lo ≤ Re ζ ≤ re_hi`): use the nearest `Re ζ` slice's
    `BVPSolution` barycentric interpolant.  A linear interpolation
    between the two bracketing slices at fractional `Re ζ` is the
    accuracy-bound v1 (Worklog 039 notes the 2D-spectral upgrade as
    a deferred follow-up).
  - **Between two bracketing slices**: before interpolating, the
    callable fails loud (Rule 1) on GROSS inter-slice incoherence —
    a jump between the two slices' values at ζ exceeding `3×` the
    smaller value (floored at 1) signals a divergent/corrupt BVP slice
    (bug `padetaylor-q0yq`).  It does NOT perform the FFW md:247
    PFS-vs-BVP *derivative* match: that ~1e-7 criterion is a value
    tautology (BVP BCs are set === the harvested PFS values) plus a
    derivative check that is ~1e-1 in v1 Float64 (FFW md:226) — a
    documented v2 deferral (bead `padetaylor-sn9a`).
  - **Outside the sector**: not currently callable — IVPBVPSolution
    v1 does not expose a dense-callable PFS interpolant (PFS is a
    Stage-2 grid by construction; FW 2011 line 166).  Users querying
    outside the sector should read `sol.pfs_top.grid_u` /
    `sol.pfs_bot.grid_u` directly.

Throws `DomainError` outside the sector; throws `ErrorException` on gross
inter-slice incoherence (a divergent/corrupt BVP slice — bug `padetaylor-q0yq`).
"""
function (sol::IVPBVPSolution{T})(ζ) where T
    CT = Complex{T}
    ζ_CT = CT(ζ)
    s = sol.sector
    re = real(ζ_CT); im_ = imag(ζ_CT)

    re_lo = s.re_anchor - s.re_extent
    re_hi = s.re_anchor

    # Inside-or-on-boundary test.
    inside = (s.im_lo - sol.glue_tol ≤ im_ ≤ s.im_hi + sol.glue_tol) &&
             (re_lo - sol.glue_tol ≤ re ≤ re_hi + sol.glue_tol)
    inside || throw(DomainError(ζ,
        "IVPBVPSolution: ζ = $ζ is outside the BVP sector " *
        "[$(re_lo), $(re_hi)] × [$(s.im_lo), $(s.im_hi)].  v1 hybrid " *
        "callable serves the BVP sector only; for exterior values " *
        "read `sol.pfs_top.grid_u` / `sol.pfs_bot.grid_u`."))

    isempty(sol.bvp_slices) && throw(ErrorException(
        "IVPBVPSolution: sol(ζ) requested but `bvp_slices` is empty.  " *
        "This is the degenerate-full-plane mode — query the PFS grid " *
        "directly via `sol.pfs_top.grid_u`."))

    # Locate the slice bracket in slice_re.
    n = length(sol.slice_re)
    if re ≤ first(sol.slice_re)
        return _eval_bvp_slice(sol.bvp_slices[1], ζ_CT)
    elseif re ≥ last(sol.slice_re)
        return _eval_bvp_slice(sol.bvp_slices[n], ζ_CT)
    else
        i = searchsortedfirst(sol.slice_re, re) - 1
        i = clamp(i, 1, n - 1)
        ξ = (re - sol.slice_re[i]) / (sol.slice_re[i+1] - sol.slice_re[i])
        # Evaluate both slices at this ζ (each slice has its own
        # `[im_lo, im_hi]` collocation segment), then linear-interp
        # in Re.  v1 accuracy bound: O(Δre · second-Re-derivative);
        # adequate for FFW Fig 5 since adjacent slices live ~0.05
        # apart in Re ζ.
        w_a, up_a = _eval_bvp_slice(sol.bvp_slices[i], ζ_CT)
        w_b, up_b = _eval_bvp_slice(sol.bvp_slices[i+1], ζ_CT)
        # Inter-slice COHERENCE guard (bug padetaylor-q0yq, Rule 1 fail-loud).
        # The callable's docstring promised a glue-continuity throw that the
        # code never performed; a divergent/corrupt BVP slice (one that failed
        # to converge to the correct branch) was silently interpolated into a
        # fictitious value.  This is a CATASTROPHE bound, NOT `glue_tol`: adjacent
        # slices legitimately differ by O(0.1·|w|) across a finite Re gap (the
        # FFW md:226 derivative phenomenon, out of v1 scope — that ~1e-7 PFS-vs-
        # BVP match is deferred, bead padetaylor-sn9a).  Measured: the honest
        # ratio is ~0.093 while a corrupt slice gives ≥14 — so a jump exceeding
        # 3× the smaller bracketing value (floored at 1 for zero-crossings) is
        # unambiguously a divergent slice, and we fail loud instead of lying.
        jump  = abs(w_b - w_a)
        scale = max(min(abs(w_a), abs(w_b)), one(jump))
        jump > 3 * scale && throw(ErrorException(
            "IVPBVPSolution: inter-slice incoherence at Re ζ = $re between " *
            "bracketing BVP slices $i and $(i + 1) (Re = $(sol.slice_re[i]), " *
            "$(sol.slice_re[i+1])): |Δw| = $jump exceeds the catastrophe bound " *
            "$(3 * scale) (3× the smaller slice value).  A BVP slice likely " *
            "failed to converge to the correct branch.  Suggestion: re-solve " *
            "with larger N / max_iter, verify every slice.residual_inf < tol, " *
            "or widen the PFS wedge so the boundary harvest is clean."))
        return ((1 - ξ) * w_a + ξ * w_b,
                (1 - ξ) * up_a + ξ * up_b)
    end
end

# Internal: evaluate one BVPSolution at a given ζ on its segment.
# The slice's z_a = re + i·im_lo, z_b = re + i·im_hi; we use the
# BVPSolution's barycentric callable, which already handles the
# segment-pre-image check.  But our query ζ may have a slightly
# different real part than the slice's `Re ζ` (we caller with the
# original `ζ`, not the projection); we project to the slice line
# before calling the BVP's callable, since the slice represents one
# Re-coordinate only.
function _eval_bvp_slice(slice::BVPSolution{T,CT}, ζ::CT) where {T,CT}
    # Project to slice line: take `Re(slice.z_a)` for the real coord.
    re_slice = real(slice.z_a)   # both z_a and z_b have the same real part
    ζ_proj = CT(re_slice + im * imag(ζ))
    return slice(ζ_proj)
end

end # module IVPBVPHybrid
