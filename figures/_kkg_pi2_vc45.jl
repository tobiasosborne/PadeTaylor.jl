# figures/_kkg_pi2_vc45.jl
#
# VC-4 and VC-5 — the per-pole validation criteria of the KKG 2015
# P_I^(2) tritronquée headline figure (ADR-0025 Phase D, beads
# `padetaylor-0ln.37.12` / `0ln.37.13`).  This helper is `include`d by
# the figure kernel `_kkg_pi2_surface_helpers.jl`; it certifies the
# genuine poles of the wedge pole field and prunes the spurious ones.
#
# ## Why a separate validation pass is needed
#
# The B3-rewired wedge walk threads a node filament through the *whole*
# pole-rich wedge out to `|x| = 20` and `extract_poles_shared_q`
# (cross-node `min_support ≥ 2`, VC-6) clusters the per-node shared-`Q`
# roots into a candidate field of ~380 poles.  VC-6 is a *cross-node*
# spurious-pole filter — it rejects single-node linear-system artefacts
# — but it cannot reject a *Froissart doublet* (a near-cancelling
# numerator/denominator pole pair that recurs across nodes) nor a
# mis-located cluster.  The A4 baseline probe (ADR-0025 Amendment 3)
# already found ≥1 spurious pole among the old 21-pole V8b field: the
# candidate at `+5.44+2.40im` fits leading coefficient `A ≈ 0`, not the
# double pole the P_I^(2) Painlevé structure demands.
#
# VC-4 and VC-5 are the per-pole certificates that close that gap.
# They are *independent of the path-network extraction algorithm under
# test* — VC-4 tests the local ODE structure, VC-5 tests a symmetry of
# the tritronquée itself — so a bug in the walk cannot mask a bad pole.
#
# ## VC-4 — dominant-balance per-pole structural certificate
#
# ADR-0025 Amendment 1 §A3 (`external/probes/pi2-pole-structure/
# REPORT.md`): every genuine movable singularity of P_I^(2) at `t = 0`
# is a **double pole** with a Laurent expansion `u = Σ aₖ ξ^{k-2}`,
# `ξ = x - x₀`, whose leading coefficient is fixed by the dominant
# balance `A(A+1)(A+3) = 0`:
#
#   * `A = -1` — the generic 4-parameter family (V₀'s wedge poles);
#   * `A = -3` — a codimension-1 special family;
#   * and in **both** families the residue vanishes exactly,
#     `a₁ = 0` (no `ξ^{-1}` term, Painlevé property, no logarithms).
#
# So for each candidate pole `x₀` we fit, by complex least squares on a
# ring of 32 points,
#
#       u(x) ≈ A·(x - x₀)^{-2} + B·(x - x₀)^{-1} + C
#
# and apply the two A3 acceptance tests:
#
#   * **VC-4a** `min(|A+1|, |A+3|) < 0.10` — the leading coefficient
#     lies in the dominant-balance family `{-1, -3}` (the ODE-structure
#     check; `0.10` is 5 % of the gap `2` between the two roots).
#   * **VC-4b** `|B| < 0.10·|A|` — the residue vanishes (the
#     zero-residue / Painlevé-property check).
#
# A candidate that fails **either** is spurious — a Froissart doublet
# (`A ≈ 0`, no double-pole structure) or an out-of-family artefact —
# and is **pruned** from the figure's pole field.  The ring radius is
# `r = min(0.05, 0.1·d_nearest)` with a floor `0.01` (`d_nearest` is
# the distance to the nearest *other* candidate): close enough that the
# `O(ξ²)` model-truncation error of the 3-term fit is ≲ 0.25 %, far
# enough that the ring stays clear of the neighbouring pole.
#
# `u` on the ring is read from the path-network solution exactly the
# way `extract_poles_shared_q` reads it: the shared-`Q` Padé of the
# visited node nearest the ring point — `t = (x - z_node)/h_node`,
# `u = P₁(t)/Q(t)`.  This is the V8b `eval_u_at` machinery
# (`external/probes/v8b-baseline/probe.jl:360`), reused verbatim so the
# fit is exercised against the same approximant the figure renders.
#
# ## VC-5 — conjugate-symmetry pole pairing
#
# The P_I^(2) ODE has real coefficients and `V₀` is real on the
# negative real axis, so `V₀(x̄) = conj(V₀(x))` everywhere: the wedge
# pole field must be **symmetric under conjugation**.  After VC-4
# pruning we match the surviving poles into conjugate pairs by greedy
# nearest-neighbour under `x ↦ x̄` and report the per-pair residual
# `|p_upper - conj(p_lower)|`.  Its median / max is itself an FW-style
# accuracy estimate of the pole field (FW 2011 `:303-311` — the
# conjugate-symmetry error estimate, manufactured because the
# tritronquée pole field has no external oracle).
#
# A pole left unpaired must sit *on* the real axis (`|Im x₀|` small) —
# its own conjugate.  A pole with a large `|Im|` and no conjugate
# partner is **flagged** as suspect: a genuine off-axis pole always has
# a mirror partner.
#
# References:
#   * `docs/adr/0025-headline-figure-re-resolution.md` — Amendment 1
#     §A3 (the VC-4 exact form), Validation Criteria Menu (VC-4, VC-5).
#   * `external/probes/pi2-pole-structure/REPORT.md` — the A3
#     dominant-balance / zero-residue derivation.
#   * `external/probes/v8b-baseline/REPORT.md` — the A4 baseline: the
#     ring-fit spot check, the one spurious `A ≈ 0` pole.
#   * `external/probes/v8b-baseline/probe.jl` — `eval_u_at` / `fit_ABC`,
#     the working reference fit this helper's `vc4_fit_pole` mirrors.

using PadeTaylor.VectorPathNetwork: VectorPathNetworkSolution

"""
    _vc_eval_poly(c, t) -> Complex

Horner evaluation of `Σ c[k+1]·tᵏ` (coefficients low-to-high) — the
convention the shared-`Q` numerators and denominator are stored in.  A
local copy of `VectorPathNetworkStage2._eval_poly` (private to that
module; copied, not imported, to respect the module boundary — exactly
as `VectorPathNetworkStage2` itself keeps a private copy of
`VectorProblems._eval_poly`, and as the V8b probe's `_eval_poly` does).
"""
function _vc_eval_poly(c::AbstractVector, t)
    s = zero(promote_type(eltype(c), typeof(t)))
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

# ======================================================================
# VC-4 / VC-5 tuning constants (ADR-0025 Amendment 1 §A3).
# ======================================================================

# VC-4 ring: 32 sample points; radius r = min(R_MAX, FRAC·d_nearest)
# floored at R_FLOOR (§A3 step 1).
const VC4_RING_NPTS  = 32
const VC4_RING_R_MAX = 0.05
const VC4_RING_FRAC  = 0.10
const VC4_RING_FLOOR = 0.01

# VC-4a / VC-4b acceptance tolerances (§A3 §7.3).
const VC4_TOL_A = 0.10      # min(|A+1|,|A+3|) < TOL_A  — the family check
const VC4_TOL_B = 0.10      # |B| < TOL_B·|A|           — the residue check

# VC-5: a pole is "on the real axis" (allowed to be unpaired) when its
# |Im| is below this.  A larger-|Im| unpaired pole is flagged suspect.
const VC5_REAL_AXIS_TOL = 0.15

# VC-5 conjugate-match acceptance: an upper / lower pole pair is a
# genuine conjugate pair only when `|p_upper - conj(p_lower)|` is below
# this.  A pole whose nearest conjugate-mirror candidate is farther
# than this has no partner — it is left unpaired (and, off-axis,
# flagged suspect) rather than force-matched into a meaningless pair.
# `0.5` is below the median pole nearest-neighbour spacing (~0.69 in
# the B3 wedge field) — far enough that a genuine mirror pair is kept,
# tight enough that two unrelated poles are not paired by accident.
const VC5_MATCH_TOL = 0.5

# ======================================================================
# VC-4 — the ring fit
# ======================================================================

"""
    vc4_eval_u(walk, z) -> ComplexF64

The first companion component `u = V₀(z)` recovered from the
path-network solution `walk` at the complex point `z`: locate the
visited node nearest `z`, rescale `t = (z - z_node)/h_node`, and
evaluate the node's shared-`Q` Padé `u = P₁(t)/Q(t)` (Horner via
`_eval_poly`).  This is exactly the approximant `extract_poles_shared_q`
roots and the figure renders — the V8b `eval_u_at` recipe
(`external/probes/v8b-baseline/probe.jl:360`).

Throws `ArgumentError` (Rule 1, fail-loud) on an empty walk — a pole
field cannot be validated against a solution with no nodes.
"""
function vc4_eval_u(walk::VectorPathNetworkSolution, z::Complex)
    N = length(walk.visited_z)
    N ≥ 1 || throw(ArgumentError(
        "vc4_eval_u: the path-network solution has no visited nodes; " *
        "cannot evaluate u for VC-4. suggestion: check the wedge walk " *
        "completed before validating its pole field."))
    idx, best = 1, abs(walk.visited_z[1] - z)
    @inbounds for i in 2:N
        d = abs(walk.visited_z[i] - z)
        d < best && ((idx, best) = (i, d))
    end
    h    = real(walk.visited_h[idx])
    t    = (z - walk.visited_z[idx]) / h
    den  = walk.visited_denominator[idx]
    q_t  = _vc_eval_poly(den, t)
    return _vc_eval_poly(walk.visited_numerators[idx][1], t) / q_t
end

"""
    vc4_ring_radius(x0, poles) -> Float64

The VC-4 ring radius for candidate `x0`: `min(0.05, 0.1·d_nearest)`
floored at `0.01`, where `d_nearest` is the distance from `x0` to the
nearest *other* pole in `poles` (ADR-0025 Amendment 1 §A3 step 1).  The
floor keeps the ring numerically resolvable when the field is locally
dense; the `0.1·d_nearest` cap keeps the ring clear of a neighbour.
"""
function vc4_ring_radius(x0::Complex, poles::AbstractVector{<:Complex})
    d_nearest = Inf
    for p in poles
        p === x0 && continue
        d = abs(p - x0)
        (d > 0 && d < d_nearest) && (d_nearest = d)
    end
    r = isfinite(d_nearest) ? min(VC4_RING_R_MAX, VC4_RING_FRAC * d_nearest) :
                              VC4_RING_R_MAX
    return max(VC4_RING_FLOOR, r)
end

"""
    vc4_fit_pole(walk, x0, poles) -> (A, B, C, r)

Fit the local Laurent model `u ≈ A·ξ⁻² + B·ξ⁻¹ + C`, `ξ = x - x0`, by
complex least squares on a ring of `VC4_RING_NPTS` points at radius
`vc4_ring_radius(x0, poles)` about `x0` (ADR-0025 Amendment 1 §A3).
`u` on the ring is read by `vc4_eval_u` (the shared-`Q` Padé of the
nearest visited node).  Returns the fitted `(A, B, C)` and the ring
radius `r`.  The `npts × 3` design matrix has rows `[ξ⁻², ξ⁻¹, 1]`;
`M \\ b` is the least-squares solve (32 equations, 3 complex unknowns).
"""
function vc4_fit_pole(walk::VectorPathNetworkSolution, x0::Complex,
                      poles::AbstractVector{<:Complex})
    r  = vc4_ring_radius(x0, poles)
    np = VC4_RING_NPTS
    M  = Matrix{ComplexF64}(undef, np, 3)
    b  = Vector{ComplexF64}(undef, np)
    for i in 1:np
        θ = 2π * (i - 1) / np
        ξ = ComplexF64(r * cis(θ))
        u = vc4_eval_u(walk, x0 + ξ)
        M[i, 1] = 1 / ξ^2
        M[i, 2] = 1 / ξ
        M[i, 3] = 1
        b[i]    = u
    end
    coef = M \ b
    return coef[1], coef[2], coef[3], r
end

"""
    vc4_validate(walk, poles) -> NamedTuple

Run VC-4 over every candidate pole in `poles`: fit `(A, B, C)` on a
ring (`vc4_fit_pole`), then apply

  * VC-4a `min(|A+1|, |A+3|) < VC4_TOL_A`  — the dominant-balance
    family check (`A ∈ {-1,-3}`, ADR-0025 A3);
  * VC-4b `|B| < VC4_TOL_B·|A|`            — the zero-residue check.

A candidate passing **both** is genuine and kept; one failing either is
spurious and pruned.  Returns
`(kept, A, B, dA, pass, family, n_pruned, prune_reason)`:

  - `kept`         : the VC-4-surviving pole field (the figure's field);
  - `A`, `B`       : the per-*candidate* fitted coefficients (diagnostic,
                     in `poles` order);
  - `dA`           : per-candidate `min(|A+1|, |A+3|)`;
  - `pass`         : per-candidate `Bool` — survived VC-4;
  - `family`       : per-candidate `:m1` (`A≈-1`) / `:m3` (`A≈-3`) /
                     `:none` (failed VC-4a) — the A-family breakdown;
  - `n_pruned`     : how many candidates were pruned;
  - `prune_reason` : per-*pruned*-candidate `(pole, :froissart | :out_of_family
                     | :nonzero_residue)` — `:froissart` is the
                     `|A| < 0.1` near-zero-leading-coefficient signature.
"""
function vc4_validate(walk::VectorPathNetworkSolution,
                      poles::AbstractVector{<:Complex})
    n = length(poles)
    A   = Vector{ComplexF64}(undef, n)
    B   = Vector{ComplexF64}(undef, n)
    dA  = Vector{Float64}(undef, n)
    pass   = Vector{Bool}(undef, n)
    family = Vector{Symbol}(undef, n)
    prune_reason = Tuple{ComplexF64, Symbol}[]
    kept = ComplexF64[]
    for k in 1:n
        x0 = ComplexF64(poles[k])
        a, bb, _, _ = vc4_fit_pole(walk, x0, poles)
        A[k], B[k]  = a, bb
        d4a   = min(abs(a + 1), abs(a + 3))
        dA[k] = d4a
        ok_a  = d4a < VC4_TOL_A
        ok_b  = abs(bb) < VC4_TOL_B * abs(a)
        family[k] = ok_a ? (abs(a + 1) ≤ abs(a + 3) ? :m1 : :m3) : :none
        pass[k]   = ok_a && ok_b
        if pass[k]
            push!(kept, x0)
        else
            # Classify the failure mode for the diagnostic report.
            reason = if abs(a) < VC4_TOL_A
                :froissart           # A ≈ 0 — no double-pole structure
            elseif !ok_a
                :out_of_family       # a double pole, but A ∉ {-1,-3}
            else
                :nonzero_residue     # right A, but |B| too large
            end
            push!(prune_reason, (x0, reason))
        end
    end
    return (kept = kept, A = A, B = B, dA = dA, pass = pass,
            family = family, n_pruned = length(prune_reason),
            prune_reason = prune_reason)
end

# ======================================================================
# VC-5 — conjugate-symmetry pole pairing
# ======================================================================

"""
    vc5_pair(poles) -> NamedTuple

Match the VC-4-surviving poles into conjugate pairs.  `V₀(x̄) =
conj(V₀(x))`, so a genuine off-axis pole `p` has a mirror partner
`conj(p)`; the pole field must be symmetric under conjugation.

**Globally-greedy assignment under `x ↦ x̄`.**  Enumerate every
`(p_upper, p_lower)` candidate pair with residual `|p_upper -
conj(p_lower)| ≤ VC5_MATCH_TOL`, sort by residual, and consume them in
ascending order — the closest mirror pair is committed first, each
pole used at most once.  This is near-optimal and strictly better than
the `|Im|`-ordered scan, which mis-assigns partners in a dense field
(it would force the *first* upper pole onto a lower pole that is a
better match for a *later* upper pole).  A candidate pair farther than
`VC5_MATCH_TOL` is **not** a conjugate pair — both poles are left
unpaired rather than force-matched, so the residual distribution stays
a meaningful accuracy estimate (FW 2011 `:303-311`).

A pole left unpaired must lie *on* the real axis (`|Im| <
VC5_REAL_AXIS_TOL` — it is its own conjugate).  An unpaired pole with a
larger `|Im|` is **flagged** suspect — a genuine off-axis tritronquée
pole always has a mirror partner.

Returns `(pairs, residuals, median_resid, max_resid, unpaired,
flagged)`:
  - `pairs`        : `Vector{Tuple{ComplexF64,ComplexF64}}` of
                     `(p_upper, p_lower)`;
  - `residuals`    : the per-pair `|p_upper - conj(p_lower)|`;
  - `median_resid`, `max_resid` : the FW-style accuracy estimate
                     (`NaN` when there are no pairs);
  - `unpaired`     : the poles with no conjugate partner;
  - `flagged`      : the subset of `unpaired` with `|Im| ≥
                     VC5_REAL_AXIS_TOL` — suspect off-axis poles.
"""
function vc5_pair(poles::AbstractVector{<:Complex})
    P     = ComplexF64.(poles)
    upper = filter(p -> imag(p) >  VC5_REAL_AXIS_TOL, P)
    lower = filter(p -> imag(p) < -VC5_REAL_AXIS_TOL, P)
    # Enumerate every admissible mirror-pair candidate (residual within
    # VC5_MATCH_TOL), then commit them globally-greedily in ascending
    # residual — each pole consumed at most once.
    cand = Tuple{Float64,Int,Int}[]
    for (iu, pu) in enumerate(upper), (il, pl) in enumerate(lower)
        d = abs(pu - conj(pl))
        d ≤ VC5_MATCH_TOL && push!(cand, (d, iu, il))
    end
    sort!(cand; by = c -> c[1])
    u_used = falses(length(upper))
    l_used = falses(length(lower))
    pairs     = Tuple{ComplexF64,ComplexF64}[]
    residuals = Float64[]
    for (d, iu, il) in cand
        (u_used[iu] || l_used[il]) && continue
        u_used[iu] = true; l_used[il] = true
        push!(pairs, (upper[iu], lower[il]))
        push!(residuals, d)
    end
    paired = Set{ComplexF64}()
    for (pu, pl) in pairs
        push!(paired, pu); push!(paired, pl)
    end
    unpaired = [p for p in P if !(p in paired)]
    flagged  = [p for p in unpaired if abs(imag(p)) ≥ VC5_REAL_AXIS_TOL]
    med = isempty(residuals) ? NaN :
          sort(residuals)[cld(length(residuals), 2)]
    mx  = isempty(residuals) ? NaN : maximum(residuals)
    return (pairs = pairs, residuals = residuals,
            median_resid = med, max_resid = mx,
            unpaired = unpaired, flagged = flagged)
end
