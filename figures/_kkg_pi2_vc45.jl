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
# pruning we match the surviving poles into conjugate pairs under
# `x ↦ x̄` and report the per-pair residual `|p_upper - conj(p_lower)|`.
# Its median / max is itself an FW-style accuracy estimate of the pole
# field (FW 2011 `:303-311` — the conjugate-symmetry error estimate,
# manufactured because the tritronquée pole field has no external
# oracle).  Crucially the matching only *pairs* poles the walk
# independently extracted; it never *constructs* a missing pole from
# its mirror — so the residual stays a genuine accuracy cross-check
# (FFW 2017 `:120-124`: FFW deliberately do not use the conjugate
# symmetry in their numerics for exactly this reason).
#
# A pole left unpaired must sit *on* the real axis (`|Im x₀|` small) —
# its own conjugate.  A pole with a large `|Im|` and no conjugate
# partner is **flagged** as suspect: a genuine off-axis pole always has
# a mirror partner.
#
# ## VC-5b — the matching is globally optimal (ADR-0025 Amendment 6)
#
# The v1 `vc5_pair` matched by a *globally-greedy* commit (tightest
# admissible mirror pair first, each pole consumed once).  Amendment 5
# reported the resulting defect: of 266 VC-4 survivors, only 93 paired
# and **72 off-axis poles were flagged unpaired**.  Bead
# `padetaylor-0ln.37.20` (VC-5b) investigated the root cause:
#
#   * the wedge *node* coverage is conjugate-symmetric (the upper- and
#     lower-half walk filaments thread near-mirror loci, median gap
#     ~0.10 even in the far wedge) — the walk does **not** miss
#     far-wedge regions;
#   * but the far wedge is at the A2 tractability ceiling: the two
#     half-trees develop different adaptive-`h` histories and resolve
#     partially-disjoint subsets of the *dense* far-wedge pole lattice,
#     so the pole field carries an intrinsic conjugate-residual
#     accuracy of **median ~0.3, p75 ~0.5** (a per-pole VC-4-residual
#     polish moves each pole only ~0.005 — every pole already sits at
#     its dominant-balance-optimal location, so this is genuine field
#     accuracy, not a correctable per-pole error);
#   * two defects in VC-5 *itself* then inflate the flag count: the
#     greedy matcher is not maximum-cardinality (it blocks pairs a
#     proper matching finds — 93 greedy vs 103 optimal on the B3
#     field), and the `VC5_MATCH_TOL = 0.5` cutoff was a guess set
#     *tighter than the measured field accuracy*.
#
# Amendment 6 fixes both: `vc5_pair` now computes a **maximum-
# cardinality bipartite matching** of the conjugate-admissibility graph
# (Kuhn augmenting paths), and `VC5_MATCH_TOL` is re-derived from the
# measured 0.69 pole spacing / ~0.3-0.5 field accuracy to **0.6**.
# Together: 72 → 52 flagged, 93 → 103 conjugate pairs.  No conjugate
# symmetry is imposed by construction — the two half-fields are still
# independently extracted and the residual stays a real accuracy
# estimate.  The residual *distribution itself widens slightly* (median
# 0.24 → 0.29) precisely because the optimal matcher no longer discards
# the harder-to-pair tail — an honest diagnostic, not a cosmetic one.
#
# References:
#   * `docs/adr/0025-headline-figure-re-resolution.md` — Amendment 1
#     §A3 (the VC-4 exact form), Amendment 5 (the 72-unpaired defect),
#     Amendment 6 (the VC-5b root cause + this optimal-matching fix),
#     Validation Criteria Menu (VC-4, VC-5).
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
#
# `VC5_MATCH_TOL` must satisfy two bounds: it must be **below** the pole
# nearest-neighbour spacing (so two *distinct* lattice poles can never
# be matched as a false conjugate pair), and **above** the pole field's
# own conjugate-residual accuracy (so a genuine mirror pair is not
# rejected for the field's intrinsic placement error).
#
# The VC-5b investigation (bead `padetaylor-0ln.37.20`, ADR-0025
# Amendment 6) *measured* both bounds on the B3 wedge field:
#
#   * pole nearest-neighbour spacing — median **0.69**;
#   * conjugate-residual accuracy — the offset between an accurately
#     extracted pole and its conjugate's extracted partner is **median
#     ~0.3, p75 ~0.5**; the far wedge is at the A2 tractability ceiling,
#     so this is the field's genuine intrinsic accuracy, not a
#     correctable extraction error (the VC-5b probe proved a per-pole
#     VC-4-residual polish moves each pole by only ~0.005 — every pole
#     already sits at its dominant-balance-optimal location; the ~0.3
#     residual is two half-trees resolving partially-disjoint subsets
#     of the dense far-wedge pole lattice).
#
# The v1 value `0.5` was a *conservative guess* — "below 0.69" — set
# before the field accuracy was measured.  It is **tighter than the
# field's own accuracy**: it bisects the conjugate-residual distribution
# and rejects ~27 % of genuine mirror pairs (the 72-unpaired VC-5
# defect of Amendment 5).  Amendment 6 re-derives it from the measured
# numbers to **0.6**: still firmly below the 0.69 spacing (87 % of it),
# and — paired with the globally-optimal matching of `vc5_pair` below
# (a pole pairs only with its *globally* best mirror) — it cannot
# false-match two distinct lattice poles.  This is a Law-1 ground-truth
# correction of a guessed parameter, not a tolerance relaxed to pass a
# test.
const VC5_MATCH_TOL = 0.6

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

# --- the optimal-matching augmenting-path core ---------------------------

"""
    _vc5_augment!(u, adj, match_l, seen) -> Bool

One Kuhn augmenting-path step of the Hopcroft–Karp / Hungarian-style
maximum-cardinality bipartite matching used by `vc5_pair`.  `adj[u]` is
the admissible lower-pole index list of upper pole `u`, *pre-sorted by
ascending conjugate residual* so the augmenting search prefers the
tightest mirror first; `match_l[il]` is the upper pole currently
matched to lower pole `il` (`0` = free); `seen` guards against revisits
within one augmentation.  Returns `true` and rewires `match_l` when an
augmenting path is found from `u`, `false` otherwise.  Standard
textbook augmenting-path recursion — no project-specific subtlety.
"""
function _vc5_augment!(u::Int, adj::Vector{Vector{Int}},
                       match_l::Vector{Int}, seen::Vector{Bool})
    for il in adj[u]
        seen[il] && continue
        seen[il] = true
        if match_l[il] == 0 ||
           _vc5_augment!(match_l[il], adj, match_l, seen)
            match_l[il] = u
            return true
        end
    end
    return false
end

"""
    vc5_pair(poles) -> NamedTuple

Match the VC-4-surviving poles into conjugate pairs.  `V₀(x̄) =
conj(V₀(x))`, so a genuine off-axis pole `p` has a mirror partner
`conj(p)`; the pole field must be symmetric under conjugation.

**Globally-optimal assignment under `x ↦ x̄` (ADR-0025 Amendment 6).**
Build the bipartite *admissibility graph* — an edge `(p_upper,
p_lower)` whenever the conjugate residual `|p_upper - conj(p_lower)| ≤
VC5_MATCH_TOL` — and compute a **maximum-cardinality matching** of it by
Kuhn augmenting paths (`_vc5_augment!`).  This is the genuine fix the
VC-5b investigation (bead `padetaylor-0ln.37.20`) identified: the v1
`vc5_pair` used a *globally-greedy* commit (tightest admissible pair
first, each pole consumed once), and greedy is **not** maximum-
cardinality — in a dense field it commits a locally-tight pair that
*blocks* two other poles from each finding their only admissible
partner.  The VC-5b probe measured the gap directly: on the 266-pole B3
field, greedy finds 93 pairs where the optimal matching finds 103.
Maximising cardinality is exactly "pair up every pole that genuinely
*has* a mirror partner", so it minimises the spurious-flag count.

Each upper pole's admissible list is sorted by ascending residual, so
among equal-cardinality matchings the augmenting search prefers the
tightest edges — the reported residual distribution stays the FW-style
accuracy estimate (FW 2011 `:303-311`), not an arbitrary one.

A pole left unpaired has **no** admissible mirror — its conjugate
partner is farther than `VC5_MATCH_TOL` (above the 0.69 pole spacing
this could not happen for a genuine pair; see the `VC5_MATCH_TOL`
docstring).  An unpaired pole must therefore lie *on* the real axis
(`|Im| < VC5_REAL_AXIS_TOL` — it is its own conjugate); an unpaired
pole with a larger `|Im|` is **flagged** suspect.

The matching is **not** a conjugate-symmetry *construction*: it pairs
poles the two independent half-fields *separately extracted* and
reports the residual — a flagged pole stays flagged, never synthesised
from its mirror.  VC-5 therefore remains a genuine accuracy cross-check
(ADR-0025 Amendment 6; FFW 2017 `:120-124`).

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
    nu, nl = length(upper), length(lower)

    # The admissibility graph: an edge wherever the conjugate residual
    # is within VC5_MATCH_TOL.  Each upper pole's adjacency list is
    # sorted by ascending residual so the augmenting search prefers the
    # tightest mirror — keeping the residual distribution honest.
    adj = [Int[] for _ in 1:nu]
    for iu in 1:nu
        pu = upper[iu]
        for il in 1:nl
            abs(pu - conj(lower[il])) ≤ VC5_MATCH_TOL && push!(adj[iu], il)
        end
        sort!(adj[iu]; by = il -> abs(pu - conj(lower[il])))
    end

    # Maximum-cardinality matching by Kuhn augmenting paths.
    match_l = zeros(Int, nl)              # lower index → matched upper, 0=free
    for u in 1:nu
        _vc5_augment!(u, adj, match_l, fill(false, nl))
    end

    pairs     = Tuple{ComplexF64,ComplexF64}[]
    residuals = Float64[]
    for il in 1:nl
        u = match_l[il]
        u == 0 && continue
        push!(pairs, (upper[u], lower[il]))
        push!(residuals, abs(upper[u] - conj(lower[il])))
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
