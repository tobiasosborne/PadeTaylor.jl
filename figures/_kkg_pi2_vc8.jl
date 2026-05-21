# figures/_kkg_pi2_vc8.jl
#
# VC-8 — the BVP endpoint higher-derivative match of the KKG 2015
# P_I^(2) tritronquée headline figure (ADR-0025 Phase D, Amendment 8;
# bead `padetaylor-0ln.37.15`).  This helper is `include`d by the figure
# kernel `_kkg_pi2_surface_helpers.jl`; it is the Rule-6 file-size
# split, sibling to `_kkg_pi2_vc45.jl` (VC-4 / VC-5).
#
# ## What FW §5.2 does, and why the companion lift is the genuine invariant
#
# FW 2011 §5.2 (`references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:192-193`) solves the *scalar*
# second-order `P_I` BVP and validates the converged Chebyshev solution
# by estimating the endpoint derivatives via the first-derivative matrix
# `D₁` — "derivatives at the endpoints can be estimated via `D₁ u`" —
# and matching them against the pole-field solver's `u'`.  "If the
# agreement is not to within some tolerance (typically `10⁻⁷` or
# `10⁻⁸`), an increase in `N` … is indicated."
#
# The figure's BVP is *not* scalar: it is the `d = 4` first-order vector
# companion system `y' = f(z, y)` with `y = (u, u', u'', u''')`
# (`surf_ray_bvp`, `surf_anchor_bvp` in the kernel).  The companion
# structure asserts the **exact** identities
#
#     y[1]' = y[2],   y[2]' = y[3],   y[3]' = y[4]
#
# — the trivial chain rows of the companion form (`src/Painleve
# Hierarchy.jl` `_rhs_PI2`).  A *converged* Chebyshev-collocation
# solution must satisfy them: the spectral derivative of component `k`
# — the Chebyshev `D₁` applied to the node-vector `(y[k])ⱼ`, rescaled
# from the `t ∈ [-1,1]` collocation variable to `z` by `1/s`
# (`s = (z_b − z_a)/2`, the affine scale `VectorBVP` introduces) — must
# reproduce component `k+1`.  This is **exactly** FW's `D₁`-derivative
# diagnostic, and it is a genuine invariant: a collocation `N` too small
# to resolve the solution leaves a residual in `D₁·y[k] − y[k+1]` that
# exceeds FW's `~10⁻⁷·10⁻⁸` band — the "increase `N`" signal.  Crucially
# it owes nothing to the asymptotic seed, and — as the PI2S.8
# mutation-proof shows — it catches an under-resolution the BVP's own
# Newton *residual* does not: a starved-`N` solve converges to the
# spectral floor on its residual yet has a companion-consistency
# residual of O(1).  It is a self-consistency test of the spectral
# solution itself.
#
# We check the identity at the two **endpoint** nodes (FW's diagnostic
# is endpoint-specific — those are where the pole-field solver supplies
# `u'`) *and* report the worst mismatch over the *whole* node set, since
# a spectral solution consistent at the endpoints but not in the
# interior is equally under-resolved.
#
# ## The asymptotic-IC cross-check
#
# `VectorBVP` pins only `y[1] = u` and `y[2] = u'` at the endpoints (the
# figure's `2+2` split BC); the higher components `y[3] = u''`,
# `y[4] = u'''` are **free** — the collocation solve determines them.
# At the deep-asymptotic endpoint `z_a` (`|z_a| = 20`) the converged
# free components should agree with the `pI2_tritronquee_ic` asymptotic
# seed's own `y[3]`, `y[4]` to the seed's truncation accuracy
# (`n_terms = 2` is `O(10⁻⁵·10⁻⁶)` at `|x| = 20`).  This is a second,
# independent VC-8 facet — it certifies the BVP converged to the
# *tritronquée* branch the asymptotic seed selects, not merely to *a*
# solution of the companion system.
#
# References:
#   * `docs/adr/0025-headline-figure-re-resolution.md` — Amendment 8
#     (VC-8); the Validation Criteria Menu.
#   * `references/markdown/FW2011_painleve_methodology_JCP230/
#     FW2011_painleve_methodology_JCP230.md:192-193` — the FW §5.2
#     endpoint-derivative diagnostic.
#   * `src/VectorBVP.jl` — the vector Chebyshev-collocation BVP solver;
#     `_chebyshev_D1` (the first-derivative matrix) and the τ-method
#     descending-node-order convention this helper reads.
#   * `src/PainleveHierarchy.jl` — `_rhs_PI2` (the companion-form chain
#     identities) and `pI2_tritronquee_ic` (the asymptotic seed).
#   * `figures/_kkg_pi2_vc45.jl` — the sibling VC-4 / VC-5 helper this
#     file mirrors in shape.

"""
    surf_vc8_companion_check(sol::VectorBVPSolution) -> NamedTuple

VC-8 — the FW 2011 §5.2 endpoint-derivative diagnostic, adapted to the
`P_I^(2)` first-order vector **companion** BVP (ADR-0025 Amendment 8,
bead `padetaylor-0ln.37.15`).  See this file's header for the full
rationale; in brief, it certifies that a converged Chebyshev-collocation
solution of the `d = 4` companion system satisfies the exact chain
identities `y[k]' = y[k+1]` — the spectral derivative `D₁·y[k]/s` must
reproduce `y[k+1]` at the endpoints (FW §5.2's check) and throughout, to
FW's `~10⁻⁷·10⁻⁸` tolerance, else the collocation `N` is too small.

Returns `(consistency_max, consistency_endpoint, ic_y3_err, ic_y4_err,
residual_inf, N)`:
  - `consistency_max`      : max over all nodes and `k ∈ {1,2,3}` of
                             `|（D₁·y[k]）ⱼ/s − y[k+1]ⱼ|` — the global
                             companion-consistency residual;
  - `consistency_endpoint` : the same max restricted to the two
                             endpoint nodes (FW §5.2's endpoint check);
  - `ic_y3_err`, `ic_y4_err` : `|y[3] − seed y[3]|`, `|y[4] − seed y[4]|`
                             at the deep-asymptotic endpoint `z_a` —
                             the asymptotic-IC cross-check;
  - `residual_inf`         : the BVP's own converged Newton residual;
  - `N`                    : the collocation `N` (for the report).
"""
function surf_vc8_companion_check(sol::VectorBVPSolution)
    N  = sol.N
    s  = (sol.z_b - sol.z_a) / 2                 # affine scale (VectorBVP)
    D1 = PadeTaylor.VectorBVP._chebyshev_D1(sol.nodes_t, Float64, N)
    d  = length(sol.Y_nodes[1])

    consistency_max      = 0.0
    consistency_endpoint = 0.0
    for k in 1:d-1
        ck  = ComplexF64[y[k]   for y in sol.Y_nodes]
        ck1 = ComplexF64[y[k+1] for y in sol.Y_nodes]
        dck = (D1 * ck) ./ s                     # d y[k]/dz, spectral
        mism = abs.(dck .- ck1)
        consistency_max      = max(consistency_max, maximum(mism))
        consistency_endpoint = max(consistency_endpoint,
                                   mism[1], mism[end])
    end

    # The deep-asymptotic endpoint is node N+1 (t = -1 ↦ z_a, the
    # descending DMSUITE node order — VectorBVP's τ-method docstring).
    seed_a = pI2_tritronquee_ic(sol.z_a; t = SURF_T, n_terms = 2)
    y_za   = sol.Y_nodes[end]
    ic_y3_err = abs(y_za[3] - seed_a[3])
    ic_y4_err = abs(y_za[4] - seed_a[4])

    return (consistency_max      = consistency_max,
            consistency_endpoint = consistency_endpoint,
            ic_y3_err            = ic_y3_err,
            ic_y4_err            = ic_y4_err,
            residual_inf         = sol.residual_inf,
            N                    = N)
end
