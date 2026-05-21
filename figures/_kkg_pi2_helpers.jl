# figures/_kkg_pi2_helpers.jl
#
# Shared pure-Julia kernel for the KKG 2015 P_I^(2) tritronquée
# pole-field figure (bead `padetaylor-0ln.18`, v0.2 plan row V8b).
# The figure script `figures/kkg_pi2_tritronquee_pole_field.jl` and the
# acceptance test `test/kkg_pi2_figure_test.jl` both `include` this file
# so they exercise the **exact same computation** — the
# `vector_bvp_solve` BVP solve followed by the
# `vector_path_network_solve` + `extract_poles_shared_q` pole-field
# march — rather than two drifting copies.  This mirrors the
# `figures/_noumi_yamada_a4_helpers.jl` pattern: the load-bearing
# numerics live in one Makie-free file the test can `include` without
# pulling CairoMakie into the package test environment.
#
# ## The equation and the solution
#
# P_I^(2) is the fourth-order member of the Painlevé-I hierarchy
# (KKG 2015 eq. (1.1); pillar C §1,
# `docs/v0p2_pillarC_painleve_hierarchy_findings.md`):
#
#     u'''' + 10 u'^2 + 20 u u'' + 40 (u^3 - 6 t u + 6 x) = 0,
#
# integrated here at the parameter `t = 0`.  Its **tritronquée**
# solution `V_0` (KKG 2015 §7) is the distinguished solution that is
# pole-free over a sector of angular width `≈ 270°` of the complex
# `x`-plane (Claeys 2011; pillar C §"Pole-free sector") and carries a
# pole field only in the remaining wedge straddling the positive real
# axis (`arg x ≈ 0`).  On the **negative real axis** it is real,
# positive and monotone, with the algebraic asymptotics
# `u → +(6|x|)^{1/3}` as `x → -∞` (KKG 2015 eq. (1.3); pillar C item 4).
#
# ## How the figure is computed — two stages
#
# **Stage A — the tritronquée on the negative real axis (BVP).**
# Forward IVP integration of a tritronquée diverges: it is a
# separatrix, one mode grows exponentially.  A global Chebyshev
# collocation BVP is therefore mandatory (`src/VectorBVP.jl`).  We
# solve the companion 4-vector system on a negative-real-axis segment
# `[x_l, x_r]` with a 2+2 split boundary condition pinning `(u, u')`
# at *both* endpoints from the KKG asymptotic seed
# `pI2_tritronquee_ic`.  This is the recipe verified by probe
# `padetaylor-0ln.29` (`external/probes/pi2-tritronquee-bvp/RECIPE.md`):
# the plain `vector_bvp_solve` converges in 2–3 Newton iterations,
# ODE residual `~1.5e-11`, matching the KKG truncated asymptotic series
# to `~2.1e-4` on `[-20,-2]`.
#
# IMPORTANT — the sign.  `pI2_tritronquee_ic` was sign-corrected (the
# tritronquée on `x<0` has `u > 0`, from the dominant balance
# `40(u^3+6x)=0 ⇒ u^3 = 6|x| > 0`).  We use the corrected library
# `pI2_tritronquee_ic` directly — never a hand-rolled seed.
#
# **Stage B — the pole-field march.**  The pole field of `V_0` lives
# in the wedge near the positive real axis.  We seed a vector
# path-network walk from a BVP-derived point near the right end of the
# negative-real segment — `y_seed = bvp_sol(z_seed)` — build a
# `VectorPadeTaylorProblem` with that `(z_seed, y_seed)` as initial
# condition, and march the 5-direction wedge walk toward targets
# covering the pole-rich wedge straddling `arg x ≈ 0`.
# `extract_poles_shared_q` then roots the per-node shared denominators
# `Q` and cross-node-clusters: because `Q` is shared across all four
# companion components, its roots are every component's poles at once
# (ADR-0019; `src/VectorPoleField.jl`).
#
# HONESTY NOTE (CLAUDE.md Rule 1, Rule 9).  Stage A is solid — the BVP
# recipe is proven.  Stage B is the uncertain part: marching the vector
# path network from a BVP point into the pole wedge is the V8b open
# risk.  `kkg_pole_field` is written to NOT throw on a march failure —
# it catches the breakdown and returns an empty pole list with a
# recorded `march_ok = false` flag, so the figure can ship Panel A
# solidly and report Panel B honestly rather than fabricate a pole
# field.  The full KKG Fig 7.4/7.5 Re/Im *surface* over `[-20,20]^2`
# is a separately-deferred v1 corner (needs Stage-2 fine-grid fill —
# bead `padetaylor-0ln.28`); this figure scatters the pole *locations*
# only, not a filled surface.
#
# References:
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015: eq. (1.1) the ODE, §7 the numerical tritronquée study,
#     Figs. 7.4–7.5 the Re/Im surface this figure's pole scatter
#     distils.
#   * `external/probes/pi2-tritronquee-bvp/RECIPE.md` — the verified
#     2+2 (u,u') BVP recipe and its numbers.
#   * `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — §1 the
#     equation, §4 the (sign-corrected) tritronquée IC, the pole-free
#     sector / pole-wedge structure.
#   * `src/PainleveHierarchy.jl` — `painleve_hierarchy`,
#     `painleve_hierarchy_jacobian`, `pI2_tritronquee_ic`,
#     `PainleveHierarchyProblem`, the `vector_bvp_solve` forwarding.
#   * `src/VectorBVP.jl`, `src/VectorPathNetwork.jl`,
#     `src/VectorPoleField.jl` — the three solver modules.

using PadeTaylor
using PadeTaylor.PainleveHierarchy: painleve_hierarchy,
                                    painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP:        vector_bvp_solve, VectorBVPSolution
using PadeTaylor.VectorProblems:   VectorPadeTaylorProblem
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution
using PadeTaylor.VectorPoleField:  extract_poles_shared_q

# ----------------------------------------------------------------------
# The figure-defining constants — fixed here, read by both call sites.
# ----------------------------------------------------------------------

# P_I^(2) parameter.  n_terms = 2 of pI2_tritronquee_ic is implemented
# only for t = 0, and the KKG tritronquée figure is the t = 0 case.
const KKG_T = 0.0

# BVP segment on the negative real axis.  [x_l, x_r] with x_l the
# large-|x| end (smaller truncation error in the asymptotic seed).
# [-20, -2] is the probe's headline segment (RECIPE.md "Converged
# solution"); the element type is complex so the BVP solution and the
# downstream path-network share one element type.
const KKG_X_L = -20.0 + 0.0im
const KKG_X_R =  -2.0 + 0.0im

# BVP discretisation.  N+1 Chebyshev collocation nodes; the probe shows
# the spectral solve is converged for any N ≥ 40, identical to ~1e-7.
const KKG_BVP_N       = 128
const KKG_BVP_TOL     = 1.0e-9
const KKG_BVP_MAXITER = 40

# Companion-system order: d = 2m = 4 for P_I^(2) (m = 2).
const KKG_D = 4

# Stage-A sample raster on the negative real axis (the Panel-A curve).
const KKG_CURVE_N = 200

# ---- Stage-B pole-field march parameters -----------------------------
# The march seeds from a BVP-derived point near the RIGHT end of the
# negative-real segment and walks toward the pole-rich wedge straddling
# arg x ≈ 0 (the positive real axis).  z_seed sits just inside the BVP
# segment so bvp_sol(z_seed) is a genuine on-solution state vector.
const KKG_Z_SEED = -3.0 + 0.0im

# Path-network IVP problem window + Taylor order + wedge step.
#
# The wedge step `h` is the load-bearing Stage-B parameter.  The
# tritronquée's pole field is *dense* in the wedge near `arg x ≈ 0`,
# so a coarse step overshoots a pole: the node lands where the state
# `‖y‖` has already blown up to `O(100…1000)`, and
# `VectorPathNetwork._canonical_pade`'s per-node shared-`Q` store then
# throws (`shared_denominator_pade`: "the Taylor jets are
# indistinguishable from zero" — the jet at a blown-up state is
# numerically degenerate).  Empirically (probe `/tmp/probe_full.jl`):
# `h = 0.3` (the FW default / the A_4 figure's value) throws;
# `h = 0.15` throws; `h = 0.1` completes — 389 visited nodes, 21 poles
# in the wedge; `h = 0.05` completes but over-resolves (each tiny
# step's local Padé barely reaches a pole, so only ~6 poles clear the
# cross-node support filter).  `h = 0.1` is the chosen value.  This is
# a Rule-9 v1 corner: a principled fix — having the wedge selector
# avoid the *shared-Q root* nearest a candidate, not just minimise
# `‖y‖` — is the deferred refinement bead `padetaylor-0ln.23`; until
# then the pole-march `h` must be hand-tuned per target region.
const KKG_PN_ORDER = 24
const KKG_PN_H     = 0.1

# Target wedge: points in the complex x-plane covering the pole-rich
# region near the positive real axis (arg x ≈ 0).  A small fan of
# radii/angles straddling the real axis — the wedge KKG Fig 7.4 shows
# the sawtooth pole structure in.
const KKG_TARGET_RADII  = (2.0, 4.0, 6.0, 8.0)
const KKG_TARGET_ANGLES = (-0.5, -0.25, 0.0, 0.25, 0.5)  # radians

# Pole-extraction filter parameters (mirror VectorPoleField defaults).
const KKG_RADIUS_T     = 5.0
const KKG_CLUSTER_ATOL = 0.2
const KKG_MIN_SUPPORT  = 2

# ----------------------------------------------------------------------
# Stage A — the tritronquée on the negative real axis via the BVP.
# ----------------------------------------------------------------------

"""
    kkg_bvp_solve() -> VectorBVPSolution

Solve the P_I^(2) tritronquée on the negative-real-axis segment
`[KKG_X_L, KKG_X_R]` via the probe-verified 2+2 (u,u') boundary-value
recipe (`external/probes/pi2-tritronquee-bvp/RECIPE.md`).

The companion 4-vector system `y = (u, u', u'', u''')` is collocated on
`KKG_BVP_N + 1` Chebyshev nodes.  The boundary condition pins `(u, u')`
at *both* endpoints from the KKG asymptotic seed
`pI2_tritronquee_ic(·; n_terms = 2)` (the sign-corrected library seed —
on `x < 0` the tritronquée has `u > 0`).  The interior initial guess is
the same asymptotic seed sampled at every node; the analytic
companion-form Jacobian is supplied.  Returns the converged
`VectorBVPSolution` (callable: `sol(z) -> Vector` of the 4 components).
"""
function kkg_bvp_solve()
    f  = painleve_hierarchy(:I, 2; t = KKG_T)
    Jf = painleve_hierarchy_jacobian(:I, 2; t = KKG_T)

    # 2+2 split boundary condition.  B_a pins (u, u') at x_l = y(z_a),
    # B_b pins (u, u') at x_r = y(z_b).  The τ-method endpoint↔node
    # pairing is handled inside vector_bvp_solve.
    CT = ComplexF64
    Ba = zeros(CT, KKG_D, KKG_D)
    Bb = zeros(CT, KKG_D, KKG_D)
    Ba[1, 1] = 1; Ba[2, 2] = 1          # rows 1-2: (u, u') at x_l
    Bb[3, 1] = 1; Bb[4, 2] = 1          # rows 3-4: (u, u') at x_r

    seed_l = pI2_tritronquee_ic(KKG_X_L; t = KKG_T, n_terms = 2)
    seed_r = pI2_tritronquee_ic(KKG_X_R; t = KKG_T, n_terms = 2)
    g = CT[seed_l[1], seed_l[2], seed_r[1], seed_r[2]]

    return vector_bvp_solve(f, KKG_X_L, KKG_X_R, Ba, Bb, g;
                            N        = KKG_BVP_N,
                            tol      = KKG_BVP_TOL,
                            maxiter  = KKG_BVP_MAXITER,
                            jacobian = Jf,
                            initial_guess =
                                z -> pI2_tritronquee_ic(z; t = KKG_T,
                                                        n_terms = 2))
end

"""
    kkg_tritronquee_curve(bvp_sol) -> (xs, us)

Sample the BVP-computed tritronquée on the negative real axis: `xs` is
the `KKG_CURVE_N`-point real raster over `[KKG_X_L, KKG_X_R]`, `us` is
`real(u(x))` (the first companion component) — the quantitative
Panel-A curve.  On the negative real axis the tritronquée is real and
positive (`u > 0`), so the imaginary part is discarded as round-off.
"""
function kkg_tritronquee_curve(bvp_sol::VectorBVPSolution)
    xs = range(real(KKG_X_L), real(KKG_X_R); length = KKG_CURVE_N)
    us = [real(bvp_sol(ComplexF64(x))[1]) for x in xs]
    return collect(xs), us
end

"""
    kkg_asymptotic(x) -> Float64

The KKG `n_terms = 2` truncated asymptotic series `u = Y + Y^{-6}` of
the tritronquée at the negative-real point `x < 0` — `Y = +(6|x|)^{1/3}`,
`Y^{-6} = |x|^{-2}/36` (KKG 2015 eq. (7.2), `c_6 = 1` at `t = 0`).  The
oracle the Panel-A curve and test `PI2F.1.2` compare the BVP solution
against.
"""
kkg_asymptotic(x) = cbrt(6.0) * abs(x)^(1 / 3) + abs(x)^(-2) / 36

"""
    kkg_ode_residual(bvp_sol, x; δ = 0.05) -> Float64

Genuine non-tautological P_I^(2) ODE residual at an interior point `x`.

The BVP solution carries the 4 companion components
`y = (u, u', u'', u''')` as *independent* barycentric interpolants — it
does NOT carry `u''''`.  So we form `u''''` two independent ways and
compare:

  - **companion / equation side** — the P_I^(2) ODE at `t = 0`,
    `u'''' = -(10 u'^2 + 20 u u'' + 40(u^3 + 6x))`, evaluated from the
    BVP state's `(u, u', u'')` components;
  - **collocation side** — a centred finite difference of the BVP
    solution's *third* component `u'''(x)`:
    `u''''(x) ≈ (u'''(x+δ) - u'''(x-δ)) / (2δ)`.

If the BVP solution is a genuine P_I^(2) solution these two `u''''`
estimates agree; their absolute difference is the residual returned.
This is non-tautological because the finite-difference side uses only
the `u'''` interpolant while the equation side uses only `(u,u',u'')`
— two disjoint slices of the state, reconciled solely *because* the
state solves the ODE.  `δ` is kept a few step-magnitudes inside the
segment so both `x ± δ` are in `[KKG_X_L, KKG_X_R]`.
"""
function kkg_ode_residual(bvp_sol::VectorBVPSolution, x; δ::Real = 0.05)
    y = bvp_sol(ComplexF64(x))
    u, up, upp, _ = y
    # Equation side: P_I^(2) at t = 0 solved for u''''.
    uxxxx_eq = -(10 * up^2 + 20 * u * upp + 40 * (u^3 + 6 * x))
    # Collocation side: centred difference of the u''' component.
    uppp_p = bvp_sol(ComplexF64(x + δ))[4]
    uppp_m = bvp_sol(ComplexF64(x - δ))[4]
    uxxxx_fd = (uppp_p - uppp_m) / (2δ)
    return abs(uxxxx_eq - uxxxx_fd)
end

# ----------------------------------------------------------------------
# Stage B — the pole-field march.
# ----------------------------------------------------------------------

"""
    kkg_target_wedge() -> Vector{ComplexF64}

The complex-`x`-plane target fan for the Stage-B pole-field walk: the
Cartesian product of `KKG_TARGET_RADII` and `KKG_TARGET_ANGLES`,
points `r·e^{iθ}` straddling the positive real axis (`arg x ≈ 0`) —
the pole-rich wedge of the P_I^(2) tritronquée (KKG Fig 7.4).
"""
function kkg_target_wedge()
    return ComplexF64[r * cis(θ) for r in KKG_TARGET_RADII
                                 for θ in KKG_TARGET_ANGLES]
end

"""
    kkg_pole_field(bvp_sol) -> NamedTuple

Stage-B pole-field march.  Seed a `VectorPadeTaylorProblem` for the
companion P_I^(2) system from the BVP-derived state
`y_seed = bvp_sol(KKG_Z_SEED)` at `z = KKG_Z_SEED` (a point inside the
negative-real BVP segment, so `y_seed` is a genuine on-tritronquée
state), run the vector path-network wedge walk toward `kkg_target_wedge`,
and `extract_poles_shared_q`.

Returns a NamedTuple `(march_ok, walk, poles, message)`:
  - `march_ok::Bool`  — `true` if the walk completed and pole
                        extraction ran;
  - `walk`            — the `VectorPathNetworkSolution` (or `nothing`
                        on failure);
  - `poles::Vector{ComplexF64}` — the extracted pole field (empty on
                        failure);
  - `message::String` — a human-readable note (the failure reason on
                        `march_ok = false`).

HONESTY (Rule 1 / Rule 9): the march is the V8b open risk — walking a
vector path network from a BVP point into the pole wedge may diverge.
This function therefore CATCHES a march breakdown and returns
`march_ok = false` with the reason, rather than throwing or fabricating
poles.  The figure script then ships Panel A solidly and reports
Panel B honestly.
"""
function kkg_pole_field(bvp_sol::VectorBVPSolution)
    f = painleve_hierarchy(:I, 2; t = KKG_T)
    y_seed = ComplexF64.(bvp_sol(ComplexF64(KKG_Z_SEED)))

    # IVP problem: companion P_I^(2), IC (z_seed, y_seed).  The zspan
    # end is a nominal point in the wedge direction; the path network
    # walks toward kkg_target_wedge, not a straight zspan integration.
    prob = VectorPadeTaylorProblem(f, y_seed,
                                   (ComplexF64(KKG_Z_SEED),
                                    ComplexF64(8.0 + 0.0im));
                                   order = KKG_PN_ORDER)

    targets = kkg_target_wedge()
    try
        walk = vector_path_network_solve(prob, targets;
                                         order = KKG_PN_ORDER,
                                         h     = KKG_PN_H)
        poles = extract_poles_shared_q(walk;
                                       radius_t     = KKG_RADIUS_T,
                                       cluster_atol = KKG_CLUSTER_ATOL,
                                       min_support  = KKG_MIN_SUPPORT)
        msg = "march OK — $(length(walk.visited_z)) visited nodes, " *
              "$(length(poles)) poles extracted"
        return (march_ok = true, walk = walk,
                poles = ComplexF64.(poles), message = msg)
    catch err
        msg = "march BLOCKED — $(typeof(err)): " *
              first(split(sprint(showerror, err), '\n'))
        return (march_ok = false, walk = nothing,
                poles = ComplexF64[], message = msg)
    end
end

# ----------------------------------------------------------------------
# Top-level driver — the full V8b computation in one struct.
# ----------------------------------------------------------------------

"""
    kkg_pi2_compute() -> NamedTuple

Run the full V8b computation and return everything the figure script
and the acceptance test need, in one NamedTuple:

  - `bvp`        — the `VectorBVPSolution` (Stage A);
  - `xs`, `us`   — the negative-real-axis tritronquée curve (Panel A);
  - `march_ok`   — whether the Stage-B march succeeded;
  - `walk`       — the `VectorPathNetworkSolution` (or `nothing`);
  - `poles`      — the extracted pole field (Panel B; empty on failure);
  - `message`    — the Stage-B status note.

Deterministic: a fixed BVP recipe + a fixed walk produce the identical
result on every call (the figure is a function of the helper constants
alone).
"""
function kkg_pi2_compute()
    bvp      = kkg_bvp_solve()
    xs, us   = kkg_tritronquee_curve(bvp)
    pf       = kkg_pole_field(bvp)
    return (bvp = bvp, xs = xs, us = us,
            march_ok = pf.march_ok, walk = pf.walk,
            poles = pf.poles, message = pf.message)
end
