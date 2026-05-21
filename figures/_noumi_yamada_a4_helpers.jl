# figures/_noumi_yamada_a4_helpers.jl
#
# Shared pure-Julia kernel for the A_4^(1) Noumi–Yamada pole-field
# figure (bead `padetaylor-0ln.17`, v0.2 plan row V8a).  The figure
# script `figures/noumi_yamada_a4_pole_field.jl` and the acceptance
# test `test/noumi_yamada_a4_figure_test.jl` both `include` this file
# so they exercise the **exact same computation** — the
# `vector_path_network_solve` + `extract_poles_shared_q` pipeline —
# rather than two drifting copies.  This mirrors the
# `figures/_ffw2017_fig_7_helpers.jl` pattern: keep the load-bearing
# numerics in one Makie-free file the test can load without pulling
# CairoMakie into the package test environment.
#
# ## What the kernel fixes
#
# Every quantity that defines "the figure" — the A_4^(1) parameters
# `α`, the all-nonzero initial condition `f0`, the integration window,
# the complex-`t`-plane target region, the Taylor order, the wedge
# step `h`, and the pole-extraction filter parameters — lives here as
# a named constant.  The figure and the test read those constants;
# neither hard-codes a duplicate.  A mutation to any constant (the
# test's mutation-proof footer corrupts `NY_ALPHA` and the RHS sign)
# is therefore seen identically by both call sites.
#
# ## The A_4^(1) system and the chosen problem
#
# `A_4^(1)` is the even-parity Noumi–Yamada system with `n = 2`,
# `l = 4`, five components `f_0,…,f_4` and five parameters
# `α_0,…,α_4` summing to `1` (pillar B §1, §1.2;
# `docs/v0p2_pillarB_noumi_yamada_findings.md`).  The trivial rational
# solutions (pillar B §5.2 — Types A/B/C) have *no* finite poles; a
# *generic* (non-rational) solution does, and that pole field is what
# this figure maps.
#
# The chosen problem:
#
#   * `α = (0.30, 0.10, 0.25, 0.15, 0.20)` — a *generic* parameter
#     set, all five distinct, summing to exactly `1` (the `k = 1`
#     normalisation).  Real ⇒ the system's RHS is conjugate-symmetric
#     (`f(t̄) = f̄(t)` solves it whenever `f(t)` does — verified
#     exactly in the test, invariant `NYF.1.5`).
#
#   * `f0 = (0.7, -0.3, 0.5, -0.55, -0.35)` — an **all-nonzero**
#     initial condition (the bead's V5b warning: a component
#     identically zero triggers the shared-`Q` degeneracy).  The
#     fifth component is set to `-Σ(first four)` so `Σf0 = 0`
#     *exactly* — the constraint `Σf_j = t` at `t0 = 0`.
#
#   * region: the square `[-2.5, 2.5]²` of the complex `t`-plane,
#     a `0.5`-spaced raster — the walk reaches it in ~190 visited
#     nodes and the extracted pole field is non-trivial (~26 poles,
#     a scattered cloud, none of them the rational-solution
#     pole-free degeneracy).
#
# These choices were arrived at empirically (probe `/tmp/probe_a4*`):
# the all-nonzero IC avoids the V5b degeneracy, and the `0.5`-raster
# over `[-2.5,2.5]²` is the largest region where the wedge walk
# completes cleanly while still showing a genuinely non-trivial
# pole cloud.
#
# ## Why the shared-Q roots ARE the system's pole field
#
# `vector_pade_step_with_pade!` fits all five components of the
# A_4^(1) state with one **shared** denominator `Q` (ADR-0019,
# type-II Hermite–Padé).  A vector meromorphic solution's components
# all blow up together at the same movable pole, so the roots of the
# single shared `Q` are *every* component's poles at once — the
# system's pole field, one consistent estimate.  `extract_poles_
# shared_q` roots each visited node's `Q`, maps roots to the
# `t`-plane, and cross-node-clusters; that is the pole field this
# figure scatters.
#
# References:
#   * `src/NoumiYamada.jl` — `NoumiYamadaProblem`, `noumi_yamada_rhs`.
#   * `src/VectorPathNetwork.jl` — `vector_path_network_solve`.
#   * `src/VectorPoleField.jl` — `extract_poles_shared_q`.
#   * `docs/v0p2_pillarB_noumi_yamada_findings.md` §1, §5.2 — the
#     A_4^(1) system + its rational (pole-free) solutions.
#   * `docs/v0p2_plan.md` row V8a — the PRD `n ≥ 4` acceptance item.

using PadeTaylor
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution
using PadeTaylor.VectorPoleField:   extract_poles_shared_q

# ----------------------------------------------------------------------
# The figure-defining constants — fixed here, read by both call sites.
# ----------------------------------------------------------------------

# A_4^(1): n = 2, d = 2n+1 = 5 components / parameters.
const NY_N = 2
const NY_D = 2 * NY_N + 1

# Generic parameter set: five distinct values summing to exactly 1.
const NY_ALPHA = ComplexF64[0.30, 0.10, 0.25, 0.15, 0.20]

# All-nonzero initial condition (V5b: no component may be identically
# zero).  The fifth entry is pinned to `-Σ(first four)` so that the
# constraint Σf0 = t0 = 0 holds *exactly* (not just to sqrt(eps)).
const NY_F0_HEAD = ComplexF64[0.7, -0.3, 0.5, -0.55]
const NY_F0 = vcat(NY_F0_HEAD, ComplexF64[-sum(NY_F0_HEAD)])

# Integration window start point t0 (= Σf0) and the nominal span end.
const NY_T0   = 0.0 + 0.0im
const NY_TEND = 3.0 + 0.0im

# Taylor truncation order + wedge step magnitude.
const NY_ORDER = 24
const NY_H     = 0.3

# Complex-t-plane target region: square [-R, R]^2, RASTER-spaced.
const NY_REGION_R    = 2.5
const NY_REGION_STEP = 0.5

# Pole-extraction filter parameters (mirror VectorPoleField defaults,
# with a slightly tighter cluster radius for the A_4 pole spacing).
const NY_RADIUS_T     = 5.0
const NY_CLUSTER_ATOL = 0.15
const NY_MIN_SUPPORT  = 2

# ----------------------------------------------------------------------
# Kernel
# ----------------------------------------------------------------------

"""
    a4_target_grid() -> Vector{ComplexF64}

The complex-`t`-plane target raster for the Stage-1 walk: a
`NY_REGION_STEP`-spaced square lattice over `[-NY_REGION_R,
NY_REGION_R]²`.  Conjugate-symmetric as a *set* (every `t` has its
conjugate `t̄` in the grid) — the symmetric sampling region the
conjugate-symmetric A_4^(1) problem deserves.
"""
function a4_target_grid()
    R, s = NY_REGION_R, NY_REGION_STEP
    return ComplexF64[x + y * im for x in -R:s:R for y in -R:s:R]
end

"""
    a4_problem() -> NoumiYamadaProblem

Build the generic A_4^(1) `NoumiYamadaProblem` with the figure's fixed
`α`, all-nonzero `f0`, window and order.  The constructor validates
`Σα = 1` and `Σf0 = t0` (Rule 1).
"""
a4_problem() = NoumiYamadaProblem(NY_N; α = NY_ALPHA, f0 = NY_F0,
                                  tspan = (NY_T0, NY_TEND),
                                  order = NY_ORDER)

"""
    a4_walk() -> VectorPathNetworkSolution

Run the Stage-1 vector path-network walk of the generic A_4^(1)
problem over the figure's target region.  Deterministic — a fixed
walk produces the same visited tree every call (PRD "reproducible
from figures/").
"""
function a4_walk()
    prob = a4_problem()
    # Explicit walk-policy pin (ADR-0025 Amendment 4, bead `padetaylor-s8z`).
    # B2 made `step_policy=:max_q_root` + `adaptive=true` the library
    # default; the V8a A_4⁽¹⁾ figure was built and validated with the V7
    # `:min_y` fixed-`h` walk, and its order-24 shared-`Q` SVD degenerates
    # at node states the `:max_q_root` path reaches.  Pinning the policy
    # reproduces the validated A_4 figure exactly and makes this figure's
    # walk independent of the mutable library default.
    return vector_path_network_solve(prob.problem, a4_target_grid();
                                     order = NY_ORDER, h = NY_H,
                                     step_policy = :min_y, adaptive = false)
end

"""
    a4_pole_field(sol) -> Vector{ComplexF64}

Extract the A_4^(1) pole field from a solved walk: the roots of the
per-node shared denominators `Q`, mapped to the `t`-plane and
cross-node-clustered.  Because `Q` is shared across all five
components, these are every component's poles at once.
"""
a4_pole_field(sol::VectorPathNetworkSolution) =
    extract_poles_shared_q(sol; radius_t = NY_RADIUS_T,
                           cluster_atol = NY_CLUSTER_ATOL,
                           min_support = NY_MIN_SUPPORT)

"""
    a4_constraint_residual(sol) -> Float64

Largest violation of the A_4^(1) flow invariant `Σf_j(t) = t` over
every visited node of the walk: `maxₖ |Σⱼ y[k][j] − z[k]|`.  The
constraint is an invariant of the flow (pillar B §1.2); a faithful
walk keeps it near zero.
"""
function a4_constraint_residual(sol::VectorPathNetworkSolution)
    worst = 0.0
    for k in eachindex(sol.visited_z)
        worst = max(worst, abs(sum(sol.visited_y[k]) - sol.visited_z[k]))
    end
    return worst
end
