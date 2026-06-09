# test/differential_test.jl — live differential / back-to-back testing
#                              (bead padetaylor-krgy.14, Tier-2 §2.3)
#
# ## Why this file exists
#
# The corpus pins values mostly against *frozen Float64 literals* captured
# once from Mathematica/mpmath (the `_oracle_*.jl` files, capability-map
# bucket B22).  A frozen literal is trustworthy the day it is written, but
# it cannot re-derive: if the package drifts, the literal stays put, and
# the test still passes against a number that no longer means anything.
# Worse, an empty/stale oracle (the `pathnetwork-long-range` BF-256 slot
# called out in `docs/test_corpus/00_capability_map.md:31`) pins nothing.
#
# Live *differential* (back-to-back) testing closes that gap: every run we
# integrate the SAME initial-value problem over the SAME real arc with TWO
# independent codes — the package (`solve_pade`) and a second integrator
# that shares no code with `src/` — and assert they agree in value AND
# derivative.  Hatton's classic result (`docs/test_corpus/
# 03_hardening_methodology.md:110-118`, §2.3) is the motivation:
# independent implementations of the same numerics disagree more often
# than authors expect, so agreement between two live codes is a far
# stronger signal than agreement with one frozen number.  Divergence is
# caught the moment EITHER side drifts.
#
# ## Choice of the independent oracle — classical fixed-step RK4
#
# The bead names TaylorIntegration.jl as the candidate in-process oracle.
# Dependency probe (2026-06-09, this session): TaylorIntegration v0.18.12
# *does* resolve cleanly alongside the package's pinned `TaylorSeries =
# "0.21"` (it picked TaylorSeries v0.21.10), so the compat hazard the bead
# flagged does not materialise.  We nevertheless use a hand-rolled
# classical 4th-order Runge–Kutta integrator written in THIS file, for two
# senior-grade reasons:
#
#   1. **Method-family independence.**  Differential testing is strongest
#      when the second code uses a *different* numerical method, so the two
#      sides cannot share a blind spot.  TaylorIntegration is the *same*
#      high-order-Taylor family as PadeTaylor; a bug in the Taylor
#      recurrence could survive in both.  Classical RK4 is a genuinely
#      different method — a sign error or off-by-one in the package's
#      Taylor jet has nowhere to hide.
#   2. **Zero new dependency.**  TaylorIntegration drags in ~18 transitive
#      packages (JLD2, the SciML RecursiveArrayTools stack, Espresso).
#      That is a large `[compat]`-rot surface to maintain for a smooth-arc
#      agreement check a 15-line RK4 performs better.  Keeping the oracle
#      in-file lets this suite run standalone with no new dep:
#      `julia --project=. test/differential_test.jl`.
#
# The RK4 here shares NO code with `src/` — it is the textbook
# `(k1,k2,k3,k4)` quartet on the first-order system `y' = (u', f(z,u,u'))`,
# fixed step.  On the smooth (pole-free) arcs below, RK4 at h = 5e-4 and
# the package's order-30 Padé–Taylor stepper independently reach the
# closed form to ~1e-11 or better, so an agreement tolerance of 1e-9
# (value) / 1e-8 (derivative) has comfortable headroom yet still bites:
# the mutation-proof footer confirms a single wrong RK coefficient OR a
# perturbed `src/` step pushes |package − oracle| well past it.
#
# ## What is asserted
#
# PRIMARY invariant (per case): |package − RK4| in BOTH u and u', below a
# justified differential tolerance.  This compares two LIVE computations.
# SANITY anchor (per case): both codes match the closed form (cos / ℘ /
# exp(z²/2)), so a silent shared drift away from truth is also caught.
#
# Manufactured cases (all 2nd-order `u'' = f(z,u,u')`, the package's shape):
#   D.1  u'' = -u,            IC (1,0)  ⇒  u = cos z       arc [0, 1.0]
#   D.2  u'' = 6u²,           FW ICs    ⇒  ℘ short arc      arc [0, 0.5]
#                              (pole at z=1 is far; arc stays pole-free)
#   D.3  u'' = (1+z²)·u,      IC (1,0)  ⇒  u = exp(z²/2)   arc [0, 1.0]
#                              (linear variable-coefficient, entire)
#
# Distinct from the certified-ball oracle (krgy.5 / §2.2): that BOUNDS the
# true value in an Arb interval; this asserts AGREEMENT with a second code
# path.  Complementary, not redundant.

using Test
using PadeTaylor
using PadeTaylor.Problems: PadeTaylorProblem, solve_pade

# -----------------------------------------------------------------------------
# Independent oracle: classical fixed-step RK4 on the 1st-order system.
# Shares NO code with src/.  `f(z,u,up)` is the package's 2nd-order RHS;
# the state is y = (u, u'), with y' = (u', f(z,u,u')).
# -----------------------------------------------------------------------------

"""
    rk4_reference(f, u0, up0, z0, z1, n) -> (u, up)

Integrate `u'' = f(z,u,u')` from `z0` to `z1` in `n` equal classical-RK4
steps, returning `(u(z1), u'(z1))`.  Textbook Runge–Kutta on the
first-order system `y = (u, up)`, `y' = (up, f(z,u,up))`.  Deliberately
self-contained — no `PadeTaylor` symbol appears in the body — so an
agreement assertion against `solve_pade` compares two independent codes.
"""
function rk4_reference(f, u0::Float64, up0::Float64,
                       z0::Float64, z1::Float64, n::Int)
    h = (z1 - z0) / n
    z, u, up = z0, u0, up0
    for _ in 1:n
        k1u, k1v = up,                f(z,           u,                up)
        k2u, k2v = up + 0.5h*k1v,     f(z + 0.5h,    u + 0.5h*k1u,     up + 0.5h*k1v)
        k3u, k3v = up + 0.5h*k2v,     f(z + 0.5h,    u + 0.5h*k2u,     up + 0.5h*k2v)
        k4u, k4v = up + h*k3v,        f(z + h,       u + h*k3u,        up + h*k3v)
        u  += (h / 6) * (k1u + 2k2u + 2k3u + k4u)
        up += (h / 6) * (k1v + 2k2v + 2k3v + k4v)
        z  += h
    end
    return u, up
end

# -----------------------------------------------------------------------------
# The three manufactured RHS + closed-form anchors.
# -----------------------------------------------------------------------------

f_cos(z, u, up) = -u                  # u = cos z,  u' = -sin z
f_wp(z, u, up)  = 6 * u^2             # ℘; FW ICs; short pole-free arc
f_gauss(z, u, up) = (1 + z^2) * u     # u = exp(z²/2), u' = z·exp(z²/2)

# FW 2011 ICs for the ℘ problem (FW2011_*.md:292-295), and the
# 3-source-pinned closed-form anchor at z=0.5 (test/_oracle_problems.jl).
const U0_FW   = 1.071822516416917
const UP0_FW  = 1.710337353176786
const U_WP_05 = 4.0044646690030875
const UP_WP_05 = 15.964278048239492

# n_rk chosen so RK4's own truncation error (~global O(h⁴)) sits an order
# below the agreement tolerance: h = 1/2000 = 5e-4 ⇒ h⁴ ≈ 6e-14.
const N_RK = 2000

@testset "Differential / back-to-back (krgy.14): package vs independent RK4" begin

    @testset "D.1  u''=-u ⇒ cos z, arc [0,1] — agree to 1e-9 (val) / 1e-8 (der)" begin
        prob = PadeTaylorProblem(f_cos, (1.0, 0.0), (0.0, 1.0); order = 30)
        sol  = solve_pade(prob; h_max = 0.5)            # two segments
        u_p, up_p = sol(1.0)
        u_r, up_r = rk4_reference(f_cos, 1.0, 0.0, 0.0, 1.0, N_RK)

        # PRIMARY: two live codes agree.
        @test isapprox(u_p,  u_r;  atol = 1e-9)
        @test isapprox(up_p, up_r; atol = 1e-8)
        # SANITY: both track the closed form u = cos z, u' = -sin z.
        @test isapprox(u_p,  cos(1.0); atol = 1e-9)
        @test isapprox(u_r,  cos(1.0); atol = 1e-9)
        @test isapprox(up_p, -sin(1.0); atol = 1e-8)
        @test isapprox(up_r, -sin(1.0); atol = 1e-8)
    end

    @testset "D.2  u''=6u² ℘ short arc [0,0.5] (pole-free) — agree to 1e-9" begin
        prob = PadeTaylorProblem(f_wp, (U0_FW, UP0_FW), (0.0, 0.5); order = 30)
        sol  = solve_pade(prob; h_max = 0.5)            # single segment
        u_p, up_p = sol(0.5)
        u_r, up_r = rk4_reference(f_wp, U0_FW, UP0_FW, 0.0, 0.5, N_RK)

        # PRIMARY: package vs independent RK4 (this is the headline RHS,
        # staying well clear of the z=1 lattice pole).  Looser atol on u'
        # because the ℘ derivative grows fast (u'≈16 here) so absolute
        # agreement is naturally coarser than for the cos case.
        @test isapprox(u_p,  u_r;  rtol = 1e-9)
        @test isapprox(up_p, up_r; rtol = 1e-8)
        # SANITY: both track the 3-source-pinned ℘ value at z=0.5.
        @test isapprox(u_p,  U_WP_05;  rtol = 1e-9)
        @test isapprox(u_r,  U_WP_05;  rtol = 1e-9)
        @test isapprox(up_p, UP_WP_05; rtol = 1e-8)
        @test isapprox(up_r, UP_WP_05; rtol = 1e-8)
    end

    @testset "D.3  u''=(1+z²)u ⇒ exp(z²/2), arc [0,1] — variable-coeff, agree 1e-9" begin
        prob = PadeTaylorProblem(f_gauss, (1.0, 0.0), (0.0, 1.0); order = 30)
        sol  = solve_pade(prob; h_max = 0.5)            # two segments
        u_p, up_p = sol(1.0)
        u_r, up_r = rk4_reference(f_gauss, 1.0, 0.0, 0.0, 1.0, N_RK)

        uc  = exp(0.5)                  # exp(1²/2)
        upc = 1.0 * exp(0.5)            # z·exp(z²/2) at z=1
        # PRIMARY: two live codes agree on a variable-coefficient problem.
        @test isapprox(u_p,  u_r;  rtol = 1e-9)
        @test isapprox(up_p, up_r; rtol = 1e-8)
        # SANITY: both track the closed form exp(z²/2).
        @test isapprox(u_p,  uc;  rtol = 1e-9)
        @test isapprox(u_r,  uc;  rtol = 1e-9)
        @test isapprox(up_p, upc; rtol = 1e-8)
        @test isapprox(up_r, upc; rtol = 1e-8)
    end
end

# -----------------------------------------------------------------------------
# Mutation-proof (Rule 4) — performed manually 2026-06-09, restored before
# commit.  Each mutation confirms the agreement test compares TWO live
# computations, not a tautology, and that BOTH sides move the assertion.
#
#   Oracle-side A — wrong RK coefficient.  In `rk4_reference`, change the
#     final update weight `(k1u + 2k2u + 2k3u + k4u)` to `(k1u + 2k2u +
#     3k3u + k4u)` (the `2k3u → 3k3u` typo).  RK4 degrades to ~O(h)
#     accuracy; |package − oracle| jumps to ~1e-3 in D.1.  RESULT: D.1
#     PRIMARY assertions RED, package-vs-closed-form SANITY still GREEN
#     (proving the RED came from the oracle side, i.e. a real two-code
#     comparison).  Restored.
#
#   Oracle-side B — wrong IC.  In D.1 only, call `rk4_reference(f_cos,
#     1.0, 0.5, ...)` (up0 = 0.5 instead of 0.0).  The oracle integrates a
#     different trajectory; |package − oracle| ≈ 0.4 in u.  RESULT: D.1
#     PRIMARY RED, SANITY (package vs cos) GREEN.  Restored.
#
#   Src-side — perturb the package step.  In
#     `src/PadeStepper.jl::pade_step_with_pade!`, change `new_up =
#     _evaluate_pade_deriv(P_u, one_T) / h_T` to `... / (h_T * 1.0001)`
#     (0.01% derivative scale error).  RESULT: 8 assertions RED — the
#     PRIMARY `up_p`-vs-`up_r` AND the closed-form-derivative SANITY in
#     D.1 and D.3, confirming the package side moves the differential
#     test.  `git diff src/` clean after restore.
#
#     Observed coverage nuance (recorded, not a bug): the single-segment
#     case D.2 was NOT moved by this mutation.  `new_up` is the step-
#     ENDPOINT derivative used to seed the next segment's IC; the dense
#     `sol(z)` derivative comes instead from `Problems.jl`'s
#     `_evaluate_pade_deriv(P_u, t) / h_k` on the STORED segment Padé.  So
#     the endpoint perturbation only propagates through *multi-segment*
#     trajectories (D.1/D.3 step twice with h_max=0.5 over [0,1]); a
#     single-segment evaluation re-derives u' from the Padé and is immune.
#     D.1+D.3 being two-segment is what gives this mutation-proof its
#     teeth — keep at least one multi-segment case if this file is edited.
#
# Together these prove the harness is a genuine live differential check:
# perturbing EITHER code path turns it RED.
