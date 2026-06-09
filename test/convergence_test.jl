# test/convergence_test.jl — convergence-order / observed-order-of-accuracy
#                            V&V tests (bead padetaylor-krgy.6, Tier-2 §2.4)
#
# ## Why this file exists — code verification, not value pinning
#
# Roy (JCP 2005; cited in `docs/test_corpus/03_hardening_methodology.md:120-133`,
# §2.4) draws the V&V line between *solution verification* ("how big is the
# error in THIS calculation?") and *code verification* ("is the discretised
# math being solved at the rate the theory predicts?").  The corpus already
# pins many single-resolution VALUES (the `_oracle_*.jl` ledgers); what it
# lacked is the code-verification complement: refine the discretisation and
# confirm the error decays at the THEORETICAL order.  That is what catches
# order-degradation bugs — a mis-scaled recurrence term, a dropped Taylor
# coefficient, a wrong differentiation-matrix power — which a single-h value
# test sails straight past because the value is still "close enough" at one
# resolution.  We use the Method of Manufactured Solutions (MMS): pick a
# KNOWN smooth (pole-free) closed form, drive the package at it, and measure
# the OBSERVED order against theory.
#
# ## CALIBRATION — the load-bearing part of this test (read before editing)
#
# A V&V order test is only meaningful inside the window where the asymptotic
# error law actually governs.  Two boundaries bracket that window, and a
# sloppy choice of either makes the fitted slope a fiction:
#
#   * SATURATION (top).  For h too large the error is O(1) and the leading
#     `h^p` term does not dominate — the log-log slope is contaminated by
#     every higher-order term at once.  We discard any refinement point with
#     error above `ERR_HI = 1e-1`.
#   * FLOOR (bottom).  For h too small the true discretisation error drops
#     below the Float64 round-off / eps floor (~1e-16) and the measured error
#     is pure cancellation noise — the slope flattens or goes to ±Inf.  This
#     is the single calibration the bead NOTES flags explicitly: at order 30
#     the leading term is `h^31`, which is below eps for any usable h, so the
#     slope is UNMEASURABLE at the production order.  We must run at a MODEST
#     order and discard any point with error below `ERR_LO = 1e-12`.
#
# We fit only over points inside `[ERR_LO, ERR_HI]` and require ≥ 3 such
# points; the test states its valid window rather than trusting a fixed h-grid.
#
# ## CALIBRATION — the order law is NOT `h^(N+1)` for this package
#
# The bead text and the §2.4 framing both reach for the textbook one-step
# law `err ~ h^(N+1)` for a Taylor method of order N.  An empirical probe
# (this session) showed that is the WRONG theory for PadeTaylor's stepper,
# and asserting it would either fail spuriously or pass with a wrong band.
# The reason is structural and lives in `src/PadeStepper.jl`
# (`pade_step_with_pade!`): the step does NOT evaluate the order-N Taylor
# polynomial — it builds a DIAGONAL Padé approximant `(m, m)` with
# `m = order ÷ 2` from the rescaled jet and evaluates the rational at t = 1.
# A diagonal `(m, m)` Padé matches the Taylor series through degree `2m`, so
# its leading truncation term is `h^(2m+1)`.  Hence
#
#       observed single-step order  =  2·(order ÷ 2) + 1.
#
# Verified empirically (sin+cos and exp anchors):
#       order 2,3  → m=1 → slope 3   (= 2·1+1)
#       order 4,5  → m=2 → slope 5   (= 2·2+1)
#       order 6,7  → m=3 → slope 7   (= 2·3+1)
# The order is a knob and the slope tracks it through `m = order÷2` — the
# strongest V&V statement.  We test three orders giving three DISTINCT
# theoretical slopes (3, 5, 7) so a single mis-scaling cannot satisfy all.
#
# ## CALIBRATION — the manufactured solution must have BOTH parities
#
# A second probe surprise: an EVEN manufactured solution (cos z, exp(z²/2) —
# both even functions, odd Taylor coefficients vanish) makes the diagonal
# Padé degenerate.  `RobustPade`'s trailing-near-zero trimmer then collapses
# the requested `(3,3)` back to `(2,2)`, so the order STOPS tracking N past
# m=2 and the slope saturates at 5 regardless of order.  That would defeat
# the "slope tracks N" assertion.  We therefore anchor on `u = exp(z)`
# (`u''=u`, IC (1,1)) and `u = sin z + cos z` (`u''=-u`, IC (1,1)) whose
# Taylor series carry BOTH even and odd terms; for these the requested
# `(m,m)` survives the trim and the order tracks `m` exactly.  This subtlety
# is exactly why a senior order-test states which solution it uses and why.
#
# ## CV.2 — BVP spectral (exponential) convergence
#
# `bvp_solve` is Chebyshev spectral collocation; for an ANALYTIC solution
# the node error decays GEOMETRICALLY, `err ~ C·ρ^(-N)`, i.e. `log(err)`
# is LINEAR in N with a clearly-negative slope — faster than ANY algebraic
# rate.  Empirically (this session, `u''=u` → cosh, `u''=-u` → sin):
#       N= 4 → 1.7e-4,  N=8 → 7e-10,  N=12 → 3e-15   (log-decay ≈ -3.0 / unit N)
# then the error FLOORS at the spectral conditioning limit ≈ cond(D₂)·eps ≈
# N²·eps ≈ 1e-15 for N ≥ ~12.  So we assert exponential decay over the
# PRE-FLOOR window N ∈ {4,6,8,10,12}: a strongly-negative log-err-vs-N slope
# AND a per-stage geometric drop, then confirm the floor (N=16 does not keep
# falling).  Mutating any D-matrix scaling destroys the geometric rate.
#
# ## Mutation-proof (Rule 4) — see the footer for the full procedure.
#
# Distinct from neighbouring corpus files: differential_test.jl asserts
# AGREEMENT with a second code at one resolution; certified_oracle_test.jl
# BOUNDS the value in an Arb ball at one resolution.  This file asserts the
# RATE at which error vanishes under refinement — the code-verification
# dimension neither of those touches.

using Test
using Statistics: mean
using PadeTaylor
using PadeTaylor.PadeStepper: PadeStepperState, pade_step!

# -----------------------------------------------------------------------------
# Shared infrastructure: observed-order fit over a calibrated error window.
# -----------------------------------------------------------------------------

# The asymptotic window: a refinement point counts toward the slope fit only
# if its error is below ERR_HI (not O(1)-saturated) AND above ERR_LO (not in
# the Float64 round-off floor where the measured "error" is cancellation
# noise).  See the "CALIBRATION" header sections.
const ERR_LO = 1.0e-12
const ERR_HI = 1.0e-1

# Single Padé-Taylor step from z=0 of size h at Taylor order N; return the
# absolute error of u(h) against the manufactured closed form.  `pade_step!`
# is the package's one-step primitive (src/PadeStepper.jl); the dense
# trajectory machinery is bypassed so we measure the STEP's own order.
function _step_error(f, exact_u, u0, up0, N::Int, h::Float64)
    st = PadeStepperState{Float64}(0.0, u0, up0)
    pade_step!(st, f, N, h)
    return abs(st.u - exact_u(h))
end

# Fit the observed order: refine h by halving from h0, keep only the points
# whose error lies inside [ERR_LO, ERR_HI], and least-squares-fit log(err)
# against log(h).  Returns (slope, n_points_in_window).
function _observed_order(f, exact_u, u0, up0, N::Int;
                         h0::Float64 = 0.4, npts::Int = 8)
    hs   = [h0 * 0.5^(k - 1) for k in 1:npts]
    errs = [_step_error(f, exact_u, u0, up0, N, h) for h in hs]
    keep = [i for i in 1:npts if ERR_LO ≤ errs[i] ≤ ERR_HI]
    lx = log.(hs[keep])
    ly = log.(errs[keep])
    mx = mean(lx)
    my = mean(ly)
    slope = sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
    return slope, length(keep)
end

# Theoretical single-step order for the diagonal (m,m) Padé, m = N÷2.
_theory_order(N::Int) = 2 * (N ÷ 2) + 1

# -----------------------------------------------------------------------------
# Manufactured smooth (pole-free, both-parity) anchors.
# -----------------------------------------------------------------------------

f_exp(z, u, up)  = u                 # u'' = u,  IC (1,1)  ⇒ u = exp z
exact_exp(z)     = exp(z)
f_trig(z, u, up) = -u                # u'' = -u, IC (1,1)  ⇒ u = sin z + cos z
exact_trig(z)    = sin(z) + cos(z)

@testset "Convergence-order / MMS (krgy.6): observed order tracks theory" begin

    # =========================================================================
    # CV.1 — IVP observed order tracks the diagonal-Padé law 2·(order÷2)+1.
    # Three orders → three DISTINCT theoretical slopes (3, 5, 7); the slope
    # tracking N is the strongest form (order is a knob the slope follows).
    # Band ±0.5 absorbs the mild positive bias from sub-leading h^(2m+2)
    # terms at the large-h end of the window (probe: slopes land 3.07, 5.13,
    # 7.22 — all within +0.3 of theory).
    # =========================================================================
    @testset "CV.1  IVP h-refinement: slope ≈ 2·(N÷2)+1, tracks N" begin
        for (name, f, exact, u0, up0) in (
                ("u=exp z (u''=u)",          f_exp,  exact_exp,  1.0, 1.0),
                ("u=sin z+cos z (u''=-u)",   f_trig, exact_trig, 1.0, 1.0))
            @testset "$name" begin
                slopes = Float64[]
                for N in (2, 4, 6)
                    slope, npts = _observed_order(f, exact, u0, up0, N)
                    theory = _theory_order(N)
                    # ≥ 3 points must survive the [ERR_LO, ERR_HI] window, else
                    # the slope is not a meaningful fit (calibration guard).
                    @test npts ≥ 3
                    # Observed order matches the diagonal-Padé theory.
                    @test theory - 0.5 ≤ slope ≤ theory + 0.5
                    push!(slopes, slope)
                end
                # The slope STRICTLY INCREASES with order: order 2 < 4 < 6
                # ⇒ slope ≈ 3 < 5 < 7.  This is "the order is a knob and the
                # observed slope follows it" — a single global mis-scaling
                # cannot satisfy three increasing, separated targets.
                @test slopes[1] < slopes[2] < slopes[3]
                @test slopes[2] - slopes[1] > 1.0    # gap ≈ 2, comfortably > 1
                @test slopes[3] - slopes[2] > 1.0
            end
        end
    end

    # =========================================================================
    # CV.2 — BVP Chebyshev spectral (exponential) convergence + floor.
    # Two analytic solutions; assert geometric decay over the pre-floor
    # window N ∈ {4,6,8,10,12}, then confirm the spectral floor at N=16.
    # =========================================================================
    @testset "CV.2  BVP spectral convergence: exponential decay then floor" begin
        f_lin(z, u)   = u                 # u'' = u   on [-1,1] ⇒ cosh(t)/cosh(1)
        ∂f_lin(z, u)  = one(u)
        f_sin(z, u)   = -u                # u'' = -u  on [0,π/2] ⇒ sin z
        ∂f_sin(z, u)  = -one(u)

        # Per-case: (label, f, ∂f, z_a, z_b, u_a, u_b, exact-in-z).
        cases = (
            ("u''=u ⇒ cosh(t)/cosh(1) on [-1,1]",
             f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0,
             z -> cosh(z) / cosh(1.0)),
            ("u''=-u ⇒ sin z on [0,π/2]",
             f_sin, ∂f_sin, 0.0, π/2, 0.0, 1.0,
             z -> sin(z)),
        )

        # Max node error of a converged BVP solve at resolution N.
        function node_error(f, ∂f, za, zb, ua, ub, exact, N)
            sol = bvp_solve(f, ∂f, za, zb, ua, ub; N = N)
            return maximum(abs(sol.u_nodes[j] - exact(real(sol.nodes_z[j])))
                           for j in 1:N+1)
        end

        for (label, f, ∂f, za, zb, ua, ub, exact) in cases
            @testset "$label" begin
                Ns   = (4, 6, 8, 10, 12)          # pre-floor convergent window
                errs = [node_error(f, ∂f, za, zb, ua, ub, exact, N) for N in Ns]

                # (a) Strongly-negative log-err-vs-N slope ⇒ exponential
                #     convergence (geometric in N).  A merely algebraic
                #     rate p (err ~ N^-p) would give a slope of log(err) vs
                #     log(N), NOT vs N; here log(err) vs N is linear and
                #     steeply negative (probe: ≈ -3.0 per unit N).  Assert
                #     clearly below an algebraic strawman: ≤ -1.5 per unit N.
                lx = Float64.(collect(Ns))
                ly = log.(errs)
                mx = mean(lx); my = mean(ly)
                logslope = sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
                @test logslope ≤ -1.5

                # (b) Per-stage geometric drop: each +2 in N shrinks the
                #     error by a large factor in the convergent regime.
                #     Probe ratios are ~1e-2 to ~1e-3 per ΔN=2; assert < 0.1.
                for k in 1:length(Ns)-1
                    @test errs[k+1] < 0.1 * errs[k]
                end

                # (c) Cumulative: N=4 → N=12 collapses by ≥ 8 decades
                #     (probe: 1.7e-4 → 3e-15 ≈ 11 decades).  This is the
                #     headline "faster than any algebraic rate" assertion.
                @test errs[end] < 1.0e-8 * errs[1]

                # (d) FLOOR: past N≈12 the error saturates at the spectral
                #     conditioning limit (≈ cond(D₂)·eps ≈ N²·eps ~ 1e-14);
                #     it does NOT keep dropping geometrically.  Confirm the
                #     N=16 error is no longer below the N=12 error by a
                #     geometric factor — it has hit the floor and is O(1e-14).
                err16 = node_error(f, ∂f, za, zb, ua, ub, exact, 16)
                @test err16 < 1.0e-12               # we are AT the floor
                @test err16 > 0.1 * errs[end]       # NOT a further geom. drop
            end
        end
    end
end

# -----------------------------------------------------------------------------
# Mutation-proof (Rule 4) — performed manually this session, src/ restored
# before commit (`git diff src/` clean).  Each mutation degrades the ORDER /
# spectral RATE and drops the observed slope below the asserted band (RED).
# A single-resolution value test would PASS every one of these.
#
#   CV.1 — IVP order degradation.  In `src/Coefficients.jl`,
#     `taylor_coefficients_2nd`, drop one Taylor term by truncating the
#     bootstrap loop early: change `for j in 2:order` to `for j in 2:order-1`
#     (the top coefficient `u_order` is left at its seeded 0).  The local
#     jet is then one degree short, the diagonal Padé loses a matched degree,
#     and the single-step order drops by ~2: order-6 observed slope falls
#     from ≈7.2 to ≈5, BELOW the `[6.5, 7.5]` band for N=6.  RESULT: CV.1
#     `theory-0.5 ≤ slope` RED for N=6 on BOTH anchors, and the strict-
#     increase `slopes[2] < slopes[3]` RED (orders 4 and 6 now share m_eff).
#     Restored.
#
#   CV.1 — alternative (step mis-scale).  In `src/PadeStepper.jl`,
#     `_rescale_by_powers`, change `h_pow = h_pow * h` to
#     `h_pow = h_pow * h * 1.0001` (a 0.01%-per-degree scale drift).  The
#     rescaled jet no longer represents `u(z + h·t)`, the Padé interpolates
#     a perturbed function, and the leading error term picks up an O(h)
#     contamination that flattens the slope by ~1 toward the lower orders.
#     RESULT: CV.1 band assertions RED.  Restored.
#
#   CV.2 — BVP spectral break.  In `src/BVP.jl`, `bvp_solve`, perturb the
#     second-derivative matrix scaling: change `D2 = D1 * D1` to
#     `D2 = D1 * D1 .* 1.01` (a 1% scale error on D₂).  The discrete
#     operator no longer matches `u''`, so the collocation solution is
#     consistently wrong by O(0.01); the node error FLOORS at ~1e-2 for
#     every N and the geometric decay vanishes.  RESULT: CV.2 (a) `logslope
#     ≤ -1.5` RED (slope ≈ 0), (b) per-stage `errs[k+1] < 0.1·errs[k]` RED,
#     and (c) the 8-decade collapse RED.  Restored.
#
#   CV.2 — alternative (affine factor).  In `src/BVP.jl`, change
#     `scale = (z_b_CT - z_a_CT)^2 / 4` to `... / 5`.  The chain-rule
#     factor on the affine map is wrong, the RHS is mis-weighted, and the
#     Newton iterate converges to the wrong function — node error floors
#     well above eps and the spectral decay breaks.  RESULT: CV.2 (a)+(c)
#     RED on both cases.  Restored.
#
# `git diff src/` is clean after all restorations.  Together these prove the
# file asserts the RATE of convergence, not just a single-resolution value:
# every order-degrading or rate-breaking src/ perturbation turns it RED.
