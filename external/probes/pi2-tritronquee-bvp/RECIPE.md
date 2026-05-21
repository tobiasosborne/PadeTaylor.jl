# RECIPE — isolating the P_I^(2) tritronquée via the vector BVP

Probe `padetaylor-0ln.29`, continuation. Steps 9–12.

## TL;DR

The prior probe's premise was **inverted by a sign error**. There was never
a "spurious +branch attractor". Once the asymptotic seed carries the
correct sign, the *plain* 2+2 (u,u') endpoint BVP — exactly what
`VectorBVP.vector_bvp_solve` already does — converges to the genuine
P_I^(2) tritronquée in **2–3 Newton iterations**, ODE residual ~1e-11,
robust across every N (40–256) and every segment tested. **No collar, no
deflation, no deviation recast, no continuation is needed.**

## The sign error (root cause — CLAUDE.md Rule 2: all bugs are deep)

P_I^(2) at t=0 in companion form has the leading dominant balance, from
`u'''' + 10u'^2 + 20uu'' + 40(u^3 + 6x) = 0`,

    40 u^3 + 240 x ≈ 0   ⇒   u^3 = -6x .

On the **negative real axis** (x < 0): `u^3 = -6x = 6|x| > 0`, so the
**real** root is

    u = +∛6 · |x|^{1/3}   (POSITIVE on x<0).

Direct check (step 11): the leading term `u = -∛6·|x|^{1/3}` (negative)
has ODE residual ≈ **-4800** at x = -10 — it is **not a P_I^(2) solution
at all**. The positive leading term has residual ≈ -0.5 (the genuine
algebraic-asymptotic decay rate).

`PainleveHierarchy.pI2_tritronquee_ic` returns `y[1] = -∛6·r^{1/3}` with
`r = |x0|` — i.e. `-∛6·|x|^{1/3}`, **negative**. KKG's asymptotic
(`tritronquee_coeff.tex:3105`, eq. `asym`) is `u ~ -∛6·x^{1/3}`. For
`x < 0` the real cube root `x^{1/3}` is **negative** (`= -|x|^{1/3}`), so
`-∛6·x^{1/3} = +∛6·|x|^{1/3}` — **positive**. The `src` builder substituted
`|x|^{1/3}` for `x^{1/3}`, flipping the sign on the negative axis.

Consequently every prior probe step pinned `u` to a *negative* endpoint
value that no real P_I^(2) solution can attain, so Newton produced
humped/oscillatory garbage that matched the (impossible) negative endpoint
data but bulged positive in the interior. That artifact was misread as a
"spurious +branch" with its own Newton basin. It has no such thing.

## The working recipe (Approach: none of A/B/C — the BVP was already correct)

Formulation: the existing `VectorBVP.vector_bvp_solve` on the companion
4-vector `y = (u,u',u'',u''')`, **unchanged**.

- **Segment:** any `[z_a, z_b]` on the negative real axis with `z_a` the
  large-|x| end; tested `[-60,-5] … [-15,-2]`. Larger `|z_a|` → smaller
  truncation error (more asymptotic terms negligible).
- **BC structure:** 2+2 (u,u') split — pin `(u,u')` at *both* endpoints:
  `B_a = diag(1,1,0,0)` rows 1–2, `B_b` rows 3–4 picking `(y_1,y_2)`;
  `g = [u(z_a), u'(z_a), u(z_b), u'(z_b)]` from the correct-sign seed.
  (This is the prior probe's split; it was always right.)
- **Correct-sign seed** (`x<0`, `r = -x`, `n_terms = 2`, t = 0):

      u  = +∛6·r^{1/3} + r^{-2}/36
      u' = -∛6·(1/3)r^{-2/3} - 2·r^{-3}/36
      u''= +∛6·(1/3)(2/3)r^{-5/3} + 6·r^{-4}/36
      u'''=-∛6·(1/3)(2/3)(5/3)r^{-8/3} - 24·r^{-5}/36

  (KKG t=0 series `u = Y + Y^{-6}`, `Y = +∛6·r^{1/3}` so `Y^3 = 6r`
  matches `u^3 = 6r`; `Y^{-6} = 6^{-2}r^{-2} = r^{-2}/36` is KKG's
  `b_1 = 1/36`, `tritronquee_coeff.tex:2688`. `d/dx = -d/dr`.)
- **Interior initial guess:** the correct-sign seed sampled at every node
  (`initial_guess = z -> tri_seed(z)`), OR a linear interpolation between
  the two endpoint seeds (KKG's documented alternative,
  `tritronquee_coeff.tex:3176`). A **zero** interior guess does NOT
  converge — the basin needs a tritronquée-leaning interior guess.
- **N:** anything ≥ 40 works identically (spectral; converged).
- **Jacobian:** analytic `painleve_hierarchy_jacobian` or `Taylor1`
  autodiff — both work, identical result.
- **tol / maxiter:** `tol = 1e-9`, `maxiter = 40` (converges in 2–3).

## Converged solution (2+2, [-20,-2], N=128, correct-sign seed)

    x      u_num        u_tri (KKG n=2)   |Δ|
   -20    4.932494     4.932494          0.0
   -16    4.578965     4.578965          5.4e-8
   -12    4.160360     4.160361          1.9e-7
   -8     3.634674     3.634675          1.3e-6
   -4     2.886202     2.886235          3.4e-5
   -2     2.296373     2.296373          0.0

- `u` is **positive and monotone** throughout (decreasing as x→0⁻).
- ODE residual `‖R‖_∞ ≈ 1.5e-11`.
- Match to the KKG truncated series: `max|u_num - u_tri|`:
  - `[-20,-2]`: 2.1e-4   `[-40,-3]`: 1.8e-5   `[-60,-5]`: 1.4e-6.
  The error shrinks with `|z_a|` because it is the **n_terms=2 truncation
  error of the reference series**, not solver error (the ODE residual is
  1e-11 regardless). At KKG's accuracy target (~1e-6) the recipe meets it
  on `[-60,-5]` and beats it everywhere with more series terms.
- Robust: identical to ~1e-7 across N ∈ {40,80,160,256} and every segment.
- BigFloat (113-bit): converges, ODE residual ~1e-29.

## Approaches A / B / C — what was tried, why they were not needed

All three were implemented and tested (steps 9, 10) BEFORE the sign error
was found, so they were all fighting a phantom:

- **A (asymptotic collar):** an over-determined weighted least-squares
  collocation pinning the (sign-wrong, negative) asymptotic on a collar of
  up to 16 nodes per end. Converged but to the humped artifact (`u` in
  `[2.3, 4.9]`, positive interior) — because the collar pinned impossible
  negative endpoint data while the collocation block (which the artifact
  *does* satisfy) dominated the least-squares. Not a real fix; moot once
  the sign is corrected.
- **B (deviation variable):** `u = u_a + v`, `u_a = -∛6·(-x)^{1/3}` (the
  sign-wrong leading order). Newton from `v=0` converged to `max|v| ≈ 13.4`
  — i.e. `v ≈ 2u_a`, the correction needed to climb from the wrong sign to
  the right one. Confirms the recast itself is sound but its `u_a` had the
  wrong sign.
- **C (deflation):** not reached — A and B already exposed that the
  "wrong branch" was the *only* genuine solution near the (impossible)
  negative BCs, i.e. there was no second genuine solution to deflate.

The KKG global initial iterate `u = -6^{1/3} x/(1+x^2)^{1/3}`
(`tritronquee_coeff.tex:3175`) was also tried (step 9). With the
sign-wrong negative-axis convention it does not help; note that KKG's own
formula, evaluated for real `x<0` with the real cube root, gives
`u → +∛6|x|^{1/3}` — i.e. KKG's iterate is itself **positive** on `x<0`,
independent confirmation of the correct sign.

## PRODUCTIONISATION RECOMMENDATION

`VectorBVP.vector_bvp_solve` needs **no change** — it solved the problem
correctly the entire time. The fix is a one-line **sign correction in
`PainleveHierarchy.pI2_tritronquee_ic`**, plus the paired doc updates.

1. **Fix `pI2_tritronquee_ic`** (`src/PainleveHierarchy.jl:452-479`).
   The companion seed on the negative real axis must be the **correct-sign
   real-branch** values. Concretely, with `r = |x0|`:

       y_1 = +∛6·r^{1/3}              (currently  -∛6·r^{1/3})
       y_2 = -∛6·(1/3)·r^{-2/3}        (sign of every derivative flips too)
       y_3 = +∛6·(2/9)·r^{-5/3}
       y_4 = -∛6·(10/27)·r^{-8/3}

   and the `n_terms=2` correction `Y^{-6}=r^{-2}/36` with `Y=+∛6·r^{1/3}`:

       Δy_1 = +r^{-2}/36   Δy_2 = -2·r^{-3}/36
       Δy_3 = +6·r^{-4}/36 Δy_4 = -24·r^{-5}/36

   The cleanest implementation is to compute `u(x)` and its three
   `x`-derivatives from the closed form `u = ∛6·(-x)^{1/3} + (-x)^{-2}/36`
   using `d/dx = -d/dr` (as `step11.jl`/`step12.jl::tri_seed` do), rather
   than the hand-transcribed component list — that removes the
   transcription-sign-slip class of bug entirely.

   IMPORTANT — this also means the pillar-C doc
   `docs/v0p2_pillarC_painleve_hierarchy_findings.md:269-281` ("Practical
   IC ... `y_1 = -∛6·|x_0|^{1/3}`") and the `PainleveHierarchy` module
   docstring (lines 76-106) are sign-wrong on the negative axis and must
   be corrected in lockstep (CLAUDE.md Law 2). The decisive ground truth
   is the ODE itself: `40u^3 = -240x ⇒ u^3 = -6x > 0` for `x<0`.

   Note the subtlety to record in the docstring: KKG write `u ~ -∛6·x^{1/3}`
   *uniformly*; the sign on the real axis depends entirely on which branch
   of `x^{1/3}` is taken. For the **real, ODE-consistent** solution on
   `x<0` the value is `+∛6·|x|^{1/3}`. On `x>0` it is `-∛6·x^{1/3}` (also
   negative-of-positive — there it genuinely is negative). The IC builder
   must pick the branch by the sign of `x0`.

2. **TDD for the production phase.** The acceptance test is exactly the
   step-11/12 result: 2+2 (u,u') BVP on `[-20,-2]`, correct-sign seed,
   assert `u(x) > 0` and monotone on the negative axis, `‖R‖_∞ < 1e-9`,
   and `max|u - u_KKG|` below the n_terms-2 truncation bound (~2.1e-4 on
   `[-20,-2]`, ~1.4e-6 on `[-60,-5]`). Mutation-proof by flipping the seed
   sign — the test must go RED (it will: KKGerr jumps to ~10, the humped
   artifact).

3. **No new module, no collar/deflation driver, no continuation wrapper.**
   The vector BVP machinery is sufficient. The `vector_bvp_solve(php, …)`
   forwarding method in `PainleveHierarchy` already wires the companion
   RHS + analytic Jacobian; once `pI2_tritronquee_ic` carries the right
   sign, a caller building `g` and the interior guess from it gets the
   tritronquée directly.

4. **Optional hardening (defer as a bead, not v0.2-blocking):** the 3+1
   and 4+0 BC splits fail to converge (ill-conditioned — the separatrix
   growing mode; step 12e). If a future target needs them, the fix is the
   over-determined least-squares collocation prototyped in step 6/10. The
   2+2 split has no such issue and is the production path.
