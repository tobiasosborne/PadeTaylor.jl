# Worklog 080 — FW Fig 4.1 step-(i) BVP at BigFloat-256 (bead `padetaylor-7zw`)

**Date**: 2026-08-23
**Scope**: Rerun the tritronquée `[-20i, +20i]` Chebyshev-Newton BVP of
worklog 016 at 256-bit precision, measure what the path achieves, pin it.
Test-only delivery (`test/fw_fig_41_test.jl`, +111 lines); no `src/`
changes.  File GREEN at 12 + **8** + 5 assertions.

> Take-home: at N=240 the BF-256 run reproduces the Float64 u(0)/u'(0)
> errors to three digits (rel 1.85e-12 / 1.75e-10) — the error is
> spectral truncation, not precision.  What 256 bits buys at this N is a
> residual floor of 1e-68 and Schwarz-imaginary floor of 1e-76 (Float64:
> 1e-8 / 1e-16), and THAT is what the new testset pins as its precision
> witness.  FW's 16-digit quoting limit is reached at N=320 (75 s).

## Ground truth

`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`:

  - :226 — eq. 4.1: `u(0) ≈ -0.1875543083404949, u'(0) ≈ 0.3049055602612289`.
  - :228 — "using an extended precision BVP solver in Maple … these two
    values above are obtained to an accuracy of better than 10-20".
  - :214-222 — Fig 4.1 caption, step (i): BVP over `[-20i, 20i]` with the
    leading term of eq. 1.2 (:54-55) as end data.

Note on bead acceptance (b): it mentions "the IC value 1.85185330…".
That number is the NIST tronquée `u'(0)` of FW Fig 4.3 (md:245; worklog
024), not part of the Fig 4.1 step-(i) recipe, whose only input data are
the leading-term BCs `-√(-z/6)` at `±20i`.  Those BCs are built directly
in `Complex{BigFloat}` (no literal needed); the two `big"…"` literals are
the eq. 4.1 targets.

## Procedure

Scratch probe (`setprecision(BigFloat, 256)`, `bvp_solve` with default
`tol = eps(BigFloat)^(3/4)`, `initial_guess = leading`, `max_iter = 20`),
relative errors against the eq. 4.1 values parsed as `big"…"`:

| N   | iters | wall   | rel-err u(0) | rel-err u'(0) | ‖R‖_∞   | max imag |
|-----|-------|--------|--------------|---------------|---------|----------|
| 120 | 8     | 3.6 s  | 3.95e-6      | 7.74e-5       | 2.3e-70 | 5e-77    |
| 160 | 8     | 9.6 s  | 2.85e-8      | 1.12e-6       | 1.5e-69 | 3e-76    |
| 200 | 7     | 17.2 s | 2.26e-10     | 1.45e-8       | 1.3e-69 | 2e-76    |
| **240** | 7 | **29.9 s** | **1.85e-12** | **1.75e-10** | 1.1e-68 | 3e-76 |
| 280 | 7     | 52.1 s | 1.54e-14     | 2.02e-12      | 1.2e-68 | 1e-76    |
| 320 | 7     | 74.3 s | 1.62e-16     | 2.25e-14      | 2.1e-68 | 2e-76    |

Julia 1.12.3, single process (Rule 7).  Wall ∝ N³ (GenericLinearAlgebra
LU on a 239×239 `Complex{BigFloat}` matrix, 7 Newton solves).

Compare worklog 016 at Float64, N=240: abs 3.5e-13 ⇔ rel 1.87e-12 for
u(0); abs 5.3e-11 ⇔ rel 1.7e-10 for u'(0).  Same numbers.  Precision is
not the bottleneck at N=240 (worklog 016 friction 3, now measured).

N=320 gives rel 1.6e-16 on u(0): the computed value
`-0.18755430834049487…` rounds to FW's 16 quoted digits, i.e. the
reference itself is exhausted.  Going beyond needs a longer reference
(FW's "10⁻²⁰" is a claim, not a printed number) and is out of scope.

## Decision: keep N=240 (bead (a)), pin what it achieves (bead (c))

  - Wall 30 s (≈60 s in-suite including BigFloat-path compilation) is
    under the 90 s ceiling; N=320 would cost 75 s + compile per run for a
    pin no test can use (reference limit).  Trade-off recorded above.
  - **FB.1.2** — `rel-err u(0) ≤ 1e-11` (measured 1.85e-12, 5.4× margin;
    tighter than the Float64 FF.1.2 pin 1e-10 abs ⇔ 5.3e-10 rel).
    `rel-err u'(0) ≤ 3e-10` (measured 1.75e-10, 1.7× margin).  The
    bead's 5× margin would be 1e-9 rel ⇔ 3.0e-10 abs, looser than the
    Float64 1e-10 abs pin; 3e-10 rel ⇔ 9.1e-11 abs keeps the BF pin no
    looser than Float64.  The thin margin is safe: MPFR and
    GenericLinearAlgebra are bit-reproducible (no LAPACK variance).
  - **FB.1.3** — `‖R‖_∞ ≤ 1e-60`, `|imag u(0)|, |imag u'(0)| ≤ 1e-60`
    (measured 1.1e-68, 3e-76).  This is the precision witness: the
    Float64 floors are ~1e-8 / ~1e-15.  Worklog 016 lesson 25 ("don't
    gate on residual") is about convergence; here the residual is gated
    as evidence the solve ran at 256 bits, eight orders below its own
    `N²·eps` estimate's Float64 analogue.
  - **FB.1.1** — converged ≤ 10 iters (7); **FB.1.4** — Dirichlet BC at
    `+20i` reproduced to 1e-70 (structural; exact at any precision).

## Mutation evidence (Rule 4)

  - M1 — `u_at_0_FW` → `big"-0.1875543083434949"` (+3e-12, 12th decimal,
    just inside the 1.9e-12-abs pin): FB.1.2 RED, `rel_u = 1.78e-11`.
    `up_at_0_FW` → `big"0.3049055603612289"` (+1e-10): RED,
    `rel_up = 5.03e-10`.
  - M2 — endpoints as `-20.0im, 20.0im` (Float64 → whole solve at 53
    bits): FB.1.3 RED on all three (`residual 3.76e-8`, imag parts
    1.2e-15 / 1.6e-15).  FB.1.2 and FB.1.4 stay GREEN, confirming they
    are not precision witnesses; FB.1.3 is.
  - Both mutants run from scratchpad copies; the in-tree file was never
    mutated.  Restored state GREEN 8/8.

## Files

  - `test/fw_fig_41_test.jl` — new testset "FW 2011 Fig 4.1 step-(i) BVP
    at BigFloat-256" + mutation footer.
  - `docs/worklog/080-fw-fig-41-bf256-pin.md` — this file.
  - `CHANGELOG.md` — [Unreleased] bullet.
