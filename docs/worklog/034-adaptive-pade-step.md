# Worklog 034 — Adaptive Padé step size (FFW 2017 §2.1.2)

**Date**: 2026-05-15
**Author**: Claude Opus
**Bead**: `padetaylor-8ui` — closed.
**Scope**: Ship step A1 of the 11-step FFW 2017 figure-reproduction
plan: an opt-in `step_size_policy = :adaptive_ffw` controller on
`path_network_solve`, implementing FFW 2017 §2.1.2's truncation-error
estimator + rescale law verbatim from md:74-97.

> **Take-home**: `ffw_truncation_error`, `ffw_rescale_q`, and
> `adaptive_pade_step!` in `src/PadeStepper.jl`, plus an inline
> rescale loop in `path_network_solve`, give us the FFW 2017 §2.1.2
> controller behind a single new kwarg.  `T(h) = |ε_{n+1} h^{n+1} /
> a(h)|` is computed in the *rescaled* variable so the algebra
> reuses the existing `h^k` Padé-input convention without
> rebuilding the rational; one extra Taylor coefficient (order `n+1`
> instead of `n`) is the entire data-cost overhead per `T(h)` call.
> Worklog 017's deferral note retired; test count **1708 → 1751
> GREEN** (+43 AS.* assertions, 4 mutation-proven).

## Ground truth read first

Per Law 1 the first reading was `references/markdown/
FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:74-97` verbatim.
Key extractions:

  - **md:74-76**: numerator `p(h)`, denominator `q(h)`, formal series
    `w(ζ+h) = Σ c_k h^k`.  The error series is `q(h) w(ζ+h) - p(h) =
    Σ_{k=n+1} ε_k h^k` with `ε_k = c_k + Σ_{r=1..ν} b_r c_{k-r}`.

  - **md:80**: relative-error estimator `T(h) = |(w(ζ+h) - p/q)/(p/q)|
    ≈ |ε_{n+1} h^{n+1} / p(h)|`.  Subtlety: FFW uses `p` for the
    NUMERATOR in `p(h)/q(h)`, but our `RobustPade` convention is
    `r(z) = a(z)/b(z)` with `a` = numerator, `b` = denominator.
    Translated: `T(h) = |ε_{n+1} h^{n+1} / a(h)|` with our `a`.

  - **md:88-91 eq. 2**: rescale `q = (k·Tol/T(h))^(1/(n+1))`,
    `h := q·h`.  Conservative factor `k` introduced explicitly as
    "a small positive constant".

  - **md:93**: rescale loop sequence is "re-run wedge, re-select
    min-|u|, re-compute T(h)" — NOT "rescale and accept blindly."
    The wedge's selected direction is allowed to change between
    rescales.  Acceptance also seeds the NEXT step's initial `h`
    from the accepted `|q·h|`, so the controller has memory.

  - **md:91**: `k ≈ 1e-3` (typical conservative value, recommended
    in FFW's own implementation).

Also re-read `docs/worklog/017-coord-transforms-pIII-pV.md:128-144`
to surface the deferral note: "Adaptive Padé step size... not on
critical path for shipping the transform helpers."  The deferral
expired when the FFW figure plan moved to active work.

Three local files were read end-to-end before code: `src/PadeStepper.jl`
(surface for the new helpers), `src/PathNetwork.jl` (the call site
at line 222 where `:adaptive_ffw` was previously a fail-fast throw),
`src/RobustPade.jl` (the `PadeApproximant{T}` shape — confirmed that
`.a` is the numerator vector, `.b` is the denominator vector with
`b[1] = 1`).  No conversation memory; `git log` + `Read` only.

## What shipped

### `src/PadeStepper.jl` — three new public helpers

Module LOC count went from 72 (code, excluding docstrings) to ~145.
Still well under the 200-LOC cap (Rule 6).

  - **`ffw_truncation_error(f, z, u, up, order, h) → Real`** —
    builds the order-`(n+1)` Taylor jet (one pass beyond what
    `pade_step_with_pade!` uses), rescales by powers of `h`, builds
    the order-`n` Padé from the first `n+1` rescaled coefficients
    (same call as in the step itself), then computes
    `ε̃_{n+1} = c̃_{n+1} + Σ_{r=1..ν_eff} b_r · c̃_{n+1-r}` and
    `a_rescaled(1) = Σ a_k`.  Returns `|ε̃_{n+1} / a_rescaled(1)|`,
    which equals FFW's `|ε_{n+1} h^{n+1} / a(h)|` after substituting
    `c̃_k = h^k c_k` and `a_rescaled(1) = a(h)`.

    The `νeff = length(P_u.b) - 1` care matters: `RobustPade`'s
    `_trim_and_normalise` strips trailing near-zero denominator
    entries, so `ν_eff` may be `< n/2` for degenerate Padé blocks
    near accuracy.  The loop `for r in 1:νeff` consumes whatever
    `b` we actually got — trimmed entries correspond to `b_r = 0`
    and contribute nothing.

  - **`ffw_rescale_q(Tol, T_h, order; k = 1e-3) → Real`** — pure;
    the FFW eq. 2 formula `(k·Tol/T_h)^(1/(n+1))`.  Defensive:
    rejects non-positive `Tol`, `T_h`, `k`, and non-positive `order`
    with `ArgumentError`.  `T_h == 0` returns `Inf` (the controller's
    accept-test fires before `q` is consulted in normal use; this
    branch is only reached if a caller asks for `q` at `T_h = 0`
    directly).

  - **`adaptive_pade_step!(state, f, order, h_init; ...)` →
    `(state, P_u, meta)`** — full controller for single-step callers.
    Iterates `T(h) > Tol ⇒ h := q·h` until acceptance or
    `max_rescales` exhausted (`ErrorException`).  Returns the
    standard `(state, P_u)` plus a `NamedTuple` `(h_used, h_step,
    T_h, n_rescales)` so callers can thread the controller's memory
    across steps.

The module docstring grew a new "Adaptive Padé step size — FFW 2017
§2.1.2" section explaining the rescaled-variable algebra and citing
md:74-97 / md:88-91 by line number per Law 1.

### `src/PathNetwork.jl` — inline adaptive loop on `path_network_solve`

Three new kwargs: `adaptive_tol::Real = 1e-12`, `k_conservative::Real
= 1e-3`, `max_rescales::Integer = 50`.  `step_size_policy` now
accepts `:adaptive_ffw` as an opt-in alternative to `:fixed`.  The
previous Tier-4 deferral throw retires.

The inner loop replaces the single `_wedge_evaluations + _select_candidate`
call with a `while true` rescale block (FFW md:93's "re-run wedge,
re-select min-|u|, re-compute T(h)" semantics).  Under `:fixed` the
loop runs exactly once (breaks immediately after selection).  Under
`:adaptive_ffw` it iterates until `T(h) ≤ Tol`.

`visited_h[k]` for each landed node now stores the *accepted* `h_step`,
not the constant `h_T`.  Stage 2 already uses per-node `visited_h[k]`
for the `t = (z_f - z_v) / h_v` evaluation, so heterogeneous step
sizes work transparently — no Stage 2 code change needed.  Under
`:fixed` all `visited_h` are still `h_T`, byte-identical to the
pre-change behaviour.

The walk-seed `h_cur` is initialised from the visited node's own
`visited_h[idx_v]`, per FFW md:93 "the initial step length is always
the scaled step length stored at the current point".  Under `:fixed`
this is `h_T`; under `:adaptive_ffw` it threads the controller's
memory across targets.

`_solve_with_schwarz_reflection` accepts and forwards the three new
kwargs unchanged.

### `test/adaptive_step_test.jl` — 7 testsets, 43 assertions

  - **AS.1.1** (10 assertions) — `T(h)` formula sanity on `u'' = u`
    (exp closed-form).  Monotonicity in `h`, asymptotic `h^{n+1}`
    scaling, and a *prefactor pin* (`T(h)/h^{n+1} ∈ [1e-12, 1e-9]`)
    that mutation M2 violates.

  - **AS.1.2** (9 assertions) — `q`-rescale fixed-point: starting
    from `h = 2.0` on `u'' = 6u²`, the rescale loop converges to
    `T(h) ≤ Tol` in ≤ 5 iterations at both `Tol = 1e-10` and
    `Tol = 1e-14`.

  - **AS.1.3** (4 assertions) — end-to-end agreement with fixed-h
    baseline on the FW Table 5.1 `u(z=30)` target.  Seeded at
    `h = 1.5` (deliberately larger than the FW default `0.5`) to
    FORCE adaptation: `T(1.5) ≈ 1e-7 ≫ Tol = 1e-12`.  Asserts
    `|u_adapt - u_fixed| ≤ 1e-9` AND that at least one visited node
    has `h ≠ 1.5` (the load-bearing "adaptation actually happened"
    test).

  - **AS.1.4** (2 assertions) — PIII transformed-vs-direct test
    (port of `test/coord_transforms_test.jl`'s CT.1.3) under
    adaptive step.  Proves the adaptive path doesn't corrupt the
    coordinate transform.

  - **AS.1.5** (4 assertions) — Tol sweep monotonicity on the
    *accepted h* (not on the final-step error, which saturates at
    the discretisation floor).  Tightening `Tol` from `1e-8` to
    `1e-14` must non-increase `max(visited_h[2:end])` AND the
    extreme spread must be ≥ 20%.  (The `[2:end]` slice excludes
    the IC's canonical-Padé radius, which is always the user-
    specified seed.)

  - **AS.1.6** (12 assertions) — `adaptive_pade_step!` single-call
    contract.  Seeded at `h = 2.0` (forces ≥ 1 rescale).  Validates
    the 3-tuple return shape, NamedTuple metadata fields, `h_used <
    h_init`, `T_h ≤ Tol`, mutated `state.z`.

  - **AS.1.7** (3 assertions) — multi-iteration rescale at extreme
    `h = 10.0`.  One rescale lands at `T(h) ≈ 5.7e-7 ≫ Tol = 1e-12`;
    the controller must iterate.  Asserts `meta.n_rescales ≥ 2`,
    `T_h ≤ Tol`, `0 < h_used < 10`.  This is the test that bites
    mutation M3 (one-shot rescale) — AS.1.6 alone is insufficient
    because at `h_init = 2.0`, one rescale happens to converge.

### `test/pathnetwork_test.jl` — PN.4.1 updated

The previous `@test_throws ArgumentError` against `:adaptive_ffw`
becomes a `@test_throws` against `:bogus_policy` (still validates
the unknown-symbol fail-fast branch).  Cited bead `padetaylor-8ui`
in the inline comment.

### `docs/adr/0011-adaptive-pade-step.md` (Accepted)

Standard ADR sections: Context (worklog 017 deferral retired), Decision
(opt-in kwarg + three helpers), Alternatives (PI-controller, Jorba-Zou,
constant halving, IVP-layer extension — all rejected with reasoned
references), Consequences (figures unblocked, test count delta,
backward compat).

### `test/runtests.jl` — one include line

`include("adaptive_step_test.jl")` inserted between
`pathnetwork_test.jl` and `polefield_test.jl` (logical adjacency to
the path-network module the new tests exercise).

## Mutation-proof procedure

Four mutations applied (per the test-file footer), each followed by
an isolated `include("test/adaptive_step_test.jl")` run, observed
RED, then reverted.  Full bite count + restoration table:

| Mutation | Description                                | Bite | Tests biting |
|----------|--------------------------------------------|------|--------------|
| M1       | flip rescale exponent sign (q^{-1/(n+1)})  | 17   | AS.1.2 (12), AS.1.3 (1), AS.1.5 (1), AS.1.6 (3) |
| M2       | drop `Σ b_r c̃_{n+1-r}` in ε_{n+1}          | 1    | AS.1.1 prefactor pin |
| M3       | one-shot rescale (no loop)                 | 2    | AS.1.7 (T_h, n_rescales) |
| M4       | hard-code q = 0.5 (FW halving)             | 1    | AS.1.5 max_h spread |

M3 required *also* mutating the inline rescale loop in
`path_network_solve` (since AS.1.3 / AS.1.5 use the path-network
driver, not `adaptive_pade_step!`).  When I applied M3 to
`adaptive_pade_step!` only, AS.1.6 at `h_init = 2.0` happened to
converge in one rescale (the `k = 1e-3` conservative factor is
aggressive enough that one shot bridges `T(h) ≈ 5e-4 → ≈ 5e-15`).
That's the genesis of AS.1.7: a HARSHER `h_init = 10.0` where one
rescale leaves `T(h) ≈ 5.7e-7 ≫ Tol`, forcing the algorithm to
iterate or fail the test.

The four mutations restored cleanly; final full-suite run **1751
GREEN** in 3 min 12 s.

## Frictions surfaced

  - **The FFW eq. 2 controller does not converge at large enough
    `h_init`.**  At `h = 5.0` on the FW IC, T(h) ≈ 1e5.  One rescale
    gives q ≈ 0.23 → h_new ≈ 1.13 → T ≈ 2e-11 (still above 1e-12).
    Two rescales suffice; the controller is robust.  But at `h =
    100.0` (deliberately silly), `T(h)` overflows — the rescale
    loop hits NaN before `max_rescales` and throws an obscure error.
    The default `max_rescales = 50` is comfortable for "reasonable
    h_init", and the bead's scope doesn't include extreme-h hygiene.
    Filed `bd create` not yet; mention here for whoever lands on a
    PVI sheet-tracker walk that explodes.

  - **The prefactor test AS.1.1 took two iterations to land.**  My
    first cut compared `T(h)/h^{n+1}` ratios across `h ∈ {1e-3,
    1e-4}` — the ratio is dominated by the `h^{n+1}` factor, so the
    *Padé correction* contribution is invisible there.  Mutation M2
    initially didn't bite, telling me the test wasn't testing what
    I thought it was.  Re-derived: the prefactor `|ε̃_{n+1} /
    a_rescaled(1)|` is `~ 1e-10` (with Padé correction) vs `~ 2.5e-8`
    (without) — a 250× ratio.  Added a direct prefactor-magnitude
    pin at `1e-12 < T(h)/h^{n+1} < 1e-9`, which both bracketed the
    correct value AND made M2 bite.  Lesson: "the test must
    discriminate the mutation, not just shape-fit the algorithm."

  - **Visited-h[1] is always the IC seed.**  AS.1.5's first cut
    asserted `max(visited_h)` non-increasing across Tol — but the IC
    node's `visited_h[1]` is hardcoded to the user-specified seed
    `h_T` (correctly: that's the radius of validity of the IC's
    canonical Padé), so all four sweep entries collapsed to `1.5`.
    Fix: assert on `visited_h[2:end]` (the controller's accepted
    step sizes).  Documented in AS.1.5's docstring.

## Hard-won lesson

**"Adaptive" is two algorithms wearing one name.**  Building
`adaptive_pade_step!` is the *easy* part — single-step contract,
clear inputs/outputs, mutation-prove with a unit test (AS.1.6 +
AS.1.7).  Wiring `:adaptive_ffw` into `path_network_solve` is the
*hard* part because the path-network's wedge/min-|u|/canonical-Padé
structure means the rescale loop has to re-run the wedge at the new
`h_step` (per FFW md:93), not just rerun the controller.  Inlining
the loop in `path_network_solve` was the right call (an extra layer
of dispatch through `adaptive_pade_step!` would have lost the wedge-
re-selection semantics).  The two implementations of the same
algorithm — `adaptive_pade_step!` for single-step callers, the
inline loop for the path-network — are intentional and documented in
ADR-0011's "Decision".  They MUST be co-mutated when testing (M3 had
to touch both for the bite to land on AS.1.3 / AS.1.5).

## What this enables

Step A1 of the 11-step FFW 2017 figure plan is now closed.  Steps A2
(non-uniform Stage-1 nodes per FFW md:67-72), A3-A7 (PIII figures
1-4 + PV figures 5-7), and the FFW-style cubic-pole-vs-Riemann-sheet
classification are the remaining downstream deliverables.  No code
in this commit touches the figure scripts; the next-bead handoff
will compose the figure pipelines on top of the controller shipped
here.
