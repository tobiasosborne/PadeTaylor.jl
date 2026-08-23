# ADR-0027 — Shared-Padé graceful degree-reduction + z^λ cancellation (the C2 axis)

**Status**: Accepted (built and shipped as beads `padetaylor-d3a`/`3p9c`, closed
2026-06-02, per worklog 066/067; live at `src/SharedPade.jl` reduction loop; status flipped 2026-08-23).
**Date**: 2026-06-02
**Beads**: `padetaylor-3p9c` (C2), depends on `padetaylor-d3a` (C1, the `+2`
window). **Supersedes** the `ρ ≥ m_cur` break, the dead `n_near > 1` guard, and
the `Q(0) ≈ 0` throw in `src/SharedPade.jl`.

## Context

ADR-0019's `shared_denominator_pade` stacks one `m×(m+1)` Toeplitz block per
component into an over-determined `dm×(m+1)` system and takes the
smallest-singular-value right vector as the shared denominator `Q`. Worklog 066
settled (unanimous, sympy-verified) that the **correct matching window is `+2`
(top-left `c_{m+1}`, the GGT diagonal `(m,m)`)**. Applying `+2` (C1) exposes a
latent defect in the `d≥2` rank handling that the buggy `+1` window masked:

- For an **entire / far-from-pole** local jet (most stepper steps, and the
  harmonic `[cos,−sin]` test) the stack is numerically **rank-deficient** at the
  stepper's degree (`m=order÷2≈15`; eff-rank ≈ 9–10 ≪ m+1).
- The current loop breaks on `ρ ≥ m_cur` (admitting the full-column-rank case),
  the `n_near > 1` isolation guard is **dead code** for `d≥2`, and after a forced
  reduction the recovered `Q` frequently has **`Q(0) ≈ 0`** — a *spurious pole at
  the expansion centre* — which the current code treats as a hard `throw`
  (Rule-1 over-throw). This breaks `VectorStepper`/`VectorProblems`/
  `calogero_moser`/`vector_pipeline_oracle` (worklog 065).

The scalar `robust_pade` never hits this: its `C̃` is `n×(n+1)` (so `ρ ≤ n`, no
full-rank case), it reduces on `ρ == n` (GGT `md:136`, "reduce `n→ρ`"), and it
**cancels the common `z^λ` factor** on a leading-zero `Q` rather than throwing
(`_trim_and_normalise:461-476`). The shared code's docstring (`:250-251`) wrongly
asserted that cancellation is "not needed here."

## Decision

Mirror the scalar `robust_pade` semantics on the vector path:

1. **Reduction criterion** — reduce until the smallest-σ null space is an
   **isolated 1-D** space, i.e. `ρ = count(σ > τ) == m_cur` (exactly one σ at/below
   `τ = tol·‖c‖`):
   - `ρ < m_cur` (rank-deficient): set `m_cur = ρ`, rebuild. *(GGT reduce-`n`-to-`ρ`.)*
   - `ρ > m_cur` (i.e. `ρ = m_cur+1`, full column rank, **no** σ ≤ τ): set
     `m_cur -= 1`, rebuild. *(No exact shared denominator at this degree.)*
   - `ρ == m_cur`: **accept** (isolated 1-D null space, Mano–Tsuda "unique iff
     rank = m", `tex:1427`).
   - `ρ == 0`: decrement `m_cur` and rebuild. If no degree `m_cur ≥ 1`
     remains, return `Q = [1]` with the full Taylor-polynomial numerators;
     only identically-zero jets fail the earlier `iszero(‖c‖)` guard
     (`docs/worklog/067-sharedpade-c1-c2-minimal-core.md:36-40`).

2. **`z^λ` cancellation** (replaces the `Q(0) ≈ 0` throw) — after the
   QR-reweighting yields `b`, cancel the common leading `z^λ` factor exactly as
   the scalar path does (`_trim_and_normalise:461-476`): `λ = ` number of leading
   `|bₖ| ≤ tol`; drop the first `λ` coefficients of `b` and of **every** numerator
   `aᵢ` (the numerator shares the leading zero, since `aᵢ[1] = c₀·b[1] ≈ 0`), then
   normalise `Q(0) = 1`. Throw only if **all** of `b` is below `tol`. This is the
   correct handling of a spurious pole at the expansion centre for a locally
   regular jet — the stepper evaluates at `t=1`, and the spurious origin
   pole–zero pair cancels.

3. **Remove** the dead `n_near > 1` guard and fix the inverted docstring at
   `SharedPade.jl:221-223`.

## Empirical validation (read-only, `external/probes/sharedpade-offbyone-confirm/c2_reduction_study.jl`)

- **Regime A** — meromorphic `d=2`, shared poles `{2,3}`, true `Q` degree 2,
  requested `m=5`: reduces `5→2`, exact null (`σ_min=0`), recovers the true `Q`
  and both components to machine precision (`eval err 0 / 4e-16`).
- **Regime B** — entire harmonic `[cos,−sin]`, `m=15`, rescaled by `h=0.7`:
  reduces `15→9` (`ρ==9` isolated), `Q(0)≈0` → **cancel `z^1`** → settles at
  degree 8, recovers `[cos,−sin]` to `2e-16 / 5e-9` with **no throw**.

## Consequences

- The `d≥2` vector stepper degrades gracefully on entire / far-from-pole jets
  instead of throwing; the meromorphic (real-workload) case recovers true shared
  poles exactly.
- **Test-accuracy consequence (Rule 5, not a band-aid)**: the correct
  construction's *honest* accuracy on the **entire** harmonic test is `~5e-9`
  (`−sin`, Float64), where the buggy `+1` incidentally hit `<1e-10`. This is the
  method-set "fixed shared-Q residual for pole-free systems" already documented
  in `padetaylor-unk` — neither window is "more correct" on an entire function
  (both place spurious poles). `VS.1.2`'s `−sin` tolerance is updated from `1e-10`
  to its honest value, with a comment citing this ADR and `padetaylor-unk`. *This
  is a correctness-driven tolerance correction, not a relaxation to pass.*
- **`reduce-vs-throw` resolved**: REDUCE (to the supported degree, with `z^λ`
  cancellation); when no degree ≥ 1 remains, return the Taylor numerator over
  `Q = [1]`. Genuinely-zero input still throws before the SVD loop.

## Relation to the deferred dispatch/Pareto layer (future ADR-0028)

This reduction is the **degree axis** of the maintainer's proposed
dual-construction + validated-Pareto dispatch: "reduce to the supported degree"
= "select the supported cell of the Padé table along the denominator-degree
axis." ADR-0028 will add the *numerator-degree axis* (the GGT `(m,m)` vs
Mano–Tsuda `(m−1,m)` cells) and a held-out-residual selection on top of this
reduction. The two compose; this ADR is the prerequisite core.

## References

- `src/RobustPade.jl:406-425` (scalar `ρ==n` reduction), `:461-476` (`z^λ`
  cancellation) — the ground truth mirrored here.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/...md:136` (reduce-`n`-to-`ρ`).
- `references/tex/hermite_pade/ManoTsuda2017_..._MathZ285/hp_arXiv_final.tex:1427`
  ("unique iff rank = m").
- `docs/worklog/066-sharedpade-window-rootcause-resolved.md` (the `+2` resolution).
- `docs/worklog/067-sharedpade-c1-c2-minimal-core.md:36-40` (the shipped `Q=1`
  Taylor fallback for `ρ==0`).
- `external/probes/sharedpade-offbyone-confirm/c2_reduction_study.jl` (validation).
- `docs/adr/0019-shared-denominator-pade.md` (the construction this refines).
