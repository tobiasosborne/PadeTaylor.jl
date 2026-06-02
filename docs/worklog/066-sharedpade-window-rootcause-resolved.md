# Worklog 066 — SharedPade window: root cause resolved (math settled)

**Date**: 2026-06-02
**Author**: Claude Opus 4.8 (1M)
**Beads**: `padetaylor-d3a` (C1, window), `padetaylor-3p9c` (C2, reduction)
**Scope**: per the maintainer's directive ("find the root cause of the spec
contradiction; deep-research it; correctness is priority; build only once the
math is clear"), settle — mathematically, before any code — whether the
SharedPade Toeplitz window should be `+1` or `+2`, and why the two project
specs disagree.

> **Resolution (unanimous: 4 derivations → adversarial reconcile → 3
> independent judges, conf 0.92–0.93; sympy/numpy-verified; cross-checked
> against the worklog-065 BigFloat probes).** The contradiction is a
> **target-approximant conflation**, not an index typo: the code implements
> Mano–Tsuda's *off-diagonal* `(m−1, m)` simultaneous-Padé window while the
> binding contract (ADR-0019 / SP.1.1 / docstring) demands the GGT *diagonal*
> `(m, m)`. The correct window for the contract is **`c_{m+1}` (the `+2`)** at
> every `d`. The over-determined stacking construction is **sound but needs
> graceful degree-reduction** for the `d≥2` entire / far-from-pole case.

## Root cause (precise, primary-source-cited)

- Mano–Tsuda's minimal type-II / simultaneous-Padé system (`hp_arXiv_final.tex`
  eq:rTm `:1060-1080`, eq:lspa `:1405-1428`, eq:spa `:414-417`) uses block
  `A^i_m(n, m+1)` with Toeplitz **subscript `j = m`** → top-left `a^i_m = c_m`
  (the `+1`), `n = m/d` rows per block (a **square** `m×(m+1)` system), and
  caps `deg P_i^{(0)} ≤ m−1` (`:411`). It is therefore the **off-diagonal
  `(m−1, m)`** approximant, matching `O(w^{2m})`.
- `docs/v0p2_pillarA_hermite_pade_findings.md:407-411` **faithfully transcribes
  that block** (`jets[i][m+r-c+1]`, top-left `c_m`) — the transcription is
  *correct*. But `§3/§7` (`:158`, `:409-411`) silently redefine `n := m`
  ("`m` rows per component"), grafting GGT's square-block **row count** onto
  Mano–Tsuda's `c_m` **window subscript**, and assert the result equals the GGT
  **diagonal `(m,m)`** `C̃` (top-left `c_{m+1}`, window `z^{m+1}…z^{2m}`,
  residual `O(z^{2m+1})`). It does not — the two differ by exactly one
  numerator degree.
- Compounding: the `SharedPade.jl:111` docstring arithmetic is itself wrong —
  it claims the slice `C = Z[m+2:2m+1,:]` gives `C[r,c] = c_{m+r-c}` (top-left
  `c_m`), but `Z[i,j]=c[i-j+1]` makes it `c_{m+r-c+1}` (top-left `c_{m+1}`).
  That 1-based→subscript mis-conversion is what made the `c_m` formula *look*
  bit-identical to `RobustPade`'s `C̃`. It is not (sympy: the `+2` block, not
  the `+1`, equals `RobustPade`'s `C`).
- The determinant **Oracle 2** (`external/probes/shared-pade-determinant-oracle/
  oracle.jl:160`) uses the **identical `+1`** and its docstring (`:92-98`)
  repeats the false `c_{m+r-c} = GGT C̃` equivalence — so it is a verbatim
  mirror of the bug and **cannot independently catch it**. This is why the
  triple-oracle suite stayed green.

## Why exact rationals are blind (and transcendentals are not)

For an exact shared-denominator rational with numerators of degree `≤ m−1`, the
true `Q` annihilates **both** windows, so `+1` and `+2` (and even the
over-determined stack) recover the same `Q` (sympy: all give `Q=[1,−1/2,1/3]`).
That is why every `SP.*` oracle (all exact rationals) is green and blind to the
bug. On a **transcendental** jet there is no exact null vector; the smallest-σ
right vector is window-dependent, and `+1` yields the `(m−1,m)` answer
(numerator degree collapses by one; 5–50 % denominator error — worklog-065 F1).

## Construction verdict — sound, needs graceful reduction (the C2 axis)

The over-determined `dm×(m+1)` stacking is a **valid** realisation: an exact
1-D null space exists iff the components share a genuine degree-`m` denominator
(Mano–Tsuda "unique iff rank = m", `:1427-1428`); otherwise the smallest-σ
vector is a legitimate least-squares "best shared `Q`". For entire / pole-free
components the stack is numerically rank-deficient at the stepper's large `m`
(numpy: eff-rank ≈ 9–10 ≪ m+1=16 at m=15). The fix is **graceful reduction**,
exactly GGT's own prescription (`GGT2013 md:136`: `σ_n=0` ⇒ reduce `n→ρ`): reduce
`m` to the honest isolated-1-D-null-space degree (a real `σ₂/σ₁` isolation test
replacing the `ρ≥m_cur` break + dead `n_near>1` guard), throw only if no degree
isolates. This axis is independent of the window but must be re-derived **after**
the `+2` fix (the window shift moves every `A_full` row).

## Math-justified fix (no code yet)

Binding spec = ADR-0019 / docstring / SP.1.1 (GGT diagonal). The findings doc is
the one to amend. The fix:
- **(A) window** `SharedPade.jl:117` `+1 → +2` (necessary at every `d`).
- **(C) docstring arithmetic** `:27,:34,:105-112` `c_{m+r-c} → c_{m+r-c+1}`.
- **(D) DO NOT touch `_upper_block`** (`:132-141`) — sympy-verified it already
  recovers the correct degree-`m` numerator under `+2`.
- **(E) Oracle 2** `oracle.jl:160` `+1 → +2` + fix its docstring `:92-98`;
  the only window-independent `d=1` oracle is the transcendental
  `robust_pade(jet,m,m;:svd)` comparison (added this session as `SP.1.7`).
- **(F) deferred / coupled C2** graceful reduction (above).
- **Jet-length guard (`:175`)** — *declined* the research's "raise to `≥2m+1`":
  `robust_pade` **zero-pads** short jets to `2m+1` (`RobustPade.jl:378`,
  padeapprox.m:62), so raising the guard would make `shared` throw where
  `robust_pade` computes, **breaking** the bit-identical `d=1` contract on
  short jets (and the `+2` block already zero-pads gracefully via
  `1 ≤ idx ≤ length(c)`). Keep `≥ m+1`. *(A Rule-3 correction of the agents'
  suggestion, verified against `robust_pade`'s actual padding.)*

## Two genuine forks the math cannot decide (project decisions)

1. **Which approximant** — GGT diagonal `(m,m)` [`+2`, contract as written] vs
   Mano–Tsuda literal `(m−1,m)` [`+1`, amend ADR-0019/SP.1.1/docstring +
   `_upper_block`]. Both are valid; the contract picks `(m,m)`.
2. **Row count for `d≥2`** — Mano–Tsuda's minimal `n=m/d` rows (square-perfect)
   vs the code's `m` rows (over-determined). ADR-0019 chose `m` rows for `d=1`
   bit-identity; whether the over-determined least-squares stack is the desired
   `d≥2` realisation (it changes the answer on transcendental/entire input) vs
   the square system is unresolved by math alone.

These two forks have an **order-vs-conditioning trade-off** (the over-determined
`(m,m)` is higher-order but rank-deficient away from poles; the square `(m−1,m)`
is lower-order but better-conditioned) — which is exactly why the maintainer's
proposed **dual-construction + dispatch + Pareto-selection** is well-posed: it
turns both forks into a per-step, data-driven choice and *unifies* them with the
graceful reduction (reduce = select a lower table cell). Evaluating that design
(an ADR) is the next step.

## References

- `references/tex/hermite_pade/ManoTsuda2017_..._MathZ285/hp_arXiv_final.tex`
  eq:rTm `:1060-1080`, eq:lspa `:1405-1428`, eq:spa `:414-417`, deg cap `:411`.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/...md:116-117` (eq.2.10),
  `:136` (reduce-`n`-to-`ρ`).
- `docs/v0p2_pillarA_hermite_pade_findings.md:407-411` (the conflation).
- `docs/adr/0019-shared-denominator-pade.md` (the binding diagonal contract).
- `src/SharedPade.jl:111,117,132-141,175`; `src/RobustPade.jl:378,413-414`.
- `external/probes/shared-pade-determinant-oracle/oracle.jl:92-98,160`.
- `docs/worklog/065-...md` (F1–F4 empirical probes); workflow `wf_f81c8ec1-ee9`.
