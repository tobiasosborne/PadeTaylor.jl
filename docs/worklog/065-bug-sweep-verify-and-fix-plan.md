# Worklog 065 — subtle-bug sweep: Verify + Synthesize phases + fix plan

**Date**: 2026-06-02
**Author**: Claude Opus 4.8 (1M)
**Beads**: created `padetaylor-d3a` (C1), `padetaylor-bez` (C3),
`padetaylor-qsj` (C4), `padetaylor-3p9c` (C2), `padetaylor-5jvd` (C6);
reconciled `padetaylor-unk`, `padetaylor-0o9`, `padetaylor-6b5`,
`padetaylor-0ln.22`.
**Scope**: run the Verify + Synthesize phases that worklog 064 stopped;
empirically confirm the leading suspect; produce a ranked, non-band-aid
fix plan and capture it in beads.

> **Take-home**: of worklog 064's 7 find-phase candidates, **5 are real**
> (C1, C2, C3, C4, C6) and **2 are refuted** (C5, C7). The CRITICAL lead
> suspect — the `SharedPade.jl:117` Toeplitz off-by-one (C1) — is now
> **empirically confirmed AND its one-character fix validated**: a probe
> shows the `d=1` shared-denominator `Q` diverges from the scalar GGT
> oracle by 5–50 % on transcendental jets (bit-identical on rationals,
> which is why the suite is green), and the `+2` correction restores
> bit-exact agreement with zero regression. Two surprises from the
> adversarial refutation: the test suite's "independent" determinant
> Oracle-2 **shares the identical off-by-one block** (so it can never
> catch C1), and the maintainer-reported *vector* zoom-figure
> discontinuity is best explained by C1 + the already-tracked unguarded
> Stage-2 fill (`padetaylor-6b5`), not by the scalar `iszero` (C5,
> refuted).

## Method

Read-only orchestrated workflow (`wf_1b29c9eb-765`), **no Julia in any
agent** (Rule 7 — the only Julia process was the serial confirmation
probe, run by the main loop): a Verify→Refute pipeline (one confirm
auditor + one adversarial skeptic per candidate, each re-deriving the
math from ground truth and tracing reachability) followed by a synthesis
lead. 13 agents, ~1.4 M subagent tokens. Two confirm agents (C1, C7) hit
a harness `StructuredOutput` transmission fault but completed full
plain-text verdicts, recovered from their transcripts.

The serial confirmation probe (`external/probes/sharedpade-offbyone-confirm/
confirm.jl`) compares `shared_denominator_pade([jet], m).Q` against the
documented oracle `robust_pade(jet, m, m; method=:svd).b` (SP.1.1 /
ADR-0019 contract) on exact-rational vs transcendental jets.

## Verdicts (7 candidates)

| ID | Bug | Verdict | Severity | Thread | Symptom |
|----|-----|---------|----------|--------|---------|
| **C1** | `SharedPade.jl:117` Toeplitz off-by-one (`+1`→`+2`) | **confirmed** (empirical + fix validated) | **critical** | vector | strong |
| **C3** | `PoleField.jl:144` far-root filter (scale-fixing heresy) | confirmed | high | scalar | strong |
| **C4** | `IVPBVPHybrid` `a₂` matched at `s⁻⁷` not `s⁻⁵` (=0 at δ=−1) | confirmed (0.95) | high | scalar | moderate |
| **C2** | `SharedPade.jl:208` rank break `ρ≥m` vs `==` | confirmed | medium | vector | moderate |
| **C6** | `PathNetwork.jl:997` wedge: linear not circular angle | confirmed | low (dormant) | scalar | moderate |
| ~~C5~~ | scalar `_evaluate_pade` bare `iszero` | **refuted** as cause | low | scalar | weak |
| ~~C7~~ | sheet-blind parent selection | **refuted** | low | scalar | weak |

### C1 — empirical confirmation (the keystone)

`robust_pade(jet,m,m;:svd).b` vs `shared_denominator_pade([jet],m).Q`,
BigFloat-256, relative ∞-norm of `ΔQ`:

| input | m=2 | m=3 | m=4 | m=5 | m=6 |
|---|---|---|---|---|---|
| rational `1/(1−z/2)` (control) | 0 | 0 | 0 | 0 | — |
| `exp(z)` | 16.7 % | 10.0 % | 7.1 % | 5.6 % | 4.5 % |
| `log(1+z)` | 50 % | 24 % | 24 % | 24 % | 21 % |

On every transcendental case the recovered numerator degree was one lower
than the oracle (`degP_shared = degP_robust − 1`) — the predicted "`a_m`
driven to ≈0" collapse. Applying `:117` `+1→+2` (validated in a
fully-reverted throwaway edit) drove **all** transcendental cases to
`rel = 0` (bit-exact) with the rational control unchanged at 0.

### Refutations that changed the picture

- **C5 refuted as the discontinuity cause.** Bare `iszero` is *required*
  by ADR-0015's dense-output graceful-degradation contract
  (`PathNetwork.jl:688-692` has no try/catch; a relative `√ε·‖Q‖` guard
  would throw on genuine large-but-finite values just inside a real pole
  and crash every `extrapolate=true` render). And the headline worklog-063
  zoom-figure discontinuity runs through the **vector** Stage-2 fill
  (`VectorPathNetworkStage2._stage2_fill`, no pole guard), not the scalar
  evaluator — already tracked as **`padetaylor-6b5`**. A narrow step-time
  fail-loud gap remains (low), deferred.
- **C7 refuted.** `_nearest_visited` (`PathNetwork.jl:882-894`) is
  sheet-blind Euclidean (true, per ADR-0013:256-261) but **deterministic
  per `rng_seed`** (default 0; `shuffle(MersenneTwister(seed), …)`).
  Labels vary across *seeds* (documented FW shuffle semantics), not
  run-to-run. B4 conflated seed-dependence with hash-order
  non-determinism. The principled improvement is Poisson-disk target
  placement (`padetaylor-pgc/zwh`), not a sheet-aware tie-break.

## Fix plan (execution order C1 → C3 → C4 → C2 → C6)

Each fix: (1) re-derive the contract from ground truth FIRST, no code;
(2) a mutation-proof test that goes RED on the bug + an independent
cross-validation oracle; (3) docs in lockstep (Rule 2). Full per-step
detail is in the bead descriptions.

1. **C1 `padetaylor-d3a`** — `:117` `+1→+2`; raise jet guard `:175`
   `≥m+1`→`≥2m+1` (corrected block reads `c_{2m}`; fail loud); fix the
   `_toeplitz_block` docstring arithmetic `:111`. Test = transcendental
   d=1 oracle (the **only** independent one — the determinant Oracle-2
   shares the bug). Then re-examine `padetaylor-unk`/`-0o9`.
2. **C3 `padetaylor-bez`** — port the landed vector S7 fix
   (`VectorPoleField.jl:250-329`, commit `19b48ba`): z-plane-distance
   filter using h-independent `h_max`, z-closest cluster representative.
   Test = varying-h Weierstrass-℘ lattice (closed-form poles).
3. **C4 `padetaylor-qsj`** — closed-form `a₂ = β²(δ+1)/(δ−2)³`; fix the
   comment algebra, `fig_5:219`, ADR-0014, worklog 040 error-budget.
   Test = FFW eq.6 published-IC oracle replacing the self-referential
   IB.1.1 (a Rule 5 violation) + sympy CAS.
4. **C2 `padetaylor-3p9c`** (depends on C1, same file) — `ρ≥m`→`==`;
   reduce branch for the full-rank case; σ₂/σ₁ isolation test replacing
   the dead `n_near>1` guard; fix inverted docstring `:221-223`.
   **Open decision (reduce-vs-throw)**: full column rank arises on
   genuinely entire/pole-free d≥2 components (`padetaylor-unk`, test
   VS.1.2 green); a naive throw regresses VS.1.2. Lean: reduce to the
   honest lower-degree `Q`, throw only if no degree gives an isolated
   1-D null space. Resolve with the maintainer before coding.
5. **C6 `padetaylor-5jvd`** (low, dormant) — `:997`
   `argmin |rem(θ_sd−off, 2π, RoundNearest)|`. Test = negative-real-goal
   `:steepest_descent` case + brute-force S¹ geodesic oracle.

## Do-not-touch (verified-correct negative space)

RobustPade + classical Padé; the Taylor/AD recurrence; all six Painlevé
RHS; the P_I^(2) companion + closed-form oracles; the Noumi–Yamada
system; the entire adjoint/transpose surface; and **`_upper_block`**
(`SharedPade.jl:132-141`) is correct — do **not** "fix" it to match the
buggy `_toeplitz_block`.

## Files this session

- `external/probes/sharedpade-offbyone-confirm/confirm.jl` — the C1
  confirmation/validation probe (committed; rerun-able).
- `docs/worklog/065-bug-sweep-verify-and-fix-plan.md` — this file.
- Beads: 5 created + 4 reconciled (above).

## References

- `external/chebfun/padeapprox.m:89,92` — the C̃ slice `C=Z(m+2:m+n+1,:)`.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/...md` eq. 2.10.
- `references/hermite_pade/ManoTsuda2017_..._MathZ285.pdf` p.12 eq. 2.6.
- `docs/adr/0019-shared-denominator-pade.md` — the bit-identical d=1 contract.
- `docs/adr/0026-vector-resilient-walk-dense-targets.md` Amd 6/7 — C3 ground truth.
- `docs/adr/0014-ivp-bvp-hybrid.md` + FFW 2017 §3 — C4 ground truth.
- `docs/bug-sweep-2026-06-01/find-*.md` — the 26 find-phase reports.
- `docs/worklog/064-bug-sweep-find-phase.md` — the find phase.
