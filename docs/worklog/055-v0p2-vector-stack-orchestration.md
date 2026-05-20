# Worklog 055 — v0.2 vector stack: orchestrated build V1–V8a

**Date**: 2026-05-20
**Author**: Claude Opus 4.7 (orchestrator) + serial Opus coding subagents
**Epic**: `padetaylor-0ln` (v0.2 — vector / multi-component Painlevé-type systems)
**Scope**: Orchestrated, serial (one subagent at a time, Rule 7), red-green-TDD
build of the v0.2 vector ODE stack — phases V1a through V8a of the 10-phase
plan in `docs/v0p2_plan.md`. **17 of 19 epic children closed this session.**

> **Take-home**: the v0.2 vector pipeline is built and validated end-to-end.
> Shared-denominator Padé → vector Taylor jets → vector stepper → vector
> driver → vector step control → vector path-network + pole extraction. The
> Noumi–Yamada A_n^(1) and Painlevé-I-hierarchy P_I^(2) systems are wired in.
> The Calogero–Moser N=2 smoke test and the A_2^(1)⇒PIV reduction both pass.
> The A_4^(1) pole-field figure (PRD `n≥4` acceptance) ships. **V8b
> (P_I^(2) tritronquée figure) was dispatched as a background agent at
> session end — its status is unknown; see "Pickup point" below.**

## Orchestration method

Parent Opus agent coordinated; each epic child was delegated to one Opus
coding subagent with a full red-green-mutation TDD brief. **Strictly serial**
— one subagent (hence one Julia process) at a time, per CLAUDE.md Rule 7.
After each: the orchestrator re-ran the new test file standalone to verify
GREEN independently (Rule 3 skepticism), fast-forwarded/confirmed `main`,
closed the bead, registered any friction beads. Full `Pkg.test()` was *not*
run (user instruction: it does not OOM but is slow — avoid for code velocity;
the memory `full-pkg-test-got-sigterm-d-twice-on` was corrected accordingly).

## What shipped — 10 new `src/` modules, 4 ADRs, 13 test files, 2 figures

| Phase | Bead | Deliverable | Commit | Tests |
|---|---|---|---|---|
| V1a | 0ln.2 | `SharedPade.jl` — shared-denominator (type-II Hermite–Padé) robust Padé; stacked `dm×(m+1)` block-Toeplitz SVD; `d=1` reduces to GGT `:svd` | `b48c4c2` | 34 |
| V1b | 0ln.3 | Calgo 766 oracle via `ccall` (FORTRAN), upgraded to double precision (~9e-16) | `d8723da`,`54d0fcd` | — (probe) |
| V1c | 0ln.4 | Independent block-Toeplitz-**determinant** shared-Q oracle (SVD-free, exact on `Rational{BigInt}`) | `24a20ff` | — (probe) |
| V1d | 0ln.5 | AAA per-component pole cross-check oracle (self-contained AAA port) | `d17581a` | — (probe) |
| V1e | 0ln.6 | Triple-oracle `SharedPade` test suite + **ADR-0019** | `1c4fcca` | 108 |
| V2 | 0ln.7 | `VectorCoefficients.jl` — first-order vector Taylor jet + **ADR-0020** | `442c976` | 158 |
| V3a | 0ln.8 | `VectorStepper.jl` — vector `pade_step!` (shared-Q) | `240a1d0` | 55 |
| V3b | 0ln.9 | `VectorProblems.jl` — `vector_solve_pade` driver + dense output | `2a1f04a` | 72 |
| V3c | 0ln.10 | `VectorStepControl.jl` — norm-based vector Jorba–Zou + **ADR-0021** | `34696b9` | 69 |
| V4 | 0ln.11 | Calogero–Moser N=2 exact-oracle end-to-end smoke test | `63c72b6` | 110 |
| V5a | 0ln.12 | `NoumiYamada.jl` — A_{2n}^(1) RHS + `NoumiYamadaProblem` + **ADR-0022** | `768d4ad` | 102 |
| V5b | 0ln.13 | A_2^(1) ⇒ scalar PIV backward-compat self-validation | `5cb4de4` | 37 |
| V5c | 0ln.14 | `NoumiYamadaSymmetry.jl` — W(A_{2n}^(1)) Bäcklund + rational oracles | `2d61345` | 175 |
| V6 | 0ln.15 | `PainleveHierarchy.jl` — `painleve_hierarchy(:I,m)` + P_I^(2) 4-vector | `d5e7ed2` | 59 |
| V7 | 0ln.16 | `VectorPathNetwork.jl` + `VectorPoleField.jl` — Stage-1 walk + shared-Q pole extraction | `45d67fa` | 77 |
| V8a | 0ln.17 | Figure: A_4^(1) Noumi–Yamada pole field (PRD `n≥4` acceptance) | `8d179ac` | 50 |

**1106 new GREEN assertions** across 13 new test files (each re-verified
standalone by the orchestrator). Every module ≤200 effective LOC (Rule 6).
Every load-bearing module mutation-proven (2–4 mutations each; all bit).

Two pillar-doc corrections were committed (`4e89216`, `be9dc3c`) — see
"Spec drift caught" below.

## PRD v0.2 acceptance checklist — status

- [x] Shared-Q `robust_pade` mutation-proof vs **3 independent oracles** (V1e).
- [x] Vector `Coefficients`/`Stepper`/`Problems` GREEN; A_2^(1)⇒PIV self-validation (V3, V5b).
- [x] An A_n^(1) pole-field figure, `n≥4`, reproducible from `figures/` (V8a).
- [ ] A P_I^(2) tritronquée pole-field figure vs KKG 2015 (**V8b — in progress**).
- [x] v0.1 tests untouched (additive architecture — no v0.1 `src/` module modified).
- [x] ADRs 0019 (shared-Q), 0020 (vector-y storage), 0021 (vector step control), 0022 (vector problem types). All **Accepted**.

## Spec drift caught (Rule 2 / Rule 3 / hard-won-lesson 9)

1. **Calogero–Moser exact solution wrong in pillar EF** — `docs/v0p2_pillarEF_*`
   §Q6(b) had two arithmetic slips (first-integral coefficient `4→8`, constant
   `4/r₀²→1/r₀²`); its formula `x=±½√(1+4t²)` contradicted its own IC. V4's
   subagent re-derived `x₁(t)=√(1+t²/2)` (for r₀=1), cross-checked to 26–30
   digits vs a from-scratch 256-bit RK4. Doc corrected (`4e89216`).
2. **P_I^(2) tritronquée IC sign slip in pillar C** — §4 printed `y_3` with a
   leading minus; differentiating `u=-∛6|x|^{1/3}` twice for `x<0` gives a
   positive `u''`. Doc corrected (`be9dc3c`).
   The P_I^(2) equation itself (§1–§2) was verified correct vs the KKG `.tex`.

## Hard-won lessons

45. **The shared-Q SVD is underdetermined for entire-function vector
    systems.** A pole-free vector ODE (e.g. harmonic `[cos,−sin]`) leaves a
    small *fixed* shared-Q fit residual (~4e-12 F64, ~3e-21 BF) that does not
    shrink with order or precision. Not a bug — the honest accuracy for
    pole-free systems; the Painlevé/CM targets are meromorphic so unaffected.
    Test tolerances on pole-free sub-cases must carry this headroom.
    (Bead registered.)
46. **The shared-Q Padé degenerates when a component jet is identically
    zero.** `shared_denominator_pade` throws `Q(0)≈0` if a component is `≡0`
    — blocks vector-solving Noumi–Yamada rational *seed* solutions with
    zero components (Type A/B). Generic non-seed solutions are unaffected.
    (Bead registered, P3.)
47. **P_I^(2) tritronquée is a separatrix — forward IVP diverges.** Seeded
    from the asymptotic series, forward integration blows up exponentially
    (KKG used a BVP). v0.2 has no vector/4th-order BVP solver. This is the
    open risk for V8b. (Bead registered, P2.)
48. **Per-node canonical Padé, not the wedge-step Padé.** V7's first cut
    stored the wedge-step approximant (centred at the *parent*, complex-`h`)
    → pole roots mapped through the wrong centre/scale. Fix mirrors v0.1
    `PathNetwork._local_pade`: recompute a shared-Q approximant centred at
    each node with the real step `h`.
49. **Calgo 766 (netlib TOMS) is precision-agnostic FORTRAN.** No `Dp/`
    tree exists, but every float is a bare `real` — `gfortran
    -fdefault-real-8 -fdefault-double-8` promotes the whole algorithm to
    double precision with zero source edits.

## Beads registered this session (friction / deferred work)

- `padetaylor-0ln.20`–`0ln.23` — V7 deferred sophistications: Stage-2 fill,
  edge detection, branch/sheet tracking, shared-Q root-distance wedge.
- `padetaylor-tji` — odd-parity `A_{2n+1}^(1)` Noumi–Yamada system (⇒PV/A_5).
- `padetaylor-qi0` — general-`m` PI-hierarchy member via Lenard recursion.
- `padetaylor-wso` — SP.1.5 `Q(0)≈0` throw lacks an executing assertion (P4).
- shared-Q zero-component degeneracy bug (P3); P_I^(2)-tritronquée-needs-BVP (P2).
- See `bd list --status=open` for the live set.

## Pickup point for the next agent

**V8b (`padetaylor-0ln.18`) was dispatched as a background subagent at
session close — its outcome is unknown.** First action:

1. `git log --oneline -5` — check whether a `V8b`/`painleve_i2_tritronquee`
   commit landed. If yes: re-verify its test standalone, confirm the bead is
   closed, read its report context (it may have shipped a bounded-region
   figure with a documented v1 corner, or STOPPED needing a vector-BVP
   solver — lesson 47). If no commit: the V8b brief is reconstructable from
   this worklog + `docs/v0p2_plan.md` V8b row; re-dispatch.
2. If V8b stopped for lack of a vector/4th-order BVP solver: that solver is
   a new module bead to add to the epic before V8b can be faithful to KKG
   Figs 7.4/7.5. Decide: build it, or ship a bounded-region figure with the
   documented v1 corner (CLAUDE.md Rule 9 — the FFW figures set the precedent).
3. **V9 (`padetaylor-0ln.19`)** — v0.2 docs + release prep — remains open:
   README/CHANGELOG/RESEARCH v0.2 sections, ADR review, HANDOFF refresh,
   the v0.2 epic close-out. This is the last planned phase.
4. The orchestration pattern that worked: serial Opus subagents, one
   red-green-mutation bead each, orchestrator re-verifies standalone +
   closes the bead. Keep it.
