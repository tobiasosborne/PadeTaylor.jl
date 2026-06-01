# Worklog 064 — orchestrated subtle-bug sweep (find phase) + beads/git sync

**Date**: 2026-06-01
**Author**: Claude Opus 4.8 (1M)
**Bead**: `padetaylor-bb2` (created this session)
**Scope**: (1) reconcile the local beads Dolt store with the git-tracked
JSONL source of truth; (2) orchestrate a read-only multi-agent audit of
the whole `src/` tree against ground-truth references, hunting a subtle,
intermittent discontinuity bug the maintainer suspects was transcribed
wrong early in development.

> **Take-home**: the maintainer's prime hypothesis (a conjugate-transpose
> `'`-vs-`.'` transcription bug in the robust-Padé QR-reweighting, à la
> the documented July-2018 Chebfun bug) is **DEAD** — verified correct
> with line-for-line citations. But the sweep surfaced a genuinely
> matching CRITICAL elsewhere: an **off-by-one in the vector
> shared-denominator Padé Toeplitz block** (`src/SharedPade.jl:117`),
> flagged INDEPENDENTLY by two agents, that makes `d=1` fail to reduce to
> the scalar GGT `(m,m)` Padé and bites intermittently on transcendental
> (non-rational) jets — exactly the "subtle transcription error at the
> heart, intermittent discontinuity" signature. **Not yet adversarially
> verified** (the Verify/Synthesize phases were stopped on user
> instruction). Verification + a careful (Rule 2, no band-aid) fix is the
> next agent's job.

## Part 1 — beads ↔ git/GitHub reconciliation

- `git remote`: `git@github.com:tobiasosborne/PadeTaylor.jl.git`. `git fetch`
  → local `main` == `origin/main` (`f249e6c`), 0 ahead / 0 behind, tree clean.
- The git-tracked beads artefact is `.beads/issues.jsonl` (the Dolt DB
  under `.beads/embeddeddolt/` is gitignored, per-machine).
- `bd import` (upsert from `.beads/issues.jsonl`) brought the local Dolt
  store from **153 → 167 issues** (+1 memory): the local DB had been
  STALE BEHIND the committed JSONL. Import reconciled it *up* to the git
  source of truth (no clobber risk — JSONL was newer). Post-import: 167
  total / 44 open / 122 closed.
- Working model recorded for future sessions: **JSONL travels through
  GitHub; the Dolt DB is a per-machine working copy** — `bd import` after
  pull, `bd export` before push.

## Part 2 — the bug sweep (orchestrated, read-only)

### Design

A 3-phase background workflow (`/tmp/sweep.workflow.js`, run
`wf_0709eaad-2bd`). All agents **read-only — no Julia** (Rule 7: one
Julia process; here, zero), each owning ONE report file in
`docs/bug-sweep-2026-06-01/` (distinct files ⇒ no write races):

1. **Find** — 26 parallel auditors: A1–A5 core arithmetic
   (RobustPade/Coefficients/StepControl/PadeStepper/SharedPade), B1–B7
   path-network + dense-output + branch/sheet + edge + pole, C1–C4
   BVP/dispatch/drivers/Laplace2D, D1–D5 equation-RHS transcription
   (Painlevé/hierarchy/Noumi–Yamada/Heun/hybrid), E1–E5 cross-cutting
   lenses (adjoint-vs-transpose, indexing, branch-cut, aliasing/RNG,
   comparison/tolerance). Each diffs code against the canonical
   references and cites by `path:line`.
2. **Verify** (NOT RUN — stopped) — every candidate hit by two
   adversarial skeptics (ground-truth re-read + code-trace reachability).
3. **Synthesize** (NOT RUN — stopped) — rank survivors by symptom match.

**The workflow was stopped on user instruction after Find completed.**
All 26 `find-*.md` reports are on disk and committed; there is no
`MASTER_REPORT.md`. Findings below are therefore from the find phase
only — **single-source except where noted as independently corroborated**.

### Headline result — the prime hypothesis is cleared

`A1` (RobustPade) and `E1` (adjoint sweep) both verified, with citations
to `external/chebfun/padeapprox.m:113` and the corrected/erroneous forms:

- `src/RobustPade.jl:443` and `src/SharedPade.jl:234` both use
  `qr(adjoint(...))` — the CONJUGATE transpose, i.e. the **correct**
  post-July-2018 Chebfun form. The buggy plain-transpose (`.'`,
  preserved in the GGT paper's Figure-1 listing `GGT…md:279`) is NOT in
  the code.
- `E1`: there is **no naked `'` adjoint, no `.'`, no `transpose(`** in
  `src/` — the entire complex-data linear-algebra surface is just those
  two `adjoint(qr(...))` sites, both correct. The hypothesis is dead.
- `A1` additionally diffed 16 RobustPade/LinAlg items (Toeplitz build,
  `C̃` 1-based slice, rank threshold `>τ`, diagonal hopping, `tol`-vs-`ts`
  split, normalisation `b₀=1`, classical FW (5.4)/(5.5)) — **all match**.

### The leading suspect (CRITICAL, independently double-flagged, UNVERIFIED)

**`src/SharedPade.jl:117` — `_toeplitz_block` builds the C̃ block one
Taylor coefficient too low.** Flagged by BOTH `A5` (SharedPade audit) and
`E2` (index sweep) independently.

- Ground truth: GGT eq. 2.10 (`GGT…md:57-58`) and `padeapprox.m:89,92`
  (`C = Z(m+2:m+n+1,:)`), faithfully realised by the scalar port
  `src/RobustPade.jl:413-414` ⇒ `C[r,c] = c_{m+r-c+1}`, top-left
  `C[1,1] = c_{m+1}`.
- Code: `idx = m + rr - cc + 1` (line 117) ⇒ `blk[1,1] = c_m`, one order
  low. At `m=1`: scalar `C̃ = [c₂, c₁]`; SharedPade `blk = [c₁, c₀]`. The
  block encodes the matching window `z^m…z^{2m-1}` (a Mano–Tsuda
  `(m-1,m)` system) instead of GGT's `z^{m+1}…z^{2m}` `(m,m)`. Fix per
  A5: index should be `m + rr - cc + 2`.
- Internal inconsistency: `_upper_block` (`:132-141`) IS correct, so the
  recovered top numerator `a_m` is driven to ≈0 by the shifted null
  vector → result collapses to numerator degree `m-1`.
- Mechanism (intermittency): invisible on EXACT-rational inputs (the true
  Q annihilates every power, both windows agree — and **every**
  `shared_pade_test.jl` oracle feeds `_ratio_jet` exact rationals, which
  is why the suite is green); on a truncated transcendental jet — the
  real workload in `VectorStepper`/`VectorProblems` — the `(m,m)` and
  realised `(m-1,m)` approximants genuinely differ, `m` is chosen
  adaptively per step, so the pole estimate (and solution value) error
  varies step-to-step ⇒ intermittent discontinuity in the vector
  trajectory. Confidence 0.9 (A5), 0.9 (E2).
- **Caveat for the next agent**: this lives in the **v0.2 vector stack**,
  not the scalar path-network of the original FW figures (whose Padé runs
  through `classical_pade_diagonal`, verified correct by A1). It plausibly
  underlies the vector P_I^(2) / zoom-figure artefacts (worklog 063), but
  is *probably not* the "appeared at the beginning of dev" scalar bug.
  Confirm scope before fixing; do NOT band-aid (Rule 2). Note ADR-0019 +
  the SharedPade docstring assert a bit-identical `d=1` reduction the code
  does not deliver on non-rational inputs — that contract is the thing to
  re-derive against Mano–Tsuda before changing the index.

### Other notable findings (single-source, UNVERIFIED — full text in reports)

HIGH:
- `B7` — scalar `PoleField` still uses the `|t*| ≤ radius_t` far-root
  filter that ADR-0026 root-caused as the cause of *intermittently
  empty/sparse pole fields*. Strong symptom match for the scalar side.
- `D5` — `IVPBVPHybrid` ships `a_2 ≈ −0.22208` where the Fig-5 family
  value is exactly 0; test IB.1.1 enshrines the wrong value so cannot
  catch it. (Adjacent to the known `padetaylor-ykg` v1 corner.)
- `E5` / `A5(MEDIUM)` — `shared_denominator_pade` rank break uses
  `ρ ≥ m_cur` instead of `ρ == m_cur`; isolation guard is dead code for
  `d ≥ 2` (accepts a full-column-rank stack → spurious denominator).

MEDIUM (intermittency-flavoured, in the walk / stepper / sheet code):
- `A4` — scalar `_evaluate_pade` uses bare `iszero` while the vector
  stepper uses a relative `√ε·‖Q‖` pole guard → near-pole steps divide by
  a tiny-but-nonzero denominator.
- `B1` / `E3` — `:steepest_descent` wedge selection uses an un-wrapped
  angular distance → picks the wrong ray near goal directions ≈ ±π.
- `B4` — sheet-blind Euclidean parent selection makes `cross_branch`
  sheet labels shuffle/order dependent (intermittent sheet labelling).

### Reassuring negative space (verified CORRECT against references)

`A1` RobustPade + classical Padé (16 items); `E1` the entire
adjoint/transpose surface; `A2` Taylor/AD recurrence (1st & 2nd order,
scalar & vector); `D1` the six Painlevé RHS term-by-term; `D2` the P_I^(2)
companion system + all closed-form oracles + tritronquée IC; `D3` the full
Noumi–Yamada RHS + Bäcklund + rotation + rational oracle. The maintainer
can trust these are not the bug.

## Files this session

- `docs/bug-sweep-2026-06-01/find-*.md` — 26 audit reports (committed).
- `docs/worklog/064-bug-sweep-find-phase.md` — this file.
- `HANDOFF.md` — new top section pointing at the leading suspect.
- `/tmp/sweep.workflow.js` — the workflow script (NOT committed; resumable
  via `Workflow({scriptPath, resumeFromRunId: "wf_0709eaad-2bd"})` to run
  the Verify/Synthesize phases, which return the 26 Find agents from
  cache and only run the new phases).

## Pickup for the next agent (in order)

1. **Resume Verify+Synthesize** from the cached find phase (cheap — the 26
   finds are journalled) to get adversarial confirmation + a ranked
   `MASTER_REPORT.md`. Or skip straight to (2) for the lead.
2. **Adversarially verify the `SharedPade.jl:117` off-by-one** with a
   single serial Julia experiment: build `robust_pade(jet, m, m)` and
   `shared_denominator_pade([jet], m)` on a TRUNCATED TRANSCENDENTAL jet
   (e.g. `exp`/`tan` Taylor at order 2m+1, NOT an exact rational) and
   compare the denominators `Q`. If they differ, the bug is real; if the
   `+2` index makes them agree, that's the fix — but first re-derive the
   intended `d≥2` matching window against Mano–Tsuda eq. 2.6 (don't break
   the genuine simultaneous-Padé case to fix the `d=1` reduction).
3. Triage the HIGH/MEDIUM list; `B7` (scalar pole-field far-root filter)
   is the best candidate for a *scalar-side* intermittent symptom.

## References

- `external/chebfun/padeapprox.m:89,92,113` — the C̃ slice + the
  documented adjoint historical bug.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/...md:57-58`
  (eq. 2.10), `:279` (Figure-1 listing with the erroneous `.'`).
- `references/hermite_pade/ManoTsuda2017_..._MathZ285.pdf` p.10,12
  (eqs. 2.2 / 2.6 — the simultaneous-Padé block the `d≥2` case targets).
- `docs/adr/0019-shared-denominator-pade.md` — the bit-identical-`d=1`
  contract the code violates on non-rational jets.
- `docs/v0p2_pillarA_hermite_pade_findings.md:410` — the spec line that
  itself carries the off-by-one (the code faithfully implemented a buggy
  spec).
