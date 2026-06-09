# Worklog 073 — Test-hardening sweep (epic `padetaylor-krgy`)

**Date:** 2026-06-09. **Outcome:** the corpus is now a multi-technique
verification suite. **Final gate:** full `Pkg.test()` GREEN at
**9219 pass / 10 broken / 0 fail** (≈25 min). The 10 broken are the intentional
`@test_broken` auto-flip markers (unchanged; investigate only FAILs).

## What this session did

Took the landed ground-truth corpus (epics `p1v0`, `25og`) and hardened it into a
rigorous benchmark + verification suite, driven by a deep-research methodology
synthesis (`docs/test_corpus/03_hardening_methodology.md`, 5-angle fan-out, 23/25
claims verified 3-0, 2 refuted). All **14 planned beads** delivered; the epic
remains open only for the one deliberately-deferred follow-on (`krgy.15`).

## The multi-technique stack now in place

| Bead | Technique | Artifact | Result |
|------|-----------|----------|--------|
| krgy.2 | Static analysis | `test/quality_test.jl` (Aqua+JET+ExplicitImports) | 8 stale imports cleaned; → `fmf8` |
| krgy.3 | Metamorphic relations | `test/metamorphic_{symmetry,consistency}_test.jl` | MM.1–5 + NY-Bäcklund; → `xhjw`, `krgy.15` |
| krgy.1 | Property-based | `test/property_test.jl` (Supposition.jl) | 5 invariants; → `jznu` |
| krgy.5 | Certified ball oracles | `test/certified_oracle_test.jl` + ADR-0029 | ℘ via Arb; CO.1 certifies the FW IC golden |
| krgy.14 | Differential testing | `test/differential_test.jl` (in-file RK4) | method-independent; no new dep |
| krgy.6 | Convergence / MMS V&V | `test/convergence_test.jl` | order law corrected to `2⌊order/2⌋+1` |
| krgy.4 | Automated mutation | `test/mutation/run_mutation_gate.jl` | 10/12 baseline; → `98pe` |
| krgy.8 | Accuracy regression | `test/accuracy_ledger_test.jl` + `accuracy_baseline.toml` | drift-catching ledger |
| krgy.9 | Snapshot/approval | `test/oracle_provenance_test.jl` + doc 04 | enforced provenance contract |
| krgy.7 | Coverage | `test/coverage/run_coverage.jl` + snapshot | 96.18%; → `ftxn` |
| krgy.11 | Perf regression | `benchmark/` (BenchmarkTools) | WSL2 slack calibrated = 0.50 |
| krgy.12 | Alloc + type stability | `test/type_stability_test.jl` (AllocCheck+JET) | hot path proven stable + alloc-free |
| krgy.13 | Formal-methods decision | ADR-0030 | DEFER-with-condition |
| krgy.10 | Quality-gate runner | `scripts/quality_gate.sh` + CLAUDE.md | fast/full/deep, strictly serial |

The in-suite gates (krgy.1/.2/.3/.5/.6/.8/.9/.12/.14) run under `Pkg.test()`; the
periodic out-of-suite tools (mutation, coverage, perf) have their own runbooks and
are NOT in `runtests.jl`. `scripts/quality_gate.sh` ties them together serially.

## Findings (8 follow-on beads — all real, none wrong-arithmetic)

The sweep confirmed the **numerical core is clean**; every finding is a missing
edge-guard, an over-broad claim, or a test-strength gap:

- **`xhjw`** (P2 bug) — `solve_pade` silently returns a degenerate one-node
  trajectory on a descending real span (Rule-1 violation). Found by metamorphic
  MR-06 (forward∘reverse).
- **`jznu`** (P2 bug) — d=1 `SharedPade` vs scalar `:svd` diverge in two degenerate
  regimes; the SP.1.1 "bit-for-bit" docstring over-claims. Found by property INV-2.
- **`fmf8`** (P3) — `Statistics` should be a weakdep (extension-only). Found by Aqua.
- **`98pe`** (P3) — mutation survivors: a real gap (`Coefficients` up-resync untested
  by `coefficients_test`) + a likely equivalent mutant (RobustPade rank threshold).
- **`ftxn`** (P3) — 7 coverage blind spots (COV-1..7), mostly degenerate/recovery
  branches; 2 may overlap `jznu`/`98pe`.
- **`krgy.15`** (P2) — deferred metamorphic relations MR-02/07/09 (need the complex
  path-network or pole-bridge arcs).
- Two non-bug **corrections folded into the work**: the diagonal-Padé order law
  `2⌊order/2⌋+1` (methodology §2.4 fixed in lockstep), and a TOML-test-env fix.

## Process lessons (recorded via `bd remember`)

1. **`static-gate-test-order-dependence`** — package-introspection gates (JET/Aqua)
   are sensitive to which extensions earlier test files loaded; assert `≤` not `==`
   on allow-listed reports.
2. **`standalone-vs-pkgtest-env-gotchas`** — three confirmed ways a test passes
   standalone (`--project=.`) but fails under `Pkg.test`: extension-load ordering;
   Pkg.test-only extras (Arblib/Makie) not standalone-runnable at all; an undeclared
   stdlib (`using TOML`) resolving via LOAD_PATH standalone but absent in the sandbox.
   **Any new in-suite test that adds a `using` MUST be gated under full `Pkg.test`.**
3. A `tail`-piped `Pkg.test()` masks Julia's real non-zero exit code; and a naive
   word-adjacency parser of Julia's two-row `Test Summary` produces a green-looking
   lie — both caught and fixed (the runner parses positionally now).

## Orchestration method

Serial, one Julia process at a time (Rule 7). Each coding bead → one Opus subagent
(claim → build → mutation-prove → verify standalone/temp-env → report); research/
recon (metamorphic-relation mining) → Sonnet read-only agent. The orchestrator
verified every report against the repo (Rule 3), ran the regression gates, raised
the finding-beads, and committed/pushed per bead. Full-suite gates were batched:
dep-adding / src-touching / in-suite beads gated under `Pkg.test`; standalone tools
(mutation, coverage, perf, runner) verified by running them directly.

## Follow-on backlog

Open, tracked: `xhjw`, `jznu` (P2 bugs); `fmf8`, `98pe`, `ftxn`, `krgy.15` (P3/P2
enhancements). The pre-existing corpus bugs (`q0yq`/`53tu`/`61um`/`v1ub`) and their
`@test_broken` markers are unchanged. The mutation gate's 2 tracked survivors
(`98pe`) keep it at exit 0 until the `Coefficients` up-resync test is strengthened.
