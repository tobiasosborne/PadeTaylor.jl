# Test-Hardening Methodology — research synthesis + adoption roadmap

> Driver for **epic `padetaylor-krgy`** ("test-hardening sweep"). Produced
> 2026-06-09 by a deep-research fan-out (5 angles → 25 sources → 117 claims →
> 25 adversarially verified, 23 confirmed 3-0, 2 refuted). This file is the
> durable record; the live tracker is the epic + its 14 children.
>
> **Skepticism note (Rule 3):** every recommendation below is tagged with its
> verification status. Two plausible claims were *refuted* in verification and
> are recorded as such — do not act on them. Where a tool choice rests on
> judgment rather than a 3-0 vote, it says so.

## 0. Where we stand (inventory, 2026-06-09)

Evidence-based gap map from grepping the 90 `_test.jl` files + 25 `_oracle_*.jl`:

**Already strong (the corpus backbone — keep, don't rebuild):**
golden/ground-truth oracles (25 `_oracle_*.jl`); fail-loud error paths
(35 files / 190 `@test_throws`); manual mutation-proofing (89 files,
narrative); cross-precision BigFloat/Arb (23 files `setprecision`);
`@test_broken` auto-flip known-bug markers (10 files).

**Absent or only partial (→ the epic's children):** see the table in §5.

The package's *numerical core is clean* (corpus-v1/v2 found edge-guard bugs,
never wrong arithmetic). So hardening ROI is in **breadth of input** (property/
fuzz/metamorphic), **strength measurement** (automated mutation, coverage),
**oracle rigor** (certified balls, live differential), and **process gates**
(static analysis, V&V convergence, perf) — not in more hand-picked value pins.

## TIER 1 — industry-standard, adopt-now

### 1.1 Static analysis + linting — `krgy.2` (P1, highest ROI)
- **What/why:** finds defects no example test asserts — method ambiguities,
  unbound type params, stale `[compat]`, type instabilities, undefined exports.
- **Tools (confirmed 3-0):** **Aqua.jl** verifies `Project.toml` hygiene incl.
  stale deps and detects method ambiguities; **JET.jl** uses Julia type
  inference — `@report_call` for call/dispatch errors, `@report_opt` for
  runtime-dispatch/type-stability. **ExplicitImports.jl** for implicit/stale
  imports. JuliaFormatter in check-mode only.
- **Wiring:** one `test/quality_test.jl`; expect first-run noise — Aqua
  `detect_ambiguities` may trip on Arblib/TaylorSeries; triage into fix-now vs
  documented allow-list. Single process, fast, no CI needed.
- Sources: github.com/JuliaTesting/Aqua.jl, github.com/aviatesk/JET.jl.

### 1.2 Property-based testing — `krgy.1` (P1)
- **What/why:** generates hundreds of inputs + shrinks failures to a minimal
  counterexample; pins INVARIANTS where goldens pin VALUES. Complement, not
  replacement.
- **Tool (confirmed 3-0):** **Supposition.jl** — Julia implementation of
  Hypothesis-style choice-sequence PBT with integrated shrinking; plugs into the
  existing corpus via `Test.AbstractTestSet`. (PropCheck.jl is the older
  alternative; Supposition.jl preferred.)
- **Wiring:** `test/property/` with 3-5 core invariants (robust-Padé reconstructs
  its Taylor coeffs; shared-Q `d=1` ≡ scalar bit-identically; step∘reverse ≈ id).
  Fixed RNG seed (honors the determinism contract in `PadeTaylor.jl`).
- Source: github.com/Seelengrab/Supposition.jl.

### 1.3 Automated mutation testing — `krgy.4` (P2)
- **What/why:** the 89 hand-mutation-proofed files are excellent but *unmeasured*
  — a refactor can silently weaken a test. Automated mutation injects src/
  perturbations and reports the **mutation score = % of mutants killed**
  (confirmed 3-0).
- **Tool:** **Gremlins.jl** (2024). ⚠️ **Vimes.jl is unmaintained** (v0.1.0,
  2019) — confirmed 3-0; Gremlins is its successor. PBMT (Property-Based Mutation
  Testing) is a research follow-on; its "coverage-quality-for-property-suites"
  positioning was **REFUTED 0-3** — do not overclaim it.
- **Wiring:** periodic deep gate (slow, not per-commit), scoped to numerical-core
  hotspots already hand-proofed. Survivors → beads.

### 1.4 Coverage as a necessary-not-sufficient floor — `krgy.7` (P2)
- **Tool:** Coverage.jl, run locally with `--code-coverage=user`; no codecov, no
  CI. One-shot measure + triage of uncovered src/ branches; pair with mutation
  score (which measures test *strength*). Do not chase 100%.

## TIER 2 — numerical / scientific-computing-specific (the high-value core)

### 2.1 Metamorphic testing — `krgy.3` (P1, single highest-leverage Tier-2 item)
- **What/why:** the principled attack on the **oracle problem** — for an ODE/
  special-function solver the exact answer is usually unknown, but a known input
  transformation forces a predictable output change. A metamorphic relation (MR)
  is an invariant between transformations; MRs are **NECESSARY-not-sufficient**
  (both confirmed 3-0).
- **LAW-1 obligation:** the Painlevé-specific MRs (scaling, Bäcklund, parity,
  conjugate-pole symmetry) **must be derived from `references/` PDFs**
  (FW2011/FFW2017/Noumi–Yamada), not paraphrased from memory — flagged open
  question. Safe oracle-free starters: forward∘reverse = id; step additivity
  (one h-step ≡ two h/2-steps); conjugate symmetry of pole fields for
  real-coefficient ODEs.
- **Tool:** plain `Test.jl` + the §1.2 generators; no new dep.
- Sources: Metamorphic-testing survey PMC7252536; en.wikipedia.org/wiki/Metamorphic_testing.

### 2.2 Certified / verified oracle values — `krgy.5` (P2)
- **What/why:** current goldens are high-precision *floats* — trustworthy but not
  *certified*. Verified numerics produce a rigorous enclosing **ball** guaranteed
  to contain the true value; the test asserts the computed value lies inside.
- **Tool (confirmed 3-0):** **Arb / ArbNumerics.jl** ball (midpoint-radius)
  arithmetic — tighter than plain interval arithmetic (which suffers the
  dependency problem / overestimation) and auto-propagates error (Arb's own
  randomized suite even exposed MPFR 3.1.3 bugs). TaylorModels.jl `(P, e)` for
  truncation capture.
- **SCOPE (reconciled):** ball arithmetic certifies the *closed-form ORACLE*
  (Weierstrass/Painlevé/Padé-of-known-fn), **not** the package's SVD-based Padé
  output — Arblib has no SVD (ADR-0002). Pattern: certify the oracle as a ball,
  assert the SVD value lands inside it.
- ⚠️ "Taylor-model integration is a turnkey single-step verified integrator" was
  **REFUTED 0-3** — treat Taylor models as truncation-capture only.
- Sources: arxiv.org/pdf/1611.02831 (Taylor models); Arb docs.

### 2.3 Live differential / back-to-back testing — `krgy.14` (P2)
- **What/why:** cross-validation today is partly *frozen Float64 literals*
  (capability map B22) — a frozen literal can rot and can't re-derive. Live
  differential testing runs an *independent* implementation in-process every run
  (Hatton: independent codes disagree more than expected).
- **Tool:** **TaylorIntegration.jl** (pure-Julia high-order Taylor IVP) in-process
  for smooth/no-pole arcs — single Julia process, Rule-7 OK; mpmath/Mathematica
  remain out-of-process capture scripts. Distinct from §2.2 (bounds the value)
  vs here (agreement with a second code path).

### 2.4 Convergence-order / MMS / observed-order-of-accuracy — `krgy.6` (P2)
- **What/why:** the V&V gold standard for solver CODE verification, distinct from
  value pinning. Method of Manufactured Solutions injects a forcing so the exact
  solution is known by construction; refine h, fit log-error vs log-h, confirm
  the **observed order** matches theory. Catches order-degradation bugs a
  single-h test passes.
- **Framing (confirmed 3-0):** Roy 2005 distinguishes **code verification**
  (math solved right — MMS/Richardson) from **solution verification** (error in
  a given calculation). We are doing *code verification*.
- **Tool:** plain `Test.jl` + refinement loop + linear fit. Assert the observed
  single-step order for the stepper; spectral (exponential) convergence for the
  Chebyshev BVP. **CORRECTION (krgy.6, as-built):** the naive textbook law
  `err ~ h^(order+1)` is WRONG for this package — `pade_step_with_pade!` builds a
  *diagonal* `(m,m)` Padé with `m = order÷2`, which matches the series to degree
  `2m`, so the single-step order is `2·(order÷2)+1` (verified m=1→3, 2→5, 3→7),
  half the naive value. The order is observable only at MODEST orders (2–6) and
  inside the error window `[~1e-12, 1e-1]`; at the production order 30 the leading
  term `h^31` is sub-eps and the slope is unmeasurable. EVEN manufactured
  solutions collapse the Padé trimmer `(m,m)→(m-1,m-1)` — use both-parity anchors.
  See `test/convergence_test.jl` header.
- Sources: Roy JCP 2005 (aoe.vt.edu cjr_jcp); Oberkampf & Roy, *V&V in Scientific
  Computing* (Cambridge); ASME MMS code-verification paper.

### 2.5 Numerical-accuracy regression tracking — `krgy.8` (P2)
- A binary pass/fail at fixed tol hides slow accuracy drift. Store per-case
  best-achieved relERR in a versioned baseline; assert `new ≤ baseline·(1+slack)`.
  Catches silent regression AND improvement. Pairs with §3.1 snapshot workflow.

### 2.6 Golden-master snapshot + approval workflow — `krgy.9` (P2)
- Formalize the existing `_oracle_*.jl` + capture-script pattern (B22): every
  golden states source (closed-form/mpmath/Mathematica), tolerance (never
  bit-exact for derived floats), and regeneration command; regenerate → review
  diff → approve. ReferenceTests.jl optional for array/image snapshots.

## TIER 3 — frontier (FLAGGED ASPIRATIONAL) — `krgy.13` (decision bead)
- Formal methods / proof assistants (Lean 4, Coq), SMT-backed verification,
  proof-carrying code, refined runtime contracts. High-cost / high-assurance;
  **no claim here reached a 3-0 vote** (lower confidence). For a small no-CI team
  this is a *conscious DEFER-with-condition* (Rule 9), not silent omission. The
  repo already has Lean 4 tooling, so the cheapest entry — IF pursued — is one
  narrow lemma behind a load-bearing recurrence. Gappa/Flocq (Melquiond) is the
  reference point for FP-rounding proofs.

## EFFICIENCY (secondary)
- **Perf-regression gate — `krgy.11` (P3):** BenchmarkTools.jl — minimum/median
  over enough samples, compare *ratio* to a committed local baseline with a wide
  slack (>1.3× = investigate) to survive WSL2 noise without CI. Open question:
  exact slack on this box.
- **Allocation + type-stability — `krgy.12` (P3):** AllocCheck.jl can *prove*
  allocation-free; `@inferred` / JET `@report_opt` assert type stability (same
  JET setup as `krgy.2`). Target the per-step inner loop; budget = O(1)/step.
- **Runner — `krgy.10` (P2):** single serial entry point (FAST tier: suite +
  static analysis; DEEP tier: mutation + full BigFloat + perf). ⚠️
  **ReTestItems.jl REJECTED** — its worker-process parallelism (confirmed 3-0)
  violates Rule 7. Stays strictly serial. Depends on `.2`, `.7`, `.11`.

## 5. Bead roadmap (epic `padetaylor-krgy`, 14 children)

Prioritized by correctness ROI (report ordering: static-analysis → metamorphic →
property → certified goldens → MMS → mutation → perf → Lean).

| Bead | Technique | Tier | Pri | Tool |
|------|-----------|------|-----|------|
| krgy.2  | Static analysis + lint gate | 1 | P1 | Aqua.jl, JET.jl, ExplicitImports.jl |
| krgy.3  | Metamorphic-relation layer | 2 | P1 | Test.jl + generators (MRs from PDFs) |
| krgy.1  | Property-based testing | 1 | P1 | Supposition.jl |
| krgy.4  | Automated mutation testing | 1 | P2 | Gremlins.jl (not Vimes) |
| krgy.5  | Certified/ball oracles | 2 | P2 | ArbNumerics.jl, TaylorModels.jl |
| krgy.14 | Live differential testing | 2 | P2 | TaylorIntegration.jl |
| krgy.6  | Convergence-order / MMS | 2 | P2 | Test.jl + Richardson fit |
| krgy.7  | Local coverage | 1 | P2 | Coverage.jl |
| krgy.8  | Accuracy regression tracking | 2 | P2 | in-suite baseline ledger |
| krgy.9  | Snapshot/approval workflow | 1 | P2 | _oracle pattern + ReferenceTests.jl |
| krgy.10 | Local quality-gate runner | — | P2 | serial script (⊥ ReTestItems) |
| krgy.11 | Perf-regression gate | eff | P3 | BenchmarkTools.jl |
| krgy.12 | Alloc + type-stability | eff | P3 | AllocCheck.jl, JET @report_opt |
| krgy.13 | Formal-methods frontier | 3 | P3 | Lean 4 (DECISION: likely defer) |

**Suggested first sprint (biggest correctness ROI, lowest cost):**
`krgy.2` (static gate — hours, finds real defects) → `krgy.3` (metamorphic —
deepest lever) → `krgy.1` (property-based). These three need no out-of-process
oracle and run inside the normal suite.

## Refuted claims (do NOT act on these)
1. *"PBMT is a coverage-quality metric for property/invariant suites rather than
   example-based suites."* — REFUTED 0-3 (arxiv 2301.13615).
2. *"Taylor-model ODE integration is a turnkey single-step verified method needing
   no separate enclosure step."* — REFUTED 0-3 (Springer 978-3-319-63501-9_1).

## Open questions (carried into the beads)
- Which MRs are provably valid for the Painlevé transcendents? → derive from
  `references/` PDFs (`krgy.3`, Law 1).
- Gremlins.jl maturity / PBMT killing strength? (`krgy.4`)
- Certified goldens vs ADR-0002 no-SVD — resolved: certify oracle side only
  (`krgy.5`).
- Perf-gate slack tolerance on this WSL2 box? (`krgy.11`)
