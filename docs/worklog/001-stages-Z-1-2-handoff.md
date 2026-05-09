# Worklog 001 — Stages 0-1 + Phases Z, 1, 2 handoff

**Date**: 2026-05-09
**Author**: Claude (single-session work)
**Scope**: Stages 0 (research) and 1 (design) plus Phases Z, 1, 2 of
the Stage 2 implementation plan.

> Worklog shards are frozen snapshots of substantive iterations. Read
> them when "git blame says I changed this line, but why?" — that's
> what these shards are for. Be honest about frictions and dead ends;
> those are the load-bearing parts.

## Context

Project goal: build PadeTaylor.jl, a general-purpose Taylor–Padé IVP
solver for analytic ODEs in the Fornberg–Weideman tradition. PRD
specifies Julia first → 100% working → then spec the TS mirror.

User explicitly mandated: "ground-truth BEFORE coding"; "red green TDD
always"; "register beads for any problems"; "Senior engineer grade
only." Discipline inherits from `scientist-workbench`'s `CLAUDE.md`,
specialised here for a Julia ODE-solver project.

User scope decisions taken during the session:
- **TIB VPN** for paper acquisition (Cloudflare blocks fall back to
  user-fetch via Windows Downloads; user dropped 5 PDFs in the
  session's first turn).
- **No author outreach** — Fornberg, Weideman, Fasondini, Reeger are
  all retired or inactive.
- **"We are alone here"** — no external code base / authors to defer
  to; cross-validation between Julia and TS impls is load-bearing.
- **Octave** installed mid-Phase-2 to capture `padeapprox.m` oracle
  values.

## What changed

### Stage 0 — RESEARCH.md (1165 lines)

- §1.1 FW 2011 read cover-to-cover (line-cited): four-layer
  decomposition, standard params `(order, h) = (30, 0.5)`,
  Toeplitz-vs-SVD contrast at §5.2, failure modes, FW Table 5.1
  cross-validation targets pinned.
- §1.2-1.5 stub citations for FW 2014, FFW 2017, RF 2014, Willers
  (Willers 1974 not acquired — ACM Cloudflare-blocked, low priority).
- §2 GGT 2013 read cover-to-cover with QR-reweighting trick from
  `padeapprox.m` lines 236-241 explicitly noted as beyond the paper.
- §3 (TIDES, TI.jl, TS.jl, ATOMFT, Jorba-Zou) populated by a sonnet
  research subagent with rigorous citations.
- §3.3 **EMPIRICALLY validated**: a sonnet probe agent ran a 5-test
  Julia program confirming `TaylorSeries.jl::Taylor1{Arb}` works at
  order 80 with 256-bit precision; no hand-rolled coefficient layer
  needed. Probe lives at `external/probes/taylorseries-arb/`.
- §4 verified-arithmetic ODE survey populated by another sonnet
  subagent (with friction recorded — this agent over-stepped its
  brief and edited RESEARCH.md directly; preserved content).
- §5 bigfloat-SVD landscape: confirmed Arblib has no SVD via source
  inspection; `GenericLinearAlgebra.svd` is the path.
- §7 resolves all 7 PRD Stage 1 design questions and all 6 open
  questions with concrete recommendations.

### Stage 1 — DESIGN.md (700 lines)

Granular plan for 9 phases of implementation. Each phase is:
- a single module + its tests
- specifies the test corpus with oracle source
- includes the mutation-proof procedure
- has a LOC budget

The plan deliberately decomposes the algorithm into the **four
layers** of FW 2011 (Coefficients / RobustPade / StepControl /
PadeStepper) so each is independently testable. This decomposition is
ADR-0001.

### ADRs

- ADR-0001: four-layer architecture.
- ADR-0002: bigfloat-SVD via `GenericLinearAlgebra` one-sided Jacobi.
  Why relative-accuracy on small SVs is load-bearing for arb-prec
  rank counting at GGT's `tol`-thresholded block boundaries.
- ADR-0003: Pkg.jl extensions for `Arblib` + `CommonSolve`.

### CLAUDE.md replaced

The generic beads-init CLAUDE.md was replaced with a 13-rule
PadeTaylor-specific discipline. Includes hallucination-risk callouts
that capture the project's specific traps:
- Don't paraphrase FW 2011 / GGT 2013 from memory.
- GGT 2013 Algorithm 2 ≠ FW 2011 §5.2.
- `TaylorSeries.jl` works for Arb empirically — don't hand-roll.
- `Arblib.jl` has NO SVD.
- `padeapprox.m` does QR-reweighting beyond GGT 2013 Algorithm 2.

### Phase Z — scaffold

- `Project.toml` with deps + weakdeps + extensions; Julia 1.10
  baseline; `[targets]` test config including `Random`.
- 6 module stubs + umbrella `src/PadeTaylor.jl`.
- 6 placeholder test files + `test/runtests.jl`.
- `bd init` with `--prefix=padetaylor --database=padetaylor`
  (workaround for the dotted-DB-name bug).

### Phase 1 — `LinAlg.pade_svd`

80 LOC + 4 tests + 1 follow-up shape test = 5 tests. Mutation-proven.

Discovered mid-Phase-2 that `padeapprox.m` needs the null-vector
column from `full=true` SVD; added `full::Bool=false` kwarg as a
follow-up commit before starting Phase 2. Did not paper over the gap
with a workaround in RobustPade.

### Phase 2 — `RobustPade.robust_pade`

210 LOC + 8 testsets + Octave oracle capture = 35 individual @test
assertions. Mutation-proven.

The Octave oracle capture is at `external/probes/padeapprox-oracle/`
and produces `oracles.txt` (copied to `test/_oracles.jl`). Captures
the 6 test cases (exp(2,2), exp(20,20), log(1.2-z)(20,20),
tan(z⁴)(20,20), 1+z²(1,1), noisy 1/(1-z)).

The (μ, ν) reductions match GGT 2013 Figs. 2, 4, 6 verbatim:
- exp(20,20) → (7, 7) — diagonal-stripe
- tan(z⁴)(20,20) → (20, 16) — 4 Froissart doublets removed
- 1+z² (1,1) → (0, 0) — defect-1 collapse
- noisy 1/(1-z) → (0, 1) — noise threshold recovery

## Why these choices

### Why decompose into four layers and not three or five

ADR-0001 §"Why decompose this way" — each FW 2011 section maps 1:1 to
a module, so a bug in one layer surfaces as a bug independently
checkable against that layer's reference. A monolithic stepper would
couple coefficient generation with Padé conversion, losing this
property.

### Why GenericLinearAlgebra and not a hand-rolled Jacobi

ADR-0002 — `GenericLinearAlgebra.jl` is the only production-grade
Julia library that provides generic SVD over `BigFloat`; it's
CI-maintained; `Arblib.jl` ships no SVD whatsoever (verified by full
source inspection); LAPACK Demmel-Kahan absolute-error guarantees are
not load-bearing-strong-enough at arb-prec for the GGT
`tol`-thresholded rank counting.

### Why port-and-verify against Octave's `padeapprox.m`

CLAUDE.md Rule 4 + the user-confirmed choice "Use Octave locally to
run padeapprox.m and capture outputs." The Octave path is the
faithful oracle. Tests assert function-value equality (load-bearing)
and only loosely per-coefficient (acknowledging GGT 2013 §7's
ill-posedness for blocks at `defect > 0`).

### Why the (20,20) BigFloat exp test does NOT reduce to (7,7)

The reduction at Float64 happens because `1/k!` for `k > 16` falls
below `1e-14 · ‖c‖₂`, so the SVD detects rank deficiency. At BigFloat
256-bit precision, `1/40! ≈ 1.2e-48` is well above the BigFloat tol
floor of `1e-74`, so no rank deficiency is detected. **This is
correct behaviour** — higher precision admits more accurate Padé
approximants — but it required updating the test from (7,7) to
(20,20) at higher precision.

## Frictions surfaced

### `bd init` rejects dotted DB names

`PadeTaylor.jl` (with the `.jl` suffix) → Dolt schema name
`PadeTaylor.jl` rejected as invalid. The init command silently
auto-detects from directory; positional args are ignored. Workaround:
`bd init --prefix=padetaylor --database=padetaylor --skip-agents
--skip-hooks`.

### `GenericLinearAlgebra` pirates into `LinearAlgebra.svd!`

When `using GenericLinearAlgebra` is in scope, `LinearAlgebra.svd(A::
Matrix{BigFloat})` works because `GenericLinearAlgebra` adds methods
to `LinearAlgebra.svd!` directly (intentional piracy in their design).

This caused the *first* mutation-proof attempt for Phase 1's BigFloat
dispatch claim to **fail to bite**: substituting
`s/GenericLinearAlgebra.svd/LinearAlgebra.svd/` did not RED any test
because both routes hit the same Jacobi impl. The correct mutation
that bites is `Matrix{Float64}(A)` downcast (loses precision required
for the 1e-65 reconstruction tolerance).

**Lesson**: when porting Julia code that touches stdlib, beware that
generic packages may be pirating into stdlib namespaces; choose
mutations that change the *arithmetic*, not just the *function name*.

### Octave `i` ≠ Julia `im` for imaginary unit

My initial `capture.m` format string was `'%.18e + %.18ei'` — the
trailing `i` was meant as MATLAB's imaginary unit suffix. Julia uses
`im`. The captured `_oracles.jl` was unparseable, with confusing
`UndefVarError: i not defined` from Julia trying to multiply by a
`i` variable.

**Fix**: changed format string to `'%.18e + %.18eim'`. Verified by
counting `im` occurrences in the rewritten oracles file.

### `padeapprox.m` needs `full=true` SVD; thin SVD doesn't expose null vector

DESIGN.md initially specified `pade_svd(A) -> (U, S, Vt)` returning
the thin SVD. While reading `padeapprox.m` line 109 (`b = V(:, n+1)`)
during Phase 2 prep, realised: for an `n × (n+1)` matrix, the thin V
is `(n+1) × n` — the (n+1)-th column doesn't exist! Need full SVD.

**Fix**: added `full::Bool=false` kwarg to `pade_svd` as a Phase 1
follow-up commit before starting Phase 2 implementation. Wrote a new
test (1.1.3b) confirming `full=true` produces the right shape. Did
NOT band-aid the API in RobustPade — fixed the abstraction properly.

### GGT 2013 §7 ill-posedness manifests empirically

For Phase 2 test 2.1.3 (log(1.2-z) (20,20) reduces to (10,10)), the
per-coefficient diff between Octave-output and Julia-output is **~1e-3
absolute** — orders of magnitude larger than I expected. The
*functional* diff (evaluate r(z) at sample points) is 1e-15 to 1e-9
(loosens approaching the branch cut at z=1.2). This is GGT 2013 §7's
ill-posedness exactly: "Padé approximation problems are sometimes
ill-posed... an arbitrarily small perturbation could fracture the
block." The (10,10) reduction sits at exactly such a fracture.

**Fix**: removed the per-coefficient match from test 2.1.3 with a
documented note; tightened the functional-evaluation match instead.
The honest test is the function-value match.

### `randn` not reproducible across language RNGs

For test 2.1.6 (noisy 1/(1-z) recovery), Octave's `randn('state',
42)` and Julia's `MersenneTwister(42)` produce different noise
samples. The captured a/b oracles depend on the noise realisation;
useless for cross-validation against an independent RNG. **Fix**:
test only the structural outcome (μ=0, ν=1) using a Julia-side RNG
sample at the same scale; document the limitation.

### Test 2.1.4 pole-position check skipped

Initial test plan included pole-position equality between our
denominator's roots and the Octave-captured pole list. Skipped at
implementation time because:
- `Polynomials.jl::roots` accuracy with high-degree real polynomials
  is itself a question (and we'll need to revisit for Phase 4).
- The 8-ray symmetry of tan(z⁴) makes the SVD/QR pole positions
  numerically delicate; small sign differences put poles on different
  rays in equally-valid implementations.

The (μ, ν) match — 4 Froissart doublets removed — is the
load-bearing structural assertion. Test passes on (μ, ν) match alone.
**Filed as a follow-up bead**: Phase 4 will need to deal with
`Polynomials.jl::roots` anyway, and once that's solved we can revisit
this test with a stricter pole-position check.

### Subagent over-stepped brief

The verified-arithmetic ODE survey subagent (Stage 0 Phase 0.3) was
told "Output: a single markdown document of ≤600 words at
/tmp/...". It wrote directly to RESEARCH.md §4 instead. Content was
good; behaviour was off-spec.

**Lesson for future subagent prompts**: be explicit that "do NOT
modify <file>" applies even when the agent thinks it's helping.
Future subagent briefs should sandbox their writes to a probe-area
temp file.

## Acceptance

Per-phase acceptance criteria (DESIGN.md §4):

| Phase | LOC | Tests | Mutation-proof | Acceptance |
|---|---|---|---|---|
| Z | scaffold | 12 placeholder | n/a | umbrella loads cleanly via `Pkg.test()` |
| 1 | 80 | 14 | ✅ Float64-downcast bites 1.1.3 + 1.1.4 | 14/14 GREEN |
| 2 | 210 | 35 | ✅ bypass QR-reweighting bites 2.1.2 + 2.1.4 | 35/35 GREEN |

Total: 61/61 GREEN at commit `45e3e90`.

## Pointers

- [`HANDOFF.md`](../../HANDOFF.md) — top-of-repo handoff for the
  next agent.
- [`RESEARCH.md`](../../RESEARCH.md) — Stage 0 deliverable.
- [`DESIGN.md`](../../DESIGN.md) — Stage 1 deliverable; §4 is the
  granular execution plan.
- [`CLAUDE.md`](../../CLAUDE.md) — discipline rules.
- [`docs/adr/0001-four-layer-architecture.md`](../adr/0001-four-layer-architecture.md)
- [`docs/adr/0002-bigfloat-svd-via-genericlinalg.md`](../adr/0002-bigfloat-svd-via-genericlinalg.md)
- [`docs/adr/0003-extensions-pattern.md`](../adr/0003-extensions-pattern.md)
- `external/probes/taylorseries-arb/SYNTHESIS.md` — empirical
  evidence the Arb-Taylor coefficient layer works.
- `external/probes/padeapprox-oracle/capture.m` — Octave oracle script.
- `references/markdown/<file>/*.md` — six papers, marker-converted,
  citable by line number.
