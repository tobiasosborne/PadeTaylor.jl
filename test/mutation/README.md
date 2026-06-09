# Mutation gate — runbook

Automated mutation-testing gate for PadeTaylor.jl's numerical-core hotspots.
Bead `padetaylor-krgy.4`; Tier-1 §1.3 of
`docs/test_corpus/03_hardening_methodology.md`.

This gate **measures** the test-strength the corpus's 89 hand-written
"MUTATION-PROOF" footers *claim* by hand. It perturbs `src/` one mutant at a
time, runs the single targeted test file that footer points at, and reports a
mutation **score** (`killed / total`). A surviving mutant is a real signal: a
test that does not assert the invariant it was supposed to.

It is a **periodic deep gate, NOT a per-commit gate**, and is deliberately
**not** included by `test/runtests.jl` (it mutates `src/`, which invalidates the
precompile cache → one package recompile per mutant).

## Run it

```bash
# Full catalogue (12 mutants, ~3–5 min: each mutant ≈ one package recompile
# (~10 s) + one targeted test file (seconds)):
julia --project=. test/mutation/run_mutation_gate.jl

# Fast subset by 1-based catalogue index (verification / spot-check):
julia --project=. test/mutation/run_mutation_gate.jl 1 4 9
```

Exit code is `0` iff every run mutant matched its catalogued `expect`: `:killed`
mutants stay killed, and KNOWN/tracked `:survived` gaps stay survived. It is
non-zero on any **UNEXPECTED** outcome — a `:killed` mutant that now SURVIVES (a
weakened/regressed test) or a tracked `:survived` mutant that now gets KILLED (a
gap was closed ⇒ flip its catalogue annotation to `:killed`). CI-free local
gate; pair with `git diff src` staying clean.

## Read the score

```
OVERALL  score = 10 / 12 killed (83.3%)
PER HOTSPOT:
  Coefficients: 2 / 3 killed
  RobustPade: 2 / 3 killed
  SharedPadeDefect: 2 / 2 killed
  SheetTracker: 4 / 4 killed
SURVIVORS:
  src/Coefficients.jl:199  [KNOWN/tracked]  — 2nd-order up-resync … (bead padetaylor-98pe)
  src/RobustPade.jl:418    [KNOWN/tracked]  — SVD rank threshold … (bead padetaylor-98pe)
```

- **KILLED** — the targeted test went RED. The test genuinely asserts the
  invariant the mutant broke. This is the healthy outcome.
- **SURVIVED** — the test stayed GREEN despite the perturbation, tagged either
  `[KNOWN/tracked]` (catalogued `expect = :survived`, an accepted gap under a
  bead) or `[NEW GAP → file a bead]` (an UNEXPECTED survivor — the gate FAILS).

### Inaugural baseline (2026-06-09)

The first full run scored **10 / 12 (83.3 %)** with two `[KNOWN/tracked]`
survivors, both filed under **`padetaylor-98pe`**:

- `Coefficients.jl:199` up-resync (`j*u[j]`→`j*u[j-1]`) — a REAL gap:
  `coefficients_test.jl` asserts the value jet `u` but not the derivative jet
  `up`. Fix = add a direct `up`-coefficient assertion, then flip the catalogue
  entry back to `:killed`.
- `RobustPade.jl:418` rank threshold (`>`→`>=`) — a LIKELY EQUIVALENT mutant
  (differs only at an exact-threshold singular value, measure-zero in F64).

## A survivor means a test-strength gap → file a bead

An **UNEXPECTED** survivor (tagged `[NEW GAP]`) is **not** a harness bug; it is a
finding. Either the perturbation is benign (an equivalent mutant) or, far more
likely, **the targeted test does not actually pin the mutated invariant** — a
footer claim that rotted, or a test weakened by a later refactor. Action:

```bash
bd create --type=task --title "mutation survivor: <hotspot> <desc> survives <test>" \
  --description "src/<file>:<line> mutated to <new> survives test/<test>.
          Strengthen the assertion (Rule 4/Rule 5), OR if equivalent, mark the
          catalogue entry expect=:survived with the bead id."
```

Do **not** weaken the catalogue to hide a NEW survivor — fix the test (or, for a
confirmed equivalent mutant, annotate it `:survived` with a bead, as 98pe did).

## Catalogue

The catalogue lives in `run_mutation_gate.jl` (`const CATALOGUE`). Each entry is
`(hotspot, src_file, targeted_test_file, old → new, description, expect)`. The
mutations reuse the ones the manual footers already claim to catch, plus new
single-perturbation operators (sign flips, threshold `<`→`<=`, index slips) on
the same load-bearing lines. To add a hotspot: add a `Mutant` with a
single-occurrence `old` string and the test file whose footer covers it.

## Safety

`src/` is never left mutated. Each mutant reads the original bytes up front,
applies the perturbation, and restores them in a `try/finally` (covers
exceptions and Ctrl-C); after restore the harness asserts `git diff --quiet`
on the file. The restore net is self-tested at startup (a simulated mid-mutant
crash) before any real mutant runs. If a final `git diff src` is ever dirty
after a run, that is a hard error — `git checkout src/<file>` and investigate.

## Why bespoke (not Gremlins.jl)

Gremlins.jl is the maintained successor to the dead Vimes.jl and is registered
(julia 1.10–1.x compat, single warm worker — Rule-7 safe). But its model is
operator-table-driven over **all** of `src/` and applies **one** `test_file`
(default `runtests.jl`, the ~22-min full suite) to **every** mutant. It cannot
map a mutant to a *per-hotspot* targeted test file (our four hotspots need four
different targeted test files), and it generates hundreds of mechanical mutants
rather than re-checking the curated, already-hand-validated catalogue. Running
its default model would cost ~22 min per mutant. The bespoke harness pairs each
curated mutant with its targeted test file → seconds of test + one recompile per
mutant. See the literate header of `run_mutation_gate.jl` for the full rationale.
