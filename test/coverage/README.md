# Coverage tool — runbook

Local code-coverage measurement for PadeTaylor.jl's `src/`.
Bead `padetaylor-krgy.7`; Tier-1 §1.4 of
`docs/test_corpus/03_hardening_methodology.md`.

Coverage tells you which `src/` lines **no test exercises** — a *necessary,
not sufficient* floor. An uncovered numerical branch is a guaranteed blind
spot. The converse does **not** hold: a covered line is **not** a tested line.
That second axis (test *strength*) is measured by the mutation gate
(`test/mutation/run_mutation_gate.jl`). Read the two together; **do not chase
100%** — defensive fail-loud guards, `verbose`-print blocks, `show` methods,
and arb-precision-only paths are legitimately cold at Float64.

It is a **periodic measurement, NOT a per-commit gate**, and is deliberately
**not** included by `test/runtests.jl`. Coverage.jl is a **tool** dependency,
installed into a throwaway temp env by the tool itself — it is **never** added
to the package `[deps]`/`[extras]`/`[targets]` (mirrors how the mutation gate
avoids a suite dep).

## Run it

Two paths; same tool. The tool **never regenerates** unless you ask it to with
`--run` — by default it only **aggregates** `.cov` already on disk.

```bash
# 1. FULL — generate fresh coverage, then aggregate (the canonical ~22-min path):
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'   # drops src/*.cov
julia test/coverage/run_coverage.jl                         # aggregate + report

# 2. AGGREGATE-ONLY (fast) — summarise src/*.cov already present (e.g. a recent
#    Pkg.test(coverage=true) left them):
julia test/coverage/run_coverage.jl

# Or let the tool drive the full regenerate itself (slow; one Julia process — Rule 7):
julia test/coverage/run_coverage.jl --run
```

`Pkg.test(coverage=true)` leaves one `src/<File>.jl.<pid>.cov` per source file.
Those are **gitignored** (`.gitignore:65 *.cov`) — never commit them. The tool
does not delete them (so you can inspect per-line `.cov` afterwards); delete
them by hand once you have the numbers (`rm src/*.cov`).

## Read the report

```
OVERALL  2543 / 2644 executable lines  (96.18%)

FILES BELOW 100% (worst first):
  Diagnostics.jl                       2/   13  15.4  %
  VectorPathNetwork.jl                86/  108  79.6  %
  ...
WORST 12 — UNCOVERED LINE NUMBERS (triage these):
  src/Diagnostics.jl  →  209-219
  ...
```

- **OVERALL** — covered / total *executable* lines (Coverage.jl's notion of
  executable; blank lines, comments, and pure docstrings are excluded).
- **FILES BELOW 100%** — every file with an uncovered executable line, worst
  first. 100%-covered files are omitted (they meet the floor; lean on the
  mutation gate for their *strength*).
- **WORST N + uncovered line numbers** — compressed `a-b` runs of the
  uncovered lines, the input to triage.

## Triage: covered-floor ≠ tested

For each uncovered region, **read the source** and classify:

- **(a) genuine BLIND SPOT** — real logic no test runs → file a bead
  (`file:lines` + what a test should assert).
- **(b) intentionally-uncovered** — a defensive fail-loud `throw`/`error`
  guard (Rule 1), a `verbose`-print / `show`-display block, an unreachable
  branch, or an arb-precision-only path not hit at Float64.

The full triage of the inaugural run lives in
`docs/test_corpus/coverage_snapshot.md` (committed). When you re-measure,
diff against that snapshot: a file that *dropped* below its snapshot %
means a test stopped exercising a branch (a regression in coverage), and a
*new* uncovered region is a new blind spot to triage.

A blind spot found here is a coverage signal; pair it with the mutation gate,
which finds the complementary failure — a line tests **do** run but do **not**
pin. Neither alone is sufficient.

## How it works (implementation)

`run_coverage.jl` (≈104 code LOC):

1. installs Coverage.jl into a throwaway temp env (`Pkg.activate(mktempdir())`
   + `Pkg.add("Coverage")`) — never the package env;
2. `Coverage.process_folder("src")` aggregates every `src/*.cov`
   (via `Base.invokelatest`, since Coverage is imported at runtime → world-age);
3. computes per-file covered/total/uncovered-lines, sorts ascending, and
   prints the report.

**Rule-1 guard:** `process_folder` does *not* return empty when no `.cov` is
present — it returns one all-`nothing` entry per source file. The tool catches
that by asserting `tc > 0` (some line was actually hit); `tc == 0` ⇒ a loud
error telling you to generate the `.cov` first, rather than a fake-looking
0/N. The tool never touches `src/`.
