# CLAUDE.md — PadeTaylor.jl

Rule numbers are stable: `src/`, `test/`, `scripts/` and the ADRs cite them.

- **Law 1 — Ground truth before code.** Open the reference
  (`references/markdown/`, `external/chebfun/padeapprox.m`), the ADR and the
  current file first. Cite `path:line`. Never paraphrase the papers from memory.
- **Law 2 / Rule 10 — Docs ship with the code.** Top-of-file docstring written
  as a chapter; an ADR in `docs/adr/` for each algorithmic decision; public API
  in `README.md`.
- **Rule 1 — Fail loud.** Throw with a suggestion; never signal failure with
  NaN, zero or `nothing`.
- **Rule 2 — Find the root cause.** No band-aids, and never relax a test
  tolerance to get to green.
- **Rule 3 — Verify, don't trust.** `git log` and the files are authoritative;
  memory, handoffs and subagent reports are not.
- **Rule 4 — Mutation-prove.** RED → GREEN, or port faithfully and then perturb
  the impl, confirm RED, restore. Cross-validate against an independent oracle.
- **Rule 5 — Tests assert known-correct values.** "Didn't throw" is not a test.
- **Rule 6 — ≤ 200 LOC per file**; split when over.
- **Rule 7 — One Julia process at a time.** Never run `julia` in parallel.
- **Rule 8 — Tracker is `bd` (beads).** Never `bd edit`. The git-tracked
  `.beads/issues.jsonl` is canonical (ADR-0035): `bd export -o
  .beads/issues.jsonl` and commit it with the code. Run
  `./scripts/git-hooks/install.sh` once per machine.
- **Rule 9 — Senior-engineer-grade only.** No "good enough for now"; a known
  limit becomes a deferred bead naming the condition that forces the work.
- **Rule 11 — Gates are local, no CI.** `scripts/quality_gate.sh fast` before
  every commit touching `src/`, `full` before push; single test files while
  iterating. Exactly 2 `@test_broken` is expected.
- **Rule 12 — No outreach** to the original authors.
- **Rule 13 — Re-read this file** at session start and after compaction.
- **Session close.** Close beads, export and commit the JSONL, then commit and
  push freely once the suite is green.
