# CLAUDE.md — guidance for AI agents working on PadeTaylor.jl

If you are an agent (Claude Code, an SDK harness, or a downstream tool)
landing in this repo, **read this top to bottom every session**. After
a context compression, re-read. The rules drift out of working memory
faster than you think; that's why they're numbered.

This project's discipline is inherited from `scientist-workbench` —
read its `CLAUDE.md` if you need extended rationale; the rules below
are the same shape, specialised for a Julia ODE-solver project.

## The Laws

Two laws. Read first, applied always.

**Law 1 — Ground truth before code.** Before writing any code: open
the ADR (or write it first), open the affected files (verify their
current shape — never trust memory or a prior conversation summary),
open the canonical reference (paper / `padeapprox.m` / textbook
section). Cite the local source by path in commit messages and ADRs
(`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:128`, not "the FW paper I
remember"). The PDFs in `references/` are the canonical mathematical
ground truth; do **not** paraphrase them from memory.

**Law 2 — Docs in lockstep with code.** Every code change ships
paired doc updates *in the same edit session*:

- new module: docstring-as-chapter at the top of the file (literate
  programming, like scientist-workbench's `tools/<name>/tool.ts`);
- new algorithmic decision: an ADR under `docs/adr/NNNN-title.md`;
- new external dependency: justification line in `Project.toml` or
  the relevant ADR;
- new feature in the public API: docstring + entry in `README.md`'s
  API table (when README exists).

Code shipped with stale docs is incomplete work, full stop.

## The Rules

Numbered, non-negotiable. Re-read after compaction.

0. **Laws 1 & 2 apply.** Ground truth before code. Docs in lockstep.

1. **Fail fast, fail loud.** Throw exceptions with `suggestion` /
   `detail` text; never return zero, NaN, or `nothing` to signal
   failure silently. Numerical breakdowns (singular `C̃`, `Q(z) = 0`
   at evaluation point, step-size convergence failure) throw with
   a suggestion. Crashes with context beat truthful-looking lies.

2. **All bugs are deep.** No band-aids, no "temporary fixes."
   Investigate root causes. A test failure that "looks like a
   precision issue" is more often a sign-flip or off-by-one in the
   port; verify via mutation-proof + cross-validation before
   patching tolerances.

3. **Skepticism.** Verify subagent output, previous-session claims,
   and your own memory against the current state of the repo.
   `git log` and `Read` are authoritative; conversation context is
   not. RESEARCH.md is a snapshot — when in doubt, re-open the PDF.

4. **TDD discipline (two valid shapes).**
   - **Spec-from-scratch:** classic RED → GREEN → refactor.
   - **Port-and-verify:** port the algorithm faithfully (GGT 2013
     Algorithm 2 verbatim from `padeapprox.m`; Jorba-Zou step
     formula verbatim from §3.2 eq. 3-8), capture invariants in
     tests, **mutation-prove** the tests catch regressions (perturb
     the impl, confirm RED, restore), cross-validate against an
     independent oracle (Chebfun's MATLAB output, closed-form
     ℘-function values, `TaylorIntegration.jl`'s step-size output)
     when available.
     Mutation-proving replaces the literal "RED first" step.
   The discipline is "tests have caught a real regression," not
   "the test was committed before the impl."

5. **"Runs without errors" is not a passing test.** Every test
   asserts an invariant against a known-correct value (FW Table 5.1
   reference numbers; closed-form `WeierstrassP` to high precision;
   pre-captured `padeapprox.m` outputs); a test that asserts only
   "didn't throw" is broken.

6. **`≤ 200 LOC per file`.** When a module exceeds this, split.
   `LOC` excludes blank lines and pure docstring blocks.

7. **No parallel Julia agents.** Julia precompile-cache contention
   makes parallel `Pkg.add` / `Pkg.precompile` brittle. **One**
   Julia process at a time across all subagents you spawn. Read-only
   subagents (literature, source reading, no `julia` invocation) are
   fine in parallel.

8. **Beads is the only tracker.** `bd create / update / claim /
   close`. No TodoWrite, no TaskCreate, no markdown TODO lists.
   Run `bd ready` at session start; `bd close <id1> <id2> ...` at
   the end. Never use `bd edit` (it opens `$EDITOR` and blocks).

9. **Senior-engineer-grade only.** No "good enough for now"; no
   "we'll fix this in v2 if it bites." If a v1-acceptable corner
   exists, document it as a deferred bead with the exact condition
   that would force v2 work. The reference impl that becomes
   "the reference" must be honest about its limits.

10. **Literate programming.** Source files are exposition.
    Top-of-file docstring expands into multiple paragraphs explaining
    *why* the code is shaped the way it is — what ground truth it
    embodies, what pitfalls motivated each defensive check, which
    references it derives from. A fresh reader should read
    `RobustPade.jl` top-to-bottom like a chapter, not piece intent
    together from scattered comments. If you find yourself writing
    terse `# add 1` comments, rewrite as prose or delete entirely.

11. **No GitHub CI.** Quality gates run locally: `julia --project=.
    -e 'using Pkg; Pkg.test()'`, mutation-proof checks, manual
    cross-validation against FW Table 5.1. Do not propose
    `.github/workflows/`; do not file "add CI" beads. Failure-noise
    from automated CI is worse than zero signal at this stage.

12. **No author outreach.** Per project decision: the original
    authors (Fornberg, Weideman, Fasondini, Reeger) are retired or
    inactive. Do not draft outreach emails. We have FW 2011 §5.2
    + GGT 2013 Figure 1 = enough to reconstruct end-to-end.

13. **Repeat rules.** Re-read this file at session start, after
    `/clear`, and after any context compression. The agent that
    re-reads catches drift the agent that doesn't ship.

## Hallucination-risk callouts

Sharp pre-emptive warnings about specific mistake categories that look
right but aren't.

- **Don't paraphrase FW 2011 / GGT 2013 from memory.** Open
  `references/markdown/<file>.md` and cite line numbers. The papers
  have specific algorithmic choices — `(order, h) = (30, 0.5)`,
  `tol = 1e-14`, `||b||₂ = 1` normalisation — that are easy to misremember.

- **GGT 2013 Algorithm 2 ≠ FW 2011 §5.2.** They solve the same
  underlying problem with different normalisations and different
  failure responses. Don't mix idioms. We adopt GGT 2013 (per
  RESEARCH.md §2 and §7.1).

- **`Taylor1{Arb}` works** (RESEARCH.md §3.3 empirical probe). Don't
  hand-roll a coefficient layer. If you find yourself writing
  generic Taylor arithmetic, you've drifted; use `TaylorSeries.jl`.

- **`Arblib.jl` has NO SVD** (RESEARCH.md §5.1 confirmed by source
  inspection). Don't search for `Arblib.svd` — it doesn't exist.
  Route `Matrix{Arb}` through `BigFloat` for SVD, accepting the
  precision-loss caveat documented in ADR-0002.

- **`padeapprox.m` does QR-reweighting beyond GGT 2013 Algorithm 2.**
  Lines 278–280 of `external/chebfun/padeapprox.m` are NOT in the
  paper. Port these explicitly; without them, near-tol blocks misbehave
  (RESEARCH.md §2.2).

- **Don't relax test tolerances to make tests pass.** If the test
  fails, find the bug, not the wrong tolerance. FW Table 5.1's
  numbers are reference values — meet them or document a precise
  reason for the gap.

- **Mutation-proof every load-bearing test.** Perturb the impl,
  confirm the test goes RED, restore. Untested tests don't catch
  regressions.

## Practical guidance

- Substrate is **Julia 1.10+**. Run tests via `julia --project=. -e
  'using Pkg; Pkg.test()'`. No build step.
- ADRs go in `docs/adr/`. Worklog (when we add it) goes in
  `docs/worklog/`.
- PDFs in `references/`; markdown extractions (via `marker_single`)
  in `references/markdown/<file>/`.
- External clones (Chebfun, Arblib.jl, etc.) go in `external/` and are
  `.gitignore`'d (reproducible via shallow clone).
- Empirical probes go in `external/probes/<probe-name>/` with their
  own `Project.toml`.
- Issue tracking: `bd` (beads). `bd ready` to find work; `bd show
  <id>` for detail; `bd close <id1> <id2> …` to close.

## Session close

When the session is winding down:

1. `bd close <id1> <id2> ...` — close completed issues.
2. Commit work with messages that cite the relevant
   `references/markdown/<file>.md:<lines>` or ADR.
3. If a non-obvious lesson surfaced, note it in a new bead or a
   docstring (don't just hold it in memory).
4. **Pushing to remote is allowed any time.** Push freely once work
   is committed and the suite is GREEN — no need to wait for an
   explicit instruction. (Superseded 2026-05-14; the earlier
   "do not auto-push" rule is retired.)

## Tool-of-last-resort

If the laws conflict with a fast path: choose the laws. "Just ship
and fix docs later" has been retired as a working mode here.
