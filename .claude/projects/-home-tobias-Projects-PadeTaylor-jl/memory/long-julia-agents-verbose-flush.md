---
name: long-julia-agents-verbose-flush
description: Long-running Julia subagents must use verbose mode + eager flush + log tailing, never poll exit status
metadata:
  type: feedback
---

When delegating a long-running Julia task (figure render, big walk, full
`Pkg.test`) to a subagent, set it up so progress is observable. The
beginner mistake — and one that already cost ~18 min of wasted polling
in the headline-figure render — is: launch `julia ... > log` in the
background and poll `BashOutput` / exit status. Julia's stdout is
line-buffered (or worse, block-buffered) under redirection, so the log
file stays *empty* until completion; every poll learns nothing.

**The principle: fast feedback, verbose mode, eager flush.**

Concrete pattern the subagent prompt should mandate:

- In the Julia script, sprinkle phase markers with explicit
  `flush(stdout)`:
  ```julia
  @info "Phase X: building Stage-1 walk (h=$h, order=$N)"
  flush(stdout); flush(stderr)
  ```
  After each significant boundary (walk built, sector solved, panel
  rendered, PNG written) — not every step.

- Wrap the invocation with `stdbuf -oL` to force line-buffered stdout
  even under redirection:
  ```bash
  stdbuf -oL -eL julia --project=figures figures/render.jl 2>&1 | tee /tmp/render.log
  ```

- The agent should **not** poll bash exit status at short intervals.
  Either (a) `wait` synchronously on the background bash (which
  blocks until it exits, with `tee` showing live output to the
  terminal), or (b) `tail -f /tmp/render.log &` with `wait` on the
  julia PID, or (c) BashOutput polls at minute+ intervals after each
  phase marker appears in the log — never seconds.

**Why:** Without flush, every short-interval poll wastes tokens and
returns no signal. With flush + verbose markers + line-buffered tee,
each poll *learns* the current phase, the agent can naturally throttle
to phase boundaries, and the user / orchestrator can `tail -f` the log
to see live progress without interrupting the agent.

**How to apply:** every subagent prompt that includes "run the figure
render" or "run the full Pkg.test" must include the
phase-marker-with-`flush` directive + the `stdbuf`/`tee` invocation
pattern. Also see [[full-pkg-test-got-sigterm-d-twice-on]] for the
related "single-file isolated tests during TDD, full Pkg.test only at
end-of-phase" rule.
