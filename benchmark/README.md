# Performance-regression gate — runbook

Local, noise-aware performance-regression gate for PadeTaylor.jl's hot paths.
Bead `padetaylor-krgy.11`; EFFICIENCY section of
`docs/test_corpus/03_hardening_methodology.md:163-167`.

The corpus pins **correctness** (values, errors, invariants) but nothing pins
**speed**. A refactor that quietly doubles the per-step Padé cost — an
accidental `O(n²)` loop, a lost `@inbounds`, a type instability that boxes the
inner jet — leaves every correctness test GREEN while the integrator crawls.
This gate catches that: it benchmarks the hot paths with **BenchmarkTools.jl**
and compares each one's current timing to a **committed baseline**
(`benchmark/baseline.toml`), flagging a regression by **ratio** with a wide,
calibrated slack band.

It is a **periodic, on-demand tool, NOT a per-commit gate**, and is deliberately
**not** included by `test/runtests.jl` (which `include`s its files explicitly —
a `benchmark/` dir is never globbed in). Wall-clock timing is too flaky to be an
always-on gate.

## Why the dependency lives here, not in the package

**BenchmarkTools.jl is a TOOL dependency, not a SUITE dependency** — exactly as
Coverage.jl is for the coverage tool and the bespoke harness is for the mutation
gate. It is declared in `benchmark/Project.toml` (which `dev`s the package), and
is **never** added to PadeTaylor's `[deps]` / `[targets].test`. Adding it there
would pull a measurement library into every `Pkg.test` run that does not need it
and blur the package's runtime contract.

Unlike the coverage tool (which installs Coverage.jl into a runtime
`mktempdir()` because it only reads `.cov` byte files and needs no stable env), a
perf **baseline** is only comparable against the *same* BenchmarkTools version
that produced it — so this tool uses a **committed, pinned `benchmark/Project.toml`**
(the idiomatic Julia / PkgBenchmark / AirspeedVelocity layout) rather than a
freshly-resolved temp env. See that file's header for the full rationale.

## Run it (single Julia process — Rule 7)

One-time env setup (the tool also self-heals if PadeTaylor is unresolved):

```bash
julia --project=benchmark -e 'import Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
```

Then, on an **idle box with threads pinned** (see *idle-box discipline* below):

```bash
# Compare current timings to the committed baseline (the gate):
JULIA_NUM_THREADS=1 julia --project=benchmark benchmark/run_perf_gate.jl

# Regenerate the baseline (EXPLICIT, reviewed-commit path — never silent):
PADETAYLOR_REGEN_PERF_BASELINE=1 JULIA_NUM_THREADS=1 \
    julia --project=benchmark benchmark/run_perf_gate.jl
```

Exit code is `0` iff every hot path is within the slack band, `1` on any flagged
regression. BLAS is pinned to one thread inside the tool; Julia thread count is
set by the `JULIA_NUM_THREADS=1` launcher and **reported** in the header.

## The benchmarks (6 hot paths)

| id | what it times | signature |
|----|---------------|-----------|
| `pade_step` | one full Padé-Taylor step (jet → rescale → SVD-Padé → eval) | `pade_step_with_pade!` |
| `robust_pade_svd_1515` | the `(15,15)` `:svd` Padé at production order — the dominant per-step cost | `robust_pade(…; method=:svd)` |
| `taylor_coefficients_2nd_o30` | the order-30 second-order Taylor jet | `taylor_coefficients_2nd` |
| `evaluate_pade` | Horner numerator/denominator evaluation of a `(15,15)` Padé | `_evaluate_pade` |
| `solve_pade_wp_multiseg` | a multi-segment IVP solve (Weierstrass-℘, 5 segments) | `solve_pade` |
| `path_network_short_walk` | a short complex-plane network walk to one target | `path_network_solve` |

All inputs are deterministic (the FW 2011 equianharmonic Weierstrass-℘ IC; no RNG
anywhere). Each `@benchmark` rebuilds fresh state per sample via `setup=` so we
time the work, not cached state.

## Read the report

```
PERFORMANCE-REGRESSION GATE  (minR = cur_min/base_min [GATE]; medR for context)
THREADS: julia=1  BLAS=1 (lbt)  [pinned]
SLACK: 50%  (GATE on the MIN ratio; flag when current_min > baseline_min·1.5)

  benchmark                      base_min   cur_min       minR    medR
  pade_step                      27.31µs    30.66µs       1.123   1.116
  ...
OK — every hot path within the 50% slack band.
```

- **The gate compares the BenchmarkTools `minimum`** (`minR`), not the median.
  The minimum is the noise-floor estimator and is the most robust to *upward*
  jitter (BenchmarkTools' own guidance) — a transient load event can only make a
  run *slower*, never faster, so the minimum filters exactly the noise we must
  not chase. The `median` (`medR`) is reported for context and is what we store
  and diff in the baseline.
- **All-flag vs one-flag.** If **every** benchmark flags at once, that is the
  signature of a whole-machine transient (load/thermal), **not** a code
  regression — re-run on a truly idle box. If **one** (or a coherent subset that
  share a callee) flags, that is a localized regression — `git diff` the `src`
  files that path touches.

## Calibrating the slack (the bead's open question, answered)

The bead asked: *what slack suits this WSL2 box with no CI?* Measured here by
running the gate back-to-back ~20 times against a warm baseline, threads pinned,
on this box (Julia 1.12.5, linux x86_64 WSL2, OpenBLAS single-thread):

- **Cheap, allocation-light kernels** (`pade_step`, `robust_pade_svd_1515`,
  `taylor_coefficients_2nd_o30`, `evaluate_pade`) hold min-ratio **≤ ~1.13**
  across clean runs — tight and well-behaved.
- **The two allocation-heavy paths** — `solve_pade_wp_multiseg` (~460 allocs) and
  `path_network_short_walk` (~725 allocs) — carry **real GC-driven min-jitter to
  ~1.25–1.30 steady-state** (verified over 6 high-sample reps: min-ratio clustered
  1.21–1.31 with occasional clean dips to ~0.99). To stabilise them the gate runs
  these two with more samples and `gcsample=true`.
- **A whole-machine transient** can push *any* path higher for a single run (we
  observed an isolated `solve_pade` min-ratio of 1.55 during a busy moment, while
  its inner kernel `pade_step` sat at 1.03 — proof it was load, not algorithm).
  The **min** estimator is markedly tighter than the median throughout (e.g.
  `path_network` medR 1.171 vs minR 1.139), which is why the gate is on the min.

**Chosen `SLACK = 0.50` (flag at 1.50× the baseline minimum).** It sits
comfortably above the worst steady-state min-jitter (~1.30, the allocation-heavy
tail) so noise never false-alarms, yet trips a genuine ≥1.5× regression and any
order-of-magnitude blow-up (the bite-proof below clears it by 19–33×). A single
transient run that flags one allocation-heavy path while its inner kernel stays
green is load, not a regression — re-run. If a future box shows wider jitter,
widen `SLACK` (in `run_perf_gate.jl`) and re-document the measured jitter here —
never narrow it to chase a noisy "regression" that is really load.

## Idle-box discipline (load-bearing)

A baseline — and any comparison against it — is only meaningful on an **idle
machine with threads pinned**:

- Launch with `JULIA_NUM_THREADS=1`; the tool pins BLAS to one thread
  (`LinearAlgebra.BLAS.set_num_threads(1)`) up front, because a multi-threaded
  OpenBLAS makes the SVD timing — `robust_pade :svd`, the dominant per-step cost
  — non-reproducible run-to-run.
- The tool **reports** the thread/BLAS fingerprint in the header and **warns**
  (and records it in the baseline `[meta]`) if threads are not pinned.
- Capture baselines with nothing else running. WSL2 in particular shares the
  Windows scheduler; close other heavy processes first.

## Regeneration / approval

Baseline updates are **explicit, reviewed commits — never silent** (mirrors the
accuracy-ledger's `PADETAYLOR_REGEN_ACCURACY_BASELINE` workflow, `krgy.8`, which
ties to the snapshot/approval bead `krgy.9`):

```bash
PADETAYLOR_REGEN_PERF_BASELINE=1 JULIA_NUM_THREADS=1 \
    julia --project=benchmark benchmark/run_perf_gate.jl
```

rewrites `benchmark/baseline.toml` (preserving its literate header) from the
current measurement, then a human reviews the `git diff` (every changed ns is
visible) and commits it as the new approved baseline. The tool **never** writes
the baseline on a normal run. Regenerate after a deliberate, understood
speed-up, or when porting to a new box (expect last-digit movement, not
order-of-magnitude change — an order-of-magnitude move is a real regression).

## Bite-proof (Rule 4 analog) — VERIFIED 2026-06-09

The gate must FLAG a genuine slowdown. Verified end-to-end: an `O(order·5000)`
busy-loop was injected into `taylor_coefficients_2nd`'s bootstrap loop in
`src/Coefficients.jl`. RESULT — the gate flagged exactly the right paths:

- `taylor_coefficients_2nd_o30` minR **33.5×** (the directly-mutated path);
- `pade_step` 28.0×, `solve_pade_wp_multiseg` 29.7×, `path_network_short_walk`
  19.5× — all three call the jet, so all three flagged;
- `robust_pade_svd_1515` (1.04×) and `evaluate_pade` (1.30×) stayed **GREEN** —
  neither calls the mutated function, correctly localizing the culprit to the
  Taylor-jet path.

All ratios vastly exceeded the 1.40× slack. The injection was then **restored**;
`git diff src/` re-confirmed byte-clean. (A subtler real regression would land
just above 1.40× rather than at 33×; the point is the gate trips, and the GREEN
paths point at the offending callee.)

## How it works (implementation)

`run_perf_gate.jl` (≈178 code LOC, ≤200):

1. self-heals the tool env (`Pkg.develop` PadeTaylor + `instantiate`) if needed;
2. pins BLAS to one thread, builds deterministic fixtures, and runs each
   `@benchmark` with modest sample counts (a few-minute total);
3. on a normal run, parses `baseline.toml` and reports the per-benchmark
   min-ratio (gate) + median-ratio (context), exit `1` on any flag;
4. on `PADETAYLOR_REGEN_PERF_BASELINE=1`, rewrites the baseline (header
   preserved) and exits — the explicit, reviewed-commit regeneration path.

The tool never touches `src/` (the bite-proof injection is a manual, restored
step). `baseline.toml` is committed and diff-friendly (sorted TOML).
