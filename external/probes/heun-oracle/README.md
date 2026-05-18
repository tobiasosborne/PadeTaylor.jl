# heun-oracle — Mathematica HeunG/HeunC ground-truth capture

## What this probe does

`capture.wl` uses Mathematica's `HeunG` and `HeunC` built-in functions to
produce 30-significant-digit reference values at a curated grid of parameter
regimes.  The output (`oracles.txt`) is consumed by downstream Julia tests in
PadeTaylor.jl that validate the project's own Heun evaluators.

## Prerequisites

- **Mathematica 14.x** with `wolframscript` on PATH.  `HeunG` and `HeunC` were
  added in Mathematica 12.0 (2019); any v12+ installation works.
- Verified locally: WolframScript 1.13.0 / Mathematica 14.3.0 (July 2025).
- To confirm availability: `wolframscript -code 'Information[HeunG]'`

## How to run

```bash
cd external/probes/heun-oracle
wolframscript -file capture.wl > oracles.txt
```

Progress lines (prefixed `[OK]`) go to stdout interleaved with data; both land
in `oracles.txt`.  To separate them:

```bash
wolframscript -file capture.wl 1>oracles.txt 2>progress.log
```

Expected wall time: **under 2 minutes** for all 42 cases at WorkingPrecision 50.

## Parameter conventions

| Function | Signature | Normalisation |
|----------|-----------|---------------|
| `HeunG`  | `HeunG[a, q, alpha, beta, gamma, delta, z]` | value = 1 at z=0 |
| `HeunC`  | `HeunC[q, alpha, gamma, delta, epsilon, z]` | value = 1 at z=0 |

Precision: all parameters supplied as `SetPrecision[x, 50]` (50-digit
arbitrary precision); results rounded to 30 significant digits.

## Test-case regimes

| Regime | Description | Cases |
|--------|-------------|-------|
| A | Easy reference: `HeunG[2,1,1,2,3,4,z]` and `HeunC[1,1,2,3,-1,z]` at z ∈ {0.1,0.5,0.9,1.5,1.9,2.5,3.0} | 14 |
| B | DLMF Ch31 special-case identity checks (degenerate & epsilon=0 params) | 9 |
| C | Hard: a=1.001 near-coincident singularities, q near apparent-singularity, large \|z\|, large \|q\| | 11 |
| D | HeunC for Teukolsky-Schwarzschild QNM (PLACEHOLDER pending teukolsky_heun_mapping.md) | 3 |
| E | Complex z: off-real upper and lower half-plane branch probes | 5 |
| **Total** | | **42** |

## Regime D note

`docs/teukolsky_heun_mapping.md` was not present when this script was written.
The D-regime cases use Fiziev (2006) approximate parameters for the l=2, n=0,
s=0 QNM mode.  Once the mapping document is produced, replace the `omegaM`,
`qT`, `alT`, `gaT`, `deT`, `epT` variables in the D section of `capture.wl`
with the exact values from that document and re-run.

## oracles.txt format

```
# header lines starting with #
HeunG|a=2|q=1|alpha=1|beta=2|gamma=3|delta=4|z=0.5|value=1.09086...-0.0i
```

Each data line: `FuncName|param=val|...|z=val|value=re±imi` where `re` and
`im` are given to 30 significant figures.

## Known fragility

- Mathematica's HeunG path-follower can return `Indeterminate` or lose
  precision near the branch cuts (real z ∈ [1, a] and [a, ∞)).  The hard
  cases in regime C deliberately probe these regimes.  If a value returns
  `Indeterminate`, the capture script will write `Indeterminate` in the value
  field; downstream tests should skip those entries.
- The `Information[HeunG]` call triggers a segfault on WolframScript 1.13 on
  Linux (known upstream issue); `HeunG` evaluation itself is unaffected.
