# Froissart-consumer calibration — FINDINGS

**Probe:** `froissart_residue.jl` (run from repo root, `--project=.`).
**Bead:** the Froissart-consumer fix (child of `padetaylor-flnr`).
**Purpose:** the ADR-0028 dispatch wiring (Phase 2a) is gated because cell B
(value-optimal, selected by the ODE defect) plants Froissart doublets that
corrupt the vector pole-root consumers. The user chose *investigate first* — this
probe calibrates the residue threshold separating genuine shared-Q poles from
cell-B doublets, on the real ℘ FW Table 5.1 walk to z=31 (the VPO.2 setup).

## Result — an ~8-order clean separation

At every one of the 342 visited nodes, the cell-A and cell-B shared-Q approximants
were rebuilt and their nearest denominator root measured (z-distance `h·|t*|`,
genuineness residue `max_c |P_c(t*)/Q'(t*)|`):

| residue at nearest root | min | median | max |
|---|---|---|---|
| cell A — **genuine ℘ poles** | **2.1e0** | 1.1e4 | 2.3e11 |
| cell B near-node poles — **doublets** | 1.6e-15 | 3.0e-13 | **5.5e-8** |

- **The defect selects cell B at all 342 nodes** (the walk threads the clear
  channel between poles, so every step is a regular stretch where B wins value).
- **143 / 342 nodes** have cell B placing a pole at `< 0.3×` cell A's z-distance to
  the genuine ℘ pole (typically z-dist ~1e-3 vs the genuine ~0.7) — the spurious
  near-node doublet that collapses the adaptive step controller.
- Every such near-node cell-B root has residue ≤ **5.5e-8**; every genuine ℘ pole
  has residue ≥ **2.1**. The gap is ~8 orders.

## Calibrated threshold

A residue filter `keep root iff max_c |P_c(t*)/Q'(t*)| ≥ min_residue` cleanly
separates the two populations for any `min_residue ∈ (5.5e-8, 2.1)`.

- **`min_residue ≈ 1e-4`** is recommended — ≥4 orders of margin above the doublet
  max (5.5e-8) and below the genuine min (2.1).
- **The scalar `PoleField.jl` default `min_residue = 1e-8` is too low here** — a
  doublet reaches 5.5e-8 and would leak through. The vector consumers need the
  calibrated (higher) threshold, not the scalar default.

## Conclusion — the consumer-filter fix is well-founded

Adding this residue Froissart filter to the two vector pole-root consumers —
`VectorWedgeStep._candidate_pole_disc` / `_adaptive_h` (the adaptive step
controller) and `VectorPoleField.extract_poles_shared_q` (the pole extractor) —
will make them skip cell B's near-node doublets while keeping the genuine ℘ poles,
so the dispatch can optimize the step value freely (defect picks B everywhere)
without corrupting the pole field. The filter is the residue criterion the scalar
`PoleField` already uses (`res = N(t*)/D'(t*)`, `:195-196`), lifted to the shared-Q
vector form and re-tuned to the calibrated `min_residue ≈ 1e-4`.

Note: filtering changes only which roots are treated as **poles** (step distance,
extraction) — never the step **value** `P(1)/Q(1)`, so the value-accuracy gains of
the dispatch are untouched.
