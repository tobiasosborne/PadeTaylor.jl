# ADR-0001 — Four-layer architecture

**Status**: Accepted (2026-05-09)
**Context**: Stage 1 design lock for PadeTaylor.jl v1.

## Decision

PadeTaylor.jl decomposes the Taylor–Padé IVP solver into **four
algorithmically independent layers**, mirrored by four tier-2 modules
plus a top-level integrator:

```
1. Coefficients   ── compute c_0, ..., c_n of the local Taylor expansion
                     of y(z + h) at the current state (z, y).
2. RobustPade    ── convert the truncated Taylor series to a
                     diagonal (n/2, n/2) Padé approximant via GGT 2013
                     Algorithm 2 + Chebfun reweighting.
3. StepControl    ── choose step length h based on either coefficient
                     decay (Jorba-Zou 2005) or denominator-root
                     distance (FW 2011 path heuristic).
4. PadeStepper    ── orchestrate one step: compute coeffs, build Padé,
                     pick step length and direction, evaluate Padé at
                     the new point.
```

Each layer is independently testable and replaceable. The interface
between adjacent layers is a typed Julia struct, not a positional
tuple — caller code does not depend on internal layout.

## Why decompose this way

The decomposition mirrors the structure of FW 2011's exposition (§2.1
"Taylor series methods", §2.2 "Padé version", §3.1 "Selection of
integration paths", §5.1 "Choice of order"). Each section maps 1:1
to a module. **This is the natural fault-line of the algorithm**:

- A bug in the coefficient layer surfaces as wrong Taylor coefficients
  before any Padé conversion, which is independently checkable
  against analytic derivatives.
- A bug in the Padé layer surfaces as wrong `(P, Q)` for known-correct
  input coefficients — checkable against `padeapprox.m`'s output.
- A bug in step control surfaces as wrong `h` for known-correct
  coefficients — checkable against Jorba-Zou's eq. 3-8 directly.
- A bug in the orchestrator surfaces as wrong end-of-step state for
  known-correct sub-layer outputs — checkable against FW Table 5.1.

The decomposition also admits clean v2 extensions:
- `step_complex(z, dz)` becomes a primitive that the v2 path-network
  layer composes.
- The BVP solver (FW 2011 §3.2) is a sibling integrator that shares
  the `Coefficients` and `RobustPade` layers but replaces
  `StepControl` and `PadeStepper`.

## Citations

- FW 2011 §2.1 (`references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:79–107`) — Taylor methods,
  efficient coefficient generation.
- FW 2011 §2.2 (lines 109–116, eq. (2.5)) — Diagonal Padé conversion
  rationale.
- FW 2011 §3.1 (lines 151–168) — Step-direction selection.
- GGT 2013 Algorithm 2 (`references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
  GGT2013_robust_pade_via_SVD_SIREV55.md:219–235`) — Robust Padé
  routine.
- Jorba & Zou 2005 §3.2 eq. 3-8 (per RESEARCH.md §3.5) — Step-size
  formula.
- RESEARCH.md §1.1, §3.5, §7.1 Q7 — synthesis underlying this ADR.

## Consequences

- Each module is a single-file `.jl` under `src/`, ≤ ~200 LOC.
- Each module ships its own test file under `test/`.
- Cross-module API is exclusively through public Julia types
  (`PadeApproximant{T}`, `StepperState{T}`, etc.).
- v2's path-network and BVP-coupled solver compose at the
  `PadeStepper` level; the inner three layers are reused unchanged.

## Alternatives considered, rejected

- **Single monolithic stepper** — would couple coefficient generation
  with Padé conversion and step control. Loses independent
  testability; impossible to A/B test step-control strategies without
  rerunning the whole pipeline.
- **Three layers** (fold step control into the stepper) — the FW path
  heuristic (`step_pade_root`) needs the *same* `PadeApproximant`
  the stepper just built; no reason to duplicate the data flow.
  Three layers would force step control to recompute SVD-derived
  data. Rejected.
- **Five layers** (separate "evaluate Padé at z" from the rest) —
  Padé evaluation is a small Horner loop; not worth its own module.
