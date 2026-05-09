# PRD: PadeTaylor — General-Purpose Taylor–Padé IVP Solver

> **Status: SKETCH (v0.1).** First draft, written from memory and a single web pass.
> Do **not** start coding from this. **Stage 0 (deep research) is the first step**;
> everything below is provisional and will be replaced by `RESEARCH.md` + `DESIGN.md`.

---

## Goal

A general-purpose Taylor–Padé IVP solver for analytic ODEs `y' = f(z, y)`,
`z ∈ ℂ`, supporting:

- Movable poles — pole-immune stepping in the Fornberg–Weideman tradition.
- Arbitrary precision and Float64.
- Complex-analytic problems beyond Painlevé (general `f` built from elementary functions).

Two reference implementations, mutually validating:

1. **Julia**: `PadeTaylor.jl` — arb-prec via `Arblib`, Float64 generic.
2. **Scientist Workbench**: TS-native package under `@cas/`, arb-prec via `@cas/cas-core` primitives.

Goal is for these to be **the reference** — well-specified, well-tested,
the thing others port from. Without 1:1 cross-validation between the two,
neither qualifies.

---

## Stage 0 — Deep research (do this BEFORE planning)

Output: `RESEARCH.md` summarising findings before any design or code.

### Primary literature (read carefully, store PDFs in `references/`)

- **Fornberg & Weideman**, *A numerical methodology for the Painlevé equations*,
  J. Comput. Phys. 230 (2011) 5957–5973. The foundational paper. Algorithmic detail in §3–4.
- **Fornberg & Weideman**, *A computational exploration of the second Painlevé equation*,
  FoCM 14 (2014) 985–1016. Refinements, PII survey.
- **Fasondini, Fornberg, Weideman**, J. Comput. Phys. 344 (2017) 36–50.
  Extension to multi-valued transcendents on Riemann surfaces.
- **Reeger & Fornberg** — PIV survey paper(s).
- **Willers (1974)** — original IVP-with-poles algorithm; FW generalises this.

### Padé robustness (essential)

- **Gonnet, Güttel, Trefethen**, *Robust Padé Approximation via SVD*, SIAM Review 55 (2013).
  This is **the** Padé routine. Read fully before implementing anything.
- Chebfun's `padeapprox` source — reference impl of GGT 2013.
- Hermite–Padé literature (Baker & Graves-Morris) — for v2 / branch-point detection.

### Adjacent Taylor-IVP packages (study designs, do not reinvent)

- **TIDES** (Barrio et al.) — mature high-order Taylor.
- **TaylorIntegration.jl** (Pérez-Hernández, Benet) — Julia; real-axis only, no Padé.
- **TaylorSeries.jl** (Benet, Sanders) — coefficient layer; assess generic-type support
  for `Arb` element type, especially at order 60+.
- **ATOMFT** (Chang & Corliss) — historical reference.
- **Jorba & Zou (2005)** — Taylor IVP package paper, design lessons.

### Validated / arb-prec ODE (adjacent, not duplicates)

- **Mezzarobba** — D-finite / linear ODEs in Arb (`ore_algebra`). Different scope, useful for verified-arithmetic patterns.
- **COSY Infinity** / Taylor models (Makino & Berz) — for verified-bound thinking in v2.

### Source code hunt

- Paper appendices — FW 2011 may include code excerpts.
- Weideman's Stellenbosch page (`appliedmaths.sun.ac.za/~weideman/`) — PFS code does not appear to be posted; confirm.
- arXiv ancillary files for FW papers.
- **Email** Weideman, Fornberg, Fasondini, Reeger. They have working MATLAB and have shared with collaborators historically.

### Ecosystem inventory

- `@cas/cas-core`: does a Taylor / formal power series class exist? Does a bigfloat SVD exist? What's missing?
- Julia: confirm `Arblib.jl` exposes SVD on `ArbMatrix`; confirm `TaylorSeries.jl` is generic enough to use `Arb` element type at order 60+ without breakage.

### Stage 0 deliverable

`RESEARCH.md` containing:

- Confirmed algorithmic spec (Taylor → Padé → step control → optional path network).
- All design choices found in literature with explicit tradeoffs.
- What pieces exist in the ecosystem; what we must build.
- Updated open questions (some below will resolve; new ones will appear).
- Recommendation: proceed to Stage 1, or further research needed.

---

## Stage 1 — Design decisions to lock (after research)

Output: `DESIGN.md` resolving at least:

1. **Padé form**: scalar per component vs Hermite–Padé with shared denominator. (v1 likely scalar; document the v2 path.)
2. **Step control**: coefficient-decay heuristic, Padé-denominator-root distance, or both switchable.
3. **Order strategy**: fixed (probably 30–60) vs adaptive.
4. **Path strategy in v1**: real-axis only with primitive `step_complex` API, vs 2D network from day one. (v1 likely real-axis only.)
5. **API surface**: standalone clean core + thin SciML `CommonSolve` wrapper for Julia.
6. **Verified arithmetic** via Arb balls — v1, v2, or never.
7. **Module boundaries** aligned with the ~200 LOC/file budget.

---

## Algorithmic core (sketch — fill in after research)

```
1. TaylorStep(f, z_n, y_n, N)   → coefficients a_0, ..., a_N
                                   (Taylor-mode AD on operator-overloaded polynomial type)
2. RobustPade(coeffs, [L/M])    → (P, Q)  via Gonnet–Güttel–Trefethen SVD
3. StepControl(coeffs, P, Q, τ) → dz_next
4. Step:  y_{n+1} = (P/Q)(z_n + dz_next)
5. (v2) PathNetwork: graph of steps with consistency checks at meeting points.
```

---

## Implementations

### Julia (`PadeTaylor.jl`)

- Deps: `TaylorSeries.jl`, `Arblib.jl`, `LinearAlgebra` (or generic SVD for `Arb`), `Polynomials.jl` (roots of `Q`).
- Modules (each ≤ ~200 LOC):
  - `RobustPade` — port of GGT 2013.
  - `PadeStepper` — coefficient generation + Padé + step control.
  - `Problems` — IVP problem and solution types.
  - `PathNetwork` (v2).
- Generic over `T <: Number`; `Float64` and `Arb` both first-class.

### Scientist Workbench (TS, monorepo)

- Package: `@cas/pade-taylor` (name TBD).
- **Mirror Julia module structure 1:1.** This is non-negotiable; mirroring is what makes them mutually a reference impl.
- Depends on `@cas/cas-core` for arb-prec primitives, polynomial type, bigfloat SVD.
- **Bigfloat SVD is the heaviest unknown.** If `@cas/cas-core` lacks it: scope as a sub-task and port a known-good algorithm (Demmel–Kahan, or one-sided Jacobi for arb-prec). **Do not hand-roll a fragile version.**

### Cross-validation

Shared test corpus; both impls run it; diff to N digits. Without 1:1 mirroring and shared corpus, neither is a reference impl.

---

## Test corpus

- PI tritronquée (Boutroux IC) — pole locations vs FW 2014 figures.
- PII Hastings–McLeod and Ablowitz–Segur — vs FW 2014 / Reeger–Fornberg.
- Riccati with closed-form meromorphic solution — analytic pole locations as ground truth.
- A Schwarzian / non-Painlevé analytic ODE — proves "general purpose."
- A bog-standard non-stiff problem on the real axis (Lorenz segment, Van der Pol away from stiff regime) — sanity check vs `DifferentialEquations.jl`.
- DLMF / Mathematica `PainleveT` tabulated values where available.

---

## Non-goals (be explicit in docs)

- Stiff systems. PFS is explicit; A-stability is not the concern here. Direct users to Rosenbrock / IMEX from `DifferentialEquations.jl`.
- DAEs.
- Essential singularities (genuinely defeats the method).
- Real-time / embedded.
- Drop-in replacement for SciML's general ODE infrastructure. **USP is correctness through pole fields and arb-prec, not breadth.**

---

## Process / philosophy constraints (per project CLAUDE.md)

- **Beads** for issue tracking. No GitHub Issues.
- **Fail-fast / fail-loud.** No silent fallbacks; numerical breakdown errors explicitly.
- **Ground truth via local PDFs.** FW papers, GGT 2013, etc. in `references/`; string-match against them rather than paraphrasing from memory.
- **No parallel agents** (Julia precompilation cache conflicts).
- **Rigorous reviewer agent** after core changes (Padé routine, stepper, AD layer).
- **~200 LOC per file.** If a module exceeds this, split.

---

## Open questions (provisional)

1. Hermite–Padé in v1 or v2?
2. Step control default: coefficient-decay or Padé-root distance?
3. Verified enclosures via Arb balls — v1, v2, or never?
4. Path-network consistency metric — value, residue, or both, to what tolerance?
5. Naming — `PadeTaylor.jl` sensible, or `MeromorphicIVP.jl` / `Pole field solver.jl` / something else?
6. Licence (MIT/BSD for ecosystem fit, GPL never).

---

## Definition of done (v1)

- [ ] `RESEARCH.md` written and reviewed.
- [ ] `DESIGN.md` written and reviewed.
- [ ] Both impls pass the shared test corpus to declared precision.
- [ ] Both impls reproduce FW 2014 PII pole-field figures qualitatively.
- [ ] Cross-impl diff ≤ 1 ULP at Float64, ≤ N digits at arb-prec for chosen N.
- [ ] Per-file LOC budget respected.
- [ ] Documentation: paper-style spec doc, API reference, minimal tutorial.

---

## Immediate next action for the follow-up agent

**Do not skip Stage 0.** Begin by producing `RESEARCH.md` per the checklist above.
Confirm or refute every provisional claim in this sketch against the literature
and the ecosystem. Return with findings and a recommendation before any code
or detailed design.
