# Unified path-network algorithm specification

> Synthesised from FW 2011 §3.1 + §3.2.2 + §4.4 + §5.3 + §5.4.1,
> FW 2014 (no algorithm change; new operational layer), FW 2015 (no
> algorithm change; new initialisation), RF 2014 §3 (codification of
> FW 2011 + pole-counting layer), FFW 2017 §2.1–§2.4 (five new
> architectural moves). Every numbered citation is to the markdown
> extractions at `references/markdown/<paper>/<file>.md:<lines>`.
> Reconciliations between conflicting paper conventions are flagged
> at the point of choice.

## §0. Scope and tier alignment

This spec covers Tiers 2–5 of `docs/figure_catalogue.md`. Tier 2 is
the FW 2011 baseline path-network. Tier 3 adds the BVP solver
(separate document; see `references/bvp_recipe.md` once acquired).
Tiers 4–5 are FFW 2017 extensions (coord transform, sheet tracking).

The spec is **layered**: Tier 2 stands alone and is the v1
deliverable. Tier 3 hooks in via a dispatcher that does not modify
Tier-2 internals. Tiers 4–5 wrap Tier-3 at the coordinate layer.
This is non-negotiable — adding sheet-tracking *into* the inner Padé
step would couple every layer below it.

## §1. Stage 1 — path-tree construction

**Goal:** for a target set `T = {z_t1, z_t2, …}` in the complex
plane, build a tree of integration paths rooted at the IC origin `z_0`,
such that each target `z_ti` is reached by some path whose every
segment is a single Padé-Taylor step bridging at most one pole.

**FW 2011 prescription** (`FW2011...md:155–166`):

  1. Initialise `visited = { z_0 ↦ (u_0, u'_0, Padé₀) }`.
  2. Choose a target order — FW use a **random shuffle** of the
     coarse-grid nodes. This avoids systematic bias toward one
     direction of the plane.
  3. For each target `z_t` in turn:
     a. Pick the **nearest already-visited point** `z_v ∈ visited`
        (Euclidean distance in ℂ) as the starting node.
     b. Step from `z_v` toward `z_t` repeatedly until the step
        endpoint is within `h` of `z_t`. Each step uses the
        5-direction wedge of §3 below.
     c. Record every intermediate landed point in `visited` with
        its full local Padé approximant (numerator + denominator
        coefficients). Paths do not cross by construction.

**Termination per target** (`FW2011...md:163`): halt when
`|z_current − z_t| ≤ h`. The final step is the Stage-2 fine-grid
extrapolation; it is not added to `visited`.

**INFERRED ambiguity:** "nearest already-visited point" is defined
by Euclidean distance in ℂ, with ties broken by lexicographic order
on (Re, Im). FW 2011 does not specify the tiebreak; we adopt
lexicographic because it gives reproducible behaviour for tests.

## §2. Stage 2 — fine-grid extrapolation

**Goal:** dense evaluation of `(u, u')` on a fine equispaced
lattice covering the same domain, without further Padé construction.

**FW 2011 prescription** (`FW2011...md:166, 397`):

  - For each visited Stage-1 node `z_v ∈ visited`, identify the
    fine-grid lattice nodes `{z_f}` within distance `h` of `z_v`.
  - For each such `z_f`: evaluate the **stored** Padé approximant
    at `z_v` at the point `t = (z_f - z_v) / h` (`z_f` is inside
    the disc of convergence of the local rational approximant).
  - The output is `(u(z_f), u'(z_f))` plus the per-lattice-node
    diagnostic `|u(z_f)|`.

**Performance contract** (`FW2011...md:166`): Stage 2 is
vectorisable — N fine-grid nodes evaluate against M stored Padé
approximants at total cost O(N) per stored approximant. RF 2014
(`RF2014...md:161–165`) makes this the *explicit* second step of a
four-step algorithm.

## §3. The 5-direction wedge

**FW 2011** (`FW2011...md:158–159`) — five candidate directions
relative to the line from `z_v` toward `z_t`:

  - **angles: 0°, ±22.5°, ±45°**. Wedge total: 90°.
  - The choice is "five Padé evaluations at negligible extra cost"
    (same Taylor jet, five direction-conditioned rescalings).

**RF 2014** (`RF2014...md:148–153`) restates the algorithm with
**±15°, ±30°** (still five directions, wedge total 60°). The
difference is a downstream choice; the codification is the same.

**Spec reconciliation:** parametrise. Public API exposes
`wedge_angles::Vector{Float64} = [-π/4, -π/8, 0, π/8, π/4]` per
FW 2011 default. RF 2014 users override to `[-π/6, -π/12, 0, π/12, π/6]`.

**FFW 2017** (`FFW2017...md:163–178`) constrains the wedge further
for PVI circumambulation: the wedge is **rotated by the cumulative
winding angle** around branch points, so the algorithm steers
deliberately clockwise or counterclockwise. See §10 below.

## §4. Step-selection rule

**FW 2011 default — min-|u|** (`FW2011...md:159, 354`):

  - For each of the 5 candidate endpoints, evaluate `|u(z_candidate)|`.
  - Select the candidate with minimum `|u|`. Rationale: minimising
    `|u|` heuristically maximises the radius of convergence of the
    local Taylor jet (because poles are where `|u| → ∞`).

**FW 2011 §5.4.1 variant — steepest descent** (`FW2011...md:362–368`):

  - Compute analytically the direction of steepest decrease of `|u|`:
    `θ = arg(−u(z_v) / u'(z_v))`.
  - If `θ` falls inside the ±45° wedge, use that single direction;
    else clip to the nearest wedge edge.
  - Saves 4 of the 5 candidate evaluations. FW report this as the
    "preferred choice when implemented in Fortran or C."

**Spec choice:** implement both. Default `step_selection = :min_u`.
`:steepest_descent` is the perf-tuned alternative. Tests for both.

## §5. Step-size policy

**FW 2011 default — fixed h** (`FW2011...md:164, 326`):

  - `h = 0.5`, `order = 30` (i.e., (15, 15) Padé). This is "a
    generally favourable combination" per FW 2011 §5.2.
  - Coarser `h = 0.3` illustrated in Fig 3.2 for visual clarity.

**FFW 2017 adaptive** (`FFW2017...md:81–92`):

  - Error estimate `T(h) = |ε_{n+1} · h^{n+1}|` from the leading
    truncation term of the local Taylor jet (we have `ε_{n+1}` as
    `c_{n+1}` from the jet at order `n`).
  - If `T(h) > Tol`, rescale `h ← q · h` where `q = (k · Tol / T(h))^{1/(n+1)}`.
  - FFW use `k = 1e-3`. The formula derives from two crude
    approximations (FW 2017 §2.4) and is *additional* to the
    direction selection, not a replacement.

**Spec choice:** parametrise `step_size_policy ∈ {:fixed, :adaptive_ffw}`.
Default `:fixed` per FW 2011 (Tier-2 scope). `:adaptive_ffw` opens
Tier-4 (PIII/PV) and is required for non-uniform pole densities.

**Composition with our existing `step_jorba_zou` / `step_pade_root`**
(Phase 4): these are Jorba-Zou §3.3.1 fixed-order step controllers,
suitable for smooth-ODE accuracy control but NOT for FW's Padé-
bridge use-case where the bigger step is exactly what bridges the
pole (per `docs/worklog/004-phase-6-pivot.md` attempt C). They are
unrelated to FFW's adaptive policy and are not used by the
path-network. They remain available for the legacy `solve_pade` driver.

## §6. Grid layout and Stage-1 node density

**FW 2011** (`FW2011...md:155–156`): **uniform rectangular grids**.
Coarse 40×40 over a 20×20 complex-plane window for Fig 3.2; fine
161×161 over the same window for Fig 3.1. No formula for scaling
with domain size or pole density — must be chosen empirically.

**FFW 2017** (`FFW2017...md:67–71`): **non-uniform Stage-1 nodes**.
For P̃III / P̃V in ζ-coordinates, pole density grows exponentially
along `Re ζ`, so a uniform grid wastes computation in the sparse
region and under-resolves the dense region. FFW use the Fornberg-
Flyer 2015 node-placement algorithm with spatial separation function
`R(ζ)` decreasing linearly with `Re ζ`.

**Spec choice:** Tier 2 ships **uniform grid only**. Non-uniform
is a Tier-4 extension and a separate bead. The path-network module
takes a `grid::Vector{Complex{T}}` parameter so callers can supply
any layout; uniform is just the default constructor.

## §7. Failure handling

| failure mode | source layer | handling |
|---|---|---|
| Singular `C̃` Toeplitz in `robust_pade` | `RobustPade` | Throw `Suggestion` (already shipped Phase 2). |
| `Q(t) = 0` at `t = 1` (pole at step endpoint) | `PadeStepper` | Detect; vault by multiplying `h` by 1.2 (FW heuristic) OR shrink `h` toward Jorba-Zou bound. Path-network catches and routes to a different candidate direction; if all 5 fail, retry with smaller `h`. (`docs/worklog/004-phase-6-pivot.md` attempts B and C are the failure-mode tests.) |
| Step-control non-convergence (FFW adaptive) | `PathNetwork` step loop | Cap rescaling at `q ∈ [1/10, 10]`; if cap binds for 3 consecutive steps, fall back to fixed-h with diagnostic message. |
| BVP edge crossing (T3) | `Dispatcher` | Edge detector flags lattice cells; path-network avoids them; smooth band solved by `BVPSmooth` module. |
| Sheet boundary (T5) | `SheetTracker` | Path tags its sheet index `s ∈ ℤ`; branch-cut crossing increments / decrements `s`; cross-sheet evaluations are *forbidden* without explicit `cross_branch=true`. |

## §8. Composition with IVP↔BVP dispatcher (Tier 3)

**FW 2011 §3.2.2 + §4.4** (`FW2011...md:204–208, 249–261`):

  - **Edge detector**: compute `|Δu|/h²` via the 5-point Laplacian
    stencil `[1, 1, −4, 1, 1]` on the equispaced lattice values
    `u(z_f)` from Stage 2.
  - Lattice cells with `|Δu|/h² > 0.001` are classified as
    pole-field (IVP-stable); cells below as smooth (IVP-unstable,
    BVP-required).
  - The threshold `0.001` is empirical for PI at `h = 0.5`. User-
    tunable.

**Architecture choice:** the dispatcher is a *post-processing*
layer over the Stage-2 fine grid, not an *online* component of
Stage 1 path construction. Path-network builds the pole-field
solution first; edge detector partitions the lattice; BVP fills
the smooth regions; the union is the final result. This avoids
coupling the path-network to a BVP that may not exist (Tier-2-only
deployments skip the dispatcher entirely).

## §9. Coordinate transformation layer (Tier 4)

**FFW 2017** (`FFW2017...md:40–48`):

  - For PIII: substitute `z = e^(ζ/2)` and `u(z) = e^(−ζ/2) w(ζ)`.
    The result `P̃III` is solved in the ζ-plane by the unmodified
    path-network. Multi-sheet z-plane output is recovered by
    partitioning the ζ-plane into 4π-wide strips and applying the
    inverse transform per strip.
  - For PV: substitute `z = e^ζ` (no factor 1/2). Strip width 2π.
  - For PIV: no transformation needed (single-valued in z).
  - For PVI: no transformation works for all sheets (two fixed
    singularities at z=0, z=1); see §10.

**Architecture choice:** the transformation is a wrapper that
   - converts the ODE coefficients to ζ-form,
   - calls the path-network in ζ,
   - splits the ζ-output by 4π-strip (or 2π for PV),
   - applies `z = e^(ζ/2)` (or `e^ζ`) to each strip,
   - assembles the multi-sheet z-plane result.

The path-network itself is unaware of ζ. This is the layering
discipline that keeps the inner solver clean.

## §10. Sheet tracking and branch-cut routing (Tier 5)

**FFW 2017** (`FFW2017...md:163–189`):

  - PVI has two fixed singularities at z=0 and z=1. No single
    coordinate transform can avoid both.
  - The path-network is modified to **tag each path with a sheet
    index**, where the sheet index is the winding number around
    each branch point accumulated along the path so far.
  - The path-selector deliberately **circumambulates** branch
    points: when the natural step direction would cross a branch
    cut, the algorithm picks the wedge direction that goes around
    instead, incrementing the sheet index by ±1 depending on
    direction.

**Architecture choice:** `SheetTracker` is a layer above
`PathNetwork`. It augments the step-selection rule with a "no
branch-cut crossing without explicit consent" constraint and
maintains the sheet-index bookkeeping. Cross-sheet evaluations
are explicit operations, not implicit.

## §11. Initialisation strategies

| problem class | strategy | source |
|---|---|---|
| PI tritronquée from origin | FW 2011 closed-form ICs from §4.1 | `FW2011...md:215–222` |
| FW Table 5.1 (℘-function test) | Closed-form `WeierstrassP` ICs | `FW2011...md:316` |
| Real PII tronquée | Asymptotic shooting from `X = 8` with 13+9 coefficients | `FW2014...md:354–412` |
| Imaginary PII tronquée | Asymptotic series at `|x| ≈ 7–10` (eqs. 10, 12, 13) | `FW2015...md:205` |
| PIII/PV asymptotic | Asymptotic series at large \|z\|, then forward into ζ-domain | `FFW2017...md:75–80` |

**Spec choice:** `init_strategy ∈ {:closed_form, :asymptotic_shoot, :user_supplied}`.
The path-network does not synthesise ICs; it takes them as input. A
separate `InitialConditions` module produces them. Tier-2 v1 only
requires `:user_supplied` (the caller hands the IC pair); the
asymptotic-shoot helpers are Tier-3 enhancements.

## §12. Open spec gaps

The following are unspecified in the source papers and require
local decisions, all documented inline:

1. **Stage-1 distance metric and tiebreak** — FW gives "nearest"; we
   adopt Euclidean + lexicographic. (§1)
2. **Wedge angle convention** — FW2011 ±22.5°/45° vs RF2014 ±15°/30°.
   Parametrised. (§3)
3. **Step-rescaling cap range** — FFW gives the formula, not the
   cap. We adopt `q ∈ [1/10, 10]` with 3-step-fallback. (§7)
4. **Edge-detector threshold portability** — 0.001 is empirical for
   PI on h=0.5 lattice. Other equations / grids need recalibration.
   User-tunable. (§8)
5. **Sheet-index initial value** — FFW's `s = 0` convention for the
   "principal sheet" is implicit, not stated. We adopt it. (§10)
6. **Cross-path agreement tolerance** — FW does not give a
   pairwise-merge criterion between independently traced paths
   reaching the same lattice node. INFERRED: no online merging;
   accept the first path's value; defer cross-validation to a
   post-pass diagnostic that flags discrepancies > 10·rtol.
7. **Maximum step count per target path** — FW does not give one;
   we adopt `max_steps_per_target = 1000` with explicit error if
   exceeded. (§1)
8. **Newton convergence tolerance for BVP** — FW only gives the
   *post-solve* derivative-match tolerance (`1e-7 / 1e-8`); the
   Newton-residual stopping rule is unspecified. Adopt `‖R‖_∞ ≤
   100 · eps(T)` per `T <: AbstractFloat`. (BVP spec, §4.)

## §13. Module decomposition preview (subject to ADR-0004 design)

```
src/
  PathNetwork.jl       — Tier 2 (§§1–7); ≤200 LOC
  EdgeDetector.jl      — Tier 2+ (§8 detector half); ≤50 LOC
  BVP.jl               — Tier 3 (separate spec doc); ≤200 LOC
  Dispatcher.jl        — Tier 3 (§8 dispatcher half); ≤100 LOC
  CoordTransforms.jl   — Tier 4 (§9); ≤100 LOC
  SheetTracker.jl      — Tier 5 (§10); ≤150 LOC
  InitialConditions.jl — Tier 2+ (§11); ≤100 LOC
```

Module dependencies fan upward: PathNetwork has no dependency on
anything Tier ≥ 3. The Dispatcher composes the Tier-2 PathNetwork
output with the BVP output. SheetTracker wraps PathNetwork with no
modification to its internals.

## §14. Pointers

  - `docs/figure_catalogue.md` — concrete acceptance criteria per
    figure, tier-aligned with this spec.
  - `references/markdown/<paper>/<file>.md` — the source ground truth.
  - `references/bvp_recipe.md` — Tier-3 BVP spec (acquired separately).
  - `docs/worklog/004-phase-6-pivot.md` — failure analysis of
    fixed-h + pole-vault + Jorba-Zou for path-less stepping; the
    motivation for the path-network.
  - Beads: `padetaylor-8cr` (umbrella P0), `padetaylor-1jf`
    (PathNetwork impl), `padetaylor-c2p` (EdgeDetector impl),
    `padetaylor-rgp` (figure catalogue), plus the BVP-module bead.
