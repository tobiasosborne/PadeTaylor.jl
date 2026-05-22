# ADR-0025 — Senior-grade re-resolution of the P_I⁽²⁾ tritronquée headline figure

**Status**: Accepted (2026-05-21); scoped work shipped 2026-05-22 —
Phases A–F done, all six v1 corners retired or justified, the
seven-criterion FW-style validation suite in place, full suite
5326/5326 GREEN. **But the figure remains inadequate** — its pole wedge
is only ~5 % filled (a sparse pole scatter, not the solution surface).
ADR-0025's scope did not include bringing the vector path-network up to
FW's full region-filling driver, which the figure actually needs. See
**Amendment 13**, worklog 059, and bead `padetaylor-0ln.40` (the
plan-first deep-dive that does that port). Amendments 1–12 record the
Phase A–F exploration outcomes and auditions.
**Bead**: `padetaylor-0ln.37` (this ADR + the child-bead plan below).
**Plan**: `docs/v0p2_plan.md`; supersedes the v1-corner ledger of
worklog 057 (`docs/worklog/057-whole-plane-kkg-surface-figure.md`).

## Decision

The headline figure of the v0.2 vector solver — the P_I⁽²⁾ tritronquée
`V₀(x,0)` over the complex x-plane `[-20,20]²`, reproducing KKG 2015
Figs 7.4/7.5 — is **re-resolved to senior-engineer grade**. Three
binding decisions:

1. **Validation-first, because there is no ground truth.** KKG 2015
   plotted the *surfaces* only; they drew **no poles**, published **no
   pole table**, and explicitly state their own figures serve "mainly
   for visualization" (KKG TeX `:3215-3217`). The pole field in the
   `arg x ≈ 0` wedge was **never numerically studied by KKG at all**.
   There is therefore nothing external to validate against. The figure
   is certified instead by a **manufactured, FW-style validation
   suite** (the *Validation Criteria Menu* below) — exactly the
   discipline Fornberg & Weideman used for their P_I pole fields, which
   had no oracle either (FW 2011 `:303-311`: conjugate-symmetry error
   estimate, double-run accuracy indicator, endpoint-derivative match).

2. **The wedge goal is the best-possible pole field — FW-grade.** Not a
   visualization-only compromise. Full *honest* coverage is pursued;
   "honest" means **no Padé approximant is ever evaluated outside its
   verified disc of validity** (the `extrapolate=true` corner is
   deleted). Bead `0ln.37`'s acceptance criterion (3) literally reads
   "full `|x|≤20` wedge coverage"; because the tritronquée pole density
   grows with `|x|`, whether honest coverage *reaches* `|x|=20` is a
   genuine numerical unknown. Criterion (3) is therefore **reframed**:
   honest coverage to the **numerically-supported frontier**, the
   frontier itself a *measured and reported* quantity. "As good as FW
   managed" = honest to the frontier, not a hard radius mandate.

3. **Every v1 corner is retired or rigorously justified** (the
   *v1-Corner Retirement Ledger* below). The triple-method majority
   vote (ADR-0024) **stays**.

## Context — the six v1 corners and the three levers

Worklog 057 §"Required follow-up" and the five-agent recon of this bead
catalogue the shipped figure's v1 corners with exact line numbers in
`figures/_kkg_pi2_surface_helpers.jl`:

- **C1** inner-arc Dirichlet datum is the `n_terms=2` KKG asymptotic
  series, `O(10⁻²)` at `|x|≈2` (`:419-432`);
- **C2** wedge Stage-2 called with `extrapolate=true` — shared-Q Padé
  approximants evaluated outside their valid disc (`:551`);
- **C3** hand-tuned wedge step `h = 0.1` (`:215`);
- **C4** `±3°` Stokes-strip NaN mask (`:224`, `:250`);
- **C5** bilinear polar ray-fan voter (`:339-354`);
- **C6** `121²` grid + `6°` ray fan — resolution-limited (`:176`,
  `:195`).

Three numerical levers drive the re-resolution:

- **Lever 1 — true-radius Stage-2 gating.** `VectorPathNetworkStage2`
  currently gates each node's validity at the *step* `h` (a disc of
  radius `h ≈ 0.1`), producing a coarse Voronoi tessellation. A node's
  shared-Q Padé is in fact valid over a far larger disc. **What that
  disc actually is is not known a priori** — recon proposed `h·min|t*|`
  (distance scaled by the nearest denominator root), but a shared-Q
  Padé *models* the poles `Q` captures and stays accurate *through*
  them, so the honest radius is bounded by the nearest singularity `Q`
  did *not* resolve, not by the nearest root. This is resolved
  empirically in Phase A1.

- **Lever 2 — the dense-wedge walk.** Tiling the wedge toward `|x|=20`
  marches through a pole field of growing density; a fixed-`h` step
  will eventually land on a pole and degenerate the shared-Q jet. The
  principled fix is bead `0ln.23` (maximise distance to the nearest
  shared-Q root) plus adaptive `h`. Absorbed here as B2.

- **Lever 3 — sector refinement.** Mechanical: voters 2/3 of ADR-0024
  are spectrally callable at any point, so a finer grid + denser rays +
  higher Laplace `N` cost only compute. The one subtlety is C5.

## The P_I⁽²⁾ local pole structure (an exact independent certificate)

The Painlevé property of P_I⁽²⁾ (KKG TeX `:3124-3125`, citing
Shimomura: "all its solutions are meromorphic functions of `x` and
`t`") guarantees every movable singularity is a **pole of fixed order**
with a Laurent expansion whose coefficients obey a recursion fixed by
the ODE. Dominant balance: near a pole `u ~ A(x-x₀)^p` in
`u'''' + 10(u')² + 20uu'' + 40(u³+6x) = 0` (`t=0`). All four of
`u''''`, `(u')²`, `uu''`, `u³` co-balance at order `(x-x₀)⁻⁶` iff
`p = -2` — a **double pole**. At that order:

```
120A + 40A² + 120A² + 40A³ = 0   ⟹   A(A+1)(A+3) = 0   ⟹   A ∈ {-1, -3}.
```

**Every genuine P_I⁽²⁾ pole has leading coefficient exactly `-1` or
`-3`** (contrast P_I, whose poles are `+1·(z-z₀)⁻²`, FW 2011 `:59-61`).
This is a per-pole structural certificate that owes **nothing** to the
path-network algorithm that located the pole — the strongest
independent validation criterion available. Its full form (which of the
two families occur for the tritronquée; the Fuchs/resonance indices;
the residue) is derived in Phase A3 and recorded by amendment.

## Validation Criteria Menu

The figure-certifying suite. `[E]` = machinery exists in the codebase;
`[N]` = new work. Each load-bearing criterion is mutation-proven
(Rule 4).

| ID | Criterion | Tier | Status |
|----|-----------|------|--------|
| VC-1 | Non-tautological ODE residual (FD of `u'''` vs equation side) | strong | `[E]` active |
| VC-2 | Sign + KKG `n=2` series cross-check, negative real axis | strong | `[E]` active |
| VC-3 | Wedge confinement + pole-free-sector purity of poles | strong | `[E]` active |
| VC-4 | Dominant-balance leading coefficient `A ∈ {-1,-3}` per pole | strong | `[N]` shipped (beads `0ln.37.12`; `figures/_kkg_pi2_vc45.jl`) |
| VC-5 | Conjugate-symmetry pole pairing (`V₀(x̄)=conj V₀(x)`) | strong | `[N]` shipped (beads `0ln.37.13`, `0ln.37.20`; `figures/_kkg_pi2_vc45.jl`; Amendment 6 optimal matching) |
| VC-6 | Cross-node support filter (`min_support ≥ 2`) | medium | `[E]` active |
| VC-7 | Loop-closure ΔP_rel certificate (`quality_diagnose`, ADR-0016) | medium | `[N]` shipped (bead `0ln.37.14`; vector adapter in `ext/PadeTaylorDiagnosticsExt.jl`; Amendment 7) |
| VC-8 | BVP endpoint higher-derivative match (FW §5.2 diagnostic) | medium | `[N]` shipped (bead `0ln.37.15`; `figures/_kkg_pi2_surface_helpers.jl` `surf_vc8_companion_check`; Amendment 8) |
| VC-9 | Weierstrass-℘ oracle for the *vector* pipeline (FW Table 5.1) | medium | `[N]` shipped (bead `0ln.37.16`; `test/vector_pipeline_oracle_test.jl`; Amendment 8) |
| VC-10 | Two-run different-path pole-disagreement accuracy indicator | medium | `[N]` shipped (bead `0ln.37.17`; `src/VectorPathNetwork.jl` `rng` kwarg + `figures/_kkg_pi2_surface_helpers.jl` `surf_vc10_two_run`; Amendment 9) |
| VC-12 | 7-fold rotational-symmetry cross-check (`V₀`→`V₁` poles) | strong | `[N]` stretch |

(VC-11 Froissart residue filter and VC-13 pole-count growth law are
recorded in the recon as lower-priority; deferred with a forcing
condition rather than scheduled.)

## v1-Corner Retirement Ledger

| Corner | Disposition | Bead |
|--------|-------------|------|
| C1 inner-arc asymptotic datum | **Retired** (Amendment 10) — each sector ray BVP extends inward to `R_INNER_BC = 1.05 < R_MIN`, so the inner arc is a BVP-*interior* datum (`~1·10⁻⁴`, vs the series' `~5·10⁻⁴`); voters 2/3 source all four `w`-rectangle edges from the ray fan, the asymptotic series gone from the voter-boundary path | C2 (`0ln.37.10`) |
| C2 `extrapolate=true` | **Retired** — Lever-1 true-radius gate; no out-of-disc evaluation | B1 |
| C3 hand-tuned `h=0.1` | **Retired** — Lever-2 adaptive / shared-Q-root step | B2 |
| C4 `±3°` Stokes mask | **Retired-or-justified** (Amendment 12) — Phase-E1 audition: the sector ODE solve is healthy *at* the 36° Stokes line and the triple-method vote stays inside the figure's honesty envelope down to a `1°` fan margin, so `SURF_SECTOR_MARGIN_DEG`/`SURF_STITCH_MASK_DEG` are narrowed `4°/3° → 1°/1°` — the grey strip shrinks `6° → 2°` (3× filled with audited-honest harmonic vote); the residual `±1°` mask is retained with a connection-formula forcing condition (deferred bead) | E1 (`0ln.37.18`) |
| C5 bilinear ray-fan voter | **Retired** (Amendment 10) — audition outcome: a harmonic voter 1 `≡ voter 2` (rejected — collapses ADR-0024 independence), plain bilinear materially too lossy (rejected); voter 1 reconstructed exact-radius (barycentric) + cubic-angular (Catmull-Rom), keeping it a direct ODE solve | C3 (`0ln.37.11`) |
| C6 `121²` grid / `6°` rays | **Retired** (Amendment 11) — C1 (`0ln.37.9`) lifts the Cartesian raster `121² → 401²` and the ray fan `6° → 2°` (≈ 141 rays); both inside the Phase-F runtime gate (figure test `2m22s`). The coupled "raise the Laplace `N` too" instruction was *not* taken speculatively — the C1 convergence probe falsified the Amendment-10 under-resolution hypothesis (the inner-arc spread is the asymptotic-seed truncation floor, irreducible by any `N`); `SURF_LAP_NX/NY` set to the measured `48/64` spectral knee | C1 (`0ln.37.9`) |

## The phased child-bead plan

Exploration spikes (Phase A) run **first** — they resolve genuine
unknowns and pre-committing the architecture before A1/A2 would violate
Rule 9. All coding is serial (Rule 7: one Julia process).

- **Phase A — exploration & auditions** (spikes, throwaway code):
  A1 Padé validity-radius study (4-gate audition → Lever-1 criterion);
  A2 wedge tractability probe (→ the honest frontier);
  A3 P_I⁽²⁾ Laurent-structure derivation (→ VC-4 exact form);
  A4 baseline-quality probe (`quality_diagnose` + conjugate pairing on
  the *current* V8b walk).
- **Phase B — wedge levers**: B1 true-radius Stage-2 gate;
  B2 shared-Q-root adaptive walk (absorbs `0ln.23`);
  B3 principled extended target fan; B4 Stage-2 interpolation audition
  (Voronoi vs in-disc blend).
- **Phase C — sector refinement**: C1 high-res grid/rays/Laplace-N;
  C2 BVP-sourced arc data; C3 ray-fan voter audition.
- **Phase D — validation suite**: VC-4, VC-5, VC-7, VC-8, VC-9, VC-10
  (one bead each; VC-12 stretch).
- **Phase E — Stokes strips**: audition the `±3°` fill.
- **Phase F — re-render, full `Pkg.test()` GREEN, worklog 058,
  ADR amendment with the A1/A2/A3 outcomes.**

## Alternatives considered, rejected

- **Crank grid + target numbers, keep `extrapolate=true`.** Rejected —
  the explicit anti-pattern of bead `0ln.37` and worklog 057: it
  produces a denser figure that is still dishonest (Padé evaluated
  outside its disc) and still uncertified.
- **Render the wedge as a pole scatter over a clamped `|u|` heatmap
  and call it KKG-faithful.** Rejected — KKG published no pole figure,
  so "faithful" has no referent; a scatter would be a *lower-ambition*
  deliverable masquerading as fidelity. The user's directive is the
  best-possible pole field.
- **Bake `h·min|t*|` as the Stage-2 gate** (the recon's first guess).
  Rejected as a *decision* — it is a plausible hypothesis but the true
  validity disc of a pole-modelling Padé is bounded by the nearest
  *uncaptured* singularity. Promoted to an audition (A1).
- **Validate against a high-resolution independent re-solve as
  "ground truth".** Rejected — any such re-solve uses the same
  Painlevé/Padé/BVP substrate; agreement would be partly circular. The
  manufactured criteria (conjugate symmetry, ODE residual,
  dominant-balance `A`, loop closure, the ℘ oracle) are genuinely
  independent of the path-network algorithm under test.

## Consequences

- This ADR governs `figures/_kkg_pi2_surface_helpers.jl`,
  `figures/kkg_pi2_tritronquee_surface.jl`, and their test
  `figures/test_kkg_pi2_surface.jl`; it may touch
  `src/VectorPathNetworkStage2.jl`, `src/VectorPathNetwork.jl`, and
  `src/VectorPoleField.jl` (Levers 1–2). Any `src/` change stays
  additive (ADR-0001) — v0.1 + v0.2 scalar tests remain bit-identical.
- New `src/` work obeys the Rule 6 ≤200-effective-LOC cap and ships a
  literate docstring + paired test with a mutation-proof block.
- The triple-method majority vote (ADR-0024) is preserved verbatim;
  C2 only changes the *source* of voters 2/3's arc boundary data.
- Bead `0ln.23` is closed by B2 (its forcing condition — "min-‖y‖
  steps a walk onto a pole in a production figure" — is now met).
- Fail-loud (Rule 1): a Stage-2 query beyond the true validity disc
  with `extrapolate=false` returns `NaN` (an honest gap), never a
  silently-extrapolated value.

## Amendment 1 — Phase A exploration outcomes (2026-05-21)

The Phase-A spikes (probes committed under `external/probes/`) resolve
the two parameters the ADR left to exploration; this amendment is the
ADR-level summary of their reports.

### A1 — Stage-2 validity gate (bead `0ln.37.1`) — RESOLVED

Five gate criteria auditioned on the 389-node V8b wedge walk:

- Gates (ii) `h·min|t*|` and (iii) nearest-uncaptured singularity:
  **rejected with evidence** — they over-estimate the honest radius
  5.7× / 30×. The shared-Q roots are genuine *far* poles; the Padé in
  fact fails by **order-24 truncation error** long before reaching any
  pole. The honest radius is *truncation-limited*, not pole-limited.
- Gate (i) `h`: unsound — over-estimates 36 % of nodes (the current
  `extrapolate=false` default is already mildly dishonest).
- Gate (iv) empirical (overlap + order-64 Taylor oracle): the
  ground-truth anchor — honest radius median ≈ `1·h`, spread `0.6h–2h`.
- **Winner — gate (v-b), per-node Jorba–Zou.** Feed each node's
  rescaled order-24 jet to the existing `vector_step_jorba_zou`; it
  honestly recovers **7–10.5× more disc area** than a global `κ·h`
  gate at equal (≤2 %) over-estimation. A global constant cannot track
  the 3.3× per-node spread; the per-node estimator collapses it to 1.5×.

Production gate for B1 (`_stage2_fill`):
```
R_gate = min( s(tol)·h_JZ·h_v ,  0.5·h_v·min|t*| )
  h_JZ    = vector_step_jorba_zou(rescaled order-24 jet, tol)
  min|t*| = min |root of the node's shared-Q|
  s(1e-6, 1e-8, 1e-10) = 0.34, 0.36 (default), 0.30
```
The truncation term binds for ~98 % of nodes; the `0.5·h·min|t*|`
pole-adjacency clamp binds for ~2 % (deep-wedge nodes with a pole
inside the disc). Honest by construction (`s<1` and the clamp both
under-cover; gaps → `NaN`). **B1 sub-task:** cache the order-24 jet in
the Stage-1 walk so the gate costs nothing extra.

**Coupling finding:** B1 alone makes the figure *holier*, not bigger —
honest discs are `med≈0.06, p90≈0.16` vs the dishonest `extrapolate=
true`. Coverage must come from a denser walk: **B1 and B2 ship
together.**

### A3 — P_I⁽²⁾ local pole structure (bead `0ln.37.3`) — RESOLVED

Dominant balance confirmed: every movable singularity is a **double
pole** (`p=-2`), leading coefficient `A ∈ {-1,-3}`. Painlevé-test
Laurent recursion `u = Σ aₖ(x-x₀)^{k-2}`:

- `A=-1`: resonances `{-1,2,5,8}` — generic 4-parameter family
  (V₀'s wedge poles expected here);
- `A=-3`: resonances `{-3,-1,8,10}` — codimension-1 special family;
- both families: **`a₁ = 0` exactly** — the residue vanishes
  universally (as for P_I); no logarithmic terms.

VC-4 acceptance form (final): per extracted pole `x₀`, fit
`u ≈ A(x-x₀)⁻² + B(x-x₀)⁻¹ + C` by complex least squares on a ring of
32 points at `r = min(0.05, 0.1·d_nearest)`:

- **VC-4a** `min(|A+1|, |A+3|) < 0.10` — the ODE-structure check;
- **VC-4b** `|B| < 0.10·|A|` — the zero-residue / Painlevé-property
  check.

Both are mutation-proof and independent of the path-network extraction
algorithm under test.

## Amendment 2 — A2 wedge-tractability outcome: the wedge panel is rescoped (2026-05-21)

### A2 — wedge tractability (bead `0ln.37.2`) — RESOLVED, with a scope consequence

The A2 spike (`external/probes/wedge-tractability/`) measured the honest
union coverage achievable by the path-network walk + per-node A1-gated
disc, extending the target fan toward `|x|=20`:

- solver tolerance (`1e-8`): honest union coverage of the `|x|≤20`
  wedge **saturates at ~8.6 %**;
- display tolerance (`1e-3`, legitimate for a clamped heatmap/3D
  surface — the colour channel resolves no finer): **~13.5 %**;
- display tol **and** Taylor order doubled to 48: **~18 %** — the
  ceiling.

Two independent obstructions bind: (i) the densest *completing*
path-network fan tops out at ~6850 visited nodes while ~13–23k honest
discs are needed to tile the area-251 wedge, and denser fans **block**
(a bridging walk lands on a pole, degenerating the shared-Q jet);
(ii) an order-48 *walk* degenerates the shared-Q Padé linear solve
outright. **A filled honest wedge surface is not viable** with the
local-Padé-tiling architecture — across three orders of tolerance and a
doubled order, coverage moves only 8.6 % → 18 %.

### Decision — the wedge panel is the validated pole FIELD

The headline figure's wedge panel is **rescoped** to its honest,
achievable, senior-grade form:

- **The wedge deliverable is the extracted, FW-validated pole field —
  the pole *locations*** — rendered in the FW 2011 Fig 4.7 idiom (a
  scatter, `figutil.jl` `pole_scatter_axis`). This is not a lesser
  deliverable: FW's own "pole-field" figures *are* pole-location plots,
  and "the pole field as good as FW managed" denotes exactly this. The
  pole field is the scientifically meaningful content of the pole-rich
  wedge — and it *is* achievable: a walk filament reaches `|x|=20` at
  any step size, so poles can be extracted (each validated by
  VC-4/VC-5/VC-6) across the whole window.
- **An honest `|u|` surface underlay** is rendered only where the
  A1-gated node discs genuinely cover (~13–18 % of the wedge); the
  remainder is `NaN`/grey. No Padé is evaluated outside its disc — the
  partial surface is honest, the gaps explicit.
- The smooth ~270° **sector** panel is unaffected — it remains the
  fully-resolved, 100 %-covered triple-method-voted harmonic surface
  (Phase C).

### Bead-scope consequences

- **B3** (`0ln.37.7`) rescopes from "tile the wedge area" to "thread
  the walk through the wedge to **extract the pole field** densely and
  accurately to `|x|=20`" — pole extraction in the FW Fig 5.1 threading
  idiom, not area tiling. Pole-field completeness/accuracy is a
  first-class deliverable.
- **B4** (`0ln.37.8`) narrows to the honest-surface *underlay* blend
  over the ~13–18 % covered region.
- **B1, B2** stand: B1 (true-radius gate) governs both the honest
  surface underlay and honest pole extraction; B2 (adaptive walk)
  threads the walk to `|x|=20` without degenerate steps.
- Bead `0ln.37` acceptance criterion (3) is finalised: *honest* wedge
  coverage = a complete validated pole field to `|x|=20` + an honest
  partial `|u|` surface underlay; no extrapolation anywhere.

### Deferred — a fully-filled honest wedge surface

A genuinely filled honest wedge surface would need a different
architecture — a 2D lattice of locally re-expanded Padé patches
(~13–23k of them) rather than walk filaments. Its error-compounding
behaviour under repeated lateral re-expansion is unproven and not
senior-grade-certain. Recorded as a deferred bead under the v0.2 epic;
forcing condition: a figure requirement for a filled wedge surface
beyond the pole field.

## Amendment 3 — A4 baseline; Phase A complete (2026-05-21)

### A4 — V8b baseline-quality probe (bead `0ln.37.4`) — RESOLVED

Baseline metrics of the *current shipped* V8b wedge walk (389 nodes,
21 poles) — the before-picture the re-resolution must beat (Phase F):

- **Loop-closure ΔP_rel** (767 non-tree Delaunay edges): median
  `9.9e-4` but p90 `0.996`, p99 `1.00`; categories well_closed 5.2 % /
  noisy 24.8 % / extrap_driven 65.6 % / depth_driven 4.4 % — **~70 % of
  loop closures catastrophic** (vs the FFW Fig 1 scalar baseline of
  6.3 %). Worst edges `extrap_max` up to 53 — midpoints 50× outside the
  canonical disc, the direct signature of the ungated `extrapolate=
  true` walk (corner C2). `bad_centroid ≈ 4.82−0.24im`.
- **Conjugate symmetry badly violated**: the 21 poles split 15 upper /
  6 lower; only 6 forced pairs, residual median ≈ `1.25` — no pole has
  a genuine conjugate partner. The VC-5 defect, confirmed.
- **A spurious pole**: VC-4 spot-check — 5 of 6 sampled poles fit
  `A=−1` to `<1.1e-5` (confirming the A3 generic family), but the pole
  at `+5.44+2.40im` fits `A≈0` — not a double pole. At least one of the
  21 V8b poles is spurious.

The shipped figure thus has three genuine, now-quantified defects;
the re-resolution is well-motivated.

### D-VC7 scope — `quality_diagnose` needs a vector adapter

`quality_diagnose` is typed for the scalar `PathNetworkSolution`; the
vector walk produces `VectorPathNetworkSolution`. Bead D-VC7
(`0ln.37.14`) must add an additive `quality_diagnose(::VectorPath
NetworkSolution)` method in `ext/PadeTaylorDiagnosticsExt.jl`: shared-Q
`Pᵢ(t)/Q(t)` Horner eval instead of `_evaluate_pade`; ΔP_rel as the
vector 2-norm; drop the sheet-0 mask (the P_I⁽²⁾ companion is
single-sheeted). The Delaunay/tree/LCA/categorise machinery ports
verbatim; the A4 probe `§2b` is a working reference implementation.

### Phase A complete

A1–A4 done; the architecture is fully determined. Phase B: B1
true-radius gate (gate v-b), B2 adaptive walk, B3 pole-field extraction
to `|x|=20`, B4 honest-underlay blend.

## Amendment 4 — B2 dense-wedge walk shipped (2026-05-21)

B2 (bead `padetaylor-0ln.37.6`, Lever 2; discharges and closes bead
`padetaylor-0ln.23`) is implemented in the new module
`src/VectorWedgeStep.jl`, `include`d before `src/VectorPathNetwork.jl`
(the Rule-6 ≤200-LOC split — `VectorPathNetwork.jl` drops to 152
code-LOC). Two coupled changes, now the **default** for vector walks:

- **Direction selector `:max_q_root`.** The V7 min-‖y‖ proxy is replaced
  by the principled criterion: of the five wedge candidates, choose the
  one whose landed node has the largest distance to its nearest
  shared-`Q` denominator root — the most pole-free disc. Exact, not a
  proxy. `:min_y` is kept as an opt-in. A candidate whose canonical
  shared-`Q` store cannot be built is *excluded* (the driver rebuilds
  that store — picking it would crash); all-excluded throws (Rule 1).
- **Adaptive `h`.** `h` is capped per-step by the parent node's
  nearest shared-`Q` pole (`POLE_SAFETY·h·min|t*|`) and grows
  geometrically toward `h_max` where the pole field is sparse. `h_max`
  is the `h` kwarg; the floor is `h_min = H_MIN_RATIO·h_max` (`1e-3`);
  a step forced below `h_min` throws (an honest wedged-walk signal).
  C3 (the hand-tuned figure `h = 0.1`) is retired.

**ADR-0021 deviation.** ADR-0021's "No vector `step_pade_root`" note
anticipated reusing the FW 2011 §3.1 *directional* `step_pade_root` on
the shared `Q` as the adaptive cap. B2 deviates deliberately: the cap
is **direction-agnostic** (`min|t*|`, the nearest pole in *any*
direction). The reason is the coupling with `:max_q_root` — the
selector picks the step *direction* to dodge a pole, so a directional
cap shrinks `h` for a pole the walk steers around (it stalled the
figure walk to `h ≈ 3·10⁻⁶` in testing). `_adaptive_h`'s docstring
records this in full.

**Measured outcome.** On the P_I⁽²⁾ tritronquée wedge: the fixed-`h=0.1`
walk BLOCKS the A2 extended target fan past `|x|≈8` with the
`shared_denominator_pade` degeneration A2 documented; the B2 adaptive
`:max_q_root` walk threads the same fan to `|x| ≈ 18–20`, extracting a
dense validated pole field. `kkg_pi2_figure_test.jl` is GREEN (the
Stage-B march now succeeds under the B2 default). Bead `0ln.23` is
closed by this amendment.

**Known cross-figure consequence.** Changing the *default* walk policy
affects every vector-walk figure call site. The A_4⁽¹⁾ Noumi–Yamada
figure (`figures/_noumi_yamada_a4_helpers.jl`, V8a) relies on the
default; under `:max_q_root` its order-24 shared-`Q` SVD degenerates at
specific node states the new walk's path reaches (`:min_y` and
`h=0.2` complete; `h=0.3` blocks) — `test/noumi_yamada_a4_figure_test.jl`
errors. B2's scope is the P_I⁽²⁾ headline figure (`figures/` files
untouched per the B2 bead); the A_4 figure needs a B3-style explicit
`step_policy = :min_y` pin or its own re-resolution — a follow-up bead
for Phase F triage.

## Amendment 5 — Phase-D VC-4 / VC-5 per-pole validation shipped (2026-05-21)

VC-4 (bead `padetaylor-0ln.37.12`) and VC-5 (bead `0ln.37.13`) — the
first two Phase-D validation criteria — are implemented in the new
figure helper `figures/_kkg_pi2_vc45.jl` (`include`d by the kernel
`figures/_kkg_pi2_surface_helpers.jl`; the Rule-6 file-size split).
The kernel's wedge `poles` field is now the **VC-4-validated** set:
`extract_poles_shared_q` produces a *candidate* field; `vc4_validate`
prunes every candidate that is not a genuine P_I⁽²⁾ double pole.

**VC-4 — measured outcome on the B3 wedge field.** The B3 extended
threading fan extracts **380 candidate poles** (`min_support ≥ 2`,
VC-6). VC-4 fits each candidate to `u ≈ A·ξ⁻² + B·ξ⁻¹ + C` on a
32-point ring and applies VC-4a (`min(|A+1|,|A+3|) < 0.10`) + VC-4b
(`|B| < 0.10·|A|`):

- **266 candidates survive** VC-4 — all in the `A = -1` generic family
  (`m1`); **zero** `A = -3` poles, confirming the A3 §5.1 genericity
  prediction (the tritronquée's wedge poles are the generic family);
- **114 pruned** — 48 `:froissart` (`|A| < 0.1`, no double-pole
  structure — the Froissart-doublet signature the A4 baseline first
  flagged at `+5.44+2.40im`), 56 `:nonzero_residue` (right `A`,
  `|B|` too large — a mis-located cluster), 10 `:out_of_family`
  (a double pole but `A ∉ {-1,-3}`).

So ~30 % of the raw VC-6 candidate field is spurious — VC-4 is
load-bearing, not cosmetic.

**VC-5 — measured outcome.** The 266 VC-4-survivors split 131 upper /
127 lower / 8 near-real — a balanced field (contrast the V8b baseline's
15 / 6 catastrophe, ADR Amendment 3). `vc5_pair` matches them by a
globally-greedy conjugate assignment (`VC5_MATCH_TOL = 0.5`): **93
conjugate pairs**, pairing-residual **median 0.24, max 0.50** — an
FW-style accuracy estimate of the pole field, ~5× better than the V8b
baseline median 1.25. The field's pole nearest-neighbour spacing is
~0.69, so a residual of 0.24 is genuine conjugate symmetry at roughly
a third of the pole spacing. VC-5 is a **diagnostic** (it reports, it
does not prune): 72 off-axis survivors have no mirror partner within
`VC5_MATCH_TOL` and are **flagged suspect** — the far-wedge
(`|x| ≳ 12`) field is extracted asymmetrically by the walk filament,
the residual Phase-F must still close.

Both criteria are **mutation-proven** (Rule 4): forcing the VC-4a
family set to `{-2}` empties the validated field (RED); disabling the
prune fails the shrink / `n_pruned > 0` / per-pole VC-4a/b assertions
(RED). A deliberately-injected fake `A ≈ 0` location is rejected as
`:froissart` while a genuine field pole through the same path is kept
— the in-test mutation-proof. `figures/test_kkg_pi2_surface.jl` PI2S.4
asserts all of the above; the suite is GREEN.

## Amendment 6 — VC-5b: the conjugate-pairing defect re-resolved (2026-05-21)

Bead `padetaylor-0ln.37.20` (Phase D, "VC-5b") investigated and fixed
the defect Amendment 5 reported — 72 of 266 VC-4-validated wedge poles
flagged as conjugate-unpaired.

### The investigation — the root cause is *not* a missed walk region

Amendment 5 hypothesised the far wedge (`|x| ≳ 12`) pole field was
extracted *asymmetrically* because the single path-network tree's
upper- and lower-half branches threaded different node sequences and
the lower branch *missed* poles the upper branch found. The VC-5b
instrumentation **falsified that hypothesis**:

- **Node coverage is conjugate-symmetric.** For every upper-half
  visited node `z`, the nearest lower-half node to `conj(z)` is a
  median `0.10` away — and `0.10` even when restricted to the far
  wedge `|x| ≥ 12`. The two half-walks thread *near-mirror loci*; the
  walk does **not** miss far-wedge regions.

- **The candidate pole field is asymmetric** (235 upper / 136 lower
  raw candidates), but a per-pole VC-4-residual *polish* — coordinate
  descent on `min(|A+1|,|A+3|)` over the walk's own shared-`Q` field —
  moves each pole only `~0.005`. **Every extracted pole already sits at
  its dominant-balance-optimal location.** The poles are individually
  accurate; what differs between the halves is *which subset of the
  dense far-wedge pole lattice each half resolves at all*.

- **The field carries an intrinsic conjugate-residual accuracy of
  median `~0.3`, p75 `~0.5`.** Fitting a Laurent model exactly *at*
  `conj(p)` of a flagged pole `p` returns `A ≈ 0` (no double pole
  there), but a VC-4-fit scan finds the genuine partner pole offset a
  median `0.51` from `conj(p)`. This is the A2 tractability ceiling
  (Amendment 2): the two half-trees develop different adaptive-`h`
  histories and the order-24 shared-`Q` approximants of one half
  resolve a slightly different subset of the lattice than the other.

The candidate fixes were tested empirically and **rejected with
evidence**: two independent half-wedge walks give *more* flagged poles
(73–81 vs 72) because each half-walk has fewer near-real nodes
informing the deep wedge; a denser per-half fan **blocks** (the A2
degeneration); angle-offset or `h`-varied second walks **block** or
give a strictly worse field (`h = 0.08`: 93 flagged; `h = 0.06`: 109);
a three-walk consensus merge only inflates the pole count. The single
whole-wedge tree at the production parameters is **measurably the
best-conditioned walk** — no walk-architecture change beats it. The
~0.3 conjugate residual is the genuine A2 ceiling, not a correctable
walk defect.

### The fix — VC-5's *own* two defects

Since the walk cannot be improved, the VC-5b fix targets the two
defects the investigation found in `vc5_pair` *itself*:

1. **Greedy → globally-optimal matching.** The v1 `vc5_pair` matched
   by a *globally-greedy* commit (tightest admissible mirror pair
   first). Greedy is **not** maximum-cardinality — it commits a
   locally-tight pair that blocks two other poles from each finding
   their *only* admissible mirror. Measured on the B3 field: greedy
   finds 93 pairs, the maximum-cardinality matching (Kuhn augmenting
   paths) finds **103**. `vc5_pair` now computes the optimal matching.

2. **`VC5_MATCH_TOL` re-derived from measured ground truth.** The v1
   `0.5` was a *guess* — "below the 0.69 pole spacing" — set before the
   field's conjugate-residual accuracy was measured. `0.5` is *tighter
   than the field's own accuracy* (~0.3 median / ~0.5 p75): it bisects
   the residual distribution and rejects ~27 % of genuine mirror pairs.
   Re-derived from the *measured* 0.69 spacing / ~0.3–0.5 accuracy to
   **`0.6`** — still firmly below the 0.69 spacing (87 % of it), and
   with the optimal matcher (a pole pairs only with its *globally* best
   mirror) it cannot false-match two distinct lattice poles. This is a
   Law-1 ground-truth correction of a guessed parameter.

**Measured outcome.** On the 266-pole B3 field: **72 → 52 flagged**
(a 28 % reduction, the VC-5b deliverable), **93 → 103 conjugate
pairs**. The pairing residual *widens slightly* — median `0.24 → 0.29`,
max `0.50 → 0.60` — precisely because the optimal matcher no longer
discards the harder-to-pair tail; this is the residual becoming a
*more honest* diagnostic, not a cosmetic improvement. No conjugate
symmetry is imposed by construction (the hard constraint): the upper
and lower poles are still extracted from disjoint node sets of one
independent walk, and `vc5_pair` only *pairs* poles — it never
*constructs* a pole from its mirror — so the residual stays a genuine
FW/FFW-style accuracy cross-check (FFW 2017 `:120-124`).

**Honest residual — the remaining 52.** The 52 still-flagged off-axis
poles are genuine poles (VC-4 dА ≈ 0) whose conjugate the field
resolves at a lattice location farther than `0.6` away. This is the
A2 far-wedge tractability ceiling, now *quantified and reported* as
the VC-5 diagnostic rather than mis-attributed to a walk defect. A
materially smaller flag count would need a different wedge
architecture (the 2D-lattice re-expansion of the Amendment-2 deferred
bead) — recorded there, not re-opened here.

**Mutation-proven** (Rule 4): reverting `vc5_pair` to the greedy
commit turns `figures/test_kkg_pi2_surface.jl` PI2S.4 RED — the
load-bearing `length(flagged) ≤ 55` assertion fails (greedy flags 58),
as do the two in-test greedy-vs-optimal comparison assertions. The
bound `55` is set strictly below the greedy result so the optimal
matcher is genuinely load-bearing, not carried by the tolerance change
alone. The suite is GREEN with the optimal matcher.

## Amendment 7 — Phase-D VC-7 loop-closure certificate: the vector adapter + a metric-scope correction (2026-05-21)

Bead `padetaylor-0ln.37.14` (Phase D, "VC-7") brings the loop-closure
quality certificate (`quality_diagnose`, ADR-0016) to the headline
figure's vector wedge walk and wires it in as a figure acceptance gate.

### The vector adapter — shipped

`quality_diagnose` was typed for the scalar `PathNetworkSolution`; the
figure's wedge walk produces a `VectorPathNetworkSolution`. The
additive method `quality_diagnose(::VectorPathNetworkSolution)` in
`ext/PadeTaylorDiagnosticsExt.jl` (Amendment 3 §"D-VC7 scope" is its
spec) differs from the scalar method in exactly three ways: it
evaluates each node's shared-`Q` approximant `Pᵢ(t)/Q(t)` by Horner
(the vector node stores `visited_numerators` + one
`visited_denominator`, not a scalar `PadeApproximant`); it generalises
`ΔP_rel` to the vector 2-norm `‖y_A − y_B‖ / (‖y_A‖ + ‖y_B‖ + ε)`; and
it drops the sheet mask (the `P_I⁽²⁾` companion is meromorphic —
single-sheeted — and `VectorPathNetworkSolution` has no `visited_sheet`
field). The Delaunay / tree-edge / LCA / categorise / aggregate
machinery ports verbatim; the A4 probe `§2b`
(`external/probes/v8b-baseline/probe.jl`) was the verified reference.
The adapter is purely additive (ADR-0001) — the scalar
`quality_diagnose` stays bit-identical (`test/diagnose_test.jl` 32/32
GREEN, unchanged). A new vector-adapter testset (`VDG.1`–`VDG.5`)
mutation-proves it: dropping the `/Q` from the shared-`Q` evaluation
turns the well-closed-lobe assertion RED. `DelaunayTriangulation` is
added to `figures/Project.toml` (a figures-project dep — it activates
the extension without becoming a hard dep of the core package).

### A4 §2b over-attribution corrected — VC-7 is *not* a B1/B2 proxy

The A4 baseline (`external/probes/v8b-baseline/REPORT.md` §2b) measured
the V8b walk's loop-closure ΔP_rel as median `9.9e-4`, p90 `0.996`,
~70 % catastrophic, and attributed the catastrophe to "the direct
signature of the ungated `extrapolate=true` walk (corner C2)". **D-VC7
measurement falsifies that attribution.** `quality_diagnose` evaluates
each node's *canonical* shared-`Q` Padé at Delaunay-edge midpoints; that
metric is governed purely by the geometry `Delaunay-edge-length vs
canonical step h`. It is **independent of the Stage-2 `extrapolate`
flag and the B1 true-radius gate** — those govern only the separate
Stage-2 grid fill, a code path `quality_diagnose` never touches. So the
B1/B2 levers cannot, and do not, move the `quality_diagnose`
distribution.

### Measured before/after — the re-resolved walk is *comparable*, not better

The vector adapter was run on both walks (the adapter bit-reproduces
the A4 §2b numbers on the old V8b walk — median `9.915e-4`, well 5.2 %,
noisy 24.8 %, extrap 65.6 %, depth 4.4 %, 767 edges — confirming it is
correct):

| metric (non-tree Delaunay edges) | old V8b (`:min_y`, fixed `h`, `|x|≤8` fan) | re-resolved B3 (`:max_q_root`, adaptive, `|x|≤20` fan) |
|----------------------------------|---------------------|-----------------------|
| visited nodes / non-tree edges   | 389 / 767           | 1878 / 3736           |
| median ΔP_rel                    | `9.9e-4`            | `0.29`                |
| p90 / p99 ΔP_rel                 | `0.996` / `1.00`    | `0.999` / `1.00`      |
| well / noisy / extrap / depth %  | 5.2 / 24.8 / 65.6 / 4.4 | 0.2 / 10.3 / 74.8 / 14.6 |
| in-disc edges (`extrap_max≤1`)   | 188 (24.5 %)        | 884 (23.7 %)          |
| in-disc median ΔP_rel            | `2.8e-8`            | `3.1e-6`              |
| in-disc tight fraction (`≤1e-6`) | 0.82                | 0.38                  |

The re-resolved B3 walk is **comparable, slightly worse on the raw
distribution**. The reason is geometric and honest: the B3 walk threads
a ~6× larger area (`|x|≤20` vs `|x|≤8`) with 1878 nodes, so its
Delaunay graph carries far more long edges whose midpoints sit outside
the small canonical discs (`h ≈ 0.1`) — `extrap_driven` by
construction. This is **not a regression** in the figure: the B1/B2
re-resolution made the *Stage-2 surface fill* honest and the *pole
field* accurate (VC-4/VC-5, Amendments 5–6); it was never claimed to
improve the canonical-Padé Delaunay-midpoint consensus, and the metric
shows it did not.

### VC-7's actual role — a structural quality certificate (the gate)

VC-7 is therefore wired in as a **structural quality certificate**, not
a before/after re-resolution proxy. The figure acceptance gate
(`figures/test_kkg_pi2_surface.jl` PI2S.7) asserts the certificate's
genuine, load-bearing invariants on the figure's wedge walk: the
certificate is non-degenerate (`n_edges > 100`, single sheet); it is
structurally sound (the five category tallies partition `n_edges`, the
quantiles are monotone — the identical contract `diagnose_test.jl` pins
for the scalar method); the **well-closed lobe genuinely closes to
machine precision** (`min ΔP_rel < 1e-8` — the load-bearing invariant
that proves the shared-`Q` `Pᵢ(t)/Q(t)` evaluation pipeline is correct
end-to-end); and a meaningful fraction of the genuine *in-disc* loop
closures close tightly (`> 25 %` of `extrap_max ≤ 1` edges have
`ΔP_rel ≤ tol_bad` — measured 0.38). The load-bearing assertions are
**mutation-proven** in-test: corrupting a visited node's shared-`Q`
store (scaling its numerators by `1e3`) poisons every incident Delaunay
edge and drives both the well-closed count and the in-disc tight
fraction down — confirming the gate has teeth. The measured
distribution is logged (`@info`), not masked — VC-7's value is the
certificate, and a certificate is only meaningful if its numbers are
surfaced (FW-style honest reporting, ADR Decision 1).

### Deferred — a loop-closure metric that *is* a re-resolution proxy

A metric that would genuinely reward the B1/B2 re-resolution would have
to evaluate each node's Padé only *inside* its B1 true-radius disc and
restrict the edge population to genuine in-disc loop closures — i.e.
fold the B1 gate into the diagnostic itself. That is a different
diagnostic (a "B1-gated loop-closure certificate") and is recorded as a
deferred bead under the v0.2 epic; forcing condition: a figure
requirement to certify the re-resolution *via* a loop-closure number
rather than via VC-4/VC-5 (the per-pole certificates that already do
certify it). VC-7 as shipped is the honest, senior-grade certificate of
the walk's loop-closure structure.

## Amendment 8 — Phase-D VC-8 + VC-9: the BVP / vector-solver chain certified (2026-05-21)

Beads `padetaylor-0ln.37.15` (VC-8) and `padetaylor-0ln.37.16` (VC-9)
— the last two medium-tier Phase-D validation criteria — certify the
BVP and the vector-solver chain that underpins the headline figure.

### VC-8 — BVP endpoint higher-derivative match (FW §5.2)

FW 2011 §5.2 (`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:192-193`) validates a converged
*scalar* Chebyshev BVP by estimating endpoint derivatives via the `D₁`
matrix and matching them to known values — a mismatch beyond
`~10⁻⁷·10⁻⁸` says "increase `N`". The figure's BVP is the `d = 4`
first-order vector **companion** system `y' = f(z,y)`,
`y = (u, u', u'', u''')`, whose companion form asserts the *exact* chain
identities `y[k]' = y[k+1]` (`k = 1,2,3`).

**The chosen invariant — companion-consistency.** A converged
collocation solution must satisfy `y[k]' = y[k+1]`; VC-8 computes the
spectral derivative `D₁·y[k]`, rescales it from the `t ∈ [-1,1]`
collocation variable to `z` by `1/s` (`s = (z_b−z_a)/2`, `VectorBVP`'s
affine scale), and checks it against component `y[k+1]` — at the two
endpoint nodes (FW §5.2's endpoint check) **and** over the whole node
set. This is FW's `D₁`-derivative diagnostic lifted to the companion
system: a genuine self-consistency invariant of the spectral solution
that owes nothing to the asymptotic seed and catches an under-resolved
`N`. The helper is `surf_vc8_companion_check` in
`figures/_kkg_pi2_surface_helpers.jl`; the gate is
`figures/test_kkg_pi2_surface.jl` PI2S.8.

**Measured agreement** on the figure's three BVP code paths (the
`SURF_BVP_N = 96` ray BVPs and the `N = 128` anchor BVP):

| BVP solve | companion-consistency max | y₃ vs asym. seed | y₄ vs asym. seed |
|-----------|---------------------------|------------------|------------------|
| ray φ=180° (neg. real axis) | `6.3e-10` | `1.1e-6` | `3.2e-6` |
| ray φ=120°                  | `1.9e-11` | `1.2e-6` | `1.7e-5` |
| anchor BVP `[-20,-2]`       | `1.7e-12` | `1.1e-6` | `3.2e-6` |

The companion-consistency residual is `1.7e-12 – 6.3e-10` — **two to
five orders below FW's `10⁻⁷·10⁻⁸` band**: the figure's BVPs are
comfortably resolved. The second VC-8 facet — the BVP's *free* higher
components `y₃ = u''`, `y₄ = u'''` at the deep-asymptotic endpoint
`|z_a| = 20` (only `y₁, y₂` are pinned by the `2+2` split BC) — agree
with the `pI2_tritronquee_ic` asymptotic seed to `~10⁻⁶·10⁻⁵`, the
seed's own `n_terms = 2` truncation accuracy at `|x| = 20`; this
certifies the BVP converged to the *tritronquée branch*, not merely to
*a* companion solution.

**Mutation-proven** (Rule 4): the same ray BVP solved at a starved
`N = 6` has companion-consistency residual **`13.0`** — far above the
FW band (RED) — while the figure `N = 96` solve is `6.3e-10` (GREEN).
The starved solve's *Newton residual* stays at the spectral floor
(`8e-12`) — proof that residual-norm convergence ≠ resolution, and that
only the `D₁` companion-consistency check catches the under-resolution:
exactly FW §5.2's insight. PI2S.8 asserts the figure-`N` GREEN and the
starved-`N` RED in-test.

### VC-9 — Weierstrass-℘ oracle for the vector pipeline

The v0.1 *scalar* pipeline is ℘-oracle-validated (FW §5.1.1
`u'' = 6u²`, the equianharmonic ℘ — `test/problems_test.jl`,
`test/polefield_test.jl`). The v0.2 **vector** pipeline
(`vector_solve_pade`, `vector_bvp_solve`, `vector_path_network_solve`,
`extract_poles_shared_q`) had no equivalent oracle — its figure-grade
exercise is the `P_I^(2)` tritronquée, which has no closed form. VC-9
closes that gap: it formulates `u'' = 6u²` as the `d = 2` companion
system `y' = (y₂, 6y₁²)` and drives it through the vector pipeline. The
new testset is `test/vector_pipeline_oracle_test.jl` (`VPO.*`),
registered in `test/runtests.jl`.

**Measured oracle results vs FW Table 5.1 / the closed-form ℘ lattice:**

- **VPO.1** — `vector_solve_pade` of the companion system on the
  pole-free segment `[0, 0.5]` reproduces the FW closed-form `u(0.5) =
  4.0044646690030875` to `rtol 1.8e-15`. The derivative component
  `u'(0.5)` is carried by the single `h = 0.5` segment's (15,15)-Padé
  and lands at `~4e-4` (the genuine single-segment vector Padé accuracy
  for the derivative on this steep ℘ segment, confirmed converging
  under finer `h`) — the tight `u'` oracle is VPO.3's job.
- **VPO.2** — `vector_path_network_solve`, tunnelling ~21 ℘ lattice
  poles along the real axis (a straight-line IVP step lands on a pole
  and degenerates the shared-`Q` jet — this is *why* the path-network
  exists), reaches the FW Table 5.1 **medium-range** value: measured
  `u(30) = 1.0950944` vs the reference `1.095098255959744` — error
  **`6.2e-6`** at walk step `h = 0.15` (264 visited nodes). The
  **long-range** `u(10⁴)` is *not* asserted: the probe measured it
  diverging (error `~5` after 58 705 tunnelled-pole nodes — the
  accumulated shared-`Q` step error compounds over thousands of poles)
  and the walk is impractically slow for a test; the achievable vector
  pole-tunnelling frontier is the medium range, recorded honestly here.
- **VPO.3** — `vector_bvp_solve` of the companion system on `[0, 0.5]`,
  given only the closed-form ℘ *values* `u(0), u(0.5)` as 2-point
  Dirichlet data, recovers the genuine ℘ *derivative* `u'(0) =
  1.710337353176786` to **`~1e-15`** — the FW §5.2 "the BVP solver
  provides `u'`" diagnostic, and a genuine oracle (the BVP is never
  told `u'(0)`). The BVP solution tracks the independent IVP ℘
  trajectory to `~1e-6`.
- **VPO.4** — `extract_poles_shared_q` on a short path-network walk
  past the nearest ℘ pole recovers the closed-form lattice pole `z = 1`
  (FW md:297) to **`~1e-6`** — the FW Fig 4.7 pole-location capability
  lifted to the shared-`Q` vector extractor, the same bar
  `polefield_test.jl` PF.1.1 pins for the scalar pipeline. (A *dense*
  2D ℘-lattice extraction is *not* claimed: the probe measured the
  vector wedge walk threading the dense doubly-periodic lattice
  erratically — worst-in-box `~1.3` — a genuine limit of a pipeline
  designed for the *sparse* `P_I^(2)` wedge, recorded honestly.)

**Mutation-proven** (Rule 4): VPO.5 perturbs the FW IC `u'(0)` by
`1e-3` — a different ℘ translate — and confirms the perturbed solve no
longer matches `u(0.5)` while the genuine IC does (RED/GREEN
discriminator in-test). Out-of-test: M1 (perturb the `u_at_30` oracle
to `1.10`) turns VPO.2 RED; M3 (corrupt the companion RHS to `5u²`)
turns VPO.1/2/3/4 RED — both verified, both reverted. The `1e-4`
tolerance on `u(30)` is the honest envelope of the vector
pole-tunnelling accuracy (tightening it to `1e-8` would be RED — the
measured error is `6.2e-6`).

### Phase-D medium tier complete

VC-4, VC-5, VC-7, VC-8, VC-9 are all shipped; the BVP and the
vector-solver chain are now certified end-to-end against an independent
oracle. VC-10 (two-run pole-disagreement) and VC-12 (7-fold symmetry,
stretch) remain.

## Amendment 9 — Phase-D VC-10: the FW/FFW two-run accuracy indicator (2026-05-21)

Bead `padetaylor-0ln.37.17` (Phase D, "VC-10") — the last Phase-D
validation criterion — implements the FW/FFW double-run pole-field
accuracy indicator for the headline figure's wedge.

### The mechanism — a controllable Stage-1 target ordering

FW 2011 §3.1 (line 156) shuffles the Stage-1 targets into a *random*
order before the path-network walk; FFW 2017
(`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:246-247`) makes that
randomness a *diagnostic*: "if we compute the same solution twice,
different paths will be run, resulting in solutions that should differ
by approximately the numerical error". The cross-ordering disagreement
is a practical, oracle-free accuracy estimate — exactly the
manufactured FW-style criterion ADR Decision 1 commits to (the
tritronquée wedge field has no external oracle).

`vector_path_network_solve` gains an additive `rng::Union{Nothing,
AbstractRNG}` kwarg (`src/VectorPathNetwork.jl`), mirroring the v0.1
scalar `PathNetwork`'s `rng_seed`/`shuffle` idiom (`src/PathNetwork.jl:
480-481`). With `rng === nothing` (the default) the targets are walked
in the *given* order — the V7/V8 behaviour, **bit-identical** (ADR-0001
additive; `test/vector_path_network_test.jl` VPN.4.1 pins it, and
v0.1 + v0.2 scalar suites are byte-unaffected). With an `AbstractRNG`
the targets are `shuffle`d by it: a different processing order builds a
different path-network tree, a different node filament, a different
independently-extracted pole field. The walk for a **fixed** `rng`
state stays fully deterministic (VPN.4.2; `test/kkg_pi2_figure_test.jl`
PI2F.1.6's exact-determinism guarantee is **per-fixed-ordering** —
unchanged, the figure passes no `rng`).

### Measured two-run disagreement

The figure helper `surf_vc10_two_run` (`figures/_kkg_pi2_surface_
helpers.jl`) runs the wedge walk twice with two fixed `MersenneTwister`
seeds, VC-4-validates each pole field, and matches the two fields
nearest-neighbour by a maximum-cardinality bipartite matching (reusing
the VC-5 `_vc5_augment!` augmenting-path core, `_kkg_pi2_vc45.jl`).
Measured on the production B3 wedge fan:

| quantity | value |
|----------|-------|
| run A / run B pole count | 240 / 227 (2102 / 2142 nodes) |
| matched pole pairs | 181 |
| matched-pole disagreement — **median / max** | **0.348 / 0.590** |
| single-run-only poles | 59 (A) + 46 (B) |

The matched-pole disagreement **median 0.348** is the VC-10 accuracy
indicator: the wedge pole field is accurate to roughly **half the pole
nearest-neighbour spacing** (~0.69, the VC-5 measurement). This is
consistent with — and an *independent* confirmation of — the VC-5b
finding (Amendment 6): the field carries an intrinsic ~0.3 median
conjugate-residual accuracy from the A2 far-wedge tractability ceiling.
VC-10 corroborates it from a *different* direction — two differently-
ordered walks, not the conjugate symmetry — so the ~0.35 figure is the
genuine accuracy of the field, not an artefact of either diagnostic.

The 59 + 46 single-run-only poles are the second FFW facet: each
ordering's adaptive-`h` history threads a slightly different subset of
the dense far-wedge pole lattice (the same A2 ceiling). A *majority* of
each field still matches (181 of 240 / 227) — the two runs genuinely
agree, they do not produce disjoint fields.

### Reconciliation with PI2F.1.6

`PI2F.1.6` asserts the V8b figure walk is exactly reproducible; it
**stays true** — the figure passes no `rng`, so it walks the fixed
default order, exactly reproducible. VC-10 adds the complementary fact:
a *different* `rng` ordering converges to the *same* pole field within
the reported ~0.35 accuracy. A one-line clarifying comment was added to
`PI2F.1.6` recording that determinism is per-fixed-ordering. All 27
PI2F.* assertions stay GREEN; the v0.1 scalar `pathnetwork_test.jl`
(107/107) is bit-identical.

### Mutation-proven

`figures/test_kkg_pi2_surface.jl` PI2S.10 asserts the indicator
(`median_disagree < 0.45`, `< 0.69` spacing; matched pairs ≥ 100;
matched + only-this-run partitions each field; per-fixed-seed
determinism; a same-seed match has *zero* disagreement). The
load-bearing mutation (M7, `test/vector_path_network_test.jl`): making
the `rng` `shuffle` inert turns the two seeds into identical trees —
`median_disagree → 0.0`, every PI2S.10 disagreement assertion and
VPN.4.3's different-tree assertion go RED; VPN.4.1 / VPN.4.2 (additive
default, fixed-seed determinism) stay GREEN, proving the shuffle is
genuinely load-bearing for VC-10 specifically. Restored to GREEN.

### Phase D complete

VC-4, VC-5, VC-7, VC-8, VC-9, VC-10 are all shipped — the entire
Phase-D validation suite is in place. Only VC-12 (7-fold rotational
symmetry, a stretch criterion) remains.

## Amendment 10 — Phase-C sector voters: C2 BVP-sourced arc data + C3 voter-1 reconstruction (2026-05-21)

Beads `padetaylor-0ln.37.10` (C2) and `padetaylor-0ln.37.11` (C3) — the
Phase-C sector-voter improvements — retire v1 corners **C1** (the
inner-arc asymptotic Dirichlet datum) and **C5** (the bilinear ray-fan
voter). The Region-1 triple-method pipeline was extracted into the
Rule-6 sibling helper `figures/_kkg_pi2_sector.jl` (mirroring the
VC-4/5/8/10 splits); the ADR-0024 median vote and the two Laplace
voters are unchanged — C2/C3 touch only the *voters' inputs* and the
voter-1 reconstruction.

### C2 — the inner arc is a BVP-interior datum, not the asymptotic seed

The recon proposed "source voters 2/3's inner/outer arc data from the
ray-fan BVP samples" (the ledger's original C2 line). The
C2-exploration probe (`external/probes/c2-arc-data/probe.jl`)
**falsified the literal form of that instruction**: the production ray
BVP's `2+2` split boundary condition *pins* `(u, u')` at the inner
endpoint `z_b = R_MIN·e^{iφ}` to the `n_terms = 2` asymptotic seed
**exactly** — the probe measured `|prod(R_MIN) − seed(R_MIN)| = 0.0`.
So `fan.U[1,:]` *is* the seed; re-routing voters 2/3 to read it would
collapse the inner-arc vote spread (all three voters then agree on the
same value) **without improving accuracy** — a gamed metric, not a fix
(CLAUDE.md Rule 2).

The genuine fix the probe found and C2 ships: each sector ray BVP is
solved over an **extended-inward segment** `[R_MAX, R_INNER_BC]` with
`R_INNER_BC = 1.05 < R_MIN = 2`, so the inner arc `|x| = R_MIN` is an
*interior* collocation point. Its value is then the genuine
ODE-constrained solution — the asymptotic-seed error is confined to the
`[R_INNER_BC, R_MIN]` boundary layer below the rendered window.
Measured (probe, vs a high-`N` deep-inward dedicated-BVP oracle
converged to `~10⁻¹⁴`):

| inner-arc datum | mean error vs oracle |
|-----------------|----------------------|
| `n_terms = 2` asymptotic seed | `~5·10⁻⁴` |
| C2 BVP-interior (`R_INNER_BC = 1.2`) | `~4.9·10⁻⁴` (barely better — `R_INNER_BC` must reach near `|x| = 1`) |
| **C2 BVP-interior (`R_INNER_BC = 1.05`)** | **`~1·10⁻⁴`** (a genuine ~5× gain) |

The probe confirmed all 48 fan-ray BVPs converge with the extended
segment, and that the **outer** arc needs no change — extending outward
moves `u(R_MAX)` only `~2·10⁻⁸` (the seed at `|x| = 20` is already that
accurate). `surf_laplace_voters` now sources *all four* `w`-rectangle
edges from the ray fan, removing the asymptotic series from the
voter-2/3 boundary path entirely (C1 fully retired). The wider
`[20, 1.05]` segment has more inner-end curvature; **`SURF_BVP_N`
was raised `96 → 128`** — the probe measured the BVP companion-
consistency (VC-8's metric) at `~10⁻⁷` (marginal vs the FW `10⁻⁷·10⁻⁸`
band) for `N = 96` on `[20, 1.05]` but `~10⁻¹¹` at `N = 128`; PI2S.8
confirms `9.9·10⁻¹²` on the production figure.

### C3 — voter-1 reconstruction: exact-radius + cubic-angular

The C3 audition (`external/probes/c2-arc-data/c3_audition.jl`) measured
the bilinear polar-raster voter-1 reconstruction's interpolation error
against dedicated through-the-point BVP oracles:

| voter-1 scheme | inter-ray error (median / max) | on-ray radial error (median / max) |
|----------------|--------------------------------|------------------------------------|
| **bilinear (the v1 corner C5 scheme)** | `7.9·10⁻⁴` / `1.3·10⁻³` | `3.6·10⁻⁴` / `1.0·10⁻³` |
| exact-radius + linear-angular | `5.9·10⁻⁴` / `7.3·10⁻⁴` | `0.0` / `2.4·10⁻¹⁵` |
| **exact-radius + cubic-angular (shipped)** | `1.5·10⁻⁷` / `6.9·10⁻⁵` | `0.0` / `2.4·10⁻¹⁵` |

The bilinear error (`1.3·10⁻³` max) is **larger than the bulk
triple-method spread** (`~5.5·10⁻⁴`, worklog 057) — at the production
`6°` fan it was the least-accurate of the three voters, feeding its
interpolation error straight into the vote. C5 is therefore genuinely
material; "retain bilinear" is *not* a valid audition outcome.

**Audition decision.** A *harmonic* voter 1 — a Laplace solve on the
`w`-rectangle with the ray fan's own edge data — would be **`≡ voter 2`**
(the in-house 2D-Chebyshev `laplace2d_solve`) and collapse the
triple-method independence ADR-0024 requires (a FEM voter, a spectral
voter, and a *third* that must stay disjoint). It is **rejected**.
Plain bilinear is **also rejected** (materially too lossy). The shipped
reconstruction keeps voter 1 a **direct ODE solve** (independent of both
Laplace voters) and only fixes the polar-raster → Cartesian sampling:
exact in the radius (the stored ray `VectorBVPSolution`s are barycentric
spectral interpolants, exact at any radius — the `40`-row raster and its
`~3.6·10⁻⁴` radial error are gone) and a Catmull-Rom cubic in the angle
across four bracketing rays (inter-ray error `1.3·10⁻³ → 6.9·10⁻⁵`,
~19×). In the two **boundary** angular cells the four-ray stencil is
degenerate (no genuine flanking ray on the Stokes-line side); there the
blend falls back to **linear** — a clamped-stencil cubic *extrapolates*
in the steep, `6°`-under-sampled near-Stokes field and was measured to
swing voter 1 `~8·10⁻³` off the bracketing-ray bracket. The three
voters stay algorithmically disjoint.

### Measured before/after — the sector spread

The triple-method `spread` map (max pairwise voter disagreement), old
pipeline vs C2+C3, over the rendered sector:

| region | OLD median | C2+C3 median | OLD max | C2+C3 max |
|--------|-----------:|-------------:|--------:|----------:|
| inner arc (`|x| ≤ 3`) | `1.8·10⁻³` | **`1.4·10⁻⁴`** (~13×) | `2.9·10⁻³` | `3.9·10⁻³` |
| bulk (`|x| ≥ 5`) | `5.4·10⁻⁴` | **`8.7·10⁻⁵`** (~6×) | `1.3·10⁻³` | `1.0·10⁻³` |

The **median** spread improves dramatically everywhere — C2/C3 make the
voters genuinely more accurate, not merely more agreeable. The
inner-arc **max** spread rises slightly (`2.9·10⁻³ → 3.9·10⁻³`). This
is **honest, expected, and in-scope**: an independent high-`N` BVP
oracle confirms that at the worst inner-arc point voter 1 is now
*accurate* (error `3.7·10⁻⁴`) while voters 2/3 — the Laplace solves —
carry a `~5·10⁻³` error near the negative-real-axis inner arc. The old
pipeline *masked* that error because its seed-contaminated bilinear
voter 1 was *also* wrong there and agreed with voters 2/3; fixing voter
1 un-masks the voter-2/3 defect, exactly as worklog-057 lesson 57
predicts a majority vote should ("a vote earns its keep by exposing
disagreement"). The voter-2/3 inner-arc inaccuracy traces to the
**`w`-rectangle Laplace resolution** near the steep negative-axis inner
arc (`SURF_LAP_NX = 40` radial nodes) — that is bead `0ln.37.9` (C1,
"high-resolution sector — finer grid, denser rays, higher Laplace N"),
a separate scheduled follow-up; C2/C3's mandate is the voter *inputs*
and the voter-1 reconstruction, not the Laplace `N`. PI2S.2's
`max < 1e-2` bound still holds (`3.9·10⁻³`).

### Mutation-proven (Rule 4)

`figures/test_kkg_pi2_surface.jl` PI2S.11 pins both decisions; PI2S.2
gains the C2 inner-arc-spread report. Both mutation-proven out-of-test:

- **C2** — reverting `surf_ray_bvp`'s `r_in` default to `SURF_R_MIN`
  (the legacy `[R_MAX, R_MIN]` pinned segment) turns **9 PI2S.11
  assertions RED**: the per-angle `err_fan < err_seed` oracle checks,
  the segment-extent checks, and the mutation-proof discriminator. It
  confirms the central probe finding — the legacy inner endpoint *is*
  the asymptotic seed exactly, and the spread-only metric does **not**
  catch it (PI2S.2's `inner_med` stays low under the mutation, because
  all three voters then agree on the same wrong value). The genuine,
  load-bearing C2 invariant is PI2S.11's *oracle* comparison, not the
  spread.
- **C3** — snapping the query radius to the `40`-row raster (the
  retired bilinear scheme's radial quantization) turns **4 PI2S.11
  assertions RED**: the exact-radius bar (`c3_radial_max` `0.0 →
  2.5·10⁻²`), the inter-ray bar (`3.8·10⁻²`), and the two
  bilinear-comparison assertions. The exact-radius reconstruction is
  genuinely load-bearing.

Both reverted; `figures/test_kkg_pi2_surface.jl` is GREEN at
**24196 / 24196** (PI2S.8 confirms `N = 128` companion-consistency
`9.9·10⁻¹²`; the figure is unaffected elsewhere).

### Phase C — C2 + C3 complete

C2 and C3 are shipped. Phase C's remaining bead is `0ln.37.9` (C1 — the
high-resolution sector: finer grid, denser rays, higher Laplace `N`),
which will additionally close the residual voter-2/3 inner-arc
under-resolution this amendment's measurement surfaced.

## Amendment 11 — C1 sector re-resolution: the lift ships, the inner-arc hypothesis is falsified (2026-05-22)

Bead `padetaylor-0ln.37.9` (C1 — the high-resolution sector) retires
v1 corner **C6** (`121²` grid / `6°` ray fan). The bead carried *two*
coupled mandates: (a) genuinely sharpen the headline figure raster,
and (b) "additionally close the residual voter-2/3 inner-arc
under-resolution" that Amendment 10 §"Measured before/after" attributed
to the `w`-rectangle Laplace resolution. Mandate (a) ships clean.
Mandate (b)'s **premise was falsified** — CLAUDE.md Rule 2, root cause
over band-aid — and that falsification is the substantive content of
this amendment.

### The lift — C6 retired

`SURF_GRID_N` `121 → 401` (Cartesian raster) and `SURF_RAY_DPHI_DEG`
`6° → 2°` (≈ 141 fan rays). Both stay odd / node-aligned: the `401`
raster lands a node exactly on `x = 0` and on `x = -15` (the PI2S.1
pin). The lift is mostly cheap: voters 2/3 are spectrally callable at
any point, and post-C3 voter 1 is exact-in-radius (barycentric) +
cubic-angular, so the finer raster is `O(N²)` interpolant evaluation,
not `O(N²)` ODE solves. The `2°` fan was measured the runtime knee:
the cubic-angular `surf_ray_eval` inner-arc reconstruction error falls
`~5·10⁻⁴` (`6°`) → `~2·10⁻⁸` (`2°`); a `1°` fan was probed and
**rejected** — ~8 s extra (one BVP per ray) for *zero* measurable gain,
because the residual error is the seed floor below, not the fan.

### The falsification — the inner-arc spread is NOT a Laplace under-resolution

Amendment 10 hypothesised the inner-arc *max* spread (`3.9·10⁻³` on the
`121²` grid) was a voter-2/3 `w`-rectangle Laplace **under-resolution**
near the steep negative-real-axis inner arc, and scheduled the fix to
this bead ("higher Laplace `N`"). The C1 convergence probe **falsified
that hypothesis** on three independent measurements:

1. **The spectral voter does not converge toward the truth with `N`.**
   Sweeping `SURF_LAP_NX` `40 → 120` does not move the inner-arc
   interior value — it stays pinned `~4.8·10⁻³` off a dedicated-BVP
   oracle. A genuine under-resolution would shrink geometrically.

2. **Perfect edge data leaves the same error.** Re-solving the
   spectral `w`-rectangle problem with *exact* dedicated-BVP edge
   traces (eliminating any ray-fan edge-data error) leaves the same
   `~4.8·10⁻³` interior error — `N`-independent. The defect is not in
   the harmonic solve.

3. **Raising `N` past the knee makes the spectral solve *worse*.** The
   seed-floored edge data is *slightly non-harmonic* (the asymptotic
   seed does not satisfy Laplace's equation exactly); a higher-`N`
   spectral solve *over-resolves* that non-harmonic contamination —
   spectral self-error rises `2.2·10⁻⁵` (`Nx = 48`) → `3.1·10⁻⁴`
   (`Nx = 64`). So a speculative `N` raise is not merely useless, it is
   actively counter-productive. `SURF_LAP_NX/NY` are therefore set to
   the measured spectral *convergence knee* `48/64` (cleanest), not
   raised; `SURF_FEM_NX/NY` stay `64/88` (the FEM voter has plateaued
   by `n_x = 64`).

### The genuine root cause — the asymptotic-seed truncation floor

Each sector ray BVP pins `(u, u')` at the inner endpoint
`R_INNER_BC = 1.05` to the KKG `n_terms = 2` asymptotic series, whose
`O(|x|^{-7/3})` truncation error at `|x| ≈ 1` propagates inward as a
boundary layer. Solving the *same* ray at several legitimate
seed-pinned inner radii `r_in ∈ [1.05, 1.8]` gives interior values
that disagree by:

| `|x|` | `r_in`-spread (the seed floor) |
|-------|-------------------------------|
| 2.0   | `~5·10⁻³` |
| 2.5   | `~4·10⁻³` |
| 5.0   | `~7·10⁻⁴` |
| 10.0  | `~1.7·10⁻⁵` |

This `r_in`-spread *is* the irreducible floor — it is `N`-independent
and fan-density-independent. `pI2_tritronquee_ic` implements only
`n_terms ∈ {1, 2}`; the seed cannot be sharpened in-place. Moving
`R_INNER_BC` *does* shrink the v1–v2 spread, but only by making the
voters agree on the *same* seed-floored value (the gameable-spread
anti-pattern C2's probe already identified) — not by improving
accuracy. `R_INNER_BC` is therefore left at the Amendment-10 /
PI2S.11-certified `1.05`.

### Measured inner-arc spread — before / after

| grid | inner-arc (`|x| ≤ 3`) median | inner-arc *max* | note |
|------|------------------------------:|----------------:|------|
| `121²` (Amendment 10) | `1.4·10⁻⁴` | `3.9·10⁻³` | coarse sampling — the worst seed-floor zone at `|x| = 2` was under-sampled |
| **`401²` (C1, this amendment)** | **`1.6·10⁻⁴`** | **`6.7·10⁻³`** | the seed floor fully exposed — the `401²` raster samples grid points right at the worst `|x| = 2` zone |

The `121²` `3.9·10⁻³` was never a real lower floor — it was a
coarse-sampling artefact. The `401²` grid does not *worsen* the
figure; it *honestly exposes* the pre-existing seed floor. The
`|x| ∈ [2.6, 3.0]` band max is `5.9·10⁻³ < 6.7·10⁻³` — the floor
**decays outward**, confirming it is a seed boundary layer and not a
global Laplace defect (a global defect would not decay with radius).
The whole-sector spread median is `8.1·10⁻⁵` — the bulk concurs
tightly. Closing the inner-arc floor needs an `n_terms ≥ 3` seed
(KKG 2015 eq. 7.2 `c₇…`), unimplemented — recorded as the deferred
bead below.

### Chosen resolutions + runtime

`SURF_GRID_N = 401`, `SURF_RAY_DPHI_DEG = 2.0`, `SURF_LAP_NX = 48`,
`SURF_LAP_NY = 64`, `SURF_FEM_NX = 64`, `SURF_FEM_NY = 88`,
`SURF_BVP_N = 128` (already converged — `N = 128` agrees with an
`N = 320` ray solve to `~10⁻¹²`). The full `figures/test_kkg_pi2_surface.jl`
suite is **GREEN at 243738 / 243738**, figure testset wall time
**`2m22s`** — comfortably inside the Phase-F runtime gate (well under
`~10 min`).

### Mutation-proven (Rule 4)

PI2S.2's C1 block pins `SURF_GRID_N == 401`, `SURF_RAY_DPHI_DEG == 2.0`,
the `401²` raster sizes + the odd-grid node-alignment invariant, and —
the load-bearing assertion — the inner-arc *max* spread bound
`inner_max < 8e-3` (the measured `6.7·10⁻³` seed floor + headroom) plus
the outward-decay check `maximum(outer_band) < inner_max`. Tightening
`inner_max`'s bound to the *falsified* Amendment-10 target (`3.9·10⁻³`)
turns **exactly one assertion RED** (line 259, the `inner_max` bound;
`243737 / 243738`) — the measured `6.7·10⁻³` genuinely binds the
assertion; it is not slack. The bound was restored and the suite is
GREEN at `243738 / 243738`.

### Deferred — the `n_terms ≥ 3` seed

The inner-arc seed-floor is closed only by an `n_terms ≥ 3` KKG
asymptotic seed (eq. 7.2 `c₇…` and beyond). `pI2_tritronquee_ic`
implements `n_terms ∈ {1, 2}` only. **Deferred bead** (forcing
condition): file/raise the `n_terms = 3` seed work *if and only if* a
downstream consumer needs inner-arc surface accuracy below the
`~5·10⁻³` seed floor at `|x| ≈ 2`. For the headline figure (KKG
Figs 7.4/7.5, "for visualization") the floor is sub-visible and
in-scope; the figure does not force the v2 seed.

### Phase C complete

C2, C3, C1 are all shipped; v1 corners C1 (Amendment 10), C5
(Amendment 10), C6 (this amendment) retired. Phase C is closed.

## Amendment 12 — Phase-E1 C4 audition: the Stokes-strip mask narrowed `±3° → ±1°` (2026-05-22)

Bead `padetaylor-0ln.37.18` (Phase E1) auditions the last open v1
corner — **C4**, the `±3°`-wide NaN mask straddling each `±36°` Stokes
line (`surf_in_mask`). The mask exists because the sector ray fan is
inset `SURF_SECTOR_MARGIN_DEG = 4°` from the Stokes line (the sector
solver covers `arg x > 40°`) and the wedge B1-gated underlay is sparse
below `|arg| ≈ 34°` — so the band `arg x ∈ (34°, 40°)` is covered by
neither region and renders grey. The audition question (Rule 9): is the
`4°` inset + `3°` mask a genuine numerical limit, or a conservative
cushion?

### The audition — evidence (`external/probes/stokes-strip-audition/`)

The probe rebuilds the production sector pipeline (`surf_ray_fan` +
`surf_laplace_voters` — the identical BVP recipe, the identical two
Laplace solvers, the C3 reconstruction) parametrised on the fan margin,
sweeping `margin ∈ {4°, 3°, 2°, 1°, 0.5°}`. For each it measures, in
1°-wide angular bands marching toward the 36° Stokes line: the
triple-method `surf_vote` spread, and — the independent oracle — the
voted value's error against a dedicated through-the-point BVP straight
along the query ray (the same oracle idiom as the PI2S.11 C2/C3
audition). Two findings:

**1. The sector ODE solve does NOT degrade at the Stokes line.** A
dedicated through-the-point BVP walked from `arg x = 30°` to `36°`:

| `arg x` | Newton residual | companion-consistency |
|---------|-----------------|------------------------|
| 30.0°   | `1.8·10⁻¹¹`     | `4.2·10⁻¹²`            |
| 33.0°   | `1.3·10⁻¹¹`     | `4.4·10⁻¹²`            |
| 35.0°   | `1.6·10⁻¹¹`     | `1.1·10⁻¹¹`            |
| 35.9°   | `1.8·10⁻¹¹`     | `1.3·10⁻¹⁰`            |
| **36.0°** | **`1.7·10⁻¹¹`** | **`1.8·10⁻¹⁰`**      |

The BVP is as healthy *at* the Stokes line as mid-sector — no blow-up.
The pole-free sector solution is analytic (harmonic) right up to the
line; the `4°` inset was a cushion, not a numerical limit. (This is
expected: the Stokes line is a feature of the *asymptotics* — where the
recessive exponential of the connection problem switches dominance —
not a singularity of `V₀`, which is analytic on the whole pole-free
sector.)

**2. The triple-method vote stays honest down to a `1°` margin.** The
near-Stokes vote spread + oracle error degrade smoothly toward 36° but
stay bounded. Measured (upper Stokes line; the lower is the mirror):

| margin | band toward 36° | vote spread med / max | oracle err med / max |
|--------|-----------------|----------------------|----------------------|
| 4° | `[40°,41°]` (fan edge) | `8.7·10⁻⁵` / `3.0·10⁻⁴` | `6.3·10⁻⁵` / `1.9·10⁻⁴` |
| 1° | `[39°,40°]` | `5.2·10⁻⁴` / `1.0·10⁻³` | `1.2·10⁻⁴` / `6.8·10⁻⁴` |
| 1° | `[38°,39°]` | `9.8·10⁻⁴` / `2.1·10⁻³` | `2.9·10⁻⁴` / `6.9·10⁻⁴` |
| 1° | `[37°,38°]` (fan edge) | `1.8·10⁻³` / `3.3·10⁻³` | `6.6·10⁻⁴` / `2.2·10⁻³` |
| 0.5° | `[36.5°,37°]` (fan edge) | `1.9·10⁻³` / `4.3·10⁻³` | `1.3·10⁻³` / `4.1·10⁻³` |

The figure's bulk honesty envelope (PI2S.2) is vote spread median
`< 3·10⁻³`, max `< 1·10⁻²`, and it already accepts an inner-arc
seed-floor max of `~6.7·10⁻³` (Amendment 11) as honest. At `margin = 1°`
the worst near-Stokes band `[37°,38°]` sits at spread max `3.3·10⁻³`,
oracle error max `2.2·10⁻³` — comfortably inside that envelope and
below the inner-arc floor. At `margin = 0.5°` the `[36.5°,37°]` band
oracle error climbs to max `4.1·10⁻³` (still inside the `< 1·10⁻²`
envelope, but the C3 voter-1 boundary-cell linear fallback is steeply
`2°`-under-sampled there and the error is climbing fast).

### Decision — narrow the mask `±3° → ±1°` (fill 3× of the grey strip)

`margin = 1°` is the honest frontier; `margin = 0.5°` is rejected (its
oracle error nearly doubles to `4.1·10⁻³` — the senior-grade call keeps
`1°` of headroom against the boundary-cell under-sampling). So:

- `SURF_SECTOR_MARGIN_DEG`: **`4° → 1°`** — the sector solver now
  covers `arg x > 37°`;
- `SURF_STITCH_MASK_DEG`: **`3° → 1°`** — the mask is `[35°, 37°]`.

The mask half-width is set **equal to the fan margin** — the tightest
honest choice: the mask covers exactly the band the sector fan does not
reach (`arg x ∈ (36°, 36°+margin)`), so no in-sector cell the fan
misses escapes the mask, and no cell the fan reaches is masked. The
headline figure's grey Stokes strip shrinks from `6°` wide
(`[33°,39°]`) to `2°` (`[35°,37°]`) — a **3× reduction** — and every
newly-filled cell is the harmonic triple-method vote, backed by the
audition's measured spread / oracle-error evidence (no Padé out of
disc, vote spread controlled).

### Retained — the residual `±1°` mask, with a forcing condition

A *full* fill (mask `→ 0`) is **not** senior-grade-reachable with the
two-region architecture. The sector ray fan is a *sector* solver — its
BVP rays are pinned by the `pI2_tritronquee_ic` asymptotic seed and the
C3 four-ray Catmull-Rom stencil degenerates in the boundary angular
cell — it cannot straddle the 36° Stokes line itself; and the wedge
region is disjoint, with the B1-gated underlay sparse near `|arg| ≈
36°`. The residual `[35°, 37°]` strip is retained.

**Forcing condition for the residual `±1°` mask:** a uniform connection
formula across the Stokes line (the asymptotics that switch dominance
*at* the line) would let the sector and wedge solutions be stitched
continuously through it, closing the last `2°`. That is a genuine new
piece of mathematics (the KKG 2015 connection problem), not a tuning
knob — recorded as a **deferred bead** under the v0.2 epic; forcing
condition: a figure requirement for a gapless Stokes-line stitch.

### Mutation-proven (Rule 4)

`figures/test_kkg_pi2_surface.jl` PI2S.12 pins the audition outcome:
the narrowed parameters (`SURF_SECTOR_MARGIN_DEG == 1.0`,
`SURF_STITCH_MASK_DEG == 1.0`, half-width `==` margin), the
`surf_in_mask` geometry (`[35°,37°]` masked, `≥37°` not), and the
load-bearing deliverable — the band `arg x ∈ (37.5°,39°)`, left grey by
the shipped `4°/3°` figure, is now `98.9 %` filled (1036 of 1048 cells)
by the harmonic triple-method vote, with band vote spread median
`9.9·10⁻⁴` / max `4.3·10⁻³` (inside the bulk envelope). Reverting
`SURF_SECTOR_MARGIN_DEG`/`SURF_STITCH_MASK_DEG` to the retired `4°/3°`
turns PI2S.12 **RED** (16 failed + 1 errored): the newly-filled band
collapses to all-`NaN`, the parameter pins and the `surf_in_mask`
geometry checks all fail. The narrowing is genuinely load-bearing.
With the narrowed values the full `figures/test_kkg_pi2_surface.jl`
suite is **GREEN at 249344 / 249344** (figure testset wall time
`2m14s`).

### C4 ledger row — updated

C4 is **Retired-or-justified**: the mask is narrowed `±3° → ±1°` (3×
of the grey strip filled with audited-honest harmonic vote), and the
residual `±1°` mask is retained with the explicit connection-formula
forcing condition above. The v1-Corner Retirement Ledger row is updated
accordingly. With C4 dispositioned, all six v1 corners are retired or
rigorously justified.

## Amendment 13 — honest reassessment: the figure is inadequate; the scope was wrong (2026-05-22)

On inspection of the rendered figure the user judged it inadequate, and
they are right. This amendment records the correction; worklog 059 is
the full diagnosis.

**The figure does not visualise the solution in the pole wedge.** It
shows a sparse ~266-pole scatter over a ~95 %-blank region. Phases A–F
of this ADR are all genuine — the v1 corners *are* retired, the
validation suite *is* built and mutation-proven, the suite *is*
5326/5326 GREEN — but they were applied to a structurally inadequate
substrate.

**The root cause is mis-scoped, not a bug.** `src/VectorPathNetwork.jl`
is a deliberately **minimal ~156-LOC port** of v0.1's full **1108-LOC**
`src/PathNetwork.jl` driver — its own docstring says so. The FW 2011
§3.1 **5-direction wedge** and the **Stage-1 tree** were ported; the
machinery that actually *fills a region* was not: a **resilient** walk
(the minimal one `throw`-aborts the entire run on the first unreachable
target), a **fine-grid** target set (the figure runs 171 coarse
targets, not a tiling grid), and the **barycentric** Stage-2 fill. With
pole-density-forced `h≈0.1` discs the wedge (area ~250) needs ~25 000
tiling nodes; the minimal brittle walk places ~2 000, hence ~5 %.

**ADR-0025's scope was wrong.** It scoped "retire the v1 corners +
certify the result" *taking the vector path-network as given*. The
prerequisite — bringing that path-network up to FW's full driver — was
not in scope, and the Phase-A2 spike that measured the ~5–18 % coverage
ceiling was wrongly recorded as a *deferred* bead (`0ln.38`) rather than
treated as a **blocker**. That deferred bead's framing ("a 2D
re-expansion lattice, a different architecture") was also wrong: it is
not a different architecture — it is the full FW driver. `0ln.38` is
superseded by **`padetaylor-0ln.40`**, the plan-first deep-dive that
ports the full driver and re-does the figure's wedge as a genuine
filled surface.

**Reusable from this ADR**: the true-radius Stage-2 gate (B1), the
adaptive `:max_q_root` walk (B2), and the entire VC-4…VC-10 validation
suite (Phase D) are all sound and carry over to the re-done figure
unchanged. This ADR stays Accepted for that work; ADR-0026 (under
`0ln.40`) will govern the full-driver port.

## References

- KKG 2015 TeX `references/tex/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41/tritronquee_coeff.tex`
  — eq. (1.1) `:128`; V₀ sector eq. (1.4) `:224-234`; Painlevé
  property `:3124-3125`; §7 figure method + "for visualization"
  `:3203-3217`; asymptotic series `:3138-3149`.
- FW 2011 `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
  — validation-without-ground-truth `:303-311`; P_I Laurent `:59-61`;
  endpoint-derivative diagnostic `:192-193`; ℘ reference values
  `:281-301`.
- FFW 2017 `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md`
  `:120-124` (conjugate-symmetry `E(ζ)`), `:246-247` (double-run).
- ADR-0016 — `Diagnostics.quality_diagnose` loop-closure certificate
  (VC-7); ADR-0019 — shared-Q Padé; ADR-0024 — the triple-method vote.
- `docs/v0p2_pillarC_painleve_hierarchy_findings.md` §4 — tritronquée
  sector geometry + the sign-corrected IC.
- `docs/worklog/057-whole-plane-kkg-surface-figure.md` — the shipped
  figure + the v1-corner diagnosis this ADR supersedes.
- CLAUDE.md Rules 1, 2, 3, 4, 5, 6, 7, 9, 10.
