# ADR-0025 — Senior-grade re-resolution of the P_I⁽²⁾ tritronquée headline figure

**Status**: Accepted (2026-05-21) — framework + validation suite locked;
two parameters (the Stage-2 validity-gate criterion and the honest wedge
frontier) are determined by the Phase-A exploration spikes and recorded
here by amendment as they resolve.
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
| VC-4 | Dominant-balance leading coefficient `A ∈ {-1,-3}` per pole | strong | `[N]` |
| VC-5 | Conjugate-symmetry pole pairing (`V₀(x̄)=conj V₀(x)`) | strong | `[N]` |
| VC-6 | Cross-node support filter (`min_support ≥ 2`) | medium | `[E]` active |
| VC-7 | Loop-closure ΔP_rel certificate (`quality_diagnose`, ADR-0016) | medium | `[E]` unused on vector walks |
| VC-8 | BVP endpoint higher-derivative match (FW §5.2 diagnostic) | medium | `[N]` |
| VC-9 | Weierstrass-℘ oracle for the *vector* pipeline (FW Table 5.1) | medium | `[N]` |
| VC-10 | Two-run different-path pole-disagreement accuracy indicator | medium | `[N]` |
| VC-12 | 7-fold rotational-symmetry cross-check (`V₀`→`V₁` poles) | strong | `[N]` stretch |

(VC-11 Froissart residue filter and VC-13 pole-count growth law are
recorded in the recon as lower-priority; deferred with a forcing
condition rather than scheduled.)

## v1-Corner Retirement Ledger

| Corner | Disposition | Bead |
|--------|-------------|------|
| C1 inner-arc asymptotic datum | **Retired** — voters 2/3 inner/outer arc data sourced from BVP ray samples, not the series | C2 |
| C2 `extrapolate=true` | **Retired** — Lever-1 true-radius gate; no out-of-disc evaluation | B1 |
| C3 hand-tuned `h=0.1` | **Retired** — Lever-2 adaptive / shared-Q-root step | B2 |
| C4 `±3°` Stokes mask | **Audition** (Phase E) — fill via overlap, or justify as masked with a forcing condition | E1 |
| C5 bilinear ray-fan voter | **Audition** (Phase C3) — harmonic reconstruction, or justify bilinear | C3 |
| C6 `121²` grid / `6°` rays | **Retired** — Lever-3 refinement | C1 |

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
