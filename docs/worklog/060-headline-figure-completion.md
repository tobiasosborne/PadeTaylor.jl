# Worklog 060 — headline figure completion: D1–S8 + tf9 arc

**Date**: 2026-05-26
**Author**: Claude Opus 4.7
**Epic**: `padetaylor-0ln` (v0.2) · **Follows**: worklog 059, ADR-0025,
ADR-0026
**Scope**: The `padetaylor-0ln.40` arc closes — the P_I⁽²⁾ tritronquée
headline figure ships with a dual-fill provenance wedge (FW-2011-style
filled with B1-certified core at full opacity, FW-style extrapolated
overlay at reduced alpha). The investigation arc — D1, D2, D3, R1,
S1–S8, tf9 — is complete. This worklog is the sober summary.

> **Take-home**: the headline figure now renders the wedge **FW-faithful
> filled with honest provenance marking**. The certified core (~29 % of
> the wedge at order 48, B1-honest, no Padé evaluated past its verified
> disc) renders at full opacity; the FW-style extrapolated overlay (the
> other ~71 %, evaluated past the verified disc the way FW 2011's own
> published figures do — FW-md:395-397) renders at reduced alpha
> (α = 0.50). The reader sees the full surface FW shows AND reads
> provenance off the alpha channel. ADR-0025's strict no-extrapolation
> rule remains the package default for every non-headline figure; the
> dual-fill is a per-figure override, justified by FW precedent and
> honest *only because* the figure marks its own provenance.

## The arc — coverage ladder

The investigation ran from worklog 059's honest reassessment ("the
figure is ~95 % blank in the wedge") to ADR-0026 Amendment 9's
converged finding ("the cap is the honesty contract, not a walk
defect"). The certified-coverage progression:

| stage | bead | certified coverage | mode |
|-------|------|---------------------|------|
| pre-resilient (V8b, order 24, sparse 171-fan) | (the inadequate figure) | **~5 %** | ADR-0025 Amendment 2 sparse |
| S2–S5 corrected stack (order 24) | `0ln.40.a-c`, S1-S5a | **~8 %** | resilient adaptive, plateau |
| **S6a (order 24 → 36)** | ADR-0026 Amendment 6 | **~22 %** | order-knee lift |
| **tf9 (order 36 → 48)** | bead `tf9`, this work | **~29 %** | the measured saturation knee |
| FW-style overlay (no validity gate) | tf9 composite | **~100 % of wedge** | reduced alpha, FW-faithful |

The certified ceiling at order 48 is the **measured order-saturation
knee** (ADR-0026 Amendment 5: order 24 → 36 → 48 → 72 yields
13 % → 23 % → 29 % → 29 %, order 72 wasted). The B1 honest disc is
uniformly truncation-limited (4737/4737 walk nodes; the pole-adjacency
clamp never binds in the order-36 walk) — Taylor order is the only
walk-side lever on certified coverage, and it is exhausted.

Beyond the certified ceiling, only the no-extrapolation honesty contract
itself (ADR-0025) caps the *visible* coverage. The contract is stricter
than FW 2011 (FW's own published pole-field figures evaluate each
Stage-2 Padé over its full step with no validity gate). The dual-fill
provenance design (this work) renders the wedge FW-faithful while
marking the contract boundary visually.

## The three scale heresies fixed

The D1–S8 arc surfaced three scale-fixing constants that broke
scale-covariance and contributed to the coverage plateau:

1. **`h_max` (fixed step size).** ADR-0025 Amendment 4's `_adaptive_h`
   was a *geometric sink* — its pole cap `POLE_SAFETY·h_prev·min|t*|`
   was multiplicative in `h_prev`, so `h` ratcheted to zero in a
   sustained pole field. ADR-0026 Amendment 4 (S2) replaced it with the
   absolute, non-ratcheting law `h = clamp(SAFETY·D_local, h_min,
   h_max)` where `D_local = h_prev·min|t*|` is the h-independent pole-
   free disc radius. The sink is gone; a node leaving a dense pocket
   recovers to `h_max` in one step.

2. **`cluster_atol = 0.2`** (absolute pole-clustering tolerance).
   ADR-0026 Amendment 3's S4 measured this was a scale heresy that
   wrongly merged/split poles. Replaced with scale-derived clustering
   (a fraction of the local pole spacing) — `cluster_atol = nothing`
   in `extract_poles_shared_q`.

3. **`radius_t = 5`** (pole-extraction filter in the rescaled
   `t`-variable). A smaller `SAFETY` made the same physical pole map
   to a larger `|t*|`, falling outside the window — the `min_support`
   cross-node filter emptied the field. ADR-0026 Amendment 7's S7
   fixed it: `extract_poles_shared_q` now accepts a shared-`Q` root by
   its z-plane distance `h_node·|t*| ≤ radius_t·h_max` (a scale-stable,
   h-independent window). This unblocked the companion `SAFETY` 0.25 →
   0.10 lowering.

All three are now scale-derived; the system is scale-covariant.

## What ships in `tf9`

- **`SURF_PN_ORDER` 36 → 48** in
  `figures/_kkg_pi2_surface_helpers.jl` (the measured saturation knee).
- **Dual-fill `surf_wedge_fill`** — one Stage-1 walk, two Stage-2 fills.
  The package call (`vector_path_network_solve` with
  `extrapolate = false`) produces the B1-certified core; the figure
  helper then re-evaluates
  `PadeTaylor.VectorPathNetworkStage2._stage2_fill` directly on the
  cached walk with `extrapolate = true` to produce the FW-style filled
  overlay. The second fill is a Horner sweep — no new walk, no rebuild.
- **Kernel API addition**: `kkg_pi2_surface()` returns `Re_u_extrap`,
  `Im_u_extrap` alongside the existing `Re_u`, `Im_u`. The certified
  matrices preserve the "covered ⟺ finite" invariant the PI2S.4 test
  pins; the extrapolated matrices are additional.
- **Figure render composite** in
  `figures/kkg_pi2_tritronquee_surface.jl`. Two layers per panel:
  extrapolated cells (the FW-style overlay where the certified layer
  leaves NaN) drawn first at α = 0.50; certified cells drawn on top at
  full opacity. The 266-pole validated field is overlaid as black
  dots, unchanged. The figure caption explains the provenance
  distinction in plain English on the figure itself.
- **No `src/` change**. The kernel reaches into private internals
  (`_stage2_fill`, `_nearest_visited`) rather than re-exporting an
  `extrapolate = true` convenience — the legitimate scope of a
  per-figure exception.

## Documentation in lockstep (Law 2)

- **ADR-0025 Amendment 14** records the headline-figure-only
  relaxation of the no-extrapolation rule, the FW-faithfulness
  argument, the provenance-marking solution, and the measured
  certified fraction. Explicit: the strict rule remains the package
  default and applies to every non-headline figure.
- **ADR-0026 Amendment 10** closes the investigation: the figure
  renders FW-style filled with B1-certified provenance overlay; the
  D1–S8 + tf9 arc is complete; ADR-0026 closes here.
- **`surf_wedge_fill` docstring** rewritten to the dual-fill provenance
  design.
- **File / script top-of-file docstrings** updated to describe the
  certified-vs-extrapolated rendering.

## Hard-won lessons

### The negative value of premature ceiling claims

Two prior agents called "ceiling" before measurement was complete.
Each claim cost hours and was falsified by the next probe:

- **Amendment 1's "density-alone premise"** (the first ceiling claim
  — the dense-target probe regressed to 1.7 %). Falsified by S1's
  diagnostic — the regression was the adaptive-`h` geometric sink, not
  a density ceiling.
- **`cfq` "the discs are too small to tile" outer-wedge mechanism**
  (Amendment 8). Falsified by S8's instrumented measurement on the
  shipped config — the "4.4× `h_v` collapse" was not reproducible on
  the full-wedge walk; the discs are not 4.5× smaller in the outer
  wedge.
- **The orchestrator's interim "architectural ceiling" worry**
  (around Amendment 5). Falsified by the order sweep (~13 % → 22 %)
  and the `SAFETY` sweep (~22 % → 31 % on the inner wedge).

Each of these was, at the time, presented with measurements and looked
load-bearing. Each turned out to be a measurement artefact, a
mis-attributed mechanism, or a premature extrapolation from a partial
probe. The value of a `cfq`-style premature ceiling claim is *negative*:
it would freeze progress at a wrong number and motivate the wrong
remediation. Per CLAUDE.md Rule 5 ("measure & report honest numbers"),
the only valid ceiling claim is one that survives the next
investigation — and only ADR-0026 Amendment 9's "the cap is the
honesty contract" was such a finding, because it was a documented
design decision (ADR-0025), not a numerical limit.

### Two of worklog 059's diagnosed gaps did not exist

The four-agent recon (ADR-0026 Context) corrected two of worklog
059's four claims: the Stage-2 fill machinery was already ported (bead
`0ln.20`); the "barycentric Stage-2 fill" does not exist in FW (FW's
Stage-2 is plain nearest-node single-Padé, FW-md:166, :395-397). The
remediation was D1 (resilient walk) + D2 (dense targets), not a
wholesale 1108-LOC port. Worklog 059's honest framing was right; two
of its specific gap diagnoses were corrected before any code shipped.

### The investigation converged in nine amendments

ADR-0026 ran nine amendments and one closeout (Amendment 10, this
work). At each stage a measured probe falsified or refined the
previous stage's hypothesis. The arc length is real: the figure was
inadequate, the cause was non-obvious, and the converged answer (the
honesty contract is the cap) was visible only after the geometric
sink, the spurious-root suspicion, the outer-wedge regression, and the
order/SAFETY trade-offs were all individually measured and ruled out.
Senior-grade is sometimes slow.

## References

- `docs/adr/0025-headline-figure-re-resolution.md` — the re-resolution
  ADR (Amendments 1–14); Amendment 14 the dual-fill provenance.
- `docs/adr/0026-vector-resilient-walk-dense-targets.md` — the
  resilient-walk + dense-targets ADR (Amendments 1–10); Amendment 10
  the closeout.
- `docs/worklog/059-headline-figure-honest-reassessment.md` — the
  honest reassessment that set up this arc.
- `figures/_kkg_pi2_surface_helpers.jl` — `surf_wedge_fill` (the
  dual-fill helper) and `kkg_pi2_surface()` (the kernel returning
  `Re_u`/`Im_u`/`Re_u_extrap`/`Im_u_extrap`).
- `figures/kkg_pi2_tritronquee_surface.jl` — the render script (the
  two-layer alpha composite).
- `figures/output/kkg_pi2_tritronquee_surface.png` — the rendered
  figure.
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
  — `:395-397` FW's Stage-2 evaluation (no validity gate).
