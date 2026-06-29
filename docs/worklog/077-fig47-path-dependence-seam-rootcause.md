# Worklog 077 — FW 2011 Fig 4.7 spurious "seam": complete root cause (path-dependence)

**Date:** 2026-06-29. **Author:** Claude Opus 4.8 (1M).
**Beads:** `padetaylor-vwgl` (P0, root cause — the bug), `padetaylor-sny7`
(P1, the seed-invariance regression gate). Related cure beads:
`padetaylor-6b5` (K-nearest consensus Stage-2), `padetaylor-fe9`/`8dg`
(smooth-region BVP fill), `padetaylor-zwh` (denser sampling).
**Probe / reproduction:** `external/probes/fig47-seam-diagnosis/`.
**Deep-research transcript:** workflow `diagnose-field-discontinuity`
(7 agents, verified against HEAD).
**Status:** **ROOT-CAUSED ONLY. No fix applied** — by explicit user
directive ("maximum critical priority; no further progress until completely
root caused"). This worklog is the root cause; the fix is gated behind it.

> **Take-home.** `path_network_solve` enforces **no path-independence**.
> A single-valued meromorphic `u(z)` reached by two different walk
> branches must agree; ours does not — by O(1), over ~35 % of the
> [-50,50]² plane. Where two branches' territories meet, the
> reconstructed field jumps and the per-node pole lattices lose
> registration: the visible **grain-boundary seam** + spurious poles
> (worst in panels (e), (a); present in all six; in (f) even inside the
> forbidden sectors). It is **not** the Stage-2 reconstruction stitch
> (refuted by measurement). It is the WALK. It has been present since
> project start (worklog 014, "y=0 row discontinuity") and was never
> pinpointed because every lens used to look at it — global statistics,
> pole-point analysis, proxy problems, "it's just Moiré" — is
> structurally blind to it.

---

## 1. How it surfaced

A request to regenerate all FW-family figures (Fig 4.7 = six PI pole
fields, `figures/fw2011_fig_4_7.jl`) was interrupted: the maintainer
recognised in the rendered figure a spurious discontinuity that has
dogged the project from the start and was never root-caused — and that
makes "all progress brittle." The pole field should be a single coherent
Boutroux quasi-lattice; instead a **seam** runs through the smooth region
of each panel, with scattered spurious poles.

This worklog also records a methodological failure: the first several
diagnostic attempts were *wrong in kind*, and understanding why is half
the value (see §6).

## 2. The symptom, precisely

- A **grain boundary**: two pole-lattice domains of slightly different
  orientation meet at a sharp line, with a band of disordered/misregistered
  poles along it. Exposed cleanly by extracting the dot centroids from the
  figure (removing render Moiré) and computing the local bond-orientation
  field ψ₆ — a *smooth* colour gradient is physical Boutroux curvature; a
  *sharp* colour discontinuity is the seam (`orient_e_UL.png` in the
  session). Worst in **(e)** (near-tritronquée) and **(a)** (`u(0)=u'(0)=0`).
- **Spurious poles**, including in panel **(f)**'s four pole-free
  (forbidden) sectors.
- It is **in the pole scatter that the figure actually plots** (not only
  in a separate field render — a correction to the deep-research report's
  one slip): the per-node denominator-root lattices inherit the seam via
  their orientation.

## 3. Root cause — path-DEPENDENCE of the walk (measured, not argued)

`path_network_solve` builds a **tree** rooted at the IC (z=0), visiting
shuffled targets (`src/PathNetwork.jl:481`), each sub-walk chaining off
the nearest already-visited node; each node stores a Padé built **solely
from its own forward-propagated state** (`src/PathNetwork.jl:605`,
`src/PadeStepper.jl`). Two spatially-adjacent regions can be reached by
*different* branches that diverged near the IC and accumulated *different*
IVP error in the smooth/pole-free sectors — which FW md:401 calls out as
the IVP-unstable regions. **Nothing enforces that the two reconstructions
agree where they meet.** That meeting line is the seam.

This was confirmed by four mutually-reinforcing measurements
(`external/probes/fig47-seam-diagnosis/`, panel (e)):

| test | result | conclusion |
|---|---|---|
| Two-seed field difference `[-50,50]²` (rng 0 vs 42) | **35 % of cells (3523/10201) differ > 1 %, max rel = 1.000** at z≈15+35i | solution is **grossly path-dependent** |
| Nearest-vs-2nd-nearest field jump `[-50,50]²` | max **1.37e-4** off-pole (in-disc 4.9e-6), **0** cells > 1e-2 | the **Stage-2 Voronoi stitch is FINE** — *refutes* the "stitching seam" hypothesis |
| Pole **count** vs seed | **8159** (seed 0) vs **7987** (seed 42) | ~172 poles are spurious / seed-dependent; the seam **moves** with the seed |
| Two-seed difference `[-25,25]²` | max **8.5e-3** | the effect **explodes with window size** (8.5e-3 → 1.0 from R≈25 to R≈50) |

The decisive falsification: the seam **moves when only the random walk
order changes** (`seed_poles.jl`). A physical feature cannot do that.
And the Voronoi-stitch — the deep-research report's leading hypothesis
H1 — is *locally consistent* (jump 1.4e-4); the disorder is global
path-dependence, not local stitching. The bug is the WALK.

What is *not* the bug (ruled out by measurement):
- **The Stage-2 reconstruction** (nearest-vs-2nd-nearest jump tiny).
- **The pole-extraction machinery** — on Weierstrass-℘ (known exact
  lattice) recall = 99/99, precision = 1.000 vs the full lattice,
  conjugate symmetry = 62/62 (`pole_groundtruth_weierstrass.jl`).
- **Real-axis symmetry** as a *fix target* — PI(a) is conjugate-symmetric
  to 0.997 and 5-fold-symmetric to 0.998–1.000 (`pi_fivefold_symmetry.jl`);
  the gross sector structure is physical. (See §5.)

## 4. THE CORE PUZZLE — "if we do the same as FW, why don't they see it?"

We are **not** doing the same as FW, and FW in fact *do* see it:

- **FW see it.** FW 2011 md:208 (verbatim): *"It is typical that low level
  'ridges' appear within the flat areas… These irregularities are caused
  by the fact that different u entries… may have been obtained following
  entirely different paths… and then been amplified by the division by h².
  Since our IVP solver is stable within and near the pole fields, these
  irregularities will not appear in regions that are of concern in the
  present context."* That
  is our seam, exactly. FW dismiss it: it lives in the smooth regions they
  don't care about, and Fig 4.7 plots **pole locations**, where a
  flat-area field ridge is invisible.
- **FW md:304–310** (the passage attached to Fig 4.7): with real ICs the
  upper/lower halves *"should be identical (up to conjugation) but… will
  not be… for the numerical solution. The difference… gives a reasonable
  estimate for the numerical error."* — FW's solution is path-dependent
  too; they *use* it as an error gauge.
- **The amplifier is window size.** FW md:147 (verbatim): *"each pole
  field shown in Fig. 4.7a–e is a **5×5 composite of 25 completely
  independent runs** (each starting from the same ICs at z = 0, and
  covering areas of the extent **20×20** in the complex plane)."* FW never
  march one monolithic tree across the 100×100 window. Our measurement shows
  why that matters: path-dependence is **8.5e-3 at `[-25,25]²` but 1.0 at
  `[-50,50]²`** (same probe script with `HALF=25` vs `50`) — it grows
  ~orders of magnitude with distance from the IC, where the IVP is unstable
  (FW md:401; the path-dependent flat-area ridges are FW md:208). FW's 20×20
  tiles stay in the harmless
  regime; our monolithic solve does not, and grows large enough to corrupt
  **pole positions**, not just flat-field values — so it surfaces in the
  pole scatter.

So the bug is the method's inherent path-dependence (which FW contain by
short-window compositing and by ignoring/not-plotting the smooth field)
left **uncontained** in our monolithic single-tree solve.

## 5. On real-axis symmetry — a deliberate non-fix

`enforce_real_axis_symmetry=true` (`src/PathNetwork.jl:64-73`) walks one
half and mirrors via conjugation, making the negative-real-axis line
vanish *by construction*. This is a **band-aid, not a fix**: it masks the
error indicator (FW md:304–310) instead of removing the error, and does
nothing for the off-axis seams. The maintainer was explicit that forcing
symmetry is the wrong track. The negative-real-axis line is just the
*most visible* locus of the general path-dependence (the conjugate halves
are reached by independent branches); H3 in the deep-research report.

## 6. Meta-lesson — why it was never found (and how to find this class of bug)

Four lenses, each structurally blind to a localized field discontinuity:

1. **Global statistics** (recall / precision / conjugate-symmetry / FFT
   spacing; and the shipped `quality_diagnose` median/p90/p99/max,
   `src/Diagnostics.jl:159-162`) **average the seam away** — a line is a
   thin tail under a ~93 %-clean median.
2. **Pole-point analysis** (NN-distance, lattice residuals, in/out-of-wedge
   counts) is blind to a coherent *line*; a grain boundary differs in
   lattice **orientation**, not spacing, so only an orientation metric
   (ψ₆) sees it — and the seam carries no pole, so recall/precision pass.
3. **Proxy testing** — Weierstrass-℘ is doubly-periodic, *closes its loops*,
   and never grows the tronquée basin that produces the seam; the proxy
   stays green while the artifact is broken.
4. **"It's just Moiré"** — dismissed without falsification. A real seam is
   fixed under grid refinement and **moves** under `rng_seed`; Moiré is the
   opposite. Extracting centroids and re-rendering removes the Moiré and
   exposes the lattice.

**The principle that replaces all four:** *interrogate the FIELD and test
the INVARIANT, not the aggregate.* The sound invariant here is
**path-independence** (`u` reached two ways must agree). The reusable
diagnostic is the **two-seed test** — re-solve with a different `rng_seed`,
require agreement — plus the ψ₆ orientation map and a field phase/magnitude
portrait. Cheap, oracle-free, and decisive (it flags 35 % of cells here).

## 7. Next steps (gated; nothing started)

1. **Regression gate first** (bead `padetaylor-sny7`): `test/field_seam_test.jl`
   — solve panel (e)/(a) at ≥2 seeds; assert two-seed field agreement,
   nearest-vs-2nd-nearest jump ≤ tol off-pole, and real-axis Im-residual
   small; fail LOUD with the offending z-locus (Rule 1). Lands **RED**.
2. **Then the cure** (bead `padetaylor-vwgl`): enforce path-independence —
   FW's own remedy of a smooth-region BVP fill (`padetaylor-fe9`/`8dg`,
   FW md:170-192) and/or a bounded-window composite (FW md:147) and/or a
   K-nearest / overlap-consensus Stage-2 (`padetaylor-6b5`). The Voronoi
   stitch itself needs no change; the WALK's path-dependence does. `sny7`
   is the acceptance test for whichever cure lands.
3. Secondary: fix the `quality_diagnose` sheet-0 mask scope (it silently
   drops ~87 % of nodes on a [-50,50] PI window — `src/Diagnostics.jl:50-58`)
   and promote it from percentiles to a localized `max ΔP_rel` gate.

## Reproduction

```
julia --project=. external/probes/fig47-seam-diagnosis/p2_jumpmap_and_seeddiff.jl   # jump map + two-seed diff
julia --project=. external/probes/fig47-seam-diagnosis/seed_poles.jl                # seed-dependent pole counts
julia --project=. external/probes/fig47-seam-diagnosis/pole_groundtruth_weierstrass.jl  # machinery is clean on ℘
julia --project=. external/probes/fig47-seam-diagnosis/pi_fivefold_symmetry.jl      # gross structure is physical
```
