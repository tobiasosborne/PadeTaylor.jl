# ADR-0026 — Vector path-network: a resilient walk + dense targets for the headline figure

**Status**: Accepted (2026-05-22). Supersedes the mis-framed deferred
bead `padetaylor-0ln.38` ("2D re-expansion lattice — a different
architecture"). **Bead**: `padetaylor-0ln.40`. **Follows**: worklog 059
(the honest reassessment), ADR-0025 Amendment 13.

## Context

Worklog 059 correctly judged the P_I⁽²⁾ tritronquée headline figure
inadequate: its pole-rich wedge is ~95 % blank. It then diagnosed the
cause as "`src/VectorPathNetwork.jl` is a minimal ~156-LOC port of
v0.1's 1108-LOC `src/PathNetwork.jl`; the region-filling machinery —
resilient walk, fine-grid targets, **barycentric Stage-2 fill** — was
never ported," and framed the fix (bead `0ln.40`) as "port the full FW
path-network driver."

A four-agent deep-dive recon (2026-05-22; the `0ln.37` plan-first
pattern) read every relevant file and the FW 2011 markdown end-to-end.
**Two of worklog 059's four diagnosed gaps do not exist.** This ADR
records the corrected diagnosis and scopes the actual work.

## The corrected diagnosis

| Worklog 059 claim | Recon finding |
|---|---|
| Missing the **barycentric Stage-2 fill** | There is no barycentric fill to port. The scalar `PathNetwork.jl:669–693` Stage-2 is *also* plain nearest-node single-Padé lookup. FW 2011 itself prescribes exactly that — "one single Padé step suffices … no tests needed" (FW-md:166, :395–397); the word *barycentric* appears in FW only for BVP Chebyshev post-processing. A multi-node distance-weighted blend *was* built and auditioned on the vector side (the "B4 audition") and **rejected on measured evidence** (`VectorPathNetworkStage2.jl:133–183`). Nothing to port; nothing missing. |
| Missing the **Stage-2 fine-grid fill** | It was ported — `VectorPathNetworkStage2.jl` exists (bead `0ln.20`, commit `1d3c32b`), wired through the `fine_grid` kwarg. Worklog 059 is stale here. |
| **Brittle walk** — `throw`-aborts on first failure | **True** — and the scalar `PathNetwork.jl` is *equally* brittle (it has no skip-and-continue either). Resilience is genuinely *new* work for both drivers, not a port. |
| **171 coarse targets**, a fill needs ~25 000 | **True.** The solver already accepts any target vector verbatim (`VectorPathNetwork.jl:428–442`, `:498`); the figure helper simply passes a sparse 9×19 fan. A denser target set is a *caller-side* change. |

The vector path-network is therefore **not a skeleton**. It already
carries machinery the scalar driver lacks — the B1 true-radius Stage-2
gate and the B2 `:max_q_root` adaptive wedge walk (ADR-0025). A
wholesale 1108-LOC port would *duplicate and in places regress* it.

The figure is blank for exactly two reasons: **(1)** the Stage-1 walk
`throw`s and aborts the entire run on the first unreachable target, and
**(2)** the figure asks for only 171 filament-endpoint targets, so the
path-network tree is 171 filaments with wide gaps between them rather
than a space-filling tiling.

## Decision

Three binding decisions. No 1108-LOC port; no Schwarz reflection (the
user has ruled it out, and FW 2011 uses conjugate symmetry only as an
accuracy *diagnostic*, FW-md:303–310, never as a work-halving shortcut);
no barycentric Stage-2 (it does not exist in FW and was already
evidence-rejected here).

### D1 — A resilient Stage-1 walk (new capability, `src/VectorPathNetwork.jl`)

`vector_path_network_solve` gains an additive kwarg
`on_target_failure::Symbol = :throw`:

- `:throw` (default) — byte-identical to today. Every existing call
  site and the ~2600 existing assertions are unaffected.
- `:skip` — the per-target walk is wrapped so that a failure (any of
  the three throw sites: unreachable `VectorPathNetwork.jl:515`,
  all-five-candidates-fail `VectorWedgeStep.jl:308`, adaptive-step
  collapse `VectorWedgeStep.jl:380`) **does not abort the run**. The
  failing target is recorded and the walk continues to the next target.

This is **not** a silent swallow — that would violate Rule 1. It is
fail-*loud-by-accounting*: every skip is captured as first-class data.
`VectorPathNetworkSolution` gains a tenth field `failed_targets`, a
concretely-typed vector of records, one per skipped target, carrying:
the target point, a `reason::Symbol`
(`:unreachable` | `:all_candidates_failed` | `:step_collapse`), the
`z` at which the walk stuck, and the residual distance `|z − target|`.
Backward-compat constructors (6-/8-/9-arg) default it to empty, exactly
as `visited_jets` was added (`VectorPathNetwork.jl:298–315`).

Partial nodes laid down by a walk that later failed are **kept** — they
are valid solution nodes (the walk pushes to `visited_*` only on an
*accepted* step) and still contribute coverage. A failed walk is "this
target was not reached," never "discard the filament."

The contract: with `on_target_failure = :skip`, the solver returns
normally; a non-empty `failed_targets` is a signal the *caller must
surface* (the figure reports the count and the achieved coverage).

### D2 — Dense, disc-spaced wedge targets (caller-side, the figure helper)

The figure's wedge target generator (`surf_wedge_targets` /
`SURF_TARGET_RADII`/`SURF_TARGET_ANGLES` in
`figures/_kkg_pi2_surface_helpers.jl`) is replaced by a generator that
**tiles the `|arg x| < 36°` wedge with targets at disc spacing**, so
adjacent path-network filaments' validity discs overlap. FW's analogue
is its 40×40 = 1600-node coarse grid feeding a 161×161 fine grid
(FW-md:153–164). The B1 true-radius gate measures honest disc radii at
median ≈ 0.06, p90 ≈ 0.16 (ADR-0025); target spacing is therefore a
*measured* parameter — start at ≈ 0.25 (≈ 2× the disc radius), render,
**measure the achieved coverage**, and tighten until coverage saturates
or runtime forces a stop. The expected envelope is a few thousand
targets; the realised count and coverage are reported, not assumed.

### D3 — Walk failures are investigated through the FW lens

A skip is a safety net, not an answer. Per CLAUDE.md Rule 2 ("all bugs
are deep"), every distinct failure cluster surfaced by `failed_targets`
is **root-caused before the figure is called done**, with the explicit
question: *how did FW solve the analogous problem?* FW's random-path
strategy empirically reached "within h of all target points" (FW-md:163)
— FW reports essentially **zero** unreached targets. So a vector walk
that skips a non-trivial fraction is exhibiting a *bug*, and FW's
near-100 % success is the benchmark. Candidate FW-side analogues to
check against: the §3.1 nearest-visited-node restart, the §5.5
"complete the pole field before stepping into smooth regions" path
ordering, the (order, h) = (30, 0.5) standard. If, after investigation,
a residual frontier is genuinely intractable, it is **measured and
reported honestly** (the ADR-0025 Amendment-13 "numerically-supported
frontier" doctrine), never quietly deferred.

## Out of scope

- The wholesale 1108-LOC port of `PathNetwork.jl` (recon shows it is
  unnecessary and partly regressive).
- Schwarz / conjugate-reflection symmetry (user-ruled-out; FW-faithful
  to omit).
- A barycentric / multi-node Stage-2 blend (absent from FW; already
  evidence-rejected by the B4 audition).
- Non-uniform `node_separation`, branch-cut tracking, the diagnostics
  layer — the other vector-port deferrals; none is needed for a filled
  honest wedge.
- The pole-free sector (the triple-method majority vote, ADR-0024) —
  sound and unchanged.

## Verification / acceptance

1. The B2 `:max_q_root` adaptive walk, the B1 true-radius Stage-2 gate,
   and the VC-4…VC-10 validation suite (`figures/test_kkg_pi2_surface.jl`)
   are **unchanged and stay GREEN**.
2. New resilience behaviour is TDD'd in `test/vector_path_network_test.jl`
   and **mutation-proven** (perturb the skip/accounting logic, confirm
   RED, restore).
3. `on_target_failure = :throw` default ⇒ the full `Pkg.test()` suite
   stays GREEN with no assertion changes (additive kwarg + additive
   struct field).
4. The headline figure re-renders with a **genuinely filled** pole
   wedge; the achieved honest coverage fraction is **measured and
   reported** in worklog 060 and in the figure's docstring.
5. Any residual unreached frontier is investigated per D3 and reported.

## Child beads

- `0ln.40.a` — resilient walk: `on_target_failure` kwarg +
  `failed_targets` field + accounting (D1). TDD + mutation-proof.
- `0ln.40.b` — dense disc-spaced wedge target generator (D2).
- `0ln.40.c` — re-render the headline figure; measure & report coverage.
- `0ln.40.d` — FW-lens investigation of walk failures + remediation (D3);
  conditional on what `.c` surfaces.
- `0ln.40.e` — re-verify (VC suite + full `Pkg.test()`), worklog 060,
  ADR finalisation, epic bookkeeping.

## Amendment 1 (2026-05-22) — the coverage probe falsifies the density-alone premise

D1 (resilient walk) shipped and is sound. D2 (dense targets) shipped its
generator — and a fast coverage probe (`figures/probe_wedge_coverage.jl`)
then **falsified D2's premise that denser targets alone fill the wedge.**

Probe, spacing 0.35, 1922 dense targets:

| metric | value |
|---|---|
| visited nodes | 135 949 (≈ 70 per target) |
| failed targets | 494 = **25.7 %**, almost all `:step_collapse` |
| honest coverage | **1.7 %** — *worse* than the legacy 171-fan's 13–18 % |

Density made coverage **regress**. The root cause is not target count: it
is that the B2 adaptive-`h` controller (`VectorWedgeStep._adaptive_h`,
ADR-0025 Amendment 4) **strangles the walk**. FW 2011's walk takes
"less than a single step per target" (FW-md:163–164) at a *fixed*
`h = 0.5`; ours takes ≈ 70 steps per target and tips 25.7 % of them into
an `h < h_min` collapse. The adaptive controller appears to ratchet `h`
down monotonically in a sustained dense pole field and cannot recover —
and the walk-termination threshold `|z − target| > visited_h[parent]`
inherits the collapsed `h`, forcing the walk to *crawl* to within
≈ 10⁻⁴ of every target. Each collapsed-`h` node also carries a tiny B1
disc, so 136 k nodes still tile almost nothing.

**Consequence for the plan.** The blank wedge has *two* causes, not one;
worklog 059 and this ADR's original framing saw only the target count.
The walk's step-size strategy is the deeper bug. Decision **D3** (the
FW-lens investigation, bead `0ln.40.d`) is therefore **promoted from a
conditional follow-up to a hard prerequisite of the re-render** — per
worklog 059's own lesson that "the deferral should have been a stop."
Bead dependency re-wired: `0ln.40.c` (re-render) now depends on
`0ln.40.d` (investigation), not the reverse. D3 must, with FW 2011 open,
root-cause the crawl/collapse and decide the walk's step-size strategy
(candidates: FW-faithful fixed `h`; a non-ratcheting adaptive law;
decoupling the termination threshold from the collapsed `h`) on
*measured* evidence, not assumption.

## Amendment 2 (2026-05-22) — D3 diagnosis: the adaptive-`h` geometric sink; decision R1

The D3 FW-lens investigation (bead `0ln.40.d`) root-caused the collapse
with instrumented experiments. Diagnosis (full detail in
`external/probes/vector-walk-collapse/`, reproducible, gitignored):

**Root cause.** `VectorWedgeStep._adaptive_h` is a *geometric sink* for
`h`. Its pole cap `POLE_SAFETY·h_prev·min|t*|` is proportional to
`h_prev`, so the step recurrence is purely multiplicative:
`h_next/h_prev = min(GROW, POLE_SAFETY·min|t*|) = min(1.5, 0.5·min|t*|)`.
`h` grows only while `min|t*| > 3`; it shrinks geometrically whenever
`min|t*| < 3` and **cannot recover** — `h = 0` is an attractor. The cap
carries no *absolute* memory of the pole distance, only a ratio.
Measured: 86.4 % of walk nodes have a collapsed `h < 0.05`. Worse, the
FW-faithful nearest-visited chaining then propagates a collapsed `h`
from a "poisoned" node into every later walk that chains off it — so
denser targets *spread* the poison (Exp 5: a target the solver flagged
`:step_collapse` walks cleanly in 111 steps at `h = 0.1` when started
fresh from the seed — the target is not hard, the parent node was
poisoned). The loop-termination test `|z−target| > visited_h[parent]`
(`VectorPathNetwork.jl`) then inherits the collapsed `h` and forces the
~70-steps-per-target crawl.

**The FW lens.** FW 2011 §3.1 (FW-md:151–168) uses a *single fixed
global `h`*; the 5-direction wedge picks only *direction*. A uniform
`h` is what makes FW's nearest-node chaining sound. B2 (ADR-0025
Amendment 4) made `visited_h` a per-node field varying 1000×, breaking
that invariant.

**The decisive finding.** The fixed-`h` "block past |x|≈8" that
motivated B2 was *mis-attributed to the step size*. It is a rare
on-pole chord; ADR-0026 **D1's `:skip` driver — which did not exist
when B2 was designed — already converts it to a single
`:all_candidates_failed` skip** and the run completes. Measured (inner
wedge, 295 targets): fixed `h = 0.1` gives **0.18 coverage with 1152
nodes**; adaptive `h_max = 0.1` gives **0.07 coverage with 5805 nodes**.
Fixed `h` is FW-faithful, simpler, and measurably better.

**Decision R1.** The headline-figure wedge walk runs at a **fixed `h`**:
`surf_wedge_fill` passes `adaptive = false`. The package default of
`vector_path_network_solve` (`adaptive = true`) is *unchanged* — the
`_adaptive_h` geometric-sink defect is real but latent (it only
catastrophises for a dense, cross-target-chained walk), so its proper
redesign (an *absolute* non-ratcheting pole cap, "R2") is filed as a
separate bead with that exact forcing condition, per Rule 9. The
arrival-threshold decouple ("R3") is automatic under fixed `h`
(`visited_h` is then uniform `= h`). The genuine open question — does a
fixed `h` thread the *outer* wedge `|x| → 20` where pole density grows,
or does the outer wedge force a finer global `h`? — is settled by
*measurement* in the R1 re-probe over the full wedge, not assumed.

## Amendment 3 (2026-05-22) — full pole-field diagnosis; the corrected plan

A four-workstream diagnosis (pole census + `src/` heresy audit + figure
heresy audit + KKG-theory query) resolved the coverage question — and
corrected Amendment 1/2's framing.

**Theory (KKG 2015, cited via `references/tex/.../tritronquee_coeff.tex`).**
The P_I⁽²⁾ tritronquée pole spacing scales as **|x|⁻¹ᐟ⁶** — a *very weak*
power (Stokes exponent `(6/7)x⁷ᐟ⁶`). Over |x| ∈ [2,20] that is only a
~1.5× scale variation. The uniform-scale coordinate is **ζ = x⁷ᐟ⁶** (the
P_I⁽²⁾ analogue of FFW 2017's `log` map for PIII/PV).

**Measurement (245 extracted poles).** Pole spacing is ~flat at **0.70**
across |x| ∈ [5,20] — consistent with the weak |x|⁻¹ᐟ⁶ law within noise.
The honest B1 disc radius is a **constant 0.106 × the local pole
spacing**, uniformly across the wedge: the system already *behaves*
scale-covariantly.

**The correction.** Amendment 1/2 (and the orchestrator's interim
"scale-mismatch ceiling" hypothesis) are **partly superseded**: the
field is *nearly uniform-scale* over |x| ≤ 20, so a fixed `h = 0.1` is
~14 % of the local pole spacing — roughly appropriate everywhere. The
~12 % coverage is **not** dominantly a scale mismatch. There is no
"architectural ceiling" either (Amendment 1's interim worry). R1 (fixed
`h`) correctly fixed the *collapse*, but fixed `h` is itself a scale
heresy and is not the coverage answer.

**The real bottlenecks (measured, none architectural):**

1. **`h` is ~2× too big.** The honest disc is ~0.044; the step `h` is
   0.1. The walk steps past where each Padé is honestly valid → gaps
   *along every filament*. `h` should be ≈ the honest disc.
2. **Coverage is node-density-limited.** Honestly tiling the wedge needs
   ~40–100 k well-spread discs of radius ~0.044; the runs lay ~8 k. Pure
   node count — feasible (the collapsed adaptive run laid 135 k).
3. **The walk re-walks covered ground.** It skips a target only on
   *exact* coincidence (`VectorPathNetwork.jl` `10·eps` test), so dense
   targets make it re-traverse covered area; the redundant chords clip
   poles → the non-monotonic `:all_candidates_failed` regression.

**The scale-fixing heresies are real** (`src/` audit: the `h = 0.5`
default in `VectorPathNetwork.jl`/`PathNetwork.jl`; figure audit:
`SURF_PN_H`, `cluster_atol = 0.2`, `radius_t = 5`). A general solver
*must* be scale-covariant, so they are fixed on principle — but for this
near-uniform field they are hygiene, not the coverage lever. The
exception: `cluster_atol = 0.2` is an absolute cluster tolerance that
wrongly merges/splits poles → it affects the pole *count*.

**Unresolved — gates the fix.** The two audits *conflict* on the
`_adaptive_h` collapse mechanism (D3 "geometric sink" vs the `src/`
audit "the pole cap is a genuine absolute z-distance, legitimate").
Working hypothesis: spurious near-origin shared-`Q` roots poison
`min|t*|` (consistent with the 37 % Froissart-doublet VC-4 prune rate).
This must be resolved first — any pole-distance-based `h` law is
vulnerable to it.

**The corrected plan (supersedes R1 as the endpoint):**

- **S1** — resolve the `_adaptive_h` / spurious-root mechanism (focused
  diagnostic) so the new `h` law rests on solid ground.
- **S2** — make `h` *scale-derived*: from the local pole distance, sized
  ≤ the honest B1 disc. Implements scale covariance, fixes bottleneck 1,
  tracks |x|⁻¹ᐟ⁶ for free. (Re-scopes bead `cqk`/"R2".)
- **S3** — *skip-if-covered*: the walk skips a target already inside a
  visited node's disc → kills the chording regression (bottleneck 3).
- **S4** — scale-derive `cluster_atol` (a fraction of the local pole
  spacing) → correct pole count.
- **S5** — dense targets + S1–S4 → re-probe, measure coverage climbing
  toward full, then render.
- ζ = x⁷ᐟ⁶ coordinate — principled, low payoff over |x| ≤ 20; deferred
  to its own bead.

## Amendment 4 (2026-05-22) — S1 verdict: the geometric sink, and the S2 step law

The S1 diagnostic (bead `padetaylor-4q6`, instrumented) resolved the
`_adaptive_h` collapse definitively.

**Verdict — a genuine geometric sink.** The recurrence is
`h_next = h_prev · min(GROW, POLE_SAFETY·min|t*|)` — *both* terms are
proportional to `h_prev`, so there is no `h`-independent floor and
`h = 0` is an attractor whenever `min|t*| < GROW/POLE_SAFETY = 3` is
sustained. The two prior audits were a *false dichotomy*: D3 was right
that it is a sink; the `src/` audit was right that the cap
`POLE_SAFETY·h_prev·min|t*|` is a genuine *scale-covariant absolute
z-distance* (S1 EXP D: `h·min|t*|` constant to 5 sig figs across a 16×
span of `h`). Both hold — a true absolute cap can still be a sink,
because the cap is `POLE_SAFETY × (the parent's disc)` and the parent's
disc tracks `h_prev` downward.

**The spurious-root hypothesis is killed.** S1 EXP B rebuilt the
cap-driving shared-`Q` root at orders 16–30 and build-steps 0.25–4×:
the root maps to a z-plane pole stable to 3 decimals — a *real* pole,
not a Froissart artefact. **`min|t*|` from the shared-`Q` is
trustworthy**; S2 may build on it directly, no separate Froissart-safe
pole indicator needed.

**Decision — the S2 step law.** `h` must be an *absolute* function of
the local pole field with **no multiplicative dependence on `h_prev`**:

```
h_next = clamp( SAFETY · D_local ,  h_min ,  h_max )
```

where `D_local = h · min|t*|` is the `h`-independent absolute pole-free
disc radius (any `h` at a fixed node yields the same `D` — S1 EXP D).
The `min(GROW·h_prev, …)` ratchet is **dropped entirely**: `h` jumps
straight to a fraction of the *true* local disc each step, with no
memory of how small `h_prev` became — a node leaving a dense pocket
recovers to `h_max` in one step. `GROW` may be retained *only* as an
optional rate-limiter on *increases* (`h_next ≤ GROW·h_prev`), never on
decreases. The `:step_collapse` throw is kept for the genuine case
`SAFETY·D_local < h_min` (the true local pole density really does
exceed what `h_min` can honestly resolve) — but it should now be rare,
firing only on genuine density, never on a self-inflicted ratchet; the
D1 `:skip` driver absorbs it.

## References

- `docs/worklog/059-headline-figure-honest-reassessment.md`
- `docs/adr/0025-headline-figure-re-resolution.md` (Amendment 13)
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
  — §3.1 path network (:151–168), §5.4.2 Stage 2 (:395–397),
  conjugate-symmetry diagnostic (:303–310)
- `src/VectorPathNetwork.jl`, `src/VectorWedgeStep.jl`,
  `src/VectorPathNetworkStage2.jl`, `src/PathNetwork.jl`
