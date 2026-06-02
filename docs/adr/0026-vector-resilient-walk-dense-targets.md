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

## Amendment 5 (2026-05-22) — S5-diag: the ~8 % plateau is NOT a ceiling

The S5a re-probe (corrected S2/S3/S4 stack) cured the collapse (~0
failures) but coverage plateaued at ~8 %; the S5a agent called it an
"architectural ceiling." The S5-diag FW-lens investigation (bead
`padetaylor-f4e`, instrumented) **falsifies that** — it is two fixable,
measured, composable geometry problems.

**H-B confirmed (hard).** The honest B1 disc is **100 %
truncation-limited** — across 4737 walk nodes, the Jorba–Zou truncation
term binds 4737/4737 times; the pole-adjacency clamp *never* binds. The
disc sits at a constant `0.106 × D_local` (the nearest-pole distance) —
it stops ~10× short of the nearest pole. `SURF_PN_ORDER = 24` is simply
too low. Order sweep (inner wedge): order 24→36→48→72 lifts coverage
13 %→23 %→29 %→29 % — the knee is order ~36–48; order 72 is wasted.
(ADR-0025 Amendment 2's "an order-48 walk is not viable" was measured on
the *pre-`:skip` fixed-h* walk and is superseded — the resilient
adaptive stack runs order 48/72 cleanly, 0 failures.)

**H-A refuted.** There is **no node pile-up** — the pile-up factor is
1.1× (4737 nodes ≈ 4453 distinct disc-cells); the walk does not re-tread
or wander (tree edge length = `h_v` exactly); 8.3 nodes/target. The
orchestrator's interim "~40× pile-up" arithmetic was wrong (it used an
over-large disc radius). The discs are genuinely too small to tile:
4737 discs of r ≈ 0.013 sum to area 12 over an inner wedge of area 36.

**Lever B (the genuine H-A residue).** The step `h` (median 0.032) is
**2.4× the honest disc** (0.013) — the walk over-strides where each Padé
is honestly valid, leaving NaN gaps along every filament (Amendment 3
bottleneck 1, still live). Sizing `h` to the disc via `SAFETY` 0.25 →
0.10 lifted coverage 23 %→31 % at fixed order 36 — measured.

**Decision — step S6 (the fix).**

- **S6a** — raise `SURF_PN_ORDER` 24 → ~36–48 (figure constant); the
  measured knee. Widens the truncation-limited honest disc.
- **S6b** — lower `SAFETY` (`src/VectorWedgeStep.jl`) from 0.25 toward
  ~0.10–0.15 so the adaptive step is sized to the honest B1 disc, not
  2.4× it. A measured sweep `SAFETY ∈ {0.25,0.15,0.10}` picks the value
  (smaller `SAFETY` ⇒ more nodes / more runtime — a real trade).
- **S6c** — re-probe the full wedge at the tuned (order, `SAFETY`),
  *sweeping target density*, to resolve the residual: after A+B, is the
  remaining shortfall a pure node-count limit (lay denser targets —
  feasible) or does it need FW's explicit coarse-Stage-1 / 16×-finer-
  Stage-2 fill split (FW-md:163-164, :395-397)? Measured, not assumed.

Neither lever is exhausted; the ~8 % is not architectural.

## Amendment 6 (2026-05-22) — S6 outcome: order 36 ships (8 %→22 %); `radius_t` heresy blocks `SAFETY`

**S6a shipped.** `SURF_PN_ORDER` 24 → 36 (the measured knee). Full-wedge
honest coverage rose from the ~8 % plateau to **~22 %** (s=0.30) — a
real 2.7× gain, every covered cell genuinely valid.

**S6b blocked.** The sweep confirmed `SAFETY` 0.25→0.10 lifts inner-wedge
coverage +9 pp — but lowering `SAFETY` *regresses the `src/` test suite*
(4 fail / 3 error: empty pole fields). Root-caused: `extract_poles_
shared_q` filters shared-`Q` roots by `radius_t` in the **rescaled
variable `t = Δz/h`** — i.e. it keeps poles within a fixed number of
*steps*. A smaller `SAFETY` ⇒ smaller `h` ⇒ the same physical pole maps
to a larger `|t*|` ⇒ it falls outside the `radius_t` window ⇒ the
`min_support ≥ 2` cross-node filter empties the field. **`radius_t` is a
scale-fixing heresy** — and the S4 review ("`radius_t` is fine post-S2")
was *reasoned, not measured*, and is wrong (Rule 5: measurement wins).
`SAFETY` stays 0.25 until `radius_t` is fixed.

**Outer-wedge regression noted.** The full-wedge target-density sweep
(order 36): coverage is 22 %/22 %/16 % at spacing 0.30/0.22/0.16 — it
*regresses* as targets densify. Banded: the inner band climbs
monotonically (31→34→37 %) but mid and outer collapse (outer 22→17→10 %).
Pure target densification is exhausted for the outer wedge.

**Decision — step S7 + the outer-wedge question.**

- **S7** — fix the `radius_t` heresy: `extract_poles_shared_q` must
  accept a shared-`Q` root by its **z-plane distance** `h·|t*|` against a
  scale-stable z-radius (tied to a fixed scale such as `h_max` / the
  local pole spacing — not the per-node shrinking `h`). TDD +
  mutation-proof. Then lower `SAFETY` toward ~0.10 (now unblocked),
  confirm the `src/` suite GREEN, re-probe.
- **The outer-wedge regression** is characterised honestly as remaining
  work — likely FW's explicit coarse-Stage-1 / finer-Stage-2 fill split
  (FW-md:163-164, :395-397) rather than denser walk targets — filed as
  its own bead, to be measured (not assumed) after S7.

## Amendment 7 (2026-05-22) — S7 shipped; the outer wedge is the remaining blocker

**S7 shipped.** The `radius_t` scale-fixing heresy is fixed — `extract_
poles_shared_q` now accepts a shared-`Q` root by its z-plane distance
`h_node·|t*| ≤ radius_t·h_max` (a scale-stable, h-independent window),
and the cluster representative is chosen by z-plane distance, not `|t*|`
(a second facet, root-caused while fixing the ℘-oracle). Both
mutation-proven. With the coupling removed, `SAFETY` is lowered 0.25 →
0.10 (sizes the step to the honest B1 disc); the `src/` suite is GREEN.

**The coverage state.** Full-wedge honest coverage: ~5 % (original) →
~8 % (S2–S5 corrected stack) → **~22 %** (S6a order 36) → **~22 %** (S7).
S7's `SAFETY` gain is real on a dense-target sub-region (inner-wedge
+8 pp) but does **not** move the full-wedge total — at the s=0.30 target
density the figure is target-limited, not step-limited, and denser
targets *regress* the outer wedge. The figure is honest at every step
(no Padé evaluated outside its verified disc) and is 4× less blank than
the start — but it is not "filled."

**The single remaining blocker — the outer wedge (bead `padetaylor-cfq`).**
Banded coverage at order 36: inner ~31 % (climbs with target density),
mid ~19 %, outer ~22 % (collapses to ~10 % as targets densify). The
total is outer-area-dominated. This is *not* a ceiling and *not* a
walk-failure (failures ≈ 0.1 %); it is an uncracked coverage-efficiency
problem in the outer wedge — the working hypothesis is FW 2011's
explicit coarse-Stage-1 / finer-Stage-2 fill split (FW-md:163-164,
:395-397), to be measured. Cracking it is an open-ended investigation;
the figure is otherwise complete and honest at ~22 %.

## Amendment 8 (2026-05-22) — outer-wedge mechanism cracked; FW coarse/fine split rejected; step S8

The outer-wedge regression (bead `padetaylor-cfq`) is **root-caused with
measurements**:

- The adaptive step `h_v` is **4.4× smaller** in the outer wedge
  (median 0.012) than the inner (0.053).
- The honest B1 disc `R_B1 = min(s·h_JZ·h_v, 0.5·h_v·min|t*|)` is
  **proportional to `h_v`** — both terms carry an `h_v` factor — so
  `R_B1` collapses 4.5× with it (outer 0.018 vs inner 0.081).
- Tiny outer discs ⇒ they pile rather than tile: measured 62 % of the
  outer disc-area is wasted on mutual overlap; outer coverage saturates
  at ~0.26 *independent of node count* (29 k nodes → 0.259; 60 k nodes →
  0.252). That is the "more targets, less coverage" paradox resolved —
  coverage is node-count-bound and the discs are too small to tile.

**This is scale non-covariance, again.** KKG theory (Amendment 3) gives
the pole spacing as ∝|x|⁻¹ᐟ⁶ — only ~1.5× across the wedge — yet `h_v`
collapses 4.4×. The adaptive law `h = SAFETY·D_local` drives `h_v` down
far steeper than the physics warrants.

**The FW coarse/fine split is NOT the fix — measured and rejected.** A
uniform coarse Stage-1 grid feeding `_stage2_fill` saturates at the same
~0.26 outer ceiling as the dense fan: coverage is a single monotone
function of node count, independent of target placement. FW's split
works for FW only because FW's Stage-2 reaches every fine point with
*no validity test*; our honest B1 gate (Rule 1, no extrapolation) caps
each node at `R_B1 ≈ 0.018` — the blocker is disc *size*, not target
architecture. Bead `cfq` closed with this verdict.

**Step S8.** Diagnose whether the 4.4× `h_v` collapse is driven by
*genuinely* denser outer poles or by *spurious near-origin shared-`Q`
roots* poisoning `min|t*|` (the vector shared-`Q` pipeline does GGT
degree-reduction but, unlike scalar `RobustPade`, no explicit
Froissart-doublet suppression — a prime suspect). Then fix the binding
constraint so the honest disc tracks the *pole field*, not the
over-collapsed step. Honest scoping (per the cfq diagnostic): S8 lifts
the realistic full-wedge coverage toward the low-to-mid 30 %s; a
*completely* filled wedge is bounded by the conservative honest B1 gate
and is likely not reachable without violating the no-extrapolation
contract — that limit, if it binds, will be measured and reported, not
assumed.

## Amendment 9 (2026-05-23) — S8: Amendment 8's mechanism withdrawn; the genuine cap is the honesty contract

The S8 diagnostic (bead `padetaylor-roe`, instrumented, live repo)
**falsifies Amendment 8's outer-wedge mechanism** and converges the
investigation.

**Amendment 8's "4.4× `h_v` collapse" is not reproducible.** S8 measured
the *shipped* config (order 36, `SAFETY` 0.10, full-wedge walk): `h_v`
median inner 0.0325 / outer 0.0350 — ratio **0.93×**; honest disc
`R_gate` inner 0.049 / outer 0.052 — uniform. The `cfq` diagnostic left
no reproducible probe; per Rule 3 (the live repo is authoritative) the
uniform measurement stands. The likely reconciliation: `cfq` probed an
*outer-targets-only* walk (threading unsupported from the seed straight
to `|x|∈[14,20]`), which is not the figure's *full-wedge* walk —
different walks, and only the full-wedge one is the figure.
**Amendment 8's mechanism and step S8's Froissart premise are
withdrawn.** S8 also established, on the code facts, that
`SharedPade.shared_denominator_pade` *does* carry GGT-2013 SVD
rank-count + degree-reduction + `padeapprox.m` QR-reweighting — i.e. it
*does* suppress Froissart doublets; the "missing Froissart suppression"
suspicion was wrong. (One minor hygiene item: `VectorStepper` builds the
canonical `Q` at `default_tol = 1e-14` rather than the walk's `1e-8`
noise floor — a frontier-only second-order effect, filed low-priority,
not a coverage lever.)

**The genuine, converged picture** (consistent across S5-diag, S6, S7,
S8 — only `cfq`'s mechanism was the outlier):

- The walk is **healthy** — ~0 % failures, no collapse, `h_v` uniform.
- The honest B1 disc is **uniformly truncation-limited**; coverage is
  ~22 % at order 36, ~29 % at order 48 (the order sweep saturates ~48).
  Taylor order is the only walk-side lever, and it is nearly exhausted.
- Denser targets do **not** lift total coverage (measured S5a, S6).

**The real cap — the no-extrapolation honesty contract.** ADR-0025
deleted the `extrapolate = true` corner for this figure: *no Padé is
evaluated outside its verified B1 disc.* That is *stricter than FW
2011* — FW's published pole-field figures FILL precisely because FW's
Stage-2 evaluates each Padé over its full step **with no validity gate**
(FW-md:395-397). Our headline figure is therefore honest-but-partial
**by design**: ~22–29 % is the B1-gate-certified coverage, and a
*completely* filled wedge is unreachable by any walk engineering while
the no-extrapolation contract holds. This is not a bug and not a walk
failure — it is a measured, principled tension between ADR-0025's
honesty contract and a visually-filled figure. The resolution is a
design decision (keep the strict gate and ship the honest ~29 %, or
render FW-style with the certified region marked) — not further walk
work. The investigation (D1–S8) is complete.

## Amendment 10 (2026-05-26) — tf9 closeout: dual-fill provenance ships; D1–S8+tf9 arc complete

Bead `padetaylor-tf9` — the final step of `padetaylor-0ln.40` — ships
the headline figure's last design decision. Amendment 9 converged the
investigation: certified coverage saturates at ~22 % (order 36) →
**~29 % at order 48** under the no-extrapolation contract; the cap is
the contract itself, not a walk defect, and a *completely* filled
wedge is unreachable while the contract holds. **Decision** (cross-
referenced from ADR-0025 Amendment 14): render FW-faithful filled,
mark provenance.

### What ships

- **Taylor order**: `SURF_PN_ORDER` 36 → **48** in
  `figures/_kkg_pi2_surface_helpers.jl` — the Amendment-5 measured
  saturation knee. Biggest certified core (~29 %) at ~2× the
  order-36 runtime; order 72 is wasted.
- **Dual-fill helper**: `surf_wedge_fill` runs **one** Stage-1 walk
  and evaluates the Stage-2 fill **twice** on the cached walk —
  `extrapolate = false` for the B1-certified core (stored as
  `walk.grid_y`, the package's normal output), then a direct call to
  `PadeTaylor.VectorPathNetworkStage2._stage2_fill` with
  `extrapolate = true` for the FW-style filled overlay. One walk, two
  fills; the second fill is a Horner sweep against the cached
  visited tree.
- **Kernel API additions**: `kkg_pi2_surface()` now returns
  `Re_u_extrap` and `Im_u_extrap` alongside `Re_u`/`Im_u`. In the
  sector the extrap matrices mirror the certified ones (the sector is
  voted, not Padé-evaluated). In the wedge they carry the FW-style
  filled rendering. `surf_wedge_fill` returns a `u_extrap` field
  alongside the existing `u`.
- **Figure render**: `figures/kkg_pi2_tritronquee_surface.jl`
  composites the wedge as a **two-layer alpha overlay** — the
  extrapolated layer at `EXTRAP_ALPHA = 0.50` (the FW-style fill where
  the certified core leaves NaN) drawn first, then the certified
  layer at full opacity on top. The wedge pole field overlay is
  unchanged.
- **No `src/` API change**: Amendment 10 ships exclusively in
  `figures/`. The kernel reaches into the existing private
  `_stage2_fill` rather than re-exporting an `extrapolate = true`
  convenience — the legitimate scope of a per-figure provenance-marked
  exception.

### Why reduced alpha (not hatch / desaturation)

Alpha was chosen over hatch overlay and desaturation because it
composes natively in Makie's `heatmap!` / `surface!`, reads cleanly on
both the 2D heatmap and the 3D surface panels in one layout pass, and
the alpha drop is a visually-monotone signal — a reader scanning the
figure sees that certain cells are "washed out" and instinctively
attributes that to lower confidence. Hatch would require a separate
overlay layer and double the render-time cost; desaturation would
flatten the `:RdBu` diverging map's sign-reading at the wedge cells
where it matters most. `EXTRAP_ALPHA = 0.50` is the chosen value: at
0.45 the overlay is too pale to read pole structure on the
colour-clamped wedge, at 0.55 the visual distinction from
full-opacity cells is too subtle on the 2D heatmap panels.

### Measured certified coverage at order 48

The exact certified fraction at order 48 is *measured at render time*
and surfaced in three places: the figure's `message` field, the
figure script's stdout, and worklog 060. The Amendment-9 ladder is
`~5 % → ~8 % → ~22 % → ~29 %`; the actual value at the production
render is the worklog 060 figure.

### Tests

- VC-4 / VC-5 / VC-7 / VC-8 / VC-10 are all unchanged — they touch
  the walk and the pole field, not the Stage-2 fill mode. PI2S.1
  through PI2S.12 stay GREEN.
- PI2S.4's load-bearing `covered ⟺ isfinite(Re_u/Im_u)` invariant
  continues to hold — `Re_u` / `Im_u` are the **certified** matrices
  (the dual-fill design's senior-grade core); the new
  `Re_u_extrap` / `Im_u_extrap` matrices are *additions*, not
  replacements. The "covered" boolean always meant "B1-certified";
  it still means exactly that.
- The full `src/`-touching test suite — `test/vector_path_network_test.jl`,
  `test/vector_path_network_stage2_test.jl`,
  `test/vector_pipeline_oracle_test.jl` — stays GREEN. No `src/`
  change ships in this amendment.

### Closes the D1–S8 + tf9 arc

The arc is complete. D1 (resilient walk) shipped Amendment 1 / Bead
`0ln.40.a`. D2 (dense targets) shipped Bead `0ln.40.b`, falsified by
Amendment 1's probe. D3 / R1 / S1–S8 (the corrected stack — scale-
derived adaptive step, skip-if-covered, scale-derived clustering,
order knee, `radius_t` fix) shipped through Amendments 2–9. The
honesty-contract cap was the converged finding (Amendment 9). The
dual-fill provenance design (this amendment) is the closeout —
ADR-0026 closes here. Worklog 060 records the full ladder, the three
scale heresies fixed (`h_max`, `radius_t`, `cluster_atol`), and the
"premature ceiling claim" lesson: two prior agents called "ceiling"
wrongly; the only thing that wasn't a ceiling was the no-extrapolation
contract — a documented design decision, not a numerical limit.

## Amendment 11 (2026-06-02) — the §S7 fix ported to the scalar `PoleField` twin (bead `padetaylor-bez`)

The §S7 scale-covariance fix (Amendment 6/7) was applied to the
*vector* extractor `VectorPoleField.extract_poles_shared_q` only; its
**scalar twin** `PoleField._extract_poles_core` (`src/PoleField.jl`)
was never ported and still carried the original two heresies:

- the far-root filter kept roots by the per-node `|t*| ≤ radius_t` (a
  window of a fixed number of *per-node* steps `h`), so a pole at a
  fixed z-distance was silently dropped whenever the per-node `h`
  shrank below the walk's ceiling — intermittent empty / sparse pole
  fields on any varying-`h` (stitched / adaptive) trajectory; and
- the greedy cluster sort key was `|t*|`, so under a varying `h` the
  cluster representative drifted to a coarse far node (smallest `|t*|`)
  rather than the z-plane-closest node.

Both halves are **one atomic fix** (this ADR's §S7 / Amendment 7): the
filter now accepts by z-plane distance `|h·t*| ≤ radius_t · h_max`
(with `h_max = maximum(|hs|)`, the `h`-independent step ceiling, and a
`zero` placeholder for an empty `hs` exactly as the vector code does),
and the cluster sort key is the z-plane distance `|h·t*|`, so the
representative is the z-plane-closest node. For a near-uniform-`h`
walk (`h ≈ h_max`) the filter reduces *exactly* to the legacy
`|t*| ≤ radius_t`, so the existing FW Table-5.1 / Fig-4.7 pole-field
tests (PF.1.*/2.*/3.*) are unchanged.

Verified by new test **PF.4.1** in `test/polefield_test.jl` against the
closed-form equianharmonic-℘ lattice pole `℘_pole(0,0) = 1`: a stitched
coarse-`h`(0.5)/fine-`h`(0.05) trajectory at `min_support = 4` drops
the genuine z=1 pole under the old per-`|t*|` filter (only 3 nodes
survive) and recovers it under the z-distance filter (all 8 nodes),
with the representative placed at 2.94e-9 (z-closest node) vs 4.57e-7
(old min-`|t*|` far node). Mutation-proved M5 (revert filter → RED) /
M6 (revert sort key → representative drifts, RED) — see the test
footer. No public-API or default-`radius_t` change. This is a pure
bug-port to bring the scalar twin into parity with §S7.

## References

- `docs/worklog/059-headline-figure-honest-reassessment.md`
- `docs/adr/0025-headline-figure-re-resolution.md` (Amendment 13)
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
  — §3.1 path network (:151–168), §5.4.2 Stage 2 (:395–397),
  conjugate-symmetry diagnostic (:303–310)
- `src/VectorPathNetwork.jl`, `src/VectorWedgeStep.jl`,
  `src/VectorPathNetworkStage2.jl`, `src/PathNetwork.jl`
