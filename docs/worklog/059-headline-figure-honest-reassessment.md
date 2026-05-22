# Worklog 059 — headline figure: honest reassessment

**Date**: 2026-05-22
**Author**: Claude Opus 4.7 (orchestrator), with the user
**Epic**: `padetaylor-0ln` (v0.2) · **Follows**: worklog 058
**Scope**: A correction. Worklog 058 reported the P_I⁽²⁾ tritronquée
headline figure "re-resolved … shipped". On inspection of the rendered
figure the user judged it inadequate — correctly. This worklog records
the honest diagnosis and points to the real remaining work.

> **Take-home**: the headline figure does **not** visualise the
> solution in the pole-rich wedge — it shows a sparse ~266-pole scatter
> over a ~95 %-blank region. The root cause is structural and was
> mis-scoped, not a bug: `src/VectorPathNetwork.jl` is a deliberately
> **minimal ~156-LOC port** of v0.1's full **1108-LOC**
> `src/PathNetwork.jl` driver. The 5-direction wedge and the Stage-1
> tree were ported; the machinery that actually *fills a region* — a
> resilient fine-grid walk and the barycentric Stage-2 fill — was not.
> The fix is bead `padetaylor-0ln.40`.

## What worklog 058 got right, and what it got wrong

Worklog 058's Phases A–F are all real: the six v1 corners genuinely
retired, the seven-criterion FW-style validation suite genuinely built
and mutation-proven, the full suite genuinely 5326/5326 GREEN. None of
that is withdrawn.

What 058 got wrong was the **frame**. It treated "re-resolve the figure"
as complete when ADR-0025's scoped work completed. But ADR-0025's scope
was *retire the v1 corners + certify the result* — taking the existing
vector path-network as given. The figure it produces is honest (no Padé
evaluated outside its disc) but **mostly empty**: in the pole wedge it
shows where ~266 poles are and nothing else. KKG Figs 7.4/7.5 show the
*solution surface*. Ours does not. "Re-resolved" was the wrong word for
"made honest and well-tested but still 95 % blank in the wedge."

## The diagnosis (precise)

The wedge is rendered by `vector_path_network_solve` (`src/Vector
PathNetwork.jl`). Reading that file end-to-end:

- **What is there.** The FW 2011 §3.1 **5-direction wedge** (`DEFAULT_
  WEDGE = [-π/4,-π/8,0,π/8,π/4]`, enforced to 5) and the **Stage-1 path
  tree** (`visited_parent`, walk-from-nearest-visited). B2 even improved
  the wedge selector and made the step adaptive. The user's belief that
  "the 5-direction path tree" is implemented is **correct**.

- **What is not there.** The module's own opening docstring states it:
  *"the vector analogue of v0.1's `PathNetwork.path_network_solve` — but
  deliberately only its Stage-1 tree construction, not the whole
  1010-LOC v0.1 driver."* `src/PathNetwork.jl` is **1108 lines**;
  `src/VectorPathNetwork.jl` is **~156**. Explicitly scoped out of the
  vector port (per its docstring): the Stage-2 **barycentric** fill,
  **non-uniform node placement**, **Schwarz-reflection symmetry**,
  branch-cut tracking, the diagnostics layer.

- **The three concrete consequences** for the figure:
  1. **Coarse targets.** The figure walks the tree to **171** target
     points (19 radial shells × 9 rays). FW's *surface* fills run the
     tree to a full fine grid (FW Fig 3.2 — a 40×40 = 1600-point grid,
     with discs at step `h≈0.5` big enough to *tile*). 171 coarse
     targets with our `h≈0.1` discs cannot tile an area-~250 wedge.
  2. **A brittle walk.** An unreachable target, an all-5-candidates-fail,
     or an adaptive-step collapse each does `throw(...)` — a *single*
     bad target aborts the *entire* run. Survivable for 171 targets;
     fatal for the ~25 000 a genuine fill needs (with no skip-and-
     continue, the chance all 25 000 are cleanly reachable is ~0).
  3. **Nearest-node-only Stage-2.** `VectorPathNetworkStage2` is a
     nearest-node lookup, not the v0.1 barycentric multi-node fill.

- **The arithmetic.** Pole density forces a small step (`h≈0.1`); the
  Padé disc is ~`h`; the wedge has area ~250. Tiling needs ~25 000+
  discs. The minimal brittle walk lays down ~2 000 then aborts → ~5 %
  honest coverage → a 95 %-blank wedge.

## Why this was not caught earlier — the orchestration lesson

The Phase-A2 exploration spike *measured* this: honest coverage by the
minimal walk ceilings at ~5–18 %, and "a genuinely filled surface needs
a different filling strategy." The orchestrator (me) recorded that as a
**deferred bead** (`0ln.38`, "2D re-expansion lattice") and proceeded to
polish the minimal-walk figure for the remaining five phases.

That was the wrong call, twice over. First, the deferral should have
been a **stop**: the figure cannot be a headline figure until the walk
can fill the wedge, and that is prerequisite work, not a follow-up.
Second, the deferred bead's framing — "a 2D re-expansion lattice, a
different architecture" — was itself wrong. It is **not** a different
architecture. It is the **full FW path-network driver** — the resilient
fine-grid walk + barycentric fill that v0.1's `PathNetwork.jl` already
contains and that the vector port simply never received. Bead `0ln.38`
is superseded by `0ln.40` accordingly.

The honest signal the orchestrator should have raised: *"the vector
solver is a skeleton of FW's method; the headline figure needs the full
driver first."*

## The work that is actually next — bead `padetaylor-0ln.40`

Port v0.1 `PathNetwork.jl`'s full driver to the vector case, then re-do
the figure's wedge as a genuine filled surface. It is a **plan-first
deep-dive** (the `0ln.37` pattern): recon → ADR-0026 → child beads →
serial implementation → re-verify. Key pieces, for the next agent:

- a **resilient** walk — skip a locally-unreachable target / wedge-fail
  / step-collapse and continue, instead of aborting the whole run;
- a **fine-grid** target set at disc spacing, not 171 coarse points;
- the **barycentric** Stage-2 multi-node fill;
- consider Schwarz-reflection symmetry and non-uniform node placement.

The deep-dive **must measure** whether the full driver actually achieves
good wedge coverage — the A2 ~18 % ceiling was measured on the *minimal
brittle* walk; a resilient fine-grid driver should do far better
(envelope: 25 000 nodes × π(0.1)² ≈ 785 ≫ 250), but this must be
verified, not assumed. If a residual fraction of the wedge is genuinely
intractable, that must be measured and reported honestly, not deferred.

`0ln.37` and ADR-0025 stay closed — their scoped work (corners +
validation) is genuinely done; `0ln.40` is the larger piece they
surfaced. The validation suite (VC-4…VC-10) and the true-radius gate
built in 058 are reusable as-is by the re-done figure.
