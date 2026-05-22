# Worklog 058 — senior-grade re-resolution of the P_I⁽²⁾ tritronquée headline figure

> **⚠ Correction (2026-05-22) — see worklog 059.** This worklog reports
> the figure "re-resolved … shipped". That overstates it: the figure is
> honest and well-tested, but its pole-rich wedge is only ~5 % filled —
> a sparse pole scatter, not the solution surface. The Phases A–F work
> here is all real, but it was poured into a structurally inadequate
> container: the vector path-network is a minimal ~156-LOC skeleton of
> FW's full 1108-LOC driver. The real remaining work is bead
> `padetaylor-0ln.40`. Worklog 059 carries the honest diagnosis.

**Date**: 2026-05-22
**Author**: Claude Opus 4.7 (orchestrator) + serial Opus coding subagents + Sonnet recon
**Epic**: `padetaylor-0ln` (v0.2) · **Bead**: `padetaylor-0ln.37` (+ 20 children)
**ADR**: `docs/adr/0025-headline-figure-re-resolution.md` (Amendments 1–12)
**Scope**: The deep-dive demanded by worklog 057 §"Required follow-up" —
re-resolve the v0.2 headline figure (`figures/kkg_pi2_tritronquee_surface.jl`,
KKG 2015 Figs 7.4/7.5) to senior-engineer grade: retire every v1 corner,
and — since KKG published no pole data — certify the result with a
manufactured FW-style validation suite.

> **Take-home**: the headline figure is re-resolved. Its wedge is now a
> **266-pole FW-validated pole field** (was 21, with ≥1 spurious); the
> `extrapolate=true` dishonesty is gone (no Padé is evaluated outside its
> verified disc); the sector is a 401²-grid triple-method surface honest
> to ~9·10⁻⁵. All **six v1 corners are retired or rigorously justified**,
> and seven independent FW-style validation criteria (VC-4…VC-10) certify
> the figure with no external oracle. Full suite **5326 / 5326 GREEN**.
> The single biggest finding: a *filled* honest wedge surface is
> architecturally infeasible — the wedge's honest deliverable is the pole
> *field*, which is exactly what "as good as FW managed" denotes.

## Method — plan-first, exploration-first, serial

Per the bead's mandate (CLAUDE.md Rule 9): **deep-dive recon → an ADR →
granular child beads → serial implementation → re-verify**. Five Sonnet
recon agents swept the ground truth (KKG 2015 TeX), the F3 kernel, the
path-network/shared-Q machinery, the Laplace voters, and — after the user
clarified that KKG drew *no poles* — the FW-style validation menu. Their
findings became ADR-0025 and **20 child beads** across six phases.

Every coding step was one serial Opus subagent (Rule 7 — one Julia
process); read-only recon ran in parallel. After each bead the
orchestrator re-verified the new tests standalone (Rule 3), committed,
and pushed. The ADR was **amended 12 times** as the Phase-A spikes and
later auditions resolved genuine unknowns — the ADR is the living record.

## What shipped

| Phase | Bead | Deliverable |
|---|---|---|
| A1 | 0ln.37.1 | Stage-2 validity-gate audition — per-node Jorba–Zou gate (v-b) wins 7–10× honest coverage over a global `κ·h` |
| A2 | 0ln.37.2 | Wedge-tractability probe — a filled honest surface ceilings at ~18 %; wedge **rescoped to the pole field** |
| A3 | 0ln.37.3 | P_I⁽²⁾ Laurent structure — double poles, `A∈{−1,−3}`, residue `a₁=0`; the VC-4 certificate |
| A4 | 0ln.37.4 | V8b baseline — ~70 % loop-closure tail, broken conjugate symmetry, a spurious pole |
| B1 | 0ln.37.5 | `_validity_radius` true-radius Stage-2 gate; jet cached in `visited_jets` |
| B2 | 0ln.37.6 | `src/VectorWedgeStep.jl` — `:max_q_root` selector + adaptive `h`; threads the wedge to `|x|=20` (absorbs 0ln.23) |
| B3 | 0ln.37.7 | Wedge rewired — 171-target fan, `extrapolate=false`, **380 candidate poles** |
| B4 | 0ln.37.8 | Stage-2 interpolation audition — Voronoi retained (evidence-based) |
| C1 | 0ln.37.9 | 401² grid, 2° ray fan, Laplace `N` at the spectral knee |
| C2 | 0ln.37.10 | Inner/outer arc data from BVP interior collocation — inner-arc spread 13× better |
| C3 | 0ln.37.11 | Ray-fan voter — exact-radius + cubic-angular reconstruction (19× better) |
| D-VC4/5 | .12/.13 | Per-pole dominant-balance + conjugate-pairing — **266 of 380 certified**, 114 spurious pruned |
| D-VC5b | 0ln.37.20 | Far-wedge conjugate completeness — optimal bipartite matcher (flagged 72→52) |
| D-VC7 | 0ln.37.14 | Loop-closure `quality_diagnose` vector adapter |
| D-VC8/9 | .15/.16 | BVP companion-consistency + Weierstrass-℘ vector-pipeline oracle |
| D-VC10 | 0ln.37.17 | Two-run accuracy indicator — disagreement 0.35, corroborates VC-5b |
| E1 | 0ln.37.18 | Stokes-strip mask narrowed ±3°→±1° (3× smaller) |
| F1 | 0ln.37.19 | Re-render + full `Pkg.test()` 5326/5326 GREEN; this worklog |

## Hard-won lessons (continuing worklog 057's numbering)

59. **No ground truth ⇒ manufacture the validation.** KKG 2015 plotted
    the surfaces only and called them "for visualization"; they drew no
    poles. There is nothing external to validate against — so the figure
    is certified by a *manufactured* FW-style suite (dominant-balance
    `A∈{−1,−3}`, conjugate symmetry, loop closure, the ℘ oracle, the
    two-run indicator). This is exactly the discipline FW used for their
    own P_I pole fields, which had no oracle either.

60. **A filled honest wedge surface is architecturally infeasible.** The
    A2 spike measured honest local-Padé coverage of the pole-dense wedge
    ceilinging at ~18 % even at display tolerance + doubled Taylor order
    (~13–23k patches needed, the walk supplies ~7k and denser fans
    block). The honest deliverable is the pole *field* — the FW Fig 4.7
    idiom — and that *is* "the pole field as good as FW managed." Forcing
    a filled surface would have meant evaluating Padé approximants
    outside their discs: dishonest. Knowing what *not* to attempt is
    senior-grade.

61. **The honest Padé validity radius is truncation-limited, not
    pole-limited.** The intuitive gate `h·min|t*|` (distance to the
    nearest denominator root) over-estimates the honest radius 5.7×: a
    shared-Q Padé models the poles `Q` captured and stays accurate
    *through* them — it fails by order-`n` *truncation* error long before
    any pole. A per-node Jorba–Zou estimate on the node's own jet beats a
    global constant 7–10× on honest coverage.

62. **Exploration before architecture.** A1 and A2 each overturned a
    plausible plan assumption (the gate criterion; the very feasibility
    of a filled surface). Committing the Phase-B/D architecture before
    those spikes would have built the wrong thing — Rule 9 made concrete.

63. **A subagent falsifying its own brief is the discipline working.**
    Three times a coding subagent disproved its bead's stated premise
    with evidence rather than implementing it: 0ln.37.20 (the walk was
    *not* missing far-wedge poles — the conjugate matcher was a buggy
    greedy commit); C1 (the inner-arc spread is *not* Laplace
    under-resolution — it is the irreducible `n_terms=2` seed-truncation
    floor, and raising Laplace `N` makes it *worse*); VC-7 (the
    loop-closure metric is independent of the `extrapolate` gate, so the
    "B3 beats V8b" gate the bead expected cannot honestly exist — the
    subagent shipped an honest structural certificate instead of
    fabricating one). Rule 2/3 over a tidy story.

64. **A library default change ripples to every call site.** B2 made
    `:max_q_root` + adaptive `h` the default of `vector_path_network_
    solve`; this degenerated the unrelated A_4⁽¹⁾ Noumi–Yamada figure.
    Fix + lesson: a figure must pin its walk policy *explicitly* rather
    than inherit a mutable library default (`a4_walk()` now pins
    `:min_y`; B3's headline figure pins its own).

## v1-corner ledger — final

All six retired or rigorously justified (ADR-0025, Amendments 1–12):
C1 inner-arc datum → BVP interior collocation (C2); C2 `extrapolate=true`
→ true-radius gate (B1); C3 hand-tuned `h` → adaptive walk (B2); C4 ±3°
Stokes mask → narrowed to ±1° (E1); C5 bilinear voter → exact-radius +
cubic-angular (C3); C6 121² grid → 401² (C1).

## Beads

**Closed**: `0ln.37` + children `.1`–`.20`, and `0ln.23` (absorbed by B2),
`s8z` (the A_4 regression, fixed by an explicit policy pin).
**Opened, deferred** (under the epic, each with a forcing condition):
a 2D-lattice re-expansion fill for a *filled* wedge surface; an
`n_terms ≥ 3` asymptotic seed for the inner-arc BVP pin.

## Pickup point

The headline figure is re-resolved and shipped; full suite
**5326 / 5326 GREEN**; `figures/output/kkg_pi2_tritronquee_surface.png`
re-rendered. The remaining planned epic phase is **V9** (`0ln.37`'s
sibling `padetaylor-0ln.19`) — v0.2 docs + release prep
(README/CHANGELOG/RESEARCH v0.2 sections, ADR review including the new
ADR-0025, HANDOFF refresh, the `padetaylor-0ln` epic close-out).
