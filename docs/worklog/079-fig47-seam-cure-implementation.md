# Worklog 079 — FW 2011 Fig 4.7 path-dependence seam: cure implementation and verification

**Date:** 2026-06-30. **Author:** Claude Opus 4.8 (1M).
**Beads:** `padetaylor-vwgl` (P0, the bug — **interior seam CURED**),
`padetaylor-xwzf` (impl — `src/WindowedComposite.jl`, DONE),
`padetaylor-sny7` (seed-invariance gate — `test/field_seam_test.jl`, DONE),
`padetaylor-ingn` (℘ verification — PASS),
`padetaylor-mro9` (per-window edge-gating — **superseded** by post-hoc filter,
closed), `padetaylor-us19` (boundary-flicker residual, open P3).
**Predecessor worklogs:** 077 (root cause), 078 (scoping + mechanism selection).
**ADR:** `docs/adr/0034-bounded-window-composite.md`.
**Status:** **IMPLEMENTATION COMPLETE AND VERIFIED.** The interior seam is
cured (97.6–97.8 %); figures rewired; PNG regeneration is the one remaining
render step.

> **Take-home.** The scoping in worklog 078 selected the cure; this session
> delivered it. Three things shipped: (1) `src/WindowedComposite.jl` — the
> bounded-window composite module, the production pole driver, and its
> supporting types; (2) `test/field_seam_test.jl` — a non-gameable,
> mutation-proven seed-invariance gate; (3) `edge_gated_windowed_poles` — the
> production figure driver, required after discovering that a plain windowed
> solve re-introduces smooth-sector bloom (7353 poles / 4475 off-wedge vs the
> edge-gated baseline's 1496 / 0). The Weierstrass-℘ oracle (probe p6, bead
> `padetaylor-ingn`) confirmed the composite's fewer poles are spurious
> seed-dependent poles, not real poles dropped — windowed recall=1.000,
> precision=1.000 > monolithic 0.998. The vwgl seam is cured in the interior
> pole lattice (97.6–97.8 % seed-invariant); residual boundary-pole flicker
> (79.9–89.6 % perimeter match) is filed as bead `padetaylor-us19` (P3, not
> v1.0-blocking).

---

## 1. What shipped

Three deliverables, all committed and gated.

**`src/WindowedComposite.jl`** (~199 eff LOC; literate chapter, under the
200-LOC ceiling). Public API exported:

- `windowed_path_network_solve(prob, xs, ys; window_extent=20.0, overlap=6.0,
  h=0.5, order, rng_seed=0, extrapolate=true, verbose=false)` — the
  bounded-window composite solve per FW 2011 md:147. Tiles the domain into
  `window_extent×window_extent` cores, solves each independently from the shared
  IC at `_window_seed(rng_seed, wi)`, and composites by Voronoi-core ownership
  (nearest core centre wins; the owning window always solved the cell, so every
  cell is written exactly once — a hard partition, no boundary seam). Defaults
  reproduce FW's recipe.
- `windowed_extract_poles(wsol; merge_atol=nothing, extract_kwargs...)` — pole-side
  counterpart of the field composite: extract poles from each window's independent
  solve, keep only those in the window's Voronoi core, union across windows.
- `edge_gated_windowed_poles(prob, xs, ys; ...)` — the production figure driver for
  Fig 4.7 / 4.8: run `edge_gated_pole_field_solve` to get the bloom-free
  pole-field mask, run `windowed_path_network_solve` for the seam-free composite,
  extract windowed poles, filter to the mask. Seam-free and bloom-free.
- `WindowedCompositeSolution{T}` — result struct carrying the composited `u` field,
  per-window solves, Voronoi assignments (`window_assign`), window extents, core
  centres, and per-window seeds.

**ADR-0034** (`docs/adr/0034-bounded-window-composite.md`) — documents the
decision, the P3/P4/P4b/P5 measurement basis, the per-window seeding design, the
alternatives considered and rejected, and the honest residuals carried forward.

**`test/field_seam_test.jl`** (bead `padetaylor-sny7`) — seed-invariance gate.
Config `[-30,30]² / order 20`. Two testsets:

- **FSEAM.1** (per-window re-randomisation): asserts `window_seeds` differ between
  global seeds for every window — proves two global seeds genuinely re-randomise
  every window's walk. Non-gameable: a frozen-seed composite passes the pole
  invariance trivially; FSEAM.1 catches it.
- **FSEAM.2** (pole-set seed-invariance): asserts the composited pole field is
  seed-invariant — bidirectional pole-set match ≥ 0.90. Composite achieves
  97.0–97.2 % GREEN; monolithic 77.6 % RED.

**`figures/fw2011_fig_4_7.jl`** and **`figures/fw2011_fig_4_8.jl`** — both
rewired from `edge_gated_pole_field_solve` to `edge_gated_windowed_poles`.

## 2. The seed-invariance gate and mutation proofs

The gate was mutation-proven on two independent mutations, each executed and
restored byte-clean:

- **M-frozen-seed** (`_window_seed` hardcoded to return a constant regardless of
  the global seed): FSEAM.1 goes RED — the per-window re-randomisation assertion
  fires immediately. The gate cannot be gamed by freezing the seed pipeline.
- **M-monolithic** (`windowed_path_network_solve` replaced by plain
  `path_network_solve`): FSEAM.2 goes RED — the seed-invariance assertion fires
  (monolithic match 77.6 %, below the 0.90 floor). The gate catches the
  regression it was written for.

Both mutations confirmed RED and restored byte-clean before the session committed.

## 3. BLOOM discovery and the edge_gated_windowed_poles fix

During figure rewiring, plain `windowed_path_network_solve` on the PI tritronquée
re-introduced smooth-sector **bloom** — IVP error inflating spurious poles off the
pole wedge. Measured on `[-50,50]²`:

| solver | total poles | off-wedge poles |
|---|---|---|
| `edge_gated_pole_field_solve` (baseline) | 1496 | 0 |
| `windowed_path_network_solve` (plain) | 7353 | **4475** |
| `edge_gated_windowed_poles` (fix) | **1452** | **0** |

The fix is composition: run the edge gate for the bloom-free pole-field footprint,
run the windowed solve for the seam-free composite, extract windowed poles, keep
only those whose nearest grid cell — or any of its 8 neighbours (matching the
gate's own 1-ring dilation) — is `true` in `field_mask`. Result: 1452 poles, 0
off-wedge, 97.6 % agreement with the edge-gated baseline in the pole field.

The `_in_field_mask` helper is inlined (not `_dilate`, per ADR-0023), avoiding a
`BitMatrix` allocation per pole. Bead `padetaylor-mro9` (per-window edge-gating, a
different approach to bloom) was superseded by this post-hoc filter driver and
closed.

## 4. INGN ℘ verification (bead padetaylor-ingn) — PASS

The P5 confirmation (worklog 078 §4) showed the composite finds ~220 fewer poles
than the monolithic walk (7223 vs 7443). Bead `padetaylor-ingn` asked whether
those missing poles are spurious or real.

Probe `p6_ingn_pole_truth.jl` (Weierstrass-℘ oracle, `u'' = 6u²`, `HALF=30`): a
known doubly-periodic lattice where every pole location is computable to high
precision, so recall and precision have a genuine truth anchor.

| metric | windowed composite | monolithic |
|---|---|---|
| recall (fraction of real poles found) | **1.000** | 1.000 |
| precision (fraction of found poles that are real) | **1.000** | 0.998 |

Windowed recall = 1.000 = monolithic: the composite drops **zero** real poles.
Windowed precision = 1.000 > monolithic 0.998: the composite additionally
suppresses spurious seed-dependent poles that the monolithic walk generates. The
composite's ~220 fewer poles are entirely the spurious ones correctly removed.
Bead `padetaylor-ingn`: **PASS**.

## 5. Boundary-flicker residual (bead padetaylor-us19)

The post-cure two-seed comparison on the full tritronquée (with
`edge_gated_windowed_poles`, over `[-50,50]²`):

- **Interior pole lattice** (≥3 rings inside the pole-field mask):
  seed-invariant at **97.6–97.8 %**. The seam is cured.
- **Perimeter poles** (boundary of the masked field):
  **79.9–89.6 %** match between seeds.

The interior seam is gone. The perimeter flicker is a separate, milder
phenomenon. The edge-gate mask is itself seed-stable (`edge_gated_pole_field_solve`
runs a deterministic seed-0 region-grow and takes no `rng_seed`, so both seeds are
filtered through the *same* mask). The flicker is in the windowed poles near that
fixed boundary: perimeter poles sit at the sparse edge of the pole field where
pole strength is marginal, so a small seed-dependent IVP-error difference flips
whether such a pole clears the extraction threshold (or lands just inside vs just
outside the 1-ring mask). The two-seed all-poles match (89.6 % / 94.4 %, Δcount 88)
is entirely this boundary-pole flicker; the interior lattice is clean.

This is benign for a single-seed figure: both seeds show the same interior pole
lattice; only the precise set of perimeter poles varies. Filed as bead
`padetaylor-us19` (P3, not v1.0-blocking).

## 6. Status — the vwgl seam is CURED (interior)

- `src/WindowedComposite.jl`: shipped, literate, ≤200 eff LOC.
- `test/field_seam_test.jl`: GREEN (composite 97.0–97.2 %, monolithic 77.6 %).
- Mutation proofs: M-frozen-seed → FSEAM.1 RED; M-monolithic → FSEAM.2 RED.
  Both executed + restored byte-clean.
- `figures/fw2011_fig_4_7.jl` and `fw2011_fig_4_8.jl` rewired to
  `edge_gated_windowed_poles`.
- INGN ℘ verification (probe p6): windowed recall=1.000, precision=1.000. **PASS**.
- Bead `padetaylor-vwgl` (P0): **CLOSED** — interior seam cured.
- Bead `padetaylor-us19` (P3): open — boundary-pole flicker, not blocking.
- **Remaining step:** regenerate `figures/output/fw2011_fig_4_7.png` and
  `fw2011_fig_4_8.png` (the figure scripts are rewired; one `julia` run needed,
  then push).

## Reproduction

```
# Full-composite acceptance (P5 probe, worklog 078):
julia --project=. external/probes/fig47-seam-diagnosis/p5_full_composite_confirm.jl

# Weierstrass-℘ oracle — INGN verification (bead ingn):
julia --project=. external/probes/fig47-seam-diagnosis/p6_ingn_pole_truth.jl

# Seed-invariance gate (standalone):
julia --project=. test/field_seam_test.jl
```
