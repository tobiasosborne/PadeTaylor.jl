# Worklog 084 — Voronoi-boundary pole double-count / drop (bead `padetaylor-580u`)

**Date**: 2026-09-29
**Scope**: Measure whether `windowed_extract_poles` double-counts or drops a
physical pole that straddles a Voronoi core line; act on the measurement.

> Take-home: the hazard is REAL and it is not a rare corner.  On a
> real-symmetric problem every real-axis pole lies *exactly* on the
> `Im z = 0` core line of an even tiling; two windows place it at
> `±1e-8 i`, and the old per-estimate ownership rule kept it twice or
> not at all.  ℘ lattice over `[-20,20]²` (2×2 windows), seed 0:
> **4 exact poles dropped, 3 duplicated (4 cross-window kept pairs) of
> 246**.  Fixed by deciding ownership once per *physical* pole
> (cluster cross-window estimates, then assign by the group mean): **0
> dropped, 0 duplicated**.  The clashing `merge_atol` kwarg now forwards
> to `PoleField.extract_poles`.

## 1. The locus (verified against the tree, Law 1)

Pre-change `src/WindowedComposite.jl` (commit b9f4132): `windowed_extract_poles`
kept window `wi`'s estimate `p` iff `_nearest_center(p, centers) == wi`
(then :363), the optional greedy `merge_atol` dedup was off by default
(:354-356, :365-366) and never enabled by `edge_gated_windowed_poles`
(:460).  `PoleField.extract_poles` has its own, different `merge_atol`
(`src/PoleField.jl:191`, default `h_max`, :331) which the composite's kwarg
shadowed.  The only in-suite consumer, `test/field_seam_test.jl` FSEAM.2,
asserts two-way `match_frac` + a count delta — duplicate-blind, as the bead says.

## 2. Measurement (the oracle: Weierstrass ℘, exact lattice)

Fixture = `external/probes/fig47-seam-diagnosis/p6_ingn_pole_truth.jl`
(`u'' = 6u²`, equianharmonic ICs; exact poles `1 + 2ωm + 2ωe^{iπ/3}n`,
ω = 1.3629079683730673), `window_extent = 20, overlap = 6, h = 0.5`.
Probe: `external/probes/voronoi-boundary-580u/p1_legacy_ownership_wp.jl`
(re-derives the legacy rule inline).  Multiplicity = number of kept poles
within 0.3 of an exact pole; "cross-window pair" = two kept poles from
different windows within `r`; `coredist` = distance to the nearest core line.

| fixture | seed | windows | kept | exact | dropped | duplicated | x-window pairs (r=0.1…1.0) | spurious |
|---|---|---|---|---|---|---|---|---|
| ℘ `[-20,20]²`, order 30 | 0  | 4 (2×2) | 310 | 246 | **4** | **3** | **4** | 0 |
| ℘ `[-20,20]²`, order 30 | 42 | 4 (2×2) | 307 | 246 | **5** | 3 (1 is same-window) | **2** | 0 |
| ℘ `[-30,30]²`, order 30 | 0  | 9 (3×3) | 645 | 550 | 0 | 1 (same-window) | 0 | 0 |
| ℘ `[-30,30]²`, order 30 | 42 | 9 (3×3) | 643 | 550 | 0 | 0 | 0 | 0 |

Every `[-20,20]²` drop had `coredist = 0.000` and was *seen* by 2–4 windows,
e.g. exact `3.726` seen by windows 1–4 at imaginary parts `+6.6e-9, +3.4e-8,
+3.7e-8, -8.5e-9` — each estimate owned by a window other than its own, so
all four were discarded.  Every duplicate was a real-axis pole (`-1.7261`,
`-15.3564`, `-20.8085`, `17.3564` at seed 0) kept from windows 1&3 or 2&4 with
imaginary parts `∓0.0`.  Cross-window estimate spread of one pole: 1.3e-7
(`[-30,30]²`), ≤ 6.4e-5 (`[-20,20]²` seed 42).  On `[-30,30]²` the 3×3 core
lines sit at `±10`, the nearest exact pole is 0.097 from one, ≫ the spread,
so nothing straddles — which is why worklog 079 §4 (INGN, recall = precision
= 1.000 at HALF=30) saw nothing: that fixture cannot exercise the hazard,
and recall/precision are duplicate-blind anyway.

PI tritronquée (the FSEAM fixture, `[-30,30]²`, order 20;
`external/probes/voronoi-boundary-580u/p2_legacy_vs_fixed_pI.jl`): legacy
2160 / 2156 poles at seeds 0 / 42, cross-window kept pairs ≤ 0.1: 1 / 0.
(At r ≥ 0.3 the PI field's own dense spacing dominates — 381 *same*-window
pairs within 0.3 — so larger radii there are not a duplicate measure.)

Decision (bead options: default dedup on / rename-forward / non-issue):
**not a non-issue**, and turning the greedy `merge_atol` dedup on is the
wrong cure — it removes duplicates but cannot restore a DROP (both
estimates already rejected).  The original 2026-08-23 bead text named the
right fix: "global support-aware cluster before assignment".

## 3. Change

- `src/WindowedTiling.jl:165-195` — new `_own_poles(ps, centers, atol)`:
  union-find links estimates from **different** windows within `atol`
  (real-part-sorted sweep; two nested loops because a `break` in Julia's
  `for a in …, b in …` form exits both), each group gets one position (the
  mean) and one owner `_nearest_center(mean)`; only the owner's estimates
  are kept, in window-major order.  Groups of one reduce exactly to the old
  rule.  Module-docstring section "Pole ownership" (:59-90).
- `src/WindowedComposite.jl:371-381` — `windowed_extract_poles` now calls
  `_own_poles`, always on; new kwarg `boundary_atol` (default: the
  forwarded `cluster_atol`, else 0.1 = `extract_poles`' default); throws
  `ArgumentError` for `boundary_atol ≤ 0`.  The composite-level `merge_atol`
  is removed, so `merge_atol` flows to `PoleField.extract_poles` through
  `extract_kwargs` — one meaning package-wide.  Docstring rewritten; module
  docstring step 3 notes the pole-side step.  No in-repo caller passed
  `merge_atol` (grep over `*.jl` excluding `external/`).
- `docs/adr/0034-bounded-window-composite.md:68` — API line amended.
- Effective LOC after change: WindowedTiling 56, WindowedComposite 148.

Result after the change: ℘ `[-20,20]²` seed 0 → kept 310, 0 dropped, 0
duplicated, 0 cross-window pairs, 0 spurious.  PI FSEAM fixture: 2159 / 2156
(the one seed-0 cross-window duplicate removed, nothing else changed);
FSEAM.2 still 97.2 % / 97.0 %, Δcount 3.  ℘ `[-30,30]²` (no straddles):
fixed rule output identical to legacy at both seeds (645 / 643 kept, same
multiplicities) — the change is inert where the hazard is absent.

## 4. Test and mutation proof

New `test/windowed_pole_ownership_test.jl` (10 assertions, ~12 s):
- **WPO.1** non-triviality: the legacy rule, re-derived inline, drops ≥1 and
  duplicates ≥1 exact pole on this fixture (measured 4 / 3).
- **WPO.2** the fix: every exact pole has multiplicity exactly 1; 0
  cross-window kept pairs within 0.5; 0 spurious; driver == `_own_poles`.
- **WPO.3** `merge_atol = 3.0` genuinely changes per-window extraction and
  reaches it through `windowed_extract_poles`; `boundary_atol = 0` throws.

Mutations (each run, observed RED, restored by copying back the saved file):
- **M1** `src/WindowedTiling.jl:194` `owner[root(a)] == win[a]` →
  `_nearest_center(z[a], centers) == win[a]` (the legacy per-estimate rule):
  WPO.2 **2 of 4 FAIL** (7 pass / 3 fail overall; the third failure was a
  since-replaced WPO.3 assertion `merge_atol = 0` changes extraction, which
  is false on ℘ and was corrected to `merge_atol = 3.0`).  Restored; GREEN 10/10.
- **M2** re-add a swallowing `merge_atol = nothing` kwarg to
  `windowed_extract_poles` (the pre-580u clash): WPO.3 **1 FAIL** at
  `test/windowed_pole_ownership_test.jl:101` (9 pass / 1 fail).  Restored.

## 5. Open

- The same-window duplicate at ℘ `[-30,30]²` seed 0 (and `[-20,20]²` seed
  42) is PoleField-level froth inside ONE window, not a Voronoi effect —
  out of this bead's scope; not fixed here.
- `boundary_atol` default 0.1 is justified by measured spread ≤ 6.4e-5 and
  PI spacing ≳ 0.3; a field denser than `cluster_atol` would need it lowered
  (the extractor's own `cluster_atol` would already be wrong there).
