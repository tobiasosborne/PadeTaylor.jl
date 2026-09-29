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

## 6. Review follow-up (after commit 5c8fe12)

Three reviewer claims, each measured before acting.  Probe:
`external/probes/voronoi-boundary-580u/p3_owner_absent_groups.jl` (seed 0,
linking radius 0.1, same fixtures as §2).

**6.1 Owner-absent groups** (the mean's nearest centre has no member in the
group, so 5c8fe12 discarded the whole group):

| fixture | groups | owner-absent singletons | owner-absent multi-window | largest group diameter |
|---|---|---|---|---|
| ℘ `[-20,20]²` | 312 | 2 (1 near a lattice point) | 0 | 1.78e-2 |
| ℘ `[-30,30]²` | 646 | 1 (1 near a lattice point) | 0 | 9.35e-2 |
| PI `[-30,30]²` | 2817 | 647 | **15** | 0.164 |

℘ oracle: no in-domain lattice pole is lost to an owner-absent singleton.
WPO.2 (`[-20,20]²`) and the p1 probe (`[-30,30]²`: mult0 = 0 at both
seeds) find every in-box exact pole kept exactly once, so the lattice
hits among the discarded singletons are poles outside the box or second
estimates of a pole that is already kept.  Singletons are therefore
**unchanged**: dropped, as under the pre-580u rule.  On PI they are
dominated by the un-gated windows' smooth-sector bloom (ADR-0034:
4475 of 7353 poles off-wedge un-gated), so keeping them would be wrong.
The 15 PI multi-window owner-absent groups are poles that two or three
independent windows agree on (< 0.1) but that the core owner did not
resolve.  Both the pre-580u rule and 5c8fe12 dropped them.

Decision: a **multi-window** owner-absent group is now kept once, from the
member window whose centre is nearest the mean (lowest index on a tie):
`src/WindowedTiling.jl` `_own_poles`, the owner loop.  Effect on the PI
FSEAM fixture: 2159 → 2174 poles (seed 0), 2156 → 2162 (seed 42).  FSEAM.2
two-way match stays 97.2 % / 97.0 %, and Δcount goes 3 → 12 (gate ≤ 10 %).
There is no PI oracle.  The evidence that the added poles are real is
indirect: had all 15 seed-0 additions been unmatched in seed 42,
match(0→42) would have fallen to about 96.5 %, and it did not move.  ℘ is
unaffected (0 multi-window owner-absent groups).

**6.2 Chaining / coupled default.**  The claim is real at the kernel level.
Synthetic case: P at -5+0.2i (core 3, seen by windows 3 and 4) and Q at
-5-0.1i (core 1, seen by windows 1 and 2) are 0.3 apart.  At
`atol = 0.4` single linkage fuses them and one real pole is lost; at 0.1
both survive.  Fix: the driver default is now a fixed `0.1`, no longer
the caller's `cluster_atol` (`src/WindowedComposite.jl`,
`windowed_extract_poles`).  I did **not** add a diameter cap that throws.
The measured largest linked group on the production PI fixture is 0.164,
which is above `atol` = 0.1, so a throw-at-`atol` cap would crash
`edge_gated_windowed_poles` on the FSEAM configuration.  A cap needs its
own measurement of what those 0.16-wide groups are (see §6.5).

**6.3 Synthetic kernel tests.**  New testset WPO.0 (10 assertions) in
`test/windowed_pole_ownership_test.jl`, with answers known by construction:
- both straddle shapes on a core line (tie-break goes to window 1);
- a 4-window corner (mean exactly 0 goes to window 1);
- an owner-absent multi-window group (kept from window 1);
- an owner-absent singleton (dropped);
- the chain case at 0.1 (both kept) and at 0.4 (one kept);
- empty input;
- BigFloat input (element type preserved, same answer).

WPO.2's self-comparison assertion is removed.  WPO.3's merge_atol routing
check is kept and labelled as a *wiring* check.  I added a
decoupling wiring check at `cluster_atol = 3.0`, together with a
non-triviality assertion that the coupled and decoupled outputs differ on
this fixture (at 0.4 the ℘ fixture cannot tell them apart).

**6.4 Mutation proof** (each observed RED, then restored by copying back
and confirmed with `cmp`; line numbers are those at run time, and a
header-comment edit afterwards moved them +1):
- M3, owner-absent fallback disabled (`length(unique(ws)) ≥ 2` changed to
  `false`): WPO.0 9/10, RED at the multi-window owner-absent assertion (:56).
- M4, fallback also applied to singletons (`≥ 2` changed to `≥ 1`): WPO.0
  9/10, RED at the singleton assertion (:58).
- M5, default re-coupled (`get(extract_kwargs, :cluster_atol, 0.1)`):
  first survived the `cluster_atol = 0.4` wiring check, which is why that
  check was replaced; with the 3.0 version, WPO.3 is RED (10/11, :146).

Final runs: `test/windowed_pole_ownership_test.jl` WPO.0 10/10 + WPO 11/11;
`test/field_seam_test.jl` 8/8, 7/7, 4/4.

**6.5 Open.**  Groups wider than `atol` on PI (0.164) mean single linkage
does chain there.  Before any cap is chosen, someone should measure whether
those groups are one pole or two.  The same-window duplicates (§5) are
unchanged.
