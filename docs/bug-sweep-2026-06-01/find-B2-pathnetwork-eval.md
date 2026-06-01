# Bug sweep B2 — PathNetwork dense-output / Stage-2 fill

## Area

`src/PathNetwork.jl` — the Stage-2 fine-grid fill (lines 660-693), the
public per-point accessors `eval_at` (1086-1102) and `eval_at_sheet`
(1026-1060), the nearest-visited assignment `_nearest_visited`
(882-894) and `_nearest_visited_on_sheet` (905-922), the
distance-`≤ h` cutoff (doc lines 30-32; code line 685 / 1095 / 1053),
and node-local Padé evaluation (`t = (z_f - z_v)/h_v`, `up = P'(t)/h_v`).

Cross-read for ground truth: `src/PadeStepper.jl` Padé build &
evaluation (`pade_step_with_pade!` 300-327, `_evaluate_pade` 367-382,
`_evaluate_pade_deriv` 397-409), the Vector analogue
`src/VectorPathNetworkStage2.jl` (the principled `_validity_radius`
gate at line 429), and `src/VectorPathNetwork.jl:921` (`_nearest_visited`
verbatim copy).

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:155-166`
  — Stage-1 + Stage-2 spec ("within a distance of h", "store ... the
  Padé coefficients", "single step right up to each nearby node").
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:368`
  — steepest-descent direction `θ = arg(−u(z₀)/u'(z₀))`, "choose the
  edge of the wedge closest to this steepest descent direction".
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:208`
  — "low-level ridges in flat areas ... different u entries ... obtained
  following entirely different paths" (the documented path-dependent
  artefact, NOT a transcription bug).
- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:60`
  — "until the target node is within a distance of |h|".
- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:62`
  — Stage-2: "Padé steps are taken from each ζ_i to the points on the
  fine grid to which it is the closest point among the ζ_i" (NO disc
  check in FFW).
- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:93,95`
  — adaptive: "stored at the point" is `|qh|`; Stage 2 fills every fine
  cell from its nearest stored point.
- `docs/unified_path_network_spec.md:48-55,62-70` — termination
  `|z_cur − z_t| ≤ h`; Stage-2 `t = (z_f − z_v)/h`; lexicographic
  tiebreak adopted as a local decision.
- `docs/adr/0004-path-network-architecture.md:54-58` — "find nearest
  visited node with coverage `|z_f - z_v| ≤ h`; evaluate stored Padé at
  `t = (z_f - z_v)/h`. ... Gap → store NaN; no silent extrapolation."
- `docs/adr/0015-extrapolate-stage-2.md` — the disc-radius gate vs
  FFW's gate-less fill; both documented and faithful to their refs.
- `docs/worklog/014-pathnetwork-symmetry-debug.md` — Bug 1
  (shuffle-induced asymmetry) and Bug 2 (y=0 row discontinuity),
  both attributed to path-dependence, not transcription.
- `docs/worklog/063-zoom-figure-iteration-arc.md:35-37,61-63,171-172`
  — the maintainer's own analysis that single-nearest-node Stage-2 fill
  produces Voronoi-cell discontinuities and perimeter spikes; a known
  structural limitation with a filed K-nearest follow-up.

## Findings

### [LOW] Steepest-descent wedge selector does not unwrap the angular
difference across the ±π branch cut

- **Location**: `src/PathNetwork.jl:994-997` (`_select_candidate`,
  `:steepest_descent` branch).
- **Ground truth**: FW 2011 md:368 — "we choose the edge of the wedge
  closest to this steepest descent direction." Closeness is on the
  circle of directions; the correct metric is the wrapped angular
  distance `min(|Δ|, 2π − |Δ|)`.
- **Code behavior**: line 997 selects
  `argmin(abs(θ_sd - off) for off in offsets)` with
  `offsets = goal_dir + wedge_angles` and `θ_sd = angle(-u/up) ∈ (-π,π]`,
  `goal_dir = angle(target - z_cur) ∈ (-π,π]` (line 531). The raw
  difference `θ_sd - off` is not reduced mod 2π, so when `θ_sd` and the
  wedge rays straddle the ±π cut (e.g. `θ_sd = +3.0`, a ray at `−3.0`),
  the true 0.28-rad-apart pair scores 6.0 and a genuinely far ray can
  win.
- **Mechanism**: only on the non-default `step_selection = :steepest_descent`.
  When a walk steers near direction ±π, an adjacent grid cell whose
  `θ_sd`/`goal_dir` lands on the other side of the cut can pick a
  different wedge ray, sending its sub-walk down a different branch →
  different accumulated path → a cell-to-cell jump in `grid_u`. Data-
  dependent (only near the cut), hence intermittent.
- **Intermittent?**: Yes (only when directions approach ±π; default
  `:min_u` path is unaffected).
- **Confidence**: 0.6 that the wrap is missing as described; ~0.3 that
  it is the maintainer's symptom, since the default selector is `:min_u`.

### [LOW] Stage-2 disc cutoff uses the step length `h_v`, not the
Padé's true validity radius — spurious near-pole spikes possible inside
`|t| ≤ 1`

- **Location**: `src/PathNetwork.jl:685` (and the identical gate in
  `eval_at` line 1095 and `eval_at_sheet` line 1053).
- **Ground truth**: ADR-0004:54-58 + spec §2 (md:65) define the gate as
  `|z_f - z_v| ≤ h` and evaluate at `t = (z_f - z_v)/h`. The Vector
  re-implementation gates instead on a *computed* validity radius:
  `src/VectorPathNetworkStage2.jl:429` `if abs(z_f - z_v) > R_gate`
  with `R_gate = _validity_radius(jet, denominator, h_v, tol)`
  (lines 419-427).
- **Code behavior**: the scalar fill admits any cell with
  `|z_f - z_v| ≤ h_v`. The stored canonical Padé was built with a real
  step `h_v` after a min-|u| selection that deliberately approaches a
  pole; its denominator `Q(t)` can have a zero at `|t| < 1` (a real
  pole nearer than `h_v`). Evaluating `_evaluate_pade` at a `t` near
  that zero returns a large finite value (or throws `DomainError` only
  on an *exact* zero, line 376) rather than NaN.
- **Mechanism**: a fine-grid cell whose nearest node's Padé has a
  denominator root just inside the unit disc gets a spike; the
  neighbouring cell assigned to a *different* node does not — a
  discontinuity at the Voronoi boundary that appears only where a node's
  local Padé happens to carry an in-disc pole. This is exactly the
  worklog-063 "perimeter / Voronoi-ray" artefact. It is a documented
  structural limitation (ADR-0015; worklog 063:171-172), not a
  transcription error, and the scalar/Vector divergence is intentional.
- **Intermittent?**: Yes (per-node, depends on whether that node's Padé
  has an in-disc denominator root).
- **Confidence**: 0.5 that this materially worsens discontinuities; 0.15
  that it is a *bug* (it is documented behaviour, kept deliberately for
  the FW-faithful scalar driver).

## Areas verified correct

- **Node-local variable shift and scale are consistent.** Stage-2 uses
  `t = (z_f - z_v)/h_v` (line 689) and `up = P'(t)/h_v` (line 691); the
  stored Padé is built by `_local_pade` with the *same* real `h_step`
  (line 605) that is pushed to `visited_h` (line 618). The Padé is
  constructed in the rescaled unit variable in
  `pade_step_with_pade!` (`PadeStepper.jl:311-320`:
  `c̃_k = h^k c_k`, evaluate at `t=1`, `up = P'(1)/h`), so `t=1` ⇔
  `z_v + h_v`. The derivative `/h_v` is the correct chain-rule factor.
  Matches ADR-0004:57 and spec §2 (md:68). No off-by-one between
  `visited_pade`/`visited_z`/`visited_h` — all five arrays are pushed in
  lockstep (root lines 470-474; per-step lines 614-618).

- **The distance cutoff is the correct strictness.** Stage-1 termination
  `while abs(z_cur - target) > term_dist` (line 523) halts at `≤`,
  matching FW md:164 / FFW md:60 "within a distance of |h|" and spec
  §1:48-50. Stage-2 `if !extrapolate && abs(z_f - z_v) > h_v → NaN`
  (line 685) keeps `≤ h_v` inside the disc, matching ADR-0004:57-58 and
  spec §2:65. `eval_at` (1095) and `eval_at_sheet` (1053) use the
  identical predicate, so the public accessors agree byte-for-byte with
  the stored `grid_u`/`grid_up` — no accessor-vs-grid seam.

- **`_nearest_visited` tie-break is deterministic and array-order
  independent.** (lines 882-894.) The loop maintains the invariant
  "`idx` = lexicographically-smallest (Re,Im) node among those at the
  current minimum distance": a strictly-closer node resets `idx`; an
  equidistant node replaces `idx` only when lexicographically smaller
  than the *current* `idx`. Worked both array orders for an equidistant
  triple — same result. Matches the spec §12 gap-1 lexicographic
  tiebreak (`unified_path_network_spec.md:52-55`). Exact float ties are
  essentially never hit; near-ties resolve at the `d < best` step
  deterministically for fixed inputs. The run-to-run intermittency lives
  in the shuffled *visited set* (worklog 014 Bug 1, the
  `shuffle(rng, ...)` at line 481), NOT in this comparator.

- **`_nearest_visited_on_sheet` matches `_nearest_visited`.** (905-922.)
  Same comparator, restricted to nodes with `visited_sheet[i] == sheet`,
  `idx=0` sentinel for "no matching-sheet node" → fail-soft NaN at the
  caller (lines 674-678). Initialisation `idx=0, best=Inf` with the
  `idx == 0 ||` first-match guard is correct; the tie-break references
  `visited_z[idx]` only after `idx` is set.

- **Schwarz-reflection mirror has the correct conjugation.** (835-844.)
  For `imag(z) < 0`, `grid_u[i] = conj(sol_upper.grid_u[j])` and
  `grid_up[i] = conj(sol_upper.grid_up[j])` with `j` indexed by the
  upper canonical rep `complex(real(z), abs(imag(z)))`. `u(z̄) = ū(z)`
  gives `u(z) = conj(u(z̄))` and `u'(z) = conj(u'(z̄))`; both correct.
  On-axis cells (`imag == 0`) take the non-conjugated else branch, also
  correct. The `< 0` strict comparison is right (on-axis is its own
  mirror). This path is opt-in and validates the IC-on-real-axis
  precondition (789-803).

- **Wedge step direction and min-|u| selection match FW.** Step
  `h·exp(i(goal_dir + θ))` (lines 944-945) with the symmetric default
  wedge `[-π/4,-π/8,0,π/8,π/4]` (line 153) = FW md:158 ±22.5°/±45°.
  `:min_u` selects `argmin(abs(e[2]))` over candidate `u_new`
  (line 992) = FW md:159 "new u value smallest in magnitude". Steepest
  descent uses `arg(-u/up)` (line 995) = FW md:368 verbatim (the
  separate wrap concern is the LOW finding above).

- **`PathNetworkSolution` field order is consistent across all three
  constructors** (struct 205-217; 10-arg fill 695-697; 11-arg
  `_attach_diagnostics` 724-728). `grid_u` and `grid_up` are never
  swapped; `visited_pade`/`visited_h`/`visited_parent` keep their slots.
