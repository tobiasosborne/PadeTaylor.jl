# Bug sweep B3 — VectorPathNetwork + Stage2 + WedgeStep

**Date**: 2026-06-01 · **Auditor**: read-only numerical-methods audit
(Opus 4.8 1M). **Scope**: shared-Q walk, single-nearest Stage-2 fill,
pole-capped step. **Mode**: Read/Grep only; no `julia`, no edits.

## Area

`src/VectorPathNetwork.jl` (1005 LOC), `src/VectorPathNetworkStage2.jl`
(461 LOC), `src/VectorWedgeStep.jl` (576 LOC). Cross-read against the
scalar twins `src/PathNetwork.jl`, `src/StepControl.jl`,
`src/VectorStepper.jl`, `src/VectorStepControl.jl`, `src/SharedPade.jl`,
and the cited references.

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:156-166`
  — the §3.1 two-stage walk: 5-direction wedge, "smallest |u|" selection
  from the *same* Padé, "within a distance of h" halt, nearest-visited
  restart, Stage 2 = single Taylor/Padé step.
- `docs/adr/0025-headline-figure-re-resolution.md:70-89` (the three
  levers), `:232-247` (the B1 gate `R_gate = min(s·h_JZ·h_v, 0.5·h_v·min|t*|)`,
  `s = 0.34/0.36/0.30`).
- `docs/adr/0026-vector-resilient-walk-dense-targets.md:195-360`
  (Amendments 2/3/4 — the geometric-sink diagnosis and the S2
  non-ratcheting step law `h = clamp(SAFETY·D_local, h_min, h_max)`),
  `:519-566` (Amendment 9 — the `default_tol` 1e-14 vs walk-1e-8 hygiene
  note; the no-extrapolation contract is the real coverage cap).
- `docs/worklog/063-zoom-figure-iteration-arc.md:13-23, 60-66, 147-191`
  — the single-nearest-node Stage-2 fill is the documented
  discontinuity/artefact source; the principled fix (K-nearest
  cross-validation) is deferred.
- `src/StepControl.jl:240-262` (`step_pade_root`, the FW §3.1 *forward
  directional* pole-distance heuristic).
- `src/VectorStepper.jl:224-270` (`vector_pade_step_with_pade!` — the
  rescale-by-h^k → shared-Q → eval-at-t=1 primitive; the `|Q(1)| ≤ √ε·‖Q‖`
  pole-hit `DomainError`).
- `src/VectorStepControl.jl:158-215` (`vector_step_jorba_zou`).
- `src/SharedPade.jl:163-268` (degree reduction + QR-reweighting,
  `default_tol`).
- `test/vector_path_network_test.jl` VPN.3.1–3.5, VPN.4.1–4.3 (the
  pinned invariants for the adaptive law, the selector, determinism).

## Findings

### [LOW] Single-nearest-node Stage-2 Voronoi fill produces intermittent seam/`Q`-zero artefacts (documented, honest-by-construction, NOT a transcription bug)

- **Location**: `src/VectorPathNetworkStage2.jl:410-439` (`_stage2_fill`
  per-grid-point loop); call site `src/VectorPathNetwork.jl:896-899`.
- **Ground truth**: FW-md:166 ("step from the end of a path to its
  nearby node points on the fine grid can be done based on a *single*
  Taylor/Padé expansion") — FW itself prescribes single-nearest-node,
  *no* validity test. ADR-0026:27 confirms "there is no barycentric fill
  to port … the scalar `PathNetwork.jl:669-693` Stage-2 is *also* plain
  nearest-node single-Padé lookup." Worklog 063:13-23, 60-66 documents
  the resulting artefacts ("Voronoi-cell rays … discontinuities
  visible"; spurious "poles" where a cell's nearest node has a `Q`-zero
  nearby) and 169-191 defers the K-nearest cross-validation fix.
- **Code behavior**: each grid point is assigned its single nearest
  visited node (`nearest_visited`) and evaluated as `Pᵢ(t)/Q(t)` of
  *that one* node's canonical shared-Q, gated by `R_gate`. Adjacent grid
  cells whose nearest node flips across a Voronoi boundary jump from one
  Padé to another; where a node's `Q(t)` has a root inside its own
  `R_gate` disc the cell renders as a spurious pole.
- **Mechanism (intermittent)**: data-dependent. The artefact appears
  only where (i) a Voronoi seam separates two nodes whose discs do not
  overlap (≈85 % of seams per the B4 audition, `VectorPathNetworkStage2.jl:155-160`),
  or (ii) a node's robust-Padé returns a spurious `Q`-zero inside its
  gated disc. Neither is always-on; both depend on the walk tree the
  RNG/target order produced. This is exactly the "intermittent
  discontinuity in computed solutions" symptom — but it is a *design
  corner FW also has*, gated honest (out-of-disc ⇒ `NaN`), and explicitly
  deferred, not a mis-transcription.
- **Intermittent?**: yes. **Confidence**: 0.9 that this is the dominant
  cause of the maintainer's observed intermittent discontinuities;
  0.95 that it is **not** a transcription/sign/index bug (it is the
  faithful FW single-Padé fill plus an honest gate, per ADR-0026:27).

### [LOW] Arrival threshold `|z_cur − target| > real(visited_h[parent])` tracks the *last adaptive step*, not a fixed `h` (FW deviation, documented, mostly benign post-S2)

- **Location**: `src/VectorPathNetwork.jl:798` (loop guard), with
  `parent` updated to the freshly-pushed node at `:836`, so on iteration
  ≥2 the threshold is `h_cur` of the most-recent accepted step.
- **Ground truth**: FW-md:163 — "Once we reach within a distance of **h**
  from the target point, this path is halted," where FW's `h` is a
  *single fixed global* step (FW-md:164, "Changing h from 0.3 to 0.5,
  the case used throughout"). ADR-0026:214-216 flags precisely this:
  "the loop-termination test `|z−target| > visited_h[parent]` … inherits
  the collapsed `h`."
- **Code behavior**: with `adaptive = true` the halt radius is the
  per-node adaptive step, which varies node-to-node. Under the *old*
  geometric-sink law this forced a ~70-step crawl to within ~1e-4 of the
  target (ADR-0026:178-180, Amendment 1 measurement). Under the current
  non-ratcheting S2 law (`VectorWedgeStep._adaptive_h`,
  `VectorWedgeStep.jl:539-574`) `h` recovers to `h_max` in one step on
  leaving a dense pocket, so the crawl pathology is largely cured.
- **Mechanism (intermittent)**: where the final approach step lands in a
  dense pole pocket, the halt radius shrinks, so the walk plants the
  terminal node much closer to the target than in a sparse region —
  giving target-to-nearest-node spacing that varies with the local pole
  field. In the single-nearest Stage-2 fill this modulates the Voronoi
  cell sizes near pole-dense regions, contributing to the seam-jump
  pattern above. Not a wrong value, but a field-dependent geometry that
  compounds finding B3-1. **Intermittent?**: yes (data-dependent).
  **Confidence**: 0.5 that this materially worsens the artefacts;
  0.85 that it is a *documented intentional* deviation (per-node
  `visited_h` is the deliberate B2 design), not a transcription bug.

### [LOW] `min|t*|` over **all** `roots()` outputs is exposed to spurious near-origin shared-Q roots in both the step cap and the selector

- **Location**: `src/VectorWedgeStep.jl:355` (`_candidate_pole_disc`:
  `min(clear, h_mag*minimum(abs, rs))`) and `:550` (`_adaptive_h`:
  `D_local = h_prev*minimum(abs, rs)`); same pattern in the Stage-2 gate
  `src/VectorPathNetworkStage2.jl:327`.
- **Ground truth**: ADR-0026 Amendment 3 (`:294-300`) raised exactly
  this — "spurious near-origin shared-`Q` roots poison `min|t*|`" — as
  the open question gating the fix. Amendment 4 (`:335-340`, S1 EXP B)
  *claims* it resolved the suspicion: the cap-driving root maps to a
  z-plane pole "stable to 3 decimals — a real pole, not a Froissart
  artefact," concluding "`min|t*|` from the shared-`Q` is trustworthy."
  `SharedPade.jl:185-265` does port GGT degree-reduction + the
  padeapprox.m QR-reweighting (so Froissart doublets *are* suppressed,
  per ADR-0026:537-540).
- **Code behavior**: both selectors take the unconditional minimum over
  every root `Polynomials.roots` returns. The `CLEAR_CAP` clamp
  (`VectorWedgeStep.jl:306,355`) only bounds the disc from *above* (a far
  pole is treated as pole-free); it does **not** filter a spurious
  *small*-`|t*|` root, which drags `D_local`/the disc score down. There
  is no Froissart/z-distance pruning at this call site (the pruning lives
  downstream in `VectorPoleField.extract_poles_shared_q`, S7, applied for
  *figure* extraction, not for the step cap).
- **Mechanism (intermittent)**: `Polynomials.roots` is a companion-matrix
  eigensolve whose spurious roots are conditioning-dependent and vary
  with the jet content and the degree-reduced `m_cur` — so a near-origin
  artefact root appears only at some nodes/orders. Where it does, the
  cap collapses `h` (or the candidate's disc score), shifting which
  wedge candidate wins and where the terminal node lands — an
  order/data-dependent perturbation of the walk tree. **Intermittent?**:
  yes. **Confidence**: 0.4. The repo's own Amendment 4 measurement
  argues against it being a live defect for the cap-driving root; I flag
  it because the resolution was measured on a *specific* config and the
  no-lower-bound exposure is structurally still present.

### [INFO] `default_tol = 1e-14` in the canonical-Q build vs the walk's `tol = 1e-8` (already filed; second-order)

- **Location**: `_canonical_pade` → `vector_pade_step_with_pade!`
  (`src/VectorPathNetwork.jl:957-958`) calls `shared_denominator_pade`
  with no `tol`, defaulting to `default_tol(T) ≈ 1e-14`
  (`src/SharedPade.jl:163-165`, `:181-183`). The Stage-2 gate, by
  contrast, takes the user `tol = 1e-8`
  (`src/VectorPathNetwork.jl:898`).
- **Ground truth**: ADR-0026:540-542 already records this — "VectorStepper
  builds the canonical `Q` at `default_tol = 1e-14` rather than the
  walk's `1e-8` noise floor — a frontier-only second-order effect, filed
  low-priority, not a coverage lever."
- **Code behavior / mechanism**: a tighter build tol keeps more
  near-noise denominator coefficients, so the canonical `Q` can carry a
  slightly higher degree (more roots) than the `1e-8`-gate would deem
  meaningful — feeding finding B3-3's `min|t*|`. Order- and
  data-dependent. **Intermittent?**: weakly yes. **Confidence**: 0.3
  that it has any visible effect; documented as known.

## Areas verified correct

- **Wedge direction & sign** — `_select_wedge` computes
  `θ = goal_dir + wedge_angles[k]`, `h_step = h_mag·(cosθ + i·sinθ)`
  (`VectorWedgeStep.jl:402-403`) and `goal_dir = angle(target − z_cur)`
  (`VectorPathNetwork.jl:808`). Byte-for-byte the scalar convention
  (`PathNetwork.jl:531, 943-945`) and FW-md:160-162 ("aiming straight at
  the goal … 22.5 and 45 to either side"). No sign flip, no
  conjugate/adjoint confusion on the complex step.

- **The S2 non-ratcheting adaptive step law** — `_adaptive_h`
  (`VectorWedgeStep.jl:539-574`) implements
  `h = clamp(SAFETY·D_local, h_min, h_max)` with
  `D_local = h_prev·min|t*|` and **no** `h_prev` multiplier on the target,
  exactly as ADR-0026 Amendment 4 (`:343-360`) prescribes. The
  `:step_collapse` throw correctly tests the *unclamped* `h_target < h_min`
  (`:565`), not the clamped `h` — pinned by VPN.3.4 (`:732-744`) and
  VPN.3.3(b) load-bearing recovery (`:666-681`). The geometric-sink
  pathology is genuinely gone by construction.

- **`step_pade_root` reuse claim is a *deliberate, documented*
  deviation, not a silent break** — the special-focus suspicion was that
  `:max_q_root` should "reuse `StepControl.step_pade_root` unchanged."
  It does **not**, and that is correct: `step_pade_root`
  (`StepControl.jl:240-262`) is a *forward, directional* projection
  (`t = Re((r − z)·conj(unit)); keep t>0`), whereas `_adaptive_h` uses a
  *direction-agnostic* `min|t*|`. `VectorWedgeStep.jl:94-99, 514-521` and
  ADR-0026:418 spell out why the directional projection mis-couples with
  the `:max_q_root` selector (which already steers around the pole). The
  divergence is principled and recorded; no transcription error.

- **B1 true-radius gate formula** — `_validity_radius`
  (`VectorPathNetworkStage2.jl:288-332`) computes
  `R_gate = min(s(tol)·h_JZ·h_v, 0.5·h_v·min|t*|)` with the jet rescaled
  *once* by `h_v^k` before `vector_step_jorba_zou`. Matches ADR-0025:240-243
  exactly; `_safety_factor` 0.34/0.36/0.30 boundary logic
  (`:244-254`) maps each calibration tol to itself with correct
  nearest-rounding at the −7/−9 midpoints; default `tol=1e-8 → 0.36`.
  No double-rescale (the cached `visited_jets` are confirmed *raw*,
  built by `vector_taylor_coefficients` at `VectorPathNetwork.jl:959`).

- **Jet rescaling and the canonical-Q centre/scale** — `_canonical_pade`
  and `_candidate_pole_disc` build the canonical store with the *same*
  real step `C(h_cur)` at the *same* landed node
  (`VectorPathNetwork.jl:957-958` vs `VectorWedgeStep.jl:345-346`); the
  Stage-2 fill rescales `t = (z_f − z_v)/h_v` with `h_v = real(visited_h)`
  (`VectorPathNetworkStage2.jl:413, 436`). Centre and scale are
  consistent end-to-end; a candidate that scores finite in
  `_candidate_pole_disc` cannot crash the driver's `_canonical_pade`
  (identical operation), so there is no selector/driver crash gap.

- **`vector_step_jorba_zou`** — `k ∈ {p−1, p}`, `(ε/nrm)^(1/k)`, the
  TI.jl second-stepsize fallback over `1≤j≤p−2`, and the `ε` resolution
  match the scalar `step_jorba_zou` (`StepControl.jl:162-204`) and
  ADR-0021 with `|c_k|→vnorm(c_k)`. The `d=1` reduction is exact
  (`VectorStepControl.jl:93`). No off-by-one in the coefficient index
  (`coefs[k+1]`/`_coef_vector(jets,k,d)`).

- **`vector_pade_step_with_pade!` pole-hit guard** — `|Q(1)| ≤ √ε·‖Q‖`
  ⇒ `DomainError` (`VectorStepper.jl:256-263`); the relative `√ε·‖Q‖`
  threshold is principled (distinguishes a true large value from
  roundoff-near-zero). The eval is at `t = 1` over the shared `Q(1)`,
  the documented shared-denominator keystone. Correct.

- **`_nearest_visited` lexicographic tiebreak** — vector
  (`VectorPathNetwork.jl:921-934`) is byte-identical to scalar
  (`PathNetwork.jl:882-894`): Euclidean min, tie on `Re` then `Im`. The
  `skip_covered`/exact-coincidence tests (`VectorPathNetwork.jl:764-770`)
  use symmetric `abs(·)` so no order asymmetry. Selection tie-break
  `(primary, −|wedge_angle|)` (`VectorWedgeStep.jl:439-443`) is strict-`>`
  hence deterministic (first index wins on a full tie). Determinism is
  pinned by VPN.4.1/4.2 (`test:780-831`).

- **`extrapolate` gate** — default `false` ⇒ `abs(z_f − z_v) > R_gate`
  yields `NaN+NaN·im` (`VectorPathNetworkStage2.jl:429-432`); `true`
  skips the gate. Matches ADR-0015 / the scalar `> h_v` gate
  (`PathNetwork.jl:685`) in semantics (strict `>`), differing only in the
  honest radius used (`R_gate` vs `h_v`) per ADR-0025 Lever 1. The
  worklog-063 zoom NaN cells are this gate firing correctly, then papered
  over by a *figure-script* BFS fill — `src/` returns honest `NaN`.

- **Resilient `:skip` driver** — only `e isa VectorWalkError` is caught;
  any other exception is `rethrow()`n (`VectorPathNetwork.jl:856`),
  honouring Rule 1; `z_stuck`/`residual_dist` are read from the live
  `z_cur`, never parsed from the message. Partial nodes are retained
  (push happens only on accepted steps). No silent swallow.
