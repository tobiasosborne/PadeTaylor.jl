# Bug sweep B4 — BranchTracker + SheetTracker (winding / branch-cut / sheet index)

Date: 2026-06-01. Auditor: read-only subagent. Scope: `src/BranchTracker.jl`,
`src/SheetTracker.jl`, plus the PathNetwork call sites that consume them.

## Area

Riemann-sheet bookkeeping for PVI in the ζ- and η-planes:

- `SheetTracker.winding_delta / accumulate_winding / sheet_index` —
  signed-angle primitives and integer sheet index.
- `SheetTracker.pVI_transformed_rhs / pVI_eta_transformed_rhs` —
  the ζ- and η-plane transformed RHS closures (FFW 2017 eq. 3 / eq. 5).
- `SheetTracker.pVI_z_to_η / pVI_η_to_z` — composed exponential
  coordinate transforms.
- `BranchTracker.segment_crosses_cut / any_cut_crossed /
  step_sheet_update / resolve_cut_angles` — walker-side cut-crossing
  predicate and per-step sheet-counter update.
- Their consumption in `PathNetwork.path_network_solve` (Stage-1 parent
  selection + sheet accumulation; Stage-2 sheet-aware lookup).

## References checked

- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:35`
  — PVI equation (z-plane), source of the transforms.
- `…/FFW2017….md:141` — `z = e^ζ` maps `z=1` to `ζ = 2πik`; fixed
  singularities of eq. (3).
- `…/FFW2017….md:144` — eq. (3), ζ-plane transformed RHS (verbatim).
- `…/FFW2017….md:149` — eq. (4), η-lattice `η = log|2πk| + i·arg(2πik)`,
  `|k| ≥ 1`.
- `…/FFW2017….md:154` — eq. (5), η-plane transformed RHS (verbatim);
  note the `-1` chain-rule artefact inside the `(dv/dη)` bracket.
- `…/FFW2017….md:157` — "the region Re η < log 2π is branch point-free".
- `…/FFW2017….md:178` — circumambulation spec: "none of the paths
  overstep the branch cuts … on a given sheet"; CCW/CW direction
  semantics.
- `…/FFW2017….md:180-189` — sheet parametrisation `(θ_k, φ_k)` index
  ranges; CCW revolution ↔ next sheet.
- `docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md:82-115` —
  cut model, predicate spec (`t ∈ (0,1)`, `s > 0`), refuse/cross mode,
  winding-wrap formula.
- `src/PathNetwork.jl:481` (shuffle), `:500-501` (Euclidean parent
  selection), `:555-558` (refuse-mode filter), `:624-633` (cross-mode
  sheet accumulation), `:669-693` (Stage-2 sheet-aware lookup),
  `:882-922` (`_nearest_visited` / `_nearest_visited_on_sheet`).

## Findings

### [MEDIUM] Sheet-blind Euclidean parent selection makes cross_branch sheet labels shuffle/order dependent

- **Location**: `src/PathNetwork.jl:500-501` (`_nearest_visited`) feeding
  `src/PathNetwork.jl:627-629` (`step_sheet_update`), with targets in
  shuffled order at `src/PathNetwork.jl:481`.
- **Ground truth**: FFW md:178 — paths reach a new sheet only by
  running "in a counterclockwise direction around ζ = 0" through the
  cut, and "none of the paths overstep the branch cuts … on a given
  sheet"; the sheet of a node is determined by the analytic-continuation
  path that physically reached it. ADR-0013:104-112 specifies the sheet
  counter is the parent counter plus `sign(Δθ_k)` for each crossed cut.
- **Code behavior**: each new Stage-1 path segment is rooted at the
  *Euclidean-nearest already-visited node* (`_nearest_visited`,
  purely `abs(visited_z[i] - target)` with a deterministic lexicographic
  tie-break; no sheet term). The new segment then inherits that parent's
  sheet vector and accumulates only the crossings *along the new
  segment*. Because targets are processed in `shuffle(rng, ...)` order
  (line 481), which visited node is nearest when a given target is
  reached is order/seed dependent. When a target sits near a cut, the
  nearest visited node can be on a different sheet than the analytic
  continuation that should reach the target; the segment then gets
  rooted on the wrong sheet and the whole sub-path is mislabelled.
- **Mechanism (intermittent discontinuity)**: Stage-2 reconstruction at
  `src/PathNetwork.jl:672` uses `_nearest_visited_on_sheet`, which
  filters visited nodes by *exact* `visited_sheet[i] == grid_sheet[i]`.
  A mislabelled Stage-1 node either (a) becomes the chosen interpolation
  node for a grid point on its (wrong) sheet, injecting a value from the
  adjacent sheet — a jump across the cut, or (b) is excluded from the
  correct sheet's pool, leaving the grid point to a farther node or NaN.
  Either way the rendered solution shows an intermittent discontinuity
  near branch cuts that moves or appears/disappears when `rng_seed` (or
  the node set / target order) changes. This is order/RNG dependent, not
  always-on — the signature the maintainer describes.
- **Intermittent?**: Yes — RNG-seed / target-order dependent.
- **Confidence**: 0.55. The mechanism is real and the code path is
  exactly as described, but ADR-0013:255-261 already lists "A5
  sheet-aware Stage-2" as in-scope and acknowledges Stage-1 parent
  selection is Euclidean; whether this rises to "the" bug versus an
  acknowledged limitation depends on whether cross_branch figures are
  actually mislabelling. Root cause lives in PathNetwork's parent
  selector, which consumes the (correct) BranchTracker primitives.

### [LOW] Knife-edge cut-crossing predicate has no epsilon band: crossings within float noise of a node/cut flip intermittently

- **Location**: `src/BranchTracker.jl:117-120` —
  `iszero(det) && return false`; `return 0 < t < 1 && s > 0`.
- **Ground truth**: ADR-0013:82-90 — predicate is open on `t` ("steps
  that *land on* the cut endpoint are not counted") and parallel
  segments (`det = 0`) return false. The module docstring
  (`src/BranchTracker.jl:41-51`) calls the endpoint cases "measure-zero
  numerical accident[s]".
- **Code behavior**: comparisons are exact float predicates with no
  tolerance. A step whose endpoint lands within rounding of the cut
  yields `t` straddling `1.0` (or `s` straddling `0.0`); a step nearly
  parallel to the cut yields `det` straddling `0.0`. The same geometric
  situation can therefore evaluate to "crosses" in one run and "does not
  cross" in another (different node set, accumulated rounding, or
  element type).
- **Mechanism (intermittent discontinuity)**: in refuse mode a knife-edge
  step is sometimes forbidden and sometimes allowed; in cross mode a
  crossing is sometimes counted and sometimes missed, shifting a node's
  sheet index by ±1 → the Stage-2 sheet filter (PathNetwork.jl:672) then
  jumps that grid region to a different sheet. Genuinely intermittent,
  but requires landing within ~eps of the cut/branch, which the FFW node
  sets are not aligned to do — so this is a latent fragility, not a
  demonstrated everyday cause.
- **Intermittent?**: Yes — float-rounding / node-set dependent at a
  measure-zero locus.
- **Confidence**: 0.2.

### [LOW] `Δθ == 0` exactly on a crossing yields no sheet bump (degenerate through-branch step)

- **Location**: `src/BranchTracker.jl:162-163` —
  `Δθ = winding_delta(...); out[k] += Δθ > 0 ? 1 : (Δθ < 0 ? -1 : 0)`.
- **Ground truth**: ADR-0013:108-111 — "the sheet counter for b_k is
  updated by sign(Δθ_k)". `sign(0) = 0` is consistent with the spec, but
  the spec assumes a genuine crossing always has a nonzero winding.
- **Code behavior**: if `segment_crosses_cut` is true yet
  `winding_delta` returns exactly `0.0`, no bump is applied. Geometrically
  this only happens when the step's endpoints are collinear with the
  branch (same ray) while the segment still intersects the cut ray — i.e.
  the segment passes through the branch point itself.
- **Mechanism**: a through-branch-point step would skip the sheet bump
  while still being a topological crossing; downstream sheet labels off
  by one. Self-limiting: the Padé stepper fails loud at `Q(z)=0` on the
  branch point (module docstring lines 44-49), so this state should not
  silently persist. Measure-zero and largely defended elsewhere.
- **Intermittent?**: Yes (data-dependent), but practically unreachable.
- **Confidence**: 0.15.

### [LOW] ADR winding-wrap formula text implies a `[-π, π)` interval; code uses `(-π, π]` — documentation/spec mismatch at the boundary

- **Location**: spec `docs/adr/0013-constrained-wedge-and-sheet-bookkeeping.md:105-106`
  (`Δθ_k = wrap((arg(z_new−b_k) − arg(z_cur−b_k)) + π) − π`) vs code
  `src/SheetTracker.jl:283-292` (normalises to `(-π, π]`:
  `Δθ ≤ -π → +2π`, `Δθ > π → -2π`).
- **Ground truth**: the SheetTracker docstring (`src/SheetTracker.jl:54-56`,
  275) and BranchTracker docstring (`src/BranchTracker.jl:60-63`) both
  state the convention is `(-π, π]`. FFW md:187-189 parametrises sheet 0
  as `-π < θ_0 ≤ π` (half-open, upper-closed) — consistent with the code.
- **Code behavior**: code is self-consistent and matches FFW's
  upper-closed convention. The ADR's `wrap(x+π)-π` idiom, if `wrap`
  targets `[0, 2π)`, would instead give `[-π, π)` (lower-closed) — the
  opposite endpoint. The two differ only at the measure-zero point
  `Δθ = ±π`, which the documented `|Δθ| < π` caller contract excludes.
- **Mechanism**: not a live bug given the caller contract; flagged only
  because endpoint-convention drift between spec and impl is exactly the
  class that breeds boundary bugs if the contract is ever relaxed.
- **Intermittent?**: No (would only matter at exact `Δθ = ±π`).
- **Confidence**: 0.2 (that it is a real spec/impl wording mismatch);
  ~0.0 that it currently produces wrong output.

## Areas verified correct

- **ζ-plane transformed RHS** (`src/SheetTracker.jl:155-174`) matches
  FFW eq. (3) `…/FFW2017….md:144` term-for-term: the `(1/2)(1/w +
  1/(w−1) + 1/(w−e^ζ))(dw/dζ)²` first term, the `−(e^ζ/(e^ζ−1) +
  e^ζ/(w−e^ζ))dw/dζ` second term, and the
  `w(w−1)(w−e^ζ)/(e^ζ−1)²·(α + βe^ζ/w² + γ(e^ζ−1)/(w−1)² +
  δe^ζ(e^ζ−1)/(w−e^ζ)²)` third term. No sign/factor error.

- **η-plane transformed RHS** (`src/SheetTracker.jl:206-226`) matches
  FFW eq. (5) `…/FFW2017….md:154` term-for-term with `E = e^(e^η)`,
  including the easily-dropped `−1` chain-rule artefact inside the
  `(dv/dη)` bracket (`src/SheetTracker.jl:216`) and the `e^(2η)`
  multiplier on the third term (`src/SheetTracker.jl:222`). No
  transcription error.

- **z↔η coordinate transforms** (`src/SheetTracker.jl:244-265`): forward
  `vp = z·ζ·up` and inverse `up = vp/(z·ζ)` are the correct chain-rule
  for `ζ = log z`, `η = log ζ` (`dζ/dη = ζ`, `dz/dζ = z`), consistent
  with FFW md:146/151's `u(z) = w(ζ) = v(η)`. Round-trip is algebraically
  exact.

- **`winding_delta` sign + wrap** (`src/SheetTracker.jl:283-292`):
  verified by hand on the canonical negative-real-axis crossing — moving
  from `−1−0.5i` to `−1+0.5i` (up through the cut, left of origin) gives
  raw `Δθ ≈ +5.36`, wrapped to `≈ −0.93 < 0`, i.e. clockwise, which is
  geometrically correct (the short way around the back of the origin is
  CW). Wrap is symmetric and lands in `(-π, π]`: `Δθ = −π → +π`, `Δθ = +π`
  stays. Matches FFW md:187-189 upper-closed sheet convention.

- **`segment_crosses_cut` linear algebra** (`src/BranchTracker.jl:112-120`):
  derived the real 2×2 system for `t·d − s·u = δ` and confirmed via
  Cramer's rule that `det = imag(d·conj(u))`, `t = imag(δ·conj(u))/det`,
  `s = imag(δ·conj(d))/det` are the exactly-correct closed forms (the
  `imag(a·conj(b))` idiom is the 2D cross product `Re a·Im b − Im a·Re
  b` with the right sign throughout). The `s > 0` half-line gate matches
  ADR-0013:88-90, and mutation B1 (dropping `s > 0`) is shown to bite 5
  assertions (`test/branch_tracker_test.jl:174-180`). No conjugate-vs-
  transpose confusion: all three uses of `conj` are correct for the 2D
  cross-product idiom, not adjoints of a matrix.

- **`step_sheet_update` direction** (`src/BranchTracker.jl:156-167`):
  CCW crossing (`Δθ > 0`) bumps `+1`, CW (`Δθ < 0`) bumps `−1`,
  consistent with FFW md:180 ("counterclockwise revolutions … the sheets
  … are parametrized by" increasing index) and ADR-0013:108-111. The
  gate `if segment_crosses_cut(...)` before the bump is present (mutation
  B3 dropping it bites 52 assertions, `test/branch_tracker_test.jl:188-194`),
  and the function returns a fresh `collect(Int, sheet_old)` rather than
  mutating in place (no buffer aliasing; `test/branch_tracker_test.jl:127`
  asserts `sheet_noop !== sheet_in`).

- **`sheet_index` rounding** (`src/SheetTracker.jl:320`): `round(Int,
  total/2π)` gives the FFW-consistent +1/-1 per revolution. Note it is
  only used in the after-the-fact analysis path (it is NOT on the
  walker's `step_sheet_update` path — confirmed by grep: no PathNetwork
  call site). Banker's-rounding asymmetry exists only at exact
  half-integer winding (`total = ±π, ±3π, …`), a measure-zero locus, and
  does not feed the solver. Not a solver-output intermittency source.

- **`resolve_cut_angles`** (`src/BranchTracker.jl:186-197`): scalar
  broadcast and length-validated tuple path are correct; mismatched
  length throws `ArgumentError` per CLAUDE.md Rule 1
  (`test/branch_tracker_test.jl:161-163`). Empty-branch no-op verified
  (`test/branch_tracker_test.jl:165-166`).

- **Refuse-mode filter** (`src/PathNetwork.jl:555-558`, `:736-745`) uses
  `z_cur` (pre-step) and the candidate `z_n` correctly; cross-mode
  accumulation (`src/PathNetwork.jl:624-633`) is invoked with `z_cur,
  z_new` BEFORE `z_cur` is advanced at line 636 — correct ordering, no
  off-by-one on the segment endpoints. Root sheet is seeded from
  `initial_sheet` (`src/PathNetwork.jl:438, 476-477`) with a length-check
  against `branch_points` (`:433-436`).

- **η branch-point-free region claim**: `src/SheetTracker.jl:204`
  ("region Re η < log(2π) is free of them") matches FFW md:157 verbatim;
  the lattice formula in the docstring (`src/SheetTracker.jl:41-42, 79`)
  matches FFW eq. (4) md:149.
