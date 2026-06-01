# Bug sweep B1 — PathNetwork Stage-1 walk

## Area

FW 2011 §3.1 path-network Stage-1 walk in `src/PathNetwork.jl`:
the 5-direction wedge, `:min_u` magnitude-minimization, the
`:steepest_descent` alternative, the seeded target shuffle, the
tree-chaining (`parent_idx` / nearest-visited start), the canonical
per-node Padé storage, and the documented default
`enforce_real_axis_symmetry = false` Schwarz-asymmetry behaviour.

Read-only audit (no `julia`, no tests run). All claims diffed against
the cited reference lines.

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:153-168`
  — Stage-1 spec: coarse grid, random target selection (line 158:
  "Initially, one of the nodes was randomly selected"; line 164: "We
  next choose, again randomly, a remaining node ... starting from the
  closest (discrete) previously visited location"), 5 directions
  ("one aiming straight at the goal and the remaining ones aiming in
  directions 22.5 and 45 to either side"), `:min_u` ("The step is then
  taken for which the new u value is smallest in magnitude among the
  five choices"), halt within distance `h`, store u/u'/Padé per point.
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:352-368`
  — §5.3/§5.4.1 step-direction discussion; line 368 is the verbatim
  `:steepest_descent` spec: "θ = arg(−u(z₀)/u'(z₀)). If this direction
  falls inside the wedge, we accept it as the direction of the next
  step. If not, we choose the edge of the wedge closest to this
  steepest descent direction as the step direction."
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:204-208`
  — FW's "low-level ridges in flat areas" caveat (path-dependence in
  smooth regions); line 399-401 §5.5 "Tuning of random path selection".
- `docs/adr/0004-path-network-architecture.md:48-77` — Stage-1/Stage-2
  decisions, wedge default, step-selection (`:min_u` default,
  `:steepest_descent` "clip to wedge"), coverage `|z_f - z_v| ≤ h`.
- `docs/worklog/008-pathnetwork-tuning.md` — the wedge-vs-canonical-Padé
  bug (found + fixed + Mutation E proven); current code at
  `PathNetwork.jl:605` is the fixed form.
- `docs/worklog/014-pathnetwork-symmetry-debug.md` — the shuffle-induced
  asymmetric-tree behaviour (Bug 1, 4-5 orders of magnitude at
  conjugate cells) and the smooth-region wedge-flip ridges (Bug 2);
  filed as `padetaylor-dtj`, cured only by the opt-in symmetry kwarg.
- `src/PadeStepper.jl:300-327, 367-409` — `pade_step_with_pade!`,
  `_evaluate_pade`, `_evaluate_pade_deriv` (the step consumed by the
  wedge + the canonical Padé).

## Findings

### [MEDIUM] Angle-wraparound in `:steepest_descent` wedge snapping picks the wrong ray near goal directions ≈ ±π

**Location**: `src/PathNetwork.jl:994-997` (`_select_candidate`,
`:steepest_descent` branch).

```julia
θ_sd = abs(up_cur) > 0 ? angle(-u_cur / up_cur) : T(goal_dir)
offsets = (T(goal_dir) + T(θ) for θ in wedge_angles)
return argmin(abs(θ_sd - off) for off in offsets)
```

**Ground truth (cited)**:
`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:368`
— "θ = arg(−u(z₀)/u'(z₀)). If this direction falls inside the wedge,
we accept it ... If not, we choose the edge of the wedge closest to
this steepest descent direction." The intended quantity is the
**angular** distance between θ_sd and each candidate ray.

**Code behavior**: `θ_sd = angle(-u/up) ∈ (-π, π]`. Each ray angle is
`off = goal_dir + θ`, with `goal_dir = angle(target - z_cur) ∈ (-π, π]`
(PathNetwork.jl:531) and `θ ∈ [-π/4, π/4]`, so `off ∈ (-5π/4, 5π/4]`.
The selection uses the raw real difference `abs(θ_sd - off)`, which is
**not** reduced modulo 2π. When the goal points near the branch cut
(`goal_dir ≈ ±π`) and θ_sd lands on the opposite side of the cut, the
naive difference is off by ~2π and `argmin` selects the *wrong* wedge
ray (often the opposite edge of the 90° wedge). Worked example:
`goal_dir = 3.0`, rays `off = {2.215, 2.607, 3.0, 3.393, 3.785}`,
`θ_sd = -3.0` (physically ≡ +3.283, i.e. essentially the +π/8 ray);
`abs(θ_sd - off)` is monotonically increasing across the 5 rays, so
`argmin` returns ray 1 (`-π/4`), the opposite edge from the correct
`+π/8`/`+π/4` ray.

The correct form snaps on the wrapped angular distance, e.g.
`abs(rem(θ_sd - off, 2π, RoundNearest))` (or `mod2pi`-normalising the
difference into `(-π, π]`). FW's "accept θ_sd if inside the wedge,
else the closest edge" is additionally not implemented (the code
always snaps to one of the 5 discrete rays) — that piece is a
documented simplification (ADR-0004:68 "clip to wedge"), but the
wraparound makes even the snapping incorrect near ±π.

**Mechanism (intermittent discontinuity)**: only triggers when (a)
`step_selection = :steepest_descent` (opt-in; `:min_u` is the FW
default and is angle-free), AND (b) the local goal direction is near
±π — i.e. walks toward the negative-real / left half-plane, common
when filling a 2-D Cartesian grid that spans negative x, AND (c) θ_sd
straddles the cut. Under those conditions a single step jumps to the
wrong wedge ray, the walk diverges from its neighbours, and the
Stage-2 value at that node (and grid cells nearest it) differs sharply
from adjacent cells — a cell-to-cell discontinuity that appears and
disappears as the goal direction sweeps through ±π across the grid.

**Intermittent?**: Yes — data/geometry-dependent (goal direction near
±π plus θ_sd cut-straddle), only in the opt-in `:steepest_descent`
mode.

**Confidence**: 0.75 that the wraparound is a real defect in the
`:steepest_descent` selection geometry; the only reason it is not
higher is that this mode is opt-in and the sole test (PN.3.1) exercises
it only with right-half-plane goals (goal_dir ≈ 0), so the wrap path is
unexercised. Demonstrated by direct construction against the FW line-368
intent; not yet observed in a running solve (tests forbidden here).

### [LOW] `:steepest_descent` never honours FW's "accept θ_sd if inside the wedge" continuous-direction rule

**Location**: `src/PathNetwork.jl:994-997` (same block).

**Ground truth (cited)**:
`FW2011_painleve_methodology_JCP230.md:368` — "If this direction falls
inside the wedge, we accept it as the direction of the next step." FW
steps in the *continuous* steepest-descent direction θ_sd when it lies
inside the 90° wedge; only when θ_sd is outside does FW fall back to
the nearest wedge **edge**.

**Code behavior**: the implementation always snaps θ_sd to the nearest
of the 5 discrete pre-evaluated rays (`argmin` over the fixed
`wedge_angles`), so a θ_sd that falls strictly between two rays is
never taken; the walker steps along a discretised approximation
instead. This is the explicitly-documented simplification in
`docs/adr/0004-path-network-architecture.md:68` ("clip to wedge"),
chosen so the candidate Padés computed in `_wedge_evaluations` can be
reused. It is a fidelity gap, not a sign/index error.

**Mechanism**: would cause `:steepest_descent` to follow a slightly
different (coarser) path than FW's RKN12† column; on a smooth-region
grid this adds the same path-dependent ridge noise FW documents
(md:208). It is bounded (≤ half the inter-ray spacing per step) and
not a discontinuity source on its own. Listed for completeness because
it is a deliberate deviation from the cited equation.

**Intermittent?**: No — systematic, bounded.

**Confidence**: 0.9 that the code deviates from FW line 368 as
described (it is documented as such in ADR-0004); low severity because
it is an intended, recorded simplification of an opt-in mode.

### [LOW] Default `enforce_real_axis_symmetry = false` reproduces the FW shuffle/path-dependence intermittency — faithful to FW, not a transcription bug

**Location**: `src/PathNetwork.jl:479-481` (seeded shuffle),
`52-73` + `751-852` (Schwarz caveat + opt-in cure), Stage-2 lookup
`669-693`.

**Ground truth (cited)**:
`FW2011_painleve_methodology_JCP230.md:158, 164` — FW *prescribes* the
random target order ("randomly selected", "again randomly"); md:204-208
documents that "different u entries ... may have been obtained
following entirely different paths through the complex plane" produce
"low-level 'ridges' in flat areas". `docs/worklog/014-...:39-45` records
the empirical 4-5-order-of-magnitude conjugate-cell divergence and that
it is **shuffle-order-dependent** ("Verified with four rng_seed values
(0, 1, 7, 42): all give non-zero asymmetry, of different magnitudes").

**Code behavior**: the RNG is *seeded deterministically* —
`rng = MersenneTwister(rng_seed)` with `rng_seed` defaulting to 0
(PathNetwork.jl:480), and `targets = shuffle(rng, ...)`
(PathNetwork.jl:481). There is **no** unseeded/global RNG in this
module (grep of `src/` shows only this seeded `MersenneTwister`; the
separate `VectorPathNetwork.jl` takes an explicit `AbstractRNG` and is
out of B1 scope). So a *single* `path_network_solve` call is fully
deterministic and reproducible. The intermittent-looking discontinuity
the maintainer observes is the documented FW path-dependence: the
shuffle creates an asymmetric visited tree, and the Stage-2
nearest-visited lookup at `z` vs `conj(z)` lands on non-conjugate
nodes. This is inherent to the faithfully-ported FW algorithm, and is
*cured* only by the opt-in `enforce_real_axis_symmetry = true`
(PathNetwork.jl:64-73, 770-852).

**Mechanism**: per-cell value depends on which tree branch reached
that cell's nearest visited node; neighbouring cells can be served by
non-conjugate / differently-rooted nodes ⇒ sharp cell-to-cell jumps
("ridges") in smooth regions and large conjugate-pair divergence. The
apparent run-to-run variability would only arise if the caller passes
*different* `rng_seed` values (or different grids / thread-unsafe `f`);
with a fixed seed it is reproducible per call.

**Intermittent?**: Yes in the spatial sense (cell-to-cell, geometry-
and shuffle-order-dependent), but deterministic per fixed-seed call.

**Confidence**: 0.9 that this is the documented FW behaviour and **not**
a transcription bug — the shuffle is seeded, the wedge is symmetric, and
the cure is correctly gated behind an opt-in kwarg. I judge this the
most likely origin of the maintainer's "intermittent discontinuities"
symptom, but the root cause is algorithmic path-dependence faithful to
FW §3.1, addressable by `enforce_real_axis_symmetry = true` (real-coeff
real-IC problems) or a fixed/consistent target order — not a sign /
off-by-one / normalization defect to be patched in the arithmetic.

## Areas verified correct

- **Wedge angle offsets** (`PathNetwork.jl:153`,
  `DEFAULT_WEDGE = [-π/4, -π/8, 0.0, π/8, π/4]`): exactly FW's "22.5°
  and 45° to either side" (`FW2011...md:158`; `π/8 = 22.5°`,
  `π/4 = 45°`) and ADR-0004:60-64. Five entries enforced
  (PathNetwork.jl:426).

- **`:min_u` selects smallest |u|** (`PathNetwork.jl:991-992`):
  `argmin(abs(e[2]) for e in evals)` where `e[2]` is `u_new` at the
  candidate's *new* point. Matches FW "the new u value is smallest in
  magnitude among the five choices" (`FW2011...md:158-159`,
  ADR-0004:66-67). Mutation-proven (test `pathnetwork_test.jl:477-483`
  Mutation A: `argmax` steers into the pole and goes RED). `argmin`
  over a `for e in evals` generator returns the linear index 1..5 of
  the minimiser, which is correctly used as `evals[idx_sel]`
  (PathNetwork.jl:561).

- **Wedge step direction** (`PathNetwork.jl:943-945`):
  `h_step = h·(cos(goal_dir+θ), sin(goal_dir+θ))` with
  `goal_dir = angle(target - z_cur)` (PathNetwork.jl:531). Correct
  rotation of the step by the wedge offset relative to the goal; the
  same `goal_dir + θ` convention is used in the steepest-descent index
  computation (PathNetwork.jl:996), so the selected index maps to the
  same physical ray that was evaluated.

- **Steepest-descent sign** (`PathNetwork.jl:995`): `angle(-u/up)`
  (descent, away from the pole), per `FW2011...md:368`
  (`θ = arg(−u(z₀)/u'(z₀))`). Mutation-proven (Mutation F: flipping to
  `angle(u/up)` steers toward the pole and goes RED,
  `pathnetwork_test.jl:515-519`, worklog 008:125-129).

- **Halt / coverage boundary** (`PathNetwork.jl:523` `> term_dist`;
  `685` `> h_v`): walk continues while strictly `> h`, halts at `≤ h`;
  Stage-2 covers `≤ h_v` and NaNs out beyond — matches "within a
  distance of h" (`FW2011...md:164`) and ADR-0004:57. Mutation-proven
  (Mutation D inverting the coverage check goes RED,
  `pathnetwork_test.jl:485-495`).

- **Tree chaining** (`PathNetwork.jl:500-513, 619, 634`): each
  per-target walk starts from `idx_v = _nearest_visited(...)` (FW
  "closest previously visited location", `FW2011...md:164`); the first
  landed node's `parent_idx = idx_v`, subsequent nodes chain off their
  immediate predecessor (`parent_idx = length(visited_z)`). Consistent
  with the `visited_parent` semantics documented at PathNetwork.jl:169-176.

- **Per-node canonical Padé** (`PathNetwork.jl:605`,
  `pade_canonical = _local_pade(..., z_new, u_new, up_new, order,
  h_step)`): stores the real-direction Padé centred at the landed node,
  not the wedge-direction `pade_sel`. This is the worklog-008 bugfix
  (Mutation E proven). Under `:fixed` (the FW default) `h_step == h_T`
  is real positive, so `t = (z_f - z_v)/h_v` in Stage 2 evaluates a
  real-direction disc — the invariant the fix restored.

- **Nearest-visited tiebreak** (`PathNetwork.jl:882-894`): lexicographic
  (Re then Im) tiebreak fires only on exact `d == best`; it yields the
  lexicographically-smallest among equal-distance ties **independent of
  insertion order** (verified by tracing both orderings), so it is
  deterministic and order-independent. Not an intermittency source
  within a fixed-seed call.

- **Threading determinism** (`PathNetwork.jl:935-955`): each `@threads`
  iteration writes only `evals[k]` for its own `k` and allocates its
  own `PadeStepperState`; no shared mutable buffer, no aliasing. Result
  ordering is by `k`, so downstream `argmin` is thread-count-independent
  (contingent on the caller's `f` being thread-safe, which is a
  documented contract, PathNetwork.jl:366-369 — not a code defect).

- **RNG seeding** (`PathNetwork.jl:480-481`): deterministically seeded
  `MersenneTwister(rng_seed)`, default seed 0; no global/unseeded RNG
  anywhere in `src/PathNetwork.jl`. A single solve is reproducible. (The
  highest-suspicion "unseeded RNG ⇒ run-to-run nondeterminism" hypothesis
  is **not** present in this module.)

- **Schwarz mirror sign** (`PathNetwork.jl:810, 837-839`): upper-canon
  representative is `complex(real(z), abs(imag(z)))`; for `imag(z) < 0`
  the output is `conj` of the upper value, i.e. `u(z) = conj(u(z̄))`,
  the correct Schwarz reflection of `u(z̄) = ū(z)`. IC-on-real-axis
  precondition checked and thrown (PathNetwork.jl:789-803).
