# Bug sweep E4 — mutable aliasing / in-place buffer reuse / RNG-and-order dependence

Auditor area: WHOLE `src/` through the **intermittency lens** — run-to-run or
order-dependent behavior that could surface as an intermittent discontinuity
in a computed solution.

Date: 2026-06-01. Read-only sweep (no `julia`, no `Pkg`, no tests run).

## Area

Every `*.jl` under `src/` (38 files, ~14.7k LOC) was scanned for the six
state-bug classes named in the assignment:

1. in-place bang functions that return / store a buffer later mutated;
2. `push!`/`append!` into a shared array that is also read;
3. default mutable kwargs that persist across calls;
4. global / unseeded RNG;
5. `sort!`/`shuffle!`/`reverse!` mutating caller-owned data;
6. closures capturing a loop variable by reference.

The headline result: **this codebase is unusually clean for this bug class.**
It is written in a near-purely-functional style — fresh `Vector{...}(undef,…)`
or array comprehensions per step, defensive `copy(state.y)` at the one place a
mutable stepper field is captured into a history, no in-place BLAS/LAPACK
(`mul!`/`ldiv!`/`lu!`/`qr!` appear **nowhere** in `src/`), no module-level
mutable globals, no `Ref`, no `global` reassignment, and the only RNG use is a
**seeded** `MersenneTwister(rng_seed)` with a deterministic default. The two
genuine order-dependence sites that exist are (a) explicitly documented as a
design limitation, and (b) confined to a diagnostic post-processing path that
cannot perturb the solution trajectory.

## References checked

- `docs/adr/0004-path-network-architecture.md:48-52` — Stage 1 spec:
  "shuffle targets, for each target find nearest visited node (Euclidean;
  **lexicographic tiebreak for reproducibility**)".
- `docs/adr/0004-path-network-architecture.md:51-52` — `max_steps_per_target
  = 1000`; exceed → `IntegrationError`.
- `docs/adr/0004-path-network-architecture.md:57` — Stage 2: "evaluate stored
  Padé at `t = (z_f - z_v)/h`. No new Taylor jets."
- `src/PathNetwork.jl:58-62` — module docstring explicitly documents the
  `shuffle(rng, targets)` → asymmetric-visited-tree consequence as a KNOWN,
  opt-out-able behavior (the `enforce_real_axis_symmetry` cure).
- `src/PathNetwork.jl:480-481` — `rng = MersenneTwister(rng_seed)` /
  `shuffle(rng, …)`; default `rng_seed = 0` (line 395) ⇒ deterministic.
- `src/PathNetwork.jl:882-894` — `_nearest_visited` lexicographic tiebreak.
- `src/PathNetwork.jl:935-955` — `_wedge_evaluations`, the `@threads` loop.
- `src/VectorPathNetwork.jl:714-723` — vector Stage-1 target ordering;
  `rng === nothing` default ⇒ no shuffle ⇒ deterministic.
- `src/VectorStepper.jl:172-173` — `VectorPadeStepperState` inner ctor copies
  `y` via `Vector{T}(y)`.
- `src/VectorStepper.jl:264-268` — `new_y = T[...]` fresh per step.
- `src/VectorProblems.jl:288` — `push!(y_vec, copy(state.y))` (defensive copy).
- `src/VectorPathNetwork.jl:828-837` — driver `push!(visited_y, y_new)` (no
  copy — analyzed below).
- `src/PoleField.jl:156` and `src/VectorPoleField.jl:320` — `sort!(candidates;
  by = c -> c[2])`.
- `src/SheetTracker.jl:283-292` — `winding_delta` principal-value normalization.
- `src/BranchTracker.jl:156-167` — `step_sheet_update` (fresh `collect(Int,…)`).
- `src/IVPBVPHybrid.jl:644-665` — problem-RHS closures (capture immutable
  scalars, not loop vars).

## Findings

### [low] `sort!(candidates; by = c -> c[2])` is unstable on equal keys — pole-representative position is tie-order-dependent

- **Location**: `src/PoleField.jl:156`, `src/VectorPoleField.jl:320`.
- **Ground truth (cited)**: `src/PoleField.jl:151-155` / `src/VectorPoleField.jl:303-312`
  document the intended semantics: greedy clustering in increasing `|t*|`
  (resp. `z_dist`), where "the first (best-placed) candidate to land in a
  cluster becomes its representative." `docs/adr/0004-path-network-architecture.md`
  does not specify a tiebreak for equal sort keys here (unlike
  `_nearest_visited`, which spec line 50-52 mandates a lexicographic tiebreak).
- **Code behavior**: Julia's default `sort!` with a non-`Base.Order.By`-on-
  primitive `by` key is **not a stable sort** (it dispatches to an introsort
  variant for general comparators). When two `candidates` entries have an
  exactly equal sort key `c[2]`, their relative post-sort order is
  implementation-defined. Because the *first* candidate to land becomes the
  cluster representative (`PoleField.jl:159-166`, `VectorPoleField.jl:324-337`),
  the reported pole's representative coordinate can depend on tie ordering.
- **Mechanism (intermittent)**: This is the *weakest* of the candidate
  intermittency vectors. Within one Julia build, `sort!` on a fixed input is
  deterministic run-to-run, so it does NOT produce true run-to-run flicker on
  one machine. It can produce a *different* reported pole list across Julia
  versions or if the upstream `candidates` build order ever changes. Crucially,
  this is a **diagnostic / pole-extraction** path (`extract_poles`), not the
  solution-evaluation path — it cannot move `grid_u` / `grid_y` values, so it
  cannot be the source of the maintainer's "intermittent discontinuities in
  *computed solutions*." Recorded for completeness; fix is `sort!(…; by, alg =
  Base.Sort.MergeSort)` plus a `(z_dist, real, imag)` tiebreak key to make the
  representative fully order-invariant.
- **Intermittent?**: No (deterministic within a Julia version); cross-version /
  cross-build only.
- **Confidence**: 0.55 that `sort!` is unstable on ties as described; 0.1 that
  this is the maintainer's bug.

### [low] Vector driver `push!(visited_y, y_new)` stores the value without the defensive `copy` its sibling `solve_vector_pade` uses

- **Location**: `src/VectorPathNetwork.jl:829` (`push!(visited_y, y_new)`),
  contrast `src/VectorProblems.jl:288` (`push!(y_vec, copy(state.y))`).
- **Ground truth (cited)**: No ADR mandates the copy; the relevant invariant is
  the general aliasing-safety rule. `src/VectorStepper.jl:264-268` shows the
  stepper builds `new_y = T[ … ]` fresh and assigns it to `state.y`, so the
  vector reaching the driver is a fresh allocation; `src/VectorWedgeStep.jl:407-409`
  stores that same `state.y` into `cand_y[k]` and returns `cand_y[best_k]` as
  `y_new`.
- **Code behavior**: `y_new` is the freshly-allocated `state.y` of a stepper
  state that is **local to and discarded by** `_select_wedge`. After the driver
  does `z_cur, y_cur = z_new, y_new` (line 837), `y_cur` and `visited_y[end]`
  alias the same vector — but neither is ever mutated element-wise afterward
  (the next step copies `y_cur` into a fresh `VectorPadeStepperState` via the
  `Vector{T}(y)` inner ctor at `VectorStepper.jl:172-173`). I grepped for any
  in-place write `visited_y[…][…] = …` / `y_cur[…] = …` / `.y[…] = …` across
  all `Vector*.jl` and found **none**.
- **Mechanism (intermittent)**: Latent only. As written, the missing `copy` is
  safe because the aliased vector is uniquely owned (the stepper state that
  produced it is unreachable) and is treated as immutable downstream. It becomes
  an aliasing bug only if a future edit introduces in-place mutation of `y_cur`
  or of a stored `visited_y[k]`. Flagged because the asymmetry with the
  `copy`-using sibling at `VectorProblems.jl:288` is a maintenance trap, not
  because it currently misbehaves.
- **Intermittent?**: No (not currently a defect).
- **Confidence**: 0.85 that no current aliasing defect exists here; ~0.1 that
  this is the maintainer's bug.

### [low] Threaded wedge evaluation (`@threads`) delegates thread-safety of the user RHS to the caller — a non-thread-safe RHS closure would produce intermittent results

- **Location**: `src/PathNetwork.jl:942-953` (`@threads for k in 1:n_w`).
- **Ground truth (cited)**: The module's own contract,
  `src/PathNetwork.jl:927-934`: "`evals[k]` is written by exactly one thread
  (its own `k`) … the user RHS must itself be thread-safe; closures over mutable
  shared state break this contract." `docs/adr/0004-path-network-architecture.md`
  does not mention threading (it predates the threaded wedge), so the
  thread-safety obligation lives only in the source docstring.
- **Code behavior**: Each wedge index `k` builds its own `PadeStepperState`
  (line 946) and writes only `evals[k]` (line 949 / 951). The result array is
  pre-allocated `undef` (line 941) and read back in `k`-order, so `argmin`
  downstream is deterministic **provided** every `f` call is pure. The internal
  arithmetic (`pade_step_with_pade!`, `_local_pade`) allocates fresh buffers and
  touches no shared state — verified clean.
- **Mechanism (intermittent)**: If a user passes an RHS `f` that closes over a
  mutable buffer (e.g. a preallocated scratch vector, a memoization `Dict`, or a
  counter) — a natural pattern users reach for — then the 5 wedge threads race
  on that buffer, and `evals[k]` becomes a function of thread interleaving.
  Different interleavings select a different wedge candidate, which forks the
  visited tree and yields **genuinely run-to-run-different solutions** — exactly
  an intermittent discontinuity. This is a *contract* exposure, not a defect in
  PadeTaylor's own code; every PadeTaylor-internal RHS (`Painleve`, `Heun`,
  `PainleveHierarchy`) is pure. The vector twin avoids the exposure entirely:
  `VectorWedgeStep._select_wedge` (`src/VectorWedgeStep.jl:401-413`) is a plain
  serial `for k in 1:n_w`, no `@threads`.
- **Intermittent?**: Yes — but only under a caller-supplied non-thread-safe RHS;
  not triggerable by any in-repo problem definition.
- **Confidence**: 0.9 that the threading is internally safe and the only
  exposure is the documented RHS contract; ~0.1 that this is the maintainer's
  bug (their problems are the pure in-repo RHSs).

## Areas verified correct

- **RNG is seeded and deterministic.** `src/PathNetwork.jl:480` constructs
  `MersenneTwister(rng_seed)` with default `rng_seed = 0`
  (`src/PathNetwork.jl:395`). No `rand`/`randn`/`shuffle` ever touches the
  global RNG. `src/VectorPathNetwork.jl:723` shuffles only when an explicit
  `rng` is passed (`rng === nothing` default ⇒ no shuffle), per its docstring
  `src/VectorPathNetwork.jl:714-723`. No run-to-run nondeterminism from RNG.

- **The `shuffle`-driven visited-tree asymmetry is a DOCUMENTED design
  limitation, not a hidden bug.** `src/PathNetwork.jl:58-62` states plainly that
  the default network is "not numerically Schwarz-symmetric" because the shuffle
  builds an asymmetric tree, and that the Stage-2 nearest-visited lookup can
  land conjugate query points on non-conjugate nodes — with the explicit
  `enforce_real_axis_symmetry = true` cure. This is order-*sensitive* but, for a
  fixed seed, fully *deterministic* and faithful to FW 2011 §3.1 line 156. It is
  not the intermittency the maintainer seeks.

- **`_nearest_visited` lexicographic tiebreak is correct and deterministic.**
  `src/PathNetwork.jl:882-894` and the identical
  `src/VectorPathNetwork.jl:921-…` compare `visited_z[i]` against the current
  best index `visited_z[idx]` on `(real, imag)` strictly, matching ADR-0004
  line 50-52 ("Euclidean; lexicographic tiebreak for reproducibility"). Ties are
  broken deterministically; no `argmin`/`findmin` (which would silently pick the
  first equal element) is used on the visited scan.

- **`_select_candidate` (scalar) and `_select_wedge` (vector) selection is
  deterministic.** `src/PathNetwork.jl:987-998` uses `argmin` over an
  index-ordered generator; `src/VectorWedgeStep.jl:432-443` uses a
  `(primary, -|wedge_angle|)` lexicographic tuple with strict `score >
  best_score` and `best_k == 0` sentinel — no float-tie ambiguity that flips
  with iteration order.

- **No in-place linear algebra.** `mul!`, `ldiv!`, `lu!`, `qr!`, `copyto!`,
  `fill!`, `broadcast!` appear **nowhere** in `src/`. The SVD path
  (`RobustPade`, routed through `BigFloat` per ADR-0002) and all polynomial
  evaluation (`_eval_poly` Horner loops, e.g. `src/VectorStepper.jl:305-311`,
  `src/VectorPathNetworkStage2.jl:453-459`) build fresh outputs. No
  buffer-reuse aliasing.

- **No module-level mutable state.** The only `const` containers are the two
  `DEFAULT_WEDGE` arrays (`src/PathNetwork.jl:153`,
  `src/VectorPathNetwork.jl:231`). Both are used solely as default kwarg values
  and are **never mutated** — no `push!`/`sort!`/element-assign/`reverse!`
  touches `wedge_angles` or `DEFAULT_WEDGE` anywhere (grepped). They are read
  via `length` and indexing only. So the classic "shared mutable default kwarg"
  bug does not occur.

- **`push!` accumulators are private to their builder.** Every `push!`/`append!`
  site (`Problems.jl:204-207`, `VectorProblems.jl:287-290`,
  `VectorPathNetwork.jl:828-834`, `Dispatcher.jl:289-380`,
  `LatticeDispatcher.jl:426`, `EdgeGatedSolve.jl:198-225`,
  `PoleField.jl:147-166`) appends to a vector that is local to the function and
  returned only after the loop finishes; none is read-during-grow in a way that
  changes results, and none is a caller-owned array being mutated.

- **`step_sheet_update` allocates fresh and is order-safe.**
  `src/BranchTracker.jl:159` `out = collect(Int, sheet_old)` — a fresh vector;
  the per-branch `out[k] += sign(Δθ)` (line 163) never aliases `sheet_old`. The
  `sign(0) → 0` increment is handled (`Δθ < 0 ? -1 : 0`). This is a
  `cross_branch = true` opt-in path only.

- **`winding_delta` principal-value normalization is the standard `(-π, π]`
  convention with a documented precondition.** `src/SheetTracker.jl:283-292`
  normalizes via `≤ -π` (adds 2π) and `> π` (subtracts 2π); the docstring
  (`src/SheetTracker.jl:275-281`) requires `|Δθ| < π` per step and proves it
  holds for `h = 0.5` and branch distance > 1 (`arcsin(0.5) ≈ 30°`), so the
  boundary is never reached in practice. No branch-cut intermittency in the
  default (no-branches) path, which pushes `Int[]` (`PathNetwork.jl:476-477`).

- **Closures capture immutable scalars, not loop variables.**
  `src/IVPBVPHybrid.jl:643-665` binds `α, β, γ, δ` once (scalars) and
  `asymptotic_ic_fn` once; the closures are not created inside a loop and
  capture nothing mutable. No "closure over loop variable by reference" hazard.

- **`_stage2_fill` per-node `radius` memo is keyed by node index and is a pure
  function of node data.** `src/VectorPathNetworkStage2.jl:408-428`: the
  `radius[idx_v]` memo caches a deterministic per-node gate; `grid_y[i]` is a
  freshly allocated comprehension or `fill(nan_C,…)` per grid point. Deterministic
  and aliasing-free. Matches ADR-0004 line 57 Stage-2 semantics
  (`t = (z_f - z_v)/h`, no new jets).

## Bottom line for this area

Through the aliasing / RNG / order-dependence lens specifically, I found **no
high- or medium-confidence defect** that would cause intermittent discontinuities
in a computed PadeTaylor solution. The three low-severity items above are either
diagnostic-only (`sort!` ties), latent-not-active (the missing `copy`), or a
documented caller contract (`@threads` RHS thread-safety). If the maintainer's
intermittency reproduces with a fixed `rng_seed` on a single machine with an
in-repo (pure) RHS, the cause is almost certainly NOT in this bug class — it is
more likely an algorithm/equation mis-transcription in the core
Padé/Taylor/step arithmetic or branch-cut/principal-value handling (the domains
of the A2/A3/A4 and B-series sweeps), surfacing as a data-dependent (not
run-to-run) discontinuity.
